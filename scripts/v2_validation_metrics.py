#!/usr/bin/env python3
"""Streaming biological validation metrics for the v2 RNA-seq workflow."""

from __future__ import annotations

from dataclasses import dataclass
import math
from pathlib import Path
from typing import Iterable

import numpy as np


@dataclass(frozen=True)
class GeneExons:
    gene_id: str
    chromosome: str
    strand: str
    intervals: tuple[tuple[int, int], ...]

    @property
    def start(self) -> int:
        return self.intervals[0][0]

    @property
    def end(self) -> int:
        return self.intervals[-1][1]


def parse_gtf_attributes(value: str) -> dict[str, str]:
    result: dict[str, str] = {}
    for part in value.strip().split(";"):
        part = part.strip()
        if not part:
            continue
        key, separator, raw = part.partition(" ")
        result[key] = raw.strip().strip('"') if separator else ""
    return result


def merge_intervals(intervals: Iterable[tuple[int, int]]) -> tuple[tuple[int, int], ...]:
    merged: list[tuple[int, int]] = []
    for start, end in sorted(intervals):
        if end <= start:
            continue
        if not merged or start > merged[-1][1]:
            merged.append((start, end))
        else:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
    return tuple(merged)


def load_gene_exons(
    path: str | Path, chromosomes: set[str] | None = None
) -> tuple[list[GeneExons], dict[str, list[int]]]:
    raw: dict[tuple[str, str, str], list[tuple[int, int]]] = {}
    with Path(path).open() as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "exon":
                continue
            chromosome, strand = fields[0], fields[6]
            if chromosomes is not None and chromosome not in chromosomes:
                continue
            attributes = parse_gtf_attributes(fields[8])
            gene_id = attributes.get("gene_id")
            if not gene_id:
                raise ValueError(f"Exon lacks gene_id: {line[:160].rstrip()}")
            raw.setdefault((chromosome, gene_id, strand), []).append(
                (int(fields[3]) - 1, int(fields[4]))
            )
    genes = [
        GeneExons(gene_id, chromosome, strand, merge_intervals(intervals))
        for (chromosome, gene_id, strand), intervals in raw.items()
    ]
    genes = [gene for gene in genes if gene.intervals]
    genes.sort(key=lambda gene: (gene.chromosome, gene.start, gene.end, gene.gene_id))
    by_chromosome: dict[str, list[int]] = {}
    for index, gene in enumerate(genes):
        by_chromosome.setdefault(gene.chromosome, []).append(index)
    return genes, by_chromosome


def interval_union_mask(
    intervals: Iterable[tuple[int, int]], start: int, end: int
) -> np.ndarray:
    mask = np.zeros(end - start, dtype=bool)
    for interval_start, interval_end in intervals:
        if interval_end <= start:
            continue
        if interval_start >= end:
            break
        local_start = max(start, interval_start) - start
        local_end = min(end, interval_end) - start
        if local_end > local_start:
            mask[local_start:local_end] = True
    return mask


class PerTrackStats:
    """Accumulate pointwise errors and Pearson sufficient statistics."""

    def __init__(self, n_tracks: int) -> None:
        self.n = np.zeros(n_tracks, dtype=np.int64)
        self.sse = np.zeros(n_tracks, dtype=np.float64)
        self.sae = np.zeros(n_tracks, dtype=np.float64)
        self.sum_x = np.zeros(n_tracks, dtype=np.float64)
        self.sum_y = np.zeros(n_tracks, dtype=np.float64)
        self.sum_xx = np.zeros(n_tracks, dtype=np.float64)
        self.sum_yy = np.zeros(n_tracks, dtype=np.float64)
        self.sum_xy = np.zeros(n_tracks, dtype=np.float64)

    def update(
        self,
        prediction: np.ndarray,
        target: np.ndarray,
        mask: np.ndarray | None = None,
    ) -> None:
        prediction = np.asarray(prediction)
        target = np.asarray(target)
        if prediction.shape != target.shape or prediction.ndim != 2:
            raise ValueError("prediction and target must have equal [tracks, positions] shape")
        if prediction.shape[0] != len(self.n):
            raise ValueError("track count mismatch")
        if mask is None and np.isfinite(prediction).all() and np.isfinite(target).all():
            error = prediction - target
            self.n += prediction.shape[1]
            self.sse += np.square(error).sum(axis=1, dtype=np.float64)
            self.sae += np.abs(error).sum(axis=1, dtype=np.float64)
            self.sum_x += prediction.sum(axis=1, dtype=np.float64)
            self.sum_y += target.sum(axis=1, dtype=np.float64)
            self.sum_xx += np.square(prediction).sum(axis=1, dtype=np.float64)
            self.sum_yy += np.square(target).sum(axis=1, dtype=np.float64)
            self.sum_xy += (prediction * target).sum(axis=1, dtype=np.float64)
            return
        valid = (
            np.ones(prediction.shape, dtype=bool)
            if mask is None
            else np.asarray(mask, dtype=bool)
        )
        if valid.ndim == 1:
            valid = np.broadcast_to(valid, prediction.shape)
        elif valid.shape != prediction.shape:
            raise ValueError("mask must have [positions] or [tracks, positions] shape")
        valid = valid & np.isfinite(prediction) & np.isfinite(target)
        for track in range(prediction.shape[0]):
            x = prediction[track, valid[track]]
            y = target[track, valid[track]]
            if not len(x):
                continue
            error = x - y
            self.n[track] += len(x)
            self.sse[track] += np.square(error).sum(dtype=np.float64)
            self.sae[track] += np.abs(error).sum(dtype=np.float64)
            self.sum_x[track] += x.sum(dtype=np.float64)
            self.sum_y[track] += y.sum(dtype=np.float64)
            self.sum_xx[track] += np.square(x).sum(dtype=np.float64)
            self.sum_yy[track] += np.square(y).sum(dtype=np.float64)
            self.sum_xy[track] += (x * y).sum(dtype=np.float64)

    def per_track(self) -> dict[str, np.ndarray]:
        n = self.n.astype(np.float64)
        mse = np.divide(self.sse, n, out=np.full_like(n, np.nan), where=n > 0)
        mae = np.divide(self.sae, n, out=np.full_like(n, np.nan), where=n > 0)
        numerator = n * self.sum_xy - self.sum_x * self.sum_y
        denominator = np.sqrt(
            np.maximum(n * self.sum_xx - np.square(self.sum_x), 0.0)
            * np.maximum(n * self.sum_yy - np.square(self.sum_y), 0.0)
        )
        pearson = np.divide(
            numerator,
            denominator,
            out=np.full_like(numerator, np.nan),
            where=denominator > 0,
        )
        return {"n": self.n.copy(), "mse": mse, "mae": mae, "pearson": pearson}


class GeneExonAccumulator:
    """Accumulate per-gene exon mean coverage without double-counting transcripts."""

    def __init__(
        self,
        genes: list[GeneExons],
        by_chromosome: dict[str, list[int]],
        n_tracks: int,
    ) -> None:
        self.genes = genes
        self.by_chromosome = by_chromosome
        self.prediction_sum = np.zeros((len(genes), n_tracks), dtype=np.float64)
        self.target_sum = np.zeros((len(genes), n_tracks), dtype=np.float64)
        self.base_count = np.zeros(len(genes), dtype=np.int64)

    def update(
        self,
        chromosome: str,
        start: int,
        prediction: np.ndarray,
        target: np.ndarray,
    ) -> None:
        prediction = np.asarray(prediction)
        target = np.asarray(target)
        if prediction.shape != target.shape or prediction.ndim != 2:
            raise ValueError("prediction and target must have equal [tracks, positions] shape")
        end = start + prediction.shape[1]
        for gene_index in self.by_chromosome.get(chromosome, []):
            gene = self.genes[gene_index]
            if gene.end <= start:
                continue
            if gene.start >= end:
                break
            for exon_start, exon_end in gene.intervals:
                local_start = max(exon_start, start) - start
                local_end = min(exon_end, end) - start
                if local_end <= local_start:
                    continue
                self.prediction_sum[gene_index] += prediction[:, local_start:local_end].sum(
                    axis=1, dtype=np.float64
                )
                self.target_sum[gene_index] += target[:, local_start:local_end].sum(
                    axis=1, dtype=np.float64
                )
                self.base_count[gene_index] += local_end - local_start

    def metrics(self, *, log1p: bool = True) -> dict[str, object]:
        valid_genes = self.base_count > 0
        if valid_genes.sum() < 2:
            raise RuntimeError("Fewer than two genes have evaluated exon bases")
        denominator = self.base_count[valid_genes, None]
        prediction = self.prediction_sum[valid_genes] / denominator
        target = self.target_sum[valid_genes] / denominator
        if log1p:
            prediction = np.log1p(np.maximum(prediction, 0.0))
            target = np.log1p(np.maximum(target, 0.0))
        stats = PerTrackStats(prediction.shape[1])
        stats.update(prediction.T, target.T)
        per_track = stats.per_track()
        finite = np.isfinite(per_track["pearson"])
        return {
            "transform": "log1p_mean_exon_coverage" if log1p else "mean_exon_coverage",
            "genes_evaluated": int(valid_genes.sum()),
            "exon_bases_evaluated": int(self.base_count.sum()),
            "mean_per_track_pearson": float(per_track["pearson"][finite].mean()),
            "finite_per_track_pearson": int(finite.sum()),
            "per_track_pearson": per_track["pearson"],
            "per_track_mse": per_track["mse"],
        }


def average_ranks(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values)
    order = np.argsort(values, kind="mergesort")
    sorted_values = values[order]
    ranks = np.empty(len(values), dtype=np.float64)
    index = 0
    while index < len(values):
        stop = index + 1
        while stop < len(values) and sorted_values[stop] == sorted_values[index]:
            stop += 1
        ranks[order[index:stop]] = (index + stop - 1) / 2.0
        index = stop
    return ranks


def pearson_arrays(x: np.ndarray, y: np.ndarray) -> float:
    x = np.asarray(x, dtype=np.float64).reshape(-1)
    y = np.asarray(y, dtype=np.float64).reshape(-1)
    valid = np.isfinite(x) & np.isfinite(y)
    x, y = x[valid], y[valid]
    if len(x) < 2:
        return math.nan
    x = x - x.mean()
    y = y - y.mean()
    denominator = math.sqrt(float(np.square(x).sum() * np.square(y).sum()))
    return float((x * y).sum() / denominator) if denominator > 0 else math.nan


def spearman_arrays(x: np.ndarray, y: np.ndarray) -> float:
    return pearson_arrays(average_ranks(x), average_ranks(y))


def distribution_metrics(
    prediction: np.ndarray,
    target: np.ndarray,
    *,
    max_spearman_points: int = 100_000,
    seed: int = 20260714,
) -> dict[str, np.ndarray | int | str]:
    prediction = np.asarray(prediction, dtype=np.float64)
    target = np.asarray(target, dtype=np.float64)
    if prediction.shape != target.shape or prediction.ndim != 2:
        raise ValueError("prediction and target must have equal [tracks, positions] shape")
    n_tracks, n_positions = prediction.shape
    rng = np.random.default_rng(seed)
    if n_positions > max_spearman_points:
        indices = np.sort(rng.choice(n_positions, max_spearman_points, replace=False))
    else:
        indices = np.arange(n_positions)
    spearman = np.array(
        [spearman_arrays(prediction[i, indices], target[i, indices]) for i in range(n_tracks)]
    )
    top1_mse = np.full(n_tracks, np.nan, dtype=np.float64)
    top1_calibration_ratio = np.full(n_tracks, np.nan, dtype=np.float64)
    top1_threshold = np.full(n_tracks, np.nan, dtype=np.float64)
    for track in range(n_tracks):
        threshold = float(np.quantile(target[track], 0.99))
        mask = target[track] >= threshold
        top1_threshold[track] = threshold
        if not mask.any():
            continue
        error = prediction[track, mask] - target[track, mask]
        top1_mse[track] = float(np.square(error).mean())
        target_mean = float(target[track, mask].mean())
        if target_mean > 0:
            top1_calibration_ratio[track] = float(prediction[track, mask].mean()) / target_mean
    return {
        "spearman_sampling": "deterministic_uniform_without_replacement_across_128bp_bins",
        "spearman_points_per_track": int(len(indices)),
        "per_track_spearman_128bp": spearman,
        "per_track_top1_mse_128bp": top1_mse,
        "per_track_top1_calibration_ratio_128bp": top1_calibration_ratio,
        "per_track_top1_target_threshold_128bp": top1_threshold,
    }


def finite_mean(values: np.ndarray) -> float:
    values = np.asarray(values, dtype=np.float64)
    finite = values[np.isfinite(values)]
    if not finite.size:
        raise RuntimeError("Metric has no finite values")
    return float(finite.mean())


def keyed(values: np.ndarray, keys: list[str]) -> dict[str, float | None]:
    if len(values) != len(keys):
        raise ValueError("value/key count mismatch")
    return {
        key: float(value) if math.isfinite(float(value)) else None
        for key, value in zip(keys, values, strict=True)
    }
