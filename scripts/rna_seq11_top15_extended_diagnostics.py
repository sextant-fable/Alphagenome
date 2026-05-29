#!/usr/bin/env python3
"""Extended valid-split diagnostics for the top RNA-seq11 checkpoints.

This script intentionally does not read the test split.  It ranks existing
run directories by valid full-resolution log1p MSE, uses Pearson as a
tie-breaker, then evaluates the selected checkpoints on the valid dataset.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
import torch
import torch.nn.functional as F

from alphagenome_pytorch import AlphaGenome
from alphagenome_rna_seq11_adapter import (
    RnaSeq11Adapter,
    load_adapter_checkpoint,
    parse_resolutions,
    pearson_from_sums,
    pick_device,
    rankdata_average,
)
from alphagenome_rna_seq11_eval import load_track_names
from torch_rna_seq_dataset import make_dataloader


FLOAT_FORMAT = ".8g"


@dataclass(frozen=True)
class Candidate:
    rank: int
    run_id: str
    run_dir: Path
    checkpoint: Path
    source: str
    best_step: int | None
    full_mse: float
    full_mae: float | None
    pearson: float
    common128_mse: float | None
    common128_mae: float | None


@dataclass(frozen=True)
class Gene:
    gene_id: str
    gene_name: str
    chrom: str
    start: int
    end: int
    strand: str
    gene_biotype: str


@dataclass(frozen=True)
class Interval:
    chrom: str
    start: int
    end: int


@dataclass(frozen=True)
class GeneOverlap:
    gene_index: int
    local_start: int
    local_end: int


class SumStats:
    """Accumulate per-track pointwise sums."""

    def __init__(self, n_tracks: int) -> None:
        self.sse = np.zeros(n_tracks, dtype=np.float64)
        self.sae = np.zeros(n_tracks, dtype=np.float64)
        self.sum_x = np.zeros(n_tracks, dtype=np.float64)
        self.sum_y = np.zeros(n_tracks, dtype=np.float64)
        self.sum_x2 = np.zeros(n_tracks, dtype=np.float64)
        self.sum_y2 = np.zeros(n_tracks, dtype=np.float64)
        self.sum_xy = np.zeros(n_tracks, dtype=np.float64)
        self.n = np.zeros(n_tracks, dtype=np.float64)

    def update_torch(
        self,
        prediction: torch.Tensor,
        target: torch.Tensor,
        mask: torch.Tensor | None = None,
    ) -> None:
        """Update from tensors shaped [B, C, S] or [C, S]."""

        if prediction.ndim == 2:
            prediction = prediction.unsqueeze(0)
            target = target.unsqueeze(0)
        if mask is None:
            mask = torch.ones(
                prediction.shape[:2] + (1,),
                dtype=torch.bool,
                device=prediction.device,
            )
        while mask.ndim < prediction.ndim:
            mask = mask.unsqueeze(-1)
        mask = mask.to(device=prediction.device, dtype=torch.bool)
        if mask.shape[-1] == 1:
            mask = mask.expand_as(prediction)

        prediction64 = prediction.to(dtype=torch.float64)
        target64 = target.to(dtype=torch.float64)
        mask64 = mask.to(dtype=torch.float64)
        error64 = prediction64 - target64

        self.sse += (
            (error64.square() * mask64).sum(dim=(0, 2)).detach().cpu().numpy()
        )
        self.sae += (
            (error64.abs() * mask64).sum(dim=(0, 2)).detach().cpu().numpy()
        )
        self.sum_x += (
            (prediction64 * mask64).sum(dim=(0, 2)).detach().cpu().numpy()
        )
        self.sum_y += (target64 * mask64).sum(dim=(0, 2)).detach().cpu().numpy()
        self.sum_x2 += (
            (prediction64.square() * mask64).sum(dim=(0, 2)).detach().cpu().numpy()
        )
        self.sum_y2 += (
            (target64.square() * mask64).sum(dim=(0, 2)).detach().cpu().numpy()
        )
        self.sum_xy += (
            (prediction64 * target64 * mask64).sum(dim=(0, 2)).detach().cpu().numpy()
        )
        self.n += mask64.sum(dim=(0, 2)).detach().cpu().numpy()

    def update_numpy(self, prediction: np.ndarray, target: np.ndarray) -> None:
        """Update from arrays shaped [C, N]."""

        if prediction.size == 0:
            return
        error = prediction.astype(np.float64, copy=False) - target.astype(
            np.float64,
            copy=False,
        )
        prediction64 = prediction.astype(np.float64, copy=False)
        target64 = target.astype(np.float64, copy=False)
        self.sse += np.square(error).sum(axis=1)
        self.sae += np.abs(error).sum(axis=1)
        self.sum_x += prediction64.sum(axis=1)
        self.sum_y += target64.sum(axis=1)
        self.sum_x2 += np.square(prediction64).sum(axis=1)
        self.sum_y2 += np.square(target64).sum(axis=1)
        self.sum_xy += (prediction64 * target64).sum(axis=1)
        self.n += prediction.shape[1]

    def metric_for_track(self, track_index: int) -> dict[str, float]:
        n = float(self.n[track_index])
        if n <= 0:
            return {"mse": math.nan, "mae": math.nan, "pearson": math.nan, "n": 0.0}
        return {
            "mse": float(self.sse[track_index] / n),
            "mae": float(self.sae[track_index] / n),
            "pearson": pearson_from_sums(
                float(self.sum_x[track_index]),
                float(self.sum_y[track_index]),
                float(self.sum_x2[track_index]),
                float(self.sum_y2[track_index]),
                float(self.sum_xy[track_index]),
                n,
            ),
            "n": n,
        }

    def overall_metric(self) -> dict[str, float]:
        n = float(self.n.sum())
        if n <= 0:
            return {"mse": math.nan, "mae": math.nan, "pearson": math.nan, "n": 0.0}
        return {
            "mse": float(self.sse.sum() / n),
            "mae": float(self.sae.sum() / n),
            "pearson": pearson_from_sums(
                float(self.sum_x.sum()),
                float(self.sum_y.sum()),
                float(self.sum_x2.sum()),
                float(self.sum_y2.sum()),
                float(self.sum_xy.sum()),
                n,
            ),
            "n": n,
        }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--runs-dir", default="runs")
    parser.add_argument(
        "--dataset-dir",
        default="alphagenome_custom/datasets/rna_seq_npz_valid",
        help="Valid NPZ dataset. Do not point this script at test.",
    )
    parser.add_argument(
        "--gtf",
        default="alphagenome_custom/reference/Caenorhabditis_elegans.WBcel235.115.gtf",
    )
    parser.add_argument(
        "--weights",
        default="weights/alphagenome_pytorch/model_all_folds.safetensors",
    )
    parser.add_argument(
        "--output-dir",
        default="runs/rna_seq11_top15_extended_diagnostics_20260529",
    )
    parser.add_argument("--top-n", type=int, default=15)
    parser.add_argument("--batch-size", type=int, default=1)
    parser.add_argument("--num-workers", type=int, default=0)
    parser.add_argument("--max-examples", type=int, default=None)
    parser.add_argument("--device", default="auto")
    parser.add_argument(
        "--rank-shard",
        type=int,
        default=0,
        help="Evaluate selected ranks where zero-based index %% num-shards equals this value.",
    )
    parser.add_argument("--num-shards", type=int, default=1)
    parser.add_argument(
        "--only-write-top15",
        action="store_true",
        help="Collect/rank models and write top15_models.tsv without loading AlphaGenome.",
    )
    parser.add_argument(
        "--spearman-sample-size",
        type=int,
        default=300000,
        help="Approximate full-resolution Spearman sample count per track.",
    )
    parser.add_argument("--spearman-seed", type=int, default=20260529)
    parser.add_argument("--promoter-radius-bp", type=int, default=1000)
    parser.add_argument(
        "--pooled-top-factors",
        default="128",
        help="Comma-separated pooling factors used for high-signal localization.",
    )
    return parser.parse_args()


def as_float(value: str | None) -> float | None:
    if value is None or value == "":
        return None
    try:
        return float(value)
    except ValueError:
        return None


def as_int(value: str | None) -> int | None:
    if value is None or value == "":
        return None
    try:
        return int(float(value))
    except ValueError:
        return None


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def select_run_metric(run_dir: Path) -> tuple[str, int | None, float, float | None, float, float | None, float | None] | None:
    metrics_path = run_dir / "metrics.tsv"
    if metrics_path.exists():
        rows = read_tsv(metrics_path)
        valid_rows = [
            row
            for row in rows
            if row.get("split") == "valid"
            and as_float(row.get("full_mse")) is not None
            and as_float(row.get("pearson")) is not None
        ]
        if valid_rows:
            best = sorted(
                valid_rows,
                key=lambda row: (
                    float(row["full_mse"]),
                    -float(row["pearson"]),
                    as_int(row.get("step")) or 0,
                ),
            )[0]
            return (
                "metrics.tsv",
                as_int(best.get("step")),
                float(best["full_mse"]),
                as_float(best.get("full_mae")),
                float(best["pearson"]),
                as_float(best.get("common128_mse")),
                as_float(best.get("common128_mae")),
            )

    point_metrics_path = run_dir / "valid_best_pointwise_metrics.tsv"
    if point_metrics_path.exists():
        for row in read_tsv(point_metrics_path):
            if row.get("metric_scope") != "overall":
                continue
            mse = as_float(row.get("mse"))
            pearson = as_float(row.get("pearson"))
            if mse is None or pearson is None:
                continue
            return (
                "valid_best_pointwise_metrics.tsv",
                None,
                mse,
                as_float(row.get("mae")),
                pearson,
                None,
                None,
            )
    return None


def collect_top_candidates(runs_dir: Path, top_n: int) -> list[Candidate]:
    skip_pattern = re.compile(r"(smoke|diagnostic|diagnostics|code_smoke)", re.IGNORECASE)
    candidates: list[Candidate] = []
    for run_dir in sorted(path for path in runs_dir.iterdir() if path.is_dir()):
        if skip_pattern.search(run_dir.name):
            continue
        checkpoint = run_dir / "adapter_head_best.pt"
        if not checkpoint.exists():
            continue
        metric = select_run_metric(run_dir)
        if metric is None:
            continue
        source, step, full_mse, full_mae, pearson, common128_mse, common128_mae = metric
        candidates.append(
            Candidate(
                rank=0,
                run_id=run_dir.name,
                run_dir=run_dir,
                checkpoint=checkpoint,
                source=source,
                best_step=step,
                full_mse=full_mse,
                full_mae=full_mae,
                pearson=pearson,
                common128_mse=common128_mse,
                common128_mae=common128_mae,
            )
        )
    candidates.sort(key=lambda item: (item.full_mse, -item.pearson, item.run_id))
    return [
        Candidate(
            rank=index + 1,
            run_id=item.run_id,
            run_dir=item.run_dir,
            checkpoint=item.checkpoint,
            source=item.source,
            best_step=item.best_step,
            full_mse=item.full_mse,
            full_mae=item.full_mae,
            pearson=item.pearson,
            common128_mse=item.common128_mse,
            common128_mae=item.common128_mae,
        )
        for index, item in enumerate(candidates[:top_n])
    ]


def write_top_models(path: Path, candidates: list[Candidate]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        fieldnames = [
            "rank",
            "run_id",
            "checkpoint",
            "source",
            "best_step",
            "full_mse",
            "full_mae",
            "pearson",
            "common128_mse",
            "common128_mae",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for candidate in candidates:
            writer.writerow(
                {
                    "rank": candidate.rank,
                    "run_id": candidate.run_id,
                    "checkpoint": str(candidate.checkpoint),
                    "source": candidate.source,
                    "best_step": "" if candidate.best_step is None else candidate.best_step,
                    "full_mse": format_float(candidate.full_mse),
                    "full_mae": format_optional_float(candidate.full_mae),
                    "pearson": format_float(candidate.pearson),
                    "common128_mse": format_optional_float(candidate.common128_mse),
                    "common128_mae": format_optional_float(candidate.common128_mae),
                }
            )


def format_float(value: float) -> str:
    if math.isnan(value):
        return "nan"
    return format(float(value), FLOAT_FORMAT)


def format_optional_float(value: float | None) -> str:
    return "" if value is None else format_float(value)


def read_manifest(dataset_dir: Path, max_examples: int | None = None) -> list[dict[str, str]]:
    rows = read_tsv(dataset_dir / "manifest.tsv")
    if max_examples is not None:
        rows = rows[:max_examples]
    return rows


def parse_attributes(value: str) -> dict[str, str]:
    attributes: dict[str, str] = {}
    for part in value.strip().split(";"):
        part = part.strip()
        if not part:
            continue
        if " " not in part:
            attributes[part] = ""
            continue
        key, raw_value = part.split(" ", maxsplit=1)
        attributes[key] = raw_value.strip().strip('"')
    return attributes


def load_gtf_annotations(
    gtf_path: Path,
    *,
    chroms: set[str],
    promoter_radius_bp: int,
    chrom_lengths: dict[str, int],
) -> tuple[list[Gene], dict[str, list[Interval]], dict[str, list[Interval]], dict[str, list[Interval]]]:
    genes: list[Gene] = []
    gene_intervals: dict[str, list[Interval]] = {chrom: [] for chrom in chroms}
    exon_intervals: dict[str, list[Interval]] = {chrom: [] for chrom in chroms}
    promoter_intervals: dict[str, list[Interval]] = {chrom: [] for chrom in chroms}

    with gtf_path.open() as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            chrom, _source, feature, start_s, end_s, _score, strand, _frame, attr_s = fields
            if chrom not in chroms:
                continue
            start = int(start_s) - 1
            end = int(end_s)
            if end <= start:
                continue
            attrs = parse_attributes(attr_s)
            if feature == "gene":
                gene_id = attrs.get("gene_id", f"{chrom}:{start}-{end}")
                gene_name = attrs.get("gene_name", gene_id)
                gene_biotype = attrs.get("gene_biotype", "")
                genes.append(
                    Gene(
                        gene_id=gene_id,
                        gene_name=gene_name,
                        chrom=chrom,
                        start=start,
                        end=end,
                        strand=strand,
                        gene_biotype=gene_biotype,
                    )
                )
                gene_intervals[chrom].append(Interval(chrom=chrom, start=start, end=end))
            elif feature == "exon":
                exon_intervals[chrom].append(Interval(chrom=chrom, start=start, end=end))
            elif feature == "transcript":
                tss = start if strand != "-" else end
                chrom_length = chrom_lengths.get(chrom, end)
                promoter_start = max(0, tss - promoter_radius_bp)
                promoter_end = min(chrom_length, tss + promoter_radius_bp)
                if promoter_end > promoter_start:
                    promoter_intervals[chrom].append(
                        Interval(chrom=chrom, start=promoter_start, end=promoter_end)
                    )

    for interval_lists in (gene_intervals, exon_intervals, promoter_intervals):
        for chrom in interval_lists:
            interval_lists[chrom] = merge_intervals(interval_lists[chrom])
    genes.sort(key=lambda gene: (gene.chrom, gene.start, gene.end, gene.gene_id))
    return genes, gene_intervals, exon_intervals, promoter_intervals


def merge_intervals(intervals: Iterable[Interval]) -> list[Interval]:
    sorted_intervals = sorted(intervals, key=lambda item: (item.chrom, item.start, item.end))
    merged: list[Interval] = []
    for interval in sorted_intervals:
        if not merged or interval.chrom != merged[-1].chrom or interval.start > merged[-1].end:
            merged.append(interval)
            continue
        previous = merged[-1]
        merged[-1] = Interval(
            chrom=previous.chrom,
            start=previous.start,
            end=max(previous.end, interval.end),
        )
    return merged


def overlapping_intervals(
    intervals: list[Interval],
    *,
    chrom: str,
    start: int,
    end: int,
) -> Iterable[Interval]:
    for interval in intervals:
        if interval.chrom != chrom:
            continue
        if interval.end <= start:
            continue
        if interval.start >= end:
            break
        yield interval


def make_interval_mask(
    intervals: list[Interval],
    *,
    chrom: str,
    window_start: int,
    window_end: int,
    local_start: int = 0,
    local_end: int | None = None,
) -> np.ndarray:
    if local_end is None:
        local_end = window_end - window_start
    length = local_end - local_start
    mask = np.zeros(length, dtype=bool)
    region_start = window_start + local_start
    region_end = window_start + local_end
    for interval in overlapping_intervals(
        intervals,
        chrom=chrom,
        start=region_start,
        end=region_end,
    ):
        start = max(interval.start, region_start) - region_start
        end = min(interval.end, region_end) - region_start
        if end > start:
            mask[start:end] = True
    return mask


def compute_core_ends(manifest_rows: list[dict[str, str]]) -> list[int]:
    core_ends: list[int] = []
    for index, row in enumerate(manifest_rows):
        start = int(row["start"])
        end = int(row["end"])
        next_start = (
            int(manifest_rows[index + 1]["start"])
            if index + 1 < len(manifest_rows)
            and manifest_rows[index + 1]["chromosome"] == row["chromosome"]
            else end
        )
        core_ends.append(min(end - start, max(0, next_start - start)))
    return core_ends


def build_gene_overlaps(
    genes: list[Gene],
    manifest_rows: list[dict[str, str]],
    core_ends: list[int],
) -> list[list[GeneOverlap]]:
    overlaps_by_window: list[list[GeneOverlap]] = []
    for row, core_end in zip(manifest_rows, core_ends):
        chrom = row["chromosome"]
        start = int(row["start"])
        region_start = start
        region_end = start + core_end
        overlaps: list[GeneOverlap] = []
        for gene_index, gene in enumerate(genes):
            if gene.chrom != chrom:
                continue
            if gene.end <= region_start:
                continue
            if gene.start >= region_end:
                break
            overlap_start = max(gene.start, region_start)
            overlap_end = min(gene.end, region_end)
            if overlap_end > overlap_start:
                overlaps.append(
                    GeneOverlap(
                        gene_index=gene_index,
                        local_start=overlap_start - start,
                        local_end=overlap_end - start,
                    )
                )
        overlaps_by_window.append(overlaps)
    return overlaps_by_window


def compute_target_quantiles(
    dataset_dir: Path,
    *,
    n_tracks: int,
    quantiles: tuple[float, ...],
    max_examples: int | None,
) -> dict[float, np.ndarray]:
    manifest_rows = read_manifest(dataset_dir, max_examples=max_examples)
    thresholds = {quantile: np.zeros(n_tracks, dtype=np.float32) for quantile in quantiles}
    for track_index in range(n_tracks):
        chunks: list[np.ndarray] = []
        for row in manifest_rows:
            with np.load(dataset_dir / row["path"], allow_pickle=False) as data:
                target = data["rna_seq"][:, track_index].astype(np.float32, copy=False)
                target = np.log1p(np.clip(target, a_min=0.0, a_max=None))
                chunks.append(np.asarray(target, dtype=np.float32))
        values = np.concatenate(chunks)
        for quantile in quantiles:
            thresholds[quantile][track_index] = float(np.quantile(values, quantile))
    return thresholds


def make_model_from_checkpoint(
    *,
    checkpoint_path: Path,
    base_model: AlphaGenome,
    n_tracks: int,
    device: torch.device,
) -> tuple[RnaSeq11Adapter, dict[str, object]]:
    checkpoint = load_adapter_checkpoint(checkpoint_path)
    embedding_resolution = int(checkpoint["embedding_resolution"])
    checkpoint_head_resolutions = tuple(
        int(resolution)
        for resolution in checkpoint.get("head_resolutions", [embedding_resolution])
    )
    head_type = str(checkpoint.get("head_type", "linear"))
    head_resolutions = parse_resolutions(
        checkpoint_head_resolutions,
        fallback_resolution=embedding_resolution,
    )
    linear_head_architecture = str(checkpoint.get("linear_head_architecture", "conv1x1"))
    linear_hidden_channels = int(checkpoint.get("linear_hidden_channels", 256))
    linear_input_bottleneck_channels = checkpoint.get("linear_input_bottleneck_channels")
    residual_scale_init = float(checkpoint.get("linear_residual_scale_init", 1.0))
    linear_dilation = int(checkpoint.get("linear_dilation", 2))
    linear_kernel_size = int(checkpoint.get("linear_kernel_size", 15))
    organism_index = int(checkpoint["organism_index"])
    residual_base_checkpoint_path = checkpoint.get("residual_base_checkpoint")
    linear_residual_fusion = str(checkpoint.get("linear_residual_fusion", "none"))
    residual_correction_scale_init = float(
        checkpoint.get("residual_correction_scale_init", 0.01)
    )
    residual_base_gate = str(checkpoint.get("residual_base_gate", "none"))
    residual_base_gate_center = float(checkpoint.get("residual_base_gate_center", 0.5))
    residual_base_gate_sharpness = float(
        checkpoint.get("residual_base_gate_sharpness", 4.0)
    )
    residual_base_gate_floor = float(checkpoint.get("residual_base_gate_floor", 0.0))
    checkpoint_head_state = checkpoint["adapter_head_state_dict"]
    track_means = checkpoint_head_state.get("track_means")

    if head_type != "linear":
        raise NotImplementedError("Extended diagnostics currently support linear heads")

    residual_base_head = None
    residual_base_input_bottleneck = None
    residual_base_resolution = 128
    if residual_base_checkpoint_path is not None:
        residual_base_checkpoint = load_adapter_checkpoint(str(residual_base_checkpoint_path))
        residual_base_resolution = int(residual_base_checkpoint["embedding_resolution"])
        residual_base_probe = RnaSeq11Adapter(
            base_model,
            n_tracks=n_tracks,
            embedding_resolution=residual_base_resolution,
            organism_index=organism_index,
            head_type="linear",
            head_resolutions=(residual_base_resolution,),
            linear_head_architecture=str(
                residual_base_checkpoint.get("linear_head_architecture", "conv1x1")
            ),
            linear_hidden_channels=int(
                residual_base_checkpoint.get("linear_hidden_channels", 256)
            ),
            linear_input_bottleneck_channels=residual_base_checkpoint.get(
                "linear_input_bottleneck_channels"
            ),
            linear_residual_scale_init=float(
                residual_base_checkpoint.get("linear_residual_scale_init", 1.0)
            ),
            linear_dilation=int(residual_base_checkpoint.get("linear_dilation", 2)),
            linear_kernel_size=int(residual_base_checkpoint.get("linear_kernel_size", 15)),
        ).to(device)
        residual_base_probe.load_adapter_head_state_dict(
            residual_base_checkpoint["adapter_head_state_dict"]
        )
        residual_base_probe.eval()
        residual_base_head = residual_base_probe.head
        residual_base_input_bottleneck = residual_base_probe.linear_input_bottleneck

    model = RnaSeq11Adapter(
        base_model,
        n_tracks=n_tracks,
        embedding_resolution=embedding_resolution,
        organism_index=organism_index,
        head_type=head_type,
        head_resolutions=head_resolutions,
        linear_head_architecture=linear_head_architecture,
        linear_hidden_channels=linear_hidden_channels,
        linear_input_bottleneck_channels=linear_input_bottleneck_channels,
        linear_track_means=checkpoint_head_state.get("track_means"),
        linear_residual_scale_init=residual_scale_init,
        linear_dilation=linear_dilation,
        linear_kernel_size=linear_kernel_size,
        linear_residual_base_head=residual_base_head,
        linear_residual_base_input_bottleneck=residual_base_input_bottleneck,
        linear_residual_base_resolution=residual_base_resolution,
        linear_residual_fusion=linear_residual_fusion,
        linear_residual_correction_scale_init=residual_correction_scale_init,
        linear_residual_base_gate=residual_base_gate,
        linear_residual_base_gate_center=residual_base_gate_center,
        linear_residual_base_gate_sharpness=residual_base_gate_sharpness,
        linear_residual_base_gate_floor=residual_base_gate_floor,
        track_means=track_means,
    ).to(device)
    model.load_adapter_head_state_dict(checkpoint_head_state)
    model.eval()
    return model, checkpoint


def metric_rows_for_stats(
    *,
    candidate: Candidate,
    scope: str,
    stats: SumStats,
    track_labels: list[str],
    extra: dict[str, object] | None = None,
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    extra = {} if extra is None else extra
    overall = stats.overall_metric()
    rows.append(
        {
            "rank": candidate.rank,
            "run_id": candidate.run_id,
            "metric_scope": scope,
            "track_index": "",
            "track_id": "",
            "track_name": "",
            **extra,
            "mse": format_float(overall["mse"]),
            "mae": format_float(overall["mae"]),
            "pearson": format_float(overall["pearson"]),
            "n_values": format_float(overall["n"]),
        }
    )
    for track_index, label in enumerate(track_labels):
        track_id, track_name = split_track_label(label)
        metrics = stats.metric_for_track(track_index)
        rows.append(
            {
                "rank": candidate.rank,
                "run_id": candidate.run_id,
                "metric_scope": scope,
                "track_index": track_index,
                "track_id": track_id,
                "track_name": track_name,
                **extra,
                "mse": format_float(metrics["mse"]),
                "mae": format_float(metrics["mae"]),
                "pearson": format_float(metrics["pearson"]),
                "n_values": format_float(metrics["n"]),
            }
        )
    return rows


def split_track_label(label: str) -> tuple[str, str]:
    if "\t" not in label:
        return label, label
    return tuple(label.split("\t", maxsplit=1))  # type: ignore[return-value]


def write_rows(path: Path, rows: list[dict[str, object]], preferred_fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames: list[str] = []
    for fieldname in preferred_fieldnames:
        if fieldname not in fieldnames:
            fieldnames.append(fieldname)
    for row in rows:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def append_rows(path: Path, rows: list[dict[str, object]], preferred_fieldnames: list[str]) -> None:
    if not rows:
        return
    if path.exists():
        with path.open() as handle:
            header = handle.readline().rstrip("\n").split("\t")
        with path.open("a", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=header, delimiter="\t")
            for row in rows:
                writer.writerow({field: row.get(field, "") for field in header})
        return
    write_rows(path, rows, preferred_fieldnames)


def update_samples(
    prediction_samples: list[list[np.ndarray]],
    target_samples: list[list[np.ndarray]],
    *,
    prediction: torch.Tensor,
    target: torch.Tensor,
    example_offset: int,
    sample_stride: int,
    sample_offset: int,
) -> None:
    batch_size = int(prediction.shape[0])
    sequence_length = int(prediction.shape[-1])
    local_positions = torch.arange(sequence_length, device=prediction.device)
    for batch_index in range(batch_size):
        global_base = (example_offset + batch_index) * sequence_length
        selected = ((global_base + local_positions + sample_offset) % sample_stride) == 0
        if not bool(selected.any()):
            continue
        pred_np = prediction[batch_index, :, selected].detach().cpu().numpy()
        target_np = target[batch_index, :, selected].detach().cpu().numpy()
        for track_index in range(prediction.shape[1]):
            prediction_samples[track_index].append(pred_np[track_index])
            target_samples[track_index].append(target_np[track_index])


def sampled_spearman_rows(
    *,
    candidate: Candidate,
    track_labels: list[str],
    prediction_samples: list[list[np.ndarray]],
    target_samples: list[list[np.ndarray]],
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for track_index, label in enumerate(track_labels):
        track_id, track_name = split_track_label(label)
        prediction = (
            np.concatenate(prediction_samples[track_index])
            if prediction_samples[track_index]
            else np.array([], dtype=np.float32)
        )
        target = (
            np.concatenate(target_samples[track_index])
            if target_samples[track_index]
            else np.array([], dtype=np.float32)
        )
        spearman = spearman_from_arrays(prediction, target)
        rows.append(
            {
                "rank": candidate.rank,
                "run_id": candidate.run_id,
                "track_index": track_index,
                "track_id": track_id,
                "track_name": track_name,
                "spearman_sampled": format_float(spearman),
                "spearman_sampled_n": len(prediction),
            }
        )
    return rows


def spearman_from_arrays(prediction: np.ndarray, target: np.ndarray) -> float:
    if len(prediction) <= 1:
        return math.nan
    prediction_ranks = rankdata_average(prediction)
    target_ranks = rankdata_average(target)
    return pearson_from_sums(
        float(prediction_ranks.sum()),
        float(target_ranks.sum()),
        float(np.square(prediction_ranks).sum()),
        float(np.square(target_ranks).sum()),
        float((prediction_ranks * target_ranks).sum()),
        float(len(prediction_ranks)),
    )


def evaluate_candidate(
    *,
    candidate: Candidate,
    base_model: AlphaGenome,
    dataloader: torch.utils.data.DataLoader,
    track_labels: list[str],
    manifest_rows: list[dict[str, str]],
    core_ends: list[int],
    genes: list[Gene],
    gene_overlaps: list[list[GeneOverlap]],
    gene_intervals: dict[str, list[Interval]],
    exon_intervals: dict[str, list[Interval]],
    promoter_intervals: dict[str, list[Interval]],
    target_quantiles: dict[float, np.ndarray],
    output_dir: Path,
    device: torch.device,
    spearman_sample_size: int,
    spearman_seed: int,
    pooled_top_factors: list[int],
) -> None:
    n_tracks = len(track_labels)
    model, checkpoint = make_model_from_checkpoint(
        checkpoint_path=candidate.checkpoint,
        base_model=base_model,
        n_tracks=n_tracks,
        device=device,
    )
    model.eval()

    overall_stats = SumStats(n_tracks)
    resolution_stats = {1: SumStats(n_tracks), 128: SumStats(n_tracks), 512: SumStats(n_tracks), 1024: SumStats(n_tracks)}
    gradient_stats = {1: SumStats(n_tracks), 128: SumStats(n_tracks), 512: SumStats(n_tracks), 1024: SumStats(n_tracks)}
    stratum_names = [
        "zero",
        "near_zero_0_0p1",
        "low_0p1_1",
        "medium_1_3",
        "high_gt3",
        "top5_true",
        "top1_true",
    ]
    stratum_stats = {name: SumStats(n_tracks) for name in stratum_names}
    region_names = [
        "promoter_tss_plusminus_1kb",
        "gene_body",
        "exon",
        "intron",
        "intergenic",
    ]
    region_stats = {name: SumStats(n_tracks) for name in region_names}
    window_rows: list[dict[str, object]] = []
    prediction_samples: list[list[np.ndarray]] = [[] for _ in range(n_tracks)]
    target_samples: list[list[np.ndarray]] = [[] for _ in range(n_tracks)]
    pooled_prediction_values: dict[int, list[list[np.ndarray]]] = {
        factor: [[] for _ in range(n_tracks)] for factor in pooled_top_factors
    }
    pooled_target_values: dict[int, list[list[np.ndarray]]] = {
        factor: [[] for _ in range(n_tracks)] for factor in pooled_top_factors
    }
    gene_prediction_sum = np.zeros((len(genes), n_tracks), dtype=np.float64)
    gene_target_sum = np.zeros((len(genes), n_tracks), dtype=np.float64)
    gene_counts = np.zeros((len(genes), n_tracks), dtype=np.float64)

    total_examples = len(dataloader.dataset)
    if dataloader.batch_size is None:
        batch_size = 1
    else:
        batch_size = int(dataloader.batch_size)
    total_positions = max(total_examples * 1_048_576, 1)
    sample_stride = max(total_positions // max(spearman_sample_size, 1), 1)
    sample_offset = spearman_seed % sample_stride
    example_offset = 0

    with torch.inference_mode():
        for batch_index, batch in enumerate(dataloader):
            dna_sequence = batch["dna_sequence"].to(device=device, dtype=torch.float32)
            target = batch["rna_seq"].to(device=device, dtype=torch.float32)
            track_mask = batch["rna_seq_mask"].to(device=device, dtype=torch.bool)
            prediction = model(dna_sequence, target_length=target.shape[-1])
            if isinstance(prediction, dict):
                raise NotImplementedError("Extended diagnostics only support linear predictions")

            overall_stats.update_torch(prediction, target, track_mask)
            update_samples(
                prediction_samples,
                target_samples,
                prediction=prediction,
                target=target,
                example_offset=example_offset,
                sample_stride=sample_stride,
                sample_offset=sample_offset,
            )

            update_window_rows(
                candidate=candidate,
                rows=window_rows,
                batch=batch,
                prediction=prediction,
                target=target,
                track_mask=track_mask,
            )
            update_resolution_and_gradient_stats(
                prediction=prediction,
                target=target,
                track_mask=track_mask,
                resolution_stats=resolution_stats,
                gradient_stats=gradient_stats,
            )
            update_strata_stats(
                prediction=prediction,
                target=target,
                stratum_stats=stratum_stats,
                target_quantiles=target_quantiles,
            )
            update_pooled_values(
                prediction=prediction,
                target=target,
                factors=pooled_top_factors,
                pooled_prediction_values=pooled_prediction_values,
                pooled_target_values=pooled_target_values,
            )
            update_region_and_gene_stats(
                batch_index=batch_index,
                prediction=prediction,
                target=target,
                batch=batch,
                manifest_rows=manifest_rows,
                core_ends=core_ends,
                gene_intervals=gene_intervals,
                exon_intervals=exon_intervals,
                promoter_intervals=promoter_intervals,
                region_stats=region_stats,
                gene_overlaps=gene_overlaps,
                gene_prediction_sum=gene_prediction_sum,
                gene_target_sum=gene_target_sum,
                gene_counts=gene_counts,
            )

            example_offset += int(dna_sequence.shape[0])

    per_track_rows = metric_rows_for_stats(
        candidate=candidate,
        scope="overall_full_resolution",
        stats=overall_stats,
        track_labels=track_labels,
    )
    spearman_rows = sampled_spearman_rows(
        candidate=candidate,
        track_labels=track_labels,
        prediction_samples=prediction_samples,
        target_samples=target_samples,
    )
    spearman_by_track = {int(row["track_index"]): row for row in spearman_rows}
    for row in per_track_rows:
        if row["track_index"] == "":
            row["spearman_sampled"] = ""
            row["spearman_sampled_n"] = ""
            continue
        extra = spearman_by_track[int(row["track_index"])]
        row["spearman_sampled"] = extra["spearman_sampled"]
        row["spearman_sampled_n"] = extra["spearman_sampled_n"]

    config_rows = checkpoint_config_rows(candidate, checkpoint)
    append_standard_outputs(
        output_dir=output_dir,
        candidate=candidate,
        track_labels=track_labels,
        per_track_rows=per_track_rows,
        stratum_stats=stratum_stats,
        region_stats=region_stats,
        resolution_stats=resolution_stats,
        gradient_stats=gradient_stats,
        window_rows=window_rows,
        gene_rows=gene_metric_rows(
            candidate=candidate,
            genes=genes,
            track_labels=track_labels,
            gene_prediction_sum=gene_prediction_sum,
            gene_target_sum=gene_target_sum,
            gene_counts=gene_counts,
        ),
        localization_rows=localization_metric_rows(
            candidate=candidate,
            track_labels=track_labels,
            pooled_prediction_values=pooled_prediction_values,
            pooled_target_values=pooled_target_values,
        ),
        config_rows=config_rows,
    )

    del model
    if device.type == "cuda":
        torch.cuda.empty_cache()


def update_window_rows(
    *,
    candidate: Candidate,
    rows: list[dict[str, object]],
    batch: dict[str, object],
    prediction: torch.Tensor,
    target: torch.Tensor,
    track_mask: torch.Tensor,
) -> None:
    mask = track_mask
    while mask.ndim < prediction.ndim:
        mask = mask.unsqueeze(-1)
    mask = mask.to(dtype=torch.bool, device=prediction.device).expand_as(prediction)
    for sample_index in range(int(prediction.shape[0])):
        sample_mask = mask[sample_index]
        pred_values = prediction[sample_index][sample_mask]
        target_values = target[sample_index][sample_mask]
        error = pred_values.to(torch.float64) - target_values.to(torch.float64)
        pred64 = pred_values.to(torch.float64)
        target64 = target_values.to(torch.float64)
        n = float(sample_mask.sum().detach().cpu())
        rows.append(
            {
                "rank": candidate.rank,
                "run_id": candidate.run_id,
                "window_index": len(rows),
                "interval_chromosome": batch["interval_chromosome"][sample_index],
                "interval_start": int(batch["interval_start"][sample_index]),
                "interval_end": int(batch["interval_end"][sample_index]),
                "mse": format_float(float(error.square().mean().detach().cpu())),
                "mae": format_float(float(error.abs().mean().detach().cpu())),
                "pearson": format_float(
                    pearson_from_sums(
                        float(pred64.sum().detach().cpu()),
                        float(target64.sum().detach().cpu()),
                        float(pred64.square().sum().detach().cpu()),
                        float(target64.square().sum().detach().cpu()),
                        float((pred64 * target64).sum().detach().cpu()),
                        n,
                    )
                ),
                "n_values": format_float(n),
            }
        )


def update_resolution_and_gradient_stats(
    *,
    prediction: torch.Tensor,
    target: torch.Tensor,
    track_mask: torch.Tensor,
    resolution_stats: dict[int, SumStats],
    gradient_stats: dict[int, SumStats],
) -> None:
    for factor, stats in resolution_stats.items():
        if factor == 1:
            pred_pooled = prediction
            target_pooled = target
        else:
            pred_pooled = F.avg_pool1d(prediction, kernel_size=factor, stride=factor)
            target_pooled = F.avg_pool1d(target, kernel_size=factor, stride=factor)
        stats.update_torch(pred_pooled, target_pooled, track_mask)
        if pred_pooled.shape[-1] > 1:
            gradient_stats[factor].update_torch(
                pred_pooled[..., 1:] - pred_pooled[..., :-1],
                target_pooled[..., 1:] - target_pooled[..., :-1],
                track_mask,
            )


def update_strata_stats(
    *,
    prediction: torch.Tensor,
    target: torch.Tensor,
    stratum_stats: dict[str, SumStats],
    target_quantiles: dict[float, np.ndarray],
) -> None:
    strata_masks = {
        "zero": target == 0,
        "near_zero_0_0p1": (target > 0) & (target <= 0.1),
        "low_0p1_1": (target > 0.1) & (target <= 1.0),
        "medium_1_3": (target > 1.0) & (target <= 3.0),
        "high_gt3": target > 3.0,
    }
    q95 = torch.as_tensor(target_quantiles[0.95], device=target.device, dtype=target.dtype).view(1, -1, 1)
    q99 = torch.as_tensor(target_quantiles[0.99], device=target.device, dtype=target.dtype).view(1, -1, 1)
    strata_masks["top5_true"] = target >= q95
    strata_masks["top1_true"] = target >= q99
    for name, mask in strata_masks.items():
        if not bool(mask.any()):
            continue
        values_prediction: list[np.ndarray] = []
        values_target: list[np.ndarray] = []
        for track_index in range(prediction.shape[1]):
            track_mask = mask[:, track_index, :]
            if not bool(track_mask.any()):
                values_prediction.append(np.array([], dtype=np.float32))
                values_target.append(np.array([], dtype=np.float32))
                continue
            values_prediction.append(
                prediction[:, track_index, :][track_mask].detach().cpu().numpy()
            )
            values_target.append(target[:, track_index, :][track_mask].detach().cpu().numpy())
        max_len = max((len(values) for values in values_prediction), default=0)
        if max_len == 0:
            continue
        # SumStats expects [C, N]; update tracks independently when stratum sizes differ.
        for track_index, (pred_values, target_values) in enumerate(
            zip(values_prediction, values_target)
        ):
            if len(pred_values) == 0:
                continue
            temp = SumStats(1)
            temp.update_numpy(
                pred_values.reshape(1, -1),
                target_values.reshape(1, -1),
            )
            copy_track_stats(stratum_stats[name], temp, dst_track=track_index)


def copy_track_stats(dst: SumStats, src: SumStats, *, dst_track: int) -> None:
    dst.sse[dst_track] += src.sse[0]
    dst.sae[dst_track] += src.sae[0]
    dst.sum_x[dst_track] += src.sum_x[0]
    dst.sum_y[dst_track] += src.sum_y[0]
    dst.sum_x2[dst_track] += src.sum_x2[0]
    dst.sum_y2[dst_track] += src.sum_y2[0]
    dst.sum_xy[dst_track] += src.sum_xy[0]
    dst.n[dst_track] += src.n[0]


def update_pooled_values(
    *,
    prediction: torch.Tensor,
    target: torch.Tensor,
    factors: list[int],
    pooled_prediction_values: dict[int, list[list[np.ndarray]]],
    pooled_target_values: dict[int, list[list[np.ndarray]]],
) -> None:
    for factor in factors:
        pred_pooled = F.avg_pool1d(prediction, kernel_size=factor, stride=factor)
        target_pooled = F.avg_pool1d(target, kernel_size=factor, stride=factor)
        pred_np = pred_pooled.detach().cpu().numpy()
        target_np = target_pooled.detach().cpu().numpy()
        for track_index in range(prediction.shape[1]):
            pooled_prediction_values[factor][track_index].append(
                pred_np[:, track_index, :].reshape(-1)
            )
            pooled_target_values[factor][track_index].append(
                target_np[:, track_index, :].reshape(-1)
            )


def update_region_and_gene_stats(
    *,
    batch_index: int,
    prediction: torch.Tensor,
    target: torch.Tensor,
    batch: dict[str, object],
    manifest_rows: list[dict[str, str]],
    core_ends: list[int],
    gene_intervals: dict[str, list[Interval]],
    exon_intervals: dict[str, list[Interval]],
    promoter_intervals: dict[str, list[Interval]],
    region_stats: dict[str, SumStats],
    gene_overlaps: list[list[GeneOverlap]],
    gene_prediction_sum: np.ndarray,
    gene_target_sum: np.ndarray,
    gene_counts: np.ndarray,
) -> None:
    pred_np = prediction.detach().cpu().numpy()
    target_np = target.detach().cpu().numpy()
    batch_size = int(pred_np.shape[0])
    for sample_index in range(batch_size):
        global_index = batch_index * batch_size + sample_index
        row = manifest_rows[global_index]
        chrom = row["chromosome"]
        start = int(row["start"])
        end = int(row["end"])
        core_end = core_ends[global_index]
        core_slice = slice(0, core_end)
        pred_core = pred_np[sample_index, :, core_slice]
        target_core = target_np[sample_index, :, core_slice]

        promoter_mask = make_interval_mask(
            promoter_intervals.get(chrom, []),
            chrom=chrom,
            window_start=start,
            window_end=end,
            local_start=0,
            local_end=core_end,
        )
        gene_mask = make_interval_mask(
            gene_intervals.get(chrom, []),
            chrom=chrom,
            window_start=start,
            window_end=end,
            local_start=0,
            local_end=core_end,
        )
        exon_mask = make_interval_mask(
            exon_intervals.get(chrom, []),
            chrom=chrom,
            window_start=start,
            window_end=end,
            local_start=0,
            local_end=core_end,
        )
        intron_mask = gene_mask & ~exon_mask
        intergenic_mask = ~gene_mask
        region_masks = {
            "promoter_tss_plusminus_1kb": promoter_mask,
            "gene_body": gene_mask,
            "exon": exon_mask,
            "intron": intron_mask,
            "intergenic": intergenic_mask,
        }
        for region_name, mask in region_masks.items():
            if not bool(mask.any()):
                continue
            region_stats[region_name].update_numpy(pred_core[:, mask], target_core[:, mask])

        for overlap in gene_overlaps[global_index]:
            local_start = overlap.local_start
            local_end = overlap.local_end
            if local_end <= local_start:
                continue
            gene_prediction_sum[overlap.gene_index, :] += pred_np[
                sample_index,
                :,
                local_start:local_end,
            ].sum(axis=1)
            gene_target_sum[overlap.gene_index, :] += target_np[
                sample_index,
                :,
                local_start:local_end,
            ].sum(axis=1)
            gene_counts[overlap.gene_index, :] += local_end - local_start


def gene_metric_rows(
    *,
    candidate: Candidate,
    genes: list[Gene],
    track_labels: list[str],
    gene_prediction_sum: np.ndarray,
    gene_target_sum: np.ndarray,
    gene_counts: np.ndarray,
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    n_tracks = len(track_labels)
    valid = gene_counts > 0
    pred_mean = np.divide(
        gene_prediction_sum,
        np.maximum(gene_counts, 1.0),
        where=valid,
        out=np.zeros_like(gene_prediction_sum),
    )
    target_mean = np.divide(
        gene_target_sum,
        np.maximum(gene_counts, 1.0),
        where=valid,
        out=np.zeros_like(gene_target_sum),
    )
    overall_pred = pred_mean[valid]
    overall_target = target_mean[valid]
    rows.append(
        {
            "rank": candidate.rank,
            "run_id": candidate.run_id,
            "metric_scope": "gene_body_mean",
            "track_index": "",
            "track_id": "",
            "track_name": "",
            "n_genes": int(valid.any(axis=1).sum()),
            "mse": format_float(float(np.mean(np.square(overall_pred - overall_target)))),
            "mae": format_float(float(np.mean(np.abs(overall_pred - overall_target)))),
            "pearson": format_float(pearson_np(overall_pred, overall_target)),
            "spearman": format_float(spearman_from_arrays(overall_pred, overall_target)),
        }
    )
    for track_index, label in enumerate(track_labels):
        track_id, track_name = split_track_label(label)
        track_valid = valid[:, track_index]
        pred = pred_mean[track_valid, track_index]
        target = target_mean[track_valid, track_index]
        rows.append(
            {
                "rank": candidate.rank,
                "run_id": candidate.run_id,
                "metric_scope": "gene_body_mean",
                "track_index": track_index,
                "track_id": track_id,
                "track_name": track_name,
                "n_genes": len(pred),
                "mse": format_float(float(np.mean(np.square(pred - target)))),
                "mae": format_float(float(np.mean(np.abs(pred - target)))),
                "pearson": format_float(pearson_np(pred, target)),
                "spearman": format_float(spearman_from_arrays(pred, target)),
            }
        )
    return rows


def pearson_np(prediction: np.ndarray, target: np.ndarray) -> float:
    if prediction.size <= 1:
        return math.nan
    prediction64 = prediction.astype(np.float64, copy=False).reshape(-1)
    target64 = target.astype(np.float64, copy=False).reshape(-1)
    return pearson_from_sums(
        float(prediction64.sum()),
        float(target64.sum()),
        float(np.square(prediction64).sum()),
        float(np.square(target64).sum()),
        float((prediction64 * target64).sum()),
        float(len(prediction64)),
    )


def localization_metric_rows(
    *,
    candidate: Candidate,
    track_labels: list[str],
    pooled_prediction_values: dict[int, list[list[np.ndarray]]],
    pooled_target_values: dict[int, list[list[np.ndarray]]],
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for factor in sorted(pooled_prediction_values):
        all_pred: list[np.ndarray] = []
        all_target: list[np.ndarray] = []
        for track_index, label in enumerate(track_labels):
            track_id, track_name = split_track_label(label)
            pred = np.concatenate(pooled_prediction_values[factor][track_index])
            target = np.concatenate(pooled_target_values[factor][track_index])
            all_pred.append(pred)
            all_target.append(target)
            for quantile, top_label in ((0.95, "top5"), (0.99, "top1")):
                rows.append(
                    {
                        "rank": candidate.rank,
                        "run_id": candidate.run_id,
                        "metric_scope": f"{factor}bp_pooled_high_signal_overlap",
                        "track_index": track_index,
                        "track_id": track_id,
                        "track_name": track_name,
                        "top_label": top_label,
                        **top_overlap_metrics(pred, target, quantile),
                    }
                )
        pred_all = np.concatenate(all_pred)
        target_all = np.concatenate(all_target)
        for quantile, top_label in ((0.95, "top5"), (0.99, "top1")):
            rows.append(
                {
                    "rank": candidate.rank,
                    "run_id": candidate.run_id,
                    "metric_scope": f"{factor}bp_pooled_high_signal_overlap",
                    "track_index": "",
                    "track_id": "",
                    "track_name": "",
                    "top_label": top_label,
                    **top_overlap_metrics(pred_all, target_all, quantile),
                }
            )
    return rows


def top_overlap_metrics(prediction: np.ndarray, target: np.ndarray, quantile: float) -> dict[str, object]:
    n = len(target)
    if n == 0:
        return {
            "k": 0,
            "true_threshold": "nan",
            "pred_threshold": "nan",
            "overlap": 0,
            "precision": "nan",
            "recall": "nan",
            "f1": "nan",
            "jaccard": "nan",
        }
    k = max(1, int(math.ceil(n * (1.0 - quantile))))
    true_indices = np.argpartition(target, n - k)[n - k :]
    pred_indices = np.argpartition(prediction, n - k)[n - k :]
    true_set = set(int(index) for index in true_indices)
    pred_set = set(int(index) for index in pred_indices)
    overlap = len(true_set & pred_set)
    precision = overlap / k
    recall = overlap / k
    f1 = 0.0 if precision + recall == 0 else 2 * precision * recall / (precision + recall)
    jaccard = overlap / max(len(true_set | pred_set), 1)
    return {
        "k": k,
        "true_threshold": format_float(float(np.min(target[true_indices]))),
        "pred_threshold": format_float(float(np.min(prediction[pred_indices]))),
        "overlap": overlap,
        "precision": format_float(precision),
        "recall": format_float(recall),
        "f1": format_float(f1),
        "jaccard": format_float(jaccard),
    }


def checkpoint_config_rows(candidate: Candidate, checkpoint: dict[str, object]) -> list[dict[str, object]]:
    fields = [
        "embedding_resolution",
        "head_type",
        "linear_head_architecture",
        "linear_hidden_channels",
        "linear_input_bottleneck_channels",
        "linear_dilation",
        "linear_kernel_size",
        "linear_target_space",
        "linear_loss_type",
        "smooth_l1_beta",
        "hybrid_loss_alpha",
        "residual_base_checkpoint",
        "linear_residual_fusion",
        "residual_correction_scale_init",
        "residual_base_gate",
        "target_transform",
        "organism_index",
        "n_tracks",
    ]
    row: dict[str, object] = {
        "rank": candidate.rank,
        "run_id": candidate.run_id,
        "checkpoint": str(candidate.checkpoint),
        "source_full_mse": format_float(candidate.full_mse),
        "source_pearson": format_float(candidate.pearson),
    }
    for field in fields:
        row[field] = checkpoint.get(field, "")
    return [row]


def append_standard_outputs(
    *,
    output_dir: Path,
    candidate: Candidate,
    track_labels: list[str],
    per_track_rows: list[dict[str, object]],
    stratum_stats: dict[str, SumStats],
    region_stats: dict[str, SumStats],
    resolution_stats: dict[int, SumStats],
    gradient_stats: dict[int, SumStats],
    window_rows: list[dict[str, object]],
    gene_rows: list[dict[str, object]],
    localization_rows: list[dict[str, object]],
    config_rows: list[dict[str, object]],
) -> None:
    append_rows(
        output_dir / "per_track_metrics.tsv",
        per_track_rows,
        [
            "rank",
            "run_id",
            "metric_scope",
            "track_index",
            "track_id",
            "track_name",
            "mse",
            "mae",
            "pearson",
            "spearman_sampled",
            "spearman_sampled_n",
            "n_values",
        ],
    )
    strata_rows: list[dict[str, object]] = []
    for stratum_name, stats in stratum_stats.items():
        strata_rows.extend(
            metric_rows_for_stats(
                candidate=candidate,
                scope="signal_stratum",
                stats=stats,
                track_labels=track_labels,
                extra={"stratum": stratum_name},
            )
        )
    append_rows(
        output_dir / "signal_strata_metrics.tsv",
        strata_rows,
        [
            "rank",
            "run_id",
            "metric_scope",
            "stratum",
            "track_index",
            "track_id",
            "track_name",
            "mse",
            "mae",
            "pearson",
            "n_values",
        ],
    )
    region_rows: list[dict[str, object]] = []
    for region_name, stats in region_stats.items():
        region_rows.extend(
            metric_rows_for_stats(
                candidate=candidate,
                scope="valid_core_unique_region",
                stats=stats,
                track_labels=track_labels,
                extra={"region": region_name},
            )
        )
    append_rows(
        output_dir / "region_metrics.tsv",
        region_rows,
        [
            "rank",
            "run_id",
            "metric_scope",
            "region",
            "track_index",
            "track_id",
            "track_name",
            "mse",
            "mae",
            "pearson",
            "n_values",
        ],
    )
    resolution_rows: list[dict[str, object]] = []
    for factor, stats in resolution_stats.items():
        resolution_rows.extend(
            metric_rows_for_stats(
                candidate=candidate,
                scope="pooled_value",
                stats=stats,
                track_labels=track_labels,
                extra={"pool_bp": factor},
            )
        )
    for factor, stats in gradient_stats.items():
        resolution_rows.extend(
            metric_rows_for_stats(
                candidate=candidate,
                scope="local_gradient",
                stats=stats,
                track_labels=track_labels,
                extra={"pool_bp": factor},
            )
        )
    append_rows(
        output_dir / "resolution_metrics.tsv",
        resolution_rows,
        [
            "rank",
            "run_id",
            "metric_scope",
            "pool_bp",
            "track_index",
            "track_id",
            "track_name",
            "mse",
            "mae",
            "pearson",
            "n_values",
        ],
    )
    append_rows(
        output_dir / "window_metrics.tsv",
        window_rows,
        [
            "rank",
            "run_id",
            "window_index",
            "interval_chromosome",
            "interval_start",
            "interval_end",
            "mse",
            "mae",
            "pearson",
            "n_values",
        ],
    )
    append_rows(
        output_dir / "gene_metrics.tsv",
        gene_rows,
        [
            "rank",
            "run_id",
            "metric_scope",
            "track_index",
            "track_id",
            "track_name",
            "n_genes",
            "mse",
            "mae",
            "pearson",
            "spearman",
        ],
    )
    append_rows(
        output_dir / "high_signal_localization_metrics.tsv",
        localization_rows,
        [
            "rank",
            "run_id",
            "metric_scope",
            "track_index",
            "track_id",
            "track_name",
            "top_label",
            "k",
            "true_threshold",
            "pred_threshold",
            "overlap",
            "precision",
            "recall",
            "f1",
            "jaccard",
        ],
    )
    append_rows(
        output_dir / "model_config.tsv",
        config_rows,
        [
            "rank",
            "run_id",
            "checkpoint",
            "source_full_mse",
            "source_pearson",
            "embedding_resolution",
            "head_type",
            "linear_head_architecture",
            "linear_hidden_channels",
            "linear_input_bottleneck_channels",
            "linear_dilation",
            "linear_kernel_size",
            "linear_target_space",
            "linear_loss_type",
            "smooth_l1_beta",
            "hybrid_loss_alpha",
            "residual_base_checkpoint",
            "linear_residual_fusion",
            "residual_correction_scale_init",
            "residual_base_gate",
            "target_transform",
            "organism_index",
            "n_tracks",
        ],
    )


def write_run_metadata(
    *,
    output_dir: Path,
    args: argparse.Namespace,
    candidates: list[Candidate],
    manifest_rows: list[dict[str, str]],
    genes: list[Gene],
) -> None:
    metadata = {
        "dataset_dir": args.dataset_dir,
        "gtf": args.gtf,
        "weights": args.weights,
        "top_n": args.top_n,
        "rank_sort": "valid full MSE ascending, Pearson descending tie-breaker",
        "test_split_used": False,
        "valid_examples": len(manifest_rows),
        "valid_chromosomes": sorted({row["chromosome"] for row in manifest_rows}),
        "gene_count_in_gtf_valid_chroms": len(genes),
        "promoter_radius_bp": args.promoter_radius_bp,
        "spearman_sample_size": args.spearman_sample_size,
        "selected_run_ids": [candidate.run_id for candidate in candidates],
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    with (output_dir / "diagnostic_metadata.json").open("w") as handle:
        json.dump(metadata, handle, indent=2, sort_keys=True)


def parse_int_list(value: str) -> list[int]:
    return [int(part.strip()) for part in value.split(",") if part.strip()]


def main() -> None:
    args = parse_args()
    dataset_dir = Path(args.dataset_dir)
    if "test" in dataset_dir.name.lower():
        raise ValueError("This diagnostics script must not be run on the test split")
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    candidates = collect_top_candidates(Path(args.runs_dir), args.top_n)
    write_top_models(output_dir / "top15_models.tsv", candidates)
    if args.only_write_top15:
        print(f"top_models\t{output_dir / 'top15_models.tsv'}")
        return

    manifest_rows = read_manifest(dataset_dir, max_examples=args.max_examples)
    core_ends = compute_core_ends(manifest_rows)
    chroms = {row["chromosome"] for row in manifest_rows}
    chrom_lengths = {
        row["chromosome"]: max(
            int(candidate_row["end"])
            for candidate_row in manifest_rows
            if candidate_row["chromosome"] == row["chromosome"]
        )
        for row in manifest_rows
    }
    genes, gene_intervals, exon_intervals, promoter_intervals = load_gtf_annotations(
        Path(args.gtf),
        chroms=chroms,
        promoter_radius_bp=args.promoter_radius_bp,
        chrom_lengths=chrom_lengths,
    )
    gene_overlaps = build_gene_overlaps(genes, manifest_rows, core_ends)
    dataloader = make_dataloader(
        dataset_dir,
        batch_size=args.batch_size,
        shuffle=False,
        num_workers=args.num_workers,
        target_transform="log1p",
        max_examples=args.max_examples,
    )
    n_tracks = int(dataloader.dataset.metadata["n_tracks"])
    track_labels = load_track_names(dataloader.dataset.metadata, n_tracks)
    target_quantiles = compute_target_quantiles(
        dataset_dir,
        n_tracks=n_tracks,
        quantiles=(0.95, 0.99),
        max_examples=args.max_examples,
    )
    write_run_metadata(
        output_dir=output_dir,
        args=args,
        candidates=candidates,
        manifest_rows=manifest_rows,
        genes=genes,
    )

    device = pick_device(args.device)
    print(f"device\t{device}")
    print(f"selected_models\t{len(candidates)}")
    print(f"valid_examples\t{len(manifest_rows)}")
    print(f"genes\t{len(genes)}")
    if device.type == "cuda":
        print(f"cuda_device\t{torch.cuda.get_device_name(device)}")
    base_model = AlphaGenome.from_pretrained(Path(args.weights), device=device)
    pooled_top_factors = parse_int_list(args.pooled_top_factors)
    selected = [
        candidate
        for index, candidate in enumerate(candidates)
        if index % args.num_shards == args.rank_shard
    ]
    print(
        "shard\t"
        f"{args.rank_shard}/{args.num_shards}\t"
        f"ranks={','.join(str(candidate.rank) for candidate in selected)}"
    )

    for candidate in selected:
        print(
            "evaluate\t"
            f"rank={candidate.rank}\t"
            f"mse={candidate.full_mse:.8f}\t"
            f"pearson={candidate.pearson:.8f}\t"
            f"run={candidate.run_id}"
        )
        evaluate_candidate(
            candidate=candidate,
            base_model=base_model,
            dataloader=dataloader,
            track_labels=track_labels,
            manifest_rows=manifest_rows,
            core_ends=core_ends,
            genes=genes,
            gene_overlaps=gene_overlaps,
            gene_intervals=gene_intervals,
            exon_intervals=exon_intervals,
            promoter_intervals=promoter_intervals,
            target_quantiles=target_quantiles,
            output_dir=output_dir,
            device=device,
            spearman_sample_size=args.spearman_sample_size,
            spearman_seed=args.spearman_seed,
            pooled_top_factors=pooled_top_factors,
        )
    print(f"output_dir\t{output_dir}")


if __name__ == "__main__":
    main()
