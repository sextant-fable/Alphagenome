#!/usr/bin/env python3
"""Dependency-light utilities for v2 biological validation applications.

This module is deliberately separate from the controlled training/evaluation
entry points.  It never resolves the chromosome-X final-test manifest and it
contains no controller state mutations.
"""

from __future__ import annotations

from collections import defaultdict
import csv
from dataclasses import dataclass
import hashlib
import json
import math
from pathlib import Path
from typing import Iterable, Mapping, Sequence

import numpy as np


AUTOSOMES = frozenset({"I", "II", "III", "IV", "V"})
DEVELOPMENT_CHROMOSOMES = AUTOSOMES | {"X"}
FORBIDDEN_SPLIT_TOKENS = ("test_locked", "final_test", "final-test")


def sha256(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        while chunk := handle.read(8 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def read_tsv(path: str | Path) -> list[dict[str, str]]:
    with Path(path).open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(
    path: str | Path,
    rows: Sequence[Mapping[str, object]],
    fieldnames: Sequence[str] | None = None,
) -> None:
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    if fieldnames is None:
        if not rows:
            raise ValueError("fieldnames are required when writing an empty TSV")
        fieldnames = list(rows[0])
    temporary = destination.with_suffix(destination.suffix + ".tmp")
    with temporary.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(fieldnames),
            delimiter="\t",
            lineterminator="\n",
            extrasaction="raise",
        )
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(destination)


def write_json(path: str | Path, payload: Mapping[str, object]) -> None:
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_suffix(destination.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    temporary.replace(destination)


def assert_no_final_test_path(path: str | Path) -> None:
    value = str(Path(path)).lower()
    if any(token in value for token in FORBIDDEN_SPLIT_TOKENS):
        raise RuntimeError(f"Final-test path is prohibited for this application: {path}")


@dataclass(frozen=True)
class ValidationSubwindow:
    interval_index: int
    chromosome: str
    start: int
    end: int
    block_id: str


def load_fold_validation_subwindows(
    path: str | Path,
    *,
    expected_fold: int,
    sequence_length: int,
) -> tuple[list[dict[str, str]], list[ValidationSubwindow]]:
    """Load only aligned fold-validation cores and reject all test-like input."""

    assert_no_final_test_path(path)
    if expected_fold not in range(1, 6):
        raise ValueError("expected_fold must be in 1..5")
    if sequence_length <= 0 or sequence_length % 128:
        raise ValueError("sequence_length must be positive and divisible by 128")
    rows = read_tsv(path)
    if not rows:
        raise ValueError(f"Empty validation manifest: {path}")
    required = {
        "chromosome",
        "start",
        "end",
        "core_start",
        "core_end",
        "fold",
        "role",
        "block_id",
    }
    missing = required - set(rows[0])
    if missing:
        raise ValueError(f"Validation manifest lacks fields: {sorted(missing)}")

    subwindows: list[ValidationSubwindow] = []
    for interval_index, row in enumerate(rows):
        if row["role"] != "valid":
            raise RuntimeError(
                f"Only role=valid is allowed; row {interval_index + 2} has {row['role']!r}"
            )
        if int(row["fold"]) != expected_fold:
            raise RuntimeError(
                f"Fold mismatch at row {interval_index + 2}: {row['fold']} != {expected_fold}"
            )
        chromosome = row["chromosome"]
        if chromosome not in DEVELOPMENT_CHROMOSOMES:
            raise RuntimeError(f"Unexpected chromosome in development validation: {chromosome}")
        window_start = int(row["start"])
        window_end = int(row["end"])
        core_start = int(row["core_start"])
        core_end = int(row["core_end"])
        if not window_start <= core_start < core_end <= window_end:
            raise RuntimeError(f"Invalid context/core bounds at row {interval_index + 2}")
        if (core_start - window_start) % 128 or (core_end - core_start) % 128:
            raise RuntimeError(f"Unaligned core at row {interval_index + 2}")
        for start in range(core_start, core_end - sequence_length + 1, sequence_length):
            subwindows.append(
                ValidationSubwindow(
                    interval_index=interval_index,
                    chromosome=chromosome,
                    start=start,
                    end=start + sequence_length,
                    block_id=row["block_id"],
                )
            )
    if not subwindows:
        raise RuntimeError("No complete core-only validation subwindows")
    by_chromosome: dict[str, list[tuple[int, int]]] = defaultdict(list)
    for window in subwindows:
        by_chromosome[window.chromosome].append((window.start, window.end))
    for chromosome, intervals in by_chromosome.items():
        intervals.sort()
        for previous, current in zip(intervals, intervals[1:]):
            if current[0] < previous[1]:
                raise RuntimeError(
                    f"Overlapping core subwindows would double-count {chromosome}: "
                    f"{previous} and {current}"
                )
    return rows, subwindows


def merge_intervals(
    intervals: Iterable[tuple[int, int]],
) -> tuple[tuple[int, int], ...]:
    merged: list[tuple[int, int]] = []
    for start, end in sorted(intervals):
        if end <= start:
            continue
        if not merged or start > merged[-1][1]:
            merged.append((start, end))
        else:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
    return tuple(merged)


def parse_gtf_attributes(value: str) -> dict[str, str]:
    attributes: dict[str, str] = {}
    for part in value.strip().split(";"):
        part = part.strip()
        if not part:
            continue
        key, separator, raw = part.partition(" ")
        attributes[key] = raw.strip().strip('"') if separator else ""
    return attributes


@dataclass(frozen=True)
class AnnotatedGene:
    gene_id: str
    gene_name: str
    chromosome: str
    strand: str
    gene_start: int
    gene_end: int
    exon_intervals: tuple[tuple[int, int], ...]

    @property
    def exon_bases(self) -> int:
        return sum(end - start for start, end in self.exon_intervals)


def load_annotated_genes(
    path: str | Path,
    chromosomes: set[str] | frozenset[str] | None = None,
) -> tuple[list[AnnotatedGene], dict[str, list[int]]]:
    gene_bounds: dict[tuple[str, str, str], tuple[int, int, str]] = {}
    exon_intervals: dict[tuple[str, str, str], list[tuple[int, int]]] = defaultdict(list)
    with Path(path).open() as handle:
        for line_number, line in enumerate(handle, 1):
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] not in {"gene", "exon"}:
                continue
            chromosome, feature, strand = fields[0], fields[2], fields[6]
            if chromosomes is not None and chromosome not in chromosomes:
                continue
            attributes = parse_gtf_attributes(fields[8])
            gene_id = attributes.get("gene_id")
            if not gene_id:
                raise ValueError(f"GTF feature lacks gene_id at line {line_number}")
            key = (chromosome, gene_id, strand)
            start, end = int(fields[3]) - 1, int(fields[4])
            if feature == "gene":
                gene_bounds[key] = (start, end, attributes.get("gene_name", gene_id))
            else:
                exon_intervals[key].append((start, end))
                if key not in gene_bounds:
                    gene_bounds[key] = (start, end, attributes.get("gene_name", gene_id))
                else:
                    old_start, old_end, old_name = gene_bounds[key]
                    gene_bounds[key] = (min(old_start, start), max(old_end, end), old_name)
    genes: list[AnnotatedGene] = []
    for key, raw_exons in exon_intervals.items():
        chromosome, gene_id, strand = key
        start, end, name = gene_bounds[key]
        merged = merge_intervals(raw_exons)
        if merged:
            genes.append(
                AnnotatedGene(gene_id, name, chromosome, strand, start, end, merged)
            )
    genes.sort(key=lambda gene: (gene.chromosome, gene.gene_start, gene.gene_end, gene.gene_id))
    by_chromosome: dict[str, list[int]] = defaultdict(list)
    for index, gene in enumerate(genes):
        by_chromosome[gene.chromosome].append(index)
    return genes, dict(by_chromosome)


class Dpy27GeneExonAccumulator:
    """Accumulate two-track predictions and observations over unique exon bases."""

    def __init__(
        self,
        genes: Sequence[AnnotatedGene],
        by_chromosome: Mapping[str, Sequence[int]],
    ) -> None:
        self.genes = list(genes)
        self.by_chromosome = {key: list(value) for key, value in by_chromosome.items()}
        self.prediction_sum = np.zeros((len(genes), 2), dtype=np.float64)
        self.observation_sum = np.zeros((len(genes), 2), dtype=np.float64)
        self.base_count = np.zeros(len(genes), dtype=np.int64)
        self.block_ids: list[set[str]] = [set() for _ in genes]

    def update(
        self,
        chromosome: str,
        start: int,
        prediction: np.ndarray,
        observation: np.ndarray,
        *,
        block_id: str,
    ) -> None:
        prediction = np.asarray(prediction)
        observation = np.asarray(observation)
        if prediction.shape != observation.shape or prediction.ndim != 2:
            raise ValueError("prediction and observation must have equal [2, positions] shape")
        if prediction.shape[0] != 2:
            raise ValueError("DPY-27 accumulation requires exactly two ordered tracks")
        if not np.isfinite(prediction).all() or not np.isfinite(observation).all():
            raise ValueError("Non-finite signal cannot enter gene-level accumulation")
        end = start + prediction.shape[1]
        for gene_index in self.by_chromosome.get(chromosome, []):
            gene = self.genes[gene_index]
            if gene.gene_end <= start:
                continue
            if gene.gene_start >= end:
                break
            for exon_start, exon_end in gene.exon_intervals:
                local_start = max(exon_start, start) - start
                local_end = min(exon_end, end) - start
                if local_end <= local_start:
                    continue
                self.prediction_sum[gene_index] += prediction[:, local_start:local_end].sum(
                    axis=1, dtype=np.float64
                )
                self.observation_sum[gene_index] += observation[:, local_start:local_end].sum(
                    axis=1, dtype=np.float64
                )
                self.base_count[gene_index] += local_end - local_start
                self.block_ids[gene_index].add(block_id)

    def rows(self, *, fold: int, seed: int) -> list[dict[str, object]]:
        rows: list[dict[str, object]] = []
        for index, gene in enumerate(self.genes):
            evaluated = int(self.base_count[index])
            if evaluated == 0:
                continue
            expected = gene.exon_bases
            if evaluated > expected:
                raise RuntimeError(
                    f"Exon bases were double-counted for {gene.gene_id}: {evaluated} > {expected}"
                )
            if len(self.block_ids[index]) != 1:
                raise RuntimeError(
                    f"Gene {gene.gene_id} overlaps multiple validation blocks: "
                    f"{sorted(self.block_ids[index])}"
                )
            predicted = self.prediction_sum[index] / evaluated
            observed = self.observation_sum[index] / evaluated
            rows.append(
                {
                    "analysis_label": "internal_genomic_block_validation",
                    "fold": fold,
                    "seed": seed,
                    "gene_id": gene.gene_id,
                    "gene_name": gene.gene_name,
                    "chromosome": gene.chromosome,
                    "strand": gene.strand,
                    "block_id": next(iter(self.block_ids[index])),
                    "exon_bases_evaluated": evaluated,
                    "exon_bases_total": expected,
                    "exon_coverage_fraction": evaluated / expected,
                    "complete_exon_coverage": evaluated == expected,
                    "predicted_dpy27_rnai_mean": float(predicted[0]),
                    "predicted_vector_rnai_mean": float(predicted[1]),
                    "observed_dpy27_rnai_mean": float(observed[0]),
                    "observed_vector_rnai_mean": float(observed[1]),
                    "predicted_log1p_contrast_g0005_minus_g0006": float(
                        np.log1p(max(predicted[0], 0.0))
                        - np.log1p(max(predicted[1], 0.0))
                    ),
                    "observed_log1p_contrast_g0005_minus_g0006": float(
                        np.log1p(max(observed[0], 0.0))
                        - np.log1p(max(observed[1], 0.0))
                    ),
                }
            )
        return rows


def average_ranks(values: Sequence[float] | np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=np.float64)
    if values.ndim != 1:
        raise ValueError("values must be one-dimensional")
    order = np.argsort(values, kind="mergesort")
    ranks = np.empty(len(values), dtype=np.float64)
    cursor = 0
    while cursor < len(values):
        stop = cursor + 1
        while stop < len(values) and values[order[stop]] == values[order[cursor]]:
            stop += 1
        ranks[order[cursor:stop]] = (cursor + stop - 1) / 2.0
        cursor = stop
    return ranks


def spearman_correlation(
    x: Sequence[float] | np.ndarray,
    y: Sequence[float] | np.ndarray,
) -> float:
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    if x.shape != y.shape or x.ndim != 1:
        raise ValueError("x and y must be equal one-dimensional arrays")
    finite = np.isfinite(x) & np.isfinite(y)
    if finite.sum() < 2:
        return math.nan
    rank_x = average_ranks(x[finite])
    rank_y = average_ranks(y[finite])
    centered_x = rank_x - rank_x.mean()
    centered_y = rank_y - rank_y.mean()
    denominator = math.sqrt(float(np.dot(centered_x, centered_x) * np.dot(centered_y, centered_y)))
    if denominator == 0:
        return math.nan
    return float(np.dot(centered_x, centered_y) / denominator)


def direction_concordance(
    predicted: Sequence[float] | np.ndarray,
    observed: Sequence[float] | np.ndarray,
    *,
    epsilon: float = 0.0,
) -> tuple[float, int]:
    predicted = np.asarray(predicted, dtype=np.float64)
    observed = np.asarray(observed, dtype=np.float64)
    if predicted.shape != observed.shape:
        raise ValueError("predicted and observed must have equal shape")
    eligible = (
        np.isfinite(predicted)
        & np.isfinite(observed)
        & (np.abs(predicted) > epsilon)
        & (np.abs(observed) > epsilon)
    )
    count = int(eligible.sum())
    if count == 0:
        return math.nan, 0
    return float(np.mean(np.sign(predicted[eligible]) == np.sign(observed[eligible]))), count


def dpy27_run_endpoints(gene_rows: Sequence[Mapping[str, object]]) -> dict[str, object]:
    primary = [row for row in gene_rows if _as_bool(row["complete_exon_coverage"])]
    x_rows = [row for row in primary if row["chromosome"] == "X"]
    autosomal = [row for row in primary if row["chromosome"] in AUTOSOMES]
    if len(x_rows) < 2 or len(autosomal) < 2:
        raise RuntimeError("At least two complete X-linked and autosomal genes are required")
    predicted = np.array(
        [float(row["predicted_log1p_contrast_g0005_minus_g0006"]) for row in primary]
    )
    observed = np.array(
        [float(row["observed_log1p_contrast_g0005_minus_g0006"]) for row in primary]
    )
    chromosomes = np.array([str(row["chromosome"]) for row in primary])
    x_mask = chromosomes == "X"
    auto_mask = np.array([chromosome in AUTOSOMES for chromosome in chromosomes])
    concordance, concordance_n = direction_concordance(predicted[x_mask], observed[x_mask])
    return {
        "n_complete_genes": len(primary),
        "n_x_genes": int(x_mask.sum()),
        "n_autosomal_genes": int(auto_mask.sum()),
        "predicted_median_contrast_x_minus_autosomes": float(
            np.median(predicted[x_mask]) - np.median(predicted[auto_mask])
        ),
        "observed_median_contrast_x_minus_autosomes": float(
            np.median(observed[x_mask]) - np.median(observed[auto_mask])
        ),
        "x_predicted_observed_direction_concordance": concordance,
        "x_direction_concordance_n": concordance_n,
        "predicted_observed_gene_spearman": spearman_correlation(predicted, observed),
        "predicted_observed_x_gene_spearman": spearman_correlation(
            predicted[x_mask], observed[x_mask]
        ),
        "predicted_x_positive_fraction": float(np.mean(predicted[x_mask] > 0)),
        "observed_x_positive_fraction": float(np.mean(observed[x_mask] > 0)),
    }


PRIMARY_ENDPOINTS = (
    "predicted_median_contrast_x_minus_autosomes",
    "observed_median_contrast_x_minus_autosomes",
    "x_predicted_observed_direction_concordance",
    "predicted_observed_gene_spearman",
)


def _as_bool(value: object) -> bool:
    if isinstance(value, bool):
        return value
    if isinstance(value, str) and value in {"True", "False"}:
        return value == "True"
    raise ValueError(f"Expected a strict boolean, received {value!r}")


def fold_seed_mean_gene_rows(
    gene_rows: Sequence[Mapping[str, object]],
) -> list[dict[str, object]]:
    """Average condition predictions across seeds before computing contrasts."""

    grouped: dict[tuple[int, str], list[Mapping[str, object]]] = defaultdict(list)
    for row in gene_rows:
        grouped[(int(row["fold"]), str(row["gene_id"]))].append(row)
    result: list[dict[str, object]] = []
    for (fold, gene_id), rows in sorted(grouped.items()):
        seeds = {int(row["seed"]) for row in rows}
        if len(rows) != 3 or len(seeds) != 3:
            raise RuntimeError(
                f"Fold {fold} gene {gene_id} must have exactly three unique seeds"
            )
        invariant_fields = (
            "gene_id",
            "gene_name",
            "chromosome",
            "strand",
            "block_id",
            "exon_bases_evaluated",
            "exon_bases_total",
            "complete_exon_coverage",
        )
        for field in invariant_fields:
            if len({str(row[field]) for row in rows}) != 1:
                raise RuntimeError(
                    f"Fold {fold} gene {gene_id} differs across seeds for {field}"
                )
        observed_dpy = np.array(
            [float(row["observed_dpy27_rnai_mean"]) for row in rows], dtype=np.float64
        )
        observed_control = np.array(
            [float(row["observed_vector_rnai_mean"]) for row in rows], dtype=np.float64
        )
        if not (
            np.allclose(observed_dpy, observed_dpy[0], rtol=0, atol=1e-10)
            and np.allclose(observed_control, observed_control[0], rtol=0, atol=1e-10)
        ):
            raise RuntimeError(
                f"Observed signals differ across seeds for fold {fold} gene {gene_id}"
            )
        predicted_dpy = float(
            np.mean([float(row["predicted_dpy27_rnai_mean"]) for row in rows])
        )
        predicted_control = float(
            np.mean([float(row["predicted_vector_rnai_mean"]) for row in rows])
        )
        template = rows[0]
        result.append(
            {
                "analysis_label": "internal_genomic_block_validation",
                "fold": fold,
                "seed": "seed_mean",
                "gene_id": gene_id,
                "gene_name": template["gene_name"],
                "chromosome": template["chromosome"],
                "strand": template["strand"],
                "block_id": template["block_id"],
                "exon_bases_evaluated": int(template["exon_bases_evaluated"]),
                "exon_bases_total": int(template["exon_bases_total"]),
                "exon_coverage_fraction": float(template["exon_coverage_fraction"]),
                "complete_exon_coverage": _as_bool(template["complete_exon_coverage"]),
                "predicted_dpy27_rnai_mean": predicted_dpy,
                "predicted_vector_rnai_mean": predicted_control,
                "observed_dpy27_rnai_mean": float(observed_dpy[0]),
                "observed_vector_rnai_mean": float(observed_control[0]),
                "predicted_log1p_contrast_g0005_minus_g0006": float(
                    np.log1p(max(predicted_dpy, 0.0))
                    - np.log1p(max(predicted_control, 0.0))
                ),
                "observed_log1p_contrast_g0005_minus_g0006": float(
                    np.log1p(max(observed_dpy[0], 0.0))
                    - np.log1p(max(observed_control[0], 0.0))
                ),
            }
        )
    return result


def fold_primary_summary(
    gene_rows: Sequence[Mapping[str, object]],
    *,
    endpoints: Sequence[str] = PRIMARY_ENDPOINTS,
    bootstrap_replicates: int = 10_000,
    bootstrap_seed: int = 20260813,
) -> tuple[
    list[dict[str, object]],
    list[dict[str, object]],
    list[dict[str, object]],
]:
    """Average per-condition predictions within gene/fold, then infer by fold."""

    if bootstrap_replicates < 100:
        raise ValueError("bootstrap_replicates must be at least 100")
    fold_gene_rows = fold_seed_mean_gene_rows(gene_rows)
    grouped: dict[int, list[Mapping[str, object]]] = defaultdict(list)
    for row in fold_gene_rows:
        grouped[int(row["fold"])].append(row)
    if sorted(grouped) != [1, 2, 3, 4, 5]:
        raise RuntimeError("Fold-primary inference requires exactly folds 1..5")
    fold_rows: list[dict[str, object]] = []
    for fold in sorted(grouped):
        rows = grouped[fold]
        endpoints_for_fold = dpy27_run_endpoints(rows)
        record: dict[str, object] = {
            "analysis_label": "internal_genomic_block_validation",
            "fold": fold,
            "n_seeds": 3,
        }
        for endpoint in endpoints:
            record[endpoint] = endpoints_for_fold[endpoint]
        fold_rows.append(record)

    rng = np.random.default_rng(bootstrap_seed)
    summaries: list[dict[str, object]] = []
    for endpoint in endpoints:
        values = np.array([float(row[endpoint]) for row in fold_rows], dtype=np.float64)
        if not np.isfinite(values).all():
            raise RuntimeError(
                f"Fold-primary endpoint {endpoint} must be finite in all five folds"
            )
        n_folds = len(values)
        estimate = float(values.mean())
        indices = rng.integers(0, n_folds, size=(bootstrap_replicates, n_folds))
        bootstrap = values[indices].mean(axis=1)
        lower, upper = (float(value) for value in np.quantile(bootstrap, [0.025, 0.975]))
        summaries.append(
            {
                "analysis_label": "internal_genomic_block_validation",
                "endpoint": endpoint,
                "estimate": estimate,
                "ci_95_low": lower,
                "ci_95_high": upper,
                "n_folds": n_folds,
                "n_runs": 15,
                "inference_unit": "fold_after_within_fold_gene_condition_seed_mean",
                "bootstrap_replicates": bootstrap_replicates,
                "bootstrap_seed": bootstrap_seed,
            }
        )
    return fold_gene_rows, fold_rows, summaries
