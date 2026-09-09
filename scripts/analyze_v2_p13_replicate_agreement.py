#!/usr/bin/env python3
"""Measure training-aggregate versus held-out-unit agreement for P13.

The analysis is read-only and CPU-only. It compares each training aggregate
with the corresponding evaluation-only held-out biological-unit aggregate over
the same P12 validation cores. This is an empirical replicate-agreement
reference, not an external benchmark or a formal upper bound on model
performance.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np

from scripts import v2_validation_metrics as full_metrics
from scripts.v2_bigwig_dataset import V2BigWigDataset


REPO_ROOT = Path(__file__).resolve().parents[1]
ROOT = REPO_ROOT / "alphagenome_custom/metadata/v2/replicate_holdout_v1"
INTERVAL_ROOT = REPO_ROOT / "alphagenome_custom/intervals/v2"
SEQUENCE_LENGTH = 131072


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=REPO_ROOT / "results/v2_p13_replicate_agreement",
        help="Ignored directory for source data and audit records.",
    )
    parser.add_argument(
        "--interval-root",
        type=Path,
        default=INTERVAL_ROOT,
        help=(
            "Directory containing fold_1 through fold_5 valid.tsv manifests. "
            "The default is intentionally rejected because it includes chromosome X."
        ),
    )
    parser.add_argument("--max-io-workers", type=int, default=16)
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def core_subwindows(
    intervals: list[dict[str, str]], sequence_length: int
) -> tuple[list[tuple[int, int]], int]:
    result: list[tuple[int, int]] = []
    eligible_bases = 0
    for index, row in enumerate(intervals):
        window_start = int(row["start"])
        window_end = int(row["end"])
        core_start = int(row["core_start"])
        core_end = int(row["core_end"])
        crop_offset = core_start - window_start
        core_length = core_end - core_start
        if not (
            window_start <= core_start < core_end <= window_end
            and crop_offset % 128 == 0
            and core_length % 128 == 0
        ):
            raise ValueError(f"Invalid 128-bp-aligned core for interval {row}")
        eligible_bases += core_length
        for offset in range(
            crop_offset,
            crop_offset + core_length - sequence_length + 1,
            sequence_length,
        ):
            result.append((index, offset))
    if not result:
        raise RuntimeError("No complete core-only subwindows")
    return result, eligible_bases


def per_track_record(
    group_ids: list[str],
    gene_metrics: dict[str, object],
    stats_128: full_metrics.PerTrackStats,
) -> tuple[dict[str, float], dict[str, float], dict[str, float]]:
    gene = np.asarray(gene_metrics["per_track_pearson"], dtype=np.float64)
    coverage = np.asarray(stats_128.per_track()["pearson"], dtype=np.float64)
    primary = 0.5 * (gene + coverage)
    return (
        {group_id: float(value) for group_id, value in zip(group_ids, gene, strict=True)},
        {group_id: float(value) for group_id, value in zip(group_ids, coverage, strict=True)},
        {group_id: float(value) for group_id, value in zip(group_ids, primary, strict=True)},
    )


def finite_mean(values: list[float]) -> float:
    array = np.asarray(values, dtype=np.float64)
    finite = array[np.isfinite(array)]
    if not finite.size:
        raise RuntimeError("No finite values")
    return float(finite.mean())


def require_non_x_development_intervals(path: Path) -> None:
    """Reject P13 inputs that would read chromosome X or locked test blocks.

    This guard intentionally runs before any :class:`V2BigWigDataset` is
    constructed, so a rejected interval manifest cannot open a BigWig handle.
    A future replicate-agreement analysis must use an explicitly registered
    I--V-only interval contract that is also comparable to its model scores.
    """

    forbidden = [
        row
        for row in read_tsv(path)
        if row["chromosome"] == "X" or row["role"] == "test_locked"
    ]
    if forbidden:
        chromosome_x_count = sum(row["chromosome"] == "X" for row in forbidden)
        locked_count = sum(row["role"] == "test_locked" for row in forbidden)
        raise PermissionError(
            f"P13 replicate agreement refuses {path}: {chromosome_x_count} "
            f"chromosome-X rows and {locked_count} locked-test rows are prohibited. "
            "No BigWig access is permitted without a separately registered I-V-only contract."
        )


def main() -> None:
    args = parse_args()
    if SEQUENCE_LENGTH % 128:
        raise ValueError("Sequence length must be divisible by 128")

    primary_index_path = ROOT / "heldout_primary_track_indices.tsv"
    training_track_path = ROOT / "training_track_manifest.tsv"
    training_group_path = ROOT / "training_group_manifest.tsv"
    heldout_track_path = ROOT / "heldout_primary_track_manifest.tsv"
    heldout_group_path = ROOT / "heldout_primary_group_manifest.tsv"
    primary_index_rows = [
        row
        for row in read_tsv(primary_index_path)
        if row["evaluation_primary"].strip().lower() == "true"
    ]
    training_indices = [int(row["track_index"]) for row in primary_index_rows]
    expected_group_ids = [row["group_id"] for row in primary_index_rows]
    if len(training_indices) != 57 or len(set(training_indices)) != len(training_indices):
        raise ValueError("Expected 57 unique evaluated primary held-out track indices")
    if len(set(expected_group_ids)) != len(expected_group_ids):
        raise ValueError("Primary held-out group IDs must be unique")

    # The training manifest retains all 241 output heads and therefore uses the
    # original global indices. The held-out manifest contains only candidates,
    # so its selected rows must be addressed by a separate local index map.
    heldout_manifest_rows = read_tsv(heldout_track_path)
    heldout_local_index = {
        row["group_id"]: index for index, row in enumerate(heldout_manifest_rows)
    }
    if len(heldout_local_index) != len(heldout_manifest_rows):
        raise ValueError("Held-out track manifest contains duplicate group IDs")
    try:
        heldout_indices = [heldout_local_index[group_id] for group_id in expected_group_ids]
    except KeyError as error:
        raise ValueError(f"Primary group missing from held-out track manifest: {error}") from error

    interval_root = args.interval_root.resolve()
    # Perform all access control before opening either aggregate manifest.
    for fold in range(1, 6):
        require_non_x_development_intervals(interval_root / f"fold_{fold}/valid.tsv")

    output = args.output_dir
    output.mkdir(parents=True, exist_ok=True)
    source_rows: list[dict[str, Any]] = []
    fold_rows: list[dict[str, Any]] = []

    for fold in range(1, 6):
        intervals_path = interval_root / f"fold_{fold}/valid.tsv"
        training = V2BigWigDataset(
            intervals_path,
            track_manifest_path=training_track_path,
            group_manifest_path=training_group_path,
            track_indices=training_indices,
            max_io_workers=args.max_io_workers,
        )
        heldout = V2BigWigDataset(
            intervals_path,
            track_manifest_path=heldout_track_path,
            group_manifest_path=heldout_group_path,
            track_indices=heldout_indices,
            max_io_workers=args.max_io_workers,
        )
        try:
            group_ids = [row["group_id"] for row in training.track_rows]
            if group_ids != expected_group_ids:
                raise ValueError(
                    f"Training group ordering differs from primary index manifest in fold {fold}"
                )
            if group_ids != [row["group_id"] for row in heldout.track_rows]:
                raise ValueError(f"Training and held-out group ordering differs in fold {fold}")
            if any(row["role"] == "test_locked" for row in training.intervals):
                raise PermissionError("P13 replicate agreement must not read locked test intervals")
            subwindows, eligible_bases = core_subwindows(training.intervals, SEQUENCE_LENGTH)
            if training.intervals != heldout.intervals:
                raise ValueError(f"Interval manifests differ in fold {fold}")
            chromosomes = {row["chromosome"] for row in training.intervals}
            genes, genes_by_chromosome = full_metrics.load_gene_exons(training.gtf_path, chromosomes)
            gene_exons = full_metrics.GeneExonAccumulator(
                genes, genes_by_chromosome, n_tracks=len(group_ids)
            )
            stats_128 = full_metrics.PerTrackStats(len(group_ids))
            for ordinal, (index, offset) in enumerate(subwindows, start=1):
                training_item = training.get_subwindow_item(
                    index,
                    shift_bp=0,
                    crop_offset_bp=offset,
                    crop_length_bp=SEQUENCE_LENGTH,
                )
                heldout_item = heldout.get_subwindow_item(
                    index,
                    shift_bp=0,
                    crop_offset_bp=offset,
                    crop_length_bp=SEQUENCE_LENGTH,
                )
                prediction_1 = training_item["target_1bp"].numpy().astype(np.float64)
                target_1 = heldout_item["target_1bp"].numpy().astype(np.float64)
                prediction_128 = training_item["target_128bp"].numpy().astype(np.float64)
                target_128 = heldout_item["target_128bp"].numpy().astype(np.float64)
                log_prediction_128 = np.log1p(np.maximum(prediction_128, 0.0))
                log_target_128 = np.log1p(np.maximum(target_128, 0.0))
                stats_128.update(log_prediction_128, log_target_128)
                gene_exons.update(
                    str(training_item["interval_chromosome"]),
                    int(training_item["interval_start"]),
                    prediction_1,
                    target_1,
                )
                print(
                    f"fold={fold}\tsubwindow={ordinal}/{len(subwindows)}\t"
                    f"{training_item['interval_chromosome']}:"
                    f"{int(training_item['interval_start'])}-{int(training_item['interval_end'])}",
                    flush=True,
                )
            gene_metrics = gene_exons.metrics(log1p=True)
            gene_values, coverage_values, primary_values = per_track_record(
                group_ids, gene_metrics, stats_128
            )
            for group_id in group_ids:
                source_rows.append(
                    {
                        "fold": fold,
                        "group_id": group_id,
                        "gene_exon_pearson_log1p": gene_values[group_id],
                        "pearson_128bp_log1p": coverage_values[group_id],
                        "primary_biological_score": primary_values[group_id],
                        "comparison": "training_aggregate_vs_heldout_biological_unit",
                    }
                )
            fold_rows.append(
                {
                    "fold": fold,
                    "validation_subwindows": len(subwindows),
                    "eligible_aligned_core_bases": eligible_bases,
                    "evaluated_bases": len(subwindows) * SEQUENCE_LENGTH,
                    "gene_exon_pearson_log1p": finite_mean(list(gene_values.values())),
                    "pearson_128bp_log1p": finite_mean(list(coverage_values.values())),
                    "primary_biological_score": finite_mean(list(primary_values.values())),
                    "comparison": "training_aggregate_vs_heldout_biological_unit",
                }
            )
        finally:
            training.close()
            heldout.close()

    summary_rows = []
    for metric in (
        "gene_exon_pearson_log1p",
        "pearson_128bp_log1p",
        "primary_biological_score",
    ):
        values = [float(row[metric]) for row in fold_rows]
        summary_rows.append(
            {
                "metric": metric,
                "equal_weight_fold_mean": finite_mean(values),
                "fold_sd": float(np.std(values, ddof=1)),
                "fold_min": min(values),
                "fold_max": max(values),
                "primary_unit": "five predefined genomic folds",
                "interpretation": (
                    "training-aggregate to held-out-biological-unit agreement; "
                    "a measurement-reliability reference, not external validation "
                    "or a formal performance ceiling"
                ),
            }
        )

    write_tsv(
        output / "p13_replicate_agreement_per_group.tsv",
        source_rows,
        [
            "fold",
            "group_id",
            "gene_exon_pearson_log1p",
            "pearson_128bp_log1p",
            "primary_biological_score",
            "comparison",
        ],
    )
    write_tsv(
        output / "p13_replicate_agreement_by_fold.tsv",
        fold_rows,
        [
            "fold",
            "validation_subwindows",
            "eligible_aligned_core_bases",
            "evaluated_bases",
            "gene_exon_pearson_log1p",
            "pearson_128bp_log1p",
            "primary_biological_score",
            "comparison",
        ],
    )
    write_tsv(
        output / "p13_replicate_agreement_summary.tsv",
        summary_rows,
        [
            "metric",
            "equal_weight_fold_mean",
            "fold_sd",
            "fold_min",
            "fold_max",
            "primary_unit",
            "interpretation",
        ],
    )
    audit = {
        "analysis": "P13 training-aggregate versus held-out biological-unit agreement",
        "scope": "read-only P11/P12 development-fold analysis",
        "prohibited": [
            "No checkpoint was loaded.",
            "No GPU was used.",
            "No training or data generation was performed.",
            "No locked final-test intervals or signals were opened.",
        ],
        "input_sha256": {
            str(path.relative_to(REPO_ROOT)): sha256(path)
            for path in (
                primary_index_path,
                training_track_path,
                training_group_path,
                heldout_track_path,
                heldout_group_path,
            )
        },
        "interval_root": str(interval_root.relative_to(REPO_ROOT)),
        "interval_sha256": {
            str((interval_root / f"fold_{fold}" / "valid.tsv").relative_to(REPO_ROOT)):
            sha256(interval_root / f"fold_{fold}" / "valid.tsv")
            for fold in range(1, 6)
        },
        "primary_group_count": len(training_indices),
        "folds": [1, 2, 3, 4, 5],
        "sequence_length": SEQUENCE_LENGTH,
        "metric_definition": (
            "Per-group primary biological score = 0.5 * gene-exon log1p Pearson "
            "+ 0.5 * 128-bp log1p Pearson."
        ),
        "interpretation": (
            "This quantifies agreement between a training-side aggregate and a "
            "held-out biological-unit aggregate. It provides a reliability "
            "reference for model-to-held-out evaluation but is not an absolute "
            "upper bound, an independent replication study, or external validation."
        ),
    }
    (output / "p13_replicate_agreement_audit.json").write_text(
        json.dumps(audit, indent=2, sort_keys=True) + "\n"
    )
    for path in sorted(output.iterdir()):
        if path.is_file():
            print(f"output\t{path.resolve().relative_to(REPO_ROOT)}\t{sha256(path)}")


if __name__ == "__main__":
    main()
