#!/usr/bin/env python3
"""Audit saved P12 primary-held-out metrics without opening signal data.

The P12 role-level table retains metric values for every configuration, fold,
seed and held-out group. This read-only P13 analysis summarizes Model B's
within-collection biological-unit-held-out outcomes by registered metadata
strata. Group distributions are descriptive: neither groups nor fold-seed
records are treated as independent biological replicates.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from collections import defaultdict
from pathlib import Path
from statistics import mean, median, stdev
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[1]
ROOT = REPO_ROOT / "alphagenome_custom/metadata/v2/replicate_holdout_v1"
REVIEW_PATH = REPO_ROOT / "alphagenome_custom/metadata/v2/audits/P12/review.json"
ROLE_METRICS_PATH = ROOT / "p12_heldout_role_metrics.tsv"
PRIMARY_INDICES_PATH = ROOT / "heldout_primary_track_indices.tsv"
PRIMARY_GROUPS_PATH = ROOT / "heldout_primary_group_manifest.tsv"

MODEL = "B_paper"
METRICS = (
    "primary_biological_score",
    "gene_exon_coverage_pearson_log1p",
    "per_track_pearson_128bp_log1p",
    "spearman_128bp",
    "top1_mse_128bp",
    "top1_calibration_ratio_128bp",
)
STRATA = (
    "sra_study_accession",
    "development_stage",
    "tissue_or_cell_type",
    "condition",
)
MISSING = "[not annotated]"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=REPO_ROOT / "results/v2_p13_heldout_descriptives",
        help="Ignored directory for P13 descriptive source data and audit records.",
    )
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise ValueError(f"Empty TSV: {path}")
    return rows


def write_tsv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def normalized_metadata(row: dict[str, str], field: str) -> str:
    value = row.get(field, "").strip()
    return value if value else MISSING


def finite(values: list[float]) -> list[float]:
    result = [value for value in values if math.isfinite(value)]
    if not result:
        raise RuntimeError("No finite metric values")
    return result


def quartiles(values: list[float]) -> tuple[float, float]:
    ordered = sorted(finite(values))
    if len(ordered) == 1:
        return ordered[0], ordered[0]
    lower = ordered[: len(ordered) // 2]
    upper = ordered[(len(ordered) + 1) // 2 :]
    return median(lower), median(upper)


def main() -> None:
    args = parse_args()
    review = json.loads(REVIEW_PATH.read_text())
    if review.get("phase") != "P12" or review.get("status") != "PASS":
        raise RuntimeError("P13 descriptive audit requires a passing P12 review")

    index_rows = [
        row
        for row in read_tsv(PRIMARY_INDICES_PATH)
        if row["evaluation_primary"].strip().lower() == "true"
    ]
    primary_groups = [row["group_id"] for row in index_rows]
    if len(primary_groups) != 57 or len(set(primary_groups)) != 57:
        raise RuntimeError("Expected exactly 57 primary evaluated groups")
    primary_group_set = set(primary_groups)

    metadata_rows = read_tsv(PRIMARY_GROUPS_PATH)
    metadata = {row["group_id"]: row for row in metadata_rows if row["group_id"] in primary_group_set}
    if set(metadata) != primary_group_set:
        raise RuntimeError("Primary group metadata does not exactly match P12 evaluated groups")

    grouped: dict[tuple[str, str], list[float]] = defaultdict(list)
    metric_rows = [
        row
        for row in read_tsv(ROLE_METRICS_PATH)
        if row["configuration"] == MODEL
        and row["role"] == "primary"
        and row["track_id"] in primary_group_set
        and row["metric"] in METRICS
    ]
    expected = len(primary_groups) * 15 * len(METRICS)
    if len(metric_rows) != expected:
        raise RuntimeError(f"Expected {expected} Model-B primary metric rows, found {len(metric_rows)}")
    for row in metric_rows:
        grouped[(row["track_id"], row["metric"])].append(float(row["value"]))
    if set(grouped) != {(group_id, metric) for group_id in primary_groups for metric in METRICS}:
        raise RuntimeError("Incomplete group-metric matrix")
    if any(len(values) != 15 for values in grouped.values()):
        raise RuntimeError("Each group-metric summary must contain five folds and three seeds")

    per_group_rows: list[dict[str, Any]] = []
    per_group_values: dict[tuple[str, str], float] = {}
    for group_id in primary_groups:
        row: dict[str, Any] = {
            "group_id": group_id,
            "configuration": MODEL,
            "role": "primary",
            "fold_seed_records": 15,
            "scope": "within-collection biological-unit holdout",
            "group_value_note": "Mean across five genomic folds and three nested seeds; descriptive only.",
        }
        for stratum in STRATA:
            row[stratum] = normalized_metadata(metadata[group_id], stratum)
        for metric in METRICS:
            values = finite(grouped[(group_id, metric)])
            value_mean = mean(values)
            per_group_values[(group_id, metric)] = value_mean
            row[f"mean_{metric}"] = value_mean
            row[f"sd_{metric}_across_fold_seed"] = stdev(values)
        calibration = per_group_values[(group_id, "top1_calibration_ratio_128bp")]
        row["abs_log2_top1_calibration_deviation"] = abs(math.log2(calibration)) if calibration > 0 else math.nan
        per_group_rows.append(row)

    stratum_rows: list[dict[str, Any]] = []
    for stratum in STRATA:
        memberships: dict[str, list[str]] = defaultdict(list)
        for group_id in primary_groups:
            memberships[normalized_metadata(metadata[group_id], stratum)].append(group_id)
        for category, group_ids in sorted(memberships.items()):
            for metric in (*METRICS, "abs_log2_top1_calibration_deviation"):
                if metric == "abs_log2_top1_calibration_deviation":
                    values = [
                        next(row[metric] for row in per_group_rows if row["group_id"] == group_id)
                        for group_id in group_ids
                    ]
                else:
                    values = [per_group_values[(group_id, metric)] for group_id in group_ids]
                values = finite(values)
                q1, q3 = quartiles(values)
                stratum_rows.append(
                    {
                        "stratum_field": stratum,
                        "stratum_value": category,
                        "metric": metric,
                        "group_count": len(group_ids),
                        "fold_seed_record_count": len(group_ids) * 15,
                        "mean_of_group_means": mean(values),
                        "median_of_group_means": median(values),
                        "q1_of_group_means": q1,
                        "q3_of_group_means": q3,
                        "minimum_group_mean": min(values),
                        "maximum_group_mean": max(values),
                        "interpretation": "Descriptive metadata stratum; groups are not independent inferential replicates.",
                    }
                )

    output = args.output_dir
    output.mkdir(parents=True, exist_ok=True)
    per_group_fields = [
        "group_id", "configuration", "role", "fold_seed_records", "scope", *STRATA,
        *[field for metric in METRICS for field in (f"mean_{metric}", f"sd_{metric}_across_fold_seed")],
        "abs_log2_top1_calibration_deviation", "group_value_note",
    ]
    stratum_fields = [
        "stratum_field", "stratum_value", "metric", "group_count", "fold_seed_record_count",
        "mean_of_group_means", "median_of_group_means", "q1_of_group_means",
        "q3_of_group_means", "minimum_group_mean", "maximum_group_mean", "interpretation",
    ]
    write_tsv(output / "p13_b_paper_primary_group_descriptives.tsv", per_group_rows, per_group_fields)
    write_tsv(output / "p13_b_paper_primary_metadata_strata.tsv", stratum_rows, stratum_fields)
    audit = {
        "analysis": "P13 saved-metric descriptive audit",
        "scope": "within-collection biological-unit-held-out Model B outcomes",
        "prohibited": [
            "No BigWig, FASTA, GTF, locked test, model checkpoint or GPU was opened.",
            "No group, study, stage, tissue or condition stratum is used as an independent inferential sample.",
            "No observed-signal-decile analysis is claimed because P12 did not retain per-bin target intensity in this source table.",
        ],
        "inputs": {
            str(path.relative_to(REPO_ROOT)): sha256(path)
            for path in (REVIEW_PATH, ROLE_METRICS_PATH, PRIMARY_INDICES_PATH, PRIMARY_GROUPS_PATH)
        },
        "counts": {"primary_groups": 57, "fold_seed_records_per_group_metric": 15, "metrics": len(METRICS)},
        "top1_calibration_definition": (
            "Mean prediction divided by mean observed held-out target among 128-bp bins at or above the observed per-track 99th percentile. "
            "A ratio of one indicates matching mean amplitude over that observed top-1% subset."
        ),
        "statistics": (
            "Group-level means summarize fold-seed outcomes descriptively. The five-fold paired configuration estimand remains the formal comparative analysis."
        ),
    }
    (output / "p13_heldout_descriptives_audit.json").write_text(json.dumps(audit, indent=2, sort_keys=True) + "\n")
    for path in sorted(output.iterdir()):
        print(f"output\t{path.relative_to(REPO_ROOT)}\t{sha256(path)}")


if __name__ == "__main__":
    main()
