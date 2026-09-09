#!/usr/bin/env python3
"""Summarize held-out per-track metrics by primary and supplementary role."""

from __future__ import annotations

import csv
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
ROOT = REPO_ROOT / "alphagenome_custom/metadata/v2/replicate_holdout_v1"
RECORDS = ROOT / "p12_replicate_holdout_job_records.tsv"
INDICES = ROOT / "heldout_track_indices.tsv"
OUTPUT = ROOT / "p12_heldout_role_metrics.tsv"
METRICS = (
    "primary_biological_score",
    "gene_exon_coverage_pearson_log1p",
    "per_track_pearson_128bp_log1p",
    "spearman_128bp",
    "top1_mse_128bp",
    "top1_calibration_ratio_128bp",
)


def read(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> None:
    index_rows = read(INDICES)
    role_by_group = {
        row["group_id"]: "primary" if row["evaluation_primary"] == "True" else "supplementary"
        for row in index_rows
    }
    output_rows = []
    for record in read(RECORDS):
        validation = json.loads((REPO_ROOT / record["heldout_validation_path"]).read_text())
        full = validation["full_metrics"]
        values = {
            "gene_exon_coverage_pearson_log1p": full["gene_exon_coverage"]["per_track_pearson"],
            "per_track_pearson_128bp_log1p": full["log1p_128bp_sum"]["per_track_pearson"],
            "spearman_128bp": full["distribution_128bp"]["per_track_spearman"],
            "top1_mse_128bp": full["distribution_128bp"]["per_track_top1_mse"],
            "top1_calibration_ratio_128bp": full["distribution_128bp"]["per_track_top1_calibration_ratio"],
        }
        values["primary_biological_score"] = {
            group_id: 0.5 * (
                float(values["gene_exon_coverage_pearson_log1p"][group_id])
                + float(values["per_track_pearson_128bp_log1p"][group_id])
            )
            for group_id in values["gene_exon_coverage_pearson_log1p"]
        }
        for metric in METRICS:
            for group_id, value in values[metric].items():
                output_rows.append({
                    "job_id": record["job_id"],
                    "configuration": record["configuration"],
                    "fold": record["fold"],
                    "seed": record["seed"],
                    "role": role_by_group[group_id],
                    "track_id": group_id,
                    "metric": metric,
                    "value": value,
                })
    with OUTPUT.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(output_rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(output_rows)
    print(f"rows={len(output_rows)}")


if __name__ == "__main__":
    main()
