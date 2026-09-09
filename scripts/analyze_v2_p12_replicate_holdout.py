#!/usr/bin/env python3
"""Compute preregistered paired effects for the P12 two-layer validation."""

from __future__ import annotations

import csv
from pathlib import Path

from scripts import v2_submission_statistics as statistics


REPO_ROOT = Path(__file__).resolve().parents[1]
ROOT = REPO_ROOT / "alphagenome_custom/metadata/v2/replicate_holdout_v1"
METRICS_PATH = ROOT / "p12_replicate_holdout_metrics.tsv"
RAW_PATH = ROOT / "p12_replicate_holdout_raw_paired_effects.tsv"
SUMMARY_PATH = ROOT / "p12_replicate_holdout_paired_effects.tsv"
FOLDS = [1, 2, 3, 4, 5]
SEEDS = [20260714, 20260715, 20260716]
METRICS = (
    "primary_biological_score",
    "gene_exon_coverage_pearson_log1p",
    "per_track_pearson_128bp_log1p",
    "paper_loss",
    "log1p_mse",
    "spearman_128bp",
    "top1_mse_128bp",
    "top1_calibration_ratio_128bp",
    "gene_body_mse_log1p_1bp",
    "exon_mse_log1p_1bp",
    "local_gradient_mse_log1p_1bp",
)
SCOPES = ("internal", "heldout")
COMPARISONS = (
    ("B_paper", "A_paper", "B_vs_A_paper"),
    ("B_paper", "C_paper", "B_vs_C_paper"),
    ("B_log1p_mse", "A_log1p_mse", "B_vs_A_log1p_mse"),
    ("B_log1p_mse", "C_log1p_mse", "B_vs_C_log1p_mse"),
    ("B_no_augmentation", "B_paper", "B_no_augmentation_vs_B"),
    ("B_no_gene_loss", "B_paper", "B_no_gene_loss_vs_B"),
    ("A_development_pool_mean", "A_paper", "A_development_pool_mean_vs_A"),
    ("B_no_lora", "B_paper", "B_no_lora_vs_B"),
    ("B_lora_only", "B_paper", "B_lora_only_vs_B"),
    ("B_1bp_head_only", "B_paper", "B_1bp_head_only_vs_B"),
    ("C_size_matched", "B_paper", "C_size_matched_vs_B"),
)


def read(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise RuntimeError(f"Refusing to write empty P12 analysis table: {path}")
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    rows = read(METRICS_PATH)
    by_key = {
        (row["configuration"], row["scope"], row["metric"], int(row["fold"]), int(row["seed"])): row
        for row in rows
    }
    raw: list[dict[str, object]] = []
    summary: list[dict[str, object]] = []
    for candidate, reference, comparison in COMPARISONS:
        repeats = [20260714] if candidate in {"B_no_augmentation", "B_no_gene_loss", "A_development_pool_mean"} else SEEDS
        for scope in SCOPES:
            for metric in METRICS:
                pairs = []
                for fold in FOLDS:
                    for seed in repeats:
                        candidate_row = by_key[(candidate, scope, metric, fold, seed)]
                        reference_row = by_key[(reference, scope, metric, fold, seed)]
                        difference = float(candidate_row["value"]) - float(reference_row["value"])
                        pair = {"fold": fold, "seed": seed, "difference": difference}
                        pairs.append(pair)
                        raw.append({
                            "comparison": comparison,
                            "scope": scope,
                            "metric": metric,
                            "fold": fold,
                            "seed": seed,
                            "candidate": candidate,
                            "reference": reference,
                            "candidate_value": candidate_row["value"],
                            "reference_value": reference_row["value"],
                            "difference": difference,
                        })
                matrix = statistics.records_to_balanced_matrix(
                    pairs,
                    value_key="difference",
                    expected_blocks=FOLDS,
                    expected_repeats=repeats,
                )
                plan = statistics.make_hierarchical_resample_plan(
                    n_blocks=len(FOLDS),
                    n_repeats=len(repeats),
                    iterations=10000,
                    seed=20260819 + len(summary),
                )
                estimate = statistics.hierarchical_block_bootstrap(matrix.values, plan, confidence_level=0.95)
                summary.append({
                    "comparison": comparison,
                    "scope": scope,
                    "metric": metric,
                    "candidate": candidate,
                    "reference": reference,
                    "paired_runs": len(pairs),
                    "folds": len(FOLDS),
                    "seeds": len(repeats),
                    "estimate": estimate.estimate,
                    "ci_95_low": estimate.ci_lower,
                    "ci_95_high": estimate.ci_upper,
                    "bootstrap_replicates": 10000,
                    "ci_method": "fold_first_seed_within_fold_hierarchical_percentile_bootstrap",
                })
    write(RAW_PATH, raw)
    write(SUMMARY_PATH, summary)
    print(f"raw_rows={len(raw)} summary_rows={len(summary)}")


if __name__ == "__main__":
    main()
