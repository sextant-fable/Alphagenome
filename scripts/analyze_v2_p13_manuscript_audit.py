#!/usr/bin/env python3
"""Create a manuscript-facing audit from completed P11/P12 artifacts.

This analysis is deliberately read-only. It never opens the locked final-test
definition or its result records; all quantitative outputs are derived from the
P11 replicate-holdout manifests and the completed P12 metrics already used for
the within-collection biological-unit evaluation.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from collections import defaultdict
from pathlib import Path
from statistics import mean, stdev
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[1]
HOLDOUT_ROOT = REPO_ROOT / "alphagenome_custom/metadata/v2/replicate_holdout_v1"
METRICS_PATH = HOLDOUT_ROOT / "p12_replicate_holdout_metrics.tsv"
RAW_EFFECTS_PATH = HOLDOUT_ROOT / "p12_replicate_holdout_raw_paired_effects.tsv"
EFFECTS_PATH = HOLDOUT_ROOT / "p12_replicate_holdout_paired_effects.tsv"
P11_SPEC_PATH = HOLDOUT_ROOT / "p11_replicate_holdout_spec.json"
P12_SPEC_PATH = REPO_ROOT / "alphagenome_custom/metadata/v2/p12_replicate_holdout_spec.json"

PRIMARY_METRIC = "primary_biological_score"
PRIMARY_FORMULA = (
    "0.5 * mean_per_track_gene_exon_coverage_pearson_log1p + "
    "0.5 * mean_per_track_pearson_128bp_log1p"
)

COMPARISONS = {
    "B_vs_A_paper": {
        "display": "Model B minus Model A",
        "sign": 1.0,
        "claim": "adaptation versus frozen-trunk baseline",
    },
    "B_vs_C_paper": {
        "display": "Model B minus protocol-matched Model C",
        "sign": 1.0,
        "claim": "adaptation versus protocol-matched from-scratch control",
    },
    "C_size_matched_vs_B": {
        "display": "Model B minus size-matched Model C",
        "sign": -1.0,
        "claim": "adaptation versus trainable-parameter-matched from-scratch control",
    },
    "B_no_lora_vs_B": {
        "display": "Full Model B minus no-LoRA configuration",
        "sign": -1.0,
        "claim": "registered LoRA-removal contrast",
    },
    "B_lora_only_vs_B": {
        "display": "Full Model B minus LoRA-only configuration",
        "sign": -1.0,
        "claim": "registered C. elegans embedding contrast conditional on LoRA",
    },
    "B_1bp_head_only_vs_B": {
        "display": "Full Model B minus 1-bp-head-only configuration",
        "sign": -1.0,
        "claim": "registered learned 128-bp-head contrast",
    },
}


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


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=REPO_ROOT / "results/v2_p13_manuscript_evidence_audit",
        help="Ignored output directory for P13 manuscript-facing audit files.",
    )
    return parser.parse_args()


def folded(rows: list[dict[str, str]], comparison: str) -> dict[int, list[float]]:
    values: dict[int, list[float]] = defaultdict(list)
    sign = COMPARISONS[comparison]["sign"]
    for row in rows:
        if (
            row["comparison"] == comparison
            and row["scope"] == "heldout"
            and row["metric"] == PRIMARY_METRIC
        ):
            values[int(row["fold"])].append(sign * float(row["difference"]))
    if sorted(values) != [1, 2, 3, 4, 5] or any(len(seed_values) != 3 for seed_values in values.values()):
        raise ValueError(f"Unexpected five-fold, three-seed P12 pairing for {comparison}: {values}")
    return values


def main() -> None:
    args = parse_args()
    output = args.output_dir
    output.mkdir(parents=True, exist_ok=True)

    for required in (METRICS_PATH, RAW_EFFECTS_PATH, EFFECTS_PATH, P11_SPEC_PATH, P12_SPEC_PATH):
        if not required.is_file():
            raise FileNotFoundError(required)

    metrics = read_tsv(METRICS_PATH)
    raw_effects = read_tsv(RAW_EFFECTS_PATH)
    effects = read_tsv(EFFECTS_PATH)
    p11_spec = json.loads(P11_SPEC_PATH.read_text())
    p12_spec = json.loads(P12_SPEC_PATH.read_text())

    fold_rows: list[dict[str, Any]] = []
    loo_rows: list[dict[str, Any]] = []
    comparison_rows: list[dict[str, Any]] = []
    for comparison, metadata in COMPARISONS.items():
        by_fold = folded(raw_effects, comparison)
        fold_means = {fold: mean(seed_values) for fold, seed_values in by_fold.items()}
        all_fold_means = list(fold_means.values())
        for fold in sorted(by_fold):
            for seed_index, value in enumerate(by_fold[fold], start=1):
                fold_rows.append(
                    {
                        "comparison": comparison,
                        "display": metadata["display"],
                        "claim": metadata["claim"],
                        "fold": fold,
                        "seed_index_within_fold": seed_index,
                        "b_advantage": value,
                        "unit": "matched fold-seed algorithmic evaluation",
                    }
                )
            leave_one_out_values = [value for heldout_fold, value in fold_means.items() if heldout_fold != fold]
            loo_rows.append(
                {
                    "comparison": comparison,
                    "display": metadata["display"],
                    "excluded_fold": fold,
                    "mean_b_advantage_remaining_four_folds": mean(leave_one_out_values),
                    "unit": "equal-weight mean of four fold means",
                }
            )
        reported = next(
            row
            for row in effects
            if row["comparison"] == comparison
            and row["scope"] == "heldout"
            and row["metric"] == PRIMARY_METRIC
        )
        sign = COMPARISONS[comparison]["sign"]
        comparison_rows.append(
            {
                "comparison": comparison,
                "display": metadata["display"],
                "claim": metadata["claim"],
                "fold_mean_b_advantage": mean(all_fold_means),
                "fold_sd": stdev(all_fold_means),
                "fold_min": min(all_fold_means),
                "fold_max": max(all_fold_means),
                "leave_one_fold_out_min": min(
                    float(row["mean_b_advantage_remaining_four_folds"])
                    for row in loo_rows
                    if row["comparison"] == comparison
                ),
                "leave_one_fold_out_max": max(
                    float(row["mean_b_advantage_remaining_four_folds"])
                    for row in loo_rows
                    if row["comparison"] == comparison
                ),
                "reported_hierarchical_bootstrap_estimate": sign * float(reported["estimate"]),
                "reported_hierarchical_bootstrap_ci_95_low": min(
                    sign * float(reported["ci_95_low"]), sign * float(reported["ci_95_high"])
                ),
                "reported_hierarchical_bootstrap_ci_95_high": max(
                    sign * float(reported["ci_95_low"]), sign * float(reported["ci_95_high"])
                ),
                "bootstrap_resamples": reported["bootstrap_replicates"],
                "ci_method": reported["ci_method"],
                "primary_inferential_unit": "five predefined genomic folds; seeds are nested algorithmic repeats",
                "track_role": "descriptive outcomes only; not independent inferential replicates",
            }
        )

    primary_metrics = [
        row
        for row in metrics
        if row["scope"] == "heldout"
        and row["configuration"] in {"A_paper", "B_paper", "C_paper", "C_size_matched"}
        and row["metric"] == PRIMARY_METRIC
    ]
    by_configuration: dict[str, list[float]] = defaultdict(list)
    for row in primary_metrics:
        by_configuration[row["configuration"]].append(float(row["value"]))
    expected = {"A_paper", "B_paper", "C_paper", "C_size_matched"}
    if set(by_configuration) != expected or any(len(values) != 15 for values in by_configuration.values()):
        raise ValueError("The formal P12 held-out configuration matrix is incomplete")
    configuration_rows = [
        {
            "configuration": configuration,
            "mean_15_fold_seed_scores": mean(values),
            "fold_seed_min": min(values),
            "fold_seed_max": max(values),
            "fold_seed_count": len(values),
            "scope": "within-collection biological-unit holdout",
        }
        for configuration, values in sorted(by_configuration.items())
    ]

    input_paths = (P11_SPEC_PATH, P12_SPEC_PATH, METRICS_PATH, RAW_EFFECTS_PATH, EFFECTS_PATH)
    audit = {
        "audit_name": "P13 manuscript evidence audit",
        "scope": "read-only audit of completed P11/P12 artifacts",
        "prohibited": [
            "No locked final-test block or final-test result record was opened.",
            "No model training, checkpoint selection, data regeneration, or controller-state modification occurred.",
        ],
        "primary_score": {
            "metric": PRIMARY_METRIC,
            "formula": PRIMARY_FORMULA,
            "direction": "higher is better",
        },
        "hierarchy": {
            "source_runs": p11_spec["counts"]["source_runs"],
            "target_groups": p11_spec["counts"]["source_groups"],
            "training_only_groups": p11_spec["counts"]["training_only_tracks"],
            "primary_heldout_groups": p11_spec["counts"]["primary_candidate_tracks"],
            "supplementary_heldout_groups": p11_spec["counts"]["supplementary_tracks"],
            "pending_qc_heldout_groups": p11_spec["counts"]["primary_pending_qc_tracks"],
            "heldout_definition": "A selected biological unit was excluded before training-label aggregation; its held-out aggregate was evaluated through the corresponding fixed output head.",
            "external_boundary": "All units remain within the audited 482-run collection; this is not study- or laboratory-isolated external validation.",
        },
        "statistics_contract": {
            "manuscript_inferential_unit": "five predefined genomic folds",
            "nested_repeats": "three training seeds within each fold",
            "track_status": "57 primary held-out groups are descriptive outcomes, not independent inferential replicates",
            "historical_spec_note": (
                "The frozen P12 spec labels its statistics unit as study source nested within track, fold and seed. "
                "The completed P12 paired-effect implementation and manuscript-facing inference use five genomic folds as the primary blocked units. "
                "Do not edit the frozen P12 spec; disclose the manuscript estimand and resampling hierarchy explicitly."
            ),
        },
        "configuration_contract": {
            "full_model": "Model B: frozen trunk, newly added C. elegans embedding, LoRA, dual-resolution RNA head",
            "lora_only": "B_lora_only: original organism embedding, LoRA, dual-resolution RNA head",
            "no_lora": "B_no_lora: newly added C. elegans embedding, no LoRA, dual-resolution RNA head",
            "head_only": "B_1bp_head_only: newly added C. elegans embedding, LoRA, learned 1-bp head with deterministic 128-bp pooling",
            "boundary": "These are registered configuration contrasts, not mechanistic estimates or evidence of unseen-condition prediction.",
        },
        "p12_matrix": p12_spec["matrix"],
        "input_sha256": {str(path.relative_to(REPO_ROOT)): sha256(path) for path in input_paths},
    }

    write_tsv(
        output / "p13_configuration_summary.tsv",
        configuration_rows,
        ["configuration", "mean_15_fold_seed_scores", "fold_seed_min", "fold_seed_max", "fold_seed_count", "scope"],
    )
    write_tsv(
        output / "p13_fold_seed_effects.tsv",
        fold_rows,
        ["comparison", "display", "claim", "fold", "seed_index_within_fold", "b_advantage", "unit"],
    )
    write_tsv(
        output / "p13_leave_one_fold_out.tsv",
        loo_rows,
        ["comparison", "display", "excluded_fold", "mean_b_advantage_remaining_four_folds", "unit"],
    )
    write_tsv(
        output / "p13_comparison_summary.tsv",
        comparison_rows,
        [
            "comparison",
            "display",
            "claim",
            "fold_mean_b_advantage",
            "fold_sd",
            "fold_min",
            "fold_max",
            "leave_one_fold_out_min",
            "leave_one_fold_out_max",
            "reported_hierarchical_bootstrap_estimate",
            "reported_hierarchical_bootstrap_ci_95_low",
            "reported_hierarchical_bootstrap_ci_95_high",
            "bootstrap_resamples",
            "ci_method",
            "primary_inferential_unit",
            "track_role",
        ],
    )
    (output / "p13_audit.json").write_text(json.dumps(audit, indent=2, sort_keys=True) + "\n")

    for path in output.iterdir():
        if path.is_file():
            print(f"{path.relative_to(REPO_ROOT)}\t{sha256(path)}")


if __name__ == "__main__":
    main()
