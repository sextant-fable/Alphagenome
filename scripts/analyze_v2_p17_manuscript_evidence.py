#!/usr/bin/env python3
"""Summarize completed P17 controls for manuscript source data.

This read-only analysis consumes saved P17 fold summaries. It reports each of
the five predefined genomic-fold contrasts and their equal-weight mean, SD,
range and leave-one-fold-out sensitivity. It does not reopen coverage, training
targets, checkpoints, chromosome X or the locked test.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
from datetime import datetime, timezone
from pathlib import Path
from statistics import fmean, stdev


REPO_ROOT = Path(__file__).resolve().parents[1]
RESULT_ROOT = REPO_ROOT / "results/v2_p17_iv_controlled_matrix"
REVIEW_PATH = REPO_ROOT / "alphagenome_custom/metadata/v2/audits/P17/review.json"
OUTPUT = REPO_ROOT / "results/v2_p17_manuscript_evidence_20260823T032000Z"
CONFIGURATIONS = ("B_iv_dual", "B_128bp_only", "Basenji2_style")
COMPARISONS = ("B_iv_dual_minus_B_128bp_only", "B_iv_dual_minus_Basenji2_style")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def finite(value: str) -> float:
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f"Non-finite P17 value: {value}")
    return result


def summarize(values: list[float]) -> dict[str, float]:
    if len(values) != 5:
        raise RuntimeError("P17 manuscript summaries require exactly five predefined folds")
    return {
        "mean": fmean(values),
        "sd": stdev(values),
        "minimum": min(values),
        "maximum": max(values),
    }


def main() -> None:
    if OUTPUT.exists():
        raise FileExistsError(f"P17 manuscript evidence output is immutable: {OUTPUT}")
    by_fold_path = RESULT_ROOT / "p17_iv_controlled_by_fold.tsv"
    paired_path = RESULT_ROOT / "p17_iv_controlled_paired_effects.tsv"
    if not by_fold_path.is_file() or not paired_path.is_file() or not REVIEW_PATH.is_file():
        raise FileNotFoundError("P17 fold summaries or review are missing")
    review = json.loads(REVIEW_PATH.read_text(encoding="utf-8"))
    if review.get("status") != "PASS":
        raise RuntimeError("P17 manuscript analysis requires the R17 pass record")
    by_fold = read_tsv(by_fold_path)
    paired = read_tsv(paired_path)
    expected_keys = {(configuration, str(fold)) for configuration in CONFIGURATIONS for fold in range(1, 6)}
    observed_keys = {(row.get("configuration", ""), row.get("fold", "")) for row in by_fold}
    if observed_keys != expected_keys or len(by_fold) != 15:
        raise RuntimeError("P17 fold table does not contain the registered 3-by-5 matrix")
    if {(row.get("comparison", ""), row.get("fold", "")) for row in paired} != {(comparison, str(fold)) for comparison in COMPARISONS for fold in range(1, 6)}:
        raise RuntimeError("P17 paired table does not contain the registered comparisons")

    configuration_rows: list[dict[str, object]] = []
    for configuration in CONFIGURATIONS:
        rows = sorted((row for row in by_fold if row["configuration"] == configuration), key=lambda row: int(row["fold"]))
        primary = [finite(row["primary_biological_score"]) for row in rows]
        gene = [finite(row["gene_exon_coverage_pearson_log1p"]) for row in rows]
        bins = [finite(row["per_track_pearson_128bp_log1p"]) for row in rows]
        for metric, values in (("primary_biological_score", primary), ("gene_exon_coverage_pearson_log1p", gene), ("per_track_pearson_128bp_log1p", bins)):
            summary = summarize(values)
            configuration_rows.append({"configuration": configuration, "metric": metric, **summary, "folds": "1,2,3,4,5", "nested_seeds": 3})

    comparison_rows: list[dict[str, object]] = []
    loo_rows: list[dict[str, object]] = []
    for comparison in COMPARISONS:
        rows = sorted((row for row in paired if row["comparison"] == comparison), key=lambda row: int(row["fold"]))
        for metric in ("primary_biological_score", "gene_exon_coverage_pearson_log1p", "per_track_pearson_128bp_log1p"):
            column = f"difference_{metric}"
            values = [finite(row[column]) for row in rows]
            summary = summarize(values)
            comparison_rows.append({"comparison": comparison, "metric": metric, **summary, "folds": "1,2,3,4,5", "inferential_unit": "predefined genomic fold", "nested_seeds": 3})
            for dropped_index, dropped in enumerate(range(1, 6)):
                retained = [value for index, value in enumerate(values, start=1) if index != dropped]
                loo_rows.append({"comparison": comparison, "metric": metric, "dropped_fold": dropped, "mean_without_dropped_fold": fmean(retained)})

    OUTPUT.mkdir(parents=True)
    config_path = OUTPUT / "p17_configuration_summary.tsv"
    paired_path_out = OUTPUT / "p17_paired_effect_summary.tsv"
    loo_path = OUTPUT / "p17_leave_one_fold_out.tsv"
    source_path = OUTPUT / "p17_figure_source_data.tsv"
    write_tsv(config_path, configuration_rows, list(configuration_rows[0]))
    write_tsv(paired_path_out, comparison_rows, list(comparison_rows[0]))
    write_tsv(loo_path, loo_rows, list(loo_rows[0]))
    source_rows = [{"table": "fold_configuration", **row} for row in by_fold] + [{"table": "fold_paired_effect", **row} for row in paired]
    source_fields = sorted({field for row in source_rows for field in row})
    write_tsv(source_path, source_rows, source_fields)
    claim = {
        "schema_version": 1,
        "created_at": utc_now(),
        "scope": "read-only P17 fold-summary analysis; no training, coverage, model checkpoint, chromosome-X or locked-test access",
        "input": {
            "by_fold": {"path": str(by_fold_path.relative_to(REPO_ROOT)), "sha256": sha256(by_fold_path)},
            "paired": {"path": str(paired_path.relative_to(REPO_ROOT)), "sha256": sha256(paired_path)},
            "review": {"path": str(REVIEW_PATH.relative_to(REPO_ROOT)), "sha256": sha256(REVIEW_PATH)},
        },
        "primary_interpretation": "B_iv_dual_minus_Basenji2_style is a task-matched architecture comparison under the strict I-V contract.",
        "resolution_boundary": "B_iv_dual_minus_B_128bp_only quantifies the registered configuration contrast only. Because the 128-bp-only head derives 1-bp values solely for metric compatibility, it does not establish a one-base biological endpoint or base-resolution utility.",
        "statistics": "Each configuration and paired effect is summarized as an equal-weight mean over five predefined genomic folds, with three nested algorithmic seeds per fold. Tracks are not independent inferential units.",
        "outputs": {
            "configuration_summary": {"path": str(config_path.relative_to(REPO_ROOT)), "sha256": sha256(config_path)},
            "paired_effect_summary": {"path": str(paired_path_out.relative_to(REPO_ROOT)), "sha256": sha256(paired_path_out)},
            "leave_one_fold_out": {"path": str(loo_path.relative_to(REPO_ROOT)), "sha256": sha256(loo_path)},
            "figure_source_data": {"path": str(source_path.relative_to(REPO_ROOT)), "sha256": sha256(source_path)},
        },
    }
    manifest = OUTPUT / "p17_manuscript_evidence_manifest.json"
    manifest.write_text(json.dumps(claim, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"status": "completed", "output": str(OUTPUT.relative_to(REPO_ROOT)), "manifest": str(manifest.relative_to(REPO_ROOT))}, indent=2))


if __name__ == "__main__":
    main()
