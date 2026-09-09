#!/usr/bin/env python3
"""Run fail-closed evidence reviews for completed v2 workflow phases."""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[1]
LEGACY_DIR = REPO_ROOT / "alphagenome_custom/metadata/legacy_v1"
AUDIT_ROOT = REPO_ROOT / "alphagenome_custom/metadata/v2/audits"
P10_IMPLEMENTATION_PATHS = (
    "scripts/run_dpy27_internal_application.py",
    "scripts/v2_biological_validation.py",
    "scripts/plot_dpy27_internal_application.py",
    "scripts/train_v2_model.py",
    "scripts/v2_training_components.py",
    "scripts/v2_bigwig_dataset.py",
)
P10_SOURCE_DATA_CONTRACT = {
    "panels": [
        "a_gene_contrasts",
        "b_fold_x_minus_autosomes",
        "c_fold_agreement",
        "c_fold_primary_summary",
    ],
    "panel_endpoints": {
        "a_gene_contrasts": ["gene_log1p_contrast"],
        "b_fold_x_minus_autosomes": [
            "observed_median_contrast_x_minus_autosomes",
            "predicted_median_contrast_x_minus_autosomes",
        ],
        "c_fold_agreement": [
            "x_predicted_observed_direction_concordance",
            "predicted_observed_gene_spearman",
        ],
        "c_fold_primary_summary": [
            "predicted_median_contrast_x_minus_autosomes",
            "observed_median_contrast_x_minus_autosomes",
            "x_predicted_observed_direction_concordance",
            "predicted_observed_gene_spearman",
        ],
    },
}


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def file_hashes(path: Path, chunk_size: int = 8 * 1024 * 1024) -> tuple[str, str]:
    sha_digest = hashlib.sha256()
    md5_digest = hashlib.md5(usedforsecurity=False)
    with path.open("rb") as handle:
        while chunk := handle.read(chunk_size):
            sha_digest.update(chunk)
            md5_digest.update(chunk)
    return sha_digest.hexdigest(), md5_digest.hexdigest()


def sha256(path: Path) -> str:
    return file_hashes(path)[0]


def repository_file(value: object) -> Path:
    raw = Path(str(value))
    path = raw.resolve() if raw.is_absolute() else (REPO_ROOT / raw).resolve()
    path.relative_to(REPO_ROOT.resolve())
    if not path.is_file():
        raise FileNotFoundError(path)
    return path


def check(check_id: str, passed: bool, evidence: str) -> dict[str, Any]:
    return {
        "check_id": check_id,
        "status": "PASS" if passed else "FAIL",
        "evidence": evidence,
    }


def _p10_implementation_sha256(repo_root: Path | None = None) -> dict[str, str]:
    repo_root = REPO_ROOT if repo_root is None else repo_root
    return {
        relative: sha256(repo_root / relative)
        for relative in P10_IMPLEMENTATION_PATHS
    }


def _p10_implementation_manifest_valid(manifest: object) -> bool:
    if not isinstance(manifest, dict) or set(manifest) != set(P10_IMPLEMENTATION_PATHS):
        return False
    return all(
        isinstance(value, str)
        and len(value) == 64
        and all(character in "0123456789abcdefABCDEF" for character in value)
        for value in manifest.values()
    )


def _p10_execution_provenance_ok(execution: dict[str, Any]) -> bool:
    python_executable = Path(str(execution.get("python_executable", "")))
    return (
        bool(execution.get("git_commit"))
        and bool(execution.get("hostname"))
        and execution.get("working_directory") == str(REPO_ROOT)
        and python_executable.is_absolute()
        and bool(python_executable.name)
        and execution.get("controller_invocation")
        == "python -m scripts.v2_phase_controller run --phase P10"
        and execution.get("registered_phase_command")
        == "python -m scripts.run_dpy27_internal_application"
    )


def _p10_close(actual: object, expected: object) -> bool:
    try:
        return math.isclose(
            float(actual), float(expected), rel_tol=1e-12, abs_tol=1e-12
        )
    except (TypeError, ValueError):
        return False


def _p10_summary_recomputation_ok(
    gene_rows: list[dict[str, Any]],
    fold_rows: list[dict[str, Any]],
    summary_rows: list[dict[str, Any]],
    statistics: dict[str, Any],
) -> bool:
    """Independently recompute all fold endpoints and deterministic bootstrap CIs."""

    try:
        from scripts import v2_biological_validation as biological_validation

        endpoints = tuple(statistics["endpoints"])
        _, expected_folds, expected_summaries = biological_validation.fold_primary_summary(
            gene_rows,
            endpoints=endpoints,
            bootstrap_replicates=int(statistics["bootstrap_replicates"]),
            bootstrap_seed=int(statistics["bootstrap_seed"]),
        )
        if len(fold_rows) != 5 or len(summary_rows) != len(endpoints):
            return False
        recorded_folds = {int(row["fold"]): row for row in fold_rows}
        recorded_summaries = {row["endpoint"]: row for row in summary_rows}
        if len(recorded_folds) != 5 or set(recorded_summaries) != set(endpoints):
            return False
        for expected in expected_folds:
            recorded = recorded_folds[int(expected["fold"])]
            if int(recorded["n_seeds"]) != 3 or any(
                not _p10_close(recorded[endpoint], expected[endpoint])
                for endpoint in endpoints
            ):
                return False
        for expected in expected_summaries:
            recorded = recorded_summaries[str(expected["endpoint"])]
            if not (
                int(recorded["n_folds"]) == int(expected["n_folds"]) == 5
                and int(recorded["n_runs"]) == int(expected["n_runs"]) == 15
                and recorded["inference_unit"] == expected["inference_unit"]
                and int(recorded["bootstrap_replicates"])
                == int(expected["bootstrap_replicates"])
                and int(recorded["bootstrap_seed"])
                == int(expected["bootstrap_seed"])
                and all(
                    _p10_close(recorded[field], expected[field])
                    for field in ("estimate", "ci_95_low", "ci_95_high")
                )
            ):
                return False
        return True
    except Exception:
        return False


def _p10_source_data_ok(
    source_rows: list[dict[str, Any]],
    fold_gene_rows: list[dict[str, Any]],
    fold_rows: list[dict[str, Any]],
    summary_rows: list[dict[str, Any]],
    contract: dict[str, Any],
) -> bool:
    """Rebuild the expected four-panel Source Data mapping and compare values."""

    try:
        if contract != P10_SOURCE_DATA_CONTRACT:
            return False
        panels = set(contract["panels"])
        if {row["panel"] for row in source_rows} != panels:
            return False
        for panel, endpoints in contract["panel_endpoints"].items():
            if {row["endpoint"] for row in source_rows if row["panel"] == panel} != set(
                endpoints
            ):
                return False
        indexed: dict[tuple[str, str, str, str], dict[str, Any]] = {}
        for row in source_rows:
            key = (
                str(row["panel"]),
                str(row["fold"]),
                str(row["gene_id"]),
                str(row["endpoint"]),
            )
            if key in indexed:
                return False
            indexed[key] = row

        expected_keys: set[tuple[str, str, str, str]] = set()
        for row in fold_gene_rows:
            complete = row["complete_exon_coverage"] is True or str(
                row["complete_exon_coverage"]
            ) == "True"
            if not complete:
                continue
            key = (
                "a_gene_contrasts",
                str(row["fold"]),
                str(row["gene_id"]),
                "gene_log1p_contrast",
            )
            expected_keys.add(key)
            actual = indexed.get(key)
            if not actual or not (
                actual["seed"] == "seed_mean"
                and actual["gene_name"] == row["gene_name"]
                and actual["chromosome"] == row["chromosome"]
                and _p10_close(
                    actual["observed"],
                    row["observed_log1p_contrast_g0005_minus_g0006"],
                )
                and _p10_close(
                    actual["predicted"],
                    row["predicted_log1p_contrast_g0005_minus_g0006"],
                )
            ):
                return False

        for row in fold_rows:
            for kind in ("observed", "predicted"):
                endpoint = f"{kind}_median_contrast_x_minus_autosomes"
                key = (
                    "b_fold_x_minus_autosomes",
                    str(row["fold"]),
                    "",
                    endpoint,
                )
                expected_keys.add(key)
                actual = indexed.get(key)
                if not actual or not (
                    actual["seed"] == "seed_mean"
                    and _p10_close(actual["estimate"], row[endpoint])
                    and _p10_close(actual[kind], row[endpoint])
                    and str(actual["predicted" if kind == "observed" else "observed"])
                    == ""
                ):
                    return False
            for endpoint in contract["panel_endpoints"]["c_fold_agreement"]:
                key = ("c_fold_agreement", str(row["fold"]), "", endpoint)
                expected_keys.add(key)
                actual = indexed.get(key)
                if not actual or not (
                    actual["seed"] == "seed_mean"
                    and _p10_close(actual["estimate"], row[endpoint])
                ):
                    return False

        for row in summary_rows:
            endpoint = str(row["endpoint"])
            key = ("c_fold_primary_summary", "all_folds", "", endpoint)
            expected_keys.add(key)
            actual = indexed.get(key)
            if not actual or not (
                actual["seed"] == "within_fold_seed_mean"
                and _p10_close(actual["estimate"], row["estimate"])
                and _p10_close(actual["ci_95_low"], row["ci_95_low"])
                and _p10_close(actual["ci_95_high"], row["ci_95_high"])
            ):
                return False
        return set(indexed) == expected_keys
    except Exception:
        return False


def review_p0() -> dict[str, Any]:
    required = [
        "raw_bigwig_inventory.tsv",
        "grouped_bigwig_inventory.tsv",
        "reference_inventory.tsv",
        "split_inventory.tsv",
        "dataset_inventory.tsv",
        "candidate_checkpoint_inventory.tsv",
        "selected_model_inventory.tsv",
        "source_file_inventory.tsv",
        "environment.tsv",
        "freeze_summary.json",
    ]
    missing = [name for name in required if not (LEGACY_DIR / name).is_file()]
    checks = [
        check("R0.01_required_outputs", not missing, f"missing={missing}"),
    ]
    if missing:
        return {"phase": "P0", "status": "FAIL", "checks": checks}

    raw = read_tsv(LEGACY_DIR / "raw_bigwig_inventory.tsv")
    grouped = read_tsv(LEGACY_DIR / "grouped_bigwig_inventory.tsv")
    references = read_tsv(LEGACY_DIR / "reference_inventory.tsv")
    splits = read_tsv(LEGACY_DIR / "split_inventory.tsv")
    datasets = read_tsv(LEGACY_DIR / "dataset_inventory.tsv")
    candidates = read_tsv(LEGACY_DIR / "candidate_checkpoint_inventory.tsv")
    selected = read_tsv(LEGACY_DIR / "selected_model_inventory.tsv")
    sources = read_tsv(LEGACY_DIR / "source_file_inventory.tsv")
    summary = json.loads((LEGACY_DIR / "freeze_summary.json").read_text())

    split_counts = {row["split"]: int(row["interval_count"]) for row in splits}
    dataset_counts = {row["split"]: int(row["example_count"]) for row in datasets}
    selected_roles = {row["role"] for row in selected}
    reference_types = {row["asset_type"] for row in references if row["exists"] == "True"}

    checks.extend(
        [
            check(
                "R0.02_raw_bigwigs",
                len(raw) == 20
                and all(row["exists"] == "True" and row["sha256"] for row in raw),
                f"count={len(raw)}",
            ),
            check(
                "R0.03_grouped_bigwigs",
                len(grouped) == 11
                and all(
                    row["exists"] == "True"
                    and row["sha256"]
                    and row["sample_ids"]
                    for row in grouped
                ),
                f"count={len(grouped)}",
            ),
            check(
                "R0.04_normalization_and_ontology_recorded",
                all(row["normalization_status"] for row in raw)
                and all(row["normalization_status"] for row in grouped)
                and all("ontology_curie" in row for row in grouped),
                "normalization status and ontology columns are present",
            ),
            check(
                "R0.05_reference_assets",
                {"FASTA", "FASTA_INDEX", "GTF"}.issubset(reference_types)
                and all(row["sha256"] for row in references if row["exists"] == "True"),
                f"types={sorted(reference_types)}",
            ),
            check(
                "R0.06_split_counts",
                split_counts
                == {"train": 116, "valid": 39, "test": 33, "all": 188},
                f"counts={split_counts}",
            ),
            check(
                "R0.07_dataset_counts",
                dataset_counts == {"train": 116, "valid": 39, "test": 33},
                f"counts={dataset_counts}; NPZ content not hashed by policy",
            ),
            check(
                "R0.08_candidate_inventory",
                len(candidates) == 171
                and len({row["run_id"] for row in candidates}) == 171
                and all(
                    row["checkpoint_exists"] == "True"
                    and row["checkpoint_sha256"]
                    and row["metrics_exists"] == "True"
                    for row in candidates
                ),
                f"count={len(candidates)} unique={len({row['run_id'] for row in candidates})}",
            ),
            check(
                "R0.09_selected_models",
                {
                    "base_weights",
                    "legacy_tested_adapter",
                    "validation_mse_winner",
                    "representation_utility_winner",
                }
                == selected_roles
                and all(row["exists"] == "True" and row["sha256"] for row in selected),
                f"roles={sorted(selected_roles)}",
            ),
            check(
                "R0.10_test_exposure",
                any(
                    row["role"] == "legacy_tested_adapter"
                    and row["chr_x_test_status"] == "read_once_recorded"
                    for row in selected
                )
                and all(
                    row["chr_x_test_status"] == "not_read_for_latest_selection"
                    for row in selected
                    if row["role"]
                    in {"validation_mse_winner", "representation_utility_winner"}
                ),
                "legacy test exposure and latest validation-only status are explicit",
            ),
            check(
                "R0.11_raw_immutability",
                summary["raw_hashes_unchanged"] is True
                and summary["protected_files_changed_during_freeze"] == [],
                "raw hashes agree before/after and protected paths were unchanged",
            ),
            check(
                "R0.12_source_evidence",
                all(row["exists"] == "True" and row["sha256"] for row in sources),
                f"source_files={len(sources)}",
            ),
            check(
                "R0.13_freeze_document",
                (REPO_ROOT / "docs/legacy_v1_freeze.md").is_file(),
                "docs/legacy_v1_freeze.md",
            ),
        ]
    )
    return {
        "schema_version": 1,
        "phase": "P0",
        "review": "R0",
        "reviewed_at": utc_now(),
        "status": "PASS"
        if all(item["status"] == "PASS" for item in checks)
        else "FAIL",
        "checks": checks,
    }


def review_p1() -> dict[str, Any]:
    manifest_path = REPO_ROOT / "alphagenome_custom/metadata/v2/rna_seq_samples_v2.tsv"
    evidence_path = (
        REPO_ROOT / "alphagenome_custom/metadata/v2/rna_seq_metadata_evidence_v2.tsv"
    )
    summary_path = REPO_ROOT / "alphagenome_custom/metadata/v2/p1_summary.json"
    full_sources_path = (
        REPO_ROOT / "alphagenome_custom/metadata/v2/p3_full_sources.tsv"
    )
    approval_path = REPO_ROOT / "docs/v2_g1_source_reprocessing_request.md"
    required = [
        manifest_path,
        evidence_path,
        summary_path,
        full_sources_path,
        approval_path,
    ]
    missing = [str(path.relative_to(REPO_ROOT)) for path in required if not path.is_file()]
    checks = [check("R1.01_required_outputs", not missing, f"missing={missing}")]
    if missing:
        return {
            "schema_version": 1,
            "phase": "P1",
            "review": "R1",
            "reviewed_at": utc_now(),
            "status": "FAIL",
            "checks": checks,
        }

    rows = read_tsv(manifest_path)
    evidence = read_tsv(evidence_path)
    full_sources = read_tsv(full_sources_path)
    summary = json.loads(summary_path.read_text())
    batch_counts: dict[str, int] = {}
    for row in rows:
        batch_counts[row["batch"]] = batch_counts.get(row["batch"], 0) + 1
    classes = {row["normalization_class"] for row in rows}
    formal = [row for row in rows if row["include_formal_v2"] == "True"]
    required_formal_fields = [
        "assay",
        "study_title",
        "development_stage",
        "tissue_or_cell_type",
        "genotype",
        "strandedness",
        "biological_replicate",
        "signal_unit",
    ]
    formal_complete = all(
        all(row[field] and row[field] != "unknown" for field in required_formal_fields)
        for row in formal
    )
    evidence_uids = {row["sample_uid"] for row in evidence}
    query_statuses = {row["ena_query_status"] for row in rows}

    checks.extend(
        [
            check(
                "R1.02_sample_scope",
                len(rows) == 485
                and len({row["run_accession"] for row in rows}) == 485
                and len({row["sample_uid"] for row in rows}) == 485,
                f"rows={len(rows)} unique_runs={len({row['run_accession'] for row in rows})}",
            ),
            check(
                "R1.03_batch_counts",
                batch_counts
                == {
                    "2026.4.20": 20,
                    "2026.5.18": 26,
                    "sharepoint_2026-07-13_new439": 439,
                },
                f"counts={batch_counts}",
            ),
            check(
                "R1.04_local_integrity",
                all(
                    row["sha256"]
                    and row["sha256_matches_prior_inventory"] == "True"
                    and row["wbcel235_exact"] == "True"
                    and (REPO_ROOT / row["local_path"]).is_file()
                    for row in rows
                )
                and summary["raw_hashes_unchanged_during_p1"] is True,
                "all local paths match prior SHA-256 and WBcel235 evidence; P1 before/after hashes agree",
            ),
            check(
                "R1.05_accession_provenance",
                all(
                    row["run_accession"]
                    and row["experiment_accession"]
                    and row["biosample_accession"]
                    and row["bioproject_accession"]
                    and row["sra_study_accession"]
                    for row in rows
                    if row["ena_query_status"] == "success"
                ),
                f"ena_success={summary['ena_success_count']} failures={len(summary['ena_failure_accessions'])}",
            ),
            check(
                "R1.06_query_failures_explicit",
                query_statuses.issubset({"success", "not_found"})
                and sum(row["ena_query_status"] == "not_found" for row in rows)
                == len(summary["ena_failure_accessions"]),
                f"statuses={sorted(query_statuses)}",
            ),
            check(
                "R1.07_normalization_classification",
                classes.issubset({"A", "B", "C", "D"})
                and all(row["normalization_class"] for row in rows)
                and sum(summary["normalization_class_counts"].values()) == 485,
                f"classes={summary['normalization_class_counts']}",
            ),
            check(
                "R1.08_no_unjustified_formal_inclusion",
                formal_complete
                and all(
                    row["normalization_class"] in {"A", "B"}
                    and row["signal_unit"] != "unknown"
                    for row in formal
                ),
                f"formal_included={len(formal)}",
            ),
            check(
                "R1.09_evidence_coverage",
                len(evidence) > 0
                and all(row["sample_uid"] in evidence_uids for row in rows)
                and all(
                    row["source_type"]
                    and row["source_locator"]
                    and row["evidence_status"]
                    for row in evidence
                ),
                f"evidence_rows={len(evidence)} samples={len(evidence_uids)}",
            ),
            check(
                "R1.10_bigwig_provenance_explicit",
                all(
                    row["bigwig_generator"]
                    and row["bigwig_reference_genome"]
                    and row["bigwig_aligner"]
                    and row["bigwig_generation_command"]
                    and row["signal_unit"]
                    for row in rows
                ),
                "unknown provenance is explicit rather than inferred from bigWig statistics",
            ),
            check(
                "R1.11_conflicts_preserved",
                all("metadata_conflict_flags" in row for row in rows)
                and any(
                    row["run_accession"] in {"SRR941651", "SRR941697"}
                    and "legacy_T38" in row["metadata_conflict_flags"]
                    for row in rows
                ),
                f"conflict_flagged={summary['conflict_flagged_count']}",
            ),
            check(
                "R1.12_g1_packet",
                approval_path.is_file()
                and "APPROVAL REQUIRED - NOT EXECUTED" in approval_path.read_text()
                and summary["fastq_total_bytes"] >= 0,
                f"fastq_bytes={summary['fastq_total_bytes']} approval_packet_recorded",
            ),
            check(
                "R1.13_full_rna_source_manifest",
                len(full_sources) == 482
                and {row["run_accession"] for row in full_sources}
                == {
                    row["run_accession"]
                    for row in rows
                    if row["assay"] == "RNA-Seq"
                }
                and all(
                    len(row["fastq_ftp"].split(";"))
                    == len(row["fastq_md5"].split(";"))
                    == len(row["fastq_bytes"].split(";"))
                    == (2 if row["library_layout"] == "PAIRED" else 1)
                    for row in full_sources
                ),
                "482 RNA-seq runs have per-file URL, MD5, byte count, and layout",
            ),
        ]
    )
    return {
        "schema_version": 1,
        "phase": "P1",
        "review": "R1",
        "reviewed_at": utc_now(),
        "status": "PASS"
        if all(item["status"] == "PASS" for item in checks)
        else "FAIL",
        "checks": checks,
    }


def review_p2() -> dict[str, Any]:
    base = REPO_ROOT / "alphagenome_custom/metadata/v2"
    paths = {
        "samples": base / "rna_seq_samples_v2.tsv",
        "sample_context": base / "rna_seq_sample_context_v2.tsv",
        "groups": base / "rna_seq_groups_v2.tsv",
        "members": base / "rna_seq_group_members_v2.tsv",
        "decisions": base / "rna_seq_sample_decisions_v2.tsv",
        "exclusions": base / "rna_seq_exclusions_v2.tsv",
        "ontology": base / "ontology_mapping_v2.tsv",
        "resolutions": base / "metadata_resolutions_v2.tsv",
        "pair_qc": base / "replicate_pair_qc_v2.tsv",
        "group_qc": base / "replicate_group_qc_v2.tsv",
        "scientific_review": base / "scientific_group_review_v2.tsv",
        "summary": base / "p2_summary.json",
    }
    missing = [name for name, path in paths.items() if not path.is_file()]
    checks = [check("R2.01_required_outputs", not missing, f"missing={missing}")]
    if missing:
        return {
            "schema_version": 1,
            "phase": "P2",
            "review": "R2",
            "reviewed_at": utc_now(),
            "status": "FAIL",
            "checks": checks,
        }

    samples = read_tsv(paths["samples"])
    sample_context = read_tsv(paths["sample_context"])
    groups = read_tsv(paths["groups"])
    members = read_tsv(paths["members"])
    decisions = read_tsv(paths["decisions"])
    exclusions = read_tsv(paths["exclusions"])
    ontology = read_tsv(paths["ontology"])
    resolutions = read_tsv(paths["resolutions"])
    pair_qc = read_tsv(paths["pair_qc"])
    group_qc = read_tsv(paths["group_qc"])
    scientific_review = read_tsv(paths["scientific_review"])
    summary = json.loads(paths["summary"].read_text())
    sample_by_uid = {row["sample_uid"]: row for row in sample_context}
    group_by_id = {row["group_id"]: row for row in groups}
    members_by_group: dict[str, list[dict[str, str]]] = {}
    membership_count: dict[str, int] = {}
    for member in members:
        members_by_group.setdefault(member["group_id"], []).append(member)
        membership_count[member["sample_uid"]] = membership_count.get(member["sample_uid"], 0) + 1
    excluded_current = {
        row["sample_uid"]
        for row in decisions
        if row["current_bigwig_group_eligible"] == "False"
    }
    ontology_ids = {row["curie"] for row in ontology if row["curie"]}
    formal_groups = [row for row in groups if row["include_formal_v2"] == "True"]
    context_fields = [
        "sra_study_accession",
        "assay",
        "strandedness",
        "tissue_or_cell_type",
        "development_stage",
        "sex",
        "strain",
        "genotype",
        "condition",
    ]
    context_consistent = True
    for group in groups:
        for member in members_by_group.get(group["group_id"], []):
            sample = sample_by_uid[member["sample_uid"]]
            group_field_names = {
                "strandedness": "strand",
            }
            for field in context_fields:
                group_field = group_field_names.get(field, field)
                if sample[field] != group[group_field]:
                    context_consistent = False
    t38_groups = {
        member["run_accession"]: member["group_id"]
        for member in members
        if member["run_accession"] in {"SRR941651", "SRR941697"}
    }
    t38_stages = {
        accession: group_by_id[group_id]["development_stage"]
        for accession, group_id in t38_groups.items()
    }
    duplicate_secondaries = {
        "rna_seq:SRR13964554",
        "rna_seq:SRR13964555",
        "rna_seq:SRR13964556",
    }
    chip_ids = {
        row["sample_uid"] for row in samples if row["assay"] != "RNA-Seq"
    }
    multi_groups = [row for row in groups if int(row["n_members"]) > 1]
    biological_multi_groups = [
        row for row in groups if int(row["n_biological_units"]) > 1
    ]
    technical_multi_groups = [
        row for row in groups if int(row["n_runs"]) > int(row["n_biological_units"])
    ]
    qc_by_group = {row["group_id"]: row for row in group_qc}
    expected_pairs = sum(
        int(row["n_members"]) * (int(row["n_members"]) - 1) // 2
        for row in multi_groups
    )
    resolution_statuses = {row["status"] for row in resolutions}

    checks.extend(
        [
            check(
                "R2.02_scope_and_decisions",
                len(samples) == 485
                and len(sample_context) == 485
                and len(decisions) == 485
                and len({row["sample_uid"] for row in decisions}) == 485,
                f"samples={len(samples)} contexts={len(sample_context)} decisions={len(decisions)}",
            ),
            check(
                "R2.03_wrong_assay_excluded",
                len(chip_ids) == 3
                and chip_ids.issubset(excluded_current)
                and not chip_ids.intersection(membership_count),
                f"chip_ids={sorted(chip_ids)}",
            ),
            check(
                "R2.04_duplicate_signals",
                duplicate_secondaries.issubset(excluded_current)
                and not duplicate_secondaries.intersection(membership_count)
                and all(
                    any(
                        row["sample_uid"] == uid
                        and row["decision"] == "exclude_duplicate_signal"
                        for row in exclusions
                    )
                    for uid in duplicate_secondaries
                ),
                "three secondary byte-identical signals excluded while provenance remains",
            ),
            check(
                "R2.05_single_membership",
                all(count == 1 for count in membership_count.values())
                and len(members) == summary["candidate_signal_members"],
                f"members={len(members)} max_memberships={max(membership_count.values(), default=0)}",
            ),
            check(
                "R2.06_group_member_counts",
                all(
                    int(group["n_members"])
                    == len(members_by_group.get(group["group_id"], []))
                    for group in groups
                ),
                f"groups={len(groups)}",
            ),
            check(
                "R2.07_context_consistency",
                context_consistent,
                "all members equal their group across the locked grouping fields",
            ),
            check(
                "R2.08_t38_split",
                set(t38_groups) == {"SRR941651", "SRR941697"}
                and len(set(t38_groups.values())) == 2
                and t38_stages
                == {"SRR941651": "L1-L3(38H)", "SRR941697": "L4(38H)"},
                f"groups={t38_groups} stages={t38_stages}",
            ),
            check(
                "R2.09_ontology",
                all(
                    (
                        group["anatomy_curie"] in ontology_ids
                        and group["life_stage_curie"] in ontology_ids
                    )
                    or bool(group["ontology_exclusion_reason"])
                    for group in groups
                )
                and all(row["ontology_sha256"] for row in ontology),
                f"mapping_rows={len(ontology)} ontology_ids={len(ontology_ids)}",
            ),
            check(
                "R2.10_replicate_relationship",
                all(
                    group["replicate_relation_status"]
                    in {
                        "source_attribute_biological_replicates",
                        "source_title_or_alias_biological_replicates",
                    }
                    and all(
                        member["biological_replicate"]
                        and member["biological_replicate_evidence"]
                        for member in members_by_group[group["group_id"]]
                    )
                    for group in biological_multi_groups
                )
                and all(
                    member["technical_unit_id"]
                    and member["technical_replicate_status"]
                    for member in members
                )
                and all(
                    group["replicate_relation_status"]
                    in {
                        "single_biological_sample_with_technical_runs",
                        "source_attribute_biological_replicates",
                        "source_title_or_alias_biological_replicates",
                    }
                    for group in technical_multi_groups
                ),
                f"biological_multi_groups={len(biological_multi_groups)} technical_multi_groups={len(technical_multi_groups)}",
            ),
            check(
                "R2.11_replicate_qc",
                len(pair_qc) == expected_pairs
                and len(group_qc) == len(multi_groups)
                and all(group["group_id"] in qc_by_group for group in multi_groups)
                and all(
                    row["log1p_1kb_pearson"]
                    and row["log1p_gene_body_pearson"]
                    and row["pair_type"]
                    and row["review_action"]
                    for row in pair_qc
                ),
                f"pairs={len(pair_qc)} expected={expected_pairs} groups={len(group_qc)}",
            ),
            check(
                "R2.12_no_automatic_qc_exclusion",
                all(
                    row["review_action"]
                    in {"none", "manual_review_no_automatic_exclusion"}
                    for row in pair_qc
                ),
                f"flagged_pairs={sum(row['review_flag'] == 'True' for row in pair_qc)}",
            ),
            check(
                "R2.13_conflict_resolutions",
                {"resolved", "resolved_current_signal_raw_reads_require_G1_review"}.issubset(
                    resolution_statuses
                )
                and all(row["authority"] for row in resolutions),
                f"resolution_rows={len(resolutions)} statuses={sorted(resolution_statuses)}",
            ),
            check(
                "R2.14_formal_fail_closed",
                len(formal_groups) == 0
                and all(
                    row["review_outcome"]
                    == "candidate_only_do_not_generate_formal_track"
                    for row in scientific_review
                ),
                f"formal_groups={len(formal_groups)} candidate_reviews={len(scientific_review)}",
            ),
            check(
                "R2.15_determinism",
                summary["deterministic_manifest_hash_match"] is True,
                "group construction repeated in memory with identical stable digest",
            ),
        ]
    )
    return {
        "schema_version": 1,
        "phase": "P2",
        "review": "R2",
        "reviewed_at": utc_now(),
        "status": "PASS"
        if all(item["status"] == "PASS" for item in checks)
        else "FAIL",
        "checks": checks,
    }


def bigwig_summary(path: Path, expected_chromosomes: dict[str, int]) -> dict[str, float]:
    import pyBigWig

    with pyBigWig.open(str(path)) as bigwig:
        if not bigwig.isBigWig():
            raise RuntimeError(f"Not a bigWig: {path}")
        if bigwig.chroms() != expected_chromosomes:
            raise RuntimeError(f"Chromosome mismatch: {path}")
        header = bigwig.header()
        decoded_total = 0.0
        for chromosome, length in expected_chromosomes.items():
            value = bigwig.stats(
                chromosome, 0, length, type="sum", exact=True
            )[0]
            if value is not None:
                decoded_total += float(value)
    summary = {
        key: float(header[key]) for key in ("minVal", "maxVal", "sumData", "sumSquared")
    }
    summary["decoded_total_signal"] = decoded_total
    if not all(math.isfinite(value) for value in summary.values()):
        raise RuntimeError(f"Non-finite bigWig header: {path}")
    if summary["minVal"] < 0 or summary["maxVal"] < 0:
        raise RuntimeError(f"Negative bigWig signal: {path}")
    return summary


def bigwig_1kb_profile(path: Path, chromosomes: dict[str, int]):
    import numpy as np
    import pyBigWig

    profiles = []
    with pyBigWig.open(str(path)) as bigwig:
        for chromosome in ("I", "II", "III", "IV", "V", "X"):
            length = chromosomes[chromosome]
            bins = math.ceil(length / 1000)
            values = bigwig.stats(
                chromosome,
                0,
                length,
                nBins=bins,
                type="mean",
                exact=True,
            )
            array = np.asarray(
                [0.0 if value is None else value for value in values],
                dtype=np.float64,
            )
            if not np.isfinite(array).all() or (array < 0).any():
                raise RuntimeError(f"Invalid 1 kb profile: {path}")
            profiles.append(array)
    return np.concatenate(profiles)


def pearson(first, second) -> float:
    import numpy as np

    if first.shape != second.shape or first.size == 0:
        return float("nan")
    if float(first.std()) == 0.0 or float(second.std()) == 0.0:
        return float("nan")
    return float(np.corrcoef(first, second)[0, 1])


def review_p3a() -> dict[str, Any]:
    import numpy as np

    metadata_dir = REPO_ROOT / "alphagenome_custom/metadata/v2"
    expected_work_rel = "shared/source_reads/v2/pilot_20260713"
    expected_output_rel = "alphagenome_custom/tracks/rna_seq_v2_normalized_pilot"
    expected_index_rel = "shared/reference_indexes/WBcel235_STAR_2.7.11b"
    approved_runs = {
        "SRR7443583",
        "SRR941632",
        "SRR10882545",
        "SRR23049895",
        "SRR36719198",
    }
    output_dir = REPO_ROOT / expected_output_rel
    pilot_manifest_path = metadata_dir / "p3_pilot_sources.tsv"
    sample_manifest_path = metadata_dir / "rna_seq_samples_v2.tsv"
    outputs_path = output_dir / "pilot_outputs.tsv"
    run_path = output_dir / "pilot_run.json"
    required = [pilot_manifest_path, sample_manifest_path, outputs_path, run_path]
    missing = [str(path.relative_to(REPO_ROOT)) for path in required if not path.is_file()]
    checks = [check("R3A.01_required_outputs", not missing, f"missing={missing}")]
    if missing:
        return {
            "schema_version": 1,
            "phase": "P3A",
            "review": "R3A",
            "reviewed_at": utc_now(),
            "status": "FAIL",
            "checks": checks,
        }

    pilot_manifest = read_tsv(pilot_manifest_path)
    samples = read_tsv(sample_manifest_path)
    outputs = read_tsv(outputs_path)
    run_metadata = json.loads(run_path.read_text())
    sample_by_run = {row["run_accession"]: row for row in samples}
    pilot_by_run = {row["run_accession"]: row for row in pilot_manifest}
    output_by_run = {row["run_accession"]: row for row in outputs}
    expected_runs = {row["run_accession"] for row in pilot_manifest}

    expected_chromosomes = {}
    with (REPO_ROOT / "alphagenome_custom/reference/genome.fa.fai").open() as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            expected_chromosomes[fields[0]] = int(fields[1])

    fastq_errors = []
    raw_errors = []
    output_errors = []
    mapping_errors = []
    bam_errors = []
    path_errors = []
    approved_work_dir = (REPO_ROOT / expected_work_rel).resolve()
    approved_output_dir = (REPO_ROOT / expected_output_rel).resolve()
    sra_manifest_path = REPO_ROOT / "alphagenome_custom/metadata/v2/p3_ncbi_sra_sources.tsv"
    sra_by_run = {
        row["run_accession"]: row for row in read_tsv(sra_manifest_path)
    }
    for accession in sorted(expected_runs):
        if accession not in output_by_run or accession not in sample_by_run:
            continue
        row = output_by_run[accession]
        sample = sample_by_run[accession]
        fastq_paths = [REPO_ROOT / value for value in row["source_fastq_paths"].split(";")]
        expected_md5 = row["source_fastq_md5"].split(";")
        expected_sha = row["source_fastq_sha256"].split(";")
        pilot_row = pilot_by_run[accession]
        backend = row.get("source_transport_backend", "ena_fastq")
        if backend == "ena_fastq":
            expected_names = [
                Path(value).name for value in pilot_row["fastq_ftp"].split(";")
            ]
            if (
                row["source_fastq_md5"] != pilot_row["fastq_md5"]
                or int(row["source_fastq_bytes"])
                != sum(int(value) for value in pilot_row["fastq_bytes"].split(";"))
                or [path.name for path in fastq_paths] != expected_names
            ):
                fastq_errors.append(f"{accession}:does not match pilot manifest")
        elif backend == "ncbi_sra" and accession == "SRR36719198":
            sra = sra_by_run[accession]
            if (
                row.get("source_reference_ena_fastq_md5") != pilot_row["fastq_md5"]
                or int(row.get("source_reference_ena_fastq_bytes", 0))
                != sum(int(value) for value in pilot_row["fastq_bytes"].split(";"))
                or [path.name for path in fastq_paths] != [f"{accession}.fastq"]
                or row.get("source_archive_url") != sra["sra_url"]
                or row.get("source_archive_md5") != sra["sra_md5"]
                or int(row.get("source_archive_bytes", 0)) != int(sra["sra_bytes"])
            ):
                fastq_errors.append(f"{accession}:SRA provenance mismatch")
            try:
                archive_path = REPO_ROOT / row["source_archive_path"]
                actual_archive_sha, actual_archive_md5 = file_hashes(archive_path)
                if (
                    actual_archive_sha != row["source_archive_sha256"]
                    or actual_archive_md5 != row["source_archive_md5"]
                    or not archive_path.resolve().is_relative_to(approved_work_dir)
                ):
                    raise ValueError("SRA archive integrity or path mismatch")
            except (KeyError, OSError, ValueError) as error:
                fastq_errors.append(f"{accession}:{error}")
        else:
            fastq_errors.append(f"{accession}:unapproved transport {backend}")
        if not (len(fastq_paths) == len(expected_md5) == len(expected_sha)):
            fastq_errors.append(f"{accession}:FASTQ field length mismatch")
        else:
            for path, md5_value, sha_value in zip(
                fastq_paths, expected_md5, expected_sha, strict=True
            ):
                if not path.is_file():
                    fastq_errors.append(f"{accession}:missing {path}")
                    continue
                actual_sha, actual_md5 = file_hashes(path)
                if actual_sha != sha_value or actual_md5 != md5_value:
                    fastq_errors.append(f"{accession}:checksum {path.name}")
                if not path.resolve().is_relative_to(approved_work_dir):
                    path_errors.append(f"FASTQ outside work dir: {path}")

        provided_path = REPO_ROOT / sample["local_path"]
        actual_provided_sha = sha256(provided_path)
        raw_values = {
            sample["sha256"],
            row["provided_bigwig_expected_sha256"],
            row["provided_bigwig_sha256_before"],
            row["provided_bigwig_sha256_after"],
            actual_provided_sha,
        }
        if len(raw_values) != 1 or row["provided_bigwig_unchanged"] != "True":
            raw_errors.append(accession)

        output_path = REPO_ROOT / row["output_path"]
        if not output_path.resolve().is_relative_to(approved_output_dir):
            path_errors.append(f"output outside output dir: {output_path}")
        try:
            summary = bigwig_summary(output_path, expected_chromosomes)
            if output_path.stat().st_size != int(row["output_size_bytes"]):
                output_errors.append(f"{accession}:size")
            if sha256(output_path) != row["output_sha256"]:
                output_errors.append(f"{accession}:sha256")
            relative_error = abs(summary["sumData"] - 100_000_000.0) / 100_000_000.0
            if relative_error > 1e-5:
                output_errors.append(f"{accession}:total={summary['sumData']}")
            if not math.isclose(
                summary["sumData"],
                float(row["output_total_signal"]),
                rel_tol=1e-12,
                abs_tol=1e-6,
            ):
                output_errors.append(f"{accession}:recorded total")
            if not math.isclose(
                summary["minVal"], float(row["output_min"]), rel_tol=1e-9, abs_tol=1e-9
            ) or not math.isclose(
                summary["maxVal"], float(row["output_max"]), rel_tol=1e-9, abs_tol=1e-9
            ):
                output_errors.append(f"{accession}:recorded extrema")
        except Exception as error:
            output_errors.append(f"{accession}:{error}")

        try:
            mapped_percent = float(row["star_uniquely_mapped_percent"].rstrip("%"))
            if not 0 <= mapped_percent <= 100:
                raise ValueError(mapped_percent)
            if int(row["star_input_reads"]) <= 0 or int(row["star_uniquely_mapped_reads"]) < 0:
                raise ValueError("invalid read counts")
        except (KeyError, ValueError) as error:
            mapping_errors.append(f"{accession}:{error}")

        try:
            primary_bam = REPO_ROOT / row["primary_bam_path"]
            if not primary_bam.resolve().is_relative_to(approved_work_dir):
                raise ValueError("BAM outside approved work directory")
            if primary_bam.stat().st_size != int(row["primary_bam_size_bytes"]):
                raise ValueError("BAM size mismatch")
            if sha256(primary_bam) != row["primary_bam_sha256"]:
                raise ValueError("BAM SHA-256 mismatch")
        except (KeyError, OSError, ValueError) as error:
            bam_errors.append(f"{accession}:{error}")

    commands_path = REPO_ROOT / run_metadata.get("commands_path", "missing")
    command_integrity = (
        commands_path.is_file()
        and sha256(commands_path) == run_metadata.get("commands_sha256")
        and len(json.loads(commands_path.read_text())) == run_metadata.get("command_count")
        and run_metadata.get("command_count", 0) > 0
    )
    tool_versions = run_metadata.get("tool_versions", {})
    expected_tool_versions = (
        "2.7.11b" in tool_versions.get("STAR", "")
        and "1.23" in tool_versions.get("samtools", "")
        and "2.31.1" in tool_versions.get("bedtools", "")
        and bool(tool_versions.get("bedGraphToBigWig"))
        and "3.4.1" in tool_versions.get("fasterq-dump", "")
        and "3.4.1" in tool_versions.get("vdb-validate", "")
    )
    partial_files = [
        str(path.relative_to(REPO_ROOT))
        for root in (
            REPO_ROOT / run_metadata["work_dir"],
            REPO_ROOT / run_metadata["output_dir"],
            REPO_ROOT / run_metadata["star_index"],
        )
        for pattern in ("*.part", "*.tmp")
        for path in root.rglob(pattern)
    ]

    comparisons = []
    comparison_errors = []
    for accession in sorted(expected_runs.intersection(output_by_run, sample_by_run)):
        output_path = REPO_ROOT / output_by_run[accession]["output_path"]
        provided_path = REPO_ROOT / sample_by_run[accession]["local_path"]
        try:
            output_profile = bigwig_1kb_profile(output_path, expected_chromosomes)
            provided_profile = bigwig_1kb_profile(provided_path, expected_chromosomes)
            raw_pearson = pearson(output_profile, provided_profile)
            log1p_pearson = pearson(
                np.log1p(output_profile), np.log1p(provided_profile)
            )
            comparisons.append(
                {
                    "run_accession": accession,
                    "bin_count": output_profile.size,
                    "raw_1kb_pearson": f"{raw_pearson:.12g}",
                    "log1p_1kb_pearson": f"{log1p_pearson:.12g}",
                    "output_zero_bin_fraction": f"{np.mean(output_profile == 0):.12g}",
                    "provided_zero_bin_fraction": f"{np.mean(provided_profile == 0):.12g}",
                    "interpretation": "diagnostic_only_unknown_provided_bigwig_unit",
                }
            )
            if not math.isfinite(raw_pearson) or not math.isfinite(log1p_pearson):
                comparison_errors.append(f"{accession}:non-finite correlation")
        except Exception as error:
            comparison_errors.append(f"{accession}:{error}")
    if comparisons:
        write_tsv(
            AUDIT_ROOT / "P3A/pilot_bigwig_comparison.tsv",
            comparisons,
            list(comparisons[0]),
        )

    formula_valid = all(
        row.get("normalization_formula")
        == "100000000/sum((end-start)*raw_coverage)"
        and float(row.get("bedgraph_scale_to_1e8", "nan")) > 0
        and float(row.get("raw_bedgraph_total", "nan")) > 0
        for row in outputs
    )
    checks.extend(
        [
            check(
                "R3A.02_bounded_scope",
                len(pilot_manifest) == 5
                and expected_runs == approved_runs
                and run_metadata.get("scope") == "five_run_reprocessing_pilot_only"
                and run_metadata.get("formal_v2_outputs") is False
                and run_metadata.get("source_fastq_bytes") == 6_221_459_671,
                f"runs={sorted(expected_runs)} formal={run_metadata.get('formal_v2_outputs')}",
            ),
            check(
                "R3A.03_manifest_integrity",
                run_metadata.get("manifest_sha256") == sha256(pilot_manifest_path)
                and run_metadata.get("sample_manifest_sha256")
                == sha256(sample_manifest_path)
                and run_metadata.get("sra_manifest_sha256")
                == sha256(sra_manifest_path)
                and run_metadata.get("sra_fallback_runs") == ["SRR36719198"],
                "pilot, sample, and locked SRA manifest SHA-256 values match the executed record",
            ),
            check(
                "R3A.04_output_membership",
                set(output_by_run) == expected_runs
                and len(list(output_dir.glob("*.bw"))) == 5,
                f"outputs={sorted(output_by_run)}",
            ),
            check(
                "R3A.05_fastq_integrity",
                not fastq_errors
                and {
                    row.get("source_transport_backend", "ena_fastq")
                    for row in outputs
                }
                == {"ena_fastq", "ncbi_sra"}
                and {
                    row["run_accession"]
                    for row in outputs
                    if row.get("source_transport_backend") == "ncbi_sra"
                }
                == {"SRR36719198"},
                f"errors={fastq_errors}",
            ),
            check(
                "R3A.06_raw_bigwig_immutability",
                not raw_errors and run_metadata.get("raw_bigwigs_unchanged") is True,
                f"errors={raw_errors}",
            ),
            check(
                "R3A.07_tool_and_command_provenance",
                command_integrity
                and expected_tool_versions,
                f"commands={run_metadata.get('command_count')} versions={tool_versions}",
            ),
            check(
                "R3A.08_mapping_metrics",
                not mapping_errors and not bam_errors,
                f"mapping_errors={mapping_errors} bam_errors={bam_errors}",
            ),
            check(
                "R3A.09_normalization_formula",
                formula_valid
                and run_metadata.get("normalization_target_total") == 100_000_000,
                "explicit 1e6 x 100 bp total-signal scale",
            ),
            check(
                "R3A.10_bigwig_integrity",
                not output_errors,
                f"errors={output_errors}",
            ),
            check(
                "R3A.11_path_and_partial_file_safety",
                not path_errors
                and not partial_files
                and run_metadata.get("partial_files") == []
                and run_metadata.get("work_dir") == expected_work_rel
                and run_metadata.get("output_dir") == expected_output_rel
                and run_metadata.get("star_index") == expected_index_rel
                and run_metadata.get("threads") == 16,
                f"path_errors={path_errors} partial_files={partial_files}",
            ),
            check(
                "R3A.12_resource_measurement",
                float(run_metadata.get("elapsed_seconds", 0)) > 0
                and int(run_metadata.get("observed_managed_peak_bytes", 0)) > 0
                and all(float(row.get("sample_elapsed_seconds", 0)) > 0 for row in outputs),
                "elapsed="
                f"{run_metadata.get('elapsed_seconds')} observed_bytes="
                f"{run_metadata.get('observed_managed_peak_bytes')}",
            ),
            check(
                "R3A.13_provided_signal_comparison",
                len(comparisons) == 5 and not comparison_errors,
                f"comparisons={len(comparisons)} errors={comparison_errors}",
            ),
        ]
    )
    return {
        "schema_version": 1,
        "phase": "P3A",
        "review": "R3A",
        "reviewed_at": utc_now(),
        "status": "PASS"
        if all(item["status"] == "PASS" for item in checks)
        else "FAIL",
        "checks": checks,
    }


def review_p3b() -> dict[str, Any]:
    metadata_dir = REPO_ROOT / "alphagenome_custom/metadata/v2"
    normalized_dir = REPO_ROOT / "alphagenome_custom/tracks/rna_seq_v2_normalized"
    grouped_dir = REPO_ROOT / "alphagenome_custom/tracks/rna_seq_v2_grouped"
    paths = {
        "samples": metadata_dir / "rna_seq_samples_v2.tsv",
        "contexts": metadata_dir / "rna_seq_sample_context_v2.tsv",
        "ledger": metadata_dir / "p3_full_outputs.tsv",
        "full_summary": metadata_dir / "p3_full_summary.json",
        "groups": metadata_dir / "rna_seq_groups_v2_final.tsv",
        "members": metadata_dir / "rna_seq_group_members_v2_final.tsv",
        "group_outputs": metadata_dir / "p3_group_outputs.tsv",
        "pair_qc": metadata_dir / "replicate_pair_qc_v2_post_reprocessing.tsv",
        "group_qc": metadata_dir / "replicate_group_qc_v2_post_reprocessing.tsv",
        "duplicate_review": metadata_dir / "duplicate_source_reuse_review_v2.tsv",
        "group_summary": metadata_dir / "p3_group_summary.json",
        "execution": metadata_dir / "p3b_execution.json",
        "sra_sources": metadata_dir / "p3_ncbi_sra_sources.tsv",
        "sra_summary": metadata_dir / "p3_ncbi_sra_summary.json",
        "audit_schema": metadata_dir / "p3_audit_schema_summary.json",
    }
    missing = [
        str(path.relative_to(REPO_ROOT))
        for path in paths.values()
        if not path.is_file()
    ]
    checks = [check("R3.01_required_outputs", not missing, f"missing={missing}")]
    if missing:
        return {
            "schema_version": 1,
            "phase": "P3B",
            "review": "R3",
            "reviewed_at": utc_now(),
            "status": "FAIL",
            "checks": checks,
        }

    samples = read_tsv(paths["samples"])
    contexts = read_tsv(paths["contexts"])
    ledger = read_tsv(paths["ledger"])
    groups = read_tsv(paths["groups"])
    members = read_tsv(paths["members"])
    group_outputs = read_tsv(paths["group_outputs"])
    pair_qc = read_tsv(paths["pair_qc"])
    group_qc = read_tsv(paths["group_qc"])
    duplicate_review = read_tsv(paths["duplicate_review"])
    full_summary = json.loads(paths["full_summary"].read_text())
    group_summary = json.loads(paths["group_summary"].read_text())
    execution = json.loads(paths["execution"].read_text())
    sra_sources = read_tsv(paths["sra_sources"])
    sra_summary = json.loads(paths["sra_summary"].read_text())
    audit_schema = json.loads(paths["audit_schema"].read_text())
    rna_runs = {row["run_accession"] for row in samples if row["assay"] == "RNA-Seq"}
    sample_by_run = {row["run_accession"]: row for row in samples}
    context_by_run = {row["run_accession"]: row for row in contexts}
    group_by_id = {row["group_id"]: row for row in groups}
    members_by_group: dict[str, list[dict[str, str]]] = {}
    for row in members:
        members_by_group.setdefault(row["group_id"], []).append(row)

    expected_chromosomes = {}
    with (REPO_ROOT / "alphagenome_custom/reference/genome.fa.fai").open() as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            expected_chromosomes[fields[0]] = int(fields[1])

    normalized_errors = []
    accepted_transports = {
        "ncbi_sra",
        "ena_fastq_fallback_after_sra_extraction_failure",
    }
    for row in ledger:
        accession = row["run_accession"]
        path = REPO_ROOT / row["output_path"]
        try:
            if path.parent.resolve() != normalized_dir.resolve():
                raise ValueError("outside normalized output directory")
            stats = bigwig_summary(path, expected_chromosomes)
            if sha256(path) != row["output_sha256"]:
                raise ValueError("SHA-256 mismatch")
            if path.stat().st_size != int(row["output_size_bytes"]):
                raise ValueError("size mismatch")
            if abs(stats["sumData"] - 100_000_000.0) / 100_000_000.0 > 1e-5:
                raise ValueError(f"total={stats['sumData']}")
            recorded_decoded_total = float(row["output_decoded_total_signal"])
            decoded_rescale = float(row["decoded_rescale_to_1e8"])
            if (
                abs(stats["decoded_total_signal"] - recorded_decoded_total)
                / recorded_decoded_total
                > 1e-9
            ):
                raise ValueError("decoded total audit mismatch")
            if (
                abs(recorded_decoded_total * decoded_rescale - 100_000_000.0)
                / 100_000_000.0
                > 1e-9
            ):
                raise ValueError("decoded rescale does not restore target total")
            if row.get("decoded_total_method") != (
                "sum_of_pyBigWig_exact_per_reference_chromosome"
            ):
                raise ValueError("decoded total method mismatch")
            if row.get("decoded_total_source_sha256") != row["output_sha256"]:
                raise ValueError("decoded total source SHA mismatch")
            if row["cleanup_status"] != "completed":
                raise ValueError("cleanup incomplete")
            if row["coverage_policy"] != "primary_unique_spliced_unstranded":
                raise ValueError("coverage policy mismatch")
            if row["output_strand"] != ".":
                raise ValueError("output strand mismatch")
            backend = row.get("source_transport_backend")
            reference_md5 = row.get("source_reference_ena_fastq_md5", "").split(";")
            source_spots = int(row.get("source_manifest_spot_count", 0))
            star_input_reads = int(row.get("star_input_reads", 0).replace(",", ""))
            common_source_valid = (
                row.get("audit_schema_version") == "2"
                and backend in accepted_transports
                and row.get("source_reference_ena_fastq_urls")
                and all(len(value) == 32 for value in reference_md5)
                and int(row.get("source_reference_ena_fastq_bytes", 0)) > 0
                and source_spots > 0
                and star_input_reads == source_spots
                and len(row.get("source_archive_md5", "")) == 32
                and len(row.get("source_archive_sha256", "")) == 64
                and row.get("sra_archive_validated") == "True"
                and any(
                    command == "STAR"
                    for command in (
                        item["argv"][0].rsplit("/", 1)[-1]
                        for item in json.loads(row["commands_json"])
                    )
                )
                and not any(
                    ambiguous in row
                    for ambiguous in (
                        "source_fastq_urls",
                        "source_fastq_md5",
                        "source_fastq_sha256",
                        "source_fastq_bytes",
                    )
                )
            )
            if backend == "ncbi_sra":
                transport_valid = (
                    bool(row.get("extracted_fastq_sha256"))
                    and int(row.get("extracted_fastq_bytes", 0)) > 0
                    and int(row.get("extracted_fastq_record_count", 0)) > 0
                    and row.get("extracted_fastq_integrity")
                    == "local_SHA256_and_layout_aware_record_count"
                )
            else:
                downloaded_sha = row.get("downloaded_ena_fastq_sha256", "").split(";")
                transport_valid = (
                    all(len(value) == 64 for value in downloaded_sha)
                    and int(row.get("downloaded_ena_fastq_bytes", 0))
                    == int(row.get("source_reference_ena_fastq_bytes", 0))
                    and int(row.get("ena_fastq_expected_record_count", 0)) > 0
                    and int(row.get("sra_extraction_return_code", 0)) != 0
                    and row.get("sra_extraction_error_type")
                    == "CalledProcessError"
                    and row.get("downloaded_ena_fastq_integrity")
                    == "ENA_MD5_and_local_SHA256_plus_STAR_spot_count"
                )
            if not common_source_valid or not transport_valid:
                raise ValueError("source audit schema/provenance mismatch")
        except Exception as error:
            normalized_errors.append(f"{accession}:{error}")

    raw_errors = []
    for row in samples:
        path = REPO_ROOT / row["local_path"]
        if not path.is_file() or sha256(path) != row["sha256"]:
            raw_errors.append(row["run_accession"])

    membership_counts: dict[str, int] = {}
    for row in members:
        membership_counts[row["run_accession"]] = (
            membership_counts.get(row["run_accession"], 0) + 1
        )
    context_errors = []
    context_fields = (
        "sra_study_accession",
        "bioproject_accession",
        "geo_accession",
        "assay",
        "tissue_or_cell_type",
        "development_stage",
        "sex",
        "strain",
        "genotype",
        "condition",
    )
    for member in members:
        group = group_by_id.get(member["group_id"])
        context = context_by_run.get(member["run_accession"])
        if group is None or context is None:
            context_errors.append(f"{member['run_accession']}:missing group/context")
            continue
        differences = [field for field in context_fields if group[field] != context[field]]
        if differences:
            context_errors.append(f"{member['run_accession']}:{','.join(differences)}")

    group_output_errors = []
    output_by_group = {row["group_id"]: row for row in group_outputs}
    for group in groups:
        group_id = group["group_id"]
        row = output_by_group.get(group_id)
        if row is None:
            group_output_errors.append(f"{group_id}:missing ledger row")
            continue
        path = REPO_ROOT / row["output_path"]
        try:
            if path.parent.resolve() != grouped_dir.resolve():
                raise ValueError("outside grouped output directory")
            stats = bigwig_summary(path, expected_chromosomes)
            if sha256(path) != row["output_sha256"]:
                raise ValueError("SHA-256 mismatch")
            if abs(stats["sumData"] - 100_000_000.0) / 100_000_000.0 > 1e-5:
                raise ValueError(f"total={stats['sumData']}")
            decoded_relative_error = (
                abs(stats["decoded_total_signal"] - 100_000_000.0)
                / 100_000_000.0
            )
            if decoded_relative_error > 1e-6:
                raise ValueError(
                    f"decoded_total={stats['decoded_total_signal']}"
                )
            if (
                abs(
                    stats["decoded_total_signal"]
                    - float(row["output_decoded_total_signal"])
                )
                / 100_000_000.0
                > 1e-9
            ):
                raise ValueError("group decoded total audit mismatch")
            if float(row["output_reconstruction_relative_error"]) > 1e-6:
                raise ValueError("group reconstruction error exceeds tolerance")
            expected_members = members_by_group.get(group_id, [])
            if int(row["n_runs"]) != len(expected_members):
                raise ValueError("run count mismatch")
            expected_policy = (
                "decoded_signal_rescaled_to_1e8_then_raw_coverage_weighted_within_"
                "biological_unit_then_equal_mean_across_units"
            )
            if row["aggregation_policy"] != expected_policy:
                raise ValueError("aggregation policy mismatch")
            if row["output_is_symlink"] != "False":
                raise ValueError("formal grouped output must be materialized")
            if (
                len(row["source_decoded_total_signal"].split(";"))
                != len(expected_members)
                or len(row["source_decoded_rescale_to_1e8"].split(";"))
                != len(expected_members)
            ):
                raise ValueError("member normalization audit count mismatch")
        except Exception as error:
            group_output_errors.append(f"{group_id}:{error}")

    expected_pairs = sum(
        len(group_members) * (len(group_members) - 1) // 2
        for group_members in members_by_group.values()
    )
    multi_groups = {
        group_id for group_id, group_members in members_by_group.items()
        if len(group_members) > 1
    }
    pair_group_errors = []
    member_group = {row["run_accession"]: row["group_id"] for row in members}
    for row in pair_qc:
        if not (
            member_group.get(row["sample_a"]) == row["group_id"]
            and member_group.get(row["sample_b"]) == row["group_id"]
        ):
            pair_group_errors.append(f"{row['sample_a']}:{row['sample_b']}")

    duplicate_pairs = {
        (row["canonical_run"], row["secondary_run"]) for row in duplicate_review
    }
    expected_duplicate_pairs = {
        ("SRR12483327", "SRR13964554"),
        ("SRR12483328", "SRR13964555"),
        ("SRR12483329", "SRR13964556"),
    }
    duplicate_valid = (
        duplicate_pairs == expected_duplicate_pairs
        and all(row["same_source_study"] == "False" for row in duplicate_review)
        and all(
            row["canonical_group_id"] != row["secondary_group_id"]
            and row["decision"]
            == "retain_distinct_study_tracks_no_cross_source_averaging"
            for row in duplicate_review
        )
    )

    failures_dir = normalized_dir / "failures"
    partial_files = [
        str(path.relative_to(REPO_ROOT))
        for root in (
            normalized_dir,
            grouped_dir,
            REPO_ROOT / "shared/source_reads/v2/full_streaming",
        )
        if root.exists()
        for pattern in ("*.part", "*.tmp")
        for path in root.rglob(pattern)
    ]
    checks.extend(
        [
            check(
                "R3.02_full_scope",
                len(rna_runs) == 482
                and set(row["run_accession"] for row in ledger) == rna_runs
                and full_summary.get("expected_runs") == 482
                and full_summary.get("completed_runs") == 482
                and full_summary.get("failures") == [],
                f"rna_runs={len(rna_runs)} ledger={len(ledger)} failures={full_summary.get('failures')}",
            ),
            check(
                "R3.03_normalized_bigwig_integrity",
                not normalized_errors and len(list(normalized_dir.glob("*.bw"))) == 482,
                f"errors={normalized_errors[:10]} count={len(list(normalized_dir.glob('*.bw')))}",
            ),
            check(
                "R3.04_raw_immutability",
                not raw_errors,
                f"changed_or_missing={raw_errors}",
            ),
            check(
                "R3.05_final_membership",
                len(groups) == 241
                and len(members) == 482
                and set(membership_counts) == rna_runs
                and all(count == 1 for count in membership_counts.values())
                and all(
                    int(group["n_members"])
                    == len(members_by_group.get(group["group_id"], []))
                    for group in groups
                ),
                f"groups={len(groups)} members={len(members)}",
            ),
            check(
                "R3.06_context_and_formal_status",
                not context_errors
                and all(
                    row["group_status"] == "formal_v2_uniform_reprocessed"
                    and row["strand"] == "."
                    and row["include_formal_v2"] == "True"
                    and row["signal_unit"] == "1e6_x_100bp_unstranded_coverage"
                    and row["anatomy_curie"]
                    and row["life_stage_curie"]
                    for row in groups
                ),
                f"context_errors={context_errors[:10]}",
            ),
            check(
                "R3.07_group_output_integrity",
                len(group_outputs) == 241 and not group_output_errors,
                f"outputs={len(group_outputs)} errors={group_output_errors[:10]}",
            ),
            check(
                "R3.08_replicate_qc",
                len(pair_qc) == expected_pairs
                and len(group_qc) == len(multi_groups)
                and {row["group_id"] for row in group_qc} == multi_groups
                and not pair_group_errors
                and all(
                    row["review_action"]
                    in {"none", "manual_review_no_automatic_exclusion"}
                    for row in pair_qc
                ),
                f"pairs={len(pair_qc)}/{expected_pairs} groups={len(group_qc)}/{len(multi_groups)} errors={pair_group_errors[:10]}",
            ),
            check(
                "R3.09_duplicate_source_reuse_resolution",
                duplicate_valid,
                f"pairs={sorted(duplicate_pairs)}",
            ),
            check(
                "R3.10_cleanup_and_partial_safety",
                not partial_files
                and (not failures_dir.exists() or not list(failures_dir.glob("*.json")))
                and all(row["cleanup_status"] == "completed" for row in ledger),
                f"partial_files={partial_files} failures={len(list(failures_dir.glob('*.json'))) if failures_dir.exists() else 0}",
            ),
            check(
                "R3.11_manifest_provenance",
                group_summary.get("normalized_runs") == 482
                and group_summary.get("formal_groups") == 241
                and group_summary.get("group_outputs") == 241
                and group_summary.get("decoded_total_relative_tolerance") == 1e-6
                and group_summary.get(
                    "group_decoded_relative_target_error_max", 1.0
                )
                <= 1e-6
                and group_summary.get("sample_manifest_sha256") == sha256(paths["samples"])
                and group_summary.get("full_ledger_sha256") == sha256(paths["ledger"])
                and group_summary.get("final_group_manifest_sha256") == sha256(paths["groups"])
                and group_summary.get("final_member_manifest_sha256") == sha256(paths["members"]),
                "tracked manifest digests agree with the P3 group summary",
            ),
            check(
                "R3.12_execution_record",
                execution.get("status") == "completed"
                and execution.get("workers") == 4
                and execution.get("threads_per_sample") == 16
                and bool(execution.get("log_path")),
                f"status={execution.get('status')} log={execution.get('log_path')}",
            ),
            check(
                "R3.13_sra_transport_provenance",
                len(sra_sources) == 482
                and len({row["run_accession"] for row in sra_sources}) == 482
                and all(len(row["sra_md5"]) == 32 for row in sra_sources)
                and sra_summary.get("run_count") == 482
                and sra_summary.get("sra_manifest_sha256")
                == sha256(paths["sra_sources"])
                and full_summary.get("sra_manifest_sha256")
                == sha256(paths["sra_sources"]),
                f"sra_sources={len(sra_sources)} total_bytes={sra_summary.get('total_sra_bytes')}",
            ),
            check(
                "R3.14_unambiguous_audit_schema",
                audit_schema.get("schema_version") == 2
                and audit_schema.get("sample_audits") == 482
                and sum(audit_schema.get("transport_counts", {}).values()) == 482
                and set(audit_schema.get("transport_counts", {}))
                <= accepted_transports
                and audit_schema.get("transport_counts")
                == dict(
                    sorted(
                        Counter(
                            row["source_transport_backend"] for row in ledger
                        ).items()
                    )
                )
                and audit_schema.get("ledger_sha256") == sha256(paths["ledger"])
                and set(audit_schema.get("ambiguous_legacy_fields_removed", []))
                == {
                    "source_fastq_urls",
                    "source_fastq_md5",
                    "source_fastq_sha256",
                    "source_fastq_bytes",
                }
                and audit_schema.get("decoded_total_method")
                == "sum_of_pyBigWig_exact_per_reference_chromosome"
                and set(audit_schema.get("signal_total_fields", []))
                == {
                    "output_header_total_signal",
                    "output_decoded_total_signal",
                    "decoded_rescale_to_1e8",
                },
                f"schema={audit_schema.get('schema_version')} audits={audit_schema.get('sample_audits')}",
            ),
        ]
    )
    return {
        "schema_version": 1,
        "phase": "P3B",
        "review": "R3",
        "reviewed_at": utc_now(),
        "status": "PASS"
        if all(item["status"] == "PASS" for item in checks)
        else "FAIL",
        "checks": checks,
    }


def review_p4() -> dict[str, Any]:
    from scripts.evaluate_v2_model import core_subwindows
    from scripts.v2_bigwig_dataset import V2BigWigDataset

    metadata_dir = REPO_ROOT / "alphagenome_custom/metadata/v2"
    interval_dir = REPO_ROOT / "alphagenome_custom/intervals/v2"
    registry_path = metadata_dir / "split_registry_v2.json"
    benchmark_path = metadata_dir / "p4_loader_benchmark.json"
    spec_path = metadata_dir / "six_chromosome_split_v2_spec.json"
    blocks_path = interval_dir / "blocks.tsv"
    migration_path = metadata_dir / "six_chromosome_revision_migration.json"
    core_migration_path = metadata_dir / "six_chromosome_core_coverage_revision.json"
    expected_paths = [
        interval_dir / f"fold_{fold}/{role}.tsv"
        for fold in range(1, 6)
        for role in ("train", "valid")
    ] + [
        interval_dir / "development_train.tsv",
        interval_dir / "test_locked.tsv",
        blocks_path,
        spec_path,
        migration_path,
        core_migration_path,
        registry_path,
        benchmark_path,
    ]
    missing = [
        str(path.relative_to(REPO_ROOT)) for path in expected_paths if not path.is_file()
    ]
    checks = [check("R4.01_required_outputs", not missing, f"missing={missing}")]
    if missing:
        return {
            "schema_version": 1,
            "phase": "P4",
            "review": "R4",
            "reviewed_at": utc_now(),
            "status": "FAIL",
            "checks": checks,
        }

    registry = json.loads(registry_path.read_text())
    benchmark = json.loads(benchmark_path.read_text())
    spec = json.loads(spec_path.read_text())
    migration = json.loads(migration_path.read_text())
    core_migration = json.loads(core_migration_path.read_text())
    blocks = read_tsv(blocks_path)
    chromosome_lengths = {}
    with (REPO_ROOT / "alphagenome_custom/reference/genome.fa.fai").open() as handle:
        for line in handle:
            name, length, *_ = line.rstrip("\n").split("\t")
            chromosome_lengths[name] = int(length)
    chromosomes = ("I", "II", "III", "IV", "V", "X")
    chromosome_set = set(chromosomes)
    expected_assignments = {
        "cv_fold_1",
        "cv_fold_2",
        "cv_fold_3",
        "cv_fold_4",
        "cv_fold_5",
        "test_locked",
    }
    window_size = 2**20
    max_shift = 1024
    buffer_bp = window_size + 2 * max_shift
    revision_id = "six_chromosome_blocks_v2"
    split_errors = []
    core_errors = []
    block_errors = []
    blocks_by_id = {row["block_id"]: row for row in blocks}
    if len(blocks) != 36 or len(blocks_by_id) != 36:
        block_errors.append(f"block_count={len(blocks)} unique={len(blocks_by_id)}")
    for chromosome in chromosomes:
        members = sorted(
            (row for row in blocks if row["chromosome"] == chromosome),
            key=lambda row: int(row["physical_order"]),
        )
        if (
            len(members) != 6
            or {row["assignment"] for row in members} != expected_assignments
            or int(members[0]["block_start"])
            != int(members[0]["outer_margin_before_bp"])
            or int(members[-1]["block_end"])
            + int(members[-1]["outer_margin_after_bp"])
            != chromosome_lengths[chromosome]
            or any(int(row["effective_bases"]) < buffer_bp for row in members)
            or any(
                int(row["effective_bases"]) % 131072 != 0 for row in members
            )
            or any(
                int(second["block_start"]) - int(first["block_end"]) != buffer_bp
                for first, second in zip(members, members[1:])
            )
        ):
            block_errors.append(f"{chromosome}:invalid assignment/size/buffer")

    def validate_rows(
        rows: list[dict[str, str]], role: str, assignments: set[str]
    ) -> list[str]:
        errors = []
        for row in rows:
            block = blocks_by_id.get(row.get("block_id"))
            if (
                block is None
                or row.get("role") != role
                or row.get("split_revision") != revision_id
                or row.get("chromosome") != block.get("chromosome")
                or block.get("assignment") not in assignments
                or int(row["end"]) - int(row["start"]) != window_size
                or int(row["start"]) < int(block["block_start"])
                or int(row["end"]) > int(block["block_end"])
            ):
                errors.append(f"{row.get('chromosome')}:{row.get('start')}-{row.get('end')}")
                continue
            if role == "train" and (
                int(row["start"]) - max_shift < int(block["block_start"])
                or int(row["end"]) + max_shift > int(block["block_end"])
            ):
                errors.append(f"{row['block_id']}:shift_escape")
            if role != "train" and not (
                int(row["start"])
                <= int(row["core_start"])
                < int(row["core_end"])
                <= int(row["end"])
            ):
                errors.append(f"{row['block_id']}:core_escape")
        return errors

    def validate_core_partitions(rows: list[dict[str, str]]) -> list[str]:
        errors = []
        grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
        for row in rows:
            grouped[row["block_id"]].append(row)
        for block_id, members in grouped.items():
            block = blocks_by_id[block_id]
            ordered = sorted(members, key=lambda row: int(row["core_start"]))
            if (
                int(ordered[0]["core_start"]) != int(block["block_start"])
                or int(ordered[-1]["core_end"]) != int(block["block_end"])
                or len(ordered)
                != int(block["effective_bases"]) // 131072
                or any(
                    int(row["core_end"]) - int(row["core_start"]) != 131072
                    or (int(row["core_start"]) - int(row["start"])) % 128 != 0
                    for row in ordered
                )
                or any(
                    int(first["core_end"]) != int(second["core_start"])
                    for first, second in zip(ordered, ordered[1:])
                )
            ):
                errors.append(f"{block_id}:core_gap_or_overlap")
        return errors

    def validate_metric_coverage(rows: list[dict[str, str]]) -> list[str]:
        try:
            subwindows, eligible_bases = core_subwindows(rows, 131072)
        except (RuntimeError, ValueError) as error:
            return [str(error)]
        core_bases = sum(
            int(row["core_end"]) - int(row["core_start"]) for row in rows
        )
        if eligible_bases != core_bases:
            return [f"eligible={eligible_bases} core={core_bases}"]
        if len(subwindows) * 131072 != core_bases:
            return [f"subwindows={len(subwindows)} core={core_bases}"]
        if len(subwindows) != len(rows):
            return [f"subwindows={len(subwindows)} rows={len(rows)}"]
        return []

    for fold in range(1, 6):
        train = read_tsv(interval_dir / f"fold_{fold}/train.tsv")
        valid = read_tsv(interval_dir / f"fold_{fold}/valid.tsv")
        train_chromosomes = {row["chromosome"] for row in train}
        valid_chromosomes = {row["chromosome"] for row in valid}
        train_assignments = {
            f"cv_fold_{value}" for value in range(1, 6) if value != fold
        }
        valid_assignments = {f"cv_fold_{fold}"}
        if train_chromosomes != chromosome_set or valid_chromosomes != chromosome_set:
            split_errors.append(
                f"fold{fold}:train={sorted(train_chromosomes)} valid={sorted(valid_chromosomes)}"
            )
        train_errors = validate_rows(train, "train", train_assignments)
        valid_errors = validate_rows(valid, "valid", valid_assignments)
        if train_errors or valid_errors:
            split_errors.append(
                f"fold{fold}:row_errors={train_errors[:3] + valid_errors[:3]}"
            )
        fold_core_errors = validate_core_partitions(valid)
        if fold_core_errors:
            core_errors.extend(f"fold{fold}:{error}" for error in fold_core_errors)
        fold_metric_errors = validate_metric_coverage(valid)
        if fold_metric_errors:
            core_errors.extend(
                f"fold{fold}:metric:{error}" for error in fold_metric_errors
            )
        train_blocks = {row["block_id"] for row in train}
        valid_blocks = {row["block_id"] for row in valid}
        test_blocks = {
            row["block_id"] for row in blocks if row["assignment"] == "test_locked"
        }
        if train_blocks & valid_blocks or train_blocks & test_blocks or valid_blocks & test_blocks:
            split_errors.append(f"fold{fold}:block_leakage")

    development = read_tsv(interval_dir / "development_train.tsv")
    development_chromosomes = {row["chromosome"] for row in development}
    if (
        development_chromosomes != chromosome_set
        or any(row["role"] != "train" for row in development)
        or validate_rows(
            development,
            "train",
            {f"cv_fold_{fold}" for fold in range(1, 6)},
        )
        or {row["block_id"] for row in development}
        & {row["block_id"] for row in blocks if row["assignment"] == "test_locked"}
        or int(registry.get("development_train_windows", -1)) != len(development)
    ):
        split_errors.append("invalid development training split")

    test_rows = read_tsv(interval_dir / "test_locked.tsv")
    test_metric_errors = validate_metric_coverage(test_rows)
    test_valid = (
        {row["chromosome"] for row in test_rows} == chromosome_set
        and all(row["role"] == "test_locked" for row in test_rows)
        and not validate_rows(test_rows, "test_locked", {"test_locked"})
        and not validate_core_partitions(test_rows)
        and not test_metric_errors
        and len({row["block_id"] for row in test_rows}) == 6
    )
    registry_hash_errors = [
        relative
        for relative, expected_hash in registry["files"].items()
        if not (REPO_ROOT / relative).is_file()
        or sha256(REPO_ROOT / relative) != expected_hash
    ]
    try:
        V2BigWigDataset(interval_dir / "test_locked.tsv")
    except PermissionError:
        test_refused = True
    else:
        test_refused = False

    expected_shapes = {
        "dna_sequence": [4, 2**20],
        "target_1bp": [241, 2**20],
        "target_128bp": [241, 8192],
        "track_mask": [241, 1],
        "track_strand": [241],
        "gene_mask": [2, 2**20],
        "core_mask": [2**20],
    }
    expected_dtypes = {
        "dna_sequence": "torch.float32",
        "target_1bp": "torch.float32",
        "target_128bp": "torch.float32",
        "track_mask": "torch.bool",
        "track_strand": "torch.int8",
        "gene_mask": "torch.bool",
        "core_mask": "torch.bool",
    }
    forbidden_npz = list(
        (REPO_ROOT / "alphagenome_custom/datasets").glob("rna_seq_v2*.npz")
    ) + list(
        (REPO_ROOT / "alphagenome_custom/datasets").glob("rna_seq_v2_npz*")
    )
    checks.extend(
        [
            check(
                "R4.02_five_fold_split",
                not split_errors
                and all(
                    set(row.get("train_chromosomes", [])) == chromosome_set
                    and set(row.get("valid_chromosomes", [])) == chromosome_set
                    for row in registry.get("folds", [])
                ),
                f"errors={split_errors} chromosomes={list(chromosomes)}",
            ),
            check(
                "R4.03_block_buffers_no_leakage_and_full_metric_coverage",
                not block_errors and not core_errors and test_valid,
                f"block_errors={block_errors} core_errors={core_errors} test_metric_errors={test_metric_errors} test_valid={test_valid}",
            ),
            check(
                "R4.04_split_lock_integrity",
                not registry_hash_errors
                and registry.get("schema_version") == 2
                and registry.get("revision_id") == revision_id
                and registry.get("window_size") == window_size
                and registry.get("stride") == 2**19
                and registry.get("max_shift_bp") == max_shift
                and registry.get("exclusion_buffer_bp") == buffer_bp
                and registry.get("evaluation_subwindow_bp") == 131072
                and registry.get("split_spec_sha256") == sha256(spec_path)
                and registry.get("blocks_sha256") == sha256(blocks_path)
                and registry.get("test_chromosomes") == list(chromosomes)
                and registry.get("test_status")
                == "embargoed_six_chromosome_locked_blocks",
                f"hash_errors={registry_hash_errors}",
            ),
            check(
                "R4.05_loader_shapes_and_dtypes",
                benchmark.get("full_window_tracks") == 241
                and benchmark.get("full_window_shapes") == expected_shapes
                and benchmark.get("full_window_dtypes") == expected_dtypes,
                f"shapes={benchmark.get('full_window_shapes')} dtypes={benchmark.get('full_window_dtypes')}",
            ),
            check(
                "R4.06_direct_bigwig_and_pooling",
                benchmark.get("direct_bigwig_errors") == []
                and benchmark.get("pooling_128bp")
                == "sum_of_128_consecutive_1bp_values",
                f"errors={benchmark.get('direct_bigwig_errors')}",
            ),
            check(
                "R4.07_worker_determinism_and_handles",
                benchmark.get("worker_deterministic") is True
                and benchmark.get("single_worker_digests")
                == benchmark.get("multi_worker_digests")
                and int(benchmark.get("open_file_descriptors_after", 10**9))
                <= int(benchmark.get("open_file_descriptors_before", 0)) + 4,
                "worker digests agree and parent file descriptors returned near baseline",
            ),
            check(
                "R4.08_test_embargo",
                test_refused and test_valid,
                f"default_loader_refused={test_refused} test_windows={len(test_rows)} chromosomes={list(chromosomes)}",
            ),
            check(
                "R4.09_io_benchmark",
                benchmark.get("device") == "cpu"
                and float(benchmark.get("full_window_seconds", 0)) > 0
                and int(benchmark.get("full_target_bytes", 0)) > 0
                and float(benchmark.get("full_fold_epoch_seconds_linear_io_estimate", 0)) > 0
                and int(benchmark.get("peak_rss_kib_after", 0)) > 0,
                f"window_seconds={benchmark.get('full_window_seconds')} target_bytes={benchmark.get('full_target_bytes')}",
            ),
            check(
                "R4.10_no_monolithic_npz",
                benchmark.get("monolithic_npz_generated") is False and not forbidden_npz,
                f"forbidden_paths={[str(path) for path in forbidden_npz]}",
            ),
            check(
                "R4.11_data_manifest_provenance",
                benchmark.get("split_registry_sha256") == sha256(registry_path)
                and benchmark.get("split_revision") == revision_id
                and benchmark.get("x_valid_role_read_verified") is True
                and benchmark.get("locked_test_block_signal_reads") == 0
                and benchmark.get("group_manifest_sha256")
                == sha256(metadata_dir / "rna_seq_groups_v2_final.tsv")
                and benchmark.get("group_outputs_sha256")
                == sha256(metadata_dir / "p3_group_outputs.tsv"),
                "benchmark digests bind loader results to splits and P3 tracks",
            ),
            check(
                "R4.12_superseded_holdout_preserved",
                migration.get("revision_id") == "six_chromosome_blocks_v1"
                and migration.get("source_test_consumed") is False
                and migration.get("old_split_disposition")
                == "superseded_for_six_chromosome_split"
                and bool(migration.get("archived_intervals"))
                and bool(migration.get("archived_metadata"))
                and bool(migration.get("archived_audits")),
                f"archived_intervals={len(migration.get('archived_intervals', []))}",
            ),
            check(
                "R4.13_incomplete_core_revision_preserved",
                core_migration.get("source_revision")
                == "six_chromosome_blocks_v1"
                and core_migration.get("target_revision") == revision_id
                and core_migration.get("source_completed_jobs") == 1
                and core_migration.get("source_locked_test_block_signal_reads") == 0
                and min(
                    (
                        float(value)
                        for value in core_migration.get(
                            "observed_validation_core_coverage_fractions", []
                        )
                    ),
                    default=1.0,
                )
                < 0.999
                and bool(core_migration.get("archived_intervals"))
                and bool(core_migration.get("archived_metadata")),
                f"source_jobs={core_migration.get('source_completed_jobs')} coverage={core_migration.get('observed_validation_core_coverage_fractions')}",
            ),
        ]
    )
    return {
        "schema_version": 1,
        "phase": "P4",
        "review": "R4",
        "reviewed_at": utc_now(),
        "status": "PASS"
        if all(item["status"] == "PASS" for item in checks)
        else "FAIL",
        "checks": checks,
    }


def review_p5() -> dict[str, Any]:
    metadata_dir = REPO_ROOT / "alphagenome_custom/metadata/v2"
    paths = {
        "means": metadata_dir / "track_nonzero_means_v2.tsv",
        "means_summary": metadata_dir / "track_nonzero_means_v2_summary.json",
        "specs": metadata_dir / "model_specs_v2.json",
        "audit": metadata_dir / "p5_component_audit.json",
    }
    missing = [
        str(path.relative_to(REPO_ROOT))
        for path in paths.values()
        if not path.is_file()
    ]
    checks = [check("R5.01_required_outputs", not missing, f"missing={missing}")]
    if missing:
        return {
            "schema_version": 1,
            "phase": "P5",
            "review": "R5",
            "reviewed_at": utc_now(),
            "status": "FAIL",
            "checks": checks,
        }
    means = read_tsv(paths["means"])
    means_summary = json.loads(paths["means_summary"].read_text())
    specs = json.loads(paths["specs"].read_text())
    audit = json.loads(paths["audit"].read_text())
    mean_fields = [
        "development_train_nonzero_mean",
        *(f"fold_{fold}_train_nonzero_mean" for fold in range(1, 6)),
    ]
    means_valid = (
        len(means) == 241
        and len({row["group_id"] for row in means}) == 241
        and all(
            math.isfinite(float(row[field])) and float(row[field]) > 0
            for row in means
            for field in mean_fields
        )
    )
    model_by_id = {row["model_id"]: row for row in specs.get("models", [])}
    expected_embedding_paths = {
        "organism_embed",
        "embedder_128bp.organism_embed",
        "embedder_1bp.organism_embed",
        "embedder_pair.organism_embed",
    }
    checks.extend(
        [
            check(
                "R5.02_leakage_safe_track_means",
                means_valid
                and means_summary.get("locked_test_block_signal_reads") == 0
                and means_summary.get("fold_mean_policy")
                == "nonzero_mean_on_registered_train_blocks_only"
                and means_summary.get("split_revision")
                == "six_chromosome_blocks_v2"
                and means_summary.get("output_sha256") == sha256(paths["means"]),
                f"tracks={len(means)} test_reads={means_summary.get('locked_test_block_signal_reads')}",
            ),
            check(
                "R5.03_model_space_scaling_and_dual_loss",
                audit.get("scale_implementation")
                == "alphagenome_pytorch.heads.targets_scaling"
                and audit.get("loss_implementation")
                == "alphagenome_pytorch.losses.multinomial_loss_on_scaled_prediction_and_target"
                and specs.get("resolutions") == [1, 128]
                and specs.get("num_segments") == 8
                and audit.get("head_loss_finite") is True
                and audit.get("head_gradients_finite") is True,
                f"head_shapes={audit.get('head_prediction_shapes')} loss={audit.get('head_loss')}",
            ),
            check(
                "R5.04_gene_loss_and_augmentation_lock",
                math.isclose(specs.get("gene_cross_track_weight", -1), 0.1)
                and specs.get("augmentation")
                == {"max_shift_bp": 1024, "reverse_complement_probability": 0.5}
                and "gene_cross_track" in audit.get("loss_metrics", {}),
                f"gene_weight={specs.get('gene_cross_track_weight')} augmentation={specs.get('augmentation')}",
            ),
            check(
                "R5.05_three_model_contract",
                set(model_by_id) == {"A", "B", "C"}
                and model_by_id.get("A", {}).get("trunk_policy") == "frozen"
                and model_by_id.get("B", {}).get("trunk_policy")
                == "worm_embeddings_and_lora_only"
                and model_by_id.get("C", {}).get("trunk_policy") == "from_scratch",
                f"models={model_by_id}",
            ),
            check(
                "R5.06_real_worm_embedding",
                audit.get("worm_embedding_rows") == 3
                and audit.get("worm_embedding_index") == 2
                and audit.get("worm_embedding_finite") is True
                and set(audit.get("worm_embedding_paths", []))
                == expected_embedding_paths
                and specs.get("model_b_organism_index") == 2,
                f"paths={audit.get('worm_embedding_paths')} index={audit.get('worm_embedding_index')}",
            ),
            check(
                "R5.07_checkpoint_freeze_and_lora_scope",
                audit.get("base_checkpoint_loaded") is True
                and len(audit.get("base_checkpoint_sha256", "")) == 64
                and audit.get("unexpected_base_trainable_parameters") == []
                and int(audit.get("lora_trainable_parameters", 0)) > 0
                and audit.get("lora_modules")
                and audit.get("lora_targets")
                == ["tower.blocks.8.mha", "tower.blocks.8.mlp"],
                f"lora_modules={len(audit.get('lora_modules', []))} unexpected={audit.get('unexpected_base_trainable_parameters')}",
            ),
            check(
                "R5.08_cpu_forward_backward",
                audit.get("device") == "cpu"
                and audit.get("head_prediction_shapes")
                == {"1": [1, 241, 1024], "128": [1, 241, 8]}
                and audit.get("baseline_prediction_shapes")
                == {"1": [1, 16, 1024], "128": [1, 16, 8]}
                and audit.get("baseline_gradients_finite") is True,
                f"head={audit.get('head_prediction_shapes')} baseline={audit.get('baseline_prediction_shapes')}",
            ),
            check(
                "R5.09_unit_and_golden_tests",
                audit.get("unit_tests_return_code") == 0
                and int(audit.get("unit_test_count", 0)) >= 31,
                f"unit_tests={audit.get('unit_test_count')} return={audit.get('unit_tests_return_code')}",
            ),
            check(
                "R5.10_test_embargo",
                specs.get("chromosome_x_access")
                == "registered_train_valid_blocks_allowed"
                and specs.get("locked_test_block_access")
                == "prohibited_until_locked_G5_P6C"
                and specs.get("final_test_scope")
                == "r6c_single_six_chromosome_block_test"
                and audit.get("locked_test_block_signal_reads") == 0,
                f"x={specs.get('chromosome_x_access')} test={specs.get('locked_test_block_access')}",
            ),
            check(
                "R5.11_manifest_hashes",
                audit.get("means_sha256") == sha256(paths["means"])
                and audit.get("model_specs_sha256") == sha256(paths["specs"])
                and audit.get("split_registry_sha256")
                == sha256(metadata_dir / "split_registry_v2.json"),
                "P5 audit hashes bind means and model specification",
            ),
        ]
    )
    return {
        "schema_version": 1,
        "phase": "P5",
        "review": "R5",
        "reviewed_at": utc_now(),
        "status": "PASS"
        if all(item["status"] == "PASS" for item in checks)
        else "FAIL",
        "checks": checks,
    }


def review_p6a() -> dict[str, Any]:
    record_path = REPO_ROOT / "alphagenome_custom/metadata/v2/p6a_execution.json"
    checks = [
        check(
            "R6A.01_execution_record",
            record_path.is_file(),
            str(record_path.relative_to(REPO_ROOT)),
        )
    ]
    if not record_path.is_file():
        return {
            "schema_version": 1,
            "phase": "P6A",
            "review": "R6A",
            "reviewed_at": utc_now(),
            "status": "FAIL",
            "checks": checks,
        }
    record = json.loads(record_path.read_text())
    run_errors = []
    runs = []
    for item in record.get("runs", []):
        try:
            run_path = REPO_ROOT / item["run_path"]
            log_path = REPO_ROOT / item["log_path"]
            checkpoint_path = REPO_ROOT / item["checkpoint_path"]
            run = json.loads(run_path.read_text())
            if sha256(checkpoint_path) != item["checkpoint_sha256"]:
                raise ValueError("checkpoint hash mismatch")
            if not log_path.is_file() or log_path.stat().st_size == 0:
                raise ValueError("missing or empty log")
            if run.get("checkpoint_reload_verified") is not True:
                raise ValueError("checkpoint reload not verified")
            runs.append(run)
        except Exception as error:
            run_errors.append(f"{item.get('model')}:{error}")
    selected = record.get("selected_physical_gpu")
    gpu_by_index = {
        row["index"]: row for row in record.get("resources_before", {}).get("gpus", [])
    }
    selected_row = gpu_by_index.get(selected, {})
    busy_uuids = {
        row["gpu_uuid"]
        for row in record.get("resources_before", {}).get("processes", [])
    }
    common_valid = all(
        run.get("fold") == 1
        and run.get("seed") == 20260714
        and run.get("loss") == "paper"
        and run.get("steps") == 2
        and run.get("sequence_length") == 131072
        and run.get("n_tracks") == 241
        and run.get("prediction_shapes")
        == {"1": [1, 241, 131072], "128": [1, 241, 1024]}
        and run.get("target_shapes")
        == {"1": [1, 241, 131072], "128": [1, 241, 1024]}
        and run.get("checkpoint_reload_verified") is True
        and run.get("claim_status") == "environment_validation"
        for run in runs
    )
    numeric_valid = all(
        run.get("metrics")
        and all(
            math.isfinite(float(metric["loss"]))
            and math.isfinite(float(metric["gradient_norm"]))
            and int(metric["cuda_max_memory_allocated_bytes"]) > 0
            for metric in run["metrics"]
        )
        for run in runs
    )
    expected_intervals = sha256(
        REPO_ROOT / "alphagenome_custom/intervals/v2/fold_1/train.tsv"
    )
    checks.extend(
        [
            check(
                "R6A.02_gpu_policy",
                selected in {2, 3}
                and int(selected_row.get("memory_free_mib", 0)) >= 70_000
                and selected_row.get("uuid") not in busy_uuids,
                f"selected={selected} free={selected_row.get('memory_free_mib')} busy={selected_row.get('uuid') in busy_uuids}",
            ),
            check(
                "R6A.03_three_model_smokes",
                record.get("status") == "completed"
                and record.get("split_revision") == "six_chromosome_blocks_v2"
                and record.get("locked_test_block_signal_reads") == 0
                and {run.get("model") for run in runs} == {"A", "B", "C"}
                and len(runs) == 3
                and not run_errors,
                f"models={[run.get('model') for run in runs]} errors={run_errors}",
            ),
            check(
                "R6A.04_common_interface_and_budget",
                common_valid,
                "A/B/C used fold 1, seed 20260714, paper loss, 131072 bp, 241 tracks, and two steps",
            ),
            check(
                "R6A.05_finite_numerics_and_memory",
                numeric_valid,
                "all losses/gradients finite and CUDA memory recorded",
            ),
            check(
                "R6A.06_checkpoint_and_log_integrity",
                not run_errors and len(runs) == 3,
                f"errors={run_errors}",
            ),
            check(
                "R6A.07_test_embargo",
                all(run.get("intervals_sha256") == expected_intervals for run in runs),
                "all smokes bind to six-chromosome fold-1 training blocks; locked test blocks were not read",
            ),
            check(
                "R6A.08_claim_scope",
                all(run.get("run_type") == "smoke" for run in runs),
                "smoke outputs are environment validation, not model results",
            ),
        ]
    )
    return {
        "schema_version": 1,
        "phase": "P6A",
        "review": "R6A",
        "reviewed_at": utc_now(),
        "status": "PASS"
        if all(item["status"] == "PASS" for item in checks)
        else "FAIL",
        "checks": checks,
    }


def review_p6b_original() -> dict[str, Any]:
    metadata_dir = REPO_ROOT / "alphagenome_custom/metadata/v2"
    paths = {
        "spec": metadata_dir / "p6b_matrix_spec.json",
        "execution": metadata_dir / "p6b_execution.json",
        "results": metadata_dir / "p6b_cv_results.tsv",
        "selection": metadata_dir / "p6b_selection.json",
        "lock": metadata_dir / "final_test_lock.json",
        "means": metadata_dir / "track_nonzero_means_v2.tsv",
    }
    missing = [
        str(path.relative_to(REPO_ROOT))
        for path in paths.values()
        if not path.is_file()
    ]
    checks = [check("R6B.01_required_outputs", not missing, f"missing={missing}")]
    if missing:
        return {
            "schema_version": 1,
            "phase": "P6B",
            "review": "R6B",
            "reviewed_at": utc_now(),
            "status": "FAIL",
            "checks": checks,
        }
    spec = json.loads(paths["spec"].read_text())
    execution = json.loads(paths["execution"].read_text())
    results = read_tsv(paths["results"])
    selection = json.loads(paths["selection"].read_text())
    lock = json.loads(paths["lock"].read_text())
    spec_sha = sha256(paths["spec"])
    expected_jobs = {
        (model, loss, fold)
        for model in ("A", "B", "C")
        for loss in ("paper", "log1p_mse")
        for fold in range(1, 6)
    }
    observed_jobs = {
        (row.get("model"), row.get("loss"), int(row.get("fold", -1)))
        for row in execution.get("jobs", [])
    }
    job_errors = []
    for job in execution.get("jobs", []):
        try:
            fold = int(job["fold"])
            run_path = REPO_ROOT / job["run_path"]
            validation_path = REPO_ROOT / job["validation_path"]
            checkpoint_path = REPO_ROOT / job["checkpoint_path"]
            log_path = REPO_ROOT / job["log_path"]
            run = json.loads(run_path.read_text())
            validation = json.loads(validation_path.read_text())
            if (
                job["spec_sha256"] != spec_sha
                or int(job["seed"]) != 20260714
                or int(job["steps"]) != 2000
                or int(job["sequence_length"]) != 131072
                or int(job["physical_gpu"]) not in {2, 3}
            ):
                raise ValueError("job registration mismatch")
            if (
                sha256(checkpoint_path) != job["checkpoint_sha256"]
                or run["checkpoint_sha256"] != job["checkpoint_sha256"]
                or not log_path.is_file()
                or log_path.stat().st_size == 0
            ):
                raise ValueError("checkpoint or log integrity mismatch")
            expected_train = sha256(
                REPO_ROOT / f"alphagenome_custom/intervals/v2/fold_{fold}/train.tsv"
            )
            expected_valid = sha256(
                REPO_ROOT / f"alphagenome_custom/intervals/v2/fold_{fold}/valid.tsv"
            )
            if (
                run.get("model") != job["model"]
                or run.get("loss") != job["loss"]
                or run.get("fold") != fold
                or run.get("seed") != 20260714
                or run.get("steps") != 2000
                or run.get("intervals_sha256") != expected_train
                or run.get("means_sha256") != sha256(paths["means"])
            ):
                raise ValueError("training contract mismatch")
            mean_metrics = validation.get("mean_metrics", {})
            if (
                validation.get("model") != job["model"]
                or validation.get("training_loss") != job["loss"]
                or validation.get("fold") != fold
                or validation.get("checkpoint_sha256") != job["checkpoint_sha256"]
                or validation.get("intervals_sha256") != expected_valid
                or validation.get("chromosome_x_read") is not False
                or float(validation.get("validation_core_coverage_fraction", 0)) < 0.98
                or int(validation.get("validation_subwindows", 0)) <= 0
                or int(validation.get("finite_per_track_pearson_128bp", 0)) <= 0
                or not math.isfinite(float(mean_metrics.get("paper_loss", "nan")))
                or not math.isfinite(float(mean_metrics.get("log1p_mse", "nan")))
                or not math.isfinite(
                    float(validation.get("mean_per_track_pearson_128bp", "nan"))
                )
            ):
                raise ValueError("validation contract mismatch")
            if (
                not math.isclose(
                    float(job["paper_loss"]),
                    float(mean_metrics["paper_loss"]),
                    rel_tol=1e-12,
                )
                or not math.isclose(
                    float(job["log1p_mse"]),
                    float(mean_metrics["log1p_mse"]),
                    rel_tol=1e-12,
                )
                or not math.isclose(
                    float(job["mean_per_track_pearson_128bp"]),
                    float(validation["mean_per_track_pearson_128bp"]),
                    rel_tol=1e-12,
                )
            ):
                raise ValueError("job summary does not match validation record")
        except Exception as error:
            job_errors.append(
                f"{job.get('model')}/{job.get('loss')}/fold{job.get('fold')}:{error}"
            )

    result_keys = {(row["model"], row["loss"]) for row in results}
    result_numeric = all(
        int(row["folds"]) == 5
        and math.isfinite(float(row["mean_five_fold_paper_loss"]))
        and math.isfinite(float(row["mean_five_fold_log1p_mse"]))
        and math.isfinite(float(row["mean_five_fold_per_track_pearson_128bp"]))
        for row in results
    )
    expected_configs = {
        (model, loss) for model in ("A", "B", "C")
        for loss in ("paper", "log1p_mse")
    }
    aggregation_errors = []
    jobs_by_config: dict[tuple[str, str], list[dict[str, Any]]] = {}
    for job in execution.get("jobs", []):
        jobs_by_config.setdefault((job["model"], job["loss"]), []).append(job)
    for row in results:
        members = jobs_by_config.get((row["model"], row["loss"]), [])
        if len(members) != 5:
            aggregation_errors.append(f"{row['model']}/{row['loss']}:members={len(members)}")
            continue
        expected_values = {
            "mean_five_fold_paper_loss": sum(
                float(member["paper_loss"]) for member in members
            )
            / 5,
            "mean_five_fold_log1p_mse": sum(
                float(member["log1p_mse"]) for member in members
            )
            / 5,
            "mean_five_fold_per_track_pearson_128bp": sum(
                float(member["mean_per_track_pearson_128bp"])
                for member in members
            )
            / 5,
        }
        for key, expected_value in expected_values.items():
            if not math.isclose(
                float(row[key]), expected_value, rel_tol=1e-12, abs_tol=1e-12
            ):
                aggregation_errors.append(f"{row['model']}/{row['loss']}:{key}")
    selected_row = min(
        results,
        key=lambda row: (
            float(row["mean_five_fold_paper_loss"]),
            float(row["mean_five_fold_log1p_mse"]),
            -float(row["mean_five_fold_per_track_pearson_128bp"]),
            row["model"],
            row["loss"],
        ),
    )
    selected = selection.get("selected", {})
    selection_valid = (
        selected.get("model") == selected_row["model"]
        and selected.get("loss") == selected_row["loss"]
        and math.isclose(
            float(selected.get("mean_five_fold_paper_loss", "nan")),
            float(selected_row["mean_five_fold_paper_loss"]),
            rel_tol=1e-12,
        )
        and selection.get("cv_results_sha256") == sha256(paths["results"])
        and selection.get("spec_sha256") == spec_sha
        and selection.get("chromosome_x_read") is False
    )
    development_errors = []
    try:
        development = execution["development"]
        run_path = REPO_ROOT / development["run_path"]
        checkpoint_path = REPO_ROOT / development["checkpoint_path"]
        run = json.loads(run_path.read_text())
        expected_development = sha256(
            REPO_ROOT / "alphagenome_custom/intervals/v2/development_train.tsv"
        )
        if (
            development["model"] != selected_row["model"]
            or development["loss"] != selected_row["loss"]
            or int(development["fold"]) != 0
            or int(development["steps"]) != 2500
            or int(development["physical_gpu"]) not in {2, 3}
            or development["spec_sha256"] != spec_sha
            or development["intervals_sha256"] != expected_development
            or sha256(checkpoint_path) != development["checkpoint_sha256"]
            or run.get("fold") != 0
            or run.get("mean_column") != "development_I_V_nonzero_mean"
            or run.get("intervals_sha256") != expected_development
        ):
            raise ValueError("development retraining mismatch")
    except Exception as error:
        development_errors.append(str(error))
    state = json.loads(
        (REPO_ROOT / "alphagenome_custom/metadata/v2/execution_state.json").read_text()
    )
    g5 = state.get("approvals", {}).get("G5_final_test", {})
    lock_valid = False
    try:
        checkpoint_path = REPO_ROOT / lock["checkpoint_path"]
        lock_valid = (
            lock["checkpoint_sha256"] == sha256(checkpoint_path)
            and lock["model"] == selected_row["model"]
            and lock["loss"] == selected_row["loss"]
            and lock["training_chromosomes"] == ["I", "II", "III", "IV", "V"]
            and lock["selection_sha256"] == sha256(paths["selection"])
            and lock["spec_sha256"] == spec_sha
            and lock["legacy_chr_x_prior_exposure_disclosed"] is True
            and lock["test_consumed"] is False
            and lock["test_consumed_at"] is None
        )
    except Exception:
        lock_valid = False
    checks.extend(
        [
            check(
                "R6B.02_preregistered_matrix",
                spec.get("locked_before_p6b") is True
                and spec.get("models") == ["A", "B", "C"]
                and spec.get("losses") == ["paper", "log1p_mse"]
                and spec.get("cv_folds") == [1, 2, 3, 4, 5]
                and spec.get("chromosome_x_access") == "prohibited",
                f"spec_sha256={spec_sha}",
            ),
            check(
                "R6B.03_complete_fair_matrix",
                execution.get("status") == "completed"
                and execution.get("jobs_expected") == 30
                and observed_jobs == expected_jobs
                and not execution.get("failures")
                and not job_errors,
                f"jobs={len(observed_jobs)}/30 errors={job_errors[:5]}",
            ),
            check(
                "R6B.04_gpu_policy",
                set(execution.get("selected_physical_gpus", [])).issubset({2, 3})
                and bool(execution.get("selected_physical_gpus")),
                f"selected={execution.get('selected_physical_gpus')}",
            ),
            check(
                "R6B.05_cv_aggregation",
                len(results) == 6
                and result_keys == expected_configs
                and result_numeric
                and not aggregation_errors,
                f"configs={sorted(result_keys)} errors={aggregation_errors}",
            ),
            check(
                "R6B.06_locked_selection_rule",
                selection_valid,
                f"selected={selected.get('model')}/{selected.get('loss')}",
            ),
            check(
                "R6B.07_development_retrain",
                not development_errors,
                f"errors={development_errors}",
            ),
            check(
                "R6B.08_final_checkpoint_lock",
                lock_valid and execution.get("lock_sha256") == sha256(paths["lock"]),
                f"checkpoint={lock.get('checkpoint_path')}",
            ),
            check(
                "R6B.09_test_embargo",
                g5.get("approved") is not True
                and all(
                    json.loads((REPO_ROOT / job["validation_path"]).read_text()).get(
                        "chromosome_x_read"
                    )
                    is False
                    for job in execution.get("jobs", [])
                ),
                "chromosome X remained blocked; G5 is unapproved",
            ),
        ]
    )
    return {
        "schema_version": 1,
        "phase": "P6B",
        "review": "R6B",
        "reviewed_at": utc_now(),
        "status": "PASS"
        if all(item["status"] == "PASS" for item in checks)
        else "FAIL",
        "checks": checks,
    }


def review_p6b_chromosome_holdout() -> dict[str, Any]:
    metadata_dir = REPO_ROOT / "alphagenome_custom/metadata/v2"
    paths = {
        "spec": metadata_dir / "p6b_amendment_spec.json",
        "execution": metadata_dir / "p6b_amendment_execution.json",
        "results": metadata_dir / "p6b_amendment_cv_results.tsv",
        "ablations": metadata_dir / "p6b_amendment_ablation_results.tsv",
        "selection": metadata_dir / "p6b_amendment_selection.json",
        "legacy": metadata_dir / "p6b_legacy_v1_disposition.json",
        "lock": metadata_dir / "final_test_lock.json",
        "old_lock": metadata_dir / "final_test_lock_original_p6b.json",
        "means": metadata_dir / "track_nonzero_means_v2.tsv",
    }
    missing = [str(path.relative_to(REPO_ROOT)) for path in paths.values() if not path.is_file()]
    checks = [check("R6B.A1_required_outputs", not missing, f"missing={missing}")]
    if missing:
        return {
            "schema_version": 2, "phase": "P6B", "review": "R6B-amended",
            "reviewed_at": utc_now(), "status": "FAIL", "checks": checks,
        }
    spec = json.loads(paths["spec"].read_text())
    execution = json.loads(paths["execution"].read_text())
    results = read_tsv(paths["results"])
    ablations = read_tsv(paths["ablations"])
    selection = json.loads(paths["selection"].read_text())
    legacy = json.loads(paths["legacy"].read_text())
    final_lock = json.loads(paths["lock"].read_text())
    old_lock = json.loads(paths["old_lock"].read_text())
    state = json.loads((metadata_dir / "execution_state.json").read_text())
    spec_sha = sha256(paths["spec"])
    formal_expected = {
        (model, loss, seed, fold)
        for model in ("A", "B", "C")
        for loss in ("paper", "log1p_mse")
        for seed in (20260714, 20260715, 20260716)
        for fold in range(1, 6)
    }
    ablation_expected = {
        (ablation, fold)
        for ablation in ("no_augmentation", "no_gene_loss", "whole_I_V_mean")
        for fold in range(1, 6)
    }
    formal_jobs = [job for job in execution.get("jobs", []) if job.get("job_kind") == "formal"]
    ablation_jobs = [job for job in execution.get("jobs", []) if job.get("job_kind") == "ablation"]
    formal_observed = {
        (job.get("model"), job.get("loss"), int(job.get("seed", -1)), int(job.get("fold", -1)))
        for job in formal_jobs
    }
    ablation_observed = {
        (job.get("ablation_id"), int(job.get("fold", -1))) for job in ablation_jobs
    }
    metric_names = (
        "primary_biological_score", "gene_exon_coverage_pearson_log1p",
        "per_track_pearson_128bp_log1p", "paper_loss", "log1p_mse",
        "spearman_128bp", "top1_mse_128bp", "top1_calibration_ratio_128bp",
        "gene_body_mse_log1p_1bp", "exon_mse_log1p_1bp",
        "local_gradient_mse_log1p_1bp",
    )
    job_errors = []
    chromosome_x_reads = 0
    for job in formal_jobs + ablation_jobs:
        try:
            checkpoint = REPO_ROOT / job["checkpoint_path"]
            validation_path = REPO_ROOT / job["validation_path"]
            job_path = REPO_ROOT / job["job_path"]
            validation = json.loads(validation_path.read_text())
            persisted = json.loads(job_path.read_text())
            fold = int(job["fold"])
            if (
                persisted.get("spec_sha256") != spec_sha
                or persisted.get("checkpoint_sha256") != sha256(checkpoint)
                or validation.get("checkpoint_sha256") != persisted.get("checkpoint_sha256")
                or validation.get("schema_version") != 2
                or validation.get("full_metrics", {}).get("metric_contract") != "v2_full_metrics_1"
                or validation.get("intervals_sha256") != sha256(
                    REPO_ROOT / f"alphagenome_custom/intervals/v2/fold_{fold}/valid.tsv"
                )
                or float(validation.get("validation_core_coverage_fraction", 0)) < 0.98
            ):
                raise ValueError("artifact or validation contract mismatch")
            chromosome_x_reads += int(validation.get("chromosome_x_read") is True)
            if not all(math.isfinite(float(job[name])) for name in metric_names):
                raise ValueError("non-finite summary metric")
            primary = validation["full_metrics"]["primary"]
            expected_score = 0.5 * float(
                primary["mean_per_track_gene_exon_coverage_pearson_log1p"]
            ) + 0.5 * float(primary["mean_per_track_pearson_128bp_log1p"])
            if not math.isclose(
                float(job["primary_biological_score"]), expected_score,
                rel_tol=1e-12, abs_tol=1e-12,
            ):
                raise ValueError("primary score mismatch")
        except Exception as error:
            job_errors.append(f"{job.get('job_id')}:{error}")

    result_keys = {(row["model"], row["loss"]) for row in results}
    expected_result_keys = {(model, loss) for model in ("A", "B", "C") for loss in ("paper", "log1p_mse")}
    result_contract = all(
        int(row["jobs"]) == 15 and int(row["folds"]) == 5 and int(row["seeds"]) == 3
        and all(math.isfinite(float(row[f"mean_{name}"])) for name in metric_names)
        for row in results
    )
    selected_row = max(
        results,
        key=lambda row: (
            float(row["mean_primary_biological_score"]),
            float(row["mean_gene_exon_coverage_pearson_log1p"]),
            float(row["mean_per_track_pearson_128bp_log1p"]),
            -float(row["mean_paper_loss"]), row["model"], row["loss"],
        ),
    )
    selected = selection.get("selected", {})
    selection_valid = (
        selected.get("model") == selected_row["model"]
        and selected.get("loss") == selected_row["loss"]
        and selection.get("amendment_spec_sha256") == spec_sha
        and selection.get("results_sha256") == sha256(paths["results"])
        and selection.get("chromosome_x_read") is False
        and math.isclose(
            float(selected.get("mean_primary_biological_score", "nan")),
            float(selected_row["mean_primary_biological_score"]), rel_tol=1e-12,
        )
    )
    ablation_contract = (
        {(row["ablation_id"], row["model"], row["loss"]) for row in ablations}
        == {
            ("no_augmentation", "B", "paper"),
            ("no_gene_loss", "B", "paper"),
            ("whole_I_V_mean", "A", "paper"),
        }
        and all(int(row["jobs"]) == 5 and int(row["folds"]) == 5 for row in ablations)
    )
    development = execution.get("development", {})
    lock_valid = False
    try:
        checkpoint = REPO_ROOT / final_lock["checkpoint_path"]
        lock_valid = (
            final_lock.get("schema_version") == 2
            and final_lock.get("superseded") is False
            and final_lock.get("test_consumed") is False
            and final_lock.get("test_consumed_at") is None
            and final_lock.get("model") == selected_row["model"]
            and final_lock.get("loss") == selected_row["loss"]
            and final_lock.get("checkpoint_sha256") == sha256(checkpoint)
            and final_lock.get("checkpoint_sha256") == development.get("checkpoint_sha256")
            and final_lock.get("seed") == 20260717
            and final_lock.get("amendment_spec_sha256") == spec_sha
            and final_lock.get("selection_sha256") == sha256(paths["selection"])
            and final_lock.get("legacy_chr_x_prior_exposure_disclosed") is True
        )
    except Exception:
        lock_valid = False
    g5 = state.get("approvals", {}).get("G5_final_test", {})
    checks.extend([
        check(
            "R6B.A2_preregistered_amendment",
            spec.get("locked_before_comparative_full_metric_backfill") is True
            and spec.get("chromosome_x_access") == "prohibited"
            and spec.get("formal_matrix", {}).get("total_jobs") == 90
            and spec.get("formal_matrix", {}).get("seeds") == [20260714, 20260715, 20260716]
            and spec.get("metrics", {}).get("contract") == "v2_full_metrics_1",
            f"spec_sha256={spec_sha}",
        ),
        check(
            "R6B.A3_complete_three_seed_matrix",
            execution.get("status") == "completed" and execution.get("jobs_expected") == 105
            and formal_observed == formal_expected and not execution.get("failures") and not job_errors,
            f"formal={len(formal_observed)}/90 errors={job_errors[:5]}",
        ),
        check(
            "R6B.A4_complete_ablations",
            ablation_observed == ablation_expected and ablation_contract,
            f"ablations={len(ablation_observed)}/15 configs={[row.get('ablation_id') for row in ablations]}",
        ),
        check(
            "R6B.A5_full_metric_aggregation",
            result_keys == expected_result_keys and len(results) == 6 and result_contract,
            f"configs={sorted(result_keys)}",
        ),
        check(
            "R6B.A6_biological_primary_selection", selection_valid,
            f"selected={selected.get('model')}/{selected.get('loss')}",
        ),
        check(
            "R6B.A7_legacy_disposition",
            legacy.get("decision") == "not_formally_comparable_as_a_label_scale_ablation"
            and legacy.get("formal_promotion_use") is False
            and legacy.get("replacement_control")
            == "Use the registered fixed-A whole-I-V-mean versus fold-train-mean ablation on the same 241 normalized v2 tracks; retain legacy_v1 results as exploratory historical evidence only.",
            legacy.get("decision", "missing"),
        ),
        check(
            "R6B.A8_development_and_final_lock", lock_valid,
            f"checkpoint={final_lock.get('checkpoint_path')}",
        ),
        check(
            "R6B.A9_old_lock_preserved",
            old_lock.get("superseded") is True and old_lock.get("test_consumed") is False,
            f"old_checkpoint={old_lock.get('checkpoint_path')}",
        ),
        check(
            "R6B.A10_test_embargo",
            chromosome_x_reads == 0 and g5.get("approved") is not True,
            f"chromosome_x_reads={chromosome_x_reads}; G5 unapproved",
        ),
    ])
    return {
        "schema_version": 2, "phase": "P6B", "review": "R6B-amended",
        "reviewed_at": utc_now(),
        "status": "PASS" if all(item["status"] == "PASS" for item in checks) else "FAIL",
        "checks": checks,
    }


def review_p6b() -> dict[str, Any]:
    metadata_dir = REPO_ROOT / "alphagenome_custom/metadata/v2"
    paths = {
        "spec": metadata_dir / "p6b_six_chromosome_v2_spec.json",
        "execution": metadata_dir / "p6b_six_chromosome_v2_execution.json",
        "results": metadata_dir / "p6b_six_chromosome_v2_cv_results.tsv",
        "ablations": metadata_dir / "p6b_six_chromosome_v2_ablation_results.tsv",
        "selection": metadata_dir / "p6b_six_chromosome_v2_selection.json",
        "lock": metadata_dir / "final_test_lock.json",
        "means": metadata_dir / "track_nonzero_means_v2.tsv",
        "registry": metadata_dir / "split_registry_v2.json",
        "test": REPO_ROOT / "alphagenome_custom/intervals/v2/test_locked.tsv",
        "migration": metadata_dir / "six_chromosome_revision_migration.json",
        "core_migration": metadata_dir / "six_chromosome_core_coverage_revision.json",
        "invalid_execution": metadata_dir
        / "six_chromosome_blocks_v1_incomplete_core/p6b_six_chromosome_execution.json",
        "old_lock": metadata_dir
        / "chromosome_holdout_original/final_test_lock.json",
    }
    missing = [
        str(path.relative_to(REPO_ROOT))
        for path in paths.values()
        if not path.is_file()
    ]
    checks = [check("R6B.S1_required_outputs", not missing, f"missing={missing}")]
    if missing:
        return {
            "schema_version": 3,
            "phase": "P6B",
            "review": "R6B-six-chromosome",
            "reviewed_at": utc_now(),
            "status": "FAIL",
            "checks": checks,
        }
    spec = json.loads(paths["spec"].read_text())
    execution = json.loads(paths["execution"].read_text())
    results = read_tsv(paths["results"])
    ablations = read_tsv(paths["ablations"])
    selection = json.loads(paths["selection"].read_text())
    lock = json.loads(paths["lock"].read_text())
    registry = json.loads(paths["registry"].read_text())
    migration = json.loads(paths["migration"].read_text())
    core_migration = json.loads(paths["core_migration"].read_text())
    invalid_execution = json.loads(paths["invalid_execution"].read_text())
    old_lock = json.loads(paths["old_lock"].read_text())
    state = json.loads((metadata_dir / "execution_state.json").read_text())
    spec_sha = sha256(paths["spec"])
    chromosomes = {"I", "II", "III", "IV", "V", "X"}
    formal_expected = {
        (model, loss, seed, fold)
        for model in ("A", "B", "C")
        for loss in ("paper", "log1p_mse")
        for seed in (20260714, 20260715, 20260716)
        for fold in range(1, 6)
    }
    ablation_expected = {
        (ablation, fold)
        for ablation in (
            "no_augmentation",
            "no_gene_loss",
            "development_pool_mean",
        )
        for fold in range(1, 6)
    }
    jobs = execution.get("jobs", [])
    formal_jobs = [job for job in jobs if job.get("job_kind") == "formal"]
    ablation_jobs = [job for job in jobs if job.get("job_kind") == "ablation"]
    formal_observed = {
        (job.get("model"), job.get("loss"), job.get("seed"), job.get("fold"))
        for job in formal_jobs
    }
    ablation_observed = {
        (job.get("ablation_id"), job.get("fold")) for job in ablation_jobs
    }
    metric_names = (
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
    job_errors = []
    locked_reads = 0
    for job in jobs:
        try:
            checkpoint = REPO_ROOT / job["checkpoint_path"]
            validation_path = REPO_ROOT / job["validation_path"]
            persisted_path = REPO_ROOT / job["job_path"]
            validation = json.loads(validation_path.read_text())
            persisted = json.loads(persisted_path.read_text())
            fold = int(job["fold"])
            locked_reads += int(
                validation.get("locked_test_block_signal_reads") is True
            )
            if (
                persisted.get("amendment_spec_sha256") != spec_sha
                or persisted.get("checkpoint_sha256") != sha256(checkpoint)
                or validation.get("checkpoint_sha256")
                != persisted.get("checkpoint_sha256")
                or validation.get("schema_version") != 2
                or validation.get("full_metrics", {}).get("metric_contract")
                != "v2_full_metrics_1"
                or validation.get("intervals_sha256")
                != sha256(
                    REPO_ROOT
                    / f"alphagenome_custom/intervals/v2/fold_{fold}/valid.tsv"
                )
                or set(validation.get("evaluated_chromosomes", [])) != chromosomes
                or validation.get("locked_test_block_signal_reads") is not False
                or validation.get("chromosome_x_train_valid_blocks_read") is not True
                or float(validation.get("validation_core_coverage_fraction", 0))
                < 0.999
                or not all(math.isfinite(float(job[name])) for name in metric_names)
            ):
                raise ValueError("artifact, metric, split, or embargo mismatch")
            run = json.loads((REPO_ROOT / job["output_dir"] / "run.json").read_text())
            if (
                run.get("model") != job.get("model")
                or run.get("loss") != job.get("loss")
                or int(run.get("fold", -1)) != fold
                or int(run.get("seed", -1)) != int(job.get("seed", -2))
                or run.get("mean_column") != job.get("mean_column")
                or run.get("intervals_sha256")
                != sha256(
                    REPO_ROOT
                    / f"alphagenome_custom/intervals/v2/fold_{fold}/train.tsv"
                )
                or run.get("means_sha256") != sha256(paths["means"])
                or not (REPO_ROOT / job["log_path"]).is_file()
                or (REPO_ROOT / job["log_path"]).stat().st_size == 0
            ):
                raise ValueError("training contract mismatch")
            primary = validation["full_metrics"]["primary"]
            expected_primary = 0.5 * float(
                primary["mean_per_track_gene_exon_coverage_pearson_log1p"]
            ) + 0.5 * float(primary["mean_per_track_pearson_128bp_log1p"])
            if not math.isclose(
                float(job["primary_biological_score"]),
                expected_primary,
                rel_tol=1e-12,
                abs_tol=1e-12,
            ):
                raise ValueError("primary score mismatch")
        except Exception as error:
            job_errors.append(f"{job.get('job_id')}:{error}")

    expected_configs = {
        (model, loss) for model in ("A", "B", "C") for loss in ("paper", "log1p_mse")
    }
    result_keys = {(row["model"], row["loss"]) for row in results}
    result_contract = all(
        int(row["jobs"]) == 15
        and int(row["folds"]) == 5
        and int(row["seeds"]) == 3
        and all(math.isfinite(float(row[f"mean_{name}"])) for name in metric_names)
        for row in results
    )
    aggregation_errors = []
    for rows, members_pool, keys in (
        (results, formal_jobs, ("model", "loss")),
        (ablations, ablation_jobs, ("ablation_id", "model", "loss")),
    ):
        for row in rows:
            members = [
                job
                for job in members_pool
                if all(job.get(key) == row[key] for key in keys)
            ]
            if not members:
                aggregation_errors.append(f"missing:{'/'.join(str(row[key]) for key in keys)}")
                continue
            for metric in metric_names:
                expected = sum(float(job[metric]) for job in members) / len(members)
                if not math.isclose(
                    float(row[f"mean_{metric}"]),
                    expected,
                    rel_tol=1e-12,
                    abs_tol=1e-12,
                ):
                    aggregation_errors.append(
                        f"{'/'.join(str(row[key]) for key in keys)}:{metric}"
                    )
    selected_row = max(
        results,
        key=lambda row: (
            float(row["mean_primary_biological_score"]),
            float(row["mean_gene_exon_coverage_pearson_log1p"]),
            float(row["mean_per_track_pearson_128bp_log1p"]),
            -float(row["mean_paper_loss"]),
            row["model"],
            row["loss"],
        ),
    )
    selected = selection.get("selected", {})
    selection_valid = (
        selected.get("model") == selected_row["model"]
        and selected.get("loss") == selected_row["loss"]
        and selection.get("spec_sha256") == spec_sha
        and selection.get("results_sha256") == sha256(paths["results"])
        and selection.get("locked_test_block_signal_reads") == 0
        and selection.get("selection_rule") == spec.get("selection")
        and math.isclose(
            float(selected.get("mean_primary_biological_score", "nan")),
            float(selected_row["mean_primary_biological_score"]),
            rel_tol=1e-12,
        )
    )
    ablation_contract = (
        {(row["ablation_id"], row["model"], row["loss"]) for row in ablations}
        == {
            ("no_augmentation", "B", "paper"),
            ("no_gene_loss", "B", "paper"),
            ("development_pool_mean", "A", "paper"),
        }
        and all(
            int(row["jobs"]) == 5 and int(row["folds"]) == 5
            for row in ablations
        )
    )
    development = execution.get("development", {})
    development_valid = False
    lock_valid = False
    try:
        development_run = json.loads(
            (REPO_ROOT / development["run_path"]).read_text()
        )
        development_checkpoint = REPO_ROOT / development["checkpoint_path"]
        expected_development = sha256(
            REPO_ROOT / "alphagenome_custom/intervals/v2/development_train.tsv"
        )
        development_valid = (
            development.get("model") == selected_row["model"]
            and development.get("loss") == selected_row["loss"]
            and development.get("seed") == 20260722
            and development.get("steps") == 2500
            and int(development.get("physical_gpu", -1)) in {2, 3}
            and development.get("spec_sha256") == spec_sha
            and development.get("intervals_sha256") == expected_development
            and development_run.get("model") == selected_row["model"]
            and development_run.get("loss") == selected_row["loss"]
            and development_run.get("fold") == 0
            and development_run.get("seed") == 20260722
            and development_run.get("mean_column")
            == "development_train_nonzero_mean"
            and development_run.get("intervals_sha256") == expected_development
            and development_run.get("means_sha256") == sha256(paths["means"])
            and sha256(development_checkpoint)
            == development.get("checkpoint_sha256")
        )
        lock_checkpoint = REPO_ROOT / lock["checkpoint_path"]
        lock_valid = (
            lock.get("schema_version") == 3
            and lock.get("superseded") is False
            and lock.get("test_consumed") is False
            and lock.get("test_consumed_at") is None
            and lock.get("model") == selected_row["model"]
            and lock.get("loss") == selected_row["loss"]
            and lock.get("checkpoint_sha256") == sha256(lock_checkpoint)
            and lock.get("checkpoint_sha256")
            == development.get("checkpoint_sha256")
            and lock.get("selection_sha256") == sha256(paths["selection"])
            and lock.get("spec_sha256") == spec_sha
            and lock.get("split_revision") == "six_chromosome_blocks_v2"
            and lock.get("split_registry_sha256") == sha256(paths["registry"])
            and lock.get("means_sha256") == sha256(paths["means"])
            and lock.get("test_intervals_sha256") == sha256(paths["test"])
            and set(lock.get("training_chromosomes", [])) == chromosomes
            and set(lock.get("final_test_chromosomes", [])) == chromosomes
            and lock.get("final_test_scope")
            == "r6c_single_six_chromosome_block_test"
        )
    except Exception:
        development_valid = False
        lock_valid = False
    g5 = state.get("approvals", {}).get("G5_final_test", {})
    checks.extend(
        [
            check(
                "R6B.S2_preregistered_revision",
                spec.get("locked_before_training") is True
                and spec.get("revision_id") == "P6B-SIXCHR-2"
                and spec.get("formal_matrix", {}).get("total_jobs") == 90
                and spec.get("metrics", {}).get("contract") == "v2_full_metrics_1"
                and spec.get("metrics", {}).get(
                    "minimum_validation_core_coverage_fraction"
                )
                == 0.999
                and spec.get("split_revision") == "six_chromosome_blocks_v2"
                and spec.get("locked_test_block_access") == "prohibited"
                and registry.get("revision_id") == spec.get("split_revision")
                and execution.get("revision_id") == spec.get("revision_id")
                and execution.get("split_registry_sha256") == sha256(paths["registry"])
                and execution.get("means_sha256") == sha256(paths["means"]),
                f"spec_sha256={spec_sha}",
            ),
            check(
                "R6B.S3_complete_matrix",
                execution.get("status") == "completed"
                and execution.get("jobs_expected") == 105
                and formal_observed == formal_expected
                and ablation_observed == ablation_expected
                and not execution.get("failures")
                and not job_errors,
                f"formal={len(formal_observed)}/90 ablations={len(ablation_observed)}/15 errors={job_errors[:5]}",
            ),
            check(
                "R6B.S4_gpu_policy",
                bool(execution.get("selected_physical_gpus"))
                and set(execution.get("selected_physical_gpus", [])).issubset({2, 3}),
                f"selected={execution.get('selected_physical_gpus')}",
            ),
            check(
                "R6B.S5_metrics_and_aggregation",
                result_keys == expected_configs
                and len(results) == 6
                and result_contract
                and not aggregation_errors
                and ablation_contract,
                f"configs={sorted(result_keys)} aggregation_errors={aggregation_errors[:5]}",
            ),
            check(
                "R6B.S6_selection_and_development",
                selection_valid and development_valid,
                f"selected={selected.get('model')}/{selected.get('loss')}",
            ),
            check(
                "R6B.S7_final_lock",
                lock_valid and execution.get("lock_sha256") == sha256(paths["lock"]),
                f"checkpoint={lock.get('checkpoint_path')}",
            ),
            check(
                "R6B.S8_test_embargo",
                locked_reads == 0
                and execution.get("locked_test_block_signal_reads") == 0
                and g5.get("approved") is not True,
                f"locked_reads={locked_reads}; G5 unapproved",
            ),
            check(
                "R6B.S9_old_evidence_preserved",
                old_lock.get("superseded") is True
                and old_lock.get("test_consumed") is False
                and old_lock.get("superseded_for_revision")
                == "six_chromosome_blocks_v1"
                and migration.get("old_split_disposition")
                == "superseded_for_six_chromosome_split",
                f"old_checkpoint={old_lock.get('checkpoint_path')}",
            ),
            check(
                "R6B.S10_incomplete_core_attempt_excluded",
                core_migration.get("source_revision")
                == "six_chromosome_blocks_v1"
                and core_migration.get("target_revision")
                == "six_chromosome_blocks_v2"
                and core_migration.get("source_completed_jobs") == 1
                and core_migration.get("source_locked_test_block_signal_reads") == 0
                and min(
                    (
                        float(value)
                        for value in core_migration.get(
                            "observed_validation_core_coverage_fractions", []
                        )
                    ),
                    default=1.0,
                )
                < 0.999
                and invalid_execution.get("status") == "running"
                and spec.get("invalid_v1_jobs_reusable") is False,
                f"source_status={core_migration.get('source_execution_status')} coverage={core_migration.get('observed_validation_core_coverage_fractions')}",
            ),
        ]
    )
    return {
        "schema_version": 3,
        "phase": "P6B",
        "review": "R6B-six-chromosome",
        "reviewed_at": utc_now(),
        "status": "PASS"
        if all(item["status"] == "PASS" for item in checks)
        else "FAIL",
        "checks": checks,
    }


def review_p6c() -> dict[str, Any]:
    metadata_dir = REPO_ROOT / "alphagenome_custom/metadata/v2"
    paths = {
        "execution": metadata_dir / "p6c_execution.json",
        "claim": metadata_dir / "final_test_claim.json",
        "report": metadata_dir / "final_test_report.json",
        "summary": metadata_dir / "v2_final_report.json",
        "document": REPO_ROOT / "docs/v2_final_report.md",
        "lock": metadata_dir / "final_test_lock.json",
        "split_registry": metadata_dir / "split_registry_v2.json",
        "means": metadata_dir / "track_nonzero_means_v2.tsv",
        "means_summary": metadata_dir / "track_nonzero_means_v2_summary.json",
        "test_intervals": REPO_ROOT
        / "alphagenome_custom/intervals/v2/test_locked.tsv",
    }
    missing = [
        str(path.relative_to(REPO_ROOT))
        for path in paths.values()
        if not path.is_file()
    ]
    checks = [check("R6C.01_required_outputs", not missing, f"missing={missing}")]
    if missing:
        return {
            "schema_version": 1,
            "phase": "P6C",
            "review": "R6C",
            "reviewed_at": utc_now(),
            "status": "FAIL",
            "checks": checks,
        }
    execution = json.loads(paths["execution"].read_text())
    claim = json.loads(paths["claim"].read_text())
    report = json.loads(paths["report"].read_text())
    summary = json.loads(paths["summary"].read_text())
    lock = json.loads(paths["lock"].read_text())
    try:
        selection_path = repository_file(lock.get("selection_path"))
    except (FileNotFoundError, ValueError):
        checks.append(
            check(
                "R6C.02_locked_selection_path",
                False,
                f"invalid={lock.get('selection_path')}",
            )
        )
        return {
            "schema_version": 1,
            "phase": "P6C",
            "review": "R6C",
            "reviewed_at": utc_now(),
            "status": "FAIL",
            "checks": checks,
        }
    selection = json.loads(selection_path.read_text())
    split_registry = json.loads(paths["split_registry"].read_text())
    means_summary = json.loads(paths["means_summary"].read_text())
    state = json.loads(
        (REPO_ROOT / "alphagenome_custom/metadata/v2/execution_state.json").read_text()
    )
    g5 = state.get("approvals", {}).get("G5_final_test", {})
    log_path = REPO_ROOT / execution.get("log_path", "missing")
    mean_metrics = report.get("mean_metrics", {})
    full_metrics = report.get("full_metrics", {})
    full_primary = full_metrics.get("primary", {})
    distribution = full_metrics.get("distribution_128bp", {})
    execution_id = execution.get("execution_id")
    relative_test = str(paths["test_intervals"].relative_to(REPO_ROOT))
    registered_test_sha = split_registry.get("files", {}).get(relative_test)
    current_test_sha = sha256(paths["test_intervals"])
    current_means_sha = sha256(paths["means"])
    split_contract = (
        split_registry.get("revision_id") == "six_chromosome_blocks_v2"
        and split_registry.get("cv_chromosomes")
        == ["I", "II", "III", "IV", "V", "X"]
        and split_registry.get("test_chromosomes")
        == ["I", "II", "III", "IV", "V", "X"]
        and split_registry.get("final_test_scope")
        == "r6c_single_six_chromosome_block_test"
    )
    full_metric_sections = {
        "raw_1bp",
        "raw_128bp_sum",
        "log1p_1bp",
        "log1p_128bp_sum",
        "gene_body_log1p_1bp",
        "exon_log1p_1bp",
        "local_gradient_log1p_1bp",
        "local_gradient_log1p_128bp",
        "gene_exon_coverage",
        "distribution_128bp",
    }
    full_metrics_valid = (
        full_metrics.get("metric_contract") == "v2_full_metrics_1"
        and full_metric_sections.issubset(full_metrics)
        and math.isfinite(
            float(
                full_primary.get(
                    "mean_per_track_gene_exon_coverage_pearson_log1p", "nan"
                )
            )
        )
        and math.isfinite(
            float(full_primary.get("mean_per_track_pearson_128bp_log1p", "nan"))
        )
        and math.isfinite(float(distribution.get("mean_per_track_spearman", "nan")))
        and math.isfinite(float(distribution.get("mean_per_track_top1_mse", "nan")))
        and math.isfinite(
            float(distribution.get("mean_per_track_top1_calibration_ratio", "nan"))
        )
    )
    summary_data = summary.get("data_lineage", {})
    summary_validation = summary.get("formal_validation", {})
    summary_final = summary.get("final_test", {})
    summary_valid = (
        summary.get("status") == "completed_one_time_final_test"
        and summary.get("execution_id") == execution_id
        and summary_data.get("input_samples") == 485
        and summary_data.get("rna_seq_runs_reprocessed") == 482
        and summary_data.get("non_rna_runs_excluded") == 3
        and summary_data.get("grouped_tracks") == 241
        and summary_data.get("monolithic_npz_generated") is False
        and summary_data.get("cv_chromosomes")
        == ["I", "II", "III", "IV", "V", "X"]
        and summary_data.get("test_chromosomes")
        == ["I", "II", "III", "IV", "V", "X"]
        and len(summary_validation.get("matrix", [])) == 6
        and len(summary_validation.get("ablations", [])) == 3
        and summary_validation.get("formal_jobs") == 90
        and summary_validation.get("ablation_jobs") == 15
        and isinstance(summary_validation.get("registered_failures"), list)
        and summary_validation.get("selected") == selection.get("selected")
        and summary_final.get("report_sha256") == sha256(paths["report"])
        and summary_final.get("checkpoint_sha256") == lock.get("checkpoint_sha256")
        and summary_final.get("test_consumed_once") is True
        and summary_final.get("no_post_test_tuning") is True
        and summary_final.get("primary_metrics") == full_primary
        and len(summary.get("limitations", [])) >= 5
        and summary.get("reproduction", {}).get("phase_command")
        == execution.get("command")
    )
    report_valid = (
        report.get("phase") == "P6C"
        and report.get("fold") == 0
        and report.get("model") == lock.get("model")
        and report.get("training_loss") == lock.get("loss")
        and report.get("checkpoint_sha256") == lock.get("checkpoint_sha256")
        and report.get("final_test_execution_id") == execution_id
        and report.get("intervals_sha256") == current_test_sha == registered_test_sha
        and split_contract
        and report.get("means_sha256") == current_means_sha
        and current_means_sha == means_summary.get("output_sha256")
        and report.get("evaluated_chromosomes")
        == ["I", "II", "III", "IV", "V", "X"]
        and report.get("locked_test_block_signal_reads") is True
        and report.get("chromosome_x_train_valid_blocks_read") is False
        and report.get("mean_column") == "development_train_nonzero_mean"
        and means_summary.get("locked_test_block_signal_reads") == 0
        and int(report.get("validation_subwindows", 0)) > 0
        and float(report.get("validation_core_coverage_fraction", 0)) >= 0.999
        and math.isfinite(float(mean_metrics.get("paper_loss", "nan")))
        and math.isfinite(float(mean_metrics.get("log1p_mse", "nan")))
        and math.isfinite(float(report.get("mean_per_track_pearson_128bp", "nan")))
        and full_metrics_valid
    )
    lock_valid = (
        lock.get("test_consumed") is True
        and bool(lock.get("test_consumed_at"))
        and lock.get("test_status") == "completed"
        and lock.get("final_report_path")
        == str(paths["report"].relative_to(REPO_ROOT))
        and lock.get("final_report_sha256") == sha256(paths["report"])
        and lock.get("selection_sha256") == sha256(selection_path)
        and lock.get("test_execution_id") == execution_id
        and int(lock.get("test_physical_gpu", -1))
        == int(execution.get("physical_gpu", -2))
        and lock.get("final_summary_path")
        == str(paths["summary"].relative_to(REPO_ROOT))
        and lock.get("final_summary_sha256") == sha256(paths["summary"])
        and lock.get("final_document_path")
        == str(paths["document"].relative_to(REPO_ROOT))
        and lock.get("final_document_sha256") == sha256(paths["document"])
        and lock.get("schema_version") == 3
        and lock.get("split_revision") == "six_chromosome_blocks_v2"
        and lock.get("split_registry_sha256") == sha256(paths["split_registry"])
        and lock.get("means_sha256") == current_means_sha
        and lock.get("test_intervals_sha256") == current_test_sha
        and lock.get("training_chromosomes") == ["I", "II", "III", "IV", "V", "X"]
        and lock.get("final_test_chromosomes")
        == ["I", "II", "III", "IV", "V", "X"]
        and lock.get("final_test_scope")
        == "r6c_single_six_chromosome_block_test"
        and lock.get("legacy_chr_x_prior_exposure_disclosed") is True
    )
    checks.extend(
        [
            check(
                "R6C.02_single_claim",
                claim.get("status") == "completed"
                and claim.get("execution_id") == execution_id
                and claim.get("checkpoint_sha256") == lock.get("checkpoint_sha256")
                and claim.get("selection_sha256") == sha256(selection_path)
                and claim.get("test_intervals_sha256") == current_test_sha
                and claim.get("means_sha256") == current_means_sha
                and int(claim.get("physical_gpu", -1)) in {2, 3},
                f"execution_id={claim.get('execution_id')} gpu={claim.get('physical_gpu')}",
            ),
            check(
                "R6C.03_execution_and_gpu",
                execution.get("status") == "completed"
                and int(execution.get("physical_gpu", -1)) in {2, 3}
                and report.get("physical_cuda_visible_devices")
                == str(execution.get("physical_gpu"))
                and "--final-test" in execution.get("command", [])
                and execution_id in execution.get("command", [])
                and log_path.is_file()
                and log_path.stat().st_size > 0,
                f"execution_id={execution_id} status={execution.get('status')} gpu={execution.get('physical_gpu')}",
            ),
            check(
                "R6C.04_locked_checkpoint_only",
                report.get("checkpoint_sha256") == lock.get("checkpoint_sha256")
                and selection.get("selected", {}).get("model") == lock.get("model")
                and selection.get("selected", {}).get("loss") == lock.get("loss")
                and lock.get("selection_path")
                == str(selection_path.relative_to(REPO_ROOT)),
                f"checkpoint={lock.get('checkpoint_path')} selection={lock.get('selection_path')}",
            ),
            check(
                "R6C.05_final_metrics",
                report_valid,
                f"paper={mean_metrics.get('paper_loss')} log1p={mean_metrics.get('log1p_mse')}",
            ),
            check(
                "R6C.06_permanent_consumption",
                lock_valid,
                f"consumed_at={lock.get('test_consumed_at')} status={lock.get('test_status')}",
            ),
            check(
                "R6C.07_scoped_final_approval",
                g5.get("approved") is True
                and g5.get("scope") == "r6c_single_six_chromosome_block_test",
                f"scope={g5.get('scope')} approved_at={g5.get('updated_at')}",
            ),
            check(
                "R6C.08_prior_exposure_disclosure",
                lock.get("legacy_chr_x_prior_exposure_disclosed") is True,
                "legacy chr X exposure remains disclosed; six-chromosome v2 test was read once after lock",
            ),
            check(
                "R6C.09_final_synthesis",
                summary_valid,
                f"samples={summary_data.get('input_samples')} tracks={summary_data.get('grouped_tracks')} formal_jobs={summary_validation.get('formal_jobs')}",
            ),
            check(
                "R6C.10_auditable_documents",
                paths["document"].stat().st_size > 0
                and claim.get("final_summary_sha256") == sha256(paths["summary"]),
                f"summary={paths['summary'].relative_to(REPO_ROOT)} document={paths['document'].relative_to(REPO_ROOT)}",
            ),
        ]
    )
    return {
        "schema_version": 3,
        "phase": "P6C",
        "review": "R6C",
        "reviewed_at": utc_now(),
        "status": "PASS"
        if all(item["status"] == "PASS" for item in checks)
        else "FAIL",
        "checks": checks,
    }


def review_p7() -> dict[str, Any]:
    """Review the bounded post-completion legacy-head architecture port."""

    metadata_dir = REPO_ROOT / "alphagenome_custom/metadata/v2"
    paths = {
        "spec": metadata_dir / "p7_legacy_residual_fold1_spec.json",
        "execution": metadata_dir / "p7_legacy_residual_fold1_execution.json",
        "comparison": metadata_dir / "p7_legacy_residual_fold1_comparison.json",
    }
    missing = [
        str(path.relative_to(REPO_ROOT))
        for path in paths.values()
        if not path.is_file()
    ]
    checks = [check("R7.01_required_outputs", not missing, f"missing={missing}")]
    if missing:
        return {
            "schema_version": 1,
            "phase": "P7",
            "review": "R7-legacy-residual-fold1",
            "reviewed_at": utc_now(),
            "status": "FAIL",
            "checks": checks,
        }
    spec = json.loads(paths["spec"].read_text())
    execution = json.loads(paths["execution"].read_text())
    comparison = json.loads(paths["comparison"].read_text())
    state = json.loads((metadata_dir / "execution_state.json").read_text())
    spec_sha = sha256(paths["spec"])
    input_errors = []
    for key, relative in spec.get("locked_inputs", {}).items():
        try:
            if sha256(repository_file(relative)) != spec["input_sha256"][key]:
                input_errors.append(key)
        except Exception:
            input_errors.append(key)
    required_execution_paths = (
        "base_run_path",
        "base_checkpoint_path",
        "run_path",
        "checkpoint_path",
        "validation_path",
        "log_path",
    )
    paths_from_execution = {
        key: repository_file(execution.get(key, ""))
        for key in required_execution_paths
        if execution.get(key)
    }
    execution_path_errors = [
        key
        for key in required_execution_paths
        if key not in paths_from_execution or not paths_from_execution[key].is_file()
    ]
    run = {}
    base_run = {}
    validation = {}
    candidate_metrics: dict[str, Any] = {}
    validation_ok = False
    try:
        base_run = json.loads(paths_from_execution["base_run_path"].read_text())
        run = json.loads(paths_from_execution["run_path"].read_text())
        validation = json.loads(paths_from_execution["validation_path"].read_text())
        primary = validation["full_metrics"]["primary"]
        candidate_metrics = comparison["candidate"]
        expected_primary = 0.5 * float(
            primary["mean_per_track_gene_exon_coverage_pearson_log1p"]
        ) + 0.5 * float(primary["mean_per_track_pearson_128bp_log1p"])
        validation_ok = (
            validation.get("schema_version") == 2
            and validation.get("model") == "D"
            and validation.get("training_loss") == "paper"
            and validation.get("fold") == 1
            and validation.get("frozen_base_checkpoint")
            == execution.get("base_checkpoint_path")
            and validation.get("frozen_base_checkpoint_sha256")
            == execution.get("base_checkpoint_sha256")
            and validation.get("locked_test_block_signal_reads") is False
            and validation.get("chromosome_x_train_valid_blocks_read") is True
            and set(validation.get("evaluated_chromosomes", []))
            == {"I", "II", "III", "IV", "V", "X"}
            and float(validation.get("validation_core_coverage_fraction", 0)) >= 0.999
            and int(validation.get("finite_per_track_pearson_128bp", 0)) == 241
            and math.isclose(
                float(candidate_metrics["primary_biological_score"]),
                expected_primary,
                rel_tol=1e-12,
                abs_tol=1e-12,
            )
        )
    except Exception:
        validation_ok = False
    training = spec.get("training", {})
    execution_ok = (
        not execution_path_errors
        and
        execution.get("phase") == "P7"
        and execution.get("status") == "completed"
        and execution.get("model") == "D"
        and execution.get("loss") == "paper"
        and execution.get("fold") == 1
        and execution.get("seed") == 20260714
        and execution.get("final_test_access") == "prohibited"
        and execution.get("locked_test_block_signal_reads") == 0
        and int(execution.get("physical_gpu", -1)) in {2, 3}
        and execution.get("spec_sha256") == spec_sha
        and execution.get("base_checkpoint_sha256")
        == sha256(paths_from_execution["base_checkpoint_path"])
        and execution.get("checkpoint_sha256") == sha256(paths_from_execution["checkpoint_path"])
        and execution.get("validation_sha256") == sha256(paths_from_execution["validation_path"])
        and execution.get("comparison_sha256") == sha256(paths["comparison"])
        and base_run.get("model") == "D_base"
        and base_run.get("loss") == "paper"
        and base_run.get("fold") == 1
        and base_run.get("seed") == 20260714
        and base_run.get("steps") == training.get("base_max_steps") == 2000
        and base_run.get("n_tracks") == 241
        and base_run.get("intervals_sha256") == spec["input_sha256"]["fold_1_train"]
        and base_run.get("means_sha256") == spec["input_sha256"]["means"]
        and run.get("model") == "D"
        and run.get("loss") == "paper"
        and run.get("fold") == 1
        and run.get("seed") == 20260714
        and run.get("steps") == training.get("residual_max_steps") == 1500
        and run.get("n_tracks") == 241
        and run.get("frozen_base_checkpoint")
        == execution.get("base_checkpoint_path")
        and run.get("frozen_base_checkpoint_sha256")
        == execution.get("base_checkpoint_sha256")
        and run.get("intervals_sha256") == spec["input_sha256"]["fold_1_train"]
        and run.get("means_sha256") == spec["input_sha256"]["means"]
        and paths_from_execution["log_path"].stat().st_size > 0
    )
    baseline = comparison.get("baseline", {})
    candidate = comparison.get("candidate", {})
    delta = comparison.get("delta", {})
    baseline_ok = (
        baseline.get("model") == "B"
        and baseline.get("loss") == "paper"
        and baseline.get("fold") == 1
        and baseline.get("seed") == 20260714
        and baseline.get("checkpoint_sha256")
        == spec["input_sha256"]["baseline_checkpoint"]
        and comparison.get("spec_sha256") == spec_sha
        and math.isclose(
            float(comparison.get("primary_biological_score_delta", "nan")),
            float(candidate.get("primary_biological_score", "nan"))
            - float(baseline.get("primary_biological_score", "nan")),
            rel_tol=1e-12,
            abs_tol=1e-12,
        )
        and "primary_biological_score" in delta
    )
    lock = repository_file(spec["locked_inputs"]["final_test_lock"])
    report = repository_file(spec["locked_inputs"]["final_test_report"])
    lock_data = json.loads(lock.read_text())
    checks.extend(
        [
            check("R7.02_locked_inputs", not input_errors, f"errors={input_errors}"),
            check(
                "R7.03_training_contract",
                execution_ok,
                f"gpu={execution.get('physical_gpu')} path_errors={execution_path_errors}",
            ),
            check("R7.04_development_validation", validation_ok, f"coverage={validation.get('validation_core_coverage_fraction')}"),
            check("R7.05_matched_baseline_comparison", baseline_ok, f"delta={comparison.get('primary_biological_score_delta')}"),
            check(
                "R7.06_final_test_preserved",
                lock_data.get("test_consumed") is True
                and lock_data.get("test_status") == "completed"
                and sha256(lock) == spec["input_sha256"]["final_test_lock"]
                and sha256(report) == spec["input_sha256"]["final_test_report"],
                "P6C lock and report hashes unchanged",
            ),
            check(
                "R7.07_controller_scope",
                state.get("current_phase") == "P7"
                and state.get("approvals", {}).get("G4_gpu_experiments", {}).get("scope")
                == "p7_single_legacy_residual_fold1",
                f"phase={state.get('current_phase')}",
            ),
        ]
    )
    return {
        "schema_version": 1,
        "phase": "P7",
        "review": "R7-legacy-residual-fold1",
        "reviewed_at": utc_now(),
        "status": "PASS"
        if all(item["status"] == "PASS" for item in checks)
        else "FAIL",
        "checks": checks,
    }


def review_p8() -> dict[str, Any]:
    """Review the bounded B-noLoRA development ablation."""

    metadata_dir = REPO_ROOT / "alphagenome_custom/metadata/v2"
    paths = {
        "spec": metadata_dir / "p8_b_no_lora_fold1_spec.json",
        "execution": metadata_dir / "p8_b_no_lora_fold1_execution.json",
        "comparison": metadata_dir / "p8_b_no_lora_fold1_comparison.json",
    }
    missing = [
        str(path.relative_to(REPO_ROOT))
        for path in paths.values()
        if not path.is_file()
    ]
    checks = [check("R8.01_required_outputs", not missing, f"missing={missing}")]
    if missing:
        return {
            "schema_version": 1,
            "phase": "P8",
            "review": "R8-b-no-lora-fold1",
            "reviewed_at": utc_now(),
            "status": "FAIL",
            "checks": checks,
        }
    spec = json.loads(paths["spec"].read_text())
    execution = json.loads(paths["execution"].read_text())
    comparison = json.loads(paths["comparison"].read_text())
    state = json.loads((metadata_dir / "execution_state.json").read_text())
    spec_sha = sha256(paths["spec"])
    input_errors = []
    for key, relative in spec.get("locked_inputs", {}).items():
        try:
            if sha256(repository_file(relative)) != spec["input_sha256"][key]:
                input_errors.append(key)
        except Exception:
            input_errors.append(key)
    required_execution_paths = ("run_path", "checkpoint_path", "validation_path", "log_path")
    paths_from_execution = {
        key: repository_file(execution.get(key, ""))
        for key in required_execution_paths
        if execution.get(key)
    }
    execution_path_errors = [
        key
        for key in required_execution_paths
        if key not in paths_from_execution or not paths_from_execution[key].is_file()
    ]
    run: dict[str, Any] = {}
    validation: dict[str, Any] = {}
    candidate_metrics: dict[str, Any] = {}
    validation_ok = False
    try:
        run = json.loads(paths_from_execution["run_path"].read_text())
        validation = json.loads(paths_from_execution["validation_path"].read_text())
        primary = validation["full_metrics"]["primary"]
        candidate_metrics = comparison["candidate"]
        expected_primary = 0.5 * float(
            primary["mean_per_track_gene_exon_coverage_pearson_log1p"]
        ) + 0.5 * float(primary["mean_per_track_pearson_128bp_log1p"])
        validation_ok = (
            validation.get("schema_version") == 2
            and validation.get("model") == "B_no_lora"
            and validation.get("training_loss") == "paper"
            and validation.get("fold") == 1
            and validation.get("locked_test_block_signal_reads") is False
            and validation.get("chromosome_x_train_valid_blocks_read") is True
            and set(validation.get("evaluated_chromosomes", []))
            == {"I", "II", "III", "IV", "V", "X"}
            and float(validation.get("validation_core_coverage_fraction", 0)) >= 0.999
            and int(validation.get("finite_per_track_pearson_128bp", 0)) == 241
            and math.isclose(
                float(candidate_metrics["primary_biological_score"]),
                expected_primary,
                rel_tol=1e-12,
                abs_tol=1e-12,
            )
        )
    except Exception:
        validation_ok = False
    adaptation = run.get("adaptation", {})
    no_lora_contract = (
        adaptation.get("base_organism_index") == 2
        and adaptation.get("c_elegans_organism_embedding") is True
        and adaptation.get("lora_enabled") is False
        and adaptation.get("lora_rank") is None
        and adaptation.get("lora_alpha") is None
        and adaptation.get("lora_target_modules") == []
        and adaptation.get("trunk_policy") == "worm_embeddings_only"
    )
    training = spec.get("training", {})
    execution_ok = (
        not execution_path_errors
        and execution.get("phase") == "P8"
        and execution.get("status") == "completed"
        and execution.get("model") == "B_no_lora"
        and execution.get("loss") == "paper"
        and execution.get("fold") == 1
        and execution.get("seed") == 20260714
        and execution.get("final_test_access") == "prohibited"
        and execution.get("locked_test_block_signal_reads") == 0
        and int(execution.get("physical_gpu", -1)) in {2, 3}
        and execution.get("spec_sha256") == spec_sha
        and execution.get("checkpoint_sha256")
        == sha256(paths_from_execution["checkpoint_path"])
        and execution.get("validation_sha256")
        == sha256(paths_from_execution["validation_path"])
        and execution.get("comparison_sha256") == sha256(paths["comparison"])
        and execution.get("candidate_adaptation") == adaptation
        and run.get("model") == "B_no_lora"
        and run.get("loss") == "paper"
        and run.get("fold") == 1
        and run.get("seed") == 20260714
        and run.get("steps") == training.get("max_steps") == 2000
        and run.get("n_tracks") == 241
        and run.get("intervals_sha256") == spec["input_sha256"]["fold_1_train"]
        and run.get("means_sha256") == spec["input_sha256"]["means"]
        and no_lora_contract
        and paths_from_execution["log_path"].stat().st_size > 0
    )
    baseline = comparison.get("baseline", {})
    candidate = comparison.get("candidate", {})
    delta = comparison.get("delta", {})
    baseline_ok = (
        baseline.get("model") == "B"
        and baseline.get("loss") == "paper"
        and baseline.get("fold") == 1
        and baseline.get("seed") == 20260714
        and baseline.get("checkpoint_sha256")
        == spec["input_sha256"]["baseline_checkpoint"]
        and candidate.get("adaptation") == adaptation
        and comparison.get("spec_sha256") == spec_sha
        and math.isclose(
            float(comparison.get("primary_biological_score_delta", "nan")),
            float(candidate.get("primary_biological_score", "nan"))
            - float(baseline.get("primary_biological_score", "nan")),
            rel_tol=1e-12,
            abs_tol=1e-12,
        )
        and "primary_biological_score" in delta
    )
    lock = repository_file(spec["locked_inputs"]["final_test_lock"])
    report = repository_file(spec["locked_inputs"]["final_test_report"])
    lock_data = json.loads(lock.read_text())
    checks.extend(
        [
            check("R8.02_locked_inputs", not input_errors, f"errors={input_errors}"),
            check(
                "R8.03_embedding_only_training_contract",
                execution_ok,
                f"gpu={execution.get('physical_gpu')} no_lora={no_lora_contract} path_errors={execution_path_errors}",
            ),
            check("R8.04_development_validation", validation_ok, f"coverage={validation.get('validation_core_coverage_fraction')}"),
            check("R8.05_matched_baseline_comparison", baseline_ok, f"delta={comparison.get('primary_biological_score_delta')}"),
            check(
                "R8.06_final_test_preserved",
                lock_data.get("test_consumed") is True
                and lock_data.get("test_status") == "completed"
                and sha256(lock) == spec["input_sha256"]["final_test_lock"]
                and sha256(report) == spec["input_sha256"]["final_test_report"],
                "P6C lock and report hashes unchanged",
            ),
            check(
                "R8.07_controller_scope",
                state.get("current_phase") == "P8"
                and state.get("approvals", {}).get("G4_gpu_experiments", {}).get("scope")
                == "p8_b_no_lora_fold1",
                f"phase={state.get('current_phase')}",
            ),
        ]
    )
    return {
        "schema_version": 1,
        "phase": "P8",
        "review": "R8-b-no-lora-fold1",
        "reviewed_at": utc_now(),
        "status": "PASS"
        if all(item["status"] == "PASS" for item in checks)
        else "FAIL",
        "checks": checks,
    }


def _p9_ordered_digest(rows: list[tuple[str, ...]]) -> str:
    digest = hashlib.sha256()
    for row in rows:
        digest.update(("\t".join(row) + "\n").encode("utf-8"))
    return digest.hexdigest()


def _review_p9_dataset_contract(spec: dict[str, Any]) -> tuple[dict[str, Any], list[str]]:
    """Re-hash the complete P9 training dataset without reading signal values."""

    errors: list[str] = []
    try:
        track_rows = read_tsv(repository_file(spec["locked_inputs"]["track_manifest"]))
        group_rows = read_tsv(repository_file(spec["locked_inputs"]["group_manifest"]))
        means_rows = read_tsv(repository_file(spec["locked_inputs"]["means"]))
        expected_count = int(spec.get("dataset_contract", {}).get("track_count", -1))
        track_ids = [row.get("group_id", "") for row in track_rows]
        group_ids = [row.get("group_id", "") for row in group_rows]
        mean_ids = [row.get("group_id", "") for row in means_rows]
        if not (
            expected_count
            == int(spec.get("dataset_contract", {}).get("group_count", -1))
            == 241
            and len(track_rows) == len(group_rows) == len(means_rows) == 241
            and len(set(track_ids)) == 241
            and track_ids == group_ids == mean_ids
            and all(row.get("include_formal_v2") == "True" for row in group_rows)
        ):
            errors.append("membership_or_order")
        registered: list[tuple[str, ...]] = []
        observed: list[tuple[str, ...]] = []
        total_bytes = 0
        for row in track_rows:
            group_id = row.get("group_id", "")
            relative = row.get("output_path", "")
            expected_sha = row.get("output_sha256", "")
            path = repository_file(relative)
            actual_sha = sha256(path)
            if len(expected_sha) != 64 or actual_sha != expected_sha:
                errors.append(f"{group_id}:hash")
            registered.append((group_id, relative, expected_sha))
            observed.append((group_id, relative, actual_sha))
            total_bytes += path.stat().st_size
        audit = {
            "status": "passed" if not errors else "failed",
            "track_count": len(track_rows),
            "group_count": len(group_rows),
            "verified_bigwig_count": len(observed),
            "verified_bigwig_bytes": total_bytes,
            "ordered_group_ids_sha256": hashlib.sha256(
                ("\n".join(track_ids) + "\n").encode("utf-8")
            ).hexdigest(),
            "registered_bigwig_manifest_digest": _p9_ordered_digest(registered),
            "observed_bigwig_manifest_digest": _p9_ordered_digest(observed),
            "bigwig_signal_reads": 0,
        }
    except Exception as error:
        errors.append(f"dataset:{error}")
        audit = {"status": "failed"}
    return audit, errors


def _review_p9_reuse_contract(spec: dict[str, Any]) -> tuple[dict[str, Any], list[str]]:
    """Independently resolve and hash every reused P6B checkpoint/run/validation."""

    errors: list[str] = []
    contract = spec.get("p6b_reuse_contract", {})
    records: list[dict[str, Any]] = []
    digest_rows: list[tuple[str, ...]] = []
    try:
        source = json.loads(repository_file(contract.get("source_execution", "")).read_text())
        selected = [
            row
            for row in source.get("formal_jobs", [])
            if row.get("model") in {"A", "B"} and row.get("loss") == "paper"
        ]
        selected.sort(
            key=lambda row: (str(row["model"]), int(row["fold"]), int(row["seed"]))
        )
        expected = {
            (model, fold, seed)
            for model in ("A", "B")
            for fold in spec["matrix"]["folds"]
            for seed in spec["matrix"]["seeds"]
        }
        actual = {
            (str(row.get("model")), int(row.get("fold", -1)), int(row.get("seed", -1)))
            for row in selected
        }
        if not (
            contract.get("source_execution") == spec["locked_inputs"]["p6b_execution"]
            and int(contract.get("record_count", -1)) == len(selected) == 30
            and int(contract.get("artifact_count", -1)) == 90
            and actual == expected
        ):
            errors.append("membership")
        for row in selected:
            output_dir = str(row["output_dir"])
            paths = {
                "checkpoint": str(row["checkpoint_path"]),
                "run": f"{output_dir}/run.json",
                "validation": str(row["validation_path"]),
            }
            hashes = {
                name: sha256(repository_file(relative))
                for name, relative in paths.items()
            }
            if hashes["checkpoint"] != row.get("checkpoint_sha256"):
                errors.append(f"{row.get('model')}:{row.get('fold')}:{row.get('seed')}:checkpoint")
            record = {
                "model": str(row["model"]),
                "loss": "paper",
                "fold": int(row["fold"]),
                "seed": int(row["seed"]),
                "checkpoint_path": paths["checkpoint"],
                "checkpoint_sha256": hashes["checkpoint"],
                "run_path": paths["run"],
                "run_sha256": hashes["run"],
                "validation_path": paths["validation"],
                "validation_sha256": hashes["validation"],
            }
            records.append(record)
            digest_rows.append(
                tuple(
                    str(record[key])
                    for key in (
                        "model",
                        "loss",
                        "fold",
                        "seed",
                        "checkpoint_path",
                        "checkpoint_sha256",
                        "run_path",
                        "run_sha256",
                        "validation_path",
                        "validation_sha256",
                    )
                )
            )
        digest = _p9_ordered_digest(digest_rows)
        if digest != contract.get("manifest_sha256"):
            errors.append("manifest_sha256")
        audit = {
            "status": "passed" if not errors else "failed",
            "record_count": len(records),
            "artifact_count": 3 * len(records),
            "manifest_sha256": digest,
            "records": records,
        }
    except Exception as error:
        errors.append(f"reuse:{error}")
        audit = {"status": "failed"}
    return audit, errors


def _p9_float_equal(first: object, second: object) -> bool:
    try:
        return math.isclose(
            float(first), float(second), rel_tol=1e-12, abs_tol=1e-12
        )
    except (TypeError, ValueError):
        return False


def _review_p9_statistics(
    records: list[dict[str, Any]],
    spec: dict[str, Any],
    paired_rows: list[dict[str, str]],
    raw_paired_rows: list[dict[str, str]],
    raw_factorial_rows: list[dict[str, str]],
    factorial_rows: list[dict[str, str]],
) -> list[str]:
    """Recompute every P9 contrast and hierarchical interval from run records."""

    from scripts import v2_submission_statistics as statistics

    errors: list[str] = []
    metrics = tuple(spec.get("analysis", {}).get("metrics", []))
    folds = list(spec["analysis"]["expected_folds"])
    seeds = list(spec["analysis"]["expected_seeds"])
    interval = spec["analysis"]["confidence_interval"]
    if len(metrics) != 11 or len(set(metrics)) != 11:
        return ["metric_contract"]
    plan = statistics.make_hierarchical_resample_plan(
        n_blocks=len(folds),
        n_repeats=len(seeds),
        iterations=int(interval["resamples"]),
        seed=int(interval["seed"]),
    )
    by_configuration: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in records:
        by_configuration[str(row.get("configuration"))].append(row)
    baseline = by_configuration.get("full_worm_lora", [])
    try:
        source = json.loads(
            repository_file(spec["baseline"]["source_execution"]).read_text()
        )
        auxiliary = [
            row
            for row in source.get("formal_jobs", [])
            if row.get("model") == "B" and row.get("loss") == "log1p_mse"
        ]
    except Exception as error:
        return [f"auxiliary:{error}"]
    comparisons = {
        configuration: rows
        for configuration, rows in by_configuration.items()
        if configuration != "full_worm_lora"
    }
    comparisons["b_log1p_mse_auxiliary"] = auxiliary

    expected_raw_paired: dict[tuple[str, str, int, int], dict[str, Any]] = {}
    expected_paired: dict[tuple[str, str], dict[str, Any]] = {}
    try:
        for configuration, candidate in sorted(comparisons.items()):
            for metric in metrics:
                differences = statistics.matched_differences(
                    candidate, baseline, value_key=metric
                )
                matrix = statistics.records_to_balanced_matrix(
                    differences,
                    value_key="difference",
                    expected_blocks=folds,
                    expected_repeats=seeds,
                )
                estimate = statistics.hierarchical_block_bootstrap(
                    matrix.values,
                    plan,
                    confidence_level=float(interval["level"]),
                )
                expected_paired[(configuration, metric)] = {
                    "estimate": estimate.estimate,
                    "low": estimate.ci_lower,
                    "high": estimate.ci_upper,
                }
                for row in differences:
                    expected_raw_paired[
                        (configuration, metric, int(row["fold"]), int(row["seed"]))
                    ] = row
    except Exception as error:
        return [f"paired_recompute:{error}"]

    observed_raw_paired = {
        (row["configuration"], row["metric"], int(row["fold"]), int(row["seed"])): row
        for row in raw_paired_rows
    }
    if len(observed_raw_paired) != len(raw_paired_rows) or set(observed_raw_paired) != set(
        expected_raw_paired
    ):
        errors.append("raw_paired_keys")
    else:
        for key, expected in expected_raw_paired.items():
            observed = observed_raw_paired[key]
            if not (
                observed.get("reference") == "B/paper"
                and _p9_float_equal(observed.get("candidate_value"), expected["candidate_value"])
                and _p9_float_equal(observed.get("reference_value"), expected["comparator_value"])
                and _p9_float_equal(observed.get("difference"), expected["difference"])
            ):
                errors.append(f"raw_paired:{key}")
                break
    observed_paired = {
        (row["configuration"], row["metric"]): row for row in paired_rows
    }
    if len(observed_paired) != len(paired_rows) or set(observed_paired) != set(expected_paired):
        errors.append("paired_keys")
    else:
        for key, expected in expected_paired.items():
            observed = observed_paired[key]
            if not (
                observed.get("reference") == "B/paper"
                and int(observed.get("paired_runs", -1)) == 15
                and int(observed.get("folds", -1)) == 5
                and int(observed.get("seeds", -1)) == 3
                and _p9_float_equal(observed.get("mean_paired_difference"), expected["estimate"])
                and _p9_float_equal(observed.get("ci_low"), expected["low"])
                and _p9_float_equal(observed.get("ci_high"), expected["high"])
                and _p9_float_equal(observed.get("ci_level"), interval["level"])
                and observed.get("ci_method") == interval["method"]
            ):
                errors.append(f"paired:{key}")
                break

    keyed = {
        (str(row.get("configuration")), int(row.get("fold", -1)), int(row.get("seed", -1))): row
        for row in records
    }
    cell_names = {
        "a": "frozen_no_worm_no_lora",
        "lora": "lora_only",
        "worm": "no_lora",
        "full": "full_worm_lora",
    }
    expected_raw_factorial: dict[tuple[str, str, int, int], float] = {}
    expected_factorial: dict[tuple[str, str], dict[str, float]] = {}
    try:
        for metric in metrics:
            values_by_effect: dict[str, list[dict[str, Any]]] = defaultdict(list)
            for fold in folds:
                for seed in seeds:
                    values = {
                        name: float(keyed[(configuration, fold, seed)][metric])
                        for name, configuration in cell_names.items()
                    }
                    effects = {
                        "worm_embedding_main": 0.5
                        * ((values["worm"] - values["a"]) + (values["full"] - values["lora"])),
                        "lora_main": 0.5
                        * ((values["lora"] - values["a"]) + (values["full"] - values["worm"])),
                        "worm_lora_interaction": values["full"]
                        - values["worm"]
                        - values["lora"]
                        + values["a"],
                    }
                    for effect, value in effects.items():
                        expected_raw_factorial[(effect, metric, fold, seed)] = value
                        values_by_effect[effect].append(
                            {"fold": fold, "seed": seed, "difference": value}
                        )
            for effect, effect_rows in values_by_effect.items():
                matrix = statistics.records_to_balanced_matrix(
                    effect_rows,
                    value_key="difference",
                    expected_blocks=folds,
                    expected_repeats=seeds,
                )
                estimate = statistics.hierarchical_block_bootstrap(
                    matrix.values,
                    plan,
                    confidence_level=float(interval["level"]),
                )
                expected_factorial[(effect, metric)] = {
                    "estimate": estimate.estimate,
                    "low": estimate.ci_lower,
                    "high": estimate.ci_upper,
                }
    except Exception as error:
        return errors + [f"factorial_recompute:{error}"]

    observed_raw_factorial = {
        (row["effect"], row["metric"], int(row["fold"]), int(row["seed"])): row
        for row in raw_factorial_rows
    }
    if len(observed_raw_factorial) != len(raw_factorial_rows) or set(
        observed_raw_factorial
    ) != set(expected_raw_factorial):
        errors.append("raw_factorial_keys")
    else:
        for key, expected in expected_raw_factorial.items():
            if not _p9_float_equal(observed_raw_factorial[key].get("value"), expected):
                errors.append(f"raw_factorial:{key}")
                break
    observed_factorial = {
        (row["effect"], row["metric"]): row for row in factorial_rows
    }
    if len(observed_factorial) != len(factorial_rows) or set(observed_factorial) != set(
        expected_factorial
    ):
        errors.append("factorial_keys")
    else:
        for key, expected in expected_factorial.items():
            observed = observed_factorial[key]
            if not (
                int(observed.get("paired_runs", -1)) == 15
                and int(observed.get("folds", -1)) == 5
                and int(observed.get("seeds", -1)) == 3
                and _p9_float_equal(observed.get("mean_effect"), expected["estimate"])
                and _p9_float_equal(observed.get("ci_low"), expected["low"])
                and _p9_float_equal(observed.get("ci_high"), expected["high"])
                and _p9_float_equal(observed.get("ci_level"), interval["level"])
                and observed.get("ci_method") == interval["method"]
            ):
                errors.append(f"factorial:{key}")
                break
    return errors


def review_p9() -> dict[str, Any]:
    """Review the development-only submission-evidence matrix."""

    metadata_dir = REPO_ROOT / "alphagenome_custom/metadata/v2"
    paths = {
        "spec": metadata_dir / "p9_submission_evidence_spec.json",
        "execution": metadata_dir / "p9_submission_evidence_execution.json",
        "results": metadata_dir / "p9_submission_evidence_results.tsv",
        "paired": metadata_dir / "p9_submission_evidence_paired_effects.tsv",
        "raw_paired": metadata_dir / "p9_submission_evidence_raw_paired_effects.tsv",
        "raw_factorial": metadata_dir / "p9_submission_evidence_raw_factorial_effects.tsv",
        "factorial": metadata_dir / "p9_submission_evidence_factorial_effects.tsv",
        "per_track": metadata_dir / "p9_submission_evidence_per_track.tsv",
    }
    missing = [
        str(path.relative_to(REPO_ROOT))
        for path in paths.values()
        if not path.is_file()
    ]
    checks = [check("R9.01_required_outputs", not missing, f"missing={missing}")]
    if missing:
        return {
            "schema_version": 1,
            "phase": "P9",
            "review": "R9-submission-evidence-matrix",
            "reviewed_at": utc_now(),
            "status": "FAIL",
            "checks": checks,
        }
    spec = json.loads(paths["spec"].read_text())
    execution = json.loads(paths["execution"].read_text())
    state = json.loads((metadata_dir / "execution_state.json").read_text())
    input_errors = []
    for key, relative in spec.get("locked_inputs", {}).items():
        try:
            if sha256(repository_file(relative)) != spec["input_sha256"][key]:
                input_errors.append(key)
        except Exception:
            input_errors.append(key)
    implementation_errors = []
    for relative, expected_hash in spec.get("implementation_sha256", {}).items():
        try:
            if sha256(repository_file(relative)) != expected_hash:
                implementation_errors.append(relative)
        except Exception:
            implementation_errors.append(relative)
    if len(spec.get("implementation_sha256", {})) != 10:
        implementation_errors.append("implementation_manifest_membership")
    dataset_audit, dataset_errors = _review_p9_dataset_contract(spec)
    reuse_audit, reuse_errors = _review_p9_reuse_contract(spec)
    dataset_preflight_ok = (
        execution.get("dataset_preflight") == dataset_audit
        and dataset_audit.get("status") == "passed"
        and dataset_audit.get("verified_bigwig_count") == 241
        and dataset_audit.get("bigwig_signal_reads") == 0
        and dataset_audit.get("registered_bigwig_manifest_digest")
        == dataset_audit.get("observed_bigwig_manifest_digest")
    )
    reuse_preflight_ok = (
        execution.get("p6b_reuse_preflight") == reuse_audit
        and reuse_audit.get("status") == "passed"
        and reuse_audit.get("record_count") == 30
        and reuse_audit.get("artifact_count") == 90
        and reuse_audit.get("manifest_sha256")
        == spec.get("p6b_reuse_contract", {}).get("manifest_sha256")
    )
    implementation_ok = (
        not implementation_errors
        and execution.get("implementation_sha256")
        == spec.get("implementation_sha256")
    )
    records = execution.get("registered_records", [])
    configurations = Counter(str(row.get("configuration")) for row in records)
    expected_configurations = {
        str(row["id"]) for row in spec["matrix"]["configurations"]
    }
    reused_sources = Counter(str(row.get("reused_source")) for row in records if row.get("base_job"))
    registered_p6b_reuse = {
        (
            str(row.get("model")),
            int(row.get("fold", -1)),
            int(row.get("seed", -1)),
            str(row.get("checkpoint_path")),
            str(row.get("checkpoint_sha256")),
            str(row.get("run_path")),
            str(row.get("run_sha256")),
            str(row.get("validation_path")),
            str(row.get("validation_sha256")),
        )
        for row in records
        if row.get("reused_source") == "P6B"
    }
    audited_p6b_reuse = {
        (
            str(row.get("model")),
            int(row.get("fold", -1)),
            int(row.get("seed", -1)),
            str(row.get("checkpoint_path")),
            str(row.get("checkpoint_sha256")),
            str(row.get("run_path")),
            str(row.get("run_sha256")),
            str(row.get("validation_path")),
            str(row.get("validation_sha256")),
        )
        for row in reuse_audit.get("records", [])
    }
    identity = {
        (str(row.get("configuration")), int(row.get("fold", -1)), int(row.get("seed", -1)))
        for row in records
    }
    matrix_ok = (
        execution.get("phase") == "P9"
        and execution.get("status") == "completed"
        and execution.get("final_test_access") == "prohibited"
        and execution.get("locked_test_block_signal_reads") == 0
        and execution.get("spec_sha256") == sha256(paths["spec"])
        and 0 < len(execution.get("selected_physical_gpus", [])) <= 2
        and set(execution.get("selected_physical_gpus", [])) <= {2, 3}
        and execution.get("registered_records_expected") == len(records) == 90
        and execution.get("new_training_jobs_expected")
        == execution.get("new_training_jobs_completed")
        == 59
        and execution.get("reused_records_expected")
        == execution.get("reused_records")
        == 31
        and len(identity) == 90
        and set(configurations) == expected_configurations
        and set(configurations.values()) == {15}
        and reused_sources == {"P6B": 30, "P8": 1}
        and len(registered_p6b_reuse) == 30
        and registered_p6b_reuse == audited_p6b_reuse
        and execution.get("controller_invocation")
        == "python -m scripts.v2_phase_controller run --phase P9"
        and execution.get("registered_phase_command")
        == "python -m scripts.run_v2_p9_submission_evidence"
        and execution.get("hostname")
        and execution.get("working_directory") == str(REPO_ROOT)
        and execution.get("python_executable")
        and not execution.get("failures")
    )
    path_errors = []
    parameter_errors = []
    validation_errors = []
    config_specs = {row["id"]: row for row in spec["matrix"]["configurations"]}
    for row in records:
        for key in ("checkpoint_path", "run_path", "validation_path"):
            try:
                artifact = repository_file(row.get(key, ""))
                expected_key = key.replace("_path", "_sha256")
                if sha256(artifact) != row.get(expected_key):
                    path_errors.append(f"{row.get('job_id')}:{key}:hash")
            except Exception:
                path_errors.append(f"{row.get('job_id')}:{key}:missing")
        try:
            validation = json.loads(
                repository_file(row.get("validation_path", "")).read_text()
            )
            fold = int(row.get("fold", -1))
            expected_interval_hash = spec["input_sha256"][f"fold_{fold}_valid"]
            if not (
                validation.get("model") == row.get("model")
                and validation.get("training_loss") == "paper"
                and validation.get("fold") == fold
                and validation.get("seed") == int(row.get("seed", -1))
                and validation.get("sequence_length")
                == spec["training"]["sequence_length"]
                and validation.get("mean_column")
                == f"fold_{fold}_train_nonzero_mean"
                and validation.get("intervals_sha256") == expected_interval_hash
                and validation.get("locked_test_block_signal_reads") is False
                and float(validation.get("validation_core_coverage_fraction", -1.0))
                == 1.0
            ):
                validation_errors.append(str(row.get("job_id")))
        except Exception:
            validation_errors.append(str(row.get("job_id")))
        registered_parameters = config_specs[str(row.get("configuration"))].get(
            "expected_trainable_parameters"
        )
        if int(row.get("trainable_parameters", -1)) != registered_parameters:
            parameter_errors.append(str(row.get("job_id")))
        if row.get("locked_test_block_signal_reads") not in {None, 0, False}:
            path_errors.append(f"{row.get('job_id')}:final_test")
    results = read_tsv(paths["results"])
    paired = read_tsv(paths["paired"])
    raw_paired = read_tsv(paths["raw_paired"])
    raw_factorial = read_tsv(paths["raw_factorial"])
    factorial = read_tsv(paths["factorial"])
    per_track = read_tsv(paths["per_track"])
    factorial_counts = Counter(row["effect"] for row in raw_factorial)
    paired_comparisons = (
        expected_configurations - {"full_worm_lora"}
    ) | {"b_log1p_mse_auxiliary"}
    raw_paired_identity = {
        (row["configuration"], row["metric"], int(row["fold"]), int(row["seed"]))
        for row in raw_paired
    }
    statistical_recompute_errors = _review_p9_statistics(
        records,
        spec,
        paired,
        raw_paired,
        raw_factorial,
        factorial,
    )
    statistical_ok = (
        len(results) == 6
        and len(paired) == 6 * 11
        and len(raw_paired) == 6 * 11 * 15
        and len(raw_paired_identity) == len(raw_paired)
        and {row["configuration"] for row in paired} == paired_comparisons
        and {row["configuration"] for row in raw_paired} == paired_comparisons
        and all(int(row["paired_runs"]) == 15 for row in paired)
        and all(int(row["folds"]) == 5 and int(row["seeds"]) == 3 for row in paired)
        and all(
            math.isfinite(float(row[column]))
            for row in raw_paired
            for column in ("candidate_value", "reference_value", "difference")
        )
        and len(raw_factorial) == 3 * 11 * 15
        and len(factorial) == 3 * 11
        and factorial_counts
        == {
            "worm_embedding_main": 11 * 15,
            "lora_main": 11 * 15,
            "worm_lora_interaction": 11 * 15,
        }
        and all(int(row["paired_runs"]) == 15 for row in factorial)
        and all(int(row["folds"]) == 5 and int(row["seeds"]) == 3 for row in factorial)
        and all(math.isfinite(float(row["mean_effect"])) for row in factorial)
        and all(float(row["ci_low"]) <= float(row["ci_high"]) for row in factorial)
        and not statistical_recompute_errors
    )
    per_track_ok = (
        len(per_track) == 90 * 241 * 3
        and {row["configuration"] for row in per_track} == set(configurations)
        and len({row["track_id"] for row in per_track}) == 241
        and {row["metric"] for row in per_track}
        == {
            "pearson_128bp_log1p",
            "gene_exon_coverage_pearson_log1p",
            "primary_biological_score",
        }
    )
    per_track_values = {
        (
            row["configuration"],
            int(row["fold"]),
            int(row["seed"]),
            row["track_id"],
            row["metric"],
        ): float(row["value"])
        for row in per_track
    }
    per_track_formula_ok = len(per_track_values) == len(per_track)
    if per_track_formula_ok:
        for configuration in configurations:
            for fold in spec["matrix"]["folds"]:
                for seed in spec["matrix"]["seeds"]:
                    for track_id in {
                        row["track_id"] for row in per_track if row["configuration"] == configuration
                    }:
                        prefix = (configuration, fold, seed, track_id)
                        observed = per_track_values.get((*prefix, "primary_biological_score"))
                        expected = 0.5 * (
                            per_track_values.get((*prefix, "pearson_128bp_log1p"), math.nan)
                            + per_track_values.get(
                                (*prefix, "gene_exon_coverage_pearson_log1p"),
                                math.nan,
                            )
                        )
                        if observed is None or not math.isclose(
                            observed, expected, rel_tol=1e-12, abs_tol=1e-12
                        ):
                            per_track_formula_ok = False
                            break
                    if not per_track_formula_ok:
                        break
                if not per_track_formula_ok:
                    break
            if not per_track_formula_ok:
                break
    per_track_ok = per_track_ok and per_track_formula_ok
    lock = repository_file(spec["locked_inputs"]["final_test_lock"])
    report = repository_file(spec["locked_inputs"]["final_test_report"])
    lock_data = json.loads(lock.read_text())
    output_hashes_ok = all(
        execution.get(key) == sha256(paths[path_key])
        for key, path_key in (
            ("results_sha256", "results"),
            ("paired_effects_sha256", "paired"),
            ("raw_paired_effects_sha256", "raw_paired"),
            ("raw_factorial_effects_sha256", "raw_factorial"),
            ("factorial_effects_sha256", "factorial"),
            ("per_track_sha256", "per_track"),
        )
    )
    checks.extend(
        [
            check(
                "R9.02_locked_inputs",
                not input_errors
                and not dataset_errors
                and not reuse_errors
                and dataset_preflight_ok
                and reuse_preflight_ok
                and implementation_ok,
                "input_errors="
                f"{input_errors} dataset_errors={dataset_errors[:3]} "
                f"reuse_errors={reuse_errors[:3]} implementation_errors={implementation_errors}",
            ),
            check("R9.03_matrix_contract", matrix_ok, f"configs={dict(configurations)} reused={dict(reused_sources)}"),
            check(
                "R9.04_artifact_hashes",
                not path_errors and not validation_errors and output_hashes_ok,
                f"path_errors={path_errors[:5]} validation_errors={validation_errors[:5]}",
            ),
            check("R9.05_parameter_contract", not parameter_errors, f"errors={parameter_errors[:5]}"),
            check(
                "R9.06_paired_factorial_statistics",
                statistical_ok,
                f"raw_paired={len(raw_paired)} raw_factorial={len(raw_factorial)} "
                f"summary={len(factorial)} recompute_errors={statistical_recompute_errors[:3]}",
            ),
            check("R9.07_per_track_source_data", per_track_ok, f"rows={len(per_track)} tracks={len({row['track_id'] for row in per_track})}"),
            check(
                "R9.08_final_test_preserved",
                lock_data.get("test_consumed") is True
                and lock_data.get("test_status") == "completed"
                and sha256(lock) == spec["input_sha256"]["final_test_lock"]
                and sha256(report) == spec["input_sha256"]["final_test_report"],
                "P6C lock and report hashes unchanged",
            ),
            check(
                "R9.09_controller_scope",
                state.get("current_phase") == "P9"
                and state.get("status") == "REVIEWING"
                and state.get("approvals", {})
                .get("G4_gpu_experiments", {})
                .get("approved")
                is True
                and state.get("approvals", {}).get("G4_gpu_experiments", {}).get("scope")
                == "p9_submission_evidence_matrix",
                f"phase={state.get('current_phase')}",
            ),
        ]
    )
    return {
        "schema_version": 1,
        "phase": "P9",
        "review": "R9-submission-evidence-matrix",
        "reviewed_at": utc_now(),
        "status": "PASS" if all(item["status"] == "PASS" for item in checks) else "FAIL",
        "checks": checks,
    }


def review_p10() -> dict[str, Any]:
    """Review the development-only DPY-27 internal application."""

    metadata_dir = REPO_ROOT / "alphagenome_custom/metadata/v2"
    spec_path = metadata_dir / "p10_dpy27_internal_application_spec.json"

    def failed_required_outputs(missing: list[str]) -> dict[str, Any]:
        return {
            "schema_version": 1,
            "phase": "P10",
            "review": "R10-dpy27-internal-application",
            "reviewed_at": utc_now(),
            "status": "FAIL",
            "checks": [
                check("R10.01_required_outputs", False, f"missing={missing}")
            ],
        }

    if not spec_path.is_file():
        return failed_required_outputs(
            [str(spec_path.relative_to(REPO_ROOT))]
        )
    spec = json.loads(spec_path.read_text())
    outputs = spec.get("output_paths", {})
    output_keys = (
        "analysis_audit",
        "preflight_audit",
        "gene_contrasts",
        "run_endpoints",
        "fold_endpoints",
        "fold_gene_contrasts",
        "summary_endpoints",
        "figure_source_data",
        "figure_pdf",
        "figure_svg",
        "figure_tiff",
        "figure_png",
        "figure_qa",
    )
    paths: dict[str, Path] = {"spec": spec_path}
    invalid_paths = []
    for key in output_keys:
        try:
            relative = Path(str(outputs[key]))
            if relative.is_absolute() or ".." in relative.parts:
                raise ValueError(relative)
            paths[key] = REPO_ROOT / relative
        except Exception:
            invalid_paths.append(key)
    missing = invalid_paths + [
        str(path.relative_to(REPO_ROOT))
        for key, path in paths.items()
        if key != "spec" and not path.is_file()
    ]
    if missing:
        return failed_required_outputs(missing)
    checks = [check("R10.01_required_outputs", True, "all frozen outputs present")]
    spec_sha = sha256(spec_path)
    state = json.loads((metadata_dir / "execution_state.json").read_text())
    preflight = json.loads(paths["preflight_audit"].read_text())
    execution = json.loads(paths["analysis_audit"].read_text())
    figure_qa = json.loads(paths["figure_qa"].read_text())

    input_errors = []
    for key, relative in spec.get("locked_inputs", {}).items():
        try:
            if sha256(repository_file(relative)) != spec["input_sha256"][key]:
                input_errors.append(f"{key}:hash")
        except Exception:
            input_errors.append(f"{key}:missing")
    expected_pairs = {
        (seed, fold)
        for seed in spec.get("expected_seeds", [])
        for fold in spec.get("expected_folds", [])
    }
    checkpoint_pairs = set()
    checkpoint_errors = []
    for row in spec.get("checkpoint_matrix", []):
        try:
            pair = (int(row["seed"]), int(row["fold"]))
            checkpoint_pairs.add(pair)
            if sha256(repository_file(row["checkpoint_path"])) != row["checkpoint_sha256"]:
                checkpoint_errors.append(f"{pair}:checkpoint")
            if sha256(repository_file(row["run_path"])) != row["run_sha256"]:
                checkpoint_errors.append(f"{pair}:run")
        except Exception:
            checkpoint_errors.append(f"{row.get('seed')}:{row.get('fold')}:missing")
    frozen_implementation = spec.get("implementation_sha256", {})
    current_implementation = _p10_implementation_sha256()
    implementation_drift = {
        relative: {
            "frozen": frozen_implementation.get(relative),
            "current": current_implementation[relative],
        }
        for relative in P10_IMPLEMENTATION_PATHS
        if not isinstance(frozen_implementation, dict)
        or frozen_implementation.get(relative) != current_implementation[relative]
    }
    static_contract_ok = (
        spec.get("schema_version") == 1
        and spec.get("phase") == "P10"
        and spec.get("analysis_contract")
        == "dpy27_internal_genomic_block_validation_v1"
        and spec.get("analysis_label") == "internal_genomic_block_validation"
        and spec.get("final_test_access") == "prohibited"
        and spec.get("expected_folds") == [1, 2, 3, 4, 5]
        and spec.get("expected_seeds") == [20260714, 20260715, 20260716]
        and spec.get("expected_subwindows_per_fold") == 83
        and spec.get("gpu_policy", {}).get("allowed_physical_indices") == [2, 3]
        and spec.get("model_contract", {}).get("model") == "B"
        and spec.get("model_contract", {}).get("loss") == "paper"
        and spec.get("statistics", {}).get("inference_unit")
        == "fold_after_within_fold_gene_condition_seed_mean"
        and spec.get("figure_source_data_contract") == P10_SOURCE_DATA_CONTRACT
        and set(spec.get("implementation_sha256", {}))
        == set(P10_IMPLEMENTATION_PATHS)
        and _p10_implementation_manifest_valid(frozen_implementation)
        and len(spec.get("checkpoint_matrix", [])) == 15
        and checkpoint_pairs == expected_pairs
        and not checkpoint_errors
    )

    interval_audit = preflight.get("dataset", {}).get("intervals", {})
    preflight_ok = (
        preflight.get("status") == "passed"
        and preflight.get("analysis_contract") == spec.get("analysis_contract")
        and preflight.get("controller_enforced") is True
        and preflight.get("spec_sha256") == spec_sha
        and preflight.get("locked_input_sha256") == spec.get("input_sha256")
        and preflight.get("checkpoint_count") == 15
        and preflight.get("dataset", {}).get("track_count") == 241
        and preflight.get("final_test_access") == "prohibited"
        and preflight.get("locked_test_block_signal_reads") == 0
        and preflight.get("model_or_bigwig_reads") == 0
        and preflight.get("implementation_sha256") == frozen_implementation
        and all(
            int(interval_audit.get(str(fold), interval_audit.get(fold, {})).get("subwindows", -1))
            == 83
            for fold in spec["expected_folds"]
        )
    )

    worker_results = execution.get("worker_results", [])
    worker_pairs = {
        (int(row.get("seed", -1)), int(row.get("fold", -1)))
        for row in worker_results
    }
    physical_gpus = execution.get("physical_gpus", [])
    execution_ok = (
        execution.get("phase") == "P10"
        and execution.get("analysis_contract") == spec.get("analysis_contract")
        and execution.get("status") == "completed"
        and execution.get("spec_sha256") == spec_sha
        and execution.get("preflight_sha256") == sha256(paths["preflight_audit"])
        and 0 < len(physical_gpus) <= 2
        and set(physical_gpus) <= {2, 3}
        and len(worker_results) == 15
        and worker_pairs == expected_pairs
        and all(int(row.get("physical_gpu", -1)) in {2, 3} for row in worker_results)
        and execution.get("run_count") == 15
        and execution.get("fold_count") == 5
        and execution.get("final_test_access") == "prohibited"
        and execution.get("locked_test_block_signal_reads") == 0
        and execution.get("preflight_model_or_bigwig_reads") == 0
        and execution.get("implementation_sha256") == frozen_implementation
        and execution.get("git_commit") == preflight.get("git_commit")
        and _p10_execution_provenance_ok(execution)
        and not execution.get("failures")
    )

    gene_fields = (
        "analysis_label",
        "fold",
        "seed",
        "gene_id",
        "gene_name",
        "chromosome",
        "strand",
        "block_id",
        "exon_bases_evaluated",
        "exon_bases_total",
        "exon_coverage_fraction",
        "complete_exon_coverage",
        "predicted_dpy27_rnai_mean",
        "predicted_vector_rnai_mean",
        "observed_dpy27_rnai_mean",
        "observed_vector_rnai_mean",
        "predicted_log1p_contrast_g0005_minus_g0006",
        "observed_log1p_contrast_g0005_minus_g0006",
    )
    endpoints = tuple(spec["statistics"]["endpoints"])
    shard_errors = []
    shard_gene_rows: list[dict[str, str]] = []
    for row in sorted(
        spec["checkpoint_matrix"], key=lambda item: (int(item["seed"]), int(item["fold"]))
    ):
        seed, fold = int(row["seed"]), int(row["fold"])
        shard_root = (
            REPO_ROOT
            / outputs["run_root"]
            / "shards"
            / f"seed_{seed}"
            / f"fold_{fold}"
        )
        gene_path = shard_root / "gene_contrasts.tsv"
        endpoint_path = shard_root / "run_endpoint.json"
        audit_path = shard_root / "shard_audit.json"
        try:
            genes = read_tsv(gene_path)
            endpoint = json.loads(endpoint_path.read_text())
            audit = json.loads(audit_path.read_text())
            if not genes or tuple(genes[0]) != gene_fields:
                raise ValueError("gene schema")
            if not (
                audit.get("status") == "completed"
                and audit.get("analysis_contract") == spec["analysis_contract"]
                and audit.get("spec_sha256") == spec_sha
                and audit.get("seed") == seed
                and audit.get("fold") == fold
                and audit.get("physical_gpu") in {2, 3}
                and audit.get("checkpoint_sha256") == row["checkpoint_sha256"]
                and audit.get("validation_sha256")
                == spec["input_sha256"][f"fold_{fold}_valid"]
                and audit.get("subwindows_evaluated") == 83
                and audit.get("gene_contrasts_sha256") == sha256(gene_path)
                and audit.get("run_endpoint_sha256") == sha256(endpoint_path)
                and audit.get("final_test_access") == "prohibited"
                and audit.get("locked_test_block_signal_reads") == 0
                and audit.get("implementation_sha256") == frozen_implementation
                and endpoint.get("seed") == seed
                and endpoint.get("fold") == fold
                and endpoint.get("checkpoint_sha256") == row["checkpoint_sha256"]
                and all(math.isfinite(float(endpoint[key])) for key in endpoints)
            ):
                raise ValueError("shard contract")
            if len({gene["gene_id"] for gene in genes}) != len(genes):
                raise ValueError("duplicate gene")
            shard_gene_rows.extend(genes)
        except Exception as error:
            shard_errors.append(f"seed={seed},fold={fold}:{error}")

    aggregate_genes = read_tsv(paths["gene_contrasts"])
    run_rows = read_tsv(paths["run_endpoints"])
    fold_rows = read_tsv(paths["fold_endpoints"])
    fold_gene_rows = read_tsv(paths["fold_gene_contrasts"])
    summary_rows = read_tsv(paths["summary_endpoints"])
    source_rows = read_tsv(paths["figure_source_data"])
    fold_gene_index = {
        (int(row["fold"]), row["gene_id"]): row for row in fold_gene_rows
    }
    fold_recomputation_ok = True
    grouped_aggregate: dict[tuple[int, str], list[dict[str, str]]] = defaultdict(list)
    for row in aggregate_genes:
        grouped_aggregate[(int(row["fold"]), row["gene_id"])].append(row)
    if set(grouped_aggregate) != set(fold_gene_index):
        fold_recomputation_ok = False
    else:
        for key, rows in grouped_aggregate.items():
            if len(rows) != 3 or {int(row["seed"]) for row in rows} != set(
                spec["expected_seeds"]
            ):
                fold_recomputation_ok = False
                break
            expected_row = fold_gene_index[key]
            predicted_case = sum(
                float(row["predicted_dpy27_rnai_mean"]) for row in rows
            ) / 3.0
            predicted_control = sum(
                float(row["predicted_vector_rnai_mean"]) for row in rows
            ) / 3.0
            expected_contrast = math.log1p(max(predicted_case, 0.0)) - math.log1p(
                max(predicted_control, 0.0)
            )
            if not (
                math.isclose(
                    float(expected_row["predicted_dpy27_rnai_mean"]),
                    predicted_case,
                    rel_tol=1e-12,
                    abs_tol=1e-12,
                )
                and math.isclose(
                    float(expected_row["predicted_vector_rnai_mean"]),
                    predicted_control,
                    rel_tol=1e-12,
                    abs_tol=1e-12,
                )
                and math.isclose(
                    float(
                        expected_row[
                            "predicted_log1p_contrast_g0005_minus_g0006"
                        ]
                    ),
                    expected_contrast,
                    rel_tol=1e-12,
                    abs_tol=1e-12,
                )
            ):
                fold_recomputation_ok = False
                break
    fold_endpoint_recomputation_ok = fold_recomputation_ok
    if fold_endpoint_recomputation_ok:
        try:
            from scripts import v2_biological_validation as biological_validation

            fold_endpoint_index = {int(row["fold"]): row for row in fold_rows}
            for fold in spec["expected_folds"]:
                typed_rows = []
                for row in fold_gene_rows:
                    if int(row["fold"]) != fold:
                        continue
                    typed = dict(row)
                    typed["complete_exon_coverage"] = (
                        row["complete_exon_coverage"] == "True"
                    )
                    typed_rows.append(typed)
                recomputed = biological_validation.dpy27_run_endpoints(typed_rows)
                recorded = fold_endpoint_index[fold]
                if any(
                    not math.isclose(
                        float(recorded[endpoint]),
                        float(recomputed[endpoint]),
                        rel_tol=1e-12,
                        abs_tol=1e-12,
                    )
                    for endpoint in endpoints
                ):
                    fold_endpoint_recomputation_ok = False
                    break
        except Exception:
            fold_endpoint_recomputation_ok = False
    summary_recomputation_ok = _p10_summary_recomputation_ok(
        aggregate_genes,
        fold_rows,
        summary_rows,
        spec["statistics"],
    )
    source_data_ok = _p10_source_data_ok(
        source_rows,
        fold_gene_rows,
        fold_rows,
        summary_rows,
        spec["figure_source_data_contract"],
    )
    aggregate_ok = (
        aggregate_genes == shard_gene_rows
        and len(aggregate_genes) == int(execution.get("gene_rows", -1))
        and fold_gene_rows
        and tuple(fold_gene_rows[0]) == gene_fields
        and {int(row["fold"]) for row in fold_gene_rows}
        == set(spec["expected_folds"])
        and {row["seed"] for row in fold_gene_rows} == {"seed_mean"}
        and len({(row["fold"], row["gene_id"]) for row in fold_gene_rows})
        == len(fold_gene_rows)
        and fold_recomputation_ok
        and fold_endpoint_recomputation_ok
        and summary_recomputation_ok
        and len(run_rows) == 15
        and {
            (int(row["seed"]), int(row["fold"])) for row in run_rows
        }
        == expected_pairs
        and all(math.isfinite(float(row[key])) for row in run_rows for key in endpoints)
        and len(fold_rows) == 5
        and {int(row["fold"]) for row in fold_rows} == set(spec["expected_folds"])
        and all(int(row["n_seeds"]) == 3 for row in fold_rows)
        and all(math.isfinite(float(row[key])) for row in fold_rows for key in endpoints)
        and len(summary_rows) == len(endpoints) == 4
        and {row["endpoint"] for row in summary_rows} == set(endpoints)
        and all(
            int(row["n_folds"]) == 5
            and int(row["n_runs"]) == 15
            and row["inference_unit"]
            == "fold_after_within_fold_gene_condition_seed_mean"
            and int(row["bootstrap_replicates"])
            == int(spec["statistics"]["bootstrap_replicates"])
            and int(row["bootstrap_seed"])
            == int(spec["statistics"]["bootstrap_seed"])
            and math.isfinite(float(row["estimate"]))
            and float(row["ci_95_low"]) <= float(row["ci_95_high"])
            for row in summary_rows
        )
        and source_data_ok
    )

    base_artifacts = {
        "preflight_audit",
        "gene_contrasts",
        "run_endpoints",
        "fold_endpoints",
        "fold_gene_contrasts",
        "summary_endpoints",
        "figure_source_data",
        "figure_pdf",
        "figure_svg",
        "figure_tiff",
        "figure_png",
        "figure_qa",
    }
    artifact_errors = []
    for key in base_artifacts:
        record = execution.get("artifacts", {}).get(key, {})
        if not (
            record.get("path") == outputs[key]
            and record.get("sha256") == sha256(paths[key])
            and int(record.get("size_bytes", -1)) == paths[key].stat().st_size
        ):
            artifact_errors.append(key)
    qa_artifact_errors = []
    for key in ("pdf", "svg", "tiff", "png"):
        output_key = f"figure_{key}"
        record = figure_qa.get("artifacts", {}).get(key, {})
        if not (
            record.get("path") == outputs[output_key]
            and record.get("sha256") == sha256(paths[output_key])
            and int(record.get("size_bytes", -1)) == paths[output_key].stat().st_size
        ):
            qa_artifact_errors.append(output_key)
    figure_ok = (
        figure_qa.get("status") == "passed"
        and figure_qa.get("spec_sha256") == spec_sha
        and figure_qa.get("source_data_path") == outputs["figure_source_data"]
        and figure_qa.get("source_data_sha256") == sha256(paths["figure_source_data"])
        and figure_qa.get("static_preflight", {}).get("summary", {}).get("ready")
        is True
        and figure_qa.get("pdf_text_audit", {}).get("auditable") is True
        and figure_qa.get("pdf_text_audit", {}).get("below_minimum_count") == 0
        and figure_qa.get("png_pixel_audit", {}).get("passed") is True
        and figure_qa.get("tiff_pixel_audit", {}).get("passed") is True
        and not qa_artifact_errors
    )

    lock = repository_file(spec["locked_inputs"]["final_test_lock"])
    final_report = repository_file(spec["locked_inputs"]["final_test_report"])
    lock_data = json.loads(lock.read_text())
    controller_ok = (
        state.get("current_phase") == "P10"
        and state.get("status") == "REVIEWING"
        and state.get("approvals", {}).get("G4_gpu_experiments", {}).get("approved")
        is True
        and state.get("approvals", {}).get("G4_gpu_experiments", {}).get("scope")
        == "p10_dpy27_internal_application"
    )
    checks.extend(
        [
            check(
                "R10.02_frozen_contract",
                static_contract_ok and not input_errors,
                f"input_errors={input_errors} checkpoint_errors={checkpoint_errors} implementation_drift={sorted(implementation_drift)}",
            ),
            check(
                "R10.03_preflight",
                preflight_ok,
                f"checkpoints={preflight.get('checkpoint_count')} tracks={preflight.get('dataset', {}).get('track_count')}",
            ),
            check(
                "R10.04_execution_matrix",
                execution_ok,
                f"status={execution.get('status')} workers={len(worker_results)} gpus={physical_gpus}",
            ),
            check(
                "R10.05_restart_shards",
                not shard_errors,
                f"shards={15 - len(shard_errors)}/15 errors={shard_errors[:3]}",
            ),
            check(
                "R10.06_fold_primary_source_data",
                aggregate_ok and not artifact_errors,
                f"runs={len(run_rows)} folds={len(fold_rows)} summaries={len(summary_rows)} artifact_errors={artifact_errors}",
            ),
            check(
                "R10.07_figure_qa",
                figure_ok,
                f"qa={figure_qa.get('status')} artifact_errors={qa_artifact_errors}",
            ),
            check(
                "R10.08_final_test_preserved",
                lock_data.get("test_consumed") is True
                and lock_data.get("test_status") == "completed"
                and sha256(lock) == spec["input_sha256"]["final_test_lock"]
                and sha256(final_report)
                == spec["input_sha256"]["final_test_report"],
                "P6C lock and report hashes unchanged",
            ),
            check(
                "R10.09_controller_scope",
                controller_ok,
                f"phase={state.get('current_phase')} status={state.get('status')}",
            ),
        ]
    )
    return {
        "schema_version": 1,
        "phase": "P10",
        "review": "R10-dpy27-internal-application",
        "reviewed_at": utc_now(),
        "status": "PASS" if all(item["status"] == "PASS" for item in checks) else "FAIL",
        "frozen_implementation_sha256": frozen_implementation,
        "current_implementation_sha256": current_implementation,
        "implementation_drift": implementation_drift,
        "checks": checks,
    }


def review_p11() -> dict[str, Any]:
    """Review the CPU-only replicate-holdout data contract."""

    root = REPO_ROOT / "alphagenome_custom/metadata/v2/replicate_holdout_v1"
    required = {
        "spec": root / "p11_replicate_holdout_spec.json",
        "execution": root / "p11_replicate_holdout_execution.json",
        "assignments": root / "holdout_assignments.tsv",
        "roles": root / "member_roles.tsv",
        "training_manifest": root / "training_track_manifest.tsv",
        "heldout_manifest": root / "heldout_track_manifest.tsv",
        "training_groups": root / "training_group_manifest.tsv",
        "heldout_groups": root / "heldout_group_manifest.tsv",
        "means": root / "track_nonzero_means.tsv",
        "means_summary": root / "track_nonzero_means_summary.json",
    }
    missing = [str(path.relative_to(REPO_ROOT)) for path in required.values() if not path.is_file()]
    if missing:
        return {
            "schema_version": 1,
            "phase": "P11",
            "review": "R11-replicate-holdout-data",
            "reviewed_at": utc_now(),
            "status": "FAIL",
            "checks": [check("R11.01_required_outputs", False, f"missing={missing}")],
        }
    assignments = list(csv.DictReader(required["assignments"].open(newline=""), delimiter="\t"))
    roles = list(csv.DictReader(required["roles"].open(newline=""), delimiter="\t"))
    training = list(csv.DictReader(required["training_manifest"].open(newline=""), delimiter="\t"))
    heldout = list(csv.DictReader(required["heldout_manifest"].open(newline=""), delimiter="\t"))
    execution = json.loads(required["execution"].read_text())
    means_summary = json.loads(required["means_summary"].read_text())
    counts = {role: sum(row.get("assignment_role") == role for row in assignments) for role in (
        "primary_candidate", "primary_pending_qc", "supplementary_candidate", "training_only"
    )}
    checks = [
        check("R11.01_required_outputs", True, "all P11 metadata and means outputs present"),
        check("R11.02_assignment_counts", counts == {"primary_candidate": 57, "primary_pending_qc": 1, "supplementary_candidate": 22, "training_only": 161}, str(counts)),
        check("R11.03_member_roles", len(roles) == 482 and len({row["run_accession"] for row in roles}) == 482, f"rows={len(roles)}"),
        check("R11.04_manifest_order", len(training) == 241 and [row["group_id"] for row in training] == sorted(row["group_id"] for row in training), f"training_rows={len(training)}"),
        check("R11.05_heldout_count", len(heldout) == 80, f"heldout_rows={len(heldout)}"),
        check("R11.06_means", means_summary.get("tracks") == 241 and means_summary.get("locked_test_block_signal_reads") == 0, str(means_summary)),
        check("R11.07_locked_test", execution.get("locked_test_access") == "prohibited" and execution.get("locked_test_signal_reads") == 0, str(execution.get("locked_test_signal_reads"))),
    ]
    status = "PASS" if all(item["status"] == "PASS" for item in checks) else "FAIL"
    return {
        "schema_version": 1,
        "phase": "P11",
        "review": "R11-replicate-holdout-data",
        "reviewed_at": utc_now(),
        "status": status,
        "checks": checks,
    }


def review_p12() -> dict[str, Any]:
    """Review P12 task completion and the two-layer validation guards."""

    root = REPO_ROOT / "alphagenome_custom/metadata/v2/replicate_holdout_v1"
    retry_spec = REPO_ROOT / "alphagenome_custom/metadata/v2/p12_replicate_holdout_retry_spec.json"
    base_spec = REPO_ROOT / "alphagenome_custom/metadata/v2/p12_replicate_holdout_spec.json"
    spec_path = retry_spec if retry_spec.is_file() else base_spec
    required = {
        "spec": spec_path,
        "execution": root / "p12_replicate_holdout_execution.json",
        "results": root / "p12_replicate_holdout_results.tsv",
        "records": root / "p12_replicate_holdout_job_records.tsv",
        "metrics": root / "p12_replicate_holdout_metrics.tsv",
        "role_metrics": root / "p12_heldout_role_metrics.tsv",
        "paired": root / "p12_replicate_holdout_paired_effects.tsv",
        "raw_paired": root / "p12_replicate_holdout_raw_paired_effects.tsv",
    }
    missing = [str(path.relative_to(REPO_ROOT)) for path in required.values() if not path.is_file()]
    if missing:
        return {
            "schema_version": 1,
            "phase": "P12",
            "review": "R12-replicate-holdout-matrix",
            "reviewed_at": utc_now(),
            "status": "FAIL",
            "checks": [check("R12.01_required_outputs", False, f"missing={missing}")],
        }
    spec = json.loads(required["spec"].read_text())
    if spec.get("base_spec"):
        base_path = REPO_ROOT / spec["base_spec"]
        spec_payload = json.loads(base_path.read_text())
        spec_payload.update({
            "retry_id": spec.get("retry_id"),
            "supersedes_spec": spec.get("supersedes_spec"),
            "correction": spec.get("correction"),
            "base_spec_sha256": spec.get("base_spec_sha256"),
        })
        spec = spec_payload
    execution = json.loads(required["execution"].read_text())
    records = list(csv.DictReader(required["records"].open(newline=""), delimiter="\t"))
    results = list(csv.DictReader(required["results"].open(newline=""), delimiter="\t"))
    metrics = list(csv.DictReader(required["metrics"].open(newline=""), delimiter="\t"))
    role_metrics = list(csv.DictReader(required["role_metrics"].open(newline=""), delimiter="\t"))
    checks = [
        check("R12.01_required_outputs", True, "all P12 outputs present"),
        check("R12.02_static_contract", spec.get("phase") == "P12" and spec.get("final_test_access") == "prohibited" and spec.get("matrix", {}).get("total_tasks") == 165, "phase/P12 task contract"),
        check("R12.03_job_count", execution.get("status") == "completed" and execution.get("registered_tasks_completed") == 165 and len(records) == 165, f"execution={execution.get('registered_tasks_completed')} records={len(records)}"),
        check("R12.04_metric_rows", len(metrics) == 165 * 2 * 11 and {row["scope"] for row in metrics} == {"internal", "heldout"}, f"metric_rows={len(metrics)}"),
        check("R12.05_role_metrics", len(role_metrics) > 0 and {row["role"] for row in role_metrics} == {"primary", "supplementary"}, f"role_metric_rows={len(role_metrics)}"),
        check("R12.06_paired_analysis", required["paired"].stat().st_size > 0 and required["raw_paired"].stat().st_size > 0, "paired effect tables present"),
        check("R12.07_configuration_counts", len(results) == 13 and {int(row["jobs"]) for row in results} >= {5, 15}, f"summary_rows={len(results)}"),
        check("R12.08_locked_test", execution.get("final_test_access") == "prohibited" and execution.get("locked_test_block_signal_reads") == 0, str(execution.get("locked_test_block_signal_reads"))),
        check("R12.09_physical_gpus", set(execution.get("selected_physical_gpus", [])) <= {2, 3} and execution.get("selected_physical_gpus"), str(execution.get("selected_physical_gpus"))),
    ]
    status = "PASS" if all(item["status"] == "PASS" for item in checks) else "FAIL"
    return {
        "schema_version": 1,
        "phase": "P12",
        "review": "R12-replicate-holdout-matrix",
        "reviewed_at": utc_now(),
        "status": status,
        "checks": checks,
    }


def review_p13() -> dict[str, Any]:
    """Review the I--V-only replicate-agreement reference without loading signal."""

    metadata_dir = REPO_ROOT / "alphagenome_custom/metadata/v2"
    contract_dir = REPO_ROOT / "results/v2_p13_iv_interval_contract"
    analysis_dir = REPO_ROOT / "results/v2_p13_replicate_agreement_iv"
    paths = {
        "execution": metadata_dir / "p13_submission_evidence_execution.json",
        "contract": contract_dir / "p13_iv_interval_contract.json",
        "analysis_audit": analysis_dir / "p13_replicate_agreement_audit.json",
        "per_group": analysis_dir / "p13_replicate_agreement_per_group.tsv",
        "by_fold": analysis_dir / "p13_replicate_agreement_by_fold.tsv",
        "summary": analysis_dir / "p13_replicate_agreement_summary.tsv",
    }
    missing = [name for name, path in paths.items() if not path.is_file()]
    checks = [check("R13.01_required_outputs", not missing, f"missing={missing}")]
    if missing:
        return {
            "schema_version": 1,
            "phase": "P13",
            "review": "R13 I--V replicate-agreement reference",
            "reviewed_at": utc_now(),
            "status": "FAIL",
            "checks": checks,
        }

    execution = json.loads(paths["execution"].read_text())
    contract = json.loads(paths["contract"].read_text())
    analysis_audit = json.loads(paths["analysis_audit"].read_text())
    per_group = read_tsv(paths["per_group"])
    by_fold = read_tsv(paths["by_fold"])
    summary = read_tsv(paths["summary"])
    contract_folds = contract.get("folds", [])
    no_x = (
        contract.get("allowed_chromosomes") == ["I", "II", "III", "IV", "V"]
        and len(contract_folds) == 5
        and all(
            set(fold.get("selected_chromosomes", [])) <= {"I", "II", "III", "IV", "V"}
            and "X" in fold.get("excluded_chromosomes", [])
            for fold in contract_folds
        )
    )
    valid_group_rows = (
        len(per_group) == 57 * 5
        and {int(row["fold"]) for row in per_group} == {1, 2, 3, 4, 5}
        and all(row["comparison"] == "training_aggregate_vs_heldout_biological_unit" for row in per_group)
    )
    valid_summary = len(by_fold) == 5 and len(summary) == 3 and {
        row["metric"] for row in summary
    } == {
        "gene_exon_pearson_log1p",
        "pearson_128bp_log1p",
        "primary_biological_score",
    }
    execution_ok = (
        execution.get("phase") == "P13"
        and execution.get("status") == "completed"
        and execution.get("controller_invocation") == "python -m scripts.v2_phase_controller run --phase P13"
        and "no gpu" in execution.get("resource_contract", "").lower()
    )
    audit_ok = (
        analysis_audit.get("folds") == [1, 2, 3, 4, 5]
        and "No GPU was used." in analysis_audit.get("prohibited", [])
        and "No locked final-test intervals or signals were opened." in analysis_audit.get("prohibited", [])
    )
    checks.extend(
        [
            check("R13.02_execution_contract", execution_ok, str(execution.get("resource_contract"))),
            check("R13.03_i_to_v_boundary", no_x, f"folds={len(contract_folds)}"),
            check("R13.04_group_coverage", valid_group_rows, f"per_group_rows={len(per_group)}"),
            check("R13.05_summary_coverage", valid_summary, f"by_fold={len(by_fold)} summary={len(summary)}"),
            check("R13.06_analysis_audit", audit_ok, str(analysis_audit.get("prohibited"))),
        ]
    )
    return {
        "schema_version": 1,
        "phase": "P13",
        "review": "R13 I--V replicate-agreement reference",
        "reviewed_at": utc_now(),
        "status": "PASS" if all(item["status"] == "PASS" for item in checks) else "FAIL",
        "checks": checks,
    }


def review_p14() -> dict[str, Any]:
    """Review P14's frozen I--V Model B reference without signal access."""

    metadata_dir = REPO_ROOT / "alphagenome_custom/metadata/v2"
    result_dir = REPO_ROOT / "results/v2_p14_iv_model_reference"
    paths = {
        "spec": metadata_dir / "p14_iv_frozen_model_reference_spec.json",
        "execution": metadata_dir / "p14_iv_frozen_model_reference_execution.json",
        "records": result_dir / "p14_iv_model_reference_records.tsv",
        "by_fold": result_dir / "p14_iv_model_reference_by_fold.tsv",
        "summary": result_dir / "p14_iv_model_reference_summary.tsv",
        "paired": result_dir / "p14_iv_model_vs_replicate_reference_by_fold.tsv",
    }
    missing = [name for name, path in paths.items() if not path.is_file()]
    checks = [check("R14.01_required_outputs", not missing, f"missing={missing}")]
    if missing:
        return {
            "schema_version": 1, "phase": "P14", "review": "R14 frozen Model B I--V reference",
            "reviewed_at": utc_now(), "status": "FAIL", "checks": checks,
        }
    spec = json.loads(paths["spec"].read_text())
    execution = json.loads(paths["execution"].read_text())
    records = read_tsv(paths["records"])
    by_fold = read_tsv(paths["by_fold"])
    summary = read_tsv(paths["summary"])
    paired = read_tsv(paths["paired"])
    expected_ids = {
        f"p14:iv_model_reference:B_paper:seed{seed}:fold{fold}"
        for fold in range(1, 6) for seed in (20260714, 20260715, 20260716)
    }
    execution_ok = (
        execution.get("phase") == "P14" and execution.get("status") == "completed"
        and execution.get("registered_tasks") == 15 and execution.get("registered_tasks_completed") == 15
        and execution.get("locked_test_block_signal_reads") == 0
        and set(execution.get("selected_physical_gpus", [])) <= {2, 3}
    )
    record_ok = (
        len(records) == 15 and {row["job_id"] for row in records} == expected_ids
        and all(row["checkpoint_sha256"] for row in records)
        and all(float(row["primary_biological_score"]) == float(row["primary_biological_score"]) for row in records)
    )
    summary_ok = (
        len(by_fold) == 5 and {int(row["fold"]) for row in by_fold} == {1, 2, 3, 4, 5}
        and len(summary) == 3 and {row["metric"] for row in summary}
        == {"primary_biological_score", "gene_exon_pearson_log1p", "pearson_128bp_log1p"}
        and len(paired) == 5
    )
    spec_ok = (
        spec.get("phase") == "P14"
        and spec.get("evaluation", {}).get("allowed_chromosomes") == ["I", "II", "III", "IV", "V"]
        and spec.get("evaluation", {}).get("final_test_access") == "prohibited"
        and len(spec.get("jobs", [])) == 15
    )
    checks.extend(
        [
            check("R14.02_preregistered_i_to_v_contract", spec_ok, str(spec.get("evaluation", {}))),
            check("R14.03_execution_contract", execution_ok, f"completed={execution.get('registered_tasks_completed')} gpus={execution.get('selected_physical_gpus')}"),
            check("R14.04_checkpoint_record_coverage", record_ok, f"records={len(records)}"),
            check("R14.05_fold_and_reference_summary", summary_ok, f"folds={len(by_fold)} summary={len(summary)} paired={len(paired)}"),
        ]
    )
    return {
        "schema_version": 1,
        "phase": "P14",
        "review": "R14 frozen Model B I--V reference",
        "reviewed_at": utc_now(),
        "status": "PASS" if all(item["status"] == "PASS" for item in checks) else "FAIL",
        "checks": checks,
    }


def review_p15() -> dict[str, Any]:
    """Review P15's metadata-only I--V training/validation contract."""

    metadata_dir = REPO_ROOT / "alphagenome_custom/metadata/v2"
    result_dir = REPO_ROOT / "results/v2_p15_iv_training_interval_contract"
    execution_path = metadata_dir / "p15_iv_training_contract_execution.json"
    contract_path = result_dir / "p15_iv_training_interval_contract.json"
    interval_paths = [
        result_dir / "intervals" / f"fold_{fold}" / role
        for fold in range(1, 6)
        for role in ("train.tsv", "valid.tsv")
    ]
    missing = [str(path.relative_to(REPO_ROOT)) for path in [execution_path, contract_path, *interval_paths] if not path.is_file()]
    checks = [check("R15.01_required_outputs", not missing, f"missing={missing}")]
    if missing:
        return {
            "schema_version": 1,
            "phase": "P15",
            "review": "R15 I--V-only training interval contract",
            "reviewed_at": utc_now(),
            "status": "FAIL",
            "checks": checks,
        }
    execution = json.loads(execution_path.read_text())
    contract = json.loads(contract_path.read_text())
    all_rows = [row for path in interval_paths for row in read_tsv(path)]
    execution_ok = (
        execution.get("phase") == "P15"
        and execution.get("status") == "completed"
        and execution.get("controller_invocation") == "python -m scripts.v2_phase_controller run --phase P15"
        and "no gpu" in execution.get("resource_contract", "").lower()
    )
    contract_ok = (
        contract.get("phase") == "P15"
        and contract.get("allowed_chromosomes") == ["I", "II", "III", "IV", "V"]
        and len(contract.get("folds", [])) == 5
    )
    interval_ok = (
        bool(all_rows)
        and all(row["chromosome"] in {"I", "II", "III", "IV", "V"} for row in all_rows)
        and all(row["role"] in {"train", "valid"} for row in all_rows)
        and all(row["role"] != "test_locked" for row in all_rows)
    )
    source_x_excluded = all(
        "X" in fold[role].get("source_chromosomes", [])
        and fold[role].get("selected_chromosomes") != ["X"]
        for fold in contract.get("folds", [])
        for role in ("train", "valid")
    )
    checks.extend(
        [
            check("R15.02_execution_contract", execution_ok, str(execution.get("resource_contract"))),
            check("R15.03_contract_structure", contract_ok, f"folds={len(contract.get('folds', []))}"),
            check("R15.04_no_x_or_locked_rows", interval_ok, f"rows={len(all_rows)}"),
            check("R15.05_source_x_exclusion", source_x_excluded, "source manifests contain X; output contract excludes it"),
        ]
    )
    return {
        "schema_version": 1,
        "phase": "P15",
        "review": "R15 I--V-only training interval contract",
        "reviewed_at": utc_now(),
        "status": "PASS" if all(item["status"] == "PASS" for item in checks) else "FAIL",
        "checks": checks,
    }


def review_p16() -> dict[str, Any]:
    """Review P16's strict I--V-only fold-normalization output."""

    metadata_dir = REPO_ROOT / "alphagenome_custom/metadata/v2"
    result_dir = REPO_ROOT / "results/v2_p16_iv_training_normalization"
    paths = {
        "spec": metadata_dir / "p16_iv_training_normalization_spec.json",
        "execution": metadata_dir / "p16_iv_training_normalization_execution.json",
        "means": result_dir / "track_nonzero_means.tsv",
        "summary": result_dir / "track_nonzero_means_summary.json",
    }
    missing = [name for name, path in paths.items() if not path.is_file()]
    checks = [check("R16.01_required_outputs", not missing, f"missing={missing}")]
    if missing:
        return {
            "schema_version": 1,
            "phase": "P16",
            "review": "R16 strict I--V training normalization",
            "reviewed_at": utc_now(),
            "status": "FAIL",
            "checks": checks,
        }
    spec = json.loads(paths["spec"].read_text())
    execution = json.loads(paths["execution"].read_text())
    summary = json.loads(paths["summary"].read_text())
    rows = read_tsv(paths["means"])
    required_columns = {"group_id", *(f"fold_{fold}_train_nonzero_mean" for fold in range(1, 6))}
    execution_ok = (
        execution.get("phase") == "P16"
        and execution.get("status") == "completed"
        and execution.get("controller_invocation") == "python -m scripts.v2_phase_controller run --phase P16"
        and "no gpu" in execution.get("resource_contract", "").lower()
    )
    contract_ok = (
        spec.get("phase") == "P16"
        and spec.get("allowed_chromosomes") == ["I", "II", "III", "IV", "V"]
        and summary.get("contract") == "strict_interval_root"
        and summary.get("allowed_chromosomes") == ["I", "II", "III", "IV", "V"]
        and summary.get("locked_test_block_signal_reads") == 0
    )
    means_ok = (
        len(rows) == 241
        and bool(rows)
        and required_columns.issubset(rows[0])
        and "development_train_nonzero_mean" not in rows[0]
        and all(float(row[column]) > 0 for row in rows for column in required_columns - {"group_id"})
    )
    checks.extend(
        [
            check("R16.02_execution_contract", execution_ok, str(execution.get("resource_contract"))),
            check("R16.03_i_to_v_source_contract", contract_ok, str(summary.get("allowed_chromosomes"))),
            check("R16.04_241_fold_means", means_ok, f"rows={len(rows)}"),
        ]
    )
    return {
        "schema_version": 1,
        "phase": "P16",
        "review": "R16 strict I--V training normalization",
        "reviewed_at": utc_now(),
        "status": "PASS" if all(item["status"] == "PASS" for item in checks) else "FAIL",
        "checks": checks,
    }


def review_p17() -> dict[str, Any]:
    """Review P17's strict I--V controlled ablation and baseline matrix."""

    metadata_dir = REPO_ROOT / "alphagenome_custom/metadata/v2"
    result_dir = REPO_ROOT / "results/v2_p17_iv_controlled_matrix"
    paths = {
        "spec": metadata_dir / "p17_iv_controlled_matrix_spec.json",
        "execution": metadata_dir / "p17_iv_controlled_matrix_execution.json",
        "records": result_dir / "p17_iv_controlled_records.tsv",
        "by_fold": result_dir / "p17_iv_controlled_by_fold.tsv",
        "paired": result_dir / "p17_iv_controlled_paired_effects.tsv",
    }


def review_p18() -> dict[str, Any]:
    """Review P18 external-source transfer and I-V-only reprocessing."""

    metadata_dir = REPO_ROOT / "alphagenome_custom/metadata/v2"
    spec_path = metadata_dir / "p18_external_continuation_spec.json"
    execution_path = metadata_dir / "p18_external_continuation_execution.json"
    checks = [check("R18.01_required_outputs", spec_path.is_file() and execution_path.is_file(), "spec and execution")]
    if not spec_path.is_file() or not execution_path.is_file():
        return {
            "schema_version": 1,
            "phase": "P18",
            "review": "R18 external RNA continuation",
            "reviewed_at": utc_now(),
            "status": "FAIL",
            "checks": checks,
        }
    spec = json.loads(spec_path.read_text())
    execution = json.loads(execution_path.read_text())
    expected_runs = {"SRR18463404", "SRR18463405", "SRR18463406", "SRR3560831"}
    output = spec.get("output", {})
    coverage_dir = REPO_ROOT / output.get("coverage_dir", "")
    summary_path = coverage_dir / "p18_reprocessing_summary.json"
    audit_dir = coverage_dir / "sample_audits"
    summary = json.loads(summary_path.read_text()) if summary_path.is_file() else {}
    audits: list[dict[str, Any]] = []
    if audit_dir.is_dir():
        for run in sorted(expected_runs):
            path = audit_dir / f"{run}.json"
            if path.is_file():
                audits.append(json.loads(path.read_text()))
    static_ok = (
        spec.get("phase") == "P18"
        and spec.get("status") == "preregistered"
        and spec.get("allowed_chromosomes") == ["I", "II", "III", "IV", "V"]
        and set(spec.get("expected_runs", [])) == expected_runs
        and spec.get("represented_head") == "RNA_V2_G0054"
        and spec.get("final_test_access") == "prohibited"
    )
    execution_ok = (
        execution.get("phase") == "P18"
        and execution.get("status") == "completed"
        and set(execution.get("expected_runs", [])) == expected_runs
        and execution.get("locked_test_block_signal_reads") == 0
        and execution.get("model_inference") == "not performed"
    )
    coverage_ok = (
        summary.get("completed_runs") == 4
        and set(summary.get("completed_accessions", [])) == expected_runs
        and len(audits) == 4
        and all(item.get("coverage_chromosomes") == ["I", "II", "III", "IV", "V"] for item in audits)
        and all((REPO_ROOT / item.get("output_path", "")).is_file() for item in audits)
        and all(item.get("output_sha256") == sha256(REPO_ROOT / item["output_path"]) for item in audits)
    )
    checks.extend(
        [
            check("R18.02_static_continuation_contract", static_ok, f"runs={sorted(spec.get('expected_runs', []))}"),
            check("R18.03_execution_boundary", execution_ok, f"inference={execution.get('model_inference')}, locked_reads={execution.get('locked_test_block_signal_reads')}"),
            check("R18.04_complete_i_to_v_source_coverage", coverage_ok, f"audits={len(audits)}"),
        ]
    )
    return {
        "schema_version": 1,
        "phase": "P18",
        "review": "R18 external RNA continuation",
        "reviewed_at": utc_now(),
        "status": "PASS" if all(item["status"] == "PASS" for item in checks) else "FAIL",
        "checks": checks,
    }


def review_p19() -> dict[str, Any]:
    """Review external represented-head inference without reopening the lock."""

    metadata_dir = REPO_ROOT / "alphagenome_custom/metadata/v2"
    spec_path = metadata_dir / "p19_external_represented_head_scoring_spec.json"
    execution_path = metadata_dir / "p19_external_represented_head_scoring_execution.json"
    result_dir = REPO_ROOT / "results/v2_p19_external_scoring"
    checks = [check("R19.01_required_outputs", spec_path.is_file() and execution_path.is_file(), "spec and execution")]
    if not spec_path.is_file() or not execution_path.is_file():
        return {"schema_version": 1, "phase": "P19", "review": "R19 external represented-head scoring", "reviewed_at": utc_now(), "status": "FAIL", "checks": checks}
    spec = json.loads(spec_path.read_text())
    execution = json.loads(execution_path.read_text())
    jobs = spec.get("jobs", [])
    records = execution.get("jobs", [])
    result_path = result_dir / "p19_external_records.tsv"
    source_path = result_dir / "p19_external_by_source.tsv"
    expected_runs = {"SRR18463404", "SRR18463405", "SRR18463406", "SRR3560831"}
    static_ok = (
        spec.get("phase") == "P19"
        and spec.get("status") == "preregistered"
        and spec.get("allowed_chromosomes") == ["I", "II", "III", "IV", "V"]
        and spec.get("final_test_access") == "prohibited"
        and spec.get("represented_head", {}).get("group_id") == "RNA_V2_G0054"
        and len(jobs) == 60
        and {job.get("run_accession") for job in jobs} == expected_runs
        and all(job.get("model") == "B" and job.get("training_loss") == "paper" for job in jobs)
    )
    execution_ok = (
        execution.get("phase") == "P19"
        and execution.get("status") == "completed"
        and execution.get("registered_tasks_completed") == 60
        and len(records) == 60
        and execution.get("locked_test_block_signal_reads") == 0
        and execution.get("final_test_access") == "prohibited"
        and execution.get("model_inference") == "external_represented_head_inference"
        and not execution.get("failures")
    )
    hashes_ok = (
        result_path.is_file()
        and source_path.is_file()
        and execution.get("records_sha256") == sha256(result_path)
        and execution.get("by_source_sha256") == sha256(source_path)
        and all((REPO_ROOT / row["result_path"]).is_file() and row.get("result_sha256") == sha256(REPO_ROOT / row["result_path"]) for row in records)
    )
    metrics_ok = all(
        row.get("evaluated_chromosomes") == "I,II,III,IV,V"
        and row.get("locked_test_block_signal_reads") is False
        and all(math.isfinite(float(row[key])) for key in ("primary_biological_score", "gene_exon_pearson_log1p", "pearson_128bp_log1p"))
        for row in records
    ) if records else False
    source_ok = all(
        (REPO_ROOT / source["coverage_path"]).is_file()
        and source.get("coverage_sha256") == sha256(REPO_ROOT / source["coverage_path"])
        and source.get("represented_group_id") == "RNA_V2_G0054"
        for source in spec.get("source", {}).get("source_records", [])
    ) and len(spec.get("source", {}).get("source_records", [])) == 4
    checks.extend([
        check("R19.02_static_inference_contract", static_ok, f"jobs={len(jobs)}, runs={sorted({job.get('run_accession') for job in jobs})}"),
        check("R19.03_execution_boundary", execution_ok, f"completed={execution.get('registered_tasks_completed')}, locked_reads={execution.get('locked_test_block_signal_reads')}"),
        check("R19.04_output_hashes", hashes_ok, f"records={len(records)}"),
        check("R19.05_metrics_and_i_to_v_only", metrics_ok, f"records={len(records)}"),
        check("R19.06_external_source_hashes", source_ok, f"sources={len(spec.get('source', {}).get('source_records', []))}"),
    ])
    return {"schema_version": 1, "phase": "P19", "review": "R19 external represented-head scoring", "reviewed_at": utc_now(), "status": "PASS" if all(item["status"] == "PASS" for item in checks) else "FAIL", "checks": checks}


def review_p20() -> dict[str, Any]:
    """Review checksum-verified public eQTL source acquisition."""

    metadata_dir = REPO_ROOT / "alphagenome_custom/metadata/v2"
    spec_path = metadata_dir / "p20_eqtl_public_source_download_spec.json"
    output_dir = REPO_ROOT / "shared/external_eqtl/p20_sources_v1"
    manifest_path = output_dir / "p20_source_manifest.json"
    checks = [check("R20.01_required_outputs", spec_path.is_file() and manifest_path.is_file(), "spec and source manifest")]
    if not spec_path.is_file() or not manifest_path.is_file():
        return {"schema_version": 1, "phase": "P20", "review": "R20 public eQTL source acquisition", "reviewed_at": utc_now(), "status": "FAIL", "checks": checks}
    spec = json.loads(spec_path.read_text())
    manifest = json.loads(manifest_path.read_text())
    records = manifest.get("source_records", [])
    required_keys = {"expression_counts", "eqtl_truth_archive", "cendr_isotype_vcf"}
    records_by_key = {row.get("key"): row for row in records}
    files_ok = (
        set(records_by_key) == required_keys
        and all(
            (REPO_ROOT / row.get("path", "")).is_file()
            and int(row.get("bytes", 0)) == (REPO_ROOT / row["path"]).stat().st_size
            and row.get("sha256") == sha256(REPO_ROOT / row["path"])
            for row in records
        )
    )
    static_ok = (
        spec.get("contract") == "p20_eqtl_public_source_download_v1"
        and spec.get("allowed_downstream_chromosomes") == ["I", "II", "III", "IV", "V"]
        and manifest.get("source_status") == "DOWNLOADED_PENDING_COHORT_SOURCE_EQUIVALENCE_AND_COORDINATE_AUDIT"
        and manifest.get("spec_sha256") == sha256(spec_path)
    )
    boundary_ok = (
        "independent eQTL validation" in manifest.get("prohibited_claims", [])
        and "regulatory-allele validation" in manifest.get("prohibited_claims", [])
        and "unseen-condition prediction" in manifest.get("prohibited_claims", [])
    )
    checks.extend([
        check("R20.02_source_inventory", set(records_by_key) == required_keys, f"keys={sorted(records_by_key)}"),
        check("R20.03_source_hashes", files_ok, f"files={len(records)}"),
        check("R20.04_static_boundary", static_ok, manifest.get("source_status")),
        check("R20.05_claim_boundary", boundary_ok, "independent/eQTL/regulatory claims remain blocked"),
    ])
    return {"schema_version": 1, "phase": "P20", "review": "R20 public eQTL source acquisition", "reviewed_at": utc_now(), "status": "PASS" if all(item["status"] == "PASS" for item in checks) else "FAIL", "checks": checks}


def review_p21() -> dict[str, Any]:
    """Review external public-candidate variant scoring on I--V only."""

    metadata_dir = REPO_ROOT / "alphagenome_custom/metadata/v2"
    execution_path = metadata_dir / "p21_eqtl_variant_scoring_execution.json"
    variant_manifest_path = REPO_ROOT / "results/v2_p20_eqtl_variant_set/p20_eqtl_variant_set_manifest.json"
    checks = [check("R21.01_required_outputs", execution_path.is_file() and variant_manifest_path.is_file(), "execution and P20 variant manifest")]
    if not execution_path.is_file() or not variant_manifest_path.is_file():
        return {"schema_version": 1, "phase": "P21", "review": "R21 external eQTL candidate variant scoring", "reviewed_at": utc_now(), "status": "FAIL", "checks": checks}

    execution = json.loads(execution_path.read_text())
    variant_manifest = json.loads(variant_manifest_path.read_text())
    jobs = execution.get("jobs", [])
    expected_keys = {
        (fold, seed)
        for fold in range(1, 6)
        for seed in (20260714, 20260715, 20260716)
    }
    observed_keys = {(int(row.get("fold", 0)), int(row.get("seed", 0))) for row in jobs}
    records_ok = (
        len(jobs) == 15
        and observed_keys == expected_keys
        and all(int(row.get("physical_gpu", 0)) in {2, 3} for row in jobs)
        and all(row.get("status") == "completed" for row in jobs)
        and all(row.get("audit_path") for row in jobs)
    )
    audit_records: list[dict[str, Any]] = []
    audit_ok = True
    for job in jobs:
        audit_path = repository_file(job.get("audit_path"))
        if not audit_path.is_file():
            audit_ok = False
            continue
        audit = json.loads(audit_path.read_text())
        audit_records.append(audit)
        if not (
            audit.get("status") == "completed"
            and int(audit.get("variant_count", 0)) > 0
            and audit.get("final_test_access") == "prohibited"
            and int(audit.get("locked_test_block_signal_reads", -1)) == 0
            and int(audit.get("model_loads", 0)) == 1
        ):
            audit_ok = False
    hash_ok = bool(audit_records) and all(
        job.get("audit_sha256") == sha256(repository_file(job["audit_path"]))
        for job in jobs
        if job.get("audit_path") and repository_file(job["audit_path"]).is_file()
    )
    candidate_count = int(variant_manifest.get("candidate_count", variant_manifest.get("unique_candidate_positions", variant_manifest.get("candidate_rows", 0))))
    chromosome_contract_ok = variant_manifest.get("allowed_downstream_chromosomes") == ["I", "II", "III", "IV", "V"] or "I--V" in str(variant_manifest.get("scope", ""))
    manifest_ok = (
        variant_manifest.get("status") == "PASS"
        and chromosome_contract_ok
        and candidate_count > 0
        and variant_manifest.get("selected_vcf_sha256") == sha256(REPO_ROOT / variant_manifest.get("selected_vcf", ""))
    )
    execution_ok = (
        execution.get("phase") == "P21"
        and execution.get("status") == "completed"
        and execution.get("registered_tasks") == 15
        and execution.get("registered_tasks_completed") == 15
        and execution.get("final_test_access") == "prohibited"
        and execution.get("locked_test_block_signal_reads") == 0
    )
    checks.extend([
        check("R21.02_variant_manifest_pass", manifest_ok, f"candidates={candidate_count}"),
        check("R21.03_execution_contract", execution_ok, f"completed={execution.get('registered_tasks_completed')}"),
        check("R21.04_complete_iv_gpu_records", records_ok, f"jobs={len(jobs)} gpus={sorted({row.get('physical_gpu') for row in jobs})}"),
        check("R21.05_scoring_audits", audit_ok, f"audits={len(audit_records)}"),
        check("R21.06_audit_hashes", hash_ok, f"hashed_audits={len(audit_records)}"),
    ])
    return {"schema_version": 1, "phase": "P21", "review": "R21 external eQTL candidate variant scoring", "reviewed_at": utc_now(), "status": "PASS" if all(item["status"] == "PASS" for item in checks) else "FAIL", "checks": checks}


def review_p22() -> dict[str, Any]:
    """Review public dosage-expression association against P21 effects."""

    output_dir = REPO_ROOT / "results/v2_p22_eqtl_effect_association"
    manifest_path = output_dir / "p22_eqtl_effect_association_manifest.json"
    summary_path = output_dir / "eqtl_variant_effects.tsv"
    model_path = output_dir / "eqtl_variant_effects_by_model.tsv"
    checks = [check("R22.01_required_outputs", all(path.is_file() for path in (manifest_path, summary_path, model_path)), "manifest and effect tables")]
    if not all(path.is_file() for path in (manifest_path, summary_path, model_path)):
        return {"schema_version": 1, "phase": "P22", "review": "R22 public eQTL effect association", "reviewed_at": utc_now(), "status": "FAIL", "checks": checks}
    manifest = json.loads(manifest_path.read_text())
    summary = read_tsv(summary_path)
    by_model = read_tsv(model_path)
    hashes_ok = (
        manifest.get("counts_sha256") == sha256(REPO_ROOT / "shared/external_eqtl/p20_sources_v1/GSE186719_Celegans_208strains_609samples_rawCounts.tsv.gz")
        and manifest.get("vcf_sha256") == sha256(REPO_ROOT / "results/v2_p20_eqtl_variant_set/published_eqtl_candidates_iv.vcf")
        and manifest.get("candidate_annotations_sha256") == sha256(REPO_ROOT / "results/v2_p20_eqtl_variant_set/published_eqtl_candidate_annotations.tsv")
    )
    status_ok = manifest.get("status") == "PASS" and int(manifest.get("analyzable_count", 0)) >= 20
    rows_ok = len(summary) == int(manifest.get("summary_count", -1)) and len(by_model) >= int(manifest.get("analyzable_count", 0)) * 15
    boundary_ok = "causal regulatory mechanism" in manifest.get("prohibited_claims", []) and "unseen-condition prediction" in manifest.get("prohibited_claims", [])
    checks.extend([
        check("R22.02_effect_association_status", status_ok, f"analyzable={manifest.get('analyzable_count')}"),
        check("R22.03_effect_table_counts", rows_ok, f"summary={len(summary)} by_model={len(by_model)}"),
        check("R22.04_input_hashes", hashes_ok, "P20 public-source hashes"),
        check("R22.05_claim_boundary", boundary_ok, "causal/unseen-condition claims remain blocked"),
    ])
    return {"schema_version": 1, "phase": "P22", "review": "R22 public eQTL effect association", "reviewed_at": utc_now(), "status": "PASS" if all(item["status"] == "PASS" for item in checks) else "FAIL", "checks": checks}


def review_p23() -> dict[str, Any]:
    """Review the I-V-only LoRA sensitivity pilot."""

    execution_path = REPO_ROOT / "alphagenome_custom/metadata/v2/p23_lora_sensitivity_execution.json"
    result_path = REPO_ROOT / "results/v2_p23_lora_sensitivity_pilot/p23_lora_sensitivity_pilot.tsv"
    spec_path = REPO_ROOT / "alphagenome_custom/metadata/v2/p23_lora_sensitivity_spec.json"
    checks = [
        check("R23.01_required_outputs", all(path.is_file() for path in (execution_path, result_path, spec_path)), "spec, execution and pilot table")
    ]
    if not all(path.is_file() for path in (execution_path, result_path, spec_path)):
        return {"schema_version": 1, "phase": "P23", "review": "R23 LoRA sensitivity pilot", "reviewed_at": utc_now(), "status": "FAIL", "checks": checks}
    execution = json.loads(execution_path.read_text())
    spec = json.loads(spec_path.read_text())
    rows = read_tsv(result_path)
    expected = {row["id"] for row in spec.get("pilot_configurations", [])}
    observed = {row.get("configuration") for row in rows}
    scope_ok = (
        execution.get("phase") == "P23"
        and execution.get("status") == "completed"
        and execution.get("registered_tasks") == len(expected)
        and execution.get("registered_tasks_completed") == len(expected)
        and execution.get("final_test_access") == "prohibited"
        and execution.get("locked_test_block_signal_reads") == 0
        and set(execution.get("selected_physical_gpus", [])) <= {2, 3}
    )
    rows_ok = (
        len(rows) == len(expected)
        and observed == expected
        and all(row.get("evaluated_chromosomes") == "I,II,III,IV,V" for row in rows)
        and all(row.get("locked_test_block_signal_reads") == "False" for row in rows)
    )
    finite_ok = all(
        all(str(row.get(field, "")).lower() not in {"", "nan", "none"} for field in ("primary_biological_score", "gene_exon_coverage_pearson_log1p", "per_track_pearson_128bp_log1p"))
        for row in rows
    )
    checks.extend([
        check("R23.02_static_scope", spec.get("allowed_chromosomes") == ["I", "II", "III", "IV", "V"] and spec.get("final_test_access") == "prohibited", "I-V-only pilot contract"),
        check("R23.03_execution_contract", scope_ok, f"completed={execution.get('registered_tasks_completed')} gpus={execution.get('selected_physical_gpus')}"),
        check("R23.04_complete_pilot_rows", rows_ok, f"rows={len(rows)} configurations={sorted(observed)}"),
        check("R23.05_finite_metrics", finite_ok, "pilot metrics are finite"),
    ])
    return {"schema_version": 1, "phase": "P23", "review": "R23 LoRA sensitivity pilot", "reviewed_at": utc_now(), "status": "PASS" if all(item["status"] == "PASS" for item in checks) else "FAIL", "checks": checks}


def review_p24() -> dict[str, Any]:
    """Review official Borzoi source feasibility without claiming an adapter."""

    execution_path = REPO_ROOT / "alphagenome_custom/metadata/v2/p24_borzoi_adapter_execution.json"
    spec_path = REPO_ROOT / "alphagenome_custom/metadata/v2/p24_borzoi_adapter_spec.json"
    checks = [check("R24.01_required_outputs", execution_path.is_file() and spec_path.is_file(), "Borzoi spec and audit execution")]
    if not execution_path.is_file() or not spec_path.is_file():
        return {"schema_version": 1, "phase": "P24", "review": "R24 Borzoi feasibility audit", "reviewed_at": utc_now(), "status": "FAIL", "checks": checks}
    execution = json.loads(execution_path.read_text())
    source_checks = execution.get("checks", {})
    source_ok = all(bool(source_checks.get(name, {}).get("pass")) for name in ("official_commit", "license_sha256", "required_source_metadata"))
    boundary_ok = execution.get("adapter_status") == "not_yet_run" and execution.get("final_test_access") == "prohibited" and execution.get("locked_test_block_signal_reads") == 0
    checks.extend([
        check("R24.02_official_source", source_ok, f"vendor_commit={execution.get('vendor_commit')} license={execution.get('license')}"),
        check("R24.03_adapter_boundary", boundary_ok, "audit is not presented as a completed C. elegans adapter"),
        check("R24.04_native_checkpoint_ineligibility", source_checks.get("native_checkpoint_fairness", {}).get("pass") is False, "native human/mouse heads are excluded from the primary comparison"),
    ])
    return {"schema_version": 1, "phase": "P24", "review": "R24 Borzoi feasibility audit", "reviewed_at": utc_now(), "status": "PASS" if all(item["status"] == "PASS" for item in checks) else "FAIL", "checks": checks}


def review_p25() -> dict[str, Any]:
    """Review the bounded Borzoi adapter pilot without upgrading its endpoint."""

    execution_path = REPO_ROOT / "alphagenome_custom/metadata/v2/p25_borzoi_adapter_pilot_execution.json"
    spec_path = REPO_ROOT / "alphagenome_custom/metadata/v2/p25_borzoi_adapter_pilot_spec.json"
    checks = [check("R25.01_required_outputs", execution_path.is_file() and spec_path.is_file(), "Borzoi pilot execution and spec")]
    if not execution_path.is_file() or not spec_path.is_file():
        return {"schema_version": 1, "phase": "P25", "review": "R25 Borzoi adapter pilot", "reviewed_at": utc_now(), "status": "FAIL", "checks": checks}
    execution = json.loads(execution_path.read_text())
    spec = json.loads(spec_path.read_text())
    output_ok = execution.get("status") == "completed_pilot" and execution.get("adapter_status") == "pilot_completed"
    boundary_ok = (
        execution.get("chromosome_x_reads") == 0
        and execution.get("locked_test_block_signal_reads") == 0
        and execution.get("native_checkpoint_head_used") is False
        and execution.get("native_trunk_weights_used") is True
        and execution.get("one_bp_result") == "not_available"
    )
    shape_ok = execution.get("trainable_parameter_count", 0) > 0 and execution.get("validation_intervals", 0) > 0
    contract_ok = (
        spec.get("input", {}).get("allowed_chromosomes") == ["I", "II", "III", "IV", "V"]
        and spec.get("target", {}).get("n_tracks") == 241
        and spec.get("target", {}).get("reported_resolution_bp") == 128
    )
    checks.extend([
        check("R25.02_pilot_output", output_ok, f"status={execution.get('status')} adapter_status={execution.get('adapter_status')}"),
        check("R25.03_task_boundary", boundary_ok, "native heads, X and locked-test reads excluded; 1-bp not claimed"),
        check("R25.04_output_shape_and_validation", shape_ok, f"trainable_parameters={execution.get('trainable_parameter_count')} validation_intervals={execution.get('validation_intervals')}"),
        check("R25.05_static_contract", contract_ok, "I-V, 241 heads, 128-bp reported endpoint"),
    ])
    return {"schema_version": 1, "phase": "P25", "review": "R25 Borzoi adapter pilot", "reviewed_at": utc_now(), "status": "PASS" if all(item["status"] == "PASS" for item in checks) else "FAIL", "checks": checks}
    missing = [name for name, path in paths.items() if not path.is_file()]
    checks = [check("R17.01_required_outputs", not missing, f"missing={missing}")]
    if missing:
        return {
            "schema_version": 1,
            "phase": "P17",
            "review": "R17 strict I--V controlled matrix",
            "reviewed_at": utc_now(),
            "status": "FAIL",
            "checks": checks,
        }
    spec = json.loads(paths["spec"].read_text())
    execution = json.loads(paths["execution"].read_text())
    records = read_tsv(paths["records"])
    by_fold = read_tsv(paths["by_fold"])
    paired = read_tsv(paths["paired"])
    expected_configurations = {"B_iv_dual", "B_128bp_only", "Basenji2_style"}
    expected_keys = {
        (configuration, str(fold), str(seed))
        for configuration in expected_configurations
        for fold in range(1, 6)
        for seed in (20260714, 20260715, 20260716)
    }
    observed_keys = {
        (row.get("configuration", ""), row.get("fold", ""), row.get("seed", ""))
        for row in records
    }
    static_ok = (
        spec.get("phase") == "P17"
        and spec.get("status") == "preregistered"
        and spec.get("allowed_chromosomes") == ["I", "II", "III", "IV", "V"]
        and spec.get("matrix", {}).get("total_jobs") == 45
        and {row.get("id") for row in spec.get("matrix", {}).get("configurations", [])}
        == expected_configurations
    )
    execution_ok = (
        execution.get("phase") == "P17"
        and execution.get("status") == "completed"
        and execution.get("registered_tasks") == 45
        and execution.get("registered_tasks_completed") == 45
        and execution.get("final_test_access") == "prohibited"
        and execution.get("locked_test_block_signal_reads") == 0
        and set(execution.get("selected_physical_gpus", [])) <= {2, 3}
    )
    records_ok = (
        len(records) == 45
        and observed_keys == expected_keys
        and all(row.get("evaluated_chromosomes") == "I,II,III,IV,V" for row in records)
        and all(row.get("locked_test_block_signal_reads") == "False" for row in records)
    )
    summaries_ok = (
        len(by_fold) == 15
        and {(row.get("configuration", ""), row.get("fold", "")) for row in by_fold}
        == {(configuration, str(fold)) for configuration in expected_configurations for fold in range(1, 6)}
        and len(paired) == 10
    )
    checks.extend(
        [
            check("R17.02_static_contract", static_ok, f"configs={sorted(expected_configurations)}"),
            check("R17.03_execution_contract", execution_ok, f"completed={execution.get('registered_tasks_completed')}"),
            check("R17.04_complete_i_to_v_records", records_ok, f"records={len(records)}"),
            check("R17.05_fold_and_paired_summaries", summaries_ok, f"by_fold={len(by_fold)} paired={len(paired)}"),
        ]
    )
    return {
        "schema_version": 1,
        "phase": "P17",
        "review": "R17 strict I--V controlled matrix",
        "reviewed_at": utc_now(),
        "status": "PASS" if all(item["status"] == "PASS" for item in checks) else "FAIL",
        "checks": checks,
    }


def write_report(report: dict[str, Any]) -> None:
    audit_dir = AUDIT_ROOT / report["phase"]
    audit_dir.mkdir(parents=True, exist_ok=True)
    (audit_dir / "review.json").write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n"
    )
    lines = [
        f"# {report.get('review', report['phase'])} Review",
        "",
        f"Status: **{report['status']}**",
        "",
        f"Reviewed at: `{report.get('reviewed_at', utc_now())}`",
        "",
        "| Check | Status | Evidence |",
        "| --- | --- | --- |",
    ]
    for item in report["checks"]:
        evidence = str(item["evidence"]).replace("|", "\\|").replace("\n", " ")
        lines.append(f"| {item['check_id']} | {item['status']} | {evidence} |")
    (audit_dir / "review.md").write_text("\n".join(lines) + "\n")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--phase",
        required=True,
        choices=[
            "P0",
            "P1",
            "P2",
            "P3A",
            "P3B",
            "P4",
            "P5",
            "P6A",
            "P6B",
            "P6C",
            "P7",
            "P8",
            "P9",
            "P10",
            "P11",
            "P12",
            "P13",
            "P14",
            "P15",
            "P16",
            "P17",
            "P18",
            "P19",
            "P20",
            "P21",
            "P22",
            "P23",
            "P24",
            "P25",
        ],
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.phase == "P0":
        report = review_p0()
    elif args.phase == "P1":
        report = review_p1()
    elif args.phase == "P2":
        report = review_p2()
    elif args.phase == "P3A":
        report = review_p3a()
    elif args.phase == "P3B":
        report = review_p3b()
    elif args.phase == "P4":
        report = review_p4()
    elif args.phase == "P5":
        report = review_p5()
    elif args.phase == "P6A":
        report = review_p6a()
    elif args.phase == "P6B":
        report = review_p6b()
    elif args.phase == "P6C":
        report = review_p6c()
    elif args.phase == "P7":
        report = review_p7()
    elif args.phase == "P8":
        report = review_p8()
    elif args.phase == "P9":
        report = review_p9()
    elif args.phase == "P10":
        report = review_p10()
    elif args.phase == "P11":
        report = review_p11()
    elif args.phase == "P12":
        report = review_p12()
    elif args.phase == "P13":
        report = review_p13()
    elif args.phase == "P14":
        report = review_p14()
    elif args.phase == "P15":
        report = review_p15()
    elif args.phase == "P16":
        report = review_p16()
    elif args.phase == "P17":
        report = review_p17()
    elif args.phase == "P18":
        report = review_p18()
    elif args.phase == "P19":
        report = review_p19()
    elif args.phase == "P20":
        report = review_p20()
    elif args.phase == "P21":
        report = review_p21()
    elif args.phase == "P22":
        report = review_p22()
    elif args.phase == "P23":
        report = review_p23()
    elif args.phase == "P24":
        report = review_p24()
    elif args.phase == "P25":
        report = review_p25()
    else:
        raise AssertionError(f"Unhandled phase: {args.phase}")
    assert report is not None
    write_report(report)
    print(json.dumps(report, indent=2, sort_keys=True))
    raise SystemExit(0 if report["status"] == "PASS" else 1)


if __name__ == "__main__":
    main()
