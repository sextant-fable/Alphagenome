#!/usr/bin/env python3
"""Run fail-closed evidence reviews for completed v2 workflow phases."""

from __future__ import annotations

import argparse
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


def check(check_id: str, passed: bool, evidence: str) -> dict[str, Any]:
    return {
        "check_id": check_id,
        "status": "PASS" if passed else "FAIL",
        "evidence": evidence,
    }


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
    summary = {
        key: float(header[key]) for key in ("minVal", "maxVal", "sumData", "sumSquared")
    }
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
    for accession in sorted(expected_runs):
        if accession not in output_by_run or accession not in sample_by_run:
            continue
        row = output_by_run[accession]
        sample = sample_by_run[accession]
        fastq_paths = [REPO_ROOT / value for value in row["source_fastq_paths"].split(";")]
        expected_md5 = row["source_fastq_md5"].split(";")
        expected_sha = row["source_fastq_sha256"].split(";")
        pilot_row = pilot_by_run[accession]
        expected_names = [Path(value).name for value in pilot_row["fastq_ftp"].split(";")]
        if (
            row["source_fastq_md5"] != pilot_row["fastq_md5"]
            or int(row["source_fastq_bytes"])
            != sum(int(value) for value in pilot_row["fastq_bytes"].split(";"))
            or [path.name for path in fastq_paths] != expected_names
        ):
            fastq_errors.append(f"{accession}:does not match pilot manifest")
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
                == sha256(sample_manifest_path),
                "pilot and 485-sample manifest SHA-256 match the executed record",
            ),
            check(
                "R3A.04_output_membership",
                set(output_by_run) == expected_runs
                and len(list(output_dir.glob("*.bw"))) == 5,
                f"outputs={sorted(output_by_run)}",
            ),
            check(
                "R3A.05_fastq_integrity",
                not fastq_errors,
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
        "--phase", required=True, choices=["P0", "P1", "P2", "P3A"]
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
    else:
        report = review_p3a()
    assert report is not None
    write_report(report)
    print(json.dumps(report, indent=2, sort_keys=True))
    raise SystemExit(0 if report["status"] == "PASS" else 1)


if __name__ == "__main__":
    main()
