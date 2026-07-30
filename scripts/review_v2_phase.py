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
    else:
        report = review_p8()
    assert report is not None
    write_report(report)
    print(json.dumps(report, indent=2, sort_keys=True))
    raise SystemExit(0 if report["status"] == "PASS" else 1)


if __name__ == "__main__":
    main()
