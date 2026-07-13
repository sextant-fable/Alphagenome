#!/usr/bin/env python3
"""Finalize the 482 reprocessed runs into audited biological v2 tracks."""

from __future__ import annotations

from collections import defaultdict
import csv
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
from typing import Any

import numpy as np
import pyBigWig

from scripts import build_rna_seq_groups_v2 as p2


REPO_ROOT = Path(__file__).resolve().parents[1]
METADATA_DIR = REPO_ROOT / "alphagenome_custom/metadata/v2"
NORMALIZED_DIR = REPO_ROOT / "alphagenome_custom/tracks/rna_seq_v2_normalized"
GROUPED_DIR = REPO_ROOT / "alphagenome_custom/tracks/rna_seq_v2_grouped"
SAMPLE_MANIFEST = METADATA_DIR / "rna_seq_samples_v2.tsv"
CONTEXT_MANIFEST = METADATA_DIR / "rna_seq_sample_context_v2.tsv"
ONTOLOGY_MANIFEST = METADATA_DIR / "ontology_mapping_v2.tsv"
CANDIDATE_GROUPS = METADATA_DIR / "rna_seq_groups_v2.tsv"
CANDIDATE_MEMBERS = METADATA_DIR / "rna_seq_group_members_v2.tsv"
FULL_LEDGER = METADATA_DIR / "p3_full_outputs.tsv"
FINAL_GROUPS = METADATA_DIR / "rna_seq_groups_v2_final.tsv"
FINAL_MEMBERS = METADATA_DIR / "rna_seq_group_members_v2_final.tsv"
GROUP_OUTPUTS = METADATA_DIR / "p3_group_outputs.tsv"
PAIR_QC = METADATA_DIR / "replicate_pair_qc_v2_post_reprocessing.tsv"
GROUP_QC = METADATA_DIR / "replicate_group_qc_v2_post_reprocessing.tsv"
DUPLICATE_REVIEW = METADATA_DIR / "duplicate_source_reuse_review_v2.tsv"
SUMMARY_PATH = METADATA_DIR / "p3_group_summary.json"
AUDIT_SCHEMA_SUMMARY = METADATA_DIR / "p3_audit_schema_summary.json"
EXPECTED_RUNS = 482
EXPECTED_GROUPS = 241
TARGET_TOTAL = 100_000_000.0
CHUNK_SIZE = 1_000_000
DUPLICATE_PAIRS = {
    "SRR13964554": "SRR12483327",
    "SRR13964555": "SRR12483328",
    "SRR13964556": "SRR12483329",
}


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise RuntimeError(f"Refusing to write empty table: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def atomic_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(8 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def audit_schema_v2(row: dict[str, Any]) -> dict[str, Any]:
    if row.get("audit_schema_version") == 2:
        return dict(row)
    if row.get("source_transport_backend") != "ncbi_sra":
        raise RuntimeError(
            f"Unexpected transport during audit migration: {row.get('run_accession')}"
        )
    extracted_sha = row.get("extracted_fastq_sha256")
    if row.get("source_fastq_sha256") != extracted_sha:
        raise RuntimeError(
            f"Extracted FASTQ SHA mismatch during audit migration: {row.get('run_accession')}"
        )
    updated = dict(row)
    updated["audit_schema_version"] = 2
    updated["source_reference_ena_fastq_urls"] = updated.pop("source_fastq_urls")
    updated["source_reference_ena_fastq_md5"] = updated.pop("source_fastq_md5")
    updated["source_reference_ena_fastq_bytes"] = updated.pop("source_fastq_bytes")
    updated.pop("source_fastq_sha256")
    updated["source_archive_integrity"] = "NCBI_SDL_MD5_and_local_SHA256"
    updated["extracted_fastq_integrity"] = "local_SHA256_and_layout_aware_record_count"
    return updated


def upgrade_full_audit_schema() -> list[dict[str, Any]]:
    audit_dir = NORMALIZED_DIR / "sample_audits"
    paths = sorted(audit_dir.glob("*.json"))
    if len(paths) != EXPECTED_RUNS:
        raise RuntimeError(
            f"Expected {EXPECTED_RUNS} sample audits before schema lock, got {len(paths)}"
        )
    rows = []
    for path in paths:
        updated = audit_schema_v2(json.loads(path.read_text()))
        atomic_json(path, updated)
        rows.append(updated)
    write_tsv(FULL_LEDGER, rows)
    summary = {
        "schema_version": 2,
        "sample_audits": len(rows),
        "transport": "ncbi_sra",
        "source_reference_fields": [
            "source_reference_ena_fastq_urls",
            "source_reference_ena_fastq_md5",
            "source_reference_ena_fastq_bytes",
        ],
        "source_archive_fields": [
            "source_archive_url",
            "source_archive_md5",
            "source_archive_sha256",
            "source_archive_bytes",
        ],
        "extracted_fastq_fields": [
            "extracted_fastq_sha256",
            "extracted_fastq_bytes",
            "extracted_fastq_record_count",
        ],
        "ambiguous_legacy_fields_removed": [
            "source_fastq_urls",
            "source_fastq_md5",
            "source_fastq_sha256",
            "source_fastq_bytes",
        ],
        "ledger_sha256": sha256(FULL_LEDGER),
        "completed_at": utc_now(),
    }
    atomic_json(AUDIT_SCHEMA_SUMMARY, summary)
    return rows


def chromosomes() -> dict[str, int]:
    result = {}
    with (REPO_ROOT / "alphagenome_custom/reference/genome.fa.fai").open() as handle:
        for line in handle:
            name, length, *_ = line.rstrip("\n").split("\t")
            result[name] = int(length)
    return result


def validate_bigwig(path: Path, expected_chromosomes: dict[str, int]) -> dict[str, float]:
    with pyBigWig.open(str(path)) as bigwig:
        if not bigwig.isBigWig() or bigwig.chroms() != expected_chromosomes:
            raise RuntimeError(f"Invalid bigWig reference: {path}")
        header = bigwig.header()
    stats = {key: float(header[key]) for key in ("minVal", "maxVal", "sumData")}
    if not all(math.isfinite(value) for value in stats.values()):
        raise RuntimeError(f"Non-finite bigWig statistics: {path}")
    if stats["minVal"] < 0 or stats["maxVal"] < 0:
        raise RuntimeError(f"Negative signal: {path}")
    stats["relative_total_error"] = abs(stats["sumData"] - TARGET_TOTAL) / TARGET_TOTAL
    if stats["relative_total_error"] > 1e-5:
        raise RuntimeError(f"Unexpected total signal {stats['sumData']}: {path}")
    return stats


def build_final_hierarchy() -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    samples = read_tsv(SAMPLE_MANIFEST)
    contexts = read_tsv(CONTEXT_MANIFEST)
    ontology = read_tsv(ONTOLOGY_MANIFEST)
    groups, members, _, _ = p2.build_groups(
        samples,
        contexts,
        ontology,
        exclude_duplicate_current_signals=False,
    )
    if len(members) != EXPECTED_RUNS or len(groups) != EXPECTED_GROUPS:
        raise RuntimeError(
            f"Expected {EXPECTED_RUNS} members/{EXPECTED_GROUPS} groups, "
            f"got {len(members)}/{len(groups)}"
        )

    existing_members: dict[str, set[str]] = defaultdict(set)
    for row in read_tsv(CANDIDATE_MEMBERS):
        existing_members[row["group_id"]].add(row["run_accession"])
    existing_by_set = {frozenset(value): key for key, value in existing_members.items()}
    temporary_members: dict[str, set[str]] = defaultdict(set)
    for row in members:
        temporary_members[row["group_id"]].add(row["run_accession"])

    new_sets = [
        frozenset(value)
        for value in temporary_members.values()
        if frozenset(value) not in existing_by_set
    ]
    if new_sets != [frozenset(DUPLICATE_PAIRS)]:
        raise RuntimeError(f"Unexpected final hierarchy changes: {new_sets}")
    next_id = max(int(value.rsplit("G", 1)[1]) for value in existing_members) + 1
    id_map = {}
    for temporary_id, accessions in temporary_members.items():
        key = frozenset(accessions)
        id_map[temporary_id] = existing_by_set.get(key, f"RNA_V2_G{next_id:04d}")

    final_groups = []
    for row in groups:
        updated = dict(row)
        updated["group_id"] = id_map[row["group_id"]]
        updated["group_status"] = "formal_v2_uniform_reprocessed"
        updated["strand"] = "."
        updated["include_formal_v2"] = "True"
        updated["formal_block_reason"] = ""
        updated["signal_unit"] = "1e6_x_100bp_unstranded_coverage"
        updated["normalization_target_total"] = "100000000"
        final_groups.append(updated)
    final_groups.sort(key=lambda row: row["group_id"])

    ledger_by_run = {row["run_accession"]: row for row in read_tsv(FULL_LEDGER)}
    final_members = []
    for row in members:
        accession = row["run_accession"]
        ledger = ledger_by_run.get(accession)
        if ledger is None:
            raise RuntimeError(f"Missing full reprocessing ledger row: {accession}")
        updated = dict(row)
        updated["group_id"] = id_map[row["group_id"]]
        updated["local_path"] = ledger["output_path"]
        updated["sha256"] = ledger["output_sha256"]
        updated["membership_status"] = "formal_v2_uniform_reprocessed"
        updated["raw_bedgraph_total"] = ledger["raw_bedgraph_total"]
        final_members.append(updated)
    final_members.sort(key=lambda row: (row["group_id"], row["run_accession"]))
    if len({row["run_accession"] for row in final_members}) != EXPECTED_RUNS:
        raise RuntimeError("Final RNA-seq membership is not one-to-one")
    return final_groups, final_members


def biological_unit(row: dict[str, Any]) -> str:
    return row["biological_replicate"] or row["technical_unit_id"]


def write_nonzero_runs(
    output: pyBigWig.pyBigWig,
    chromosome: str,
    offset: int,
    values: np.ndarray,
) -> None:
    if values.size == 0:
        return
    boundaries = np.flatnonzero(values[1:] != values[:-1]) + 1
    starts = np.concatenate((np.array([0]), boundaries))
    ends = np.concatenate((boundaries, np.array([values.size])))
    run_values = values[starts]
    keep = run_values != 0
    if not np.any(keep):
        return
    absolute_starts = (starts[keep] + offset).astype(np.int64).tolist()
    absolute_ends = (ends[keep] + offset).astype(np.int64).tolist()
    output.addEntries(
        [chromosome] * len(absolute_starts),
        absolute_starts,
        ends=absolute_ends,
        values=run_values[keep].astype(float).tolist(),
    )


def aggregate_group(
    group: dict[str, Any],
    members: list[dict[str, Any]],
    expected_chromosomes: dict[str, int],
) -> dict[str, Any]:
    group_id = group["group_id"]
    output_path = GROUPED_DIR / f"{group_id}.bw"
    audit_path = GROUPED_DIR / "audits" / f"{group_id}.json"
    source_signature = [
        {
            "run_accession": row["run_accession"],
            "sha256": row["sha256"],
            "raw_bedgraph_total": row["raw_bedgraph_total"],
            "biological_unit": biological_unit(row),
        }
        for row in members
    ]
    signature = hashlib.sha256(
        json.dumps(source_signature, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    if output_path.exists() and audit_path.exists():
        audit = json.loads(audit_path.read_text())
        if (
            audit.get("source_signature_sha256") == signature
            and audit.get("output_sha256") == sha256(output_path)
        ):
            validate_bigwig(output_path, expected_chromosomes)
            return audit

    GROUPED_DIR.mkdir(parents=True, exist_ok=True)
    audit_path.parent.mkdir(parents=True, exist_ok=True)
    temporary = output_path.with_suffix(".bw.tmp")
    temporary.unlink(missing_ok=True)
    handles = {
        row["run_accession"]: pyBigWig.open(str(REPO_ROOT / row["local_path"]))
        for row in members
    }
    output = None
    try:
        if any(handle.chroms() != expected_chromosomes for handle in handles.values()):
            raise RuntimeError(f"Input chromosome mismatch in {group_id}")
        units: dict[str, list[dict[str, Any]]] = defaultdict(list)
        for row in members:
            units[biological_unit(row)].append(row)
        if len(members) == 1:
            source_path = REPO_ROOT / members[0]["local_path"]
            temporary.symlink_to(Path(os.path.relpath(source_path, temporary.parent)))
        else:
            output = pyBigWig.open(str(temporary), "w")
            output.addHeader(list(expected_chromosomes.items()))
            for chromosome, length in expected_chromosomes.items():
                for start in range(0, length, CHUNK_SIZE):
                    end = min(start + CHUNK_SIZE, length)
                    group_values = np.zeros(end - start, dtype=np.float64)
                    for unit_members in units.values():
                        unit_values = np.zeros(end - start, dtype=np.float64)
                        total_weight = 0.0
                        for member in unit_members:
                            weight = float(member["raw_bedgraph_total"])
                            values = handles[member["run_accession"]].values(
                                chromosome, start, end, numpy=True
                            )
                            values = np.nan_to_num(
                                values, nan=0.0, posinf=0.0, neginf=0.0
                            ).astype(np.float64, copy=False)
                            unit_values += values * weight
                            total_weight += weight
                        if total_weight <= 0:
                            raise RuntimeError(f"Non-positive unit weight in {group_id}")
                        group_values += unit_values / total_weight
                    group_values /= len(units)
                    write_nonzero_runs(
                        output,
                        chromosome,
                        start,
                        group_values.astype(np.float32),
                    )
            output.close()
            output = None
        temporary.replace(output_path)
    finally:
        if output is not None:
            output.close()
        for handle in handles.values():
            handle.close()
        if temporary.exists():
            temporary.unlink()

    stats = validate_bigwig(output_path, expected_chromosomes)
    audit = {
        "group_id": group_id,
        "n_runs": len(members),
        "n_biological_units": len({biological_unit(row) for row in members}),
        "run_accessions": ";".join(row["run_accession"] for row in members),
        "source_signature_sha256": signature,
        "aggregation_policy": (
            "identity_symlink_singleton"
            if len(members) == 1
            else "raw_coverage_weighted_within_biological_unit_then_equal_mean_across_units"
        ),
        "output_is_symlink": output_path.is_symlink(),
        "output_path": str(output_path.relative_to(REPO_ROOT)),
        "output_size_bytes": output_path.stat().st_size,
        "output_sha256": sha256(output_path),
        "output_total_signal": f"{stats['sumData']:.12g}",
        "output_relative_total_error": f"{stats['relative_total_error']:.12g}",
        "completed_at": utc_now(),
    }
    atomic_json(audit_path, audit)
    return audit


def duplicate_source_review(
    final_members: list[dict[str, Any]], expected_chromosomes: dict[str, int]
) -> list[dict[str, Any]]:
    by_run = {row["run_accession"]: row for row in final_members}
    rows = []
    for secondary, canonical in DUPLICATE_PAIRS.items():
        first = by_run[canonical]
        second = by_run[secondary]
        first_path = REPO_ROOT / first["local_path"]
        second_path = REPO_ROOT / second["local_path"]
        first_profile = p2.signal_profile(
            {"run_accession": canonical, "local_path": first["local_path"], "sha256": first["sha256"]},
            {},
        )["bin_means"]
        second_profile = p2.signal_profile(
            {"run_accession": secondary, "local_path": second["local_path"], "sha256": second["sha256"]},
            {},
        )["bin_means"]
        correlation = p2.safe_pearson(np.log1p(first_profile), np.log1p(second_profile))
        rows.append(
            {
                "canonical_run": canonical,
                "secondary_run": secondary,
                "provided_bigwig_was_byte_identical": "True",
                "reprocessed_bigwig_sha256_equal": str(sha256(first_path) == sha256(second_path)),
                "reprocessed_log1p_1kb_pearson": f"{correlation:.12g}",
                "canonical_group_id": first["group_id"],
                "secondary_group_id": second["group_id"],
                "same_source_study": str(
                    next(row for row in read_tsv(SAMPLE_MANIFEST) if row["run_accession"] == canonical)["sra_study_accession"]
                    == next(row for row in read_tsv(SAMPLE_MANIFEST) if row["run_accession"] == secondary)["sra_study_accession"]
                ),
                "decision": "retain_distinct_study_tracks_no_cross_source_averaging",
            }
        )
    return rows


def main() -> None:
    upgrade_full_audit_schema()
    expected_chromosomes = chromosomes()
    groups, members = build_final_hierarchy()
    write_tsv(FINAL_GROUPS, groups)
    write_tsv(FINAL_MEMBERS, members)
    members_by_group: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in members:
        members_by_group[row["group_id"]].append(row)

    audits = []
    for index, group in enumerate(groups, start=1):
        print(f"group\t{index}/{len(groups)}\t{group['group_id']}", flush=True)
        audits.append(
            aggregate_group(
                group, members_by_group[group["group_id"]], expected_chromosomes
            )
        )
        write_tsv(GROUP_OUTPUTS, audits)

    samples = read_tsv(SAMPLE_MANIFEST)
    normalized_by_run = {row["run_accession"]: row for row in members}
    normalized_samples = []
    for sample in samples:
        if sample["assay"] != "RNA-Seq":
            continue
        normalized = dict(sample)
        member = normalized_by_run[sample["run_accession"]]
        normalized["local_path"] = member["local_path"]
        normalized["sha256"] = member["sha256"]
        normalized_samples.append(normalized)
    pair_qc, group_qc = p2.build_replicate_qc(
        normalized_samples, groups, members
    )
    write_tsv(PAIR_QC, pair_qc)
    write_tsv(GROUP_QC, group_qc)
    duplicate_rows = duplicate_source_review(members, expected_chromosomes)
    write_tsv(DUPLICATE_REVIEW, duplicate_rows)

    summary = {
        "schema_version": 1,
        "phase": "P3B",
        "created_at": utc_now(),
        "normalized_runs": len(members),
        "formal_groups": len(groups),
        "multi_run_groups": sum(int(row["n_runs"]) > 1 for row in groups),
        "biological_replicate_groups": sum(
            int(row["n_biological_units"]) > 1 for row in groups
        ),
        "group_outputs": len(audits),
        "replicate_pair_qc_rows": len(pair_qc),
        "replicate_group_qc_rows": len(group_qc),
        "replicate_flagged_pairs": sum(row["review_flag"] == "True" for row in pair_qc),
        "duplicate_source_pairs_resolved": len(duplicate_rows),
        "aggregation_policy": (
            "raw_coverage_weighted_within_biological_unit_then_equal_mean_across_units"
        ),
        "normalization_target_total": int(TARGET_TOTAL),
        "sample_manifest_sha256": sha256(SAMPLE_MANIFEST),
        "full_ledger_sha256": sha256(FULL_LEDGER),
        "final_group_manifest_sha256": sha256(FINAL_GROUPS),
        "final_member_manifest_sha256": sha256(FINAL_MEMBERS),
    }
    atomic_json(SUMMARY_PATH, summary)
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
