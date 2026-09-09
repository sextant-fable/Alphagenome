#!/usr/bin/env python3
"""Fail-closed metadata collision screen for a proposed P14 source.

This utility never downloads external data, reads BigWigs, starts inference or
accesses chromosome X. It compares a caller-supplied external metadata registry
against the tracked v2 source ledgers and writes a machine-auditable screen for
manual provenance review. A clear screen is deliberately not a Gate 0 pass.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[1]
V2_ROOT = REPO_ROOT / "alphagenome_custom" / "metadata" / "v2"
LOCAL_SAMPLES = V2_ROOT / "rna_seq_samples_v2.tsv"
LOCAL_CONTEXT = V2_ROOT / "rna_seq_sample_context_v2.tsv"
LOCAL_MEMBERS = V2_ROOT / "rna_seq_group_members_v2_final.tsv"

IDENTITY_FIELDS = (
    "run_accession",
    "experiment_accession",
    "biosample_accession",
    "sample_alias",
    "sra_study_accession",
    "bioproject_accession",
)
OPTIONAL_MATCH_FIELDS = ("geo_accession", "fastq_md5", "processed_sha256")
REGISTRY_COLUMNS = IDENTITY_FIELDS + OPTIONAL_MATCH_FIELDS + ("source_equivalence_audit_status",)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--external-registry", type=Path, required=True, help="TSV with the P14 candidate metadata fields.")
    parser.add_argument("--output-dir", type=Path, required=True, help="Explicit new directory for the audit TSV and JSON.")
    return parser.parse_args()


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def normalize(value: str | None) -> str:
    return (value or "").strip().casefold()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def require_columns(rows: list[dict[str, str]], columns: tuple[str, ...], path: Path) -> None:
    if not rows:
        raise ValueError(f"Registry is empty: {path}")
    missing = [column for column in columns if column not in rows[0]]
    if missing:
        raise ValueError(f"Registry lacks required columns: {', '.join(missing)}")


def checksum_values(value: str | None) -> set[str]:
    """Return individual checksum values from ENA's multi-file checksum fields."""

    return {
        normalized
        for item in re.split(r"[;,\s]+", value or "")
        if (normalized := normalize(item))
    }


def add_index(index: dict[str, dict[str, set[str]]], field: str, value: str | None, run_accession: str | None) -> None:
    values = checksum_values(value) if field in {"fastq_md5", "processed_sha256"} else {normalize(value)}
    for normalized in values:
        if normalized:
            index[field][normalized].add((run_accession or "unknown").strip() or "unknown")


def build_local_index() -> dict[str, dict[str, set[str]]]:
    for path in (LOCAL_SAMPLES, LOCAL_CONTEXT, LOCAL_MEMBERS):
        if not path.is_file():
            raise FileNotFoundError(f"Required local v2 ledger is missing: {path}")

    index: dict[str, dict[str, set[str]]] = defaultdict(lambda: defaultdict(set))
    for row in read_tsv(LOCAL_SAMPLES):
        run = row.get("run_accession")
        for field in ("run_accession", "experiment_accession", "biosample_accession", "sra_study_accession", "bioproject_accession", "geo_accession", "fastq_md5"):
            add_index(index, field, row.get(field), run)
    for row in read_tsv(LOCAL_CONTEXT):
        run = row.get("run_accession")
        for field in ("run_accession", "experiment_accession", "biosample_accession", "sample_alias", "sra_study_accession", "bioproject_accession", "geo_accession"):
            add_index(index, field, row.get(field), run)
    for row in read_tsv(LOCAL_MEMBERS):
        add_index(index, "processed_sha256", row.get("sha256"), row.get("run_accession"))
    return index


def collision_details(row: dict[str, str], index: dict[str, dict[str, set[str]]]) -> dict[str, list[str]]:
    matches: dict[str, list[str]] = {}
    for field in IDENTITY_FIELDS + OPTIONAL_MATCH_FIELDS:
        values = checksum_values(row.get(field)) if field in {"fastq_md5", "processed_sha256"} else {normalize(row.get(field))}
        found = set()
        for value in values:
            if value:
                found.update(index[field].get(value, set()))
        if found:
            matches[field] = sorted(found)
    return matches


def row_status(row: dict[str, str], matches: dict[str, list[str]]) -> tuple[str, str, str, str, list[str]]:
    missing_identity = [field for field in IDENTITY_FIELDS if not normalize(row.get(field))]
    identity_matches = {field: values for field, values in matches.items() if field in IDENTITY_FIELDS}
    checksum_matches = {field: values for field, values in matches.items() if field in ("fastq_md5", "processed_sha256")}

    if identity_matches:
        metadata_status = "REJECT_COLLISION"
    elif missing_identity:
        metadata_status = "PENDING_IDENTITY_METADATA"
    else:
        metadata_status = "CLEAR_MACHINE_SCREEN"

    raw_checksum_match = "fastq_md5" in checksum_matches
    processed_checksum_match = "processed_sha256" in checksum_matches
    if raw_checksum_match:
        raw_checksum_status = "REJECT_COLLISION"
    elif checksum_values(row.get("fastq_md5")):
        raw_checksum_status = "CLEAR_MACHINE_SCREEN"
    else:
        raw_checksum_status = "PENDING_RAW_CHECKSUM_METADATA"

    if processed_checksum_match:
        processed_checksum_status = "REJECT_COLLISION"
    elif checksum_values(row.get("processed_sha256")):
        processed_checksum_status = "CLEAR_MACHINE_SCREEN"
    else:
        processed_checksum_status = "PENDING_POSTPROCESSING_CHECKSUM"

    source_equivalence = normalize(row.get("source_equivalence_audit_status"))
    if source_equivalence != "cleared_no_equivalent_source":
        equivalence_status = "PENDING_MANUAL_EQUIVALENCE_AUDIT"
    else:
        equivalence_status = "ATTESTED_PENDING_REVIEW"

    if "REJECT_COLLISION" in (metadata_status, raw_checksum_status, processed_checksum_status):
        overall = "REJECT"
    elif metadata_status == "CLEAR_MACHINE_SCREEN" and raw_checksum_status == "CLEAR_MACHINE_SCREEN" and equivalence_status == "ATTESTED_PENDING_REVIEW":
        overall = "CLEAR_FOR_MANUAL_GATE0_REVIEW"
    else:
        overall = "PENDING"
    return overall, metadata_status, raw_checksum_status, processed_checksum_status, missing_identity


def write_tsv(path: Path, rows: list[dict[str, Any]]) -> None:
    fields = [
        "candidate_row",
        "run_accession",
        "overall_status",
        "metadata_status",
        "raw_checksum_status",
        "processed_checksum_status",
        "source_equivalence_status",
        "missing_identity_fields",
        "collision_fields",
        "collision_local_runs",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    args = parse_args()
    registry = args.external_registry.resolve()
    output = args.output_dir.resolve()
    if not registry.is_file():
        raise FileNotFoundError(f"External registry does not exist: {registry}")
    if output.exists():
        raise FileExistsError(f"Output directory already exists; choose a new immutable audit location: {output}")

    rows = read_tsv(registry)
    require_columns(rows, REGISTRY_COLUMNS, registry)
    index = build_local_index()
    report_rows: list[dict[str, Any]] = []
    counts: dict[str, int] = defaultdict(int)
    for candidate_row, row in enumerate(rows, start=2):
        matches = collision_details(row, index)
        overall, metadata_status, raw_checksum_status, processed_checksum_status, missing_identity = row_status(row, matches)
        equivalence_status = "ATTESTED_PENDING_REVIEW" if normalize(row.get("source_equivalence_audit_status")) == "cleared_no_equivalent_source" else "PENDING_MANUAL_EQUIVALENCE_AUDIT"
        counts[overall] += 1
        report_rows.append(
            {
                "candidate_row": candidate_row,
                "run_accession": row.get("run_accession", "").strip(),
                "overall_status": overall,
                "metadata_status": metadata_status,
                "raw_checksum_status": raw_checksum_status,
                "processed_checksum_status": processed_checksum_status,
                "source_equivalence_status": equivalence_status,
                "missing_identity_fields": ";".join(missing_identity),
                "collision_fields": ";".join(sorted(matches)),
                "collision_local_runs": json.dumps(matches, sort_keys=True),
            }
        )

    output.mkdir(parents=True)
    table_path = output / "p14_collision_audit.tsv"
    json_path = output / "p14_collision_audit.json"
    write_tsv(table_path, report_rows)
    payload = {
        "schema_version": 2,
        "created_at": utc_now(),
        "scope": "metadata-only collision screen; no external data download, signal read, inference or locked-test access",
        "gate0_status": "MANUAL_REVIEW_REQUIRED",
        "machine_screen_status_counts": dict(sorted(counts.items())),
        "clearance_rule": "CLEAR_FOR_MANUAL_GATE0_REVIEW requires no direct identity or raw FASTQ checksum match and a supplied equivalence attestation. A processed-data checksum is unavailable before reprocessing and remains a separate post-processing Gate 1 requirement.",
        "required_registry_columns": list(REGISTRY_COLUMNS),
        "input_sha256": {
            str(registry.relative_to(REPO_ROOT)) if registry.is_relative_to(REPO_ROOT) else str(registry): sha256(registry),
            str(LOCAL_SAMPLES.relative_to(REPO_ROOT)): sha256(LOCAL_SAMPLES),
            str(LOCAL_CONTEXT.relative_to(REPO_ROOT)): sha256(LOCAL_CONTEXT),
            str(LOCAL_MEMBERS.relative_to(REPO_ROOT)): sha256(LOCAL_MEMBERS),
        },
        "table": str(table_path),
    }
    json_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"status": "completed_metadata_screen", "table": str(table_path), "audit": str(json_path), "counts": dict(sorted(counts.items())), "gate0_status": "MANUAL_REVIEW_REQUIRED"}, indent=2))


if __name__ == "__main__":
    main()
