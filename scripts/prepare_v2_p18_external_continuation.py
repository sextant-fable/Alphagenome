#!/usr/bin/env python3
"""Freeze an external-RNA continuation after replacement-source Gate 0 review.

The original P18 transport contract remains immutable: its L4 source failed
post-transfer checksum validation and is excluded here. This preparation step
creates a new four-source contract consisting of the three previously verified
young-adult sources and one separately attested young-adult replacement. It
does not download RNA, create coverage, or score a model.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[1]
METADATA_ROOT = REPO_ROOT / "alphagenome_custom/metadata/v2"
ORIGINAL_TRANSPORT = REPO_ROOT / "results/v2_p18_external_source_transport_20260822T185343Z/p18_external_source_transport.tsv"
ORIGINAL_AUDIT_DIR = REPO_ROOT / "results/v2_p18_external_download"
REPLACEMENT_ATTESTED = REPO_ROOT / "results/v2_p18_replacement_source_attestation_20260823T030000Z/p18_attested_external_registry.tsv"
REPLACEMENT_TRANSPORT = REPO_ROOT / "results/v2_p18_replacement_source_transport_20260823T030000Z/p18_external_source_transport.tsv"
SPEC_PATH = METADATA_ROOT / "p18_external_continuation_spec.json"
RESULT_DIR = REPO_ROOT / "results/v2_p18_external_continuation_contract_20260823T030000Z"
ORIGINAL_RUNS = {"SRR18463404", "SRR18463405", "SRR18463406"}
REPLACEMENT_RUN = "SRR3560831"
EXPECTED_RUNS = ORIGINAL_RUNS | {REPLACEMENT_RUN}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=RESULT_DIR)
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, str]], fields: list[str]) -> None:
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


def artifact(path: Path) -> dict[str, str]:
    if not path.is_file():
        raise FileNotFoundError(path)
    return {"path": str(path.relative_to(REPO_ROOT)), "sha256": sha256(path)}


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def validate_original_rows(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    by_run = {row.get("run_accession", ""): row for row in rows}
    if set(by_run) != {"SRR18463404", "SRR18463405", "SRR18463406", "SRR2005820"}:
        raise RuntimeError("The original P18 transport manifest is not the immutable four-run contract")
    selected = [by_run[run] for run in sorted(ORIGINAL_RUNS)]
    if any(row.get("proposed_group_id") != "RNA_V2_G0054" for row in selected):
        raise RuntimeError("The retained original sources no longer map to RNA_V2_G0054")
    for row in selected:
        audit_path = ORIGINAL_AUDIT_DIR / f"{row['run_accession']}.json"
        if not audit_path.is_file():
            raise FileNotFoundError(f"Missing verified original source audit: {audit_path}")
        audit = json.loads(audit_path.read_text(encoding="utf-8"))
        if audit.get("source_fastq_actual_md5") != row.get("fastq_md5"):
            raise RuntimeError(f"Original download checksum audit does not match transport metadata: {row['run_accession']}")
    return selected


def validate_replacement() -> dict[str, str]:
    rows = read_tsv(REPLACEMENT_ATTESTED)
    if len(rows) != 1 or rows[0].get("run_accession") != REPLACEMENT_RUN:
        raise RuntimeError("Replacement provenance attestation must contain only SRR3560831")
    if rows[0].get("source_equivalence_audit_status") != "cleared_no_equivalent_source":
        raise RuntimeError("Replacement source-equivalence attestation is incomplete")
    transport_rows = read_tsv(REPLACEMENT_TRANSPORT)
    if len(transport_rows) != 1 or transport_rows[0].get("run_accession") != REPLACEMENT_RUN:
        raise RuntimeError("Replacement transport manifest must contain only SRR3560831")
    row = transport_rows[0]
    if row.get("proposed_group_id") != "RNA_V2_G0054":
        raise RuntimeError("Replacement does not map to the registered G0054 head")
    if row.get("source_equivalence_audit_status") != "cleared_no_equivalent_source":
        raise RuntimeError("Replacement transport has no cleared equivalence attestation")
    if not row.get("fastq_md5") or not row.get("fastq_ftp") or not row.get("fastq_bytes"):
        raise RuntimeError("Replacement transport metadata is incomplete")
    return row


def main() -> None:
    args = parse_args()
    output = args.output_dir.resolve()
    if output.exists() or SPEC_PATH.exists():
        raise FileExistsError("P18 continuation contract already exists; immutable preregistration is not overwritten")
    originals = validate_original_rows(read_tsv(ORIGINAL_TRANSPORT))
    replacement = validate_replacement()
    rows = sorted([*originals, replacement], key=lambda row: row["run_accession"])
    if {row["run_accession"] for row in rows} != EXPECTED_RUNS:
        raise RuntimeError("P18 continuation sources are incomplete")
    fields = list(rows[0])
    output.mkdir(parents=True)
    transport = output / "p18_external_continuation_transport.tsv"
    write_tsv(transport, rows, fields)
    exclusion = {
        "excluded_run": "SRR2005820",
        "reason": "The completed local transfer did not match the ENA-published MD5 after a full-byte recheck. Its retained .part file and failed transfer evidence remain excluded from all continuation inputs.",
        "claim_boundary": "Exclusion is a transport-integrity decision, not a biological quality judgment.",
    }
    exclusion_path = output / "p18_excluded_l4_source.json"
    exclusion_path.write_text(json.dumps(exclusion, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    spec = {
        "schema_version": 1,
        "phase": "P18",
        "status": "preregistered",
        "created_at": utc_now(),
        "purpose": "I-V-only uniform reprocessing of a provenance-cleared replacement external-RNA source set for later represented-context evaluation.",
        "allowed_chromosomes": ["I", "II", "III", "IV", "V"],
        "expected_runs": sorted(EXPECTED_RUNS),
        "represented_head": "RNA_V2_G0054",
        "source": {
            "original_transport": artifact(ORIGINAL_TRANSPORT),
            "original_download_audits": {run: artifact(ORIGINAL_AUDIT_DIR / f"{run}.json") for run in sorted(ORIGINAL_RUNS)},
            "replacement_attestation": artifact(REPLACEMENT_ATTESTED),
            "replacement_transport": artifact(REPLACEMENT_TRANSPORT),
            "continuation_transport": artifact(transport),
            "excluded_l4_source": artifact(exclusion_path),
        },
        "output": {
            "download_audit_dir": "results/v2_p18_external_continuation_download",
            "fastq_root": "shared/source_reads/v2/p18_external_continuation_fastqs",
            "work_dir": "shared/source_reads/v2/p18_external_continuation_iv_work",
            "coverage_dir": "alphagenome_custom/tracks/rna_seq_v2_external_p18_continuation",
        },
        "resource_contract": "CPU-only source download and reprocessing; no model training or inference; no chromosome-X or locked-test signal access.",
        "claim_boundary": [
            "The four sources measure one represented broad young-adult target context and do not establish unseen-target or unseen-condition prediction.",
            "The sources span two external studies but are insufficient for study-population inference.",
            "This phase produces source coverage only; model evaluation requires a later frozen scoring contract.",
        ],
        "final_test_access": "prohibited",
    }
    SPEC_PATH.write_text(json.dumps(spec, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"status": "preregistered", "spec": str(SPEC_PATH.relative_to(REPO_ROOT)), "transport": str(transport.relative_to(REPO_ROOT)), "runs": len(rows)}, indent=2))


if __name__ == "__main__":
    main()
