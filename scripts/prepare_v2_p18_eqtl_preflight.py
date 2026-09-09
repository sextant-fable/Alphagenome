#!/usr/bin/env python3
"""Summarize the fail-closed eQTL source and coordinate gate.

This script accepts only the metadata registry and its direct-collision screen.
It makes no provenance attestation and deliberately cannot download expression,
VCF or eQTL data. Its output distinguishes a clear direct-collision screen from
the still-required source-equivalence, cohort-context and coordinate gates.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--registry", type=Path, required=True)
    parser.add_argument("--collision-table", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(REPO_ROOT.resolve()))


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def main() -> None:
    args = parse_args()
    registry = args.registry.resolve()
    collision = args.collision_table.resolve()
    output = args.output_dir.resolve()
    if output.exists():
        raise FileExistsError(f"Output directory is immutable: {output}")
    if not registry.is_file() or not collision.is_file():
        raise FileNotFoundError("eQTL registry or collision table is missing")
    registry_rows = read_tsv(registry)
    collision_rows = read_tsv(collision)
    registry_runs = {row.get("run_accession", "") for row in registry_rows}
    collision_runs = {row.get("run_accession", "") for row in collision_rows}
    if not registry_runs or registry_runs != collision_runs:
        raise RuntimeError("eQTL registry and collision screen do not contain the same non-empty run set")
    direct_collision_rows = [row for row in collision_rows if row.get("collision_fields") or row.get("overall_status") == "REJECT"]
    clear_direct_identity = all(
        row.get("metadata_status") == "CLEAR_MACHINE_SCREEN"
        and row.get("raw_checksum_status") == "CLEAR_MACHINE_SCREEN"
        and not row.get("collision_fields")
        for row in collision_rows
    )
    output.mkdir(parents=True)
    payload = {
        "schema_version": 1,
        "created_at": utc_now(),
        "scope": "metadata-only eQTL preflight; no expression, genotype, eQTL truth, coverage, model, chromosome-X or locked-test data were read",
        "registry": relative(registry),
        "registry_sha256": sha256(registry),
        "collision_table": relative(collision),
        "collision_table_sha256": sha256(collision),
        "candidate_runs": len(registry_runs),
        "direct_identity_and_fastq_md5_collision_rows": len(direct_collision_rows),
        "direct_collision_screen": "CLEAR" if clear_direct_identity else "NOT_CLEAR",
        "gate_status": "BLOCKED_PENDING_PROVENANCE_COHORT_AND_COORDINATE_AUDITS",
        "required_before_download_or_scoring": [
            "Attest source equivalence for the complete wild-strain cohort, not only absent direct identifiers.",
            "Obtain the expression-to-strain and eQTL truth manifest and check its sample/accession overlap.",
            "Obtain the exact CeNDR VCF release used by the cohort and its checksum.",
            "Demonstrate allele and coordinate compatibility from the source reference release to WBcel235 on the intended I-V loci.",
            "Freeze a context/head mapping and LD-aware locus-level statistical estimand before model scoring.",
        ],
        "prohibited_claims": [
            "No independent eQTL or regulatory-allele validation exists.",
            "No expression or genotype payload was downloaded.",
            "No fixed output head has been assigned to the wild-strain cohort.",
        ],
    }
    report = output / "p18_eqtl_preflight.json"
    report.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"status": "blocked_preflight_recorded", "runs": len(registry_runs), "direct_collision_screen": payload["direct_collision_screen"], "report": relative(report)}, indent=2))


if __name__ == "__main__":
    main()
