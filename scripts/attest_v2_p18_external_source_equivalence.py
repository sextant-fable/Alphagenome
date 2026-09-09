#!/usr/bin/env python3
"""Record a bounded provenance attestation after a clean Gate 0 metadata screen.

This is intentionally narrower than an external benchmark pass. It requires
that each selected source has no direct local identity or FASTQ-MD5 collision,
then records that the public source is not source-equivalent to a member of the
482-run collection under the audited identifiers. It neither downloads data nor
establishes post-reprocessing equivalence, RNA quality or model performance.
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
    parser.add_argument("--selected-registry", type=Path, required=True)
    parser.add_argument("--collision-table", type=Path, required=True)
    parser.add_argument("--fixed-head-mapping", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
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


def repo_relative(path: Path) -> str:
    return str(path.resolve().relative_to(REPO_ROOT.resolve()))


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def main() -> None:
    args = parse_args()
    selected = args.selected_registry.resolve()
    collision = args.collision_table.resolve()
    mapping = args.fixed_head_mapping.resolve()
    output = args.output_dir.resolve()
    if output.exists():
        raise FileExistsError(f"Output directory is immutable: {output}")
    if not all(path.is_file() for path in (selected, collision, mapping)):
        raise FileNotFoundError("Selected registry, collision table, or fixed-head mapping is missing")

    selected_rows = read_tsv(selected)
    collision_rows = read_tsv(collision)
    mapping_rows = read_tsv(mapping)
    if not selected_rows or not collision_rows or not mapping_rows:
        raise RuntimeError("Gate 0 input table is empty")
    selected_runs = {row.get("run_accession", "") for row in selected_rows}
    collision_by_run = {row.get("run_accession", ""): row for row in collision_rows}
    mapping_by_run = {row.get("run_accession", ""): row for row in mapping_rows}
    if set(collision_by_run) != selected_runs or set(mapping_by_run) != selected_runs:
        raise RuntimeError("Gate 0 tables do not refer to exactly the same selected source runs")

    attestation_rows: list[dict[str, str]] = []
    cleared_rows: list[dict[str, str]] = []
    for row in selected_rows:
        run = row["run_accession"]
        screen = collision_by_run[run]
        head = mapping_by_run[run]
        expected_screen = {
            "metadata_status": "CLEAR_MACHINE_SCREEN",
            "raw_checksum_status": "CLEAR_MACHINE_SCREEN",
            "collision_fields": "",
            "collision_local_runs": "{}",
        }
        violations = [field for field, expected in expected_screen.items() if screen.get(field, "") != expected]
        if violations:
            raise RuntimeError(f"Cannot attest {run}; collision screen violates {violations}")
        if head.get("mapping_status") != "CONTEXT_SELECTED_PENDING_COLLISION_AND_MANUAL_GATE0":
            raise RuntimeError(f"Cannot attest {run}; fixed-head mapping is not the expected pending state")
        cleared = dict(row)
        cleared["source_equivalence_audit_status"] = "cleared_no_equivalent_source"
        cleared_rows.append(cleared)
        attestation_rows.append(
            {
                "run_accession": run,
                "sra_study_accession": row.get("sra_study_accession", ""),
                "bioproject_accession": row.get("bioproject_accession", ""),
                "geo_accession": row.get("geo_accession", ""),
                "proposed_group_id": head.get("proposed_group_id", ""),
                "metadata_screen": "no direct run, experiment, BioSample, alias, study or BioProject collision",
                "raw_checksum_screen": "no FASTQ MD5 collision in the local 482-run ledger",
                "source_equivalence_conclusion": "cleared_no_equivalent_source_under_audited_identifiers",
                "remaining_requirement": "reprocess after a separate G1 approval and compare post-processing source hashes before scoring",
                "claim_boundary": "This provenance attestation does not establish RNA quality, external prediction accuracy, study-population generalization, or unseen-condition prediction.",
            }
        )

    output.mkdir(parents=True)
    registry_path = output / "p18_attested_external_registry.tsv"
    attestation_path = output / "p18_source_equivalence_attestation.tsv"
    write_tsv(registry_path, sorted(cleared_rows, key=lambda row: row["run_accession"]), list(selected_rows[0]))
    write_tsv(attestation_path, sorted(attestation_rows, key=lambda row: row["run_accession"]), list(attestation_rows[0]))
    manifest = {
        "schema_version": 1,
        "created_at": utc_now(),
        "scope": "Gate 0 provenance attestation from metadata and public FASTQ MD5 identifiers only; no external signal, chromosome-X or locked-test data were read",
        "selected_registry": repo_relative(selected),
        "selected_registry_sha256": sha256(selected),
        "collision_table": repo_relative(collision),
        "collision_table_sha256": sha256(collision),
        "fixed_head_mapping": repo_relative(mapping),
        "fixed_head_mapping_sha256": sha256(mapping),
        "attested_registry": repo_relative(registry_path),
        "attested_registry_sha256": sha256(registry_path),
        "attestation": repo_relative(attestation_path),
        "attestation_sha256": sha256(attestation_path),
        "attested_runs": len(cleared_rows),
        "required_next_gate": "Run the collision screen again using the attested registry. Its clear result remains subject to a separate G1 approval before download or reprocessing.",
        "prohibited_claims": [
            "This is not an external RNA performance result.",
            "This is not a post-reprocessing checksum comparison.",
            "This is not an unseen-target or unseen-condition validation.",
        ],
    }
    manifest_path = output / "p18_source_equivalence_attestation_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"status": "attested", "runs": len(cleared_rows), "registry": repo_relative(registry_path), "attestation": repo_relative(attestation_path)}, indent=2))


if __name__ == "__main__":
    main()
