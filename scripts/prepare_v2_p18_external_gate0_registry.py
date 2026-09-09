#!/usr/bin/env python3
"""Create an immutable, metadata-only Gate 0 registry for external RNA sources.

The registry is deliberately limited to samples whose public metadata supports
mapping to an existing fixed output head. It does not attest source equivalence,
download reads, inspect coverage or access chromosome X. A separate collision
screen and manual provenance decision are required before data acquisition.
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
DEFAULT_SELECTION = REPO_ROOT / "alphagenome_custom/metadata/v2/p18_external_context_selection.tsv"
DEFAULT_HEADS = REPO_ROOT / "alphagenome_custom/metadata/v2/replicate_holdout_v1/training_group_manifest.tsv"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidate-registry", type=Path, required=True)
    parser.add_argument("--selection", type=Path, default=DEFAULT_SELECTION)
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


def repo_relative(path: Path) -> str:
    return str(path.resolve().relative_to(REPO_ROOT.resolve()))


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def write_tsv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    args = parse_args()
    candidate_registry = args.candidate_registry.resolve()
    selection = args.selection.resolve()
    output = args.output_dir.resolve()
    if output.exists():
        raise FileExistsError(f"Output directory is immutable: {output}")
    if not candidate_registry.is_file() or not selection.is_file() or not DEFAULT_HEADS.is_file():
        raise FileNotFoundError("Candidate registry, context selection, or represented-head manifest is missing")

    candidates = read_tsv(candidate_registry)
    selected = read_tsv(selection)
    heads = read_tsv(DEFAULT_HEADS)
    if not candidates or not selected or len(heads) != 241:
        raise RuntimeError("Expected non-empty metadata inputs and exactly 241 represented heads")
    candidate_by_run = {row.get("run_accession", ""): row for row in candidates}
    head_by_id = {row.get("group_id", ""): row for row in heads}
    required_selection = (
        "run_accession", "geo_accession", "bioproject_accession", "proposed_group_id", "selection_status",
        "assay_match", "anatomy_match", "stage_match", "strain_genotype_match", "condition_match",
        "reference_coordinate_plan", "normalization_plan", "context_evidence_url", "mapping_rationale", "claim_boundary",
    )
    missing = sorted(set(required_selection).difference(selected[0]))
    if missing:
        raise RuntimeError(f"Selection file lacks required columns: {missing}")

    registry_rows: list[dict[str, str]] = []
    mapping_rows: list[dict[str, str]] = []
    seen_runs: set[str] = set()
    for row in selected:
        run = row["run_accession"].strip()
        head = row["proposed_group_id"].strip()
        if row["selection_status"].strip() != "selected":
            raise RuntimeError(f"Only selected rows are allowed in Gate 0 registry: {run}")
        if run in seen_runs or run not in candidate_by_run:
            raise RuntimeError(f"Selected run is duplicate or absent from candidate registry: {run}")
        if head not in head_by_id:
            raise RuntimeError(f"Selected represented head is absent: {head}")
        source = dict(candidate_by_run[run])
        if source.get("geo_accession", "") not in ("", row["geo_accession"]):
            raise RuntimeError(f"GEO mismatch for {run}")
        if source.get("bioproject_accession", "") != row["bioproject_accession"]:
            raise RuntimeError(f"BioProject mismatch for {run}")
        source["geo_accession"] = row["geo_accession"]
        source["source_equivalence_audit_status"] = "pending_manual_equivalence_audit"
        registry_rows.append(source)
        head_row = head_by_id[head]
        mapping_rows.append(
            {
                "run_accession": run,
                "proposed_group_id": head,
                "head_assay": head_row.get("assay", ""),
                "head_anatomy": head_row.get("tissue_or_cell_type", ""),
                "head_stage": head_row.get("development_stage", ""),
                "head_strain": head_row.get("strain", ""),
                "head_genotype": head_row.get("genotype", ""),
                **{field: row[field] for field in required_selection if field not in {"run_accession", "proposed_group_id"}},
                "mapping_status": "CONTEXT_SELECTED_PENDING_COLLISION_AND_MANUAL_GATE0",
            }
        )
        seen_runs.add(run)

    output.mkdir(parents=True)
    registry_path = output / "p18_selected_external_registry.tsv"
    mapping_path = output / "p18_fixed_head_mapping.tsv"
    registry_fields = list(candidates[0])
    mapping_fields = list(mapping_rows[0])
    write_tsv(registry_path, sorted(registry_rows, key=lambda row: row["run_accession"]), registry_fields)
    write_tsv(mapping_path, sorted(mapping_rows, key=lambda row: row["run_accession"]), mapping_fields)
    manifest = {
        "schema_version": 1,
        "created_at": utc_now(),
        "scope": "metadata-only external RNA Gate 0 candidate selection and fixed-head context mapping; no external signal, chromosome-X or locked-test data were read",
        "candidate_registry": repo_relative(candidate_registry),
        "candidate_registry_sha256": sha256(candidate_registry),
        "selection": repo_relative(selection),
        "selection_sha256": sha256(selection),
        "represented_head_manifest": repo_relative(DEFAULT_HEADS),
        "represented_head_manifest_sha256": sha256(DEFAULT_HEADS),
        "selected_registry": repo_relative(registry_path),
        "selected_registry_sha256": sha256(registry_path),
        "fixed_head_mapping": repo_relative(mapping_path),
        "fixed_head_mapping_sha256": sha256(mapping_path),
        "selected_runs": len(registry_rows),
        "selected_studies": len({row["bioproject_accession"] for row in registry_rows}),
        "next_required_actions": [
            "Run the fail-closed direct identity and FASTQ-MD5 collision screen.",
            "Complete a manual source-equivalence attestation after reviewing the collision screen.",
            "Approve a separate G1 scope before any raw external RNA download or reprocessing.",
        ],
        "prohibited_claims": [
            "This registry does not establish external prediction performance.",
            "This registry does not establish unseen-target or unseen-condition prediction.",
            "No raw RNA reads have been downloaded or reprocessed.",
        ],
    }
    manifest_path = output / "p18_external_gate0_selection_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"status": "prepared", "selected_runs": len(registry_rows), "selected_studies": manifest["selected_studies"], "registry": repo_relative(registry_path), "mapping": repo_relative(mapping_path)}, indent=2))


if __name__ == "__main__":
    main()
