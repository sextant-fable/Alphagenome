#!/usr/bin/env python3
"""Compile the selected replacement source into auditable P18 Gate 0 inputs.

This compiler turns an immutable public-ENA candidate registry plus a reviewed
represented-head decision into the exact registry and mapping table consumed by
the existing collision and source-equivalence tools. It never opens an RNA
signal file or downloads source reads.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
CANDIDATE_REGISTRY = REPO_ROOT / "results/v2_p18_replacement_candidate_metadata_20260822T215000Z/p14_external_candidate_registry.tsv"
SELECTION_CONTRACT = REPO_ROOT / "alphagenome_custom/metadata/v2/p18_external_replacement_context_selection.tsv"
GROUP_MANIFEST = REPO_ROOT / "alphagenome_custom/metadata/v2/replicate_holdout_v1/training_group_manifest.tsv"
DEFAULT_OUTPUT = REPO_ROOT / "results/v2_p18_replacement_gate0_inputs_20260823T030000Z"
EXPECTED_RUN = "SRR3560831"
EXPECTED_GROUP = "RNA_V2_G0054"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
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


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def main() -> None:
    args = parse_args()
    output = args.output_dir.resolve()
    if output.exists():
        raise FileExistsError(f"P18 replacement Gate 0 inputs are immutable: {output}")
    candidates = read_tsv(CANDIDATE_REGISTRY)
    selections = read_tsv(SELECTION_CONTRACT)
    groups = read_tsv(GROUP_MANIFEST)
    if len(candidates) != 1 or candidates[0].get("run_accession") != EXPECTED_RUN:
        raise RuntimeError("Replacement candidate registry must contain only SRR3560831")
    if len(selections) != 1 or selections[0].get("run_accession") != EXPECTED_RUN:
        raise RuntimeError("Replacement selection contract must contain only SRR3560831")
    selection = selections[0]
    candidate = candidates[0]
    if selection.get("selection_status") != "selected" or selection.get("proposed_group_id") != EXPECTED_GROUP:
        raise RuntimeError("Replacement selection is not approved for RNA_V2_G0054")
    for field in ("geo_accession", "bioproject_accession"):
        if candidate.get(field) and selection.get(field) != candidate.get(field):
            raise RuntimeError(f"Replacement selection does not match candidate {field}")
    group = next((row for row in groups if row.get("group_id") == EXPECTED_GROUP), None)
    if group is None:
        raise RuntimeError("RNA_V2_G0054 is absent from the 241-head training manifest")
    expected_context = {
        "assay": "RNA-Seq",
        "tissue_or_cell_type": "whole animal",
        "development_stage": "Young Adult",
        "strain": "N2",
        "genotype": "wild type",
    }
    if any(group.get(field) != value for field, value in expected_context.items()):
        raise RuntimeError("RNA_V2_G0054 no longer has the registered broad young-adult context")

    selected = dict(candidate)
    selected["geo_accession"] = selection["geo_accession"]
    selected["source_equivalence_audit_status"] = "pending_manual_equivalence_audit"
    mapping = {
        "run_accession": EXPECTED_RUN,
        "proposed_group_id": EXPECTED_GROUP,
        "mapping_status": "CONTEXT_SELECTED_PENDING_COLLISION_AND_MANUAL_GATE0",
        "assay_match": selection["assay_match"],
        "anatomy_match": selection["anatomy_match"],
        "stage_match": selection["stage_match"],
        "strain_genotype_match": selection["strain_genotype_match"],
        "condition_match": selection["condition_match"],
        "reference_coordinate_match": selection["reference_coordinate_plan"],
        "normalization_compatibility": selection["normalization_plan"],
        "mapping_rationale": selection["mapping_rationale"],
        "claim_boundary": selection["claim_boundary"],
    }
    output.mkdir(parents=True)
    registry_path = output / "p18_replacement_selected_registry.tsv"
    mapping_path = output / "p18_replacement_fixed_head_mapping.tsv"
    write_tsv(registry_path, [selected], list(candidate))
    write_tsv(mapping_path, [mapping], list(mapping))
    manifest = {
        "schema_version": 1,
        "created_at": utc_now(),
        "scope": "P18 replacement Gate 0 input compilation from metadata and a represented-head context decision; no raw RNA, coverage, checkpoint, chromosome-X or locked-test signal was read",
        "candidate_registry": str(CANDIDATE_REGISTRY.relative_to(REPO_ROOT)),
        "candidate_registry_sha256": sha256(CANDIDATE_REGISTRY),
        "selection_contract": str(SELECTION_CONTRACT.relative_to(REPO_ROOT)),
        "selection_contract_sha256": sha256(SELECTION_CONTRACT),
        "group_manifest": str(GROUP_MANIFEST.relative_to(REPO_ROOT)),
        "group_manifest_sha256": sha256(GROUP_MANIFEST),
        "selected_registry": str(registry_path.relative_to(REPO_ROOT)),
        "selected_registry_sha256": sha256(registry_path),
        "fixed_head_mapping": str(mapping_path.relative_to(REPO_ROOT)),
        "fixed_head_mapping_sha256": sha256(mapping_path),
        "required_next_gate": "Run the fail-closed local collision screen and then source-equivalence attestation before any G1 download/reprocessing approval.",
        "prohibited_claims": [
            "The source is not an unseen target or condition.",
            "No external RNA signal or model result exists at this step.",
            "A selected head mapping alone does not establish external performance.",
        ],
    }
    manifest_path = output / "p18_replacement_gate0_input_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"status": "prepared", "registry": str(registry_path.relative_to(REPO_ROOT)), "mapping": str(mapping_path.relative_to(REPO_ROOT))}, indent=2))


if __name__ == "__main__":
    main()
