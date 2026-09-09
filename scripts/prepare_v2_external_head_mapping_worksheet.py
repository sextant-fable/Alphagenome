#!/usr/bin/env python3
"""Create a metadata-only worksheet for external fixed-output-head mapping.

The current model has 241 fixed RNA output heads.  An external candidate cannot
be treated as an unseen condition merely because its accession is new.  This
utility writes two immutable, signal-free tables: a candidate review worksheet
and the represented-head catalogue needed to make each mapping decision
auditable before any raw external RNA download.
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
DEFAULT_REGISTRY = REPO_ROOT / "results/v2_p14_metadata_collection_20260822T174811Z/p14_external_candidate_registry.tsv"
GROUP_CATALOGUE = REPO_ROOT / "alphagenome_custom/metadata/v2/replicate_holdout_v1/training_group_manifest.tsv"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--external-registry", type=Path, default=DEFAULT_REGISTRY)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def repo_relative(path: Path) -> str:
    return str(path.resolve().relative_to(REPO_ROOT.resolve()))


def main() -> None:
    args = parse_args()
    registry = args.external_registry.resolve()
    output = args.output_dir.resolve()
    if output.exists():
        raise FileExistsError(f"Output directory is immutable: {output}")
    if not registry.is_file() or not GROUP_CATALOGUE.is_file():
        raise FileNotFoundError("External registry or 241-head catalogue is missing")
    candidates = read_tsv(registry)
    groups = read_tsv(GROUP_CATALOGUE)
    if not candidates or not groups or len(groups) != 241:
        raise RuntimeError("Expected non-empty candidate registry and exactly 241 represented heads")

    review_rows: list[dict[str, Any]] = []
    for index, row in enumerate(candidates, start=2):
        review_rows.append(
            {
                "candidate_row": index,
                "candidate_set": row.get("candidate_set", ""),
                "candidate_label": row.get("candidate_label", ""),
                "run_accession": row.get("run_accession", ""),
                "experiment_accession": row.get("experiment_accession", ""),
                "biosample_accession": row.get("biosample_accession", ""),
                "sample_alias": row.get("sample_alias", ""),
                "sra_study_accession": row.get("sra_study_accession", ""),
                "bioproject_accession": row.get("bioproject_accession", ""),
                "geo_accession": row.get("geo_accession", ""),
                "library_source": row.get("library_source", ""),
                "library_strategy": row.get("library_strategy", ""),
                "library_layout": row.get("library_layout", ""),
                "scientific_name": row.get("scientific_name", ""),
                "sample_title": row.get("sample_title", ""),
                "experiment_title": row.get("experiment_title", ""),
                "study_title": row.get("study_title", ""),
                "study_alias": row.get("study_alias", ""),
                "first_public": row.get("first_public", ""),
                "collision_screen_status": "PENDING_MANUAL_GATE0_REVIEW",
                "proposed_group_id": "",
                "head_mapping_status": "UNASSIGNED",
                "assay_match": "",
                "anatomy_match": "",
                "stage_match": "",
                "strain_genotype_match": "",
                "condition_match": "",
                "reference_coordinate_match": "",
                "normalization_compatibility": "",
                "mapping_rationale": "",
                "manual_reviewer": "",
                "manual_reviewed_at": "",
                "final_disposition": "PENDING",
            }
        )
    head_fields = [
        "group_id", "group_status", "assay", "strand", "tissue_or_cell_type", "anatomy_curie",
        "anatomy_label", "development_stage", "life_stage_curie", "life_stage_label", "sex", "strain",
        "genotype", "condition", "sra_study_accession", "bioproject_accession", "geo_accession",
        "n_biological_units", "replicate_holdout_role", "signal_unit", "normalization_target_total",
    ]
    missing = [field for field in head_fields if field not in groups[0]]
    if missing:
        raise RuntimeError(f"Head catalogue lacks required context fields: {missing}")
    output.mkdir(parents=True)
    review_path = output / "external_candidate_head_mapping_worksheet.tsv"
    head_path = output / "represented_241_head_catalogue.tsv"
    write_tsv(review_path, review_rows, list(review_rows[0]))
    write_tsv(head_path, [{field: row[field] for field in head_fields} for row in groups], head_fields)
    manifest = {
        "schema_version": 1,
        "created_at": utc_now(),
        "scope": "metadata-only fixed-output-head mapping worksheet; no external raw RNA, coverage, checkpoint, chromosome-X or locked-test signal was read",
        "candidate_registry": repo_relative(registry),
        "candidate_registry_sha256": sha256(registry),
        "represented_head_catalogue": repo_relative(GROUP_CATALOGUE),
        "represented_head_catalogue_sha256": sha256(GROUP_CATALOGUE),
        "candidate_rows": len(review_rows),
        "represented_heads": len(groups),
        "gate_a2_rule": "A candidate may proceed only after a reviewer records one compatible represented head, complete assay/context/reference/normalization rationale, and a final disposition. Empty proposed_group_id is not a mapping.",
        "prohibited_claims": ["No candidate is currently mapped.", "No source is cleared for external RNA download.", "This worksheet does not establish unseen-condition prediction."],
    }
    manifest_path = output / "external_head_mapping_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"status": "prepared", "worksheet": repo_relative(review_path), "head_catalogue": repo_relative(head_path), "candidates": len(review_rows), "heads": len(groups)}, indent=2))


if __name__ == "__main__":
    main()
