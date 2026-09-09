#!/usr/bin/env python3
"""Derive primary and supplementary held-out manifest views without reading signal."""

from __future__ import annotations

import csv
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
ROOT = REPO_ROOT / "alphagenome_custom/metadata/v2/replicate_holdout_v1"


def read(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write(path: Path, rows: list[dict[str, str]]) -> None:
    if not rows:
        raise RuntimeError(f"Refusing to write empty role manifest: {path}")
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    heldout = read(ROOT / "heldout_track_manifest.tsv")
    groups = {row["group_id"]: row for row in read(ROOT / "heldout_group_manifest.tsv")}
    indices = read(ROOT / "heldout_track_indices.tsv")
    by_group = {row["group_id"]: row for row in indices}
    roles = {
        "primary": {"primary_candidate", "primary_pending_qc"},
        "supplementary": {"supplementary_candidate"},
    }
    for label, allowed in roles.items():
        selected = [row for row in heldout if by_group[row["group_id"]]["assignment_role"] in allowed]
        selected_groups = [groups[row["group_id"]] for row in selected]
        selected_indices = [by_group[row["group_id"]] for row in selected]
        write(ROOT / f"heldout_{label}_track_manifest.tsv", selected)
        write(ROOT / f"heldout_{label}_group_manifest.tsv", selected_groups)
        write(ROOT / f"heldout_{label}_track_indices.tsv", selected_indices)
    payload = json.loads((ROOT / "p11_replicate_holdout_execution.json").read_text())
    payload["role_manifest_views"] = {
        "primary": {
            "track_manifest": "alphagenome_custom/metadata/v2/replicate_holdout_v1/heldout_primary_track_manifest.tsv",
            "group_manifest": "alphagenome_custom/metadata/v2/replicate_holdout_v1/heldout_primary_group_manifest.tsv",
            "track_indices": "alphagenome_custom/metadata/v2/replicate_holdout_v1/heldout_primary_track_indices.tsv",
            "track_count": len([row for row in heldout if by_group[row["group_id"]]["assignment_role"] in roles["primary"]]),
        },
        "supplementary": {
            "track_manifest": "alphagenome_custom/metadata/v2/replicate_holdout_v1/heldout_supplementary_track_manifest.tsv",
            "group_manifest": "alphagenome_custom/metadata/v2/replicate_holdout_v1/heldout_supplementary_group_manifest.tsv",
            "track_indices": "alphagenome_custom/metadata/v2/replicate_holdout_v1/heldout_supplementary_track_indices.tsv",
            "track_count": len([row for row in heldout if by_group[row["group_id"]]["assignment_role"] in roles["supplementary"]]),
        },
    }
    (ROOT / "p11_replicate_holdout_execution.json").write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps(payload["role_manifest_views"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
