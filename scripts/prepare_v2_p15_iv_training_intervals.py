#!/usr/bin/env python3
"""Create a fail-closed I--V-only train/validation interval contract for P15.

The completed P12 training manifests contain chromosome-X rows. P15 derives
new train and validation manifests before a future ablation or public-baseline
training process can open a sequence or BigWig file. This script is metadata
only and never instantiates a dataset.
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
SOURCE_ROOT = REPO_ROOT / "alphagenome_custom" / "intervals" / "v2"
ALLOWED_CHROMOSOMES = ("I", "II", "III", "IV", "V")
ALLOWED_SET = set(ALLOWED_CHROMOSOMES)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def read_tsv(path: Path) -> tuple[list[dict[str, str]], list[str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader), list(reader.fieldnames or [])


def write_tsv(path: Path, rows: list[dict[str, str]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def select_rows(
    rows: list[dict[str, str]], *, expected_role: str, source: Path
) -> list[dict[str, str]]:
    selected = [
        row
        for row in rows
        if row["chromosome"] in ALLOWED_SET and row["role"] == expected_role
    ]
    if not selected:
        raise ValueError(f"No I--V {expected_role} rows in {source}")
    if any(row["chromosome"] not in ALLOWED_SET for row in selected):
        raise AssertionError("Chromosome-X row survived the P15 filter")
    if any(row["role"] != expected_role or row["role"] == "test_locked" for row in selected):
        raise AssertionError("Unexpected interval role survived the P15 filter")
    return selected


def main() -> None:
    args = parse_args()
    output = args.output_dir.resolve()
    if output.exists():
        raise FileExistsError(f"Output directory already exists: {output}")
    output.relative_to(REPO_ROOT.resolve())
    output.mkdir(parents=True)

    folds: list[dict[str, Any]] = []
    required = {"chromosome", "role", "fold", "start", "end", "core_start", "core_end"}
    for fold in range(1, 6):
        fold_records: dict[str, Any] = {"fold": fold}
        for role in ("train", "valid"):
            source = SOURCE_ROOT / f"fold_{fold}" / f"{role}.tsv"
            rows, fields = read_tsv(source)
            if not required.issubset(fields):
                raise ValueError(f"Missing required columns in {source}: {sorted(required - set(fields))}")
            selected = select_rows(rows, expected_role=role, source=source)
            destination = output / "intervals" / f"fold_{fold}" / f"{role}.tsv"
            write_tsv(destination, selected, fields)
            fold_records[role] = {
                "source": str(source.relative_to(REPO_ROOT)),
                "source_sha256": sha256(source),
                "source_rows": len(rows),
                "source_chromosomes": sorted({row["chromosome"] for row in rows}),
                "output": str(destination.relative_to(REPO_ROOT)),
                "output_sha256": sha256(destination),
                "selected_i_to_v_rows": len(selected),
                "selected_chromosomes": sorted({row["chromosome"] for row in selected}),
                "excluded_rows": len(rows) - len(selected),
            }
        folds.append(fold_records)

    contract = {
        "schema_version": 1,
        "phase": "P15",
        "created_at": utc_now(),
        "scope": "I--V-only training and validation interval contract for future ablations and public baselines.",
        "prohibited": [
            "No DNA, BigWig, checkpoint or model file was opened.",
            "Chromosome X is excluded from every output interval manifest.",
            "Locked-test intervals are excluded from every output interval manifest.",
            "This contract does not itself establish external validation or out-of-distribution generalization.",
        ],
        "allowed_chromosomes": list(ALLOWED_CHROMOSOMES),
        "source_split_revision": "six_chromosome_blocks_v2",
        "folds": folds,
    }
    contract_path = output / "p15_iv_training_interval_contract.json"
    contract_path.write_text(json.dumps(contract, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"status": "completed", "contract": str(contract_path.relative_to(REPO_ROOT))}, indent=2))


if __name__ == "__main__":
    main()
