#!/usr/bin/env python3
"""Create the fail-closed I--V-only interval contract for P13.

The completed P12 validation manifests contain chromosome-X development cores.
This utility copies only non-X, non-locked validation rows into a new immutable
output directory before any BigWig data are opened. It never reads signal data.
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
ALLOWED_CHROMOSOMES = {"I", "II", "III", "IV", "V"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


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
        raise FileExistsError(f"Output directory already exists: {output}")
    output.relative_to(REPO_ROOT.resolve())

    output.mkdir(parents=True)
    folds: list[dict[str, Any]] = []
    for fold in range(1, 6):
        source = SOURCE_ROOT / f"fold_{fold}" / "valid.tsv"
        rows, fields = read_tsv(source)
        required = {"chromosome", "role", "fold", "start", "end", "core_start", "core_end"}
        if not required.issubset(fields):
            raise ValueError(f"Missing required columns in {source}: {sorted(required - set(fields))}")
        selected = [
            row
            for row in rows
            if row["chromosome"] in ALLOWED_CHROMOSOMES and row["role"] == "valid"
        ]
        rejected = [row for row in rows if row not in selected]
        if not selected:
            raise ValueError(f"No I--V validation rows in {source}")
        if any(row["chromosome"] not in ALLOWED_CHROMOSOMES for row in selected):
            raise AssertionError("Chromosome-X row survived the P13 filter")
        if any(row["role"] == "test_locked" for row in selected):
            raise AssertionError("Locked-test row survived the P13 filter")
        destination = output / "intervals" / f"fold_{fold}" / "valid.tsv"
        write_tsv(destination, selected, fields)
        folds.append(
            {
                "fold": fold,
                "source": str(source.relative_to(REPO_ROOT)),
                "source_sha256": sha256(source),
                "output": str(destination.relative_to(REPO_ROOT)),
                "output_sha256": sha256(destination),
                "source_rows": len(rows),
                "selected_i_to_v_rows": len(selected),
                "excluded_rows": len(rejected),
                "selected_chromosomes": sorted({row["chromosome"] for row in selected}),
                "excluded_chromosomes": sorted({row["chromosome"] for row in rejected}),
            }
        )

    contract = {
        "schema_version": 1,
        "phase": "P13",
        "created_at": utc_now(),
        "scope": "I--V-only replicate-agreement reference; no BigWig signal was read while creating this contract",
        "source_interval_root": str(SOURCE_ROOT.relative_to(REPO_ROOT)),
        "allowed_chromosomes": sorted(ALLOWED_CHROMOSOMES),
        "prohibited": ["chromosome X", "test_locked rows", "locked-test signal", "model selection"],
        "folds": folds,
    }
    contract_path = output / "p13_iv_interval_contract.json"
    contract_path.write_text(json.dumps(contract, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"status": "prepared", "contract": str(contract_path), "folds": folds}, indent=2))


if __name__ == "__main__":
    main()
