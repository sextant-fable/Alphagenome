#!/usr/bin/env python3
"""Compute leakage-safe nonzero means for each I-V training fold."""

from __future__ import annotations

import csv
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import pyBigWig


REPO_ROOT = Path(__file__).resolve().parents[1]
METADATA_DIR = REPO_ROOT / "alphagenome_custom/metadata/v2"
GROUP_OUTPUTS = METADATA_DIR / "p3_group_outputs.tsv"
SPLIT_REGISTRY = METADATA_DIR / "split_registry_v2.json"
OUTPUT_PATH = METADATA_DIR / "track_nonzero_means_v2.tsv"
SUMMARY_PATH = METADATA_DIR / "track_nonzero_means_v2_summary.json"
CV_CHROMOSOMES = ("I", "II", "III", "IV", "V")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(8 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def nonzero_totals(
    path: Path, chromosomes: tuple[str, ...] = CV_CHROMOSOMES
) -> dict[str, tuple[float, int]]:
    result = {}
    with pyBigWig.open(str(path)) as bigwig:
        for chromosome in chromosomes:
            total = 0.0
            count = 0
            for start, end, value in bigwig.intervals(chromosome) or ():
                if not math.isfinite(value) or value < 0:
                    raise RuntimeError(f"Invalid signal in {path} at {chromosome}:{start}-{end}")
                if value > 0:
                    width = end - start
                    total += width * value
                    count += width
            result[chromosome] = (total, count)
    return result


def mean_for_chromosomes(
    totals: dict[str, tuple[float, int]], chromosomes: list[str] | tuple[str, ...]
) -> float:
    total = sum(totals[chromosome][0] for chromosome in chromosomes)
    count = sum(totals[chromosome][1] for chromosome in chromosomes)
    if total <= 0 or count <= 0:
        raise RuntimeError(f"Non-positive nonzero signal for chromosomes {chromosomes}")
    return total / count


def write_tsv(path: Path, rows: list[dict[str, Any]]) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def main() -> None:
    outputs = read_tsv(GROUP_OUTPUTS)
    registry = json.loads(SPLIT_REGISTRY.read_text())
    if len(outputs) != 241:
        raise RuntimeError(f"Expected 241 group outputs, got {len(outputs)}")
    fold_train = {
        int(row["fold"]): row["train_chromosomes"] for row in registry["folds"]
    }
    rows = []
    for index, output in enumerate(outputs, start=1):
        print(f"track_mean\t{index}/{len(outputs)}\t{output['group_id']}", flush=True)
        path = REPO_ROOT / output["output_path"]
        totals = nonzero_totals(path)
        row: dict[str, Any] = {
            "group_id": output["group_id"],
            "development_I_V_nonzero_mean": f"{mean_for_chromosomes(totals, CV_CHROMOSOMES):.12g}",
        }
        for fold in range(1, 6):
            row[f"fold_{fold}_train_nonzero_mean"] = (
                f"{mean_for_chromosomes(totals, fold_train[fold]):.12g}"
            )
        rows.append(row)
    write_tsv(OUTPUT_PATH, rows)
    summary = {
        "schema_version": 1,
        "tracks": len(rows),
        "development_mean_chromosomes": list(CV_CHROMOSOMES),
        "fold_mean_policy": "nonzero_mean_on_four_training_chromosomes_only",
        "controlled_comparison": "development_I_V_mean_vs_fold_train_only_mean",
        "chromosome_x_read": False,
        "group_outputs_sha256": sha256(GROUP_OUTPUTS),
        "split_registry_sha256": sha256(SPLIT_REGISTRY),
        "output_sha256": sha256(OUTPUT_PATH),
    }
    SUMMARY_PATH.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
