#!/usr/bin/env python3
"""Compute leakage-safe nonzero means from registered training blocks only."""

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
INTERVAL_DIR = REPO_ROOT / "alphagenome_custom/intervals/v2"
CHROMOSOMES = ("I", "II", "III", "IV", "V", "X")


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
    path: Path, chromosomes: tuple[str, ...] = CHROMOSOMES
) -> dict[str, tuple[float, int]]:
    """Compatibility helper for per-chromosome diagnostics and unit tests."""
    result = {}
    with pyBigWig.open(str(path)) as bigwig:
        for chromosome in chromosomes:
            total = 0.0
            count = 0
            for start, end, value in bigwig.intervals(chromosome) or ():
                if not math.isfinite(value) or value < 0:
                    raise RuntimeError(
                        f"Invalid signal in {path} at {chromosome}:{start}-{end}"
                    )
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


def training_blocks(path: Path) -> list[tuple[str, int, int, str]]:
    rows = read_tsv(path)
    blocks = {
        (
            row["chromosome"],
            int(row["block_start"]),
            int(row["block_end"]),
            row["block_id"],
        )
        for row in rows
    }
    if not blocks or any(row["role"] != "train" for row in rows):
        raise RuntimeError(f"Expected non-empty train-only manifest: {path}")
    return sorted(blocks)


def nonzero_mean_for_blocks(
    path: Path, blocks: list[tuple[str, int, int, str]]
) -> float:
    return mean_from_block_totals(nonzero_totals_for_blocks(path, blocks), blocks)


def nonzero_totals_for_blocks(
    path: Path, blocks: list[tuple[str, int, int, str]]
) -> dict[tuple[str, int, int, str], tuple[float, int]]:
    result = {}
    with pyBigWig.open(str(path)) as bigwig:
        for block in blocks:
            chromosome, block_start, block_end, _ = block
            total = 0.0
            count = 0
            for start, end, value in (
                bigwig.intervals(chromosome, block_start, block_end) or ()
            ):
                if not math.isfinite(value) or value < 0:
                    raise RuntimeError(f"Invalid signal in {path} at {chromosome}:{start}-{end}")
                if value > 0:
                    width = max(0, min(end, block_end) - max(start, block_start))
                    total += width * value
                    count += width
            result[block] = (total, count)
    return result


def mean_from_block_totals(
    totals: dict[tuple[str, int, int, str], tuple[float, int]],
    blocks: list[tuple[str, int, int, str]],
) -> float:
    total = sum(totals[block][0] for block in blocks)
    count = sum(totals[block][1] for block in blocks)
    if total <= 0 or count <= 0:
        raise RuntimeError("Non-positive nonzero signal for registered blocks")
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
    if (
        registry.get("revision_id") != "six_chromosome_blocks_v1"
        or registry.get("test_chromosomes") != list(CHROMOSOMES)
    ):
        raise RuntimeError("Track means require the six-chromosome block split")
    development_blocks = training_blocks(INTERVAL_DIR / "development_train.tsv")
    fold_blocks = {
        fold: training_blocks(INTERVAL_DIR / f"fold_{fold}/train.tsv")
        for fold in range(1, 6)
    }
    rows = []
    for index, output in enumerate(outputs, start=1):
        print(f"track_mean\t{index}/{len(outputs)}\t{output['group_id']}", flush=True)
        path = REPO_ROOT / output["output_path"]
        block_totals = nonzero_totals_for_blocks(path, development_blocks)
        row: dict[str, Any] = {
            "group_id": output["group_id"],
            "development_train_nonzero_mean": (
                f"{mean_from_block_totals(block_totals, development_blocks):.12g}"
            ),
        }
        for fold in range(1, 6):
            row[f"fold_{fold}_train_nonzero_mean"] = (
                f"{mean_from_block_totals(block_totals, fold_blocks[fold]):.12g}"
            )
        rows.append(row)
    write_tsv(OUTPUT_PATH, rows)
    summary = {
        "schema_version": 2,
        "split_revision": registry["revision_id"],
        "tracks": len(rows),
        "development_mean_chromosomes": list(CHROMOSOMES),
        "fold_mean_policy": "nonzero_mean_on_registered_train_blocks_only",
        "controlled_comparison": "development_pool_train_blocks_vs_fold_train_blocks",
        "locked_test_block_signal_reads": 0,
        "development_block_count": len(development_blocks),
        "fold_train_block_counts": {
            str(fold): len(fold_blocks[fold]) for fold in range(1, 6)
        },
        "group_outputs_sha256": sha256(GROUP_OUTPUTS),
        "split_registry_sha256": sha256(SPLIT_REGISTRY),
        "output_sha256": sha256(OUTPUT_PATH),
    }
    SUMMARY_PATH.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
