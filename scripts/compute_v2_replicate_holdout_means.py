#!/usr/bin/env python3
"""Compute fold-train means for the replicate-heldout training labels."""

from __future__ import annotations

import argparse
import csv
from concurrent.futures import ThreadPoolExecutor, as_completed
import hashlib
import json
from pathlib import Path
from typing import Any

from scripts import compute_v2_track_means as base


REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INTERVAL_DIR = REPO_ROOT / "alphagenome_custom/intervals/v2"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(8 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def write_tsv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--track-manifest",
        default="alphagenome_custom/metadata/v2/replicate_holdout_v1/training_track_manifest.tsv",
    )
    parser.add_argument(
        "--interval-dir",
        default="alphagenome_custom/intervals/v2",
    )
    parser.add_argument(
        "--interval-root",
        help="Optional root containing fold_<n>/train.tsv manifests. Omits the development-wide mean.",
    )
    parser.add_argument(
        "--allowed-chromosomes",
        nargs="+",
        help="Required with --interval-root. Reject every train row outside this exact set.",
    )
    parser.add_argument(
        "--output",
        default="alphagenome_custom/metadata/v2/replicate_holdout_v1/track_nonzero_means.tsv",
    )
    parser.add_argument(
        "--summary",
        default="alphagenome_custom/metadata/v2/replicate_holdout_v1/track_nonzero_means_summary.json",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    track_manifest = REPO_ROOT / args.track_manifest
    interval_dir = REPO_ROOT / args.interval_dir
    output_path = REPO_ROOT / args.output
    summary_path = REPO_ROOT / args.summary
    outputs = read_tsv(track_manifest)
    if len(outputs) != 241:
        raise RuntimeError(f"Expected 241 training track rows, got {len(outputs)}")
    group_ids = [row["group_id"] for row in outputs]
    if group_ids != sorted(group_ids):
        raise RuntimeError("Replicate-heldout training track manifest is not sorted")
    strict_interval_root = REPO_ROOT / args.interval_root if args.interval_root else None
    if strict_interval_root is not None and not args.allowed_chromosomes:
        raise ValueError("--allowed-chromosomes is required with --interval-root")
    allowed_chromosomes = tuple(args.allowed_chromosomes or ())
    allowed_set = set(allowed_chromosomes)
    if strict_interval_root is not None and len(allowed_set) != len(allowed_chromosomes):
        raise ValueError("--allowed-chromosomes must not contain duplicates")

    def strict_train_blocks(path: Path) -> list[tuple[str, int, int, str]]:
        rows = read_tsv(path)
        if not rows:
            raise RuntimeError(f"Expected a non-empty train manifest: {path}")
        if any(row["role"] != "train" for row in rows):
            raise RuntimeError(f"Expected train-only rows: {path}")
        observed = {row["chromosome"] for row in rows}
        if not observed <= allowed_set or "X" in observed or any(row["role"] == "test_locked" for row in rows):
            raise RuntimeError(f"Strict interval root contains prohibited rows: {path}")
        return base.training_blocks(path)

    development_blocks = (
        None
        if strict_interval_root is not None
        else base.training_blocks(interval_dir / "development_train.tsv")
    )
    fold_blocks = {
        fold: (
            strict_train_blocks(strict_interval_root / f"fold_{fold}" / "train.tsv")
            if strict_interval_root is not None
            else base.training_blocks(interval_dir / f"fold_{fold}/train.tsv")
        )
        for fold in range(1, 6)
    }
    all_blocks = sorted(
        set(() if development_blocks is None else development_blocks).union(
            *(set(blocks) for blocks in fold_blocks.values())
        )
    )

    def compute_one(output: dict[str, str]) -> dict[str, Any]:
        path = REPO_ROOT / output["output_path"]
        totals = base.nonzero_totals_for_blocks(path, all_blocks)
        row: dict[str, Any] = {"group_id": output["group_id"]}
        if development_blocks is not None:
            row["development_train_nonzero_mean"] = f"{base.mean_from_block_totals(totals, development_blocks):.12g}"
        for fold in range(1, 6):
            row[f"fold_{fold}_train_nonzero_mean"] = f"{base.mean_from_block_totals(totals, fold_blocks[fold]):.12g}"
        return row

    rows_by_group: dict[str, dict[str, Any]] = {}
    with ThreadPoolExecutor(max_workers=8) as executor:
        futures = {executor.submit(compute_one, output): output["group_id"] for output in outputs}
        for index, future in enumerate(as_completed(futures), start=1):
            group_id = futures[future]
            rows_by_group[group_id] = future.result()
            print(f"track_mean\t{index}/{len(outputs)}\t{group_id}", flush=True)
    rows = [rows_by_group[output["group_id"]] for output in outputs]
    write_tsv(output_path, rows)
    summary = {
        "schema_version": 1,
        "contract": "strict_interval_root" if strict_interval_root is not None else "replicate_holdout_v1",
        "tracks": len(rows),
        "fold_mean_policy": "nonzero_mean_on_replicate_holdout_training_labels_and_registered_train_blocks_only",
        "locked_test_block_signal_reads": 0,
        "track_manifest": str(track_manifest.relative_to(REPO_ROOT)),
        "track_manifest_sha256": sha256(track_manifest),
        "interval_dir": str((strict_interval_root or interval_dir).relative_to(REPO_ROOT)),
        "allowed_chromosomes": list(allowed_chromosomes) if strict_interval_root is not None else None,
        "output_sha256": sha256(output_path),
        "development_block_count": None if development_blocks is None else len(development_blocks),
        "fold_train_block_counts": {str(fold): len(fold_blocks[fold]) for fold in range(1, 6)},
    }
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
