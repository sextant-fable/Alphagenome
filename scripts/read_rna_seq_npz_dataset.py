#!/usr/bin/env python3
"""Read back AlphaGenome-like RNA_SEQ NPZ examples and validate batch shapes."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "dataset_dir",
        nargs="?",
        default="alphagenome_custom/datasets/rna_seq_npz_pilot_train",
    )
    parser.add_argument("--batch-size", type=int, default=2)
    parser.add_argument(
        "--summary-output",
        default=None,
        help="Optional TSV summary output path.",
    )
    parser.add_argument(
        "--scan-all",
        action="store_true",
        help="Scan all examples for global RNA-seq min/max/mean and non-finite values.",
    )
    return parser.parse_args()


def read_manifest(path: Path) -> list[dict[str, str]]:
    return list(csv.DictReader(path.open(), delimiter="\t"))


def load_example(dataset_dir: Path, rel_path: str) -> dict[str, np.ndarray]:
    with np.load(dataset_dir / rel_path, allow_pickle=False) as data:
        return {key: data[key] for key in data.files}


def batch_examples(examples: list[dict[str, np.ndarray]]) -> dict[str, np.ndarray]:
    return {
        "dna_sequence": np.stack(
            [ex["dna_sequence"].astype(np.float32, copy=False) for ex in examples],
            axis=0,
        ),
        "rna_seq": np.stack(
            [ex["rna_seq"].astype(np.float32, copy=False) for ex in examples],
            axis=0,
        ),
        "rna_seq_mask": np.stack([ex["rna_seq_mask"] for ex in examples], axis=0),
        "rna_seq_strand": np.stack([ex["rna_seq_strand"] for ex in examples], axis=0),
        "interval_chromosome": np.stack([ex["interval_chromosome"] for ex in examples]),
        "interval_start": np.stack([ex["interval_start"] for ex in examples]),
        "interval_end": np.stack([ex["interval_end"] for ex in examples]),
    }


def scan_dataset(dataset_dir: Path, manifest: list[dict[str, str]]) -> dict[str, str]:
    rna_min = np.inf
    rna_max = -np.inf
    rna_sum = 0.0
    rna_count = 0
    nonfinite_count = 0

    for row in manifest:
        example = load_example(dataset_dir, row["path"])
        rna_seq = example["rna_seq"].astype(np.float32, copy=False)
        finite = np.isfinite(rna_seq)
        nonfinite_count += int((~finite).sum())
        finite_values = rna_seq if bool(finite.all()) else rna_seq[finite]
        if finite_values.size:
            rna_min = min(rna_min, float(finite_values.min()))
            rna_max = max(rna_max, float(finite_values.max()))
            rna_sum += float(finite_values.sum(dtype=np.float64))
            rna_count += int(finite_values.size)

    rna_mean = rna_sum / rna_count if rna_count else float("nan")
    return {
        "scan_all_examples": "true",
        "examples_scanned": str(len(manifest)),
        "rna_seq_all_min": f"{rna_min:.6g}",
        "rna_seq_all_max": f"{rna_max:.6g}",
        "rna_seq_all_mean": f"{rna_mean:.6g}",
        "rna_seq_nonfinite_count": str(nonfinite_count),
    }


def main() -> None:
    args = parse_args()
    dataset_dir = Path(args.dataset_dir)
    manifest = read_manifest(dataset_dir / "manifest.tsv")
    with (dataset_dir / "dataset_metadata.json").open() as handle:
        metadata = json.load(handle)
    selected = manifest[: args.batch_size]
    examples = [load_example(dataset_dir, row["path"]) for row in selected]
    batch = batch_examples(examples)

    summary = {
        "dataset_dir": str(dataset_dir),
        "n_examples_manifest": len(manifest),
        "batch_size": len(selected),
        "dna_sequence_shape": "x".join(map(str, batch["dna_sequence"].shape)),
        "rna_seq_shape": "x".join(map(str, batch["rna_seq"].shape)),
        "rna_seq_mask_shape": "x".join(map(str, batch["rna_seq_mask"].shape)),
        "rna_seq_strand_shape": "x".join(map(str, batch["rna_seq_strand"].shape)),
        "dna_dtype_after_load": str(batch["dna_sequence"].dtype),
        "rna_seq_dtype_after_load": str(batch["rna_seq"].dtype),
        "rna_seq_min": f"{float(batch['rna_seq'].min()):.6g}",
        "rna_seq_max": f"{float(batch['rna_seq'].max()):.6g}",
        "rna_seq_mean": f"{float(batch['rna_seq'].mean()):.6g}",
        "first_interval": (
            f"{str(batch['interval_chromosome'][0])}:"
            f"{int(batch['interval_start'][0])}-{int(batch['interval_end'][0])}"
        ),
        "format": metadata["format"],
    }
    if args.scan_all:
        summary.update(scan_dataset(dataset_dir, manifest))

    for key, value in summary.items():
        print(f"{key}\t{value}")

    if args.summary_output:
        output = Path(args.summary_output)
        output.parent.mkdir(parents=True, exist_ok=True)
        with output.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(summary.keys()), delimiter="\t")
            writer.writeheader()
            writer.writerow(summary)


if __name__ == "__main__":
    main()
