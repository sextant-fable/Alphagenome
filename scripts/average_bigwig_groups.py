#!/usr/bin/env python3
"""Average replicate bigWigs by AlphaGenome-style biological context groups."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
import shutil

import numpy as np
import pyBigWig


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--groups",
        default="alphagenome_custom/metadata/track_groups.tsv",
        help="Replicate group table.",
    )
    parser.add_argument(
        "--base-dir",
        default="alphagenome_custom",
        help="Base directory used by relative file paths in group table.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Overwrite existing grouped bigWigs.",
    )
    return parser.parse_args()


def read_tsv(path: str | Path) -> list[dict[str, str]]:
    return list(csv.DictReader(Path(path).open(), delimiter="\t"))


def open_bigwigs(paths: list[Path]) -> list[pyBigWig.pyBigWig]:
    return [pyBigWig.open(str(path)) for path in paths]


def validate_chromosomes(handles: list[pyBigWig.pyBigWig]) -> list[tuple[str, int]]:
    chroms = handles[0].chroms()
    for handle in handles[1:]:
        if handle.chroms() != chroms:
            raise ValueError("Input bigWigs have different chromosome definitions.")
    return list(chroms.items())


def write_average_bigwig(input_paths: list[Path], output_path: Path) -> None:
    handles = open_bigwigs(input_paths)
    try:
        header = validate_chromosomes(handles)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        tmp_path = output_path.with_suffix(output_path.suffix + ".tmp")
        if tmp_path.exists():
            tmp_path.unlink()
        out = pyBigWig.open(str(tmp_path), "w")
        out.addHeader(header)
        for chrom, length in header:
            arrays = [
                np.nan_to_num(
                    handle.values(chrom, 0, length, numpy=True),
                    nan=0.0,
                    posinf=0.0,
                    neginf=0.0,
                ).astype(np.float32, copy=False)
                for handle in handles
            ]
            mean_values = np.mean(arrays, axis=0, dtype=np.float32)
            out.addEntries(
                chrom,
                0,
                values=mean_values.tolist(),
                span=1,
                step=1,
            )
        out.close()
        tmp_path.replace(output_path)
    finally:
        for handle in handles:
            handle.close()


def main() -> None:
    args = parse_args()
    base_dir = Path(args.base_dir)
    groups = read_tsv(args.groups)
    for index, group in enumerate(groups, start=1):
        input_paths = [base_dir / p for p in group["input_files"].split(",") if p]
        output_path = base_dir / group["aggregated_file_path"]
        if output_path.exists() and not args.force:
            print(f"[{index}/{len(groups)}] exists {output_path}")
            continue
        print(
            f"[{index}/{len(groups)}] {group['group_id']} "
            f"n={len(input_paths)} -> {output_path}"
        )
        if len(input_paths) == 1:
            output_path.parent.mkdir(parents=True, exist_ok=True)
            if output_path.exists():
                output_path.unlink()
            shutil.copy2(input_paths[0], output_path)
        else:
            write_average_bigwig(input_paths, output_path)


if __name__ == "__main__":
    main()
