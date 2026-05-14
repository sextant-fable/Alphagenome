#!/usr/bin/env python3
"""Build an AlphaGenome-like RNA_SEQ dataset in per-interval NPZ files.

This is a lightweight custom format for C. elegans experiments. It keeps field
names aligned with AlphaGenome's RNA_SEQ training batch schema:

- dna_sequence: [S, 4], uint8 one-hot on disk
- rna_seq: [S, C], float32 on disk by default; float16 is optional
- rna_seq_mask: [1, C], bool
- rna_seq_strand: [1, C], int32
- interval/chromosome, interval/start, interval/end
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import time

import numpy as np
import pyBigWig


DNA_TO_INDEX = {"A": 0, "C": 1, "G": 2, "T": 3}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--base-dir", default="alphagenome_custom")
    parser.add_argument("--split", choices=["train", "valid", "test"], default="train")
    parser.add_argument(
        "--output-dir",
        default=None,
        help="Dataset output directory. Defaults to alphagenome_custom/datasets/rna_seq_npz_<split>.",
    )
    parser.add_argument(
        "--max-intervals",
        type=int,
        default=None,
        help="Build only the first N intervals for smoke/pilot datasets.",
    )
    parser.add_argument(
        "--target-dtype",
        choices=["float16", "float32"],
        default="float32",
        help="On-disk dtype for rna_seq targets.",
    )
    parser.add_argument(
        "--compress",
        action="store_true",
        help="Use compressed NPZ. Smaller but slower.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Overwrite existing examples.",
    )
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    return list(csv.DictReader(path.open(), delimiter="\t"))


def read_fai(path: Path) -> dict[str, tuple[int, int, int, int]]:
    result = {}
    with path.open() as handle:
        for line in handle:
            name, length, offset, line_bases, line_width = line.rstrip("\n").split("\t")
            result[name] = (int(length), int(offset), int(line_bases), int(line_width))
    return result


def fetch_fasta_sequence(
    fasta_path: Path,
    fai: dict[str, tuple[int, int, int, int]],
    chrom: str,
    start: int,
    end: int,
) -> str:
    length, offset, line_bases, line_width = fai[chrom]
    if start < 0 or end > length:
        raise ValueError(f"Interval out of bounds: {chrom}:{start}-{end}")
    seq_len = end - start
    start_line = start // line_bases
    start_in_line = start % line_bases
    byte_start = offset + start_line * line_width + start_in_line
    n_lines = ((start_in_line + seq_len - 1) // line_bases) + 1
    byte_count = seq_len + n_lines * (line_width - line_bases)
    with fasta_path.open("rb") as handle:
        handle.seek(byte_start)
        raw = handle.read(byte_count)
    return raw.replace(b"\n", b"").replace(b"\r", b"")[:seq_len].decode("ascii").upper()


def one_hot_encode_uint8(seq: str) -> np.ndarray:
    arr = np.zeros((len(seq), 4), dtype=np.uint8)
    for index, base in enumerate(seq):
        if base in DNA_TO_INDEX:
            arr[index, DNA_TO_INDEX[base]] = 1
    return arr


def read_intervals(path: Path, max_intervals: int | None) -> list[tuple[str, int, int]]:
    intervals = []
    with path.open() as handle:
        for line in handle:
            if not line.strip():
                continue
            chrom, start, end = line.rstrip("\n").split("\t")[:3]
            intervals.append((chrom, int(start), int(end)))
            if max_intervals is not None and len(intervals) >= max_intervals:
                break
    return intervals


def write_manifest(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "example_id",
        "split",
        "chromosome",
        "start",
        "end",
        "width",
        "path",
        "n_tracks",
        "dna_dtype",
        "rna_seq_dtype",
    ]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def clean_and_cast_signal(
    values: np.ndarray,
    target_dtype: type[np.floating],
    context: str,
) -> np.ndarray:
    cleaned = np.nan_to_num(values, nan=0.0, posinf=0.0, neginf=0.0)
    if target_dtype == np.float16:
        max_abs = float(np.max(np.abs(cleaned), initial=0.0))
        float16_max = float(np.finfo(np.float16).max)
        if max_abs > float16_max:
            raise ValueError(
                f"{context} has signal {max_abs:g}, above float16 max "
                f"{float16_max:g}; rerun with --target-dtype float32"
            )
    return cleaned.astype(target_dtype, copy=False)


def main() -> None:
    args = parse_args()
    base_dir = Path(args.base_dir)
    output_dir = (
        Path(args.output_dir)
        if args.output_dir is not None
        else base_dir / "datasets" / f"rna_seq_npz_{args.split}"
    )
    examples_dir = output_dir / "examples"
    examples_dir.mkdir(parents=True, exist_ok=True)

    fasta_path = base_dir / "reference" / "genome.fa"
    fai = read_fai(base_dir / "reference" / "genome.fa.fai")
    intervals = read_intervals(base_dir / "intervals" / f"{args.split}.bed", args.max_intervals)
    track_rows = read_tsv(base_dir / "metadata" / "track_metadata_grouped.tsv")
    bw_handles = [pyBigWig.open(str(base_dir / row["file_path"])) for row in track_rows]

    target_dtype = np.float16 if args.target_dtype == "float16" else np.float32
    save_fn = np.savez_compressed if args.compress else np.savez
    manifest_rows: list[dict[str, object]] = []
    start_time = time.time()
    try:
        for idx, (chrom, start, end) in enumerate(intervals):
            example_id = f"{args.split}_{idx:06d}_{chrom}_{start}_{end}"
            output_path = examples_dir / f"{example_id}.npz"
            if output_path.exists() and not args.force:
                print(f"[{idx + 1}/{len(intervals)}] exists {output_path.name}")
            else:
                print(f"[{idx + 1}/{len(intervals)}] building {output_path.name}")
                seq = fetch_fasta_sequence(fasta_path, fai, chrom, start, end)
                dna_sequence = one_hot_encode_uint8(seq)
                target_arrays = []
                for row, handle in zip(track_rows, bw_handles, strict=True):
                    values = handle.values(chrom, start, end, numpy=True)
                    track_label = (
                        row.get("track_id")
                        or row.get("group_id")
                        or Path(row["file_path"]).stem
                    )
                    context = f"{track_label} {chrom}:{start}-{end}"
                    target_arrays.append(clean_and_cast_signal(values, target_dtype, context))
                rna_seq = np.stack(target_arrays, axis=1)
                rna_seq_mask = np.ones((1, len(track_rows)), dtype=bool)
                # 0 means unstranded in this custom RNA_SEQ-only dataset.
                rna_seq_strand = np.zeros((1, len(track_rows)), dtype=np.int32)
                save_fn(
                    output_path,
                    dna_sequence=dna_sequence,
                    rna_seq=rna_seq,
                    rna_seq_mask=rna_seq_mask,
                    rna_seq_strand=rna_seq_strand,
                    interval_chromosome=np.array(chrom),
                    interval_start=np.array(start, dtype=np.int64),
                    interval_end=np.array(end, dtype=np.int64),
                )

            manifest_rows.append(
                {
                    "example_id": example_id,
                    "split": args.split,
                    "chromosome": chrom,
                    "start": start,
                    "end": end,
                    "width": end - start,
                    "path": str(output_path.relative_to(output_dir)),
                    "n_tracks": len(track_rows),
                    "dna_dtype": "uint8",
                    "rna_seq_dtype": args.target_dtype,
                }
            )
    finally:
        for handle in bw_handles:
            handle.close()

    write_manifest(output_dir / "manifest.tsv", manifest_rows)
    metadata = {
        "format": "alphagenome_custom_rna_seq_npz",
        "split": args.split,
        "base_dir": str(base_dir),
        "output_dir": str(output_dir),
        "n_examples": len(manifest_rows),
        "sequence_length": intervals[0][2] - intervals[0][1] if intervals else None,
        "n_tracks": len(track_rows),
        "dna_sequence_disk_dtype": "uint8",
        "rna_seq_disk_dtype": args.target_dtype,
        "rna_seq_mask_shape": [1, len(track_rows)],
        "rna_seq_strand_shape": [1, len(track_rows)],
        "rna_seq_overflow_policy": "error_on_float16_overflow",
        "track_metadata": str(base_dir / "metadata" / "track_metadata_grouped.tsv"),
        "intervals": str(base_dir / "intervals" / f"{args.split}.bed"),
        "compressed_npz": bool(args.compress),
        "elapsed_seconds": round(time.time() - start_time, 3),
    }
    with (output_dir / "dataset_metadata.json").open("w") as handle:
        json.dump(metadata, handle, indent=2)
        handle.write("\n")
    print(output_dir / "manifest.tsv")
    print(output_dir / "dataset_metadata.json")


if __name__ == "__main__":
    main()
