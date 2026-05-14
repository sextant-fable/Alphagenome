#!/usr/bin/env python3
"""Smoke-test extraction of AlphaGenome-style DNA and RNA-seq targets."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np
import pyBigWig


DNA_TO_INDEX = {"A": 0, "C": 1, "G": 2, "T": 3}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--base-dir", default="alphagenome_custom")
    parser.add_argument("--interval-bed", default="alphagenome_custom/intervals/train.bed")
    parser.add_argument("--output", default="alphagenome_custom/metadata/training_input_smoke_test.tsv")
    return parser.parse_args()


def read_fai(path: Path) -> dict[str, tuple[int, int, int, int]]:
    result = {}
    with path.open() as handle:
        for line in handle:
            name, length, offset, line_bases, line_width = line.rstrip("\n").split("\t")
            result[name] = (int(length), int(offset), int(line_bases), int(line_width))
    return result


def fetch_fasta_sequence(fasta_path: Path, fai: dict[str, tuple[int, int, int, int]], chrom: str, start: int, end: int) -> str:
    length, offset, line_bases, line_width = fai[chrom]
    if start < 0 or end > length:
        raise ValueError(f"Interval out of bounds: {chrom}:{start}-{end}")
    seq_len = end - start
    start_line = start // line_bases
    start_in_line = start % line_bases
    byte_start = offset + start_line * line_width + start_in_line
    # Read enough bytes to span line breaks, then remove newlines.
    n_lines = ((start_in_line + seq_len - 1) // line_bases) + 1
    byte_count = seq_len + n_lines * (line_width - line_bases)
    with fasta_path.open("rb") as handle:
        handle.seek(byte_start)
        raw = handle.read(byte_count)
    return raw.replace(b"\n", b"").replace(b"\r", b"")[:seq_len].decode("ascii").upper()


def one_hot_encode(seq: str) -> np.ndarray:
    arr = np.zeros((len(seq), 4), dtype=np.float32)
    for index, base in enumerate(seq):
        if base in DNA_TO_INDEX:
            arr[index, DNA_TO_INDEX[base]] = 1.0
    return arr


def main() -> None:
    args = parse_args()
    base_dir = Path(args.base_dir)
    fasta_path = base_dir / "reference/genome.fa"
    fai = read_fai(base_dir / "reference/genome.fa.fai")

    with Path(args.interval_bed).open() as handle:
        chrom, start_s, end_s = handle.readline().rstrip("\n").split("\t")[:3]
    start = int(start_s)
    end = int(end_s)

    sequence = fetch_fasta_sequence(fasta_path, fai, chrom, start, end)
    dna = one_hot_encode(sequence)

    track_rows = list(
        csv.DictReader((base_dir / "metadata/track_metadata_grouped.tsv").open(), delimiter="\t")
    )
    target_arrays = []
    for row in track_rows:
        bw_path = base_dir / row["file_path"]
        bw = pyBigWig.open(str(bw_path))
        values = bw.values(chrom, start, end, numpy=True)
        bw.close()
        target_arrays.append(np.nan_to_num(values, nan=0.0).astype(np.float32))
    targets = np.stack(target_arrays, axis=1)

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "chromosome",
                "start",
                "end",
                "sequence_length",
                "dna_shape",
                "rna_seq_shape",
                "n_tracks",
                "dna_non_acgt_bases",
                "rna_seq_min",
                "rna_seq_max",
                "rna_seq_mean",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerow(
            {
                "chromosome": chrom,
                "start": start,
                "end": end,
                "sequence_length": len(sequence),
                "dna_shape": "x".join(map(str, dna.shape)),
                "rna_seq_shape": "x".join(map(str, targets.shape)),
                "n_tracks": len(track_rows),
                "dna_non_acgt_bases": int((dna.sum(axis=1) == 0).sum()),
                "rna_seq_min": f"{float(targets.min()):.6g}",
                "rna_seq_max": f"{float(targets.max()):.6g}",
                "rna_seq_mean": f"{float(targets.mean()):.6g}",
            }
        )
    print(out)


if __name__ == "__main__":
    main()
