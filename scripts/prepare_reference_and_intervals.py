#!/usr/bin/env python3
"""Prepare C. elegans reference metadata and AlphaGenome-style intervals."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


DEFAULT_SPLITS = {
    "I": "train",
    "II": "train",
    "III": "train",
    "IV": "train",
    "V": "valid",
    "X": "test",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--base-dir", default="alphagenome_custom")
    parser.add_argument(
        "--fasta",
        default="alphagenome_custom/reference/genome.fa",
        help="Uncompressed FASTA file.",
    )
    parser.add_argument(
        "--gtf",
        default="alphagenome_custom/reference/annotation.gtf",
        help="GTF annotation file.",
    )
    parser.add_argument("--window-size", type=int, default=2**20)
    parser.add_argument("--stride", type=int, default=2**19)
    return parser.parse_args()


def fasta_index_rows(fasta_path: Path) -> list[dict[str, int | str]]:
    rows: list[dict[str, int | str]] = []
    with fasta_path.open("rb") as handle:
        name: str | None = None
        length = 0
        offset = 0
        line_bases: int | None = None
        line_width: int | None = None

        while True:
            pos = handle.tell()
            line = handle.readline()
            if not line:
                break
            if line.startswith(b">"):
                if name is not None:
                    rows.append(
                        {
                            "name": name,
                            "length": length,
                            "offset": offset,
                            "line_bases": line_bases or 0,
                            "line_width": line_width or 0,
                        }
                    )
                name = line[1:].split()[0].decode("utf-8")
                length = 0
                offset = handle.tell()
                line_bases = None
                line_width = None
            else:
                stripped = line.rstrip(b"\r\n")
                if not stripped:
                    continue
                if line_bases is None:
                    line_bases = len(stripped)
                    line_width = len(line)
                length += len(stripped)

        if name is not None:
            rows.append(
                {
                    "name": name,
                    "length": length,
                    "offset": offset,
                    "line_bases": line_bases or 0,
                    "line_width": line_width or 0,
                }
            )
    return rows


def write_fai(fasta_path: Path, rows: list[dict[str, int | str]]) -> None:
    with Path(f"{fasta_path}.fai").open("w") as handle:
        for row in rows:
            handle.write(
                "\t".join(
                    str(row[key])
                    for key in ["name", "length", "offset", "line_bases", "line_width"]
                )
                + "\n"
            )


def read_bigwig_chroms(qc_path: Path) -> dict[str, int]:
    with qc_path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if row.get("is_bigwig") == "True" and row.get("chromosomes"):
                return {
                    item.split(":")[0]: int(item.split(":")[1])
                    for item in row["chromosomes"].split(",")
                }
    return {}


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def make_windows(chrom: str, length: int, window_size: int, stride: int) -> list[tuple[str, int, int]]:
    if length < window_size:
        return []
    starts = list(range(0, length - window_size + 1, stride))
    final_start = length - window_size
    if starts[-1] != final_start:
        starts.append(final_start)
    return [(chrom, start, start + window_size) for start in sorted(set(starts))]


def summarize_gtf(gtf_path: Path, output_path: Path) -> None:
    counts: dict[tuple[str, str], int] = {}
    with gtf_path.open() as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "gene":
                continue
            chrom = fields[0]
            biotype = "unknown"
            for attr in fields[8].split(";"):
                attr = attr.strip()
                if attr.startswith("gene_biotype "):
                    biotype = attr.split('"')[1]
                    break
            counts[(chrom, biotype)] = counts.get((chrom, biotype), 0) + 1

    rows = [
        {"chromosome": chrom, "gene_biotype": biotype, "n_genes": count}
        for (chrom, biotype), count in sorted(counts.items())
    ]
    write_tsv(output_path, rows, ["chromosome", "gene_biotype", "n_genes"])


def main() -> None:
    args = parse_args()
    base_dir = Path(args.base_dir)
    fasta_path = Path(args.fasta)
    gtf_path = Path(args.gtf)
    metadata_dir = base_dir / "metadata"
    intervals_dir = base_dir / "intervals"

    fai_rows = fasta_index_rows(fasta_path)
    write_fai(fasta_path, fai_rows)

    chrom_size_rows = [
        {"chromosome": row["name"], "length": row["length"]} for row in fai_rows
    ]
    write_tsv(metadata_dir / "chrom_sizes.tsv", chrom_size_rows, ["chromosome", "length"])

    bigwig_chroms = read_bigwig_chroms(metadata_dir / "bigwig_qc.tsv")
    validation_rows = []
    fasta_chroms = {str(row["name"]): int(row["length"]) for row in fai_rows}
    for chrom in sorted(set(fasta_chroms) | set(bigwig_chroms)):
        validation_rows.append(
            {
                "chromosome": chrom,
                "fasta_length": fasta_chroms.get(chrom, ""),
                "bigwig_length": bigwig_chroms.get(chrom, ""),
                "length_match": fasta_chroms.get(chrom) == bigwig_chroms.get(chrom),
            }
        )
    write_tsv(
        metadata_dir / "reference_bigwig_validation.tsv",
        validation_rows,
        ["chromosome", "fasta_length", "bigwig_length", "length_match"],
    )

    split_rows = []
    interval_rows_by_split: dict[str, list[tuple[str, int, int]]] = {
        "train": [],
        "valid": [],
        "test": [],
        "excluded": [],
    }
    for chrom, length in fasta_chroms.items():
        split = DEFAULT_SPLITS.get(chrom, "excluded")
        windows = make_windows(chrom, length, args.window_size, args.stride)
        interval_rows_by_split[split].extend(windows)
        split_rows.append(
            {
                "chromosome": chrom,
                "length": length,
                "split": split,
                "window_size": args.window_size,
                "stride": args.stride,
                "n_intervals": len(windows),
            }
        )

    write_tsv(
        metadata_dir / "interval_split_summary.tsv",
        split_rows,
        ["chromosome", "length", "split", "window_size", "stride", "n_intervals"],
    )

    for split in ["train", "valid", "test"]:
        with (intervals_dir / f"{split}.bed").open("w") as handle:
            for chrom, start, end in interval_rows_by_split[split]:
                handle.write(f"{chrom}\t{start}\t{end}\n")

    with (intervals_dir / "all_intervals.bed").open("w") as handle:
        for split in ["train", "valid", "test"]:
            for chrom, start, end in interval_rows_by_split[split]:
                handle.write(f"{chrom}\t{start}\t{end}\t{split}\n")

    summarize_gtf(gtf_path, metadata_dir / "gtf_gene_summary.tsv")

    print(f"{fasta_path}.fai")
    print(metadata_dir / "chrom_sizes.tsv")
    print(metadata_dir / "reference_bigwig_validation.tsv")
    print(metadata_dir / "interval_split_summary.tsv")
    print(intervals_dir / "train.bed")
    print(intervals_dir / "valid.bed")
    print(intervals_dir / "test.bed")


if __name__ == "__main__":
    main()
