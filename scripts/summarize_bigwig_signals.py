#!/usr/bin/env python3
"""Summarize bigWig signal scales and replicate groups."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path
import statistics

import pyBigWig


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--manifest",
        default="remote_inventory/bigwig_metadata_manifest.tsv",
        help="Per-file metadata manifest.",
    )
    parser.add_argument(
        "--groups",
        default="alphagenome_custom/metadata/track_groups.tsv",
        help="Replicate group table.",
    )
    parser.add_argument(
        "--tracks-dir",
        default="alphagenome_custom/tracks/rna_seq",
        help="Directory containing input bigWigs.",
    )
    parser.add_argument(
        "--sample-output",
        default="alphagenome_custom/metadata/bigwig_signal_summary.tsv",
    )
    parser.add_argument(
        "--group-output",
        default="alphagenome_custom/metadata/group_signal_summary.tsv",
    )
    return parser.parse_args()


def read_tsv(path: str | Path) -> list[dict[str, str]]:
    return list(csv.DictReader(Path(path).open(), delimiter="\t"))


def write_tsv(path: str | Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    args = parse_args()
    tracks_dir = Path(args.tracks_dir)
    manifest = read_tsv(args.manifest)
    groups = read_tsv(args.groups)

    sample_rows: list[dict[str, object]] = []
    sample_by_id: dict[str, dict[str, object]] = {}
    for row in manifest:
        bw_path = tracks_dir / row["filename"]
        bw = pyBigWig.open(str(bw_path))
        header = bw.header()
        chroms = bw.chroms()
        bw.close()
        n_bases = int(header["nBasesCovered"])
        sum_data = float(header["sumData"])
        mean_signal = sum_data / n_bases if n_bases else 0.0
        out = {
            "sample_id": row["sample_id"],
            "filename": row["filename"],
            "output_type": row["output_type"],
            "tissue": row["tissue"],
            "stage": row["stage"],
            "n_chromosomes": len(chroms),
            "n_bases_covered": n_bases,
            "sum_data": int(sum_data),
            "mean_signal": f"{mean_signal:.6g}",
            "min_val": header["minVal"],
            "max_val": header["maxVal"],
            "file_size_mb": f"{bw_path.stat().st_size / 1024 / 1024:.3f}",
        }
        sample_rows.append(out)
        sample_by_id[row["sample_id"]] = out

    sample_fields = [
        "sample_id",
        "filename",
        "output_type",
        "tissue",
        "stage",
        "n_chromosomes",
        "n_bases_covered",
        "sum_data",
        "mean_signal",
        "min_val",
        "max_val",
        "file_size_mb",
    ]
    write_tsv(args.sample_output, sample_rows, sample_fields)

    group_rows: list[dict[str, object]] = []
    for group in groups:
        sample_ids = [s for s in group["sample_ids"].split(",") if s]
        means = [float(sample_by_id[s]["mean_signal"]) for s in sample_ids]
        sums = [int(sample_by_id[s]["sum_data"]) for s in sample_ids]
        mean_of_means = statistics.mean(means)
        cv = (
            statistics.stdev(means) / mean_of_means
            if len(means) > 1 and mean_of_means
            else 0.0
        )
        group_rows.append(
            {
                "group_id": group["group_id"],
                "output_type": group["output_type"],
                "tissue": group["tissue"],
                "stage": group["stage"],
                "n_replicates": len(sample_ids),
                "sample_ids": ",".join(sample_ids),
                "sample_mean_signals": ",".join(f"{m:.6g}" for m in means),
                "mean_signal_after_average": f"{mean_of_means:.6g}",
                "replicate_mean_cv": f"{cv:.6g}",
                "sample_sum_data": ",".join(str(s) for s in sums),
            }
        )

    group_fields = [
        "group_id",
        "output_type",
        "tissue",
        "stage",
        "n_replicates",
        "sample_ids",
        "sample_mean_signals",
        "mean_signal_after_average",
        "replicate_mean_cv",
        "sample_sum_data",
    ]
    write_tsv(args.group_output, group_rows, group_fields)
    print(args.sample_output)
    print(args.group_output)


if __name__ == "__main__":
    main()
