#!/usr/bin/env python3
"""Run lightweight QC on downloaded bigWig files."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import pyBigWig


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--manifest",
        default="remote_inventory/bigwig_metadata_manifest.tsv",
        help="TSV with filename and sample metadata.",
    )
    parser.add_argument(
        "--tracks-dir",
        default="alphagenome_custom/tracks/rna_seq",
        help="Directory containing downloaded bigWig files.",
    )
    parser.add_argument(
        "--output",
        default="alphagenome_custom/metadata/bigwig_qc.tsv",
        help="Output QC TSV.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    tracks_dir = Path(args.tracks_dir)
    rows = list(csv.DictReader(Path(args.manifest).open(), delimiter="\t"))

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "sample_id",
        "filename",
        "exists",
        "is_bigwig",
        "file_size_mb",
        "chromosomes",
        "n_bases_covered",
        "min_val",
        "max_val",
        "sum_data",
        "error",
    ]
    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            path = tracks_dir / row["filename"]
            result = {
                "sample_id": row["sample_id"],
                "filename": row["filename"],
                "exists": path.exists(),
                "is_bigwig": "",
                "file_size_mb": f"{path.stat().st_size / 1024 / 1024:.3f}"
                if path.exists()
                else "",
                "chromosomes": "",
                "n_bases_covered": "",
                "min_val": "",
                "max_val": "",
                "sum_data": "",
                "error": "",
            }
            if path.exists():
                try:
                    bw = pyBigWig.open(str(path))
                    result["is_bigwig"] = bw.isBigWig()
                    result["chromosomes"] = ",".join(
                        f"{name}:{length}" for name, length in bw.chroms().items()
                    )
                    header = bw.header()
                    result["n_bases_covered"] = header.get("nBasesCovered", "")
                    result["min_val"] = header.get("minVal", "")
                    result["max_val"] = header.get("maxVal", "")
                    result["sum_data"] = header.get("sumData", "")
                    bw.close()
                except Exception as exc:  # noqa: BLE001
                    result["error"] = f"{type(exc).__name__}: {exc}"
            writer.writerow(result)

    print(output)


if __name__ == "__main__":
    main()
