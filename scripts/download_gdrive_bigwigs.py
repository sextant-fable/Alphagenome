#!/usr/bin/env python3
"""Download Google Drive bigWig files listed in a manifest.

The manifest is expected to be a TSV with at least:
sample_id, filename, drive_id.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
import time

import gdown


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--manifest",
        default="remote_inventory/bigwig_metadata_manifest.tsv",
        help="TSV produced from the Google Drive inventory.",
    )
    parser.add_argument(
        "--output-dir",
        default="alphagenome_custom/tracks/rna_seq",
        help="Directory for downloaded bigWig files.",
    )
    parser.add_argument("--retries", type=int, default=5)
    parser.add_argument("--sleep", type=float, default=10.0)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    manifest = Path(args.manifest)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    rows = list(csv.DictReader(manifest.open(), delimiter="\t"))
    for index, row in enumerate(rows, start=1):
        output_path = output_dir / row["filename"]
        if output_path.exists() and output_path.stat().st_size > 0:
            size_mb = output_path.stat().st_size / 1024 / 1024
            print(f"[{index}/{len(rows)}] exists {output_path.name} {size_mb:.1f} MB")
            continue

        url = f"https://drive.google.com/uc?id={row['drive_id']}"
        for attempt in range(1, args.retries + 1):
            print(
                f"[{index}/{len(rows)}] downloading {output_path.name} "
                f"(attempt {attempt}/{args.retries})"
            )
            try:
                result = gdown.download(
                    url=url,
                    output=str(output_path),
                    quiet=False,
                    use_cookies=False,
                    resume=True,
                )
                if result is not None and output_path.exists():
                    size_mb = output_path.stat().st_size / 1024 / 1024
                    print(f"  done {output_path.name} {size_mb:.1f} MB")
                    break
            except Exception as exc:  # noqa: BLE001
                print(f"  error: {type(exc).__name__}: {exc}")

            if attempt == args.retries:
                raise RuntimeError(f"Failed to download {output_path.name}")
            time.sleep(args.sleep)


if __name__ == "__main__":
    main()
