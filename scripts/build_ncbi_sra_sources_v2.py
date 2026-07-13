#!/usr/bin/env python3
"""Resolve and lock NCBI full-quality SRA transport objects for P3B."""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import subprocess

try:
    from scripts.run_v2_reprocessing_full import DOH_ARGS, select_sra_archive
except ModuleNotFoundError:
    from run_v2_reprocessing_full import DOH_ARGS, select_sra_archive


REPO_ROOT = Path(__file__).resolve().parents[1]
SOURCE_MANIFEST = REPO_ROOT / "alphagenome_custom/metadata/v2/p3_full_sources.tsv"
OUTPUT_PATH = REPO_ROOT / "alphagenome_custom/metadata/v2/p3_ncbi_sra_sources.tsv"
SUMMARY_PATH = REPO_ROOT / "alphagenome_custom/metadata/v2/p3_ncbi_sra_summary.json"
EXPECTED_RUNS = 482


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(8 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def read_accessions() -> list[str]:
    with SOURCE_MANIFEST.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    accessions = [row["run_accession"] for row in rows]
    if len(accessions) != EXPECTED_RUNS or len(set(accessions)) != EXPECTED_RUNS:
        raise RuntimeError("P3 full source manifest is not the exact 482-run scope")
    return accessions


def resolve(accession: str) -> dict[str, str | int]:
    url = (
        "https://locate.ncbi.nlm.nih.gov/sdl/2/retrieve"
        f"?acc={accession}&accept-alternate-locations=yes"
    )
    completed = subprocess.run(
        [
            "curl",
            "-fsSL",
            "--connect-timeout",
            "30",
            *DOH_ARGS,
            url,
        ],
        cwd=REPO_ROOT,
        capture_output=True,
        check=True,
        text=True,
    )
    payload = json.loads(completed.stdout)
    archive = select_sra_archive(payload, accession)
    return {
        "run_accession": accession,
        "sra_url": archive["url"],
        "sra_md5": archive["md5"],
        "sra_bytes": archive["size"],
        "sdl_query_url": url,
    }


def write_tsv(path: Path, rows: list[dict[str, str | int]]) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=8)
    args = parser.parse_args()
    accessions = read_accessions()
    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        rows = list(executor.map(resolve, accessions))
    rows.sort(key=lambda row: str(row["run_accession"]))
    if len(rows) != EXPECTED_RUNS:
        raise RuntimeError(f"Expected {EXPECTED_RUNS} resolved SRA objects")
    write_tsv(OUTPUT_PATH, rows)
    summary = {
        "schema_version": 1,
        "resolved_at": datetime.now(timezone.utc).replace(microsecond=0).isoformat(),
        "source_manifest": str(SOURCE_MANIFEST.relative_to(REPO_ROOT)),
        "source_manifest_sha256": sha256(SOURCE_MANIFEST),
        "sra_manifest": str(OUTPUT_PATH.relative_to(REPO_ROOT)),
        "sra_manifest_sha256": sha256(OUTPUT_PATH),
        "run_count": len(rows),
        "total_sra_bytes": sum(int(row["sra_bytes"]) for row in rows),
        "transport": "NCBI SDL full-quality SRA public ODP S3",
    }
    SUMMARY_PATH.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
