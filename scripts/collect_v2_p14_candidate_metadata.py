#!/usr/bin/env python3
"""Collect only public ENA metadata for candidate P14 sources.

This utility issues small ENA Portal metadata requests. It does not download
FASTQ, BAM, coverage, model weights or any biological signal. Its output is a
candidate registry for the separate fail-closed collision screen.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen


REPO_ROOT = Path(__file__).resolve().parents[1]
ENA_URL = "https://www.ebi.ac.uk/ena/portal/api/filereport"
ENA_FIELDS = (
    "run_accession,experiment_accession,sample_accession,study_accession,"
    "secondary_study_accession,sample_alias,fastq_md5,submitted_md5,"
    "library_source,library_strategy,library_layout,scientific_name,"
    "sample_title,experiment_title,study_title,study_alias,first_public"
)
REGISTRY_FIELDS = [
    "candidate_set",
    "candidate_label",
    "query_accession",
    "run_accession",
    "experiment_accession",
    "biosample_accession",
    "sample_alias",
    "sra_study_accession",
    "bioproject_accession",
    "geo_accession",
    "fastq_md5",
    "processed_sha256",
    "source_equivalence_audit_status",
    "library_source",
    "library_strategy",
    "library_layout",
    "scientific_name",
    "sample_title",
    "experiment_title",
    "study_title",
    "study_alias",
    "first_public",
]


def parse_candidate(value: str) -> tuple[str, str, str]:
    try:
        candidate_set, label, accession = value.split(":", 2)
    except ValueError as error:
        raise argparse.ArgumentTypeError(
            "--candidate must be SET:LABEL:ACCESSION"
        ) from error
    if not all((candidate_set, label, accession)):
        raise argparse.ArgumentTypeError("Candidate fields must be non-empty")
    return candidate_set, label, accession


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidate", type=parse_candidate, action="append", required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def fetch_rows(accession: str) -> tuple[list[dict[str, str]], str]:
    query = urlencode(
        {
            "accession": accession,
            "result": "read_run",
            "fields": ENA_FIELDS,
            "format": "tsv",
            "download": "false",
            "limit": "0",
        }
    )
    url = f"{ENA_URL}?{query}"
    with urlopen(url, timeout=60) as response:
        payload = response.read().decode("utf-8")
    rows = list(csv.DictReader(io.StringIO(payload), delimiter="\t"))
    if not rows:
        raise RuntimeError(f"ENA returned no read-run rows for {accession}")
    return rows, url


def main() -> None:
    args = parse_args()
    output = args.output_dir.resolve()
    if output.exists():
        raise FileExistsError(f"Output directory already exists: {output}")
    output.relative_to(REPO_ROOT.resolve())

    output.mkdir(parents=True)
    registry_rows: list[dict[str, str]] = []
    queries: list[dict[str, object]] = []
    for candidate_set, label, accession in args.candidate:
        rows, url = fetch_rows(accession)
        eligible = [
            row
            for row in rows
            if row.get("scientific_name", "").casefold() == "caenorhabditis elegans"
            and row.get("library_source", "").casefold() == "transcriptomic"
            and row.get("library_strategy", "").casefold() == "rna-seq"
        ]
        if not eligible:
            raise RuntimeError(f"No C. elegans transcriptomic RNA-seq rows for {accession}")
        for row in eligible:
            registry_rows.append(
                {
                    "candidate_set": candidate_set,
                    "candidate_label": label,
                    "query_accession": accession,
                    "run_accession": row.get("run_accession", ""),
                    "experiment_accession": row.get("experiment_accession", ""),
                    "biosample_accession": row.get("sample_accession", ""),
                    "sample_alias": row.get("sample_alias", ""),
                    "sra_study_accession": row.get("secondary_study_accession", ""),
                    "bioproject_accession": row.get("study_accession", ""),
                    "geo_accession": "",
                    "fastq_md5": row.get("fastq_md5", "") or row.get("submitted_md5", ""),
                    "processed_sha256": "",
                    "source_equivalence_audit_status": "pending_manual_equivalence_audit",
                    "library_source": row.get("library_source", ""),
                    "library_strategy": row.get("library_strategy", ""),
                    "library_layout": row.get("library_layout", ""),
                    "scientific_name": row.get("scientific_name", ""),
                    "sample_title": row.get("sample_title", ""),
                    "experiment_title": row.get("experiment_title", ""),
                    "study_title": row.get("study_title", ""),
                    "study_alias": row.get("study_alias", ""),
                    "first_public": row.get("first_public", ""),
                }
            )
        queries.append(
            {
                "candidate_set": candidate_set,
                "candidate_label": label,
                "accession": accession,
                "query_url": url,
                "all_read_run_rows": len(rows),
                "eligible_rna_seq_rows": len(eligible),
            }
        )

    registry = output / "p14_external_candidate_registry.tsv"
    with registry.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=REGISTRY_FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(sorted(registry_rows, key=lambda row: (row["candidate_set"], row["candidate_label"], row["run_accession"])))
    manifest = {
        "schema_version": 1,
        "created_at": utc_now(),
        "scope": "public accession and checksum metadata only; no raw or processed biological signal was downloaded",
        "queries": queries,
        "registry": str(registry.relative_to(REPO_ROOT)),
        "registry_sha256": sha256(registry),
        "rows": len(registry_rows),
        "manual_next_step": "Complete source-equivalence adjudication and collision review before any raw RNA download.",
    }
    manifest_path = output / "p14_metadata_collection_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"status": "metadata_collected", "registry": str(registry), "manifest": str(manifest_path), "rows": len(registry_rows)}, indent=2))


if __name__ == "__main__":
    main()
