#!/usr/bin/env python3
"""Resolve immutable ENA FASTQ transport metadata for attested P18 sources.

This metadata-only step resolves the exact FASTQ URL, checksum, byte count and
read count for each already-attested external run. It intentionally does not
download a byte of source RNA data. The resulting manifest is the only allowed
input to the later bounded download/reprocessing sidecar.
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
    "secondary_study_accession,library_layout,library_selection,library_source,"
    "library_strategy,read_count,fastq_ftp,fastq_md5,fastq_bytes"
)
TRANSPORT_FIELDS = [
    "run_accession",
    "experiment_accession",
    "biosample_accession",
    "sra_study_accession",
    "bioproject_accession",
    "geo_accession",
    "proposed_group_id",
    "library_layout",
    "library_selection",
    "read_count",
    "fastq_ftp",
    "fastq_md5",
    "fastq_bytes",
    "source_equivalence_audit_status",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--attested-registry", type=Path, required=True)
    parser.add_argument("--fixed-head-mapping", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=TRANSPORT_FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def repo_relative(path: Path) -> str:
    return str(path.resolve().relative_to(REPO_ROOT.resolve()))


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def fetch_run(run_accession: str) -> tuple[dict[str, str], str]:
    query = urlencode({"accession": run_accession, "result": "read_run", "fields": ENA_FIELDS, "format": "tsv", "download": "false"})
    url = f"{ENA_URL}?{query}"
    with urlopen(url, timeout=60) as response:
        rows = list(csv.DictReader(io.StringIO(response.read().decode("utf-8")), delimiter="\t"))
    if len(rows) != 1 or rows[0].get("run_accession") != run_accession:
        raise RuntimeError(f"Expected exactly one ENA read-run record for {run_accession}")
    return rows[0], url


def validate_transport(row: dict[str, str]) -> None:
    layout = row["library_layout"]
    expected_files = 2 if layout == "PAIRED" else 1 if layout == "SINGLE" else 0
    if expected_files == 0:
        raise RuntimeError(f"Unsupported library layout for {row['run_accession']}: {layout}")
    urls = row["fastq_ftp"].split(";")
    md5s = row["fastq_md5"].split(";")
    bytes_values = row["fastq_bytes"].split(";")
    if not (len(urls) == len(md5s) == len(bytes_values) == expected_files):
        raise RuntimeError(f"ENA FASTQ list is inconsistent for {row['run_accession']}")
    if any(not value.startswith("ftp.sra.ebi.ac.uk/") for value in urls):
        raise RuntimeError(f"Unexpected ENA FASTQ host for {row['run_accession']}")
    if any(len(value) != 32 for value in md5s) or any(int(value) <= 0 for value in bytes_values):
        raise RuntimeError(f"Invalid ENA FASTQ checksum or byte count for {row['run_accession']}")
    if int(row["read_count"]) <= 0:
        raise RuntimeError(f"Invalid ENA read count for {row['run_accession']}")


def main() -> None:
    args = parse_args()
    registry = args.attested_registry.resolve()
    mapping = args.fixed_head_mapping.resolve()
    output = args.output_dir.resolve()
    if output.exists():
        raise FileExistsError(f"Output directory is immutable: {output}")
    if not registry.is_file() or not mapping.is_file():
        raise FileNotFoundError("Attested registry or fixed-head mapping is missing")
    registry_rows = read_tsv(registry)
    mapping_rows = read_tsv(mapping)
    mapping_by_run = {row.get("run_accession", ""): row for row in mapping_rows}
    if not registry_rows or set(mapping_by_run) != {row.get("run_accession", "") for row in registry_rows}:
        raise RuntimeError("Registry and fixed-head mapping must contain the same non-empty run set")
    output.mkdir(parents=True)
    rows: list[dict[str, str]] = []
    queries: list[dict[str, str]] = []
    for source in sorted(registry_rows, key=lambda row: row["run_accession"]):
        if source.get("source_equivalence_audit_status") != "cleared_no_equivalent_source":
            raise RuntimeError(f"Source-equivalence attestation is missing: {source['run_accession']}")
        ena, url = fetch_run(source["run_accession"])
        if ena.get("experiment_accession") != source.get("experiment_accession") or ena.get("sample_accession") != source.get("biosample_accession"):
            raise RuntimeError(f"ENA identity changed since Gate 0 for {source['run_accession']}")
        transport = {
            "run_accession": source["run_accession"],
            "experiment_accession": ena["experiment_accession"],
            "biosample_accession": ena["sample_accession"],
            "sra_study_accession": ena["secondary_study_accession"],
            "bioproject_accession": ena["study_accession"],
            "geo_accession": source["geo_accession"],
            "proposed_group_id": mapping_by_run[source["run_accession"]]["proposed_group_id"],
            "library_layout": ena["library_layout"],
            "library_selection": ena["library_selection"],
            "read_count": ena["read_count"],
            "fastq_ftp": ena["fastq_ftp"],
            "fastq_md5": ena["fastq_md5"],
            "fastq_bytes": ena["fastq_bytes"],
            "source_equivalence_audit_status": source["source_equivalence_audit_status"],
        }
        validate_transport(transport)
        rows.append(transport)
        queries.append({"run_accession": source["run_accession"], "query_url": url})
    manifest_path = output / "p18_external_source_transport_manifest.json"
    table_path = output / "p18_external_source_transport.tsv"
    write_tsv(table_path, rows)
    total_bytes = sum(sum(int(value) for value in row["fastq_bytes"].split(";")) for row in rows)
    payload = {
        "schema_version": 1,
        "created_at": utc_now(),
        "scope": "ENA transport metadata only; no FASTQ, coverage, checkpoint, chromosome-X or locked-test signal was read",
        "attested_registry": repo_relative(registry),
        "attested_registry_sha256": sha256(registry),
        "fixed_head_mapping": repo_relative(mapping),
        "fixed_head_mapping_sha256": sha256(mapping),
        "transport_manifest": repo_relative(table_path),
        "transport_manifest_sha256": sha256(table_path),
        "source_runs": len(rows),
        "source_fastq_bytes": total_bytes,
        "queries": queries,
        "required_next_gate": "A separate G1 scope is required before any source FASTQ download or uniform reprocessing.",
        "prohibited_claims": ["No external RNA data have been downloaded.", "No external performance result exists.", "No coordinate compatibility has been empirically verified."],
    }
    manifest_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"status": "prepared", "runs": len(rows), "fastq_bytes": total_bytes, "transport_manifest": repo_relative(table_path)}, indent=2))


if __name__ == "__main__":
    main()
