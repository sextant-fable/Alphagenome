#!/usr/bin/env python3
"""Resume-safe, checksum-verified download for the four P18 external RNA sources.

This sidecar consumes only the immutable P18 transport manifest. It is guarded
by a distinct G1 approval recorded while P17 remains the active serial matrix.
It never opens a bigWig, model, chromosome-X interval or locked-test asset.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import shutil
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[1]
P18_SCOPE = "p18_external_rna_download_and_reprocessing"
EXPECTED_RUNS = {"SRR18463404", "SRR18463405", "SRR18463406", "SRR2005820"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--transport-manifest", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument("--audit-dir", type=Path, required=True)
    parser.add_argument("--expected-run", action="append", default=[])
    parser.add_argument("--controller-phase", default="P17")
    parser.add_argument("--approval-scope", default=P18_SCOPE)
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def require_scoped_approval(phase: str, scope: str) -> None:
    from scripts import v2_phase_controller as controller

    state = controller.load_state()
    controller.validate_approval_scope("G1", scope, phase)
    if state.get("current_phase") != phase or state.get("status") != "RUNNING":
        raise RuntimeError(f"External download requires the active {phase} RUNNING period")
    if not controller.approval_satisfies(state, "G1", scope):
        raise RuntimeError(f"External download requires G1:{scope}")


def valid_row(row: dict[str, str], expected_runs: set[str]) -> None:
    run = row.get("run_accession", "")
    if run not in expected_runs:
        raise RuntimeError(f"Unexpected P18 run: {run}")
    if row.get("source_equivalence_audit_status") != "cleared_no_equivalent_source":
        raise RuntimeError(f"Source-equivalence status is not clear for {run}")
    layout = row.get("library_layout")
    expected_files = 1 if layout == "SINGLE" else 2 if layout == "PAIRED" else 0
    urls = row.get("fastq_ftp", "").split(";")
    md5s = row.get("fastq_md5", "").split(";")
    bytes_values = row.get("fastq_bytes", "").split(";")
    if expected_files == 0 or not (len(urls) == len(md5s) == len(bytes_values) == expected_files):
        raise RuntimeError(f"Invalid source transport fields for {run}")
    if any(not url.startswith("ftp.sra.ebi.ac.uk/") for url in urls):
        raise RuntimeError(f"Unexpected remote host for {run}")
    if any(len(value) != 32 for value in md5s) or any(int(value) <= 0 for value in bytes_values):
        raise RuntimeError(f"Invalid checksum or byte metadata for {run}")


def main() -> None:
    args = parse_args()
    manifest = args.transport_manifest.resolve()
    work_dir = args.work_dir.resolve()
    audit_dir = args.audit_dir.resolve()
    if not manifest.is_file():
        raise FileNotFoundError(f"P18 transport manifest is missing: {manifest}")
    for path in (work_dir, audit_dir):
        path.relative_to(REPO_ROOT.resolve())
    if shutil.disk_usage(REPO_ROOT).free < 25 * 1024**3:
        raise RuntimeError("P18 external source download requires at least 25 GiB free disk")
    expected_runs = set(args.expected_run) if args.expected_run else EXPECTED_RUNS
    require_scoped_approval(args.controller_phase, args.approval_scope)
    rows = read_tsv(manifest)
    if {row.get("run_accession", "") for row in rows} != expected_runs or len(rows) != len(expected_runs):
        raise RuntimeError("External transport manifest does not contain exactly the registered runs")
    for row in rows:
        valid_row(row, expected_runs)
    work_dir.mkdir(parents=True, exist_ok=True)
    audit_dir.mkdir(parents=True, exist_ok=True)

    from scripts.run_v2_reprocessing_full import SampleRunner, download_fastq

    completed: list[dict[str, Any]] = []
    for row in sorted(rows, key=lambda value: value["run_accession"]):
        require_scoped_approval(args.controller_phase, args.approval_scope)
        run = row["run_accession"]
        sample_dir = work_dir / run
        sample_dir.mkdir(parents=True, exist_ok=True)
        audit_path = audit_dir / f"{run}.json"
        if audit_path.is_file():
            existing = json.loads(audit_path.read_text(encoding="utf-8"))
            if existing.get("transport_manifest_sha256") == sha256(manifest):
                files = [Path(path) for path in existing.get("local_fastq_paths", [])]
                if files and all(path.is_file() for path in files):
                    completed.append(existing)
                    continue
        runner = SampleRunner(sample_dir)
        started = utc_now()
        local_paths: list[str] = []
        source_sha256: list[str] = []
        source_md5: list[str] = []
        for index, (url, expected_md5, expected_bytes) in enumerate(
            zip(row["fastq_ftp"].split(";"), row["fastq_md5"].split(";"), row["fastq_bytes"].split(";")), start=1
        ):
            filename = Path(url).name
            local_path = sample_dir / filename
            actual_sha, actual_md5 = download_fastq(runner, url, expected_md5, int(expected_bytes), local_path)
            if actual_md5 != expected_md5:
                raise RuntimeError(f"Unexpected post-download MD5 for {run} file {index}")
            local_paths.append(str(local_path.relative_to(REPO_ROOT)))
            source_sha256.append(actual_sha)
            source_md5.append(actual_md5)
        audit = {
            "schema_version": 1,
            "run_accession": run,
            "proposed_group_id": row["proposed_group_id"],
            "transport_manifest": str(manifest.relative_to(REPO_ROOT)),
            "transport_manifest_sha256": sha256(manifest),
            "source_fastq_urls": row["fastq_ftp"],
            "source_fastq_expected_md5": row["fastq_md5"],
            "source_fastq_actual_md5": ";".join(source_md5),
            "source_fastq_sha256": ";".join(source_sha256),
            "source_fastq_bytes": sum(int(value) for value in row["fastq_bytes"].split(";")),
            "local_fastq_paths": local_paths,
            "commands": runner.commands,
            "peak_work_bytes": runner.peak_bytes,
            "started_at": started,
            "completed_at": utc_now(),
            "claim_boundary": "Source files were integrity-checked for an external represented-context evaluation candidate. No coverage, model inference or external performance result is created by this transfer.",
        }
        temporary = audit_path.with_suffix(".json.tmp")
        temporary.write_text(json.dumps(audit, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        temporary.replace(audit_path)
        completed.append(audit)
        progress = {
            "schema_version": 1,
            "scope": f"{args.controller_phase} checksum-verified external FASTQ transfer only",
            "expected_runs": sorted(expected_runs),
            "completed_runs": len(completed),
            "completed_accessions": [item["run_accession"] for item in completed],
            "updated_at": utc_now(),
        }
        (audit_dir / "p18_download_progress.json").write_text(json.dumps(progress, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        print(json.dumps({"run_accession": run, "status": "download_verified", "completed": len(completed), "expected": len(rows)}, sort_keys=True), flush=True)
    summary = {
        "schema_version": 1,
        "scope": f"{args.controller_phase} checksum-verified external FASTQ transfer only",
        "transport_manifest": str(manifest.relative_to(REPO_ROOT)),
        "transport_manifest_sha256": sha256(manifest),
        "completed_runs": len(completed),
        "completed_accessions": sorted(item["run_accession"] for item in completed),
        "completed_at": utc_now(),
        "next_required_gate": "Uniform reprocessing and post-processing checksum audit remain separate steps.",
    }
    (audit_dir / "p18_download_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"status": "complete", "runs": len(completed), "audit_dir": str(audit_dir.relative_to(REPO_ROOT))}, sort_keys=True))


if __name__ == "__main__":
    main()
