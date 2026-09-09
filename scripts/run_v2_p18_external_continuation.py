#!/usr/bin/env python3
"""Execute the preregistered P18 external-RNA continuation.

P18 reuses three source FASTQs only after fresh MD5 verification against the
new immutable continuation manifest, downloads the attested replacement, then
uniformly reprocesses all four sources to I-V-only WBcel235 bigWigs. This phase
does not load a model, perform inference, or access chromosome X.
"""

from __future__ import annotations

import csv
import hashlib
import json
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[1]
METADATA_ROOT = REPO_ROOT / "alphagenome_custom/metadata/v2"
SPEC_PATH = METADATA_ROOT / "p18_external_continuation_spec.json"
EXECUTION_PATH = METADATA_ROOT / "p18_external_continuation_execution.json"
SCOPE = "p18_external_replacement_download_and_reprocessing"
ORIGINAL_AUDIT_DIR = REPO_ROOT / "results/v2_p18_external_download"
ORIGINAL_RUNS = {"SRR18463404", "SRR18463405", "SRR18463406"}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def md5(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def require_contract() -> tuple[dict[str, Any], Path, list[dict[str, str]]]:
    from scripts import v2_phase_controller as controller

    state = controller.load_state()
    if state.get("current_phase") != "P18" or state.get("status") != "RUNNING":
        raise RuntimeError("P18 continuation requires controller state P18/RUNNING")
    if not controller.approval_satisfies(state, "G1", SCOPE):
        raise RuntimeError(f"P18 continuation requires G1:{SCOPE}")
    if not SPEC_PATH.is_file():
        raise FileNotFoundError(f"P18 preregistration is missing: {SPEC_PATH}")
    spec = json.loads(SPEC_PATH.read_text(encoding="utf-8"))
    if spec.get("phase") != "P18" or spec.get("status") != "preregistered":
        raise RuntimeError("P18 preregistration is invalid")
    if spec.get("allowed_chromosomes") != ["I", "II", "III", "IV", "V"] or spec.get("final_test_access") != "prohibited":
        raise RuntimeError("P18 does not satisfy the I-V/no-final-test contract")
    transport = REPO_ROOT / spec["source"]["continuation_transport"]["path"]
    if not transport.is_file() or sha256(transport) != spec["source"]["continuation_transport"]["sha256"]:
        raise RuntimeError("P18 continuation transport manifest is missing or altered")
    rows = read_tsv(transport)
    expected = set(spec.get("expected_runs", []))
    if not expected or {row.get("run_accession", "") for row in rows} != expected or len(rows) != len(expected):
        raise RuntimeError("P18 continuation transport does not match its preregistered source set")
    if any(row.get("source_equivalence_audit_status") != "cleared_no_equivalent_source" for row in rows):
        raise RuntimeError("P18 continuation source-equivalence evidence is incomplete")
    return spec, transport, rows


def rebind_verified_originals(rows: list[dict[str, str]], transport: Path, audit_dir: Path) -> None:
    """Bind previously verified source bytes to the new immutable manifest."""

    transport_hash = sha256(transport)
    by_run = {row["run_accession"]: row for row in rows}
    audit_dir.mkdir(parents=True, exist_ok=True)
    for run in sorted(ORIGINAL_RUNS):
        row = by_run[run]
        original_path = ORIGINAL_AUDIT_DIR / f"{run}.json"
        if not original_path.is_file():
            raise FileNotFoundError(f"Missing original verified-download audit: {original_path}")
        original = json.loads(original_path.read_text(encoding="utf-8"))
        paths = [REPO_ROOT / value for value in original.get("local_fastq_paths", [])]
        if len(paths) != 1 or not paths[0].is_file():
            raise RuntimeError(f"Original verified FASTQ is unavailable for {run}")
        if md5(paths[0]) != row.get("fastq_md5"):
            raise RuntimeError(f"Fresh MD5 verification failed for reused source {run}")
        audit_path = audit_dir / f"{run}.json"
        rebound = {
            "schema_version": 1,
            "run_accession": run,
            "proposed_group_id": row["proposed_group_id"],
            "transport_manifest": str(transport.relative_to(REPO_ROOT)),
            "transport_manifest_sha256": transport_hash,
            "source_fastq_urls": row["fastq_ftp"],
            "source_fastq_expected_md5": row["fastq_md5"],
            "source_fastq_actual_md5": row["fastq_md5"],
            "source_fastq_sha256": sha256(paths[0]),
            "source_fastq_bytes": int(row["fastq_bytes"]),
            "local_fastq_paths": [str(paths[0].relative_to(REPO_ROOT))],
            "source_origin": str(original_path.relative_to(REPO_ROOT)),
            "verification": "P18 fresh MD5 verification before reuse under the continuation manifest",
            "completed_at": utc_now(),
            "claim_boundary": "A prior checksum-verified source was rebound to a new external-RNA manifest; no coverage or model result was produced.",
        }
        audit_path.write_text(json.dumps(rebound, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def run(command: list[str]) -> None:
    print("command\t" + " ".join(command), flush=True)
    subprocess.run(command, cwd=REPO_ROOT, check=True)


def main() -> None:
    spec, transport, rows = require_contract()
    output = spec["output"]
    audit_dir = REPO_ROOT / output["download_audit_dir"]
    coverage_dir = REPO_ROOT / output["coverage_dir"]
    if EXECUTION_PATH.exists():
        previous = json.loads(EXECUTION_PATH.read_text(encoding="utf-8"))
        if previous.get("transport_manifest_sha256") == sha256(transport) and previous.get("status") == "completed":
            print(json.dumps({"status": "already_completed", "execution": str(EXECUTION_PATH.relative_to(REPO_ROOT))}, indent=2))
            return
        raise RuntimeError("P18 execution record exists but is not the same completed contract")
    started = utc_now()
    rebind_verified_originals(rows, transport, audit_dir)
    common = [item for run in sorted(row["run_accession"] for row in rows) for item in ("--expected-run", run)]
    run([
        sys.executable, "-m", "scripts.run_v2_p18_external_rna_download",
        "--transport-manifest", str(transport.relative_to(REPO_ROOT)),
        "--work-dir", output["fastq_root"],
        "--audit-dir", output["download_audit_dir"],
        "--controller-phase", "P18", "--approval-scope", SCOPE,
        *common,
    ])
    run([
        sys.executable, "-m", "scripts.run_v2_p18_external_iv_reprocessing",
        "--transport-manifest", str(transport.relative_to(REPO_ROOT)),
        "--download-audit-dir", output["download_audit_dir"],
        "--fastq-root", output["fastq_root"],
        "--work-dir", output["work_dir"],
        "--output-dir", output["coverage_dir"],
        "--threads", "8", "--controller-phase", "P18", "--approval-scope", SCOPE,
        *common,
    ])
    summary_path = coverage_dir / "p18_reprocessing_summary.json"
    if not summary_path.is_file():
        raise RuntimeError("P18 reprocessing did not produce its required summary")
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    expected = sorted(row["run_accession"] for row in rows)
    if summary.get("completed_runs") != len(expected) or summary.get("completed_accessions") != expected:
        raise RuntimeError("P18 reprocessing summary is incomplete")
    execution = {
        "schema_version": 1,
        "phase": "P18",
        "status": "completed",
        "controller_invocation": "python -m scripts.v2_phase_controller run --phase P18",
        "registered_phase_command": "python -m scripts.run_v2_p18_external_continuation",
        "started_at": started,
        "completed_at": utc_now(),
        "transport_manifest": str(transport.relative_to(REPO_ROOT)),
        "transport_manifest_sha256": sha256(transport),
        "expected_runs": expected,
        "coverage_dir": output["coverage_dir"],
        "coverage_summary": str(summary_path.relative_to(REPO_ROOT)),
        "coverage_summary_sha256": sha256(summary_path),
        "resource_contract": spec["resource_contract"],
        "final_test_access": "prohibited",
        "locked_test_block_signal_reads": 0,
        "model_inference": "not performed",
    }
    EXECUTION_PATH.write_text(json.dumps(execution, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"status": "completed", "execution": str(EXECUTION_PATH.relative_to(REPO_ROOT))}, indent=2))


if __name__ == "__main__":
    main()
