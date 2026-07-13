#!/usr/bin/env python3
"""Run the approved 482-run P3B streaming reprocessing phase."""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
import shutil
import subprocess


REPO_ROOT = Path(__file__).resolve().parents[1]
STATE_PATH = REPO_ROOT / "alphagenome_custom/metadata/v2/execution_state.json"
ENV_NAME = "alphagenome"
MIN_FREE_BYTES = 200 * 1024**3
RUN_RECORD = REPO_ROOT / "alphagenome_custom/metadata/v2/p3b_execution.json"
REQUIRED_APPROVALS = {
    "G1_large_source_download_or_realignment": "p3_full_rna_streaming_482",
    "G2_v2_manifest_scientific_review": "v2_manifest_candidate_hierarchy",
    "G3_dataset_or_cache_generation": "p3_full_outputs_and_p4_loader",
}


def require_scoped_approvals() -> None:
    state = json.loads(STATE_PATH.read_text())
    if state.get("current_phase") != "P3B":
        raise RuntimeError(
            f"P3B may run only while current_phase=P3B, got {state.get('current_phase')}"
        )
    missing = []
    for key, scope in REQUIRED_APPROVALS.items():
        record = state.get("approvals", {}).get(key, {})
        if record.get("approved") is not True or record.get("scope") != scope:
            missing.append(f"{key}:{scope}")
    if missing:
        raise RuntimeError("Missing scoped approvals: " + ",".join(missing))


def main() -> None:
    require_scoped_approvals()
    conda = shutil.which("conda")
    if conda is None:
        raise RuntimeError("conda is required")
    free_bytes = shutil.disk_usage(REPO_ROOT).free
    if free_bytes < MIN_FREE_BYTES:
        raise RuntimeError(
            f"P3B requires at least {MIN_FREE_BYTES} free bytes; found {free_bytes}"
        )
    commands = [
        [
            conda,
            "run",
            "--no-capture-output",
            "-n",
            ENV_NAME,
            "python",
            "scripts/build_ncbi_sra_sources_v2.py",
            "--workers",
            "8",
        ],
        [
            conda,
            "run",
            "--no-capture-output",
            "-n",
            ENV_NAME,
            "python",
            "scripts/run_v2_reprocessing_full.py",
            "--threads-per-sample",
            "16",
            "--workers",
            "4",
            "--work-dir",
            "shared/source_reads/v2/full_streaming",
            "--output-dir",
            "alphagenome_custom/tracks/rna_seq_v2_normalized",
            "--star-index",
            "shared/reference_indexes/WBcel235_STAR_2.7.11b",
            "--download-backend",
            "ncbi_sra",
        ],
        [
            conda,
            "run",
            "--no-capture-output",
            "-n",
            ENV_NAME,
            "python",
            "scripts/finalize_v2_reprocessed_groups.py",
        ],
    ]
    log_path = REPO_ROOT / "logs/v2_p3b_20260714/p3b_full.log"
    log_path.parent.mkdir(parents=True, exist_ok=True)
    started_at = datetime.now(timezone.utc).replace(microsecond=0).isoformat()
    record = {
        "phase": "P3B",
        "status": "running",
        "started_at": started_at,
        "completed_at": None,
        "commands": commands,
        "log_path": str(log_path.relative_to(REPO_ROOT)),
        "workers": 4,
        "threads_per_sample": 16,
        "download_backend": "ncbi_sra",
        "sra_toolkit_version": "3.4.1",
        "sra_toolkit_archive_sha256": "b950362c054765a4184af41947f022f040e94e964862017c0ecb0b0273db3596",
        "free_bytes_at_start": free_bytes,
    }
    RUN_RECORD.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
    try:
        with log_path.open("a") as log:
            for command in commands:
                print("command\t" + " ".join(command), flush=True)
                log.write("command\t" + " ".join(command) + "\n")
                log.flush()
                subprocess.run(
                    command,
                    cwd=REPO_ROOT,
                    stdout=log,
                    stderr=subprocess.STDOUT,
                    check=True,
                )
    except Exception:
        record["status"] = "failed"
        raise
    else:
        record["status"] = "completed"
    finally:
        record["completed_at"] = (
            datetime.now(timezone.utc).replace(microsecond=0).isoformat()
        )
        RUN_RECORD.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
