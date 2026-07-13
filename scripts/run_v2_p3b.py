#!/usr/bin/env python3
"""Run the approved 482-run P3B streaming reprocessing phase."""

from __future__ import annotations

import json
from pathlib import Path
import shutil
import subprocess


REPO_ROOT = Path(__file__).resolve().parents[1]
STATE_PATH = REPO_ROOT / "alphagenome_custom/metadata/v2/execution_state.json"
ENV_NAME = "alphagenome"
MIN_FREE_BYTES = 200 * 1024**3
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
    command = [
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
    ]
    print("command\t" + " ".join(command), flush=True)
    subprocess.run(command, cwd=REPO_ROOT, check=True)


if __name__ == "__main__":
    main()
