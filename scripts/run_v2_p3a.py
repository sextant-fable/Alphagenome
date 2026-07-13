#!/usr/bin/env python3
"""Prepare the isolated environment and run the approved five-run P3A pilot."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import shutil
import subprocess
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[1]
STATE_PATH = REPO_ROOT / "alphagenome_custom/metadata/v2/execution_state.json"
ENV_NAME = "alphagenome-rnaseq-v2"
MIN_FREE_BYTES = 50 * 1024**3
REQUIRED_APPROVALS = {
    "G1_large_source_download_or_realignment": "p3a_five_run_pilot",
    "G2_v2_manifest_scientific_review": "v2_manifest_candidate_hierarchy",
    "G3_dataset_or_cache_generation": "p3a_five_run_pilot_outputs",
}
ENV_PACKAGES = [
    "python=3.12",
    "star=2.7.11b",
    "samtools=1.24",
    "bedtools=2.31.1",
    "ucsc-bedgraphtobigwig",
    "pybigwig=0.3.25",
    "numpy",
    "curl",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--plan-only",
        action="store_true",
        help="Print the immutable command plan without checking approvals or writing data.",
    )
    return parser.parse_args()


def command_plan(conda: str) -> dict[str, Any]:
    create_command = [
        conda,
        "create",
        "-y",
        "-n",
        ENV_NAME,
        "-c",
        "conda-forge",
        "-c",
        "bioconda",
        *ENV_PACKAGES,
    ]
    pilot_command = [
        conda,
        "run",
        "--no-capture-output",
        "-n",
        ENV_NAME,
        "python",
        "scripts/run_v2_reprocessing_pilot.py",
        "--manifest",
        "alphagenome_custom/metadata/v2/p3_pilot_sources.tsv",
        "--threads",
        "16",
        "--work-dir",
        "shared/source_reads/v2/pilot_20260713",
        "--output-dir",
        "alphagenome_custom/tracks/rna_seq_v2_normalized_pilot",
        "--star-index",
        "shared/reference_indexes/WBcel235_STAR_2.7.11b",
    ]
    return {
        "phase": "P3A",
        "environment": ENV_NAME,
        "minimum_free_bytes": MIN_FREE_BYTES,
        "create_environment_command": create_command,
        "pilot_command": pilot_command,
    }


def require_scoped_approvals() -> None:
    if not STATE_PATH.is_file():
        raise RuntimeError(f"Missing workflow state: {STATE_PATH}")
    state = json.loads(STATE_PATH.read_text())
    if state.get("current_phase") != "P3A":
        raise RuntimeError(
            f"P3A may run only while current_phase=P3A, got {state.get('current_phase')}"
        )
    missing = []
    for key, scope in REQUIRED_APPROVALS.items():
        record = state.get("approvals", {}).get(key, {})
        if record.get("approved") is not True or record.get("scope") != scope:
            missing.append(f"{key}:{scope}")
    if missing:
        raise RuntimeError("Missing scoped approvals: " + ",".join(missing))


def environment_exists(conda: str) -> bool:
    completed = subprocess.run(
        [conda, "env", "list", "--json"],
        cwd=REPO_ROOT,
        check=True,
        capture_output=True,
        text=True,
    )
    environments = json.loads(completed.stdout).get("envs", [])
    return any(Path(path).name == ENV_NAME for path in environments)


def run(command: list[str]) -> None:
    print("command\t" + " ".join(command), flush=True)
    subprocess.run(command, cwd=REPO_ROOT, check=True)


def main() -> None:
    args = parse_args()
    conda = shutil.which("conda")
    if conda is None:
        raise RuntimeError("conda is required for the isolated P3A environment")
    plan = command_plan(conda)
    if args.plan_only:
        print(json.dumps(plan, indent=2, sort_keys=True))
        return

    require_scoped_approvals()
    free_bytes = shutil.disk_usage(REPO_ROOT).free
    if free_bytes < MIN_FREE_BYTES:
        raise RuntimeError(
            f"P3A requires at least {MIN_FREE_BYTES} free bytes; found {free_bytes}"
        )
    if not environment_exists(conda):
        run(plan["create_environment_command"])
    run(plan["pilot_command"])


if __name__ == "__main__":
    main()
