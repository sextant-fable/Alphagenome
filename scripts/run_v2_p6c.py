#!/usr/bin/env python3
"""Run the separately approved one-time chromosome-X evaluation."""

from __future__ import annotations

from datetime import datetime, timezone
import json
import os
from pathlib import Path
import subprocess
import sys

from scripts import run_v2_p6b
from scripts import v2_gpu_resources


REPO_ROOT = Path(__file__).resolve().parents[1]
METADATA_DIR = REPO_ROOT / "alphagenome_custom/metadata/v2"
STATE_PATH = METADATA_DIR / "execution_state.json"
LOCK_PATH = METADATA_DIR / "final_test_lock.json"
CLAIM_PATH = METADATA_DIR / "final_test_claim.json"
REPORT_PATH = METADATA_DIR / "final_test_report.json"
EXECUTION_PATH = METADATA_DIR / "p6c_execution.json"


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def atomic_json(path: Path, payload: dict[str, object]) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)


def require_approvals() -> dict[str, object]:
    state = json.loads(STATE_PATH.read_text())
    g4 = state.get("approvals", {}).get("G4_gpu_experiments", {})
    g5 = state.get("approvals", {}).get("G5_final_test", {})
    if not (
        state.get("current_phase") == "P6C"
        and g4.get("approved") is True
        and g4.get("scope") == "r6_gpu_auto_available_2_3"
        and g5.get("approved") is True
        and g5.get("scope") == "r6c_single_chr_x_test"
    ):
        raise RuntimeError("P6C requires current phase P6C and exact scoped G4/G5")
    lock = json.loads(LOCK_PATH.read_text())
    if lock.get("test_consumed") is True:
        raise RuntimeError("Chromosome-X final test has already been consumed")
    return lock


def main() -> None:
    lock = require_approvals()
    selected, snapshots = run_v2_p6b.wait_for_gpus(70_000)
    gpu = selected[0]
    claim = {
        "schema_version": 1,
        "claimed_at": utc_now(),
        "checkpoint_path": lock["checkpoint_path"],
        "checkpoint_sha256": lock["checkpoint_sha256"],
        "physical_gpu": gpu,
        "status": "claimed",
    }
    try:
        descriptor = os.open(CLAIM_PATH, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o644)
    except FileExistsError as error:
        raise RuntimeError("Final-test claim already exists; refusing a second entry") from error
    with os.fdopen(descriptor, "w") as handle:
        handle.write(json.dumps(claim, indent=2, sort_keys=True) + "\n")
    command = [
        sys.executable,
        "scripts/evaluate_v2_model.py",
        "--model",
        str(lock["model"]),
        "--training-loss",
        str(lock["loss"]),
        "--fold",
        "0",
        "--checkpoint",
        str(lock["checkpoint_path"]),
        "--sequence-length",
        "131072",
        "--hidden-channels",
        "64",
        "--output",
        str(REPORT_PATH.relative_to(REPO_ROOT)),
        "--device",
        "cuda",
        "--final-test",
    ]
    log_path = REPO_ROOT / "logs/v2_p6c_20260714/final_test.log"
    record = {
        "schema_version": 1,
        "phase": "P6C",
        "status": "running",
        "started_at": utc_now(),
        "completed_at": None,
        "physical_gpu": gpu,
        "gpu_wait_snapshots": snapshots,
        "command": command,
        "claim_path": str(CLAIM_PATH.relative_to(REPO_ROOT)),
        "report_path": str(REPORT_PATH.relative_to(REPO_ROOT)),
        "log_path": str(log_path.relative_to(REPO_ROOT)),
    }
    atomic_json(EXECUTION_PATH, record)
    environment = os.environ.copy()
    environment["CUDA_VISIBLE_DEVICES"] = str(gpu)
    log_path.parent.mkdir(parents=True, exist_ok=True)
    try:
        with log_path.open("w") as log:
            subprocess.run(
                command,
                cwd=REPO_ROOT,
                env=environment,
                stdout=log,
                stderr=subprocess.STDOUT,
                check=True,
            )
    except Exception:
        record["status"] = "failed"
        current_lock = json.loads(LOCK_PATH.read_text())
        if current_lock.get("test_consumed") is True:
            current_lock["test_status"] = "failed_after_consumption"
            atomic_json(LOCK_PATH, current_lock)
        claim["status"] = "failed"
        atomic_json(CLAIM_PATH, claim)
        raise
    else:
        record["status"] = "completed"
        claim["status"] = "completed"
        atomic_json(CLAIM_PATH, claim)
    finally:
        record["completed_at"] = utc_now()
        record["resources_after"] = v2_gpu_resources.snapshot()
        atomic_json(EXECUTION_PATH, record)


if __name__ == "__main__":
    main()
