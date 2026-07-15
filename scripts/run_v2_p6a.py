#!/usr/bin/env python3
"""Run the approved A/B/C GPU smoke checkpoint on an idle GPU 2 or 3."""

from __future__ import annotations

from datetime import datetime, timezone
import json
import os
from pathlib import Path
import subprocess

from scripts import v2_gpu_resources
from scripts import v2_subprocess


REPO_ROOT = Path(__file__).resolve().parents[1]
STATE_PATH = REPO_ROOT / "alphagenome_custom/metadata/v2/execution_state.json"
RECORD_PATH = REPO_ROOT / "alphagenome_custom/metadata/v2/p6a_execution.json"
REQUIRED_SCOPE = "r6_gpu_auto_available_2_3"


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def require_approval() -> None:
    state = json.loads(STATE_PATH.read_text())
    gate = state.get("approvals", {}).get("G4_gpu_experiments", {})
    if not (
        state.get("current_phase") == "P6A"
        and gate.get("approved") is True
        and gate.get("scope") == REQUIRED_SCOPE
    ):
        raise RuntimeError("P6A requires current phase P6A and standing scoped G4")


def main() -> None:
    require_approval()
    resources_before = v2_gpu_resources.snapshot()
    selected = v2_gpu_resources.select_available(
        resources_before, count=1, minimum_free_mib=70_000
    )[0]
    run_root = REPO_ROOT / "runs/v2_p6a_smoke_20260714"
    log_root = REPO_ROOT / "logs/v2_p6a_smoke_20260714"
    run_root.mkdir(parents=True, exist_ok=True)
    log_root.mkdir(parents=True, exist_ok=True)
    commands = []
    for model in ("A", "B", "C"):
        output_dir = run_root / model
        commands.append(
            v2_subprocess.module_command(
                "train_v2_model",
                "--model",
                model,
                "--fold",
                "1",
                "--seed",
                "20260714",
                "--loss",
                "paper",
                "--max-steps",
                "2",
                "--sequence-length",
                "131072",
                "--output-dir",
                str(output_dir.relative_to(REPO_ROOT)),
                "--device",
                "cuda",
            )
        )
    record = {
        "schema_version": 1,
        "phase": "P6A",
        "status": "running",
        "started_at": utc_now(),
        "completed_at": None,
        "selected_physical_gpu": selected,
        "resources_before": resources_before,
        "commands": commands,
        "runs": [],
    }
    RECORD_PATH.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
    environment = os.environ.copy()
    environment["CUDA_VISIBLE_DEVICES"] = str(selected)
    try:
        for model, command in zip(("A", "B", "C"), commands, strict=True):
            log_path = log_root / f"{model}.log"
            print("command\t" + " ".join(command), flush=True)
            with log_path.open("w") as log:
                subprocess.run(
                    command,
                    cwd=REPO_ROOT,
                    env=environment,
                    stdout=log,
                    stderr=subprocess.STDOUT,
                    check=True,
                )
            run_path = run_root / model / "run.json"
            run = json.loads(run_path.read_text())
            record["runs"].append(
                {
                    "model": model,
                    "run_path": str(run_path.relative_to(REPO_ROOT)),
                    "log_path": str(log_path.relative_to(REPO_ROOT)),
                    "checkpoint_path": run["checkpoint_path"],
                    "checkpoint_sha256": run["checkpoint_sha256"],
                }
            )
            RECORD_PATH.write_text(
                json.dumps(record, indent=2, sort_keys=True) + "\n"
            )
    except Exception:
        record["status"] = "failed"
        raise
    else:
        record["status"] = "completed"
    finally:
        record["completed_at"] = utc_now()
        record["resources_after"] = v2_gpu_resources.snapshot()
        RECORD_PATH.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
