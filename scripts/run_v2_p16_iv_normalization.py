#!/usr/bin/env python3
"""Run P16's preregistered I--V-only normalization calculation."""

from __future__ import annotations

import csv
import hashlib
import json
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
METADATA_ROOT = REPO_ROOT / "alphagenome_custom" / "metadata" / "v2"
SPEC_PATH = METADATA_ROOT / "p16_iv_training_normalization_spec.json"
EXECUTION_PATH = METADATA_ROOT / "p16_iv_training_normalization_execution.json"
STATE_PATH = METADATA_ROOT / "execution_state.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def require_controller_scope() -> None:
    state = json.loads(STATE_PATH.read_text())
    approval = state.get("approvals", {}).get("G3_dataset_or_cache_generation", {})
    if not (
        state.get("current_phase") == "P16"
        and state.get("status") == "RUNNING"
        and approval.get("approved") is True
        and approval.get("scope") == "p16_iv_training_normalization"
    ):
        raise RuntimeError("P16 requires controller state P16/RUNNING and G3:p16_iv_training_normalization")


def validate_outputs(spec: dict[str, object]) -> None:
    outputs = spec["outputs"]
    means_path = REPO_ROOT / str(outputs["means"])
    summary_path = REPO_ROOT / str(outputs["summary"])
    if not means_path.is_file() or not summary_path.is_file():
        raise RuntimeError("P16 normalization outputs are incomplete")
    rows = read_tsv(means_path)
    required_columns = {"group_id", *(f"fold_{fold}_train_nonzero_mean" for fold in range(1, 6))}
    if len(rows) != 241 or not rows or not required_columns.issubset(rows[0]):
        raise RuntimeError("P16 normalization means do not cover the 241 training heads")
    if "development_train_nonzero_mean" in rows[0]:
        raise RuntimeError("P16 strict I--V normalization must not retain a development-wide mean")
    summary = json.loads(summary_path.read_text())
    if not (
        summary.get("contract") == "strict_interval_root"
        and summary.get("allowed_chromosomes") == ["I", "II", "III", "IV", "V"]
        and summary.get("locked_test_block_signal_reads") == 0
    ):
        raise RuntimeError("P16 normalization summary violates its I--V contract")


def main() -> None:
    if EXECUTION_PATH.exists():
        raise FileExistsError(f"P16 execution record already exists: {EXECUTION_PATH}")
    if not SPEC_PATH.is_file():
        raise FileNotFoundError("P16 requires an immutable preregistered specification")
    require_controller_scope()
    spec = json.loads(SPEC_PATH.read_text())
    if spec.get("phase") != "P16" or spec.get("status") != "preregistered":
        raise RuntimeError("P16 specification is invalid")
    if spec.get("allowed_chromosomes") != ["I", "II", "III", "IV", "V"]:
        raise RuntimeError("P16 specification does not enforce I--V-only normalization")
    outputs = spec["outputs"]
    means_path = REPO_ROOT / str(outputs["means"])
    summary_path = REPO_ROOT / str(outputs["summary"])
    if means_path.exists() or summary_path.exists():
        raise FileExistsError("P16 output path already exists; preserve and review the prior attempt")
    command = [
        str(Path(sys.executable).resolve()), "-m", "scripts.compute_v2_replicate_holdout_means",
        "--track-manifest", str(spec["source"]["training_track_manifest"]),
        "--interval-root", "results/v2_p15_iv_training_interval_contract/intervals",
        "--allowed-chromosomes", "I", "II", "III", "IV", "V",
        "--output", str(outputs["means"]), "--summary", str(outputs["summary"]),
    ]
    print("command\t" + " ".join(command), flush=True)
    subprocess.run(command, cwd=REPO_ROOT, check=True)
    validate_outputs(spec)
    execution = {
        "schema_version": 1,
        "phase": "P16",
        "status": "completed",
        "started_and_completed_at": utc_now(),
        "git_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=REPO_ROOT, text=True).strip(),
        "controller_invocation": "python -m scripts.v2_phase_controller run --phase P16",
        "registered_phase_command": "python -m scripts.run_v2_p16_iv_normalization",
        "resource_contract": "CPU-only I-V training-aggregate BigWig means; no GPU, DNA, checkpoint, external source, chromosome-X or locked-test access",
        "spec": str(SPEC_PATH.relative_to(REPO_ROOT)),
        "spec_sha256": sha256(SPEC_PATH),
        "means": str(means_path.relative_to(REPO_ROOT)),
        "means_sha256": sha256(means_path),
        "summary": str(summary_path.relative_to(REPO_ROOT)),
        "summary_sha256": sha256(summary_path),
    }
    EXECUTION_PATH.write_text(json.dumps(execution, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"status": "completed", "execution": str(EXECUTION_PATH.relative_to(REPO_ROOT))}, indent=2))


if __name__ == "__main__":
    main()
