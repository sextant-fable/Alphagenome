#!/usr/bin/env python3
"""Run P15's metadata-only I--V training/validation contract construction."""

from __future__ import annotations

import hashlib
import json
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
METADATA_ROOT = REPO_ROOT / "alphagenome_custom" / "metadata" / "v2"
CONTRACT_DIR = REPO_ROOT / "results" / "v2_p15_iv_training_interval_contract"
CONTRACT_PATH = CONTRACT_DIR / "p15_iv_training_interval_contract.json"
EXECUTION_PATH = METADATA_ROOT / "p15_iv_training_contract_execution.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def main() -> None:
    if EXECUTION_PATH.exists():
        raise FileExistsError(f"P15 execution record already exists: {EXECUTION_PATH}")
    if CONTRACT_DIR.exists():
        raise FileExistsError(f"P15 output directory already exists: {CONTRACT_DIR}")

    python = str(Path(sys.executable).resolve())
    command = [
        python,
        "-m",
        "scripts.prepare_v2_p15_iv_training_intervals",
        "--output-dir",
        str(CONTRACT_DIR.relative_to(REPO_ROOT)),
    ]
    print("command\t" + " ".join(command), flush=True)
    subprocess.run(command, cwd=REPO_ROOT, check=True)
    if not CONTRACT_PATH.is_file():
        raise RuntimeError(f"P15 contract was not created: {CONTRACT_PATH}")
    contract = json.loads(CONTRACT_PATH.read_text())
    if contract.get("allowed_chromosomes") != ["I", "II", "III", "IV", "V"]:
        raise RuntimeError("P15 contract has an unexpected chromosome set")

    execution = {
        "schema_version": 1,
        "phase": "P15",
        "status": "completed",
        "started_and_completed_at": utc_now(),
        "git_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=REPO_ROOT, text=True).strip(),
        "working_directory": str(REPO_ROOT),
        "python_executable": python,
        "controller_invocation": "python -m scripts.v2_phase_controller run --phase P15",
        "registered_phase_command": "python -m scripts.run_v2_p15_iv_training_contract",
        "commands": [command],
        "resource_contract": "metadata-only interval filtering; no GPU, DNA, BigWig, checkpoint, training, chromosome-X or locked-test access",
        "contract": str(CONTRACT_PATH.relative_to(REPO_ROOT)),
        "contract_sha256": sha256(CONTRACT_PATH),
        "output_sha256": {
            str(path.relative_to(REPO_ROOT)): sha256(path)
            for fold in range(1, 6)
            for path in (
                CONTRACT_DIR / "intervals" / f"fold_{fold}" / "train.tsv",
                CONTRACT_DIR / "intervals" / f"fold_{fold}" / "valid.tsv",
            )
        },
    }
    EXECUTION_PATH.write_text(json.dumps(execution, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"status": "completed", "execution": str(EXECUTION_PATH.relative_to(REPO_ROOT))}, indent=2))


if __name__ == "__main__":
    main()
