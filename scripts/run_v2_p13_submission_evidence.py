#!/usr/bin/env python3
"""Run the bounded P13 I--V replicate-agreement evidence task.

P13 deliberately creates a new I--V-only interval contract before opening a
BigWig. It is CPU-only, reads no checkpoint, and cannot access chromosome X or
the locked-test block. The phase controller owns its execution and review.
"""

from __future__ import annotations

import hashlib
import json
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
METADATA_ROOT = REPO_ROOT / "alphagenome_custom" / "metadata" / "v2"
CONTRACT_DIR = REPO_ROOT / "results" / "v2_p13_iv_interval_contract"
ANALYSIS_DIR = REPO_ROOT / "results" / "v2_p13_replicate_agreement_iv"
EXECUTION_PATH = METADATA_ROOT / "p13_submission_evidence_execution.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def run(command: list[str]) -> None:
    print("command\t" + " ".join(command), flush=True)
    subprocess.run(command, cwd=REPO_ROOT, check=True)


def git_commit() -> str:
    return subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=REPO_ROOT, text=True).strip()


def main() -> None:
    if EXECUTION_PATH.exists():
        raise FileExistsError(f"P13 execution record already exists: {EXECUTION_PATH}")

    python = str(Path(sys.executable).resolve())
    prepare_command = [
        python,
        "-m",
        "scripts.prepare_v2_p13_iv_intervals",
        "--output-dir",
        str(CONTRACT_DIR.relative_to(REPO_ROOT)),
    ]
    analysis_command = [
        python,
        "-m",
        "scripts.analyze_v2_p13_replicate_agreement",
        "--interval-root",
        str((CONTRACT_DIR / "intervals").relative_to(REPO_ROOT)),
        "--output-dir",
        str(ANALYSIS_DIR.relative_to(REPO_ROOT)),
        "--max-io-workers",
        "4",
    ]
    contract = CONTRACT_DIR / "p13_iv_interval_contract.json"
    analysis_audit = ANALYSIS_DIR / "p13_replicate_agreement_audit.json"
    required = [
        contract,
        analysis_audit,
        ANALYSIS_DIR / "p13_replicate_agreement_per_group.tsv",
        ANALYSIS_DIR / "p13_replicate_agreement_by_fold.tsv",
        ANALYSIS_DIR / "p13_replicate_agreement_summary.tsv",
    ]
    if CONTRACT_DIR.exists():
        if not contract.is_file():
            raise RuntimeError(f"P13 interval contract directory is incomplete: {CONTRACT_DIR}")
    else:
        run(prepare_command)
    analysis_complete = all(path.is_file() for path in required[1:])
    if ANALYSIS_DIR.exists() and not analysis_complete:
        raise RuntimeError(f"P13 analysis directory is incomplete: {ANALYSIS_DIR}")
    if not ANALYSIS_DIR.exists():
        run(analysis_command)

    missing = [str(path.relative_to(REPO_ROOT)) for path in required if not path.is_file()]
    if missing:
        raise RuntimeError(f"P13 required artifacts missing after analysis: {missing}")
    contract_payload = json.loads(contract.read_text())
    if contract_payload.get("allowed_chromosomes") != ["I", "II", "III", "IV", "V"]:
        raise RuntimeError("P13 interval contract has an unexpected chromosome set")

    execution = {
        "schema_version": 1,
        "phase": "P13",
        "status": "completed",
        "started_and_completed_at": utc_now(),
        "git_commit": git_commit(),
        "working_directory": str(REPO_ROOT),
        "python_executable": python,
        "controller_invocation": "python -m scripts.v2_phase_controller run --phase P13",
        "registered_phase_command": "python -m scripts.run_v2_p13_submission_evidence",
        "commands": [prepare_command, analysis_command],
        "resume_policy": (
            "The registered execution may reuse a complete immutable interval contract and complete analysis artifacts "
            "only after a post-processing/provenance failure; no BigWig reread or metric recomputation occurs during resume."
        ),
        "resource_contract": "CPU-only BigWig reads; no GPU, checkpoint, training, chromosome-X or locked-test access",
        "contract": str(contract.relative_to(REPO_ROOT)),
        "contract_sha256": sha256(contract),
        "analysis_audit": str(analysis_audit.relative_to(REPO_ROOT)),
        "analysis_audit_sha256": sha256(analysis_audit),
        "output_sha256": {str(path.relative_to(REPO_ROOT)): sha256(path) for path in required},
    }
    EXECUTION_PATH.write_text(json.dumps(execution, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"status": "completed", "execution": str(EXECUTION_PATH.relative_to(REPO_ROOT))}, indent=2))


if __name__ == "__main__":
    main()
