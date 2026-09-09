#!/usr/bin/env python3
"""Audit the official Borzoi source and adapter prerequisites.

This phase deliberately stops before C. elegans training. A successful audit is
not a completed Borzoi adapter result; it authorizes a separate adapter run
only when the task interface and checkpoint terms are documented.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import socket
import subprocess
from datetime import datetime, timezone


REPO_ROOT = Path(__file__).resolve().parents[1]
METADATA = REPO_ROOT / "alphagenome_custom/metadata/v2"
SPEC_PATH = METADATA / "p24_borzoi_adapter_spec.json"
EXECUTION_PATH = METADATA / "p24_borzoi_adapter_execution.json"
VENDOR_PATH = Path("/home/zelinli6/vendor/borzoi")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def require_scope() -> None:
    state = json.loads((METADATA / "execution_state.json").read_text())
    gate = state.get("approvals", {}).get("G4_gpu_experiments", {})
    if not (state.get("current_phase") == "P24" and state.get("status") == "RUNNING" and gate.get("approved") is True and gate.get("scope") == "p24_borzoi_adapter"):
        raise RuntimeError("P24 requires controller P24/RUNNING and G4:p24_borzoi_adapter")


def main() -> None:
    require_scope()
    spec = json.loads(SPEC_PATH.read_text())
    if not VENDOR_PATH.is_dir():
        raise FileNotFoundError(VENDOR_PATH)
    commit = subprocess.check_output(["git", "-C", str(VENDOR_PATH), "rev-parse", "HEAD"], text=True).strip()
    license_path = VENDOR_PATH / "LICENSE"
    readme = (VENDOR_PATH / "README.md").read_text()
    pyproject = (VENDOR_PATH / "pyproject.toml").read_text()
    required_urls = {
        "human_weight_url": spec["human_weight_url"] in readme,
        "mouse_or_human_model_availability": "model0_best.h5" in readme,
        "apache_license_text": "Apache License" in license_path.read_text(),
        "tensorflow_dependency": "tensorflow~=2.15.0" in pyproject,
        "python_3_10_documented": "Python == 3.10" in readme,
        "native_context_documented": "524kb input" in readme.lower(),
        "native_32bp_output_documented": "32bp resolution" in readme.lower(),
    }
    checks = {
        "official_commit": {"expected": spec["official_commit"], "observed": commit, "pass": commit == spec["official_commit"]},
        "license_sha256": {"observed": sha256(license_path), "pass": license_path.is_file()},
        "required_source_metadata": {"observed": required_urls, "pass": all(required_urls.values())},
        "c_elegans_adapter_defined": {"observed": False, "pass": False, "note": "No C. elegans 241-head adapter has been trained or registered yet."},
        "native_checkpoint_fairness": {"observed": False, "pass": False, "note": "Native human/mouse output heads are not task-faithful for WBcel235."},
    }
    execution = {
        "schema_version": 1,
        "phase": "P24",
        "status": "completed_feasibility_audit",
        "audited_at": utc_now(),
        "hostname": socket.gethostname(),
        "working_directory": str(REPO_ROOT),
        "controller_invocation": "python -m scripts.v2_phase_controller run --phase P24",
        "registered_phase_command": "python -m scripts.run_v2_p24_borzoi_adapter",
        "spec_path": str(SPEC_PATH.relative_to(REPO_ROOT)),
        "spec_sha256": sha256(SPEC_PATH),
        "vendor_path": str(VENDOR_PATH),
        "vendor_commit": commit,
        "license": "Apache-2.0",
        "native_model_boundary": "human/mouse 524kb input and 32bp output; not directly comparable to WBcel235 241-head task",
        "adapter_status": "not_yet_run",
        "final_test_access": "prohibited",
        "chromosome_x_reads": 0,
        "locked_test_block_signal_reads": 0,
        "checks": checks,
        "next_gate": "define and register a task-faithful C. elegans adapter before any Borzoi training",
    }
    EXECUTION_PATH.write_text(json.dumps(execution, indent=2, sort_keys=True) + "\n")
    print(json.dumps(execution, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
