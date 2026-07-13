#!/usr/bin/env python3
"""Fail-closed controller for the C. elegans RNA-seq v2 workflow."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import subprocess
import sys
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[1]
STATE_PATH = REPO_ROOT / "alphagenome_custom/metadata/v2/execution_state.json"
PHASES = ("P0", "P1", "P2", "P3", "P4", "P5", "P6")
GATE_KEYS = {
    "G1": "G1_large_source_download_or_realignment",
    "G2": "G2_v2_manifest_scientific_review",
    "G3": "G3_dataset_or_cache_generation",
    "G4": "G4_gpu_experiments",
    "G5": "G5_final_test",
}
PHASE_REQUIRED_APPROVALS = {
    "P3": ("G1", "G2", "G3"),
    "P5": ("G4",),
    "P6": ("G4", "G5"),
}
PHASE_COMMANDS = {
    "P0": [sys.executable, "scripts/freeze_legacy_v1.py"],
    "P1": [sys.executable, "scripts/build_rna_seq_samples_v2.py"],
    "P2": [sys.executable, "scripts/build_rna_seq_groups_v2.py"],
}
REVIEW_COMMANDS = {
    "P0": [sys.executable, "scripts/review_v2_phase.py", "--phase", "P0"],
    "P1": [sys.executable, "scripts/review_v2_phase.py", "--phase", "P1"],
    "P2": [sys.executable, "scripts/review_v2_phase.py", "--phase", "P2"],
}


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def initial_state() -> dict[str, Any]:
    return {
        "schema_version": 1,
        "workflow": "c_elegans_rna_seq_v2",
        "current_phase": "P0",
        "status": "PENDING",
        "updated_at": utc_now(),
        "approvals": {
            "G1_large_source_download_or_realignment": False,
            "G2_v2_manifest_scientific_review": False,
            "G3_dataset_or_cache_generation": False,
            "G4_gpu_experiments": False,
            "G5_final_test": False,
        },
        "history": [],
    }


def atomic_write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = path.with_suffix(path.suffix + ".tmp")
    tmp_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    tmp_path.replace(path)


def load_state() -> dict[str, Any]:
    if not STATE_PATH.exists():
        state = initial_state()
        atomic_write_json(STATE_PATH, state)
        return state
    return json.loads(STATE_PATH.read_text())


def save_state(state: dict[str, Any]) -> None:
    state["updated_at"] = utc_now()
    atomic_write_json(STATE_PATH, state)


def append_history(
    state: dict[str, Any], phase: str, status: str, detail: str
) -> None:
    state["history"].append(
        {
            "phase": phase,
            "status": status,
            "detail": detail,
            "timestamp": utc_now(),
        }
    )


def run_command(command: list[str]) -> int:
    print("command\t" + " ".join(command), flush=True)
    completed = subprocess.run(command, cwd=REPO_ROOT, check=False)
    return completed.returncode


def next_phase(phase: str) -> str | None:
    index = PHASES.index(phase)
    return PHASES[index + 1] if index + 1 < len(PHASES) else None


def run_phase(state: dict[str, Any], phase: str) -> int:
    if phase != state["current_phase"]:
        raise SystemExit(
            f"Requested {phase}, but current phase is {state['current_phase']}"
        )
    missing_approvals = [
        gate
        for gate in PHASE_REQUIRED_APPROVALS.get(phase, ())
        if not state["approvals"][GATE_KEYS[gate]]
    ]
    if missing_approvals:
        state["status"] = "APPROVAL_REQUIRED"
        append_history(
            state,
            phase,
            "APPROVAL_REQUIRED",
            f"Missing approvals: {','.join(missing_approvals)}.",
        )
        save_state(state)
        return 2
    if phase not in PHASE_COMMANDS or phase not in REVIEW_COMMANDS:
        state["status"] = "APPROVAL_REQUIRED"
        append_history(
            state,
            phase,
            "APPROVAL_REQUIRED",
            "Phase implementation has not been registered in the controller.",
        )
        save_state(state)
        return 2

    state["status"] = "RUNNING"
    append_history(state, phase, "RUNNING", "Phase command started.")
    save_state(state)
    phase_code = run_command(PHASE_COMMANDS[phase])
    if phase_code != 0:
        state["status"] = "FAIL"
        append_history(
            state, phase, "FAIL", f"Phase command exited with code {phase_code}."
        )
        save_state(state)
        return 1

    state["status"] = "REVIEWING"
    append_history(state, phase, "REVIEWING", "Phase review started.")
    save_state(state)
    review_code = run_command(REVIEW_COMMANDS[phase])
    if review_code != 0:
        state["status"] = "FAIL" if review_code == 1 else "APPROVAL_REQUIRED"
        append_history(
            state,
            phase,
            state["status"],
            f"Phase review exited with code {review_code}.",
        )
        save_state(state)
        return review_code

    append_history(state, phase, "PASS", "All phase review checks passed.")
    following = next_phase(phase)
    if following is None:
        state["status"] = "COMPLETE"
    else:
        state["current_phase"] = following
        state["status"] = "PENDING"
        append_history(
            state, following, "PENDING", f"Advanced automatically after {phase} PASS."
        )
    save_state(state)
    return 0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command", required=True)
    subparsers.add_parser("status")
    run_parser = subparsers.add_parser("run")
    run_parser.add_argument("--phase", choices=PHASES)
    reopen_parser = subparsers.add_parser("reopen")
    reopen_parser.add_argument("--phase", required=True, choices=PHASES)
    reopen_parser.add_argument("--reason", required=True)
    approve_parser = subparsers.add_parser("approve")
    approve_parser.add_argument("--gate", required=True, choices=sorted(GATE_KEYS))
    approve_parser.add_argument("--note", required=True)
    revoke_parser = subparsers.add_parser("revoke")
    revoke_parser.add_argument("--gate", required=True, choices=sorted(GATE_KEYS))
    revoke_parser.add_argument("--note", required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    state = load_state()
    if args.command == "status":
        print(json.dumps(state, indent=2, sort_keys=True))
        return
    if args.command == "reopen":
        target_index = PHASES.index(args.phase)
        current_index = PHASES.index(state["current_phase"])
        if target_index > current_index:
            raise SystemExit("Cannot reopen a phase that has not been reached")
        append_history(state, args.phase, "REOPENED", args.reason)
        state["current_phase"] = args.phase
        state["status"] = "PENDING"
        save_state(state)
        print(json.dumps(state, indent=2, sort_keys=True))
        return
    if args.command in {"approve", "revoke"}:
        approved = args.command == "approve"
        state["approvals"][GATE_KEYS[args.gate]] = approved
        append_history(
            state,
            state["current_phase"],
            f"{args.gate}_{'APPROVED' if approved else 'REVOKED'}",
            args.note,
        )
        if approved and state["status"] == "APPROVAL_REQUIRED":
            state["status"] = "PENDING"
        save_state(state)
        print(json.dumps(state, indent=2, sort_keys=True))
        return
    phase = args.phase or state["current_phase"]
    raise SystemExit(run_phase(state, phase))


if __name__ == "__main__":
    main()
