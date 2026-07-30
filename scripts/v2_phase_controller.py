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
SCHEMA_VERSION = 2
PHASES = (
    "P0",
    "P1",
    "P2",
    "P3A",
    "P3B",
    "P4",
    "P5",
    "P6A",
    "P6B",
    "P6C",
    "P7",
    "P8",
)
GATE_KEYS = {
    "G1": "G1_large_source_download_or_realignment",
    "G2": "G2_v2_manifest_scientific_review",
    "G3": "G3_dataset_or_cache_generation",
    "G4": "G4_gpu_experiments",
    "G5": "G5_final_test",
}
APPROVAL_SCOPES = {
    "G1": {
        "p3a_five_run_pilot",
        "p3b_full_reprocessing",
        "p3_full_rna_streaming_482",
    },
    "G2": {"v2_manifest_candidate_hierarchy"},
    "G3": {
        "p3a_five_run_pilot_outputs",
        "p3b_full_normalized_outputs",
        "p4_dynamic_loader_split_pilot",
        "p3_full_outputs_and_p4_loader",
        "p4_six_chromosome_block_split",
    },
    "G4": {
        "r6_gpu_experiment_matrix",
        "r6_gpu_auto_available_2_3",
        "p7_single_legacy_residual_fold1",
        "p8_b_no_lora_fold1",
    },
    "G5": {
        "r6c_single_chr_x_test",
        "r6c_single_six_chromosome_block_test",
    },
}
APPROVAL_SCOPE_PHASES = {
    "p3a_five_run_pilot": {"P3A"},
    "v2_manifest_candidate_hierarchy": {"P3A"},
    "p3a_five_run_pilot_outputs": {"P3A"},
    "p3b_full_reprocessing": {"P3B"},
    "p3b_full_normalized_outputs": {"P3B"},
    "p4_dynamic_loader_split_pilot": {"P4"},
    "r6_gpu_experiment_matrix": {"P6A"},
    "p3_full_rna_streaming_482": {"P3A"},
    "p3_full_outputs_and_p4_loader": {"P3A"},
    # The approval is recorded while the superseded P6C state is still active,
    # then remains bound to P4 after the controlled archival transition.
    "p4_six_chromosome_block_split": {"P4", "P6C"},
    "r6_gpu_auto_available_2_3": {"P3A"},
    "r6c_single_chr_x_test": {"P6C"},
    "r6c_single_six_chromosome_block_test": {"P6C"},
    "p7_single_legacy_residual_fold1": {"P7"},
    "p8_b_no_lora_fold1": {"P8"},
}
PHASE_REQUIRED_APPROVALS = {
    "P3A": (
        ("G1", "p3_full_rna_streaming_482"),
        ("G2", "v2_manifest_candidate_hierarchy"),
        ("G3", "p3_full_outputs_and_p4_loader"),
    ),
    "P3B": (
        ("G1", "p3_full_rna_streaming_482"),
        ("G2", "v2_manifest_candidate_hierarchy"),
        ("G3", "p3_full_outputs_and_p4_loader"),
    ),
    "P4": (("G3", "p4_six_chromosome_block_split"),),
    "P6A": (("G4", "r6_gpu_auto_available_2_3"),),
    "P6B": (("G4", "r6_gpu_auto_available_2_3"),),
    "P6C": (
        ("G4", "r6_gpu_auto_available_2_3"),
        ("G5", "r6c_single_six_chromosome_block_test"),
    ),
    "P7": (("G4", "p7_single_legacy_residual_fold1"),),
    "P8": (("G4", "p8_b_no_lora_fold1"),),
}
def module_command(module: str, *arguments: str) -> list[str]:
    return [sys.executable, "-m", f"scripts.{module}", *arguments]


PHASE_COMMANDS = {
    "P0": module_command("freeze_legacy_v1"),
    "P1": module_command("build_rna_seq_samples_v2"),
    "P2": module_command("build_rna_seq_groups_v2"),
    "P3A": module_command("run_v2_p3a"),
    "P3B": module_command("run_v2_p3b"),
    "P4": module_command("run_v2_p4"),
    "P5": module_command("run_v2_p5"),
    "P6A": module_command("run_v2_p6a"),
    "P6B": module_command("run_v2_p6b_six_chromosome"),
    "P6C": module_command("run_v2_p6c"),
    "P7": module_command("run_v2_p7_legacy_residual"),
    "P8": module_command("run_v2_p8_b_no_lora"),
}
REVIEW_COMMANDS = {
    phase: module_command("review_v2_phase", "--phase", phase)
    for phase in PHASES
}
SIX_CHROMOSOME_MIGRATION_COMMAND = module_command(
    "prepare_v2_six_chromosome_revision"
)


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def initial_state() -> dict[str, Any]:
    return {
        "schema_version": SCHEMA_VERSION,
        "workflow": "c_elegans_rna_seq_v2",
        "current_phase": "P0",
        "status": "PENDING",
        "updated_at": utc_now(),
        "approvals": {
            key: {
                "approved": False,
                "scope": None,
                "note": None,
                "updated_at": None,
            }
            for key in GATE_KEYS.values()
        },
        "history": [],
    }


def migrate_state(state: dict[str, Any]) -> bool:
    version = int(state.get("schema_version", 1))
    if version > SCHEMA_VERSION:
        raise RuntimeError(
            f"State schema {version} is newer than supported schema {SCHEMA_VERSION}"
        )
    changed = False
    if version < 2:
        migrated_approvals = {}
        for key in GATE_KEYS.values():
            legacy_value = state.get("approvals", {}).get(key, False)
            approved = bool(legacy_value) if isinstance(legacy_value, bool) else False
            migrated_approvals[key] = {
                "approved": approved,
                "scope": "legacy_unscoped" if approved else None,
                "note": "Migrated from schema v1" if approved else None,
                "updated_at": None,
            }
        state["approvals"] = migrated_approvals
        if state.get("current_phase") == "P3":
            state["current_phase"] = "P3A"
            append_history(
                state,
                "P3A",
                "MIGRATED",
                "Schema v1 P3 migrated to the bounded five-run P3A pilot.",
            )
        state["schema_version"] = 2
        changed = True
    return changed


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
    state = json.loads(STATE_PATH.read_text())
    if migrate_state(state):
        save_state(state)
    return state


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


def approval_satisfies(state: dict[str, Any], gate: str, scope: str) -> bool:
    record = state["approvals"].get(GATE_KEYS[gate], {})
    return record.get("approved") is True and record.get("scope") == scope


def missing_approvals(state: dict[str, Any], phase: str) -> list[str]:
    return [
        f"{gate}:{scope}"
        for gate, scope in PHASE_REQUIRED_APPROVALS.get(phase, ())
        if not approval_satisfies(state, gate, scope)
    ]


def validate_approval_scope(gate: str, scope: str, phase: str) -> None:
    if scope not in APPROVAL_SCOPES[gate]:
        allowed = ",".join(sorted(APPROVAL_SCOPES[gate]))
        raise ValueError(f"Invalid {gate} scope {scope!r}; allowed: {allowed}")
    if phase not in APPROVAL_SCOPE_PHASES[scope]:
        allowed_phases = ",".join(sorted(APPROVAL_SCOPE_PHASES[scope]))
        raise ValueError(
            f"Scope {scope!r} may be approved only in {allowed_phases}; "
            f"current phase is {phase}"
        )


def run_phase(state: dict[str, Any], phase: str) -> int:
    if phase != state["current_phase"]:
        raise SystemExit(
            f"Requested {phase}, but current phase is {state['current_phase']}"
        )
    if state["status"] in {"FAIL", "RUNNING", "REVIEWING", "COMPLETE"}:
        raise SystemExit(
            f"Cannot run {phase} while status={state['status']}; use reopen after FAIL"
        )
    missing = missing_approvals(state, phase)
    if missing:
        state["status"] = "APPROVAL_REQUIRED"
        append_history(
            state,
            phase,
            "APPROVAL_REQUIRED",
            f"Missing scoped approvals: {','.join(missing)}.",
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


def prepare_six_chromosome_revision(state: dict[str, Any]) -> int:
    """Archive the old holdout lineage and reopen the controlled P4 revision."""
    scope = "p4_six_chromosome_block_split"
    if state.get("current_phase") != "P6C":
        raise SystemExit(
            "Six-chromosome revision preparation is only available from P6C"
        )
    if state.get("approvals", {}).get(GATE_KEYS["G5"], {}).get("approved") is True:
        raise SystemExit(
            "Cannot revise the split after a final-test approval has been recorded"
        )
    if not approval_satisfies(state, "G3", scope):
        state["status"] = "APPROVAL_REQUIRED"
        append_history(
            state,
            "P6C",
            "APPROVAL_REQUIRED",
            f"Missing scoped approval for split revision: G3:{scope}.",
        )
        save_state(state)
        return 2
    if state.get("status") in {"RUNNING", "REVIEWING", "COMPLETE"}:
        raise SystemExit(
            f"Cannot prepare six-chromosome revision while status={state['status']}"
        )

    state["status"] = "RUNNING"
    append_history(
        state,
        "P6C",
        "SPLIT_REVISION_RUNNING",
        "Archiving chromosome-holdout evidence before six-chromosome P4 reopen.",
    )
    save_state(state)
    code = run_command(SIX_CHROMOSOME_MIGRATION_COMMAND)
    if code != 0:
        state["status"] = "FAIL"
        append_history(
            state,
            "P6C",
            "SPLIT_REVISION_FAIL",
            f"Revision preparation exited with code {code}.",
        )
        save_state(state)
        return 1

    append_history(
        state,
        "P6C",
        "SPLIT_REVISION_PREPARED",
        "Chromosome-holdout evidence archived and final lock superseded.",
    )
    state["current_phase"] = "P4"
    state["status"] = "PENDING"
    append_history(
        state,
        "P4",
        "REOPENED",
        "Six-chromosome split revision prepared; regenerate manifests and means.",
    )
    save_state(state)
    return 0


def start_p7_legacy_residual(state: dict[str, Any]) -> int:
    """Open one bounded post-completion development-only architecture check."""

    if state.get("current_phase") != "P6C" or state.get("status") != "COMPLETE":
        raise SystemExit("P7 may start only after the completed P6C workflow")
    lock_path = REPO_ROOT / "alphagenome_custom/metadata/v2/final_test_lock.json"
    report_path = REPO_ROOT / "alphagenome_custom/metadata/v2/final_test_report.json"
    if not lock_path.is_file() or not report_path.is_file():
        raise SystemExit("P7 requires the completed P6C lock and report")
    lock = json.loads(lock_path.read_text())
    if not (
        lock.get("test_consumed") is True
        and lock.get("test_status") == "completed"
    ):
        raise SystemExit("P7 requires the already-consumed completed P6C test lock")
    state["current_phase"] = "P7"
    state["status"] = "PENDING"
    append_history(
        state,
        "P7",
        "PENDING",
        "User requested one fold-1 legacy 1bp-B architecture port on all 241 v2 tracks; development validation only and no final-test access.",
    )
    save_state(state)
    return 0


def start_p8_b_no_lora(state: dict[str, Any]) -> int:
    """Open one bounded post-completion B-noLoRA development ablation."""

    if state.get("current_phase") != "P7" or state.get("status") != "COMPLETE":
        raise SystemExit("P8 may start only after the completed P7 architecture check")
    lock_path = REPO_ROOT / "alphagenome_custom/metadata/v2/final_test_lock.json"
    report_path = REPO_ROOT / "alphagenome_custom/metadata/v2/final_test_report.json"
    if not lock_path.is_file() or not report_path.is_file():
        raise SystemExit("P8 requires the completed P6C lock and report")
    lock = json.loads(lock_path.read_text())
    if not (
        lock.get("test_consumed") is True
        and lock.get("test_status") == "completed"
    ):
        raise SystemExit("P8 requires the already-consumed completed P6C test lock")
    state["current_phase"] = "P8"
    state["status"] = "PENDING"
    append_history(
        state,
        "P8",
        "PENDING",
        "User requested one fold-1 B-noLoRA ablation retaining the C. elegans embedding and dual RNA head; development validation only and no final-test access.",
    )
    save_state(state)
    return 0


def run_until_boundary(state: dict[str, Any], first_phase: str) -> int:
    phase = first_phase
    while True:
        code = run_phase(state, phase)
        if code != 0 or state["status"] == "COMPLETE":
            return code
        phase = state["current_phase"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command", required=True)
    subparsers.add_parser("status")
    run_parser = subparsers.add_parser("run")
    run_parser.add_argument("--phase", choices=PHASES)
    run_parser.add_argument(
        "--auto",
        action="store_true",
        help="Continue through PASS phases until FAIL, an approval boundary, or completion.",
    )
    reopen_parser = subparsers.add_parser("reopen")
    reopen_parser.add_argument("--phase", required=True, choices=PHASES)
    reopen_parser.add_argument("--reason", required=True)
    subparsers.add_parser("prepare-six-chromosome-revision")
    subparsers.add_parser("start-p7-legacy-residual")
    subparsers.add_parser("start-p8-b-no-lora")
    approve_parser = subparsers.add_parser("approve")
    approve_parser.add_argument("--gate", required=True, choices=sorted(GATE_KEYS))
    approve_parser.add_argument("--scope", required=True)
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
    if args.command == "prepare-six-chromosome-revision":
        raise SystemExit(prepare_six_chromosome_revision(state))
    if args.command == "start-p7-legacy-residual":
        raise SystemExit(start_p7_legacy_residual(state))
    if args.command == "start-p8-b-no-lora":
        raise SystemExit(start_p8_b_no_lora(state))
    if args.command in {"approve", "revoke"}:
        approved = args.command == "approve"
        if approved:
            try:
                validate_approval_scope(args.gate, args.scope, state["current_phase"])
            except ValueError as error:
                raise SystemExit(str(error)) from error
        state["approvals"][GATE_KEYS[args.gate]] = {
            "approved": approved,
            "scope": args.scope if approved else None,
            "note": args.note,
            "updated_at": utc_now(),
        }
        append_history(
            state,
            state["current_phase"],
            f"{args.gate}_{'APPROVED' if approved else 'REVOKED'}",
            f"scope={args.scope if approved else 'all'}; {args.note}",
        )
        if approved and state["status"] == "APPROVAL_REQUIRED":
            state["status"] = (
                "PENDING"
                if not missing_approvals(state, state["current_phase"])
                else "APPROVAL_REQUIRED"
            )
        save_state(state)
        print(json.dumps(state, indent=2, sort_keys=True))
        return
    phase = args.phase or state["current_phase"]
    if args.auto:
        raise SystemExit(run_until_boundary(state, phase))
    raise SystemExit(run_phase(state, phase))


if __name__ == "__main__":
    main()
