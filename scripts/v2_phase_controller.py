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
    "P9",
    "P10",
    "P11",
    "P12",
)
# P13 is an explicitly requested post-completion manuscript-evidence phase.
# Keeping it outside PHASES preserves P12 as the terminal state of the
# original controlled v2 workflow.
POST_COMPLETION_PHASES = ("P13", "P14", "P15", "P16", "P17", "P18", "P19", "P20", "P21", "P22", "P23", "P24", "P25")
KNOWN_PHASES = PHASES + POST_COMPLETION_PHASES
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
        # The P18 external-data sidecar is bound to P17 so transfer can run
        # concurrently without changing the serial GPU-matrix phase state.
        "p18_external_rna_download_and_reprocessing",
        "p18_external_replacement_download_and_reprocessing",
        "p20_eqtl_public_source_download_and_audit",
    },
    "G2": {"v2_manifest_candidate_hierarchy"},
    "G3": {
        "p3a_five_run_pilot_outputs",
        "p3b_full_normalized_outputs",
        "p4_dynamic_loader_split_pilot",
        "p3_full_outputs_and_p4_loader",
        "p4_six_chromosome_block_split",
        "p11_replicate_holdout_data",
        "p13_iv_replicate_agreement",
        "p15_iv_training_interval_contract",
        "p16_iv_training_normalization",
        "p22_eqtl_effect_association",
    },
    "G4": {
        "r6_gpu_experiment_matrix",
        "r6_gpu_auto_available_2_3",
        "p7_single_legacy_residual_fold1",
        "p8_b_no_lora_fold1",
        "p9_submission_evidence_matrix",
        "p10_dpy27_internal_application",
        "p12_replicate_holdout_matrix",
        "p14_iv_frozen_model_reference",
        "p17_iv_controlled_ablation_and_basenji2",
        "p19_external_represented_head_scoring",
        "p21_eqtl_variant_scoring",
        "p23_lora_sensitivity",
        "p24_borzoi_adapter",
        "p25_borzoi_adapter_pilot",
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
    "p18_external_rna_download_and_reprocessing": {"P17"},
    "p18_external_replacement_download_and_reprocessing": {"P18"},
    "p3_full_outputs_and_p4_loader": {"P3A"},
    # The approval is recorded while the superseded P6C state is still active,
    # then remains bound to P4 after the controlled archival transition.
    "p4_six_chromosome_block_split": {"P4", "P6C"},
    "r6_gpu_auto_available_2_3": {"P3A"},
    "r6c_single_chr_x_test": {"P6C"},
    "r6c_single_six_chromosome_block_test": {"P6C"},
    "p7_single_legacy_residual_fold1": {"P7"},
    "p8_b_no_lora_fold1": {"P8"},
    "p9_submission_evidence_matrix": {"P9"},
    "p10_dpy27_internal_application": {"P10"},
    "p11_replicate_holdout_data": {"P11"},
    "p12_replicate_holdout_matrix": {"P12"},
    "p13_iv_replicate_agreement": {"P13"},
    "p14_iv_frozen_model_reference": {"P14"},
    "p15_iv_training_interval_contract": {"P15"},
    "p16_iv_training_normalization": {"P16"},
    "p17_iv_controlled_ablation_and_basenji2": {"P17"},
    "p19_external_represented_head_scoring": {"P19"},
    "p20_eqtl_public_source_download_and_audit": {"P20"},
    "p21_eqtl_variant_scoring": {"P21"},
    "p22_eqtl_effect_association": {"P22"},
    "p23_lora_sensitivity": {"P23"},
    "p24_borzoi_adapter": {"P24"},
    "p25_borzoi_adapter_pilot": {"P25"},
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
    "P9": (("G4", "p9_submission_evidence_matrix"),),
    "P10": (("G4", "p10_dpy27_internal_application"),),
    "P11": (("G3", "p11_replicate_holdout_data"),),
    "P12": (("G4", "p12_replicate_holdout_matrix"),),
    "P13": (("G3", "p13_iv_replicate_agreement"),),
    "P14": (("G4", "p14_iv_frozen_model_reference"),),
    "P15": (("G3", "p15_iv_training_interval_contract"),),
    "P16": (("G3", "p16_iv_training_normalization"),),
    "P17": (("G4", "p17_iv_controlled_ablation_and_basenji2"),),
    "P18": (("G1", "p18_external_replacement_download_and_reprocessing"),),
    "P19": (("G4", "p19_external_represented_head_scoring"),),
    "P20": (("G1", "p20_eqtl_public_source_download_and_audit"),),
    "P21": (("G4", "p21_eqtl_variant_scoring"),),
    "P22": (("G3", "p22_eqtl_effect_association"),),
    "P23": (("G4", "p23_lora_sensitivity"),),
    "P24": (("G4", "p24_borzoi_adapter"),),
    "P25": (("G4", "p25_borzoi_adapter_pilot"),),
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
    "P9": module_command("run_v2_p9_submission_evidence"),
    "P10": module_command("run_dpy27_internal_application"),
    "P11": module_command("run_v2_p11_replicate_holdout"),
    "P12": module_command("run_v2_p12_replicate_holdout"),
    "P13": module_command("run_v2_p13_submission_evidence"),
    "P14": module_command("run_v2_p14_iv_model_reference"),
    "P15": module_command("run_v2_p15_iv_training_contract"),
    "P16": module_command("run_v2_p16_iv_normalization"),
    "P17": module_command("run_v2_p17_iv_controlled_matrix"),
    "P18": module_command("run_v2_p18_external_continuation"),
    "P19": module_command("run_v2_p19_external_scoring"),
    "P20": module_command("run_v2_p20_eqtl_source_download"),
    "P21": module_command("run_v2_p21_eqtl_variant_scoring"),
    "P22": module_command("run_v2_p22_eqtl_effect_association"),
    "P23": module_command("run_v2_p23_lora_sensitivity"),
    "P24": module_command("run_v2_p24_borzoi_adapter"),
    "P25": module_command("run_v2_p25_borzoi_adapter_pilot"),
}
REVIEW_COMMANDS = {
    phase: module_command("review_v2_phase", "--phase", phase)
    for phase in KNOWN_PHASES
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
    if phase in POST_COMPLETION_PHASES:
        return None
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


def start_p9_submission_evidence(state: dict[str, Any]) -> int:
    """Open the registered post-completion submission-evidence matrix."""

    if state.get("current_phase") != "P8" or state.get("status") != "COMPLETE":
        raise SystemExit("P9 may start only after the completed P8 ablation")
    lock_path = REPO_ROOT / "alphagenome_custom/metadata/v2/final_test_lock.json"
    report_path = REPO_ROOT / "alphagenome_custom/metadata/v2/final_test_report.json"
    if not lock_path.is_file() or not report_path.is_file():
        raise SystemExit("P9 requires the completed P6C lock and report")
    lock = json.loads(lock_path.read_text())
    if not (lock.get("test_consumed") is True and lock.get("test_status") == "completed"):
        raise SystemExit("P9 requires the preserved completed P6C test lock")
    state["current_phase"] = "P9"
    state["status"] = "PENDING"
    append_history(
        state,
        "P9",
        "PENDING",
        "User requested the 59-job submission-evidence matrix; development validation only and no final-test access.",
    )
    save_state(state)
    return 0


def start_p11_replicate_holdout(state: dict[str, Any]) -> int:
    """Open the CPU-only biological-replicate holdout data phase."""

    if state.get("current_phase") != "P10" or state.get("status") != "COMPLETE":
        raise SystemExit("P11 may start only after the completed P10 workflow")
    lock_path = REPO_ROOT / "alphagenome_custom/metadata/v2/final_test_lock.json"
    report_path = REPO_ROOT / "alphagenome_custom/metadata/v2/final_test_report.json"
    if not lock_path.is_file() or not report_path.is_file():
        raise SystemExit("P11 requires the preserved completed P6C lock and report")
    lock = json.loads(lock_path.read_text())
    if lock.get("test_consumed") is not True or lock.get("test_status") != "completed":
        raise SystemExit("P11 requires the already-consumed completed P6C test lock")
    state["current_phase"] = "P11"
    state["status"] = "PENDING"
    append_history(
        state,
        "P11",
        "PENDING",
        "User requested a fresh biological-replicate holdout data layer; CPU-only generation and no locked-test access.",
    )
    save_state(state)
    return 0


def start_p13_submission_evidence(state: dict[str, Any]) -> int:
    """Open the bounded I--V-only measurement-reliability analysis."""

    if state.get("current_phase") != "P12" or state.get("status") != "COMPLETE":
        raise SystemExit("P13 may start only after the completed P12 workflow")
    state["current_phase"] = "P13"
    state["status"] = "PENDING"
    append_history(
        state,
        "P13",
        "PENDING",
        "User requested I--V-only replicate agreement for manuscript completion; chromosome X, locked-test access, GPU use and model selection are prohibited.",
    )
    save_state(state)
    return 0


def start_p14_iv_model_reference(state: dict[str, Any]) -> int:
    """Open bounded frozen-checkpoint inference after the P13 reference passes."""

    if state.get("current_phase") != "P13" or state.get("status") != "COMPLETE":
        raise SystemExit("P14 may start only after the completed P13 replicate-agreement phase")
    spec_path = REPO_ROOT / "alphagenome_custom/metadata/v2/p14_iv_frozen_model_reference_spec.json"
    if not spec_path.is_file():
        raise SystemExit("P14 requires its immutable preregistered frozen-model specification")
    state["current_phase"] = "P14"
    state["status"] = "PENDING"
    append_history(
        state,
        "P14",
        "PENDING",
        "User authorized frozen P12 Model B inference on the P13 I--V-only held-out contract; GPU 2/3 only, no training, chromosome-X or locked-test access.",
    )
    save_state(state)
    return 0


def start_p15_iv_training_contract(state: dict[str, Any]) -> int:
    """Open a metadata-only I--V contract before any new model training."""

    if state.get("current_phase") != "P14" or state.get("status") != "COMPLETE":
        raise SystemExit("P15 may start only after the completed P14 frozen-reference phase")
    review_path = REPO_ROOT / "alphagenome_custom/metadata/v2/audits/P14/review.json"
    if not review_path.is_file() or json.loads(review_path.read_text()).get("status") != "PASS":
        raise SystemExit("P15 requires the completed R14 pass record")
    state["current_phase"] = "P15"
    state["status"] = "PENDING"
    append_history(
        state,
        "P15",
        "PENDING",
        "User authorized an I--V-only train/validation contract before new ablations or public-baseline training; no signal, checkpoint or GPU access.",
    )
    save_state(state)
    return 0


def start_p16_iv_training_normalization(state: dict[str, Any]) -> int:
    """Open CPU-only I--V normalization after the P15 interval contract passes."""

    if state.get("current_phase") != "P15" or state.get("status") != "COMPLETE":
        raise SystemExit("P16 may start only after the completed P15 training-interval contract")
    spec_path = REPO_ROOT / "alphagenome_custom/metadata/v2/p16_iv_training_normalization_spec.json"
    review_path = REPO_ROOT / "alphagenome_custom/metadata/v2/audits/P15/review.json"
    if not spec_path.is_file() or not review_path.is_file():
        raise SystemExit("P16 requires a preregistered normalization specification and R15 pass record")
    if json.loads(review_path.read_text()).get("status") != "PASS":
        raise SystemExit("P16 requires the completed R15 pass record")
    state["current_phase"] = "P16"
    state["status"] = "PENDING"
    append_history(
        state,
        "P16",
        "PENDING",
        "User authorized strict I--V fold-normalization means before new ablations or public-baseline training; CPU-only with chromosome X and locked test prohibited.",
    )
    save_state(state)
    return 0


def start_p17_iv_controlled_matrix(state: dict[str, Any]) -> int:
    """Open the bounded I--V-only ablation and public-baseline matrix."""

    if state.get("current_phase") != "P16" or state.get("status") != "COMPLETE":
        raise SystemExit("P17 may start only after the completed P16 normalization phase")
    spec_path = REPO_ROOT / "alphagenome_custom/metadata/v2/p17_iv_controlled_matrix_spec.json"
    review_path = REPO_ROOT / "alphagenome_custom/metadata/v2/audits/P16/review.json"
    if not spec_path.is_file() or not review_path.is_file():
        raise SystemExit("P17 requires a preregistered matrix specification and R16 pass record")
    if json.loads(review_path.read_text()).get("status") != "PASS":
        raise SystemExit("P17 requires the completed R16 pass record")
    state["current_phase"] = "P17"
    state["status"] = "PENDING"
    append_history(
        state,
        "P17",
        "PENDING",
        "User authorized a minimal 45-job strict I--V controlled matrix: retrained full Model B, 128-bp-only Model B and a Basenji2-style from-scratch public-architecture baseline; GPU 2/3 only with chromosome-X and locked-test access prohibited.",
    )
    save_state(state)
    return 0


def start_p18_external_continuation(state: dict[str, Any]) -> int:
    """Open a separately scoped external-RNA continuation after P17."""

    if state.get("current_phase") != "P17" or state.get("status") != "COMPLETE":
        raise SystemExit("P18 may start only after the completed P17 controlled matrix")
    spec_path = REPO_ROOT / "alphagenome_custom/metadata/v2/p18_external_continuation_spec.json"
    review_path = REPO_ROOT / "alphagenome_custom/metadata/v2/audits/P17/review.json"
    if not spec_path.is_file() or not review_path.is_file():
        raise SystemExit("P18 requires a preregistered continuation contract and R17 pass record")
    if json.loads(review_path.read_text()).get("status") != "PASS":
        raise SystemExit("P18 requires the completed R17 pass record")
    state["current_phase"] = "P18"
    state["status"] = "PENDING"
    append_history(
        state,
        "P18",
        "PENDING",
        "User directed continuation after P17: checksum-verified replacement external RNA transfer and I-V-only reprocessing; model inference, chromosome X and locked test remain prohibited.",
    )
    save_state(state)
    return 0


def start_p19_external_scoring(state: dict[str, Any]) -> int:
    """Open the separately scoped external represented-head inference phase."""

    if state.get("current_phase") != "P18" or state.get("status") != "COMPLETE":
        raise SystemExit("P19 may start only after the completed P18 external reprocessing phase")
    spec_path = REPO_ROOT / "alphagenome_custom/metadata/v2/p19_external_represented_head_scoring_spec.json"
    review_path = REPO_ROOT / "alphagenome_custom/metadata/v2/audits/P18/review.json"
    if not spec_path.is_file() or not review_path.is_file():
        raise SystemExit("P19 requires a preregistered scoring contract and R18 pass record")
    if json.loads(review_path.read_text()).get("status") != "PASS":
        raise SystemExit("P19 requires the completed R18 pass record")
    state["current_phase"] = "P19"
    state["status"] = "PENDING"
    append_history(
        state,
        "P19",
        "PENDING",
        "User authorized external represented-head inference on checksum-verified P18 RNA tracks using frozen P17 Model B I--V checkpoints; no training, chromosome-X or locked-test signal access is allowed.",
    )
    save_state(state)
    return 0


def start_p20_eqtl_public_source_download(state: dict[str, Any]) -> int:
    """Open the public eQTL source acquisition and audit phase."""

    if state.get("current_phase") != "P19" or state.get("status") != "COMPLETE":
        raise SystemExit("P20 may start only after the completed P19 external scoring phase")
    spec_path = REPO_ROOT / "alphagenome_custom/metadata/v2/p20_eqtl_public_source_download_spec.json"
    review_path = REPO_ROOT / "alphagenome_custom/metadata/v2/audits/P19/review.json"
    if not spec_path.is_file() or not review_path.is_file():
        raise SystemExit("P20 requires a preregistered source specification and R19 pass record")
    if json.loads(review_path.read_text()).get("status") != "PASS":
        raise SystemExit("P20 requires the completed R19 pass record")
    state["current_phase"] = "P20"
    state["status"] = "PENDING"
    append_history(
        state,
        "P20",
        "PENDING",
        "User authorized public eQTL source acquisition: download the prespecified GEO expression matrix, published eQTL truth archive and CeNDR release; no model inference, chromosome-X or locked-test access is allowed.",
    )
    save_state(state)
    return 0


def start_p21_eqtl_variant_scoring(state: dict[str, Any]) -> int:
    """Open the registered external eQTL variant-scoring phase."""

    if state.get("current_phase") != "P20" or state.get("status") != "COMPLETE":
        raise SystemExit("P21 may start only after the completed P20 public-source phase")
    spec_path = REPO_ROOT / "results/v2_p20_eqtl_variant_set/p20_eqtl_variant_set_manifest.json"
    review_path = REPO_ROOT / "alphagenome_custom/metadata/v2/audits/P20/review.json"
    if not spec_path.is_file() or not review_path.is_file():
        raise SystemExit("P21 requires P20 source review and an I--V variant-set manifest")
    if json.loads(review_path.read_text()).get("status") != "PASS":
        raise SystemExit("P21 requires the completed R20 pass record")
    if json.loads(spec_path.read_text()).get("status") != "PASS":
        raise SystemExit("P21 requires a passing WBcel235 I--V coordinate audit")
    state["current_phase"] = "P21"
    state["status"] = "PENDING"
    append_history(
        state,
        "P21",
        "PENDING",
        "User authorized external I--V eQTL candidate scoring with the frozen P17 Model B checkpoint matrix; no training, chromosome-X signal or locked-test access is allowed.",
    )
    save_state(state)
    return 0


def start_p22_eqtl_effect_association(state: dict[str, Any]) -> int:
    """Open the CPU-only public eQTL effect-association phase."""

    if state.get("current_phase") != "P21" or state.get("status") != "COMPLETE":
        raise SystemExit("P22 may start only after the completed P21 variant-scoring phase")
    review_path = REPO_ROOT / "alphagenome_custom/metadata/v2/audits/P21/review.json"
    execution_path = REPO_ROOT / "alphagenome_custom/metadata/v2/p21_eqtl_variant_scoring_execution.json"
    if not review_path.is_file() or not execution_path.is_file():
        raise SystemExit("P22 requires the completed R21 review and P21 execution record")
    if json.loads(review_path.read_text()).get("status") != "PASS":
        raise SystemExit("P22 requires a passing R21 review")
    state["current_phase"] = "P22"
    state["status"] = "PENDING"
    append_history(
        state,
        "P22",
        "PENDING",
        "User authorized public wild-strain dosage-expression association analysis against frozen P21 represented-head variant effects; no X or locked-test access.",
    )
    save_state(state)
    return 0


def start_p23_lora_sensitivity(state: dict[str, Any]) -> int:
    """Open the I-V-only LoRA rank/location sensitivity phase."""

    if state.get("current_phase") != "P22" or state.get("status") != "COMPLETE":
        raise SystemExit("P23 may start only after the completed P22 association phase")
    review_path = REPO_ROOT / "alphagenome_custom/metadata/v2/audits/P22/review.json"
    spec_path = REPO_ROOT / "alphagenome_custom/metadata/v2/p23_lora_sensitivity_spec.json"
    if not review_path.is_file() or not spec_path.is_file():
        raise SystemExit("P23 requires the completed R22 review and a preregistered sensitivity spec")
    if json.loads(review_path.read_text()).get("status") != "PASS":
        raise SystemExit("P23 requires a passing R22 review")
    state["current_phase"] = "P23"
    state["status"] = "PENDING"
    append_history(
        state,
        "P23",
        "PENDING",
        "User authorized I-V-only LoRA rank/location sensitivity; physical GPUs 2/3 only, chromosome X and locked-test access prohibited.",
    )
    save_state(state)
    return 0


def start_p24_borzoi_adapter(state: dict[str, Any]) -> int:
    """Open the conditional Borzoi feasibility/adapter phase."""

    if state.get("current_phase") != "P23" or state.get("status") != "COMPLETE":
        raise SystemExit("P24 may start only after the completed P23 sensitivity phase")
    review_path = REPO_ROOT / "alphagenome_custom/metadata/v2/audits/P23/review.json"
    spec_path = REPO_ROOT / "alphagenome_custom/metadata/v2/p24_borzoi_adapter_spec.json"
    if not review_path.is_file() or not spec_path.is_file():
        raise SystemExit("P24 requires the completed R23 review and a preregistered Borzoi spec")
    if json.loads(review_path.read_text()).get("status") != "PASS":
        raise SystemExit("P24 requires a passing R23 review")
    state["current_phase"] = "P24"
    state["status"] = "PENDING"
    append_history(
        state,
        "P24",
        "PENDING",
        "User authorized the official Borzoi feasibility audit and conditional I-V-only adapter; no native human/mouse checkpoint leaderboard, chromosome X or locked-test access.",
    )
    save_state(state)
    return 0


def start_p25_borzoi_adapter_pilot(state: dict[str, Any]) -> int:
    """Open the task-faithful I-V-only Borzoi adapter pilot."""

    if state.get("current_phase") != "P24" or state.get("status") != "COMPLETE":
        raise SystemExit("P25 may start only after the completed P24 Borzoi source audit")
    review_path = REPO_ROOT / "alphagenome_custom/metadata/v2/audits/P24/review.json"
    spec_path = REPO_ROOT / "alphagenome_custom/metadata/v2/p25_borzoi_adapter_pilot_spec.json"
    if not review_path.is_file() or not spec_path.is_file():
        raise SystemExit("P25 requires the completed R24 review and a preregistered adapter pilot spec")
    if json.loads(review_path.read_text()).get("status") != "PASS":
        raise SystemExit("P25 requires a passing R24 review")
    state["current_phase"] = "P25"
    state["status"] = "PENDING"
    append_history(
        state,
        "P25",
        "PENDING",
        "User authorized the I-V-only task-faithful Borzoi C. elegans adapter pilot; native human/mouse heads, chromosome X and locked-test access remain prohibited.",
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
    run_parser.add_argument("--phase", choices=KNOWN_PHASES)
    run_parser.add_argument(
        "--auto",
        action="store_true",
        help="Continue through PASS phases until FAIL, an approval boundary, or completion.",
    )
    reopen_parser = subparsers.add_parser("reopen")
    reopen_parser.add_argument("--phase", required=True, choices=KNOWN_PHASES)
    reopen_parser.add_argument("--reason", required=True)
    subparsers.add_parser("prepare-six-chromosome-revision")
    subparsers.add_parser("start-p7-legacy-residual")
    subparsers.add_parser("start-p8-b-no-lora")
    subparsers.add_parser("start-p9-submission-evidence")
    subparsers.add_parser("start-p11-replicate-holdout")
    subparsers.add_parser("start-p13-submission-evidence")
    subparsers.add_parser("start-p14-iv-model-reference")
    subparsers.add_parser("start-p15-iv-training-contract")
    subparsers.add_parser("start-p16-iv-training-normalization")
    subparsers.add_parser("start-p17-iv-controlled-matrix")
    subparsers.add_parser("start-p18-external-continuation")
    subparsers.add_parser("start-p19-external-scoring")
    subparsers.add_parser("start-p20-eqtl-public-source-download")
    subparsers.add_parser("start-p21-eqtl-variant-scoring")
    subparsers.add_parser("start-p22-eqtl-effect-association")
    subparsers.add_parser("start-p23-lora-sensitivity")
    subparsers.add_parser("start-p24-borzoi-adapter")
    subparsers.add_parser("start-p25-borzoi-adapter-pilot")
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
        if args.phase in POST_COMPLETION_PHASES:
            if state.get("current_phase") not in {"P12", *POST_COMPLETION_PHASES}:
                raise SystemExit("Post-completion evidence phases may be reopened only after the completed P12 workflow")
        else:
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
    if args.command == "start-p9-submission-evidence":
        raise SystemExit(start_p9_submission_evidence(state))
    if args.command == "start-p11-replicate-holdout":
        raise SystemExit(start_p11_replicate_holdout(state))
    if args.command == "start-p13-submission-evidence":
        raise SystemExit(start_p13_submission_evidence(state))
    if args.command == "start-p14-iv-model-reference":
        raise SystemExit(start_p14_iv_model_reference(state))
    if args.command == "start-p15-iv-training-contract":
        raise SystemExit(start_p15_iv_training_contract(state))
    if args.command == "start-p16-iv-training-normalization":
        raise SystemExit(start_p16_iv_training_normalization(state))
    if args.command == "start-p17-iv-controlled-matrix":
        raise SystemExit(start_p17_iv_controlled_matrix(state))
    if args.command == "start-p18-external-continuation":
        raise SystemExit(start_p18_external_continuation(state))
    if args.command == "start-p19-external-scoring":
        raise SystemExit(start_p19_external_scoring(state))
    if args.command == "start-p20-eqtl-public-source-download":
        raise SystemExit(start_p20_eqtl_public_source_download(state))
    if args.command == "start-p21-eqtl-variant-scoring":
        raise SystemExit(start_p21_eqtl_variant_scoring(state))
    if args.command == "start-p22-eqtl-effect-association":
        raise SystemExit(start_p22_eqtl_effect_association(state))
    if args.command == "start-p23-lora-sensitivity":
        raise SystemExit(start_p23_lora_sensitivity(state))
    if args.command == "start-p24-borzoi-adapter":
        raise SystemExit(start_p24_borzoi_adapter(state))
    if args.command == "start-p25-borzoi-adapter-pilot":
        raise SystemExit(start_p25_borzoi_adapter_pilot(state))
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
