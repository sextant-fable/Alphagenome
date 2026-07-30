#!/usr/bin/env python3
"""Run one bounded development-only port of the legacy 1bp-B architecture."""

from __future__ import annotations

from datetime import datetime, timezone
import json
from pathlib import Path
import subprocess
import time
from typing import Any

from scripts import run_v2_p6b as base
from scripts import run_v2_p6b_amendment as shared
from scripts import v2_gpu_resources
from scripts import v2_subprocess


REPO_ROOT = Path(__file__).resolve().parents[1]
METADATA_DIR = REPO_ROOT / "alphagenome_custom/metadata/v2"
SPEC_PATH = METADATA_DIR / "p7_legacy_residual_fold1_spec.json"
EXECUTION_PATH = METADATA_DIR / "p7_legacy_residual_fold1_execution.json"
COMPARISON_PATH = METADATA_DIR / "p7_legacy_residual_fold1_comparison.json"
RUN_ROOT = Path("runs/v2_p7_legacy_residual_fold1_20260730")
LOG_ROOT = Path("logs/v2_p7_legacy_residual_fold1_20260730")


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def git_commit() -> str:
    return subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=REPO_ROOT,
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()


def require_contract(spec: dict[str, Any]) -> None:
    state = shared.read_json(METADATA_DIR / "execution_state.json")
    g4 = state.get("approvals", {}).get("G4_gpu_experiments", {})
    if not (
        state.get("current_phase") == "P7"
        and g4.get("approved") is True
        and g4.get("scope") == "p7_single_legacy_residual_fold1"
    ):
        raise RuntimeError("P7 requires its exact scoped G4 approval")
    if spec.get("phase") != "P7" or spec.get("model") != "D":
        raise RuntimeError("P7 spec must register model D")
    training = spec.get("training", {})
    if not (
        training.get("fold") == 1
        and training.get("seed") == 20260714
        and training.get("loss") == "paper"
        and training.get("base_max_steps") == 2000
        and training.get("residual_max_steps") == 1500
        and training.get("sequence_length") == 131072
    ):
        raise RuntimeError("P7 spec must retain the bounded fold-1 training contract")
    if spec.get("final_test_access") != "prohibited":
        raise RuntimeError("P7 must prohibit final-test access")
    for key, relative in spec["locked_inputs"].items():
        expected = spec["input_sha256"][key]
        if base.sha256(REPO_ROOT / relative) != expected:
            raise RuntimeError(f"P7 input hash mismatch: {key}")
    lock = shared.read_json(METADATA_DIR / "final_test_lock.json")
    if not (
        lock.get("test_consumed") is True
        and lock.get("test_status") == "completed"
    ):
        raise RuntimeError("P7 requires the preserved completed P6C test lock")


def commands(spec: dict[str, Any]) -> tuple[list[str], list[str], list[str]]:
    training = spec["training"]
    base_root = RUN_ROOT / "base_128bp" / "fold_1"
    residual_root = RUN_ROOT / "model_D_paper" / "fold_1"
    base_train = v2_subprocess.module_command(
        "train_v2_model",
        "--model", "D_base",
        "--fold", str(training["fold"]),
        "--seed", str(training["seed"]),
        "--loss", training["loss"],
        "--max-steps", str(training["base_max_steps"]),
        "--sequence-length", str(training["sequence_length"]),
        "--learning-rate", str(training["learning_rate"]),
        "--weight-decay", str(training["weight_decay"]),
        "--gene-weight", str(training["gene_weight"]),
        "--max-shift-bp", str(training["max_shift_bp"]),
        "--reverse-complement-probability",
        str(training["reverse_complement_probability"]),
        "--hidden-channels", str(training["hidden_channels"]),
        "--mean-column", training["mean_column"],
        "--output-dir", str(base_root),
        "--device", "cuda",
    )
    residual_train = v2_subprocess.module_command(
        "train_v2_model",
        "--model", "D",
        "--fold", str(training["fold"]),
        "--seed", str(training["seed"]),
        "--loss", training["loss"],
        "--max-steps", str(training["residual_max_steps"]),
        "--sequence-length", str(training["sequence_length"]),
        "--learning-rate", str(training["learning_rate"]),
        "--weight-decay", str(training["weight_decay"]),
        "--gene-weight", str(training["gene_weight"]),
        "--max-shift-bp", str(training["max_shift_bp"]),
        "--reverse-complement-probability",
        str(training["reverse_complement_probability"]),
        "--hidden-channels", str(training["hidden_channels"]),
        "--mean-column", training["mean_column"],
        "--frozen-base-checkpoint", str(base_root / "checkpoint.pt"),
        "--output-dir", str(residual_root),
        "--device", "cuda",
    )
    evaluate = v2_subprocess.module_command(
        "evaluate_v2_model",
        "--model", "D",
        "--training-loss", training["loss"],
        "--fold", str(training["fold"]),
        "--checkpoint", str(residual_root / "checkpoint.pt"),
        "--sequence-length", str(training["sequence_length"]),
        "--hidden-channels", str(training["hidden_channels"]),
        "--mean-column", training["mean_column"],
        "--frozen-base-checkpoint", str(base_root / "checkpoint.pt"),
        "--output", str(residual_root / "validation_full_metrics.json"),
        "--device", "cuda",
    )
    if any(
        "--final-test" in command
        for command in (base_train, residual_train, evaluate)
    ):
        raise RuntimeError("P7 command unexpectedly requests final-test access")
    return base_train, residual_train, evaluate


def main() -> None:
    spec = shared.read_json(SPEC_PATH)
    require_contract(spec)
    spec_sha = base.sha256(SPEC_PATH)
    selected_gpus, snapshots = base.wait_for_gpus(
        int(spec["gpu_policy"]["minimum_free_mib"])
    )
    gpu = selected_gpus[0]
    base_root = RUN_ROOT / "base_128bp" / "fold_1"
    root = RUN_ROOT / "model_D_paper" / "fold_1"
    base_output_dir = REPO_ROOT / base_root
    output_dir = REPO_ROOT / root
    log_path = REPO_ROOT / LOG_ROOT / "model_D_paper" / "fold_1.log"
    checkpoint_path = output_dir / "checkpoint.pt"
    base_checkpoint_path = base_output_dir / "checkpoint.pt"
    validation_path = output_dir / "validation_full_metrics.json"
    execution: dict[str, Any] = {
        "schema_version": 1,
        "phase": "P7",
        "status": "running",
        "started_at": utc_now(),
        "completed_at": None,
        "git_commit": git_commit(),
        "spec_path": str(SPEC_PATH.relative_to(REPO_ROOT)),
        "spec_sha256": spec_sha,
        "model": "D",
        "loss": "paper",
        "fold": 1,
        "seed": 20260714,
        "physical_gpu": gpu,
        "gpu_wait_snapshots": snapshots,
        "locked_test_block_signal_reads": 0,
        "final_test_access": "prohibited",
        "output_dir": str(root),
        "log_path": str(log_path.relative_to(REPO_ROOT)),
        "failures": [],
    }
    base.atomic_json(EXECUTION_PATH, execution)
    started = time.monotonic()
    try:
        base_train_command, residual_train_command, evaluate_command = commands(spec)
        base.run_command(base_train_command, log_path, gpu)
        base.run_command(residual_train_command, log_path, gpu)
        base.run_command(evaluate_command, log_path, gpu)
        base_run = shared.read_json(base_output_dir / "run.json")
        run = shared.read_json(output_dir / "run.json")
        validation = shared.read_json(validation_path)
        candidate = shared.validation_summary(validation)
        baseline_path = REPO_ROOT / spec["baseline"]["validation_path"]
        baseline = shared.validation_summary(shared.read_json(baseline_path))
        primary_delta = (
            candidate["primary_biological_score"]
            - baseline["primary_biological_score"]
        )
        comparison = {
            "schema_version": 1,
            "phase": "P7",
            "comparison_type": "single_fold_exploratory_architecture_port",
            "spec_sha256": spec_sha,
            "candidate": {
                "model": "D",
                "loss": "paper",
                "fold": 1,
                "seed": 20260714,
                "checkpoint_path": str(checkpoint_path.relative_to(REPO_ROOT)),
                "checkpoint_sha256": base.sha256(checkpoint_path),
                **candidate,
            },
            "baseline": {
                "model": "B",
                "loss": "paper",
                "fold": 1,
                "seed": 20260714,
                "checkpoint_path": spec["baseline"]["checkpoint_path"],
                "checkpoint_sha256": spec["input_sha256"]["baseline_checkpoint"],
                **baseline,
            },
            "delta": {
                key: candidate[key] - baseline[key]
                for key in candidate
                if key in baseline
            },
            "primary_biological_score_delta": primary_delta,
            "interpretation": "One fold and one seed only; exploratory architecture check, not a new model selection or final-test result.",
            "created_at": utc_now(),
        }
        base.atomic_json(COMPARISON_PATH, comparison)
        execution.update(
            {
                "status": "completed",
                "completed_at": utc_now(),
                "elapsed_seconds": time.monotonic() - started,
                "run_path": str((output_dir / "run.json").relative_to(REPO_ROOT)),
                "base_run_path": str(
                    (base_output_dir / "run.json").relative_to(REPO_ROOT)
                ),
                "base_checkpoint_path": str(base_checkpoint_path.relative_to(REPO_ROOT)),
                "base_checkpoint_sha256": base.sha256(base_checkpoint_path),
                "checkpoint_path": str(checkpoint_path.relative_to(REPO_ROOT)),
                "checkpoint_sha256": base.sha256(checkpoint_path),
                "validation_path": str(validation_path.relative_to(REPO_ROOT)),
                "validation_sha256": base.sha256(validation_path),
                "comparison_path": str(COMPARISON_PATH.relative_to(REPO_ROOT)),
                "comparison_sha256": base.sha256(COMPARISON_PATH),
                "candidate_metrics": candidate,
                "baseline_metrics": baseline,
                "primary_biological_score_delta": primary_delta,
                "run_checkpoint_sha256": run["checkpoint_sha256"],
                "base_run_checkpoint_sha256": base_run["checkpoint_sha256"],
                "resources_after": v2_gpu_resources.snapshot(),
            }
        )
        base.atomic_json(EXECUTION_PATH, execution)
    except Exception as error:
        execution["status"] = "failed"
        execution["completed_at"] = utc_now()
        execution["failures"] = [str(error)]
        base.atomic_json(EXECUTION_PATH, execution)
        raise


if __name__ == "__main__":
    main()
