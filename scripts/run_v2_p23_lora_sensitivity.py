#!/usr/bin/env python3
"""Run the registered I-V-only LoRA sensitivity pilot."""

from __future__ import annotations

import csv
import hashlib
import json
import os
from pathlib import Path
import socket
import sys
import time

from scripts import run_v2_p6b as base
from scripts import run_v2_p6b_amendment as shared


REPO_ROOT = Path(__file__).resolve().parents[1]
METADATA = REPO_ROOT / "alphagenome_custom/metadata/v2"
SPEC_PATH = METADATA / "p23_lora_sensitivity_spec.json"
EXECUTION_PATH = METADATA / "p23_lora_sensitivity_execution.json"
RESULT_ROOT = REPO_ROOT / "results/v2_p23_lora_sensitivity_pilot"
RUN_ROOT = Path("runs/v2_p23_lora_sensitivity_pilot")
LOG_ROOT = Path("logs/v2_p23_lora_sensitivity_pilot")
ALLOWED_CHROMOSOMES = ["I", "II", "III", "IV", "V"]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_json(path: Path) -> dict:
    return json.loads(path.read_text())


def write_tsv(path: Path, rows: list[dict]) -> None:
    if not rows:
        raise RuntimeError(f"Refusing to write empty result: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def require_scope() -> None:
    state = read_json(METADATA / "execution_state.json")
    gate = state.get("approvals", {}).get("G4_gpu_experiments", {})
    if not (
        state.get("current_phase") == "P23"
        and state.get("status") == "RUNNING"
        and gate.get("approved") is True
        and gate.get("scope") == "p23_lora_sensitivity"
    ):
        raise RuntimeError("P23 requires controller P23/RUNNING and G4:p23_lora_sensitivity")


def command_pair(config: dict, spec: dict) -> tuple[list[str], list[str]]:
    training = spec["training"]
    fold = 1
    seed = 20260714
    train_interval = REPO_ROOT / spec["source_contract"]["interval_root"] / f"fold_{fold}/train.tsv"
    valid_interval = REPO_ROOT / spec["source_contract"]["interval_root"] / f"fold_{fold}/valid.tsv"
    root = RUN_ROOT / config["id"] / f"seed_{seed}" / f"fold_{fold}"
    evaluation = RESULT_ROOT / "evaluations" / config["id"] / f"seed_{seed}" / f"fold_{fold}.json"
    log_path = REPO_ROOT / LOG_ROOT / config["id"] / f"seed_{seed}" / f"fold_{fold}.log"
    common_train = [
        str(Path(sys.executable).resolve()), "-m", "scripts.train_v2_model",
        "--model", "B", "--fold", str(fold), "--seed", str(seed), "--loss", "paper",
        "--max-steps", str(training["max_steps"]), "--sequence-length", str(training["sequence_length"]),
        "--learning-rate", str(training["learning_rate"]), "--weight-decay", str(training["weight_decay"]),
        "--gene-weight", str(training["gene_weight"]), "--loss-1bp-weight", str(training["loss_1bp_weight"]),
        "--loss-128bp-weight", str(training["loss_128bp_weight"]), "--max-shift-bp", str(training["max_shift_bp"]),
        "--reverse-complement-probability", str(training["reverse_complement_probability"]),
        "--hidden-channels", str(training["hidden_channels"]), "--mean-column", "fold_1_train_nonzero_mean",
        "--means-path", training["means_path"],
        "--track-manifest", spec["source_contract"]["track_manifest"],
        "--group-manifest", spec["source_contract"]["group_manifest"],
        "--train-intervals-path", str(train_interval.relative_to(REPO_ROOT)),
        "--output-dir", str(root), "--device", "cuda",
        "--lora-rank", str(config["lora_rank"]), "--lora-alpha", str(config["lora_alpha"]),
        "--lora-target-set", config["lora_target_set"],
    ]
    evaluate = [
        str(Path(sys.executable).resolve()), "-m", "scripts.evaluate_v2_model",
        "--model", "B", "--training-loss", "paper", "--fold", str(fold),
        "--checkpoint", str(root / "checkpoint.pt"), "--sequence-length", str(training["sequence_length"]),
        "--hidden-channels", str(training["hidden_channels"]), "--mean-column", "fold_1_train_nonzero_mean",
        "--means-path", training["means_path"],
        "--track-manifest", spec["source_contract"]["evaluation_track_manifest"],
        "--group-manifest", spec["source_contract"]["evaluation_group_manifest"],
        "--prediction-track-indices", spec["source_contract"]["prediction_track_indices"],
        "--valid-intervals-path", str(valid_interval.relative_to(REPO_ROOT)),
        "--output", str(evaluation), "--device", "cuda",
        "--lora-rank", str(config["lora_rank"]), "--lora-alpha", str(config["lora_alpha"]),
        "--lora-target-set", config["lora_target_set"],
    ]
    serialized = " ".join([*common_train, *evaluate])
    if any(token in serialized for token in ("test_locked", "--final-test", "/fold_0/", "chromosome_X")):
        raise RuntimeError("P23 command violates the I-V/no-final-test contract")
    return common_train, evaluate, log_path, root, evaluation


def main() -> None:
    require_scope()
    spec = read_json(SPEC_PATH)
    if spec.get("status") != "preregistered" or spec.get("allowed_chromosomes") != ALLOWED_CHROMOSOMES:
        raise RuntimeError("P23 static specification mismatch")
    spec_sha = sha256(SPEC_PATH)
    selected_gpus, snapshots = base.wait_for_gpus(spec["gpu_policy"]["minimum_free_mib"])
    if not selected_gpus or not set(selected_gpus) <= {2, 3}:
        raise RuntimeError(f"P23 selected prohibited GPUs: {selected_gpus}")
    execution = {
        "schema_version": 1, "phase": "P23", "status": "running", "started_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "hostname": socket.gethostname(), "working_directory": str(REPO_ROOT),
        "controller_invocation": "python -m scripts.v2_phase_controller run --phase P23",
        "registered_phase_command": "python -m scripts.run_v2_p23_lora_sensitivity",
        "spec_path": str(SPEC_PATH.relative_to(REPO_ROOT)), "spec_sha256": spec_sha,
        "selected_physical_gpus": selected_gpus, "gpu_wait_snapshots": snapshots,
        "registered_tasks": len(spec["pilot_configurations"]), "jobs": [], "failures": [],
        "allowed_chromosomes": ALLOWED_CHROMOSOMES, "final_test_access": "prohibited",
        "locked_test_block_signal_reads": 0,
    }
    base.atomic_json(EXECUTION_PATH, execution)
    try:
        for index, config in enumerate(spec["pilot_configurations"]):
            gpu = selected_gpus[index % len(selected_gpus)]
            train, evaluate, log_path, root, evaluation = command_pair(config, spec)
            started = time.monotonic()
            base.run_command(train, log_path, gpu)
            base.run_command(evaluate, log_path, gpu)
            if not (REPO_ROOT / root / "checkpoint.pt").is_file() or not evaluation.is_file():
                raise RuntimeError(f"P23 missing output for {config['id']}")
            eval_record = read_json(evaluation)
            if eval_record.get("evaluated_chromosomes") != ALLOWED_CHROMOSOMES or eval_record.get("locked_test_block_signal_reads") is not False:
                raise RuntimeError(f"P23 evaluation scope failed for {config['id']}")
            run_record = read_json(REPO_ROOT / root / "run.json")
            metrics = shared.validation_summary(eval_record)
            execution["jobs"].append({
                "configuration": config["id"], "fold": 1, "seed": 20260714,
                "physical_gpu": gpu, "elapsed_seconds": time.monotonic() - started,
                "checkpoint_path": str((REPO_ROOT / root / "checkpoint.pt").relative_to(REPO_ROOT)),
                "checkpoint_sha256": sha256(REPO_ROOT / root / "checkpoint.pt"),
                "run_sha256": sha256(REPO_ROOT / root / "run.json"),
                "evaluation_path": str(evaluation.relative_to(REPO_ROOT)),
                "evaluation_sha256": sha256(evaluation),
                "trainable_parameters": int(run_record["trainable_parameters"]),
                "primary_biological_score": float(metrics["primary_biological_score"]),
                "gene_exon_coverage_pearson_log1p": float(metrics["gene_exon_coverage_pearson_log1p"]),
                "per_track_pearson_128bp_log1p": float(metrics["per_track_pearson_128bp_log1p"]),
                "evaluated_chromosomes": ",".join(ALLOWED_CHROMOSOMES),
                "locked_test_block_signal_reads": False,
            })
            execution["registered_tasks_completed"] = len(execution["jobs"])
            base.atomic_json(EXECUTION_PATH, execution)
        rows = execution["jobs"]
        write_tsv(RESULT_ROOT / "p23_lora_sensitivity_pilot.tsv", rows)
        execution.update({"status": "completed", "completed_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), "records_sha256": sha256(RESULT_ROOT / "p23_lora_sensitivity_pilot.tsv")})
        base.atomic_json(EXECUTION_PATH, execution)
    except Exception:
        execution.update({"status": "failed", "completed_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())})
        base.atomic_json(EXECUTION_PATH, execution)
        raise
    print(json.dumps({"status": "completed", "registered_tasks_completed": len(execution["jobs"])}, indent=2))


if __name__ == "__main__":
    main()
