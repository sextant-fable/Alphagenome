#!/usr/bin/env python3
"""Run the pre-registered five-fold A/B/C and loss comparison on GPU 2/3."""

from __future__ import annotations

from concurrent.futures import Future, ThreadPoolExecutor, as_completed
import csv
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import time
from typing import Any

from scripts import v2_gpu_resources
from scripts import v2_subprocess


REPO_ROOT = Path(__file__).resolve().parents[1]
METADATA_DIR = REPO_ROOT / "alphagenome_custom/metadata/v2"
STATE_PATH = METADATA_DIR / "execution_state.json"
SPEC_PATH = METADATA_DIR / "p6b_matrix_spec.json"
EXECUTION_PATH = METADATA_DIR / "p6b_execution.json"
RESULTS_PATH = METADATA_DIR / "p6b_cv_results.tsv"
SELECTION_PATH = METADATA_DIR / "p6b_selection.json"
LOCK_PATH = METADATA_DIR / "final_test_lock.json"
RUN_ROOT = REPO_ROOT / "runs/v2_p6b_20260714"
LOG_ROOT = REPO_ROOT / "logs/v2_p6b_20260714"


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(8 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def atomic_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)


def require_approval() -> None:
    state = json.loads(STATE_PATH.read_text())
    gate = state.get("approvals", {}).get("G4_gpu_experiments", {})
    if not (
        state.get("current_phase") == "P6B"
        and gate.get("approved") is True
        and gate.get("scope") == "r6_gpu_auto_available_2_3"
    ):
        raise RuntimeError("P6B requires current phase P6B and standing scoped G4")


def wait_for_gpus(
    minimum_free_mib: int,
    *,
    timeout_seconds: int = 24 * 60 * 60,
    poll_seconds: int = 60,
) -> tuple[list[int], list[dict[str, Any]]]:
    started = time.monotonic()
    snapshots = []
    while True:
        snapshot = v2_gpu_resources.snapshot()
        snapshot["observed_at"] = utc_now()
        snapshots.append(snapshot)
        selected = []
        for count in (2, 1):
            try:
                selected = v2_gpu_resources.select_available(
                    snapshot, count=count, minimum_free_mib=minimum_free_mib
                )
                break
            except RuntimeError:
                continue
        if selected:
            return selected, snapshots
        if time.monotonic() - started >= timeout_seconds:
            raise RuntimeError("No allowed GPU 2/3 became available within 24 hours")
        print("gpu_wait\tno idle GPU 2/3; retrying", flush=True)
        time.sleep(poll_seconds)


def run_command(command: list[str], log_path: Path, gpu: int) -> None:
    environment = os.environ.copy()
    environment["CUDA_VISIBLE_DEVICES"] = str(gpu)
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("a") as log:
        log.write("command\t" + " ".join(command) + "\n")
        log.flush()
        subprocess.run(
            command,
            cwd=REPO_ROOT,
            env=environment,
            stdout=log,
            stderr=subprocess.STDOUT,
            check=True,
        )


def run_cv_job(
    job: dict[str, Any],
    gpu: int,
    spec: dict[str, Any],
    spec_sha256: str,
) -> dict[str, Any]:
    model = job["model"]
    loss = job["loss"]
    fold = int(job["fold"])
    relative_root = Path("runs/v2_p6b_20260714") / model / loss / f"fold_{fold}"
    output_dir = REPO_ROOT / relative_root
    log_path = LOG_ROOT / model / loss / f"fold_{fold}.log"
    job_path = output_dir / "job.json"
    if job_path.is_file():
        existing = json.loads(job_path.read_text())
        validation_path = REPO_ROOT / existing.get("validation_path", "missing")
        checkpoint_path = REPO_ROOT / existing.get("checkpoint_path", "missing")
        if (
            existing.get("status") == "completed"
            and existing.get("spec_sha256") == spec_sha256
            and validation_path.is_file()
            and checkpoint_path.is_file()
            and sha256(checkpoint_path) == existing.get("checkpoint_sha256")
        ):
            existing["reused"] = True
            return existing
    output_dir.mkdir(parents=True, exist_ok=True)
    training = spec["training"]
    train_command = v2_subprocess.module_command(
        "train_v2_model",
        "--model",
        model,
        "--fold",
        str(fold),
        "--seed",
        str(training["seed"]),
        "--loss",
        loss,
        "--max-steps",
        str(training["max_steps"]),
        "--sequence-length",
        str(training["sequence_length"]),
        "--learning-rate",
        str(training["learning_rate"]),
        "--weight-decay",
        str(training["weight_decay"]),
        "--gene-weight",
        str(training["gene_weight"]),
        "--max-shift-bp",
        str(spec["augmentation"]["max_shift_bp"]),
        "--reverse-complement-probability",
        str(spec["augmentation"]["reverse_complement_probability"]),
        "--hidden-channels",
        str(training["hidden_channels"]),
        "--output-dir",
        str(relative_root),
        "--device",
        "cuda",
    )
    validation_relative = relative_root / "validation.json"
    evaluate_command = v2_subprocess.module_command(
        "evaluate_v2_model",
        "--model",
        model,
        "--training-loss",
        loss,
        "--fold",
        str(fold),
        "--checkpoint",
        str(relative_root / "checkpoint.pt"),
        "--sequence-length",
        str(spec["validation"]["sequence_length"]),
        "--hidden-channels",
        str(training["hidden_channels"]),
        "--output",
        str(validation_relative),
        "--device",
        "cuda",
    )
    started = time.monotonic()
    run_command(train_command, log_path, gpu)
    run_command(evaluate_command, log_path, gpu)
    run_path = output_dir / "run.json"
    checkpoint_path = output_dir / "checkpoint.pt"
    validation_path = REPO_ROOT / validation_relative
    run = json.loads(run_path.read_text())
    validation = json.loads(validation_path.read_text())
    record = {
        "schema_version": 1,
        "status": "completed",
        "model": model,
        "loss": loss,
        "fold": fold,
        "seed": training["seed"],
        "steps": training["max_steps"],
        "sequence_length": training["sequence_length"],
        "physical_gpu": gpu,
        "spec_sha256": spec_sha256,
        "run_path": str(run_path.relative_to(REPO_ROOT)),
        "validation_path": str(validation_path.relative_to(REPO_ROOT)),
        "checkpoint_path": str(checkpoint_path.relative_to(REPO_ROOT)),
        "checkpoint_sha256": sha256(checkpoint_path),
        "log_path": str(log_path.relative_to(REPO_ROOT)),
        "paper_loss": validation["mean_metrics"]["paper_loss"],
        "log1p_mse": validation["mean_metrics"]["log1p_mse"],
        "mean_per_track_pearson_128bp": validation[
            "mean_per_track_pearson_128bp"
        ],
        "elapsed_seconds": time.monotonic() - started,
        "reused": False,
        "run_checkpoint_sha256": run["checkpoint_sha256"],
    }
    atomic_json(job_path, record)
    return record


def write_results(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    grouped: dict[tuple[str, str], list[dict[str, Any]]] = {}
    for row in rows:
        grouped.setdefault((row["model"], row["loss"]), []).append(row)
    results = []
    for (model, loss), members in sorted(grouped.items()):
        if len(members) != 5 or {int(row["fold"]) for row in members} != set(range(1, 6)):
            raise RuntimeError(f"Incomplete five-fold result for {model}/{loss}")
        results.append(
            {
                "model": model,
                "loss": loss,
                "folds": 5,
                "mean_five_fold_paper_loss": sum(
                    float(row["paper_loss"]) for row in members
                )
                / 5,
                "mean_five_fold_log1p_mse": sum(
                    float(row["log1p_mse"]) for row in members
                )
                / 5,
                "mean_five_fold_per_track_pearson_128bp": sum(
                    float(row["mean_per_track_pearson_128bp"]) for row in members
                )
                / 5,
            }
        )
    temporary = RESULTS_PATH.with_suffix(".tsv.tmp")
    with temporary.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(results[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(results)
    temporary.replace(RESULTS_PATH)
    return results


def run_development(
    selected: dict[str, Any],
    gpu: int,
    spec: dict[str, Any],
    spec_sha256: str,
) -> dict[str, Any]:
    model = selected["model"]
    loss = selected["loss"]
    relative_root = Path("runs/v2_p6b_20260714/development_selected")
    output_dir = REPO_ROOT / relative_root
    log_path = LOG_ROOT / "development_selected.log"
    output_dir.mkdir(parents=True, exist_ok=True)
    training = spec["training"]
    max_steps = spec["development_retrain"]["max_steps"]
    command = v2_subprocess.module_command(
        "train_v2_model",
        "--model",
        model,
        "--fold",
        "0",
        "--seed",
        str(training["seed"]),
        "--loss",
        loss,
        "--max-steps",
        str(max_steps),
        "--sequence-length",
        str(training["sequence_length"]),
        "--learning-rate",
        str(training["learning_rate"]),
        "--weight-decay",
        str(training["weight_decay"]),
        "--gene-weight",
        str(training["gene_weight"]),
        "--max-shift-bp",
        str(spec["augmentation"]["max_shift_bp"]),
        "--reverse-complement-probability",
        str(spec["augmentation"]["reverse_complement_probability"]),
        "--hidden-channels",
        str(training["hidden_channels"]),
        "--output-dir",
        str(relative_root),
        "--device",
        "cuda",
    )
    run_command(command, log_path, gpu)
    run_path = output_dir / "run.json"
    checkpoint_path = output_dir / "checkpoint.pt"
    run = json.loads(run_path.read_text())
    return {
        "model": model,
        "loss": loss,
        "fold": 0,
        "seed": training["seed"],
        "steps": max_steps,
        "spec_sha256": spec_sha256,
        "physical_gpu": gpu,
        "run_path": str(run_path.relative_to(REPO_ROOT)),
        "checkpoint_path": str(checkpoint_path.relative_to(REPO_ROOT)),
        "checkpoint_sha256": sha256(checkpoint_path),
        "intervals_sha256": run["intervals_sha256"],
        "means_sha256": run["means_sha256"],
        "log_path": str(log_path.relative_to(REPO_ROOT)),
    }


def main() -> None:
    require_approval()
    spec = json.loads(SPEC_PATH.read_text())
    if spec.get("chromosome_x_access") != "prohibited":
        raise RuntimeError("P6B specification must prohibit chromosome X")
    spec_sha = sha256(SPEC_PATH)
    selected_gpus, waiting_snapshots = wait_for_gpus(
        int(spec["gpu_policy"]["minimum_free_mib"])
    )
    RUN_ROOT.mkdir(parents=True, exist_ok=True)
    LOG_ROOT.mkdir(parents=True, exist_ok=True)
    jobs = [
        {"fold": fold, "model": model, "loss": loss}
        for fold in spec["cv_folds"]
        for model in spec["models"]
        for loss in spec["losses"]
    ]
    execution = {
        "schema_version": 1,
        "phase": "P6B",
        "status": "running",
        "started_at": utc_now(),
        "completed_at": None,
        "spec_path": str(SPEC_PATH.relative_to(REPO_ROOT)),
        "spec_sha256": spec_sha,
        "selected_physical_gpus": selected_gpus,
        "gpu_wait_snapshots": waiting_snapshots,
        "jobs_expected": len(jobs),
        "jobs": [],
    }
    atomic_json(EXECUTION_PATH, execution)
    failures = []
    job_iterator = iter(jobs)
    with ThreadPoolExecutor(max_workers=len(selected_gpus)) as executor:
        futures: dict[Future[dict[str, Any]], int] = {}
        for gpu in selected_gpus:
            try:
                job = next(job_iterator)
            except StopIteration:
                break
            futures[executor.submit(run_cv_job, job, gpu, spec, spec_sha)] = gpu
        while futures:
            future = next(as_completed(futures))
            gpu = futures.pop(future)
            try:
                record = future.result()
                execution["jobs"].append(record)
            except Exception as error:
                failures.append(f"GPU{gpu}:{type(error).__name__}:{error}")
            execution["failures"] = failures
            atomic_json(EXECUTION_PATH, execution)
            if not failures:
                try:
                    job = next(job_iterator)
                except StopIteration:
                    continue
                futures[executor.submit(run_cv_job, job, gpu, spec, spec_sha)] = gpu
    if failures or len(execution["jobs"]) != len(jobs):
        execution["status"] = "failed"
        execution["completed_at"] = utc_now()
        atomic_json(EXECUTION_PATH, execution)
        raise RuntimeError(f"P6B matrix failed: {failures}")

    results = write_results(execution["jobs"])
    selected = min(
        results,
        key=lambda row: (
            float(row["mean_five_fold_paper_loss"]),
            float(row["mean_five_fold_log1p_mse"]),
            -float(row["mean_five_fold_per_track_pearson_128bp"]),
            row["model"],
            row["loss"],
        ),
    )
    selection = {
        "schema_version": 1,
        "selected": selected,
        "selection_rule": spec["selection"],
        "cv_results_path": str(RESULTS_PATH.relative_to(REPO_ROOT)),
        "cv_results_sha256": sha256(RESULTS_PATH),
        "spec_sha256": spec_sha,
        "chromosome_x_read": False,
        "selected_at": utc_now(),
    }
    atomic_json(SELECTION_PATH, selection)
    development = run_development(
        selected, selected_gpus[0], spec, spec_sha
    )
    lock = {
        "schema_version": 1,
        "checkpoint_path": development["checkpoint_path"],
        "checkpoint_sha256": development["checkpoint_sha256"],
        "model": development["model"],
        "loss": development["loss"],
        "seed": development["seed"],
        "steps": development["steps"],
        "training_chromosomes": ["I", "II", "III", "IV", "V"],
        "selection_path": str(SELECTION_PATH.relative_to(REPO_ROOT)),
        "selection_sha256": sha256(SELECTION_PATH),
        "spec_sha256": spec_sha,
        "legacy_chr_x_prior_exposure_disclosed": True,
        "test_consumed": False,
        "test_consumed_at": None,
        "locked_at": utc_now(),
    }
    atomic_json(LOCK_PATH, lock)
    execution["status"] = "completed"
    execution["completed_at"] = utc_now()
    execution["failures"] = []
    execution["cv_results_sha256"] = sha256(RESULTS_PATH)
    execution["selection_sha256"] = sha256(SELECTION_PATH)
    execution["development"] = development
    execution["lock_sha256"] = sha256(LOCK_PATH)
    execution["resources_after"] = v2_gpu_resources.snapshot()
    atomic_json(EXECUTION_PATH, execution)
    print(json.dumps(execution, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
