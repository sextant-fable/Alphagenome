#!/usr/bin/env python3
"""Complete the preregistered P6B amendment without accessing chromosome X."""

from __future__ import annotations

import csv
from datetime import datetime, timezone
import json
import math
import os
from pathlib import Path
import queue
import shutil
import threading
import time
from typing import Any

from scripts import run_v2_p6b as base
from scripts import v2_gpu_resources
from scripts import v2_subprocess


REPO_ROOT = Path(__file__).resolve().parents[1]
METADATA_DIR = REPO_ROOT / "alphagenome_custom/metadata/v2"
SPEC_PATH = METADATA_DIR / "p6b_amendment_spec.json"
BASE_EXECUTION_PATH = METADATA_DIR / "p6b_execution.json"
EXECUTION_PATH = METADATA_DIR / "p6b_amendment_execution.json"
RESULTS_PATH = METADATA_DIR / "p6b_amendment_cv_results.tsv"
ABLATION_RESULTS_PATH = METADATA_DIR / "p6b_amendment_ablation_results.tsv"
SELECTION_PATH = METADATA_DIR / "p6b_amendment_selection.json"
LOCK_PATH = METADATA_DIR / "final_test_lock.json"
OLD_LOCK_PATH = METADATA_DIR / "final_test_lock_original_p6b.json"
RUN_ROOT = Path("runs/v2_p6b_amendment_20260720")
LOG_ROOT = Path("logs/v2_p6b_amendment_20260720")


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text())


def add_option(command: list[str], name: str, value: Any) -> None:
    command.extend((name, str(value)))


def validation_summary(validation: dict[str, Any]) -> dict[str, float]:
    metrics = validation.get("full_metrics", {})
    if metrics.get("metric_contract") != "v2_full_metrics_1":
        raise RuntimeError("Validation lacks v2_full_metrics_1")
    primary = metrics["primary"]
    gene = float(primary["mean_per_track_gene_exon_coverage_pearson_log1p"])
    track = float(primary["mean_per_track_pearson_128bp_log1p"])
    values = {
        "primary_biological_score": 0.5 * gene + 0.5 * track,
        "gene_exon_coverage_pearson_log1p": gene,
        "per_track_pearson_128bp_log1p": track,
        "paper_loss": float(validation["mean_metrics"]["paper_loss"]),
        "log1p_mse": float(validation["mean_metrics"]["log1p_mse"]),
        "spearman_128bp": float(
            metrics["distribution_128bp"]["mean_per_track_spearman"]
        ),
        "top1_mse_128bp": float(
            metrics["distribution_128bp"]["mean_per_track_top1_mse"]
        ),
        "top1_calibration_ratio_128bp": float(
            metrics["distribution_128bp"][
                "mean_per_track_top1_calibration_ratio"
            ]
        ),
        "gene_body_mse_log1p_1bp": float(
            metrics["gene_body_log1p_1bp"]["mean_per_track_mse"]
        ),
        "exon_mse_log1p_1bp": float(
            metrics["exon_log1p_1bp"]["mean_per_track_mse"]
        ),
        "local_gradient_mse_log1p_1bp": float(
            metrics["local_gradient_log1p_1bp"]["mean_per_track_mse"]
        ),
    }
    if not all(math.isfinite(value) for value in values.values()):
        raise RuntimeError("Validation contains a non-finite summary metric")
    return values


def commands_for_job(job: dict[str, Any], spec: dict[str, Any]) -> tuple[list[str] | None, list[str]]:
    train_command: list[str] | None = None
    if not job.get("base_job"):
        training = spec["training"]
        train_command = v2_subprocess.module_command(
            "train_v2_model",
            "--model", job["model"],
            "--fold", str(job["fold"]),
            "--seed", str(job["seed"]),
            "--loss", job["loss"],
            "--max-steps", str(training["max_steps"]),
            "--sequence-length", str(training["sequence_length"]),
            "--learning-rate", str(training["learning_rate"]),
            "--weight-decay", str(training["weight_decay"]),
            "--gene-weight", str(job["gene_weight"]),
            "--max-shift-bp", str(job["max_shift_bp"]),
            "--reverse-complement-probability", str(job["reverse_complement_probability"]),
            "--hidden-channels", str(training["hidden_channels"]),
            "--output-dir", job["output_dir"],
            "--device", "cuda",
        )
        if job.get("mean_column"):
            add_option(train_command, "--mean-column", job["mean_column"])
    evaluate_command = v2_subprocess.module_command(
        "evaluate_v2_model",
        "--model", job["model"],
        "--training-loss", job["loss"],
        "--fold", str(job["fold"]),
        "--checkpoint", job["checkpoint_path"],
        "--sequence-length", str(spec["validation"]["sequence_length"]),
        "--hidden-channels", str(spec["training"]["hidden_channels"]),
        "--output", job["validation_path"],
        "--device", "cuda",
    )
    if job.get("mean_column"):
        add_option(evaluate_command, "--mean-column", job["mean_column"])
    return train_command, evaluate_command


def completed_job(job: dict[str, Any], spec_sha: str) -> dict[str, Any] | None:
    record_path = REPO_ROOT / job["job_path"]
    checkpoint_path = REPO_ROOT / job["checkpoint_path"]
    validation_path = REPO_ROOT / job["validation_path"]
    if not (record_path.is_file() and checkpoint_path.is_file() and validation_path.is_file()):
        return None
    try:
        record = read_json(record_path)
        validation = read_json(validation_path)
        if (
            record.get("status") != "completed"
            or record.get("amendment_spec_sha256") != spec_sha
            or record.get("checkpoint_sha256") != base.sha256(checkpoint_path)
            or validation.get("checkpoint_sha256") != record.get("checkpoint_sha256")
            or validation.get("chromosome_x_read") is not False
        ):
            return None
        validation_summary(validation)
    except Exception:
        return None
    record["reused"] = True
    return record


def run_job(job: dict[str, Any], gpu: int, spec: dict[str, Any], spec_sha: str) -> dict[str, Any]:
    existing = completed_job(job, spec_sha)
    if existing is not None:
        return existing
    checkpoint_path = REPO_ROOT / job["checkpoint_path"]
    validation_path = REPO_ROOT / job["validation_path"]
    record_path = REPO_ROOT / job["job_path"]
    log_path = REPO_ROOT / job["log_path"]
    record_path.parent.mkdir(parents=True, exist_ok=True)
    train_command, evaluate_command = commands_for_job(job, spec)
    started = time.monotonic()
    if train_command is not None:
        base.run_command(train_command, log_path, gpu)
    if not checkpoint_path.is_file():
        raise RuntimeError(f"Missing checkpoint after training: {checkpoint_path}")
    base.run_command(evaluate_command, log_path, gpu)
    validation = read_json(validation_path)
    summary = validation_summary(validation)
    record = {
        **job,
        "schema_version": 1,
        "status": "completed",
        "physical_gpu": gpu,
        "amendment_spec_sha256": spec_sha,
        "checkpoint_sha256": base.sha256(checkpoint_path),
        "elapsed_seconds": time.monotonic() - started,
        "reused": False,
        **summary,
    }
    base.atomic_json(record_path, record)
    return record


def formal_jobs(spec: dict[str, Any]) -> list[dict[str, Any]]:
    base_execution = read_json(BASE_EXECUTION_PATH)
    base_jobs = {
        (row["model"], row["loss"], int(row["fold"]), int(row["seed"])): row
        for row in base_execution["jobs"]
    }
    jobs = []
    for seed in spec["formal_matrix"]["seeds"]:
        for fold in spec["cv_folds"]:
            for model in spec["formal_matrix"]["models"]:
                for loss in spec["formal_matrix"]["losses"]:
                    key = (model, loss, fold, seed)
                    old = base_jobs.get(key)
                    if old:
                        output_dir = str(Path(old["checkpoint_path"]).parent)
                        checkpoint_path = old["checkpoint_path"]
                        log_path = str(Path(old["log_path"]).with_name(Path(old["log_path"]).stem + "_full_metrics.log"))
                    else:
                        root = RUN_ROOT / "formal" / model / loss / f"seed_{seed}" / f"fold_{fold}"
                        output_dir = str(root)
                        checkpoint_path = str(root / "checkpoint.pt")
                        log_path = str(LOG_ROOT / "formal" / model / loss / f"seed_{seed}" / f"fold_{fold}.log")
                    jobs.append({
                        "job_kind": "formal",
                        "job_id": f"formal:{model}:{loss}:seed{seed}:fold{fold}",
                        "model": model,
                        "loss": loss,
                        "seed": seed,
                        "fold": fold,
                        "gene_weight": spec["training"]["gene_weight"],
                        "max_shift_bp": spec["training"]["max_shift_bp"],
                        "reverse_complement_probability": spec["training"]["reverse_complement_probability"],
                        "mean_column": f"fold_{fold}_train_nonzero_mean",
                        "base_job": bool(old),
                        "output_dir": output_dir,
                        "checkpoint_path": checkpoint_path,
                        "validation_path": str(Path(output_dir) / "validation_full_metrics.json"),
                        "job_path": str(RUN_ROOT / "jobs" / "formal" / model / loss / f"seed_{seed}" / f"fold_{fold}.json"),
                        "log_path": log_path,
                    })
    return jobs


def ablation_jobs(spec: dict[str, Any]) -> list[dict[str, Any]]:
    configurations = list(spec["augmentation_ablations"])
    mean_ablation = spec["mean_scope_ablation"]
    configurations.append({
        "id": mean_ablation["id"],
        "model": mean_ablation["model"],
        "loss": mean_ablation["loss"],
        "gene_weight": mean_ablation["gene_weight"],
        "max_shift_bp": mean_ablation["max_shift_bp"],
        "reverse_complement_probability": mean_ablation["reverse_complement_probability"],
        "mean_column": mean_ablation["whole_mean_column"],
    })
    jobs = []
    seed = spec["formal_matrix"]["seeds"][0]
    for config in configurations:
        for fold in spec["cv_folds"]:
            root = RUN_ROOT / "ablations" / config["id"] / f"fold_{fold}"
            jobs.append({
                "job_kind": "ablation",
                "job_id": f"ablation:{config['id']}:fold{fold}",
                "ablation_id": config["id"],
                "model": config["model"],
                "loss": config["loss"],
                "seed": seed,
                "fold": fold,
                "gene_weight": config["gene_weight"],
                "max_shift_bp": config["max_shift_bp"],
                "reverse_complement_probability": config["reverse_complement_probability"],
                "mean_column": config.get("mean_column", f"fold_{fold}_train_nonzero_mean"),
                "base_job": False,
                "output_dir": str(root),
                "checkpoint_path": str(root / "checkpoint.pt"),
                "validation_path": str(root / "validation_full_metrics.json"),
                "job_path": str(RUN_ROOT / "jobs" / "ablations" / config["id"] / f"fold_{fold}.json"),
                "log_path": str(LOG_ROOT / "ablations" / config["id"] / f"fold_{fold}.log"),
            })
    return jobs


def run_queue(jobs: list[dict[str, Any]], gpus: list[int], spec: dict[str, Any], spec_sha: str, execution: dict[str, Any]) -> list[dict[str, Any]]:
    pending: queue.Queue[dict[str, Any]] = queue.Queue()
    for job in jobs:
        pending.put(job)
    records: list[dict[str, Any]] = []
    failures: list[dict[str, Any]] = []
    lock = threading.Lock()

    def worker(gpu: int) -> None:
        while not failures:
            try:
                job = pending.get_nowait()
            except queue.Empty:
                return
            try:
                record = run_job(job, gpu, spec, spec_sha)
            except Exception as error:
                with lock:
                    failures.append({"job_id": job["job_id"], "physical_gpu": gpu, "error": repr(error)})
                    execution["failures"] = failures
                    base.atomic_json(EXECUTION_PATH, execution)
                return
            else:
                with lock:
                    records.append(record)
                    execution["jobs"] = sorted(
                        execution.get("jobs", []) + [record], key=lambda row: row["job_id"]
                    )
                    base.atomic_json(EXECUTION_PATH, execution)
            finally:
                pending.task_done()

    threads = [threading.Thread(target=worker, args=(gpu,), daemon=False) for gpu in gpus]
    for thread in threads:
        thread.start()
    for thread in threads:
        thread.join()
    if failures:
        raise RuntimeError(f"P6B amendment job failed: {failures[0]}")
    return records


def write_tsv(path: Path, rows: list[dict[str, Any]]) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def aggregate(records: list[dict[str, Any]], keys: tuple[str, ...]) -> list[dict[str, Any]]:
    metric_names = [
        "primary_biological_score", "gene_exon_coverage_pearson_log1p",
        "per_track_pearson_128bp_log1p", "paper_loss", "log1p_mse",
        "spearman_128bp", "top1_mse_128bp", "top1_calibration_ratio_128bp",
        "gene_body_mse_log1p_1bp", "exon_mse_log1p_1bp",
        "local_gradient_mse_log1p_1bp",
    ]
    groups: dict[tuple[Any, ...], list[dict[str, Any]]] = {}
    for record in records:
        groups.setdefault(tuple(record[key] for key in keys), []).append(record)
    rows = []
    for group_key, members in sorted(groups.items()):
        row = {key: value for key, value in zip(keys, group_key, strict=True)}
        row["jobs"] = len(members)
        row["folds"] = len({member["fold"] for member in members})
        row["seeds"] = len({member["seed"] for member in members})
        for metric in metric_names:
            row[f"mean_{metric}"] = sum(float(member[metric]) for member in members) / len(members)
        rows.append(row)
    return rows


def run_development(selected: dict[str, Any], gpu: int, spec: dict[str, Any], spec_sha: str) -> dict[str, Any]:
    root = RUN_ROOT / "development_selected"
    checkpoint = REPO_ROOT / root / "checkpoint.pt"
    run_path = REPO_ROOT / root / "run.json"
    expected_seed = spec["development_retrain"]["seed"]
    if run_path.is_file() and checkpoint.is_file():
        run = read_json(run_path)
        if run.get("model") == selected["model"] and run.get("loss") == selected["loss"] and run.get("seed") == expected_seed:
            return {"model": selected["model"], "loss": selected["loss"], "seed": expected_seed, "steps": run["steps"], "checkpoint_path": str(root / "checkpoint.pt"), "checkpoint_sha256": base.sha256(checkpoint), "run_path": str(root / "run.json"), "physical_gpu": gpu}
    command = v2_subprocess.module_command(
        "train_v2_model", "--model", selected["model"], "--fold", "0",
        "--seed", str(expected_seed), "--loss", selected["loss"],
        "--max-steps", str(spec["development_retrain"]["max_steps"]),
        "--sequence-length", str(spec["training"]["sequence_length"]),
        "--learning-rate", str(spec["training"]["learning_rate"]),
        "--weight-decay", str(spec["training"]["weight_decay"]),
        "--gene-weight", str(spec["training"]["gene_weight"]),
        "--max-shift-bp", str(spec["training"]["max_shift_bp"]),
        "--reverse-complement-probability", str(spec["training"]["reverse_complement_probability"]),
        "--hidden-channels", str(spec["training"]["hidden_channels"]),
        "--output-dir", str(root), "--device", "cuda",
    )
    base.run_command(command, REPO_ROOT / LOG_ROOT / "development_selected.log", gpu)
    run = read_json(run_path)
    return {"model": selected["model"], "loss": selected["loss"], "seed": expected_seed, "steps": run["steps"], "checkpoint_path": str(root / "checkpoint.pt"), "checkpoint_sha256": base.sha256(checkpoint), "run_path": str(root / "run.json"), "physical_gpu": gpu, "amendment_spec_sha256": spec_sha}


def main() -> None:
    base.require_approval()
    spec = read_json(SPEC_PATH)
    spec_sha = base.sha256(SPEC_PATH)
    selected_gpus, snapshots = base.wait_for_gpus(spec["gpu_policy"]["minimum_free_mib"])
    formal = formal_jobs(spec)
    ablations = ablation_jobs(spec)
    execution: dict[str, Any] = {
        "schema_version": 1, "phase": "P6B", "amendment_id": spec["amendment_id"],
        "status": "running", "started_at": utc_now(), "completed_at": None,
        "amendment_spec_path": str(SPEC_PATH.relative_to(REPO_ROOT)),
        "amendment_spec_sha256": spec_sha, "selected_physical_gpus": selected_gpus,
        "gpu_wait_snapshots": snapshots, "jobs_expected": len(formal) + len(ablations),
        "jobs": [], "failures": [],
    }
    base.atomic_json(EXECUTION_PATH, execution)
    try:
        formal_records = run_queue(formal, selected_gpus, spec, spec_sha, execution)
        ablation_records = run_queue(ablations, selected_gpus, spec, spec_sha, execution)
        formal_rows = aggregate(formal_records, ("model", "loss"))
        ablation_rows = aggregate(ablation_records, ("ablation_id", "model", "loss"))
        write_tsv(RESULTS_PATH, formal_rows)
        write_tsv(ABLATION_RESULTS_PATH, ablation_rows)
        selected = max(
            formal_rows,
            key=lambda row: (
                float(row["mean_primary_biological_score"]),
                float(row["mean_gene_exon_coverage_pearson_log1p"]),
                float(row["mean_per_track_pearson_128bp_log1p"]),
                -float(row["mean_paper_loss"]),
                row["model"], row["loss"],
            ),
        )
        selection = {
            "schema_version": 1, "amendment_id": spec["amendment_id"],
            "selected_at": utc_now(), "selected": selected,
            "amendment_spec_sha256": spec_sha,
            "results_path": str(RESULTS_PATH.relative_to(REPO_ROOT)),
            "results_sha256": base.sha256(RESULTS_PATH), "chromosome_x_read": False,
            "selection_rule": spec["selection"],
        }
        base.atomic_json(SELECTION_PATH, selection)
        development = run_development(selected, selected_gpus[0], spec, spec_sha)
        if LOCK_PATH.is_file() and not OLD_LOCK_PATH.exists():
            shutil.copy2(LOCK_PATH, OLD_LOCK_PATH)
        final_lock = {
            "schema_version": 2, "model": selected["model"], "loss": selected["loss"],
            "seed": development["seed"], "steps": development["steps"],
            "checkpoint_path": development["checkpoint_path"],
            "checkpoint_sha256": development["checkpoint_sha256"],
            "selection_path": str(SELECTION_PATH.relative_to(REPO_ROOT)),
            "selection_sha256": base.sha256(SELECTION_PATH),
            "amendment_spec_sha256": spec_sha, "training_chromosomes": ["I", "II", "III", "IV", "V"],
            "legacy_chr_x_prior_exposure_disclosed": True, "locked_at": utc_now(),
            "superseded": False, "test_consumed": False, "test_consumed_at": None,
        }
        base.atomic_json(LOCK_PATH, final_lock)
        execution.update({
            "status": "completed", "completed_at": utc_now(), "formal_jobs": formal_records,
            "ablation_jobs": ablation_records, "development": development,
            "results_sha256": base.sha256(RESULTS_PATH),
            "ablation_results_sha256": base.sha256(ABLATION_RESULTS_PATH),
            "selection_sha256": base.sha256(SELECTION_PATH), "lock_sha256": base.sha256(LOCK_PATH),
            "resources_after": v2_gpu_resources.snapshot(),
        })
    except Exception:
        execution["status"] = "failed"
        execution["completed_at"] = utc_now()
        base.atomic_json(EXECUTION_PATH, execution)
        raise
    base.atomic_json(EXECUTION_PATH, execution)


if __name__ == "__main__":
    main()
