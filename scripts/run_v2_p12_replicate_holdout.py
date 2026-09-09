#!/usr/bin/env python3
"""Run the fully fresh P12 replicate-holdout training and validation matrix."""

from __future__ import annotations

import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import queue
import socket
import threading
import time
from typing import Any
import subprocess
import sys

from scripts import run_v2_p6b as base
from scripts import run_v2_p6b_amendment as shared


REPO_ROOT = Path(__file__).resolve().parents[1]
METADATA_ROOT = REPO_ROOT / "alphagenome_custom/metadata/v2"
HOLDOUT_ROOT = METADATA_ROOT / "replicate_holdout_v1"
BASE_SPEC_PATH = METADATA_ROOT / "p12_replicate_holdout_spec.json"
SPEC_PATH = METADATA_ROOT / "p12_replicate_holdout_retry_spec.json"
EXECUTION_PATH = HOLDOUT_ROOT / "p12_replicate_holdout_execution.json"
RUN_ROOT = Path("runs/v2_p12_replicate_holdout_retry_20260819")
LOG_ROOT = Path("logs/v2_p12_replicate_holdout_retry_20260819")
METRIC_NAMES = (
    "primary_biological_score",
    "gene_exon_coverage_pearson_log1p",
    "per_track_pearson_128bp_log1p",
    "paper_loss",
    "log1p_mse",
    "spearman_128bp",
    "top1_mse_128bp",
    "top1_calibration_ratio_128bp",
    "gene_body_mse_log1p_1bp",
    "exon_mse_log1p_1bp",
    "local_gradient_mse_log1p_1bp",
)


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text())


def load_spec() -> dict[str, Any]:
    retry = read_json(SPEC_PATH)
    base_spec = read_json(BASE_SPEC_PATH)
    expected_base = retry.get("base_spec_sha256")
    observed_base = sha256(BASE_SPEC_PATH)
    if expected_base != observed_base:
        raise RuntimeError(
            f"P12 retry base spec hash mismatch: expected={expected_base} observed={observed_base}"
        )
    base_spec.update(
        {
            "retry_id": retry["retry_id"],
            "supersedes_spec": retry["supersedes_spec"],
            "correction": retry["correction"],
            "base_spec_sha256": observed_base,
            "run_root": retry["run_root"],
            "log_root": retry["log_root"],
            "failed_attempt_execution": retry["failed_attempt_execution"],
        }
    )
    return base_spec


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(8 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def write_tsv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise RuntimeError(f"Refusing to write empty P12 table: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def implementation_sha256() -> dict[str, str]:
    paths = (
        "scripts/run_v2_p12_replicate_holdout.py",
        "scripts/run_v2_p11_replicate_holdout.py",
        "scripts/build_v2_replicate_holdout.py",
        "scripts/compute_v2_replicate_holdout_means.py",
        "scripts/train_v2_model.py",
        "scripts/evaluate_v2_model.py",
        "scripts/v2_bigwig_dataset.py",
        "scripts/v2_training_components.py",
    )
    return {path: sha256(REPO_ROOT / path) for path in paths}


def require_contract(spec: dict[str, Any]) -> dict[str, Any]:
    if spec.get("phase") != "P12" or spec.get("final_test_access") != "prohibited":
        raise RuntimeError("P12 spec is not a development-only locked matrix")
    matrix = spec["matrix"]
    expected = int(matrix["formal_tasks"]) + int(matrix["ablation_tasks"]) + int(matrix["component_tasks"])
    if expected != int(matrix["total_tasks"]) or expected != 165:
        raise RuntimeError("P12 task counts do not equal 165")
    p11_spec = REPO_ROOT / spec["source_p11_spec"]
    p11_execution = REPO_ROOT / spec["source_p11_execution"]
    if not p11_spec.is_file() or not p11_execution.is_file():
        raise FileNotFoundError("P12 requires completed P11 metadata and execution")
    p11 = read_json(p11_execution)
    if p11.get("phase") != "P11" or p11.get("locked_test_signal_reads") != 0:
        raise RuntimeError("P11 execution is incomplete or accessed locked test data")
    lock_path = METADATA_ROOT / "final_test_lock.json"
    report_path = METADATA_ROOT / "final_test_report.json"
    p11_spec_payload = read_json(p11_spec)
    if p11_spec_payload.get("p6c_lock_sha256") != sha256(lock_path):
        raise RuntimeError("P6C lock changed since P11 freeze")
    if p11_spec_payload.get("p6c_report_sha256") != sha256(report_path):
        raise RuntimeError("P6C report changed since P11 freeze")
    training = spec["training"]
    required = (
        training["means_path"],
        training["track_manifest"],
        training["group_manifest"],
        training["heldout_track_manifest"],
        training["heldout_group_manifest"],
        training["heldout_track_indices"],
        training["heldout_primary_track_manifest"],
        training["heldout_primary_group_manifest"],
        training["heldout_primary_track_indices"],
        training["heldout_supplementary_track_manifest"],
        training["heldout_supplementary_group_manifest"],
        training["heldout_supplementary_track_indices"],
    )
    for relative in required:
        if not (REPO_ROOT / relative).is_file():
            raise FileNotFoundError(REPO_ROOT / relative)
    with (REPO_ROOT / training["track_manifest"]).open(newline="") as handle:
        tracks = list(csv.DictReader(handle, delimiter="\t"))
    with (REPO_ROOT / training["group_manifest"]).open(newline="") as handle:
        groups = list(csv.DictReader(handle, delimiter="\t"))
    if not (len(tracks) == len(groups) == 241):
        raise RuntimeError("P12 training manifests must contain 241 rows in the same order")
    if any("test_locked" in str(path) for path in required):
        raise RuntimeError("P12 contract includes a prohibited locked-test path")
    return {"p11_execution": p11, "implementation_sha256": implementation_sha256()}


def matrix_jobs(spec: dict[str, Any]) -> list[dict[str, Any]]:
    matrix = spec["matrix"]
    training = spec["training"]
    jobs: list[dict[str, Any]] = []

    def add(configuration: dict[str, Any], fold: int, seed: int, kind: str) -> None:
        config_id = str(configuration["id"])
        root = RUN_ROOT / config_id / f"seed_{seed}" / f"fold_{fold}"
        jobs.append(
            {
                "job_id": f"p12:{kind}:{config_id}:seed{seed}:fold{fold}",
                "job_kind": kind,
                "configuration": config_id,
                "model": configuration["model"],
                "loss": configuration["loss"],
                "fold": fold,
                "seed": seed,
                "gene_weight": configuration.get("gene_weight", training["gene_weight"]),
                "max_shift_bp": configuration.get("max_shift_bp", training["max_shift_bp"]),
                "reverse_complement_probability": configuration.get(
                    "reverse_complement_probability", training["reverse_complement_probability"]
                ),
                "mean_column": configuration.get("mean_column", f"fold_{fold}_train_nonzero_mean"),
                "hidden_channels": configuration.get("hidden_channels", training["hidden_channels_default"]),
                "output_dir": str(root),
                "checkpoint_path": str(root / "checkpoint.pt"),
                "internal_validation_path": str(root / "validation_internal.json"),
        "heldout_validation_path": str(root / "validation_heldout.json"),
                "job_path": str(HOLDOUT_ROOT / "jobs" / config_id / f"seed_{seed}" / f"fold_{fold}.json"),
                "log_path": str(LOG_ROOT / config_id / f"seed_{seed}" / f"fold_{fold}.log"),
            }
        )

    for configuration in matrix["formal_configurations"]:
        for seed in matrix["seeds"]:
            for fold in matrix["folds"]:
                add(configuration, fold, seed, "formal")
    for configuration in matrix["ablation_configurations"]:
        for fold in matrix["folds"]:
            add(configuration, fold, int(configuration["seed"]), "ablation")
    for configuration in matrix["component_configurations"]:
        for seed in matrix["seeds"]:
            for fold in matrix["folds"]:
                add(configuration, fold, seed, "component")
    if len(jobs) != int(matrix["total_tasks"]):
        raise RuntimeError(f"P12 matrix generated {len(jobs)} jobs, expected 165")
    if len({job["job_id"] for job in jobs}) != len(jobs):
        raise RuntimeError("P12 matrix contains duplicate job IDs")
    return jobs


def commands_for_job(job: dict[str, Any], spec: dict[str, Any]) -> tuple[list[str], list[str], list[str]]:
    training = spec["training"]
    train = shared.v2_subprocess.module_command(
        "train_v2_model",
        "--model", job["model"], "--fold", str(job["fold"]), "--seed", str(job["seed"]),
        "--loss", job["loss"], "--max-steps", str(training["max_steps"]),
        "--sequence-length", str(training["sequence_length"]), "--learning-rate", str(training["learning_rate"]),
        "--weight-decay", str(training["weight_decay"]), "--gene-weight", str(job["gene_weight"]),
        "--max-shift-bp", str(job["max_shift_bp"]),
        "--reverse-complement-probability", str(job["reverse_complement_probability"]),
        "--hidden-channels", str(job["hidden_channels"]), "--mean-column", job["mean_column"],
        "--means-path", training["means_path"], "--track-manifest", training["track_manifest"],
        "--group-manifest", training["group_manifest"], "--output-dir", job["output_dir"], "--device", "cuda",
    )
    internal = shared.v2_subprocess.module_command(
        "evaluate_v2_model", "--model", job["model"], "--training-loss", job["loss"],
        "--fold", str(job["fold"]), "--checkpoint", job["checkpoint_path"],
        "--sequence-length", str(training["sequence_length"]), "--hidden-channels", str(job["hidden_channels"]),
        "--mean-column", job["mean_column"], "--means-path", training["means_path"],
        "--track-manifest", training["track_manifest"], "--group-manifest", training["group_manifest"],
        "--output", job["internal_validation_path"], "--device", "cuda",
    )
    heldout = shared.v2_subprocess.module_command(
        "evaluate_v2_model", "--model", job["model"], "--training-loss", job["loss"],
        "--fold", str(job["fold"]), "--checkpoint", job["checkpoint_path"],
        "--sequence-length", str(training["sequence_length"]), "--hidden-channels", str(job["hidden_channels"]),
        "--mean-column", job["mean_column"], "--means-path", training["means_path"],
        "--track-manifest", training["heldout_track_manifest"], "--group-manifest", training["heldout_group_manifest"],
        "--prediction-track-indices", training["heldout_track_indices"],
        "--valid-intervals-path", f"alphagenome_custom/intervals/v2/fold_{job['fold']}/valid.tsv",
        "--output", job["heldout_validation_path"], "--device", "cuda",
    )
    serialized = " ".join([*train, *internal, *heldout])
    if "--final-test" in serialized or "test_locked.tsv" in serialized:
        raise RuntimeError("P12 command unexpectedly requests final-test access")
    return train, internal, heldout


def record_job(job: dict[str, Any], spec_sha: str, gpu: int, elapsed: float) -> dict[str, Any]:
    checkpoint = REPO_ROOT / job["checkpoint_path"]
    run_path = REPO_ROOT / job["output_dir"] / "run.json"
    internal_path = REPO_ROOT / job["internal_validation_path"]
    heldout_path = REPO_ROOT / job["heldout_validation_path"]
    if not all(path.is_file() for path in (checkpoint, run_path, internal_path, heldout_path)):
        raise RuntimeError(f"P12 output missing for {job['job_id']}")
    run = read_json(run_path)
    internal = read_json(internal_path)
    heldout = read_json(heldout_path)
    if run.get("model") != job["model"] or run.get("seed") != job["seed"] or run.get("fold") != job["fold"]:
        raise RuntimeError(f"P12 run contract mismatch for {job['job_id']}")
    for validation in (internal, heldout):
        if validation.get("locked_test_block_signal_reads") is not False:
            raise RuntimeError(f"P12 validation accessed locked test for {job['job_id']}")
    record = {
        **job,
        "schema_version": 1,
        "status": "completed",
        "spec_sha256": spec_sha,
        "physical_gpu": gpu,
        "elapsed_seconds": elapsed,
        "checkpoint_sha256": sha256(checkpoint),
        "run_sha256": sha256(run_path),
        "internal_validation_sha256": sha256(internal_path),
        "heldout_validation_sha256": sha256(heldout_path),
        "trainable_parameters": int(run["trainable_parameters"]),
        "internal_metrics": shared.validation_summary(internal),
        "heldout_metrics": shared.validation_summary(heldout),
    }
    return record


def validate_resume_records(
    jobs: list[dict[str, Any]], execution: dict[str, Any], spec_sha: str
) -> list[dict[str, Any]]:
    """Accept only completed records whose four immutable outputs still hash-match."""

    by_id = {job["job_id"]: job for job in jobs}
    records: list[dict[str, Any]] = []
    seen: set[str] = set()
    for record in execution.get("jobs", []):
        job_id = record.get("job_id")
        if job_id in seen:
            raise RuntimeError(f"P12 resume has duplicate job record: {job_id}")
        if job_id not in by_id:
            raise RuntimeError(f"P12 resume has an unregistered job record: {job_id}")
        if record.get("status") != "completed" or record.get("spec_sha256") != spec_sha:
            raise RuntimeError(f"P12 resume record is not complete under the locked spec: {job_id}")
        job = by_id[job_id]
        checkpoint = REPO_ROOT / job["checkpoint_path"]
        run_path = REPO_ROOT / job["output_dir"] / "run.json"
        internal_path = REPO_ROOT / job["internal_validation_path"]
        heldout_path = REPO_ROOT / job["heldout_validation_path"]
        paths = (checkpoint, run_path, internal_path, heldout_path)
        if not all(path.is_file() for path in paths):
            raise RuntimeError(f"P12 resume output is missing: {job_id}")
        recorded_hashes = {
            checkpoint: record.get("checkpoint_sha256"),
            run_path: record.get("run_sha256"),
            internal_path: record.get("internal_validation_sha256"),
            heldout_path: record.get("heldout_validation_sha256"),
        }
        for path, expected in recorded_hashes.items():
            observed = sha256(path)
            if expected != observed:
                raise RuntimeError(
                    f"P12 resume hash mismatch for {job_id}: {path} expected={expected} observed={observed}"
                )
        for validation_path in (internal_path, heldout_path):
            validation = read_json(validation_path)
            if validation.get("locked_test_block_signal_reads") is not False:
                raise RuntimeError(f"P12 resume validation accessed locked test: {job_id}")
        records.append(record)
        seen.add(job_id)
    return records


def run_queue(
    jobs: list[dict[str, Any]],
    gpus: list[int],
    spec: dict[str, Any],
    spec_sha: str,
    execution: dict[str, Any],
    initial_records: list[dict[str, Any]] | None = None,
) -> list[dict[str, Any]]:
    pending: queue.Queue[dict[str, Any]] = queue.Queue()
    for job in jobs:
        pending.put(job)
    records: list[dict[str, Any]] = list(initial_records or [])
    failures: list[dict[str, Any]] = []
    lock = threading.Lock()

    def worker(gpu: int) -> None:
        while not failures:
            try:
                job = pending.get_nowait()
            except queue.Empty:
                return
            started = time.monotonic()
            try:
                train, internal, heldout = commands_for_job(job, spec)
                log_path = REPO_ROOT / job["log_path"]
                base.run_command(train, log_path, gpu)
                base.run_command(internal, log_path, gpu)
                base.run_command(heldout, log_path, gpu)
                record = record_job(job, spec_sha, gpu, time.monotonic() - started)
            except Exception as error:
                with lock:
                    failures.append({"job_id": job["job_id"], "physical_gpu": gpu, "error": repr(error)})
                    execution["failures"] = failures
                    base.atomic_json(EXECUTION_PATH, execution)
                return
            else:
                with lock:
                    records.append(record)
                    execution["jobs"] = sorted(execution.get("jobs", []) + [record], key=lambda row: row["job_id"])
                    execution["registered_tasks_completed"] = len(execution["jobs"])
                    base.atomic_json(EXECUTION_PATH, execution)
            finally:
                pending.task_done()

    threads = [threading.Thread(target=worker, args=(gpu,), daemon=False) for gpu in gpus]
    for thread in threads:
        thread.start()
    for thread in threads:
        thread.join()
    if failures:
        raise RuntimeError(f"P12 job failed: {failures[0]}")
    return records


def aggregate_records(records: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    by_config: dict[str, list[dict[str, Any]]] = {}
    for record in records:
        by_config.setdefault(record["configuration"], []).append(record)
    for configuration, members in sorted(by_config.items()):
        row = {
            "configuration": configuration,
            "model": members[0]["model"],
            "loss": members[0]["loss"],
            "job_kind": members[0]["job_kind"],
            "jobs": len(members),
            "folds": len({member["fold"] for member in members}),
            "seeds": len({member["seed"] for member in members}),
        }
        for scope in ("internal", "heldout"):
            for metric in METRIC_NAMES:
                row[f"mean_{scope}_{metric}"] = sum(float(member[f"{scope}_metrics"][metric]) for member in members) / len(members)
        rows.append(row)
    return rows


def main() -> None:
    spec = load_spec()
    preflight = require_contract(spec)
    jobs = matrix_jobs(spec)
    spec_sha = sha256(SPEC_PATH)
    selected_gpus, snapshots = base.wait_for_gpus(spec["gpu_policy"]["minimum_free_mib"])
    allowed = set(spec["gpu_policy"]["allowed_physical_indices"])
    if not selected_gpus or not set(selected_gpus) <= allowed:
        raise RuntimeError(f"P12 selected disallowed physical GPUs: {selected_gpus}")
    fresh_execution: dict[str, Any] = {
        "schema_version": 1,
        "phase": "P12",
        "contract": "replicate_holdout_v1_model_matrix",
        "status": "running",
        "started_at": utc_now(),
        "completed_at": None,
        "hostname": socket.gethostname(),
        "working_directory": str(REPO_ROOT),
        "spec_path": str(SPEC_PATH.relative_to(REPO_ROOT)),
        "spec_sha256": spec_sha,
        "implementation_sha256": preflight["implementation_sha256"],
        "selected_physical_gpus": selected_gpus,
        "gpu_wait_snapshots": snapshots,
        "registered_tasks": len(jobs),
        "jobs": [],
        "failures": [],
        "locked_test_block_signal_reads": 0,
        "final_test_access": "prohibited",
    }
    existing_execution: dict[str, Any] | None = None
    existing_records: list[dict[str, Any]] = []
    if EXECUTION_PATH.is_file():
        candidate = read_json(EXECUTION_PATH)
        if (
            candidate.get("phase") == "P12"
            and candidate.get("spec_path") == str(SPEC_PATH.relative_to(REPO_ROOT))
            and candidate.get("spec_sha256") == spec_sha
            and int(candidate.get("registered_tasks", 0)) == len(jobs)
            and candidate.get("final_test_access") == "prohibited"
            and not candidate.get("failures")
        ):
            existing_execution = candidate
            existing_records = validate_resume_records(jobs, candidate, spec_sha)
    execution = existing_execution or fresh_execution
    if existing_execution is not None:
        completed_ids = {record["job_id"] for record in existing_records}
        jobs = [job for job in jobs if job["job_id"] not in completed_ids]
        execution.update(
            {
                "status": "running",
                "completed_at": None,
                "resume_count": int(execution.get("resume_count", 0)) + 1,
                "resumed_at": utc_now(),
                "resume_implementation_sha256": preflight["implementation_sha256"],
                "selected_physical_gpus": selected_gpus,
                "gpu_wait_snapshots": snapshots,
            }
        )
    else:
        base.atomic_json(EXECUTION_PATH, execution)
    try:
        records = run_queue(jobs, selected_gpus, spec, spec_sha, execution, existing_records)
        records.sort(key=lambda row: row["job_id"])
        write_tsv(HOLDOUT_ROOT / "p12_replicate_holdout_results.tsv", aggregate_records(records))
        metric_rows = []
        for record in records:
            for scope in ("internal", "heldout"):
                for metric in METRIC_NAMES:
                    metric_rows.append(
                        {
                            "job_id": record["job_id"],
                            "configuration": record["configuration"],
                            "job_kind": record["job_kind"],
                            "model": record["model"],
                            "loss": record["loss"],
                            "fold": record["fold"],
                            "seed": record["seed"],
                            "scope": scope,
                            "metric": metric,
                            "value": record[f"{scope}_metrics"][metric],
                        }
                    )
        write_tsv(HOLDOUT_ROOT / "p12_replicate_holdout_metrics.tsv", metric_rows)
        write_tsv(
            HOLDOUT_ROOT / "p12_replicate_holdout_job_records.tsv",
            [{key: value for key, value in record.items() if key not in {"internal_metrics", "heldout_metrics"}} for record in records],
        )
        role_analysis = subprocess.run(
            [sys.executable, "-m", "scripts.analyze_v2_p12_heldout_roles"],
            cwd=REPO_ROOT,
            check=False,
        )
        if role_analysis.returncode != 0:
            raise RuntimeError("P12 heldout role analysis failed")
        paired_analysis = subprocess.run(
            [sys.executable, "-m", "scripts.analyze_v2_p12_replicate_holdout"],
            cwd=REPO_ROOT,
            check=False,
        )
        if paired_analysis.returncode != 0:
            raise RuntimeError("P12 paired analysis failed")
        execution.update({
            "status": "completed",
            "completed_at": utc_now(),
            "registered_tasks_completed": len(records),
            "locked_test_block_signal_reads": 0,
            "result_sha256": sha256(HOLDOUT_ROOT / "p12_replicate_holdout_results.tsv"),
            "metric_records_sha256": sha256(HOLDOUT_ROOT / "p12_replicate_holdout_metrics.tsv"),
            "heldout_role_metrics_sha256": sha256(HOLDOUT_ROOT / "p12_heldout_role_metrics.tsv"),
            "paired_effects_sha256": sha256(HOLDOUT_ROOT / "p12_replicate_holdout_paired_effects.tsv"),
            "job_records_sha256": sha256(HOLDOUT_ROOT / "p12_replicate_holdout_job_records.tsv"),
        })
        base.atomic_json(EXECUTION_PATH, execution)
    except Exception:
        execution["status"] = "failed"
        execution["completed_at"] = utc_now()
        base.atomic_json(EXECUTION_PATH, execution)
        raise
    print(json.dumps(execution, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
