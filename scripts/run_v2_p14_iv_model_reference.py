#!/usr/bin/env python3
"""Run the preregistered P14 frozen Model B I--V evaluation matrix.

This is inference only. It reuses P12 checkpoints without training, evaluates
only P13's I--V-only intervals and records one immutable result per fold/seed.
"""

from __future__ import annotations

import csv
import hashlib
import json
import queue
import socket
import threading
import time
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from scripts import run_v2_p6b as base


REPO_ROOT = Path(__file__).resolve().parents[1]
METADATA_ROOT = REPO_ROOT / "alphagenome_custom" / "metadata" / "v2"
SPEC_PATH = METADATA_ROOT / "p14_iv_frozen_model_reference_spec.json"
EXECUTION_PATH = METADATA_ROOT / "p14_iv_frozen_model_reference_execution.json"
STATE_PATH = METADATA_ROOT / "execution_state.json"
RESULT_ROOT = REPO_ROOT / "results/v2_p14_iv_model_reference"
LOG_ROOT = REPO_ROOT / "logs/v2_p14_iv_model_reference"


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    temporary.replace(path)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    if not rows:
        raise RuntimeError(f"Refusing to write empty P14 table: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def load_spec() -> dict[str, Any]:
    if not SPEC_PATH.is_file():
        raise FileNotFoundError(f"P14 requires preregistration: {SPEC_PATH}")
    spec = json.loads(SPEC_PATH.read_text())
    if spec.get("phase") != "P14" or spec.get("status") != "preregistered":
        raise RuntimeError("P14 specification is invalid")
    if spec.get("evaluation", {}).get("allowed_chromosomes") != ["I", "II", "III", "IV", "V"]:
        raise RuntimeError("P14 specification does not enforce I--V-only evaluation")
    jobs = spec.get("jobs", [])
    if len(jobs) != 15 or len({job["job_id"] for job in jobs}) != 15:
        raise RuntimeError("P14 requires exactly 15 unique B_paper inference jobs")
    for job in jobs:
        checkpoint = REPO_ROOT / job["checkpoint_path"]
        interval = REPO_ROOT / job["interval_path"]
        if not checkpoint.is_file() or sha256(checkpoint) != job["checkpoint_sha256"]:
            raise RuntimeError(f"P14 checkpoint hash mismatch: {job['job_id']}")
        if not interval.is_file() or sha256(interval) != job["interval_sha256"]:
            raise RuntimeError(f"P14 interval hash mismatch: {job['job_id']}")
        rows = read_tsv(interval)
        if not rows or any(row["chromosome"] == "X" or row["role"] == "test_locked" for row in rows):
            raise RuntimeError(f"P14 prohibited interval row: {job['job_id']}")
    return spec


def require_controller_scope() -> None:
    if not STATE_PATH.is_file():
        raise FileNotFoundError(STATE_PATH)
    state = json.loads(STATE_PATH.read_text())
    approval = state.get("approvals", {}).get("G4_gpu_experiments", {})
    if not (
        state.get("current_phase") == "P14"
        and state.get("status") == "RUNNING"
        and approval.get("approved") is True
        and approval.get("scope") == "p14_iv_frozen_model_reference"
    ):
        raise RuntimeError("P14 GPU inference requires controller state P14/RUNNING and G4:p14_iv_frozen_model_reference")


def output_path(job: dict[str, Any]) -> Path:
    return RESULT_ROOT / "evaluations" / f"seed_{job['seed']}" / f"fold_{job['fold']}.json"


def log_path(job: dict[str, Any]) -> Path:
    return LOG_ROOT / f"seed_{job['seed']}" / f"fold_{job['fold']}.log"


def command_for(job: dict[str, Any], spec: dict[str, Any]) -> list[str]:
    evaluation = spec["evaluation"]
    return [
        str(Path(sys.executable).resolve()), "-m", "scripts.evaluate_v2_model",
        "--model", job["model"], "--training-loss", job["training_loss"],
        "--fold", str(job["fold"]), "--checkpoint", job["checkpoint_path"],
        "--sequence-length", "131072", "--hidden-channels", str(job["hidden_channels"]),
        "--mean-column", job["mean_column"], "--means-path", evaluation["means_path"],
        "--track-manifest", evaluation["track_manifest"], "--group-manifest", evaluation["group_manifest"],
        "--prediction-track-indices", evaluation["prediction_track_indices"],
        "--valid-intervals-path", job["interval_path"],
        "--output", str(output_path(job).relative_to(REPO_ROOT)), "--device", "cuda",
    ]


def load_completed_record(job: dict[str, Any], spec_sha: str) -> dict[str, Any] | None:
    path = output_path(job)
    if not path.is_file():
        return None
    payload = json.loads(path.read_text())
    required_chromosomes = ["I", "II", "III", "IV", "V"]
    if not (
        payload.get("model") == "B"
        and payload.get("training_loss") == "paper"
        and int(payload.get("fold", -1)) == int(job["fold"])
        and int(payload.get("seed", -1)) == int(job["seed"])
        and payload.get("checkpoint_sha256") == job["checkpoint_sha256"]
        and payload.get("intervals_sha256") == job["interval_sha256"]
        and payload.get("evaluated_chromosomes") == required_chromosomes
        and payload.get("locked_test_block_signal_reads") is False
        and float(payload.get("validation_core_coverage_fraction", 0.0)) == 1.0
    ):
        return None
    primary = payload["full_metrics"]["primary"]
    return {
        "job_id": job["job_id"], "source_p12_job_id": job["source_p12_job_id"],
        "fold": int(job["fold"]), "seed": int(job["seed"]), "physical_gpu": None,
        "checkpoint_path": job["checkpoint_path"], "checkpoint_sha256": job["checkpoint_sha256"],
        "result_path": str(path.relative_to(REPO_ROOT)), "result_sha256": sha256(path),
        "primary_biological_score": 0.5 * (
            float(primary["mean_per_track_gene_exon_coverage_pearson_log1p"])
            + float(primary["mean_per_track_pearson_128bp_log1p"])
        ),
        "gene_exon_pearson_log1p": float(primary["mean_per_track_gene_exon_coverage_pearson_log1p"]),
        "pearson_128bp_log1p": float(primary["mean_per_track_pearson_128bp_log1p"]),
        "spec_sha256": spec_sha,
    }


def run_queue(jobs: list[dict[str, Any]], spec: dict[str, Any], spec_sha: str, execution: dict[str, Any], gpus: list[int]) -> list[dict[str, Any]]:
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
            started = time.monotonic()
            try:
                base.run_command(command_for(job, spec), log_path(job), gpu)
                record = load_completed_record(job, spec_sha)
                if record is None:
                    raise RuntimeError(f"P14 evaluation contract did not validate: {job['job_id']}")
                record["physical_gpu"] = gpu
                record["elapsed_seconds"] = time.monotonic() - started
            except Exception as error:
                with lock:
                    failures.append({"job_id": job["job_id"], "physical_gpu": gpu, "error": repr(error)})
                    execution["failures"] = failures
                    atomic_json(EXECUTION_PATH, execution)
                return
            else:
                with lock:
                    records.append(record)
                    execution["jobs"] = sorted(execution.get("jobs", []) + [record], key=lambda row: row["job_id"])
                    execution["registered_tasks_completed"] = len(execution["jobs"])
                    atomic_json(EXECUTION_PATH, execution)
            finally:
                pending.task_done()

    threads = [threading.Thread(target=worker, args=(gpu,), daemon=False) for gpu in gpus]
    for thread in threads:
        thread.start()
    for thread in threads:
        thread.join()
    if failures:
        raise RuntimeError(f"P14 job failed: {failures[0]}")
    return records


def summarize(records: list[dict[str, Any]], spec: dict[str, Any]) -> None:
    if len(records) != 15:
        raise RuntimeError(f"P14 requires 15 completed records, observed {len(records)}")
    metrics = ("primary_biological_score", "gene_exon_pearson_log1p", "pearson_128bp_log1p")
    by_fold: list[dict[str, Any]] = []
    for fold in range(1, 6):
        members = [row for row in records if row["fold"] == fold]
        if len(members) != 3:
            raise RuntimeError(f"P14 fold {fold} has {len(members)} records")
        by_fold.append({"fold": fold, **{metric: sum(float(row[metric]) for row in members) / len(members) for metric in metrics}})
    summary = []
    for metric in metrics:
        values = [float(row[metric]) for row in by_fold]
        mean = sum(values) / len(values)
        sd = (sum((value - mean) ** 2 for value in values) / (len(values) - 1)) ** 0.5
        summary.append({"metric": metric, "equal_weight_fold_mean": mean, "fold_sd": sd, "fold_min": min(values), "fold_max": max(values), "primary_unit": "five predefined genomic folds"})
    write_tsv(RESULT_ROOT / "p14_iv_model_reference_records.tsv", sorted(records, key=lambda row: (row["fold"], row["seed"])), list(records[0]))
    write_tsv(RESULT_ROOT / "p14_iv_model_reference_by_fold.tsv", by_fold, list(by_fold[0]))
    write_tsv(RESULT_ROOT / "p14_iv_model_reference_summary.tsv", summary, list(summary[0]))

    p13 = {
        int(row["fold"]): row
        for row in read_tsv(REPO_ROOT / "results/v2_p13_replicate_agreement_iv/p13_replicate_agreement_by_fold.tsv")
    }
    paired = []
    for row in by_fold:
        replicate = float(p13[row["fold"]]["primary_biological_score"])
        paired.append({"fold": row["fold"], "model_b_primary_biological_score": row["primary_biological_score"], "training_aggregate_to_heldout_reference": replicate, "model_minus_replicate_reference": row["primary_biological_score"] - replicate, "interpretation": "Same I--V-only held-out biological-unit interval contract; reference is not an absolute performance ceiling or external validation."})
    write_tsv(RESULT_ROOT / "p14_iv_model_vs_replicate_reference_by_fold.tsv", paired, list(paired[0]))


def main() -> None:
    require_controller_scope()
    spec = load_spec()
    spec_sha = sha256(SPEC_PATH)
    policy = spec["gpu_policy"]
    selected_gpus, snapshots = base.wait_for_gpus(int(policy["minimum_free_mib"]))
    if not selected_gpus or not set(selected_gpus) <= {2, 3}:
        raise RuntimeError(f"P14 selected disallowed GPUs: {selected_gpus}")
    existing = json.loads(EXECUTION_PATH.read_text()) if EXECUTION_PATH.is_file() else None
    records = []
    if existing:
        if existing.get("phase") != "P14" or existing.get("spec_sha256") != spec_sha or existing.get("failures"):
            raise RuntimeError("P14 existing execution cannot be resumed")
        for job in spec["jobs"]:
            record = load_completed_record(job, spec_sha)
            if record is not None:
                records.append(record)
    execution: dict[str, Any] = existing or {
        "schema_version": 1, "phase": "P14", "status": "running", "started_at": utc_now(),
        "hostname": socket.gethostname(), "working_directory": str(REPO_ROOT),
        "spec_path": str(SPEC_PATH.relative_to(REPO_ROOT)), "spec_sha256": spec_sha,
        "registered_tasks": 15, "jobs": [], "failures": [],
        "locked_test_block_signal_reads": 0, "final_test_access": "prohibited",
    }
    execution.update({"status": "running", "selected_physical_gpus": selected_gpus, "gpu_wait_snapshots": snapshots, "resumed_at": utc_now() if existing else None, "registered_tasks_completed": len(records)})
    execution["jobs"] = sorted(records, key=lambda row: row["job_id"])
    atomic_json(EXECUTION_PATH, execution)
    completed_ids = {record["job_id"] for record in records}
    pending = [job for job in spec["jobs"] if job["job_id"] not in completed_ids]
    try:
        new_records = run_queue(pending, spec, spec_sha, execution, selected_gpus)
        records = sorted(records + new_records, key=lambda row: row["job_id"])
        summarize(records, spec)
        execution.update({"status": "completed", "completed_at": utc_now(), "registered_tasks_completed": len(records), "jobs": records})
        atomic_json(EXECUTION_PATH, execution)
    except Exception:
        execution.update({"status": "failed", "completed_at": utc_now()})
        atomic_json(EXECUTION_PATH, execution)
        raise
    print(json.dumps({"status": "completed", "records": len(records), "output": str(RESULT_ROOT.relative_to(REPO_ROOT))}, indent=2))


if __name__ == "__main__":
    main()
