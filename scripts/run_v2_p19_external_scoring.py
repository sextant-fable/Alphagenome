#!/usr/bin/env python3
"""Run P19 external represented-head inference on physical GPUs 2 and 3."""

from __future__ import annotations

import csv
import hashlib
import json
import queue
import socket
import sys
import threading
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from scripts import run_v2_p6b as base


REPO_ROOT = Path(__file__).resolve().parents[1]
METADATA_ROOT = REPO_ROOT / "alphagenome_custom/metadata/v2"
SPEC_PATH = METADATA_ROOT / "p19_external_represented_head_scoring_spec.json"
EXECUTION_PATH = METADATA_ROOT / "p19_external_represented_head_scoring_execution.json"
STATE_PATH = METADATA_ROOT / "execution_state.json"
RESULT_ROOT = REPO_ROOT / "results/v2_p19_external_scoring"
LOG_ROOT = REPO_ROOT / "logs/v2_p19_external_scoring"
ALLOWED_CHROMOSOMES = ["I", "II", "III", "IV", "V"]


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


def require_controller_scope() -> None:
    state = json.loads(STATE_PATH.read_text())
    approval = state.get("approvals", {}).get("G4_gpu_experiments", {})
    if not (
        state.get("current_phase") == "P19"
        and state.get("status") == "RUNNING"
        and approval.get("approved") is True
        and approval.get("scope") == "p19_external_represented_head_scoring"
    ):
        raise RuntimeError("P19 requires controller state P19/RUNNING and G4:p19_external_represented_head_scoring")


def load_spec() -> dict[str, Any]:
    if not SPEC_PATH.is_file():
        raise FileNotFoundError(SPEC_PATH)
    spec = json.loads(SPEC_PATH.read_text())
    if spec.get("phase") != "P19" or spec.get("status") != "preregistered":
        raise RuntimeError("P19 specification is invalid")
    if spec.get("allowed_chromosomes") != ALLOWED_CHROMOSOMES or spec.get("final_test_access") != "prohibited":
        raise RuntimeError("P19 specification violates the I--V-only boundary")
    jobs = spec.get("jobs", [])
    if len(jobs) != 60 or len({job["job_id"] for job in jobs}) != 60:
        raise RuntimeError("P19 requires exactly 60 unique source-by-fold-by-seed jobs")
    for job in jobs:
        checkpoint = REPO_ROOT / job["checkpoint_path"]
        interval = REPO_ROOT / job["interval_path"]
        track = REPO_ROOT / job["track_manifest"]
        group = REPO_ROOT / job["group_manifest"]
        indices = REPO_ROOT / job["prediction_track_indices"]
        if not checkpoint.is_file() or sha256(checkpoint) != job["checkpoint_sha256"]:
            raise RuntimeError(f"P19 checkpoint hash mismatch: {job['job_id']}")
        if not interval.is_file() or sha256(interval) != job["interval_sha256"]:
            raise RuntimeError(f"P19 interval hash mismatch: {job['job_id']}")
        if not all(path.is_file() for path in (track, group, indices)):
            raise RuntimeError(f"P19 manifest missing: {job['job_id']}")
        interval_rows = read_tsv(interval)
        if not interval_rows or any(row.get("chromosome") not in ALLOWED_CHROMOSOMES for row in interval_rows) or any(row.get("role") == "test_locked" for row in interval_rows):
            raise RuntimeError(f"P19 prohibited interval row: {job['job_id']}")
        bw = REPO_ROOT / read_tsv(track)[0]["output_path"]
        if not bw.is_file():
            raise RuntimeError(f"P19 external BigWig missing: {job['job_id']}")
    return spec


def output_path(job: dict[str, Any]) -> Path:
    return RESULT_ROOT / "evaluations" / job["run_accession"] / f"seed_{job['seed']}" / f"fold_{job['fold']}.json"


def log_path(job: dict[str, Any]) -> Path:
    return LOG_ROOT / job["run_accession"] / f"seed_{job['seed']}" / f"fold_{job['fold']}.log"


def command_for(job: dict[str, Any]) -> list[str]:
    return [
        str(Path(sys.executable).resolve()), "-m", "scripts.evaluate_v2_model",
        "--model", "B", "--training-loss", "paper", "--fold", str(job["fold"]),
        "--checkpoint", job["checkpoint_path"], "--sequence-length", "131072",
        "--hidden-channels", "64", "--mean-column", job["mean_column"],
        "--means-path", "results/v2_p16_iv_training_normalization/track_nonzero_means.tsv",
        "--track-manifest", job["track_manifest"], "--group-manifest", job["group_manifest"],
        "--prediction-track-indices", job["prediction_track_indices"],
        "--valid-intervals-path", job["interval_path"], "--output", str(output_path(job).relative_to(REPO_ROOT)),
        "--device", "cuda",
    ]


def load_completed_record(job: dict[str, Any], spec_sha: str) -> dict[str, Any] | None:
    path = output_path(job)
    if not path.is_file():
        return None
    try:
        payload = json.loads(path.read_text())
        primary = payload["full_metrics"]["primary"]
        metrics = {
            "primary_biological_score": 0.5 * (float(primary["mean_per_track_gene_exon_coverage_pearson_log1p"]) + float(primary["mean_per_track_pearson_128bp_log1p"])),
            "gene_exon_pearson_log1p": float(primary["mean_per_track_gene_exon_coverage_pearson_log1p"]),
            "pearson_128bp_log1p": float(primary["mean_per_track_pearson_128bp_log1p"]),
        }
        if not (
            payload.get("model") == "B" and payload.get("training_loss") == "paper"
            and int(payload.get("fold", -1)) == int(job["fold"])
            and int(payload.get("seed", -1)) == int(job["seed"])
            and payload.get("checkpoint_sha256") == job["checkpoint_sha256"]
            and payload.get("intervals_sha256") == job["interval_sha256"]
            and payload.get("evaluated_chromosomes") == ALLOWED_CHROMOSOMES
            and payload.get("locked_test_block_signal_reads") is False
            and float(payload.get("validation_core_coverage_fraction", 0.0)) == 1.0
            and all(value == value and abs(value) < float("inf") for value in metrics.values())
        ):
            return None
    except Exception:
        return None
    return {
        "job_id": job["job_id"], "run_accession": job["run_accession"], "study_accession": job["study_accession"],
        "source_p17_job_id": job["source_p17_job_id"], "fold": int(job["fold"]), "seed": int(job["seed"]),
        "physical_gpu": None, "checkpoint_path": job["checkpoint_path"], "checkpoint_sha256": job["checkpoint_sha256"],
        "result_path": str(path.relative_to(REPO_ROOT)), "result_sha256": sha256(path), "spec_sha256": spec_sha,
        "locked_test_block_signal_reads": False, "evaluated_chromosomes": ",".join(ALLOWED_CHROMOSOMES), **metrics,
    }


def run_queue(jobs: list[dict[str, Any]], spec_sha: str, execution: dict[str, Any], gpus: list[int]) -> list[dict[str, Any]]:
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
                base.run_command(command_for(job), log_path(job), gpu)
                record = load_completed_record(job, spec_sha)
                if record is None:
                    raise RuntimeError(f"P19 evaluation contract did not validate: {job['job_id']}")
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
        raise RuntimeError(f"P19 job failed: {failures[0]}")
    return records


def write_summary(records: list[dict[str, Any]]) -> None:
    fields = ["run_accession", "study_accession", "fold", "seed", "primary_biological_score", "gene_exon_pearson_log1p", "pearson_128bp_log1p"]
    rows = [{key: row[key] for key in fields} for row in sorted(records, key=lambda row: (row["run_accession"], row["fold"], row["seed"]))]
    RESULT_ROOT.mkdir(parents=True, exist_ok=True)
    with (RESULT_ROOT / "p19_external_records.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)
    by_source = []
    for run in sorted({row["run_accession"] for row in records}):
        members = [row for row in records if row["run_accession"] == run]
        by_source.append({"run_accession": run, "study_accession": members[0]["study_accession"], "n_jobs": len(members), "primary_mean": sum(float(row["primary_biological_score"]) for row in members) / len(members), "gene_exon_mean": sum(float(row["gene_exon_pearson_log1p"]) for row in members) / len(members), "pearson_128bp_mean": sum(float(row["pearson_128bp_log1p"]) for row in members) / len(members), "primary_unit": "15 frozen P17 checkpoints"})
    with (RESULT_ROOT / "p19_external_by_source.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(by_source[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(by_source)


def main() -> None:
    require_controller_scope()
    spec = load_spec()
    spec_sha = sha256(SPEC_PATH)
    selected_gpus, snapshots = base.wait_for_gpus(70000)
    if not selected_gpus or not set(selected_gpus) <= {2, 3}:
        raise RuntimeError(f"P19 selected disallowed GPUs: {selected_gpus}")
    existing = json.loads(EXECUTION_PATH.read_text()) if EXECUTION_PATH.is_file() else None
    if existing and (existing.get("phase") != "P19" or existing.get("spec_sha256") != spec_sha or existing.get("failures")):
        raise RuntimeError("P19 existing execution cannot be resumed")
    execution: dict[str, Any] = existing or {"schema_version": 1, "phase": "P19", "status": "running", "started_at": utc_now(), "hostname": socket.gethostname(), "working_directory": str(REPO_ROOT), "spec_path": str(SPEC_PATH.relative_to(REPO_ROOT)), "spec_sha256": spec_sha, "registered_tasks": 60, "jobs": [], "failures": [], "locked_test_block_signal_reads": 0, "final_test_access": "prohibited", "model_inference": "external_represented_head_inference"}
    records = []
    for job in spec["jobs"]:
        record = load_completed_record(job, spec_sha)
        if record is not None:
            records.append(record)
    execution.update({"status": "running", "selected_physical_gpus": selected_gpus, "gpu_wait_snapshots": snapshots, "resumed_at": utc_now() if existing else None, "registered_tasks_completed": len(records), "jobs": sorted(records, key=lambda row: row["job_id"])})
    atomic_json(EXECUTION_PATH, execution)
    completed_ids = {row["job_id"] for row in records}
    pending = [job for job in spec["jobs"] if job["job_id"] not in completed_ids]
    try:
        records = sorted(records + run_queue(pending, spec_sha, execution, selected_gpus), key=lambda row: row["job_id"])
        if len(records) != 60:
            raise RuntimeError(f"P19 requires 60 completed jobs, observed {len(records)}")
        write_summary(records)
        execution.update({"status": "completed", "completed_at": utc_now(), "registered_tasks_completed": len(records), "jobs": records, "records_sha256": sha256(RESULT_ROOT / "p19_external_records.tsv"), "by_source_sha256": sha256(RESULT_ROOT / "p19_external_by_source.tsv")})
        atomic_json(EXECUTION_PATH, execution)
    except Exception:
        execution.update({"status": "failed", "completed_at": utc_now()})
        atomic_json(EXECUTION_PATH, execution)
        raise
    print(json.dumps({"status": "completed", "phase": "P19", "records": len(records), "output": str(RESULT_ROOT.relative_to(REPO_ROOT))}, indent=2))


if __name__ == "__main__":
    main()
