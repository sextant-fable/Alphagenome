#!/usr/bin/env python3
"""Run P17's bounded strict I--V ablation and Basenji2-style matrix."""

from __future__ import annotations

import csv
import hashlib
import json
import queue
import sys
import socket
import threading
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from scripts import run_v2_p6b as base
from scripts import run_v2_p6b_amendment as shared


REPO_ROOT = Path(__file__).resolve().parents[1]
METADATA_ROOT = REPO_ROOT / "alphagenome_custom" / "metadata" / "v2"
SPEC_PATH = METADATA_ROOT / "p17_iv_controlled_matrix_spec.json"
EXECUTION_PATH = METADATA_ROOT / "p17_iv_controlled_matrix_execution.json"
STATE_PATH = METADATA_ROOT / "execution_state.json"
P15_ROOT = REPO_ROOT / "results/v2_p15_iv_training_interval_contract"
RESULT_ROOT = REPO_ROOT / "results/v2_p17_iv_controlled_matrix"
RUN_ROOT = Path("runs/v2_p17_iv_controlled_matrix")
LOG_ROOT = Path("logs/v2_p17_iv_controlled_matrix")
ALLOWED_CHROMOSOMES = ["I", "II", "III", "IV", "V"]
PRIMARY_METRICS = (
    "primary_biological_score",
    "gene_exon_coverage_pearson_log1p",
    "per_track_pearson_128bp_log1p",
)


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text())


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise RuntimeError(f"Refusing to write empty P17 table: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def implementation_sha256() -> dict[str, str]:
    paths = (
        "scripts/run_v2_p17_iv_controlled_matrix.py",
        "scripts/prepare_v2_p17_iv_controlled_matrix.py",
        "scripts/train_v2_model.py",
        "scripts/evaluate_v2_model.py",
        "scripts/v2_training_components.py",
        "scripts/v2_bigwig_dataset.py",
    )
    return {path: sha256(REPO_ROOT / path) for path in paths}


def require_controller_scope() -> None:
    state = read_json(STATE_PATH)
    gate = state.get("approvals", {}).get("G4_gpu_experiments", {})
    if not (
        state.get("current_phase") == "P17"
        and state.get("status") == "RUNNING"
        and gate.get("approved") is True
        and gate.get("scope") == "p17_iv_controlled_ablation_and_basenji2"
    ):
        raise RuntimeError("P17 requires controller state P17/RUNNING and G4:p17_iv_controlled_ablation_and_basenji2")


def load_spec() -> dict[str, Any]:
    if not SPEC_PATH.is_file():
        raise FileNotFoundError(f"P17 preregistration missing: {SPEC_PATH}")
    spec = read_json(SPEC_PATH)
    matrix = spec.get("matrix", {})
    configurations = matrix.get("configurations", [])
    expected = {"B_iv_dual", "B_128bp_only", "Basenji2_style"}
    if not (
        spec.get("phase") == "P17"
        and spec.get("status") == "preregistered"
        and spec.get("allowed_chromosomes") == ALLOWED_CHROMOSOMES
        and matrix.get("folds") == [1, 2, 3, 4, 5]
        and matrix.get("seeds") == [20260714, 20260715, 20260716]
        and matrix.get("total_jobs") == 45
        and {row.get("id") for row in configurations} == expected
        and spec.get("final_test_access") == "prohibited"
    ):
        raise RuntimeError("P17 static matrix contract mismatch")
    p16_execution = METADATA_ROOT / "p16_iv_training_normalization_execution.json"
    p16_review = METADATA_ROOT / "audits/P16/review.json"
    if read_json(p16_execution).get("status") != "completed" or read_json(p16_review).get("status") != "PASS":
        raise RuntimeError("P17 requires completed P16 and R16")
    for source in spec["source"].values():
        entries = source if isinstance(source, list) else [source]
        for entry in entries:
            if not isinstance(entry, dict) or "path" not in entry:
                continue
            path = REPO_ROOT / entry["path"]
            if not path.is_file() or sha256(path) != entry.get("sha256"):
                raise RuntimeError(f"P17 source hash mismatch: {path}")
    for fold in matrix["folds"]:
        for role in ("train", "valid"):
            path = P15_ROOT / "intervals" / f"fold_{fold}" / f"{role}.tsv"
            rows = read_tsv(path)
            if not rows or any(row["chromosome"] not in ALLOWED_CHROMOSOMES or row["role"] != role for row in rows):
                raise RuntimeError(f"P17 prohibited P15 interval content: {path}")
    return spec


def jobs_from(spec: dict[str, Any]) -> list[dict[str, Any]]:
    jobs: list[dict[str, Any]] = []
    for config in spec["matrix"]["configurations"]:
        for fold in spec["matrix"]["folds"]:
            for seed in spec["matrix"]["seeds"]:
                root = RUN_ROOT / config["id"] / f"seed_{seed}" / f"fold_{fold}"
                evaluation = RESULT_ROOT / "evaluations" / config["id"] / f"seed_{seed}" / f"fold_{fold}.json"
                jobs.append(
                    {
                        "job_id": f"p17:{config['id']}:seed{seed}:fold{fold}",
                        "configuration": config["id"],
                        "model": config["model"],
                        "loss": config["loss"],
                        "fold": fold,
                        "seed": seed,
                        "gene_weight": config["gene_weight"],
                        "loss_1bp_weight": config["loss_1bp_weight"],
                        "loss_128bp_weight": config["loss_128bp_weight"],
                        "mean_column": f"fold_{fold}_train_nonzero_mean",
                        "output_dir": str(root),
                        "checkpoint_path": str(root / "checkpoint.pt"),
                        "run_path": str(root / "run.json"),
                        "evaluation_path": str(evaluation),
                        "log_path": str(LOG_ROOT / config["id"] / f"seed_{seed}" / f"fold_{fold}.log"),
                    }
                )
    if len(jobs) != 45 or len({job["job_id"] for job in jobs}) != 45:
        raise RuntimeError("P17 generated an invalid job matrix")
    return jobs


def command_pair(job: dict[str, Any], spec: dict[str, Any]) -> tuple[list[str], list[str]]:
    training = spec["training"]
    fold = int(job["fold"])
    train_interval = P15_ROOT / "intervals" / f"fold_{fold}" / "train.tsv"
    valid_interval = P15_ROOT / "intervals" / f"fold_{fold}" / "valid.tsv"
    train = [
        str(Path(sys.executable).resolve()), "-m", "scripts.train_v2_model",
        "--model", str(job["model"]), "--fold", str(fold), "--seed", str(job["seed"]),
        "--loss", str(job["loss"]), "--max-steps", str(training["max_steps"]),
        "--sequence-length", str(training["sequence_length"]), "--learning-rate", str(training["learning_rate"]),
        "--weight-decay", str(training["weight_decay"]), "--gene-weight", str(job["gene_weight"]),
        "--loss-1bp-weight", str(job["loss_1bp_weight"]), "--loss-128bp-weight", str(job["loss_128bp_weight"]),
        "--max-shift-bp", str(training["max_shift_bp"]),
        "--reverse-complement-probability", str(training["reverse_complement_probability"]),
        "--hidden-channels", str(training["hidden_channels"]), "--mean-column", str(job["mean_column"]),
        "--means-path", str(training["means_path"]),
        "--track-manifest", "alphagenome_custom/metadata/v2/replicate_holdout_v1/training_track_manifest.tsv",
        "--group-manifest", "alphagenome_custom/metadata/v2/replicate_holdout_v1/training_group_manifest.tsv",
        "--train-intervals-path", str(train_interval.relative_to(REPO_ROOT)),
        "--output-dir", str(job["output_dir"]), "--device", "cuda",
    ]
    evaluate = [
        str(Path(sys.executable).resolve()), "-m", "scripts.evaluate_v2_model",
        "--model", str(job["model"]), "--training-loss", str(job["loss"]), "--fold", str(fold),
        "--checkpoint", str(job["checkpoint_path"]), "--sequence-length", str(training["sequence_length"]),
        "--hidden-channels", str(training["hidden_channels"]), "--mean-column", str(job["mean_column"]),
        "--means-path", str(training["means_path"]),
        "--track-manifest", str(spec["evaluation"]["track_manifest"]),
        "--group-manifest", str(spec["evaluation"]["group_manifest"]),
        "--prediction-track-indices", str(spec["evaluation"]["prediction_track_indices"]),
        "--valid-intervals-path", str(valid_interval.relative_to(REPO_ROOT)),
        "--output", str(job["evaluation_path"]), "--device", "cuda",
    ]
    serialized = " ".join([*train, *evaluate])
    if "test_locked" in serialized or "--final-test" in serialized or "/fold_0/" in serialized:
        raise RuntimeError("P17 command violates the no-final-test I--V contract")
    return train, evaluate


def completed_record(job: dict[str, Any], spec_sha: str) -> dict[str, Any] | None:
    checkpoint = REPO_ROOT / job["checkpoint_path"]
    run_path = REPO_ROOT / job["run_path"]
    evaluation_path = REPO_ROOT / job["evaluation_path"]
    if not all(path.is_file() for path in (checkpoint, run_path, evaluation_path)):
        return None
    try:
        run = read_json(run_path)
        evaluation = read_json(evaluation_path)
        if not (
            run.get("model") == job["model"]
            and int(run.get("fold", -1)) == int(job["fold"])
            and int(run.get("seed", -1)) == int(job["seed"])
            and run.get("mean_column") == job["mean_column"]
            and evaluation.get("model") == job["model"]
            and int(evaluation.get("fold", -1)) == int(job["fold"])
            and evaluation.get("checkpoint_sha256") == sha256(checkpoint)
            and evaluation.get("evaluated_chromosomes") == ALLOWED_CHROMOSOMES
            and evaluation.get("locked_test_block_signal_reads") is False
            and float(evaluation.get("validation_core_coverage_fraction", 0.0)) == 1.0
        ):
            return None
        metrics = shared.validation_summary(evaluation)
    except Exception:
        return None
    return {
        **job,
        "status": "completed",
        "spec_sha256": spec_sha,
        "physical_gpu": None,
        "elapsed_seconds": None,
        "checkpoint_sha256": sha256(checkpoint),
        "run_sha256": sha256(run_path),
        "evaluation_sha256": sha256(evaluation_path),
        "trainable_parameters": int(run["trainable_parameters"]),
        "evaluated_chromosomes": ",".join(ALLOWED_CHROMOSOMES),
        "locked_test_block_signal_reads": False,
        **{key: metrics[key] for key in PRIMARY_METRICS},
    }


def run_queue(
    jobs: list[dict[str, Any]], spec: dict[str, Any], spec_sha: str, execution: dict[str, Any], gpus: list[int]
) -> list[dict[str, Any]]:
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
                train, evaluate = command_pair(job, spec)
                base.run_command(train, REPO_ROOT / job["log_path"], gpu)
                base.run_command(evaluate, REPO_ROOT / job["log_path"], gpu)
                record = completed_record(job, spec_sha)
                if record is None:
                    raise RuntimeError(f"P17 output failed contract validation: {job['job_id']}")
                record["physical_gpu"] = gpu
                record["elapsed_seconds"] = time.monotonic() - started
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
        raise RuntimeError(f"P17 job failed: {failures[0]}")
    return records


def summarize(records: list[dict[str, Any]]) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    by_fold: list[dict[str, Any]] = []
    for configuration in ("B_iv_dual", "B_128bp_only", "Basenji2_style"):
        for fold in range(1, 6):
            members = [row for row in records if row["configuration"] == configuration and int(row["fold"]) == fold]
            if len(members) != 3:
                raise RuntimeError(f"P17 incomplete {configuration} fold {fold}")
            by_fold.append(
                {
                    "configuration": configuration,
                    "fold": fold,
                    "seeds": 3,
                    **{metric: sum(float(row[metric]) for row in members) / len(members) for metric in PRIMARY_METRICS},
                }
            )
    indexed = {(row["configuration"], int(row["fold"])): row for row in by_fold}
    paired: list[dict[str, Any]] = []
    for comparison in ("B_128bp_only", "Basenji2_style"):
        for fold in range(1, 6):
            full = indexed[("B_iv_dual", fold)]
            other = indexed[(comparison, fold)]
            paired.append(
                {
                    "comparison": f"B_iv_dual_minus_{comparison}",
                    "fold": fold,
                    **{f"difference_{metric}": float(full[metric]) - float(other[metric]) for metric in PRIMARY_METRICS},
                    "interpretation": "Paired across matched fold and three nested seeds under the same strict I--V training and held-out interval contract.",
                }
            )
    return by_fold, paired


def main() -> None:
    require_controller_scope()
    spec = load_spec()
    spec_sha = sha256(SPEC_PATH)
    jobs = jobs_from(spec)
    selected_gpus, snapshots = base.wait_for_gpus(spec["gpu_policy"]["minimum_free_mib"])
    if not selected_gpus or not set(selected_gpus) <= {2, 3}:
        raise RuntimeError(f"P17 selected prohibited GPU(s): {selected_gpus}")
    existing_records: list[dict[str, Any]] = []
    if EXECUTION_PATH.is_file():
        existing = read_json(EXECUTION_PATH)
        if existing.get("spec_sha256") != spec_sha or existing.get("final_test_access") != "prohibited":
            raise RuntimeError("P17 prior execution is not resumable under this exact spec")
        for job in jobs:
            record = completed_record(job, spec_sha)
            if record is not None:
                existing_records.append(record)
        execution = existing
        execution.update({
            "status": "running", "completed_at": None,
            "resume_count": int(execution.get("resume_count", 0)) + 1,
            "resumed_at": utc_now(), "selected_physical_gpus": selected_gpus,
            "gpu_wait_snapshots": snapshots, "failures": [],
        })
    else:
        execution = {
            "schema_version": 1,
            "phase": "P17",
            "status": "running",
            "started_at": utc_now(),
            "completed_at": None,
            "hostname": socket.gethostname(),
            "working_directory": str(REPO_ROOT),
            "controller_invocation": "python -m scripts.v2_phase_controller run --phase P17",
            "registered_phase_command": "python -m scripts.run_v2_p17_iv_controlled_matrix",
            "spec_path": str(SPEC_PATH.relative_to(REPO_ROOT)),
            "spec_sha256": spec_sha,
            "implementation_sha256": implementation_sha256(),
            "registered_tasks": 45,
            "registered_tasks_completed": 0,
            "selected_physical_gpus": selected_gpus,
            "gpu_wait_snapshots": snapshots,
            "jobs": [],
            "failures": [],
            "final_test_access": "prohibited",
            "locked_test_block_signal_reads": 0,
        }
    execution["jobs"] = sorted(existing_records, key=lambda row: row["job_id"])
    execution["registered_tasks_completed"] = len(existing_records)
    base.atomic_json(EXECUTION_PATH, execution)
    remaining = [job for job in jobs if job["job_id"] not in {row["job_id"] for row in existing_records}]
    try:
        new_records = run_queue(remaining, spec, spec_sha, execution, selected_gpus)
        records = sorted(existing_records + new_records, key=lambda row: row["job_id"])
        if len(records) != 45:
            raise RuntimeError(f"P17 requires 45 complete records, observed {len(records)}")
        write_tsv(RESULT_ROOT / "p17_iv_controlled_records.tsv", records)
        by_fold, paired = summarize(records)
        write_tsv(RESULT_ROOT / "p17_iv_controlled_by_fold.tsv", by_fold)
        write_tsv(RESULT_ROOT / "p17_iv_controlled_paired_effects.tsv", paired)
        execution.update(
            {
                "status": "completed",
                "completed_at": utc_now(),
                "registered_tasks_completed": 45,
                "records_sha256": sha256(RESULT_ROOT / "p17_iv_controlled_records.tsv"),
                "by_fold_sha256": sha256(RESULT_ROOT / "p17_iv_controlled_by_fold.tsv"),
                "paired_effects_sha256": sha256(RESULT_ROOT / "p17_iv_controlled_paired_effects.tsv"),
            }
        )
        base.atomic_json(EXECUTION_PATH, execution)
    except Exception:
        execution.update({"status": "failed", "completed_at": utc_now()})
        base.atomic_json(EXECUTION_PATH, execution)
        raise
    print(json.dumps({"status": "completed", "registered_tasks_completed": 45}, indent=2))


if __name__ == "__main__":
    main()
