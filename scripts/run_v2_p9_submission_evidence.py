#!/usr/bin/env python3
"""Run the registered development-only P9 submission-evidence matrix."""

from __future__ import annotations

import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import socket
import subprocess
import sys
import time
from typing import Any
import threading
import queue
import torch

from scripts import run_v2_p6b as base
from scripts import run_v2_p6b_amendment as shared
from scripts import v2_gpu_resources
from scripts import v2_submission_statistics as submission_statistics


REPO_ROOT = Path(__file__).resolve().parents[1]
METADATA_DIR = REPO_ROOT / "alphagenome_custom/metadata/v2"
SPEC_PATH = METADATA_DIR / "p9_submission_evidence_spec.json"
EXECUTION_PATH = METADATA_DIR / "p9_submission_evidence_execution.json"
RESULTS_PATH = METADATA_DIR / "p9_submission_evidence_results.tsv"
EFFECTS_PATH = METADATA_DIR / "p9_submission_evidence_paired_effects.tsv"
RAW_EFFECTS_PATH = METADATA_DIR / "p9_submission_evidence_raw_paired_effects.tsv"
RAW_FACTORIAL_EFFECTS_PATH = METADATA_DIR / "p9_submission_evidence_raw_factorial_effects.tsv"
FACTORIAL_EFFECTS_PATH = METADATA_DIR / "p9_submission_evidence_factorial_effects.tsv"
TRACK_RESULTS_PATH = METADATA_DIR / "p9_submission_evidence_per_track.tsv"
RUN_ROOT = Path("runs/v2_p9_submission_evidence_20260813")
LOG_ROOT = Path("logs/v2_p9_submission_evidence_20260813")
IMPLEMENTATION_PATHS = (
    "scripts/run_v2_p9_submission_evidence.py",
    "scripts/train_v2_model.py",
    "scripts/evaluate_v2_model.py",
    "scripts/v2_training_components.py",
    "scripts/v2_submission_statistics.py",
    "scripts/v2_bigwig_dataset.py",
    "scripts/run_v2_p6b.py",
    "scripts/run_v2_p6b_amendment.py",
    "scripts/v2_gpu_resources.py",
    "scripts/v2_subprocess.py",
)
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


def git_commit() -> str:
    return subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=REPO_ROOT,
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text())


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _repository_file(relative: object) -> Path:
    raw = Path(str(relative))
    path = raw.resolve() if raw.is_absolute() else (REPO_ROOT / raw).resolve()
    path.relative_to(REPO_ROOT.resolve())
    if not path.is_file():
        raise FileNotFoundError(path)
    return path


def _ordered_digest(rows: list[tuple[str, ...]]) -> str:
    digest = hashlib.sha256()
    for row in rows:
        digest.update(("\t".join(row) + "\n").encode("utf-8"))
    return digest.hexdigest()


def implementation_sha256() -> dict[str, str]:
    return {
        relative: base.sha256(REPO_ROOT / relative)
        for relative in IMPLEMENTATION_PATHS
    }


def audit_p6b_reuse_manifest(spec: dict[str, Any]) -> dict[str, Any]:
    """Resolve and hash the 30 frozen P6B A/B paper-loss reuse artifacts."""

    contract = spec.get("p6b_reuse_contract", {})
    source_path = _repository_file(spec["locked_inputs"]["p6b_execution"])
    source = read_json(source_path)
    selected = [
        row
        for row in source.get("formal_jobs", [])
        if row.get("model") in {"A", "B"} and row.get("loss") == "paper"
    ]
    selected.sort(key=lambda row: (str(row["model"]), int(row["fold"]), int(row["seed"])))
    expected_pairs = {
        (model, fold, seed)
        for model in ("A", "B")
        for fold in spec["matrix"]["folds"]
        for seed in spec["matrix"]["seeds"]
    }
    actual_pairs = {
        (str(row.get("model")), int(row.get("fold", -1)), int(row.get("seed", -1)))
        for row in selected
    }
    if not (
        int(contract.get("record_count", -1)) == len(selected) == 30
        and actual_pairs == expected_pairs
        and contract.get("source_execution") == spec["locked_inputs"]["p6b_execution"]
    ):
        raise RuntimeError("P9 frozen P6B reuse membership mismatch")

    records: list[dict[str, Any]] = []
    digest_rows: list[tuple[str, ...]] = []
    for row in selected:
        output_dir = str(row["output_dir"])
        paths = {
            "checkpoint": str(row["checkpoint_path"]),
            "run": f"{output_dir}/run.json",
            "validation": str(row["validation_path"]),
        }
        hashes = {
            name: base.sha256(_repository_file(relative))
            for name, relative in paths.items()
        }
        if hashes["checkpoint"] != row.get("checkpoint_sha256"):
            raise RuntimeError(
                f"P9 P6B checkpoint changed: {row['model']}/fold{row['fold']}/seed{row['seed']}"
            )
        record = {
            "model": str(row["model"]),
            "loss": "paper",
            "fold": int(row["fold"]),
            "seed": int(row["seed"]),
            "checkpoint_path": paths["checkpoint"],
            "checkpoint_sha256": hashes["checkpoint"],
            "run_path": paths["run"],
            "run_sha256": hashes["run"],
            "validation_path": paths["validation"],
            "validation_sha256": hashes["validation"],
        }
        records.append(record)
        digest_rows.append(
            tuple(str(record[key]) for key in (
                "model",
                "loss",
                "fold",
                "seed",
                "checkpoint_path",
                "checkpoint_sha256",
                "run_path",
                "run_sha256",
                "validation_path",
                "validation_sha256",
            ))
        )
    digest = _ordered_digest(digest_rows)
    if digest != contract.get("manifest_sha256"):
        raise RuntimeError("P9 frozen P6B reuse manifest hash mismatch")
    return {
        "status": "passed",
        "record_count": len(records),
        "artifact_count": 3 * len(records),
        "manifest_sha256": digest,
        "records": records,
    }


def audit_frozen_dataset(spec: dict[str, Any]) -> dict[str, Any]:
    """Hash every registered training track once without reading BigWig signal."""

    contract = spec.get("dataset_contract", {})
    track_manifest = _repository_file(spec["locked_inputs"]["track_manifest"])
    group_manifest = _repository_file(spec["locked_inputs"]["group_manifest"])
    means_path = _repository_file(spec["locked_inputs"]["means"])
    track_rows = read_tsv(track_manifest)
    group_rows = read_tsv(group_manifest)
    means_rows = read_tsv(means_path)
    expected_count = int(contract.get("track_count", -1))
    track_ids = [row.get("group_id", "") for row in track_rows]
    group_ids = [row.get("group_id", "") for row in group_rows]
    mean_ids = [row.get("group_id", "") for row in means_rows]
    if not (
        expected_count == int(contract.get("group_count", -1)) == 241
        and len(track_rows) == len(group_rows) == len(means_rows) == expected_count
        and len(set(track_ids)) == expected_count
        and track_ids == group_ids == mean_ids
        and all(row.get("include_formal_v2") == "True" for row in group_rows)
    ):
        raise RuntimeError("P9 frozen 241-track order or membership mismatch")

    registered: list[tuple[str, str, str]] = []
    observed: list[tuple[str, str, str]] = []
    total_bytes = 0
    for row in track_rows:
        group_id = row.get("group_id", "")
        relative = row.get("output_path", "")
        expected_sha = row.get("output_sha256", "")
        if len(expected_sha) != 64:
            raise RuntimeError(f"P9 missing registered BigWig hash: {group_id}")
        path = _repository_file(relative)
        actual_sha = base.sha256(path)
        if actual_sha != expected_sha:
            raise RuntimeError(f"P9 BigWig hash mismatch: {group_id}")
        registered.append((group_id, relative, expected_sha))
        observed.append((group_id, relative, actual_sha))
        total_bytes += path.stat().st_size
    return {
        "status": "passed",
        "track_count": len(track_rows),
        "group_count": len(group_rows),
        "verified_bigwig_count": len(observed),
        "verified_bigwig_bytes": total_bytes,
        "ordered_group_ids_sha256": hashlib.sha256(
            ("\n".join(track_ids) + "\n").encode("utf-8")
        ).hexdigest(),
        "registered_bigwig_manifest_digest": _ordered_digest(registered),
        "observed_bigwig_manifest_digest": _ordered_digest(observed),
        "bigwig_signal_reads": 0,
    }


def require_contract(spec: dict[str, Any]) -> dict[str, Any]:
    state = read_json(METADATA_DIR / "execution_state.json")
    gate = state.get("approvals", {}).get("G4_gpu_experiments", {})
    matrix = spec.get("matrix", {})
    analysis = spec.get("analysis", {})
    if not (
        state.get("current_phase") == "P9"
        and state.get("status") == "RUNNING"
        and gate.get("approved") is True
        and gate.get("scope") == "p9_submission_evidence_matrix"
        and spec.get("phase") == "P9"
        and spec.get("final_test_access") == "prohibited"
        and matrix.get("folds") == [1, 2, 3, 4, 5]
        and matrix.get("seeds") == [20260714, 20260715, 20260716]
        and matrix.get("loss") == "paper"
        and matrix.get("new_training_jobs") == 59
        and matrix.get("registered_records") == 90
        and matrix.get("reused_records") == 31
        and tuple(analysis.get("metrics", [])) == METRIC_NAMES
    ):
        raise RuntimeError("P9 controller, approval, or frozen matrix contract mismatch")
    for key, relative in spec["locked_inputs"].items():
        if base.sha256(REPO_ROOT / relative) != spec["input_sha256"][key]:
            raise RuntimeError(f"P9 locked input hash mismatch: {key}")
    observed_implementation = implementation_sha256()
    if observed_implementation != spec.get("implementation_sha256"):
        raise RuntimeError("P9 frozen implementation hash mismatch")
    lock = read_json(METADATA_DIR / "final_test_lock.json")
    if not (lock.get("test_consumed") is True and lock.get("test_status") == "completed"):
        raise RuntimeError("P9 requires the preserved completed P6C test lock")
    return {
        "dataset": audit_frozen_dataset(spec),
        "p6b_reuse": audit_p6b_reuse_manifest(spec),
        "implementation_sha256": observed_implementation,
    }


def matrix_jobs(spec: dict[str, Any]) -> list[dict[str, Any]]:
    jobs = []
    reused = spec["matrix"]["reused_record"]
    for seed in spec["matrix"]["seeds"]:
        for fold in spec["matrix"]["folds"]:
            for config in spec["matrix"]["configurations"]:
                reuse_p8 = (
                    config["id"] == reused["configuration"]
                    and seed == reused["seed"]
                    and fold == reused["fold"]
                )
                reuse_p6b = config.get("reuse_source") == "P6B"
                root = (
                    Path("runs/v2_p8_b_no_lora_fold1_20260730")
                    / "model_B_no_lora_paper/fold_1"
                    if reuse_p8
                    else Path("runs/v2_p6b_six_chromosome_corefix_20260722")
                    / "formal"
                    / config["model"]
                    / "paper"
                    / f"seed_{seed}"
                    / f"fold_{fold}"
                    if reuse_p6b
                    else RUN_ROOT
                    / config["id"]
                    / f"seed_{seed}"
                    / f"fold_{fold}"
                )
                jobs.append(
                    {
                        "job_kind": "submission_evidence",
                        "job_id": f"p9:{config['id']}:seed{seed}:fold{fold}",
                        "configuration": config["id"],
                        "semantics": config["semantics"],
                        "model": config["model"],
                        "loss": "paper",
                        "seed": seed,
                        "fold": fold,
                        "gene_weight": spec["training"]["gene_weight"],
                        "max_shift_bp": spec["training"]["max_shift_bp"],
                        "reverse_complement_probability": spec["training"][
                            "reverse_complement_probability"
                        ],
                        "hidden_channels": config.get(
                            "hidden_channels", spec["training"]["hidden_channels_default"]
                        ),
                        "mean_column": f"fold_{fold}_train_nonzero_mean",
                        "base_job": reuse_p8 or reuse_p6b,
                        "reused_source": (
                            "P8" if reuse_p8 else "P6B" if reuse_p6b else None
                        ),
                        "output_dir": str(root),
                        "checkpoint_path": str(root / "checkpoint.pt"),
                        "validation_path": str(root / "validation_full_metrics.json"),
                        "run_path": str(root / "run.json"),
                        "job_path": str(
                            RUN_ROOT
                            / "jobs"
                            / config["id"]
                            / f"seed_{seed}"
                            / f"fold_{fold}.json"
                        ),
                        "log_path": str(
                            LOG_ROOT
                            / config["id"]
                            / f"seed_{seed}"
                            / f"fold_{fold}.log"
                        ),
                    }
                )
    return jobs


def commands_for_job(
    job: dict[str, Any], spec: dict[str, Any]
) -> tuple[list[str] | None, list[str] | None]:
    if job["base_job"]:
        return None, None
    training = spec["training"]
    train = shared.v2_subprocess.module_command(
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
        "--hidden-channels", str(job["hidden_channels"]),
        "--mean-column", job["mean_column"],
        "--output-dir", job["output_dir"],
        "--device", "cuda",
    )
    evaluate = shared.v2_subprocess.module_command(
        "evaluate_v2_model",
        "--model", job["model"],
        "--training-loss", job["loss"],
        "--fold", str(job["fold"]),
        "--checkpoint", job["checkpoint_path"],
        "--sequence-length", str(training["sequence_length"]),
        "--hidden-channels", str(job["hidden_channels"]),
        "--mean-column", job["mean_column"],
        "--output", job["validation_path"],
        "--device", "cuda",
    )
    serialized = " ".join([*train, *evaluate])
    if "--final-test" in serialized or "test_locked.tsv" in serialized:
        raise RuntimeError("P9 command unexpectedly requests final-test access")
    return train, evaluate


def _record_from_validation(job: dict[str, Any], spec_sha: str) -> dict[str, Any]:
    checkpoint = REPO_ROOT / job["checkpoint_path"]
    run_path = REPO_ROOT / job["run_path"]
    validation_path = REPO_ROOT / job["validation_path"]
    if not (checkpoint.is_file() and run_path.is_file() and validation_path.is_file()):
        raise RuntimeError(f"P9 job output missing: {job['job_id']}")
    run = read_json(run_path)
    validation = read_json(validation_path)
    summary = shared.validation_summary(validation)
    if not (
        run.get("model") == job["model"]
        and run.get("loss") == job["loss"]
        and run.get("fold") == job["fold"]
        and run.get("seed") == job["seed"]
        and run.get("hidden_channels") == job["hidden_channels"]
        and validation.get("checkpoint_sha256") == base.sha256(checkpoint)
        and validation.get("locked_test_block_signal_reads") is False
        and validation.get("fold") == job["fold"]
    ):
        raise RuntimeError(f"P9 run/validation contract mismatch: {job['job_id']}")
    trainable_parameters = run.get("trainable_parameters")
    if trainable_parameters is None:
        checkpoint_payload = torch.load(
            checkpoint, map_location="cpu", weights_only=False
        )
        trainable_parameters = sum(
            value.numel()
            for value in checkpoint_payload["trainable_model_state"].values()
        )
    return {
        **job,
        "schema_version": 1,
        "status": "completed",
        "spec_sha256": spec_sha,
        "checkpoint_sha256": base.sha256(checkpoint),
        "run_sha256": base.sha256(run_path),
        "validation_sha256": base.sha256(validation_path),
        "trainable_parameters": int(trainable_parameters),
        "adaptation": run.get("adaptation"),
        "reused": bool(job["base_job"]),
        **summary,
    }


def completed_job(job: dict[str, Any], spec_sha: str) -> dict[str, Any] | None:
    if job["base_job"]:
        return _record_from_validation(job, spec_sha)
    record_path = REPO_ROOT / job["job_path"]
    if not record_path.is_file():
        return None
    try:
        record = read_json(record_path)
        if record.get("spec_sha256") != spec_sha or record.get("status") != "completed":
            return None
        verified = _record_from_validation(job, spec_sha)
        if any(
            record.get(key) != verified[key]
            for key in ("checkpoint_sha256", "run_sha256", "validation_sha256")
        ):
            return None
    except Exception:
        return None
    verified["physical_gpu"] = record.get("physical_gpu")
    verified["elapsed_seconds"] = record.get("elapsed_seconds")
    verified["reused"] = True
    return verified


def run_job(
    job: dict[str, Any], gpu: int, spec: dict[str, Any], spec_sha: str
) -> dict[str, Any]:
    existing = completed_job(job, spec_sha)
    if existing is not None:
        return existing
    train, evaluate = commands_for_job(job, spec)
    if train is None or evaluate is None:
        raise RuntimeError("Missing P9 commands for a new job")
    started = time.monotonic()
    base.run_command(train, REPO_ROOT / job["log_path"], gpu)
    base.run_command(evaluate, REPO_ROOT / job["log_path"], gpu)
    record = _record_from_validation(job, spec_sha)
    record["physical_gpu"] = gpu
    record["elapsed_seconds"] = time.monotonic() - started
    base.atomic_json(REPO_ROOT / job["job_path"], record)
    return record


def run_queue(
    jobs: list[dict[str, Any]],
    gpus: list[int],
    spec: dict[str, Any],
    spec_sha: str,
    execution: dict[str, Any],
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
            try:
                record = run_job(job, gpu, spec, spec_sha)
            except Exception as error:
                with lock:
                    failures.append(
                        {"job_id": job["job_id"], "physical_gpu": gpu, "error": repr(error)}
                    )
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
        raise RuntimeError(f"P9 job failed: {failures[0]}")
    return records


def baseline_records(spec: dict[str, Any]) -> list[dict[str, Any]]:
    execution = read_json(REPO_ROOT / spec["baseline"]["source_execution"])
    records = [
        row
        for row in execution["formal_jobs"]
        if row["model"] == spec["baseline"]["model"]
        and row["loss"] == spec["baseline"]["loss"]
    ]
    if len(records) != 15:
        raise RuntimeError("P9 baseline must contain exactly 15 B/paper records")
    return records


def loss_records(spec: dict[str, Any]) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    execution = read_json(REPO_ROOT / spec["baseline"]["source_execution"])
    paper = [row for row in execution["formal_jobs"] if row["model"] == "B" and row["loss"] == "paper"]
    log1p = [row for row in execution["formal_jobs"] if row["model"] == "B" and row["loss"] == "log1p_mse"]
    if len(paper) != 15 or len(log1p) != 15:
        raise RuntimeError("P9 loss comparison requires 15 paired records per loss")
    return paper, log1p


def paired_bootstrap_interval(
    pairs: list[dict[str, Any]], metric: str, analysis: dict[str, Any]
) -> tuple[float, float, float]:
    try:
        matrix = submission_statistics.records_to_balanced_matrix(
            pairs,
            value_key="difference",
            expected_blocks=analysis["expected_folds"],
            expected_repeats=analysis["expected_seeds"],
        )
    except ValueError as error:
        raise RuntimeError(
            f"Incomplete paired fold/seed matrix for {metric}: {error}"
        ) from error
    interval = analysis["confidence_interval"]
    plan = submission_statistics.make_hierarchical_resample_plan(
        n_blocks=len(matrix.blocks),
        n_repeats=len(matrix.repeats),
        iterations=int(interval["resamples"]),
        seed=int(interval["seed"]),
    )
    estimate = submission_statistics.hierarchical_block_bootstrap(
        matrix.values,
        plan,
        confidence_level=float(interval["level"]),
    )
    return estimate.estimate, estimate.ci_lower, estimate.ci_upper


def paired_effect_rows(
    candidate_records: list[dict[str, Any]], spec: dict[str, Any]
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    baseline = baseline_records(spec)
    configurations: dict[str, list[dict[str, Any]]] = {}
    for row in candidate_records:
        configurations.setdefault(str(row["configuration"]), []).append(row)
    configurations.pop("full_worm_lora", None)
    paper, log1p = loss_records(spec)
    configurations["b_log1p_mse_auxiliary"] = log1p
    raw: list[dict[str, Any]] = []
    result: list[dict[str, Any]] = []
    for configuration, records in sorted(configurations.items()):
        comparison_baseline = (
            paper if configuration == "b_log1p_mse_auxiliary" else baseline
        )
        for metric in METRIC_NAMES:
            pairs = submission_statistics.matched_differences(
                records,
                comparison_baseline,
                value_key=metric,
            )
            raw.extend(
                {
                    "configuration": configuration,
                    "reference": "B/paper",
                    "metric": metric,
                    "fold": row["fold"],
                    "seed": row["seed"],
                    "candidate_value": row["candidate_value"],
                    "reference_value": row["comparator_value"],
                    "difference": row["difference"],
                }
                for row in pairs
            )
            effect, low, high = paired_bootstrap_interval(pairs, metric, spec["analysis"])
            result.append(
                {
                    "configuration": configuration,
                    "reference": "B/paper",
                    "metric": metric,
                    "paired_runs": len(pairs),
                    "folds": len({row["fold"] for row in pairs}),
                    "seeds": len({row["seed"] for row in pairs}),
                    "mean_paired_difference": effect,
                    "ci_level": spec["analysis"]["confidence_interval"]["level"],
                    "ci_low": low,
                    "ci_high": high,
                    "ci_method": spec["analysis"]["confidence_interval"]["method"],
                }
            )
    return raw, result


def factorial_effect_rows(
    records: list[dict[str, Any]], spec: dict[str, Any]
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    keyed = {
        (str(row["configuration"]), int(row["fold"]), int(row["seed"])): row
        for row in records
    }
    cells = {
        "a": "frozen_no_worm_no_lora",
        "lora": "lora_only",
        "worm": "no_lora",
        "full": "full_worm_lora",
    }
    raw: list[dict[str, Any]] = []
    summaries: list[dict[str, Any]] = []
    for metric in METRIC_NAMES:
        for fold in spec["matrix"]["folds"]:
            for seed in spec["matrix"]["seeds"]:
                values = {
                    name: float(keyed[(configuration, fold, seed)][metric])
                    for name, configuration in cells.items()
                }
                effects = {
                    "worm_embedding_main": 0.5
                    * ((values["worm"] - values["a"]) + (values["full"] - values["lora"])),
                    "lora_main": 0.5
                    * ((values["lora"] - values["a"]) + (values["full"] - values["worm"])),
                    "worm_lora_interaction": values["full"]
                    - values["worm"]
                    - values["lora"]
                    + values["a"],
                }
                raw.extend(
                    {
                        "effect": effect,
                        "metric": metric,
                        "fold": fold,
                        "seed": seed,
                        "value": value,
                    }
                    for effect, value in effects.items()
                )
        for effect in ("worm_embedding_main", "lora_main", "worm_lora_interaction"):
            pairs = [
                {"fold": row["fold"], "seed": row["seed"], "difference": row["value"]}
                for row in raw
                if row["metric"] == metric and row["effect"] == effect
            ]
            estimate, low, high = paired_bootstrap_interval(pairs, metric, spec["analysis"])
            summaries.append(
                {
                    "effect": effect,
                    "metric": metric,
                    "paired_runs": 15,
                    "folds": 5,
                    "seeds": 3,
                    "mean_effect": estimate,
                    "ci_level": spec["analysis"]["confidence_interval"]["level"],
                    "ci_low": low,
                    "ci_high": high,
                    "ci_method": spec["analysis"]["confidence_interval"]["method"],
                }
            )
    return raw, summaries


def per_track_rows(records: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rows = []
    for record in records:
        validation = read_json(REPO_ROOT / record["validation_path"])
        full = validation["full_metrics"]
        components = {
            "pearson_128bp_log1p": full["log1p_128bp_sum"]["per_track_pearson"],
            "gene_exon_coverage_pearson_log1p": full["gene_exon_coverage"]["per_track_pearson"],
        }
        track_ids = set(components["pearson_128bp_log1p"])
        if set(components["gene_exon_coverage_pearson_log1p"]) != track_ids:
            raise RuntimeError(
                f"P9 per-track metric identifiers differ: {record['job_id']}"
            )
        sources = {
            **components,
            "primary_biological_score": {
                track_id: 0.5
                * (
                    float(components["pearson_128bp_log1p"][track_id])
                    + float(
                        components["gene_exon_coverage_pearson_log1p"][track_id]
                    )
                )
                for track_id in sorted(track_ids)
            },
        }
        for metric, values in sources.items():
            for track_id in sorted(values):
                rows.append(
                    {
                        "configuration": record["configuration"],
                        "fold": record["fold"],
                        "seed": record["seed"],
                        "track_id": track_id,
                        "metric": metric,
                        "value": values[track_id],
                    }
                )
    return rows


def aggregate_rows(records: list[dict[str, Any]]) -> list[dict[str, Any]]:
    groups: dict[str, list[dict[str, Any]]] = {}
    for row in records:
        groups.setdefault(str(row["configuration"]), []).append(row)
    result = []
    for configuration, members in sorted(groups.items()):
        row: dict[str, Any] = {
            "configuration": configuration,
            "model": members[0]["model"],
            "loss": members[0]["loss"],
            "jobs": len(members),
            "folds": len({member["fold"] for member in members}),
            "seeds": len({member["seed"] for member in members}),
            "mean_trainable_parameters": sum(member["trainable_parameters"] for member in members) / len(members),
        }
        for metric in METRIC_NAMES:
            row[f"mean_{metric}"] = sum(float(member[metric]) for member in members) / len(members)
        result.append(row)
    return result


def write_tsv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise ValueError(f"Cannot write empty P9 table: {path}")
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def main() -> None:
    spec = read_json(SPEC_PATH)
    preflight = require_contract(spec)
    spec_sha = base.sha256(SPEC_PATH)
    jobs = matrix_jobs(spec)
    new_jobs = [job for job in jobs if not job["base_job"]]
    if len(jobs) != 90 or len(new_jobs) != 59:
        raise RuntimeError("P9 must register 90 records and exactly 59 new jobs")
    selected_gpus, snapshots = base.wait_for_gpus(spec["gpu_policy"]["minimum_free_mib"])
    if not selected_gpus or not set(selected_gpus) <= set(
        spec["gpu_policy"]["allowed_physical_indices"]
    ):
        raise RuntimeError(f"P9 selected disallowed physical GPUs: {selected_gpus}")
    execution: dict[str, Any] = {
        "schema_version": 1,
        "phase": "P9",
        "status": "running",
        "started_at": utc_now(),
        "completed_at": None,
        "git_commit": git_commit(),
        "hostname": socket.gethostname(),
        "working_directory": str(REPO_ROOT),
        "python_executable": sys.executable,
        "controller_invocation": "python -m scripts.v2_phase_controller run --phase P9",
        "registered_phase_command": "python -m scripts.run_v2_p9_submission_evidence",
        "spec_path": str(SPEC_PATH.relative_to(REPO_ROOT)),
        "spec_sha256": spec_sha,
        "implementation_sha256": preflight["implementation_sha256"],
        "dataset_preflight": preflight["dataset"],
        "p6b_reuse_preflight": preflight["p6b_reuse"],
        "selected_physical_gpus": selected_gpus,
        "gpu_wait_snapshots": snapshots,
        "registered_records_expected": 90,
        "new_training_jobs_expected": 59,
        "reused_records_expected": 31,
        "jobs": [],
        "failures": [],
        "locked_test_block_signal_reads": 0,
        "final_test_access": "prohibited",
    }
    base.atomic_json(EXECUTION_PATH, execution)
    try:
        records = run_queue(jobs, selected_gpus, spec, spec_sha, execution)
        records.sort(
            key=lambda row: (
                str(row["configuration"]),
                int(row["fold"]),
                int(row["seed"]),
            )
        )
        result_rows = aggregate_rows(records)
        raw_effect_rows, effect_rows = paired_effect_rows(records, spec)
        raw_factorial_rows, factorial_rows = factorial_effect_rows(records, spec)
        track_rows = per_track_rows(records)
        write_tsv(RESULTS_PATH, result_rows)
        write_tsv(EFFECTS_PATH, effect_rows)
        write_tsv(RAW_EFFECTS_PATH, raw_effect_rows)
        write_tsv(RAW_FACTORIAL_EFFECTS_PATH, raw_factorial_rows)
        write_tsv(FACTORIAL_EFFECTS_PATH, factorial_rows)
        write_tsv(TRACK_RESULTS_PATH, track_rows)
        execution.update(
            {
                "status": "completed",
                "completed_at": utc_now(),
                "registered_records": records,
                "new_training_jobs_completed": sum(not row["base_job"] for row in records),
                "reused_records": sum(bool(row["base_job"]) for row in records),
                "results_path": str(RESULTS_PATH.relative_to(REPO_ROOT)),
                "results_sha256": base.sha256(RESULTS_PATH),
                "paired_effects_path": str(EFFECTS_PATH.relative_to(REPO_ROOT)),
                "paired_effects_sha256": base.sha256(EFFECTS_PATH),
                "raw_paired_effects_path": str(RAW_EFFECTS_PATH.relative_to(REPO_ROOT)),
                "raw_paired_effects_sha256": base.sha256(RAW_EFFECTS_PATH),
                "raw_factorial_effects_path": str(
                    RAW_FACTORIAL_EFFECTS_PATH.relative_to(REPO_ROOT)
                ),
                "raw_factorial_effects_sha256": base.sha256(RAW_FACTORIAL_EFFECTS_PATH),
                "factorial_effects_path": str(FACTORIAL_EFFECTS_PATH.relative_to(REPO_ROOT)),
                "factorial_effects_sha256": base.sha256(FACTORIAL_EFFECTS_PATH),
                "per_track_path": str(TRACK_RESULTS_PATH.relative_to(REPO_ROOT)),
                "per_track_sha256": base.sha256(TRACK_RESULTS_PATH),
                "resources_after": v2_gpu_resources.snapshot(),
            }
        )
    except Exception as error:
        execution["status"] = "failed"
        execution["completed_at"] = utc_now()
        execution["failures"] = execution.get("failures", []) + [str(error)]
        base.atomic_json(EXECUTION_PATH, execution)
        raise
    base.atomic_json(EXECUTION_PATH, execution)


if __name__ == "__main__":
    main()
