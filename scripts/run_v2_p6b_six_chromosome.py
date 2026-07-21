#!/usr/bin/env python3
"""Run the preregistered six-chromosome P6B matrix on GPU 2/3."""

from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from scripts import run_v2_p6b as base
from scripts import run_v2_p6b_amendment as shared
from scripts import v2_gpu_resources
from scripts import v2_subprocess


REPO_ROOT = Path(__file__).resolve().parents[1]
METADATA_DIR = REPO_ROOT / "alphagenome_custom/metadata/v2"
SPEC_PATH = METADATA_DIR / "p6b_six_chromosome_spec.json"
EXECUTION_PATH = METADATA_DIR / "p6b_six_chromosome_execution.json"
RESULTS_PATH = METADATA_DIR / "p6b_six_chromosome_cv_results.tsv"
ABLATION_RESULTS_PATH = METADATA_DIR / "p6b_six_chromosome_ablation_results.tsv"
SELECTION_PATH = METADATA_DIR / "p6b_six_chromosome_selection.json"
LOCK_PATH = METADATA_DIR / "final_test_lock.json"
SPLIT_REGISTRY_PATH = METADATA_DIR / "split_registry_v2.json"
MEANS_PATH = METADATA_DIR / "track_nonzero_means_v2.tsv"
TEST_INTERVALS_PATH = REPO_ROOT / "alphagenome_custom/intervals/v2/test_locked.tsv"
RUN_ROOT = Path("runs/v2_p6b_six_chromosome_20260721")
LOG_ROOT = Path("logs/v2_p6b_six_chromosome_20260721")


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def formal_jobs(spec: dict[str, Any]) -> list[dict[str, Any]]:
    jobs = []
    for seed in spec["formal_matrix"]["seeds"]:
        for fold in spec["cv_folds"]:
            for model in spec["formal_matrix"]["models"]:
                for loss in spec["formal_matrix"]["losses"]:
                    root = RUN_ROOT / "formal" / model / loss / f"seed_{seed}" / f"fold_{fold}"
                    jobs.append(
                        {
                            "job_kind": "formal",
                            "job_id": f"formal:{model}:{loss}:seed{seed}:fold{fold}",
                            "model": model,
                            "loss": loss,
                            "seed": seed,
                            "fold": fold,
                            "gene_weight": spec["training"]["gene_weight"],
                            "max_shift_bp": spec["training"]["max_shift_bp"],
                            "reverse_complement_probability": spec["training"][
                                "reverse_complement_probability"
                            ],
                            "mean_column": f"fold_{fold}_train_nonzero_mean",
                            "base_job": False,
                            "output_dir": str(root),
                            "checkpoint_path": str(root / "checkpoint.pt"),
                            "validation_path": str(root / "validation_full_metrics.json"),
                            "job_path": str(
                                RUN_ROOT
                                / "jobs/formal"
                                / model
                                / loss
                                / f"seed_{seed}"
                                / f"fold_{fold}.json"
                            ),
                            "log_path": str(
                                LOG_ROOT
                                / "formal"
                                / model
                                / loss
                                / f"seed_{seed}"
                                / f"fold_{fold}.log"
                            ),
                        }
                    )
    return jobs


def ablation_jobs(spec: dict[str, Any]) -> list[dict[str, Any]]:
    seed = spec["formal_matrix"]["seeds"][0]
    jobs = []
    for config in spec["ablation_matrix"]:
        for fold in spec["cv_folds"]:
            root = RUN_ROOT / "ablations" / config["id"] / f"fold_{fold}"
            jobs.append(
                {
                    "job_kind": "ablation",
                    "job_id": f"ablation:{config['id']}:fold{fold}",
                    "ablation_id": config["id"],
                    "model": config["model"],
                    "loss": config["loss"],
                    "seed": seed,
                    "fold": fold,
                    "gene_weight": config["gene_weight"],
                    "max_shift_bp": config["max_shift_bp"],
                    "reverse_complement_probability": config[
                        "reverse_complement_probability"
                    ],
                    "mean_column": config.get(
                        "mean_column", f"fold_{fold}_train_nonzero_mean"
                    ),
                    "base_job": False,
                    "output_dir": str(root),
                    "checkpoint_path": str(root / "checkpoint.pt"),
                    "validation_path": str(root / "validation_full_metrics.json"),
                    "job_path": str(
                        RUN_ROOT / "jobs/ablations" / config["id"] / f"fold_{fold}.json"
                    ),
                    "log_path": str(
                        LOG_ROOT / "ablations" / config["id"] / f"fold_{fold}.log"
                    ),
                }
            )
    return jobs


def run_development(
    selected: dict[str, Any], gpu: int, spec: dict[str, Any], spec_sha: str
) -> dict[str, Any]:
    root = RUN_ROOT / "development_selected"
    command = v2_subprocess.module_command(
        "train_v2_model",
        "--model", selected["model"],
        "--fold", "0",
        "--seed", str(spec["development_retrain"]["seed"]),
        "--loss", selected["loss"],
        "--max-steps", str(spec["development_retrain"]["max_steps"]),
        "--sequence-length", str(spec["training"]["sequence_length"]),
        "--learning-rate", str(spec["training"]["learning_rate"]),
        "--weight-decay", str(spec["training"]["weight_decay"]),
        "--gene-weight", str(spec["training"]["gene_weight"]),
        "--max-shift-bp", str(spec["training"]["max_shift_bp"]),
        "--reverse-complement-probability",
        str(spec["training"]["reverse_complement_probability"]),
        "--hidden-channels", str(spec["training"]["hidden_channels"]),
        "--mean-column", spec["development_retrain"]["mean_column"],
        "--output-dir", str(root),
        "--device", "cuda",
    )
    log_path = REPO_ROOT / LOG_ROOT / "development_selected.log"
    base.run_command(command, log_path, gpu)
    run_path = REPO_ROOT / root / "run.json"
    checkpoint_path = REPO_ROOT / root / "checkpoint.pt"
    run = shared.read_json(run_path)
    return {
        "model": selected["model"],
        "loss": selected["loss"],
        "fold": 0,
        "seed": spec["development_retrain"]["seed"],
        "steps": run["steps"],
        "spec_sha256": spec_sha,
        "physical_gpu": gpu,
        "run_path": str(run_path.relative_to(REPO_ROOT)),
        "checkpoint_path": str(checkpoint_path.relative_to(REPO_ROOT)),
        "checkpoint_sha256": base.sha256(checkpoint_path),
        "intervals_sha256": run["intervals_sha256"],
        "means_sha256": run["means_sha256"],
        "log_path": str(log_path.relative_to(REPO_ROOT)),
    }


def main() -> None:
    base.require_approval()
    spec = shared.read_json(SPEC_PATH)
    split_registry = shared.read_json(SPLIT_REGISTRY_PATH)
    state = shared.read_json(METADATA_DIR / "execution_state.json")
    g5 = state.get("approvals", {}).get("G5_final_test", {})
    if not (
        spec.get("locked_before_training") is True
        and spec.get("split_revision") == "six_chromosome_blocks_v1"
        and spec.get("locked_test_block_access") == "prohibited"
        and split_registry.get("revision_id") == spec.get("split_revision")
        and g5.get("approved") is not True
    ):
        raise RuntimeError("P6B six-chromosome preregistration or test embargo mismatch")
    spec_sha = base.sha256(SPEC_PATH)
    split_sha = base.sha256(SPLIT_REGISTRY_PATH)
    selected_gpus, snapshots = base.wait_for_gpus(
        spec["gpu_policy"]["minimum_free_mib"]
    )
    formal = formal_jobs(spec)
    ablations = ablation_jobs(spec)
    execution: dict[str, Any] = {
        "schema_version": 1,
        "phase": "P6B",
        "revision_id": spec["revision_id"],
        "split_revision": spec["split_revision"],
        "status": "running",
        "started_at": utc_now(),
        "completed_at": None,
        "spec_path": str(SPEC_PATH.relative_to(REPO_ROOT)),
        "spec_sha256": spec_sha,
        "split_registry_sha256": split_sha,
        "means_sha256": base.sha256(MEANS_PATH),
        "selected_physical_gpus": selected_gpus,
        "gpu_wait_snapshots": snapshots,
        "jobs_expected": len(formal) + len(ablations),
        "jobs": [],
        "failures": [],
        "locked_test_block_signal_reads": 0,
    }
    base.atomic_json(EXECUTION_PATH, execution)
    try:
        formal_records = shared.run_queue(
            formal, selected_gpus, spec, spec_sha, execution, EXECUTION_PATH
        )
        ablation_records = shared.run_queue(
            ablations, selected_gpus, spec, spec_sha, execution, EXECUTION_PATH
        )
        formal_rows = shared.aggregate(formal_records, ("model", "loss"))
        ablation_rows = shared.aggregate(
            ablation_records, ("ablation_id", "model", "loss")
        )
        shared.write_tsv(RESULTS_PATH, formal_rows)
        shared.write_tsv(ABLATION_RESULTS_PATH, ablation_rows)
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
            "schema_version": 1,
            "revision_id": spec["revision_id"],
            "split_revision": spec["split_revision"],
            "selected_at": utc_now(),
            "selected": selected,
            "spec_sha256": spec_sha,
            "results_path": str(RESULTS_PATH.relative_to(REPO_ROOT)),
            "results_sha256": base.sha256(RESULTS_PATH),
            "selection_rule": spec["selection"],
            "locked_test_block_signal_reads": 0,
        }
        base.atomic_json(SELECTION_PATH, selection)
        development = run_development(selected, selected_gpus[0], spec, spec_sha)
        final_lock = {
            "schema_version": 3,
            "model": selected["model"],
            "loss": selected["loss"],
            "seed": development["seed"],
            "steps": development["steps"],
            "checkpoint_path": development["checkpoint_path"],
            "checkpoint_sha256": development["checkpoint_sha256"],
            "selection_path": str(SELECTION_PATH.relative_to(REPO_ROOT)),
            "selection_sha256": base.sha256(SELECTION_PATH),
            "spec_sha256": spec_sha,
            "split_revision": spec["split_revision"],
            "split_registry_sha256": split_sha,
            "means_sha256": base.sha256(MEANS_PATH),
            "test_intervals_sha256": base.sha256(TEST_INTERVALS_PATH),
            "training_chromosomes": ["I", "II", "III", "IV", "V", "X"],
            "final_test_chromosomes": ["I", "II", "III", "IV", "V", "X"],
            "final_test_scope": "r6c_single_six_chromosome_block_test",
            "legacy_chr_x_prior_exposure_disclosed": True,
            "locked_at": utc_now(),
            "superseded": False,
            "test_consumed": False,
            "test_consumed_at": None,
        }
        base.atomic_json(LOCK_PATH, final_lock)
        execution.update(
            {
                "status": "completed",
                "completed_at": utc_now(),
                "formal_jobs": formal_records,
                "ablation_jobs": ablation_records,
                "development": development,
                "results_sha256": base.sha256(RESULTS_PATH),
                "ablation_results_sha256": base.sha256(ABLATION_RESULTS_PATH),
                "selection_sha256": base.sha256(SELECTION_PATH),
                "lock_sha256": base.sha256(LOCK_PATH),
                "resources_after": v2_gpu_resources.snapshot(),
            }
        )
    except Exception:
        execution["status"] = "failed"
        execution["completed_at"] = utc_now()
        base.atomic_json(EXECUTION_PATH, execution)
        raise
    base.atomic_json(EXECUTION_PATH, execution)


if __name__ == "__main__":
    main()
