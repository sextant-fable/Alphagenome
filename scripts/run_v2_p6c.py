#!/usr/bin/env python3
"""Run the separately approved one-time six-chromosome block evaluation."""

from __future__ import annotations

import csv
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import subprocess
import uuid

from scripts import run_v2_p6b
from scripts import v2_gpu_resources
from scripts import v2_subprocess


REPO_ROOT = Path(__file__).resolve().parents[1]
METADATA_DIR = REPO_ROOT / "alphagenome_custom/metadata/v2"
STATE_PATH = METADATA_DIR / "execution_state.json"
LOCK_PATH = METADATA_DIR / "final_test_lock.json"
CLAIM_PATH = METADATA_DIR / "final_test_claim.json"
REPORT_PATH = METADATA_DIR / "final_test_report.json"
EXECUTION_PATH = METADATA_DIR / "p6c_execution.json"
FINAL_SUMMARY_PATH = METADATA_DIR / "v2_final_report.json"
FINAL_DOCUMENT_PATH = REPO_ROOT / "docs/v2_final_report.md"
SPLIT_REGISTRY_PATH = METADATA_DIR / "split_registry_v2.json"
MEANS_PATH = METADATA_DIR / "track_nonzero_means_v2.tsv"
MEANS_SUMMARY_PATH = METADATA_DIR / "track_nonzero_means_v2_summary.json"
TEST_INTERVALS_PATH = REPO_ROOT / "alphagenome_custom/intervals/v2/test_locked.tsv"


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def atomic_json(path: Path, payload: dict[str, object]) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(8 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def repository_path(value: object, label: str) -> Path:
    raw = Path(str(value))
    path = raw.resolve() if raw.is_absolute() else (REPO_ROOT / raw).resolve()
    try:
        path.relative_to(REPO_ROOT.resolve())
    except ValueError as error:
        raise RuntimeError(f"Locked {label} escapes the repository: {raw}") from error
    if not path.is_file():
        raise RuntimeError(f"Locked {label} is missing: {raw}")
    return path


def validate_locked_inputs(lock: dict[str, object]) -> dict[str, str]:
    checkpoint = repository_path(lock.get("checkpoint_path"), "checkpoint")
    selection = repository_path(lock.get("selection_path"), "selection")
    observed_checkpoint = sha256(checkpoint)
    observed_selection = sha256(selection)
    if observed_checkpoint != lock.get("checkpoint_sha256"):
        raise RuntimeError("Locked checkpoint SHA-256 mismatch before final-test claim")
    if observed_selection != lock.get("selection_sha256"):
        raise RuntimeError("Locked selection SHA-256 mismatch before final-test claim")

    split_registry = json.loads(SPLIT_REGISTRY_PATH.read_text())
    if not (
        split_registry.get("revision_id") == "six_chromosome_blocks_v1"
        and split_registry.get("test_chromosomes") == ["I", "II", "III", "IV", "V", "X"]
        and split_registry.get("final_test_scope")
        == "r6c_single_six_chromosome_block_test"
    ):
        raise RuntimeError("Final test requires the locked six-chromosome split")
    relative_test = str(TEST_INTERVALS_PATH.relative_to(REPO_ROOT))
    observed_test = sha256(TEST_INTERVALS_PATH)
    if split_registry.get("files", {}).get(relative_test) != observed_test:
        raise RuntimeError("Final-test intervals no longer match the preregistered split")
    means_summary = json.loads(MEANS_SUMMARY_PATH.read_text())
    observed_means = sha256(MEANS_PATH)
    if means_summary.get("output_sha256") != observed_means:
        raise RuntimeError("Development means no longer match their preregistered summary")
    if means_summary.get("split_registry_sha256") != sha256(SPLIT_REGISTRY_PATH):
        raise RuntimeError("Development means were not derived from the current split registry")
    if not (
        lock.get("schema_version") == 3
        and lock.get("split_revision") == split_registry.get("revision_id")
        and lock.get("split_registry_sha256") == sha256(SPLIT_REGISTRY_PATH)
        and lock.get("means_sha256") == observed_means
        and lock.get("test_intervals_sha256") == observed_test
        and lock.get("training_chromosomes")
        == ["I", "II", "III", "IV", "V", "X"]
        and lock.get("final_test_chromosomes")
        == ["I", "II", "III", "IV", "V", "X"]
        and lock.get("final_test_scope")
        == "r6c_single_six_chromosome_block_test"
    ):
        raise RuntimeError("Final-test lock does not bind the six-chromosome inputs")
    return {
        "checkpoint_sha256": observed_checkpoint,
        "selection_sha256": observed_selection,
        "test_intervals_sha256": observed_test,
        "means_sha256": observed_means,
        "split_registry_sha256": sha256(SPLIT_REGISTRY_PATH),
    }


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def tracked_artifact(path: Path) -> dict[str, object]:
    return {
        "path": str(path.relative_to(REPO_ROOT)),
        "sha256": sha256(path),
        "bytes": path.stat().st_size,
    }


def build_final_summary(
    lock: dict[str, object],
    execution_id: str,
    command: list[str],
    report: dict[str, object],
) -> dict[str, object]:
    p1_path = METADATA_DIR / "p1_summary.json"
    p2_path = METADATA_DIR / "p2_summary.json"
    p3_full_path = METADATA_DIR / "p3_full_summary.json"
    p3_group_path = METADATA_DIR / "p3_group_summary.json"
    p4_path = METADATA_DIR / "p4_loader_benchmark.json"
    p5_path = METADATA_DIR / "p5_component_audit.json"
    p6a_path = METADATA_DIR / "p6a_execution.json"
    p6b_execution_path = METADATA_DIR / "p6b_six_chromosome_execution.json"
    p6b_results_path = METADATA_DIR / "p6b_six_chromosome_cv_results.tsv"
    p6b_ablations_path = METADATA_DIR / "p6b_six_chromosome_ablation_results.tsv"
    p6b_review_path = METADATA_DIR / "audits/P6B/review.json"
    selection_path = repository_path(lock["selection_path"], "selection")
    sample_manifest_path = METADATA_DIR / "rna_seq_samples_v2.tsv"
    group_manifest_path = METADATA_DIR / "rna_seq_groups_v2_final.tsv"
    group_outputs_path = METADATA_DIR / "p3_group_outputs.tsv"
    artifacts = [
        p1_path,
        p2_path,
        p3_full_path,
        p3_group_path,
        p4_path,
        p5_path,
        p6a_path,
        p6b_execution_path,
        p6b_results_path,
        p6b_ablations_path,
        p6b_review_path,
        selection_path,
        sample_manifest_path,
        group_manifest_path,
        group_outputs_path,
        SPLIT_REGISTRY_PATH,
        MEANS_PATH,
        TEST_INTERVALS_PATH,
        REPORT_PATH,
    ]
    if not all(path.is_file() for path in artifacts):
        missing = [str(path) for path in artifacts if not path.is_file()]
        raise RuntimeError(f"Cannot build final summary; missing artifacts: {missing}")

    p1 = json.loads(p1_path.read_text())
    p2 = json.loads(p2_path.read_text())
    p3_full = json.loads(p3_full_path.read_text())
    p3_group = json.loads(p3_group_path.read_text())
    p4 = json.loads(p4_path.read_text())
    p5 = json.loads(p5_path.read_text())
    p6b_execution = json.loads(p6b_execution_path.read_text())
    selection = json.loads(selection_path.read_text())
    cv_results = read_tsv(p6b_results_path)
    ablations = read_tsv(p6b_ablations_path)
    split_registry = json.loads(SPLIT_REGISTRY_PATH.read_text())
    git_commit = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=REPO_ROOT,
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()
    full_metrics = report["full_metrics"]
    return {
        "schema_version": 1,
        "workflow": "c_elegans_rna_seq_v2",
        "status": "completed_one_time_final_test",
        "generated_at": utc_now(),
        "execution_id": execution_id,
        "code_commit": git_commit,
        "host": os.uname().nodename,
        "data_lineage": {
            "reference_genome": "WBcel235",
            "input_samples": p1["sample_count"],
            "batch_counts": p1["batch_counts"],
            "source_normalization_classes": p1["normalization_class_counts"],
            "rna_seq_runs_reprocessed": p3_full["completed_runs"],
            "non_rna_runs_excluded": p2["non_rna_excluded"],
            "source_fastq_bytes": p3_full["source_fastq_bytes"],
            "normalized_runs": p3_group["normalized_runs"],
            "formal_groups": p3_group["formal_groups"],
            "grouped_tracks": p3_group["group_outputs"],
            "biological_replicate_groups": p3_group[
                "biological_replicate_groups"
            ],
            "normalization_target_total": p3_group["normalization_target_total"],
            "normalization_interpretation": "1e6 times 100-bp coverage",
            "aggregation_policy": p3_group["aggregation_policy"],
            "dynamic_loader_tracks": p4["full_window_tracks"],
            "monolithic_npz_generated": p4["monolithic_npz_generated"],
            "cv_chromosomes": split_registry["cv_chromosomes"],
            "test_chromosomes": split_registry["test_chromosomes"],
            "test_windows": split_registry["test_windows"],
        },
        "method_contract": {
            "loss_implementation": p5["loss_implementation"],
            "scale_implementation": p5["scale_implementation"],
            "worm_embedding_rows": p5["worm_embedding_rows"],
            "lora_trainable_parameters": p5["lora_trainable_parameters"],
            "blocked_cv": "five-fold within-chromosome block validation over I-V,X",
            "model_classes": ["A", "B", "C"],
            "formal_seeds": 3,
        },
        "formal_validation": {
            "matrix": cv_results,
            "ablations": ablations,
            "selected": selection["selected"],
            "selection_rule": selection["selection_rule"],
            "formal_jobs": len(p6b_execution["formal_jobs"]),
            "ablation_jobs": len(p6b_execution["ablation_jobs"]),
            "registered_failures": p6b_execution.get("failures", []),
            "development_checkpoint": p6b_execution["development"],
        },
        "final_test": {
            "checkpoint_path": lock["checkpoint_path"],
            "checkpoint_sha256": lock["checkpoint_sha256"],
            "model": lock["model"],
            "loss": lock["loss"],
            "seed": lock["seed"],
            "steps": lock["steps"],
            "primary_metrics": full_metrics["primary"],
            "mean_losses": report["mean_metrics"],
            "validation_subwindows": report["validation_subwindows"],
            "validation_bases": report["validation_bases"],
            "core_coverage_fraction": report["validation_core_coverage_fraction"],
            "report_path": str(REPORT_PATH.relative_to(REPO_ROOT)),
            "report_sha256": sha256(REPORT_PATH),
            "test_consumed_once": True,
            "no_post_test_tuning": True,
        },
        "limitations": [
            "Chromosome X was viewed in legacy exploratory work, so its locked test block is not fully pristine; every final test block remains one-time-use within the revised v2 workflow.",
            "All 485 source runs were class C before reprocessing; 482 RNA-seq runs were uniformly rebuilt and three non-RNA runs were excluded.",
            "The augmentation and gene-loss ablations used one seed across five folds and are directional rather than high-power estimates.",
            f"The selected {selection['selected']['model']}/{selection['selected']['loss']} configuration won the preregistered biological primary score but need not dominate every secondary metric.",
            "The one-time six-chromosome locked-block result must be reported regardless of outcome and cannot be used for further tuning in this study cycle.",
        ],
        "reproduction": {
            "phase_command": command,
            "controller_command": [
                str(command[0]),
                "-m",
                "scripts.v2_phase_controller",
                "run",
                "--phase",
                "P6C",
                "--auto",
            ],
            "cuda_visible_devices": report["physical_cuda_visible_devices"],
            "python_environment": str(command[0]),
        },
        "artifacts": {
            str(path.relative_to(REPO_ROOT)): tracked_artifact(path)
            for path in artifacts
        },
    }


def write_final_document(summary: dict[str, object]) -> None:
    data = summary["data_lineage"]
    selected = summary["formal_validation"]["selected"]
    final = summary["final_test"]
    primary = final["primary_metrics"]
    lines = [
        "# C. elegans RNA-seq v2 Final Report",
        "",
        f"Status: `{summary['status']}`",
        "",
        f"Execution ID: `{summary['execution_id']}`",
        "",
        "## Data",
        "",
        f"The auditable input set contains {data['input_samples']} runs: "
        f"{data['batch_counts']}. Uniform reprocessing retained "
        f"{data['rna_seq_runs_reprocessed']} RNA-seq runs, excluded "
        f"{data['non_rna_runs_excluded']} non-RNA runs, and produced "
        f"{data['grouped_tracks']} grouped tracks at a target total of "
        f"{data['normalization_target_total']} ({data['normalization_interpretation']}).",
        "",
        "The model reads these tracks on demand. No monolithic v2 NPZ was generated. "
        "Each of I, II, III, IV, V, and X contributes leakage-buffered train, "
        "validation, and locked final-test blocks. The final blocks were consumed once "
        "by the locked final execution.",
        "",
        "## Model Selection",
        "",
        f"The formal matrix selected model {selected['model']} with "
        f"{selected['loss']} loss over {selected['folds']} folds, "
        f"{selected['seeds']} seeds, and {selected['jobs']} jobs. Its development "
        f"primary biological score was {selected['mean_primary_biological_score']:.6f}.",
        "",
        "## Final Test",
        "",
        f"Gene-exon coverage Pearson (log1p): "
        f"{primary['mean_per_track_gene_exon_coverage_pearson_log1p']:.6f}",
        "",
        f"Per-track 128-bp Pearson (log1p): "
        f"{primary['mean_per_track_pearson_128bp_log1p']:.6f}",
        "",
        f"Paper loss: {final['mean_losses']['paper_loss']:.6f}; "
        f"log1p MSE: {final['mean_losses']['log1p_mse']:.6f}.",
        "",
        "## Limitations",
        "",
    ]
    lines.extend(f"- {item}" for item in summary["limitations"])
    lines.extend(
        [
            "",
            "## Reproduction",
            "",
            "The exact command, code commit, input hashes, complete metric report, "
            "formal matrix, ablations, and failure record are stored in "
            "`alphagenome_custom/metadata/v2/v2_final_report.json`.",
            "",
        ]
    )
    temporary = FINAL_DOCUMENT_PATH.with_suffix(FINAL_DOCUMENT_PATH.suffix + ".tmp")
    temporary.write_text("\n".join(lines))
    temporary.replace(FINAL_DOCUMENT_PATH)


def require_approvals() -> tuple[dict[str, object], dict[str, object]]:
    state = json.loads(STATE_PATH.read_text())
    g4 = state.get("approvals", {}).get("G4_gpu_experiments", {})
    g5 = state.get("approvals", {}).get("G5_final_test", {})
    if not (
        state.get("current_phase") == "P6C"
        and g4.get("approved") is True
        and g4.get("scope") == "r6_gpu_auto_available_2_3"
        and g5.get("approved") is True
        and g5.get("scope") == "r6c_single_six_chromosome_block_test"
    ):
        raise RuntimeError("P6C requires current phase P6C and exact scoped G4/G5")
    lock = json.loads(LOCK_PATH.read_text())
    if lock.get("superseded") is True:
        raise RuntimeError("Final-test lock was superseded; amended R6B must issue a new lock")
    if lock.get("test_consumed") is True:
        raise RuntimeError("Six-chromosome final test has already been consumed")
    return lock, state


def main() -> None:
    lock, state = require_approvals()
    locked_hashes = validate_locked_inputs(lock)
    selected, snapshots = run_v2_p6b.wait_for_gpus(70_000)
    gpu = selected[0]
    execution_id = (
        "p6c-"
        + datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
        + "-"
        + uuid.uuid4().hex[:12]
    )
    g5 = state["approvals"]["G5_final_test"]
    claim = {
        "schema_version": 1,
        "execution_id": execution_id,
        "claimed_at": utc_now(),
        "g5_approved_at": g5.get("updated_at"),
        "g5_scope": g5.get("scope"),
        "checkpoint_path": lock["checkpoint_path"],
        "checkpoint_sha256": lock["checkpoint_sha256"],
        "selection_path": lock["selection_path"],
        **locked_hashes,
        "physical_gpu": gpu,
        "status": "claimed",
    }
    try:
        descriptor = os.open(CLAIM_PATH, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o644)
    except FileExistsError as error:
        raise RuntimeError("Final-test claim already exists; refusing a second entry") from error
    with os.fdopen(descriptor, "w") as handle:
        handle.write(json.dumps(claim, indent=2, sort_keys=True) + "\n")
    command = v2_subprocess.module_command(
        "evaluate_v2_model",
        "--model",
        str(lock["model"]),
        "--training-loss",
        str(lock["loss"]),
        "--fold",
        "0",
        "--checkpoint",
        str(lock["checkpoint_path"]),
        "--sequence-length",
        "131072",
        "--hidden-channels",
        "64",
        "--output",
        str(REPORT_PATH.relative_to(REPO_ROOT)),
        "--device",
        "cuda",
        "--final-test",
        "--final-test-execution-id",
        execution_id,
    )
    log_path = REPO_ROOT / "logs" / execution_id / "final_test.log"
    record = {
        "schema_version": 1,
        "phase": "P6C",
        "execution_id": execution_id,
        "status": "running",
        "started_at": utc_now(),
        "completed_at": None,
        "g5_approved_at": g5.get("updated_at"),
        "physical_gpu": gpu,
        "checkpoint_path": lock["checkpoint_path"],
        "checkpoint_sha256": lock["checkpoint_sha256"],
        "selection_path": lock["selection_path"],
        "locked_input_hashes": locked_hashes,
        "gpu_wait_snapshots": snapshots,
        "command": command,
        "claim_path": str(CLAIM_PATH.relative_to(REPO_ROOT)),
        "report_path": str(REPORT_PATH.relative_to(REPO_ROOT)),
        "log_path": str(log_path.relative_to(REPO_ROOT)),
    }
    atomic_json(EXECUTION_PATH, record)
    environment = os.environ.copy()
    environment["CUDA_VISIBLE_DEVICES"] = str(gpu)
    log_path.parent.mkdir(parents=True, exist_ok=True)
    try:
        with log_path.open("w") as log:
            subprocess.run(
                command,
                cwd=REPO_ROOT,
                env=environment,
                stdout=log,
                stderr=subprocess.STDOUT,
                check=True,
            )
        report = json.loads(REPORT_PATH.read_text())
        if not (
            report.get("final_test_execution_id") == execution_id
            and report.get("checkpoint_sha256") == lock["checkpoint_sha256"]
        ):
            raise RuntimeError("Final-test report does not match its exclusive claim")
        summary = build_final_summary(lock, execution_id, command, report)
        atomic_json(FINAL_SUMMARY_PATH, summary)
        write_final_document(summary)
        current_lock = json.loads(LOCK_PATH.read_text())
        if not (
            current_lock.get("test_status") == "completed"
            and current_lock.get("test_execution_id") == execution_id
        ):
            raise RuntimeError("Final-test lock was not finalized before summary creation")
        current_lock["final_summary_path"] = str(
            FINAL_SUMMARY_PATH.relative_to(REPO_ROOT)
        )
        current_lock["final_summary_sha256"] = sha256(FINAL_SUMMARY_PATH)
        current_lock["final_document_path"] = str(
            FINAL_DOCUMENT_PATH.relative_to(REPO_ROOT)
        )
        current_lock["final_document_sha256"] = sha256(FINAL_DOCUMENT_PATH)
        atomic_json(LOCK_PATH, current_lock)
        record["final_summary_path"] = current_lock["final_summary_path"]
        record["final_summary_sha256"] = current_lock["final_summary_sha256"]
        record["final_document_path"] = current_lock["final_document_path"]
        record["final_document_sha256"] = current_lock["final_document_sha256"]
    except Exception:
        record["status"] = "failed"
        current_lock = json.loads(LOCK_PATH.read_text())
        if current_lock.get("test_consumed") is True:
            current_lock["test_status"] = "failed_after_consumption"
            atomic_json(LOCK_PATH, current_lock)
        claim["status"] = "failed"
        atomic_json(CLAIM_PATH, claim)
        raise
    else:
        record["status"] = "completed"
        claim["status"] = "completed"
        claim["final_summary_sha256"] = record["final_summary_sha256"]
        atomic_json(CLAIM_PATH, claim)
    finally:
        record["completed_at"] = utc_now()
        record["resources_after"] = v2_gpu_resources.snapshot()
        atomic_json(EXECUTION_PATH, record)


if __name__ == "__main__":
    main()
