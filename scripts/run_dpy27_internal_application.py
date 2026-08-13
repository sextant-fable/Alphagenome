#!/usr/bin/env python3
"""Run the P10 DPY-27 internal genomic-block validation application.

The no-argument entry point is controlled by P10.  It validates every frozen
input before selecting a GPU, runs the 15 existing B/paper CV checkpoints on
their matching validation cores, and supports checksum-verified restart shards.
It never resolves or reads the final-test interval manifest.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
import json
import math
import os
from pathlib import Path
import socket
import subprocess
import sys
import time
from typing import Any, Mapping

import numpy as np

from scripts import v2_biological_validation as biology
from scripts import v2_gpu_resources


REPO_ROOT = Path(__file__).resolve().parents[1]
METADATA_DIR = REPO_ROOT / "alphagenome_custom/metadata/v2"
STATE_PATH = METADATA_DIR / "execution_state.json"
SPEC_PATH = METADATA_DIR / "p10_dpy27_internal_application_spec.json"
ANALYSIS_CONTRACT = "dpy27_internal_genomic_block_validation_v1"
ANALYSIS_LABEL = "internal_genomic_block_validation"
EXPECTED_MODEL = "B"
EXPECTED_LOSS = "paper"
EXPECTED_SEEDS = (20260714, 20260715, 20260716)
EXPECTED_FOLDS = (1, 2, 3, 4, 5)
CASE_GROUP_ID = "RNA_V2_G0005"
CONTROL_GROUP_ID = "RNA_V2_G0006"
DEFAULT_RUN_ROOT = "runs/v2_p6b_six_chromosome_corefix_20260722/formal/B/paper"
CONTROLLER_INVOCATION = "python -m scripts.v2_phase_controller run --phase P10"
REGISTERED_PHASE_COMMAND = "python -m scripts.run_dpy27_internal_application"
IMPLEMENTATION_PATHS = (
    "scripts/run_dpy27_internal_application.py",
    "scripts/v2_biological_validation.py",
    "scripts/plot_dpy27_internal_application.py",
    "scripts/train_v2_model.py",
    "scripts/v2_training_components.py",
    "scripts/v2_bigwig_dataset.py",
)
GENE_FIELDS = (
    "analysis_label",
    "fold",
    "seed",
    "gene_id",
    "gene_name",
    "chromosome",
    "strand",
    "block_id",
    "exon_bases_evaluated",
    "exon_bases_total",
    "exon_coverage_fraction",
    "complete_exon_coverage",
    "predicted_dpy27_rnai_mean",
    "predicted_vector_rnai_mean",
    "observed_dpy27_rnai_mean",
    "observed_vector_rnai_mean",
    "predicted_log1p_contrast_g0005_minus_g0006",
    "observed_log1p_contrast_g0005_minus_g0006",
)


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def read_json(path: str | Path) -> dict[str, Any]:
    payload = json.loads(Path(path).read_text())
    if not isinstance(payload, dict):
        raise ValueError(f"Expected a JSON object: {path}")
    return payload


def implementation_sha256() -> dict[str, str]:
    return {
        relative: biology.sha256(REPO_ROOT / relative)
        for relative in IMPLEMENTATION_PATHS
    }


def git_commit() -> str:
    return subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=REPO_ROOT,
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()


def execution_provenance() -> dict[str, str]:
    """Return the stable controller/host fields required in the P10 audit."""

    return {
        "git_commit": git_commit(),
        "hostname": socket.gethostname(),
        "working_directory": str(REPO_ROOT),
        "python_executable": sys.executable,
        "controller_invocation": CONTROLLER_INVOCATION,
        "registered_phase_command": REGISTERED_PHASE_COMMAND,
    }


def repo_path(relative: str) -> Path:
    path = Path(relative)
    if path.is_absolute() or ".." in path.parts:
        raise ValueError(f"Spec paths must be repository-relative: {relative}")
    return REPO_ROOT / path


def _validate_static_contract(spec: Mapping[str, Any]) -> None:
    model = spec.get("model_contract", {})
    statistics = spec.get("statistics", {})
    figure_source = spec.get("figure_source_data_contract", {})
    if not (
        spec.get("schema_version") == 1
        and spec.get("phase") == "P10"
        and spec.get("analysis_contract") == ANALYSIS_CONTRACT
        and spec.get("analysis_label") == ANALYSIS_LABEL
        and spec.get("final_test_access") == "prohibited"
        and tuple(spec.get("expected_folds", [])) == EXPECTED_FOLDS
        and tuple(spec.get("expected_seeds", [])) == EXPECTED_SEEDS
        and model.get("model") == EXPECTED_MODEL
        and model.get("loss") == EXPECTED_LOSS
        and model.get("sequence_length") == 131072
        and model.get("hidden_channels") == 64
        and statistics.get("case_group_id") == CASE_GROUP_ID
        and statistics.get("control_group_id") == CONTROL_GROUP_ID
        and statistics.get("inference_unit")
        == "fold_after_within_fold_gene_condition_seed_mean"
        and statistics.get("complete_exon_coverage_required") is True
    ):
        raise RuntimeError("P10 spec does not match the frozen DPY-27 contract")
    gpu_policy = spec.get("gpu_policy", {})
    if gpu_policy.get("allowed_physical_indices") != [2, 3]:
        raise RuntimeError("P10 must restrict physical GPUs to 2 and 3")
    if gpu_policy.get("maximum_concurrent_jobs_per_gpu") != 1:
        raise RuntimeError("P10 must allow at most one worker per physical GPU")
    expected_panels = [
        "a_gene_contrasts",
        "b_fold_x_minus_autosomes",
        "c_fold_agreement",
        "c_fold_primary_summary",
    ]
    expected_panel_endpoints = {
        "a_gene_contrasts": ["gene_log1p_contrast"],
        "b_fold_x_minus_autosomes": [
            "observed_median_contrast_x_minus_autosomes",
            "predicted_median_contrast_x_minus_autosomes",
        ],
        "c_fold_agreement": [
            "x_predicted_observed_direction_concordance",
            "predicted_observed_gene_spearman",
        ],
        "c_fold_primary_summary": list(biology.PRIMARY_ENDPOINTS),
    }
    if (
        figure_source.get("panels") != expected_panels
        or figure_source.get("panel_endpoints") != expected_panel_endpoints
    ):
        raise RuntimeError("P10 Source Data panel contract changed")
    outputs = spec.get("output_paths", {})
    required_outputs = {
        "analysis_audit",
        "preflight_audit",
        "run_root",
        "gene_contrasts",
        "fold_gene_contrasts",
        "run_endpoints",
        "fold_endpoints",
        "summary_endpoints",
        "figure_source_data",
        "figure_pdf",
        "figure_png",
        "figure_qa",
        "figure_svg",
        "figure_tiff",
    }
    if required_outputs - set(outputs):
        raise RuntimeError("P10 output path contract is incomplete")
    for path in outputs.values():
        biology.assert_no_final_test_path(path)
        repo_path(str(path))
    locked = spec.get("locked_inputs", {})
    for key, path in locked.items():
        lowered = str(path).lower()
        if key not in {"final_test_lock", "final_test_report"} and any(
            token in lowered for token in biology.FORBIDDEN_SPLIT_TOKENS
        ):
            raise RuntimeError(f"P10 data input cannot name final-test content: {key}")
    registered_implementation = spec.get("implementation_sha256", {})
    if (
        set(registered_implementation) != set(IMPLEMENTATION_PATHS)
        or registered_implementation != implementation_sha256()
    ):
        raise RuntimeError("P10 frozen implementation hash mismatch")


def _validate_controller_state() -> None:
    state = read_json(STATE_PATH)
    approval = state.get("approvals", {}).get("G4_gpu_experiments", {})
    if not (
        state.get("current_phase") == "P10"
        and state.get("status") in {"PENDING", "RUNNING"}
        and approval.get("approved") is True
        and approval.get("scope") == "p10_dpy27_internal_application"
    ):
        raise RuntimeError("P10 requires its exact active controller phase and G4 scope")


def _validate_locked_inputs(spec: Mapping[str, Any]) -> dict[str, str]:
    locked = spec.get("locked_inputs", {})
    expected = spec.get("input_sha256", {})
    if set(locked) != set(expected):
        raise RuntimeError("P10 locked_inputs and input_sha256 keys differ")
    verified: dict[str, str] = {}
    for key, relative in locked.items():
        path = repo_path(str(relative))
        if not path.is_file():
            raise FileNotFoundError(f"Missing P10 input {key}: {path}")
        actual = biology.sha256(path)
        if actual != expected[key]:
            raise RuntimeError(f"P10 input hash mismatch: {key}")
        verified[key] = actual
    lock = read_json(repo_path(str(locked["final_test_lock"])))
    if not (lock.get("test_consumed") is True and lock.get("test_status") == "completed"):
        raise RuntimeError("P10 must preserve the already-completed final-test lock")
    return verified


def _validate_checkpoint_matrix(spec: Mapping[str, Any]) -> list[dict[str, Any]]:
    matrix = spec.get("checkpoint_matrix")
    if not isinstance(matrix, list) or len(matrix) != 15:
        raise RuntimeError("P10 requires exactly 15 frozen checkpoints")
    expected_pairs = {(seed, fold) for seed in EXPECTED_SEEDS for fold in EXPECTED_FOLDS}
    observed_pairs: set[tuple[int, int]] = set()
    validated: list[dict[str, Any]] = []
    model = spec["model_contract"]
    for raw in matrix:
        row = dict(raw)
        seed, fold = int(row["seed"]), int(row["fold"])
        pair = (seed, fold)
        if pair in observed_pairs:
            raise RuntimeError(f"Duplicate P10 checkpoint pair: {pair}")
        observed_pairs.add(pair)
        checkpoint_path = repo_path(str(row["checkpoint_path"]))
        run_path = repo_path(str(row["run_path"]))
        if biology.sha256(checkpoint_path) != row["checkpoint_sha256"]:
            raise RuntimeError(f"Checkpoint hash mismatch for seed={seed}, fold={fold}")
        if biology.sha256(run_path) != row["run_sha256"]:
            raise RuntimeError(f"run.json hash mismatch for seed={seed}, fold={fold}")
        run = read_json(run_path)
        expected_run = {
            "fold": fold,
            "seed": seed,
            "loss": EXPECTED_LOSS,
            "sequence_length": model["sequence_length"],
            "hidden_channels": model["hidden_channels"],
            "mean_column": model["mean_column_template"].format(fold=fold),
            "checkpoint_sha256": row["checkpoint_sha256"],
            "checkpoint_reload_verified": True,
        }
        mismatches = {
            key: (run.get(key), value)
            for key, value in expected_run.items()
            if run.get(key) != value
        }
        if run.get("model") != EXPECTED_MODEL and run.get("model_id") != EXPECTED_MODEL:
            mismatches["model"] = (run.get("model", run.get("model_id")), EXPECTED_MODEL)
        if mismatches:
            raise RuntimeError(f"Frozen run contract mismatch for {pair}: {mismatches}")
        validated.append(row)
    if observed_pairs != expected_pairs:
        raise RuntimeError("P10 checkpoint matrix is not the exact 3 seeds x 5 folds product")
    return sorted(validated, key=lambda row: (int(row["seed"]), int(row["fold"])))


def _validate_tracks_and_intervals(spec: Mapping[str, Any]) -> dict[str, Any]:
    locked = spec["locked_inputs"]
    group_rows = biology.read_tsv(repo_path(locked["group_manifest"]))
    by_group = {row["group_id"]: row for row in group_rows}
    if CASE_GROUP_ID not in by_group or CONTROL_GROUP_ID not in by_group:
        raise RuntimeError("DPY-27 case/control groups are missing from the formal manifest")
    case = by_group[CASE_GROUP_ID]
    control = by_group[CONTROL_GROUP_ID]
    if not (
        case.get("include_formal_v2") == "True"
        and control.get("include_formal_v2") == "True"
        and case.get("sra_study_accession") == control.get("sra_study_accession")
        and case.get("development_stage") == control.get("development_stage")
        and case.get("normalization_target_total") == "100000000"
        and control.get("normalization_target_total") == "100000000"
    ):
        raise RuntimeError("DPY-27 case/control metadata no longer match the frozen contrast")
    track_rows = biology.read_tsv(repo_path(locked["track_manifest"]))
    track_ids = [row["group_id"] for row in track_rows]
    if len(track_ids) != 241 or len(set(track_ids)) != 241:
        raise RuntimeError("Expected exactly 241 unique formal v2 tracks")
    case_index, control_index = track_ids.index(CASE_GROUP_ID), track_ids.index(CONTROL_GROUP_ID)
    for group_id, key, index in (
        (CASE_GROUP_ID, "case_bigwig", case_index),
        (CONTROL_GROUP_ID, "control_bigwig", control_index),
    ):
        row = track_rows[index]
        if repo_path(row["output_path"]) != repo_path(locked[key]):
            raise RuntimeError(f"Frozen BigWig path mismatch for {group_id}")
        if row.get("output_sha256") != spec["input_sha256"][key]:
            raise RuntimeError(f"Manifest BigWig hash mismatch for {group_id}")
    means_rows = biology.read_tsv(repo_path(locked["means"]))
    if [row["group_id"] for row in means_rows] != track_ids:
        raise RuntimeError("Track means order differs from model output track order")
    intervals: dict[int, dict[str, Any]] = {}
    for fold in EXPECTED_FOLDS:
        relative = locked[f"fold_{fold}_valid"]
        rows, subwindows = biology.load_fold_validation_subwindows(
            repo_path(relative),
            expected_fold=fold,
            sequence_length=int(spec["model_contract"]["sequence_length"]),
        )
        if len(subwindows) != int(spec["expected_subwindows_per_fold"]):
            raise RuntimeError(f"Fold {fold} does not have the frozen subwindow count")
        intervals[fold] = {
            "path": relative,
            "sha256": spec["input_sha256"][f"fold_{fold}_valid"],
            "rows": len(rows),
            "subwindows": len(subwindows),
            "chromosome_counts": {
                chromosome: sum(window.chromosome == chromosome for window in subwindows)
                for chromosome in sorted({window.chromosome for window in subwindows})
            },
        }
    return {
        "case_track_index": case_index,
        "control_track_index": control_index,
        "case_metadata": case,
        "control_metadata": control,
        "track_count": len(track_ids),
        "intervals": intervals,
    }


def preflight_contract(
    spec_path: str | Path = SPEC_PATH,
    *,
    enforce_controller: bool = True,
) -> dict[str, Any]:
    """Validate the complete frozen P10 scope without model or BigWig reads."""

    spec_path = Path(spec_path)
    spec = read_json(spec_path)
    _validate_static_contract(spec)
    if enforce_controller:
        _validate_controller_state()
    locked_hashes = _validate_locked_inputs(spec)
    matrix = _validate_checkpoint_matrix(spec)
    dataset = _validate_tracks_and_intervals(spec)
    return {
        "schema_version": 1,
        "analysis_contract": ANALYSIS_CONTRACT,
        "analysis_label": ANALYSIS_LABEL,
        "status": "passed",
        "checked_at": utc_now(),
        "controller_enforced": enforce_controller,
        "spec_path": str(spec_path.relative_to(REPO_ROOT)),
        "spec_sha256": biology.sha256(spec_path),
        "locked_input_sha256": locked_hashes,
        "checkpoint_count": len(matrix),
        "checkpoint_matrix": matrix,
        "dataset": dataset,
        "final_test_access": "prohibited",
        "locked_test_block_signal_reads": 0,
        "model_or_bigwig_reads": 0,
        "git_commit": git_commit(),
        "implementation_sha256": implementation_sha256(),
    }


def _require_worker_gpu(physical_gpu: int) -> "Any":
    if physical_gpu not in {2, 3}:
        raise RuntimeError("P10 worker physical GPU must be 2 or 3")
    visible = os.environ.get("CUDA_VISIBLE_DEVICES")
    if visible != str(physical_gpu):
        raise RuntimeError(
            f"Worker CUDA_VISIBLE_DEVICES={visible!r} does not match physical GPU {physical_gpu}"
        )
    import torch

    if not torch.cuda.is_available() or torch.cuda.device_count() != 1:
        raise RuntimeError("P10 worker requires exactly one visible CUDA device")
    return torch


def _load_checkpoint(
    model: "Any",
    checkpoint_path: Path,
    *,
    fold: int,
    seed: int,
    spec: Mapping[str, Any],
) -> dict[str, Any]:
    import torch

    checkpoint = torch.load(checkpoint_path, map_location="cpu", weights_only=False)
    contract = spec["model_contract"]
    expected = {
        "model_id": EXPECTED_MODEL,
        "loss": EXPECTED_LOSS,
        "fold": fold,
        "seed": seed,
        "sequence_length": contract["sequence_length"],
        "hidden_channels": contract["hidden_channels"],
        "mean_column": contract["mean_column_template"].format(fold=fold),
    }
    mismatches = {
        key: (checkpoint.get(key), value)
        for key, value in expected.items()
        if checkpoint.get(key) != value
    }
    if mismatches:
        raise RuntimeError(f"Checkpoint tensor metadata mismatch: {mismatches}")
    trainable_names = {
        name for name, parameter in model.named_parameters() if parameter.requires_grad
    }
    trainable_state = checkpoint.get("trainable_model_state", {})
    if set(trainable_state) != trainable_names:
        raise RuntimeError("Checkpoint trainable parameter set mismatch")
    state = model.state_dict()
    for name, value in trainable_state.items():
        if name not in state or state[name].shape != value.shape:
            raise RuntimeError(f"Checkpoint tensor mismatch: {name}")
        state[name] = value
    model.load_state_dict(state, strict=True)
    return checkpoint


def _shard_paths(spec: Mapping[str, Any], *, seed: int, fold: int) -> dict[str, Path]:
    root = repo_path(spec["output_paths"]["run_root"]) / "shards" / f"seed_{seed}" / f"fold_{fold}"
    return {
        "root": root,
        "genes": root / "gene_contrasts.tsv",
        "endpoint": root / "run_endpoint.json",
        "audit": root / "shard_audit.json",
    }


def _completed_shard(
    spec: Mapping[str, Any],
    job: Mapping[str, Any],
    spec_sha256: str,
) -> tuple[list[dict[str, str]], dict[str, Any]] | None:
    paths = _shard_paths(spec, seed=int(job["seed"]), fold=int(job["fold"]))
    if not any(path.exists() for key, path in paths.items() if key != "root"):
        return None
    if not all(path.is_file() for key, path in paths.items() if key != "root"):
        raise RuntimeError(f"Incomplete restart shard: {paths['root']}")
    audit = read_json(paths["audit"])
    expected = {
        "status": "completed",
        "analysis_contract": ANALYSIS_CONTRACT,
        "spec_sha256": spec_sha256,
        "checkpoint_sha256": job["checkpoint_sha256"],
        "fold": int(job["fold"]),
        "seed": int(job["seed"]),
        "validation_sha256": spec["input_sha256"][f"fold_{int(job['fold'])}_valid"],
        "gene_contrasts_sha256": biology.sha256(paths["genes"]),
        "run_endpoint_sha256": biology.sha256(paths["endpoint"]),
        "final_test_access": "prohibited",
        "locked_test_block_signal_reads": 0,
        "implementation_sha256": implementation_sha256(),
    }
    mismatches = {
        key: (audit.get(key), value)
        for key, value in expected.items()
        if audit.get(key) != value
    }
    if mismatches:
        raise RuntimeError(f"Restart shard contract mismatch at {paths['root']}: {mismatches}")
    return biology.read_tsv(paths["genes"]), read_json(paths["endpoint"])


def run_one_checkpoint(
    spec: Mapping[str, Any],
    job: Mapping[str, Any],
    *,
    spec_sha256: str,
    physical_gpu: int,
) -> dict[str, Any]:
    """Run one frozen checkpoint; intended for the isolated worker process."""

    torch = _require_worker_gpu(physical_gpu)
    _validate_controller_state()
    if biology.sha256(SPEC_PATH) != spec_sha256:
        raise RuntimeError("P10 spec changed before worker execution")
    existing = _completed_shard(spec, job, spec_sha256)
    if existing is not None:
        genes, endpoint = existing
        return {"status": "reused", "gene_rows": len(genes), **endpoint}

    fold, seed = int(job["fold"]), int(job["seed"])
    checkpoint_path = repo_path(str(job["checkpoint_path"]))
    valid_path = repo_path(spec["locked_inputs"][f"fold_{fold}_valid"])
    # Revalidate every job-specific and model input before importing model/data readers.
    if biology.sha256(checkpoint_path) != job["checkpoint_sha256"]:
        raise RuntimeError("Worker checkpoint changed after main preflight")
    if biology.sha256(repo_path(str(job["run_path"]))) != job["run_sha256"]:
        raise RuntimeError("Worker run record changed after main preflight")
    if biology.sha256(valid_path) != spec["input_sha256"][f"fold_{fold}_valid"]:
        raise RuntimeError("Worker validation manifest changed after main preflight")
    for key in ("model_weights", "means", "track_manifest", "group_manifest", "fasta", "fai", "gtf", "case_bigwig", "control_bigwig"):
        if biology.sha256(repo_path(spec["locked_inputs"][key])) != spec["input_sha256"][key]:
            raise RuntimeError(f"Worker input changed after main preflight: {key}")

    from scripts import train_v2_model
    from scripts import v2_training_components as components
    from scripts.v2_bigwig_dataset import V2BigWigDataset

    means_rows = biology.read_tsv(repo_path(spec["locked_inputs"]["means"]))
    mean_column = spec["model_contract"]["mean_column_template"].format(fold=fold)
    fold_means = torch.tensor(
        [float(row[mean_column]) for row in means_rows], dtype=torch.float32
    )
    device = torch.device("cuda:0")
    model = train_v2_model.build_model(
        EXPECTED_MODEL,
        fold_means,
        device,
        int(spec["model_contract"]["hidden_channels"]),
    )
    checkpoint = _load_checkpoint(
        model,
        checkpoint_path,
        fold=fold,
        seed=seed,
        spec=spec,
    )
    track_ids = [row["group_id"] for row in means_rows]
    track_indices = [track_ids.index(CASE_GROUP_ID), track_ids.index(CONTROL_GROUP_ID)]
    interval_rows, subwindows = biology.load_fold_validation_subwindows(
        valid_path,
        expected_fold=fold,
        sequence_length=int(spec["model_contract"]["sequence_length"]),
    )
    genes, genes_by_chromosome = biology.load_annotated_genes(
        repo_path(spec["locked_inputs"]["gtf"]), biology.DEVELOPMENT_CHROMOSOMES
    )
    accumulator = biology.Dpy27GeneExonAccumulator(genes, genes_by_chromosome)
    dataset = V2BigWigDataset(
        valid_path,
        track_manifest_path=repo_path(spec["locked_inputs"]["track_manifest"]),
        group_manifest_path=repo_path(spec["locked_inputs"]["group_manifest"]),
        fasta_path=repo_path(spec["locked_inputs"]["fasta"]),
        fai_path=repo_path(spec["locked_inputs"]["fai"]),
        gtf_path=repo_path(spec["locked_inputs"]["gtf"]),
        track_indices=track_indices,
        max_io_workers=2,
    )
    started = time.monotonic()
    model.eval()
    try:
        with torch.no_grad():
            for ordinal, window in enumerate(subwindows, 1):
                context_start = int(interval_rows[window.interval_index]["start"])
                item = dataset.get_subwindow_item(
                    window.interval_index,
                    shift_bp=0,
                    crop_offset_bp=window.start - context_start,
                    crop_length_bp=int(spec["model_contract"]["sequence_length"]),
                )
                if int(item["interval_start"]) != window.start:
                    raise RuntimeError("Dataset subwindow coordinate mismatch")
                dna = item["dna_sequence"].unsqueeze(0).to(device)
                with torch.autocast(device_type="cuda", dtype=torch.bfloat16):
                    model_space = model(dna)[1]
                prediction = components.unscale_predictions_experimental_space(
                    model_space.float(), fold_means.to(device), 1
                )[0, track_indices].clamp_min(0).cpu().numpy()
                observation = item["target_1bp"].clamp_min(0).numpy()
                accumulator.update(
                    window.chromosome,
                    window.start,
                    prediction,
                    observation,
                    block_id=window.block_id,
                )
                print(
                    f"p10_subwindow\tseed={seed}\tfold={fold}\t{ordinal}/{len(subwindows)}\t"
                    f"{window.chromosome}:{window.start}-{window.end}",
                    flush=True,
                )
    finally:
        dataset.close()
    gene_rows = accumulator.rows(fold=fold, seed=seed)
    endpoint = {
        "analysis_label": ANALYSIS_LABEL,
        "fold": fold,
        "seed": seed,
        "checkpoint_path": str(job["checkpoint_path"]),
        "checkpoint_sha256": str(job["checkpoint_sha256"]),
        **biology.dpy27_run_endpoints(gene_rows),
    }
    paths = _shard_paths(spec, seed=seed, fold=fold)
    biology.write_tsv(paths["genes"], gene_rows, GENE_FIELDS)
    biology.write_json(paths["endpoint"], endpoint)
    shard_audit = {
        "schema_version": 1,
        "analysis_contract": ANALYSIS_CONTRACT,
        "analysis_label": ANALYSIS_LABEL,
        "git_commit": git_commit(),
        "implementation_sha256": implementation_sha256(),
        "status": "completed",
        "completed_at": utc_now(),
        "elapsed_seconds": time.monotonic() - started,
        "spec_sha256": spec_sha256,
        "fold": fold,
        "seed": seed,
        "physical_gpu": physical_gpu,
        "checkpoint_sha256": job["checkpoint_sha256"],
        "checkpoint_step": checkpoint["step"],
        "validation_sha256": spec["input_sha256"][f"fold_{fold}_valid"],
        "subwindows_evaluated": len(subwindows),
        "genes_with_any_exon_coverage": len(gene_rows),
        "genes_with_complete_exon_coverage": sum(
            bool(row["complete_exon_coverage"]) for row in gene_rows
        ),
        "gene_contrasts_path": str(paths["genes"].relative_to(REPO_ROOT)),
        "gene_contrasts_sha256": biology.sha256(paths["genes"]),
        "run_endpoint_path": str(paths["endpoint"].relative_to(REPO_ROOT)),
        "run_endpoint_sha256": biology.sha256(paths["endpoint"]),
        "normalization_boundary": spec["normalization_boundary"],
        "final_test_access": "prohibited",
        "locked_test_block_signal_reads": 0,
    }
    biology.write_json(paths["audit"], shard_audit)
    return {"status": "completed", "gene_rows": len(gene_rows), **endpoint}


def _run_worker(args: argparse.Namespace) -> None:
    spec = read_json(SPEC_PATH)
    _validate_static_contract(spec)
    matrix = {
        (int(row["seed"]), int(row["fold"])): row for row in spec["checkpoint_matrix"]
    }
    job = matrix.get((args.worker_seed, args.worker_fold))
    if job is None:
        raise RuntimeError("Worker job is outside the frozen checkpoint matrix")
    result = run_one_checkpoint(
        spec,
        job,
        spec_sha256=args.spec_sha256,
        physical_gpu=args.physical_gpu,
    )
    print(json.dumps(result, sort_keys=True))


def _launch_worker(
    spec_sha256: str,
    job: Mapping[str, Any],
    physical_gpu: int,
    log_path: Path,
) -> dict[str, Any]:
    command = [
        sys.executable,
        "-m",
        "scripts.run_dpy27_internal_application",
        "--worker-seed",
        str(job["seed"]),
        "--worker-fold",
        str(job["fold"]),
        "--physical-gpu",
        str(physical_gpu),
        "--spec-sha256",
        spec_sha256,
    ]
    environment = os.environ.copy()
    environment["CUDA_VISIBLE_DEVICES"] = str(physical_gpu)
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("a") as handle:
        handle.write("command\t" + " ".join(command) + "\n")
        handle.flush()
        subprocess.run(
            command,
            cwd=REPO_ROOT,
            env=environment,
            stdout=handle,
            stderr=subprocess.STDOUT,
            check=True,
        )
    return {
        "seed": int(job["seed"]),
        "fold": int(job["fold"]),
        "physical_gpu": physical_gpu,
        "log_path": str(log_path.relative_to(REPO_ROOT)),
    }


def _run_gpu_queue(
    spec_sha256: str,
    jobs: list[Mapping[str, Any]],
    physical_gpu: int,
    log_root: Path,
) -> list[dict[str, Any]]:
    """Run one serial queue so a physical GPU is never oversubscribed."""

    results: list[dict[str, Any]] = []
    for job in jobs:
        log_path = log_root / f"seed_{job['seed']}_fold_{job['fold']}.log"
        results.append(_launch_worker(spec_sha256, job, physical_gpu, log_path))
    return results


def _load_all_shards(
    spec: Mapping[str, Any], spec_sha256: str
) -> tuple[list[dict[str, str]], list[dict[str, Any]]]:
    gene_rows: list[dict[str, str]] = []
    endpoints: list[dict[str, Any]] = []
    for job in sorted(
        spec["checkpoint_matrix"], key=lambda row: (int(row["seed"]), int(row["fold"]))
    ):
        completed = _completed_shard(spec, job, spec_sha256)
        if completed is None:
            raise RuntimeError(f"Missing completed P10 shard: seed={job['seed']}, fold={job['fold']}")
        genes, endpoint = completed
        gene_rows.extend(genes)
        endpoints.append(endpoint)
    return gene_rows, endpoints


def _figure_source_rows(
    fold_gene_rows: list[dict[str, object]],
    fold_rows: list[dict[str, object]],
    summaries: list[dict[str, object]],
) -> list[dict[str, object]]:
    source: list[dict[str, object]] = []
    complete = [row for row in fold_gene_rows if bool(row["complete_exon_coverage"])]
    for row in sorted(complete, key=lambda value: (int(value["fold"]), str(value["gene_id"]))):
        source.append(
            {
                "panel": "a_gene_contrasts",
                "fold": row["fold"],
                "seed": "seed_mean",
                "gene_id": row["gene_id"],
                "gene_name": row["gene_name"],
                "chromosome": row["chromosome"],
                "endpoint": "gene_log1p_contrast",
                "observed": row["observed_log1p_contrast_g0005_minus_g0006"],
                "predicted": row["predicted_log1p_contrast_g0005_minus_g0006"],
                "estimate": "",
                "ci_95_low": "",
                "ci_95_high": "",
            }
        )
    for row in fold_rows:
        for kind in ("observed", "predicted"):
            endpoint = f"{kind}_median_contrast_x_minus_autosomes"
            source.append(
                {
                    "panel": "b_fold_x_minus_autosomes",
                    "fold": row["fold"],
                    "seed": "seed_mean",
                    "gene_id": "",
                    "gene_name": "",
                    "chromosome": "X_minus_autosomes",
                    "endpoint": endpoint,
                    "observed": row[endpoint] if kind == "observed" else "",
                    "predicted": row[endpoint] if kind == "predicted" else "",
                    "estimate": row[endpoint],
                    "ci_95_low": "",
                    "ci_95_high": "",
                }
            )
        for endpoint in (
            "x_predicted_observed_direction_concordance",
            "predicted_observed_gene_spearman",
        ):
            source.append(
                {
                    "panel": "c_fold_agreement",
                    "fold": row["fold"],
                    "seed": "seed_mean",
                    "gene_id": "",
                    "gene_name": "",
                    "chromosome": "",
                    "endpoint": endpoint,
                    "observed": "",
                    "predicted": "",
                    "estimate": row[endpoint],
                    "ci_95_low": "",
                    "ci_95_high": "",
                }
            )
    for row in summaries:
        source.append(
            {
                "panel": "c_fold_primary_summary",
                "fold": "all_folds",
                "seed": "within_fold_seed_mean",
                "gene_id": "",
                "gene_name": "",
                "chromosome": "",
                "endpoint": row["endpoint"],
                "observed": "",
                "predicted": "",
                "estimate": row["estimate"],
                "ci_95_low": row["ci_95_low"],
                "ci_95_high": row["ci_95_high"],
            }
        )
    return source


def _finalize_outputs(
    spec: Mapping[str, Any],
    spec_sha256: str,
    preflight: Mapping[str, Any],
    execution: dict[str, Any],
) -> None:
    gene_rows, run_rows = _load_all_shards(spec, spec_sha256)
    fold_gene_rows, fold_rows, summaries = biology.fold_primary_summary(
        gene_rows,
        endpoints=tuple(spec["statistics"]["endpoints"]),
        bootstrap_replicates=int(spec["statistics"]["bootstrap_replicates"]),
        bootstrap_seed=int(spec["statistics"]["bootstrap_seed"]),
    )
    outputs = spec["output_paths"]
    gene_path = repo_path(outputs["gene_contrasts"])
    run_path = repo_path(outputs["run_endpoints"])
    fold_path = repo_path(outputs["fold_endpoints"])
    fold_gene_path = repo_path(outputs["fold_gene_contrasts"])
    summary_path = repo_path(outputs["summary_endpoints"])
    source_path = repo_path(outputs["figure_source_data"])
    biology.write_tsv(gene_path, gene_rows, GENE_FIELDS)
    biology.write_tsv(run_path, run_rows)
    biology.write_tsv(fold_path, fold_rows)
    biology.write_tsv(fold_gene_path, fold_gene_rows, GENE_FIELDS)
    biology.write_tsv(summary_path, summaries)
    source_rows = _figure_source_rows(fold_gene_rows, fold_rows, summaries)
    biology.write_tsv(
        source_path,
        source_rows,
        (
            "panel",
            "fold",
            "seed",
            "gene_id",
            "gene_name",
            "chromosome",
            "endpoint",
            "observed",
            "predicted",
            "estimate",
            "ci_95_low",
            "ci_95_high",
        ),
    )
    artifact_paths = {
        "gene_contrasts": gene_path,
        "run_endpoints": run_path,
        "fold_endpoints": fold_path,
        "fold_gene_contrasts": fold_gene_path,
        "summary_endpoints": summary_path,
        "figure_source_data": source_path,
        "preflight_audit": repo_path(outputs["preflight_audit"]),
    }
    execution.update(
        {
            "status": "analysis_completed_figure_pending",
            "completed_at": utc_now(),
            "artifacts": {
                key: {
                    "path": str(path.relative_to(REPO_ROOT)),
                    "sha256": biology.sha256(path),
                    "size_bytes": path.stat().st_size,
                }
                for key, path in artifact_paths.items()
            },
            "run_count": len(run_rows),
            "fold_count": len(fold_rows),
            "gene_rows": len(gene_rows),
            "complete_gene_rows": sum(row["complete_exon_coverage"] == "True" for row in gene_rows),
            "summary": summaries,
            "normalization_boundary": spec["normalization_boundary"],
            "biological_claim_boundary": spec["biological_claim_boundary"],
            "final_test_access": "prohibited",
            "locked_test_block_signal_reads": 0,
            "preflight_model_or_bigwig_reads": preflight["model_or_bigwig_reads"],
        }
    )
    biology.write_json(repo_path(outputs["analysis_audit"]), execution)


def _render_figure_and_complete(
    spec: Mapping[str, Any], spec_sha256: str, execution: dict[str, Any]
) -> None:
    outputs = spec["output_paths"]
    runtime = spec["plot_runtime"]
    python_path = Path(runtime["python_executable"])
    if not python_path.is_file() or biology.sha256(python_path) != runtime["python_executable_sha256"]:
        raise RuntimeError("Frozen P10 plot interpreter is missing or changed")
    for key in ("validate_figure", "audit_pdf_text"):
        path = Path(runtime[f"{key}_path"])
        if not path.is_file() or biology.sha256(path) != runtime[f"{key}_sha256"]:
            raise RuntimeError(f"Frozen P10 plot QA tool is missing or changed: {key}")
    version_command = [
        str(python_path),
        "-c",
        (
            "import json,sys,numpy,matplotlib,PIL; "
            "print(json.dumps({'python':sys.version.split()[0],"
            "'numpy':numpy.__version__,'matplotlib':matplotlib.__version__,"
            "'pillow':PIL.__version__},sort_keys=True))"
        ),
    ]
    versions = json.loads(
        subprocess.run(
            version_command,
            cwd=REPO_ROOT,
            check=True,
            capture_output=True,
            text=True,
        ).stdout
    )
    expected_versions = {
        "python": runtime["python_version"],
        "numpy": runtime["numpy_version"],
        "matplotlib": runtime["matplotlib_version"],
        "pillow": runtime["pillow_version"],
    }
    if versions != expected_versions:
        raise RuntimeError(f"P10 plot runtime version mismatch: {versions} != {expected_versions}")
    command = [
        str(python_path),
        "-m",
        "scripts.plot_dpy27_internal_application",
        "--spec",
        str(SPEC_PATH),
        "--spec-sha256",
        spec_sha256,
    ]
    completed = subprocess.run(
        command,
        cwd=REPO_ROOT,
        check=False,
        capture_output=True,
        text=True,
    )
    if completed.returncode != 0:
        raise RuntimeError(
            f"P10 plot finalizer failed ({completed.returncode}):\n"
            f"{completed.stdout}\n{completed.stderr}"
        )
    qa_path = repo_path(outputs["figure_qa"])
    qa = read_json(qa_path)
    if qa.get("status") != "passed" or qa.get("spec_sha256") != spec_sha256:
        raise RuntimeError("P10 figure QA contract did not pass")
    figure_artifacts: dict[str, dict[str, object]] = {}
    for key in ("figure_pdf", "figure_svg", "figure_tiff", "figure_png", "figure_qa"):
        path = repo_path(outputs[key])
        if not path.is_file() or path.stat().st_size <= 0:
            raise RuntimeError(f"P10 figure artifact is missing: {key}")
        figure_artifacts[key] = {
            "path": str(path.relative_to(REPO_ROOT)),
            "sha256": biology.sha256(path),
            "size_bytes": path.stat().st_size,
        }
    execution["artifacts"].update(figure_artifacts)
    execution.update(
        {
            "status": "completed",
            "completed_at": utc_now(),
            "plot_command": command,
            "plot_runtime": versions,
            "figure_qa": qa,
            "final_test_access": "prohibited",
            "locked_test_block_signal_reads": 0,
        }
    )
    biology.write_json(repo_path(outputs["analysis_audit"]), execution)


def _run_controller_entry() -> None:
    started = time.monotonic()
    spec = read_json(SPEC_PATH)
    preflight = preflight_contract(SPEC_PATH, enforce_controller=True)
    spec_sha256 = preflight["spec_sha256"]
    outputs = spec["output_paths"]
    preflight_path = repo_path(outputs["preflight_audit"])
    biology.write_json(preflight_path, preflight)
    execution_path = repo_path(outputs["analysis_audit"])
    execution: dict[str, Any] = {
        "schema_version": 1,
        "phase": "P10",
        "analysis_contract": ANALYSIS_CONTRACT,
        "analysis_label": ANALYSIS_LABEL,
        "status": "running",
        "started_at": utc_now(),
        **execution_provenance(),
        "implementation_sha256": implementation_sha256(),
        "completed_at": None,
        "spec_path": str(SPEC_PATH.relative_to(REPO_ROOT)),
        "spec_sha256": spec_sha256,
        "preflight_path": str(preflight_path.relative_to(REPO_ROOT)),
        "preflight_sha256": biology.sha256(preflight_path),
        "physical_gpus": [],
        "worker_results": [],
        "failures": [],
        "final_test_access": "prohibited",
        "locked_test_block_signal_reads": 0,
    }
    biology.write_json(execution_path, execution)
    try:
        snapshot = v2_gpu_resources.snapshot()
        minimum_free = int(spec["gpu_policy"]["minimum_free_mib"])
        try:
            gpus = v2_gpu_resources.select_available(snapshot, count=2, minimum_free_mib=minimum_free)
        except RuntimeError:
            gpus = v2_gpu_resources.select_available(snapshot, count=1, minimum_free_mib=minimum_free)
        execution["physical_gpus"] = gpus
        execution["gpu_snapshot_before"] = snapshot
        biology.write_json(execution_path, execution)
        jobs = list(spec["checkpoint_matrix"])
        log_root = repo_path(outputs["run_root"]) / "logs"
        queues = {
            gpu: [job for ordinal, job in enumerate(jobs) if gpus[ordinal % len(gpus)] == gpu]
            for gpu in gpus
        }
        with ThreadPoolExecutor(max_workers=len(gpus)) as pool:
            futures = {
                pool.submit(_run_gpu_queue, spec_sha256, queue, gpu, log_root): gpu
                for gpu, queue in queues.items()
            }
            for future in as_completed(futures):
                execution["worker_results"].extend(future.result())
                biology.write_json(execution_path, execution)
        execution["elapsed_seconds"] = time.monotonic() - started
        execution["gpu_snapshot_after"] = v2_gpu_resources.snapshot()
        _finalize_outputs(spec, spec_sha256, preflight, execution)
        _render_figure_and_complete(spec, spec_sha256, execution)
    except Exception as error:
        execution.update(
            {
                "status": "failed",
                "completed_at": utc_now(),
                "elapsed_seconds": time.monotonic() - started,
                "failures": [str(error)],
            }
        )
        biology.write_json(execution_path, execution)
        raise


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--worker-seed", type=int, choices=EXPECTED_SEEDS)
    parser.add_argument("--worker-fold", type=int, choices=EXPECTED_FOLDS)
    parser.add_argument("--physical-gpu", type=int, choices=(2, 3))
    parser.add_argument("--spec-sha256")
    return parser


def main(argv: list[str] | None = None) -> None:
    args = build_parser().parse_args(argv)
    worker_values = (
        args.worker_seed,
        args.worker_fold,
        args.physical_gpu,
        args.spec_sha256,
    )
    if any(value is not None for value in worker_values):
        if not all(value is not None for value in worker_values):
            raise SystemExit("All worker arguments must be provided together")
        _run_worker(args)
        return
    _run_controller_entry()


if __name__ == "__main__":
    main()
