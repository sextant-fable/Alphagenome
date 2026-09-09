#!/usr/bin/env python3
"""Extract deterministic held-out coverage examples from one frozen P12 model.

This is a read-only inference utility for figure Source Data. It uses fold 1
development-validation intervals, the P12 B/paper checkpoint and only the
held-out biological-unit labels. It cannot read the locked final-test block.
"""

from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import sys
from typing import Any

import numpy as np
import torch

# Keep repository-local imports valid when this utility is invoked by file path.
REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts import evaluate_v2_model
from scripts import train_v2_model
from scripts import v2_training_components as components
from scripts.v2_bigwig_dataset import V2BigWigDataset


P12_ROOT = REPO_ROOT / "alphagenome_custom/metadata/v2/replicate_holdout_v1"
DEFAULT_CHECKPOINT = REPO_ROOT / "runs/v2_p12_replicate_holdout_retry_20260819/B_paper/seed_20260714/fold_1/checkpoint.pt"
DEFAULT_OUTPUT = REPO_ROOT / "results/v2_p12_replicate_holdout_figure/coverage_examples"
SEED = 20260822


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(8 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def select_tracks(validation: dict[str, Any], group_rows: list[dict[str, str]], assignments: list[dict[str, str]]) -> list[dict[str, Any]]:
    primary_ids = {
        row["group_id"]
        for row in assignments
        if row["assignment_role"] == "primary_candidate"
    }
    per_track = validation["per_track_pearson_128bp"]
    candidates = [(group_id, float(value)) for group_id, value in per_track.items() if group_id in primary_ids and value is not None]
    candidates.sort(key=lambda item: (item[1], item[0]))
    if len(candidates) < 3:
        raise RuntimeError("Fewer than three finite primary held-out tracks")
    quantile_positions = [int(round((len(candidates) - 1) * q)) for q in (0.25, 0.50, 0.75)]
    group_by_id = {row["group_id"]: row for row in group_rows}
    selected = []
    for quantile, position in zip(("lower", "median", "upper"), quantile_positions):
        group_id, pearson = candidates[position]
        selected.append({
            "selection_quantile": quantile,
            "group_id": group_id,
            "heldout_pearson_128bp": pearson,
            "tissue_or_cell_type": group_by_id[group_id].get("tissue_or_cell_type", ""),
            "development_stage": group_by_id[group_id].get("development_stage", ""),
            "sra_study_accession": group_by_id[group_id].get("sra_study_accession", ""),
        })
    return selected


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--checkpoint", default=str(DEFAULT_CHECKPOINT))
    parser.add_argument("--output-dir", default=str(DEFAULT_OUTPUT))
    parser.add_argument("--device", default="cuda")
    args = parser.parse_args()

    checkpoint_path = Path(args.checkpoint).resolve()
    output_dir = Path(args.output_dir).resolve()
    if not checkpoint_path.is_file():
        raise FileNotFoundError(checkpoint_path)
    if args.device not in {"cuda", "cpu"}:
        raise ValueError("--device must be cuda or cpu")
    visible = os.environ.get("CUDA_VISIBLE_DEVICES", "")
    if args.device == "cuda":
        if not torch.cuda.is_available():
            raise RuntimeError("CUDA was requested but the NVIDIA driver is unavailable")
        if visible not in {"2", "3"}:
            raise RuntimeError("Coverage extraction must use exactly one physical GPU 2 or 3")
    else:
        visible = "cpu"
    spec_path = REPO_ROOT / "alphagenome_custom/metadata/v2/p12_replicate_holdout_spec.json"
    if spec_path.read_text().find('"final_test_access": "prohibited"') < 0:
        raise RuntimeError("P12 spec contract is missing final-test prohibition")

    validation_path = REPO_ROOT / "runs/v2_p12_replicate_holdout_retry_20260819/B_paper/seed_20260714/fold_1/validation_heldout.json"
    validation = json.loads(validation_path.read_text())
    if validation.get("locked_test_block_signal_reads") is not False:
        raise RuntimeError("Validation record indicates locked-test access")
    if validation.get("fold") != 1 or validation.get("seed") != 20260714:
        raise RuntimeError("Coverage extraction checkpoint selection changed")

    track_indices_rows = read_tsv(P12_ROOT / "heldout_primary_track_indices.tsv")
    group_rows = read_tsv(P12_ROOT / "heldout_primary_group_manifest.tsv")
    assignments = read_tsv(P12_ROOT / "holdout_assignments.tsv")
    selected = select_tracks(validation, group_rows, assignments)
    heldout_index_by_group = {row["group_id"]: int(row["track_index"]) for row in track_indices_rows}
    selected_positions = [heldout_index_by_group[row["group_id"]] for row in selected]

    device = torch.device(args.device)
    means_path = P12_ROOT / "track_nonzero_means.tsv"
    mean_rows = read_tsv(means_path)
    mean_column = "fold_1_train_nonzero_mean"
    model_fold_means = torch.tensor([float(row[mean_column]) for row in mean_rows], dtype=torch.float32, device=device)
    model = train_v2_model.build_model("B", model_fold_means.detach().cpu(), device, 64, None)
    checkpoint = evaluate_v2_model.load_checkpoint(
        model,
        checkpoint_path,
        expected_model="B",
        expected_loss="paper",
        expected_fold=1,
        expected_sequence_length=131072,
        expected_hidden_channels=64,
        expected_mean_column=mean_column,
    )
    checkpoint_sha = sha256(checkpoint_path)
    if checkpoint_sha != validation.get("checkpoint_sha256"):
        raise RuntimeError("Coverage checkpoint hash differs from the audited P12 validation record")

    prediction_indices = [int(row["track_index"]) for row in track_indices_rows]
    eval_fold_means = model_fold_means[prediction_indices]
    intervals_path = REPO_ROOT / "alphagenome_custom/intervals/v2/fold_1/valid.tsv"
    dataset = V2BigWigDataset(
        intervals_path,
        track_manifest_path=P12_ROOT / "heldout_primary_track_manifest.tsv",
        group_manifest_path=P12_ROOT / "heldout_primary_group_manifest.tsv",
        track_indices=list(range(len(prediction_indices))),
        max_io_workers=16,
    )
    subwindows, _ = evaluate_v2_model.core_subwindows(dataset.intervals, 131072)
    best: dict[str, dict[str, Any]] = {row["group_id"]: {"target_sum": -1.0} for row in selected}
    model.eval()
    with torch.no_grad():
        for ordinal, (index, offset) in enumerate(subwindows, start=1):
            item = dataset.get_subwindow_item(index, shift_bp=0, crop_offset_bp=offset, crop_length_bp=131072)
            dna = item["dna_sequence"].unsqueeze(0).to(device)
            if device.type == "cuda":
                with torch.autocast(device_type="cuda", dtype=torch.bfloat16):
                    predictions = model(dna)
            else:
                predictions = model(dna)
            prediction_128 = components.unscale_predictions_experimental_space(
                predictions[128][:, prediction_indices, ...].float(), eval_fold_means, 128
            ).clamp_min(0)[0].cpu().numpy()
            target_128 = item["target_128bp"].clamp_min(0).cpu().numpy()
            chromosome = str(item["interval_chromosome"])
            interval_start = int(item["interval_start"])
            interval_end = int(item["interval_end"])
            for output_index, selected_row in zip(selected_positions, selected):
                group_id = selected_row["group_id"]
                observed = target_128[output_index]
                predicted = prediction_128[output_index]
                target_sum = float(np.sum(observed))
                if target_sum > best[group_id]["target_sum"]:
                    best[group_id] = {
                        "target_sum": target_sum,
                        "predicted": predicted.astype(float).tolist(),
                        "observed": observed.astype(float).tolist(),
                        "chromosome": chromosome,
                        "interval_start": interval_start,
                        "interval_end": interval_end,
                        "subwindow_ordinal": ordinal,
                    }
            if ordinal == 1 or ordinal % 10 == 0 or ordinal == len(subwindows):
                print(f"coverage extraction progress: {ordinal}/{len(subwindows)} subwindows", flush=True)
    dataset.close()
    if any(value["target_sum"] < 0 for value in best.values()):
        raise RuntimeError("No selected coverage window was extracted")

    output_dir.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, Any]] = []
    selected_by_id = {row["group_id"]: row for row in selected}
    for group_id, result in best.items():
        meta = selected_by_id[group_id]
        for bin_index, (observed, predicted) in enumerate(zip(result["observed"], result["predicted"])):
            rows.append({
                "group_id": group_id,
                "selection_quantile": meta["selection_quantile"],
                "heldout_pearson_128bp": meta["heldout_pearson_128bp"],
                "sra_study_accession": meta["sra_study_accession"],
                "development_stage": meta["development_stage"],
                "tissue_or_cell_type": meta["tissue_or_cell_type"],
                "chromosome": result["chromosome"],
                "interval_start": result["interval_start"],
                "interval_end": result["interval_end"],
                "subwindow_ordinal": result["subwindow_ordinal"],
                "bin_index": bin_index,
                "genomic_start": result["interval_start"] + 128 * bin_index,
                "observed_128bp": observed,
                "predicted_128bp": predicted,
                "selection_rule": "fold1 seed20260714 primary track Pearson quantile; validation window with largest observed 128-bp signal",
            })
    source_path = output_dir / "coverage_examples_source_data.tsv"
    fields = list(rows[0])
    with source_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    metadata = {
        "schema_version": 1,
        "created_at": datetime.now(timezone.utc).replace(microsecond=0).isoformat(),
        "checkpoint": str(checkpoint_path.relative_to(REPO_ROOT)),
        "checkpoint_sha256": checkpoint_sha,
        "fold": 1,
        "seed": 20260714,
        "physical_cuda_visible_devices": visible,
        "validation_intervals": str(intervals_path.relative_to(REPO_ROOT)),
        "validation_intervals_sha256": sha256(intervals_path),
        "means_path": str(means_path.relative_to(REPO_ROOT)),
        "means_sha256": sha256(means_path),
        "heldout_track_manifest": str((P12_ROOT / "heldout_primary_track_manifest.tsv").relative_to(REPO_ROOT)),
        "selected_tracks": selected,
        "subwindows_evaluated": len(subwindows),
        "locked_test_block_signal_reads": 0,
        "source_data": str(source_path.relative_to(REPO_ROOT)),
        "source_data_sha256": sha256(source_path),
        "selection_rule": "three deterministic primary-track quantiles from audited B/paper fold1/seed20260714 heldout Pearson; one window per track selected by largest observed 128-bp signal",
    }
    (output_dir / "coverage_examples_manifest.json").write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n")
    print(json.dumps(metadata, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
