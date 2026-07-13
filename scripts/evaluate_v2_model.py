#!/usr/bin/env python3
"""Evaluate one v2 checkpoint on fixed core-only blocked-CV subwindows."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path

import numpy as np
import torch

from scripts import train_v2_model
from scripts import v2_training_components as components
from scripts.v2_bigwig_dataset import V2BigWigDataset


REPO_ROOT = Path(__file__).resolve().parents[1]
METADATA_DIR = REPO_ROOT / "alphagenome_custom/metadata/v2"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--model", choices=("A", "B", "C"), required=True)
    parser.add_argument("--training-loss", choices=("paper", "log1p_mse"), required=True)
    parser.add_argument("--fold", type=int, choices=range(1, 6), required=True)
    parser.add_argument("--checkpoint", required=True)
    parser.add_argument("--sequence-length", type=int, default=131072)
    parser.add_argument("--hidden-channels", type=int, default=64)
    parser.add_argument("--output", required=True)
    parser.add_argument("--device", default="cuda")
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(8 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def core_subwindows(
    intervals: list[dict[str, str]], sequence_length: int
) -> tuple[list[tuple[int, int]], int]:
    if sequence_length <= 0 or sequence_length % 128:
        raise ValueError("sequence_length must be positive and divisible by 128")
    result = []
    eligible_bases = 0
    for index, row in enumerate(intervals):
        core_start = ((int(row["core_start"]) + 127) // 128) * 128
        core_end = (int(row["core_end"]) // 128) * 128
        eligible_bases += max(0, core_end - core_start)
        for start in range(core_start, core_end - sequence_length + 1, sequence_length):
            result.append((index, start - int(row["start"])))
    if not result:
        raise RuntimeError("No complete core-only validation subwindows")
    return result, eligible_bases


def load_checkpoint(
    model: torch.nn.Module,
    checkpoint_path: Path,
    *,
    expected_model: str,
    expected_loss: str,
    expected_fold: int,
    expected_sequence_length: int,
    expected_hidden_channels: int,
) -> dict[str, object]:
    checkpoint = torch.load(checkpoint_path, map_location="cpu", weights_only=False)
    expected = {
        "model_id": expected_model,
        "loss": expected_loss,
        "fold": expected_fold,
        "sequence_length": expected_sequence_length,
        "hidden_channels": expected_hidden_channels,
    }
    mismatches = {
        key: (checkpoint.get(key), value)
        for key, value in expected.items()
        if checkpoint.get(key) != value
    }
    if mismatches:
        raise RuntimeError(f"Checkpoint configuration mismatch: {mismatches}")
    trainable_names = {
        name for name, parameter in model.named_parameters() if parameter.requires_grad
    }
    checkpoint_state = checkpoint.get("trainable_model_state", {})
    if set(checkpoint_state) != trainable_names:
        raise RuntimeError("Checkpoint trainable parameter set mismatch")
    state = model.state_dict()
    for name, value in checkpoint_state.items():
        if name not in state or state[name].shape != value.shape:
            raise RuntimeError(f"Checkpoint tensor mismatch: {name}")
        state[name] = value
    model.load_state_dict(state, strict=True)
    return checkpoint


def pearson_from_sums(
    count: int,
    sum_x: np.ndarray,
    sum_y: np.ndarray,
    sum_xx: np.ndarray,
    sum_yy: np.ndarray,
    sum_xy: np.ndarray,
) -> np.ndarray:
    numerator = count * sum_xy - sum_x * sum_y
    denominator = np.sqrt(
        np.maximum(count * sum_xx - np.square(sum_x), 0)
        * np.maximum(count * sum_yy - np.square(sum_y), 0)
    )
    result = np.full(sum_x.shape, np.nan, dtype=np.float64)
    valid = denominator > 0
    result[valid] = numerator[valid] / denominator[valid]
    return result


def main() -> None:
    args = parse_args()
    if args.sequence_length % 128:
        raise ValueError("sequence-length must be divisible by 128")
    device = torch.device(args.device)
    visible = os.environ.get("CUDA_VISIBLE_DEVICES", "")
    if device.type != "cuda" or not torch.cuda.is_available():
        raise RuntimeError("Registered P6B evaluation requires CUDA")
    if visible not in {"2", "3"}:
        raise RuntimeError("Evaluation must use exactly one physical GPU 2 or 3")

    means_path = METADATA_DIR / "track_nonzero_means_v2.tsv"
    means_rows = read_tsv(means_path)
    mean_column = f"fold_{args.fold}_train_nonzero_mean"
    fold_means = torch.tensor(
        [float(row[mean_column]) for row in means_rows],
        dtype=torch.float32,
        device=device,
    )
    intervals_path = (
        REPO_ROOT / f"alphagenome_custom/intervals/v2/fold_{args.fold}/valid.tsv"
    )
    dataset = V2BigWigDataset(intervals_path, max_io_workers=16)
    subwindows, eligible_bases = core_subwindows(
        dataset.intervals, args.sequence_length
    )
    model = train_v2_model.build_model(
        args.model,
        fold_means.detach().cpu(),
        device,
        args.hidden_channels,
    )
    checkpoint_path = REPO_ROOT / args.checkpoint
    checkpoint = load_checkpoint(
        model,
        checkpoint_path,
        expected_model=args.model,
        expected_loss=args.training_loss,
        expected_fold=args.fold,
        expected_sequence_length=args.sequence_length,
        expected_hidden_channels=args.hidden_channels,
    )
    model.eval()
    metric_sums: dict[str, float] = {}
    n_tracks = len(means_rows)
    count_128 = 0
    sum_x = np.zeros(n_tracks, dtype=np.float64)
    sum_y = np.zeros(n_tracks, dtype=np.float64)
    sum_xx = np.zeros(n_tracks, dtype=np.float64)
    sum_yy = np.zeros(n_tracks, dtype=np.float64)
    sum_xy = np.zeros(n_tracks, dtype=np.float64)
    with torch.no_grad():
        for ordinal, (index, offset) in enumerate(subwindows, start=1):
            item = dataset.get_subwindow_item(
                index,
                shift_bp=0,
                crop_offset_bp=offset,
                crop_length_bp=args.sequence_length,
            )
            dna = item["dna_sequence"].unsqueeze(0).to(device)
            targets = {
                1: item["target_1bp"].unsqueeze(0).to(device),
                128: item["target_128bp"].unsqueeze(0).to(device),
            }
            track_mask = item["track_mask"].unsqueeze(0).to(device)
            track_strand = item["track_strand"].unsqueeze(0).to(device)
            gene_mask = item["gene_mask"].unsqueeze(0).to(device)
            with torch.autocast(device_type="cuda", dtype=torch.bfloat16):
                predictions = model(dna)
                paper_loss, paper_metrics = components.dual_resolution_paper_loss(
                    predictions,
                    targets,
                    track_means=fold_means,
                    track_mask=track_mask,
                    track_strand=track_strand,
                    gene_mask=gene_mask,
                )
                log_loss, log_metrics = components.log1p_mse_loss(
                    predictions,
                    targets,
                    track_means=fold_means,
                    track_mask=track_mask,
                )
            current = {
                "paper_loss": paper_loss,
                "log1p_mse": log_loss,
                **{f"paper_{key}": value for key, value in paper_metrics.items()},
                **{f"log1p_{key}": value for key, value in log_metrics.items()},
            }
            for key, value in current.items():
                numeric = float(value.detach())
                if not math.isfinite(numeric):
                    raise RuntimeError(f"Non-finite validation metric {key}")
                metric_sums[key] = metric_sums.get(key, 0.0) + numeric

            prediction_128 = components.unscale_predictions_experimental_space(
                predictions[128].float(), fold_means, 128
            )
            x = torch.log1p(prediction_128.clamp_min(0))[0].double().cpu().numpy()
            y = torch.log1p(targets[128].clamp_min(0))[0].double().cpu().numpy()
            count_128 += x.shape[1]
            sum_x += x.sum(axis=1)
            sum_y += y.sum(axis=1)
            sum_xx += np.square(x).sum(axis=1)
            sum_yy += np.square(y).sum(axis=1)
            sum_xy += (x * y).sum(axis=1)
            print(
                f"validation_subwindow\t{ordinal}/{len(subwindows)}\t"
                f"{item['interval_chromosome']}:{int(item['interval_start'])}-{int(item['interval_end'])}",
                flush=True,
            )
    dataset.close()
    means = {key: value / len(subwindows) for key, value in metric_sums.items()}
    per_track_pearson = pearson_from_sums(
        count_128, sum_x, sum_y, sum_xx, sum_yy, sum_xy
    )
    finite_pearson = per_track_pearson[np.isfinite(per_track_pearson)]
    if not finite_pearson.size:
        raise RuntimeError("No finite per-track validation Pearson values")
    record = {
        "schema_version": 1,
        "phase": "P6B",
        "model": args.model,
        "training_loss": args.training_loss,
        "fold": args.fold,
        "seed": checkpoint["seed"],
        "checkpoint_path": args.checkpoint,
        "checkpoint_sha256": sha256(checkpoint_path),
        "sequence_length": args.sequence_length,
        "validation_policy": "all_complete_nonoverlapping_core_only_131072bp_subwindows",
        "validation_subwindows": len(subwindows),
        "validation_bases": len(subwindows) * args.sequence_length,
        "eligible_aligned_core_bases": eligible_bases,
        "validation_core_coverage_fraction": (
            len(subwindows) * args.sequence_length / eligible_bases
        ),
        "mean_metrics": means,
        "mean_per_track_pearson_128bp": float(finite_pearson.mean()),
        "finite_per_track_pearson_128bp": int(finite_pearson.size),
        "min_per_track_pearson_128bp": float(finite_pearson.min()),
        "max_per_track_pearson_128bp": float(finite_pearson.max()),
        "per_track_pearson_128bp": {
            row["group_id"]: float(value) if math.isfinite(value) else None
            for row, value in zip(means_rows, per_track_pearson, strict=True)
        },
        "intervals_sha256": sha256(intervals_path),
        "means_sha256": sha256(means_path),
        "physical_cuda_visible_devices": visible,
        "cuda_device_name": torch.cuda.get_device_name(0),
        "chromosome_x_read": False,
    }
    output_path = REPO_ROOT / args.output
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
    print(json.dumps(record, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
