#!/usr/bin/env python3
"""Diagnose RNA-seq prediction scale for adapter checkpoints."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import torch

from alphagenome_pytorch import AlphaGenome
from alphagenome_rna_seq11_adapter import (
    RnaSeq11Adapter,
    bin_rna_seq_target,
    count_parameters,
    load_adapter_checkpoint,
    parse_resolutions,
    pearson_from_sums,
    pick_device,
)
from alphagenome_rna_seq11_eval import load_track_names
from torch_rna_seq_dataset import make_dataloader


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--dataset-dir",
        default="alphagenome_custom/datasets/rna_seq_npz_valid",
        help="Raw-target NPZ dataset directory to diagnose.",
    )
    parser.add_argument(
        "--weights",
        default="weights/alphagenome_pytorch/model_all_folds.safetensors",
        help="Local alphagenome-pytorch safetensors checkpoint.",
    )
    parser.add_argument("--checkpoint", required=True)
    parser.add_argument("--metrics-output", required=True)
    parser.add_argument("--batch-size", type=int, default=1)
    parser.add_argument("--num-workers", type=int, default=0)
    parser.add_argument("--max-examples", type=int, default=None)
    parser.add_argument("--embedding-resolution", type=int, choices=[1, 128], default=None)
    parser.add_argument("--head-type", choices=["linear", "genome-tracks"], default=None)
    parser.add_argument("--head-resolutions", default=None)
    parser.add_argument("--organism-index", type=int, default=None)
    parser.add_argument(
        "--device",
        default="auto",
        help="auto, cpu, cuda, cuda:0, mps, etc.",
    )
    return parser.parse_args()


def expand_mask(mask: torch.Tensor, tensor: torch.Tensor) -> torch.Tensor:
    while mask.ndim < tensor.ndim:
        mask = mask.unsqueeze(-1)
    return mask.to(device=tensor.device, dtype=torch.bool)


def raw_prediction_128bp(
    model: RnaSeq11Adapter,
    prediction: torch.Tensor | dict[int, torch.Tensor],
    *,
    checkpoint_target_transform: str,
) -> torch.Tensor:
    if isinstance(prediction, dict):
        if 128 in prediction:
            return prediction[128].clamp_min(0.0)
        if 1 in prediction:
            return bin_rna_seq_target(prediction[1], 128)
        raise ValueError("Prediction dict contains neither 1 bp nor 128 bp output")

    raw_prediction = prediction
    if checkpoint_target_transform == "log1p":
        raw_prediction = torch.expm1(raw_prediction).clamp_min(0.0)
    elif checkpoint_target_transform == "none":
        raw_prediction = raw_prediction.clamp_min(0.0)
    else:
        raise ValueError(
            "checkpoint target_transform must be 'log1p' or 'none', got "
            f"{checkpoint_target_transform}"
        )
    return bin_rna_seq_target(raw_prediction, 128)


def safe_ratio(numerator: torch.Tensor, denominator: torch.Tensor) -> torch.Tensor:
    return numerator / denominator.clamp_min(1e-12)


def main() -> None:
    args = parse_args()
    dataset_dir = Path(args.dataset_dir)
    weights = Path(args.weights)
    checkpoint_path = Path(args.checkpoint)
    metrics_output = Path(args.metrics_output)
    device = pick_device(args.device)

    checkpoint = load_adapter_checkpoint(checkpoint_path)
    embedding_resolution = (
        args.embedding_resolution
        if args.embedding_resolution is not None
        else int(checkpoint["embedding_resolution"])
    )
    checkpoint_head_resolutions = tuple(
        int(resolution)
        for resolution in checkpoint.get("head_resolutions", [embedding_resolution])
    )
    head_type = (
        args.head_type
        if args.head_type is not None
        else str(checkpoint.get("head_type", "linear"))
    )
    head_resolutions = (
        parse_resolutions(args.head_resolutions, fallback_resolution=embedding_resolution)
        if args.head_resolutions is not None
        else checkpoint_head_resolutions
    )
    organism_index = (
        args.organism_index
        if args.organism_index is not None
        else int(checkpoint["organism_index"])
    )
    checkpoint_target_transform = str(checkpoint.get("target_transform", "log1p"))

    dataloader = make_dataloader(
        dataset_dir,
        batch_size=args.batch_size,
        shuffle=False,
        num_workers=args.num_workers,
        target_transform="none",
        max_examples=args.max_examples,
    )
    n_tracks = int(dataloader.dataset.metadata["n_tracks"])
    if int(checkpoint["n_tracks"]) != n_tracks:
        raise ValueError(
            f"Checkpoint n_tracks={checkpoint['n_tracks']}, dataset n_tracks={n_tracks}"
        )

    print(f"dataset_dir\t{dataset_dir}")
    print(f"weights\t{weights}")
    print(f"checkpoint\t{checkpoint_path}")
    print(f"n_examples\t{len(dataloader.dataset)}")
    print(f"n_tracks\t{n_tracks}")
    print(f"checkpoint_target_transform\t{checkpoint_target_transform}")
    print(f"head_type\t{head_type}")
    print(f"head_resolutions\t{','.join(map(str, head_resolutions))}")
    print(f"embedding_resolution\t{embedding_resolution}")
    print(f"organism_index\t{organism_index}")
    print(f"device\t{device}")
    if device.type == "cuda":
        print(f"cuda_device\t{torch.cuda.get_device_name(device)}")

    base_model = AlphaGenome.from_pretrained(weights, device=device)
    checkpoint_head_state = checkpoint["adapter_head_state_dict"]
    track_means = checkpoint_head_state.get("track_means")
    model = RnaSeq11Adapter(
        base_model,
        n_tracks=n_tracks,
        embedding_resolution=embedding_resolution,
        organism_index=organism_index,
        head_type=head_type,
        head_resolutions=head_resolutions,
        track_means=track_means,
    ).to(device)
    model.head.load_state_dict(checkpoint_head_state)
    model.eval()

    print(f"base_parameters\t{count_parameters(model.base_model)}")
    print(f"head_parameters\t{count_parameters(model.head)}")

    pred_mean_sum = torch.zeros(n_tracks, dtype=torch.float64)
    target_mean_sum = torch.zeros(n_tracks, dtype=torch.float64)
    pred_bin_sum_sum = torch.zeros(n_tracks, dtype=torch.float64)
    target_bin_sum_sum = torch.zeros(n_tracks, dtype=torch.float64)
    target_nonzero_mean_sum = torch.zeros(n_tracks, dtype=torch.float64)
    target_nonzero_bins = torch.zeros(n_tracks, dtype=torch.float64)
    sse = torch.zeros(n_tracks, dtype=torch.float64)
    sae = torch.zeros(n_tracks, dtype=torch.float64)
    sum_x = torch.zeros(n_tracks, dtype=torch.float64)
    sum_y = torch.zeros(n_tracks, dtype=torch.float64)
    sum_x2 = torch.zeros(n_tracks, dtype=torch.float64)
    sum_y2 = torch.zeros(n_tracks, dtype=torch.float64)
    sum_xy = torch.zeros(n_tracks, dtype=torch.float64)
    denominator = torch.zeros(n_tracks, dtype=torch.float64)
    batches = 0
    examples = 0

    with torch.no_grad():
        for batch in dataloader:
            dna_sequence = batch["dna_sequence"].to(device=device, dtype=torch.float32)
            rna_seq = batch["rna_seq"].to(device=device, dtype=torch.float32)
            rna_seq_mask = batch["rna_seq_mask"].to(device=device)

            prediction = model(
                dna_sequence,
                target_length=rna_seq.shape[-1],
                return_scaled=False,
            )
            prediction_128 = raw_prediction_128bp(
                model,
                prediction,
                checkpoint_target_transform=checkpoint_target_transform,
            )
            target_128 = bin_rna_seq_target(rna_seq, 128)
            mask = expand_mask(rna_seq_mask, prediction_128)
            mask64 = mask.to(dtype=torch.float64)

            prediction_mean = prediction_128 / 128.0
            target_mean = target_128 / 128.0
            prediction_metric = torch.log1p(prediction_mean)
            target_metric = torch.log1p(target_mean)
            error = prediction_metric.to(dtype=torch.float64) - target_metric.to(
                dtype=torch.float64
            )

            pred_mean64 = prediction_mean.to(dtype=torch.float64)
            target_mean64 = target_mean.to(dtype=torch.float64)
            pred_bin64 = prediction_128.to(dtype=torch.float64)
            target_bin64 = target_128.to(dtype=torch.float64)
            prediction_metric64 = prediction_metric.to(dtype=torch.float64)
            target_metric64 = target_metric.to(dtype=torch.float64)

            pred_mean_sum += (pred_mean64 * mask64).sum(dim=(0, 2)).detach().cpu()
            target_mean_sum += (target_mean64 * mask64).sum(dim=(0, 2)).detach().cpu()
            pred_bin_sum_sum += (pred_bin64 * mask64).sum(dim=(0, 2)).detach().cpu()
            target_bin_sum_sum += (target_bin64 * mask64).sum(dim=(0, 2)).detach().cpu()
            nonzero = (target_bin64 > 0) & mask
            target_nonzero_mean_sum += (
                target_mean64 * nonzero.to(dtype=torch.float64)
            ).sum(dim=(0, 2)).detach().cpu()
            target_nonzero_bins += nonzero.sum(dim=(0, 2)).to(dtype=torch.float64).cpu()
            sse += (error.square() * mask64).sum(dim=(0, 2)).detach().cpu()
            sae += (error.abs() * mask64).sum(dim=(0, 2)).detach().cpu()
            sum_x += (prediction_metric64 * mask64).sum(dim=(0, 2)).detach().cpu()
            sum_y += (target_metric64 * mask64).sum(dim=(0, 2)).detach().cpu()
            sum_x2 += (
                prediction_metric64.square() * mask64
            ).sum(dim=(0, 2)).detach().cpu()
            sum_y2 += (target_metric64.square() * mask64).sum(dim=(0, 2)).detach().cpu()
            sum_xy += (
                prediction_metric64 * target_metric64 * mask64
            ).sum(dim=(0, 2)).detach().cpu()
            denominator += mask64.expand_as(prediction_metric64).sum(dim=(0, 2)).cpu()

            batches += 1
            examples += int(dna_sequence.shape[0])

    mse = sse / denominator.clamp_min(1.0)
    mae = sae / denominator.clamp_min(1.0)
    pred_mean = pred_mean_sum / denominator.clamp_min(1.0)
    target_mean = target_mean_sum / denominator.clamp_min(1.0)
    pred_bin_mean = pred_bin_sum_sum / denominator.clamp_min(1.0)
    target_bin_mean = target_bin_sum_sum / denominator.clamp_min(1.0)
    target_nonzero_fraction = target_nonzero_bins / denominator.clamp_min(1.0)
    target_nonzero_mean = target_nonzero_mean_sum / target_nonzero_bins.clamp_min(1.0)
    pred_over_target = safe_ratio(pred_mean, target_mean)
    per_track_pearson = [
        pearson_from_sums(
            float(sum_x[index]),
            float(sum_y[index]),
            float(sum_x2[index]),
            float(sum_y2[index]),
            float(sum_xy[index]),
            float(denominator[index]),
        )
        for index in range(n_tracks)
    ]

    overall_denominator = denominator.sum().clamp_min(1.0)
    overall_mse = float(sse.sum() / overall_denominator)
    overall_mae = float(sae.sum() / overall_denominator)
    overall_pred_mean = float(pred_mean_sum.sum() / overall_denominator)
    overall_target_mean = float(target_mean_sum.sum() / overall_denominator)
    overall_pred_bin_mean = float(pred_bin_sum_sum.sum() / overall_denominator)
    overall_target_bin_mean = float(target_bin_sum_sum.sum() / overall_denominator)
    overall_nonzero_fraction = float(target_nonzero_bins.sum() / overall_denominator)
    overall_nonzero_mean = float(
        target_nonzero_mean_sum.sum() / target_nonzero_bins.sum().clamp_min(1.0)
    )
    overall_pearson = pearson_from_sums(
        float(sum_x.sum()),
        float(sum_y.sum()),
        float(sum_x2.sum()),
        float(sum_y2.sum()),
        float(sum_xy.sum()),
        float(denominator.sum()),
    )

    track_labels = load_track_names(dataloader.dataset.metadata, n_tracks)
    metrics_output.parent.mkdir(parents=True, exist_ok=True)
    with metrics_output.open("w", newline="") as handle:
        fieldnames = [
            "metric_scope",
            "track_index",
            "track_id",
            "track_name",
            "bins",
            "prediction_mean_signal",
            "target_mean_signal",
            "prediction_over_target_mean_signal",
            "prediction_128bp_bin_sum_mean",
            "target_128bp_bin_sum_mean",
            "target_nonzero_bin_fraction",
            "target_nonzero_mean_signal",
            "common128_log1pmean_mse",
            "common128_log1pmean_mae",
            "common128_log1pmean_pearson",
            "batches",
            "examples",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerow(
            {
                "metric_scope": "overall",
                "track_index": "",
                "track_id": "",
                "track_name": "",
                "bins": int(denominator.sum().item()),
                "prediction_mean_signal": f"{overall_pred_mean:.8g}",
                "target_mean_signal": f"{overall_target_mean:.8g}",
                "prediction_over_target_mean_signal": (
                    f"{overall_pred_mean / max(overall_target_mean, 1e-12):.8g}"
                ),
                "prediction_128bp_bin_sum_mean": f"{overall_pred_bin_mean:.8g}",
                "target_128bp_bin_sum_mean": f"{overall_target_bin_mean:.8g}",
                "target_nonzero_bin_fraction": f"{overall_nonzero_fraction:.8g}",
                "target_nonzero_mean_signal": f"{overall_nonzero_mean:.8g}",
                "common128_log1pmean_mse": f"{overall_mse:.8g}",
                "common128_log1pmean_mae": f"{overall_mae:.8g}",
                "common128_log1pmean_pearson": f"{overall_pearson:.8g}",
                "batches": batches,
                "examples": examples,
            }
        )
        for index, track_label in enumerate(track_labels):
            track_id, track_name = track_label.split("\t", maxsplit=1)
            writer.writerow(
                {
                    "metric_scope": "track",
                    "track_index": index,
                    "track_id": track_id,
                    "track_name": track_name,
                    "bins": int(denominator[index].item()),
                    "prediction_mean_signal": f"{float(pred_mean[index]):.8g}",
                    "target_mean_signal": f"{float(target_mean[index]):.8g}",
                    "prediction_over_target_mean_signal": (
                        f"{float(pred_over_target[index]):.8g}"
                    ),
                    "prediction_128bp_bin_sum_mean": (
                        f"{float(pred_bin_mean[index]):.8g}"
                    ),
                    "target_128bp_bin_sum_mean": f"{float(target_bin_mean[index]):.8g}",
                    "target_nonzero_bin_fraction": (
                        f"{float(target_nonzero_fraction[index]):.8g}"
                    ),
                    "target_nonzero_mean_signal": (
                        f"{float(target_nonzero_mean[index]):.8g}"
                    ),
                    "common128_log1pmean_mse": f"{float(mse[index]):.8g}",
                    "common128_log1pmean_mae": f"{float(mae[index]):.8g}",
                    "common128_log1pmean_pearson": (
                        f"{float(per_track_pearson[index]):.8g}"
                    ),
                    "batches": batches,
                    "examples": examples,
                }
            )

    print(f"batches\t{batches}")
    print(f"examples\t{examples}")
    print(f"prediction_mean_signal\t{overall_pred_mean:.8g}")
    print(f"target_mean_signal\t{overall_target_mean:.8g}")
    print(
        "prediction_over_target_mean_signal\t"
        f"{overall_pred_mean / max(overall_target_mean, 1e-12):.8g}"
    )
    print(f"target_nonzero_bin_fraction\t{overall_nonzero_fraction:.8g}")
    print(f"target_nonzero_mean_signal\t{overall_nonzero_mean:.8g}")
    print(f"common128_log1pmean_mse\t{overall_mse:.8g}")
    print(f"common128_log1pmean_mae\t{overall_mae:.8g}")
    print(f"common128_log1pmean_pearson\t{overall_pearson:.8g}")
    print(f"metrics_output\t{metrics_output}")
    if device.type == "cuda":
        print(f"cuda_max_memory_allocated_mb\t{torch.cuda.max_memory_allocated() / 1024**2:.1f}")


if __name__ == "__main__":
    main()
