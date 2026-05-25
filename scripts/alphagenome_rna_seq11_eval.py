#!/usr/bin/env python3
"""Reload and evaluate an 11-track C. elegans RNA-seq adapter checkpoint."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import torch

from alphagenome_pytorch import AlphaGenome
from alphagenome_rna_seq11_adapter import (
    RnaSeq11Adapter,
    count_parameters,
    evaluate_128bp_binned_pointwise_metrics,
    evaluate_masked_mse,
    evaluate_masked_mse_by_track,
    evaluate_pointwise_metrics,
    load_adapter_checkpoint,
    parse_resolutions,
    pearson_from_sums,
    pick_device,
)
from torch_rna_seq_dataset import make_dataloader


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--dataset-dir",
        default="alphagenome_custom/datasets/rna_seq_npz_valid",
        help="NPZ dataset directory to evaluate.",
    )
    parser.add_argument(
        "--weights",
        default="weights/alphagenome_pytorch/model_all_folds.safetensors",
        help="Local alphagenome-pytorch safetensors checkpoint.",
    )
    parser.add_argument(
        "--checkpoint",
        required=True,
        help="Adapter-only checkpoint produced by alphagenome_rna_seq11_finetune.py.",
    )
    parser.add_argument("--batch-size", type=int, default=1)
    parser.add_argument("--num-workers", type=int, default=0)
    parser.add_argument("--max-examples", type=int, default=None)
    parser.add_argument(
        "--per-track",
        action="store_true",
        help="Also compute and print per-track masked MSE.",
    )
    parser.add_argument(
        "--metrics-output",
        default=None,
        help="Optional TSV path for overall and per-track metrics.",
    )
    parser.add_argument(
        "--diagnostic-output",
        default=None,
        help="Optional combined TSV with overall, track, window, and stratum metrics.",
    )
    parser.add_argument(
        "--point-metrics",
        action="store_true",
        help=(
            "Compute exact pointwise MSE/MAE/Pearson and optional sampled "
            "Spearman per track."
        ),
    )
    parser.add_argument(
        "--common-128bp-metrics",
        choices=["raw-sum", "log1p-mean"],
        default=None,
        help=(
            "Compute common 128 bp binned metrics for any head. raw-sum compares "
            "128 bp summed raw signal; log1p-mean compares log1p of 128 bp mean "
            "raw signal. Uses raw NPZ targets regardless of checkpoint target "
            "transform."
        ),
    )
    parser.add_argument(
        "--spearman-sample-size",
        type=int,
        default=0,
        help="Approximate number of sampled positions per track for Spearman.",
    )
    parser.add_argument(
        "--spearman-seed",
        type=int,
        default=0,
        help="Deterministic offset seed for sampled Spearman positions.",
    )
    parser.add_argument(
        "--embedding-resolution",
        type=int,
        choices=[1, 128],
        default=None,
        help="Override checkpoint embedding resolution.",
    )
    parser.add_argument(
        "--head-type",
        choices=["linear", "genome-tracks"],
        default=None,
        help="Override checkpoint head type.",
    )
    parser.add_argument(
        "--head-resolutions",
        default=None,
        help="Override checkpoint head resolutions, e.g. 128 or 1,128.",
    )
    parser.add_argument(
        "--linear-head-architecture",
        choices=[
            "conv1x1",
            "mlp1x1",
            "conv3",
            "conv3x2",
            "conv5",
            "conv7",
            "dilated-conv3",
            "residual-conv1x1",
            "residual-conv3",
        ],
        default=None,
        help="Override checkpoint linear head architecture.",
    )
    parser.add_argument(
        "--linear-hidden-channels",
        type=int,
        default=None,
        help="Override checkpoint linear hidden channels.",
    )
    parser.add_argument(
        "--linear-input-bottleneck-channels",
        type=int,
        default=None,
        help="Override checkpoint linear input bottleneck channels.",
    )
    parser.add_argument(
        "--residual-scale-init",
        type=float,
        default=None,
        help="Override residual scale init before checkpoint load.",
    )
    parser.add_argument(
        "--linear-dilation",
        type=int,
        default=None,
        help="Override checkpoint linear dilation for dilated-conv3.",
    )
    parser.add_argument(
        "--residual-base-checkpoint",
        default=None,
        help="Override checkpoint frozen residual-base checkpoint path.",
    )
    parser.add_argument(
        "--residual-correction-scale-init",
        type=float,
        default=None,
        help="Override residual correction scale init before checkpoint load.",
    )
    parser.add_argument(
        "--linear-target-space",
        choices=["full-log1p", "binned128-log1p-mean"],
        default=None,
        help="Override checkpoint linear target space.",
    )
    parser.add_argument(
        "--linear-loss-type",
        choices=["mse", "smooth-l1", "hybrid", "hybrid-mse-smooth-l1"],
        default=None,
        help="Override checkpoint linear loss type for scalar valid_loss.",
    )
    parser.add_argument("--smooth-l1-beta", type=float, default=None)
    parser.add_argument("--hybrid-loss-alpha", type=float, default=None)
    parser.add_argument(
        "--organism-index",
        type=int,
        default=None,
        help="Override checkpoint organism index.",
    )
    parser.add_argument(
        "--target-transform",
        choices=["log1p", "none"],
        default=None,
        help="Override checkpoint target transform.",
    )
    parser.add_argument(
        "--loss-type",
        choices=["mse", "poisson-multinomial", "hybrid-mse-poisson"],
        default=None,
        help="Override checkpoint loss type for scalar valid_loss evaluation.",
    )
    parser.add_argument("--multinomial-num-segments", type=int, default=None)
    parser.add_argument("--positional-weight", type=float, default=None)
    parser.add_argument("--count-weight", type=float, default=None)
    parser.add_argument("--mse-weight", type=float, default=None)
    parser.add_argument("--poisson-weight", type=float, default=None)
    parser.add_argument(
        "--device",
        default="auto",
        help="auto, cpu, cuda, cuda:0, mps, etc.",
    )
    return parser.parse_args()


def load_track_names(dataset_metadata: dict[str, object], n_tracks: int) -> list[str]:
    metadata_path = dataset_metadata.get("track_metadata")
    if not metadata_path:
        return [f"track_{index}\ttrack_{index}" for index in range(n_tracks)]

    rows: list[tuple[int, str, str]] = []
    with Path(str(metadata_path)).open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if row.get("output_type") != "RNA_SEQ":
                continue
            track_index = int(row["track_index"])
            group_id = row.get("group_id", f"track_{track_index}")
            name = row.get("name", group_id)
            rows.append((track_index, group_id, name))

    rows.sort(key=lambda item: item[0])
    if len(rows) != n_tracks:
        return [f"track_{index}\ttrack_{index}" for index in range(n_tracks)]
    return [f"{group_id}\t{name}" for _track_index, group_id, name in rows]


def write_metrics_tsv(
    path: Path,
    *,
    overall_loss: float,
    per_track_losses: list[float] | None,
    track_labels: list[str],
    batches: int,
    examples: int,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "metric_scope",
                "track_index",
                "track_id",
                "track_name",
                "loss",
                "batches",
                "examples",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerow(
            {
                "metric_scope": "overall",
                "track_index": "",
                "track_id": "",
                "track_name": "",
                "loss": f"{overall_loss:.8g}",
                "batches": batches,
                "examples": examples,
            }
        )
        if per_track_losses is None:
            return
        for index, loss in enumerate(per_track_losses):
            track_id, track_name = track_labels[index].split("\t", maxsplit=1)
            writer.writerow(
                {
                    "metric_scope": "track",
                    "track_index": index,
                    "track_id": track_id,
                    "track_name": track_name,
                    "loss": f"{loss:.8g}",
                    "batches": batches,
                    "examples": examples,
                }
            )


def write_point_metrics_tsv(
    path: Path,
    *,
    metrics: dict[str, object],
    track_labels: list[str],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    spearman_values = metrics["per_track_spearman_sampled"]
    spearman_ns = metrics["per_track_spearman_sampled_n"]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "metric_scope",
                "track_index",
                "track_id",
                "track_name",
                "mse",
                "mae",
                "pearson",
                "spearman_sampled",
                "spearman_sampled_n",
                "batches",
                "examples",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerow(
            {
                "metric_scope": "overall",
                "track_index": "",
                "track_id": "",
                "track_name": "",
                "mse": f"{float(metrics['overall_mse']):.8g}",
                "mae": f"{float(metrics['overall_mae']):.8g}",
                "pearson": f"{float(metrics['overall_pearson']):.8g}",
                "spearman_sampled": "",
                "spearman_sampled_n": "",
                "batches": metrics["batches"],
                "examples": metrics["examples"],
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
                    "mse": f"{metrics['per_track_mse'][index]:.8g}",
                    "mae": f"{metrics['per_track_mae'][index]:.8g}",
                    "pearson": f"{metrics['per_track_pearson'][index]:.8g}",
                    "spearman_sampled": (
                        "" if spearman_values is None
                        else f"{spearman_values[index]:.8g}"
                    ),
                    "spearman_sampled_n": (
                        "" if spearman_ns is None else spearman_ns[index]
                    ),
                    "batches": metrics["batches"],
                    "examples": metrics["examples"],
                }
            )


def metrics_from_sums(
    *,
    sse: float,
    sae: float,
    sum_x: float,
    sum_y: float,
    sum_x2: float,
    sum_y2: float,
    sum_xy: float,
    n: float,
) -> dict[str, float]:
    denominator = max(n, 1.0)
    return {
        "mse": sse / denominator,
        "mae": sae / denominator,
        "pearson": pearson_from_sums(sum_x, sum_y, sum_x2, sum_y2, sum_xy, n),
        "n_values": n,
    }


def evaluate_window_and_strata_metrics(
    model: RnaSeq11Adapter,
    dataloader: torch.utils.data.DataLoader,
    *,
    device: torch.device,
) -> tuple[list[dict[str, object]], dict[str, dict[str, float]]]:
    if model.head_type != "linear":
        raise NotImplementedError("Window/stratified diagnostics currently support linear heads")

    strata = {
        "zero": {
            "predicate": lambda target: target == 0,
            "sse": 0.0,
            "sae": 0.0,
            "sum_x": 0.0,
            "sum_y": 0.0,
            "sum_x2": 0.0,
            "sum_y2": 0.0,
            "sum_xy": 0.0,
            "n": 0.0,
        },
        "0<log1p<=1": {
            "predicate": lambda target: (target > 0) & (target <= 1),
            "sse": 0.0,
            "sae": 0.0,
            "sum_x": 0.0,
            "sum_y": 0.0,
            "sum_x2": 0.0,
            "sum_y2": 0.0,
            "sum_xy": 0.0,
            "n": 0.0,
        },
        "1<log1p<=3": {
            "predicate": lambda target: (target > 1) & (target <= 3),
            "sse": 0.0,
            "sae": 0.0,
            "sum_x": 0.0,
            "sum_y": 0.0,
            "sum_x2": 0.0,
            "sum_y2": 0.0,
            "sum_xy": 0.0,
            "n": 0.0,
        },
        "log1p>3": {
            "predicate": lambda target: target > 3,
            "sse": 0.0,
            "sae": 0.0,
            "sum_x": 0.0,
            "sum_y": 0.0,
            "sum_x2": 0.0,
            "sum_y2": 0.0,
            "sum_xy": 0.0,
            "n": 0.0,
        },
    }
    window_rows: list[dict[str, object]] = []
    model.eval()
    with torch.no_grad():
        for batch in dataloader:
            dna_sequence = batch["dna_sequence"].to(device=device, dtype=torch.float32)
            rna_seq = batch["rna_seq"].to(device=device, dtype=torch.float32)
            rna_seq_mask = batch["rna_seq_mask"].to(device=device)
            prediction = model(dna_sequence, target_length=rna_seq.shape[-1])
            if isinstance(prediction, dict):
                raise NotImplementedError(
                    "Window/stratified diagnostics currently support linear heads"
                )

            mask = rna_seq_mask
            while mask.ndim < prediction.ndim:
                mask = mask.unsqueeze(-1)
            mask = mask.to(dtype=torch.bool, device=prediction.device).expand_as(prediction)
            prediction64 = prediction.to(dtype=torch.float64)
            target64 = rna_seq.to(dtype=torch.float64)
            error64 = prediction64 - target64

            for batch_index in range(int(prediction.shape[0])):
                sample_mask = mask[batch_index]
                n = float(sample_mask.sum().detach().cpu())
                pred_values = prediction64[batch_index][sample_mask]
                target_values = target64[batch_index][sample_mask]
                error_values = pred_values - target_values
                stats = metrics_from_sums(
                    sse=float(error_values.square().sum().detach().cpu()),
                    sae=float(error_values.abs().sum().detach().cpu()),
                    sum_x=float(pred_values.sum().detach().cpu()),
                    sum_y=float(target_values.sum().detach().cpu()),
                    sum_x2=float(pred_values.square().sum().detach().cpu()),
                    sum_y2=float(target_values.square().sum().detach().cpu()),
                    sum_xy=float((pred_values * target_values).sum().detach().cpu()),
                    n=n,
                )
                window_rows.append(
                    {
                        "interval_chromosome": batch["interval_chromosome"][batch_index],
                        "interval_start": int(batch["interval_start"][batch_index]),
                        "interval_end": int(batch["interval_end"][batch_index]),
                        **stats,
                    }
                )

            for stratum_name, stratum in strata.items():
                stratum_mask = stratum["predicate"](target64) & mask
                if not bool(stratum_mask.any()):
                    continue
                pred_values = prediction64[stratum_mask]
                target_values = target64[stratum_mask]
                error_values = pred_values - target_values
                stratum["sse"] += float(error_values.square().sum().detach().cpu())
                stratum["sae"] += float(error_values.abs().sum().detach().cpu())
                stratum["sum_x"] += float(pred_values.sum().detach().cpu())
                stratum["sum_y"] += float(target_values.sum().detach().cpu())
                stratum["sum_x2"] += float(pred_values.square().sum().detach().cpu())
                stratum["sum_y2"] += float(target_values.square().sum().detach().cpu())
                stratum["sum_xy"] += float((pred_values * target_values).sum().detach().cpu())
                stratum["n"] += float(stratum_mask.sum().detach().cpu())

    model.train()
    model.base_model.eval()
    stratum_metrics = {
        name: metrics_from_sums(
            sse=float(stats["sse"]),
            sae=float(stats["sae"]),
            sum_x=float(stats["sum_x"]),
            sum_y=float(stats["sum_y"]),
            sum_x2=float(stats["sum_x2"]),
            sum_y2=float(stats["sum_y2"]),
            sum_xy=float(stats["sum_xy"]),
            n=float(stats["n"]),
        )
        for name, stats in strata.items()
    }
    return window_rows, stratum_metrics


def write_diagnostic_tsv(
    path: Path,
    *,
    point_metrics: dict[str, object],
    window_rows: list[dict[str, object]],
    stratum_metrics: dict[str, dict[str, float]],
    track_labels: list[str],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    spearman_values = point_metrics["per_track_spearman_sampled"]
    spearman_ns = point_metrics["per_track_spearman_sampled_n"]
    fieldnames = [
        "metric_scope",
        "track_index",
        "track_id",
        "track_name",
        "window_index",
        "interval_chromosome",
        "interval_start",
        "interval_end",
        "stratum",
        "mse",
        "mae",
        "pearson",
        "spearman_sampled",
        "spearman_sampled_n",
        "n_values",
    ]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerow(
            {
                "metric_scope": "overall",
                "mse": f"{float(point_metrics['overall_mse']):.8g}",
                "mae": f"{float(point_metrics['overall_mae']):.8g}",
                "pearson": f"{float(point_metrics['overall_pearson']):.8g}",
                "n_values": "",
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
                    "mse": f"{point_metrics['per_track_mse'][index]:.8g}",
                    "mae": f"{point_metrics['per_track_mae'][index]:.8g}",
                    "pearson": f"{point_metrics['per_track_pearson'][index]:.8g}",
                    "spearman_sampled": (
                        "" if spearman_values is None
                        else f"{spearman_values[index]:.8g}"
                    ),
                    "spearman_sampled_n": (
                        "" if spearman_ns is None else spearman_ns[index]
                    ),
                }
            )
        for index, row in enumerate(window_rows):
            writer.writerow(
                {
                    "metric_scope": "window",
                    "window_index": index,
                    "interval_chromosome": row["interval_chromosome"],
                    "interval_start": row["interval_start"],
                    "interval_end": row["interval_end"],
                    "mse": f"{float(row['mse']):.8g}",
                    "mae": f"{float(row['mae']):.8g}",
                    "pearson": f"{float(row['pearson']):.8g}",
                    "n_values": f"{float(row['n_values']):.0f}",
                }
            )
        for stratum_name, stats in stratum_metrics.items():
            writer.writerow(
                {
                    "metric_scope": "stratum",
                    "stratum": stratum_name,
                    "mse": f"{float(stats['mse']):.8g}",
                    "mae": f"{float(stats['mae']):.8g}",
                    "pearson": f"{float(stats['pearson']):.8g}",
                    "n_values": f"{float(stats['n_values']):.0f}",
                }
            )


def main() -> None:
    args = parse_args()
    dataset_dir = Path(args.dataset_dir)
    weights = Path(args.weights)
    checkpoint_path = Path(args.checkpoint)
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
        parse_resolutions(
            args.head_resolutions,
            fallback_resolution=embedding_resolution,
        )
        if args.head_resolutions is not None
        else checkpoint_head_resolutions
    )
    linear_head_architecture = (
        args.linear_head_architecture
        if args.linear_head_architecture is not None
        else str(checkpoint.get("linear_head_architecture", "conv1x1"))
    )
    linear_hidden_channels = (
        args.linear_hidden_channels
        if args.linear_hidden_channels is not None
        else int(checkpoint.get("linear_hidden_channels", 256))
    )
    linear_input_bottleneck_channels = (
        args.linear_input_bottleneck_channels
        if args.linear_input_bottleneck_channels is not None
        else checkpoint.get("linear_input_bottleneck_channels")
    )
    residual_scale_init = (
        args.residual_scale_init
        if args.residual_scale_init is not None
        else float(checkpoint.get("linear_residual_scale_init", 1.0))
    )
    linear_dilation = (
        args.linear_dilation
        if args.linear_dilation is not None
        else int(checkpoint.get("linear_dilation", 2))
    )
    residual_base_checkpoint_path = (
        args.residual_base_checkpoint
        if args.residual_base_checkpoint is not None
        else checkpoint.get("residual_base_checkpoint")
    )
    residual_correction_scale_init = (
        args.residual_correction_scale_init
        if args.residual_correction_scale_init is not None
        else float(checkpoint.get("residual_correction_scale_init", 0.01))
    )
    linear_target_space = (
        args.linear_target_space
        if args.linear_target_space is not None
        else str(checkpoint.get("linear_target_space", "full-log1p"))
    )
    linear_loss_type = (
        args.linear_loss_type
        if args.linear_loss_type is not None
        else str(checkpoint.get("linear_loss_type", "mse"))
    )
    smooth_l1_beta = (
        args.smooth_l1_beta
        if args.smooth_l1_beta is not None
        else float(checkpoint.get("smooth_l1_beta", 1.0))
    )
    hybrid_loss_alpha = (
        args.hybrid_loss_alpha
        if args.hybrid_loss_alpha is not None
        else float(checkpoint.get("hybrid_loss_alpha", 0.5))
    )
    organism_index = (
        args.organism_index
        if args.organism_index is not None
        else int(checkpoint["organism_index"])
    )
    target_transform = (
        args.target_transform
        if args.target_transform is not None
        else str(checkpoint.get("target_transform", "log1p"))
    )
    checkpoint_target_transform = str(checkpoint.get("target_transform", "log1p"))
    loss_type = (
        args.loss_type
        if args.loss_type is not None
        else str(checkpoint.get("loss_type", "mse"))
    )
    multinomial_num_segments = (
        args.multinomial_num_segments
        if args.multinomial_num_segments is not None
        else int(checkpoint.get("multinomial_num_segments", 8))
    )
    positional_weight = (
        args.positional_weight
        if args.positional_weight is not None
        else float(checkpoint.get("positional_weight", 5.0))
    )
    count_weight = (
        args.count_weight
        if args.count_weight is not None
        else float(checkpoint.get("count_weight", 1.0))
    )
    mse_weight = (
        args.mse_weight
        if args.mse_weight is not None
        else float(checkpoint.get("mse_weight", 1.0))
    )
    poisson_weight = (
        args.poisson_weight
        if args.poisson_weight is not None
        else float(checkpoint.get("poisson_weight", 1.0))
    )
    if head_type == "genome-tracks" and target_transform != "none":
        raise ValueError(
            "GenomeTracksHead checkpoints must be evaluated with raw targets: "
            "--target-transform none"
        )
    if loss_type in {"poisson-multinomial", "hybrid-mse-poisson"} and head_type != "genome-tracks":
        raise ValueError(f"{loss_type} evaluation requires genome-tracks")
    if linear_dilation < 1:
        raise ValueError("--linear-dilation must be >= 1")
    if not 0.0 <= hybrid_loss_alpha <= 1.0:
        raise ValueError("--hybrid-loss-alpha must be between 0 and 1")
    if args.diagnostic_output is not None and args.common_128bp_metrics is not None:
        raise ValueError("--diagnostic-output cannot be combined with --common-128bp-metrics")

    if args.common_128bp_metrics is not None:
        dataset_target_transform = "none"
    elif head_type == "linear" and linear_target_space == "binned128-log1p-mean" and not (
        args.point_metrics or args.diagnostic_output is not None
    ):
        dataset_target_transform = "none"
    else:
        dataset_target_transform = target_transform

    dataloader = make_dataloader(
        dataset_dir,
        batch_size=args.batch_size,
        shuffle=False,
        num_workers=args.num_workers,
        target_transform=dataset_target_transform,
        max_examples=args.max_examples,
    )
    n_tracks = int(dataloader.dataset.metadata["n_tracks"])
    checkpoint_n_tracks = int(checkpoint["n_tracks"])
    if checkpoint_n_tracks != n_tracks:
        raise ValueError(
            f"Checkpoint n_tracks={checkpoint_n_tracks}, dataset n_tracks={n_tracks}"
        )

    print(f"dataset_dir\t{dataset_dir}")
    print(f"weights\t{weights}")
    print(f"checkpoint\t{checkpoint_path}")
    print(f"n_examples\t{len(dataloader.dataset)}")
    print(f"n_tracks\t{n_tracks}")
    print(f"target_transform\t{target_transform}")
    print(f"dataset_target_transform\t{dataset_target_transform}")
    print(f"head_type\t{head_type}")
    print(f"head_resolutions\t{','.join(map(str, head_resolutions))}")
    print(f"linear_head_architecture\t{linear_head_architecture}")
    print(f"linear_hidden_channels\t{linear_hidden_channels}")
    print(f"linear_input_bottleneck_channels\t{linear_input_bottleneck_channels}")
    print(f"residual_scale_init\t{residual_scale_init}")
    print(f"linear_dilation\t{linear_dilation}")
    print(f"residual_base_checkpoint\t{residual_base_checkpoint_path}")
    print(f"residual_correction_scale_init\t{residual_correction_scale_init}")
    print(f"linear_target_space\t{linear_target_space}")
    print(f"linear_loss_type\t{linear_loss_type}")
    print(f"smooth_l1_beta\t{smooth_l1_beta}")
    print(f"hybrid_loss_alpha\t{hybrid_loss_alpha}")
    print(f"loss_space\t{checkpoint.get('loss_space', 'target_transform')}")
    print(f"loss_type\t{loss_type}")
    print(f"multinomial_num_segments\t{multinomial_num_segments}")
    print(f"positional_weight\t{positional_weight}")
    print(f"count_weight\t{count_weight}")
    print(f"mse_weight\t{mse_weight}")
    print(f"poisson_weight\t{poisson_weight}")
    if args.common_128bp_metrics is not None:
        print(f"common_128bp_metrics\t{args.common_128bp_metrics}")
    print(f"embedding_resolution\t{embedding_resolution}")
    print(f"organism_index\t{organism_index}")
    print(f"device\t{device}")
    if device.type == "cuda":
        print(f"cuda_device\t{torch.cuda.get_device_name(device)}")

    base_model = AlphaGenome.from_pretrained(weights, device=device)
    checkpoint_head_state = checkpoint["adapter_head_state_dict"]
    track_means = checkpoint_head_state.get("track_means")
    residual_base_head = None
    residual_base_input_bottleneck = None
    residual_base_resolution = 128
    if residual_base_checkpoint_path is not None:
        residual_base_checkpoint = load_adapter_checkpoint(
            str(residual_base_checkpoint_path)
        )
        residual_base_resolution = int(residual_base_checkpoint["embedding_resolution"])
        residual_base_probe = RnaSeq11Adapter(
            base_model,
            n_tracks=n_tracks,
            embedding_resolution=residual_base_resolution,
            organism_index=organism_index,
            head_type="linear",
            head_resolutions=(residual_base_resolution,),
            linear_head_architecture=str(
                residual_base_checkpoint.get("linear_head_architecture", "conv1x1")
            ),
            linear_hidden_channels=int(
                residual_base_checkpoint.get("linear_hidden_channels", 256)
            ),
            linear_input_bottleneck_channels=residual_base_checkpoint.get(
                "linear_input_bottleneck_channels"
            ),
            linear_residual_scale_init=float(
                residual_base_checkpoint.get("linear_residual_scale_init", 1.0)
            ),
            linear_dilation=int(residual_base_checkpoint.get("linear_dilation", 2)),
        ).to(device)
        residual_base_probe.load_adapter_head_state_dict(
            residual_base_checkpoint["adapter_head_state_dict"]
        )
        residual_base_head = residual_base_probe.head
        residual_base_input_bottleneck = residual_base_probe.linear_input_bottleneck
    model = RnaSeq11Adapter(
        base_model,
        n_tracks=n_tracks,
        embedding_resolution=embedding_resolution,
        organism_index=organism_index,
        head_type=head_type,
        head_resolutions=head_resolutions,
        linear_head_architecture=linear_head_architecture,
        linear_hidden_channels=linear_hidden_channels,
        linear_input_bottleneck_channels=linear_input_bottleneck_channels,
        linear_track_means=checkpoint_head_state.get("track_means"),
        linear_residual_scale_init=residual_scale_init,
        linear_dilation=linear_dilation,
        linear_residual_base_head=residual_base_head,
        linear_residual_base_input_bottleneck=residual_base_input_bottleneck,
        linear_residual_base_resolution=residual_base_resolution,
        linear_residual_correction_scale_init=residual_correction_scale_init,
        track_means=track_means,
    ).to(device)
    model.load_adapter_head_state_dict(checkpoint_head_state)

    print(f"base_parameters\t{count_parameters(model.base_model)}")
    print(f"total_parameters\t{count_parameters(model)}")
    print(f"trainable_parameters\t{count_parameters(model, trainable_only=True)}")
    print(f"head_parameters\t{count_parameters(model.head)}")

    per_track_losses = None
    point_metrics = None
    track_labels = load_track_names(dataloader.dataset.metadata, n_tracks)
    if args.common_128bp_metrics is not None:
        point_metrics = evaluate_128bp_binned_pointwise_metrics(
            model,
            dataloader,
            device=device,
            linear_prediction_transform=checkpoint_target_transform,
            linear_target_space=linear_target_space,
            metric_transform=args.common_128bp_metrics,
        )
        loss = float(point_metrics["overall_mse"])
        batches = int(point_metrics["batches"])
        examples = int(point_metrics["examples"])
    elif args.point_metrics or args.diagnostic_output is not None:
        point_metrics = evaluate_pointwise_metrics(
            model,
            dataloader,
            device=device,
            spearman_sample_size=args.spearman_sample_size,
            spearman_seed=args.spearman_seed,
        )
        loss = float(point_metrics["overall_mse"])
        batches = int(point_metrics["batches"])
        examples = int(point_metrics["examples"])
    elif args.per_track:
        loss, per_track_losses, batches, examples = evaluate_masked_mse_by_track(
            model,
            dataloader,
            device=device,
        )
    else:
        loss, batches, examples = evaluate_masked_mse(
            model,
            dataloader,
            device=device,
            loss_type=loss_type,
            multinomial_num_segments=multinomial_num_segments,
            positional_weight=positional_weight,
            count_weight=count_weight,
            mse_weight=mse_weight,
            poisson_weight=poisson_weight,
            linear_loss_type=linear_loss_type,
            linear_target_space=linear_target_space,
            smooth_l1_beta=smooth_l1_beta,
            hybrid_loss_alpha=hybrid_loss_alpha,
        )
    print(f"valid_batches\t{batches}")
    print(f"valid_examples\t{examples}")
    print(f"valid_loss\t{loss:.8g}")
    if point_metrics is not None:
        print(f"valid_mae\t{float(point_metrics['overall_mae']):.8g}")
        print(f"valid_pearson\t{float(point_metrics['overall_pearson']):.8g}")
        if point_metrics["per_track_spearman_sampled"] is not None:
            print(f"spearman_sample_size\t{args.spearman_sample_size}")
            print(f"spearman_seed\t{args.spearman_seed}")
        print(
            "point_metric_header\ttrack_index\ttrack_id\ttrack_name\t"
            "mse\tmae\tpearson\tspearman_sampled\tspearman_sampled_n"
        )
        spearman_values = point_metrics["per_track_spearman_sampled"]
        spearman_ns = point_metrics["per_track_spearman_sampled_n"]
        for index, track_label in enumerate(track_labels):
            track_id, track_name = track_label.split("\t", maxsplit=1)
            spearman_value = "" if spearman_values is None else f"{spearman_values[index]:.8g}"
            spearman_n = "" if spearman_ns is None else str(spearman_ns[index])
            print(
                "point_metric\t"
                f"{index}\t{track_id}\t{track_name}\t"
                f"{point_metrics['per_track_mse'][index]:.8g}\t"
                f"{point_metrics['per_track_mae'][index]:.8g}\t"
                f"{point_metrics['per_track_pearson'][index]:.8g}\t"
                f"{spearman_value}\t{spearman_n}"
            )
    elif per_track_losses is not None:
        print("per_track_loss_header\ttrack_index\ttrack_id\ttrack_name\tloss")
        for index, track_loss in enumerate(per_track_losses):
            track_id, track_name = track_labels[index].split("\t", maxsplit=1)
            print(
                "per_track_loss\t"
                f"{index}\t{track_id}\t{track_name}\t{track_loss:.8g}"
            )
    if args.metrics_output is not None:
        if point_metrics is not None:
            write_point_metrics_tsv(
                Path(args.metrics_output),
                metrics=point_metrics,
                track_labels=track_labels,
            )
        else:
            write_metrics_tsv(
                Path(args.metrics_output),
                overall_loss=loss,
                per_track_losses=per_track_losses,
                track_labels=track_labels,
                batches=batches,
                examples=examples,
            )
        print(f"metrics_output\t{args.metrics_output}")
    if args.diagnostic_output is not None:
        if point_metrics is None:
            raise RuntimeError("diagnostic output requires point metrics")
        window_rows, stratum_metrics = evaluate_window_and_strata_metrics(
            model,
            dataloader,
            device=device,
        )
        write_diagnostic_tsv(
            Path(args.diagnostic_output),
            point_metrics=point_metrics,
            window_rows=window_rows,
            stratum_metrics=stratum_metrics,
            track_labels=track_labels,
        )
        print(f"diagnostic_output\t{args.diagnostic_output}")
        print(f"diagnostic_track_rows\t{len(track_labels)}")
        print(f"diagnostic_window_rows\t{len(window_rows)}")
        print(f"diagnostic_stratum_rows\t{len(stratum_metrics)}")
    if device.type == "cuda":
        print(f"cuda_max_memory_allocated_mb\t{torch.cuda.max_memory_allocated() / 1024**2:.1f}")


if __name__ == "__main__":
    main()
