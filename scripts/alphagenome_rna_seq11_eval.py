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
        choices=["mse", "poisson-multinomial"],
        default=None,
        help="Override checkpoint loss type for scalar valid_loss evaluation.",
    )
    parser.add_argument("--multinomial-num-segments", type=int, default=None)
    parser.add_argument("--positional-weight", type=float, default=None)
    parser.add_argument("--count-weight", type=float, default=None)
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
    if head_type == "genome-tracks" and target_transform != "none":
        raise ValueError(
            "GenomeTracksHead checkpoints must be evaluated with raw targets: "
            "--target-transform none"
        )
    if loss_type == "poisson-multinomial" and head_type != "genome-tracks":
        raise ValueError("poisson-multinomial evaluation requires genome-tracks")

    dataset_target_transform = (
        "none" if args.common_128bp_metrics is not None else target_transform
    )

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
    print(f"loss_space\t{checkpoint.get('loss_space', 'target_transform')}")
    print(f"loss_type\t{loss_type}")
    print(f"multinomial_num_segments\t{multinomial_num_segments}")
    print(f"positional_weight\t{positional_weight}")
    print(f"count_weight\t{count_weight}")
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
            metric_transform=args.common_128bp_metrics,
        )
        loss = float(point_metrics["overall_mse"])
        batches = int(point_metrics["batches"])
        examples = int(point_metrics["examples"])
    elif args.point_metrics:
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
    if device.type == "cuda":
        print(f"cuda_max_memory_allocated_mb\t{torch.cuda.max_memory_allocated() / 1024**2:.1f}")


if __name__ == "__main__":
    main()
