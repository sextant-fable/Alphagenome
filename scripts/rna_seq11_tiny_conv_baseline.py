#!/usr/bin/env python3
"""Train a tiny Conv1d learned baseline on C. elegans RNA-seq NPZ data."""

from __future__ import annotations

import argparse
import csv
import itertools
import json
import random
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import torch

from torch_rna_seq_dataset import (
    TinyConvRnaSeqModel,
    make_dataloader,
    masked_mse_loss,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--train-dataset-dir",
        default="alphagenome_custom/datasets/rna_seq_npz_train",
        help="Training NPZ dataset directory.",
    )
    parser.add_argument(
        "--valid-dataset-dir",
        default="alphagenome_custom/datasets/rna_seq_npz_valid",
        help="Validation NPZ dataset directory.",
    )
    parser.add_argument(
        "--output-dir",
        required=True,
        help="Ignored run directory for config, metrics, and checkpoint.",
    )
    parser.add_argument("--batch-size", type=int, default=1)
    parser.add_argument("--num-workers", type=int, default=0)
    parser.add_argument("--max-steps", type=int, default=500)
    parser.add_argument("--eval-every", type=int, default=100)
    parser.add_argument("--hidden-channels", type=int, default=16)
    parser.add_argument("--max-train-examples", type=int, default=None)
    parser.add_argument("--max-valid-examples", type=int, default=None)
    parser.add_argument("--learning-rate", type=float, default=3e-4)
    parser.add_argument("--weight-decay", type=float, default=0.0)
    parser.add_argument("--seed", type=int, default=20260515)
    parser.add_argument(
        "--target-transform",
        choices=["log1p", "none"],
        default="log1p",
        help="Target transform to match adapter training/evaluation.",
    )
    parser.add_argument(
        "--spearman-sample-size",
        type=int,
        default=200000,
        help="Approximate sampled positions per track for final Spearman.",
    )
    parser.add_argument(
        "--spearman-seed",
        type=int,
        default=20260515,
        help="Deterministic offset seed for sampled Spearman positions.",
    )
    parser.add_argument(
        "--device",
        default="auto",
        help="auto, cpu, cuda, cuda:0, mps, etc.",
    )
    parser.add_argument(
        "--no-save-checkpoint",
        action="store_true",
        help="Write config and metrics only; skip model checkpoint.",
    )
    return parser.parse_args()


def pick_device(requested: str) -> torch.device:
    if requested != "auto":
        return torch.device(requested)
    if torch.cuda.is_available():
        return torch.device("cuda")
    if hasattr(torch.backends, "mps") and torch.backends.mps.is_available():
        return torch.device("mps")
    return torch.device("cpu")


def set_seed(seed: int) -> None:
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)


def count_parameters(module: torch.nn.Module, *, trainable_only: bool = False) -> int:
    parameters = module.parameters()
    if trainable_only:
        parameters = (p for p in parameters if p.requires_grad)
    return sum(p.numel() for p in parameters)


def pearson_from_sums(
    sum_x: float,
    sum_y: float,
    sum_x2: float,
    sum_y2: float,
    sum_xy: float,
    n: float,
) -> float:
    if n <= 1:
        return float("nan")
    numerator = n * sum_xy - sum_x * sum_y
    denominator_x = n * sum_x2 - sum_x * sum_x
    denominator_y = n * sum_y2 - sum_y * sum_y
    denominator = float(np.sqrt(max(denominator_x, 0.0) * max(denominator_y, 0.0)))
    if denominator == 0.0:
        return float("nan")
    return float(numerator / denominator)


def rankdata_average(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values)
    sorter = np.argsort(values, kind="mergesort")
    sorted_values = values[sorter]
    ranks = np.empty(len(values), dtype=np.float64)
    if len(values) == 0:
        return ranks

    boundaries = np.concatenate(
        (
            np.array([0]),
            np.flatnonzero(sorted_values[1:] != sorted_values[:-1]) + 1,
            np.array([len(values)]),
        )
    )
    for start, end in zip(boundaries[:-1], boundaries[1:]):
        average_rank = 0.5 * (start + end - 1) + 1.0
        ranks[sorter[start:end]] = average_rank
    return ranks


def sampled_spearman(prediction: np.ndarray, target: np.ndarray) -> float:
    if len(prediction) <= 1:
        return float("nan")
    prediction_ranks = rankdata_average(prediction)
    target_ranks = rankdata_average(target)
    return pearson_from_sums(
        float(prediction_ranks.sum()),
        float(target_ranks.sum()),
        float(np.square(prediction_ranks).sum()),
        float(np.square(target_ranks).sum()),
        float((prediction_ranks * target_ranks).sum()),
        float(len(prediction_ranks)),
    )


def load_track_labels(
    dataset_metadata: dict[str, object],
    n_tracks: int,
) -> list[tuple[str, str]]:
    metadata_path = dataset_metadata.get("track_metadata")
    if not metadata_path:
        return [(f"track_{index}", f"track_{index}") for index in range(n_tracks)]

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
        return [(f"track_{index}", f"track_{index}") for index in range(n_tracks)]
    return [(group_id, name) for _track_index, group_id, name in rows]


def write_config(
    path: Path,
    *,
    args: argparse.Namespace,
    device: torch.device,
    n_train_examples: int,
    n_valid_examples: int,
    n_tracks: int,
) -> None:
    config = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "args": vars(args),
        "device": str(device),
        "torch_version": torch.__version__,
        "torch_cuda": torch.version.cuda,
        "cuda_device": (
            torch.cuda.get_device_name(device) if device.type == "cuda" else None
        ),
        "n_train_examples": n_train_examples,
        "n_valid_examples": n_valid_examples,
        "n_tracks": n_tracks,
        "model": "TinyConvRnaSeqModel",
    }
    with path.open("w") as handle:
        json.dump(config, handle, indent=2, sort_keys=True)
        handle.write("\n")


def append_metric(
    path: Path,
    *,
    split: str,
    step: int,
    loss: float,
    batches: int,
    examples: int,
    lr: float,
    cuda_max_memory_allocated_mb: float | None,
) -> None:
    fieldnames = [
        "split",
        "step",
        "loss",
        "batches",
        "examples",
        "learning_rate",
        "cuda_max_memory_allocated_mb",
    ]
    write_header = not path.exists()
    with path.open("a", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        if write_header:
            writer.writeheader()
        writer.writerow(
            {
                "split": split,
                "step": step,
                "loss": f"{loss:.8g}",
                "batches": batches,
                "examples": examples,
                "learning_rate": f"{lr:.8g}",
                "cuda_max_memory_allocated_mb": (
                    "" if cuda_max_memory_allocated_mb is None
                    else f"{cuda_max_memory_allocated_mb:.1f}"
                ),
            }
        )


def evaluate_masked_mse(
    model: torch.nn.Module,
    dataloader: torch.utils.data.DataLoader,
    *,
    device: torch.device,
) -> tuple[float, int, int]:
    model.eval()
    weighted_loss_sum = 0.0
    n_batches = 0
    n_examples = 0
    with torch.no_grad():
        for batch in dataloader:
            dna_sequence = batch["dna_sequence"].to(device=device, dtype=torch.float32)
            rna_seq = batch["rna_seq"].to(device=device, dtype=torch.float32)
            rna_seq_mask = batch["rna_seq_mask"].to(device=device)

            prediction = model(dna_sequence)
            loss = masked_mse_loss(prediction, rna_seq, rna_seq_mask)
            batch_size = int(dna_sequence.shape[0])
            weighted_loss_sum += float(loss.detach().cpu()) * batch_size
            n_batches += 1
            n_examples += batch_size

    model.train()
    return weighted_loss_sum / max(n_examples, 1), n_batches, n_examples


def evaluate_pointwise_metrics(
    model: torch.nn.Module,
    dataloader: torch.utils.data.DataLoader,
    *,
    device: torch.device,
    spearman_sample_size: int,
    spearman_seed: int,
) -> dict[str, object]:
    model.eval()
    n_tracks: int | None = None
    sse: torch.Tensor | None = None
    sae: torch.Tensor | None = None
    sum_x: torch.Tensor | None = None
    sum_y: torch.Tensor | None = None
    sum_x2: torch.Tensor | None = None
    sum_y2: torch.Tensor | None = None
    sum_xy: torch.Tensor | None = None
    denominator: torch.Tensor | None = None
    prediction_samples: list[list[np.ndarray]] | None = None
    target_samples: list[list[np.ndarray]] | None = None
    sample_stride: int | None = None
    sample_offset = 0
    example_offset = 0
    n_batches = 0
    n_examples = 0

    with torch.no_grad():
        for batch in dataloader:
            dna_sequence = batch["dna_sequence"].to(device=device, dtype=torch.float32)
            rna_seq = batch["rna_seq"].to(device=device, dtype=torch.float32)
            rna_seq_mask = batch["rna_seq_mask"].to(device=device)
            prediction = model(dna_sequence)

            if n_tracks is None:
                n_tracks = int(prediction.shape[1])
                sse = torch.zeros(n_tracks, dtype=torch.float64)
                sae = torch.zeros(n_tracks, dtype=torch.float64)
                sum_x = torch.zeros(n_tracks, dtype=torch.float64)
                sum_y = torch.zeros(n_tracks, dtype=torch.float64)
                sum_x2 = torch.zeros(n_tracks, dtype=torch.float64)
                sum_y2 = torch.zeros(n_tracks, dtype=torch.float64)
                sum_xy = torch.zeros(n_tracks, dtype=torch.float64)
                denominator = torch.zeros(n_tracks, dtype=torch.float64)
                if spearman_sample_size > 0:
                    prediction_samples = [[] for _ in range(n_tracks)]
                    target_samples = [[] for _ in range(n_tracks)]

            mask = rna_seq_mask
            while mask.ndim < prediction.ndim:
                mask = mask.unsqueeze(-1)
            mask = mask.to(dtype=torch.bool, device=prediction.device)
            mask64 = mask.to(dtype=torch.float64)
            prediction64 = prediction.to(dtype=torch.float64)
            target64 = rna_seq.to(dtype=torch.float64)
            error64 = prediction64 - target64

            assert sse is not None
            assert sae is not None
            assert sum_x is not None
            assert sum_y is not None
            assert sum_x2 is not None
            assert sum_y2 is not None
            assert sum_xy is not None
            assert denominator is not None
            sse += (error64.square() * mask64).sum(dim=(0, 2)).detach().cpu()
            sae += (error64.abs() * mask64).sum(dim=(0, 2)).detach().cpu()
            sum_x += (prediction64 * mask64).sum(dim=(0, 2)).detach().cpu()
            sum_y += (target64 * mask64).sum(dim=(0, 2)).detach().cpu()
            sum_x2 += (prediction64.square() * mask64).sum(dim=(0, 2)).detach().cpu()
            sum_y2 += (target64.square() * mask64).sum(dim=(0, 2)).detach().cpu()
            sum_xy += (prediction64 * target64 * mask64).sum(dim=(0, 2)).detach().cpu()
            denominator += mask64.expand_as(prediction64).sum(dim=(0, 2)).detach().cpu()

            if spearman_sample_size > 0:
                assert prediction_samples is not None
                assert target_samples is not None
                batch_size = int(prediction.shape[0])
                sequence_length = int(prediction.shape[-1])
                if sample_stride is None:
                    total_positions = max(len(dataloader.dataset) * sequence_length, 1)
                    sample_stride = max(total_positions // spearman_sample_size, 1)
                    sample_offset = spearman_seed % sample_stride

                local_positions = torch.arange(sequence_length, device=device)
                example_indices = torch.arange(
                    example_offset,
                    example_offset + batch_size,
                    device=device,
                ).unsqueeze(1)
                global_positions = example_indices * sequence_length + local_positions
                sample_position_mask = (
                    (global_positions + sample_offset) % sample_stride == 0
                )
                for track_index in range(n_tracks):
                    track_mask = sample_position_mask & mask[:, track_index, :]
                    if not bool(track_mask.any()):
                        continue
                    prediction_samples[track_index].append(
                        prediction[:, track_index, :][track_mask]
                        .detach()
                        .cpu()
                        .numpy()
                        .astype(np.float64, copy=False)
                    )
                    target_samples[track_index].append(
                        rna_seq[:, track_index, :][track_mask]
                        .detach()
                        .cpu()
                        .numpy()
                        .astype(np.float64, copy=False)
                    )

            batch_size = int(dna_sequence.shape[0])
            example_offset += batch_size
            n_batches += 1
            n_examples += batch_size

    if denominator is None:
        raise ValueError("Cannot evaluate an empty dataloader")

    per_track_mse = sse / denominator.clamp_min(1.0)
    per_track_mae = sae / denominator.clamp_min(1.0)
    per_track_pearson = [
        pearson_from_sums(
            float(sum_x[index]),
            float(sum_y[index]),
            float(sum_x2[index]),
            float(sum_y2[index]),
            float(sum_xy[index]),
            float(denominator[index]),
        )
        for index in range(len(denominator))
    ]

    per_track_spearman: list[float] | None = None
    per_track_spearman_n: list[int] | None = None
    if prediction_samples is not None and target_samples is not None:
        per_track_spearman = []
        per_track_spearman_n = []
        for track_prediction_samples, track_target_samples in zip(
            prediction_samples,
            target_samples,
        ):
            prediction_sample = np.concatenate(track_prediction_samples)
            target_sample = np.concatenate(track_target_samples)
            per_track_spearman_n.append(int(len(prediction_sample)))
            per_track_spearman.append(sampled_spearman(prediction_sample, target_sample))

    overall_n = float(denominator.sum())
    overall_mse = float(sse.sum() / denominator.sum().clamp_min(1.0))
    overall_mae = float(sae.sum() / denominator.sum().clamp_min(1.0))
    overall_pearson = pearson_from_sums(
        float(sum_x.sum()),
        float(sum_y.sum()),
        float(sum_x2.sum()),
        float(sum_y2.sum()),
        float(sum_xy.sum()),
        overall_n,
    )

    model.train()
    return {
        "overall_mse": overall_mse,
        "overall_mae": overall_mae,
        "overall_pearson": overall_pearson,
        "per_track_mse": [float(value) for value in per_track_mse],
        "per_track_mae": [float(value) for value in per_track_mae],
        "per_track_pearson": per_track_pearson,
        "per_track_spearman_sampled": per_track_spearman,
        "per_track_spearman_sampled_n": per_track_spearman_n,
        "batches": n_batches,
        "examples": n_examples,
    }


def write_point_metrics_tsv(
    path: Path,
    *,
    metrics: dict[str, object],
    track_labels: list[tuple[str, str]],
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
        for index, (track_id, track_name) in enumerate(track_labels):
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
    set_seed(args.seed)

    train_dataset_dir = Path(args.train_dataset_dir)
    valid_dataset_dir = Path(args.valid_dataset_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    config_path = output_dir / "config.json"
    metrics_path = output_dir / "metrics.tsv"
    point_metrics_path = output_dir / "valid_pointwise_metrics.tsv"
    checkpoint_path = output_dir / "tiny_conv.pt"
    device = pick_device(args.device)

    train_loader = make_dataloader(
        train_dataset_dir,
        batch_size=args.batch_size,
        shuffle=True,
        num_workers=args.num_workers,
        target_transform=args.target_transform,
        max_examples=args.max_train_examples,
    )
    valid_loader = make_dataloader(
        valid_dataset_dir,
        batch_size=args.batch_size,
        shuffle=False,
        num_workers=args.num_workers,
        target_transform=args.target_transform,
        max_examples=args.max_valid_examples,
    )
    n_tracks = int(train_loader.dataset.metadata["n_tracks"])
    valid_n_tracks = int(valid_loader.dataset.metadata["n_tracks"])
    if valid_n_tracks != n_tracks:
        raise ValueError(f"Train n_tracks={n_tracks}, valid n_tracks={valid_n_tracks}")

    write_config(
        config_path,
        args=args,
        device=device,
        n_train_examples=len(train_loader.dataset),
        n_valid_examples=len(valid_loader.dataset),
        n_tracks=n_tracks,
    )

    print(f"train_dataset_dir\t{train_dataset_dir}")
    print(f"valid_dataset_dir\t{valid_dataset_dir}")
    print(f"output_dir\t{output_dir}")
    print(f"config_path\t{config_path}")
    print(f"metrics_path\t{metrics_path}")
    print(f"point_metrics_path\t{point_metrics_path}")
    print(f"checkpoint_path\t{checkpoint_path if not args.no_save_checkpoint else 'disabled'}")
    print(f"n_train_examples\t{len(train_loader.dataset)}")
    print(f"n_valid_examples\t{len(valid_loader.dataset)}")
    print(f"n_tracks\t{n_tracks}")
    print(f"target_transform\t{args.target_transform}")
    print(f"hidden_channels\t{args.hidden_channels}")
    print(f"learning_rate\t{args.learning_rate}")
    print(f"weight_decay\t{args.weight_decay}")
    print(f"seed\t{args.seed}")
    print(f"device\t{device}")
    if device.type == "cuda":
        print(f"cuda_device\t{torch.cuda.get_device_name(device)}")

    model = TinyConvRnaSeqModel(
        n_tracks=n_tracks,
        hidden_channels=args.hidden_channels,
    ).to(device)
    optimizer = torch.optim.AdamW(
        model.parameters(),
        lr=args.learning_rate,
        weight_decay=args.weight_decay,
    )

    print(f"total_parameters\t{count_parameters(model)}")
    print(f"trainable_parameters\t{count_parameters(model, trainable_only=True)}")

    model.train()
    for step, batch in enumerate(itertools.cycle(train_loader), start=1):
        dna_sequence = batch["dna_sequence"].to(device=device, dtype=torch.float32)
        rna_seq = batch["rna_seq"].to(device=device, dtype=torch.float32)
        rna_seq_mask = batch["rna_seq_mask"].to(device=device)

        prediction = model(dna_sequence)
        loss = masked_mse_loss(prediction, rna_seq, rna_seq_mask)

        optimizer.zero_grad(set_to_none=True)
        loss.backward()
        optimizer.step()

        train_loss = float(loss.detach().cpu())
        cuda_memory = (
            torch.cuda.max_memory_allocated() / 1024**2
            if device.type == "cuda" else None
        )
        append_metric(
            metrics_path,
            split="train",
            step=step,
            loss=train_loss,
            batches=1,
            examples=int(dna_sequence.shape[0]),
            lr=optimizer.param_groups[0]["lr"],
            cuda_max_memory_allocated_mb=cuda_memory,
        )
        print(f"train_step\t{step}\tloss\t{train_loss:.8g}")

        if step % args.eval_every == 0 or step == args.max_steps:
            valid_loss, valid_batches, valid_examples = evaluate_masked_mse(
                model,
                valid_loader,
                device=device,
            )
            append_metric(
                metrics_path,
                split="valid",
                step=step,
                loss=valid_loss,
                batches=valid_batches,
                examples=valid_examples,
                lr=optimizer.param_groups[0]["lr"],
                cuda_max_memory_allocated_mb=cuda_memory,
            )
            print(
                "valid_eval\t"
                f"step\t{step}\tloss\t{valid_loss:.8g}\t"
                f"batches\t{valid_batches}\texamples\t{valid_examples}"
            )
            model.train()

        if step >= args.max_steps:
            break

    if not args.no_save_checkpoint:
        torch.save(
            {
                "model_state_dict": model.state_dict(),
                "n_tracks": n_tracks,
                "hidden_channels": args.hidden_channels,
                "target_transform": args.target_transform,
                "max_steps": args.max_steps,
                "learning_rate": args.learning_rate,
                "weight_decay": args.weight_decay,
                "seed": args.seed,
            },
            checkpoint_path,
        )
        print(f"checkpoint_saved\t{checkpoint_path}")

    point_metrics = evaluate_pointwise_metrics(
        model,
        valid_loader,
        device=device,
        spearman_sample_size=args.spearman_sample_size,
        spearman_seed=args.spearman_seed,
    )
    track_labels = load_track_labels(valid_loader.dataset.metadata, n_tracks)
    write_point_metrics_tsv(
        point_metrics_path,
        metrics=point_metrics,
        track_labels=track_labels,
    )
    print(f"valid_mse\t{float(point_metrics['overall_mse']):.8g}")
    print(f"valid_mae\t{float(point_metrics['overall_mae']):.8g}")
    print(f"valid_pearson\t{float(point_metrics['overall_pearson']):.8g}")
    print(f"point_metrics_saved\t{point_metrics_path}")


if __name__ == "__main__":
    main()
