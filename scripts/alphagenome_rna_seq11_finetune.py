#!/usr/bin/env python3
"""Fine-tune an 11-track C. elegans RNA-seq adapter on frozen AlphaGenome."""

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

from alphagenome_pytorch import AlphaGenome
from alphagenome_rna_seq11_adapter import (
    RnaSeq11Adapter,
    count_parameters,
    evaluate_masked_mse,
    pick_device,
)
from torch_rna_seq_dataset import make_dataloader, masked_mse_loss


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
        "--weights",
        default="weights/alphagenome_pytorch/model_all_folds.safetensors",
        help="Local alphagenome-pytorch safetensors checkpoint.",
    )
    parser.add_argument(
        "--output-dir",
        required=True,
        help="Ignored run directory for config, metrics, and adapter checkpoint.",
    )
    parser.add_argument("--batch-size", type=int, default=1)
    parser.add_argument("--num-workers", type=int, default=0)
    parser.add_argument("--max-steps", type=int, default=20)
    parser.add_argument("--eval-every", type=int, default=5)
    parser.add_argument("--max-train-examples", type=int, default=None)
    parser.add_argument("--max-valid-examples", type=int, default=None)
    parser.add_argument("--learning-rate", type=float, default=1e-3)
    parser.add_argument("--weight-decay", type=float, default=0.0)
    parser.add_argument("--seed", type=int, default=20260514)
    parser.add_argument("--organism-index", type=int, default=0)
    parser.add_argument(
        "--embedding-resolution",
        type=int,
        choices=[1, 128],
        default=128,
        help="Frozen AlphaGenome embedding resolution used by the adapter head.",
    )
    parser.add_argument(
        "--target-transform",
        choices=["log1p", "none"],
        default="log1p",
        help="Use log1p for stable training on high dynamic range targets.",
    )
    parser.add_argument(
        "--device",
        default="auto",
        help="auto, cpu, cuda, cuda:0, mps, etc.",
    )
    parser.add_argument(
        "--no-save-checkpoint",
        action="store_true",
        help="Write config and metrics only; skip adapter checkpoint.",
    )
    return parser.parse_args()


def set_seed(seed: int) -> None:
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)


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


def main() -> None:
    args = parse_args()
    set_seed(args.seed)

    train_dataset_dir = Path(args.train_dataset_dir)
    valid_dataset_dir = Path(args.valid_dataset_dir)
    weights = Path(args.weights)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    config_path = output_dir / "config.json"
    metrics_path = output_dir / "metrics.tsv"
    checkpoint_path = output_dir / "adapter_head.pt"
    best_checkpoint_path = output_dir / "adapter_head_best.pt"
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
    print(f"weights\t{weights}")
    print(f"output_dir\t{output_dir}")
    print(f"config_path\t{config_path}")
    print(f"metrics_path\t{metrics_path}")
    print(f"checkpoint_path\t{checkpoint_path if not args.no_save_checkpoint else 'disabled'}")
    print(
        "best_checkpoint_path\t"
        f"{best_checkpoint_path if not args.no_save_checkpoint else 'disabled'}"
    )
    print(f"n_train_examples\t{len(train_loader.dataset)}")
    print(f"n_valid_examples\t{len(valid_loader.dataset)}")
    print(f"n_tracks\t{n_tracks}")
    print(f"target_transform\t{args.target_transform}")
    print(f"embedding_resolution\t{args.embedding_resolution}")
    print(f"organism_index\t{args.organism_index}")
    print(f"learning_rate\t{args.learning_rate}")
    print(f"weight_decay\t{args.weight_decay}")
    print(f"seed\t{args.seed}")
    print(f"device\t{device}")
    if device.type == "cuda":
        print(f"cuda_device\t{torch.cuda.get_device_name(device)}")

    base_model = AlphaGenome.from_pretrained(weights, device=device)
    model = RnaSeq11Adapter(
        base_model,
        n_tracks=n_tracks,
        embedding_resolution=args.embedding_resolution,
        organism_index=args.organism_index,
    ).to(device)
    optimizer = torch.optim.AdamW(
        model.head.parameters(),
        lr=args.learning_rate,
        weight_decay=args.weight_decay,
    )

    print(f"base_parameters\t{count_parameters(model.base_model)}")
    print(f"total_parameters\t{count_parameters(model)}")
    print(f"trainable_parameters\t{count_parameters(model, trainable_only=True)}")
    print(f"head_parameters\t{count_parameters(model.head)}")

    model.train()
    model.base_model.eval()
    best_valid_loss: float | None = None
    best_valid_step: int | None = None
    for step, batch in enumerate(itertools.cycle(train_loader), start=1):
        dna_sequence = batch["dna_sequence"].to(device=device, dtype=torch.float32)
        rna_seq = batch["rna_seq"].to(device=device, dtype=torch.float32)
        rna_seq_mask = batch["rna_seq_mask"].to(device=device)

        prediction = model(dna_sequence, target_length=rna_seq.shape[-1])
        loss = masked_mse_loss(prediction, rna_seq, rna_seq_mask)

        optimizer.zero_grad(set_to_none=True)
        loss.backward()
        optimizer.step()

        train_loss = float(loss.detach().cpu())
        base_has_grad = any(p.grad is not None for p in model.base_model.parameters())
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
            lr=args.learning_rate,
            cuda_max_memory_allocated_mb=cuda_memory,
        )

        print(f"step\t{step}")
        print(f"train_loss\t{train_loss:.6g}")
        print(f"dna_sequence_shape\t{'x'.join(map(str, dna_sequence.shape))}")
        print(f"rna_seq_shape\t{'x'.join(map(str, rna_seq.shape))}")
        print(f"prediction_shape\t{'x'.join(map(str, prediction.shape))}")
        print(f"base_has_grad\t{base_has_grad}")
        print(f"head_weight_grad_norm\t{float(model.head.weight.grad.norm().detach().cpu()):.6g}")
        if cuda_memory is not None:
            print(f"cuda_max_memory_allocated_mb\t{cuda_memory:.1f}")

        if step % args.eval_every == 0 or step == args.max_steps:
            valid_loss, valid_batches, valid_examples = evaluate_masked_mse(
                model,
                valid_loader,
                device=device,
            )
            cuda_memory = (
                torch.cuda.max_memory_allocated() / 1024**2
                if device.type == "cuda" else None
            )
            append_metric(
                metrics_path,
                split="valid",
                step=step,
                loss=valid_loss,
                batches=valid_batches,
                examples=valid_examples,
                lr=args.learning_rate,
                cuda_max_memory_allocated_mb=cuda_memory,
            )
            print(f"valid_step\t{step}")
            print(f"valid_batches\t{valid_batches}")
            print(f"valid_examples\t{valid_examples}")
            print(f"valid_loss\t{valid_loss:.6g}")
            if (
                not args.no_save_checkpoint
                and (best_valid_loss is None or valid_loss < best_valid_loss)
            ):
                best_valid_loss = valid_loss
                best_valid_step = step
                torch.save(
                    {
                        "adapter_head_state_dict": model.head.state_dict(),
                        "n_tracks": n_tracks,
                        "embedding_resolution": args.embedding_resolution,
                        "organism_index": args.organism_index,
                        "target_transform": args.target_transform,
                        "base_weights": str(weights),
                        "best_valid_loss": best_valid_loss,
                        "best_valid_step": best_valid_step,
                    },
                    best_checkpoint_path,
                )
                print(f"saved_best_checkpoint\t{best_checkpoint_path}")
                print(f"best_valid_step\t{best_valid_step}")
                print(f"best_valid_loss\t{best_valid_loss:.6g}")

        if step >= args.max_steps:
            break

    if not args.no_save_checkpoint:
        torch.save(
            {
                "adapter_head_state_dict": model.head.state_dict(),
                "n_tracks": n_tracks,
                "embedding_resolution": args.embedding_resolution,
                "organism_index": args.organism_index,
                "target_transform": args.target_transform,
                "base_weights": str(weights),
                "best_valid_loss": best_valid_loss,
                "best_valid_step": best_valid_step,
            },
            checkpoint_path,
        )
        print(f"saved_checkpoint\t{checkpoint_path}")


if __name__ == "__main__":
    main()
