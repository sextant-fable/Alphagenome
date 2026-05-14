#!/usr/bin/env python3
"""Smoke-train a tiny PyTorch model on the custom AlphaGenome-like NPZ dataset."""

from __future__ import annotations

import argparse
from pathlib import Path

import torch

from torch_rna_seq_dataset import (
    TinyConvRnaSeqModel,
    make_dataloader,
    masked_mse_loss,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--dataset-dir",
        default="alphagenome_custom/datasets/rna_seq_npz_pilot_train",
        help="NPZ dataset directory containing manifest.tsv and examples/.",
    )
    parser.add_argument("--batch-size", type=int, default=1)
    parser.add_argument("--num-workers", type=int, default=0)
    parser.add_argument("--max-steps", type=int, default=2)
    parser.add_argument("--hidden-channels", type=int, default=16)
    parser.add_argument(
        "--target-transform",
        choices=["log1p", "none"],
        default="log1p",
        help="Use log1p for stable smoke training on high dynamic range bigWig targets.",
    )
    parser.add_argument(
        "--device",
        default="auto",
        help="auto, cpu, cuda, cuda:0, mps, etc.",
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


def main() -> None:
    args = parse_args()
    dataset_dir = Path(args.dataset_dir)
    device = pick_device(args.device)

    dataloader = make_dataloader(
        dataset_dir,
        batch_size=args.batch_size,
        shuffle=True,
        num_workers=args.num_workers,
        target_transform=args.target_transform,
    )
    n_tracks = int(dataloader.dataset.metadata["n_tracks"])
    model = TinyConvRnaSeqModel(
        n_tracks=n_tracks,
        hidden_channels=args.hidden_channels,
    ).to(device)
    optimizer = torch.optim.AdamW(model.parameters(), lr=1e-3)

    print(f"dataset_dir\t{dataset_dir}")
    print(f"n_examples\t{len(dataloader.dataset)}")
    print(f"n_tracks\t{n_tracks}")
    print(f"target_transform\t{args.target_transform}")
    print(f"device\t{device}")

    model.train()
    for step, batch in enumerate(dataloader, start=1):
        dna_sequence = batch["dna_sequence"].to(device=device, dtype=torch.float32)
        rna_seq = batch["rna_seq"].to(device=device, dtype=torch.float32)
        rna_seq_mask = batch["rna_seq_mask"].to(device=device)

        prediction = model(dna_sequence)
        loss = masked_mse_loss(prediction, rna_seq, rna_seq_mask)

        optimizer.zero_grad(set_to_none=True)
        loss.backward()
        optimizer.step()

        print(f"step\t{step}")
        print(f"dna_sequence_shape\t{'x'.join(map(str, dna_sequence.shape))}")
        print(f"rna_seq_shape\t{'x'.join(map(str, rna_seq.shape))}")
        print(f"prediction_shape\t{'x'.join(map(str, prediction.shape))}")
        print(f"loss\t{float(loss.detach().cpu()):.6g}")

        if step >= args.max_steps:
            break


if __name__ == "__main__":
    main()
