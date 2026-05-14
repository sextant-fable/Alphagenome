#!/usr/bin/env python3
"""PyTorch Dataset utilities for the custom AlphaGenome-like RNA_SEQ NPZ files."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Literal

import numpy as np
import torch
from torch.utils.data import DataLoader, Dataset


TargetTransform = Literal["none", "log1p"]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


class AlphaGenomeRnaSeqNpzDataset(Dataset):
    """Loads one 1 Mb interval per NPZ file.

    Returned tensors are channel-first by default:

    - dna_sequence: [4, S]
    - rna_seq: [C, S]
    - rna_seq_mask: [C, 1]
    - rna_seq_strand: [C]
    """

    def __init__(
        self,
        dataset_dir: str | Path,
        *,
        channels_first: bool = True,
        target_transform: TargetTransform = "log1p",
        max_examples: int | None = None,
    ) -> None:
        self.dataset_dir = Path(dataset_dir)
        self.channels_first = channels_first
        self.target_transform = target_transform
        self.manifest = read_tsv(self.dataset_dir / "manifest.tsv")
        if max_examples is not None:
            self.manifest = self.manifest[:max_examples]
        with (self.dataset_dir / "dataset_metadata.json").open() as handle:
            self.metadata = json.load(handle)

    def __len__(self) -> int:
        return len(self.manifest)

    def __getitem__(self, index: int) -> dict[str, object]:
        row = self.manifest[index]
        with np.load(self.dataset_dir / row["path"], allow_pickle=False) as data:
            dna_sequence = data["dna_sequence"].astype(np.float32, copy=False)
            rna_seq = data["rna_seq"].astype(np.float32, copy=False)
            rna_seq_mask = data["rna_seq_mask"].astype(bool, copy=False)
            rna_seq_strand = data["rna_seq_strand"].astype(np.int64, copy=False)
            chrom = str(data["interval_chromosome"])
            start = int(data["interval_start"])
            end = int(data["interval_end"])

        if self.target_transform == "log1p":
            rna_seq = np.log1p(np.clip(rna_seq, a_min=0.0, a_max=None))
        elif self.target_transform != "none":
            raise ValueError(f"Unknown target_transform: {self.target_transform}")

        if self.channels_first:
            dna_sequence = np.ascontiguousarray(dna_sequence.T)
            rna_seq = np.ascontiguousarray(rna_seq.T)
            rna_seq_mask = np.ascontiguousarray(rna_seq_mask.reshape(-1, 1))
            rna_seq_strand = np.ascontiguousarray(rna_seq_strand.reshape(-1))
        else:
            dna_sequence = np.ascontiguousarray(dna_sequence)
            rna_seq = np.ascontiguousarray(rna_seq)

        return {
            "dna_sequence": torch.from_numpy(dna_sequence),
            "rna_seq": torch.from_numpy(rna_seq),
            "rna_seq_mask": torch.from_numpy(rna_seq_mask),
            "rna_seq_strand": torch.from_numpy(rna_seq_strand),
            "interval_chromosome": chrom,
            "interval_start": torch.tensor(start, dtype=torch.long),
            "interval_end": torch.tensor(end, dtype=torch.long),
        }


def make_dataloader(
    dataset_dir: str | Path,
    *,
    batch_size: int = 1,
    shuffle: bool = False,
    num_workers: int = 0,
    target_transform: TargetTransform = "log1p",
    max_examples: int | None = None,
) -> DataLoader:
    dataset = AlphaGenomeRnaSeqNpzDataset(
        dataset_dir,
        target_transform=target_transform,
        max_examples=max_examples,
    )
    return DataLoader(
        dataset,
        batch_size=batch_size,
        shuffle=shuffle,
        num_workers=num_workers,
        pin_memory=torch.cuda.is_available(),
    )


def masked_mse_loss(
    prediction: torch.Tensor,
    target: torch.Tensor,
    mask: torch.Tensor,
) -> torch.Tensor:
    """Mean squared error over valid RNA_SEQ tracks.

    prediction/target: [B, C, S]
    mask: [B, C, 1] or [B, C]
    """

    while mask.ndim < prediction.ndim:
        mask = mask.unsqueeze(-1)
    mask = mask.to(dtype=prediction.dtype, device=prediction.device)
    squared_error = (prediction - target).square() * mask
    denominator = mask.sum().clamp_min(1.0) * prediction.shape[-1]
    return squared_error.sum() / denominator


class TinyConvRnaSeqModel(torch.nn.Module):
    """Tiny same-length Conv1d baseline for dataloader smoke tests."""

    def __init__(self, n_tracks: int, hidden_channels: int = 16) -> None:
        super().__init__()
        self.net = torch.nn.Sequential(
            torch.nn.Conv1d(4, hidden_channels, kernel_size=15, padding=7),
            torch.nn.GELU(),
            torch.nn.Conv1d(hidden_channels, hidden_channels, kernel_size=15, padding=7),
            torch.nn.GELU(),
            torch.nn.Conv1d(hidden_channels, n_tracks, kernel_size=1),
        )

    def forward(self, dna_sequence: torch.Tensor) -> torch.Tensor:
        return self.net(dna_sequence)
