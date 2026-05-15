#!/usr/bin/env python3
"""Shared frozen AlphaGenome adapter utilities for C. elegans RNA-seq."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import torch
import torch.nn.functional as F

from alphagenome_pytorch import AlphaGenome
from torch_rna_seq_dataset import masked_mse_loss


class RnaSeq11Adapter(torch.nn.Module):
    """Project frozen AlphaGenome sequence embeddings to custom RNA-seq tracks."""

    def __init__(
        self,
        base_model: AlphaGenome,
        *,
        n_tracks: int,
        embedding_resolution: int = 128,
        organism_index: int = 0,
    ) -> None:
        super().__init__()
        if embedding_resolution not in {1, 128}:
            raise ValueError("embedding_resolution must be 1 or 128")

        self.base_model = base_model
        self.embedding_resolution = embedding_resolution
        self.organism_index = organism_index
        in_channels = 1536 if embedding_resolution == 1 else 3072
        self.head = torch.nn.Conv1d(in_channels, n_tracks, kernel_size=1)

        for parameter in self.base_model.parameters():
            parameter.requires_grad = False
        self.base_model.eval()

    def forward(self, dna_sequence: torch.Tensor, target_length: int) -> torch.Tensor:
        """Return RNA-seq predictions in [B, C, S] format."""

        organism_index = torch.full(
            (dna_sequence.shape[0],),
            self.organism_index,
            dtype=torch.long,
            device=dna_sequence.device,
        )
        dna_sequence_nlc = dna_sequence.transpose(1, 2).contiguous()

        with torch.no_grad():
            embeddings = self.base_model.encode(
                dna_sequence_nlc,
                organism_index,
                resolutions=(self.embedding_resolution,),
                channels_last=False,
            )

        key = f"embeddings_{self.embedding_resolution}bp"
        prediction = self.head(embeddings[key].to(dtype=self.head.weight.dtype))
        if prediction.shape[-1] != target_length:
            prediction = F.interpolate(
                prediction,
                size=target_length,
                mode="linear",
                align_corners=False,
            )
        return prediction


def pick_device(requested: str) -> torch.device:
    if requested != "auto":
        return torch.device(requested)
    if torch.cuda.is_available():
        return torch.device("cuda")
    if hasattr(torch.backends, "mps") and torch.backends.mps.is_available():
        return torch.device("mps")
    return torch.device("cpu")


def count_parameters(module: torch.nn.Module, *, trainable_only: bool = False) -> int:
    parameters = module.parameters()
    if trainable_only:
        parameters = (p for p in parameters if p.requires_grad)
    return sum(p.numel() for p in parameters)


def evaluate_masked_mse(
    model: RnaSeq11Adapter,
    dataloader: torch.utils.data.DataLoader,
    *,
    device: torch.device,
    max_batches: int | None = None,
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

            prediction = model(dna_sequence, target_length=rna_seq.shape[-1])
            loss = masked_mse_loss(prediction, rna_seq, rna_seq_mask)
            batch_size = int(dna_sequence.shape[0])
            weighted_loss_sum += float(loss.detach().cpu()) * batch_size
            n_batches += 1
            n_examples += batch_size
            if max_batches is not None and n_batches >= max_batches:
                break

    model.train()
    model.base_model.eval()
    return weighted_loss_sum / max(n_examples, 1), n_batches, n_examples


def evaluate_masked_mse_by_track(
    model: RnaSeq11Adapter,
    dataloader: torch.utils.data.DataLoader,
    *,
    device: torch.device,
    max_batches: int | None = None,
) -> tuple[float, list[float], int, int]:
    model.eval()
    per_track_sse: torch.Tensor | None = None
    per_track_denominator: torch.Tensor | None = None
    n_batches = 0
    n_examples = 0

    with torch.no_grad():
        for batch in dataloader:
            dna_sequence = batch["dna_sequence"].to(device=device, dtype=torch.float32)
            rna_seq = batch["rna_seq"].to(device=device, dtype=torch.float32)
            rna_seq_mask = batch["rna_seq_mask"].to(device=device)

            prediction = model(dna_sequence, target_length=rna_seq.shape[-1])
            mask = rna_seq_mask
            while mask.ndim < prediction.ndim:
                mask = mask.unsqueeze(-1)
            mask = mask.to(dtype=prediction.dtype, device=prediction.device)

            squared_error = (prediction - rna_seq).square() * mask
            batch_sse = squared_error.sum(dim=(0, 2)).detach().cpu()
            batch_denominator = mask.expand_as(prediction).sum(dim=(0, 2)).detach().cpu()

            if per_track_sse is None:
                per_track_sse = batch_sse
                per_track_denominator = batch_denominator
            else:
                per_track_sse += batch_sse
                per_track_denominator += batch_denominator

            n_batches += 1
            n_examples += int(dna_sequence.shape[0])
            if max_batches is not None and n_batches >= max_batches:
                break

    if per_track_sse is None or per_track_denominator is None:
        raise ValueError("Cannot evaluate an empty dataloader")

    per_track_loss = per_track_sse / per_track_denominator.clamp_min(1.0)
    overall_loss = per_track_sse.sum() / per_track_denominator.sum().clamp_min(1.0)
    model.train()
    model.base_model.eval()
    return (
        float(overall_loss),
        [float(value) for value in per_track_loss],
        n_batches,
        n_examples,
    )


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


def evaluate_pointwise_metrics(
    model: RnaSeq11Adapter,
    dataloader: torch.utils.data.DataLoader,
    *,
    device: torch.device,
    max_batches: int | None = None,
    spearman_sample_size: int = 0,
    spearman_seed: int = 0,
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

            prediction = model(dna_sequence, target_length=rna_seq.shape[-1])
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
                    total_examples = len(dataloader.dataset)
                    if max_batches is not None:
                        total_examples = min(
                            total_examples,
                            max_batches * int(dataloader.batch_size or batch_size),
                        )
                    total_positions = max(total_examples * sequence_length, 1)
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
            if max_batches is not None and n_batches >= max_batches:
                break

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
    model.base_model.eval()
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


def load_adapter_checkpoint(path: str | Path) -> dict[str, Any]:
    checkpoint = torch.load(Path(path), map_location="cpu")
    if "adapter_head_state_dict" not in checkpoint:
        raise ValueError(f"Missing adapter_head_state_dict in checkpoint: {path}")
    return checkpoint
