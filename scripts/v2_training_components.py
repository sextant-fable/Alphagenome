#!/usr/bin/env python3
"""Audited v2 loss, augmentation, and A/B/C RNA-seq model components."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Mapping, Sequence

import torch
import torch.nn.functional as F

from alphagenome_pytorch.heads import (
    GenomeTracksHead,
    predictions_scaling,
    targets_scaling,
)
from alphagenome_pytorch.losses import multinomial_loss


RESOLUTIONS = (1, 128)
NUM_SEGMENTS = 8
POSITIONAL_WEIGHT = 5.0
GENE_CROSS_TRACK_WEIGHT = 0.1


def normalize_track_means(
    track_means: torch.Tensor | Sequence[float],
    batch_size: int,
    device: torch.device,
) -> torch.Tensor:
    means = torch.as_tensor(track_means, dtype=torch.float32, device=device)
    if means.ndim == 1:
        means = means.unsqueeze(0)
    if means.ndim != 2:
        raise ValueError("track_means must have shape [C] or [B,C]")
    if means.shape[0] == 1 and batch_size > 1:
        means = means.expand(batch_size, -1)
    if means.shape[0] != batch_size or torch.any(means <= 0):
        raise ValueError("track_means batch mismatch or non-positive value")
    return means


def scale_targets_model_space(
    target: torch.Tensor,
    track_means: torch.Tensor | Sequence[float],
    resolution: int,
) -> torch.Tensor:
    if resolution not in RESOLUTIONS:
        raise ValueError(f"Unsupported resolution: {resolution}")
    means = normalize_track_means(track_means, target.shape[0], target.device)
    return targets_scaling(
        target,
        track_means=means,
        resolution=resolution,
        apply_squashing=True,
        channels_last=False,
    )


def unscale_predictions_experimental_space(
    prediction: torch.Tensor,
    track_means: torch.Tensor | Sequence[float],
    resolution: int,
) -> torch.Tensor:
    means = normalize_track_means(
        track_means, prediction.shape[0], prediction.device
    )
    return predictions_scaling(
        prediction,
        track_means=means,
        resolution=resolution,
        apply_squashing=True,
        channels_last=False,
    )


def gene_cross_track_loss(
    prediction: torch.Tensor,
    target: torch.Tensor,
    track_mask: torch.Tensor,
    track_strand: torch.Tensor,
    gene_mask: torch.Tensor,
) -> torch.Tensor:
    """KL across tracks after strand-aware gene-body aggregation."""

    if prediction.shape != target.shape or prediction.ndim != 3:
        raise ValueError("prediction and target must have shape [B,C,L]")
    if gene_mask.shape != (prediction.shape[0], 2, prediction.shape[-1]):
        raise ValueError("gene_mask must have shape [B,2,L]")
    if track_strand.shape != prediction.shape[:2]:
        raise ValueError("track_strand must have shape [B,C]")
    valid_tracks = track_mask.reshape(prediction.shape[0], prediction.shape[1]).bool()
    losses = []
    for batch_index in range(prediction.shape[0]):
        for gene_strand_index, gene_strand in enumerate((1, -1)):
            positions = gene_mask[batch_index, gene_strand_index].bool()
            compatible = valid_tracks[batch_index] & (
                (track_strand[batch_index] == 0)
                | (track_strand[batch_index] == gene_strand)
            )
            if not positions.any() or compatible.sum() < 2:
                continue
            pred_counts = prediction[batch_index, compatible][:, positions].sum(dim=1)
            true_counts = target[batch_index, compatible][:, positions].sum(dim=1)
            total_true = true_counts.sum()
            if total_true <= 0:
                continue
            true_probability = true_counts / total_true
            pred_probability = pred_counts.clamp_min(1e-7)
            pred_probability = pred_probability / pred_probability.sum()
            kl = (
                true_probability
                * (
                    torch.log(true_probability.clamp_min(1e-7))
                    - torch.log(pred_probability.clamp_min(1e-7))
                )
            ).sum()
            losses.append(kl)
    if not losses:
        return prediction.sum() * 0.0
    return torch.stack(losses).mean()


def dual_resolution_paper_loss(
    predictions_model_space: Mapping[int, torch.Tensor],
    targets_experimental_space: Mapping[int, torch.Tensor],
    *,
    track_means: torch.Tensor | Sequence[float],
    track_mask: torch.Tensor,
    track_strand: torch.Tensor,
    gene_mask: torch.Tensor,
    resolution_weights: Mapping[int, float] | None = None,
    positional_weight: float = POSITIONAL_WEIGHT,
    gene_weight: float = GENE_CROSS_TRACK_WEIGHT,
) -> tuple[torch.Tensor, dict[str, torch.Tensor]]:
    weights = dict(resolution_weights or {1: 1.0, 128: 1.0})
    metrics: dict[str, torch.Tensor] = {}
    total = next(iter(predictions_model_space.values())).sum() * 0.0
    scaled_targets = {}
    for resolution in RESOLUTIONS:
        if resolution not in predictions_model_space or resolution not in targets_experimental_space:
            raise ValueError(f"Missing {resolution} bp prediction or target")
        prediction = predictions_model_space[resolution]
        target = targets_experimental_space[resolution]
        if prediction.shape != target.shape:
            raise ValueError(
                f"Resolution {resolution} shape mismatch: {prediction.shape} vs {target.shape}"
            )
        if prediction.shape[-1] % NUM_SEGMENTS:
            raise ValueError(
                f"Resolution {resolution} length must divide into {NUM_SEGMENTS} segments"
            )
        if not torch.isfinite(prediction).all() or torch.any(prediction < 0):
            raise ValueError("Model-space predictions must be finite and non-negative")
        scaled_target = scale_targets_model_space(target, track_means, resolution)
        scaled_targets[resolution] = scaled_target
        loss_values = multinomial_loss(
            y_true=scaled_target,
            y_pred=prediction,
            mask=track_mask.bool(),
            multinomial_resolution=prediction.shape[-1] // NUM_SEGMENTS,
            positional_weight=positional_weight,
            channels_last=False,
        )
        weight = float(weights.get(resolution, 0.0))
        total = total + weight * loss_values["loss"]
        metrics[f"loss_{resolution}bp"] = loss_values["loss"]
        metrics[f"count_{resolution}bp"] = loss_values["loss_total"]
        metrics[f"position_{resolution}bp"] = loss_values["loss_positional"]
    gene_loss = gene_cross_track_loss(
        predictions_model_space[1],
        scaled_targets[1],
        track_mask,
        track_strand,
        gene_mask,
    )
    total = total + gene_weight * gene_loss
    metrics["gene_cross_track"] = gene_loss
    metrics["loss"] = total
    return total, metrics


def log1p_mse_loss(
    predictions_model_space: Mapping[int, torch.Tensor],
    targets_experimental_space: Mapping[int, torch.Tensor],
    *,
    track_means: torch.Tensor | Sequence[float],
    track_mask: torch.Tensor,
    resolution_weights: Mapping[int, float] | None = None,
) -> tuple[torch.Tensor, dict[str, torch.Tensor]]:
    weights = dict(resolution_weights or {1: 1.0, 128: 1.0})
    total = next(iter(predictions_model_space.values())).sum() * 0.0
    metrics = {}
    mask = track_mask.bool()
    for resolution in RESOLUTIONS:
        prediction = unscale_predictions_experimental_space(
            predictions_model_space[resolution], track_means, resolution
        )
        target = targets_experimental_space[resolution]
        squared = (torch.log1p(prediction) - torch.log1p(target.clamp_min(0))).square()
        expanded_mask = mask.expand_as(squared)
        value = squared[expanded_mask].mean() if expanded_mask.any() else squared.sum() * 0
        total = total + float(weights.get(resolution, 0.0)) * value
        metrics[f"loss_{resolution}bp"] = value
    metrics["loss"] = total
    return total, metrics


def reverse_complement_item(
    item: Mapping[str, Any],
    strand_pair_index: torch.Tensor | Sequence[int] | None = None,
) -> dict[str, Any]:
    result = {
        key: value.clone() if isinstance(value, torch.Tensor) else value
        for key, value in item.items()
    }
    n_tracks = int(result["target_1bp"].shape[-2])
    if strand_pair_index is None:
        pair_index = torch.arange(n_tracks, dtype=torch.long)
    else:
        pair_index = torch.as_tensor(strand_pair_index, dtype=torch.long)
    if pair_index.shape != (n_tracks,) or not torch.equal(
        pair_index[pair_index], torch.arange(n_tracks)
    ):
        raise ValueError("strand_pair_index must be an involution over tracks")
    result["dna_sequence"] = result["dna_sequence"][[3, 2, 1, 0]].flip(-1)
    result["target_1bp"] = result["target_1bp"][pair_index].flip(-1)
    result["target_128bp"] = result["target_128bp"][pair_index].flip(-1)
    result["track_mask"] = result["track_mask"][pair_index]
    result["track_strand"] = -result["track_strand"][pair_index]
    result["gene_mask"] = result["gene_mask"][[1, 0]].flip(-1)
    result["core_mask"] = result["core_mask"].flip(-1)
    result["reverse_complemented"] = not bool(item.get("reverse_complemented", False))
    return result


def synchronized_crop(
    item: Mapping[str, Any], *, start: int, length: int
) -> dict[str, Any]:
    """Crop a flank-extended item while keeping all coordinate tensors aligned."""

    if start < 0 or length <= 0:
        raise ValueError("Invalid crop")
    end = start + length
    if end > item["dna_sequence"].shape[-1] or length % 128:
        raise ValueError("Crop escapes item or is not divisible by 128")
    result = {
        key: value.clone() if isinstance(value, torch.Tensor) else value
        for key, value in item.items()
    }
    result["dna_sequence"] = item["dna_sequence"][..., start:end]
    result["target_1bp"] = item["target_1bp"][..., start:end]
    result["gene_mask"] = item["gene_mask"][..., start:end]
    result["core_mask"] = item["core_mask"][..., start:end]
    result["target_128bp"] = result["target_1bp"].reshape(
        result["target_1bp"].shape[-2], length // 128, 128
    ).sum(dim=-1)
    result["shift_bp"] = start
    return result


def _expand_embedding(embedding: torch.nn.Embedding) -> torch.nn.Embedding:
    if embedding.num_embeddings != 2:
        raise ValueError("Expected a two-organism embedding")
    expanded = torch.nn.Embedding(
        3,
        embedding.embedding_dim,
        device=embedding.weight.device,
        dtype=embedding.weight.dtype,
    )
    with torch.no_grad():
        expanded.weight[:2].copy_(embedding.weight)
        expanded.weight[2].copy_(embedding.weight.mean(dim=0))
    gradient_mask = torch.zeros_like(expanded.weight)
    gradient_mask[2] = 1
    expanded.weight.register_hook(lambda gradient: gradient * gradient_mask)
    expanded.weight.requires_grad_(True)
    return expanded


def add_c_elegans_organism_embeddings(model: torch.nn.Module) -> list[str]:
    """Add actual third-row organism embeddings used by AlphaGenome encode."""

    paths = (
        "organism_embed",
        "embedder_128bp.organism_embed",
        "embedder_1bp.organism_embed",
        "embedder_pair.organism_embed",
    )
    expanded_paths = []
    for path in paths:
        parent = model
        parts = path.split(".")
        for part in parts[:-1]:
            parent = getattr(parent, part)
        name = parts[-1]
        setattr(parent, name, _expand_embedding(getattr(parent, name)))
        if hasattr(parent, "num_organisms"):
            parent.num_organisms = 3
        expanded_paths.append(path)
    model.num_organisms = 3
    return expanded_paths


class DualResolutionRnaHead(torch.nn.Module):
    def __init__(self, n_tracks: int, track_means: torch.Tensor) -> None:
        super().__init__()
        means = torch.as_tensor(track_means, dtype=torch.float32).reshape(1, n_tracks)
        self.head = GenomeTracksHead(
            in_channels={1: 1536, 128: 3072},
            num_tracks=n_tracks,
            resolutions=RESOLUTIONS,
            num_organisms=1,
            track_means=means,
            apply_squashing=True,
        )

    def forward(self, embeddings: Mapping[int, torch.Tensor]) -> dict[int, torch.Tensor]:
        batch_size = next(iter(embeddings.values())).shape[0]
        organism = torch.zeros(
            batch_size,
            dtype=torch.long,
            device=next(iter(embeddings.values())).device,
        )
        return self.head(
            dict(embeddings), organism, return_scaled=True, channels_last=False
        )


class AlphaGenomeRnaModel(torch.nn.Module):
    """Model A/B wrapper with separate trunk and one-organism worm RNA head."""

    def __init__(
        self,
        base_model: torch.nn.Module,
        n_tracks: int,
        track_means: torch.Tensor,
        *,
        base_organism_index: int,
        encode_requires_grad: bool,
    ) -> None:
        super().__init__()
        self.base_model = base_model
        self.rna_head = DualResolutionRnaHead(n_tracks, track_means)
        self.base_organism_index = int(base_organism_index)
        self.encode_requires_grad = bool(encode_requires_grad)

    def train(self, mode: bool = True) -> "AlphaGenomeRnaModel":
        super().train(mode)
        if not self.encode_requires_grad:
            self.base_model.eval()
        return self

    def forward(self, dna: torch.Tensor) -> dict[int, torch.Tensor]:
        organism = torch.full(
            (dna.shape[0],),
            self.base_organism_index,
            dtype=torch.long,
            device=dna.device,
        )
        context = torch.enable_grad() if self.encode_requires_grad else torch.no_grad()
        with context:
            encoded = self.base_model.encode(
                dna.transpose(1, 2).contiguous(),
                organism,
                resolutions=RESOLUTIONS,
                channels_last=False,
            )
        embeddings = {
            1: encoded["embeddings_1bp"],
            128: encoded["embeddings_128bp"],
        }
        return self.rna_head(embeddings)


class ResidualBlock(torch.nn.Module):
    def __init__(self, channels: int, dilation: int) -> None:
        super().__init__()
        self.norm = torch.nn.GroupNorm(8, channels)
        self.conv = torch.nn.Conv1d(
            channels, channels, kernel_size=7, padding=3 * dilation, dilation=dilation
        )
        self.pointwise = torch.nn.Conv1d(channels, channels, kernel_size=1)

    def forward(self, value: torch.Tensor) -> torch.Tensor:
        update = self.pointwise(F.gelu(self.conv(F.gelu(self.norm(value)))))
        return value + update


class WormSequenceBaseline(torch.nn.Module):
    """Model C: a residual sequence baseline with matched dual outputs."""

    def __init__(self, n_tracks: int, hidden_channels: int = 64) -> None:
        super().__init__()
        if hidden_channels % 8:
            raise ValueError("hidden_channels must be divisible by 8")
        self.stem = torch.nn.Conv1d(4, hidden_channels, kernel_size=15, padding=7)
        self.blocks = torch.nn.Sequential(
            *(ResidualBlock(hidden_channels, dilation) for dilation in (1, 2, 4, 8, 16, 32))
        )
        self.head_1bp = torch.nn.Conv1d(hidden_channels, n_tracks, kernel_size=1)
        self.head_128bp = torch.nn.Sequential(
            torch.nn.Conv1d(hidden_channels, hidden_channels, kernel_size=3, padding=1),
            torch.nn.GELU(),
            torch.nn.Conv1d(hidden_channels, n_tracks, kernel_size=1),
        )

    def forward(self, dna: torch.Tensor) -> dict[int, torch.Tensor]:
        hidden = self.blocks(self.stem(dna))
        pooled = F.avg_pool1d(hidden, kernel_size=128, stride=128)
        return {
            1: F.softplus(self.head_1bp(hidden)),
            128: F.softplus(self.head_128bp(pooled)),
        }


@dataclass(frozen=True)
class ModelSpec:
    model_id: str
    description: str
    base_organism_index: int | None
    trunk_policy: str


MODEL_SPECS = (
    ModelSpec("A", "frozen AlphaGenome trunk with new dual RNA head", 0, "frozen"),
    ModelSpec(
        "B",
        "new C. elegans embedding plus final-tower LoRA and dual RNA head",
        2,
        "worm_embeddings_and_lora_only",
    ),
    ModelSpec("C", "residual sequence baseline trained from scratch", None, "from_scratch"),
)


class AugmentedV2Dataset(torch.utils.data.Dataset):
    """Deterministic per-epoch +/- shift and 50% reverse-complement wrapper."""

    def __init__(
        self,
        dataset: torch.utils.data.Dataset,
        *,
        max_shift_bp: int = 1024,
        reverse_complement_probability: float = 0.5,
        seed: int = 0,
        strand_pair_index: Sequence[int] | None = None,
        sequence_length: int | None = None,
    ) -> None:
        if max_shift_bp < 0:
            raise ValueError("max_shift_bp must be non-negative")
        if not 0 <= reverse_complement_probability <= 1:
            raise ValueError("reverse_complement_probability must be in [0,1]")
        if not hasattr(dataset, "get_subwindow_item"):
            raise TypeError("dataset must implement get_subwindow_item")
        self.dataset = dataset
        self.max_shift_bp = int(max_shift_bp)
        self.reverse_complement_probability = float(reverse_complement_probability)
        self.seed = int(seed)
        self.epoch = 0
        self.strand_pair_index = strand_pair_index
        self.sequence_length = sequence_length

    def set_epoch(self, epoch: int) -> None:
        self.epoch = int(epoch)

    def __len__(self) -> int:
        return len(self.dataset)

    def __getitem__(self, index: int) -> dict[str, Any]:
        row = self.dataset.intervals[index]
        chromosome_length = self.dataset.fai[row["chromosome"]][0]
        block_start = int(row.get("block_start", 0))
        block_end = int(row.get("block_end", chromosome_length))
        minimum = max(
            -self.max_shift_bp,
            block_start - int(row["start"]),
        )
        maximum = min(
            self.max_shift_bp,
            block_end - int(row["end"]),
        )
        generator = torch.Generator()
        generator.manual_seed(
            (self.seed + 1_000_003 * self.epoch + 9_176 * int(index)) % (2**63 - 1)
        )
        shift = int(
            torch.randint(minimum, maximum + 1, (1,), generator=generator).item()
        )
        full_width = int(row["end"]) - int(row["start"])
        sequence_length = full_width if self.sequence_length is None else int(self.sequence_length)
        if sequence_length <= 0 or sequence_length > full_width or sequence_length % 128:
            raise ValueError("sequence_length must divide into 128 bp bins within the window")
        crop_choices = (full_width - sequence_length) // 128 + 1
        crop_offset = 128 * int(
            torch.randint(0, crop_choices, (1,), generator=generator).item()
        )
        item = self.dataset.get_subwindow_item(
            index,
            shift_bp=shift,
            crop_offset_bp=crop_offset,
            crop_length_bp=sequence_length,
        )
        apply_reverse = bool(
            torch.rand((), generator=generator).item()
            < self.reverse_complement_probability
        )
        if apply_reverse:
            item = reverse_complement_item(item, self.strand_pair_index)
        item["shift_bp"] = shift
        item["crop_offset_bp"] = crop_offset
        return item
