#!/usr/bin/env python3
"""Shared frozen AlphaGenome adapter utilities for C. elegans RNA-seq."""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Any, Literal, Sequence

import numpy as np
import torch
import torch.nn.functional as F

from alphagenome_pytorch import AlphaGenome
from alphagenome_pytorch.heads import GenomeTracksHead
from alphagenome_pytorch.losses import multinomial_loss as alphagenome_multinomial_loss
from torch_rna_seq_dataset import masked_mse_loss


HeadType = Literal["linear", "genome-tracks"]
LossType = Literal["mse", "poisson-multinomial", "hybrid-mse-poisson"]


def parse_resolutions(
    value: str | Sequence[int] | None,
    *,
    fallback_resolution: int,
) -> tuple[int, ...]:
    """Parse 1 bp / 128 bp head resolutions."""

    if value is None:
        resolutions = (fallback_resolution,)
    elif isinstance(value, str):
        resolutions = tuple(int(part.strip()) for part in value.split(",") if part.strip())
    else:
        resolutions = tuple(int(resolution) for resolution in value)

    if not resolutions:
        raise ValueError("At least one resolution is required")
    invalid = sorted(set(resolutions) - {1, 128})
    if invalid:
        raise ValueError(f"Unsupported resolution(s): {invalid}; use 1 and/or 128")
    return tuple(sorted(set(resolutions)))


def embedding_channels_for_resolution(resolution: int) -> int:
    if resolution == 1:
        return 1536
    if resolution == 128:
        return 3072
    raise ValueError(f"Unsupported resolution: {resolution}")


def normalize_track_means(
    track_means: torch.Tensor | Sequence[float] | None,
    *,
    n_tracks: int,
    num_organisms: int,
) -> torch.Tensor | None:
    """Return track means as [num_organisms, n_tracks] or None."""

    if track_means is None:
        return None

    tensor = torch.as_tensor(track_means, dtype=torch.float32)
    if tensor.ndim == 1:
        if int(tensor.shape[0]) != n_tracks:
            raise ValueError(
                f"track_means has {tensor.shape[0]} tracks; expected {n_tracks}"
            )
        tensor = tensor.unsqueeze(0)
    elif tensor.ndim != 2:
        raise ValueError("track_means must be a 1D or 2D tensor/sequence")

    if int(tensor.shape[1]) != n_tracks:
        raise ValueError(f"track_means has {tensor.shape[1]} tracks; expected {n_tracks}")
    if int(tensor.shape[0]) == 1 and num_organisms > 1:
        tensor = tensor.repeat(num_organisms, 1)
    if int(tensor.shape[0]) < num_organisms:
        raise ValueError(
            f"track_means has {tensor.shape[0]} organisms; expected {num_organisms}"
        )
    return tensor[:num_organisms].contiguous()


def load_track_means_from_grouped_qc(path: str | Path, *, n_tracks: int) -> torch.Tensor:
    """Load per-track mean_signal values from grouped_bigwig_qc.tsv."""

    rows: list[tuple[int, float]] = []
    with Path(path).open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            group_id = row.get("group_id", "")
            if not group_id.startswith("RNA_SEQ_"):
                continue
            try:
                track_index = int(group_id.rsplit("_", maxsplit=1)[1]) - 1
                mean_signal = float(row["mean_signal"])
            except (KeyError, ValueError) as error:
                raise ValueError(f"Cannot parse track mean row in {path}: {row}") from error
            rows.append((track_index, mean_signal))

    rows.sort(key=lambda item: item[0])
    if len(rows) != n_tracks:
        raise ValueError(f"Loaded {len(rows)} track means from {path}; expected {n_tracks}")
    expected_indices = list(range(n_tracks))
    observed_indices = [index for index, _mean in rows]
    if observed_indices != expected_indices:
        raise ValueError(
            f"Track mean indices in {path} are {observed_indices}; "
            f"expected {expected_indices}"
        )
    return torch.tensor([mean for _index, mean in rows], dtype=torch.float32)


def compute_nonzero_track_means_from_dataset(dataset) -> torch.Tensor:
    """Compute AlphaGenome-style nonzero means from a raw-target NPZ dataset."""

    n_tracks = int(dataset.metadata["n_tracks"])
    sums = torch.zeros(n_tracks, dtype=torch.float64)
    counts = torch.zeros(n_tracks, dtype=torch.float64)
    for index in range(len(dataset)):
        rna_seq = dataset[index]["rna_seq"].to(dtype=torch.float32)
        if rna_seq.ndim != 2:
            raise ValueError(f"Expected rna_seq [tracks, sequence], got {rna_seq.shape}")
        rna_seq = rna_seq.clamp_min(0.0)
        nonzero = rna_seq > 0
        sums += (rna_seq * nonzero).sum(dim=1, dtype=torch.float64)
        counts += nonzero.sum(dim=1).to(dtype=torch.float64)

    means = torch.where(counts > 0, sums / counts.clamp_min(1.0), torch.ones_like(sums))
    return means.to(dtype=torch.float32)


def bin_rna_seq_target(target: torch.Tensor, resolution: int) -> torch.Tensor:
    """Bin raw RNA-seq targets by summing bases into native head resolution."""

    target = target.clamp_min(0.0)
    if resolution == 1:
        return target

    usable_length = (target.shape[-1] // resolution) * resolution
    if usable_length == 0:
        raise ValueError(
            f"Target length {target.shape[-1]} is too short for resolution {resolution}"
        )
    target = target[..., :usable_length]
    return target.reshape(*target.shape[:-1], usable_length // resolution, resolution).sum(
        dim=-1
    )


def format_prediction_shape(prediction: torch.Tensor | dict[int, torch.Tensor]) -> str:
    if isinstance(prediction, dict):
        return ",".join(
            f"{resolution}:{'x'.join(map(str, tensor.shape))}"
            for resolution, tensor in sorted(prediction.items())
        )
    return "x".join(map(str, prediction.shape))


class RnaSeq11Adapter(torch.nn.Module):
    """Project frozen AlphaGenome sequence embeddings to custom RNA-seq tracks."""

    def __init__(
        self,
        base_model: AlphaGenome,
        *,
        n_tracks: int,
        embedding_resolution: int = 128,
        organism_index: int = 0,
        encode_requires_grad: bool = False,
        head_type: HeadType = "linear",
        head_resolutions: Sequence[int] | None = None,
        track_means: torch.Tensor | Sequence[float] | None = None,
        num_head_organisms: int | None = None,
    ) -> None:
        super().__init__()
        if embedding_resolution not in {1, 128}:
            raise ValueError("embedding_resolution must be 1 or 128")
        if head_type not in {"linear", "genome-tracks"}:
            raise ValueError("head_type must be 'linear' or 'genome-tracks'")

        self.base_model = base_model
        self.embedding_resolution = embedding_resolution
        self.organism_index = organism_index
        self.encode_requires_grad = encode_requires_grad
        self.head_type = head_type
        self.head_resolutions = parse_resolutions(
            head_resolutions,
            fallback_resolution=embedding_resolution,
        )
        self.n_tracks = n_tracks

        if self.head_type == "linear":
            if len(self.head_resolutions) != 1:
                raise ValueError("linear head supports exactly one embedding resolution")
            if self.head_resolutions[0] != embedding_resolution:
                raise ValueError(
                    "linear head_resolutions must match embedding_resolution"
                )
            in_channels = embedding_channels_for_resolution(embedding_resolution)
            self.head = torch.nn.Conv1d(in_channels, n_tracks, kernel_size=1)
        else:
            num_organisms = (
                max(organism_index + 1, 1)
                if num_head_organisms is None
                else num_head_organisms
            )
            normalized_means = normalize_track_means(
                track_means,
                n_tracks=n_tracks,
                num_organisms=num_organisms,
            )
            self.head = GenomeTracksHead(
                in_channels=None,
                num_tracks=n_tracks,
                resolutions=self.head_resolutions,
                num_organisms=num_organisms,
                apply_squashing=True,
                track_means=normalized_means,
            )

        for parameter in self.base_model.parameters():
            parameter.requires_grad = False
        self.base_model.eval()

    @property
    def uses_scaled_targets(self) -> bool:
        return self.head_type == "genome-tracks"

    def organism_index_tensor(self, batch_size: int, device: torch.device) -> torch.Tensor:
        return torch.full(
            (batch_size,),
            self.organism_index,
            dtype=torch.long,
            device=device,
        )

    def scale_targets_for_loss(
        self,
        target: torch.Tensor,
        *,
        resolution: int,
    ) -> torch.Tensor:
        if self.head_type != "genome-tracks":
            return target

        organism_index = self.organism_index_tensor(target.shape[0], target.device)
        target_at_resolution = bin_rna_seq_target(target, resolution)
        return self.head.scale(
            target_at_resolution,
            organism_index,
            resolution,
            channels_last=False,
        )

    def forward(
        self,
        dna_sequence: torch.Tensor,
        target_length: int | None = None,
        *,
        return_scaled: bool = False,
    ) -> torch.Tensor | dict[int, torch.Tensor]:
        """Return RNA-seq predictions in [B, C, S] format."""

        organism_index = self.organism_index_tensor(
            dna_sequence.shape[0],
            dna_sequence.device,
        )
        dna_sequence_nlc = dna_sequence.transpose(1, 2).contiguous()

        grad_context = (
            torch.enable_grad()
            if self.encode_requires_grad and torch.is_grad_enabled()
            else torch.no_grad()
        )
        with grad_context:
            embeddings = self.base_model.encode(
                dna_sequence_nlc,
                organism_index,
                resolutions=self.head_resolutions,
                channels_last=False,
            )

        if self.head_type == "genome-tracks":
            head_dtype = next(self.head.parameters()).dtype
            embeddings_by_resolution = {
                resolution: embeddings[f"embeddings_{resolution}bp"].to(dtype=head_dtype)
                for resolution in self.head_resolutions
            }
            return self.head(
                embeddings_by_resolution,
                organism_index,
                return_scaled=return_scaled,
                channels_last=False,
            )

        key = f"embeddings_{self.embedding_resolution}bp"
        prediction = self.head(embeddings[key].to(dtype=self.head.weight.dtype))
        if target_length is not None and prediction.shape[-1] != target_length:
            prediction = F.interpolate(
                prediction,
                size=target_length,
                mode="linear",
                align_corners=False,
            )
        return prediction


def prediction_return_scaled_for_loss(model: RnaSeq11Adapter, loss_type: LossType) -> bool:
    """Return whether the head should emit model-space predictions for loss."""

    return model.head_type == "genome-tracks" and loss_type in {
        "mse",
        "hybrid-mse-poisson",
    }


def _mask_for_prediction(mask: torch.Tensor, prediction: torch.Tensor) -> torch.Tensor:
    while mask.ndim < prediction.ndim:
        mask = mask.unsqueeze(-1)
    return mask.to(device=prediction.device, dtype=torch.bool)


def _genome_tracks_mse_component(
    model: RnaSeq11Adapter,
    prediction: dict[int, torch.Tensor],
    target: torch.Tensor,
    mask: torch.Tensor,
) -> torch.Tensor:
    losses = []
    for resolution, resolution_prediction in sorted(prediction.items()):
        scaled_target = model.scale_targets_for_loss(target, resolution=resolution)
        losses.append(
            masked_mse_loss(
                resolution_prediction,
                scaled_target.to(dtype=resolution_prediction.dtype),
                mask,
            )
        )
    if not losses:
        raise ValueError("No prediction resolutions were returned")
    return torch.stack(losses).mean()


def _poisson_multinomial_component(
    model: RnaSeq11Adapter,
    prediction: dict[int, torch.Tensor],
    target: torch.Tensor,
    mask: torch.Tensor,
    *,
    prediction_is_scaled: bool,
    multinomial_num_segments: int,
    positional_weight: float,
    count_weight: float,
) -> torch.Tensor:
    if multinomial_num_segments < 1:
        raise ValueError("multinomial_num_segments must be >= 1")

    organism_index = model.organism_index_tensor(target.shape[0], target.device)
    losses = []
    for resolution, resolution_prediction in sorted(prediction.items()):
        target_at_resolution = bin_rna_seq_target(target, resolution)
        raw_prediction = resolution_prediction
        if prediction_is_scaled:
            raw_prediction = model.head.unscale(
                resolution_prediction,
                organism_index,
                resolution,
                channels_last=False,
            )
        if target_at_resolution.shape != raw_prediction.shape:
            raise ValueError(
                "Prediction/target shape mismatch for "
                f"{resolution} bp: {raw_prediction.shape} vs "
                f"{target_at_resolution.shape}"
            )
        sequence_length = int(raw_prediction.shape[-1])
        if sequence_length % multinomial_num_segments != 0:
            raise ValueError(
                f"{sequence_length=} is not divisible by "
                f"{multinomial_num_segments=}"
            )
        segment_length = sequence_length // multinomial_num_segments
        loss_dict = alphagenome_multinomial_loss(
            y_true=target_at_resolution.to(dtype=raw_prediction.dtype),
            y_pred=raw_prediction.clamp_min(1e-7),
            mask=_mask_for_prediction(mask, raw_prediction),
            multinomial_resolution=segment_length,
            positional_weight=positional_weight,
            count_weight=count_weight,
            channels_last=False,
        )
        losses.append(loss_dict["loss"])
    if not losses:
        raise ValueError("No prediction resolutions were returned")
    return torch.stack(losses).mean()


def masked_adapter_loss(
    model: RnaSeq11Adapter,
    prediction: torch.Tensor | dict[int, torch.Tensor],
    target: torch.Tensor,
    mask: torch.Tensor,
    *,
    loss_type: LossType = "mse",
    multinomial_num_segments: int = 8,
    positional_weight: float = 5.0,
    count_weight: float = 1.0,
    mse_weight: float = 1.0,
    poisson_weight: float = 1.0,
) -> torch.Tensor:
    """Compute legacy MSE or AlphaGenome-style RNA-seq count/position loss."""

    if loss_type == "poisson-multinomial":
        if model.head_type != "genome-tracks" or not isinstance(prediction, dict):
            raise ValueError(
                "poisson-multinomial loss requires a GenomeTracksHead prediction dict"
            )
        return _poisson_multinomial_component(
            model,
            prediction,
            target,
            mask,
            prediction_is_scaled=False,
            multinomial_num_segments=multinomial_num_segments,
            positional_weight=positional_weight,
            count_weight=count_weight,
        )

    if loss_type == "hybrid-mse-poisson":
        if model.head_type != "genome-tracks" or not isinstance(prediction, dict):
            raise ValueError(
                "hybrid-mse-poisson loss requires a GenomeTracksHead prediction dict"
            )
        mse_loss = _genome_tracks_mse_component(model, prediction, target, mask)
        poisson_loss = _poisson_multinomial_component(
            model,
            prediction,
            target,
            mask,
            prediction_is_scaled=True,
            multinomial_num_segments=multinomial_num_segments,
            positional_weight=positional_weight,
            count_weight=count_weight,
        )
        return mse_weight * mse_loss + poisson_weight * poisson_loss

    if loss_type != "mse":
        raise ValueError(f"Unsupported loss_type: {loss_type}")

    if isinstance(prediction, dict):
        return _genome_tracks_mse_component(model, prediction, target, mask)

    return masked_mse_loss(prediction, target, mask)


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
    loss_type: LossType = "mse",
    multinomial_num_segments: int = 8,
    positional_weight: float = 5.0,
    count_weight: float = 1.0,
    mse_weight: float = 1.0,
    poisson_weight: float = 1.0,
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

            prediction = model(
                dna_sequence,
                target_length=rna_seq.shape[-1],
                return_scaled=prediction_return_scaled_for_loss(model, loss_type),
            )
            loss = masked_adapter_loss(
                model,
                prediction,
                rna_seq,
                rna_seq_mask,
                loss_type=loss_type,
                multinomial_num_segments=multinomial_num_segments,
                positional_weight=positional_weight,
                count_weight=count_weight,
                mse_weight=mse_weight,
                poisson_weight=poisson_weight,
            )
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

            prediction = model(
                dna_sequence,
                target_length=rna_seq.shape[-1],
                return_scaled=model.uses_scaled_targets,
            )
            if isinstance(prediction, dict):
                prediction_items = sorted(prediction.items())
            else:
                prediction_items = [(model.embedding_resolution, prediction)]
            mask = rna_seq_mask
            while mask.ndim < 3:
                mask = mask.unsqueeze(-1)

            batch_sse = None
            batch_denominator = None
            for resolution, resolution_prediction in prediction_items:
                if isinstance(prediction, dict):
                    target_for_loss = model.scale_targets_for_loss(
                        rna_seq,
                        resolution=resolution,
                    )
                else:
                    target_for_loss = rna_seq
                resolution_mask = mask.to(
                    dtype=resolution_prediction.dtype,
                    device=resolution_prediction.device,
                )
                squared_error = (
                    resolution_prediction - target_for_loss
                ).square() * resolution_mask
                resolution_sse = squared_error.sum(dim=(0, 2)).detach().cpu()
                resolution_denominator = (
                    resolution_mask.expand_as(resolution_prediction)
                    .sum(dim=(0, 2))
                    .detach()
                    .cpu()
                )
                if batch_sse is None:
                    batch_sse = resolution_sse
                    batch_denominator = resolution_denominator
                else:
                    batch_sse += resolution_sse
                    batch_denominator += resolution_denominator

            assert batch_sse is not None
            assert batch_denominator is not None

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
    if model.head_type != "linear":
        raise NotImplementedError(
            "Pointwise metrics are currently implemented for the legacy linear "
            "head only. Use evaluate_masked_mse for GenomeTracksHead model-space "
            "sanity checks."
        )

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


def evaluate_128bp_binned_pointwise_metrics(
    model: RnaSeq11Adapter,
    dataloader: torch.utils.data.DataLoader,
    *,
    device: torch.device,
    linear_prediction_transform: str,
    metric_transform: Literal["raw-sum", "log1p-mean"] = "log1p-mean",
    max_batches: int | None = None,
) -> dict[str, object]:
    """Evaluate any RNA-seq head on common 128 bp binned raw-target metrics."""

    if metric_transform not in {"raw-sum", "log1p-mean"}:
        raise ValueError(f"Unsupported metric_transform: {metric_transform}")

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
    n_batches = 0
    n_examples = 0

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
            if isinstance(prediction, dict):
                if 128 in prediction:
                    prediction_128 = prediction[128].clamp_min(0.0)
                elif 1 in prediction:
                    prediction_128 = bin_rna_seq_target(prediction[1], 128)
                else:
                    raise ValueError("GenomeTracksHead returned neither 1 bp nor 128 bp")
            else:
                prediction_raw = prediction
                if linear_prediction_transform == "log1p":
                    prediction_raw = torch.expm1(prediction_raw).clamp_min(0.0)
                elif linear_prediction_transform == "none":
                    prediction_raw = prediction_raw.clamp_min(0.0)
                else:
                    raise ValueError(
                        "linear_prediction_transform must be 'log1p' or 'none'"
                    )
                prediction_128 = bin_rna_seq_target(prediction_raw, 128)

            target_128 = bin_rna_seq_target(rna_seq, 128)
            if metric_transform == "log1p-mean":
                prediction_metric = torch.log1p(prediction_128 / 128.0)
                target_metric = torch.log1p(target_128 / 128.0)
            else:
                prediction_metric = prediction_128
                target_metric = target_128

            if n_tracks is None:
                n_tracks = int(prediction_metric.shape[1])
                sse = torch.zeros(n_tracks, dtype=torch.float64)
                sae = torch.zeros(n_tracks, dtype=torch.float64)
                sum_x = torch.zeros(n_tracks, dtype=torch.float64)
                sum_y = torch.zeros(n_tracks, dtype=torch.float64)
                sum_x2 = torch.zeros(n_tracks, dtype=torch.float64)
                sum_y2 = torch.zeros(n_tracks, dtype=torch.float64)
                sum_xy = torch.zeros(n_tracks, dtype=torch.float64)
                denominator = torch.zeros(n_tracks, dtype=torch.float64)

            mask = _mask_for_prediction(rna_seq_mask, prediction_metric)
            mask64 = mask.to(dtype=torch.float64)
            prediction64 = prediction_metric.to(dtype=torch.float64)
            target64 = target_metric.to(dtype=torch.float64)
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

            n_batches += 1
            n_examples += int(dna_sequence.shape[0])
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
        "per_track_spearman_sampled": None,
        "per_track_spearman_sampled_n": None,
        "batches": n_batches,
        "examples": n_examples,
    }


def load_adapter_checkpoint(path: str | Path) -> dict[str, Any]:
    checkpoint = torch.load(Path(path), map_location="cpu")
    if "adapter_head_state_dict" not in checkpoint:
        raise ValueError(f"Missing adapter_head_state_dict in checkpoint: {path}")
    return checkpoint
