#!/usr/bin/env python3
"""Fine-tune an 11-track C. elegans RNA-seq adapter on frozen AlphaGenome."""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import random
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import torch
import torch.distributed as dist
import torch.nn.functional as F
from torch.nn.parallel import DistributedDataParallel
from torch.utils.data.distributed import DistributedSampler

from alphagenome_pytorch import AlphaGenome
from alphagenome_pytorch.extensions.finetuning.adapters import (
    LoRA,
    apply_lora,
    get_adapter_params,
)
from alphagenome_rna_seq11_adapter import (
    RnaSeq11Adapter,
    bin_rna_seq_target,
    compute_linear_track_means_from_dataset,
    compute_nonzero_track_means_from_dataset,
    count_parameters,
    format_prediction_shape,
    is_residual_linear_architecture,
    load_track_means_from_grouped_qc,
    load_adapter_checkpoint,
    masked_adapter_loss,
    parse_resolutions,
    pick_device,
    prepare_linear_target,
    prediction_return_scaled_for_loss,
)
from torch_rna_seq_dataset import AlphaGenomeRnaSeqNpzDataset, make_dataloader


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
    parser.add_argument(
        "--init-adapter-checkpoint",
        default=None,
        help="Optional adapter checkpoint used to initialize the 11-track head.",
    )
    parser.add_argument(
        "--init-adapter-checkpoint-mode",
        choices=["strict", "matching-shapes"],
        default="strict",
        help=(
            "How to load --init-adapter-checkpoint. strict requires an exact "
            "adapter architecture match; matching-shapes copies only adapter "
            "state_dict entries with identical names and shapes."
        ),
    )
    parser.add_argument("--batch-size", type=int, default=1)
    parser.add_argument("--num-workers", type=int, default=0)
    parser.add_argument("--max-steps", type=int, default=20)
    parser.add_argument("--eval-every", type=int, default=5)
    parser.add_argument("--max-train-examples", type=int, default=None)
    parser.add_argument("--max-valid-examples", type=int, default=None)
    parser.add_argument("--learning-rate", type=float, default=1e-3)
    parser.add_argument("--weight-decay", type=float, default=0.0)
    parser.add_argument(
        "--grad-accum-steps",
        type=int,
        default=1,
        help="Number of microbatches to accumulate per optimizer step.",
    )
    parser.add_argument(
        "--grad-clip-norm",
        type=float,
        default=None,
        help="Clip adapter-head gradients to this global norm before optimizer step.",
    )
    parser.add_argument(
        "--warmup-steps",
        type=int,
        default=0,
        help="Linearly warm the learning rate from 0 to --learning-rate.",
    )
    parser.add_argument(
        "--lr-schedule",
        choices=["constant", "cosine", "step"],
        default="constant",
        help="Learning-rate schedule after optional warmup.",
    )
    parser.add_argument(
        "--step-decay-steps",
        default=None,
        help=(
            "Comma-separated optimizer steps for --lr-schedule step, e.g. 3000,4000."
        ),
    )
    parser.add_argument(
        "--step-decay-learning-rates",
        default=None,
        help=(
            "Comma-separated learning rates after each step boundary. A single value, "
            "e.g. 3e-4, is reused for all boundaries."
        ),
    )
    parser.add_argument(
        "--selection-metric",
        choices=["loss", "full-mse", "pearson", "common128-mse"],
        default="full-mse",
        help="Validation metric used for best checkpoint and early stopping.",
    )
    parser.add_argument(
        "--enable-last-block-lora",
        action="store_true",
        help="Add LoRA to the final transformer tower block MHA/MLP Linear layers.",
    )
    parser.add_argument("--lora-rank", type=int, default=4)
    parser.add_argument("--lora-alpha", type=int, default=8)
    parser.add_argument(
        "--freeze-head",
        action="store_true",
        help="Freeze the RNA-seq adapter head and train only enabled LoRA parameters.",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=20260514,
        help="Random seed. Use -1 to draw and record a non-deterministic seed.",
    )
    parser.add_argument("--organism-index", type=int, default=0)
    parser.add_argument(
        "--embedding-resolution",
        type=int,
        choices=[1, 128],
        default=128,
        help="Frozen AlphaGenome embedding resolution used by the adapter head.",
    )
    parser.add_argument(
        "--head-type",
        choices=["linear", "genome-tracks"],
        default="linear",
        help=(
            "linear keeps the legacy 1x1 Conv1d head; genome-tracks uses an "
            "AlphaGenome-style RNA-seq GenomeTracksHead."
        ),
    )
    parser.add_argument(
        "--head-resolutions",
        default=None,
        help=(
            "Comma-separated GenomeTracksHead resolutions, e.g. 128 or 1,128. "
            "Defaults to --embedding-resolution."
        ),
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
            "depthwise-separable-conv",
            "residual-conv1x1",
            "residual-conv3",
        ],
        default="conv1x1",
        help="Trainable simple-head architecture used when --head-type linear.",
    )
    parser.add_argument(
        "--linear-hidden-channels",
        type=int,
        default=256,
        help="Hidden channels for non-legacy linear head variants.",
    )
    parser.add_argument(
        "--linear-input-bottleneck-channels",
        type=int,
        default=None,
        help=(
            "Optional 1x1 bottleneck before pooling/interpolation for linear heads, "
            "for example 1536 -> 64 before 1 bp to 128 bp pooling."
        ),
    )
    parser.add_argument(
        "--residual-scale-init",
        type=float,
        default=1.0,
        help="Initial residual scale for residual linear heads.",
    )
    parser.add_argument(
        "--linear-dilation",
        type=int,
        default=2,
        help="Dilation for --linear-head-architecture dilated-conv3.",
    )
    parser.add_argument(
        "--linear-kernel-size",
        type=int,
        default=15,
        help="Kernel size for --linear-head-architecture depthwise-separable-conv.",
    )
    parser.add_argument(
        "--residual-base-checkpoint",
        default=None,
        help=(
            "Optional frozen linear checkpoint whose prediction is used as the "
            "base for 1 bp residual correction."
        ),
    )
    parser.add_argument(
        "--linear-residual-fusion",
        choices=["none", "base-prediction"],
        default="none",
        help=(
            "Optional fusion input for 1 bp residual correction. base-prediction "
            "concatenates the frozen residual-base prediction to the correction "
            "head input before predicting the residual update."
        ),
    )
    parser.add_argument(
        "--residual-correction-scale-init",
        type=float,
        default=0.01,
        help="Initial trainable output scale for --residual-base-checkpoint.",
    )
    parser.add_argument(
        "--residual-correction-l2",
        type=float,
        default=0.0,
        help=(
            "Training-only L2 penalty weight on the 1 bp residual update after "
            "scale/gating. Use with --residual-base-checkpoint."
        ),
    )
    parser.add_argument(
        "--residual-base-gate",
        choices=["none", "sigmoid"],
        default="none",
        help=(
            "Optional non-parametric gate on 1 bp residual correction based on "
            "the frozen base prediction."
        ),
    )
    parser.add_argument(
        "--residual-base-gate-center",
        type=float,
        default=0.5,
        help="Center used by --residual-base-gate sigmoid in log1p prediction space.",
    )
    parser.add_argument(
        "--residual-base-gate-sharpness",
        type=float,
        default=4.0,
        help="Positive sharpness used by --residual-base-gate sigmoid.",
    )
    parser.add_argument(
        "--residual-base-gate-floor",
        type=float,
        default=0.0,
        help="Minimum gate value for --residual-base-gate sigmoid.",
    )
    parser.add_argument(
        "--linear-output-calibration",
        choices=["none", "track-affine", "high-signal-gated"],
        default="none",
        help="Optional lightweight calibration applied after the linear head output.",
    )
    parser.add_argument(
        "--calibration-gate-center",
        type=float,
        default=3.0,
        help="Center for --linear-output-calibration high-signal-gated.",
    )
    parser.add_argument(
        "--calibration-gate-sharpness",
        type=float,
        default=3.0,
        help="Sharpness for --linear-output-calibration high-signal-gated.",
    )
    parser.add_argument(
        "--train-only-calibration",
        action="store_true",
        help="Freeze existing head modules and update only calibration parameters.",
    )
    parser.add_argument(
        "--linear-loss-type",
        choices=["mse", "smooth-l1", "hybrid", "hybrid-mse-smooth-l1"],
        default="mse",
        help="Loss for linear heads; GenomeTracks loss selection uses --loss-type.",
    )
    parser.add_argument(
        "--smooth-l1-beta",
        type=float,
        default=1.0,
        help="Beta parameter for --linear-loss-type smooth-l1.",
    )
    parser.add_argument(
        "--hybrid-loss-alpha",
        type=float,
        default=0.5,
        help=(
            "MSE weight alpha for --linear-loss-type hybrid; "
            "loss = alpha*MSE + (1-alpha)*SmoothL1."
        ),
    )
    parser.add_argument(
        "--track-loss-weight-preset",
        choices=[
            "uniform",
            "intestine-boost-1.25",
            "intestine-boost-1.5",
            "intestine-boost-2.0",
            "intestine-hard-boost-1.5",
            "intestine-hard-boost-2.0",
        ],
        default="uniform",
        help=(
            "Optional per-track training loss weighting preset for linear heads. "
            "Validation full MSE/Pearson/common128 metrics remain unweighted."
        ),
    )
    parser.add_argument(
        "--track-loss-weights",
        default=None,
        help=(
            "Optional comma-separated per-track training loss weights. Overrides "
            "--track-loss-weight-preset. Validation full metrics remain unweighted."
        ),
    )
    parser.add_argument(
        "--signal-loss-weight-preset",
        choices=[
            "none",
            "high-gt3-x1.5",
            "high-gt3-x2",
            "high-gt3-x3",
            "top5-x2",
            "top5-x3",
            "top1-x5",
            "intestine-high-x3",
        ],
        default="none",
        help=(
            "Optional elementwise training loss weighting from train-batch signal. "
            "Validation metrics remain unweighted."
        ),
    )
    parser.add_argument(
        "--signal-high-threshold",
        type=float,
        default=3.0,
        help="Loss-space threshold for high-gt3 signal weighting presets.",
    )
    parser.add_argument(
        "--signal-top-quantile",
        type=float,
        default=None,
        help=(
            "Optional dynamic train-batch quantile for top-signal presets. "
            "Defaults are 0.95 for top5 and 0.99 for top1."
        ),
    )
    parser.add_argument(
        "--region-loss-weight-preset",
        choices=["none", "exon-gene-v1", "exon-gene-v2"],
        default="none",
        help=(
            "Optional GTF-derived train-batch region weighting. Validation metrics "
            "remain unweighted."
        ),
    )
    parser.add_argument(
        "--gtf",
        default="alphagenome_custom/reference/Caenorhabditis_elegans.WBcel235.115.gtf",
        help="GTF used for train-only region weighting and gene auxiliary loss.",
    )
    parser.add_argument(
        "--promoter-radius-bp",
        type=int,
        default=1000,
        help="Promoter half-width around transcript TSS for region weighting.",
    )
    parser.add_argument(
        "--gene-aux-loss-weight",
        type=float,
        default=0.0,
        help="Training-only auxiliary MSE on per-gene gene-body mean predictions.",
    )
    parser.add_argument(
        "--gene-aux-min-bases",
        type=int,
        default=128,
        help="Minimum overlapped bases/bins for a gene to contribute aux loss.",
    )
    parser.add_argument(
        "--linear-target-space",
        choices=["full-log1p", "binned128-log1p-mean"],
        default="full-log1p",
        help=(
            "full-log1p keeps legacy full-resolution log1p targets; "
            "binned128-log1p-mean trains against 128 bp log1p(mean raw) bins."
        ),
    )
    parser.add_argument(
        "--track-means-source",
        choices=["ones", "grouped-qc", "train-nonzero"],
        default="grouped-qc",
        help=(
            "Track means for genome-tracks scaling. grouped-qc uses "
            "alphagenome_custom/metadata/grouped_bigwig_qc.tsv; train-nonzero "
            "streams raw train NPZ targets and computes nonzero means."
        ),
    )
    parser.add_argument(
        "--track-means-tsv",
        default="alphagenome_custom/metadata/grouped_bigwig_qc.tsv",
        help="TSV used when --track-means-source grouped-qc.",
    )
    parser.add_argument(
        "--track-means-max-examples",
        type=int,
        default=None,
        help="Optional train NPZ example limit for --track-means-source train-nonzero.",
    )
    parser.add_argument(
        "--target-transform",
        choices=["log1p", "none"],
        default="log1p",
        help="Use log1p for stable training on high dynamic range targets.",
    )
    parser.add_argument(
        "--loss-type",
        choices=["mse", "poisson-multinomial", "hybrid-mse-poisson"],
        default="mse",
        help=(
            "mse keeps the current masked MSE objective; poisson-multinomial "
            "uses AlphaGenome-style total-count plus positional count loss for "
            "GenomeTracksHead; hybrid-mse-poisson combines both with explicit "
            "weights."
        ),
    )
    parser.add_argument(
        "--multinomial-num-segments",
        type=int,
        default=8,
        help="Number of equal sequence segments for poisson-multinomial loss.",
    )
    parser.add_argument(
        "--positional-weight",
        type=float,
        default=5.0,
        help="Positional multinomial loss weight.",
    )
    parser.add_argument(
        "--count-weight",
        type=float,
        default=1.0,
        help="Poisson total-count loss weight.",
    )
    parser.add_argument(
        "--mse-weight",
        type=float,
        default=1.0,
        help="Model-space MSE weight for hybrid-mse-poisson loss.",
    )
    parser.add_argument(
        "--poisson-weight",
        type=float,
        default=1.0,
        help="Poisson-multinomial component weight for hybrid-mse-poisson loss.",
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
    parser.add_argument(
        "--early-stopping-patience",
        type=int,
        default=None,
        help="Stop after this many validations without improvement.",
    )
    parser.add_argument(
        "--early-stopping-min-steps",
        type=int,
        default=0,
        help="Do not early-stop before this optimizer step.",
    )
    return parser.parse_args()


def setup_distributed(requested_device: str) -> tuple[bool, int, int, int, torch.device]:
    world_size = int(os.environ.get("WORLD_SIZE", "1"))
    if world_size <= 1:
        return False, 0, 0, 1, pick_device(requested_device)

    if not torch.cuda.is_available():
        raise RuntimeError("Distributed training requires CUDA in this script")

    local_rank = int(os.environ["LOCAL_RANK"])
    rank = int(os.environ["RANK"])
    torch.cuda.set_device(local_rank)
    dist.init_process_group(backend="nccl")
    return True, rank, local_rank, world_size, torch.device(f"cuda:{local_rank}")


def cleanup_distributed(is_distributed: bool) -> None:
    if is_distributed and dist.is_initialized():
        dist.destroy_process_group()


def is_main_process(rank: int) -> bool:
    return rank == 0


def print_main(rank: int, *values: object) -> None:
    if is_main_process(rank):
        print(*values)


def reduce_mean(value: float, device: torch.device, world_size: int) -> float:
    if world_size <= 1:
        return value
    tensor = torch.tensor(value, dtype=torch.float64, device=device)
    dist.all_reduce(tensor, op=dist.ReduceOp.SUM)
    tensor /= world_size
    return float(tensor.detach().cpu())


def grad_norm(parameters) -> float:
    norms = [
        parameter.grad.detach().float().norm()
        for parameter in parameters
        if parameter.grad is not None
    ]
    if not norms:
        return 0.0
    return float(torch.linalg.vector_norm(torch.stack(norms)).cpu())


def lora_state_dict(model: torch.nn.Module) -> dict[str, torch.Tensor]:
    return {
        name: parameter.detach().cpu()
        for name, parameter in model.named_parameters()
        if parameter.requires_grad and ".lora_" in name
    }


def parse_int_list(value: str | None, *, option_name: str) -> list[int]:
    if value is None:
        return []
    try:
        values = [int(part.strip()) for part in value.split(",") if part.strip()]
    except ValueError as error:
        raise ValueError(f"{option_name} must be a comma-separated integer list") from error
    if any(item < 1 for item in values):
        raise ValueError(f"{option_name} values must be >= 1")
    if values != sorted(values):
        raise ValueError(f"{option_name} values must be sorted ascending")
    return values


def parse_float_list(value: str | None, *, option_name: str) -> list[float]:
    if value is None:
        return []
    try:
        values = [float(part.strip()) for part in value.split(",") if part.strip()]
    except ValueError as error:
        raise ValueError(f"{option_name} must be a comma-separated float list") from error
    if any(item < 0.0 for item in values):
        raise ValueError(f"{option_name} values must be non-negative")
    return values


def resolve_track_loss_weights(
    *,
    weights_text: str | None,
    preset: str,
    n_tracks: int,
) -> list[float] | None:
    """Return optional per-track training weights; None means uniform."""

    if weights_text:
        weights = parse_float_list(weights_text, option_name="--track-loss-weights")
        if len(weights) != n_tracks:
            raise ValueError(
                f"--track-loss-weights has {len(weights)} values; expected {n_tracks}"
            )
    else:
        weights = [1.0] * n_tracks
        if preset == "uniform":
            pass
        elif preset.startswith("intestine-boost-"):
            if n_tracks != 11:
                raise ValueError(f"{preset} expects 11 RNA-seq tracks")
            boost = float(preset.rsplit("-", maxsplit=1)[1])
            for index in range(5):
                weights[index] = boost
        elif preset.startswith("intestine-hard-boost-"):
            if n_tracks != 11:
                raise ValueError(f"{preset} expects 11 RNA-seq tracks")
            boost = float(preset.rsplit("-", maxsplit=1)[1])
            for index in (1, 3, 4):
                weights[index] = boost
        else:
            raise ValueError(f"Unsupported --track-loss-weight-preset: {preset}")

    if any(weight <= 0.0 for weight in weights):
        raise ValueError("Track loss weights must be > 0")
    if all(abs(weight - 1.0) < 1e-12 for weight in weights):
        return None
    return weights


def format_track_loss_weights(weights: list[float] | None) -> str:
    if weights is None:
        return "uniform"
    return ",".join(f"{weight:.8g}" for weight in weights)


@dataclass(frozen=True)
class GenomicInterval:
    chrom: str
    start: int
    end: int


@dataclass(frozen=True)
class GtfAnnotations:
    genes: dict[str, list[GenomicInterval]]
    exons: dict[str, list[GenomicInterval]]
    promoters: dict[str, list[GenomicInterval]]


def parse_gtf_attributes(value: str) -> dict[str, str]:
    attributes: dict[str, str] = {}
    for part in value.strip().split(";"):
        part = part.strip()
        if not part:
            continue
        if " " not in part:
            attributes[part] = ""
            continue
        key, raw_value = part.split(" ", maxsplit=1)
        attributes[key] = raw_value.strip().strip('"')
    return attributes


def merge_genomic_intervals(intervals: list[GenomicInterval]) -> list[GenomicInterval]:
    merged: list[GenomicInterval] = []
    for interval in sorted(intervals, key=lambda item: (item.chrom, item.start, item.end)):
        if (
            not merged
            or interval.chrom != merged[-1].chrom
            or interval.start > merged[-1].end
        ):
            merged.append(interval)
            continue
        previous = merged[-1]
        merged[-1] = GenomicInterval(
            chrom=previous.chrom,
            start=previous.start,
            end=max(previous.end, interval.end),
        )
    return merged


def load_gtf_annotations_for_training(
    path: Path,
    *,
    promoter_radius_bp: int,
) -> GtfAnnotations:
    if promoter_radius_bp < 0:
        raise ValueError("--promoter-radius-bp must be >= 0")
    genes: dict[str, list[GenomicInterval]] = {}
    exons: dict[str, list[GenomicInterval]] = {}
    promoters: dict[str, list[GenomicInterval]] = {}

    with path.open() as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            chrom, _source, feature, start_s, end_s, _score, strand, _frame, attr_s = fields
            try:
                start = int(start_s) - 1
                end = int(end_s)
            except ValueError:
                continue
            if end <= start:
                continue
            if feature == "gene":
                genes.setdefault(chrom, []).append(
                    GenomicInterval(chrom=chrom, start=start, end=end)
                )
            elif feature == "exon":
                exons.setdefault(chrom, []).append(
                    GenomicInterval(chrom=chrom, start=start, end=end)
                )
            elif feature == "transcript":
                _attrs = parse_gtf_attributes(attr_s)
                tss = start if strand != "-" else end
                promoter_start = max(0, tss - promoter_radius_bp)
                promoter_end = tss + promoter_radius_bp
                if promoter_end > promoter_start:
                    promoters.setdefault(chrom, []).append(
                        GenomicInterval(
                            chrom=chrom,
                            start=promoter_start,
                            end=promoter_end,
                        )
                    )

    for container in (genes, exons, promoters):
        for chrom, intervals in list(container.items()):
            container[chrom] = merge_genomic_intervals(intervals)
    return GtfAnnotations(genes=genes, exons=exons, promoters=promoters)


def batch_chromosome(batch: dict[str, object], index: int) -> str:
    chromosomes = batch["interval_chromosome"]
    if isinstance(chromosomes, (list, tuple)):
        return str(chromosomes[index])
    return str(chromosomes)


def batch_interval_start_end(batch: dict[str, object], index: int) -> tuple[int, int]:
    starts = batch["interval_start"]
    ends = batch["interval_end"]
    if torch.is_tensor(starts):
        start = int(starts[index].item())
    elif isinstance(starts, (list, tuple)):
        start = int(starts[index])
    else:
        start = int(starts)
    if torch.is_tensor(ends):
        end = int(ends[index].item())
    elif isinstance(ends, (list, tuple)):
        end = int(ends[index])
    else:
        end = int(ends)
    return start, end


def interval_mask_for_window(
    intervals: dict[str, list[GenomicInterval]],
    *,
    chrom: str,
    window_start: int,
    window_end: int,
    target_length: int,
) -> np.ndarray:
    mask = np.zeros(target_length, dtype=bool)
    window_width = window_end - window_start
    if target_length <= 0 or window_width <= 0:
        return mask
    scale = window_width / float(target_length)
    for interval in intervals.get(chrom, []):
        if interval.end <= window_start:
            continue
        if interval.start >= window_end:
            break
        overlap_start = max(interval.start, window_start)
        overlap_end = min(interval.end, window_end)
        if overlap_end <= overlap_start:
            continue
        local_start = overlap_start - window_start
        local_end = overlap_end - window_start
        if target_length == window_width:
            start_index = int(local_start)
            end_index = int(local_end)
        else:
            start_index = int(math.floor(local_start / scale))
            end_index = int(math.ceil(local_end / scale))
        start_index = max(0, min(target_length, start_index))
        end_index = max(0, min(target_length, end_index))
        if end_index > start_index:
            mask[start_index:end_index] = True
    return mask


def region_preset_weights(preset: str) -> dict[str, float]:
    if preset == "exon-gene-v1":
        return {
            "exon": 2.0,
            "gene_body": 1.5,
            "promoter": 1.25,
            "intron": 1.0,
            "intergenic": 0.75,
        }
    if preset == "exon-gene-v2":
        return {
            "exon": 1.5,
            "gene_body": 1.25,
            "promoter": 1.25,
            "intron": 1.0,
            "intergenic": 1.0,
        }
    raise ValueError(f"Unsupported --region-loss-weight-preset: {preset}")


def build_region_loss_weights(
    *,
    target: torch.Tensor,
    batch: dict[str, object],
    annotations: GtfAnnotations,
    preset: str,
) -> torch.Tensor | None:
    if preset == "none":
        return None
    preset_weights = region_preset_weights(preset)
    batch_size = int(target.shape[0])
    target_length = int(target.shape[-1])
    sample_weights: list[torch.Tensor] = []
    for batch_index in range(batch_size):
        chrom = batch_chromosome(batch, batch_index)
        window_start, window_end = batch_interval_start_end(batch, batch_index)
        gene_mask = interval_mask_for_window(
            annotations.genes,
            chrom=chrom,
            window_start=window_start,
            window_end=window_end,
            target_length=target_length,
        )
        exon_mask = interval_mask_for_window(
            annotations.exons,
            chrom=chrom,
            window_start=window_start,
            window_end=window_end,
            target_length=target_length,
        )
        promoter_mask = interval_mask_for_window(
            annotations.promoters,
            chrom=chrom,
            window_start=window_start,
            window_end=window_end,
            target_length=target_length,
        )
        intron_mask = gene_mask & ~exon_mask
        weights = np.full(
            target_length,
            preset_weights["intergenic"],
            dtype=np.float32,
        )
        weights[intron_mask] = np.maximum(weights[intron_mask], preset_weights["intron"])
        weights[promoter_mask] = np.maximum(
            weights[promoter_mask],
            preset_weights["promoter"],
        )
        weights[gene_mask] = np.maximum(weights[gene_mask], preset_weights["gene_body"])
        weights[exon_mask] = np.maximum(weights[exon_mask], preset_weights["exon"])
        sample_weights.append(
            torch.from_numpy(weights).to(device=target.device, dtype=target.dtype)
        )
    return torch.stack(sample_weights, dim=0).unsqueeze(1)


def signal_preset_factor(preset: str) -> float:
    factors = {
        "high-gt3-x1.5": 1.5,
        "high-gt3-x2": 2.0,
        "high-gt3-x3": 3.0,
        "top5-x2": 2.0,
        "top5-x3": 3.0,
        "top1-x5": 5.0,
        "intestine-high-x3": 3.0,
    }
    try:
        return factors[preset]
    except KeyError as error:
        raise ValueError(f"Unsupported --signal-loss-weight-preset: {preset}") from error


def build_signal_loss_weights(
    *,
    target: torch.Tensor,
    preset: str,
    high_threshold: float,
    top_quantile: float | None,
) -> torch.Tensor | None:
    if preset == "none":
        return None
    if high_threshold < 0.0:
        raise ValueError("--signal-high-threshold must be >= 0")
    factor = signal_preset_factor(preset)
    target_for_mask = target.detach()

    if preset.startswith("high-gt3-"):
        selected = target_for_mask >= high_threshold
    elif preset == "intestine-high-x3":
        if int(target.shape[1]) != 11:
            raise ValueError("intestine-high-x3 expects 11 RNA-seq tracks")
        hard_tracks = torch.zeros(
            int(target.shape[1]),
            dtype=torch.bool,
            device=target.device,
        )
        hard_tracks[torch.tensor([1, 3, 4], device=target.device)] = True
        selected = (target_for_mask >= high_threshold) & hard_tracks.view(1, -1, 1)
    elif preset.startswith("top5-") or preset.startswith("top1-"):
        quantile = top_quantile
        if quantile is None:
            quantile = 0.99 if preset.startswith("top1-") else 0.95
        if not 0.0 < quantile < 1.0:
            raise ValueError("--signal-top-quantile must be between 0 and 1")
        thresholds = torch.quantile(
            target_for_mask.to(dtype=torch.float32),
            quantile,
            dim=-1,
            keepdim=True,
        ).to(dtype=target.dtype)
        selected = target_for_mask >= thresholds
    else:
        raise ValueError(f"Unsupported --signal-loss-weight-preset: {preset}")

    return torch.where(
        selected,
        torch.as_tensor(factor, dtype=target.dtype, device=target.device),
        torch.ones((), dtype=target.dtype, device=target.device),
    )


def combine_element_loss_weights(
    *weights: torch.Tensor | None,
) -> torch.Tensor | None:
    combined: torch.Tensor | None = None
    for weight in weights:
        if weight is None:
            continue
        combined = weight if combined is None else combined * weight
    return combined


def gene_body_overlaps_for_window(
    annotations: GtfAnnotations,
    *,
    chrom: str,
    window_start: int,
    window_end: int,
    target_length: int,
    min_bases: int,
) -> list[tuple[int, int]]:
    overlaps: list[tuple[int, int]] = []
    window_width = window_end - window_start
    if target_length <= 0 or window_width <= 0:
        return overlaps
    scale = window_width / float(target_length)
    for interval in annotations.genes.get(chrom, []):
        if interval.end <= window_start:
            continue
        if interval.start >= window_end:
            break
        overlap_start = max(interval.start, window_start)
        overlap_end = min(interval.end, window_end)
        if overlap_end <= overlap_start:
            continue
        local_start = overlap_start - window_start
        local_end = overlap_end - window_start
        if target_length == window_width:
            start_index = int(local_start)
            end_index = int(local_end)
        else:
            start_index = int(math.floor(local_start / scale))
            end_index = int(math.ceil(local_end / scale))
        start_index = max(0, min(target_length, start_index))
        end_index = max(0, min(target_length, end_index))
        if end_index - start_index >= min_bases:
            overlaps.append((start_index, end_index))
    return overlaps


def gene_body_mean_aux_loss(
    prediction: torch.Tensor,
    target: torch.Tensor,
    mask: torch.Tensor,
    *,
    batch: dict[str, object],
    annotations: GtfAnnotations,
    min_bases: int,
) -> torch.Tensor:
    if min_bases < 1:
        raise ValueError("--gene-aux-min-bases must be >= 1")
    while mask.ndim < prediction.ndim:
        mask = mask.unsqueeze(-1)
    mask = mask.to(device=prediction.device, dtype=prediction.dtype)

    total = prediction.new_tensor(0.0)
    count = prediction.new_tensor(0.0)
    target_length = int(prediction.shape[-1])
    for batch_index in range(int(prediction.shape[0])):
        chrom = batch_chromosome(batch, batch_index)
        window_start, window_end = batch_interval_start_end(batch, batch_index)
        overlaps = gene_body_overlaps_for_window(
            annotations,
            chrom=chrom,
            window_start=window_start,
            window_end=window_end,
            target_length=target_length,
            min_bases=min_bases,
        )
        if not overlaps:
            continue
        track_mask = mask[batch_index, :, 0]
        for start_index, end_index in overlaps:
            prediction_mean = prediction[batch_index, :, start_index:end_index].mean(dim=-1)
            target_mean = target[batch_index, :, start_index:end_index].mean(dim=-1)
            error = (prediction_mean - target_mean.to(dtype=prediction.dtype)).square()
            total = total + (error * track_mask).sum()
            count = count + track_mask.sum()
    return total / count.clamp_min(1.0)


def learning_rate_for_step(
    *,
    step: int,
    base_lr: float,
    max_steps: int,
    warmup_steps: int,
    schedule: str,
    step_decay_steps: list[int],
    step_decay_learning_rates: list[float],
) -> float:
    if warmup_steps > 0 and step <= warmup_steps:
        return base_lr * (step / warmup_steps)

    if schedule == "step":
        lr = base_lr
        for boundary, boundary_lr in zip(step_decay_steps, step_decay_learning_rates):
            if step >= boundary:
                lr = boundary_lr
        return lr

    if schedule == "cosine":
        decay_steps = max(max_steps - warmup_steps, 1)
        decay_step = min(max(step - warmup_steps, 0), decay_steps)
        cosine = 0.5 * (1.0 + math.cos(math.pi * decay_step / decay_steps))
        return base_lr * cosine

    return base_lr


def masked_mse_loss_local(
    prediction: torch.Tensor,
    target: torch.Tensor,
    mask: torch.Tensor,
    track_loss_weights: torch.Tensor | None = None,
    loss_weights: torch.Tensor | None = None,
) -> torch.Tensor:
    while mask.ndim < prediction.ndim:
        mask = mask.unsqueeze(-1)
    mask = mask.to(device=prediction.device, dtype=prediction.dtype)
    error = (prediction - target.to(dtype=prediction.dtype)).square()
    denominator = mask.expand_as(error)
    if loss_weights is not None:
        weights = loss_weights.to(device=prediction.device, dtype=prediction.dtype)
        if weights.shape != prediction.shape:
            weights = weights.expand_as(prediction)
        error = error * weights
        denominator = denominator * weights
    if track_loss_weights is not None:
        if prediction.ndim != 3:
            raise ValueError("Track loss weights require [batch, tracks, length] tensors")
        if int(track_loss_weights.numel()) != int(prediction.shape[1]):
            raise ValueError(
                "Track loss weights have "
                f"{track_loss_weights.numel()} tracks; expected {prediction.shape[1]}"
            )
        weights = track_loss_weights.to(device=prediction.device, dtype=prediction.dtype)
        weights = weights.view(1, -1, 1)
        error = error * weights
        denominator = denominator * weights
    return (error * mask).sum() / denominator.sum().clamp_min(1.0)


def masked_smooth_l1_loss_local(
    prediction: torch.Tensor,
    target: torch.Tensor,
    mask: torch.Tensor,
    *,
    beta: float,
    track_loss_weights: torch.Tensor | None = None,
    loss_weights: torch.Tensor | None = None,
) -> torch.Tensor:
    while mask.ndim < prediction.ndim:
        mask = mask.unsqueeze(-1)
    mask = mask.to(device=prediction.device, dtype=prediction.dtype)
    loss = F.smooth_l1_loss(
        prediction,
        target.to(dtype=prediction.dtype),
        reduction="none",
        beta=beta,
    )
    denominator = mask.expand_as(loss)
    if loss_weights is not None:
        weights = loss_weights.to(device=prediction.device, dtype=prediction.dtype)
        if weights.shape != prediction.shape:
            weights = weights.expand_as(prediction)
        loss = loss * weights
        denominator = denominator * weights
    if track_loss_weights is not None:
        if prediction.ndim != 3:
            raise ValueError("Track loss weights require [batch, tracks, length] tensors")
        if int(track_loss_weights.numel()) != int(prediction.shape[1]):
            raise ValueError(
                "Track loss weights have "
                f"{track_loss_weights.numel()} tracks; expected {prediction.shape[1]}"
            )
        weights = track_loss_weights.to(device=prediction.device, dtype=prediction.dtype)
        weights = weights.view(1, -1, 1)
        loss = loss * weights
        denominator = denominator * weights
    return (loss * mask).sum() / denominator.sum().clamp_min(1.0)


def masked_residual_update_l2(
    residual_update: torch.Tensor,
    mask: torch.Tensor,
) -> torch.Tensor:
    while mask.ndim < residual_update.ndim:
        mask = mask.unsqueeze(-1)
    mask = mask.to(device=residual_update.device, dtype=residual_update.dtype)
    penalty = residual_update.square()
    return (penalty * mask).sum() / mask.expand_as(penalty).sum().clamp_min(1.0)


def adapter_loss_for_training(
    model: RnaSeq11Adapter,
    prediction: torch.Tensor | dict[int, torch.Tensor],
    target: torch.Tensor,
    mask: torch.Tensor,
    *,
    loss_type: str,
    multinomial_num_segments: int,
    positional_weight: float,
    count_weight: float,
    mse_weight: float,
    poisson_weight: float,
    linear_loss_type: str,
    smooth_l1_beta: float,
    hybrid_loss_alpha: float,
    track_loss_weights: torch.Tensor | None = None,
    loss_weights: torch.Tensor | None = None,
) -> torch.Tensor:
    if isinstance(prediction, dict):
        if track_loss_weights is not None or loss_weights is not None:
            raise ValueError("Weighted training losses are only supported for linear heads")
        return masked_adapter_loss(
            model,
            prediction,
            target,
            mask,
            loss_type=loss_type,
            multinomial_num_segments=multinomial_num_segments,
            positional_weight=positional_weight,
            count_weight=count_weight,
            mse_weight=mse_weight,
            poisson_weight=poisson_weight,
            linear_loss_type=linear_loss_type,
            smooth_l1_beta=smooth_l1_beta,
            hybrid_loss_alpha=hybrid_loss_alpha,
        )

    if (
        track_loss_weights is None
        and loss_weights is None
        and linear_loss_type != "hybrid-mse-smooth-l1"
    ):
        return masked_adapter_loss(
            model,
            prediction,
            target,
            mask,
            loss_type=loss_type,
            multinomial_num_segments=multinomial_num_segments,
            positional_weight=positional_weight,
            count_weight=count_weight,
            mse_weight=mse_weight,
            poisson_weight=poisson_weight,
            linear_loss_type=linear_loss_type,
            smooth_l1_beta=smooth_l1_beta,
            hybrid_loss_alpha=hybrid_loss_alpha,
        )

    if linear_loss_type == "mse":
        return masked_mse_loss_local(
            prediction,
            target,
            mask,
            track_loss_weights=track_loss_weights,
            loss_weights=loss_weights,
        )
    if linear_loss_type == "smooth-l1":
        return masked_smooth_l1_loss_local(
            prediction,
            target,
            mask,
            beta=smooth_l1_beta,
            track_loss_weights=track_loss_weights,
            loss_weights=loss_weights,
        )
    if linear_loss_type not in {"hybrid", "hybrid-mse-smooth-l1"}:
        raise ValueError(f"Unsupported linear loss type: {linear_loss_type}")

    mse_loss = masked_mse_loss_local(
        prediction,
        target,
        mask,
        track_loss_weights=track_loss_weights,
        loss_weights=loss_weights,
    )
    smooth_l1_loss = masked_smooth_l1_loss_local(
        prediction,
        target,
        mask,
        beta=smooth_l1_beta,
        track_loss_weights=track_loss_weights,
        loss_weights=loss_weights,
    )
    return hybrid_loss_alpha * mse_loss + (1.0 - hybrid_loss_alpha) * smooth_l1_loss


def masked_metric_sums(
    prediction: torch.Tensor,
    target: torch.Tensor,
    mask: torch.Tensor,
) -> dict[str, float]:
    if prediction.shape != target.shape:
        raise ValueError(f"Prediction/target shape mismatch: {prediction.shape} vs {target.shape}")
    mask = metric_mask_for_prediction(mask, prediction)
    values_x = prediction.to(dtype=torch.float64)[mask]
    values_y = target.to(device=prediction.device, dtype=torch.float64)[mask]
    error = values_x - values_y
    return {
        "sse": float(error.square().sum().detach().cpu()),
        "sae": float(error.abs().sum().detach().cpu()),
        "count": float(values_x.numel()),
        "sum_x": float(values_x.sum().detach().cpu()),
        "sum_y": float(values_y.sum().detach().cpu()),
        "sum_x2": float(values_x.square().sum().detach().cpu()),
        "sum_y2": float(values_y.square().sum().detach().cpu()),
        "sum_xy": float((values_x * values_y).sum().detach().cpu()),
    }


def metric_mask_for_prediction(mask: torch.Tensor, prediction: torch.Tensor) -> torch.Tensor:
    while mask.ndim < prediction.ndim:
        mask = mask.unsqueeze(-1)
    mask = mask.to(device=prediction.device, dtype=torch.bool)
    if mask.shape[-1] == 1:
        return mask.expand_as(prediction)
    if mask.shape[-1] != prediction.shape[-1]:
        if mask.shape[-1] % prediction.shape[-1] != 0:
            raise ValueError(
                f"Mask/prediction length mismatch: {mask.shape[-1]} vs {prediction.shape[-1]}"
            )
        factor = mask.shape[-1] // prediction.shape[-1]
        usable_length = prediction.shape[-1] * factor
        mask = mask[..., :usable_length]
        mask = mask.reshape(*mask.shape[:-1], prediction.shape[-1], factor).any(dim=-1)
    return mask.expand_as(prediction)


def add_metric_sums(total: dict[str, float], update: dict[str, float]) -> None:
    for key, value in update.items():
        total[key] = total.get(key, 0.0) + value


def primary_prediction_tensor(
    model: RnaSeq11Adapter,
    prediction: torch.Tensor | dict[int, torch.Tensor],
    *,
    target_length: int,
) -> torch.Tensor:
    if not isinstance(prediction, dict):
        if prediction.shape[-1] == target_length:
            return prediction
        return F.interpolate(prediction, size=target_length, mode="linear", align_corners=False)

    if 1 in prediction:
        return prediction[1]
    if 128 not in prediction:
        raise ValueError(f"No supported prediction resolution in {sorted(prediction)}")
    return F.interpolate(
        prediction[128] / 128.0,
        size=target_length,
        mode="linear",
        align_corners=False,
    )


def bin_mean_128(tensor: torch.Tensor) -> torch.Tensor:
    usable_length = (tensor.shape[-1] // 128) * 128
    if usable_length == 0:
        raise ValueError(f"Sequence length {tensor.shape[-1]} is shorter than 128 bp")
    tensor = tensor[..., :usable_length]
    return tensor.reshape(*tensor.shape[:-1], usable_length // 128, 128).mean(dim=-1)


def transform_full_target(raw_target: torch.Tensor, target_transform: str) -> torch.Tensor:
    if target_transform == "log1p":
        return torch.log1p(raw_target.clamp_min(0.0))
    if target_transform == "none":
        return raw_target
    raise ValueError(f"Unsupported target transform: {target_transform}")


def prediction_common128_tensor(
    model: RnaSeq11Adapter,
    prediction: torch.Tensor | dict[int, torch.Tensor],
    *,
    linear_prediction_transform: str,
    linear_target_space: str,
) -> torch.Tensor:
    if isinstance(prediction, dict):
        if 128 in prediction:
            prediction_128 = prediction[128].clamp_min(0.0)
        elif 1 in prediction:
            prediction_128 = bin_rna_seq_target(prediction[1], 128)
        else:
            raise ValueError(f"No supported prediction resolution in {sorted(prediction)}")
        return torch.log1p(prediction_128 / 128.0)

    if linear_target_space == "binned128-log1p-mean":
        return prediction

    prediction_raw = prediction
    if linear_prediction_transform == "log1p":
        prediction_raw = torch.expm1(prediction_raw).clamp_min(0.0)
    elif linear_prediction_transform == "none":
        prediction_raw = prediction_raw.clamp_min(0.0)
    else:
        raise ValueError(
            "linear_prediction_transform must be 'log1p' or 'none'"
        )
    return torch.log1p(bin_rna_seq_target(prediction_raw, 128) / 128.0)


def selection_value(metrics: dict[str, float], selection_metric: str) -> float:
    if selection_metric == "loss":
        return metrics["loss"]
    if selection_metric == "full-mse":
        return metrics["full_mse"]
    if selection_metric == "pearson":
        return metrics["pearson"]
    if selection_metric == "common128-mse":
        return metrics["common128_mse"]
    raise ValueError(f"Unsupported selection metric: {selection_metric}")


def selection_improved(
    current: float,
    best: float | None,
    *,
    selection_metric: str,
) -> bool:
    if best is None:
        return True
    if selection_metric == "pearson":
        return current > best
    return current < best


def evaluate_validation_metrics(
    model: RnaSeq11Adapter,
    dataloader: torch.utils.data.DataLoader,
    *,
    device: torch.device,
    loss_type: str,
    multinomial_num_segments: int,
    positional_weight: float,
    count_weight: float,
    mse_weight: float,
    poisson_weight: float,
    linear_loss_type: str,
    linear_target_space: str,
    target_transform: str,
    smooth_l1_beta: float,
    hybrid_loss_alpha: float,
    track_loss_weights: torch.Tensor | None = None,
) -> dict[str, float | int]:
    model.eval()
    weighted_loss_sum = 0.0
    full_sums: dict[str, float] = {}
    common128_sums: dict[str, float] = {}
    n_batches = 0
    n_examples = 0

    with torch.no_grad():
        for batch in dataloader:
            dna_sequence = batch["dna_sequence"].to(device=device, dtype=torch.float32)
            raw_rna_seq = batch["rna_seq"].to(device=device, dtype=torch.float32)
            rna_seq_mask = batch["rna_seq_mask"].to(device=device)
            full_target = transform_full_target(raw_rna_seq, target_transform)
            if model.head_type == "linear":
                if linear_target_space == "binned128-log1p-mean":
                    loss_target = prepare_linear_target(
                        raw_rna_seq,
                        target_space=linear_target_space,
                    )
                else:
                    loss_target = full_target
            else:
                loss_target = raw_rna_seq

            return_scaled_for_loss = prediction_return_scaled_for_loss(model, loss_type)
            loss_prediction = model(
                dna_sequence,
                target_length=loss_target.shape[-1],
                return_scaled=return_scaled_for_loss,
            )
            loss = adapter_loss_for_training(
                model,
                loss_prediction,
                loss_target if model.head_type == "linear" else raw_rna_seq,
                rna_seq_mask,
                loss_type=loss_type,
                multinomial_num_segments=multinomial_num_segments,
                positional_weight=positional_weight,
                count_weight=count_weight,
                mse_weight=mse_weight,
                poisson_weight=poisson_weight,
                linear_loss_type=linear_loss_type,
                smooth_l1_beta=smooth_l1_beta,
                hybrid_loss_alpha=hybrid_loss_alpha,
                track_loss_weights=track_loss_weights,
            )

            metric_prediction = (
                model(
                    dna_sequence,
                    target_length=loss_target.shape[-1],
                    return_scaled=False,
                )
                if return_scaled_for_loss else loss_prediction
            )
            full_prediction = primary_prediction_tensor(
                model,
                metric_prediction,
                target_length=raw_rna_seq.shape[-1],
            )
            add_metric_sums(
                full_sums,
                masked_metric_sums(full_prediction, full_target, rna_seq_mask),
            )

            common128_prediction = prediction_common128_tensor(
                model,
                metric_prediction,
                linear_prediction_transform=target_transform,
                linear_target_space=linear_target_space,
            )
            common128_target = torch.log1p(bin_rna_seq_target(raw_rna_seq, 128) / 128.0)
            add_metric_sums(
                common128_sums,
                masked_metric_sums(
                    common128_prediction,
                    common128_target,
                    rna_seq_mask,
                ),
            )

            batch_size = int(dna_sequence.shape[0])
            weighted_loss_sum += float(loss.detach().cpu()) * batch_size
            n_batches += 1
            n_examples += batch_size

    model.train()
    model.base_model.eval()

    full_count = max(full_sums.get("count", 0.0), 1.0)
    common_count = max(common128_sums.get("count", 0.0), 1.0)
    pearson_denominator = math.sqrt(
        max(full_sums.get("sum_x2", 0.0) - full_sums.get("sum_x", 0.0) ** 2 / full_count, 0.0)
        * max(full_sums.get("sum_y2", 0.0) - full_sums.get("sum_y", 0.0) ** 2 / full_count, 0.0)
    )
    pearson_numerator = (
        full_sums.get("sum_xy", 0.0)
        - full_sums.get("sum_x", 0.0) * full_sums.get("sum_y", 0.0) / full_count
    )
    pearson = pearson_numerator / pearson_denominator if pearson_denominator > 0.0 else float("nan")

    return {
        "loss": weighted_loss_sum / max(n_examples, 1),
        "full_mse": full_sums.get("sse", 0.0) / full_count,
        "full_mae": full_sums.get("sae", 0.0) / full_count,
        "pearson": pearson,
        "common128_mse": common128_sums.get("sse", 0.0) / common_count,
        "common128_mae": common128_sums.get("sae", 0.0) / common_count,
        "batches": n_batches,
        "examples": n_examples,
    }


def loss_space_for_config(
    head_type: str,
    loss_type: str,
    linear_target_space: str = "full-log1p",
    linear_loss_type: str = "mse",
) -> str:
    if head_type == "genome-tracks" and loss_type == "mse":
        return "genome_tracks_model_scaled"
    if loss_type == "poisson-multinomial":
        return "experimental_counts"
    if loss_type == "hybrid-mse-poisson":
        return "hybrid_genome_tracks_mse_plus_counts"
    if head_type == "linear":
        return f"{linear_target_space}_{linear_loss_type}"
    return "target_transform"


def infinite_batches(
    dataloader: torch.utils.data.DataLoader,
    sampler: DistributedSampler | None = None,
):
    epoch = 0
    while True:
        if sampler is not None:
            sampler.set_epoch(epoch)
        yield from dataloader
        epoch += 1


def set_seed(seed: int) -> None:
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)


def resolve_seed(seed: int) -> int:
    if seed >= 0:
        return seed
    return int.from_bytes(os.urandom(4), byteorder="big") % (2**31 - 1)


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
        "distributed": dist.is_available() and dist.is_initialized(),
        "distributed_world_size": (
            dist.get_world_size() if dist.is_available() and dist.is_initialized() else 1
        ),
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


def load_matching_adapter_head_state_dict(
    model: RnaSeq11Adapter,
    state_dict: dict[str, torch.Tensor],
) -> tuple[list[str], list[str]]:
    """Load adapter parameters whose names and shapes match the current model."""

    target_state = model.adapter_head_state_dict()
    merged_state: dict[str, torch.Tensor] = {}
    loaded_keys: list[str] = []
    skipped_keys: list[str] = []
    for key, target_value in target_state.items():
        source_value = state_dict.get(key)
        if source_value is not None and tuple(source_value.shape) == tuple(
            target_value.shape
        ):
            merged_state[key] = source_value
            loaded_keys.append(key)
        else:
            merged_state[key] = target_value
            skipped_keys.append(key)
    model.load_adapter_head_state_dict(merged_state)
    return loaded_keys, skipped_keys


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
    full_mse: float | None = None,
    full_mae: float | None = None,
    pearson: float | None = None,
    common128_mse: float | None = None,
    common128_mae: float | None = None,
    selection_metric: str | None = None,
    selection_value_for_row: float | None = None,
) -> None:
    fieldnames = [
        "split",
        "step",
        "loss",
        "full_mse",
        "full_mae",
        "pearson",
        "common128_mse",
        "common128_mae",
        "selection_metric",
        "selection_value",
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
                "full_mse": "" if full_mse is None else f"{full_mse:.8g}",
                "full_mae": "" if full_mae is None else f"{full_mae:.8g}",
                "pearson": "" if pearson is None else f"{pearson:.8g}",
                "common128_mse": (
                    "" if common128_mse is None else f"{common128_mse:.8g}"
                ),
                "common128_mae": (
                    "" if common128_mae is None else f"{common128_mae:.8g}"
                ),
                "selection_metric": "" if selection_metric is None else selection_metric,
                "selection_value": (
                    "" if selection_value_for_row is None
                    else f"{selection_value_for_row:.8g}"
                ),
                "batches": batches,
                "examples": examples,
                "learning_rate": f"{lr:.8g}",
                "cuda_max_memory_allocated_mb": (
                    "" if cuda_max_memory_allocated_mb is None
                    else f"{cuda_max_memory_allocated_mb:.1f}"
                ),
            }
        )


def save_training_checkpoint(
    path: Path,
    *,
    model: RnaSeq11Adapter,
    n_tracks: int,
    embedding_resolution: int,
    organism_index: int,
    target_transform: str,
    weights: Path,
    best_valid_loss: float | None,
    best_valid_step: int | None,
    best_selection_metric: str,
    best_selection_value: float | None,
    lora_enabled: bool,
    lora_target_modules: list[str],
    lora_rank: int,
    lora_alpha: int,
    freeze_head: bool,
    init_adapter_checkpoint: str | None,
    init_adapter_checkpoint_mode: str,
    head_type: str,
    head_resolutions: tuple[int, ...],
    linear_head_architecture: str,
    linear_hidden_channels: int,
    linear_input_bottleneck_channels: int | None,
    linear_residual_scale_init: float,
    linear_dilation: int,
    linear_kernel_size: int,
    residual_base_checkpoint: str | None,
    linear_residual_fusion: str,
    residual_correction_scale_init: float,
    residual_correction_l2: float,
    residual_base_gate: str,
    residual_base_gate_center: float,
    residual_base_gate_sharpness: float,
    residual_base_gate_floor: float,
    linear_output_calibration: str,
    calibration_gate_center: float,
    calibration_gate_sharpness: float,
    train_only_calibration: bool,
    linear_loss_type: str,
    linear_target_space: str,
    smooth_l1_beta: float,
    hybrid_loss_alpha: float,
    track_loss_weight_preset: str,
    track_loss_weights: list[float] | None,
    signal_loss_weight_preset: str,
    signal_high_threshold: float,
    signal_top_quantile: float | None,
    region_loss_weight_preset: str,
    gtf: str,
    promoter_radius_bp: int,
    gene_aux_loss_weight: float,
    gene_aux_min_bases: int,
    loss_space: str,
    loss_type: str,
    multinomial_num_segments: int,
    positional_weight: float,
    count_weight: float,
    mse_weight: float,
    poisson_weight: float,
    track_means_source: str,
    track_means_tsv: str,
    track_means_max_examples: int | None,
) -> None:
    checkpoint = {
        "adapter_head_state_dict": model.adapter_head_state_dict(),
        "n_tracks": n_tracks,
        "embedding_resolution": embedding_resolution,
        "head_type": head_type,
        "head_resolutions": list(head_resolutions),
        "linear_head_architecture": linear_head_architecture,
        "linear_hidden_channels": linear_hidden_channels,
        "linear_input_bottleneck_channels": linear_input_bottleneck_channels,
        "linear_residual_scale_init": linear_residual_scale_init,
        "linear_dilation": linear_dilation,
        "linear_kernel_size": linear_kernel_size,
        "residual_base_checkpoint": residual_base_checkpoint,
        "linear_residual_fusion": linear_residual_fusion,
        "residual_correction_scale_init": residual_correction_scale_init,
        "residual_correction_l2": residual_correction_l2,
        "residual_base_gate": residual_base_gate,
        "residual_base_gate_center": residual_base_gate_center,
        "residual_base_gate_sharpness": residual_base_gate_sharpness,
        "residual_base_gate_floor": residual_base_gate_floor,
        "linear_output_calibration": linear_output_calibration,
        "calibration_gate_center": calibration_gate_center,
        "calibration_gate_sharpness": calibration_gate_sharpness,
        "train_only_calibration": train_only_calibration,
        "linear_loss_type": linear_loss_type,
        "linear_target_space": linear_target_space,
        "smooth_l1_beta": smooth_l1_beta,
        "hybrid_loss_alpha": hybrid_loss_alpha,
        "track_loss_weight_preset": track_loss_weight_preset,
        "track_loss_weights": track_loss_weights,
        "signal_loss_weight_preset": signal_loss_weight_preset,
        "signal_high_threshold": signal_high_threshold,
        "signal_top_quantile": signal_top_quantile,
        "region_loss_weight_preset": region_loss_weight_preset,
        "gtf": gtf,
        "promoter_radius_bp": promoter_radius_bp,
        "gene_aux_loss_weight": gene_aux_loss_weight,
        "gene_aux_min_bases": gene_aux_min_bases,
        "organism_index": organism_index,
        "target_transform": target_transform,
        "loss_space": loss_space,
        "loss_type": loss_type,
        "multinomial_num_segments": multinomial_num_segments,
        "positional_weight": positional_weight,
        "count_weight": count_weight,
        "mse_weight": mse_weight,
        "poisson_weight": poisson_weight,
        "base_weights": str(weights),
        "best_valid_loss": best_valid_loss,
        "best_valid_step": best_valid_step,
        "best_selection_metric": best_selection_metric,
        "best_selection_value": best_selection_value,
        "lora_enabled": lora_enabled,
        "lora_target_modules": lora_target_modules,
        "lora_rank": lora_rank,
        "lora_alpha": lora_alpha,
        "lora_state_dict": (
            lora_state_dict(model.base_model) if lora_enabled else {}
        ),
        "freeze_head": freeze_head,
        "init_adapter_checkpoint": init_adapter_checkpoint,
        "init_adapter_checkpoint_mode": init_adapter_checkpoint_mode,
        "track_means_source": track_means_source,
        "track_means_tsv": track_means_tsv,
        "track_means_max_examples": track_means_max_examples,
    }
    torch.save(checkpoint, path)


def main() -> None:
    args = parse_args()
    requested_seed = args.seed
    args.seed = resolve_seed(args.seed)
    args.seed_request = requested_seed
    set_seed(args.seed)
    is_distributed, rank, local_rank, world_size, device = setup_distributed(args.device)
    if args.grad_accum_steps < 1:
        raise ValueError("--grad-accum-steps must be >= 1")
    if args.warmup_steps < 0:
        raise ValueError("--warmup-steps must be >= 0")
    if args.freeze_head and not args.enable_last_block_lora:
        raise ValueError("--freeze-head requires --enable-last-block-lora")
    if args.linear_hidden_channels < 1:
        raise ValueError("--linear-hidden-channels must be >= 1")
    if (
        args.linear_input_bottleneck_channels is not None
        and args.linear_input_bottleneck_channels < 1
    ):
        raise ValueError("--linear-input-bottleneck-channels must be >= 1")
    if args.linear_dilation < 1:
        raise ValueError("--linear-dilation must be >= 1")
    if args.linear_kernel_size < 1:
        raise ValueError("--linear-kernel-size must be >= 1")
    if args.smooth_l1_beta <= 0.0:
        raise ValueError("--smooth-l1-beta must be > 0")
    if not 0.0 <= args.hybrid_loss_alpha <= 1.0:
        raise ValueError("--hybrid-loss-alpha must be between 0 and 1")
    if args.residual_correction_l2 < 0.0:
        raise ValueError("--residual-correction-l2 must be >= 0")
    if args.residual_base_gate_sharpness <= 0.0:
        raise ValueError("--residual-base-gate-sharpness must be > 0")
    if not 0.0 <= args.residual_base_gate_floor <= 1.0:
        raise ValueError("--residual-base-gate-floor must be between 0 and 1")
    if args.calibration_gate_sharpness <= 0.0:
        raise ValueError("--calibration-gate-sharpness must be > 0")
    if args.signal_high_threshold < 0.0:
        raise ValueError("--signal-high-threshold must be >= 0")
    if args.signal_top_quantile is not None and not 0.0 < args.signal_top_quantile < 1.0:
        raise ValueError("--signal-top-quantile must be between 0 and 1")
    if args.gene_aux_loss_weight < 0.0:
        raise ValueError("--gene-aux-loss-weight must be >= 0")
    if args.gene_aux_min_bases < 1:
        raise ValueError("--gene-aux-min-bases must be >= 1")
    step_decay_steps = parse_int_list(
        args.step_decay_steps,
        option_name="--step-decay-steps",
    )
    step_decay_learning_rates = parse_float_list(
        args.step_decay_learning_rates,
        option_name="--step-decay-learning-rates",
    )
    if args.lr_schedule == "step":
        if not step_decay_steps:
            raise ValueError("--lr-schedule step requires --step-decay-steps")
        if not step_decay_learning_rates:
            raise ValueError(
                "--lr-schedule step requires --step-decay-learning-rates"
            )
        if len(step_decay_learning_rates) == 1:
            step_decay_learning_rates = step_decay_learning_rates * len(step_decay_steps)
        if len(step_decay_learning_rates) != len(step_decay_steps):
            raise ValueError(
                "--step-decay-learning-rates must contain one value or one value "
                "per --step-decay-steps boundary"
            )
    elif step_decay_steps or step_decay_learning_rates:
        raise ValueError("step decay options require --lr-schedule step")
    if args.early_stopping_patience is not None and args.early_stopping_patience < 1:
        raise ValueError("--early-stopping-patience must be >= 1")
    if args.early_stopping_min_steps < 0:
        raise ValueError("--early-stopping-min-steps must be >= 0")
    head_resolutions = parse_resolutions(
        args.head_resolutions,
        fallback_resolution=args.embedding_resolution,
    )
    if args.head_type == "linear" and head_resolutions != (args.embedding_resolution,):
        raise ValueError("linear head requires --head-resolutions to match --embedding-resolution")
    if args.head_type == "genome-tracks" and args.target_transform != "none":
        raise ValueError(
            "--head-type genome-tracks requires --target-transform none so raw "
            "RNA-seq targets can be scaled into model space."
        )
    if args.head_type == "genome-tracks" and args.linear_target_space != "full-log1p":
        raise ValueError("--linear-target-space is only supported for linear heads")
    if args.head_type != "linear" and args.linear_output_calibration != "none":
        raise ValueError("--linear-output-calibration is only supported for linear heads")
    if args.train_only_calibration and args.linear_output_calibration == "none":
        raise ValueError("--train-only-calibration requires --linear-output-calibration")
    if (
        args.head_type != "linear"
        and (
            args.signal_loss_weight_preset != "none"
            or args.region_loss_weight_preset != "none"
            or args.gene_aux_loss_weight > 0.0
        )
    ):
        raise ValueError("Signal/region/gene auxiliary losses are only supported for linear heads")
    if args.residual_base_checkpoint is not None:
        if args.head_type != "linear":
            raise ValueError("--residual-base-checkpoint requires --head-type linear")
        if args.linear_target_space != "full-log1p":
            raise ValueError(
                "--residual-base-checkpoint currently requires full-log1p targets"
            )
    elif args.residual_correction_l2 > 0.0 or args.residual_base_gate != "none":
        raise ValueError(
            "--residual-correction-l2 and --residual-base-gate require "
            "--residual-base-checkpoint"
        )
    if (
        args.linear_residual_fusion != "none"
        and args.residual_base_checkpoint is None
    ):
        raise ValueError("--linear-residual-fusion requires --residual-base-checkpoint")
    if args.loss_type in {"poisson-multinomial", "hybrid-mse-poisson"}:
        if args.head_type != "genome-tracks":
            raise ValueError(f"--loss-type {args.loss_type} requires genome-tracks")
        if args.target_transform != "none":
            raise ValueError(f"--loss-type {args.loss_type} requires raw targets")
        if args.multinomial_num_segments < 1:
            raise ValueError("--multinomial-num-segments must be >= 1")
    if args.loss_type == "hybrid-mse-poisson":
        if args.mse_weight < 0.0 or args.poisson_weight < 0.0:
            raise ValueError("Hybrid loss weights must be non-negative")
        if args.mse_weight == 0.0 and args.poisson_weight == 0.0:
            raise ValueError("At least one hybrid loss weight must be positive")
    dataset_target_transform = (
        "none"
        if args.head_type == "linear"
        and args.linear_target_space == "binned128-log1p-mean"
        else args.target_transform
    )

    train_dataset_dir = Path(args.train_dataset_dir)
    valid_dataset_dir = Path(args.valid_dataset_dir)
    weights = Path(args.weights)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    config_path = output_dir / "config.json"
    metrics_path = output_dir / "metrics.tsv"
    checkpoint_path = output_dir / "adapter_head.pt"
    best_checkpoint_path = output_dir / "adapter_head_best.pt"

    train_sampler = None
    if is_distributed:
        train_dataset_for_sampler = make_dataloader(
            train_dataset_dir,
            batch_size=args.batch_size,
            shuffle=False,
            num_workers=0,
            target_transform=dataset_target_transform,
            max_examples=args.max_train_examples,
        ).dataset
        train_sampler = DistributedSampler(
            train_dataset_for_sampler,
            num_replicas=world_size,
            rank=rank,
            shuffle=True,
            seed=args.seed,
            drop_last=False,
        )
        train_loader = torch.utils.data.DataLoader(
            train_dataset_for_sampler,
            batch_size=args.batch_size,
            sampler=train_sampler,
            num_workers=args.num_workers,
            pin_memory=torch.cuda.is_available(),
        )
    else:
        train_loader = make_dataloader(
            train_dataset_dir,
            batch_size=args.batch_size,
            shuffle=True,
            num_workers=args.num_workers,
            target_transform=dataset_target_transform,
            max_examples=args.max_train_examples,
        )
    valid_loader = make_dataloader(
        valid_dataset_dir,
        batch_size=args.batch_size,
        shuffle=False,
        num_workers=args.num_workers,
        target_transform="none",
        max_examples=args.max_valid_examples,
    )
    n_tracks = int(train_loader.dataset.metadata["n_tracks"])
    track_loss_weights = resolve_track_loss_weights(
        weights_text=args.track_loss_weights,
        preset=args.track_loss_weight_preset,
        n_tracks=n_tracks,
    )
    track_loss_weights_tensor = (
        None
        if track_loss_weights is None else torch.tensor(track_loss_weights, dtype=torch.float32)
    )
    if track_loss_weights_tensor is not None and args.head_type != "linear":
        raise ValueError("--track-loss-weights are only supported for linear heads")

    gtf_annotations = None
    if args.region_loss_weight_preset != "none" or args.gene_aux_loss_weight > 0.0:
        gtf_path = Path(args.gtf)
        if not gtf_path.exists():
            raise FileNotFoundError(f"GTF not found: {gtf_path}")
        gtf_annotations = load_gtf_annotations_for_training(
            gtf_path,
            promoter_radius_bp=args.promoter_radius_bp,
        )

    track_means = None
    if args.head_type == "genome-tracks":
        if args.track_means_source == "ones":
            track_means = torch.ones(n_tracks, dtype=torch.float32)
        elif args.track_means_source == "grouped-qc":
            track_means = load_track_means_from_grouped_qc(
                args.track_means_tsv,
                n_tracks=n_tracks,
            )
        elif args.track_means_source == "train-nonzero":
            track_mean_dataset = AlphaGenomeRnaSeqNpzDataset(
                train_dataset_dir,
                target_transform="none",
                max_examples=args.track_means_max_examples,
            )
            track_means = compute_nonzero_track_means_from_dataset(track_mean_dataset)
        else:
            raise ValueError(f"Unknown track means source: {args.track_means_source}")

    linear_track_means = None
    if (
        args.head_type == "linear"
        and is_residual_linear_architecture(args.linear_head_architecture)
    ):
        linear_mean_dataset = AlphaGenomeRnaSeqNpzDataset(
            train_dataset_dir,
            target_transform=dataset_target_transform,
            max_examples=args.max_train_examples,
        )
        linear_track_means = compute_linear_track_means_from_dataset(
            linear_mean_dataset,
            target_space=args.linear_target_space,
        )

    if is_main_process(rank):
        write_config(
            config_path,
            args=args,
            device=device,
            n_train_examples=len(train_loader.dataset),
            n_valid_examples=len(valid_loader.dataset),
            n_tracks=n_tracks,
        )

    print_main(rank, f"train_dataset_dir\t{train_dataset_dir}")
    print_main(rank, f"valid_dataset_dir\t{valid_dataset_dir}")
    print_main(rank, f"weights\t{weights}")
    print_main(rank, f"output_dir\t{output_dir}")
    print_main(rank, f"config_path\t{config_path}")
    print_main(rank, f"metrics_path\t{metrics_path}")
    print_main(rank, f"init_adapter_checkpoint\t{args.init_adapter_checkpoint}")
    print_main(
        rank,
        f"init_adapter_checkpoint_mode\t{args.init_adapter_checkpoint_mode}",
    )
    print_main(
        rank,
        f"checkpoint_path\t{checkpoint_path if not args.no_save_checkpoint else 'disabled'}",
    )
    print_main(
        rank,
        "best_checkpoint_path\t"
        f"{best_checkpoint_path if not args.no_save_checkpoint else 'disabled'}",
    )
    print_main(rank, f"n_train_examples\t{len(train_loader.dataset)}")
    print_main(rank, f"n_valid_examples\t{len(valid_loader.dataset)}")
    print_main(rank, f"n_tracks\t{n_tracks}")
    print_main(rank, f"target_transform\t{args.target_transform}")
    print_main(rank, f"dataset_target_transform\t{dataset_target_transform}")
    print_main(rank, "valid_dataset_target_transform\tnone")
    print_main(rank, f"head_type\t{args.head_type}")
    print_main(rank, f"head_resolutions\t{','.join(map(str, head_resolutions))}")
    print_main(rank, f"linear_head_architecture\t{args.linear_head_architecture}")
    print_main(rank, f"linear_hidden_channels\t{args.linear_hidden_channels}")
    print_main(
        rank,
        "linear_input_bottleneck_channels\t"
        f"{args.linear_input_bottleneck_channels}",
    )
    print_main(rank, f"residual_scale_init\t{args.residual_scale_init}")
    print_main(rank, f"linear_dilation\t{args.linear_dilation}")
    print_main(rank, f"linear_kernel_size\t{args.linear_kernel_size}")
    print_main(rank, f"residual_base_checkpoint\t{args.residual_base_checkpoint}")
    print_main(rank, f"linear_residual_fusion\t{args.linear_residual_fusion}")
    print_main(
        rank,
        f"residual_correction_scale_init\t{args.residual_correction_scale_init}",
    )
    print_main(rank, f"residual_correction_l2\t{args.residual_correction_l2}")
    print_main(rank, f"residual_base_gate\t{args.residual_base_gate}")
    print_main(
        rank,
        f"residual_base_gate_center\t{args.residual_base_gate_center}",
    )
    print_main(
        rank,
        f"residual_base_gate_sharpness\t{args.residual_base_gate_sharpness}",
    )
    print_main(rank, f"residual_base_gate_floor\t{args.residual_base_gate_floor}")
    print_main(rank, f"linear_output_calibration\t{args.linear_output_calibration}")
    print_main(rank, f"calibration_gate_center\t{args.calibration_gate_center}")
    print_main(rank, f"calibration_gate_sharpness\t{args.calibration_gate_sharpness}")
    print_main(rank, f"train_only_calibration\t{args.train_only_calibration}")
    print_main(rank, f"linear_loss_type\t{args.linear_loss_type}")
    print_main(rank, f"smooth_l1_beta\t{args.smooth_l1_beta}")
    print_main(rank, f"hybrid_loss_alpha\t{args.hybrid_loss_alpha}")
    print_main(rank, f"track_loss_weight_preset\t{args.track_loss_weight_preset}")
    print_main(
        rank,
        f"track_loss_weights\t{format_track_loss_weights(track_loss_weights)}",
    )
    print_main(rank, f"signal_loss_weight_preset\t{args.signal_loss_weight_preset}")
    print_main(rank, f"signal_high_threshold\t{args.signal_high_threshold}")
    print_main(rank, f"signal_top_quantile\t{args.signal_top_quantile}")
    print_main(rank, f"region_loss_weight_preset\t{args.region_loss_weight_preset}")
    print_main(rank, f"gtf\t{args.gtf}")
    print_main(rank, f"promoter_radius_bp\t{args.promoter_radius_bp}")
    print_main(rank, f"gene_aux_loss_weight\t{args.gene_aux_loss_weight}")
    print_main(rank, f"gene_aux_min_bases\t{args.gene_aux_min_bases}")
    print_main(rank, f"linear_target_space\t{args.linear_target_space}")
    if linear_track_means is not None:
        print_main(
            rank,
            "linear_track_means\t"
            + ",".join(f"{float(value):.8g}" for value in linear_track_means),
        )
    print_main(
        rank,
        "loss_space\t"
        f"{loss_space_for_config(args.head_type, args.loss_type, args.linear_target_space, args.linear_loss_type)}",
    )
    print_main(rank, f"loss_type\t{args.loss_type}")
    print_main(rank, f"multinomial_num_segments\t{args.multinomial_num_segments}")
    print_main(rank, f"positional_weight\t{args.positional_weight}")
    print_main(rank, f"count_weight\t{args.count_weight}")
    print_main(rank, f"mse_weight\t{args.mse_weight}")
    print_main(rank, f"poisson_weight\t{args.poisson_weight}")
    print_main(rank, f"track_means_source\t{args.track_means_source}")
    if track_means is not None:
        print_main(
            rank,
            "track_means\t"
            + ",".join(f"{float(value):.8g}" for value in track_means),
        )
    print_main(rank, f"embedding_resolution\t{args.embedding_resolution}")
    print_main(rank, f"organism_index\t{args.organism_index}")
    print_main(rank, f"learning_rate\t{args.learning_rate}")
    print_main(rank, f"weight_decay\t{args.weight_decay}")
    print_main(rank, f"grad_accum_steps\t{args.grad_accum_steps}")
    print_main(rank, f"grad_clip_norm\t{args.grad_clip_norm}")
    print_main(rank, f"warmup_steps\t{args.warmup_steps}")
    print_main(rank, f"lr_schedule\t{args.lr_schedule}")
    print_main(rank, f"step_decay_steps\t{','.join(map(str, step_decay_steps))}")
    print_main(
        rank,
        "step_decay_learning_rates\t"
        f"{','.join(f'{value:.8g}' for value in step_decay_learning_rates)}",
    )
    print_main(rank, f"selection_metric\t{args.selection_metric}")
    print_main(rank, f"init_adapter_checkpoint\t{args.init_adapter_checkpoint}")
    print_main(rank, f"enable_last_block_lora\t{args.enable_last_block_lora}")
    print_main(rank, f"lora_rank\t{args.lora_rank}")
    print_main(rank, f"lora_alpha\t{args.lora_alpha}")
    print_main(rank, f"freeze_head\t{args.freeze_head}")
    print_main(rank, f"per_process_batch_size\t{args.batch_size}")
    print_main(
        rank,
        f"global_effective_batch_size\t{args.batch_size * args.grad_accum_steps * world_size}",
    )
    print_main(rank, f"seed_request\t{args.seed_request}")
    print_main(rank, f"seed\t{args.seed}")
    print_main(rank, f"device\t{device}")
    print_main(rank, f"distributed\t{is_distributed}")
    print_main(rank, f"rank\t{rank}")
    print_main(rank, f"local_rank\t{local_rank}")
    print_main(rank, f"world_size\t{world_size}")
    if device.type == "cuda":
        print_main(rank, f"cuda_device\t{torch.cuda.get_device_name(device)}")

    base_model = AlphaGenome.from_pretrained(weights, device=device)
    residual_base_head = None
    residual_base_input_bottleneck = None
    residual_base_resolution = 128
    if args.residual_base_checkpoint is not None:
        residual_base_checkpoint = load_adapter_checkpoint(args.residual_base_checkpoint)
        residual_base_n_tracks = int(residual_base_checkpoint["n_tracks"])
        residual_base_head_type = str(
            residual_base_checkpoint.get("head_type", "linear")
        )
        residual_base_target_space = str(
            residual_base_checkpoint.get("linear_target_space", "full-log1p")
        )
        residual_base_resolution = int(residual_base_checkpoint["embedding_resolution"])
        if residual_base_n_tracks != n_tracks:
            raise ValueError(
                "Residual base checkpoint n_tracks="
                f"{residual_base_n_tracks}, dataset n_tracks={n_tracks}"
            )
        if residual_base_head_type != "linear":
            raise ValueError("--residual-base-checkpoint must be a linear head")
        if residual_base_target_space != "full-log1p":
            raise ValueError(
                "--residual-base-checkpoint currently must use full-log1p"
            )
        residual_base_probe = RnaSeq11Adapter(
            base_model,
            n_tracks=n_tracks,
            embedding_resolution=residual_base_resolution,
            organism_index=args.organism_index,
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
            linear_kernel_size=int(
                residual_base_checkpoint.get("linear_kernel_size", 15)
            ),
            linear_output_calibration=str(
                residual_base_checkpoint.get("linear_output_calibration", "none")
            ),
            linear_output_calibration_gate_center=float(
                residual_base_checkpoint.get("calibration_gate_center", 3.0)
            ),
            linear_output_calibration_gate_sharpness=float(
                residual_base_checkpoint.get("calibration_gate_sharpness", 3.0)
            ),
        ).to(device)
        residual_base_probe.load_adapter_head_state_dict(
            residual_base_checkpoint["adapter_head_state_dict"]
        )
        residual_base_head = residual_base_probe.head
        residual_base_input_bottleneck = residual_base_probe.linear_input_bottleneck

    model = RnaSeq11Adapter(
        base_model,
        n_tracks=n_tracks,
        embedding_resolution=args.embedding_resolution,
        organism_index=args.organism_index,
        encode_requires_grad=args.enable_last_block_lora,
        head_type=args.head_type,
        head_resolutions=head_resolutions,
        linear_head_architecture=args.linear_head_architecture,
        linear_hidden_channels=args.linear_hidden_channels,
        linear_input_bottleneck_channels=args.linear_input_bottleneck_channels,
        linear_track_means=linear_track_means,
        linear_residual_scale_init=args.residual_scale_init,
        linear_dilation=args.linear_dilation,
        linear_kernel_size=args.linear_kernel_size,
        linear_residual_base_head=residual_base_head,
        linear_residual_base_input_bottleneck=residual_base_input_bottleneck,
        linear_residual_base_resolution=residual_base_resolution,
        linear_residual_fusion=args.linear_residual_fusion,
        linear_residual_correction_scale_init=args.residual_correction_scale_init,
        linear_residual_base_gate=args.residual_base_gate,
        linear_residual_base_gate_center=args.residual_base_gate_center,
        linear_residual_base_gate_sharpness=args.residual_base_gate_sharpness,
        linear_residual_base_gate_floor=args.residual_base_gate_floor,
        linear_output_calibration=args.linear_output_calibration,
        linear_output_calibration_gate_center=args.calibration_gate_center,
        linear_output_calibration_gate_sharpness=args.calibration_gate_sharpness,
        track_means=track_means,
    ).to(device)

    if args.init_adapter_checkpoint is not None:
        init_checkpoint = load_adapter_checkpoint(args.init_adapter_checkpoint)
        checkpoint_n_tracks = int(init_checkpoint["n_tracks"])
        checkpoint_resolution = int(init_checkpoint["embedding_resolution"])
        checkpoint_head_type = str(init_checkpoint.get("head_type", "linear"))
        checkpoint_head_resolutions = tuple(
            int(resolution)
            for resolution in init_checkpoint.get(
                "head_resolutions",
                [checkpoint_resolution],
            )
        )
        checkpoint_linear_head_architecture = str(
            init_checkpoint.get("linear_head_architecture", "conv1x1")
        )
        checkpoint_linear_hidden_channels = int(
            init_checkpoint.get("linear_hidden_channels", args.linear_hidden_channels)
        )
        checkpoint_linear_input_bottleneck_channels = init_checkpoint.get(
            "linear_input_bottleneck_channels"
        )
        checkpoint_linear_dilation = int(
            init_checkpoint.get("linear_dilation", args.linear_dilation)
        )
        checkpoint_linear_kernel_size = int(
            init_checkpoint.get("linear_kernel_size", args.linear_kernel_size)
        )
        checkpoint_linear_residual_fusion = str(
            init_checkpoint.get("linear_residual_fusion", "none")
        )
        checkpoint_linear_target_space = str(
            init_checkpoint.get("linear_target_space", "full-log1p")
        )
        checkpoint_linear_output_calibration = str(
            init_checkpoint.get("linear_output_calibration", "none")
        )
        if checkpoint_n_tracks != n_tracks:
            raise ValueError(
                f"Init checkpoint n_tracks={checkpoint_n_tracks}, dataset n_tracks={n_tracks}"
            )
        if checkpoint_head_type != args.head_type:
            raise ValueError(
                f"Init checkpoint head_type={checkpoint_head_type}, "
                f"requested={args.head_type}"
            )
        if checkpoint_head_resolutions != head_resolutions:
            raise ValueError(
                "Init checkpoint head_resolutions="
                f"{checkpoint_head_resolutions}, requested={head_resolutions}"
            )
        if (
            args.head_type == "linear"
            and checkpoint_resolution != args.embedding_resolution
        ):
            raise ValueError(
                "Init checkpoint embedding_resolution="
                f"{checkpoint_resolution}, requested={args.embedding_resolution}"
            )
        if args.head_type == "linear":
            if checkpoint_linear_target_space != args.linear_target_space:
                raise ValueError(
                    "Init checkpoint linear_target_space="
                    f"{checkpoint_linear_target_space}, requested={args.linear_target_space}"
                )
        if args.head_type == "linear" and args.init_adapter_checkpoint_mode == "strict":
            if checkpoint_linear_head_architecture != args.linear_head_architecture:
                raise ValueError(
                    "Init checkpoint linear_head_architecture="
                    f"{checkpoint_linear_head_architecture}, "
                    f"requested={args.linear_head_architecture}"
                )
            if checkpoint_linear_head_architecture != "conv1x1" and (
                checkpoint_linear_hidden_channels != args.linear_hidden_channels
            ):
                raise ValueError(
                    "Init checkpoint linear_hidden_channels="
                    f"{checkpoint_linear_hidden_channels}, "
                    f"requested={args.linear_hidden_channels}"
                )
            if (
                checkpoint_linear_input_bottleneck_channels
                != args.linear_input_bottleneck_channels
            ):
                raise ValueError(
                    "Init checkpoint linear_input_bottleneck_channels="
                    f"{checkpoint_linear_input_bottleneck_channels}, "
                    "requested="
                    f"{args.linear_input_bottleneck_channels}"
                )
            if (
                checkpoint_linear_head_architecture == "dilated-conv3"
                and checkpoint_linear_dilation != args.linear_dilation
            ):
                raise ValueError(
                    "Init checkpoint linear_dilation="
                    f"{checkpoint_linear_dilation}, requested={args.linear_dilation}"
                )
            if (
                checkpoint_linear_head_architecture == "depthwise-separable-conv"
                and checkpoint_linear_kernel_size != args.linear_kernel_size
            ):
                raise ValueError(
                    "Init checkpoint linear_kernel_size="
                    f"{checkpoint_linear_kernel_size}, "
                    f"requested={args.linear_kernel_size}"
                )
            if checkpoint_linear_residual_fusion != args.linear_residual_fusion:
                raise ValueError(
                    "Init checkpoint linear_residual_fusion="
                    f"{checkpoint_linear_residual_fusion}, "
                    f"requested={args.linear_residual_fusion}"
                )
            if (
                checkpoint_linear_output_calibration != args.linear_output_calibration
                and not (
                    checkpoint_linear_output_calibration == "none"
                    and args.linear_output_calibration != "none"
                )
            ):
                raise ValueError(
                    "Init checkpoint linear_output_calibration="
                    f"{checkpoint_linear_output_calibration}, "
                    f"requested={args.linear_output_calibration}"
                )
        if args.init_adapter_checkpoint_mode == "strict":
            model.load_adapter_head_state_dict(
                init_checkpoint["adapter_head_state_dict"]
            )
        else:
            loaded_keys, skipped_keys = load_matching_adapter_head_state_dict(
                model,
                init_checkpoint["adapter_head_state_dict"],
            )
            print_main(rank, f"init_matching_loaded_keys\t{len(loaded_keys)}")
            print_main(rank, f"init_matching_skipped_keys\t{len(skipped_keys)}")
            if loaded_keys:
                print_main(
                    rank,
                    "init_matching_loaded_key_list\t" + ",".join(loaded_keys),
                )
            if skipped_keys:
                print_main(
                    rank,
                    "init_matching_skipped_key_list\t" + ",".join(skipped_keys),
                )

    lora_target_modules: list[str] = []
    lora_parameters: list[torch.nn.Parameter] = []
    lora_parameter_ids: set[int] = set()
    lora_module_names: list[str] = []
    if args.enable_last_block_lora:
        lora_target_modules = ["tower.blocks.8.mha", "tower.blocks.8.mlp"]
        apply_lora(
            model.base_model,
            target_modules=lora_target_modules,
            rank=args.lora_rank,
            alpha=args.lora_alpha,
        )
        model.to(device)
        lora_module_names = [
            name for name, module in model.base_model.named_modules()
            if isinstance(module, LoRA)
        ]
        if not lora_module_names:
            raise RuntimeError(
                f"No LoRA modules were applied for targets: {lora_target_modules}"
            )
        lora_parameters = get_adapter_params(model.base_model)
        lora_parameter_ids = {id(parameter) for parameter in lora_parameters}

    if args.freeze_head:
        for parameter in model.head.parameters():
            parameter.requires_grad = False

    if args.train_only_calibration:
        if model.linear_output_calibration is None:
            raise ValueError("--train-only-calibration requires a calibration module")
        for parameter in model.parameters():
            parameter.requires_grad = False
        for parameter in model.linear_output_calibration.parameters():
            parameter.requires_grad = True

    trainable_parameters = [
        parameter for parameter in model.parameters() if parameter.requires_grad
    ]
    if not trainable_parameters:
        raise ValueError("No trainable parameters selected")

    train_model: torch.nn.Module = model
    if is_distributed:
        train_model = DistributedDataParallel(
            model,
            device_ids=[local_rank],
            output_device=local_rank,
            find_unused_parameters=False,
        )
    optimizer = torch.optim.AdamW(
        trainable_parameters,
        lr=args.learning_rate,
        weight_decay=args.weight_decay,
    )

    base_non_lora_trainable_parameters = sum(
        parameter.numel()
        for parameter in model.base_model.parameters()
        if parameter.requires_grad and id(parameter) not in lora_parameter_ids
    )
    print_main(rank, f"base_parameters\t{count_parameters(model.base_model)}")
    print_main(rank, f"total_parameters\t{count_parameters(model)}")
    print_main(
        rank,
        f"trainable_parameters\t{count_parameters(model, trainable_only=True)}",
    )
    print_main(rank, f"head_parameters\t{count_parameters(model.head)}")
    print_main(
        rank,
        f"head_trainable_parameters\t"
        f"{count_parameters(model.head, trainable_only=True)}",
    )
    print_main(
        rank,
        "calibration_trainable_parameters\t"
        f"{count_parameters(model.linear_output_calibration, trainable_only=True) if model.linear_output_calibration is not None else 0}",
    )
    print_main(
        rank,
        f"lora_trainable_parameters\t{sum(p.numel() for p in lora_parameters)}",
    )
    print_main(
        rank,
        f"base_non_lora_trainable_parameters\t{base_non_lora_trainable_parameters}",
    )
    print_main(
        rank,
        f"optimizer_trainable_parameters\t"
        f"{sum(p.numel() for p in trainable_parameters)}",
    )
    if args.enable_last_block_lora:
        print_main(rank, f"lora_target_modules\t{','.join(lora_target_modules)}")
        print_main(rank, f"lora_modules_applied\t{len(lora_module_names)}")
        for name in lora_module_names:
            print_main(rank, f"lora_module\t{name}")

    train_model.train()
    model.base_model.eval()
    best_valid_loss: float | None = None
    best_selection_value: float | None = None
    best_valid_step: int | None = None
    validations_without_improvement = 0
    train_batches = infinite_batches(train_loader, train_sampler)
    for step in range(1, args.max_steps + 1):
        lr = learning_rate_for_step(
            step=step,
            base_lr=args.learning_rate,
            max_steps=args.max_steps,
            warmup_steps=args.warmup_steps,
            schedule=args.lr_schedule,
            step_decay_steps=step_decay_steps,
            step_decay_learning_rates=step_decay_learning_rates,
        )
        for group in optimizer.param_groups:
            group["lr"] = lr
        optimizer.zero_grad(set_to_none=True)
        accumulated_train_loss = 0.0
        accumulated_examples = 0
        last_dna_shape: tuple[int, ...] | None = None
        last_rna_seq_shape: tuple[int, ...] | None = None
        last_prediction_shape: str | None = None

        for _ in range(args.grad_accum_steps):
            batch = next(train_batches)
            dna_sequence = batch["dna_sequence"].to(device=device, dtype=torch.float32)
            rna_seq = batch["rna_seq"].to(device=device, dtype=torch.float32)
            rna_seq_mask = batch["rna_seq_mask"].to(device=device)
            loss_target = (
                prepare_linear_target(rna_seq, target_space=args.linear_target_space)
                if args.head_type == "linear" else rna_seq
            )

            prediction_result = train_model(
                dna_sequence,
                target_length=loss_target.shape[-1],
                return_scaled=prediction_return_scaled_for_loss(model, args.loss_type),
                return_linear_residual_update=args.residual_correction_l2 > 0.0,
            )
            residual_update = None
            if isinstance(prediction_result, tuple):
                prediction, residual_update = prediction_result
            else:
                prediction = prediction_result
            signal_loss_weights = build_signal_loss_weights(
                target=loss_target,
                preset=args.signal_loss_weight_preset,
                high_threshold=args.signal_high_threshold,
                top_quantile=args.signal_top_quantile,
            )
            region_loss_weights = None
            if args.region_loss_weight_preset != "none":
                if gtf_annotations is None:
                    raise RuntimeError("GTF annotations were not loaded")
                region_loss_weights = build_region_loss_weights(
                    target=loss_target,
                    batch=batch,
                    annotations=gtf_annotations,
                    preset=args.region_loss_weight_preset,
                )
            element_loss_weights = combine_element_loss_weights(
                signal_loss_weights,
                region_loss_weights,
            )
            loss = adapter_loss_for_training(
                model,
                prediction,
                loss_target if args.head_type == "linear" else rna_seq,
                rna_seq_mask,
                loss_type=args.loss_type,
                multinomial_num_segments=args.multinomial_num_segments,
                positional_weight=args.positional_weight,
                count_weight=args.count_weight,
                mse_weight=args.mse_weight,
                poisson_weight=args.poisson_weight,
                linear_loss_type=args.linear_loss_type,
                smooth_l1_beta=args.smooth_l1_beta,
                hybrid_loss_alpha=args.hybrid_loss_alpha,
                track_loss_weights=track_loss_weights_tensor,
                loss_weights=element_loss_weights,
            )
            if args.gene_aux_loss_weight > 0.0:
                if gtf_annotations is None:
                    raise RuntimeError("GTF annotations were not loaded")
                if not isinstance(prediction, torch.Tensor):
                    raise RuntimeError("Gene auxiliary loss requires a linear prediction tensor")
                loss = loss + args.gene_aux_loss_weight * gene_body_mean_aux_loss(
                    prediction,
                    loss_target,
                    rna_seq_mask,
                    batch=batch,
                    annotations=gtf_annotations,
                    min_bases=args.gene_aux_min_bases,
                )
            if args.residual_correction_l2 > 0.0:
                if residual_update is None:
                    raise RuntimeError(
                        "--residual-correction-l2 requires residual update output"
                    )
                loss = loss + args.residual_correction_l2 * masked_residual_update_l2(
                    residual_update,
                    rna_seq_mask,
                )
            (loss / args.grad_accum_steps).backward()

            accumulated_train_loss += float(loss.detach().cpu())
            accumulated_examples += int(dna_sequence.shape[0])
            last_dna_shape = tuple(dna_sequence.shape)
            last_rna_seq_shape = tuple(rna_seq.shape)
            last_prediction_shape = format_prediction_shape(prediction)

        if args.grad_clip_norm is not None:
            torch.nn.utils.clip_grad_norm_(trainable_parameters, args.grad_clip_norm)
        optimizer.step()

        train_loss = accumulated_train_loss / args.grad_accum_steps
        train_loss = reduce_mean(train_loss, device, world_size)
        if args.enable_last_block_lora:
            base_has_grad = any(
                p.grad is not None for p in model.base_model.parameters()
                if id(p) not in lora_parameter_ids
            )
        else:
            base_has_grad = any(p.grad is not None for p in model.base_model.parameters())
        base_has_grad = bool(int(reduce_mean(float(base_has_grad), device, world_size)))
        head_weight_grad_norm = grad_norm(list(model.head.parameters()))
        head_weight_grad_norm = reduce_mean(head_weight_grad_norm, device, world_size)
        lora_has_grad = any(parameter.grad is not None for parameter in lora_parameters)
        lora_has_grad = bool(int(reduce_mean(float(lora_has_grad), device, world_size)))
        lora_weight_grad_norm = reduce_mean(
            grad_norm(lora_parameters),
            device,
            world_size,
        )
        cuda_memory = (
            torch.cuda.max_memory_allocated() / 1024**2
            if device.type == "cuda" else None
        )
        if is_main_process(rank):
            append_metric(
                metrics_path,
                split="train",
                step=step,
                loss=train_loss,
                batches=args.grad_accum_steps * world_size,
                examples=accumulated_examples * world_size,
                lr=lr,
                cuda_max_memory_allocated_mb=cuda_memory,
            )

            print(f"step\t{step}")
            print(f"train_loss\t{train_loss:.6g}")
            if last_dna_shape is not None:
                print(f"dna_sequence_shape\t{'x'.join(map(str, last_dna_shape))}")
            if last_rna_seq_shape is not None:
                print(f"rna_seq_shape\t{'x'.join(map(str, last_rna_seq_shape))}")
            if last_prediction_shape is not None:
                print(f"prediction_shape\t{last_prediction_shape}")
            print(f"base_has_grad\t{base_has_grad}")
            if args.enable_last_block_lora:
                print(f"base_non_lora_has_grad\t{base_has_grad}")
                print(f"lora_has_grad\t{lora_has_grad}")
                print(f"lora_grad_norm\t{lora_weight_grad_norm:.6g}")
            print(f"head_weight_grad_norm\t{head_weight_grad_norm:.6g}")
            if cuda_memory is not None:
                print(f"cuda_max_memory_allocated_mb\t{cuda_memory:.1f}")

        if step % args.eval_every == 0 or step == args.max_steps:
            stop_training = False
            if is_main_process(rank):
                valid_metrics = evaluate_validation_metrics(
                    model,
                    valid_loader,
                    device=device,
                    loss_type=args.loss_type,
                    multinomial_num_segments=args.multinomial_num_segments,
                    positional_weight=args.positional_weight,
                    count_weight=args.count_weight,
                    mse_weight=args.mse_weight,
                    poisson_weight=args.poisson_weight,
                    linear_loss_type=args.linear_loss_type,
                    linear_target_space=args.linear_target_space,
                    target_transform=args.target_transform,
                    smooth_l1_beta=args.smooth_l1_beta,
                    hybrid_loss_alpha=args.hybrid_loss_alpha,
                )
                valid_loss = float(valid_metrics["loss"])
                valid_selection_value = selection_value(
                    valid_metrics,
                    args.selection_metric,
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
                    full_mse=float(valid_metrics["full_mse"]),
                    full_mae=float(valid_metrics["full_mae"]),
                    pearson=float(valid_metrics["pearson"]),
                    common128_mse=float(valid_metrics["common128_mse"]),
                    common128_mae=float(valid_metrics["common128_mae"]),
                    selection_metric=args.selection_metric,
                    selection_value_for_row=valid_selection_value,
                    batches=int(valid_metrics["batches"]),
                    examples=int(valid_metrics["examples"]),
                    lr=lr,
                    cuda_max_memory_allocated_mb=cuda_memory,
                )
                print(f"valid_step\t{step}")
                print(f"valid_batches\t{valid_metrics['batches']}")
                print(f"valid_examples\t{valid_metrics['examples']}")
                print(f"valid_loss\t{valid_loss:.6g}")
                print(f"valid_full_mse\t{float(valid_metrics['full_mse']):.6g}")
                print(f"valid_full_mae\t{float(valid_metrics['full_mae']):.6g}")
                print(f"valid_pearson\t{float(valid_metrics['pearson']):.6g}")
                print(f"valid_common128_mse\t{float(valid_metrics['common128_mse']):.6g}")
                print(f"valid_common128_mae\t{float(valid_metrics['common128_mae']):.6g}")
                print(f"selection_metric\t{args.selection_metric}")
                print(f"selection_value\t{valid_selection_value:.6g}")
                improved = selection_improved(
                    valid_selection_value,
                    best_selection_value,
                    selection_metric=args.selection_metric,
                )
                if improved:
                    best_valid_loss = valid_loss
                    best_selection_value = valid_selection_value
                    best_valid_step = step
                    validations_without_improvement = 0
                else:
                    validations_without_improvement += 1
                if not args.no_save_checkpoint and improved:
                    save_training_checkpoint(
                        best_checkpoint_path,
                        model=model,
                        n_tracks=n_tracks,
                        embedding_resolution=args.embedding_resolution,
                        organism_index=args.organism_index,
                        target_transform=args.target_transform,
                        weights=weights,
                        best_valid_loss=best_valid_loss,
                        best_valid_step=best_valid_step,
                        best_selection_metric=args.selection_metric,
                        best_selection_value=best_selection_value,
                        lora_enabled=args.enable_last_block_lora,
                        lora_target_modules=lora_target_modules,
                        lora_rank=args.lora_rank,
                        lora_alpha=args.lora_alpha,
                        freeze_head=args.freeze_head,
                        init_adapter_checkpoint=args.init_adapter_checkpoint,
                        init_adapter_checkpoint_mode=args.init_adapter_checkpoint_mode,
                        head_type=args.head_type,
                        head_resolutions=head_resolutions,
                        linear_head_architecture=args.linear_head_architecture,
                        linear_hidden_channels=args.linear_hidden_channels,
                        linear_input_bottleneck_channels=(
                            args.linear_input_bottleneck_channels
                        ),
                        linear_residual_scale_init=args.residual_scale_init,
                        linear_dilation=args.linear_dilation,
                        linear_kernel_size=args.linear_kernel_size,
                        residual_base_checkpoint=args.residual_base_checkpoint,
                        linear_residual_fusion=args.linear_residual_fusion,
                        residual_correction_scale_init=(
                            args.residual_correction_scale_init
                        ),
                        residual_correction_l2=args.residual_correction_l2,
                        residual_base_gate=args.residual_base_gate,
                        residual_base_gate_center=args.residual_base_gate_center,
                        residual_base_gate_sharpness=(
                            args.residual_base_gate_sharpness
                        ),
                        residual_base_gate_floor=args.residual_base_gate_floor,
                        linear_output_calibration=args.linear_output_calibration,
                        calibration_gate_center=args.calibration_gate_center,
                        calibration_gate_sharpness=args.calibration_gate_sharpness,
                        train_only_calibration=args.train_only_calibration,
                        linear_loss_type=args.linear_loss_type,
                        linear_target_space=args.linear_target_space,
                        smooth_l1_beta=args.smooth_l1_beta,
                        hybrid_loss_alpha=args.hybrid_loss_alpha,
                        track_loss_weight_preset=args.track_loss_weight_preset,
                        track_loss_weights=track_loss_weights,
                        signal_loss_weight_preset=args.signal_loss_weight_preset,
                        signal_high_threshold=args.signal_high_threshold,
                        signal_top_quantile=args.signal_top_quantile,
                        region_loss_weight_preset=args.region_loss_weight_preset,
                        gtf=args.gtf,
                        promoter_radius_bp=args.promoter_radius_bp,
                        gene_aux_loss_weight=args.gene_aux_loss_weight,
                        gene_aux_min_bases=args.gene_aux_min_bases,
                        loss_space=loss_space_for_config(
                            args.head_type,
                            args.loss_type,
                            args.linear_target_space,
                            args.linear_loss_type,
                        ),
                        loss_type=args.loss_type,
                        multinomial_num_segments=args.multinomial_num_segments,
                        positional_weight=args.positional_weight,
                        count_weight=args.count_weight,
                        mse_weight=args.mse_weight,
                        poisson_weight=args.poisson_weight,
                        track_means_source=args.track_means_source,
                        track_means_tsv=args.track_means_tsv,
                        track_means_max_examples=args.track_means_max_examples,
                    )
                    print(f"saved_best_checkpoint\t{best_checkpoint_path}")
                    print(f"best_valid_step\t{best_valid_step}")
                    print(f"best_valid_loss\t{best_valid_loss:.6g}")
                    print(f"best_selection_value\t{best_selection_value:.6g}")
                if (
                    args.early_stopping_patience is not None
                    and step >= args.early_stopping_min_steps
                    and validations_without_improvement >= args.early_stopping_patience
                ):
                    stop_training = True
                    print(f"early_stopping_step\t{step}")
                    print(
                        "early_stopping_validations_without_improvement\t"
                        f"{validations_without_improvement}"
                    )
            if is_distributed:
                stop_tensor = torch.tensor(
                    int(stop_training),
                    dtype=torch.int64,
                    device=device,
                )
                dist.broadcast(stop_tensor, src=0)
                stop_training = bool(int(stop_tensor.item()))
            if is_distributed:
                dist.barrier()
            if stop_training:
                break

    if is_main_process(rank) and not args.no_save_checkpoint:
        save_training_checkpoint(
            checkpoint_path,
            model=model,
            n_tracks=n_tracks,
            embedding_resolution=args.embedding_resolution,
            organism_index=args.organism_index,
            target_transform=args.target_transform,
            weights=weights,
            best_valid_loss=best_valid_loss,
            best_valid_step=best_valid_step,
            best_selection_metric=args.selection_metric,
            best_selection_value=best_selection_value,
            lora_enabled=args.enable_last_block_lora,
            lora_target_modules=lora_target_modules,
            lora_rank=args.lora_rank,
            lora_alpha=args.lora_alpha,
            freeze_head=args.freeze_head,
            init_adapter_checkpoint=args.init_adapter_checkpoint,
            init_adapter_checkpoint_mode=args.init_adapter_checkpoint_mode,
            head_type=args.head_type,
            head_resolutions=head_resolutions,
            linear_head_architecture=args.linear_head_architecture,
            linear_hidden_channels=args.linear_hidden_channels,
            linear_input_bottleneck_channels=args.linear_input_bottleneck_channels,
            linear_residual_scale_init=args.residual_scale_init,
            linear_dilation=args.linear_dilation,
            linear_kernel_size=args.linear_kernel_size,
            residual_base_checkpoint=args.residual_base_checkpoint,
            linear_residual_fusion=args.linear_residual_fusion,
            residual_correction_scale_init=args.residual_correction_scale_init,
            residual_correction_l2=args.residual_correction_l2,
            residual_base_gate=args.residual_base_gate,
            residual_base_gate_center=args.residual_base_gate_center,
            residual_base_gate_sharpness=args.residual_base_gate_sharpness,
            residual_base_gate_floor=args.residual_base_gate_floor,
            linear_output_calibration=args.linear_output_calibration,
            calibration_gate_center=args.calibration_gate_center,
            calibration_gate_sharpness=args.calibration_gate_sharpness,
            train_only_calibration=args.train_only_calibration,
            linear_loss_type=args.linear_loss_type,
            linear_target_space=args.linear_target_space,
            smooth_l1_beta=args.smooth_l1_beta,
            hybrid_loss_alpha=args.hybrid_loss_alpha,
            track_loss_weight_preset=args.track_loss_weight_preset,
            track_loss_weights=track_loss_weights,
            signal_loss_weight_preset=args.signal_loss_weight_preset,
            signal_high_threshold=args.signal_high_threshold,
            signal_top_quantile=args.signal_top_quantile,
            region_loss_weight_preset=args.region_loss_weight_preset,
            gtf=args.gtf,
            promoter_radius_bp=args.promoter_radius_bp,
            gene_aux_loss_weight=args.gene_aux_loss_weight,
            gene_aux_min_bases=args.gene_aux_min_bases,
            loss_space=loss_space_for_config(
                args.head_type,
                args.loss_type,
                args.linear_target_space,
                args.linear_loss_type,
            ),
            loss_type=args.loss_type,
            multinomial_num_segments=args.multinomial_num_segments,
            positional_weight=args.positional_weight,
            count_weight=args.count_weight,
            mse_weight=args.mse_weight,
            poisson_weight=args.poisson_weight,
            track_means_source=args.track_means_source,
            track_means_tsv=args.track_means_tsv,
            track_means_max_examples=args.track_means_max_examples,
        )
        print(f"saved_checkpoint\t{checkpoint_path}")
    cleanup_distributed(is_distributed)


if __name__ == "__main__":
    main()
