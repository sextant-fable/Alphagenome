#!/usr/bin/env python3
"""Fine-tune an 11-track C. elegans RNA-seq adapter on frozen AlphaGenome."""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import random
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import torch
import torch.distributed as dist
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
    compute_nonzero_track_means_from_dataset,
    count_parameters,
    evaluate_masked_mse,
    format_prediction_shape,
    load_track_means_from_grouped_qc,
    load_adapter_checkpoint,
    masked_adapter_loss,
    parse_resolutions,
    pick_device,
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
        choices=["constant", "cosine"],
        default="constant",
        help="Learning-rate schedule after optional warmup.",
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


def learning_rate_for_step(
    *,
    step: int,
    base_lr: float,
    max_steps: int,
    warmup_steps: int,
    schedule: str,
) -> float:
    if warmup_steps > 0 and step <= warmup_steps:
        return base_lr * (step / warmup_steps)

    if schedule == "cosine":
        decay_steps = max(max_steps - warmup_steps, 1)
        decay_step = min(max(step - warmup_steps, 0), decay_steps)
        cosine = 0.5 * (1.0 + math.cos(math.pi * decay_step / decay_steps))
        return base_lr * cosine

    return base_lr


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
    lora_enabled: bool,
    lora_target_modules: list[str],
    lora_rank: int,
    lora_alpha: int,
    freeze_head: bool,
    init_adapter_checkpoint: str | None,
    head_type: str,
    head_resolutions: tuple[int, ...],
    loss_space: str,
    track_means_source: str,
    track_means_tsv: str,
    track_means_max_examples: int | None,
) -> None:
    checkpoint = {
        "adapter_head_state_dict": model.head.state_dict(),
        "n_tracks": n_tracks,
        "embedding_resolution": embedding_resolution,
        "head_type": head_type,
        "head_resolutions": list(head_resolutions),
        "organism_index": organism_index,
        "target_transform": target_transform,
        "loss_space": loss_space,
        "base_weights": str(weights),
        "best_valid_loss": best_valid_loss,
        "best_valid_step": best_valid_step,
        "lora_enabled": lora_enabled,
        "lora_target_modules": lora_target_modules,
        "lora_rank": lora_rank,
        "lora_alpha": lora_alpha,
        "lora_state_dict": (
            lora_state_dict(model.base_model) if lora_enabled else {}
        ),
        "freeze_head": freeze_head,
        "init_adapter_checkpoint": init_adapter_checkpoint,
        "track_means_source": track_means_source,
        "track_means_tsv": track_means_tsv,
        "track_means_max_examples": track_means_max_examples,
    }
    torch.save(checkpoint, path)


def main() -> None:
    args = parse_args()
    set_seed(args.seed)
    is_distributed, rank, local_rank, world_size, device = setup_distributed(args.device)
    if args.grad_accum_steps < 1:
        raise ValueError("--grad-accum-steps must be >= 1")
    if args.warmup_steps < 0:
        raise ValueError("--warmup-steps must be >= 0")
    if args.freeze_head and not args.enable_last_block_lora:
        raise ValueError("--freeze-head requires --enable-last-block-lora")
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
            target_transform=args.target_transform,
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
    print_main(rank, f"head_type\t{args.head_type}")
    print_main(rank, f"head_resolutions\t{','.join(map(str, head_resolutions))}")
    print_main(
        rank,
        "loss_space\t"
        f"{'genome_tracks_model_scaled' if args.head_type == 'genome-tracks' else 'target_transform'}",
    )
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
    print_main(rank, f"seed\t{args.seed}")
    print_main(rank, f"device\t{device}")
    print_main(rank, f"distributed\t{is_distributed}")
    print_main(rank, f"rank\t{rank}")
    print_main(rank, f"local_rank\t{local_rank}")
    print_main(rank, f"world_size\t{world_size}")
    if device.type == "cuda":
        print_main(rank, f"cuda_device\t{torch.cuda.get_device_name(device)}")

    base_model = AlphaGenome.from_pretrained(weights, device=device)
    model = RnaSeq11Adapter(
        base_model,
        n_tracks=n_tracks,
        embedding_resolution=args.embedding_resolution,
        organism_index=args.organism_index,
        encode_requires_grad=args.enable_last_block_lora,
        head_type=args.head_type,
        head_resolutions=head_resolutions,
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
        if args.head_type == "linear" and checkpoint_resolution != args.embedding_resolution:
            raise ValueError(
                "Init checkpoint embedding_resolution="
                f"{checkpoint_resolution}, requested={args.embedding_resolution}"
            )
        model.head.load_state_dict(init_checkpoint["adapter_head_state_dict"])

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
    best_valid_step: int | None = None
    train_batches = infinite_batches(train_loader, train_sampler)
    for step in range(1, args.max_steps + 1):
        lr = learning_rate_for_step(
            step=step,
            base_lr=args.learning_rate,
            max_steps=args.max_steps,
            warmup_steps=args.warmup_steps,
            schedule=args.lr_schedule,
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

            prediction = train_model(
                dna_sequence,
                target_length=rna_seq.shape[-1],
                return_scaled=model.uses_scaled_targets,
            )
            loss = masked_adapter_loss(model, prediction, rna_seq, rna_seq_mask)
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
            if is_main_process(rank):
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
                    lr=lr,
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
                        lora_enabled=args.enable_last_block_lora,
                        lora_target_modules=lora_target_modules,
                        lora_rank=args.lora_rank,
                        lora_alpha=args.lora_alpha,
                        freeze_head=args.freeze_head,
                        init_adapter_checkpoint=args.init_adapter_checkpoint,
                        head_type=args.head_type,
                        head_resolutions=head_resolutions,
                        loss_space=(
                            "genome_tracks_model_scaled"
                            if args.head_type == "genome-tracks"
                            else "target_transform"
                        ),
                        track_means_source=args.track_means_source,
                        track_means_tsv=args.track_means_tsv,
                        track_means_max_examples=args.track_means_max_examples,
                    )
                    print(f"saved_best_checkpoint\t{best_checkpoint_path}")
                    print(f"best_valid_step\t{best_valid_step}")
                    print(f"best_valid_loss\t{best_valid_loss:.6g}")
            if is_distributed:
                dist.barrier()

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
            lora_enabled=args.enable_last_block_lora,
            lora_target_modules=lora_target_modules,
            lora_rank=args.lora_rank,
            lora_alpha=args.lora_alpha,
            freeze_head=args.freeze_head,
            init_adapter_checkpoint=args.init_adapter_checkpoint,
            head_type=args.head_type,
            head_resolutions=head_resolutions,
            loss_space=(
                "genome_tracks_model_scaled"
                if args.head_type == "genome-tracks"
                else "target_transform"
            ),
            track_means_source=args.track_means_source,
            track_means_tsv=args.track_means_tsv,
            track_means_max_examples=args.track_means_max_examples,
        )
        print(f"saved_checkpoint\t{checkpoint_path}")
    cleanup_distributed(is_distributed)


if __name__ == "__main__":
    main()
