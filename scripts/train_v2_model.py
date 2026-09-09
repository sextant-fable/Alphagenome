#!/usr/bin/env python3
"""Train one registered v2 A/B/C configuration with auditable inputs."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import random
import time

import numpy as np
import torch
from torch.utils.data import DataLoader

from alphagenome_pytorch import AlphaGenome
from alphagenome_pytorch.extensions.finetuning.adapters import apply_lora
from scripts import v2_bigwig_dataset
from scripts.v2_bigwig_dataset import V2BigWigDataset
from scripts import v2_training_components as components


REPO_ROOT = Path(__file__).resolve().parents[1]
METADATA_DIR = REPO_ROOT / "alphagenome_custom/metadata/v2"
WEIGHTS_PATH = REPO_ROOT / "weights/alphagenome_pytorch/model_all_folds.safetensors"
B_MODEL_IDS = {"B", "B_no_lora", "B_lora_only", "B_1bp_head_only", "B_128bp_head_only"}
LORA_TARGET_MODULES = ["tower.blocks.8.mha", "tower.blocks.8.mlp"]
LORA_TARGET_MODULE_SETS = {
    "final_block": LORA_TARGET_MODULES,
    "previous_block": ["tower.blocks.7.mha", "tower.blocks.7.mlp"],
}


def adaptation_metadata(
    model_id: str,
    *,
    lora_rank: int = 8,
    lora_alpha: int = 16,
    lora_target_modules: list[str] | None = None,
) -> dict[str, object]:
    target_modules = list(
        LORA_TARGET_MODULES if lora_target_modules is None else lora_target_modules
    )
    if model_id == "B":
        return {
            "base_organism_index": 2,
            "c_elegans_organism_embedding": True,
            "lora_enabled": True,
            "lora_rank": lora_rank,
            "lora_alpha": lora_alpha,
            "lora_target_modules": target_modules,
            "trunk_policy": "worm_embeddings_and_lora_only",
        }
    if model_id == "B_no_lora":
        return {
            "base_organism_index": 2,
            "c_elegans_organism_embedding": True,
            "lora_enabled": False,
            "lora_rank": None,
            "lora_alpha": None,
            "lora_target_modules": [],
            "trunk_policy": "worm_embeddings_only",
        }
    if model_id == "B_lora_only":
        return {
            "base_organism_index": 0,
            "c_elegans_organism_embedding": False,
            "lora_enabled": True,
            "lora_rank": lora_rank,
            "lora_alpha": lora_alpha,
            "lora_target_modules": target_modules,
            "trunk_policy": "lora_only",
        }
    if model_id == "B_1bp_head_only":
        return {
            "base_organism_index": 2,
            "c_elegans_organism_embedding": True,
            "lora_enabled": True,
            "lora_rank": lora_rank,
            "lora_alpha": lora_alpha,
            "lora_target_modules": target_modules,
            "trunk_policy": "worm_embeddings_and_lora_only",
            "learned_head_resolutions": [1],
            "derived_output_resolutions": [128],
            "derived_128bp_policy": "sum_pool_unscaled_1bp_then_rescale_128bp",
            "output_head_policy": "learned_1bp_branch_only_with_full_b_adaptation",
        }
    if model_id == "B_128bp_head_only":
        return {
            "base_organism_index": 2,
            "c_elegans_organism_embedding": True,
            "lora_enabled": True,
            "lora_rank": lora_rank,
            "lora_alpha": lora_alpha,
            "lora_target_modules": target_modules,
            "trunk_policy": "worm_embeddings_and_lora_only",
            "learned_head_resolutions": [128],
            "derived_output_resolutions": [1],
            "derived_1bp_policy": "uniformly_distribute_unscaled_128bp_then_rescale_1bp",
            "output_head_policy": "learned_128bp_branch_only_with_full_b_adaptation",
        }
    if model_id == "Basenji2_style":
        return {
            "base_organism_index": None,
            "c_elegans_organism_embedding": False,
            "lora_enabled": False,
            "lora_rank": None,
            "lora_alpha": None,
            "lora_target_modules": [],
            "trunk_policy": "from_scratch_basenji2_style",
            "learned_head_resolutions": [128],
            "derived_output_resolutions": [1],
            "derived_1bp_policy": "uniformly_distribute_unscaled_128bp_then_rescale_1bp",
            "public_baseline_source": "calico/basenji commit 06ce5d387e20b47184d05433b3983163c5f923cd",
        }
    if model_id == "C_size_matched":
        return {
            "base_organism_index": None,
            "c_elegans_organism_embedding": False,
            "lora_enabled": False,
            "lora_rank": None,
            "lora_alpha": None,
            "lora_target_modules": [],
            "trunk_policy": "from_scratch_size_matched",
        }
    return {
        "base_organism_index": 0 if model_id in {"A", "D_base", "D"} else None,
        "c_elegans_organism_embedding": False,
        "lora_enabled": False,
        "lora_rank": None,
        "lora_alpha": None,
        "lora_target_modules": [],
        "trunk_policy": "frozen" if model_id != "C" else "from_scratch",
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--model",
        choices=[
            "A", "B", "B_no_lora", "B_lora_only", "B_1bp_head_only", "B_128bp_head_only",
            "Basenji2_style",
            "C", "C_size_matched", "D_base", "D",
        ],
        required=True,
    )
    parser.add_argument("--fold", type=int, choices=range(0, 6), default=1)
    parser.add_argument("--seed", type=int, default=20260714)
    parser.add_argument("--loss", choices=["paper", "log1p_mse"], default="paper")
    parser.add_argument("--max-steps", type=int, default=2)
    parser.add_argument("--sequence-length", type=int, default=131072)
    parser.add_argument("--learning-rate", type=float, default=1e-4)
    parser.add_argument("--weight-decay", type=float, default=1e-4)
    parser.add_argument("--gene-weight", type=float, default=0.1)
    parser.add_argument("--loss-1bp-weight", type=float, default=1.0)
    parser.add_argument("--loss-128bp-weight", type=float, default=1.0)
    parser.add_argument("--max-shift-bp", type=int, default=1024)
    parser.add_argument("--reverse-complement-probability", type=float, default=0.5)
    parser.add_argument("--hidden-channels", type=int, default=64)
    parser.add_argument("--mean-column")
    parser.add_argument("--means-path")
    parser.add_argument("--track-manifest")
    parser.add_argument("--group-manifest")
    parser.add_argument("--train-intervals-path")
    parser.add_argument(
        "--frozen-base-checkpoint",
        help="Required for model D; a D_base checkpoint whose Conv5 head is frozen.",
    )
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--device", default="cuda")
    parser.add_argument("--lora-rank", type=int, default=8)
    parser.add_argument("--lora-alpha", type=int, default=16)
    parser.add_argument(
        "--lora-target-set",
        choices=sorted(LORA_TARGET_MODULE_SETS),
        default="final_block",
    )
    return parser.parse_args()


def trainable_parameter_count(model: torch.nn.Module) -> int:
    return sum(parameter.numel() for parameter in model.parameters() if parameter.requires_grad)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(8 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def set_seed(seed: int) -> None:
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    torch.use_deterministic_algorithms(True, warn_only=True)


def build_optimizer(
    model: torch.nn.Module,
    *,
    learning_rate: float,
    weight_decay: float,
) -> tuple[torch.optim.Optimizer, list[torch.nn.Parameter]]:
    trainable_named = [
        (name, parameter)
        for name, parameter in model.named_parameters()
        if parameter.requires_grad
    ]
    if not trainable_named:
        raise RuntimeError("Model has no trainable parameters")
    organism_embeddings = [
        parameter
        for name, parameter in trainable_named
        if name.endswith("organism_embed.weight")
    ]
    organism_embedding_ids = {id(parameter) for parameter in organism_embeddings}
    regular_parameters = [
        parameter
        for _, parameter in trainable_named
        if id(parameter) not in organism_embedding_ids
    ]
    parameter_groups = []
    if regular_parameters:
        parameter_groups.append(
            {"params": regular_parameters, "weight_decay": weight_decay}
        )
    if organism_embeddings:
        # The gradient hook freezes rows 0/1; zero decay prevents AdamW from
        # changing those rows despite their zero gradients.
        parameter_groups.append({"params": organism_embeddings, "weight_decay": 0.0})
    optimizer = torch.optim.AdamW(parameter_groups, lr=learning_rate)
    return optimizer, [parameter for _, parameter in trainable_named]


def build_model(
    model_id: str,
    fold_means: torch.Tensor,
    device: torch.device,
    hidden_channels: int,
    frozen_base_checkpoint: Path | None = None,
    lora_rank: int = 8,
    lora_alpha: int = 16,
    lora_target_modules: list[str] | None = None,
) -> torch.nn.Module:
    if model_id in {"C", "C_size_matched"}:
        if model_id == "C_size_matched" and hidden_channels != 152:
            raise ValueError("C_size_matched requires hidden_channels=152")
        baseline_class = (
            components.SizeMatchedWormSequenceBaseline
            if model_id == "C_size_matched"
            else components.WormSequenceBaseline
        )
        return baseline_class(
            n_tracks=fold_means.numel(), hidden_channels=hidden_channels
        ).to(device)
    if model_id == "Basenji2_style":
        return components.Basenji2StyleRnaBaseline(
            n_tracks=fold_means.numel(),
            track_means=fold_means,
            hidden_channels=hidden_channels,
        ).to(device)
    base_model = AlphaGenome.from_pretrained(WEIGHTS_PATH, device=device)
    for parameter in base_model.parameters():
        parameter.requires_grad_(False)
    if model_id == "A":
        return components.AlphaGenomeRnaModel(
            base_model,
            n_tracks=fold_means.numel(),
            track_means=fold_means,
            base_organism_index=0,
            encode_requires_grad=False,
        ).to(device)
    if model_id == "D_base":
        return components.Legacy128bpBaseRnaModel(
            base_model, n_tracks=fold_means.numel()
        ).to(device)
    if model_id == "D":
        if frozen_base_checkpoint is None:
            raise ValueError("Model D requires --frozen-base-checkpoint")
        checkpoint = torch.load(
            frozen_base_checkpoint, map_location="cpu", weights_only=False
        )
        if checkpoint.get("model_id") != "D_base":
            raise ValueError("--frozen-base-checkpoint must be from model D_base")
        model = components.LegacyResidualRnaModel(
            base_model, n_tracks=fold_means.numel()
        ).to(device)
        model.load_frozen_base_state(checkpoint.get("trainable_model_state", {}))
        return model
    if model_id not in B_MODEL_IDS:
        raise ValueError(f"Unsupported model {model_id}")
    if model_id != "B_lora_only":
        components.add_c_elegans_organism_embeddings(base_model)
    base_model.encoder.gradient_checkpointing = True
    base_model.tower.gradient_checkpointing = True
    base_model.decoder.gradient_checkpointing = True
    if model_id != "B_no_lora":
        apply_lora(
            base_model,
            target_modules=(
                LORA_TARGET_MODULES
                if lora_target_modules is None
                else lora_target_modules
            ),
            rank=lora_rank,
            alpha=lora_alpha,
        )
    return components.AlphaGenomeRnaModel(
        base_model,
        n_tracks=fold_means.numel(),
        track_means=fold_means,
        base_organism_index=0 if model_id == "B_lora_only" else 2,
        encode_requires_grad=True,
        head_resolutions=(1,) if model_id == "B_1bp_head_only" else (128,) if model_id == "B_128bp_head_only" else components.RESOLUTIONS,
        derive_128bp_from_1bp=model_id == "B_1bp_head_only",
        derive_1bp_from_128bp=model_id == "B_128bp_head_only",
    ).to(device)


def main() -> None:
    args = parse_args()
    if args.sequence_length % 128 or args.sequence_length > 2**20:
        raise ValueError("sequence-length must be <= 1048576 and divisible by 128")
    if args.loss_1bp_weight < 0 or args.loss_128bp_weight < 0:
        raise ValueError("resolution loss weights must be non-negative")
    if args.loss_1bp_weight == 0 and args.loss_128bp_weight == 0:
        raise ValueError("at least one resolution loss weight must be positive")
    output_dir = REPO_ROOT / args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    frozen_base_checkpoint = (
        REPO_ROOT / args.frozen_base_checkpoint
        if args.frozen_base_checkpoint is not None
        else None
    )
    if args.model == "D" and (
        frozen_base_checkpoint is None or not frozen_base_checkpoint.is_file()
    ):
        raise ValueError("Model D requires an existing --frozen-base-checkpoint")
    if args.model != "D" and frozen_base_checkpoint is not None:
        raise ValueError("--frozen-base-checkpoint is valid only for model D")
    device = torch.device(args.device)
    if device.type != "cuda" or not torch.cuda.is_available():
        raise RuntimeError("Registered P6 training requires CUDA")
    visible = os.environ.get("CUDA_VISIBLE_DEVICES", "")
    if not visible or any(index not in {"2", "3"} for index in visible.split(",")):
        raise RuntimeError("CUDA_VISIBLE_DEVICES must restrict P6 to physical GPU 2 and/or 3")
    set_seed(args.seed)

    means_path = REPO_ROOT / (
        args.means_path
        if args.means_path is not None
        else "alphagenome_custom/metadata/v2/track_nonzero_means_v2.tsv"
    )
    means_rows = read_tsv(means_path)
    default_mean_column = (
        "development_train_nonzero_mean"
        if args.fold == 0
        else f"fold_{args.fold}_train_nonzero_mean"
    )
    mean_column = args.mean_column or default_mean_column
    if args.fold == 0 and mean_column != "development_train_nonzero_mean":
        raise ValueError("fold 0 requires development_train_nonzero_mean")
    if mean_column not in means_rows[0]:
        raise ValueError(f"Unknown mean column: {mean_column}")
    fold_means = torch.tensor(
        [float(row[mean_column]) for row in means_rows],
        dtype=torch.float32,
        device=device,
    )
    intervals = (
        REPO_ROOT / args.train_intervals_path
        if args.train_intervals_path is not None
        else (
            REPO_ROOT / "alphagenome_custom/intervals/v2/development_train.tsv"
            if args.fold == 0
            else REPO_ROOT / f"alphagenome_custom/intervals/v2/fold_{args.fold}/train.tsv"
        )
    )
    track_manifest = (
        REPO_ROOT / args.track_manifest
        if args.track_manifest is not None
        else v2_bigwig_dataset.DEFAULT_TRACKS
    )
    group_manifest = (
        REPO_ROOT / args.group_manifest
        if args.group_manifest is not None
        else v2_bigwig_dataset.DEFAULT_GROUPS
    )
    base_dataset = V2BigWigDataset(
        intervals,
        track_manifest_path=track_manifest,
        group_manifest_path=group_manifest,
        max_io_workers=16,
    )
    dataset = components.AugmentedV2Dataset(
        base_dataset,
        max_shift_bp=args.max_shift_bp,
        reverse_complement_probability=args.reverse_complement_probability,
        seed=args.seed,
        sequence_length=args.sequence_length,
    )
    generator = torch.Generator().manual_seed(args.seed)
    dataloader = DataLoader(
        dataset,
        batch_size=1,
        shuffle=True,
        num_workers=0,
        generator=generator,
    )
    lora_target_modules = LORA_TARGET_MODULE_SETS[args.lora_target_set]
    model = build_model(
        args.model,
        fold_means.detach().cpu(),
        device,
        args.hidden_channels,
        frozen_base_checkpoint,
        lora_rank=args.lora_rank,
        lora_alpha=args.lora_alpha,
        lora_target_modules=lora_target_modules,
    )
    optimizer, trainable = build_optimizer(
        model,
        learning_rate=args.learning_rate,
        weight_decay=args.weight_decay,
    )
    scaler_enabled = False
    metrics_rows = []
    started = time.monotonic()
    model.train()
    epoch = 0
    while len(metrics_rows) < args.max_steps:
        dataset.set_epoch(epoch)
        for batch in dataloader:
            step = len(metrics_rows) + 1
            dna = batch["dna_sequence"].to(device=device, dtype=torch.float32)
            targets = {
                1: batch["target_1bp"].to(device=device, dtype=torch.float32),
                128: batch["target_128bp"].to(device=device, dtype=torch.float32),
            }
            track_mask = batch["track_mask"].to(device)
            track_strand = batch["track_strand"].to(device)
            gene_mask = batch["gene_mask"].to(device)
            optimizer.zero_grad(set_to_none=True)
            with torch.autocast(device_type="cuda", dtype=torch.bfloat16):
                predictions = model(dna)
                if args.loss == "paper":
                    loss, loss_metrics = components.dual_resolution_paper_loss(
                        predictions,
                        targets,
                        track_means=fold_means,
                        track_mask=track_mask,
                        track_strand=track_strand,
                        gene_mask=gene_mask,
                        gene_weight=args.gene_weight,
                        resolution_weights={1: args.loss_1bp_weight, 128: args.loss_128bp_weight},
                    )
                else:
                    loss, loss_metrics = components.log1p_mse_loss(
                        predictions,
                        targets,
                        track_means=fold_means,
                        track_mask=track_mask,
                        resolution_weights={1: args.loss_1bp_weight, 128: args.loss_128bp_weight},
                    )
            if not torch.isfinite(loss):
                raise RuntimeError(f"Non-finite loss at step {step}: {loss}")
            loss.backward()
            gradient_norm = torch.nn.utils.clip_grad_norm_(trainable, max_norm=1.0)
            if not torch.isfinite(gradient_norm):
                raise RuntimeError(f"Non-finite gradient at step {step}")
            optimizer.step()
            metrics_rows.append(
                {
                    "step": step,
                    "epoch": epoch,
                    "loss": float(loss.detach()),
                    "gradient_norm": float(gradient_norm.detach()),
                    "shift_bp": int(batch["shift_bp"].item()),
                    "crop_offset_bp": int(batch["crop_offset_bp"].item()),
                    "reverse_complemented": bool(batch["reverse_complemented"].item()),
                    "cuda_memory_allocated_bytes": torch.cuda.memory_allocated(),
                    "cuda_max_memory_allocated_bytes": torch.cuda.max_memory_allocated(),
                    **{
                        key: float(value.detach())
                        for key, value in loss_metrics.items()
                    },
                }
            )
            print(json.dumps(metrics_rows[-1], sort_keys=True), flush=True)
            if step >= args.max_steps:
                break
        epoch += 1
    checkpoint_path = output_dir / "checkpoint.pt"
    trainable_names = {
        name for name, parameter in model.named_parameters() if parameter.requires_grad
    }
    trainable_state = {
        name: value.detach().cpu()
        for name, value in model.state_dict().items()
        if name in trainable_names
    }
    torch.save(
        {
            "trainable_model_state": trainable_state,
            "trainable_parameter_names": sorted(trainable_names),
            "optimizer": optimizer.state_dict(),
            "model_id": args.model,
            "fold": args.fold,
            "seed": args.seed,
            "loss": args.loss,
            "step": len(metrics_rows),
            "sequence_length": args.sequence_length,
            "hidden_channels": args.hidden_channels,
            "mean_column": mean_column,
            "gene_weight": args.gene_weight,
            "loss_1bp_weight": args.loss_1bp_weight,
            "loss_128bp_weight": args.loss_128bp_weight,
            "max_shift_bp": args.max_shift_bp,
            "reverse_complement_probability": args.reverse_complement_probability,
            "adaptation": adaptation_metadata(
                args.model,
                lora_rank=args.lora_rank,
                lora_alpha=args.lora_alpha,
                lora_target_modules=lora_target_modules,
            ),
            "frozen_base_checkpoint": args.frozen_base_checkpoint,
            "frozen_base_checkpoint_sha256": (
                sha256(frozen_base_checkpoint)
                if frozen_base_checkpoint is not None
                else None
            ),
        },
        checkpoint_path,
    )
    reloaded = torch.load(checkpoint_path, map_location="cpu", weights_only=False)
    if (
        reloaded["step"] != len(metrics_rows)
        or set(reloaded["trainable_model_state"]) != trainable_names
    ):
        raise RuntimeError("Checkpoint reload mismatch")
    run_record = {
        "schema_version": 1,
        "run_type": "smoke" if args.max_steps <= 2 else "training",
        "model": args.model,
        "fold": args.fold,
        "mean_column": mean_column,
        "seed": args.seed,
        "loss": args.loss,
        "steps": len(metrics_rows),
        "epochs_started": epoch,
        "sequence_length": args.sequence_length,
        "n_tracks": len(means_rows),
        "prediction_shapes": {
            str(resolution): list(value.shape)
            for resolution, value in predictions.items()
        },
        "target_shapes": {
            str(resolution): list(value.shape)
            for resolution, value in targets.items()
        },
        "learning_rate": args.learning_rate,
        "gene_weight": args.gene_weight,
        "loss_1bp_weight": args.loss_1bp_weight,
        "loss_128bp_weight": args.loss_128bp_weight,
        "max_shift_bp": args.max_shift_bp,
        "reverse_complement_probability": args.reverse_complement_probability,
        "adaptation": adaptation_metadata(
            args.model,
            lora_rank=args.lora_rank,
            lora_alpha=args.lora_alpha,
            lora_target_modules=lora_target_modules,
        ),
        "frozen_base_checkpoint": args.frozen_base_checkpoint,
        "frozen_base_checkpoint_sha256": (
            sha256(frozen_base_checkpoint)
            if frozen_base_checkpoint is not None
            else None
        ),
        "physical_cuda_visible_devices": visible,
        "cuda_device_name": torch.cuda.get_device_name(0),
        "trainable_parameters": trainable_parameter_count(model),
        "hidden_channels": args.hidden_channels,
        "metrics": metrics_rows,
        "checkpoint_path": str(checkpoint_path.relative_to(REPO_ROOT)),
        "checkpoint_sha256": sha256(checkpoint_path),
        "checkpoint_reload_verified": True,
        "means_sha256": sha256(means_path),
        "intervals_sha256": sha256(intervals),
        "track_manifest": str(track_manifest.relative_to(REPO_ROOT)),
        "track_manifest_sha256": sha256(track_manifest),
        "group_manifest": str(group_manifest.relative_to(REPO_ROOT)),
        "group_manifest_sha256": sha256(group_manifest),
        "elapsed_seconds": time.monotonic() - started,
        "amp_dtype": "bfloat16",
        "gradient_scaler_enabled": scaler_enabled,
        "claim_status": "environment_validation"
        if args.max_steps <= 2
        else "provisional_training",
    }
    (output_dir / "run.json").write_text(
        json.dumps(run_record, indent=2, sort_keys=True) + "\n"
    )
    base_dataset.close()
    print(json.dumps(run_record, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
