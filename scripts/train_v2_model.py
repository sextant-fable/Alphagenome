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
from scripts.v2_bigwig_dataset import V2BigWigDataset
from scripts import v2_training_components as components


REPO_ROOT = Path(__file__).resolve().parents[1]
METADATA_DIR = REPO_ROOT / "alphagenome_custom/metadata/v2"
WEIGHTS_PATH = REPO_ROOT / "weights/alphagenome_pytorch/model_all_folds.safetensors"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--model", choices=["A", "B", "C"], required=True)
    parser.add_argument("--fold", type=int, choices=range(0, 6), default=1)
    parser.add_argument("--seed", type=int, default=20260714)
    parser.add_argument("--loss", choices=["paper", "log1p_mse"], default="paper")
    parser.add_argument("--max-steps", type=int, default=2)
    parser.add_argument("--sequence-length", type=int, default=131072)
    parser.add_argument("--learning-rate", type=float, default=1e-4)
    parser.add_argument("--weight-decay", type=float, default=1e-4)
    parser.add_argument("--gene-weight", type=float, default=0.1)
    parser.add_argument("--max-shift-bp", type=int, default=1024)
    parser.add_argument("--reverse-complement-probability", type=float, default=0.5)
    parser.add_argument("--hidden-channels", type=int, default=64)
    parser.add_argument("--mean-column")
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--device", default="cuda")
    return parser.parse_args()


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
) -> torch.nn.Module:
    if model_id == "C":
        return components.WormSequenceBaseline(
            n_tracks=fold_means.numel(), hidden_channels=hidden_channels
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
    components.add_c_elegans_organism_embeddings(base_model)
    base_model.encoder.gradient_checkpointing = True
    base_model.tower.gradient_checkpointing = True
    base_model.decoder.gradient_checkpointing = True
    apply_lora(
        base_model,
        target_modules=["tower.blocks.8.mha", "tower.blocks.8.mlp"],
        rank=8,
        alpha=16,
    )
    return components.AlphaGenomeRnaModel(
        base_model,
        n_tracks=fold_means.numel(),
        track_means=fold_means,
        base_organism_index=2,
        encode_requires_grad=True,
    ).to(device)


def main() -> None:
    args = parse_args()
    if args.sequence_length % 128 or args.sequence_length > 2**20:
        raise ValueError("sequence-length must be <= 1048576 and divisible by 128")
    output_dir = REPO_ROOT / args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    device = torch.device(args.device)
    if device.type != "cuda" or not torch.cuda.is_available():
        raise RuntimeError("Registered P6 training requires CUDA")
    visible = os.environ.get("CUDA_VISIBLE_DEVICES", "")
    if not visible or any(index not in {"2", "3"} for index in visible.split(",")):
        raise RuntimeError("CUDA_VISIBLE_DEVICES must restrict P6 to physical GPU 2 and/or 3")
    set_seed(args.seed)

    means_path = METADATA_DIR / "track_nonzero_means_v2.tsv"
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
        REPO_ROOT / "alphagenome_custom/intervals/v2/development_train.tsv"
        if args.fold == 0
        else REPO_ROOT / f"alphagenome_custom/intervals/v2/fold_{args.fold}/train.tsv"
    )
    base_dataset = V2BigWigDataset(intervals, max_io_workers=16)
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
    model = build_model(
        args.model, fold_means.detach().cpu(), device, args.hidden_channels
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
                    )
                else:
                    loss, loss_metrics = components.log1p_mse_loss(
                        predictions,
                        targets,
                        track_means=fold_means,
                        track_mask=track_mask,
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
            "max_shift_bp": args.max_shift_bp,
            "reverse_complement_probability": args.reverse_complement_probability,
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
        "max_shift_bp": args.max_shift_bp,
        "reverse_complement_probability": args.reverse_complement_probability,
        "physical_cuda_visible_devices": visible,
        "cuda_device_name": torch.cuda.get_device_name(0),
        "trainable_parameters": sum(parameter.numel() for parameter in trainable),
        "hidden_channels": args.hidden_channels,
        "metrics": metrics_rows,
        "checkpoint_path": str(checkpoint_path.relative_to(REPO_ROOT)),
        "checkpoint_sha256": sha256(checkpoint_path),
        "checkpoint_reload_verified": True,
        "means_sha256": sha256(means_path),
        "intervals_sha256": sha256(intervals),
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
