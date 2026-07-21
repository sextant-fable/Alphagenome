#!/usr/bin/env python3
"""Audit v2 scaling, loss, augmentation, and A/B/C model construction on CPU."""

from __future__ import annotations

import csv
import gc
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import torch

from alphagenome_pytorch import AlphaGenome
from alphagenome_pytorch.extensions.finetuning.adapters import (
    LoRA,
    apply_lora,
    get_adapter_params,
)
from scripts import v2_training_components as components
from scripts import v2_subprocess


REPO_ROOT = Path(__file__).resolve().parents[1]
METADATA_DIR = REPO_ROOT / "alphagenome_custom/metadata/v2"
WEIGHTS_PATH = REPO_ROOT / "weights/alphagenome_pytorch/model_all_folds.safetensors"
MEANS_PATH = METADATA_DIR / "track_nonzero_means_v2.tsv"
AUDIT_PATH = METADATA_DIR / "p5_component_audit.json"
MODEL_SPEC_PATH = METADATA_DIR / "model_specs_v2.json"
SPLIT_REGISTRY_PATH = METADATA_DIR / "split_registry_v2.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(8 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def resolve_module(root: torch.nn.Module, path: str) -> torch.nn.Module:
    module = root
    for part in path.split("."):
        module = getattr(module, part)
    return module


def main() -> None:
    subprocess.run(
        v2_subprocess.module_command("compute_v2_track_means"),
        cwd=REPO_ROOT,
        check=True,
    )
    test_run = subprocess.run(
        [sys.executable, "-m", "unittest", "discover", "-s", "tests", "-v"],
        cwd=REPO_ROOT,
        check=False,
        capture_output=True,
        text=True,
    )
    if test_run.returncode != 0:
        print(test_run.stdout)
        print(test_run.stderr, file=sys.stderr)
        raise RuntimeError("P5 unit tests failed")

    means_rows = read_tsv(MEANS_PATH)
    if len(means_rows) != 241:
        raise RuntimeError(f"Expected 241 track means, got {len(means_rows)}")
    fold_means = torch.tensor(
        [float(row["fold_1_train_nonzero_mean"]) for row in means_rows],
        dtype=torch.float32,
    )
    if not WEIGHTS_PATH.is_file():
        raise RuntimeError(f"Missing base checkpoint: {WEIGHTS_PATH}")

    base_model = AlphaGenome.from_pretrained(WEIGHTS_PATH, device="cpu")
    for parameter in base_model.parameters():
        parameter.requires_grad_(False)
    expanded_paths = components.add_c_elegans_organism_embeddings(base_model)
    lora_targets = ["tower.blocks.8.mha", "tower.blocks.8.mlp"]
    apply_lora(base_model, target_modules=lora_targets, rank=8, alpha=16)
    lora_modules = [
        name for name, module in base_model.named_modules() if isinstance(module, LoRA)
    ]
    lora_parameters = list(get_adapter_params(base_model))
    embedding_parameters = [
        resolve_module(base_model, path).weight for path in expanded_paths
    ]
    allowed_parameter_ids = {
        id(parameter) for parameter in [*lora_parameters, *embedding_parameters]
    }
    unexpected_trainable = [
        name
        for name, parameter in base_model.named_parameters()
        if parameter.requires_grad and id(parameter) not in allowed_parameter_ids
    ]
    if unexpected_trainable or not lora_modules:
        raise RuntimeError(
            f"Invalid B trainable scope: unexpected={unexpected_trainable} lora={lora_modules}"
        )
    worm_embedding_finite = all(
        torch.isfinite(resolve_module(base_model, path)(torch.tensor([2]))).all()
        for path in expanded_paths
    )

    torch.manual_seed(20260714)
    head = components.DualResolutionRnaHead(241, fold_means)
    embeddings = {
        1: torch.randn(1, 1536, 1024),
        128: torch.randn(1, 3072, 8),
    }
    predictions = head(embeddings)
    target_1 = torch.rand(1, 241, 1024) * 4
    target_128 = target_1.reshape(1, 241, 8, 128).sum(dim=-1)
    track_mask = torch.ones(1, 241, 1, dtype=torch.bool)
    track_strand = torch.zeros(1, 241, dtype=torch.int8)
    gene_mask = torch.ones(1, 2, 1024, dtype=torch.bool)
    loss, loss_metrics = components.dual_resolution_paper_loss(
        predictions,
        {1: target_1, 128: target_128},
        track_means=fold_means,
        track_mask=track_mask,
        track_strand=track_strand,
        gene_mask=gene_mask,
    )
    loss.backward()
    head_gradients_finite = all(
        parameter.grad is not None and torch.isfinite(parameter.grad).all()
        for parameter in head.parameters()
    )

    baseline = components.WormSequenceBaseline(n_tracks=16, hidden_channels=32)
    baseline_outputs = baseline(torch.rand(1, 4, 1024))
    baseline_loss = baseline_outputs[1].mean() + baseline_outputs[128].mean()
    baseline_loss.backward()
    baseline_gradients_finite = all(
        parameter.grad is not None and torch.isfinite(parameter.grad).all()
        for parameter in baseline.parameters()
    )

    model_specs = {
        "schema_version": 1,
        "models": [spec.__dict__ for spec in components.MODEL_SPECS],
        "resolutions": list(components.RESOLUTIONS),
        "n_tracks": 241,
        "losses": ["log1p_mse", "scaled_poisson_multinomial_gene"],
        "num_segments": components.NUM_SEGMENTS,
        "positional_weight": components.POSITIONAL_WEIGHT,
        "gene_cross_track_weight": components.GENE_CROSS_TRACK_WEIGHT,
        "augmentation": {
            "max_shift_bp": 1024,
            "reverse_complement_probability": 0.5,
        },
        "model_b_lora_targets": lora_targets,
        "model_b_lora_rank": 8,
        "model_b_lora_alpha": 16,
        "model_b_organism_index": 2,
        "promotion_metric": "mean_five_fold_paper_loss_then_log1p_mse_then_mean_per_track_pearson",
        "split_revision": "six_chromosome_blocks_v1",
        "chromosome_x_access": "registered_train_valid_blocks_allowed",
        "locked_test_block_access": "prohibited_until_locked_G5_P6C",
        "final_test_scope": "r6c_single_six_chromosome_block_test",
    }
    MODEL_SPEC_PATH.write_text(
        json.dumps(model_specs, indent=2, sort_keys=True) + "\n"
    )
    audit = {
        "schema_version": 1,
        "phase": "P5",
        "device": "cpu",
        "unit_tests_return_code": test_run.returncode,
        "unit_test_count": test_run.stderr.count(" ... ok"),
        "base_checkpoint_path": str(WEIGHTS_PATH.relative_to(REPO_ROOT)),
        "base_checkpoint_sha256": sha256(WEIGHTS_PATH),
        "base_checkpoint_loaded": True,
        "worm_embedding_paths": expanded_paths,
        "worm_embedding_rows": 3,
        "worm_embedding_index": 2,
        "worm_embedding_finite": worm_embedding_finite,
        "lora_targets": lora_targets,
        "lora_modules": lora_modules,
        "lora_trainable_parameters": sum(p.numel() for p in lora_parameters),
        "worm_embedding_parameter_elements": sum(
            p.numel() for p in embedding_parameters
        ),
        "unexpected_base_trainable_parameters": unexpected_trainable,
        "head_prediction_shapes": {
            str(resolution): list(value.shape)
            for resolution, value in predictions.items()
        },
        "head_loss": float(loss.detach()),
        "head_loss_finite": bool(torch.isfinite(loss).item()),
        "head_gradients_finite": head_gradients_finite,
        "baseline_prediction_shapes": {
            str(resolution): list(value.shape)
            for resolution, value in baseline_outputs.items()
        },
        "baseline_loss": float(baseline_loss.detach()),
        "baseline_gradients_finite": baseline_gradients_finite,
        "loss_metrics": {
            key: float(value.detach()) for key, value in loss_metrics.items()
        },
        "means_path": str(MEANS_PATH.relative_to(REPO_ROOT)),
        "means_sha256": sha256(MEANS_PATH),
        "model_specs_sha256": sha256(MODEL_SPEC_PATH),
        "split_registry_sha256": sha256(SPLIT_REGISTRY_PATH),
        "locked_test_block_signal_reads": 0,
        "scale_implementation": "alphagenome_pytorch.heads.targets_scaling",
        "loss_implementation": "alphagenome_pytorch.losses.multinomial_loss_on_scaled_prediction_and_target",
    }
    AUDIT_PATH.write_text(json.dumps(audit, indent=2, sort_keys=True) + "\n")
    print(json.dumps(audit, indent=2, sort_keys=True))
    del base_model, head, baseline
    gc.collect()
    if not all(
        (
            worm_embedding_finite,
            head_gradients_finite,
            baseline_gradients_finite,
            bool(torch.isfinite(loss).item()),
        )
    ):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
