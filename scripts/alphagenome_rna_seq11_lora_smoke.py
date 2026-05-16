#!/usr/bin/env python3
"""Smoke-test last-block LoRA adaptation for the RNA-seq 11-track adapter."""

from __future__ import annotations

import argparse
import itertools
from pathlib import Path

import torch

from alphagenome_pytorch import AlphaGenome
from alphagenome_pytorch.extensions.finetuning.adapters import (
    LoRA,
    apply_lora,
    get_adapter_params,
)
from alphagenome_rna_seq11_adapter import (
    RnaSeq11Adapter,
    count_parameters,
    evaluate_masked_mse,
    pick_device,
)
from torch_rna_seq_dataset import make_dataloader, masked_mse_loss


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--dataset-dir",
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
    parser.add_argument("--batch-size", type=int, default=1)
    parser.add_argument("--num-workers", type=int, default=0)
    parser.add_argument("--max-steps", type=int, default=1)
    parser.add_argument("--max-train-examples", type=int, default=1)
    parser.add_argument("--max-valid-examples", type=int, default=1)
    parser.add_argument("--max-valid-batches", type=int, default=1)
    parser.add_argument("--learning-rate", type=float, default=1e-4)
    parser.add_argument("--weight-decay", type=float, default=0.0)
    parser.add_argument("--lora-rank", type=int, default=4)
    parser.add_argument("--lora-alpha", type=int, default=8)
    parser.add_argument("--organism-index", type=int, default=0)
    parser.add_argument(
        "--target-transform",
        choices=["log1p", "none"],
        default="log1p",
        help="Use log1p for stable smoke training on high dynamic range targets.",
    )
    parser.add_argument(
        "--device",
        default="auto",
        help="auto, cpu, cuda, cuda:0, mps, etc.",
    )
    return parser.parse_args()


def grad_norm(parameters: list[torch.nn.Parameter]) -> float:
    norms = [
        parameter.grad.detach().float().norm()
        for parameter in parameters
        if parameter.grad is not None
    ]
    if not norms:
        return 0.0
    return float(torch.linalg.vector_norm(torch.stack(norms)).cpu())


def main() -> None:
    args = parse_args()
    dataset_dir = Path(args.dataset_dir)
    valid_dataset_dir = Path(args.valid_dataset_dir)
    weights = Path(args.weights)
    device = pick_device(args.device)

    dataloader = make_dataloader(
        dataset_dir,
        batch_size=args.batch_size,
        shuffle=True,
        num_workers=args.num_workers,
        target_transform=args.target_transform,
        max_examples=args.max_train_examples,
    )
    valid_dataloader = make_dataloader(
        valid_dataset_dir,
        batch_size=args.batch_size,
        shuffle=False,
        num_workers=args.num_workers,
        target_transform=args.target_transform,
        max_examples=args.max_valid_examples,
    )
    n_tracks = int(dataloader.dataset.metadata["n_tracks"])

    print(f"dataset_dir\t{dataset_dir}")
    print(f"valid_dataset_dir\t{valid_dataset_dir}")
    print(f"weights\t{weights}")
    print(f"n_examples\t{len(dataloader.dataset)}")
    print(f"n_valid_examples\t{len(valid_dataloader.dataset)}")
    print(f"n_tracks\t{n_tracks}")
    print(f"target_transform\t{args.target_transform}")
    print("embedding_resolution\t128")
    print(f"organism_index\t{args.organism_index}")
    print(f"learning_rate\t{args.learning_rate}")
    print(f"weight_decay\t{args.weight_decay}")
    print(f"lora_rank\t{args.lora_rank}")
    print(f"lora_alpha\t{args.lora_alpha}")
    print(f"device\t{device}")
    if device.type == "cuda":
        print(f"cuda_device\t{torch.cuda.get_device_name(device)}")

    base_model = AlphaGenome.from_pretrained(weights, device=device)
    model = RnaSeq11Adapter(
        base_model,
        n_tracks=n_tracks,
        embedding_resolution=128,
        organism_index=args.organism_index,
        encode_requires_grad=True,
    ).to(device)

    target_modules = ["tower.blocks.8.mha", "tower.blocks.8.mlp"]
    apply_lora(
        model.base_model,
        target_modules=target_modules,
        rank=args.lora_rank,
        alpha=args.lora_alpha,
    )
    model.to(device)
    lora_module_names = [
        name for name, module in model.base_model.named_modules()
        if isinstance(module, LoRA)
    ]
    if not lora_module_names:
        raise RuntimeError(f"No LoRA modules were applied for targets: {target_modules}")

    lora_parameters = get_adapter_params(model.base_model)
    lora_parameter_ids = {id(parameter) for parameter in lora_parameters}
    base_non_lora_trainable = sum(
        parameter.numel()
        for parameter in model.base_model.parameters()
        if parameter.requires_grad and id(parameter) not in lora_parameter_ids
    )
    trainable_parameters = list(model.head.parameters()) + lora_parameters
    optimizer = torch.optim.AdamW(
        trainable_parameters,
        lr=args.learning_rate,
        weight_decay=args.weight_decay,
    )

    print(f"lora_target_modules\t{','.join(target_modules)}")
    print(f"lora_modules_applied\t{len(lora_module_names)}")
    for name in lora_module_names:
        print(f"lora_module\t{name}")
    print(f"base_parameters\t{count_parameters(model.base_model)}")
    print(f"head_parameters\t{count_parameters(model.head)}")
    print(f"lora_trainable_parameters\t{sum(p.numel() for p in lora_parameters)}")
    print(f"base_non_lora_trainable_parameters\t{base_non_lora_trainable}")
    print(f"total_trainable_parameters\t{sum(p.numel() for p in trainable_parameters)}")

    model.train()
    model.base_model.eval()
    for step, batch in enumerate(itertools.cycle(dataloader), start=1):
        dna_sequence = batch["dna_sequence"].to(device=device, dtype=torch.float32)
        rna_seq = batch["rna_seq"].to(device=device, dtype=torch.float32)
        rna_seq_mask = batch["rna_seq_mask"].to(device=device)

        prediction = model(dna_sequence, target_length=rna_seq.shape[-1])
        loss = masked_mse_loss(prediction, rna_seq, rna_seq_mask)

        optimizer.zero_grad(set_to_none=True)
        loss.backward()
        optimizer.step()

        base_non_lora_has_grad = any(
            parameter.grad is not None
            for parameter in model.base_model.parameters()
            if id(parameter) not in lora_parameter_ids
        )
        lora_has_grad = any(parameter.grad is not None for parameter in lora_parameters)
        head_has_grad = any(parameter.grad is not None for parameter in model.head.parameters())

        print(f"step\t{step}")
        print(f"dna_sequence_shape\t{'x'.join(map(str, dna_sequence.shape))}")
        print(f"rna_seq_shape\t{'x'.join(map(str, rna_seq.shape))}")
        print(f"prediction_shape\t{'x'.join(map(str, prediction.shape))}")
        print(f"loss\t{float(loss.detach().cpu()):.6g}")
        print(f"base_non_lora_has_grad\t{base_non_lora_has_grad}")
        print(f"lora_has_grad\t{lora_has_grad}")
        print(f"head_has_grad\t{head_has_grad}")
        print(f"lora_grad_norm\t{grad_norm(lora_parameters):.6g}")
        print(f"head_grad_norm\t{grad_norm(list(model.head.parameters())):.6g}")
        if device.type == "cuda":
            print(
                "cuda_max_memory_allocated_mb\t"
                f"{torch.cuda.max_memory_allocated() / 1024**2:.1f}"
            )

        if step >= args.max_steps:
            break

    valid_loss, valid_batches, valid_examples = evaluate_masked_mse(
        model,
        valid_dataloader,
        device=device,
        max_batches=args.max_valid_batches,
    )
    print(f"valid_batches\t{valid_batches}")
    print(f"valid_examples\t{valid_examples}")
    print(f"valid_loss\t{valid_loss:.6g}")


if __name__ == "__main__":
    main()
