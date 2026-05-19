#!/usr/bin/env python3
"""Smoke-train an 11-track C. elegans RNA-seq head on frozen AlphaGenome."""

from __future__ import annotations

import argparse
import itertools
from pathlib import Path

import torch

from alphagenome_pytorch import AlphaGenome
from alphagenome_rna_seq11_adapter import (
    RnaSeq11Adapter,
    count_parameters,
    evaluate_masked_mse,
    format_prediction_shape,
    load_track_means_from_grouped_qc,
    masked_adapter_loss,
    parse_resolutions,
    pick_device,
)
from torch_rna_seq_dataset import make_dataloader


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--dataset-dir",
        default="alphagenome_custom/datasets/rna_seq_npz_pilot_train",
        help="NPZ dataset directory containing manifest.tsv and examples/.",
    )
    parser.add_argument(
        "--valid-dataset-dir",
        default=None,
        help="Optional validation NPZ dataset directory.",
    )
    parser.add_argument(
        "--weights",
        default="weights/alphagenome_pytorch/model_all_folds.safetensors",
        help="Local alphagenome-pytorch safetensors checkpoint.",
    )
    parser.add_argument("--batch-size", type=int, default=1)
    parser.add_argument("--num-workers", type=int, default=0)
    parser.add_argument("--max-steps", type=int, default=1)
    parser.add_argument("--max-train-examples", type=int, default=None)
    parser.add_argument("--max-valid-examples", type=int, default=None)
    parser.add_argument("--max-valid-batches", type=int, default=2)
    parser.add_argument("--learning-rate", type=float, default=1e-3)
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
        help="Head implementation to smoke-train.",
    )
    parser.add_argument(
        "--head-resolutions",
        default=None,
        help="Comma-separated head resolutions, e.g. 128 or 1,128.",
    )
    parser.add_argument(
        "--track-means-source",
        choices=["ones", "grouped-qc"],
        default="grouped-qc",
        help="Track means for genome-tracks scaling.",
    )
    parser.add_argument(
        "--track-means-tsv",
        default="alphagenome_custom/metadata/grouped_bigwig_qc.tsv",
        help="TSV used when --track-means-source grouped-qc.",
    )
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


def grad_norm(module: torch.nn.Module) -> float:
    norms = [
        parameter.grad.detach().float().norm()
        for parameter in module.parameters()
        if parameter.grad is not None
    ]
    if not norms:
        return 0.0
    return float(torch.linalg.vector_norm(torch.stack(norms)).cpu())


def main() -> None:
    args = parse_args()
    head_resolutions = parse_resolutions(
        args.head_resolutions,
        fallback_resolution=args.embedding_resolution,
    )
    if args.head_type == "linear" and head_resolutions != (args.embedding_resolution,):
        raise ValueError("linear head requires --head-resolutions to match --embedding-resolution")
    if args.head_type == "genome-tracks" and args.target_transform != "none":
        raise ValueError(
            "--head-type genome-tracks requires --target-transform none for raw target scaling"
        )
    dataset_dir = Path(args.dataset_dir)
    valid_dataset_dir = (
        Path(args.valid_dataset_dir) if args.valid_dataset_dir is not None else None
    )
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
    valid_dataloader = None
    if valid_dataset_dir is not None:
        valid_dataloader = make_dataloader(
            valid_dataset_dir,
            batch_size=args.batch_size,
            shuffle=False,
            num_workers=args.num_workers,
            target_transform=args.target_transform,
            max_examples=args.max_valid_examples,
        )
    n_tracks = int(dataloader.dataset.metadata["n_tracks"])
    track_means = None
    if args.head_type == "genome-tracks":
        if args.track_means_source == "ones":
            track_means = torch.ones(n_tracks, dtype=torch.float32)
        else:
            track_means = load_track_means_from_grouped_qc(
                args.track_means_tsv,
                n_tracks=n_tracks,
            )

    print(f"dataset_dir\t{dataset_dir}")
    if valid_dataset_dir is not None:
        print(f"valid_dataset_dir\t{valid_dataset_dir}")
    print(f"weights\t{weights}")
    print(f"n_examples\t{len(dataloader.dataset)}")
    if valid_dataloader is not None:
        print(f"n_valid_examples\t{len(valid_dataloader.dataset)}")
    print(f"n_tracks\t{n_tracks}")
    print(f"target_transform\t{args.target_transform}")
    print(f"head_type\t{args.head_type}")
    print(f"head_resolutions\t{','.join(map(str, head_resolutions))}")
    print(
        "loss_space\t"
        f"{'genome_tracks_model_scaled' if args.head_type == 'genome-tracks' else 'target_transform'}"
    )
    print(f"embedding_resolution\t{args.embedding_resolution}")
    print(f"organism_index\t{args.organism_index}")
    print(f"device\t{device}")
    if device.type == "cuda":
        print(f"cuda_device\t{torch.cuda.get_device_name(device)}")

    base_model = AlphaGenome.from_pretrained(weights, device=device)
    model = RnaSeq11Adapter(
        base_model,
        n_tracks=n_tracks,
        embedding_resolution=args.embedding_resolution,
        organism_index=args.organism_index,
        head_type=args.head_type,
        head_resolutions=head_resolutions,
        track_means=track_means,
    ).to(device)
    optimizer = torch.optim.AdamW(model.head.parameters(), lr=args.learning_rate)

    print(f"base_parameters\t{count_parameters(model.base_model)}")
    print(f"total_parameters\t{count_parameters(model)}")
    print(f"trainable_parameters\t{count_parameters(model, trainable_only=True)}")
    print(f"head_parameters\t{count_parameters(model.head)}")

    model.train()
    model.base_model.eval()
    for step, batch in enumerate(itertools.cycle(dataloader), start=1):
        dna_sequence = batch["dna_sequence"].to(device=device, dtype=torch.float32)
        rna_seq = batch["rna_seq"].to(device=device, dtype=torch.float32)
        rna_seq_mask = batch["rna_seq_mask"].to(device=device)

        prediction = model(
            dna_sequence,
            target_length=rna_seq.shape[-1],
            return_scaled=model.uses_scaled_targets,
        )
        loss = masked_adapter_loss(model, prediction, rna_seq, rna_seq_mask)

        optimizer.zero_grad(set_to_none=True)
        loss.backward()
        optimizer.step()

        base_has_grad = any(
            p.grad is not None for p in model.base_model.parameters()
        )
        print(f"step\t{step}")
        print(f"dna_sequence_shape\t{'x'.join(map(str, dna_sequence.shape))}")
        print(f"rna_seq_shape\t{'x'.join(map(str, rna_seq.shape))}")
        print(f"prediction_shape\t{format_prediction_shape(prediction)}")
        print(f"loss\t{float(loss.detach().cpu()):.6g}")
        print(f"base_has_grad\t{base_has_grad}")
        print(f"head_weight_grad_norm\t{grad_norm(model.head):.6g}")
        if device.type == "cuda":
            print(f"cuda_max_memory_allocated_mb\t{torch.cuda.max_memory_allocated() / 1024**2:.1f}")

        if step >= args.max_steps:
            break

    if valid_dataloader is not None:
        valid_loss, valid_batches, _valid_examples = evaluate_masked_mse(
            model,
            valid_dataloader,
            device=device,
            max_batches=args.max_valid_batches,
        )
        print(f"valid_batches\t{valid_batches}")
        print(f"valid_loss\t{valid_loss:.6g}")


if __name__ == "__main__":
    main()
