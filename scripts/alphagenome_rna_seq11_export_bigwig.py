#!/usr/bin/env python3
"""Export RNA-seq11 adapter predictions as per-track bigWig files."""

from __future__ import annotations

import argparse
import csv
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pyBigWig
import torch

from alphagenome_pytorch import AlphaGenome
from alphagenome_rna_seq11_adapter import (
    RnaSeq11Adapter,
    count_parameters,
    load_adapter_checkpoint,
    parse_resolutions,
    pick_device,
)
from torch_rna_seq_dataset import make_dataloader


DEFAULT_CHECKPOINT = (
    "runs/rna_seq11_1bpB_broad_w12_track_hard15_beta05_from_w10best1000_"
    "lr1e-5_seedrand_1500steps/adapter_head_best.pt"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--checkpoint",
        default=DEFAULT_CHECKPOINT,
        help="Adapter checkpoint to export.",
    )
    parser.add_argument(
        "--weights",
        default="weights/alphagenome_pytorch/model_all_folds.safetensors",
        help="Local alphagenome-pytorch safetensors checkpoint.",
    )
    parser.add_argument(
        "--valid-dataset-dir",
        default="alphagenome_custom/datasets/rna_seq_npz_valid",
        help="Validation NPZ dataset directory.",
    )
    parser.add_argument(
        "--test-dataset-dir",
        default="alphagenome_custom/datasets/rna_seq_npz_test",
        help="Test NPZ dataset directory.",
    )
    parser.add_argument(
        "--splits",
        default="valid,test",
        help="Comma-separated splits to export from {valid,test}.",
    )
    parser.add_argument(
        "--output-dir",
        required=True,
        help="Output directory for split subdirectories containing bigWigs.",
    )
    parser.add_argument(
        "--chrom-sizes",
        default="alphagenome_custom/metadata/chrom_sizes.tsv",
        help="Chromosome sizes TSV.",
    )
    parser.add_argument("--batch-size", type=int, default=1)
    parser.add_argument("--num-workers", type=int, default=0)
    parser.add_argument("--max-examples", type=int, default=None)
    parser.add_argument(
        "--prediction-transform",
        choices=["expm1", "as-is"],
        default="expm1",
        help=(
            "Transform model output before writing. expm1 converts log1p-space "
            "linear predictions back to raw signal scale."
        ),
    )
    parser.add_argument(
        "--clip-min",
        type=float,
        default=0.0,
        help="Clip exported values to this minimum. Use nan to disable.",
    )
    parser.add_argument(
        "--value-chunk-size",
        type=int,
        default=100_000,
        help="Number of bases per pyBigWig addEntries call.",
    )
    parser.add_argument(
        "--device",
        default="auto",
        help="auto, cpu, cuda, cuda:0, mps, etc.",
    )
    return parser.parse_args()


def read_chrom_sizes(path: str | Path) -> list[tuple[str, int]]:
    rows: list[tuple[str, int]] = []
    with Path(path).open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            chrom = row.get("chromosome") or row.get("chrom")
            length = row.get("length") or row.get("size")
            if chrom is None or length is None:
                raise ValueError(f"Cannot parse chromosome sizes from {path}")
            rows.append((chrom, int(length)))
    if not rows:
        raise ValueError(f"No chromosome sizes found in {path}")
    return rows


def read_track_labels(track_metadata: str | Path, n_tracks: int) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    with Path(track_metadata).open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if row.get("output_type") != "RNA_SEQ":
                continue
            rows.append(row)
    rows.sort(key=lambda row: int(row["track_index"]))
    if len(rows) != n_tracks:
        raise ValueError(
            f"Loaded {len(rows)} RNA_SEQ rows from {track_metadata}; expected {n_tracks}"
        )
    return rows


def optional_int(value: object) -> int | None:
    if value is None or value == "":
        return None
    return int(value)


def checkpoint_value(
    checkpoint: dict[str, Any],
    key: str,
    default: Any,
    override: Any = None,
) -> Any:
    if override is not None:
        return override
    return checkpoint.get(key, default)


def build_model_from_checkpoint(
    *,
    checkpoint_path: Path,
    weights: Path,
    n_tracks: int,
    device: torch.device,
) -> tuple[RnaSeq11Adapter, dict[str, Any]]:
    checkpoint = load_adapter_checkpoint(checkpoint_path)
    embedding_resolution = int(checkpoint["embedding_resolution"])
    head_type = str(checkpoint.get("head_type", "linear"))
    head_resolutions = tuple(
        int(resolution)
        for resolution in checkpoint.get("head_resolutions", [embedding_resolution])
    )
    linear_head_architecture = str(
        checkpoint.get("linear_head_architecture", "conv1x1")
    )
    linear_hidden_channels = int(checkpoint.get("linear_hidden_channels", 256))
    linear_input_bottleneck_channels = optional_int(
        checkpoint.get("linear_input_bottleneck_channels")
    )
    residual_scale_init = float(checkpoint.get("linear_residual_scale_init", 1.0))
    linear_dilation = int(checkpoint.get("linear_dilation", 2))
    linear_kernel_size = int(checkpoint.get("linear_kernel_size", 15))
    residual_base_checkpoint_path = checkpoint.get("residual_base_checkpoint")
    linear_residual_fusion = str(checkpoint.get("linear_residual_fusion", "none"))
    residual_correction_scale_init = float(
        checkpoint.get("residual_correction_scale_init", 0.01)
    )
    residual_base_gate = str(checkpoint.get("residual_base_gate", "none"))
    residual_base_gate_center = float(checkpoint.get("residual_base_gate_center", 0.5))
    residual_base_gate_sharpness = float(
        checkpoint.get("residual_base_gate_sharpness", 4.0)
    )
    residual_base_gate_floor = float(checkpoint.get("residual_base_gate_floor", 0.0))
    organism_index = int(checkpoint["organism_index"])

    base_model = AlphaGenome.from_pretrained(weights, device=device)
    checkpoint_head_state = checkpoint["adapter_head_state_dict"]
    track_means = checkpoint_head_state.get("track_means")

    residual_base_head = None
    residual_base_input_bottleneck = None
    residual_base_resolution = 128
    if residual_base_checkpoint_path is not None:
        residual_base_checkpoint = load_adapter_checkpoint(residual_base_checkpoint_path)
        residual_base_resolution = int(residual_base_checkpoint["embedding_resolution"])
        residual_base_probe = RnaSeq11Adapter(
            base_model,
            n_tracks=n_tracks,
            embedding_resolution=residual_base_resolution,
            organism_index=organism_index,
            head_type="linear",
            head_resolutions=(residual_base_resolution,),
            linear_head_architecture=str(
                residual_base_checkpoint.get("linear_head_architecture", "conv1x1")
            ),
            linear_hidden_channels=int(
                residual_base_checkpoint.get("linear_hidden_channels", 256)
            ),
            linear_input_bottleneck_channels=optional_int(
                residual_base_checkpoint.get("linear_input_bottleneck_channels")
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
        embedding_resolution=embedding_resolution,
        organism_index=organism_index,
        head_type=head_type,
        head_resolutions=parse_resolutions(
            head_resolutions,
            fallback_resolution=embedding_resolution,
        ),
        linear_head_architecture=linear_head_architecture,
        linear_hidden_channels=linear_hidden_channels,
        linear_input_bottleneck_channels=linear_input_bottleneck_channels,
        linear_track_means=checkpoint_head_state.get("track_means"),
        linear_residual_scale_init=residual_scale_init,
        linear_dilation=linear_dilation,
        linear_kernel_size=linear_kernel_size,
        linear_residual_base_head=residual_base_head,
        linear_residual_base_input_bottleneck=residual_base_input_bottleneck,
        linear_residual_base_resolution=residual_base_resolution,
        linear_residual_fusion=linear_residual_fusion,
        linear_residual_correction_scale_init=residual_correction_scale_init,
        linear_residual_base_gate=residual_base_gate,
        linear_residual_base_gate_center=residual_base_gate_center,
        linear_residual_base_gate_sharpness=residual_base_gate_sharpness,
        linear_residual_base_gate_floor=residual_base_gate_floor,
        linear_output_calibration=str(checkpoint.get("linear_output_calibration", "none")),
        linear_output_calibration_gate_center=float(
            checkpoint.get("calibration_gate_center", 3.0)
        ),
        linear_output_calibration_gate_sharpness=float(
            checkpoint.get("calibration_gate_sharpness", 3.0)
        ),
        track_means=track_means,
    ).to(device)
    model.load_adapter_head_state_dict(checkpoint_head_state)
    model.eval()
    model.base_model.eval()
    return model, checkpoint


def transform_prediction(values: np.ndarray, *, transform: str, clip_min: float) -> np.ndarray:
    values = values.astype(np.float32, copy=False)
    if transform == "expm1":
        values = np.expm1(values).astype(np.float32, copy=False)
    elif transform != "as-is":
        raise ValueError(f"Unsupported prediction transform: {transform}")
    values = np.nan_to_num(values, nan=0.0, posinf=0.0, neginf=0.0)
    if not np.isnan(clip_min):
        values = np.maximum(values, np.float32(clip_min))
    return values.astype(np.float32, copy=False)


def ensure_accumulator(
    accumulators: dict[str, dict[str, np.ndarray]],
    *,
    chrom: str,
    chrom_sizes: dict[str, int],
    n_tracks: int,
) -> dict[str, np.ndarray]:
    if chrom in accumulators:
        return accumulators[chrom]
    length = chrom_sizes[chrom]
    accumulators[chrom] = {
        "sum": np.zeros((n_tracks, length), dtype=np.float32),
        "count": np.zeros((n_tracks, length), dtype=np.uint16),
    }
    return accumulators[chrom]


def write_bigwig(
    path: Path,
    *,
    chrom_header: list[tuple[str, int]],
    chrom_values: dict[str, np.ndarray],
    chunk_size: int,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = path.with_suffix(path.suffix + ".tmp")
    if tmp_path.exists():
        tmp_path.unlink()
    out = pyBigWig.open(str(tmp_path), "w")
    try:
        out.addHeader(chrom_header)
        for chrom, _length in chrom_header:
            values = chrom_values.get(chrom)
            if values is None:
                continue
            for start in range(0, len(values), chunk_size):
                chunk = values[start:start + chunk_size]
                out.addEntries(
                    chrom,
                    start,
                    values=chunk.tolist(),
                    span=1,
                    step=1,
                )
    finally:
        out.close()
    tmp_path.replace(path)


def export_split(
    *,
    split: str,
    dataset_dir: Path,
    output_dir: Path,
    model: RnaSeq11Adapter,
    device: torch.device,
    chrom_header: list[tuple[str, int]],
    track_rows: list[dict[str, str]],
    batch_size: int,
    num_workers: int,
    max_examples: int | None,
    prediction_transform: str,
    clip_min: float,
    value_chunk_size: int,
) -> dict[str, object]:
    if batch_size != 1:
        raise ValueError("BigWig export currently requires --batch-size 1")

    dataloader = make_dataloader(
        dataset_dir,
        batch_size=batch_size,
        shuffle=False,
        num_workers=num_workers,
        target_transform="none",
        max_examples=max_examples,
    )
    n_tracks = int(dataloader.dataset.metadata["n_tracks"])
    chrom_sizes = dict(chrom_header)
    accumulators: dict[str, dict[str, np.ndarray]] = {}
    n_batches = 0
    n_examples = 0

    with torch.inference_mode():
        for batch in dataloader:
            dna_sequence = batch["dna_sequence"].to(device=device, dtype=torch.float32)
            rna_seq = batch["rna_seq"]
            mask = batch["rna_seq_mask"].numpy().astype(bool)
            prediction = model(dna_sequence, target_length=rna_seq.shape[-1])
            if isinstance(prediction, dict):
                raise NotImplementedError(
                    "BigWig export currently supports linear full-length predictions."
                )
            prediction_np = prediction.detach().cpu().numpy()
            batch_size_observed = int(prediction_np.shape[0])
            for batch_index in range(batch_size_observed):
                chrom = str(batch["interval_chromosome"][batch_index])
                start = int(batch["interval_start"][batch_index])
                end = int(batch["interval_end"][batch_index])
                if chrom not in chrom_sizes:
                    raise ValueError(f"Chromosome {chrom!r} missing from chrom sizes")
                if end > chrom_sizes[chrom]:
                    raise ValueError(
                        f"Interval {chrom}:{start}-{end} exceeds chromosome length "
                        f"{chrom_sizes[chrom]}"
                    )
                span = end - start
                values = transform_prediction(
                    prediction_np[batch_index, :, :span],
                    transform=prediction_transform,
                    clip_min=clip_min,
                )
                accumulator = ensure_accumulator(
                    accumulators,
                    chrom=chrom,
                    chrom_sizes=chrom_sizes,
                    n_tracks=n_tracks,
                )
                for track_index in range(n_tracks):
                    if not bool(mask[batch_index, track_index].all()):
                        continue
                    accumulator["sum"][track_index, start:end] += values[track_index]
                    accumulator["count"][track_index, start:end] += 1
            n_batches += 1
            n_examples += batch_size_observed
            print(
                f"export_progress\t{split}\tbatch={n_batches}\texamples={n_examples}",
                flush=True,
            )

    split_dir = output_dir / split
    split_dir.mkdir(parents=True, exist_ok=True)
    output_rows = []
    for track_index, track_row in enumerate(track_rows):
        group_id = track_row.get("group_id", f"RNA_SEQ_{track_index + 1:03d}")
        track_name = track_row.get("name", group_id)
        chrom_values: dict[str, np.ndarray] = {}
        covered_bases = 0
        for chrom, accumulator in accumulators.items():
            counts = accumulator["count"][track_index]
            covered = counts > 0
            averaged = np.zeros_like(accumulator["sum"][track_index])
            averaged[covered] = accumulator["sum"][track_index][covered] / counts[covered]
            chrom_values[chrom] = averaged
            covered_bases += int(covered.sum())
        bw_path = split_dir / f"{group_id}.pred.{prediction_transform}.bw"
        write_bigwig(
            bw_path,
            chrom_header=chrom_header,
            chrom_values=chrom_values,
            chunk_size=value_chunk_size,
        )
        output_rows.append(
            {
                "split": split,
                "track_index": track_index,
                "group_id": group_id,
                "track_name": track_name,
                "output_path": str(bw_path),
                "covered_bases": covered_bases,
            }
        )
        print(f"wrote_bigwig\t{bw_path}", flush=True)

    manifest_path = split_dir / "bigwig_manifest.tsv"
    with manifest_path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "split",
                "track_index",
                "group_id",
                "track_name",
                "output_path",
                "covered_bases",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(output_rows)

    return {
        "split": split,
        "dataset_dir": str(dataset_dir),
        "n_batches": n_batches,
        "n_examples": n_examples,
        "chromosomes": sorted(accumulators),
        "bigwigs": [row["output_path"] for row in output_rows],
        "manifest": str(manifest_path),
    }


def main() -> None:
    args = parse_args()
    checkpoint_path = Path(args.checkpoint)
    weights = Path(args.weights)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    chrom_header = read_chrom_sizes(args.chrom_sizes)
    device = pick_device(args.device)

    split_to_dataset = {
        "valid": Path(args.valid_dataset_dir),
        "test": Path(args.test_dataset_dir),
    }
    splits = [part.strip() for part in args.splits.split(",") if part.strip()]
    invalid = sorted(set(splits) - set(split_to_dataset))
    if invalid:
        raise ValueError(f"Unsupported split(s): {invalid}")
    if args.value_chunk_size < 1:
        raise ValueError("--value-chunk-size must be positive")

    first_dataset = make_dataloader(
        split_to_dataset[splits[0]],
        batch_size=1,
        shuffle=False,
        num_workers=0,
        target_transform="none",
        max_examples=1,
    ).dataset
    n_tracks = int(first_dataset.metadata["n_tracks"])
    track_metadata = first_dataset.metadata["track_metadata"]
    track_rows = read_track_labels(track_metadata, n_tracks)

    print(f"checkpoint\t{checkpoint_path}")
    print(f"weights\t{weights}")
    print(f"output_dir\t{output_dir}")
    print(f"device\t{device}")
    print(f"prediction_transform\t{args.prediction_transform}")
    print(f"clip_min\t{args.clip_min}")
    print(f"splits\t{','.join(splits)}")

    model, checkpoint = build_model_from_checkpoint(
        checkpoint_path=checkpoint_path,
        weights=weights,
        n_tracks=n_tracks,
        device=device,
    )
    print(f"checkpoint_target_transform\t{checkpoint.get('target_transform', 'log1p')}")
    print(f"checkpoint_linear_target_space\t{checkpoint.get('linear_target_space', '')}")
    print(f"head_type\t{checkpoint.get('head_type', 'linear')}")
    print(f"embedding_resolution\t{checkpoint['embedding_resolution']}")
    print(f"base_parameters\t{count_parameters(model.base_model)}")
    print(f"total_parameters\t{count_parameters(model)}")
    print(f"head_parameters\t{count_parameters(model.head)}")

    split_summaries = []
    for split in splits:
        split_summaries.append(
            export_split(
                split=split,
                dataset_dir=split_to_dataset[split],
                output_dir=output_dir,
                model=model,
                device=device,
                chrom_header=chrom_header,
                track_rows=track_rows,
                batch_size=args.batch_size,
                num_workers=args.num_workers,
                max_examples=args.max_examples,
                prediction_transform=args.prediction_transform,
                clip_min=float(args.clip_min),
                value_chunk_size=args.value_chunk_size,
            )
        )

    summary = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "checkpoint": str(checkpoint_path),
        "weights": str(weights),
        "output_dir": str(output_dir),
        "device": str(device),
        "prediction_transform": args.prediction_transform,
        "clip_min": args.clip_min,
        "merge_method": "mean_over_overlapping_windows",
        "chrom_sizes": str(args.chrom_sizes),
        "track_metadata": str(track_metadata),
        "splits": split_summaries,
    }
    summary_path = output_dir / "export_summary.json"
    with summary_path.open("w") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)
        handle.write("\n")
    print(f"summary\t{summary_path}")


if __name__ == "__main__":
    main()
