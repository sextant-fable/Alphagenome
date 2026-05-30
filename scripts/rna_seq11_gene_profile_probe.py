#!/usr/bin/env python3
"""Gene-profile representation probes for RNA-seq11 checkpoints.

This script intentionally uses only train and valid splits. It never reads the
held-out test split. For each candidate checkpoint, it computes gene-level
predicted and true 11-track expression profiles on train chromosomes and valid
chromosome V, then evaluates representation-oriented utility metrics:

- gene-level Pearson/Spearman/MSE
- top expressed gene retrieval
- intestine T1/T3/T4 retrieval and gene metrics
- cross-track correlation structure preservation
- lightweight linear probes trained on train chromosomes and evaluated on valid
  genes for family, time, and dominant-track labels.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import torch

from alphagenome_pytorch import AlphaGenome
from alphagenome_rna_seq11_adapter import (
    RnaSeq11Adapter,
    pearson_from_sums,
    pick_device,
    rankdata_average,
)
from alphagenome_rna_seq11_eval import load_track_names
from rna_seq11_top15_extended_diagnostics import (
    Candidate,
    Gene,
    GeneOverlap,
    build_gene_overlaps,
    collect_top_candidates,
    compute_core_ends,
    load_gtf_annotations,
    make_model_from_checkpoint,
    read_manifest,
    read_tsv,
)
from torch_rna_seq_dataset import make_dataloader


FLOAT_FORMAT = ".8g"
FAMILY_NAMES = ("intestine", "muscle", "pharynx")
TIME_NAMES = ("T0", "T1", "T2", "T3", "T4")


@dataclass(frozen=True)
class SplitProfiles:
    genes: list[Gene]
    prediction_mean: np.ndarray
    target_mean: np.ndarray
    counts: np.ndarray


@dataclass(frozen=True)
class SharedSplitProfiles:
    genes: list[Gene]
    prediction_mean_by_model: np.ndarray
    target_mean: np.ndarray
    counts: np.ndarray


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--runs-dir", default="runs")
    parser.add_argument(
        "--candidate-tsv",
        default=None,
        help=(
            "Optional candidate TSV written by rna_seq11_top15_extended_diagnostics.py. "
            "When omitted, candidates are collected from --runs-dir."
        ),
    )
    parser.add_argument("--top-n", type=int, default=999)
    parser.add_argument(
        "--train-dataset-dir",
        default="alphagenome_custom/datasets/rna_seq_npz_train",
    )
    parser.add_argument(
        "--valid-dataset-dir",
        default="alphagenome_custom/datasets/rna_seq_npz_valid",
    )
    parser.add_argument(
        "--gtf",
        default="alphagenome_custom/reference/Caenorhabditis_elegans.WBcel235.115.gtf",
    )
    parser.add_argument(
        "--weights",
        default="weights/alphagenome_pytorch/model_all_folds.safetensors",
    )
    parser.add_argument(
        "--output-dir",
        default="runs/rna_seq11_gene_profile_probe_20260530",
    )
    parser.add_argument("--batch-size", type=int, default=1)
    parser.add_argument("--num-workers", type=int, default=0)
    parser.add_argument("--max-train-examples", type=int, default=None)
    parser.add_argument("--max-valid-examples", type=int, default=None)
    parser.add_argument("--gene-min-bases", type=int, default=128)
    parser.add_argument("--top-gene-fractions", default="0.01,0.05,0.10")
    parser.add_argument("--probe-steps", type=int, default=300)
    parser.add_argument("--probe-lr", type=float, default=0.05)
    parser.add_argument("--probe-weight-decay", type=float, default=1e-4)
    parser.add_argument("--probe-seed", type=int, default=20260530)
    parser.add_argument("--device", default="auto")
    parser.add_argument("--rank-shard", type=int, default=0)
    parser.add_argument("--num-shards", type=int, default=1)
    parser.add_argument(
        "--eval-mode",
        choices=["per-model", "shared-embeddings"],
        default="per-model",
        help=(
            "per-model preserves the original simple implementation. "
            "shared-embeddings computes each split's frozen AlphaGenome embeddings "
            "once per batch and evaluates all selected linear heads from them."
        ),
    )
    parser.add_argument(
        "--skip-existing",
        action="store_true",
        help="Skip run_ids already present in gene_profile_probe_metrics.tsv.",
    )
    return parser.parse_args()


def parse_float_list(value: str) -> list[float]:
    return [float(part.strip()) for part in value.split(",") if part.strip()]


def format_float(value: float) -> str:
    if math.isnan(value):
        return "nan"
    return format(float(value), FLOAT_FORMAT)


def candidate_from_row(row: dict[str, str]) -> Candidate:
    return Candidate(
        rank=int(row["rank"]),
        run_id=row["run_id"],
        run_dir=Path(row["checkpoint"]).parent,
        checkpoint=Path(row["checkpoint"]),
        source=row.get("source", ""),
        best_step=int(float(row["best_step"])) if row.get("best_step") else None,
        full_mse=float(row["full_mse"]),
        full_mae=float(row["full_mae"]) if row.get("full_mae") else None,
        pearson=float(row["pearson"]),
        common128_mse=float(row["common128_mse"]) if row.get("common128_mse") else None,
        common128_mae=float(row["common128_mae"]) if row.get("common128_mae") else None,
    )


def load_candidates(args: argparse.Namespace) -> list[Candidate]:
    if args.candidate_tsv:
        rows = read_tsv(Path(args.candidate_tsv))
        return [candidate_from_row(row) for row in rows[: args.top_n]]
    return collect_top_candidates(Path(args.runs_dir), args.top_n)


def chrom_lengths_from_manifest(rows: list[dict[str, str]]) -> dict[str, int]:
    chrom_lengths: dict[str, int] = {}
    for row in rows:
        chrom = row["chromosome"]
        chrom_lengths[chrom] = max(chrom_lengths.get(chrom, 0), int(row["end"]))
    return chrom_lengths


def load_split_context(
    dataset_dir: Path,
    *,
    gtf_path: Path,
    promoter_radius_bp: int,
    max_examples: int | None,
) -> tuple[list[dict[str, str]], list[int], list[Gene], list[list[GeneOverlap]]]:
    if "test" in dataset_dir.name.lower():
        raise ValueError("This script must not be run on the test split")
    manifest_rows = read_manifest(dataset_dir, max_examples=max_examples)
    core_ends = compute_core_ends(manifest_rows)
    chroms = {row["chromosome"] for row in manifest_rows}
    genes, _gene_intervals, _exon_intervals, _promoter_intervals = load_gtf_annotations(
        gtf_path,
        chroms=chroms,
        promoter_radius_bp=promoter_radius_bp,
        chrom_lengths=chrom_lengths_from_manifest(manifest_rows),
    )
    gene_overlaps = build_gene_overlaps(genes, manifest_rows, core_ends)
    return manifest_rows, core_ends, genes, gene_overlaps


def compute_gene_profiles(
    *,
    model: torch.nn.Module,
    dataloader: torch.utils.data.DataLoader,
    genes: list[Gene],
    gene_overlaps: list[list[GeneOverlap]],
    device: torch.device,
) -> SplitProfiles:
    n_tracks = int(dataloader.dataset.metadata["n_tracks"])
    pred_sum = np.zeros((len(genes), n_tracks), dtype=np.float64)
    target_sum = np.zeros((len(genes), n_tracks), dtype=np.float64)
    counts = np.zeros((len(genes), n_tracks), dtype=np.float64)

    if dataloader.batch_size is None:
        batch_size = 1
    else:
        batch_size = int(dataloader.batch_size)

    model.eval()
    with torch.inference_mode():
        for batch_index, batch in enumerate(dataloader):
            dna_sequence = batch["dna_sequence"].to(device=device, dtype=torch.float32)
            target = batch["rna_seq"].to(device=device, dtype=torch.float32)
            track_mask = batch["rna_seq_mask"].to(device=device, dtype=torch.bool)
            prediction = model(dna_sequence, target_length=target.shape[-1])
            if isinstance(prediction, dict):
                raise NotImplementedError("Gene-profile probes support linear predictions")
            pred_np = prediction.detach().cpu().numpy()
            target_np = target.detach().cpu().numpy()
            mask_np = track_mask.detach().cpu().numpy().astype(bool)
            for sample_index in range(int(prediction.shape[0])):
                window_index = batch_index * batch_size + sample_index
                if window_index >= len(gene_overlaps):
                    continue
                sample_mask = mask_np[sample_index].reshape(n_tracks)
                for overlap in gene_overlaps[window_index]:
                    if overlap.local_end <= overlap.local_start:
                        continue
                    segment = slice(overlap.local_start, overlap.local_end)
                    length = float(overlap.local_end - overlap.local_start)
                    if length <= 0:
                        continue
                    pred_values = pred_np[sample_index, :, segment].sum(axis=1)
                    target_values = target_np[sample_index, :, segment].sum(axis=1)
                    pred_sum[overlap.gene_index, sample_mask] += pred_values[sample_mask]
                    target_sum[overlap.gene_index, sample_mask] += target_values[sample_mask]
                    counts[overlap.gene_index, sample_mask] += length

    pred_mean = np.divide(
        pred_sum,
        counts,
        out=np.full_like(pred_sum, np.nan, dtype=np.float64),
        where=counts > 0,
    )
    target_mean = np.divide(
        target_sum,
        counts,
        out=np.full_like(target_sum, np.nan, dtype=np.float64),
        where=counts > 0,
    )
    return SplitProfiles(
        genes=genes,
        prediction_mean=pred_mean,
        target_mean=target_mean,
        counts=counts,
    )


def predict_from_shared_embeddings(
    model: RnaSeq11Adapter,
    embeddings: dict[str, torch.Tensor],
    *,
    target_length: int,
) -> torch.Tensor:
    """Match RnaSeq11Adapter.forward after its shared base_model.encode call."""
    if model.head_type != "linear":
        raise NotImplementedError("Shared embedding probes support linear heads only")
    base_prediction = None
    if model.linear_residual_base_head is not None:
        base_prediction = model._linear_prediction_from_embeddings(
            embeddings,
            resolution=model.linear_residual_base_resolution,
            prediction_head=model.linear_residual_base_head,
            input_bottleneck=model.linear_residual_base_input_bottleneck,
            target_length=target_length,
        )
    prediction = model._linear_prediction_from_embeddings(
        embeddings,
        resolution=model.embedding_resolution,
        prediction_head=model.head,
        input_bottleneck=model.linear_input_bottleneck,
        target_length=target_length,
        extra_head_input=(
            None
            if model.linear_residual_fusion == "none"
            else base_prediction.detach()
        ),
    )
    if base_prediction is not None:
        assert model.linear_residual_correction_scale is not None
        prediction = base_prediction.detach() + (
            model.linear_residual_correction_scale.to(dtype=prediction.dtype)
            * model._linear_residual_gate(base_prediction, prediction.dtype)
            * prediction
        )
    if model.linear_output_calibration is not None:
        prediction = model.linear_output_calibration(prediction)
    return prediction


def shared_encode_resolutions(models: list[RnaSeq11Adapter]) -> tuple[int, ...]:
    resolutions: set[int] = set()
    for model in models:
        resolutions.update(int(resolution) for resolution in model.head_resolutions)
        if model.linear_residual_base_head is not None:
            resolutions.add(int(model.linear_residual_base_resolution))
    return tuple(sorted(resolutions))


def accumulate_gene_sums(
    *,
    output_sum: np.ndarray,
    values_np: np.ndarray,
    mask_np: np.ndarray,
    gene_overlaps: list[list[GeneOverlap]],
    batch_index: int,
    batch_size: int,
) -> None:
    n_tracks = output_sum.shape[1]
    for sample_index in range(int(values_np.shape[0])):
        window_index = batch_index * batch_size + sample_index
        if window_index >= len(gene_overlaps):
            continue
        sample_mask = mask_np[sample_index].reshape(n_tracks)
        for overlap in gene_overlaps[window_index]:
            if overlap.local_end <= overlap.local_start:
                continue
            segment = slice(overlap.local_start, overlap.local_end)
            values = values_np[sample_index, :, segment].sum(axis=1)
            output_sum[overlap.gene_index, sample_mask] += values[sample_mask]


def accumulate_gene_counts(
    *,
    counts: np.ndarray,
    mask_np: np.ndarray,
    gene_overlaps: list[list[GeneOverlap]],
    batch_index: int,
    batch_size: int,
) -> None:
    n_tracks = counts.shape[1]
    for sample_index in range(int(mask_np.shape[0])):
        window_index = batch_index * batch_size + sample_index
        if window_index >= len(gene_overlaps):
            continue
        sample_mask = mask_np[sample_index].reshape(n_tracks)
        for overlap in gene_overlaps[window_index]:
            if overlap.local_end <= overlap.local_start:
                continue
            length = float(overlap.local_end - overlap.local_start)
            if length <= 0:
                continue
            counts[overlap.gene_index, sample_mask] += length


def compute_gene_profiles_shared(
    *,
    models: list[RnaSeq11Adapter],
    dataloader: torch.utils.data.DataLoader,
    genes: list[Gene],
    gene_overlaps: list[list[GeneOverlap]],
    device: torch.device,
) -> SharedSplitProfiles:
    if not models:
        raise ValueError("compute_gene_profiles_shared requires at least one model")
    organism_indices = {int(model.organism_index) for model in models}
    if len(organism_indices) != 1:
        raise ValueError("Shared embedding mode requires one organism_index")
    n_tracks = int(dataloader.dataset.metadata["n_tracks"])
    pred_sum = np.zeros((len(models), len(genes), n_tracks), dtype=np.float64)
    target_sum = np.zeros((len(genes), n_tracks), dtype=np.float64)
    counts = np.zeros((len(genes), n_tracks), dtype=np.float64)
    batch_size = 1 if dataloader.batch_size is None else int(dataloader.batch_size)
    encode_resolutions = shared_encode_resolutions(models)
    base_model = models[0].base_model
    base_model.eval()
    for model in models:
        model.eval()

    with torch.inference_mode():
        for batch_index, batch in enumerate(dataloader):
            dna_sequence = batch["dna_sequence"].to(device=device, dtype=torch.float32)
            target = batch["rna_seq"].to(device=device, dtype=torch.float32)
            track_mask = batch["rna_seq_mask"].to(device=device, dtype=torch.bool)
            mask_np = track_mask.detach().cpu().numpy().astype(bool)
            target_np = target.detach().cpu().numpy()
            accumulate_gene_sums(
                output_sum=target_sum,
                values_np=target_np,
                mask_np=mask_np,
                gene_overlaps=gene_overlaps,
                batch_index=batch_index,
                batch_size=batch_size,
            )
            accumulate_gene_counts(
                counts=counts,
                mask_np=mask_np,
                gene_overlaps=gene_overlaps,
                batch_index=batch_index,
                batch_size=batch_size,
            )
            organism_index = models[0].organism_index_tensor(
                dna_sequence.shape[0],
                device,
            )
            embeddings = base_model.encode(
                dna_sequence.transpose(1, 2).contiguous(),
                organism_index,
                resolutions=encode_resolutions,
                channels_last=False,
            )
            for model_index, model in enumerate(models):
                prediction = predict_from_shared_embeddings(
                    model,
                    embeddings,
                    target_length=target.shape[-1],
                )
                pred_np = prediction.detach().cpu().numpy()
                accumulate_gene_sums(
                    output_sum=pred_sum[model_index],
                    values_np=pred_np,
                    mask_np=mask_np,
                    gene_overlaps=gene_overlaps,
                    batch_index=batch_index,
                    batch_size=batch_size,
                )
                del prediction
            del embeddings

    pred_mean = np.divide(
        pred_sum,
        counts[None, :, :],
        out=np.full_like(pred_sum, np.nan, dtype=np.float64),
        where=counts[None, :, :] > 0,
    )
    target_mean = np.divide(
        target_sum,
        counts,
        out=np.full_like(target_sum, np.nan, dtype=np.float64),
        where=counts > 0,
    )
    return SharedSplitProfiles(
        genes=genes,
        prediction_mean_by_model=pred_mean,
        target_mean=target_mean,
        counts=counts,
    )


def finite_gene_mask(profiles: SplitProfiles, min_bases: int) -> np.ndarray:
    return finite_gene_mask_arrays(
        profiles.prediction_mean,
        profiles.target_mean,
        profiles.counts,
        min_bases,
    )


def finite_gene_mask_arrays(
    prediction_mean: np.ndarray,
    target_mean: np.ndarray,
    counts: np.ndarray,
    min_bases: int,
) -> np.ndarray:
    return (
        np.all(np.isfinite(prediction_mean), axis=1)
        & np.all(np.isfinite(target_mean), axis=1)
        & np.all(counts >= float(min_bases), axis=1)
    )


def spearman(x: np.ndarray, y: np.ndarray) -> float:
    mask = np.isfinite(x) & np.isfinite(y)
    if int(mask.sum()) <= 1:
        return math.nan
    xr = rankdata_average(x[mask])
    yr = rankdata_average(y[mask])
    return pearson_from_sums(
        float(xr.sum()),
        float(yr.sum()),
        float(np.square(xr).sum()),
        float(np.square(yr).sum()),
        float((xr * yr).sum()),
        float(len(xr)),
    )


def pearson(x: np.ndarray, y: np.ndarray) -> float:
    mask = np.isfinite(x) & np.isfinite(y)
    if int(mask.sum()) <= 1:
        return math.nan
    x64 = x[mask].astype(np.float64, copy=False)
    y64 = y[mask].astype(np.float64, copy=False)
    return pearson_from_sums(
        float(x64.sum()),
        float(y64.sum()),
        float(np.square(x64).sum()),
        float(np.square(y64).sum()),
        float((x64 * y64).sum()),
        float(len(x64)),
    )


def mse(x: np.ndarray, y: np.ndarray) -> float:
    mask = np.isfinite(x) & np.isfinite(y)
    if int(mask.sum()) == 0:
        return math.nan
    return float(np.mean(np.square(x[mask] - y[mask])))


def profile_labels(values: np.ndarray) -> dict[str, np.ndarray]:
    intestine = np.nanmean(values[:, 0:5], axis=1)
    muscle = np.nanmean(values[:, 5:10], axis=1)
    pharynx = values[:, 10]
    family = np.nanargmax(np.stack([intestine, muscle, pharynx], axis=1), axis=1)
    time_scores = np.stack(
        [
            np.nanmax(values[:, [0, 5]], axis=1),
            np.nanmax(values[:, [1, 6]], axis=1),
            np.nanmax(values[:, [2, 7]], axis=1),
            np.nanmax(values[:, [3, 8]], axis=1),
            np.nanmax(values[:, [4, 9, 10]], axis=1),
        ],
        axis=1,
    )
    time = np.nanargmax(time_scores, axis=1)
    dominant_track = np.nanargmax(values, axis=1)
    return {"family": family, "time": time, "track": dominant_track}


def balanced_accuracy(y_true: np.ndarray, y_pred: np.ndarray, n_classes: int) -> float:
    recalls: list[float] = []
    for label in range(n_classes):
        mask = y_true == label
        if not bool(mask.any()):
            continue
        recalls.append(float(np.mean(y_pred[mask] == label)))
    if not recalls:
        return math.nan
    return float(np.mean(recalls))


def standardize_train_valid(
    train_x: np.ndarray,
    valid_x: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    mean = np.nanmean(train_x, axis=0)
    std = np.nanstd(train_x, axis=0)
    std = np.where(std < 1e-6, 1.0, std)
    return (train_x - mean) / std, (valid_x - mean) / std


def fit_linear_probe(
    *,
    train_x: np.ndarray,
    train_y: np.ndarray,
    valid_x: np.ndarray,
    valid_y: np.ndarray,
    n_classes: int,
    steps: int,
    lr: float,
    weight_decay: float,
    seed: int,
) -> dict[str, float]:
    train_mask = np.all(np.isfinite(train_x), axis=1) & np.isfinite(train_y)
    valid_mask = np.all(np.isfinite(valid_x), axis=1) & np.isfinite(valid_y)
    train_x = train_x[train_mask]
    valid_x = valid_x[valid_mask]
    train_y = train_y[train_mask].astype(np.int64)
    valid_y = valid_y[valid_mask].astype(np.int64)
    present = sorted(set(int(value) for value in train_y.tolist()))
    if len(train_y) == 0 or len(valid_y) == 0 or len(present) < 2:
        return {
            "accuracy": math.nan,
            "balanced_accuracy": math.nan,
            "n_train": float(len(train_y)),
            "n_valid": float(len(valid_y)),
        }
    train_x, valid_x = standardize_train_valid(train_x, valid_x)
    torch.manual_seed(seed)
    model = torch.nn.Linear(train_x.shape[1], n_classes)
    optimizer = torch.optim.AdamW(model.parameters(), lr=lr, weight_decay=weight_decay)
    x_train_t = torch.as_tensor(train_x, dtype=torch.float32)
    y_train_t = torch.as_tensor(train_y, dtype=torch.long)
    for _step in range(steps):
        optimizer.zero_grad(set_to_none=True)
        loss = torch.nn.functional.cross_entropy(model(x_train_t), y_train_t)
        loss.backward()
        optimizer.step()
    with torch.inference_mode():
        pred = model(torch.as_tensor(valid_x, dtype=torch.float32)).argmax(dim=1).numpy()
    return {
        "accuracy": float(np.mean(pred == valid_y)),
        "balanced_accuracy": balanced_accuracy(valid_y, pred, n_classes),
        "n_train": float(len(train_y)),
        "n_valid": float(len(valid_y)),
    }


def top_retrieval_metrics(
    pred: np.ndarray,
    target: np.ndarray,
    *,
    fractions: list[float],
) -> list[dict[str, float]]:
    rows: list[dict[str, float]] = []
    for track_index in range(pred.shape[1]):
        mask = np.isfinite(pred[:, track_index]) & np.isfinite(target[:, track_index])
        pred_values = pred[mask, track_index]
        target_values = target[mask, track_index]
        n = len(pred_values)
        if n == 0:
            continue
        for fraction in fractions:
            k = max(1, int(math.ceil(n * fraction)))
            true_top = set(np.argpartition(target_values, -k)[-k:].tolist())
            pred_top = set(np.argpartition(pred_values, -k)[-k:].tolist())
            overlap = len(true_top & pred_top)
            rows.append(
                {
                    "track_index": float(track_index),
                    "fraction": float(fraction),
                    "k": float(k),
                    "recall": float(overlap / k),
                    "precision": float(overlap / k),
                    "spearman": spearman(pred_values, target_values),
                    "pearson": pearson(pred_values, target_values),
                    "mse": mse(pred_values, target_values),
                }
            )
    return rows


def summarize_retrieval(rows: list[dict[str, float]], fraction: float) -> dict[str, float]:
    selected = [row for row in rows if abs(row["fraction"] - fraction) < 1e-12]
    if not selected:
        return {"recall": math.nan, "spearman": math.nan}
    return {
        "recall": float(np.mean([row["recall"] for row in selected])),
        "spearman": float(np.mean([row["spearman"] for row in selected])),
    }


def summarize_intestine_t134(rows: list[dict[str, float]], fraction: float) -> dict[str, float]:
    selected = [
        row
        for row in rows
        if abs(row["fraction"] - fraction) < 1e-12
        and int(row["track_index"]) in {1, 3, 4}
    ]
    if not selected:
        return {"recall": math.nan, "spearman": math.nan, "mse": math.nan}
    return {
        "recall": float(np.mean([row["recall"] for row in selected])),
        "spearman": float(np.mean([row["spearman"] for row in selected])),
        "mse": float(np.mean([row["mse"] for row in selected])),
    }


def track_structure_metrics(pred: np.ndarray, target: np.ndarray) -> dict[str, float]:
    mask = np.all(np.isfinite(pred), axis=1) & np.all(np.isfinite(target), axis=1)
    if int(mask.sum()) <= 2:
        return {"corr_pearson": math.nan, "corr_spearman": math.nan, "corr_mse": math.nan}
    pred_corr = np.corrcoef(pred[mask].T)
    target_corr = np.corrcoef(target[mask].T)
    tri = np.triu_indices(pred_corr.shape[0], k=1)
    pred_vec = pred_corr[tri]
    target_vec = target_corr[tri]
    return {
        "corr_pearson": pearson(pred_vec, target_vec),
        "corr_spearman": spearman(pred_vec, target_vec),
        "corr_mse": mse(pred_vec, target_vec),
    }


def write_rows(path: Path, fieldnames: list[str], rows: list[dict[str, object]]) -> None:
    exists = path.exists()
    with path.open("a", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        if not exists:
            writer.writeheader()
        for row in rows:
            writer.writerow(row)


def read_existing_run_ids(path: Path) -> set[str]:
    if not path.exists():
        return set()
    with path.open() as handle:
        return {row["run_id"] for row in csv.DictReader(handle, delimiter="\t")}


def evaluate_candidate_profiles(
    *,
    candidate: Candidate,
    train_profiles: SplitProfiles,
    valid_profiles: SplitProfiles,
    track_labels: list[str],
    args: argparse.Namespace,
    fractions: list[float],
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    train_mask = finite_gene_mask(train_profiles, args.gene_min_bases)
    valid_mask = finite_gene_mask(valid_profiles, args.gene_min_bases)
    train_pred = train_profiles.prediction_mean[train_mask]
    train_target = train_profiles.target_mean[train_mask]
    valid_pred = valid_profiles.prediction_mean[valid_mask]
    valid_target = valid_profiles.target_mean[valid_mask]

    retrieval_rows_raw = top_retrieval_metrics(
        valid_pred,
        valid_target,
        fractions=fractions,
    )
    per_track_rows: list[dict[str, object]] = []
    for row in retrieval_rows_raw:
        track_index = int(row["track_index"])
        per_track_rows.append(
            {
                "rank": candidate.rank,
                "run_id": candidate.run_id,
                "track_index": track_index,
                "track_name": track_labels[track_index],
                "top_fraction": format_float(row["fraction"]),
                "k": int(row["k"]),
                "top_gene_recall": format_float(row["recall"]),
                "top_gene_precision": format_float(row["precision"]),
                "gene_pearson": format_float(row["pearson"]),
                "gene_spearman": format_float(row["spearman"]),
                "gene_mse": format_float(row["mse"]),
            }
        )

    labels_train = profile_labels(train_target)
    labels_valid = profile_labels(valid_target)
    labels_pred_valid = profile_labels(valid_pred)
    direct_family_acc = float(np.mean(labels_pred_valid["family"] == labels_valid["family"]))
    direct_time_acc = float(np.mean(labels_pred_valid["time"] == labels_valid["time"]))
    direct_track_acc = float(np.mean(labels_pred_valid["track"] == labels_valid["track"]))
    family_probe = fit_linear_probe(
        train_x=train_pred,
        train_y=labels_train["family"],
        valid_x=valid_pred,
        valid_y=labels_valid["family"],
        n_classes=3,
        steps=args.probe_steps,
        lr=args.probe_lr,
        weight_decay=args.probe_weight_decay,
        seed=args.probe_seed + candidate.rank,
    )
    time_probe = fit_linear_probe(
        train_x=train_pred,
        train_y=labels_train["time"],
        valid_x=valid_pred,
        valid_y=labels_valid["time"],
        n_classes=5,
        steps=args.probe_steps,
        lr=args.probe_lr,
        weight_decay=args.probe_weight_decay,
        seed=args.probe_seed + 1000 + candidate.rank,
    )
    track_probe = fit_linear_probe(
        train_x=train_pred,
        train_y=labels_train["track"],
        valid_x=valid_pred,
        valid_y=labels_valid["track"],
        n_classes=11,
        steps=args.probe_steps,
        lr=args.probe_lr,
        weight_decay=args.probe_weight_decay,
        seed=args.probe_seed + 2000 + candidate.rank,
    )
    structure = track_structure_metrics(valid_pred, valid_target)
    top01 = summarize_retrieval(retrieval_rows_raw, 0.01)
    top05 = summarize_retrieval(retrieval_rows_raw, 0.05)
    top10 = summarize_retrieval(retrieval_rows_raw, 0.10)
    intestine_top05 = summarize_intestine_t134(retrieval_rows_raw, 0.05)

    summary_rows = [
        {
            "rank": candidate.rank,
            "run_id": candidate.run_id,
            "checkpoint": str(candidate.checkpoint),
            "source": candidate.source,
            "best_step": "" if candidate.best_step is None else candidate.best_step,
            "full_mse": format_float(candidate.full_mse),
            "full_mae": "" if candidate.full_mae is None else format_float(candidate.full_mae),
            "pearson": format_float(candidate.pearson),
            "common128_mse": ""
            if candidate.common128_mse is None
            else format_float(candidate.common128_mse),
            "train_genes": int(train_mask.sum()),
            "valid_genes": int(valid_mask.sum()),
            "valid_top01_recall": format_float(top01["recall"]),
            "valid_top05_recall": format_float(top05["recall"]),
            "valid_top10_recall": format_float(top10["recall"]),
            "valid_top05_gene_spearman": format_float(top05["spearman"]),
            "intestine_t134_top05_recall": format_float(intestine_top05["recall"]),
            "intestine_t134_gene_spearman": format_float(intestine_top05["spearman"]),
            "intestine_t134_gene_mse": format_float(intestine_top05["mse"]),
            "track_corr_pearson": format_float(structure["corr_pearson"]),
            "track_corr_spearman": format_float(structure["corr_spearman"]),
            "track_corr_mse": format_float(structure["corr_mse"]),
            "direct_family_accuracy": format_float(direct_family_acc),
            "direct_time_accuracy": format_float(direct_time_acc),
            "direct_track_accuracy": format_float(direct_track_acc),
            "family_probe_accuracy": format_float(family_probe["accuracy"]),
            "family_probe_balanced_accuracy": format_float(
                family_probe["balanced_accuracy"]
            ),
            "time_probe_accuracy": format_float(time_probe["accuracy"]),
            "time_probe_balanced_accuracy": format_float(time_probe["balanced_accuracy"]),
            "track_probe_accuracy": format_float(track_probe["accuracy"]),
            "track_probe_balanced_accuracy": format_float(
                track_probe["balanced_accuracy"]
            ),
        }
    ]
    return summary_rows, per_track_rows


def evaluate_candidate_probe(
    *,
    candidate: Candidate,
    base_model: AlphaGenome,
    train_loader: torch.utils.data.DataLoader,
    valid_loader: torch.utils.data.DataLoader,
    train_genes: list[Gene],
    valid_genes: list[Gene],
    train_overlaps: list[list[GeneOverlap]],
    valid_overlaps: list[list[GeneOverlap]],
    track_labels: list[str],
    args: argparse.Namespace,
    device: torch.device,
    fractions: list[float],
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    model, _checkpoint = make_model_from_checkpoint(
        checkpoint_path=candidate.checkpoint,
        base_model=base_model,
        n_tracks=len(track_labels),
        device=device,
    )
    train_profiles = compute_gene_profiles(
        model=model,
        dataloader=train_loader,
        genes=train_genes,
        gene_overlaps=train_overlaps,
        device=device,
    )
    valid_profiles = compute_gene_profiles(
        model=model,
        dataloader=valid_loader,
        genes=valid_genes,
        gene_overlaps=valid_overlaps,
        device=device,
    )
    del model
    if device.type == "cuda":
        torch.cuda.empty_cache()
    return evaluate_candidate_profiles(
        candidate=candidate,
        train_profiles=train_profiles,
        valid_profiles=valid_profiles,
        track_labels=track_labels,
        args=args,
        fractions=fractions,
    )


def main() -> None:
    args = parse_args()
    train_dataset_dir = Path(args.train_dataset_dir)
    valid_dataset_dir = Path(args.valid_dataset_dir)
    if "test" in train_dataset_dir.name.lower() or "test" in valid_dataset_dir.name.lower():
        raise ValueError("This script must not read the test split")
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    candidates = load_candidates(args)
    selected = [
        candidate
        for index, candidate in enumerate(candidates)
        if index % args.num_shards == args.rank_shard
    ]
    if args.skip_existing:
        existing_run_ids = read_existing_run_ids(output_dir / "gene_profile_probe_metrics.tsv")
        selected = [
            candidate for candidate in selected if candidate.run_id not in existing_run_ids
        ]
    fractions = parse_float_list(args.top_gene_fractions)

    train_manifest, _train_core_ends, train_genes, train_overlaps = load_split_context(
        train_dataset_dir,
        gtf_path=Path(args.gtf),
        promoter_radius_bp=1000,
        max_examples=args.max_train_examples,
    )
    valid_manifest, _valid_core_ends, valid_genes, valid_overlaps = load_split_context(
        valid_dataset_dir,
        gtf_path=Path(args.gtf),
        promoter_radius_bp=1000,
        max_examples=args.max_valid_examples,
    )
    train_loader = make_dataloader(
        train_dataset_dir,
        batch_size=args.batch_size,
        shuffle=False,
        num_workers=args.num_workers,
        target_transform="log1p",
        max_examples=args.max_train_examples,
    )
    valid_loader = make_dataloader(
        valid_dataset_dir,
        batch_size=args.batch_size,
        shuffle=False,
        num_workers=args.num_workers,
        target_transform="log1p",
        max_examples=args.max_valid_examples,
    )
    n_tracks = int(valid_loader.dataset.metadata["n_tracks"])
    track_labels = load_track_names(valid_loader.dataset.metadata, n_tracks)
    device = pick_device(args.device)
    print(f"device\t{device}")
    print(f"selected_models\t{len(selected)}")
    print(f"train_examples\t{len(train_manifest)}")
    print(f"valid_examples\t{len(valid_manifest)}")
    print(f"train_genes\t{len(train_genes)}")
    print(f"valid_genes\t{len(valid_genes)}")
    if device.type == "cuda":
        print(f"cuda_device\t{torch.cuda.get_device_name(device)}")

    metadata = {
        "candidate_tsv": args.candidate_tsv,
        "runs_dir": args.runs_dir,
        "top_n": args.top_n,
        "rank_shard": args.rank_shard,
        "num_shards": args.num_shards,
        "train_dataset_dir": str(train_dataset_dir),
        "valid_dataset_dir": str(valid_dataset_dir),
        "test_split_used": False,
        "eval_mode": args.eval_mode,
        "skip_existing": bool(args.skip_existing),
        "top_gene_fractions": fractions,
        "probe_steps": args.probe_steps,
        "probe_lr": args.probe_lr,
        "probe_weight_decay": args.probe_weight_decay,
    }
    with (output_dir / "probe_metadata.json").open("w") as handle:
        json.dump(metadata, handle, indent=2, sort_keys=True)

    base_model = AlphaGenome.from_pretrained(Path(args.weights), device=device)
    summary_fieldnames = [
        "rank",
        "run_id",
        "checkpoint",
        "source",
        "best_step",
        "full_mse",
        "full_mae",
        "pearson",
        "common128_mse",
        "train_genes",
        "valid_genes",
        "valid_top01_recall",
        "valid_top05_recall",
        "valid_top10_recall",
        "valid_top05_gene_spearman",
        "intestine_t134_top05_recall",
        "intestine_t134_gene_spearman",
        "intestine_t134_gene_mse",
        "track_corr_pearson",
        "track_corr_spearman",
        "track_corr_mse",
        "direct_family_accuracy",
        "direct_time_accuracy",
        "direct_track_accuracy",
        "family_probe_accuracy",
        "family_probe_balanced_accuracy",
        "time_probe_accuracy",
        "time_probe_balanced_accuracy",
        "track_probe_accuracy",
        "track_probe_balanced_accuracy",
    ]
    per_track_fieldnames = [
        "rank",
        "run_id",
        "track_index",
        "track_name",
        "top_fraction",
        "k",
        "top_gene_recall",
        "top_gene_precision",
        "gene_pearson",
        "gene_spearman",
        "gene_mse",
    ]
    if args.eval_mode == "shared-embeddings" and selected:
        models: list[RnaSeq11Adapter] = []
        for candidate in selected:
            model, _checkpoint = make_model_from_checkpoint(
                checkpoint_path=candidate.checkpoint,
                base_model=base_model,
                n_tracks=len(track_labels),
                device=device,
            )
            models.append(model)
        train_shared = compute_gene_profiles_shared(
            models=models,
            dataloader=train_loader,
            genes=train_genes,
            gene_overlaps=train_overlaps,
            device=device,
        )
        valid_shared = compute_gene_profiles_shared(
            models=models,
            dataloader=valid_loader,
            genes=valid_genes,
            gene_overlaps=valid_overlaps,
            device=device,
        )
        for model_index, candidate in enumerate(selected):
            print(
                "evaluate\t"
                f"rank={candidate.rank}\t"
                f"mse={candidate.full_mse:.8f}\t"
                f"pearson={candidate.pearson:.8f}\t"
                f"run={candidate.run_id}"
            )
            train_profiles = SplitProfiles(
                genes=train_shared.genes,
                prediction_mean=train_shared.prediction_mean_by_model[model_index],
                target_mean=train_shared.target_mean,
                counts=train_shared.counts,
            )
            valid_profiles = SplitProfiles(
                genes=valid_shared.genes,
                prediction_mean=valid_shared.prediction_mean_by_model[model_index],
                target_mean=valid_shared.target_mean,
                counts=valid_shared.counts,
            )
            summary_rows, per_track_rows = evaluate_candidate_profiles(
                candidate=candidate,
                train_profiles=train_profiles,
                valid_profiles=valid_profiles,
                track_labels=track_labels,
                args=args,
                fractions=fractions,
            )
            write_rows(
                output_dir / "gene_profile_probe_metrics.tsv",
                summary_fieldnames,
                summary_rows,
            )
            write_rows(
                output_dir / "gene_profile_probe_per_track.tsv",
                per_track_fieldnames,
                per_track_rows,
            )
        print(f"output_dir\t{output_dir}")
        return

    for candidate in selected:
        print(
            "evaluate\t"
            f"rank={candidate.rank}\t"
            f"mse={candidate.full_mse:.8f}\t"
            f"pearson={candidate.pearson:.8f}\t"
            f"run={candidate.run_id}"
        )
        summary_rows, per_track_rows = evaluate_candidate_probe(
            candidate=candidate,
            base_model=base_model,
            train_loader=train_loader,
            valid_loader=valid_loader,
            train_genes=train_genes,
            valid_genes=valid_genes,
            train_overlaps=train_overlaps,
            valid_overlaps=valid_overlaps,
            track_labels=track_labels,
            args=args,
            device=device,
            fractions=fractions,
        )
        write_rows(output_dir / "gene_profile_probe_metrics.tsv", summary_fieldnames, summary_rows)
        write_rows(
            output_dir / "gene_profile_probe_per_track.tsv",
            per_track_fieldnames,
            per_track_rows,
        )
    print(f"output_dir\t{output_dir}")


if __name__ == "__main__":
    main()
