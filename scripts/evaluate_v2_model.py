#!/usr/bin/env python3
"""Evaluate one v2 checkpoint on fixed core-only blocked-CV subwindows."""

from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import fcntl
import hashlib
import json
import math
import os
from pathlib import Path

import numpy as np
import torch

from scripts import train_v2_model
from scripts import v2_training_components as components
from scripts import v2_validation_metrics as full_metrics
from scripts.v2_bigwig_dataset import V2BigWigDataset


REPO_ROOT = Path(__file__).resolve().parents[1]
METADATA_DIR = REPO_ROOT / "alphagenome_custom/metadata/v2"
FINAL_CLAIM_PATH = METADATA_DIR / "final_test_claim.json"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--model", choices=("A", "B", "C"), required=True)
    parser.add_argument("--training-loss", choices=("paper", "log1p_mse"), required=True)
    parser.add_argument("--fold", type=int, choices=range(0, 6), required=True)
    parser.add_argument("--checkpoint", required=True)
    parser.add_argument("--sequence-length", type=int, default=131072)
    parser.add_argument("--hidden-channels", type=int, default=64)
    parser.add_argument("--mean-column")
    parser.add_argument("--output", required=True)
    parser.add_argument("--device", default="cuda")
    parser.add_argument("--final-test", action="store_true")
    parser.add_argument("--final-test-execution-id")
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


def atomic_json(path: Path, payload: dict[str, object]) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)


def consume_final_test(
    checkpoint_sha256: str, execution_id: str, physical_gpu: int
) -> None:
    lock_path = METADATA_DIR / "final_test_lock.json"
    with FINAL_CLAIM_PATH.open() as claim_handle:
        fcntl.flock(claim_handle.fileno(), fcntl.LOCK_EX)
        claim = json.load(claim_handle)
        if not (
            claim.get("status") == "claimed"
            and claim.get("execution_id") == execution_id
            and claim.get("checkpoint_sha256") == checkpoint_sha256
            and claim.get("physical_gpu") == physical_gpu
        ):
            raise RuntimeError("Final-test claim is mismatched or inactive")
        lock = json.loads(lock_path.read_text())
        if (
            lock.get("superseded") is True
            or lock.get("checkpoint_sha256") != checkpoint_sha256
            or lock.get("test_consumed") is True
        ):
            raise RuntimeError("Final-test lock is mismatched or already consumed")
        lock["test_consumed"] = True
        lock["test_consumed_at"] = datetime.now(timezone.utc).replace(
            microsecond=0
        ).isoformat()
        lock["test_status"] = "running"
        lock["test_execution_id"] = execution_id
        lock["test_physical_gpu"] = physical_gpu
        atomic_json(lock_path, lock)


def finalize_final_test(
    checkpoint_sha256: str, execution_id: str, report_path: Path
) -> None:
    lock_path = METADATA_DIR / "final_test_lock.json"
    lock = json.loads(lock_path.read_text())
    if (
        lock.get("superseded") is True
        or
        lock.get("checkpoint_sha256") != checkpoint_sha256
        or lock.get("test_consumed") is not True
        or lock.get("test_status") != "running"
        or lock.get("test_execution_id") != execution_id
    ):
        raise RuntimeError("Final-test lock changed during evaluation")
    lock["test_status"] = "completed"
    lock["final_report_path"] = str(report_path.relative_to(REPO_ROOT))
    lock["final_report_sha256"] = sha256(report_path)
    atomic_json(lock_path, lock)


def core_subwindows(
    intervals: list[dict[str, str]], sequence_length: int
) -> tuple[list[tuple[int, int]], int]:
    if sequence_length <= 0 or sequence_length % 128:
        raise ValueError("sequence_length must be positive and divisible by 128")
    result = []
    eligible_bases = 0
    for index, row in enumerate(intervals):
        core_start = ((int(row["core_start"]) + 127) // 128) * 128
        core_end = (int(row["core_end"]) // 128) * 128
        eligible_bases += max(0, core_end - core_start)
        for start in range(core_start, core_end - sequence_length + 1, sequence_length):
            result.append((index, start - int(row["start"])))
    if not result:
        raise RuntimeError("No complete core-only validation subwindows")
    return result, eligible_bases


def load_checkpoint(
    model: torch.nn.Module,
    checkpoint_path: Path,
    *,
    expected_model: str,
    expected_loss: str,
    expected_fold: int,
    expected_sequence_length: int,
    expected_hidden_channels: int,
    expected_mean_column: str,
) -> dict[str, object]:
    checkpoint = torch.load(checkpoint_path, map_location="cpu", weights_only=False)
    expected = {
        "model_id": expected_model,
        "loss": expected_loss,
        "fold": expected_fold,
        "sequence_length": expected_sequence_length,
        "hidden_channels": expected_hidden_channels,
        "mean_column": expected_mean_column,
    }
    mismatches = {
        key: (checkpoint.get(key), value)
        for key, value in expected.items()
        if checkpoint.get(key) != value
    }
    if mismatches:
        raise RuntimeError(f"Checkpoint configuration mismatch: {mismatches}")
    trainable_names = {
        name for name, parameter in model.named_parameters() if parameter.requires_grad
    }
    checkpoint_state = checkpoint.get("trainable_model_state", {})
    if set(checkpoint_state) != trainable_names:
        raise RuntimeError("Checkpoint trainable parameter set mismatch")
    state = model.state_dict()
    for name, value in checkpoint_state.items():
        if name not in state or state[name].shape != value.shape:
            raise RuntimeError(f"Checkpoint tensor mismatch: {name}")
        state[name] = value
    model.load_state_dict(state, strict=True)
    return checkpoint


def pearson_from_sums(
    count: int,
    sum_x: np.ndarray,
    sum_y: np.ndarray,
    sum_xx: np.ndarray,
    sum_yy: np.ndarray,
    sum_xy: np.ndarray,
) -> np.ndarray:
    numerator = count * sum_xy - sum_x * sum_y
    denominator = np.sqrt(
        np.maximum(count * sum_xx - np.square(sum_x), 0)
        * np.maximum(count * sum_yy - np.square(sum_y), 0)
    )
    result = np.full(sum_x.shape, np.nan, dtype=np.float64)
    valid = denominator > 0
    result[valid] = numerator[valid] / denominator[valid]
    return result


def stats_record(
    stats: full_metrics.PerTrackStats, group_ids: list[str]
) -> dict[str, object]:
    values = stats.per_track()
    return {
        "positions_per_track": full_metrics.keyed(values["n"], group_ids),
        "mean_per_track_mse": full_metrics.finite_mean(values["mse"]),
        "mean_per_track_mae": full_metrics.finite_mean(values["mae"]),
        "mean_per_track_pearson": full_metrics.finite_mean(values["pearson"]),
        "finite_per_track_pearson": int(np.isfinite(values["pearson"]).sum()),
        "per_track_mse": full_metrics.keyed(values["mse"], group_ids),
        "per_track_mae": full_metrics.keyed(values["mae"], group_ids),
        "per_track_pearson": full_metrics.keyed(values["pearson"], group_ids),
    }


def main() -> None:
    args = parse_args()
    if args.sequence_length % 128:
        raise ValueError("sequence-length must be divisible by 128")
    device = torch.device(args.device)
    visible = os.environ.get("CUDA_VISIBLE_DEVICES", "")
    if device.type != "cuda" or not torch.cuda.is_available():
        raise RuntimeError("Registered P6B evaluation requires CUDA")
    if visible not in {"2", "3"}:
        raise RuntimeError("Evaluation must use exactly one physical GPU 2 or 3")
    if (args.fold == 0) != args.final_test:
        raise RuntimeError("fold 0 is reserved exclusively for the one-time final test")
    if args.final_test != bool(args.final_test_execution_id):
        raise RuntimeError(
            "--final-test and --final-test-execution-id must be provided together"
        )

    means_path = METADATA_DIR / "track_nonzero_means_v2.tsv"
    means_rows = read_tsv(means_path)
    default_mean_column = (
        "development_train_nonzero_mean"
        if args.final_test
        else f"fold_{args.fold}_train_nonzero_mean"
    )
    mean_column = args.mean_column or default_mean_column
    if args.final_test and mean_column != "development_train_nonzero_mean":
        raise ValueError("final test requires development_train_nonzero_mean")
    if mean_column not in means_rows[0]:
        raise ValueError(f"Unknown mean column: {mean_column}")
    fold_means = torch.tensor(
        [float(row[mean_column]) for row in means_rows],
        dtype=torch.float32,
        device=device,
    )
    intervals_path = (
        REPO_ROOT / "alphagenome_custom/intervals/v2/test_locked.tsv"
        if args.final_test
        else REPO_ROOT / f"alphagenome_custom/intervals/v2/fold_{args.fold}/valid.tsv"
    )
    model = train_v2_model.build_model(
        args.model,
        fold_means.detach().cpu(),
        device,
        args.hidden_channels,
    )
    checkpoint_path = REPO_ROOT / args.checkpoint
    checkpoint = load_checkpoint(
        model,
        checkpoint_path,
        expected_model=args.model,
        expected_loss=args.training_loss,
        expected_fold=args.fold,
        expected_sequence_length=args.sequence_length,
        expected_hidden_channels=args.hidden_channels,
        expected_mean_column=mean_column,
    )
    checkpoint_sha = sha256(checkpoint_path)
    dataset = V2BigWigDataset(
        intervals_path,
        max_io_workers=16,
        final_test_checkpoint_sha256=checkpoint_sha if args.final_test else None,
        final_test_execution_id=(
            args.final_test_execution_id if args.final_test else None
        ),
    )
    subwindows, eligible_bases = core_subwindows(
        dataset.intervals, args.sequence_length
    )
    n_tracks = len(means_rows)
    group_ids = [row["group_id"] for row in means_rows]
    chromosomes = {row["chromosome"] for row in dataset.intervals}
    genes, genes_by_chromosome = full_metrics.load_gene_exons(
        dataset.gtf_path, chromosomes
    )
    gene_exons = full_metrics.GeneExonAccumulator(
        genes, genes_by_chromosome, n_tracks=n_tracks
    )
    exon_intervals: dict[str, tuple[tuple[int, int], ...]] = {}
    for chromosome in chromosomes:
        exon_intervals[chromosome] = full_metrics.merge_intervals(
            interval
            for gene_index in genes_by_chromosome.get(chromosome, [])
            for interval in genes[gene_index].intervals
        )
    raw_1bp_stats = full_metrics.PerTrackStats(n_tracks)
    raw_128bp_stats = full_metrics.PerTrackStats(n_tracks)
    log1p_1bp_stats = full_metrics.PerTrackStats(n_tracks)
    log1p_128bp_stats = full_metrics.PerTrackStats(n_tracks)
    gene_body_stats = full_metrics.PerTrackStats(n_tracks)
    exon_stats = full_metrics.PerTrackStats(n_tracks)
    gradient_1bp_stats = full_metrics.PerTrackStats(n_tracks)
    gradient_128bp_stats = full_metrics.PerTrackStats(n_tracks)
    prediction_128_chunks: list[np.ndarray] = []
    target_128_chunks: list[np.ndarray] = []
    if args.final_test:
        consume_final_test(
            checkpoint_sha, str(args.final_test_execution_id), int(visible)
        )
    model.eval()
    metric_sums: dict[str, float] = {}
    count_128 = 0
    sum_x = np.zeros(n_tracks, dtype=np.float64)
    sum_y = np.zeros(n_tracks, dtype=np.float64)
    sum_xx = np.zeros(n_tracks, dtype=np.float64)
    sum_yy = np.zeros(n_tracks, dtype=np.float64)
    sum_xy = np.zeros(n_tracks, dtype=np.float64)
    with torch.no_grad():
        for ordinal, (index, offset) in enumerate(subwindows, start=1):
            item = dataset.get_subwindow_item(
                index,
                shift_bp=0,
                crop_offset_bp=offset,
                crop_length_bp=args.sequence_length,
            )
            dna = item["dna_sequence"].unsqueeze(0).to(device)
            targets = {
                1: item["target_1bp"].unsqueeze(0).to(device),
                128: item["target_128bp"].unsqueeze(0).to(device),
            }
            track_mask = item["track_mask"].unsqueeze(0).to(device)
            track_strand = item["track_strand"].unsqueeze(0).to(device)
            gene_mask = item["gene_mask"].unsqueeze(0).to(device)
            with torch.autocast(device_type="cuda", dtype=torch.bfloat16):
                predictions = model(dna)
                paper_loss, paper_metrics = components.dual_resolution_paper_loss(
                    predictions,
                    targets,
                    track_means=fold_means,
                    track_mask=track_mask,
                    track_strand=track_strand,
                    gene_mask=gene_mask,
                )
                log_loss, log_metrics = components.log1p_mse_loss(
                    predictions,
                    targets,
                    track_means=fold_means,
                    track_mask=track_mask,
                )
            current = {
                "paper_loss": paper_loss,
                "log1p_mse": log_loss,
                **{f"paper_{key}": value for key, value in paper_metrics.items()},
                **{f"log1p_{key}": value for key, value in log_metrics.items()},
            }
            for key, value in current.items():
                numeric = float(value.detach())
                if not math.isfinite(numeric):
                    raise RuntimeError(f"Non-finite validation metric {key}")
                metric_sums[key] = metric_sums.get(key, 0.0) + numeric

            prediction_128 = components.unscale_predictions_experimental_space(
                predictions[128].float(), fold_means, 128
            )
            prediction_1 = components.unscale_predictions_experimental_space(
                predictions[1].float(), fold_means, 1
            )
            raw_prediction_1 = prediction_1.clamp_min(0)[0].float().cpu().numpy()
            raw_target_1 = targets[1].clamp_min(0)[0].float().cpu().numpy()
            raw_prediction_128 = prediction_128.clamp_min(0)[0].float().cpu().numpy()
            raw_target_128 = targets[128].clamp_min(0)[0].float().cpu().numpy()
            log_prediction_1 = np.log1p(raw_prediction_1)
            log_target_1 = np.log1p(raw_target_1)
            x = np.log1p(raw_prediction_128)
            y = np.log1p(raw_target_128)
            raw_1bp_stats.update(raw_prediction_1, raw_target_1)
            raw_128bp_stats.update(raw_prediction_128, raw_target_128)
            log1p_1bp_stats.update(log_prediction_1, log_target_1)
            log1p_128bp_stats.update(x, y)
            gradient_1bp_stats.update(
                np.diff(log_prediction_1, axis=1), np.diff(log_target_1, axis=1)
            )
            gradient_128bp_stats.update(np.diff(x, axis=1), np.diff(y, axis=1))
            chromosome = str(item["interval_chromosome"])
            interval_start = int(item["interval_start"])
            interval_end = int(item["interval_end"])
            gene_body_mask = item["gene_mask"].any(dim=0).cpu().numpy()
            exon_mask = full_metrics.interval_union_mask(
                exon_intervals.get(chromosome, ()), interval_start, interval_end
            )
            gene_body_stats.update(log_prediction_1, log_target_1, gene_body_mask)
            exon_stats.update(log_prediction_1, log_target_1, exon_mask)
            gene_exons.update(
                chromosome,
                interval_start,
                raw_prediction_1,
                raw_target_1,
            )
            prediction_128_chunks.append(raw_prediction_128)
            target_128_chunks.append(raw_target_128)
            count_128 += x.shape[1]
            sum_x += x.sum(axis=1)
            sum_y += y.sum(axis=1)
            sum_xx += np.square(x).sum(axis=1)
            sum_yy += np.square(y).sum(axis=1)
            sum_xy += (x * y).sum(axis=1)
            print(
                f"validation_subwindow\t{ordinal}/{len(subwindows)}\t"
                f"{item['interval_chromosome']}:{int(item['interval_start'])}-{int(item['interval_end'])}",
                flush=True,
            )
    dataset.close()
    means = {key: value / len(subwindows) for key, value in metric_sums.items()}
    per_track_pearson = pearson_from_sums(
        count_128, sum_x, sum_y, sum_xx, sum_yy, sum_xy
    )
    finite_pearson = per_track_pearson[np.isfinite(per_track_pearson)]
    if not finite_pearson.size:
        raise RuntimeError("No finite per-track validation Pearson values")
    gene_exon_metrics = gene_exons.metrics(log1p=True)
    distribution = full_metrics.distribution_metrics(
        np.concatenate(prediction_128_chunks, axis=1),
        np.concatenate(target_128_chunks, axis=1),
        seed=int(checkpoint["seed"]),
    )
    full_metric_record = {
        "metric_contract": "v2_full_metrics_1",
        "primary": {
            "mean_per_track_gene_exon_coverage_pearson_log1p": gene_exon_metrics[
                "mean_per_track_pearson"
            ],
            "mean_per_track_pearson_128bp_log1p": float(finite_pearson.mean()),
        },
        "raw_1bp": stats_record(raw_1bp_stats, group_ids),
        "raw_128bp_sum": stats_record(raw_128bp_stats, group_ids),
        "log1p_1bp": stats_record(log1p_1bp_stats, group_ids),
        "log1p_128bp_sum": stats_record(log1p_128bp_stats, group_ids),
        "gene_body_log1p_1bp": stats_record(gene_body_stats, group_ids),
        "exon_log1p_1bp": stats_record(exon_stats, group_ids),
        "local_gradient_log1p_1bp": stats_record(gradient_1bp_stats, group_ids),
        "local_gradient_log1p_128bp": stats_record(gradient_128bp_stats, group_ids),
        "gene_exon_coverage": {
            key: value
            for key, value in gene_exon_metrics.items()
            if not isinstance(value, np.ndarray)
        }
        | {
            "per_track_pearson": full_metrics.keyed(
                gene_exon_metrics["per_track_pearson"], group_ids
            ),
            "per_track_mse": full_metrics.keyed(
                gene_exon_metrics["per_track_mse"], group_ids
            ),
        },
        "distribution_128bp": {
            "spearman_sampling": distribution["spearman_sampling"],
            "spearman_points_per_track": distribution["spearman_points_per_track"],
            "mean_per_track_spearman": full_metrics.finite_mean(
                distribution["per_track_spearman_128bp"]
            ),
            "mean_per_track_top1_mse": full_metrics.finite_mean(
                distribution["per_track_top1_mse_128bp"]
            ),
            "mean_per_track_top1_calibration_ratio": full_metrics.finite_mean(
                distribution["per_track_top1_calibration_ratio_128bp"]
            ),
            "per_track_spearman": full_metrics.keyed(
                distribution["per_track_spearman_128bp"], group_ids
            ),
            "per_track_top1_mse": full_metrics.keyed(
                distribution["per_track_top1_mse_128bp"], group_ids
            ),
            "per_track_top1_calibration_ratio": full_metrics.keyed(
                distribution["per_track_top1_calibration_ratio_128bp"], group_ids
            ),
            "per_track_top1_target_threshold": full_metrics.keyed(
                distribution["per_track_top1_target_threshold_128bp"], group_ids
            ),
        },
    }
    record = {
        "schema_version": 2,
        "phase": "P6C" if args.final_test else "P6B",
        "model": args.model,
        "training_loss": args.training_loss,
        "fold": args.fold,
        "mean_column": mean_column,
        "seed": checkpoint["seed"],
        "checkpoint_path": args.checkpoint,
        "checkpoint_sha256": checkpoint_sha,
        "sequence_length": args.sequence_length,
        "validation_policy": "all_complete_nonoverlapping_core_only_131072bp_subwindows",
        "validation_subwindows": len(subwindows),
        "validation_bases": len(subwindows) * args.sequence_length,
        "eligible_aligned_core_bases": eligible_bases,
        "validation_core_coverage_fraction": (
            len(subwindows) * args.sequence_length / eligible_bases
        ),
        "mean_metrics": means,
        "mean_per_track_pearson_128bp": float(finite_pearson.mean()),
        "finite_per_track_pearson_128bp": int(finite_pearson.size),
        "min_per_track_pearson_128bp": float(finite_pearson.min()),
        "max_per_track_pearson_128bp": float(finite_pearson.max()),
        "per_track_pearson_128bp": {
            row["group_id"]: float(value) if math.isfinite(value) else None
            for row, value in zip(means_rows, per_track_pearson, strict=True)
        },
        "full_metrics": full_metric_record,
        "intervals_sha256": sha256(intervals_path),
        "means_sha256": sha256(means_path),
        "physical_cuda_visible_devices": visible,
        "cuda_device_name": torch.cuda.get_device_name(0),
        "evaluated_chromosomes": sorted(chromosomes),
        "locked_test_block_signal_reads": args.final_test,
        "chromosome_x_train_valid_blocks_read": (
            "X" in chromosomes and not args.final_test
        ),
        "final_test_execution_id": (
            args.final_test_execution_id if args.final_test else None
        ),
    }
    output_path = REPO_ROOT / args.output
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
    if args.final_test:
        finalize_final_test(
            checkpoint_sha, str(args.final_test_execution_id), output_path
        )
    print(json.dumps(record, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
