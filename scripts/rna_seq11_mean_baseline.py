#!/usr/bin/env python3
"""Evaluate a per-track mean RNA-seq baseline on NPZ datasets."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--train-dataset-dir",
        default="alphagenome_custom/datasets/rna_seq_npz_train",
        help="Training NPZ dataset directory used to estimate per-track means.",
    )
    parser.add_argument(
        "--eval-dataset-dir",
        default="alphagenome_custom/datasets/rna_seq_npz_valid",
        help="Evaluation NPZ dataset directory.",
    )
    parser.add_argument(
        "--target-transform",
        choices=["log1p", "none"],
        default="log1p",
        help="Target transform to match model training/evaluation.",
    )
    parser.add_argument("--max-train-examples", type=int, default=None)
    parser.add_argument("--max-eval-examples", type=int, default=None)
    parser.add_argument(
        "--metrics-output",
        default=None,
        help="Optional TSV path for overall and per-track baseline metrics.",
    )
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_metadata(dataset_dir: Path) -> dict[str, object]:
    with (dataset_dir / "dataset_metadata.json").open() as handle:
        return json.load(handle)


def load_track_labels(metadata: dict[str, object], n_tracks: int) -> list[tuple[str, str]]:
    metadata_path = metadata.get("track_metadata")
    if not metadata_path:
        return [(f"track_{index}", f"track_{index}") for index in range(n_tracks)]

    rows: list[tuple[int, str, str]] = []
    with Path(str(metadata_path)).open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if row.get("output_type") != "RNA_SEQ":
                continue
            rows.append(
                (
                    int(row["track_index"]),
                    row.get("group_id", f"track_{row['track_index']}"),
                    row.get("name", row.get("group_id", f"track_{row['track_index']}")),
                )
            )
    rows.sort(key=lambda item: item[0])
    if len(rows) != n_tracks:
        return [(f"track_{index}", f"track_{index}") for index in range(n_tracks)]
    return [(group_id, name) for _index, group_id, name in rows]


def iter_npz_paths(dataset_dir: Path, max_examples: int | None = None) -> list[Path]:
    rows = read_tsv(dataset_dir / "manifest.tsv")
    if max_examples is not None:
        rows = rows[:max_examples]
    return [dataset_dir / row["path"] for row in rows]


def transform_target(rna_seq: np.ndarray, target_transform: str) -> np.ndarray:
    rna_seq = rna_seq.astype(np.float64, copy=False)
    if target_transform == "log1p":
        return np.log1p(np.clip(rna_seq, a_min=0.0, a_max=None))
    if target_transform == "none":
        return rna_seq
    raise ValueError(f"Unknown target_transform: {target_transform}")


def compute_track_mean(
    dataset_dir: Path,
    *,
    target_transform: str,
    max_examples: int | None,
) -> tuple[np.ndarray, int]:
    metadata = read_metadata(dataset_dir)
    n_tracks = int(metadata["n_tracks"])
    total = np.zeros(n_tracks, dtype=np.float64)
    denominator = np.zeros(n_tracks, dtype=np.float64)
    paths = iter_npz_paths(dataset_dir, max_examples=max_examples)

    for path in paths:
        with np.load(path, allow_pickle=False) as data:
            rna_seq = transform_target(data["rna_seq"], target_transform)
            mask = data["rna_seq_mask"].astype(bool, copy=False).reshape(-1)
        if rna_seq.shape[-1] != n_tracks:
            raise ValueError(f"Unexpected rna_seq shape in {path}: {rna_seq.shape}")
        mask_float = mask.astype(np.float64)
        total += (rna_seq * mask_float.reshape(1, -1)).sum(axis=0)
        denominator += mask_float * rna_seq.shape[0]

    return total / np.clip(denominator, a_min=1.0, a_max=None), len(paths)


def evaluate_track_mean(
    dataset_dir: Path,
    *,
    track_mean: np.ndarray,
    target_transform: str,
    max_examples: int | None,
) -> tuple[float, np.ndarray, int]:
    metadata = read_metadata(dataset_dir)
    n_tracks = int(metadata["n_tracks"])
    if len(track_mean) != n_tracks:
        raise ValueError(f"track_mean has {len(track_mean)} tracks; expected {n_tracks}")

    sse = np.zeros(n_tracks, dtype=np.float64)
    denominator = np.zeros(n_tracks, dtype=np.float64)
    paths = iter_npz_paths(dataset_dir, max_examples=max_examples)

    for path in paths:
        with np.load(path, allow_pickle=False) as data:
            rna_seq = transform_target(data["rna_seq"], target_transform)
            mask = data["rna_seq_mask"].astype(bool, copy=False).reshape(-1)
        mask_float = mask.astype(np.float64)
        error = rna_seq - track_mean.reshape(1, -1)
        sse += (np.square(error) * mask_float.reshape(1, -1)).sum(axis=0)
        denominator += mask_float * rna_seq.shape[0]

    per_track_loss = sse / np.clip(denominator, a_min=1.0, a_max=None)
    overall_loss = float(sse.sum() / max(float(denominator.sum()), 1.0))
    return overall_loss, per_track_loss, len(paths)


def write_metrics(
    path: Path,
    *,
    overall_loss: float,
    per_track_loss: np.ndarray,
    track_mean: np.ndarray,
    track_labels: list[tuple[str, str]],
    train_examples: int,
    eval_examples: int,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "metric_scope",
                "track_index",
                "track_id",
                "track_name",
                "loss",
                "train_mean",
                "train_examples",
                "eval_examples",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerow(
            {
                "metric_scope": "overall",
                "track_index": "",
                "track_id": "",
                "track_name": "",
                "loss": f"{overall_loss:.8g}",
                "train_mean": "",
                "train_examples": train_examples,
                "eval_examples": eval_examples,
            }
        )
        for index, (track_id, track_name) in enumerate(track_labels):
            writer.writerow(
                {
                    "metric_scope": "track",
                    "track_index": index,
                    "track_id": track_id,
                    "track_name": track_name,
                    "loss": f"{per_track_loss[index]:.8g}",
                    "train_mean": f"{track_mean[index]:.8g}",
                    "train_examples": train_examples,
                    "eval_examples": eval_examples,
                }
            )


def main() -> None:
    args = parse_args()
    train_dataset_dir = Path(args.train_dataset_dir)
    eval_dataset_dir = Path(args.eval_dataset_dir)

    track_mean, train_examples = compute_track_mean(
        train_dataset_dir,
        target_transform=args.target_transform,
        max_examples=args.max_train_examples,
    )
    overall_loss, per_track_loss, eval_examples = evaluate_track_mean(
        eval_dataset_dir,
        track_mean=track_mean,
        target_transform=args.target_transform,
        max_examples=args.max_eval_examples,
    )
    eval_metadata = read_metadata(eval_dataset_dir)
    track_labels = load_track_labels(eval_metadata, len(track_mean))

    print(f"train_dataset_dir\t{train_dataset_dir}")
    print(f"eval_dataset_dir\t{eval_dataset_dir}")
    print(f"target_transform\t{args.target_transform}")
    print(f"train_examples\t{train_examples}")
    print(f"eval_examples\t{eval_examples}")
    print(f"n_tracks\t{len(track_mean)}")
    print(f"valid_loss\t{overall_loss:.8g}")
    print("per_track_loss_header\ttrack_index\ttrack_id\ttrack_name\tloss\ttrain_mean")
    for index, (track_id, track_name) in enumerate(track_labels):
        print(
            "per_track_loss\t"
            f"{index}\t{track_id}\t{track_name}\t"
            f"{per_track_loss[index]:.8g}\t{track_mean[index]:.8g}"
        )

    if args.metrics_output is not None:
        write_metrics(
            Path(args.metrics_output),
            overall_loss=overall_loss,
            per_track_loss=per_track_loss,
            track_mean=track_mean,
            track_labels=track_labels,
            train_examples=train_examples,
            eval_examples=eval_examples,
        )
        print(f"metrics_output\t{args.metrics_output}")


if __name__ == "__main__":
    main()
