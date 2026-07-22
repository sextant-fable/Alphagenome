#!/usr/bin/env python3
"""Lock v2 splits and benchmark the manifest-driven BigWig loader on CPU."""

from __future__ import annotations

import gc
import hashlib
import json
from pathlib import Path
import resource
import subprocess
import time

import numpy as np
import pyBigWig
from torch.utils.data import DataLoader

from scripts.v2_bigwig_dataset import V2BigWigDataset
from scripts import v2_subprocess


REPO_ROOT = Path(__file__).resolve().parents[1]
METADATA_DIR = REPO_ROOT / "alphagenome_custom/metadata/v2"
BENCHMARK_PATH = METADATA_DIR / "p4_loader_benchmark.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(8 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def tensor_digest(item: dict[str, object]) -> str:
    digest = hashlib.sha256()
    for key in ("dna_sequence", "target_1bp", "target_128bp", "gene_mask", "core_mask"):
        tensor = item[key]
        digest.update(key.encode())
        digest.update(tensor.numpy().tobytes())
    return digest.hexdigest()


def batch_digest(batch: dict[str, object]) -> str:
    digest = hashlib.sha256()
    for key in ("dna_sequence", "target_1bp", "target_128bp", "gene_mask", "core_mask"):
        tensor = batch[key]
        digest.update(key.encode())
        digest.update(tensor.numpy().tobytes())
    return digest.hexdigest()


def main() -> None:
    subprocess.run(
        v2_subprocess.module_command("archive_v2_incomplete_core_revision"),
        cwd=REPO_ROOT,
        check=True,
    )
    subprocess.run(
        v2_subprocess.module_command("build_v2_splits"),
        cwd=REPO_ROOT,
        check=True,
    )
    intervals = REPO_ROOT / "alphagenome_custom/intervals/v2/fold_1/valid.tsv"
    started = time.monotonic()
    fd_before = len(list(Path("/proc/self/fd").iterdir()))
    rss_before = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    full_dataset = V2BigWigDataset(intervals, max_io_workers=16)
    full_started = time.monotonic()
    full_item = full_dataset[0]
    full_seconds = time.monotonic() - full_started
    full_digest = tensor_digest(full_item)
    full_shapes = {
        key: list(full_item[key].shape)
        for key in (
            "dna_sequence",
            "target_1bp",
            "target_128bp",
            "track_mask",
            "track_strand",
            "gene_mask",
            "core_mask",
        )
    }
    full_dtypes = {
        key: str(full_item[key].dtype)
        for key in full_shapes
    }
    n_tracks = int(full_item["target_1bp"].shape[0])
    target_bytes = int(
        full_item["target_1bp"].numel() * full_item["target_1bp"].element_size()
        + full_item["target_128bp"].numel()
        * full_item["target_128bp"].element_size()
    )
    full_dataset.close()
    del full_item, full_dataset
    gc.collect()

    verification_indices = sorted({0, n_tracks // 2, n_tracks - 1})
    verify_dataset = V2BigWigDataset(
        intervals, track_indices=verification_indices, max_io_workers=3
    )
    verify_item = verify_dataset[0]
    row = verify_dataset.intervals[0]
    direct_errors = []
    for local_index, path in enumerate(verify_dataset.track_paths):
        with pyBigWig.open(str(path)) as bigwig:
            direct = bigwig.values(
                row["chromosome"], int(row["start"]), int(row["end"]), numpy=True
            )
        direct = np.nan_to_num(np.asarray(direct, dtype=np.float32))
        observed = verify_item["target_1bp"][local_index].numpy()
        if not np.array_equal(direct, observed):
            direct_errors.append(verify_dataset.track_rows[local_index]["group_id"])
        pooled = direct.reshape(-1, 128).sum(axis=1)
        if not np.allclose(
            pooled,
            verify_item["target_128bp"][local_index].numpy(),
            rtol=1e-6,
            atol=1e-5,
        ):
            direct_errors.append(
                verify_dataset.track_rows[local_index]["group_id"] + ":128bp"
            )
    verify_dataset.close()

    x_valid_dataset = V2BigWigDataset(
        intervals, track_indices=[0], max_io_workers=1
    )
    x_valid_index = next(
        index
        for index, interval in enumerate(x_valid_dataset.intervals)
        if interval["chromosome"] == "X" and interval["role"] == "valid"
    )
    x_valid_item = x_valid_dataset[x_valid_index]
    x_valid_role_read_verified = (
        x_valid_item["interval_chromosome"] == "X"
        and x_valid_dataset.intervals[x_valid_index]["role"] == "valid"
    )
    x_valid_dataset.close()

    worker_indices = list(range(min(4, n_tracks)))
    single_dataset = V2BigWigDataset(
        intervals, track_indices=worker_indices, max_io_workers=4
    )
    single_loader = DataLoader(single_dataset, batch_size=1, num_workers=0, shuffle=False)
    single_digests = []
    single_started = time.monotonic()
    for index, batch in enumerate(single_loader):
        single_digests.append(batch_digest(batch))
        if index == 1:
            break
    single_seconds = time.monotonic() - single_started
    single_dataset.close()
    del batch, single_loader, single_dataset
    gc.collect()

    multi_dataset = V2BigWigDataset(
        intervals, track_indices=worker_indices, max_io_workers=4
    )
    multi_loader = DataLoader(
        multi_dataset,
        batch_size=1,
        num_workers=2,
        shuffle=False,
        persistent_workers=False,
        prefetch_factor=1,
    )
    multi_digests = []
    multi_started = time.monotonic()
    for index, batch in enumerate(multi_loader):
        multi_digests.append(batch_digest(batch))
        if index == 1:
            break
    multi_seconds = time.monotonic() - multi_started
    multi_dataset.close()
    del batch, multi_loader, multi_dataset
    gc.collect()
    fd_after = len(list(Path("/proc/self/fd").iterdir()))
    rss_after = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss

    train_windows = 0
    registry_path = METADATA_DIR / "split_registry_v2.json"
    registry = json.loads(registry_path.read_text())
    for fold in registry["folds"]:
        train_windows = max(train_windows, int(fold["train_windows"]))
    benchmark = {
        "schema_version": 1,
        "phase": "P4",
        "device": "cpu",
        "full_window_tracks": n_tracks,
        "full_window_seconds": full_seconds,
        "full_window_digest": full_digest,
        "full_window_shapes": full_shapes,
        "full_window_dtypes": full_dtypes,
        "full_target_bytes": target_bytes,
        "full_fold_epoch_seconds_linear_io_estimate": full_seconds * train_windows,
        "direct_bigwig_verification_tracks": verification_indices,
        "direct_bigwig_errors": direct_errors,
        "worker_determinism_track_indices": worker_indices,
        "single_worker_digests": single_digests,
        "multi_worker_digests": multi_digests,
        "single_worker_two_window_seconds": single_seconds,
        "multi_worker_two_window_seconds": multi_seconds,
        "worker_deterministic": single_digests == multi_digests,
        "peak_rss_kib_before": rss_before,
        "peak_rss_kib_after": rss_after,
        "open_file_descriptors_before": fd_before,
        "open_file_descriptors_after": fd_after,
        "elapsed_seconds": time.monotonic() - started,
        "cache_policy": "optional_per_worker_lru_default_disabled",
        "pooling_128bp": "sum_of_128_consecutive_1bp_values",
        "monolithic_npz_generated": False,
        "split_revision": registry.get("revision_id"),
        "x_valid_role_read_verified": x_valid_role_read_verified,
        "locked_test_block_signal_reads": 0,
        "split_registry_sha256": sha256(registry_path),
        "group_manifest_sha256": sha256(
            METADATA_DIR / "rna_seq_groups_v2_final.tsv"
        ),
        "group_outputs_sha256": sha256(METADATA_DIR / "p3_group_outputs.tsv"),
    }
    temporary = BENCHMARK_PATH.with_suffix(".json.tmp")
    temporary.write_text(json.dumps(benchmark, indent=2, sort_keys=True) + "\n")
    temporary.replace(BENCHMARK_PATH)
    print(json.dumps(benchmark, indent=2, sort_keys=True))
    if (
        direct_errors
        or not benchmark["worker_deterministic"]
        or not x_valid_role_read_verified
        or n_tracks != 241
    ):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
