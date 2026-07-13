#!/usr/bin/env python3
"""Build locked five-fold chromosome splits and non-overlapping eval cores."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[1]
OUTPUT_DIR = REPO_ROOT / "alphagenome_custom/intervals/v2"
METADATA_DIR = REPO_ROOT / "alphagenome_custom/metadata/v2"
FAI_PATH = REPO_ROOT / "alphagenome_custom/reference/genome.fa.fai"
GROUP_MANIFEST = METADATA_DIR / "rna_seq_groups_v2_final.tsv"
GROUP_OUTPUTS = METADATA_DIR / "p3_group_outputs.tsv"
WINDOW_SIZE = 2**20
STRIDE = 2**19
CV_CHROMOSOMES = ("I", "II", "III", "IV", "V")
TEST_CHROMOSOME = "X"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(8 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def read_chromosomes(path: Path = FAI_PATH) -> dict[str, int]:
    chromosomes = {}
    with path.open() as handle:
        for line in handle:
            name, length, *_ = line.rstrip("\n").split("\t")
            chromosomes[name] = int(length)
    return chromosomes


def make_windows(length: int) -> list[tuple[int, int]]:
    if length < WINDOW_SIZE:
        raise ValueError(f"Chromosome length {length} is shorter than {WINDOW_SIZE}")
    starts = list(range(0, length - WINDOW_SIZE + 1, STRIDE))
    final_start = length - WINDOW_SIZE
    if starts[-1] != final_start:
        starts.append(final_start)
    return [(start, start + WINDOW_SIZE) for start in sorted(set(starts))]


def eval_cores(
    windows: list[tuple[int, int]], chromosome_length: int
) -> list[tuple[int, int]]:
    centers = [(start + end) // 2 for start, end in windows]
    boundaries = [0]
    boundaries.extend((first + second) // 2 for first, second in zip(centers, centers[1:]))
    boundaries.append(chromosome_length)
    cores = list(zip(boundaries[:-1], boundaries[1:], strict=True))
    for (window_start, window_end), (core_start, core_end) in zip(
        windows, cores, strict=True
    ):
        if not (window_start <= core_start < core_end <= window_end):
            raise RuntimeError(
                f"Evaluation core {core_start}-{core_end} escapes "
                f"window {window_start}-{window_end}"
            )
    return cores


def write_tsv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def interval_rows(
    chromosome: str,
    length: int,
    fold: int | str,
    role: str,
) -> list[dict[str, Any]]:
    windows = make_windows(length)
    cores = eval_cores(windows, length) if role in {"valid", "test_locked"} else windows
    return [
        {
            "chromosome": chromosome,
            "start": start,
            "end": end,
            "core_start": core_start,
            "core_end": core_end,
            "fold": fold,
            "role": role,
        }
        for (start, end), (core_start, core_end) in zip(windows, cores, strict=True)
    ]


def main() -> None:
    if not GROUP_MANIFEST.is_file() or not GROUP_OUTPUTS.is_file():
        raise RuntimeError("P3B final group artifacts are required before locking P4")
    chromosomes = read_chromosomes()
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    files = []
    fold_summaries = []
    for fold, valid_chromosome in enumerate(CV_CHROMOSOMES, start=1):
        train_rows = []
        for chromosome in CV_CHROMOSOMES:
            if chromosome != valid_chromosome:
                train_rows.extend(
                    interval_rows(chromosome, chromosomes[chromosome], fold, "train")
                )
        valid_rows = interval_rows(
            valid_chromosome, chromosomes[valid_chromosome], fold, "valid"
        )
        fold_dir = OUTPUT_DIR / f"fold_{fold}"
        train_path = fold_dir / "train.tsv"
        valid_path = fold_dir / "valid.tsv"
        write_tsv(train_path, train_rows)
        write_tsv(valid_path, valid_rows)
        files.extend((train_path, valid_path))
        fold_summaries.append(
            {
                "fold": fold,
                "valid_chromosome": valid_chromosome,
                "train_chromosomes": [
                    chromosome for chromosome in CV_CHROMOSOMES
                    if chromosome != valid_chromosome
                ],
                "train_windows": len(train_rows),
                "valid_windows": len(valid_rows),
            }
        )

    development_rows = []
    for chromosome in CV_CHROMOSOMES:
        development_rows.extend(
            interval_rows(
                chromosome,
                chromosomes[chromosome],
                "development",
                "train",
            )
        )
    development_path = OUTPUT_DIR / "development_train.tsv"
    write_tsv(development_path, development_rows)
    files.append(development_path)

    test_rows = interval_rows(
        TEST_CHROMOSOME, chromosomes[TEST_CHROMOSOME], "final", "test_locked"
    )
    test_path = OUTPUT_DIR / "test_locked.tsv"
    write_tsv(test_path, test_rows)
    files.append(test_path)
    registry = {
        "schema_version": 1,
        "window_size": WINDOW_SIZE,
        "stride": STRIDE,
        "cv_chromosomes": list(CV_CHROMOSOMES),
        "test_chromosome": TEST_CHROMOSOME,
        "test_status": "embargoed_prior_exposure_disclosed",
        "evaluation_core_policy": "nearest_window_center_partition_no_overlap",
        "folds": fold_summaries,
        "development_train_chromosomes": list(CV_CHROMOSOMES),
        "development_train_windows": len(development_rows),
        "test_windows": len(test_rows),
        "files": {
            str(path.relative_to(REPO_ROOT)): sha256(path) for path in files
        },
        "reference_fai_sha256": sha256(FAI_PATH),
        "group_manifest_sha256": sha256(GROUP_MANIFEST),
        "group_outputs_sha256": sha256(GROUP_OUTPUTS),
    }
    registry_path = METADATA_DIR / "split_registry_v2.json"
    registry_path.write_text(json.dumps(registry, indent=2, sort_keys=True) + "\n")
    print(json.dumps(registry, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
