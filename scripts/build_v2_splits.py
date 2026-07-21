#!/usr/bin/env python3
"""Build locked six-chromosome block splits and non-overlapping eval cores."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
import random
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[1]
OUTPUT_DIR = REPO_ROOT / "alphagenome_custom/intervals/v2"
METADATA_DIR = REPO_ROOT / "alphagenome_custom/metadata/v2"
FAI_PATH = REPO_ROOT / "alphagenome_custom/reference/genome.fa.fai"
GROUP_MANIFEST = METADATA_DIR / "rna_seq_groups_v2_final.tsv"
GROUP_OUTPUTS = METADATA_DIR / "p3_group_outputs.tsv"
SPEC_PATH = METADATA_DIR / "six_chromosome_split_spec.json"
BLOCKS_PATH = OUTPUT_DIR / "blocks.tsv"
CHROMOSOMES = ("I", "II", "III", "IV", "V", "X")
CV_FOLDS = tuple(range(1, 6))
WINDOW_SIZE = 2**20
STRIDE = 2**19
MAX_SHIFT_BP = 1024
EXCLUSION_BUFFER_BP = WINDOW_SIZE + 2 * MAX_SHIFT_BP
EVALUATION_SUBWINDOW_BP = 2**17


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


def make_windows(
    length: int,
    *,
    interval_start: int = 0,
    interval_end: int | None = None,
    shift_margin: int = 0,
) -> list[tuple[int, int]]:
    interval_end = length if interval_end is None else interval_end
    usable_start = interval_start + shift_margin
    usable_end = interval_end - shift_margin
    if usable_start < 0 or usable_end > length or usable_end - usable_start < WINDOW_SIZE:
        raise ValueError(
            f"Interval {interval_start}-{interval_end} on length {length} cannot fit "
            f"a {WINDOW_SIZE}-bp window with {shift_margin}-bp shift margins"
        )
    starts = list(range(usable_start, usable_end - WINDOW_SIZE + 1, STRIDE))
    final_start = usable_end - WINDOW_SIZE
    if starts[-1] != final_start:
        starts.append(final_start)
    return [(start, start + WINDOW_SIZE) for start in sorted(set(starts))]


def eval_cores(
    windows: list[tuple[int, int]],
    chromosome_length: int,
    *,
    interval_start: int = 0,
    interval_end: int | None = None,
) -> list[tuple[int, int]]:
    interval_end = chromosome_length if interval_end is None else interval_end
    centers = [(start + end) // 2 for start, end in windows]
    boundaries = [interval_start]
    boundaries.extend((first + second) // 2 for first, second in zip(centers, centers[1:]))
    boundaries.append(interval_end)
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


def split_blocks(
    chromosomes: dict[str, int], spec: dict[str, Any]
) -> list[dict[str, Any]]:
    labels = list(spec["block_labels"])
    buffer_bp = int(spec["exclusion_buffer_bp"])
    minimum_bp = int(spec["minimum_effective_block_bp"])
    seed = int(spec["seed"])
    rows = []
    for chromosome in CHROMOSOMES:
        length = chromosomes[chromosome]
        effective_bases = length - (len(labels) - 1) * buffer_bp
        aligned_minimum_bp = (
            (minimum_bp + EVALUATION_SUBWINDOW_BP - 1)
            // EVALUATION_SUBWINDOW_BP
            * EVALUATION_SUBWINDOW_BP
        )
        if effective_bases < len(labels) * aligned_minimum_bp:
            raise RuntimeError(
                f"{chromosome} cannot fit {len(labels)} aligned blocks, buffers, "
                "and minimum size"
            )
        block_size = (
            effective_bases // len(labels) // EVALUATION_SUBWINDOW_BP
        ) * EVALUATION_SUBWINDOW_BP
        sizes = [block_size] * len(labels)
        outer_remainder = effective_bases - block_size * len(labels)
        outer_before = outer_remainder // 2
        outer_after = outer_remainder - outer_before
        stable_seed = int.from_bytes(
            hashlib.sha256(f"{seed}:{chromosome}".encode()).digest()[:8], "big"
        )
        assigned = labels.copy()
        random.Random(stable_seed).shuffle(assigned)
        cursor = outer_before
        for physical_order, (label, size) in enumerate(zip(assigned, sizes, strict=True)):
            start = cursor
            end = start + size
            rows.append(
                {
                    "chromosome": chromosome,
                    "block_id": f"{chromosome}_block_{physical_order + 1}",
                    "physical_order": physical_order + 1,
                    "assignment": label,
                    "block_start": start,
                    "block_end": end,
                    "effective_bases": size,
                    "buffer_before_bp": 0 if physical_order == 0 else buffer_bp,
                    "buffer_after_bp": 0
                    if physical_order == len(labels) - 1
                    else buffer_bp,
                    "outer_margin_before_bp": outer_before
                    if physical_order == 0
                    else 0,
                    "outer_margin_after_bp": outer_after
                    if physical_order == len(labels) - 1
                    else 0,
                }
            )
            cursor = end + (buffer_bp if physical_order < len(labels) - 1 else 0)
        if cursor + outer_after != length:
            raise RuntimeError(f"Block construction did not consume {chromosome} length")
    return rows


def interval_rows(
    chromosome: str,
    length: int,
    fold: int | str,
    role: str,
    block: dict[str, Any],
    revision_id: str,
) -> list[dict[str, Any]]:
    block_start = int(block["block_start"])
    block_end = int(block["block_end"])
    shift_margin = MAX_SHIFT_BP if role == "train" else 0
    windows = make_windows(
        length,
        interval_start=block_start,
        interval_end=block_end,
        shift_margin=shift_margin,
    )
    cores = (
        eval_cores(
            windows,
            length,
            interval_start=block_start,
            interval_end=block_end,
        )
        if role in {"valid", "test_locked"}
        else windows
    )
    return [
        {
            "chromosome": chromosome,
            "start": start,
            "end": end,
            "core_start": core_start,
            "core_end": core_end,
            "fold": fold,
            "role": role,
            "block_id": block["block_id"],
            "block_start": block_start,
            "block_end": block_end,
            "split_revision": revision_id,
        }
        for (start, end), (core_start, core_end) in zip(windows, cores, strict=True)
    ]


def main() -> None:
    if not GROUP_MANIFEST.is_file() or not GROUP_OUTPUTS.is_file():
        raise RuntimeError("P3B final group artifacts are required before locking P4")
    if not SPEC_PATH.is_file():
        raise RuntimeError(f"Missing preregistered split specification: {SPEC_PATH}")
    spec = json.loads(SPEC_PATH.read_text())
    expected = {
        "revision_id": "six_chromosome_blocks_v1",
        "chromosomes": list(CHROMOSOMES),
        "window_size_bp": WINDOW_SIZE,
        "stride_bp": STRIDE,
        "max_shift_bp": MAX_SHIFT_BP,
        "exclusion_buffer_bp": EXCLUSION_BUFFER_BP,
        "evaluation_subwindow_bp": EVALUATION_SUBWINDOW_BP,
        "effective_block_alignment_bp": EVALUATION_SUBWINDOW_BP,
    }
    mismatches = {
        key: (spec.get(key), value)
        for key, value in expected.items()
        if spec.get(key) != value
    }
    if mismatches:
        raise RuntimeError(f"Split specification mismatch: {mismatches}")
    chromosomes = read_chromosomes()
    if tuple(name for name in chromosomes if name != "MtDNA") != CHROMOSOMES:
        raise RuntimeError("WBcel235 nuclear chromosome names/order do not match I-V,X")
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    blocks = split_blocks(chromosomes, spec)
    write_tsv(BLOCKS_PATH, blocks)
    blocks_by_assignment = {
        (row["chromosome"], row["assignment"]): row for row in blocks
    }
    files = []
    fold_summaries = []
    for fold in CV_FOLDS:
        train_rows = []
        valid_rows = []
        for chromosome in CHROMOSOMES:
            for train_fold in CV_FOLDS:
                if train_fold == fold:
                    continue
                block = blocks_by_assignment[(chromosome, f"cv_fold_{train_fold}")]
                train_rows.extend(
                    interval_rows(
                        chromosome,
                        chromosomes[chromosome],
                        fold,
                        "train",
                        block,
                        spec["revision_id"],
                    )
                )
            valid_block = blocks_by_assignment[(chromosome, f"cv_fold_{fold}")]
            valid_rows.extend(
                interval_rows(
                    chromosome,
                    chromosomes[chromosome],
                    fold,
                    "valid",
                    valid_block,
                    spec["revision_id"],
                )
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
                "train_chromosomes": list(CHROMOSOMES),
                "valid_chromosomes": list(CHROMOSOMES),
                "train_block_assignments": [
                    f"cv_fold_{value}" for value in CV_FOLDS if value != fold
                ],
                "valid_block_assignment": f"cv_fold_{fold}",
                "train_windows": len(train_rows),
                "valid_windows": len(valid_rows),
                "train_effective_bases": sum(
                    int(row["effective_bases"])
                    for row in blocks
                    if row["assignment"] in {
                        f"cv_fold_{value}" for value in CV_FOLDS if value != fold
                    }
                ),
                "valid_effective_bases": sum(
                    int(row["effective_bases"])
                    for row in blocks
                    if row["assignment"] == f"cv_fold_{fold}"
                ),
            }
        )

    development_rows = []
    for chromosome in CHROMOSOMES:
        for fold in CV_FOLDS:
            development_rows.extend(
                interval_rows(
                    chromosome,
                    chromosomes[chromosome],
                    "development",
                    "train",
                    blocks_by_assignment[(chromosome, f"cv_fold_{fold}")],
                    spec["revision_id"],
                )
            )
    development_path = OUTPUT_DIR / "development_train.tsv"
    write_tsv(development_path, development_rows)
    files.append(development_path)

    test_rows = []
    for chromosome in CHROMOSOMES:
        test_rows.extend(
            interval_rows(
                chromosome,
                chromosomes[chromosome],
                "final",
                "test_locked",
                blocks_by_assignment[(chromosome, "test_locked")],
                spec["revision_id"],
            )
        )
    test_path = OUTPUT_DIR / "test_locked.tsv"
    write_tsv(test_path, test_rows)
    files.append(test_path)
    registry = {
        "schema_version": 2,
        "revision_id": spec["revision_id"],
        "split_spec_path": str(SPEC_PATH.relative_to(REPO_ROOT)),
        "split_spec_sha256": sha256(SPEC_PATH),
        "window_size": WINDOW_SIZE,
        "stride": STRIDE,
        "max_shift_bp": MAX_SHIFT_BP,
        "exclusion_buffer_bp": EXCLUSION_BUFFER_BP,
        "evaluation_subwindow_bp": EVALUATION_SUBWINDOW_BP,
        "cv_chromosomes": list(CHROMOSOMES),
        "test_chromosomes": list(CHROMOSOMES),
        "test_status": "embargoed_six_chromosome_locked_blocks",
        "final_test_scope": spec["final_test_scope"],
        "evaluation_core_policy": spec["evaluation_core_policy"],
        "folds": fold_summaries,
        "development_train_chromosomes": list(CHROMOSOMES),
        "development_train_windows": len(development_rows),
        "test_windows": len(test_rows),
        "test_effective_bases": sum(
            int(row["effective_bases"])
            for row in blocks
            if row["assignment"] == "test_locked"
        ),
        "blocks_path": str(BLOCKS_PATH.relative_to(REPO_ROOT)),
        "blocks_sha256": sha256(BLOCKS_PATH),
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
