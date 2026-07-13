#!/usr/bin/env python3
"""Manifest-driven WBcel235 DNA and grouped RNA-seq BigWig dataset."""

from __future__ import annotations

from collections import OrderedDict, defaultdict
from concurrent.futures import ThreadPoolExecutor
import csv
import json
import os
from pathlib import Path
from typing import Any, Sequence

import numpy as np
import pyBigWig
import torch
from torch.utils.data import Dataset


REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_FASTA = (
    REPO_ROOT
    / "alphagenome_custom/reference/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa"
)
DEFAULT_FAI = REPO_ROOT / "alphagenome_custom/reference/genome.fa.fai"
DEFAULT_GTF = (
    REPO_ROOT / "alphagenome_custom/reference/Caenorhabditis_elegans.WBcel235.115.gtf"
)
DEFAULT_TRACKS = REPO_ROOT / "alphagenome_custom/metadata/v2/p3_group_outputs.tsv"
DEFAULT_GROUPS = REPO_ROOT / "alphagenome_custom/metadata/v2/rna_seq_groups_v2_final.tsv"
STATE_PATH = REPO_ROOT / "alphagenome_custom/metadata/v2/execution_state.json"
FINAL_LOCK_PATH = REPO_ROOT / "alphagenome_custom/metadata/v2/final_test_lock.json"
DNA_TO_INDEX = {ord("A"): 0, ord("C"): 1, ord("G"): 2, ord("T"): 3}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_fai(path: Path) -> dict[str, tuple[int, int, int, int]]:
    result = {}
    with path.open() as handle:
        for line in handle:
            name, length, offset, line_bases, line_width = line.rstrip("\n").split("\t")
            result[name] = (
                int(length), int(offset), int(line_bases), int(line_width)
            )
    return result


def fetch_fasta(
    fasta_path: Path,
    fai: dict[str, tuple[int, int, int, int]],
    chromosome: str,
    start: int,
    end: int,
) -> bytes:
    length, offset, line_bases, line_width = fai[chromosome]
    if start < 0 or end > length or end <= start:
        raise ValueError(f"Out-of-bounds FASTA interval {chromosome}:{start}-{end}")
    sequence_length = end - start
    start_line = start // line_bases
    start_in_line = start % line_bases
    byte_start = offset + start_line * line_width + start_in_line
    line_count = ((start_in_line + sequence_length - 1) // line_bases) + 1
    byte_count = sequence_length + line_count * (line_width - line_bases)
    with fasta_path.open("rb") as handle:
        handle.seek(byte_start)
        raw = handle.read(byte_count)
    sequence = raw.replace(b"\n", b"").replace(b"\r", b"")[:sequence_length].upper()
    if len(sequence) != sequence_length:
        raise RuntimeError(f"Short FASTA read at {chromosome}:{start}-{end}")
    return sequence


def one_hot(sequence: bytes) -> np.ndarray:
    result = np.zeros((4, len(sequence)), dtype=np.float32)
    for index, base in enumerate(sequence):
        channel = DNA_TO_INDEX.get(base)
        if channel is not None:
            result[channel, index] = 1.0
    return result


def parse_gene_intervals(path: Path) -> dict[str, dict[str, list[tuple[int, int]]]]:
    intervals: dict[str, dict[str, list[tuple[int, int]]]] = defaultdict(
        lambda: {"+": [], "-": []}
    )
    with path.open() as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 7 and fields[2] == "gene" and fields[6] in {"+", "-"}:
                intervals[fields[0]][fields[6]].append(
                    (int(fields[3]) - 1, int(fields[4]))
                )
    for chromosome in intervals.values():
        for values in chromosome.values():
            values.sort()
    return dict(intervals)


def interval_mask(
    intervals: list[tuple[int, int]], start: int, end: int
) -> np.ndarray:
    mask = np.zeros(end - start, dtype=bool)
    for interval_start, interval_end in intervals:
        if interval_end <= start:
            continue
        if interval_start >= end:
            break
        local_start = max(interval_start, start) - start
        local_end = min(interval_end, end) - start
        if local_end > local_start:
            mask[local_start:local_end] = True
    return mask


def require_final_test_access(checkpoint_sha256: str | None) -> None:
    state = json.loads(STATE_PATH.read_text())
    gate = state.get("approvals", {}).get("G5_final_test", {})
    if not (
        state.get("current_phase") == "P6C"
        and gate.get("approved") is True
        and gate.get("scope") == "r6c_single_chr_x_test"
        and checkpoint_sha256
        and FINAL_LOCK_PATH.is_file()
    ):
        raise PermissionError(
            "Chromosome X is embargoed until P6C, scoped G5 approval, and a locked checkpoint"
        )
    lock = json.loads(FINAL_LOCK_PATH.read_text())
    if lock.get("checkpoint_sha256") != checkpoint_sha256:
        raise PermissionError("Checkpoint SHA-256 does not match the final-test lock")
    if lock.get("test_consumed") is True:
        raise PermissionError("The one-time chromosome-X test entry has already been consumed")


class V2BigWigDataset(Dataset):
    """Read one DNA window and selected formal v2 tracks on demand."""

    def __init__(
        self,
        intervals_path: str | Path,
        *,
        track_manifest_path: str | Path = DEFAULT_TRACKS,
        group_manifest_path: str | Path = DEFAULT_GROUPS,
        fasta_path: str | Path = DEFAULT_FASTA,
        fai_path: str | Path = DEFAULT_FAI,
        gtf_path: str | Path = DEFAULT_GTF,
        track_indices: Sequence[int] | None = None,
        max_io_workers: int = 16,
        cache_size: int = 0,
        final_test_checkpoint_sha256: str | None = None,
    ) -> None:
        self.intervals_path = Path(intervals_path)
        self.track_manifest_path = Path(track_manifest_path)
        self.group_manifest_path = Path(group_manifest_path)
        self.fasta_path = Path(fasta_path)
        self.fai_path = Path(fai_path)
        self.gtf_path = Path(gtf_path)
        self.intervals = read_tsv(self.intervals_path)
        if not self.intervals:
            raise ValueError("Interval manifest is empty")
        if any(row["chromosome"] == "X" for row in self.intervals):
            require_final_test_access(final_test_checkpoint_sha256)
        outputs = read_tsv(self.track_manifest_path)
        groups = {row["group_id"]: row for row in read_tsv(self.group_manifest_path)}
        if track_indices is None:
            selected = outputs
        else:
            if len(set(track_indices)) != len(track_indices):
                raise ValueError("track_indices must be unique")
            selected = [outputs[index] for index in track_indices]
        self.track_rows = selected
        self.group_rows = [groups[row["group_id"]] for row in selected]
        self.track_paths = [
            path if path.is_absolute() else REPO_ROOT / path
            for path in (Path(row["output_path"]) for row in selected)
        ]
        self.track_strand = np.zeros(len(selected), dtype=np.int8)
        self.track_mask = np.ones((len(selected), 1), dtype=bool)
        self.fai = read_fai(self.fai_path)
        self.genes = parse_gene_intervals(self.gtf_path)
        self.max_io_workers = max(1, int(max_io_workers))
        self.cache_size = max(0, int(cache_size))
        self._cache: OrderedDict[int, dict[str, Any]] = OrderedDict()
        self._owner_pid: int | None = None
        self._handles: list[pyBigWig.pyBigWig] | None = None
        self._executor: ThreadPoolExecutor | None = None

    def __len__(self) -> int:
        return len(self.intervals)

    def _close_handles(self) -> None:
        if self._executor is not None:
            self._executor.shutdown(wait=True)
            self._executor = None
        if self._handles is not None:
            for handle in self._handles:
                handle.close()
            self._handles = None

    def _ensure_handles(self) -> None:
        pid = os.getpid()
        if self._owner_pid != pid:
            self._close_handles()
            self._owner_pid = pid
        if self._handles is None:
            self._handles = [pyBigWig.open(str(path)) for path in self.track_paths]
            if len(self._handles) > 1:
                self._executor = ThreadPoolExecutor(
                    max_workers=min(self.max_io_workers, len(self._handles))
                )

    def _read_signal(self, index: int, chromosome: str, start: int, end: int) -> np.ndarray:
        assert self._handles is not None
        values = self._handles[index].values(chromosome, start, end, numpy=True)
        values = np.asarray(values, dtype=np.float32)
        values = np.nan_to_num(values, nan=0.0, posinf=0.0, neginf=0.0)
        if (values < 0).any():
            raise RuntimeError(f"Negative signal in track {self.track_rows[index]['group_id']}")
        return values

    def _all_signals(self, chromosome: str, start: int, end: int) -> np.ndarray:
        if self._executor is None:
            values = [
                self._read_signal(index, chromosome, start, end)
                for index in range(len(self.track_rows))
            ]
        else:
            futures = [
                self._executor.submit(self._read_signal, index, chromosome, start, end)
                for index in range(len(self.track_rows))
            ]
            values = [future.result() for future in futures]
        return np.stack(values, axis=0)

    def _load_item(self, index: int, shift_bp: int = 0) -> dict[str, Any]:
        if shift_bp == 0 and index in self._cache:
            item = self._cache.pop(index)
            self._cache[index] = item
            return item
        self._ensure_handles()
        row = self.intervals[index]
        chromosome = row["chromosome"]
        if shift_bp and row["role"] != "train":
            raise ValueError("Random shifts are allowed only for training intervals")
        start = int(row["start"]) + shift_bp
        end = int(row["end"]) + shift_bp
        chromosome_length = self.fai[chromosome][0]
        if start < 0 or end > chromosome_length:
            raise ValueError(
                f"Shift {shift_bp} moves {chromosome}:{start}-{end} out of bounds"
            )
        width = end - start
        if width % 128:
            raise ValueError(f"Window width {width} is not divisible by 128")
        dna = one_hot(fetch_fasta(self.fasta_path, self.fai, chromosome, start, end))
        target_1bp = self._all_signals(chromosome, start, end)
        target_128bp = target_1bp.reshape(len(self.track_rows), width // 128, 128).sum(axis=2)
        chromosome_genes = self.genes.get(chromosome, {"+": [], "-": []})
        gene_mask = np.stack(
            [
                interval_mask(chromosome_genes["+"], start, end),
                interval_mask(chromosome_genes["-"], start, end),
            ],
            axis=0,
        )
        core_mask = np.zeros(width, dtype=bool)
        core_start = int(row["core_start"]) + shift_bp - start
        core_end = int(row["core_end"]) + shift_bp - start
        core_mask[core_start:core_end] = True
        item = {
            "dna_sequence": torch.from_numpy(np.ascontiguousarray(dna)),
            "target_1bp": torch.from_numpy(np.ascontiguousarray(target_1bp)),
            "target_128bp": torch.from_numpy(np.ascontiguousarray(target_128bp)),
            "track_mask": torch.from_numpy(self.track_mask.copy()),
            "track_strand": torch.from_numpy(self.track_strand.copy()),
            "gene_mask": torch.from_numpy(gene_mask),
            "core_mask": torch.from_numpy(core_mask),
            "interval_chromosome": chromosome,
            "interval_start": torch.tensor(start, dtype=torch.long),
            "interval_end": torch.tensor(end, dtype=torch.long),
            "group_ids": tuple(row["group_id"] for row in self.track_rows),
            "shift_bp": shift_bp,
            "reverse_complemented": False,
        }
        if self.cache_size and shift_bp == 0:
            self._cache[index] = item
            while len(self._cache) > self.cache_size:
                self._cache.popitem(last=False)
        return item

    def __getitem__(self, index: int) -> dict[str, Any]:
        return self._load_item(index, shift_bp=0)

    def get_shifted_item(self, index: int, shift_bp: int) -> dict[str, Any]:
        return self._load_item(index, shift_bp=int(shift_bp))

    def __getstate__(self) -> dict[str, Any]:
        state = self.__dict__.copy()
        state["_owner_pid"] = None
        state["_handles"] = None
        state["_executor"] = None
        state["_cache"] = OrderedDict()
        return state

    def close(self) -> None:
        if self._owner_pid == os.getpid():
            self._close_handles()

    def __del__(self) -> None:
        try:
            self.close()
        except Exception:
            pass
