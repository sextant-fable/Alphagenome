#!/usr/bin/env python3
"""Auditable, dependency-light primitives for v2 sequence-variant scoring."""

from __future__ import annotations

import csv
from dataclasses import dataclass
import gzip
import math
from pathlib import Path
from typing import Callable, Iterable, Iterator, Mapping, Sequence

import numpy as np

from scripts.v2_biological_validation import AnnotatedGene


VARIANT_SCORING_CONTRACT = "v2_variant_scoring_v1"
EXPECTED_PREDICTION_RESOLUTIONS = frozenset({1, 128})
DNA_ALPHABET = frozenset("ACGTN")
CANONICAL_BASES = "ACGT"
DNA_TO_INDEX = {base: index for index, base in enumerate(CANONICAL_BASES)}
COMPLEMENT = str.maketrans("ACGTN", "TGCAN")


@dataclass(frozen=True)
class VcfAllele:
    chromosome: str
    position_1based: int
    identifier: str
    reference: str
    alternate: str
    quality: str
    filter_value: str
    info: str
    source_line: int
    alternate_index: int

    @property
    def position_0based(self) -> int:
        return self.position_1based - 1

    @property
    def key(self) -> str:
        return (
            f"{self.chromosome}:{self.position_1based}:"
            f"{self.reference}>{self.alternate}:alt{self.alternate_index}"
        )


@dataclass(frozen=True)
class AlleleWindow:
    allele: VcfAllele
    start: int
    end: int
    anchor_index: int
    reference_sequence: str
    alternate_sequence: str
    alignment_policy: str


def supports_coordinate_gene_aggregation(allele: VcfAllele) -> bool:
    """Gene/exon coordinates are comparable only for equal-length alleles."""

    return len(allele.reference) == len(allele.alternate)


def _open_text(path: str | Path):
    path = Path(path)
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open()


def parse_vcf(
    path: str | Path,
    *,
    require_pass: bool = False,
) -> Iterator[VcfAllele]:
    """Yield one record per concrete ALT allele from a VCF 4.2-style file."""

    saw_header = False
    with _open_text(path) as handle:
        for line_number, line in enumerate(handle, 1):
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                saw_header = True
                continue
            if line.startswith("#"):
                continue
            if not saw_header:
                raise ValueError("VCF #CHROM header must precede data records")
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                raise ValueError(f"VCF line {line_number} has fewer than eight fields")
            chromosome, raw_position, identifier, reference, raw_alts, quality, filter_value, info = fields[:8]
            try:
                position = int(raw_position)
            except ValueError as error:
                raise ValueError(f"Invalid VCF POS at line {line_number}: {raw_position}") from error
            if position < 1:
                raise ValueError(f"VCF POS must be one-based and positive at line {line_number}")
            reference = reference.upper()
            if not reference or set(reference) - set(CANONICAL_BASES):
                raise ValueError(f"Unsupported REF bases at line {line_number}: {reference}")
            if require_pass and filter_value not in {"PASS", "."}:
                continue
            for alternate_index, alternate in enumerate(raw_alts.split(","), 1):
                alternate = alternate.upper()
                if (
                    not alternate
                    or alternate == "*"
                    or alternate.startswith("<")
                    or "[" in alternate
                    or "]" in alternate
                ):
                    raise ValueError(
                        f"Symbolic/spanning/breakend ALT is unsupported at line {line_number}: {alternate}"
                    )
                if set(alternate) - set(CANONICAL_BASES):
                    raise ValueError(f"Unsupported ALT bases at line {line_number}: {alternate}")
                if alternate == reference:
                    raise ValueError(f"ALT equals REF at line {line_number}: {alternate}")
                yield VcfAllele(
                    chromosome=chromosome,
                    position_1based=position,
                    identifier=identifier,
                    reference=reference,
                    alternate=alternate,
                    quality=quality,
                    filter_value=filter_value,
                    info=info,
                    source_line=line_number,
                    alternate_index=alternate_index,
                )
    if not saw_header:
        raise ValueError("VCF lacks a #CHROM header")


def read_fai(path: str | Path) -> dict[str, tuple[int, int, int, int]]:
    records: dict[str, tuple[int, int, int, int]] = {}
    with Path(path).open() as handle:
        for line_number, line in enumerate(handle, 1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 5:
                raise ValueError(f"Invalid FAI row at line {line_number}")
            name, length, offset, line_bases, line_width = fields
            if name in records:
                raise ValueError(f"Duplicate FAI contig: {name}")
            records[name] = tuple(map(int, (length, offset, line_bases, line_width)))
    if not records:
        raise ValueError("FAI is empty")
    return records


def fetch_fasta(
    fasta_path: str | Path,
    fai: Mapping[str, tuple[int, int, int, int]],
    chromosome: str,
    start: int,
    end: int,
) -> str:
    if chromosome not in fai:
        raise ValueError(f"VCF chromosome is absent from FAI: {chromosome}")
    length, offset, line_bases, line_width = fai[chromosome]
    if not 0 <= start < end <= length:
        raise ValueError(f"Out-of-bounds reference interval {chromosome}:{start}-{end}")
    sequence_length = end - start
    start_line = start // line_bases
    start_in_line = start % line_bases
    byte_start = offset + start_line * line_width + start_in_line
    line_count = ((start_in_line + sequence_length - 1) // line_bases) + 1
    byte_count = sequence_length + line_count * (line_width - line_bases)
    with Path(fasta_path).open("rb") as handle:
        handle.seek(byte_start)
        raw = handle.read(byte_count)
    sequence = raw.replace(b"\n", b"").replace(b"\r", b"")[:sequence_length].decode("ascii").upper()
    if len(sequence) != sequence_length:
        raise RuntimeError(f"Short FASTA read at {chromosome}:{start}-{end}")
    if set(sequence) - DNA_ALPHABET:
        raise ValueError(f"Reference contains unsupported bases at {chromosome}:{start}-{end}")
    return sequence


def build_allele_window(
    allele: VcfAllele,
    *,
    fasta_path: str | Path,
    fai: Mapping[str, tuple[int, int, int, int]],
    sequence_length: int,
    anchor_index: int | None = None,
) -> AlleleWindow:
    """Build equal-length ref/alt windows with the variant start at one index.

    For an indel, downstream reference context is shifted in the alternate
    haplotype and the right edge is cropped or extended to retain fixed length.
    This is an anchored window-index comparison, not a post-indel genomic
    coordinate realignment.
    """

    if sequence_length <= 0:
        raise ValueError("sequence_length must be positive")
    anchor = sequence_length // 2 if anchor_index is None else int(anchor_index)
    if not 0 <= anchor < sequence_length:
        raise ValueError("anchor_index must lie inside the sequence window")
    if len(allele.reference) > sequence_length - anchor:
        raise ValueError("REF allele does not fit to the right of the fixed anchor")
    position = allele.position_0based
    start = position - anchor
    extra_right = max(0, len(allele.reference) - len(allele.alternate))
    extended_end = start + sequence_length + extra_right
    extended = fetch_fasta(fasta_path, fai, allele.chromosome, start, extended_end)
    observed_reference = extended[anchor : anchor + len(allele.reference)]
    if observed_reference != allele.reference:
        raise ValueError(
            f"Reference allele mismatch for {allele.key}: "
            f"VCF={allele.reference}, FASTA={observed_reference}"
        )
    reference_sequence = extended[:sequence_length]
    alternate_haplotype = (
        extended[:anchor]
        + allele.alternate
        + extended[anchor + len(allele.reference) :]
    )
    alternate_sequence = alternate_haplotype[:sequence_length]
    if len(reference_sequence) != sequence_length or len(alternate_sequence) != sequence_length:
        raise RuntimeError("Failed to construct fixed-length allele windows")
    return AlleleWindow(
        allele=allele,
        start=start,
        end=start + sequence_length,
        anchor_index=anchor,
        reference_sequence=reference_sequence,
        alternate_sequence=alternate_sequence,
        alignment_policy=(
            "coordinate_aligned_substitution"
            if supports_coordinate_gene_aggregation(allele)
            else "fixed_window_index_aligned_diagnostic"
        ),
    )


def one_hot(sequence: str) -> np.ndarray:
    sequence = sequence.upper()
    if set(sequence) - DNA_ALPHABET:
        raise ValueError("Sequence contains unsupported bases")
    encoded = np.zeros((4, len(sequence)), dtype=np.float32)
    for position, base in enumerate(sequence):
        index = DNA_TO_INDEX.get(base)
        if index is not None:
            encoded[index, position] = 1.0
    return encoded


def reverse_complement(sequence: str) -> str:
    sequence = sequence.upper()
    if set(sequence) - DNA_ALPHABET:
        raise ValueError("Sequence contains unsupported bases")
    return sequence.translate(COMPLEMENT)[::-1]


PredictionMap = Mapping[int, np.ndarray]
Predictor = Callable[[np.ndarray], PredictionMap]


def _as_numpy(value: object) -> np.ndarray:
    if hasattr(value, "detach"):
        value = value.detach()
    if hasattr(value, "cpu"):
        value = value.cpu()
    return np.asarray(value)


def validate_prediction_map(
    predictions: PredictionMap,
    *,
    batch_size: int,
) -> dict[int, np.ndarray]:
    result: dict[int, np.ndarray] = {}
    if not predictions:
        raise ValueError("Predictor returned no resolutions")
    for raw_resolution, raw_values in predictions.items():
        resolution = int(raw_resolution)
        values = _as_numpy(raw_values).astype(np.float32, copy=False)
        if resolution <= 0 or values.ndim != 3 or values.shape[0] != batch_size:
            raise ValueError("Predictions must map positive resolution to [batch, tracks, positions]")
        if not np.isfinite(values).all():
            raise ValueError("Predictor returned non-finite values")
        result[resolution] = values
    if set(result) != EXPECTED_PREDICTION_RESOLUTIONS:
        raise ValueError(
            "Predictor must return exactly the 1-bp and 128-bp heads"
        )
    return result


def forward_reverse_complement_ensemble(
    predictor: Predictor,
    sequences: Sequence[str],
) -> dict[int, np.ndarray]:
    """Average forward and orientation-restored RC predictions."""

    if not sequences or len({len(sequence) for sequence in sequences}) != 1:
        raise ValueError("sequences must be non-empty and equal length")
    forward_input = np.stack([one_hot(sequence) for sequence in sequences], axis=0)
    reverse_input = np.stack(
        [one_hot(reverse_complement(sequence)) for sequence in sequences], axis=0
    )
    forward = validate_prediction_map(predictor(forward_input), batch_size=len(sequences))
    reverse = validate_prediction_map(predictor(reverse_input), batch_size=len(sequences))
    if set(forward) != set(reverse):
        raise ValueError("Forward and reverse predictions expose different resolutions")
    result: dict[int, np.ndarray] = {}
    for resolution in forward:
        if forward[resolution].shape != reverse[resolution].shape:
            raise ValueError("Forward and reverse prediction shapes differ")
        result[resolution] = (
            forward[resolution] + reverse[resolution][..., ::-1]
        ) / 2.0
    return result


def score_allele_window(
    predictor: Predictor,
    window: AlleleWindow,
) -> tuple[dict[int, np.ndarray], dict[int, np.ndarray], dict[int, np.ndarray]]:
    predictions = forward_reverse_complement_ensemble(
        predictor, [window.reference_sequence, window.alternate_sequence]
    )
    reference = {resolution: values[0] for resolution, values in predictions.items()}
    alternate = {resolution: values[1] for resolution, values in predictions.items()}
    delta = {
        resolution: alternate[resolution] - reference[resolution]
        for resolution in predictions
    }
    return reference, alternate, delta


def _overlap_bases(intervals: Iterable[tuple[int, int]], start: int, end: int) -> int:
    return sum(max(0, min(end, right) - max(start, left)) for left, right in intervals)


def _interval_delta_sum(
    delta: np.ndarray,
    *,
    interval_start: int,
    interval_end: int,
    window_start: int,
    resolution: int,
) -> tuple[np.ndarray, int]:
    window_end = window_start + delta.shape[1] * resolution
    start = max(interval_start, window_start)
    end = min(interval_end, window_end)
    if end <= start:
        return np.zeros(delta.shape[0], dtype=np.float64), 0
    first_bin = (start - window_start) // resolution
    stop_bin = (end - window_start + resolution - 1) // resolution
    indices = np.arange(first_bin, stop_bin, dtype=np.int64)
    bin_starts = window_start + indices * resolution
    overlaps = np.maximum(
        0,
        np.minimum(bin_starts + resolution, end) - np.maximum(bin_starts, start),
    ).astype(np.float64)
    weighted = delta[:, indices] @ (overlaps / resolution)
    return weighted, int(end - start)


def aggregate_gene_exon_delta(
    delta: np.ndarray,
    *,
    chromosome: str,
    window_start: int,
    genes: Sequence[AnnotatedGene],
    track_ids: Sequence[str],
    resolution: int = 1,
) -> list[dict[str, object]]:
    """Aggregate positional deltas over gene bodies and unique exon unions."""

    delta = np.asarray(delta, dtype=np.float64)
    if delta.ndim != 2 or delta.shape[0] != len(track_ids):
        raise ValueError("delta must have [tracks, positions] matching track_ids")
    if resolution <= 0:
        raise ValueError("resolution must be positive")
    window_end = window_start + delta.shape[1] * resolution
    rows: list[dict[str, object]] = []
    for gene in genes:
        if gene.chromosome != chromosome or gene.gene_end <= window_start or gene.gene_start >= window_end:
            continue
        body_sum, body_bases = _interval_delta_sum(
            delta,
            interval_start=gene.gene_start,
            interval_end=gene.gene_end,
            window_start=window_start,
            resolution=resolution,
        )
        exon_sum = np.zeros(len(track_ids), dtype=np.float64)
        exon_bases = 0
        for exon_start, exon_end in gene.exon_intervals:
            current_sum, current_bases = _interval_delta_sum(
                delta,
                interval_start=exon_start,
                interval_end=exon_end,
                window_start=window_start,
                resolution=resolution,
            )
            exon_sum += current_sum
            exon_bases += current_bases
        for track_index, track_id in enumerate(track_ids):
            rows.append(
                {
                    "gene_id": gene.gene_id,
                    "gene_name": gene.gene_name,
                    "chromosome": gene.chromosome,
                    "strand": gene.strand,
                    "track_id": track_id,
                    "resolution": resolution,
                    "gene_body_bases": body_bases,
                    "exon_bases": exon_bases,
                    "gene_body_delta_sum": float(body_sum[track_index]),
                    "gene_body_delta_mean": (
                        float(body_sum[track_index] / body_bases) if body_bases else math.nan
                    ),
                    "exon_delta_sum": float(exon_sum[track_index]),
                    "exon_delta_mean": (
                        float(exon_sum[track_index] / exon_bases) if exon_bases else math.nan
                    ),
                }
            )
    return rows


def aggregate_variant_gene_exon_delta(
    window: AlleleWindow,
    delta: np.ndarray,
    *,
    genes: Sequence[AnnotatedGene],
    track_ids: Sequence[str],
    resolution: int = 1,
) -> list[dict[str, object]]:
    """Return coordinate-aligned gene effects; indels deliberately return none."""

    if not supports_coordinate_gene_aggregation(window.allele):
        return []
    return aggregate_gene_exon_delta(
        delta,
        chromosome=window.allele.chromosome,
        window_start=window.start,
        genes=genes,
        track_ids=track_ids,
        resolution=resolution,
    )


@dataclass(frozen=True)
class IsmMutation:
    sequence_index: int
    genomic_position_1based: int
    reference: str
    alternate: str
    sequence: str


def enumerate_ism_mutations(
    sequence: str,
    *,
    window_start: int,
    positions: Iterable[int] | None = None,
) -> Iterator[IsmMutation]:
    """Enumerate all three canonical substitutions at selected window indices."""

    sequence = sequence.upper()
    selected = range(len(sequence)) if positions is None else positions
    seen: set[int] = set()
    for raw_index in selected:
        index = int(raw_index)
        if index in seen:
            continue
        seen.add(index)
        if not 0 <= index < len(sequence):
            raise ValueError(f"ISM index outside sequence: {index}")
        reference = sequence[index]
        if reference not in CANONICAL_BASES:
            continue
        for alternate in CANONICAL_BASES:
            if alternate == reference:
                continue
            yield IsmMutation(
                sequence_index=index,
                genomic_position_1based=window_start + index + 1,
                reference=reference,
                alternate=alternate,
                sequence=sequence[:index] + alternate + sequence[index + 1 :],
            )


def score_ism(
    predictor: Predictor,
    *,
    reference_sequence: str,
    window_start: int,
    positions: Iterable[int] | None = None,
    max_mutants: int = 10_000,
) -> Iterator[tuple[IsmMutation, dict[int, np.ndarray]]]:
    if max_mutants <= 0:
        raise ValueError("max_mutants must be positive")
    reference_prediction = forward_reverse_complement_ensemble(
        predictor, [reference_sequence]
    )
    mutations = enumerate_ism_mutations(
        reference_sequence, window_start=window_start, positions=positions
    )
    for ordinal, mutation in enumerate(mutations, 1):
        if ordinal > max_mutants:
            raise RuntimeError(f"ISM exceeds max_mutants={max_mutants}")
        alternate_prediction = forward_reverse_complement_ensemble(
            predictor, [mutation.sequence]
        )
        yield mutation, {
            resolution: alternate_prediction[resolution][0] - values[0]
            for resolution, values in reference_prediction.items()
        }


def summarize_track_deltas(
    delta: np.ndarray,
    track_ids: Sequence[str],
    *,
    resolution: int,
) -> list[dict[str, object]]:
    delta = np.asarray(delta, dtype=np.float64)
    if delta.ndim != 2 or delta.shape[0] != len(track_ids):
        raise ValueError("delta must have [tracks, positions] matching track_ids")
    rows = []
    for index, track_id in enumerate(track_ids):
        values = delta[index]
        maximum_index = int(np.argmax(np.abs(values)))
        rows.append(
            {
                "track_id": track_id,
                "resolution": resolution,
                "delta_sum": float(values.sum()),
                "delta_mean": float(values.mean()),
                "delta_l1": float(np.abs(values).sum()),
                "max_abs_delta": float(abs(values[maximum_index])),
                "max_abs_delta_index": maximum_index,
            }
        )
    return rows
