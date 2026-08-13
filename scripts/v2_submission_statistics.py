"""Dependency-light statistics helpers for matched v2 model comparisons.

The helpers deliberately separate the experimental hierarchy from plotting.
Folds are the primary blocks and training seeds are algorithmic repeats within
each block.  They are generic enough to accept future P9 run-level records.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Iterable, Mapping, Sequence

import numpy as np


@dataclass(frozen=True)
class BalancedMatrix:
    """A complete block-by-repeat matrix and its ordered identifiers."""

    values: np.ndarray
    blocks: tuple[Any, ...]
    repeats: tuple[Any, ...]


@dataclass(frozen=True)
class HierarchicalResamplePlan:
    """Indices for fold-first, repeat-within-fold bootstrap resampling."""

    block_indices: np.ndarray
    repeat_indices: np.ndarray
    iterations: int
    seed: int


@dataclass(frozen=True)
class BootstrapEstimate:
    """Fold-primary estimate and percentile confidence interval."""

    estimate: float
    ci_lower: float
    ci_upper: float
    confidence_level: float
    iterations: int
    seed: int


def records_to_balanced_matrix(
    records: Iterable[Mapping[str, Any]],
    *,
    value_key: str,
    block_key: str = "fold",
    repeat_key: str = "seed",
    expected_blocks: Sequence[Any] | None = None,
    expected_repeats: Sequence[Any] | None = None,
) -> BalancedMatrix:
    """Validate and reshape records into a complete block-by-repeat matrix."""

    indexed: dict[tuple[Any, Any], float] = {}
    for record in records:
        key = (record[block_key], record[repeat_key])
        if key in indexed:
            raise ValueError(f"duplicate block/repeat record: {key}")
        value = float(record[value_key])
        if not np.isfinite(value):
            raise ValueError(f"non-finite value for {value_key} at {key}: {value}")
        indexed[key] = value

    if not indexed:
        raise ValueError("no records supplied")

    blocks = tuple(
        expected_blocks
        if expected_blocks is not None
        else sorted({key[0] for key in indexed})
    )
    repeats = tuple(
        expected_repeats
        if expected_repeats is not None
        else sorted({key[1] for key in indexed})
    )
    expected = {(block, repeat) for block in blocks for repeat in repeats}
    actual = set(indexed)
    if actual != expected:
        missing = sorted(expected - actual)
        extra = sorted(actual - expected)
        raise ValueError(f"unbalanced block/repeat matrix; missing={missing}, extra={extra}")

    values = np.asarray(
        [[indexed[(block, repeat)] for repeat in repeats] for block in blocks],
        dtype=np.float64,
    )
    return BalancedMatrix(values=values, blocks=blocks, repeats=repeats)


def matched_differences(
    candidate_records: Iterable[Mapping[str, Any]],
    comparator_records: Iterable[Mapping[str, Any]],
    *,
    value_key: str,
    block_key: str = "fold",
    repeat_key: str = "seed",
) -> list[dict[str, Any]]:
    """Return candidate-minus-comparator values at identical block/repeat keys."""

    def index_records(
        records: Iterable[Mapping[str, Any]], label: str
    ) -> dict[tuple[Any, Any], Mapping[str, Any]]:
        indexed: dict[tuple[Any, Any], Mapping[str, Any]] = {}
        for record in records:
            key = (record[block_key], record[repeat_key])
            if key in indexed:
                raise ValueError(f"duplicate {label} block/repeat record: {key}")
            indexed[key] = record
        return indexed

    candidate = index_records(candidate_records, "candidate")
    comparator = index_records(comparator_records, "comparator")
    if set(candidate) != set(comparator):
        candidate_only = sorted(set(candidate) - set(comparator))
        comparator_only = sorted(set(comparator) - set(candidate))
        raise ValueError(
            "paired keys differ; "
            f"candidate_only={candidate_only}, comparator_only={comparator_only}"
        )

    output = []
    for block, repeat in sorted(candidate):
        candidate_value = float(candidate[(block, repeat)][value_key])
        comparator_value = float(comparator[(block, repeat)][value_key])
        if not np.isfinite(candidate_value) or not np.isfinite(comparator_value):
            raise ValueError(f"non-finite paired value at {(block, repeat)}")
        output.append(
            {
                block_key: block,
                repeat_key: repeat,
                "candidate_value": candidate_value,
                "comparator_value": comparator_value,
                "difference": candidate_value - comparator_value,
            }
        )
    return output


def fold_primary_estimate(values: np.ndarray) -> float:
    """Average repeat means within blocks, then average blocks equally."""

    array = np.asarray(values, dtype=np.float64)
    if array.ndim != 2 or min(array.shape) < 1:
        raise ValueError("values must be a non-empty two-dimensional array")
    if not np.isfinite(array).all():
        raise ValueError("values contain non-finite entries")
    return float(array.mean(axis=1).mean())


def make_hierarchical_resample_plan(
    *,
    n_blocks: int,
    n_repeats: int,
    iterations: int,
    seed: int,
) -> HierarchicalResamplePlan:
    """Create fold-first and seed-within-sampled-fold bootstrap indices."""

    if n_blocks < 1 or n_repeats < 1 or iterations < 1:
        raise ValueError("n_blocks, n_repeats and iterations must be positive")
    rng = np.random.default_rng(seed)
    block_indices = rng.integers(
        0, n_blocks, size=(iterations, n_blocks), endpoint=False
    )
    repeat_indices = rng.integers(
        0,
        n_repeats,
        size=(iterations, n_blocks, n_repeats),
        endpoint=False,
    )
    return HierarchicalResamplePlan(
        block_indices=block_indices,
        repeat_indices=repeat_indices,
        iterations=iterations,
        seed=seed,
    )


def hierarchical_block_bootstrap(
    values: np.ndarray,
    plan: HierarchicalResamplePlan,
    *,
    confidence_level: float = 0.95,
) -> BootstrapEstimate:
    """Estimate a fold-primary mean with a hierarchical percentile interval.

    Each bootstrap replicate samples blocks with replacement.  For every
    sampled block occurrence, repeats are independently sampled with
    replacement.  Repeats are averaged within occurrence before the sampled
    block means are averaged, retaining the block as the primary unit.
    """

    array = np.asarray(values, dtype=np.float64)
    if array.ndim != 2 or min(array.shape) < 1:
        raise ValueError("values must be a non-empty two-dimensional array")
    if not np.isfinite(array).all():
        raise ValueError("values contain non-finite entries")
    if plan.block_indices.shape != (plan.iterations, array.shape[0]):
        raise ValueError("bootstrap block plan does not match values")
    if plan.repeat_indices.shape != (
        plan.iterations,
        array.shape[0],
        array.shape[1],
    ):
        raise ValueError("bootstrap repeat plan does not match values")
    if not 0.0 < confidence_level < 1.0:
        raise ValueError("confidence_level must be between zero and one")

    sampled = array[plan.block_indices[:, :, None], plan.repeat_indices]
    bootstrap_means = sampled.mean(axis=2).mean(axis=1)
    alpha = (1.0 - confidence_level) / 2.0
    lower, upper = np.quantile(bootstrap_means, [alpha, 1.0 - alpha])
    return BootstrapEstimate(
        estimate=fold_primary_estimate(array),
        ci_lower=float(lower),
        ci_upper=float(upper),
        confidence_level=confidence_level,
        iterations=plan.iterations,
        seed=plan.seed,
    )
