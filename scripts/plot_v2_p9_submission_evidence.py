#!/usr/bin/env python3
"""Render the R9-PASS P9 component-attribution submission figure on CPU."""

from __future__ import annotations

import argparse
from collections import defaultdict
import csv
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import socket
import subprocess
import sys
from typing import Any, Mapping, Sequence

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image


REPO_ROOT = Path(__file__).resolve().parents[1]
METADATA_DIR = REPO_ROOT / "alphagenome_custom/metadata/v2"
DEFAULT_REVIEW = METADATA_DIR / "audits/P9/review.json"
DEFAULT_PAIRED = METADATA_DIR / "p9_submission_evidence_paired_effects.tsv"
DEFAULT_RAW_PAIRED = METADATA_DIR / "p9_submission_evidence_raw_paired_effects.tsv"
DEFAULT_FACTORIAL = METADATA_DIR / "p9_submission_evidence_factorial_effects.tsv"
DEFAULT_RAW_FACTORIAL = METADATA_DIR / "p9_submission_evidence_raw_factorial_effects.tsv"
DEFAULT_PER_TRACK = METADATA_DIR / "p9_submission_evidence_per_track.tsv"
DEFAULT_OUTPUT_DIR = REPO_ROOT / "results/v2_p9_submission_evidence_figure"

FIGURE_STEM = "figure_p9_component_attribution"
FIG_WIDTH_MM = 183
FIG_HEIGHT_MM = 112
RASTER_DPI = 600
FONT_SIZE_PT = 7
PRIMARY_METRIC = "primary_biological_score"
EXPECTED_FOLDS = (1, 2, 3, 4, 5)
EXPECTED_SEEDS = (20260714, 20260715, 20260716)
REFERENCE_CONFIGURATION = "full_worm_lora"
REFERENCE_DISPLAY = "Full B"
CONFIGURATION_ORDER = (
    "size_matched_scratch",
    "frozen_no_worm_no_lora",
    "no_lora",
    "lora_only",
    "learned_1bp_head_only",
)
CONFIGURATION_LABELS = {
    "size_matched_scratch": "Size-matched C",
    "frozen_no_worm_no_lora": "Frozen A",
    "no_lora": "B without LoRA",
    "lora_only": "LoRA only",
    "learned_1bp_head_only": "1-bp head only",
}
EFFECT_ORDER = ("lora_main", "worm_embedding_main", "worm_lora_interaction")
EFFECT_LABELS = {
    "lora_main": "LoRA main effect",
    "worm_embedding_main": "Worm embedding main effect",
    "worm_lora_interaction": "Worm x LoRA interaction",
}
COLORS = {
    "size_matched_scratch": "#C4473F",
    "frozen_no_worm_no_lora": "#6B6B6B",
    "no_lora": "#D08A3C",
    "lora_only": "#177E89",
    "learned_1bp_head_only": "#6D6AA8",
    "lora_main": "#C4473F",
    "worm_embedding_main": "#177E89",
    "worm_lora_interaction": "#6B6B6B",
    "fold": "#B8B8B8",
    "reference": "#303030",
    "grid": "#D8D8D8",
}
SOURCE_FIELDS = (
    "panel",
    "row_type",
    "configuration",
    "effect",
    "display_label",
    "metric",
    "fold",
    "seed",
    "track_id",
    "estimate",
    "ci_low",
    "ci_high",
    "n_folds",
    "seeds_per_fold",
    "inference_role",
    "notes",
)


mpl.rcParams.update(
    {
        "font.family": "Liberation Sans",
        "font.sans-serif": ["Liberation Sans", "Arial", "DejaVu Sans"],
        "font.size": FONT_SIZE_PT,
        "axes.labelsize": FONT_SIZE_PT,
        "axes.titlesize": 8,
        "xtick.labelsize": 6,
        "ytick.labelsize": 6,
        "legend.fontsize": 6,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
        "axes.linewidth": 0.7,
        "xtick.major.width": 0.7,
        "ytick.major.width": 0.7,
    }
)


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def sha256(path: Path, chunk_size: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def read_json(path: Path) -> dict[str, Any]:
    payload = json.loads(path.read_text())
    if not isinstance(payload, dict):
        raise ValueError(f"Expected a JSON object: {path}")
    return payload


def read_tsv(path: Path, required: set[str]) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows or required - set(rows[0]):
        raise ValueError(f"Missing required columns in {path.name}: {required - set(rows[0] if rows else ())}")
    return rows


def finite_float(value: object, label: str) -> float:
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f"Non-finite value for {label}")
    return result


def _require_exact_cells(
    rows: Sequence[Mapping[str, str]],
    *,
    group_key: str,
    expected_groups: Sequence[str],
) -> None:
    observed = {
        (str(row[group_key]), int(row["fold"]), int(row["seed"])) for row in rows
    }
    expected = {
        (group, fold, seed)
        for group in expected_groups
        for fold in EXPECTED_FOLDS
        for seed in EXPECTED_SEEDS
    }
    if len(observed) != len(rows) or observed != expected:
        raise ValueError(f"Unbalanced {group_key} x fold x seed contract")


def load_figure_inputs(
    *,
    review_path: Path = DEFAULT_REVIEW,
    paired_path: Path = DEFAULT_PAIRED,
    raw_paired_path: Path = DEFAULT_RAW_PAIRED,
    factorial_path: Path = DEFAULT_FACTORIAL,
    raw_factorial_path: Path = DEFAULT_RAW_FACTORIAL,
    per_track_path: Path = DEFAULT_PER_TRACK,
) -> dict[str, Any]:
    """Load and fail closed on the formal R9-PASS figure inputs."""

    paths = {
        "review": review_path,
        "paired": paired_path,
        "raw_paired": raw_paired_path,
        "factorial": factorial_path,
        "raw_factorial": raw_factorial_path,
        "per_track": per_track_path,
    }
    for path in paths.values():
        if not path.is_file():
            raise FileNotFoundError(path)
        lowered = path.name.lower()
        if "final_test" in lowered or "test_locked" in lowered:
            raise ValueError("The P9 figure must not consume final-test artifacts")

    review = read_json(review_path)
    checks = review.get("checks", [])
    if not (
        review.get("phase") == "P9"
        and review.get("status") == "PASS"
        and len(checks) == 9
        and all(item.get("status") == "PASS" for item in checks)
    ):
        raise RuntimeError("The component-attribution figure requires a formal R9 PASS")

    paired = read_tsv(
        paired_path,
        {
            "configuration",
            "reference",
            "metric",
            "paired_runs",
            "folds",
            "seeds",
            "mean_paired_difference",
            "ci_low",
            "ci_high",
            "ci_method",
        },
    )
    raw_paired = read_tsv(
        raw_paired_path,
        {
            "configuration",
            "reference",
            "metric",
            "fold",
            "seed",
            "candidate_value",
            "reference_value",
            "difference",
        },
    )
    factorial = read_tsv(
        factorial_path,
        {
            "effect",
            "metric",
            "paired_runs",
            "folds",
            "seeds",
            "mean_effect",
            "ci_low",
            "ci_high",
            "ci_method",
        },
    )
    raw_factorial = read_tsv(
        raw_factorial_path, {"effect", "metric", "fold", "seed", "value"}
    )
    per_track = read_tsv(
        per_track_path,
        {"configuration", "fold", "seed", "track_id", "metric", "value"},
    )
    if not (
        len(paired) == 66
        and len(raw_paired) == 990
        and len(factorial) == 33
        and len(raw_factorial) == 495
        and len(per_track) == 65070
    ):
        raise ValueError("Formal P9 table row counts changed")
    return {
        "review": review,
        "paired": paired,
        "raw_paired": raw_paired,
        "factorial": factorial,
        "raw_factorial": raw_factorial,
        "per_track": per_track,
        "paths": paths,
    }


def prepare_figure_data(inputs: Mapping[str, Any]) -> dict[str, Any]:
    """Derive all plotted values without treating tracks or seeds as independent n."""

    paired_summary = [
        row
        for row in inputs["paired"]
        if row["metric"] == PRIMARY_METRIC
        and row["configuration"] in CONFIGURATION_ORDER
    ]
    if {row["configuration"] for row in paired_summary} != set(CONFIGURATION_ORDER):
        raise ValueError("Panel a lacks one or more registered configurations")
    paired_by_configuration = {row["configuration"]: row for row in paired_summary}
    for configuration, row in paired_by_configuration.items():
        if not (
            row["reference"] == "B/paper"
            and int(row["paired_runs"]) == 15
            and int(row["folds"]) == 5
            and int(row["seeds"]) == 3
            and finite_float(row["ci_low"], "paired ci_low")
            <= finite_float(row["mean_paired_difference"], "paired estimate")
            <= finite_float(row["ci_high"], "paired ci_high")
        ):
            raise ValueError(f"Invalid paired summary for {configuration}")

    paired_raw = [
        row
        for row in inputs["raw_paired"]
        if row["metric"] == PRIMARY_METRIC
        and row["configuration"] in CONFIGURATION_ORDER
    ]
    _require_exact_cells(
        paired_raw,
        group_key="configuration",
        expected_groups=CONFIGURATION_ORDER,
    )
    paired_fold_means: list[dict[str, Any]] = []
    for configuration in CONFIGURATION_ORDER:
        for fold in EXPECTED_FOLDS:
            members = [
                row
                for row in paired_raw
                if row["configuration"] == configuration and int(row["fold"]) == fold
            ]
            paired_fold_means.append(
                {
                    "configuration": configuration,
                    "fold": fold,
                    "estimate": float(
                        np.mean(
                            [finite_float(row["difference"], "paired difference") for row in members]
                        )
                    ),
                }
            )

    factorial_summary = [
        row
        for row in inputs["factorial"]
        if row["metric"] == PRIMARY_METRIC and row["effect"] in EFFECT_ORDER
    ]
    if {row["effect"] for row in factorial_summary} != set(EFFECT_ORDER):
        raise ValueError("Panel b lacks one or more registered factorial effects")
    factorial_by_effect = {row["effect"]: row for row in factorial_summary}
    factorial_raw = [
        row
        for row in inputs["raw_factorial"]
        if row["metric"] == PRIMARY_METRIC and row["effect"] in EFFECT_ORDER
    ]
    _require_exact_cells(
        factorial_raw,
        group_key="effect",
        expected_groups=EFFECT_ORDER,
    )
    factorial_fold_means: list[dict[str, Any]] = []
    for effect in EFFECT_ORDER:
        for fold in EXPECTED_FOLDS:
            members = [
                row
                for row in factorial_raw
                if row["effect"] == effect and int(row["fold"]) == fold
            ]
            factorial_fold_means.append(
                {
                    "effect": effect,
                    "fold": fold,
                    "estimate": float(
                        np.mean([finite_float(row["value"], "factorial value") for row in members])
                    ),
                }
            )

    primary_track_rows = [
        row for row in inputs["per_track"] if row["metric"] == PRIMARY_METRIC
    ]
    expected_configurations = (*CONFIGURATION_ORDER, REFERENCE_CONFIGURATION)
    tracks = sorted({row["track_id"] for row in primary_track_rows})
    if len(tracks) != 241:
        raise ValueError(f"Expected 241 tracks, observed {len(tracks)}")
    keyed: dict[tuple[str, int, int, str], float] = {}
    for row in primary_track_rows:
        configuration = row["configuration"]
        if configuration not in expected_configurations:
            continue
        key = (
            configuration,
            int(row["fold"]),
            int(row["seed"]),
            row["track_id"],
        )
        if key in keyed:
            raise ValueError(f"Duplicate per-track cell: {key}")
        keyed[key] = finite_float(row["value"], "per-track value")
    expected_keys = {
        (configuration, fold, seed, track)
        for configuration in expected_configurations
        for fold in EXPECTED_FOLDS
        for seed in EXPECTED_SEEDS
        for track in tracks
    }
    if set(keyed) != expected_keys:
        raise ValueError("Per-track primary-score matrix is incomplete")

    track_effects: list[dict[str, Any]] = []
    for configuration in CONFIGURATION_ORDER:
        for track in tracks:
            fold_means = []
            for fold in EXPECTED_FOLDS:
                differences = [
                    keyed[(configuration, fold, seed, track)]
                    - keyed[(REFERENCE_CONFIGURATION, fold, seed, track)]
                    for seed in EXPECTED_SEEDS
                ]
                fold_means.append(float(np.mean(differences)))
            track_effects.append(
                {
                    "configuration": configuration,
                    "track_id": track,
                    "estimate": float(np.mean(fold_means)),
                }
            )
    return {
        "paired_summary": paired_by_configuration,
        "paired_raw": paired_raw,
        "paired_fold_means": paired_fold_means,
        "factorial_summary": factorial_by_effect,
        "factorial_raw": factorial_raw,
        "factorial_fold_means": factorial_fold_means,
        "track_effects": track_effects,
        "track_count": len(tracks),
    }


def build_source_rows(data: Mapping[str, Any]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    common = {
        "metric": PRIMARY_METRIC,
        "n_folds": 5,
        "seeds_per_fold": 3,
    }
    for row in data["paired_raw"]:
        configuration = row["configuration"]
        rows.append(
            {
                "panel": "a_paired_primary",
                "row_type": "seed_pair",
                "configuration": configuration,
                "effect": "",
                "display_label": CONFIGURATION_LABELS[configuration],
                "fold": row["fold"],
                "seed": row["seed"],
                "track_id": "",
                "estimate": row["difference"],
                "ci_low": "",
                "ci_high": "",
                "inference_role": "nested_algorithmic_repeat",
                "notes": "candidate minus B/paper",
                **common,
            }
        )
    for row in data["paired_fold_means"]:
        configuration = row["configuration"]
        rows.append(
            {
                "panel": "a_paired_primary",
                "row_type": "fold_seed_mean",
                "configuration": configuration,
                "effect": "",
                "display_label": CONFIGURATION_LABELS[configuration],
                "fold": row["fold"],
                "seed": "seed_mean",
                "track_id": "",
                "estimate": row["estimate"],
                "ci_low": "",
                "ci_high": "",
                "inference_role": "independent_block",
                "notes": "mean of three paired seed effects within fold",
                **common,
            }
        )
    for configuration in CONFIGURATION_ORDER:
        source = data["paired_summary"][configuration]
        rows.append(
            {
                "panel": "a_paired_primary",
                "row_type": "hierarchical_summary",
                "configuration": configuration,
                "effect": "",
                "display_label": CONFIGURATION_LABELS[configuration],
                "fold": "all_folds",
                "seed": "nested_within_fold",
                "track_id": "",
                "estimate": source["mean_paired_difference"],
                "ci_low": source["ci_low"],
                "ci_high": source["ci_high"],
                "inference_role": "fold_primary_95_ci",
                "notes": source["ci_method"],
                **common,
            }
        )

    for row in data["factorial_raw"]:
        effect = row["effect"]
        rows.append(
            {
                "panel": "b_factorial_primary",
                "row_type": "seed_effect",
                "configuration": "",
                "effect": effect,
                "display_label": EFFECT_LABELS[effect],
                "fold": row["fold"],
                "seed": row["seed"],
                "track_id": "",
                "estimate": row["value"],
                "ci_low": "",
                "ci_high": "",
                "inference_role": "nested_algorithmic_repeat",
                "notes": "registered 2 x 2 factorial contrast",
                **common,
            }
        )
    for row in data["factorial_fold_means"]:
        effect = row["effect"]
        rows.append(
            {
                "panel": "b_factorial_primary",
                "row_type": "fold_seed_mean",
                "configuration": "",
                "effect": effect,
                "display_label": EFFECT_LABELS[effect],
                "fold": row["fold"],
                "seed": "seed_mean",
                "track_id": "",
                "estimate": row["estimate"],
                "ci_low": "",
                "ci_high": "",
                "inference_role": "independent_block",
                "notes": "mean of three seed effects within fold",
                **common,
            }
        )
    for effect in EFFECT_ORDER:
        source = data["factorial_summary"][effect]
        rows.append(
            {
                "panel": "b_factorial_primary",
                "row_type": "hierarchical_summary",
                "configuration": "",
                "effect": effect,
                "display_label": EFFECT_LABELS[effect],
                "fold": "all_folds",
                "seed": "nested_within_fold",
                "track_id": "",
                "estimate": source["mean_effect"],
                "ci_low": source["ci_low"],
                "ci_high": source["ci_high"],
                "inference_role": "fold_primary_95_ci",
                "notes": source["ci_method"],
                **common,
            }
        )

    for row in data["track_effects"]:
        configuration = row["configuration"]
        rows.append(
            {
                "panel": "c_track_distribution",
                "row_type": "track_fold_primary_effect",
                "configuration": configuration,
                "effect": "",
                "display_label": CONFIGURATION_LABELS[configuration],
                "fold": "all_folds",
                "seed": "nested_within_fold",
                "track_id": row["track_id"],
                "estimate": row["estimate"],
                "ci_low": "",
                "ci_high": "",
                "inference_role": "descriptive_track_outcome",
                "notes": "tracks are outcomes, not independent replicates",
                **common,
            }
        )
    expected_rows = 75 + 25 + 5 + 45 + 15 + 3 + 5 * 241
    if len(rows) != expected_rows:
        raise RuntimeError(f"Unexpected Figure Source Data size: {len(rows)}")
    return rows


def write_tsv(path: Path, rows: Sequence[Mapping[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=SOURCE_FIELDS, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def _panel_label(axis: plt.Axes, label: str) -> None:
    axis.text(
        -0.15,
        1.08,
        label,
        transform=axis.transAxes,
        fontsize=9,
        fontweight="bold",
        va="top",
    )


def _stable_jitter(identifier: str, half_width: float = 0.19) -> float:
    code = int(hashlib.sha256(identifier.encode("utf-8")).hexdigest()[:8], 16)
    return ((code / float(0xFFFFFFFF)) - 0.5) * 2.0 * half_width


def _forest_panel(
    axis: plt.Axes,
    *,
    order: Sequence[str],
    labels: Mapping[str, str],
    summaries: Mapping[str, Mapping[str, str]],
    fold_rows: Sequence[Mapping[str, Any]],
    group_key: str,
    estimate_key: str,
    title: str,
    xlabel: str,
) -> None:
    positions = np.arange(len(order), dtype=float)
    fold_offsets = np.linspace(-0.13, 0.13, len(EXPECTED_FOLDS))
    all_extents: list[float] = [0.0]
    for position, group in zip(positions, order):
        color = COLORS[group]
        members = sorted(
            [row for row in fold_rows if row[group_key] == group],
            key=lambda row: int(row["fold"]),
        )
        if len(members) != 5:
            raise ValueError(f"Forest panel requires five fold means for {group}")
        fold_values = np.array([finite_float(row["estimate"], "fold effect") for row in members])
        axis.scatter(
            fold_values,
            position + fold_offsets,
            s=10,
            facecolors="white",
            edgecolors=color,
            linewidths=0.6,
            alpha=0.86,
            zorder=2,
        )
        summary = summaries[group]
        estimate = finite_float(summary[estimate_key], "summary estimate")
        lower = finite_float(summary["ci_low"], "summary ci_low")
        upper = finite_float(summary["ci_high"], "summary ci_high")
        all_extents.extend([*fold_values.tolist(), lower, upper])
        open_marker = group == "learned_1bp_head_only"
        axis.errorbar(
            estimate,
            position,
            xerr=np.array([[estimate - lower], [upper - estimate]]),
            fmt="D",
            markersize=4.6,
            markerfacecolor="white" if open_marker else color,
            markeredgecolor=color,
            markeredgewidth=0.9,
            color=color,
            capsize=2.2,
            linewidth=1.3,
            zorder=4,
        )
    data_min, data_max = min(all_extents), max(all_extents)
    span = max(data_max - data_min, 0.02)
    axis.set_xlim(data_min - 0.08 * span, data_max + 0.16 * span)
    axis.axvline(0, color=COLORS["reference"], linewidth=0.75, linestyle="--", zorder=1)
    axis.set_yticks(positions, [labels[group] for group in order])
    axis.set_ylim(len(order) - 0.55, -0.55)
    axis.set_xlabel(xlabel)
    axis.set_title(title, loc="left", fontweight="bold", pad=4)
    axis.grid(axis="x", color=COLORS["grid"], linewidth=0.45, alpha=0.65)
    axis.set_axisbelow(True)
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    axis.text(
        0.0,
        -0.31,
        "Circles: fold means; diamonds: estimate and 95% CI",
        transform=axis.transAxes,
        fontsize=5.5,
        color=COLORS["reference"],
        va="top",
    )


def build_figure(data: Mapping[str, Any]) -> plt.Figure:
    fig = plt.figure(
        figsize=(FIG_WIDTH_MM / 25.4, FIG_HEIGHT_MM / 25.4),
        facecolor="white",
    )
    grid = fig.add_gridspec(
        2,
        2,
        height_ratios=(1.0, 1.18),
        width_ratios=(1.22, 0.88),
        hspace=0.62,
        wspace=0.48,
    )
    ax_a = fig.add_subplot(grid[0, 0])
    ax_b = fig.add_subplot(grid[0, 1])
    ax_c = fig.add_subplot(grid[1, :])

    _forest_panel(
        ax_a,
        order=CONFIGURATION_ORDER,
        labels=CONFIGURATION_LABELS,
        summaries=data["paired_summary"],
        fold_rows=data["paired_fold_means"],
        group_key="configuration",
        estimate_key="mean_paired_difference",
        title="Matched model and component comparisons",
        xlabel="Paired primary-score effect (candidate - B; higher is better)",
    )
    head_position = CONFIGURATION_ORDER.index("learned_1bp_head_only")
    head_summary = data["paired_summary"]["learned_1bp_head_only"]
    head_upper = finite_float(head_summary["ci_high"], "1-bp head ci_high")
    if not (
        finite_float(head_summary["ci_low"], "1-bp head ci_low")
        <= 0.0
        <= head_upper
    ):
        raise ValueError("The registered 1-bp-head-only primary CI no longer includes zero")
    ax_a.annotate(
        "95% CI includes 0",
        xy=(head_upper, head_position),
        xytext=(8, 0),
        textcoords="offset points",
        fontsize=5.5,
        color=COLORS["learned_1bp_head_only"],
        ha="left",
        va="center",
    )
    _panel_label(ax_a, "a")

    _forest_panel(
        ax_b,
        order=EFFECT_ORDER,
        labels=EFFECT_LABELS,
        summaries=data["factorial_summary"],
        fold_rows=data["factorial_fold_means"],
        group_key="effect",
        estimate_key="mean_effect",
        title="2 x 2 factorial attribution",
        xlabel="Factorial effect on primary score",
    )
    ax_b.text(
        0.0,
        -0.45,
        "Interaction is statistical non-additivity;\nnot a mechanistic interaction.",
        transform=ax_b.transAxes,
        fontsize=5.5,
        color=COLORS["reference"],
        va="top",
        wrap=True,
    )
    _panel_label(ax_b, "b")

    positions = np.arange(len(CONFIGURATION_ORDER), dtype=float)
    all_track_values: list[float] = []
    for position, configuration in zip(positions, CONFIGURATION_ORDER):
        members = [
            row for row in data["track_effects"] if row["configuration"] == configuration
        ]
        if len(members) != 241:
            raise ValueError(f"Panel c requires 241 track effects for {configuration}")
        members.sort(key=lambda row: row["track_id"])
        values = np.array([finite_float(row["estimate"], "track effect") for row in members])
        all_track_values.extend(values.tolist())
        violin = ax_c.violinplot(
            [values],
            positions=[position],
            vert=False,
            widths=0.68,
            showmeans=False,
            showmedians=False,
            showextrema=False,
            points=100,
        )
        for body in violin["bodies"]:
            body.set_facecolor(COLORS[configuration])
            body.set_edgecolor(COLORS[configuration])
            body.set_linewidth(0.7)
            body.set_alpha(0.18)
        jitter = np.array(
            [_stable_jitter(f"{configuration}:{row['track_id']}") for row in members]
        )
        ax_c.scatter(
            values,
            position + jitter,
            s=4,
            color=COLORS[configuration],
            alpha=0.30,
            edgecolors="none",
            rasterized=True,
            zorder=2,
        )
        lower, median, upper = np.quantile(values, [0.25, 0.5, 0.75])
        ax_c.plot([lower, upper], [position, position], color=COLORS["reference"], linewidth=1.8, zorder=4)
        ax_c.scatter(
            median,
            position,
            marker="|",
            s=55,
            linewidths=1.2,
            color=COLORS["reference"],
            zorder=5,
        )
    span = max(max(all_track_values) - min(all_track_values), 0.02)
    ax_c.set_xlim(
        min(min(all_track_values), 0.0) - 0.04 * span,
        max(max(all_track_values), 0.0) + 0.04 * span,
    )
    ax_c.axvline(0, color=COLORS["reference"], linewidth=0.75, linestyle="--", zorder=1)
    ax_c.set_yticks(positions, [CONFIGURATION_LABELS[item] for item in CONFIGURATION_ORDER])
    ax_c.set_ylim(len(CONFIGURATION_ORDER) - 0.55, -0.55)
    ax_c.set_xlabel("Track-wise paired primary-score effect (candidate - B; higher is better)")
    ax_c.set_title("Across-track heterogeneity", loc="left", fontweight="bold", pad=4)
    ax_c.grid(axis="x", color=COLORS["grid"], linewidth=0.45, alpha=0.65)
    ax_c.set_axisbelow(True)
    ax_c.spines["top"].set_visible(False)
    ax_c.spines["right"].set_visible(False)
    ax_c.text(
        0.0,
        -0.31,
        "Each point is one of 241 tracks after fold/seed aggregation; distributions are descriptive, not independent n.",
        transform=ax_c.transAxes,
        fontsize=5.5,
        color=COLORS["reference"],
        va="top",
    )
    _panel_label(ax_c, "c")

    fig.suptitle(
        "Component attribution under matched genomic-block validation",
        x=0.055,
        y=0.985,
        ha="left",
        fontsize=9,
        fontweight="bold",
    )
    fig.text(
        0.055,
        0.012,
        "n = 5 genomic folds; three seeds are nested within fold. Hierarchical percentile bootstrap: 10,000 resamples.",
        ha="left",
        va="bottom",
        fontsize=5.7,
        color=COLORS["reference"],
    )
    fig.subplots_adjust(left=0.17, right=0.985, bottom=0.19, top=0.91)
    return fig


def export_figure(
    figure: plt.Figure,
    output_dir: Path,
    *,
    raster_dpi: int = RASTER_DPI,
) -> dict[str, Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    paths = {
        "pdf": output_dir / f"{FIGURE_STEM}.pdf",
        "svg": output_dir / f"{FIGURE_STEM}.svg",
        "png": output_dir / f"{FIGURE_STEM}.png",
        "tiff": output_dir / f"{FIGURE_STEM}.tiff",
    }
    figure.savefig(paths["pdf"], facecolor="white")
    figure.savefig(paths["svg"], facecolor="white")
    figure.savefig(paths["png"], dpi=raster_dpi, facecolor="white")
    figure.savefig(
        paths["tiff"],
        dpi=raster_dpi,
        facecolor="white",
        pil_kwargs={"compression": "tiff_lzw"},
    )
    for path in paths.values():
        if not path.is_file() or path.stat().st_size <= 0:
            raise RuntimeError(f"Figure export is missing or empty: {path.name}")
    return paths


def raster_pixel_audit(path: Path, expected_dpi: int) -> dict[str, Any]:
    with Image.open(path) as image:
        rgb = np.asarray(image.convert("RGB"), dtype=np.uint8)
        dpi = image.info.get("dpi", (None, None))
        width, height = image.size
    expected_width = round(FIG_WIDTH_MM / 25.4 * expected_dpi)
    expected_height = round(FIG_HEIGHT_MM / 25.4 * expected_dpi)
    nonwhite_fraction = float(np.mean(np.any(rgb < 248, axis=2)))
    channel_std = [float(value) for value in rgb.reshape(-1, 3).std(axis=0)]
    border = np.concatenate(
        (rgb[0, :, :], rgb[-1, :, :], rgb[:, 0, :], rgb[:, -1, :]), axis=0
    )
    clean_border_fraction = float(np.mean(np.all(border >= 245, axis=1)))
    passed = (
        abs(width - expected_width) <= 2
        and abs(height - expected_height) <= 2
        and nonwhite_fraction > 0.015
        and max(channel_std) > 12.0
        and clean_border_fraction > 0.97
    )
    return {
        "passed": passed,
        "width_pixels": width,
        "height_pixels": height,
        "expected_width_pixels": expected_width,
        "expected_height_pixels": expected_height,
        "dpi": [float(value) if value is not None else None for value in dpi],
        "nonwhite_fraction": nonwhite_fraction,
        "channel_std": channel_std,
        "clean_border_fraction": clean_border_fraction,
    }


def _run_json(command: list[str]) -> dict[str, Any]:
    completed = subprocess.run(command, check=False, capture_output=True, text=True)
    if completed.returncode != 0:
        raise RuntimeError(
            f"QA command failed ({completed.returncode}): {' '.join(command)}\n"
            f"{completed.stdout}\n{completed.stderr}"
        )
    payload = json.loads(completed.stdout)
    if not isinstance(payload, dict):
        raise RuntimeError("QA command did not return a JSON object")
    return payload


def _display_path(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(REPO_ROOT.resolve()))
    except ValueError:
        return path.name


def artifact_record(path: Path) -> dict[str, Any]:
    return {
        "path": _display_path(path),
        "sha256": sha256(path),
        "size_bytes": path.stat().st_size,
    }


def git_commit() -> str:
    return subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=REPO_ROOT,
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()


def write_json(path: Path, payload: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)


def render_bundle(
    inputs: Mapping[str, Any],
    *,
    output_dir: Path,
    skill_dir: Path,
    manual_visual_pass: bool,
    visual_review_note: str,
) -> dict[str, Any]:
    data = prepare_figure_data(inputs)
    source_rows = build_source_rows(data)
    source_path = output_dir / "figure_source_data.tsv"
    write_tsv(source_path, source_rows)
    figure = build_figure(data)
    try:
        figure_paths = export_figure(figure, output_dir)
    finally:
        plt.close(figure)

    validator = skill_dir / "scripts/validate_figure.py"
    pdf_auditor = skill_dir / "scripts/audit_pdf_text.py"
    if not validator.is_file() or not pdf_auditor.is_file():
        raise FileNotFoundError("Nature figure QA tools are unavailable")
    static_preflight = _run_json(
        [
            sys.executable,
            str(validator),
            str(Path(__file__).resolve()),
            "--json",
            "--strict",
        ]
    )
    pdf_text_audit = _run_json(
        [
            sys.executable,
            str(pdf_auditor),
            str(figure_paths["pdf"]),
            "--min-pt",
            "5",
            "--json",
        ]
    )
    png_audit = raster_pixel_audit(figure_paths["png"], RASTER_DPI)
    tiff_audit = raster_pixel_audit(figure_paths["tiff"], RASTER_DPI)
    panel_audit = [
        {
            "panel": "a",
            "unique_claim": "same-fold/same-seed primary effects versus B",
            "center": "fold-primary paired mean",
            "interval": "95% hierarchical percentile bootstrap CI",
            "replicate_unit": "5 genomic folds; 3 seeds nested within fold",
            "boundary": "1-bp-head-only CI includes zero; no equivalence claim",
            "visual_review": "passed" if manual_visual_pass else "pending",
        },
        {
            "panel": "b",
            "unique_claim": "registered 2 x 2 component main and interaction effects",
            "center": "fold-primary factorial effect",
            "interval": "95% hierarchical percentile bootstrap CI",
            "replicate_unit": "5 genomic folds; 3 seeds nested within fold",
            "boundary": "interaction is statistical non-additivity, not mechanism",
            "visual_review": "passed" if manual_visual_pass else "pending",
        },
        {
            "panel": "c",
            "unique_claim": "heterogeneity across all 241 registered tracks",
            "center": "track-wise fold-primary paired effect; median and IQR shown",
            "interval": "none; descriptive distribution",
            "replicate_unit": "tracks are outcomes, not independent replicates",
            "boundary": "no track-level inference",
            "visual_review": "passed" if manual_visual_pass else "pending",
        },
    ]
    static_ready = static_preflight.get("summary", {}).get("ready") is True
    pdf_ready = (
        pdf_text_audit.get("auditable") is True
        and int(pdf_text_audit.get("below_minimum_count", -1)) == 0
    )
    automated_pass = static_ready and pdf_ready and png_audit["passed"] and tiff_audit["passed"]
    qa_status = "passed" if automated_pass and manual_visual_pass else "pending_visual_review"
    artifacts = {key: artifact_record(path) for key, path in figure_paths.items()}
    artifacts["source_data"] = artifact_record(source_path)
    qa = {
        "schema_version": 1,
        "status": qa_status,
        "created_at": utc_now(),
        "backend": "python_matplotlib",
        "backend_exclusive": True,
        "cpu_only": True,
        "final_size_mm": [FIG_WIDTH_MM, FIG_HEIGHT_MM],
        "raster_dpi": RASTER_DPI,
        "minimum_pdf_glyph_pt": 5,
        "source_data_rows": len(source_rows),
        "source_data_sha256": sha256(source_path),
        "source_counts": {
            "panel_a_seed_pairs": 75,
            "panel_a_fold_means": 25,
            "panel_a_summaries": 5,
            "panel_b_seed_effects": 45,
            "panel_b_fold_means": 15,
            "panel_b_summaries": 3,
            "panel_c_tracks": 1205,
            "unique_tracks": data["track_count"],
        },
        "static_preflight": static_preflight,
        "pdf_text_audit": pdf_text_audit,
        "png_pixel_audit": png_audit,
        "tiff_pixel_audit": tiff_audit,
        "panel_audit": panel_audit,
        "manual_visual_review": {
            "status": "passed" if manual_visual_pass else "pending",
            "note": visual_review_note,
        },
        "artifacts": artifacts,
    }
    qa_path = output_dir / "figure_qa.json"
    write_json(qa_path, qa)

    runtime = {
        "python": sys.version.split()[0],
        "numpy": np.__version__,
        "matplotlib": mpl.__version__,
        "pillow": Image.__version__,
    }
    manifest = {
        "schema_version": 1,
        "status": qa_status,
        "created_at": utc_now(),
        "git_commit": git_commit(),
        "hostname": socket.gethostname(),
        "runtime": runtime,
        "command": [
            _display_path(Path(__file__)),
            "--output-dir",
            _display_path(output_dir),
            *( ["--manual-visual-pass"] if manual_visual_pass else [] ),
        ],
        "figure_contract": {
            "core_conclusion": "Full B exceeds matched scratch and major component removals; LoRA is the dominant component effect, while the 1-bp-head-only primary CI includes zero.",
            "archetype": "quantitative_grid",
            "inference_unit": "5 genomic folds with 3 seeds nested within fold",
            "track_boundary": "241 tracks are descriptive outcomes, not independent n",
            "interaction_boundary": "statistical non-additivity only",
            "final_test_access": "none",
        },
        "inputs": {
            key: artifact_record(path) for key, path in inputs["paths"].items()
        },
        "implementation": artifact_record(Path(__file__)),
        "qa_tools": {
            "validate_figure_sha256": sha256(validator),
            "audit_pdf_text_sha256": sha256(pdf_auditor),
        },
        "artifacts": {**artifacts, "qa": artifact_record(qa_path)},
        "source_data_schema": list(SOURCE_FIELDS),
        "source_data_rows": len(source_rows),
        "data_selection": {
            "rule": "prespecified primary_biological_score for the five requested candidates versus full B; all folds, seeds, and 241 tracks retained",
            "paired_summary_rows": {"input": 66, "selected": 5},
            "paired_raw_rows": {"input": 990, "selected": 75},
            "factorial_summary_rows": {"input": 33, "selected": 3},
            "factorial_raw_rows": {"input": 495, "selected": 45},
            "per_track_rows": {
                "input": 65070,
                "primary_rows_for_candidates_and_reference": 21690,
                "plotted_track_effects": 1205,
            },
        },
        "no_exclusions_within_requested_contract": True,
        "no_model_split_bigwig_or_final_test_reads": True,
    }
    manifest_path = output_dir / "manifest.json"
    write_json(manifest_path, manifest)
    return {"qa": qa, "manifest": manifest, "output_dir": output_dir}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--review", type=Path, default=DEFAULT_REVIEW)
    parser.add_argument("--paired", type=Path, default=DEFAULT_PAIRED)
    parser.add_argument("--raw-paired", type=Path, default=DEFAULT_RAW_PAIRED)
    parser.add_argument("--factorial", type=Path, default=DEFAULT_FACTORIAL)
    parser.add_argument("--raw-factorial", type=Path, default=DEFAULT_RAW_FACTORIAL)
    parser.add_argument("--per-track", type=Path, default=DEFAULT_PER_TRACK)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    default_skill_dir = Path(
        os.environ.get(
            "NATURE_FIGURE_SKILL_DIR",
            str(Path.home() / ".codex/skills/nature-figure"),
        )
    )
    parser.add_argument("--nature-figure-skill-dir", type=Path, default=default_skill_dir)
    parser.add_argument("--manual-visual-pass", action="store_true")
    parser.add_argument("--visual-review-note", default="")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    inputs = load_figure_inputs(
        review_path=args.review,
        paired_path=args.paired,
        raw_paired_path=args.raw_paired,
        factorial_path=args.factorial,
        raw_factorial_path=args.raw_factorial,
        per_track_path=args.per_track,
    )
    result = render_bundle(
        inputs,
        output_dir=args.output_dir,
        skill_dir=args.nature_figure_skill_dir,
        manual_visual_pass=args.manual_visual_pass,
        visual_review_note=args.visual_review_note,
    )
    print(result["output_dir"])
    print(result["qa"]["status"])


if __name__ == "__main__":
    main()
