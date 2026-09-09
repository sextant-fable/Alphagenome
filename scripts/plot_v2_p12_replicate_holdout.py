#!/usr/bin/env python3
"""Render the P12 biological-replicate holdout figures from frozen Source Data.

This script never trains a model and never reads the locked final-test block.
The optional coverage-example panel is rendered by a separate read-only
extraction script after its Source Data contract has been satisfied.
"""

from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys
from typing import Any, Iterable

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import numpy as np
from PIL import Image


REPO_ROOT = Path(__file__).resolve().parents[1]
METADATA = REPO_ROOT / "alphagenome_custom/metadata/v2"
P12_ROOT = METADATA / "replicate_holdout_v1"
DEFAULT_OUTPUT = REPO_ROOT / "results/v2_p12_replicate_holdout_figure"
SKILL_ROOT = Path("/home/zelinli6/.codex/skills/nature-figure")
WIDTH_MM = 183
RASTER_DPI = 600
EXPECTED_FOLDS = (1, 2, 3, 4, 5)
EXPECTED_SEEDS = (20260714, 20260715, 20260716)

COLORS = {
    "b": "#B23A48",
    "a": "#666666",
    "c": "#D08A3C",
    "c_size": "#177E89",
    "primary": "#B23A48",
    "supplementary": "#A9A9A9",
    "positive": "#177E89",
    "negative": "#B23A48",
    "ink": "#2F2F2F",
    "grid": "#D8D8D8",
    "soft": "#F2F2F2",
}

mpl.rcParams.update(
    {
        "font.family": "Liberation Sans",
        "font.sans-serif": ["Liberation Sans", "Arial", "DejaVu Sans"],
        "font.size": 7,
        "axes.labelsize": 7,
        "axes.titlesize": 8,
        "xtick.labelsize": 6,
        "ytick.labelsize": 6,
        "legend.fontsize": 6,
        "axes.linewidth": 0.7,
        "xtick.major.width": 0.7,
        "ytick.major.width": 0.7,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
    }
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(8 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def read_json(path: Path) -> dict[str, Any]:
    payload = json.loads(path.read_text())
    if not isinstance(payload, dict):
        raise ValueError(f"Expected JSON object: {path}")
    return payload


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise ValueError(f"Empty TSV: {path}")
    return rows


def write_tsv(path: Path, rows: Iterable[dict[str, Any]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def finite(value: object) -> float:
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f"Non-finite plotted value: {value}")
    return result


def setup_axis(axis: plt.Axes) -> None:
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    axis.grid(axis="y", color=COLORS["grid"], linewidth=0.45, zorder=0)
    axis.set_axisbelow(True)


def panel_label(axis: plt.Axes, label: str) -> None:
    axis.text(-0.15, 1.08, label, transform=axis.transAxes, fontsize=9, fontweight="bold", va="top")


def export_figure(fig: plt.Figure, output_dir: Path, stem: str, width_mm: float, height_mm: float) -> dict[str, str]:
    output_dir.mkdir(parents=True, exist_ok=True)
    fig.set_size_inches(width_mm / 25.4, height_mm / 25.4)
    paths = {
        "pdf": output_dir / f"{stem}.pdf",
        "svg": output_dir / f"{stem}.svg",
        "png": output_dir / f"{stem}.png",
        "tiff": output_dir / f"{stem}.tiff",
    }
    fig.savefig(paths["pdf"], bbox_inches="tight", facecolor="white")
    fig.savefig(paths["svg"], bbox_inches="tight", facecolor="white")
    fig.savefig(paths["png"], dpi=RASTER_DPI, bbox_inches="tight", facecolor="white")
    fig.savefig(paths["tiff"], dpi=RASTER_DPI, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return {key: str(path.relative_to(REPO_ROOT)) for key, path in paths.items()}


def validate_inputs() -> dict[str, Any]:
    review_path = METADATA / "audits/P12/review.json"
    review = read_json(review_path)
    if review.get("status") != "PASS" or review.get("phase") != "P12":
        raise RuntimeError("P12 figure generation requires a passing R12 review")
    if any(item.get("status") != "PASS" for item in review.get("checks", [])):
        raise RuntimeError("P12 review contains a non-passing check")
    execution = read_json(P12_ROOT / "p12_replicate_holdout_execution.json")
    if execution.get("status") != "completed" or execution.get("registered_tasks_completed") != 165:
        raise RuntimeError("P12 execution is incomplete")
    if execution.get("failures"):
        raise RuntimeError("P12 execution contains failures")
    if execution.get("locked_test_block_signal_reads") != 0:
        raise RuntimeError("P12 figure generation detected locked-test reads")
    summary = read_tsv(P12_ROOT / "p12_replicate_holdout_results.tsv")
    paired = read_tsv(P12_ROOT / "p12_replicate_holdout_paired_effects.tsv")
    metrics = read_tsv(P12_ROOT / "p12_replicate_holdout_metrics.tsv")
    roles = read_tsv(P12_ROOT / "p12_heldout_role_metrics.tsv")
    if len(summary) != 13 or len(paired) != 242 or len(metrics) != 3630 or len(roles) != 79200:
        raise RuntimeError("P12 figure input row counts do not match the R12 contract")
    if {int(row["fold"]) for row in metrics} != set(EXPECTED_FOLDS):
        raise RuntimeError("P12 metrics do not contain all five folds")
    return {
        "review_path": review_path,
        "review": review,
        "execution_path": P12_ROOT / "p12_replicate_holdout_execution.json",
        "execution": execution,
        "summary": summary,
        "paired": paired,
        "metrics": metrics,
        "roles": roles,
        "input_paths": [
            review_path,
            P12_ROOT / "p12_replicate_holdout_execution.json",
            P12_ROOT / "p12_replicate_holdout_results.tsv",
            P12_ROOT / "p12_replicate_holdout_paired_effects.tsv",
            P12_ROOT / "p12_replicate_holdout_metrics.tsv",
            P12_ROOT / "p12_heldout_role_metrics.tsv",
            P12_ROOT / "heldout_primary_group_manifest.tsv",
            P12_ROOT / "heldout_supplementary_group_manifest.tsv",
            P12_ROOT / "holdout_assignments.tsv",
        ],
    }


def make_design_figure(data: dict[str, Any], output: Path) -> dict[str, str]:
    assignments = read_tsv(P12_ROOT / "holdout_assignments.tsv")
    counts = {}
    for row in assignments:
        counts[row["assignment_role"]] = counts.get(row["assignment_role"], 0) + 1
    fig, axes = plt.subplots(1, 3, figsize=(WIDTH_MM / 25.4, 102 / 25.4), gridspec_kw={"wspace": 0.38})

    ax = axes[0]
    ax.axis("off")
    panel_label(ax, "a")
    ax.text(0.02, 0.95, "Audited source inventory", transform=ax.transAxes, fontsize=8, fontweight="bold", va="top")
    y = 0.78
    inventory = [("482", "verified RNA-seq runs", COLORS["ink"]), ("241", "grouped RNA tracks", COLORS["b"]), ("27", "source studies", COLORS["c_size"])]
    for number, label, color in inventory:
        ax.text(0.08, y, number, transform=ax.transAxes, fontsize=17, fontweight="bold", color=color, va="center")
        ax.text(0.38, y, label, transform=ax.transAxes, fontsize=7.5, va="center")
        y -= 0.18
    ax.text(0.08, 0.20, "WBcel235 reference genome", transform=ax.transAxes, fontsize=7, color=COLORS["ink"])
    ax.text(0.08, 0.12, "241 output heads retained in every model", transform=ax.transAxes, fontsize=7, color=COLORS["ink"])

    ax = axes[1]
    ax.axis("off")
    panel_label(ax, "b")
    ax.text(0.02, 0.95, "Biological-unit holdout", transform=ax.transAxes, fontsize=8, fontweight="bold", va="top")
    boxes = [
        (0.05, 0.70, 0.90, 0.16, "Member ledger\n482 runs", COLORS["soft"]),
        (0.04, 0.40, 0.42, 0.16, "Training\nunit excluded", "#E5F0F1"),
        (0.54, 0.40, 0.42, 0.16, "Held-out\nevaluation only", "#F5E4E5"),
    ]
    for x, y0, w, h, text, face in boxes:
        patch = FancyBboxPatch((x, y0), w, h, boxstyle="round,pad=0.012", transform=ax.transAxes, facecolor=face, edgecolor=COLORS["ink"], linewidth=0.7)
        ax.add_patch(patch)
        ax.text(x + w / 2, y0 + h / 2, text, transform=ax.transAxes, ha="center", va="center", fontsize=6.4)
    for x1, y1, x2, y2 in [(0.50, 0.70, 0.25, 0.56), (0.50, 0.70, 0.75, 0.56)]:
        ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), transform=ax.transAxes, arrowstyle="-|>", mutation_scale=8, linewidth=0.7, color=COLORS["ink"]))
    ax.text(0.05, 0.22, "Duplicate-source equivalence classes stay on one side", transform=ax.transAxes, fontsize=6.5)
    ax.text(0.05, 0.13, f"Primary candidates: {counts.get('primary_candidate', 0)} | pending QC: {counts.get('primary_pending_qc', 0)}", transform=ax.transAxes, fontsize=6.5)
    ax.text(0.05, 0.06, f"Supplementary candidates: {counts.get('supplementary_candidate', 0)} | training-only: {counts.get('training_only', 0)}", transform=ax.transAxes, fontsize=6.5)

    ax = axes[2]
    ax.axis("off")
    panel_label(ax, "c")
    ax.text(0.02, 0.95, "Frozen evaluation contract", transform=ax.transAxes, fontsize=8, fontweight="bold", va="top")
    y = 0.78
    lines = [
        ("Input", "131,072-bp DNA window"),
        ("Output", "241 RNA tracks at 1-bp and 128-bp"),
        ("Blocks", "5 chromosome folds (I-V)"),
        ("Repeats", "3 seeds nested within fold"),
        ("Scopes", "internal labels and held-out labels"),
        ("Lock", "0 reads from locked final-test block"),
    ]
    for key, value in lines:
        ax.text(0.03, y, key, transform=ax.transAxes, fontsize=7, fontweight="bold", color=COLORS["b"])
        ax.text(0.40, y, value, transform=ax.transAxes, fontsize=7)
        y -= 0.105
    ax.text(0.03, 0.10, "All 165 P12 jobs used new training-label aggregates.", transform=ax.transAxes, fontsize=6.5, color=COLORS["ink"])
    return export_figure(fig, output, "figure_p12_holdout_design", WIDTH_MM, 102)


def summary_row(rows: list[dict[str, str]], configuration: str) -> dict[str, str]:
    matches = [row for row in rows if row["configuration"] == configuration]
    if len(matches) != 1:
        raise ValueError(f"Expected one summary row for {configuration}")
    return matches[0]


def heldout_effect(rows: list[dict[str, str]], comparison: str, metric: str, *, invert: bool = False) -> tuple[float, float, float]:
    matches = [row for row in rows if row["scope"] == "heldout" and row["comparison"] == comparison and row["metric"] == metric]
    if len(matches) != 1:
        raise ValueError(f"Expected one heldout effect for {comparison}/{metric}")
    row = matches[0]
    values = (finite(row["estimate"]), finite(row["ci_95_low"]), finite(row["ci_95_high"]))
    if invert:
        return (-values[0], -values[2], -values[1])
    return values


def fold_means(metrics: list[dict[str, str]], configuration: str, scope: str, metric: str) -> list[tuple[int, float]]:
    selected = [row for row in metrics if row["configuration"] == configuration and row["scope"] == scope and row["metric"] == metric]
    grouped: dict[int, list[float]] = {}
    for row in selected:
        grouped.setdefault(int(row["fold"]), []).append(finite(row["value"]))
    if set(grouped) != set(EXPECTED_FOLDS):
        raise ValueError(f"Missing folds for {configuration}/{scope}/{metric}")
    return [(fold, float(np.mean(grouped[fold]))) for fold in EXPECTED_FOLDS]


def make_performance_figure(data: dict[str, Any], output: Path) -> dict[str, str]:
    summary = data["summary"]
    paired = data["paired"]
    metrics = data["metrics"]
    roles = data["roles"]
    metric_keys = ["primary_biological_score", "gene_exon_coverage_pearson_log1p", "per_track_pearson_128bp_log1p"]
    metric_labels = ["Primary score", "Gene-exon Pearson", "128-bp Pearson"]
    configs = [("A_paper", "Model A", COLORS["a"]), ("B_paper", "Model B", COLORS["b"]), ("C_paper", "Model C", COLORS["c"]), ("C_size_matched", "Size-matched C", COLORS["c_size"])]
    fig, axes = plt.subplots(2, 2, figsize=(WIDTH_MM / 25.4, 124 / 25.4), gridspec_kw={"wspace": 0.35, "hspace": 0.42})

    ax = axes[0, 0]
    setup_axis(ax)
    x = np.arange(len(metric_keys))
    offsets = np.linspace(-0.27, 0.27, len(configs))
    for idx, (configuration, label, color) in enumerate(configs):
        means = []
        lows = []
        highs = []
        for metric in metric_keys:
            row = summary_row(summary, configuration)
            means.append(finite(row[f"mean_heldout_{metric}"]))
            # The configuration summary is a point estimate. Fold-level points
            # expose the uncertainty without inventing a configuration CI.
            values = np.array([value for _, value in fold_means(metrics, configuration, "heldout", metric)])
            lows.append(float(values.min()))
            highs.append(float(values.max()))
        for xi, value, low, high in zip(x + offsets[idx], means, lows, highs):
            ax.plot([xi, xi], [low, high], color=color, alpha=0.25, linewidth=0.8, zorder=1)
            ax.scatter([xi] * len(EXPECTED_FOLDS), [value for _, value in fold_means(metrics, configuration, "heldout", metric)], s=8, color=color, alpha=0.55, edgecolors="white", linewidths=0.25, zorder=2)
            ax.scatter([xi], [value], s=28 if configuration == "B_paper" else 20, color=color, marker="D", edgecolors="white", linewidths=0.4, zorder=3, label=label if xi == x[0] + offsets[idx] else None)
    ax.set_xticks(x, metric_labels)
    ax.set_ylabel("Held-out score")
    ax.set_title("Model ordering persists on unseen biological units", loc="left", fontweight="bold")
    ax.legend(frameon=False, ncol=2, loc="lower right", handletextpad=0.4, columnspacing=0.8)
    panel_label(ax, "a")

    ax = axes[0, 1]
    setup_axis(ax)
    comparisons = [("B_vs_A_paper", "B - A", False), ("B_vs_C_paper", "B - C", False), ("C_size_matched_vs_B", "B - C matched", True)]
    y_positions = np.arange(len(comparisons) * len(metric_keys))[::-1]
    yticks = []
    ypos = 0
    for comparison, label, invert in comparisons:
        for metric, metric_label in zip(metric_keys, metric_labels):
            estimate, low, high = heldout_effect(paired, comparison, metric, invert=invert)
            y = y_positions[ypos]
            ax.errorbar(estimate, y, xerr=[[estimate - low], [high - estimate]], fmt="o", color=COLORS["positive"] if estimate >= 0 else COLORS["negative"], markersize=4.5, capsize=2.2, linewidth=0.9)
            yticks.append(f"{label}: {metric_label}")
            ypos += 1
    ax.axvline(0, color=COLORS["ink"], linewidth=0.7)
    ax.set_yticks(y_positions, yticks, fontsize=5.7)
    ax.set_xlabel("Paired effect (positive favors Model B)")
    ax.set_title("Matched held-out effects", loc="left", fontweight="bold")
    panel_label(ax, "b")

    ax = axes[1, 0]
    setup_axis(ax)
    role_values = {}
    for role in ("primary", "supplementary"):
        selected = [row for row in roles if row["configuration"] == "B_paper" and row["metric"] == "primary_biological_score" and row["role"] == role]
        by_track: dict[str, list[float]] = {}
        for row in selected:
            by_track.setdefault(row["track_id"], []).append(finite(row["value"]))
        role_values[role] = np.array([np.mean(values) for values in by_track.values()])
    bp = ax.boxplot([role_values["primary"], role_values["supplementary"]], positions=[0, 1], widths=0.42, patch_artist=True, showfliers=False, medianprops={"color": COLORS["ink"], "linewidth": 1.0}, whiskerprops={"color": COLORS["ink"], "linewidth": 0.7}, capprops={"color": COLORS["ink"], "linewidth": 0.7})
    for patch, color in zip(bp["boxes"], [COLORS["primary"], COLORS["supplementary"]]):
        patch.set_facecolor(color)
        patch.set_alpha(0.45)
    rng = np.random.default_rng(20260822)
    for idx, role in enumerate(("primary", "supplementary")):
        values = role_values[role]
        ax.scatter(np.full(values.size, idx) + rng.uniform(-0.13, 0.13, values.size), values, s=7, color=COLORS[role], alpha=0.65, edgecolors="none", rasterized=True)
    ax.set_xticks([0, 1], [f"Primary\n(n={role_values['primary'].size})", f"Supplementary\n(n={role_values['supplementary'].size})"])
    ax.set_ylabel("Held-out primary score")
    ax.set_title("Track-level heterogeneity", loc="left", fontweight="bold")
    ax.text(0.02, 0.03, "Descriptive tracks; not independent replicates", transform=ax.transAxes, fontsize=5.6, color=COLORS["ink"])
    panel_label(ax, "c")

    ax = axes[1, 1]
    setup_axis(ax)
    for metric, label, color in zip(metric_keys, metric_labels, [COLORS["b"], COLORS["c_size"], COLORS["positive"]]):
        internal = np.array([value for _, value in fold_means(metrics, "B_paper", "internal", metric)])
        heldout = np.array([value for _, value in fold_means(metrics, "B_paper", "heldout", metric)])
        ax.plot(np.arange(5), internal, color=color, alpha=0.45, linewidth=0.8)
        ax.scatter(np.arange(5), internal, color=color, marker="o", s=12, alpha=0.55)
        ax.plot(np.arange(5), heldout, color=color, linewidth=1.2, marker="D", markersize=4, label=label)
    ax.set_xticks(np.arange(5), ["I", "II", "III", "IV", "V"])
    ax.set_xlabel("Validation chromosome fold")
    ax.set_ylabel("Model B score")
    ax.set_title("Performance is retained after replicate holdout", loc="left", fontweight="bold")
    ax.legend(frameon=False, fontsize=5.5, loc="upper right", bbox_to_anchor=(0.98, 0.98), ncol=1)
    panel_label(ax, "d")
    return export_figure(fig, output, "figure_p12_heldout_performance", WIDTH_MM, 124)


def make_component_figure(data: dict[str, Any], output: Path) -> dict[str, str]:
    paired = data["paired"]
    roles = data["roles"]
    metrics = ["primary_biological_score", "gene_exon_coverage_pearson_log1p", "per_track_pearson_128bp_log1p"]
    metric_labels = ["Primary score", "Gene-exon Pearson", "128-bp Pearson"]
    components = [("B_no_lora_vs_B", "no LoRA"), ("B_lora_only_vs_B", "LoRA only"), ("B_1bp_head_only_vs_B", "1-bp head"), ("B_no_augmentation_vs_B", "no aug."), ("B_no_gene_loss_vs_B", "no gene loss"), ("C_size_matched_vs_B", "size-match C")]
    fig, axes = plt.subplots(1, 3, figsize=(WIDTH_MM / 25.4, 106 / 25.4), gridspec_kw={"wspace": 0.50})

    ax = axes[0]
    setup_axis(ax)
    y = np.arange(len(components))[::-1]
    for yi, (comparison, label) in zip(y, components):
        estimate, low, high = heldout_effect(paired, comparison, "primary_biological_score", invert=True)
        color = COLORS["positive"] if low > 0 else COLORS["negative"] if high < 0 else COLORS["c_size"]
        ax.errorbar(estimate, yi, xerr=[[estimate - low], [high - estimate]], fmt="o", color=color, markersize=4.5, capsize=2.3, linewidth=0.9)
    ax.axvline(0, color=COLORS["ink"], linewidth=0.7)
    ax.set_yticks(y, [label for _, label in components], fontsize=6)
    ax.set_xlabel("B advantage on held-out primary score")
    ax.set_title("Measured component effects", loc="left", fontweight="bold")
    panel_label(ax, "a")

    ax = axes[1]
    setup_axis(ax)
    matrix = []
    for comparison, _ in components:
        matrix.append([heldout_effect(paired, comparison, metric, invert=True)[0] for metric in metrics])
    matrix_array = np.asarray(matrix)
    vmax = max(abs(matrix_array.min()), abs(matrix_array.max()))
    image = ax.imshow(matrix_array, cmap="PuOr_r", vmin=-vmax, vmax=vmax, aspect="auto")
    ax.set_xticks(np.arange(3), metric_labels, rotation=25, ha="right")
    ax.set_yticks(np.arange(len(components)), [label for _, label in components], fontsize=5.6)
    for i in range(matrix_array.shape[0]):
        for j in range(matrix_array.shape[1]):
            ax.text(j, i, f"{matrix_array[i, j]:+.3f}", ha="center", va="center", fontsize=5.5, color=COLORS["ink"])
    ax.set_title("Effect consistency across metrics", loc="left", fontweight="bold")
    panel_label(ax, "b")
    cbar = fig.colorbar(image, ax=ax, fraction=0.045, pad=0.03)
    cbar.ax.tick_params(labelsize=5.5, width=0.5)
    cbar.set_label("B advantage", fontsize=6)

    ax = axes[2]
    setup_axis(ax)
    primary_track_values: dict[str, np.ndarray] = {}
    for comparison, label in components:
        candidate = comparison.split("_vs_")[0]
        b_rows = [row for row in roles if row["configuration"] == "B_paper" and row["role"] == "primary" and row["metric"] == "primary_biological_score"]
        c_rows = [row for row in roles if row["configuration"] == candidate and row["role"] == "primary" and row["metric"] == "primary_biological_score"]
        b_by = {}
        c_by = {}
        for row in b_rows:
            b_by.setdefault(row["track_id"], []).append(finite(row["value"]))
        for row in c_rows:
            c_by.setdefault(row["track_id"], []).append(finite(row["value"]))
        common = sorted(set(b_by) & set(c_by))
        primary_track_values[label] = np.array([np.mean(b_by[k]) - np.mean(c_by[k]) for k in common])
    positions = np.arange(len(components))
    bp = ax.boxplot([primary_track_values[label] for _, label in components], positions=positions, widths=0.45, patch_artist=True, showfliers=False, medianprops={"color": COLORS["ink"], "linewidth": 0.9}, whiskerprops={"color": COLORS["ink"], "linewidth": 0.6}, capprops={"color": COLORS["ink"], "linewidth": 0.6})
    for patch in bp["boxes"]:
        patch.set_facecolor(COLORS["soft"])
        patch.set_edgecolor(COLORS["ink"])
    rng = np.random.default_rng(20260822)
    for pos, (_, label) in enumerate(components):
        values = primary_track_values[label]
        ax.scatter(np.full(values.size, pos) + rng.uniform(-0.14, 0.14, values.size), values, s=5, alpha=0.42, color=COLORS["b"], edgecolors="none", rasterized=True)
    ax.axhline(0, color=COLORS["ink"], linewidth=0.7)
    ax.set_xticks(positions, [label.replace(" ", "\n", 1) for _, label in components], rotation=35, ha="right", fontsize=5.1)
    ax.set_ylabel("Per-track B advantage")
    ax.set_title("Descriptive track heterogeneity", loc="left", fontweight="bold")
    ax.text(0.02, 0.03, "n = 57 primary tracks; no track-level inference", transform=ax.transAxes, fontsize=5.6, color=COLORS["ink"])
    panel_label(ax, "c")
    return export_figure(fig, output, "figure_p12_component_attribution", WIDTH_MM, 106)


def write_source_data(data: dict[str, Any], output: Path) -> dict[str, str]:
    output.mkdir(parents=True, exist_ok=True)
    performance_fields = ["panel", "configuration", "comparison", "scope", "metric", "fold", "seed", "role", "track_id", "value", "estimate", "ci_95_low", "ci_95_high", "notes"]
    performance_rows: list[dict[str, Any]] = []
    for row in data["summary"]:
        for scope in ("internal", "heldout"):
            for metric in ("primary_biological_score", "gene_exon_coverage_pearson_log1p", "per_track_pearson_128bp_log1p"):
                performance_rows.append({"panel": "2a_summary", "configuration": row["configuration"], "scope": scope, "metric": metric, "value": row[f"mean_{scope}_{metric}"], "notes": "15-job mean; fold-level values are separately provided in panel 2d source rows."})
    for row in data["paired"]:
        if row["scope"] == "heldout" and row["metric"] in {"primary_biological_score", "gene_exon_coverage_pearson_log1p", "per_track_pearson_128bp_log1p"}:
            performance_rows.append({"panel": "2b_paired_effect", "comparison": row["comparison"], "scope": row["scope"], "metric": row["metric"], "estimate": row["estimate"], "ci_95_low": row["ci_95_low"], "ci_95_high": row["ci_95_high"], "notes": row["ci_method"]})
    for row in data["roles"]:
        if row["configuration"] == "B_paper" and row["metric"] == "primary_biological_score":
            performance_rows.append({"panel": "2c_track_distribution", "configuration": row["configuration"], "role": row["role"], "fold": row["fold"], "seed": row["seed"], "track_id": row["track_id"], "value": row["value"], "notes": "Track-level descriptive value; not an inferential replicate."})
    for metric in ("primary_biological_score", "gene_exon_coverage_pearson_log1p", "per_track_pearson_128bp_log1p"):
        for scope in ("internal", "heldout"):
            for fold, value in fold_means(data["metrics"], "B_paper", scope, metric):
                performance_rows.append({"panel": "2d_scope_fold", "configuration": "B_paper", "scope": scope, "metric": metric, "fold": fold, "value": value, "notes": "Mean across three seeds within fold."})
    perf_path = output / "figure_source_data_heldout_performance.tsv"
    write_tsv(perf_path, performance_rows, performance_fields)

    component_fields = ["panel", "comparison", "configuration", "metric", "estimate", "ci_95_low", "ci_95_high", "role", "fold", "seed", "track_id", "value", "notes"]
    component_rows: list[dict[str, Any]] = []
    component_comparisons = {"B_no_lora_vs_B", "B_lora_only_vs_B", "B_1bp_head_only_vs_B", "B_no_augmentation_vs_B", "B_no_gene_loss_vs_B", "C_size_matched_vs_B"}
    for row in data["paired"]:
        if row["scope"] == "heldout" and row["comparison"] in component_comparisons and row["metric"] in {"primary_biological_score", "gene_exon_coverage_pearson_log1p", "per_track_pearson_128bp_log1p"}:
            component_rows.append({"panel": "3a_3b_paired_effect", "comparison": row["comparison"], "metric": row["metric"], "estimate": row["estimate"], "ci_95_low": row["ci_95_low"], "ci_95_high": row["ci_95_high"], "notes": "Candidate-minus-B source effect; plotting reverses sign to display B advantage."})
    for row in data["roles"]:
        if row["role"] == "primary" and row["metric"] == "primary_biological_score" and (row["configuration"] == "B_paper" or row["configuration"] in {comparison.split("_vs_")[0] for comparison in component_comparisons}):
            component_rows.append({"panel": "3c_track_distribution", "configuration": row["configuration"], "fold": row["fold"], "seed": row["seed"], "track_id": row["track_id"], "value": row["value"], "notes": "Track-level descriptive value; 57 primary tracks."})
    comp_path = output / "figure_source_data_component_attribution.tsv"
    write_tsv(comp_path, component_rows, component_fields)
    return {"performance": str(perf_path.relative_to(REPO_ROOT)), "component": str(comp_path.relative_to(REPO_ROOT))}


def run_qa(output: Path, source_script: Path) -> dict[str, Any]:
    validator = SKILL_ROOT / "scripts/validate_figure.py"
    pdf_auditor = SKILL_ROOT / "scripts/audit_pdf_text.py"
    qa: dict[str, Any] = {"status": "passed", "created_at": utc_now(), "backend": "python_matplotlib", "figures": {}}
    for pdf in sorted(output.glob("*.pdf")):
        stem = pdf.stem
        validation = subprocess.run([sys.executable, str(validator), str(source_script), "--backend", "python", "--json"], capture_output=True, text=True, check=False)
        pdf_audit = subprocess.run([sys.executable, str(pdf_auditor), str(pdf), "--min-pt", "5", "--json"], capture_output=True, text=True, check=False)
        try:
            validation_json = json.loads(validation.stdout)
        except json.JSONDecodeError:
            validation_json = {"status": "invalid", "stdout": validation.stdout, "stderr": validation.stderr}
        try:
            pdf_json = json.loads(pdf_audit.stdout)
        except json.JSONDecodeError:
            pdf_json = {"status": "invalid", "stdout": pdf_audit.stdout, "stderr": pdf_audit.stderr}
        png = output / f"{stem}.png"
        tiff = output / f"{stem}.tiff"
        raster = {}
        for path in (png, tiff):
            with Image.open(path) as image:
                raster[path.suffix.lstrip(".")] = {"size": list(image.size), "dpi": [float(v) if v is not None else None for v in image.info.get("dpi", (None, None))], "mode": image.mode}
        if validation.returncode != 0 or pdf_audit.returncode != 0:
            qa["status"] = "failed"
        qa["figures"][stem] = {"pdf": str(pdf.relative_to(REPO_ROOT)), "pdf_sha256": sha256(pdf), "svg_sha256": sha256(output / f"{stem}.svg"), "png_sha256": sha256(png), "tiff_sha256": sha256(tiff), "source_validation": validation_json, "pdf_text_audit": pdf_json, "raster": raster}
    return qa


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", default=str(DEFAULT_OUTPUT))
    args = parser.parse_args()
    output = Path(args.output_dir).resolve()
    data = validate_inputs()
    output.mkdir(parents=True, exist_ok=True)
    design_paths = make_design_figure(data, output)
    performance_paths = make_performance_figure(data, output)
    component_paths = make_component_figure(data, output)
    source_paths = write_source_data(data, output)
    manifest = {
        "schema_version": 1,
        "figure_contract": "nature_methods_manuscript/plans/figure_contracts/fig_p12_replicate_holdout.md",
        "created_at": utc_now(),
        "backend": "python_matplotlib",
        "figure_width_mm": WIDTH_MM,
        "raster_dpi": RASTER_DPI,
        "figures": {"design": design_paths, "performance": performance_paths, "component": component_paths},
        "source_data": source_paths,
        "input_sha256": {str(path.relative_to(REPO_ROOT)): sha256(path) for path in data["input_paths"] if path.is_file()},
        "counts": {"p12_jobs": 165, "p12_metric_rows": 3630, "p12_role_metric_rows": 79200, "paired_effect_rows": 242, "heldout_primary_candidate_groups": 57, "heldout_pending_qc_groups": 1, "heldout_supplementary_candidate_groups": 22},
        "statistics": {"primary_unit": "five genomic folds", "nested_repeat": "three seeds within fold", "bootstrap": "fold-first, seed-within-fold hierarchical percentile bootstrap", "bootstrap_replicates": 10000, "track_role": "descriptive outcomes, not independent inferential replicates"},
        "coverage_examples": "pending read-only frozen-checkpoint extraction",
    }
    manifest_path = output / "manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    qa = run_qa(output, Path(__file__))
    qa["manifest_sha256"] = sha256(manifest_path)
    qa_path = output / "figure_qa.json"
    qa_path.write_text(json.dumps(qa, indent=2, sort_keys=True) + "\n")
    if qa["status"] != "passed":
        raise SystemExit("P12 figure QA failed")
    print(json.dumps({"status": "passed", "output": str(output.relative_to(REPO_ROOT)), "figures": manifest["figures"], "qa": str(qa_path.relative_to(REPO_ROOT))}, indent=2))


if __name__ == "__main__":
    main()
