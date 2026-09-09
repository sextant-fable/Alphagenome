#!/usr/bin/env python3
"""Render P13 Nature-style figures from audited, saved P11/P12 artifacts.

This module is intentionally render-only: it does not open BigWigs, FASTA,
GTF, checkpoints or the locked test. It turns completed P11/P12 metadata and
saved P12 metrics into four evidence-bound manuscript figures. Independent
study/laboratory and regulatory-variant figures are deliberately excluded until
their data contracts and result artifacts exist.
"""

from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
import math
from collections import Counter, defaultdict
from pathlib import Path
import subprocess
import sys
from typing import Any, Iterable

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Rectangle
import numpy as np
from PIL import Image


REPO_ROOT = Path(__file__).resolve().parents[1]
P12_ROOT = REPO_ROOT / "alphagenome_custom/metadata/v2/replicate_holdout_v1"
P12_REVIEW = REPO_ROOT / "alphagenome_custom/metadata/v2/audits/P12/review.json"
P12_SPEC = REPO_ROOT / "alphagenome_custom/metadata/v2/p12_replicate_holdout_spec.json"
P11_SPEC = P12_ROOT / "p11_replicate_holdout_spec.json"
P13_AUDIT = REPO_ROOT / "results/v2_p13_manuscript_evidence_audit"
P13_DESCRIPTIVES = REPO_ROOT / "results/v2_p13_heldout_descriptives"
P13_REPLICATE_REFERENCE = REPO_ROOT / "results/v2_p13_replicate_agreement_iv"
P13_REVIEW = REPO_ROOT / "alphagenome_custom/metadata/v2/audits/P13/review.json"
P14_MODEL_REFERENCE = REPO_ROOT / "results/v2_p14_iv_model_reference"
P14_REVIEW = REPO_ROOT / "alphagenome_custom/metadata/v2/audits/P14/review.json"
COVERAGE_ROOT = REPO_ROOT / "results/v2_p12_replicate_holdout_figure/coverage_examples"
SKILL_ROOT = Path("/home/zelinli6/.codex/skills/nature-figure")
DEFAULT_OUTPUT = REPO_ROOT / "results/v2_p13_manuscript_figures"
WIDTH_MM = 183
RASTER_DPI = 600
FOLDS = (1, 2, 3, 4, 5)

COLORS = {
    "ink": "#1F2933",
    "muted": "#52606D",
    "grid": "#E6E8EB",
    "paper": "#FFFFFF",
    "soft": "#F7F8FA",
    "train": "#D9EAF2",
    "heldout": "#F6E3C5",
    "alert": "#C44E52",
    "b": "#007C91",
    "a": "#6B7280",
    "c": "#B7791F",
    "c_size": "#374151",
    "positive": "#007C91",
    "neutral": "#9A6700",
    "stage": ["#007C91", "#4C78A8", "#F58518", "#54A24B", "#B279A2", "#9D9D9D"],
}

mpl.rcParams.update(
    {
        "font.family": "Liberation Sans",
        "font.sans-serif": ["Liberation Sans", "Arial", "DejaVu Sans"],
        "font.size": 7.2,
        "axes.labelsize": 7,
        "axes.titlesize": 8,
        "xtick.labelsize": 6,
        "ytick.labelsize": 6,
        "legend.fontsize": 6,
        "axes.linewidth": 0.65,
        "xtick.major.width": 0.7,
        "ytick.major.width": 0.7,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
        "figure.facecolor": "white",
        "axes.facecolor": "white",
    }
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    return parser.parse_args()


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


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
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def finite(value: str | float | int) -> float:
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f"Non-finite plotted value: {value}")
    return result


def panel_label(axis: plt.Axes, label: str) -> None:
    axis.text(-0.16, 1.08, label, transform=axis.transAxes, fontsize=9, fontweight="bold", va="top")


def setup_axis(axis: plt.Axes, *, grid: bool = True) -> None:
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    if grid:
        axis.grid(axis="y", color=COLORS["grid"], linewidth=0.45, zorder=0)
        axis.set_axisbelow(True)


def export_figure(fig: plt.Figure, output: Path, stem: str, height_mm: float) -> dict[str, str]:
    output.mkdir(parents=True, exist_ok=True)
    fig.set_size_inches(WIDTH_MM / 25.4, height_mm / 25.4)
    paths = {kind: output / f"{stem}.{extension}" for kind, extension in {"pdf": "pdf", "svg": "svg", "png": "png", "tiff": "tiff"}.items()}
    fig.savefig(paths["pdf"], bbox_inches="tight", facecolor="white")
    fig.savefig(paths["svg"], bbox_inches="tight", facecolor="white")
    fig.savefig(paths["png"], dpi=RASTER_DPI, bbox_inches="tight", facecolor="white")
    fig.savefig(paths["tiff"], dpi=RASTER_DPI, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return {kind: str(path.relative_to(REPO_ROOT)) for kind, path in paths.items()}


def validate_inputs() -> dict[str, Any]:
    required = {
        "p12_review": P12_REVIEW,
        "p12_spec": P12_SPEC,
        "p11_spec": P11_SPEC,
        "assignments": P12_ROOT / "holdout_assignments.tsv",
        "training_groups": P12_ROOT / "training_group_manifest.tsv",
        "metrics": P12_ROOT / "p12_replicate_holdout_metrics.tsv",
        "paired": P12_ROOT / "p12_replicate_holdout_paired_effects.tsv",
        "job_records": P12_ROOT / "p12_replicate_holdout_job_records.tsv",
        "configuration_summary": P13_AUDIT / "p13_configuration_summary.tsv",
        "fold_seed_effects": P13_AUDIT / "p13_fold_seed_effects.tsv",
        "loo": P13_AUDIT / "p13_leave_one_fold_out.tsv",
        "comparison_summary": P13_AUDIT / "p13_comparison_summary.tsv",
        "group_descriptives": P13_DESCRIPTIVES / "p13_b_paper_primary_group_descriptives.tsv",
        "strata": P13_DESCRIPTIVES / "p13_b_paper_primary_metadata_strata.tsv",
        "p13_replicate_by_fold": P13_REPLICATE_REFERENCE / "p13_replicate_agreement_by_fold.tsv",
        "p13_replicate_summary": P13_REPLICATE_REFERENCE / "p13_replicate_agreement_summary.tsv",
        "p13_review": P13_REVIEW,
        "p14_model_by_fold": P14_MODEL_REFERENCE / "p14_iv_model_reference_by_fold.tsv",
        "p14_model_summary": P14_MODEL_REFERENCE / "p14_iv_model_reference_summary.tsv",
        "p14_model_vs_replicate": P14_MODEL_REFERENCE / "p14_iv_model_vs_replicate_reference_by_fold.tsv",
        "p14_review": P14_REVIEW,
        "coverage_manifest": COVERAGE_ROOT / "coverage_examples_manifest.json",
        "coverage_source": COVERAGE_ROOT / "coverage_examples_source_data.tsv",
    }
    for path in required.values():
        if not path.is_file():
            raise FileNotFoundError(path)
    review = read_json(required["p12_review"])
    if review.get("phase") != "P12" or review.get("status") != "PASS":
        raise RuntimeError("P13 figures require a passing P12 review")
    if any(check.get("status") != "PASS" for check in review.get("checks", [])):
        raise RuntimeError("P12 review contains a non-passing check")
    for name in ("p13_review", "p14_review"):
        phase_review = read_json(required[name])
        if phase_review.get("status") != "PASS":
            raise RuntimeError(f"{name} must pass before rendering the matched I--V reference")
    spec = read_json(required["p12_spec"])
    p11_spec = read_json(required["p11_spec"])
    execution = read_json(P12_ROOT / "p12_replicate_holdout_execution.json")
    if execution.get("status") != "completed" or execution.get("registered_tasks_completed") != 165:
        raise RuntimeError("P12 execution is incomplete")
    if execution.get("locked_test_block_signal_reads") != 0:
        raise RuntimeError("P12 execution reports locked-test reads")
    coverage = read_tsv(required["coverage_source"])
    if any(row["chromosome"] == "X" for row in coverage):
        raise RuntimeError("P13 figures refuse coverage source rows on chromosome X")
    return {"paths": required, "review": review, "spec": spec, "p11_spec": p11_spec, "execution": execution}


def load_data(validated: dict[str, Any]) -> dict[str, Any]:
    paths = validated["paths"]
    return {
        **validated,
        "assignments": read_tsv(paths["assignments"]),
        "training_groups": read_tsv(paths["training_groups"]),
        "metrics": read_tsv(paths["metrics"]),
        "paired": read_tsv(paths["paired"]),
        "job_records": read_tsv(paths["job_records"]),
        "configuration_summary": read_tsv(paths["configuration_summary"]),
        "fold_seed_effects": read_tsv(paths["fold_seed_effects"]),
        "loo": read_tsv(paths["loo"]),
        "comparison_summary": read_tsv(paths["comparison_summary"]),
        "group_descriptives": read_tsv(paths["group_descriptives"]),
        "strata": read_tsv(paths["strata"]),
        "p13_replicate_by_fold": read_tsv(paths["p13_replicate_by_fold"]),
        "p13_replicate_summary": read_tsv(paths["p13_replicate_summary"]),
        "p14_model_by_fold": read_tsv(paths["p14_model_by_fold"]),
        "p14_model_summary": read_tsv(paths["p14_model_summary"]),
        "p14_model_vs_replicate": read_tsv(paths["p14_model_vs_replicate"]),
        "coverage_manifest": read_json(paths["coverage_manifest"]),
        "coverage_source": read_tsv(paths["coverage_source"]),
    }


def grouped_mean(rows: list[dict[str, str]], keys: tuple[str, ...], value: str) -> dict[tuple[str, ...], float]:
    grouped: dict[tuple[str, ...], list[float]] = defaultdict(list)
    for row in rows:
        grouped[tuple(row[key] for key in keys)].append(finite(row[value]))
    return {key: float(np.mean(values)) for key, values in grouped.items()}


def paired_effect(row: dict[str, str], *, invert: bool) -> tuple[float, float, float]:
    estimate, low, high = (finite(row[key]) for key in ("estimate", "ci_95_low", "ci_95_high"))
    return (-estimate, -high, -low) if invert else (estimate, low, high)


def parameter_counts(job_records: list[dict[str, str]]) -> dict[str, int]:
    result: dict[str, set[int]] = defaultdict(set)
    for row in job_records:
        result[row["configuration"]].add(int(row["trainable_parameters"]))
    inconsistent = {config: values for config, values in result.items() if len(values) != 1}
    if inconsistent:
        raise RuntimeError(f"Configuration has inconsistent trainable parameter count: {inconsistent}")
    return {config: next(iter(values)) for config, values in result.items()}


def make_fig1(data: dict[str, Any], output: Path) -> tuple[dict[str, str], list[dict[str, Any]]]:
    assignments = data["assignments"]
    counts = Counter(row["assignment_role"] for row in assignments)
    source_studies = len({row["sra_study_accession"] for row in data["training_groups"] if row["sra_study_accession"]})
    source_runs = data["p11_spec"]["counts"]["source_runs"]
    target_groups = data["p11_spec"]["counts"]["source_groups"]
    fig = plt.figure(figsize=(WIDTH_MM / 25.4, 132 / 25.4))
    grid = fig.add_gridspec(
        2,
        2,
        left=0.075,
        right=0.985,
        bottom=0.075,
        top=0.94,
        width_ratios=(0.84, 1.16),
        height_ratios=(1.0, 1.0),
        wspace=0.34,
        hspace=0.38,
    )
    source_rows: list[dict[str, Any]] = []

    ax = fig.add_subplot(grid[0, 0])
    ax.axis("off")
    panel_label(ax, "a")
    ax.text(0.0, 0.99, "Audited RNA-seq resource", transform=ax.transAxes, fontsize=8, fontweight="bold", va="top")
    inventory = [(str(source_runs), "verified RNA-seq runs", COLORS["ink"]), (str(source_studies), "source studies", COLORS["c_size"]), (str(target_groups), "registered RNA target groups", COLORS["b"])]
    for index, (count, label, color) in enumerate(inventory):
        y = 0.73 - index * 0.22
        ax.text(0.05, y, count, transform=ax.transAxes, fontsize=14, fontweight="bold", color=color, va="center")
        ax.text(0.31, y, label, transform=ax.transAxes, fontsize=6.6, va="center")
        source_rows.append({"panel": "1a", "item": label, "value": count, "note": "Descriptive inventory."})
    ax.text(0.05, 0.08, "WBcel235 | 131,072-bp input | 241 fixed output heads", transform=ax.transAxes, fontsize=6.6, color=COLORS["ink"])

    ax = fig.add_subplot(grid[0, 1])
    ax.axis("off")
    panel_label(ax, "b")
    ax.text(0.0, 0.99, "Target aggregation must follow biological-unit assignment", transform=ax.transAxes, fontsize=8, fontweight="bold", va="top")
    columns = [(0.02, "Invalid: aggregate first", "Held-out unit already\ncontributes to training label", "#F6E6E8"), (0.53, "Valid: assign first", "Held-out unit is excluded\nbefore training aggregation", "#E4F0ED")]
    for x, title, detail, face in columns:
        ax.add_patch(Rectangle((x, 0.12), 0.42, 0.72, transform=ax.transAxes, facecolor=face, edgecolor=COLORS["ink"], linewidth=0.7))
        ax.text(x + 0.02, 0.77, title, transform=ax.transAxes, fontsize=6.0, fontweight="bold")
        members = ["unit 1", "unit 2", "unit 3"]
        for i, member in enumerate(members):
            y = 0.64 - i * 0.13
            fill = COLORS["alert"] if i == 2 else COLORS["train"]
            ax.add_patch(Rectangle((x + 0.04, y), 0.105, 0.07, transform=ax.transAxes, facecolor=fill, edgecolor="none"))
            ax.text(x + 0.165, y + 0.035, member if i < 2 else "held-out unit", transform=ax.transAxes, va="center", fontsize=5.4)
        if x < 0.5:
            for y in (0.675, 0.545, 0.415):
                ax.add_patch(FancyArrowPatch((x + 0.27, y), (x + 0.27, 0.28), transform=ax.transAxes, arrowstyle="-|>", mutation_scale=7, linewidth=0.55, color=COLORS["ink"]))
            ax.add_patch(Rectangle((x + 0.06, 0.22), 0.30, 0.08, transform=ax.transAxes, facecolor=COLORS["paper"], edgecolor=COLORS["ink"], linewidth=0.55))
            ax.text(x + 0.21, 0.26, "training aggregate", transform=ax.transAxes, ha="center", va="center", fontsize=5.2)
        else:
            for y in (0.675, 0.545):
                ax.add_patch(FancyArrowPatch((x + 0.27, y), (x + 0.27, 0.32), transform=ax.transAxes, arrowstyle="-|>", mutation_scale=7, linewidth=0.55, color=COLORS["ink"]))
            ax.add_patch(FancyArrowPatch((x + 0.27, 0.415), (x + 0.35, 0.22), transform=ax.transAxes, arrowstyle="-|>", mutation_scale=7, linewidth=0.55, color=COLORS["alert"]))
            ax.add_patch(Rectangle((x + 0.03, 0.26), 0.24, 0.08, transform=ax.transAxes, facecolor=COLORS["paper"], edgecolor=COLORS["ink"], linewidth=0.55))
            ax.add_patch(Rectangle((x + 0.25, 0.15), 0.15, 0.08, transform=ax.transAxes, facecolor=COLORS["paper"], edgecolor=COLORS["alert"], linewidth=0.55))
            ax.text(x + 0.15, 0.30, "training", transform=ax.transAxes, ha="center", va="center", fontsize=5.0)
            ax.text(x + 0.325, 0.19, "held-out", transform=ax.transAxes, ha="center", va="center", fontsize=5.0)
        ax.text(x + 0.02, 0.045, detail, transform=ax.transAxes, fontsize=5.1, va="top")
    source_rows.extend([
        {"panel": "1b", "item": "split-after-aggregation", "value": "dependent", "note": "Counterfactual schematic."},
        {"panel": "1b", "item": "assignment-before-aggregation", "value": "label-protected", "note": "Completed biological-unit-aware design."},
    ])

    ax = fig.add_subplot(grid[1, 0])
    ax.axis("off")
    panel_label(ax, "c")
    ax.text(0.0, 0.99, "Objects in the evaluation contract", transform=ax.transAxes, fontsize=8, fontweight="bold", va="top")
    levels = ["source\nrun/member", "biological\nunit", "target\ngroup", "training or\nheld-out aggregate", "fixed output\nhead"]
    x_positions = np.linspace(0.02, 0.80, len(levels))
    for index, (x, text) in enumerate(zip(x_positions, levels)):
        ax.add_patch(Rectangle((x, 0.46), 0.15, 0.18, transform=ax.transAxes, facecolor=COLORS["soft"], edgecolor=COLORS["ink"], linewidth=0.6))
        ax.text(x + 0.075, 0.55, text, transform=ax.transAxes, ha="center", va="center", fontsize=5.7)
        if index < len(levels) - 1:
            ax.add_patch(FancyArrowPatch((x + 0.15, 0.55), (x_positions[index + 1], 0.55), transform=ax.transAxes, arrowstyle="-|>", mutation_scale=7, linewidth=0.6, color=COLORS["ink"]))
    ax.text(0.02, 0.22, "The 57 primary aggregates are represented output targets with a biological unit held out of label construction.", transform=ax.transAxes, fontsize=6.1, wrap=True)
    source_rows.append({"panel": "1c", "item": "hierarchy", "value": "run -> unit -> group -> aggregate -> fixed head", "note": "Terminology contract."})

    ax = fig.add_subplot(grid[1, 1])
    ax.axis("off")
    panel_label(ax, "d")
    ax.text(0.0, 0.99, "Completed within-collection evaluation", transform=ax.transAxes, fontsize=8, fontweight="bold", va="top")
    primary_evaluated = sum(
        row["assignment_role"] == "primary_candidate"
        and row["group_qc_decision"] == "candidate_keep"
        for row in assignments
    )
    categories = [("training only", counts.get("training_only", 0), COLORS["train"]), ("primary evaluated", primary_evaluated, COLORS["heldout"]), ("pending QC", counts.get("primary_pending_qc", 0), "#E5D9A9"), ("supplementary", counts.get("supplementary_candidate", 0), "#D9D9D9")]
    max_count = max(count for _, count, _ in categories)
    for index, (label, count, color) in enumerate(categories):
        y = 0.76 - index * 0.16
        width = 0.66 * count / max_count
        ax.add_patch(Rectangle((0.03, y), width, 0.09, transform=ax.transAxes, facecolor=color, edgecolor="none"))
        ax.text(0.72, y + 0.045, f"{count:>3}  {label}", transform=ax.transAxes, va="center", fontsize=6.6)
        source_rows.append({"panel": "1d", "item": label, "value": count, "note": "Assignment category count."})
    ax.text(0.03, 0.15, "Primary comparison: five predefined genomic folds; three seeds nested within fold.\nAll 165 registered jobs used reconstructed training aggregates.\nThis is not study/laboratory-isolated external validation.", transform=ax.transAxes, fontsize=6.1, va="top")
    return export_figure(fig, output, "figure_p13_biological_unit_design", 132), source_rows


def make_fig2(data: dict[str, Any], output: Path) -> tuple[dict[str, str], list[dict[str, Any]]]:
    metrics = data["metrics"]
    group_rows = data["group_descriptives"]
    summary = {row["configuration"]: row for row in data["configuration_summary"]}
    comparison = {row["comparison"]: row for row in data["comparison_summary"]}
    raw_effects = data["fold_seed_effects"]
    loo = data["loo"]
    formal = [("A_paper", "Model A", COLORS["a"]), ("B_paper", "Model B", COLORS["b"]), ("C_paper", "Protocol C", COLORS["c"]), ("C_size_matched", "Size-matched C", COLORS["c_size"])]
    fold_scores = grouped_mean([row for row in metrics if row["scope"] == "heldout" and row["metric"] == "primary_biological_score" and row["configuration"] in {item[0] for item in formal}], ("configuration", "fold"), "value")
    fig, axes = plt.subplots(2, 2, constrained_layout=True)
    source_rows: list[dict[str, Any]] = []

    ax = axes[0, 0]
    setup_axis(ax)
    for x, (config, label, color) in enumerate(formal):
        values = [fold_scores[(config, str(fold))] for fold in FOLDS]
        ax.scatter(np.full(len(values), x), values, color=color, alpha=0.65, s=18, edgecolors="white", linewidths=0.35, zorder=2)
        ax.scatter(x, finite(summary[config]["mean_15_fold_seed_scores"]), color=color, marker="D", s=34, edgecolors="white", linewidths=0.45, zorder=3)
        for fold, value in zip(FOLDS, values):
            source_rows.append({"panel": "2a", "configuration": config, "fold": fold, "value": value, "note": "Mean across three nested seeds."})
    ax.set_xticks(range(len(formal)), [label for _, label, _ in formal])
    ax.set_ylabel("Held-out primary score")
    ax.set_ylim(0.43, 0.70)
    ax.set_title("Model ordering across five genomic folds", loc="left", fontweight="bold")
    ax.text(0.02, 0.02, "Circles: fold means; diamonds: equal-weight mean", transform=ax.transAxes, fontsize=5.6)
    panel_label(ax, "a")

    ax = axes[0, 1]
    setup_axis(ax)
    comparison_ids = [("B_vs_A_paper", "B - A"), ("B_vs_C_paper", "B - protocol C"), ("C_size_matched_vs_B", "B - size-matched C")]
    for y, (comparison_id, label) in enumerate(comparison_ids[::-1]):
        row = comparison[comparison_id]
        raw_by_fold = grouped_mean([effect for effect in raw_effects if effect["comparison"] == comparison_id], ("fold",), "b_advantage")
        values = [raw_by_fold[(str(fold),)] for fold in FOLDS]
        ax.scatter(values, np.full(5, y), color=COLORS["positive"], s=16, alpha=0.72, edgecolors="white", linewidths=0.3, zorder=3)
        estimate, low, high = (finite(row[key]) for key in ("reported_hierarchical_bootstrap_estimate", "reported_hierarchical_bootstrap_ci_95_low", "reported_hierarchical_bootstrap_ci_95_high"))
        ax.errorbar(estimate, y, xerr=[[estimate - low], [high - estimate]], fmt="D", color=COLORS["b"], markersize=4.4, capsize=2.2, linewidth=1.0, zorder=4)
        for fold, value in zip(FOLDS, values):
            source_rows.append({"panel": "2b", "comparison": comparison_id, "fold": fold, "value": value, "note": "Fold mean of paired seed differences."})
        source_rows.append({"panel": "2b", "comparison": comparison_id, "estimate": estimate, "ci_95_low": low, "ci_95_high": high, "note": "Registered fold-first, seed-within-fold bootstrap interval."})
    ax.axvline(0, color=COLORS["ink"], linewidth=0.7)
    ax.set_yticks(range(3), [label for _, label in comparison_ids[::-1]])
    ax.set_xlabel("Paired B advantage")
    ax.set_title("Matched improvements survive label protection", loc="left", fontweight="bold")
    panel_label(ax, "b")

    ax = axes[1, 0]
    setup_axis(ax)
    for index, (comparison_id, label) in enumerate(comparison_ids):
        rows = sorted((row for row in loo if row["comparison"] == comparison_id), key=lambda row: int(row["excluded_fold"]))
        values = [finite(row["mean_b_advantage_remaining_four_folds"]) for row in rows]
        ax.plot(range(1, 6), values, marker="o", markersize=3.6, linewidth=1.0, label=label)
        for row, value in zip(rows, values):
            source_rows.append({"panel": "2c", "comparison": comparison_id, "excluded_fold": row["excluded_fold"], "value": value, "note": "Equal-weight mean of remaining four fold means."})
    ax.set_xticks(FOLDS)
    ax.set_xlabel("Fold excluded")
    ax.set_ylabel("Four-fold B advantage")
    ax.set_title("Leave-one-fold-out sensitivity", loc="left", fontweight="bold")
    ax.legend(frameon=False, fontsize=5.5, loc="lower right")
    panel_label(ax, "c")

    ax = axes[1, 1]
    setup_axis(ax)
    model_by_fold = {int(row["fold"]): row for row in data["p14_model_by_fold"]}
    replicate_by_fold = {int(row["fold"]): row for row in data["p13_replicate_by_fold"]}
    if set(model_by_fold) != set(FOLDS) or set(replicate_by_fold) != set(FOLDS):
        raise RuntimeError("P13/P14 I--V reference must contain exactly five matched folds")
    model_values = [finite(model_by_fold[fold]["primary_biological_score"]) for fold in FOLDS]
    replicate_values = [finite(replicate_by_fold[fold]["primary_biological_score"]) for fold in FOLDS]
    for fold, model_value, replicate_value in zip(FOLDS, model_values, replicate_values):
        ax.plot([0, 1], [model_value, replicate_value], color=COLORS["grid"], linewidth=0.9, zorder=1)
        ax.scatter(0, model_value, color=COLORS["b"], s=22, edgecolors="white", linewidths=0.35, zorder=3)
        ax.scatter(1, replicate_value, color=COLORS["ink"], s=22, edgecolors="white", linewidths=0.35, zorder=3)
        source_rows.append({"panel": "2d", "fold": fold, "model_b_primary_score": model_value, "training_aggregate_to_heldout_reference": replicate_value, "model_minus_reference": model_value - replicate_value, "note": "Matched I--V-only contract; measurement-reliability reference, not an external validation or absolute performance ceiling."})
    ax.scatter(0, float(np.mean(model_values)), color=COLORS["b"], marker="D", s=34, edgecolors="white", linewidths=0.45, zorder=4)
    ax.scatter(1, float(np.mean(replicate_values)), color=COLORS["ink"], marker="D", s=34, edgecolors="white", linewidths=0.45, zorder=4)
    ax.set_xticks([0, 1], ["Model B", "Training aggregate\nvs held-out unit"])
    ax.set_xlim(-0.35, 1.35)
    ax.set_ylim(0.45, 1.02)
    ax.set_ylabel("Primary score")
    ax.set_title("Reference measurement agreement exceeds Model B", loc="left", fontweight="bold")
    ax.text(0.02, 0.02, "Circles: folds; diamonds: means\nNot external validation or an absolute ceiling", transform=ax.transAxes, fontsize=5.3, va="bottom")
    panel_label(ax, "d")
    return export_figure(fig, output, "figure_p13_heldout_performance", 128), source_rows


def make_fig3(data: dict[str, Any], output: Path) -> tuple[dict[str, str], list[dict[str, Any]]]:
    counts = parameter_counts(data["job_records"])
    comparison = {row["comparison"]: row for row in data["comparison_summary"]}
    raw_effects = data["fold_seed_effects"]
    paired = data["paired"]
    fig, axes = plt.subplots(2, 2, constrained_layout=True)
    source_rows: list[dict[str, Any]] = []

    ax = axes[0, 0]
    setup_axis(ax)
    controls = [("B_vs_A_paper", "frozen Model A", "A_paper"), ("B_vs_C_paper", "protocol-matched C", "C_paper"), ("C_size_matched_vs_B", "size-matched C", "C_size_matched")]
    for y, (comparison_id, label, config) in enumerate(controls[::-1]):
        row = comparison[comparison_id]
        raw_by_fold = grouped_mean([effect for effect in raw_effects if effect["comparison"] == comparison_id], ("fold",), "b_advantage")
        values = [raw_by_fold[(str(fold),)] for fold in FOLDS]
        estimate, low, high = (finite(row[key]) for key in ("reported_hierarchical_bootstrap_estimate", "reported_hierarchical_bootstrap_ci_95_low", "reported_hierarchical_bootstrap_ci_95_high"))
        ax.scatter(values, np.full(5, y), color=COLORS["positive"], s=14, alpha=0.7, edgecolors="white", linewidths=0.3)
        ax.errorbar(estimate, y, xerr=[[estimate - low], [high - estimate]], fmt="D", color=COLORS["b"], markersize=4.4, capsize=2.2, linewidth=1.0)
        source_rows.append({"panel": "3a", "comparison": comparison_id, "trainable_parameters": counts[config], "estimate": estimate, "ci_95_low": low, "ci_95_high": high, "note": "Positive effect favors Model B."})
    ax.axvline(0, color=COLORS["ink"], linewidth=0.7)
    ax.set_yticks(range(3), [f"{label}\n({counts[config] / 1e6:.2f}M trainable)" for _, label, config in controls[::-1]])
    ax.set_xlabel("B advantage on held-out primary score")
    ax.set_title("Adaptation advantage exceeds parameter count", loc="left", fontweight="bold")
    panel_label(ax, "a")

    ax = axes[0, 1]
    setup_axis(ax)
    components = [("B_no_lora_vs_B", "LoRA removed"), ("B_lora_only_vs_B", "Added C. elegans\nembedding (conditional)"), ("B_1bp_head_only_vs_B", "Learned 128-bp\nhead")]
    for y, (comparison_id, label) in enumerate(components[::-1]):
        row = comparison[comparison_id]
        raw_by_fold = grouped_mean([effect for effect in raw_effects if effect["comparison"] == comparison_id], ("fold",), "b_advantage")
        values = [raw_by_fold[(str(fold),)] for fold in FOLDS]
        estimate, low, high = (finite(row[key]) for key in ("reported_hierarchical_bootstrap_estimate", "reported_hierarchical_bootstrap_ci_95_low", "reported_hierarchical_bootstrap_ci_95_high"))
        color = COLORS["positive"] if low > 0 else COLORS["neutral"]
        ax.scatter(values, np.full(5, y), color=color, s=14, alpha=0.7, edgecolors="white", linewidths=0.3)
        ax.errorbar(estimate, y, xerr=[[estimate - low], [high - estimate]], fmt="D", color=color, markersize=4.4, capsize=2.2, linewidth=1.0)
        source_rows.append({"panel": "3b", "comparison": comparison_id, "estimate": estimate, "ci_95_low": low, "ci_95_high": high, "note": "Registered configuration contrast; not a biological mechanism."})
    ax.axvline(0, color=COLORS["ink"], linewidth=0.7)
    ax.set_yticks(range(3), [label for _, label in components[::-1]])
    ax.set_xlabel("Full Model B advantage")
    ax.set_title("LoRA is the largest registered contrast", loc="left", fontweight="bold")
    panel_label(ax, "b")

    ax = axes[1, 0]
    setup_axis(ax)
    for comparison_id, label in components:
        raw_by_fold = grouped_mean([effect for effect in raw_effects if effect["comparison"] == comparison_id], ("fold",), "b_advantage")
        values = [raw_by_fold[(str(fold),)] for fold in FOLDS]
        ax.plot(FOLDS, values, marker="o", markersize=3.5, linewidth=1.0, label=label.replace("\n", " "))
        for fold, value in zip(FOLDS, values):
            source_rows.append({"panel": "3c", "comparison": comparison_id, "fold": fold, "value": value, "note": "Fold mean of paired seed differences."})
    ax.axhline(0, color=COLORS["ink"], linewidth=0.7)
    ax.set_xticks(FOLDS)
    ax.set_xlabel("Genomic fold")
    ax.set_ylabel("Full Model B advantage")
    ax.set_title("Configuration contrasts across folds", loc="left", fontweight="bold")
    ax.legend(frameon=False, fontsize=5.3, loc="upper right")
    panel_label(ax, "c")

    ax = axes[1, 1]
    setup_axis(ax)
    endpoint_rows: list[tuple[str, str, float, float, float]] = []
    endpoint_labels = [("gene_exon_coverage_pearson_log1p", "Gene-exon Pearson"), ("per_track_pearson_128bp_log1p", "128-bp Pearson")]
    for component_index, (comparison_id, label) in enumerate(components):
        invert = True
        for endpoint_index, (metric, metric_label) in enumerate(endpoint_labels):
            matches = [row for row in paired if row["scope"] == "heldout" and row["comparison"] == comparison_id and row["metric"] == metric]
            if len(matches) != 1:
                raise RuntimeError(f"Missing P12 paired endpoint {comparison_id}/{metric}")
            estimate, low, high = paired_effect(matches[0], invert=invert)
            y = component_index + (-0.13 if endpoint_index == 0 else 0.13)
            color = COLORS["b"] if endpoint_index == 0 else COLORS["c_size"]
            ax.errorbar(estimate, y, xerr=[[estimate - low], [high - estimate]], fmt="o", color=color, markersize=3.8, capsize=2.0, linewidth=0.9, label=metric_label if component_index == 0 else None)
            endpoint_rows.append((comparison_id, metric, estimate, low, high))
    ax.axvline(0, color=COLORS["ink"], linewidth=0.7)
    ax.set_yticks(range(len(components)), [label for _, label in components])
    ax.set_xlabel("Full Model B advantage")
    ax.set_title("Component effects across primary endpoints", loc="left", fontweight="bold")
    ax.legend(frameon=False, fontsize=5.3, loc="lower right")
    panel_label(ax, "d")
    source_rows.extend({"panel": "3d", "comparison": comparison_id, "metric": metric, "estimate": estimate, "ci_95_low": low, "ci_95_high": high, "note": "Positive effect favors full Model B."} for comparison_id, metric, estimate, low, high in endpoint_rows)
    return export_figure(fig, output, "figure_p13_configuration_attribution", 128), source_rows


def window_pearson(rows: list[dict[str, str]]) -> float:
    observed = np.log1p(np.asarray([finite(row["observed_128bp"]) for row in rows]))
    predicted = np.log1p(np.asarray([finite(row["predicted_128bp"]) for row in rows]))
    return float(np.corrcoef(observed, predicted)[0, 1])


def make_fig4(data: dict[str, Any], output: Path) -> tuple[dict[str, str], list[dict[str, Any]]]:
    group_rows = data["group_descriptives"]
    coverage = data["coverage_source"]
    selection = data["coverage_manifest"]["selected_tracks"]
    fig = plt.figure(constrained_layout=True)
    grid = fig.add_gridspec(3, 2, width_ratios=(1.55, 0.85), hspace=0.15, wspace=0.30)
    source_rows: list[dict[str, Any]] = []

    for index, selected in enumerate(selection):
        ax = fig.add_subplot(grid[index, 0])
        setup_axis(ax)
        rows = [row for row in coverage if row["group_id"] == selected["group_id"]]
        rows.sort(key=lambda row: int(row["genomic_start"]))
        x = np.asarray([(int(row["genomic_start"]) - int(rows[0]["genomic_start"])) / 1000 for row in rows])
        observed = np.log1p(np.asarray([finite(row["observed_128bp"]) for row in rows]))
        predicted = np.log1p(np.asarray([finite(row["predicted_128bp"]) for row in rows]))
        ax.plot(x, observed, color="#4B4B4B", linewidth=0.75, label="Observed" if index == 0 else None)
        ax.plot(x, predicted, color=COLORS["b"], linewidth=0.75, label="Predicted" if index == 0 else None)
        local_r = window_pearson(rows)
        ax.set_ylabel("log1p coverage")
        if index == len(selection) - 1:
            ax.set_xlabel("Position within selected 131-kb window (kb)")
        else:
            ax.tick_params(labelbottom=False)
        ax.set_title(f"{selected['selection_quantile'].title()} global track quantile | {selected['group_id']} | window r = {local_r:.3f}", loc="left", fontsize=6.7, fontweight="bold")
        if index == 0:
            ax.legend(frameon=False, ncol=2, loc="upper right")
            panel_label(ax, "a")
        for row in rows:
            source_rows.append({"panel": "4a", **row, "window_pearson_log1p": local_r, "note": "Selected-window Pearson differs from the global track-quantile selection metric."})

    ax = fig.add_subplot(grid[0, 1])
    setup_axis(ax)
    stage_order = ["Embryo", "L1", "L2", "L3", "L4", "Young Adult"]
    stage_values = []
    labels = []
    for stage in stage_order:
        values = [finite(row["mean_primary_biological_score"]) for row in group_rows if row["development_stage"] == stage]
        if values:
            stage_values.append(values)
            labels.append(f"{stage}\n(n={len(values)})")
    boxes = ax.boxplot(stage_values, positions=np.arange(len(stage_values)), widths=0.48, patch_artist=True, showfliers=False, medianprops={"color": COLORS["ink"], "linewidth": 0.8}, whiskerprops={"color": COLORS["ink"], "linewidth": 0.55}, capprops={"color": COLORS["ink"], "linewidth": 0.55})
    for index, patch in enumerate(boxes["boxes"]):
        patch.set_facecolor(COLORS["stage"][index])
        patch.set_alpha(0.32)
    for index, values in enumerate(stage_values):
        offsets = np.linspace(-0.12, 0.12, len(values))
        ax.scatter(np.full(len(values), index) + offsets, values, s=8, color=COLORS["stage"][index], alpha=0.65, edgecolors="none")
    ax.set_xticks(range(len(labels)), labels, fontsize=5.0)
    ax.set_ylabel("Primary score")
    ax.set_title("Per-group score heterogeneity", loc="left", fontweight="bold")
    ax.text(0.01, 0.02, "Descriptive groups", transform=ax.transAxes, fontsize=5.2)
    panel_label(ax, "b")
    for row in group_rows:
        source_rows.append({"panel": "4b", "group_id": row["group_id"], "development_stage": row["development_stage"], "value": row["mean_primary_biological_score"], "note": "Descriptive mean across 15 fold-seed records."})

    ax = fig.add_subplot(grid[1, 1])
    setup_axis(ax)
    ratio_values = [finite(row["mean_top1_calibration_ratio_128bp"]) for row in group_rows]
    if np.any(np.asarray(ratio_values) <= 0):
        raise ValueError("Non-positive top-1% calibration ratios are prohibited before logarithmic plotting")
    ax.scatter(np.arange(len(ratio_values)), sorted(ratio_values), color=COLORS["b"], s=10, alpha=0.75, edgecolors="none")
    ax.axhline(1.0, color=COLORS["ink"], linewidth=0.8, linestyle="--")
    ax.set_yscale("log")
    ax.set_ylim(0.045, 1.5)
    ax.set_yticks([0.05, 0.1, 0.2, 1.0], ["0.05", "0.1", "0.2", "1.0"])
    ax.set_xlabel("Primary groups, ordered")
    ax.set_ylabel("Predicted / observed")
    ax.set_title("Top-1% observed-bin amplitude", loc="left", fontweight="bold")
    ax.text(0.02, 0.05, "Reference = 1", transform=ax.transAxes, fontsize=5.3)
    panel_label(ax, "c")
    for row in group_rows:
        source_rows.append({"panel": "4c", "group_id": row["group_id"], "value": row["mean_top1_calibration_ratio_128bp"], "note": "Mean prediction / mean observed target on observed top-1% 128-bp bins; descriptive."})

    ax = fig.add_subplot(grid[2, 1])
    setup_axis(ax)
    stage_color = {stage: COLORS["stage"][index] for index, stage in enumerate(stage_order)}
    for stage in stage_order:
        rows = [row for row in group_rows if row["development_stage"] == stage]
        if not rows:
            continue
        x = [finite(row["abs_log2_top1_calibration_deviation"]) for row in rows]
        y = [finite(row["mean_primary_biological_score"]) for row in rows]
        ax.scatter(x, y, s=15, color=stage_color[stage], alpha=0.75, edgecolors="white", linewidths=0.25, label=stage)
        for row, x_value, y_value in zip(rows, x, y):
            source_rows.append({"panel": "4d", "group_id": row["group_id"], "development_stage": row["development_stage"], "abs_log2_top1_calibration_deviation": x_value, "primary_score": y_value, "note": "Descriptive failure-mode map."})
    ax.set_xlabel("|log2 top-1% calibration ratio|")
    ax.set_ylabel("Primary score")
    ax.set_title("Correlation does not remove amplitude error", loc="left", fontweight="bold")
    ax.legend(frameon=False, fontsize=5.0, ncol=2, loc="lower left")
    panel_label(ax, "d")
    return export_figure(fig, output, "figure_p13_fidelity_failure_modes", 152), source_rows


def run_qa(output: Path, source_script: Path) -> dict[str, Any]:
    validator = SKILL_ROOT / "scripts/validate_figure.py"
    pdf_auditor = SKILL_ROOT / "scripts/audit_pdf_text.py"
    result: dict[str, Any] = {"status": "passed", "backend": "python_matplotlib", "created_at": utc_now(), "figures": {}}
    source_validation = subprocess.run([sys.executable, str(validator), str(source_script), "--backend", "python", "--json"], capture_output=True, text=True, check=False)
    try:
        source_validation_payload = json.loads(source_validation.stdout)
    except json.JSONDecodeError:
        source_validation_payload = {"status": "invalid", "stdout": source_validation.stdout, "stderr": source_validation.stderr}
    if source_validation.returncode != 0:
        result["status"] = "failed"
    for pdf in sorted(output.glob("*.pdf")):
        stem = pdf.stem
        audit = subprocess.run([sys.executable, str(pdf_auditor), str(pdf), "--min-pt", "5", "--json"], capture_output=True, text=True, check=False)
        try:
            audit_payload = json.loads(audit.stdout)
        except json.JSONDecodeError:
            audit_payload = {"status": "invalid", "stdout": audit.stdout, "stderr": audit.stderr}
        png = output / f"{stem}.png"
        tiff = output / f"{stem}.tiff"
        raster: dict[str, Any] = {}
        for path in (png, tiff):
            with Image.open(path) as image:
                raster[path.suffix.lstrip(".")] = {"size": list(image.size), "dpi": [float(value) if value is not None else None for value in image.info.get("dpi", (None, None))], "mode": image.mode}
        if audit.returncode != 0:
            result["status"] = "failed"
        result["figures"][stem] = {"pdf": str(pdf.relative_to(REPO_ROOT)), "pdf_sha256": sha256(pdf), "svg_sha256": sha256(output / f"{stem}.svg"), "png_sha256": sha256(png), "tiff_sha256": sha256(tiff), "source_validation": source_validation_payload, "pdf_text_audit": audit_payload, "raster": raster}
    return result


def main() -> None:
    args = parse_args()
    output = args.output_dir.resolve()
    data = load_data(validate_inputs())
    output.mkdir(parents=True, exist_ok=True)
    figures: dict[str, dict[str, str]] = {}
    all_source_rows: list[dict[str, Any]] = []
    for name, maker in (("fig1", make_fig1), ("fig2", make_fig2), ("fig3", make_fig3), ("fig4", make_fig4)):
        paths, source_rows = maker(data, output)
        figures[name] = paths
        all_source_rows.extend(source_rows)
    fields = sorted({field for row in all_source_rows for field in row})
    source_path = output / "figure_source_data.tsv"
    write_tsv(source_path, all_source_rows, fields)
    inputs = [*data["paths"].values(), P12_ROOT / "p12_replicate_holdout_execution.json"]
    manifest = {
        "schema_version": 1,
        "created_at": utc_now(),
        "backend": "python_matplotlib",
        "figure_width_mm": WIDTH_MM,
        "raster_dpi": RASTER_DPI,
        "figures": figures,
        "source_data": str(source_path.relative_to(REPO_ROOT)),
        "scope": "completed within-collection P11/P12 evidence plus completed P13/P14 I--V matched measurement-reliability reference",
        "prohibited_claims": ["No external study/laboratory result is plotted.", "No regulatory-variant result is plotted.", "The replicate-agreement reference is not an absolute performance ceiling.", "No new signal, checkpoint or chromosome-X read occurred during rendering."],
        "statistics": {"primary_estimand": "equal-weight paired difference across five predefined genomic folds", "nested_repeats": "three seeds within fold", "track_status": "57 primary groups are descriptive outcomes, not independent inferential replicates"},
        "input_sha256": {str(path.relative_to(REPO_ROOT)): sha256(path) for path in inputs},
    }
    manifest_path = output / "manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    qa = run_qa(output, Path(__file__))
    qa["manifest_sha256"] = sha256(manifest_path)
    qa_path = output / "figure_qa.json"
    qa_path.write_text(json.dumps(qa, indent=2, sort_keys=True) + "\n")
    if qa["status"] != "passed":
        raise SystemExit("P13 figure QA failed")
    print(json.dumps({"status": "passed", "output": str(output.relative_to(REPO_ROOT)), "figures": figures, "source_data": str(source_path.relative_to(REPO_ROOT)), "qa": str(qa_path.relative_to(REPO_ROOT))}, indent=2))


if __name__ == "__main__":
    main()
