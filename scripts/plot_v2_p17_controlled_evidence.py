#!/usr/bin/env python3
"""Render the P17 controlled-comparison manuscript addendum with Python."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image


REPO_ROOT = Path(__file__).resolve().parents[1]
DATA_ROOT = REPO_ROOT / "results/v2_p17_manuscript_evidence_20260823T032000Z"
REVIEW_PATH = REPO_ROOT / "alphagenome_custom/metadata/v2/audits/P17/review.json"
DEFAULT_OUTPUT = REPO_ROOT / "results/v2_p17_manuscript_figure_20260823T032000Z"
WIDTH_MM = 183
DPI = 600
CONFIGS = ("B_iv_dual", "B_128bp_only", "Basenji2_style")
LABELS = {"B_iv_dual": "Model B, dual resolution", "B_128bp_only": "Model B, 128-bp only", "Basenji2_style": "Basenji2-style scratch"}
COLORS = {"B_iv_dual": "#007C91", "B_128bp_only": "#B7791F", "Basenji2_style": "#6B7280"}

mpl.rcParams.update({"font.family": "Liberation Sans", "font.sans-serif": ["Liberation Sans", "Arial", "DejaVu Sans"], "font.size": 7.2, "axes.labelsize": 7, "axes.titlesize": 8, "xtick.labelsize": 6.5, "ytick.labelsize": 6.5, "axes.linewidth": 0.65, "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none", "figure.facecolor": "white", "axes.facecolor": "white"})


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def setup(axis: plt.Axes) -> None:
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    axis.grid(axis="y", color="#D9D9D9", linewidth=0.45)
    axis.set_axisbelow(True)


def panel_label(axis: plt.Axes, label: str) -> None:
    axis.text(-0.14, 1.08, label, transform=axis.transAxes, fontsize=9, fontweight="bold", va="top")


def main() -> None:
    args = parse_args()
    output = args.output_dir.resolve()
    if output.exists():
        raise FileExistsError(f"P17 figure output is immutable: {output}")
    source_path = DATA_ROOT / "p17_figure_source_data.tsv"
    manifest_path = DATA_ROOT / "p17_manuscript_evidence_manifest.json"
    if not source_path.is_file() or not manifest_path.is_file() or not REVIEW_PATH.is_file():
        raise FileNotFoundError("P17 source data, manifest, or review is missing")
    if json.loads(REVIEW_PATH.read_text()).get("status") != "PASS":
        raise RuntimeError("P17 figure requires R17 PASS")
    source = read_tsv(source_path)
    configuration = sorted((row for row in source if row["table"] == "fold_configuration" and row["configuration"] in CONFIGS), key=lambda row: (int(row["fold"]), CONFIGS.index(row["configuration"])))
    paired = sorted((row for row in source if row["table"] == "fold_paired_effect" and row["comparison"] in {"B_iv_dual_minus_B_128bp_only", "B_iv_dual_minus_Basenji2_style"}), key=lambda row: (row["comparison"], int(row["fold"])))
    if len(configuration) != 15 or len(paired) != 10:
        raise RuntimeError("P17 source data does not contain the registered full fold matrix")
    fig = plt.figure(figsize=(WIDTH_MM / 25.4, 62 / 25.4))
    grid = fig.add_gridspec(1, 2, left=0.09, right=0.985, bottom=0.21, top=0.86, width_ratios=(1.25, 1.0), wspace=0.36)
    axis = fig.add_subplot(grid[0, 0])
    panel_label(axis, "a")
    axis.set_title("Strict I-V controlled comparison", loc="left", fontweight="bold")
    by_config = {config: {int(row["fold"]): float(row["primary_biological_score"]) for row in configuration if row["configuration"] == config} for config in CONFIGS}
    folds = np.arange(1, 6)
    offsets = {"B_iv_dual": -0.18, "B_128bp_only": 0.0, "Basenji2_style": 0.18}
    for config in CONFIGS:
        values = [by_config[config][fold] for fold in folds]
        axis.plot(folds + offsets[config], values, color=COLORS[config], linewidth=1.0, marker="o", markersize=4.2, label=LABELS[config], zorder=3)
    for fold in folds:
        axis.plot([fold + offsets["B_iv_dual"], fold + offsets["B_128bp_only"]], [by_config["B_iv_dual"][fold], by_config["B_128bp_only"][fold]], color="#B7B7B7", linewidth=0.55, zorder=1)
        axis.plot([fold + offsets["B_iv_dual"], fold + offsets["Basenji2_style"]], [by_config["B_iv_dual"][fold], by_config["Basenji2_style"][fold]], color="#D7D7D7", linewidth=0.55, zorder=1)
    setup(axis)
    axis.set_xlim(0.55, 5.45)
    axis.set_xticks(folds)
    axis.set_xlabel("Predefined genomic fold")
    axis.set_ylabel("Primary biological score")
    axis.set_ylim(0.50, 0.69)
    axis.legend(loc="lower left", frameon=False, fontsize=5.4, handlelength=1.7)
    axis.text(0.02, -0.34, "Each point is the mean of three nested seeds.\nI-V only; within-collection biological-unit holdout.", transform=axis.transAxes, fontsize=5.5, va="top")

    axis = fig.add_subplot(grid[0, 1])
    panel_label(axis, "b")
    axis.set_title("Paired effect across five folds", loc="left", fontweight="bold")
    comparisons = [("B_iv_dual_minus_Basenji2_style", "Dual B - Basenji2-style", "#B23A48"), ("B_iv_dual_minus_B_128bp_only", "Dual B - 128-bp-only B", "#8E7A3F")]
    for index, (comparison, label, color) in enumerate(comparisons):
        values = [float(row["difference_primary_biological_score"]) for row in paired if row["comparison"] == comparison]
        jitter = np.linspace(-0.07, 0.07, len(values))
        axis.scatter(np.full(len(values), index) + jitter, values, s=23, color=color, edgecolor="white", linewidth=0.45, zorder=3)
        mean = float(np.mean(values))
        sd = float(np.std(values, ddof=1))
        axis.errorbar(index, mean, yerr=sd, color="#252525", linewidth=0.8, capsize=2.0, zorder=4)
        axis.hlines(mean, index - 0.19, index + 0.19, color="#252525", linewidth=1.1, zorder=5)
        axis.text(index, mean + sd + 0.004, f"mean {mean:+.3f}\nSD {sd:.3f}", ha="center", va="bottom", fontsize=5.3)
    setup(axis)
    axis.axhline(0.0, color="#252525", linewidth=0.7)
    axis.set_xlim(-0.45, 1.45)
    axis.set_xticks([0, 1], ["vs scratch", "vs 128-bp-only"])
    axis.set_ylabel("Paired difference in primary score")
    axis.set_ylim(-0.035, 0.10)
    axis.text(0.02, -0.34, "Positive values favour dual-resolution Model B.\nPoints: folds; bars: mean +/- fold SD (n=5).", transform=axis.transAxes, fontsize=5.5, va="top")

    output.mkdir(parents=True)
    paths = {
        "pdf": output / "figure_p17_controlled_comparison.pdf",
        "svg": output / "figure_p17_controlled_comparison.svg",
        "png": output / "figure_p17_controlled_comparison.png",
        "tiff": output / "figure_p17_controlled_comparison.tiff",
    }
    fig.savefig(output / "figure_p17_controlled_comparison.pdf", bbox_inches="tight", facecolor="white")
    fig.savefig(output / "figure_p17_controlled_comparison.svg", bbox_inches="tight", facecolor="white")
    fig.savefig(output / "figure_p17_controlled_comparison.png", dpi=DPI, bbox_inches="tight", facecolor="white")
    fig.savefig(output / "figure_p17_controlled_comparison.tiff", dpi=DPI, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    with Image.open(paths["png"]) as image:
        if image.width < 1000 or image.height < 500:
            raise RuntimeError("P17 figure raster export is unexpectedly small")
    qa = {
        "schema_version": 1,
        "backend": "python/matplotlib",
        "figure_contract": {
            "core_conclusion": "The strict I-V dual-resolution Model B outperforms the task-matched Basenji2-style scratch baseline, while the 128-bp-only Model B has a higher registered primary score than dual-resolution B.",
            "archetype": "quantitative grid",
            "hero_evidence": "All five paired fold-level primary scores for the three registered configurations.",
            "boundary": "The 128-bp-only comparison does not establish biological one-base utility because its 1-bp values were metric-compatibility expansions.",
            "statistics": "Five predefined genomic folds are the top-level units; three seeds are nested means. No track-level pseudo-replication or population CI is displayed.",
        },
        "input_source_data": {"path": str(source_path.relative_to(REPO_ROOT)), "sha256": sha256(source_path)},
        "output": {name: str(path.relative_to(REPO_ROOT)) for name, path in paths.items()},
        "exclusions": "None; all registered configurations and five folds are displayed.",
        "image_integrity": "No biological raster images; all marks are rendered directly from tabular P17 source data.",
    }
    (output / "figure_p17_qa.json").write_text(json.dumps(qa, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"status": "rendered", "output": str(output.relative_to(REPO_ROOT)), "figure": str(paths["pdf"].relative_to(REPO_ROOT))}, indent=2))


if __name__ == "__main__":
    main()
