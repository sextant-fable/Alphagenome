"""Render descriptive calibration and metadata-stratification evidence."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = ROOT / "results/v2_p13_heldout_descriptives/p13_b_paper_primary_group_descriptives.tsv"


def export(output: Path, source: pd.DataFrame, strata: pd.DataFrame) -> None:
    output.mkdir(parents=True, exist_ok=True)
    source.to_csv(output / "ed6_group_source_data.tsv", sep="\t", index=False)
    strata.to_csv(output / "ed6_strata_source_data.tsv", sep="\t", index=False)

    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 7.2,
        "axes.titlesize": 8.2,
        "axes.labelsize": 7.2,
        "xtick.labelsize": 6.2,
        "ytick.labelsize": 6.5,
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
    })
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 3.05), dpi=600)
    fig.subplots_adjust(left=0.08, right=0.96, bottom=0.25, top=0.78, wspace=0.58)
    stage_colors = {
        "Embryo": "#3b6ea8", "L1": "#4c956c", "L2": "#d77a2f",
        "L3": "#8b5fbf", "L4": "#b4475d", "Young Adult": "#444444",
    }
    colors = [stage_colors.get(str(v), "#777777") for v in source["development_stage"]]

    ax = axes[0]
    ax.scatter(source["mean_primary_biological_score"], source["mean_top1_calibration_ratio_128bp"],
               c=colors, s=20, alpha=0.85, edgecolors="white", linewidths=0.3)
    ax.axhline(1.0, color="#777777", lw=0.8, ls="--")
    ax.set_xlabel("Primary score")
    ax.set_ylabel("Top-1% predicted / observed")
    ax.set_title("Group-level calibration")
    ax.text(0.03, 0.96, "n=57 groups", transform=ax.transAxes, va="top", fontsize=6.5)

    ax = axes[1]
    stages = [s for s in ["Embryo", "L1", "L2", "L3", "L4", "Young Adult"] if s in set(source["development_stage"])]
    values = [source.loc[source["development_stage"] == s, "mean_top1_calibration_ratio_128bp"].to_numpy() for s in stages]
    positions = np.arange(1, len(stages) + 1)
    ax.boxplot(values, positions=positions, widths=0.55, showfliers=False,
               medianprops={"color": "#111111", "lw": 1.2},
               boxprops={"color": "#555555"}, whiskerprops={"color": "#555555"},
               capprops={"color": "#555555"})
    for pos, vals, stage in zip(positions, values, stages):
        ax.scatter(np.full(len(vals), pos), vals, s=12, color=stage_colors[stage], alpha=0.75, zorder=3)
    ax.axhline(1.0, color="#777777", lw=0.8, ls="--")
    ax.set_xticks(positions, stages)
    for tick in ax.get_xticklabels():
        tick.set_rotation(35)
        tick.set_ha("right")
        tick.set_rotation_mode("anchor")
    ax.set_ylabel("Top-1% ratio")
    ax.set_title("Developmental stage")

    ax = axes[2]
    field = "sra_study_accession"
    counts = source[field].value_counts()
    selected = counts[counts >= 3].index.tolist()
    selected = sorted(selected, key=lambda v: (-counts[v], v))
    med = [source.loc[source[field] == s, "mean_primary_biological_score"].median() for s in selected]
    ratio = [source.loc[source[field] == s, "mean_top1_calibration_ratio_128bp"].median() for s in selected]
    y = np.arange(len(selected))
    ax.scatter(med, y, s=24, color="#2f718e", label="Score median")
    ax2 = ax.twiny()
    ax2.scatter(ratio, y, s=24, color="#c47f18", label="Calibration median")
    ax.set_yticks(y, selected)
    ax.set_xlabel("Study median primary score", color="#2f718e")
    ax2.set_xlabel("Study median top-1% ratio", color="#c47f18")
    ax.set_title("Study strata (n >= 3)")
    ax.grid(axis="y", alpha=0.15)
    fig.suptitle("Calibration and metadata strata in biological-unit holdout", x=0.01, ha="left", fontsize=10.5, fontweight="bold")
    fig.text(0.08, 0.06, "Ratios below one indicate underestimation in the observed highest 1% of 128-bp bins; no new biological replicate is inferred.", fontsize=6.2, color="#555555")
    fig.savefig(output / "figure_ed6_calibration_strata.pdf", bbox_inches="tight", pad_inches=0.04)
    fig.savefig(output / "figure_ed6_calibration_strata.svg", bbox_inches="tight", pad_inches=0.04)
    fig.savefig(output / "figure_ed6_calibration_strata.png", dpi=600, bbox_inches="tight", pad_inches=0.04)
    fig.savefig(output / "figure_ed6_calibration_strata.tiff", dpi=600, bbox_inches="tight", pad_inches=0.04)
    plt.close(fig)
    (output / "figure_ed6_manifest.json").write_text(json.dumps({
        "figure": "ed6_calibration_strata",
        "source": str(DEFAULT_INPUT.relative_to(ROOT)),
        "rows": int(len(source)),
        "study_strata_plotted": int(len(selected)),
        "status": "descriptive_existing_p13_data",
        "claim_boundary": "not population-level inference; no new training or locked-test reads",
    }, indent=2) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    source = pd.read_csv(args.input, sep="\t")
    required = {
        "sra_study_accession", "development_stage", "tissue_or_cell_type", "condition",
        "mean_primary_biological_score", "mean_top1_calibration_ratio_128bp",
    }
    missing = required - set(source.columns)
    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}")
    strata_rows = []
    for field in ["sra_study_accession", "development_stage", "tissue_or_cell_type", "condition"]:
        for value, group in source.groupby(field, dropna=False):
            strata_rows.append({
                "stratum_field": field,
                "stratum_value": str(value),
                "group_count": len(group),
                "median_primary_score": group["mean_primary_biological_score"].median(),
                "median_top1_calibration_ratio": group["mean_top1_calibration_ratio_128bp"].median(),
                "mean_top1_calibration_ratio": group["mean_top1_calibration_ratio_128bp"].mean(),
            })
    export(args.output, source, pd.DataFrame(strata_rows))


if __name__ == "__main__":
    main()
