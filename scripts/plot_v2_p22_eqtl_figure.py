#!/usr/bin/env python3
"""Plot the public eQTL candidate effect-association figure."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[1]
SUMMARY = REPO_ROOT / "results/v2_p22_eqtl_effect_association/eqtl_variant_effects.tsv"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def finite(row: dict[str, str], key: str) -> float:
    try:
        value = float(row[key])
    except (KeyError, TypeError, ValueError):
        return float("nan")
    return value if np.isfinite(value) else float("nan")


def save_figure(fig: plt.Figure, stem: Path) -> None:
    stem.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(f"{stem}.svg", bbox_inches="tight")
    fig.savefig(f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(f"{stem}.tiff", dpi=600, bbox_inches="tight")
    fig.savefig(f"{stem}.png", dpi=300, bbox_inches="tight")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--summary", type=Path, default=SUMMARY)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    rows = read_tsv(args.summary)
    usable = [
        row for row in rows
        if int(row.get("strain_pairs", 0)) >= 20
        and int(row.get("model_jobs_mapped", 0)) == 15
        and np.isfinite(finite(row, "observed_signed_slope"))
        and np.isfinite(finite(row, "model_gene_body_delta_mean"))
    ]
    if not usable:
        raise RuntimeError("No analyzable public eQTL candidate rows")

    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
        "font.size": 7.2,
        "axes.spines.right": False,
        "axes.spines.top": False,
        "axes.linewidth": 0.65,
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "figure.facecolor": "white",
        "axes.facecolor": "white",
    })
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.55), gridspec_kw={"width_ratios": [1.25, 1.0, 1.0]})
    ax = axes[0]
    x = np.array([finite(row, "observed_signed_slope") for row in usable])
    y = np.array([finite(row, "model_gene_body_delta_mean") for row in usable])
    colors = ["#2f6f8f" if row.get("eQTL_classification") == "Local eQTL" else "#8c6bb1" for row in usable]
    ax.scatter(x, y, s=22, c=colors, edgecolor="white", linewidth=0.35, alpha=0.9, zorder=3)
    lim = max(np.max(np.abs(x)), np.max(np.abs(y)), 1e-6) * 1.12
    ax.axhline(0, color="#9b9b9b", lw=0.7, zorder=1)
    ax.axvline(0, color="#9b9b9b", lw=0.7, zorder=1)
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
    ax.set_xlabel("Observed dosage-expression slope\n(log1p CPM per allele)")
    ax.set_ylabel("Model gene-body delta\n(ALT minus REF)")
    ax.set_title("Observed versus predicted sign", loc="left", fontweight="bold")
    ax.text(0.02, 0.98, f"n={len(usable)} candidates", transform=ax.transAxes, va="top", fontsize=6.5)
    ax.text(0.03, 0.05, "same sign", transform=ax.transAxes, color="#555555", fontsize=6)
    ax.legend(handles=[
        mpl.lines.Line2D([], [], marker="o", linestyle="", color="#2f6f8f", label="Local eQTL", markersize=4),
        mpl.lines.Line2D([], [], marker="o", linestyle="", color="#8c6bb1", label="Distant eQTL", markersize=4),
    ], loc="lower right", fontsize=6, handletextpad=0.3, borderpad=0.2)

    ax = axes[1]
    order = np.argsort([finite(row, "model_observed_sign_concordance") for row in usable])
    concordance = np.array([finite(usable[index], "model_observed_sign_concordance") for index in order])
    ax.scatter(np.arange(len(concordance)), concordance, s=15, color="#2f6f8f", alpha=0.9)
    ax.axhline(0.5, color="#9b9b9b", ls="--", lw=0.8)
    ax.set_ylim(-0.05, 1.05)
    ax.set_xlabel("Candidate rank")
    ax.set_ylabel("15-checkpoint sign concordance")
    ax.set_title("Checkpoint stability", loc="left", fontweight="bold")
    ax.text(0.02, 0.98, "0.5 = chance-level sign split", transform=ax.transAxes, va="top", fontsize=6.5)

    ax = axes[2]
    pairs = np.array([int(row["strain_pairs"]) for row in usable])
    ax.hist(pairs, bins=np.arange(pairs.min() - 0.5, pairs.max() + 1.5, 1), color="#c5a14a", edgecolor="white", linewidth=0.4)
    ax.set_xlabel("Wild strains per candidate")
    ax.set_ylabel("Candidate count")
    ax.set_title("Public-cohort support", loc="left", fontweight="bold")
    ax.text(0.02, 0.98, f"median={int(np.median(pairs))}", transform=ax.transAxes, va="top", fontsize=6.5)

    fig.text(0.01, 1.01, "Public eQTL candidate allele analysis", fontsize=9, fontweight="bold", ha="left")
    fig.text(0.01, -0.03, "Candidate set is I--V only; finemapping-gene mapping and dosage-expression association are descriptive, not causal validation.", fontsize=6.5, color="#555555", ha="left")
    fig.tight_layout(w_pad=1.6)
    save_figure(fig, args.output)
    plt.close(fig)


if __name__ == "__main__":
    main()
