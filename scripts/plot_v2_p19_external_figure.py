#!/usr/bin/env python3
"""Render the R19 external represented-head benchmark for manuscript review."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image


REPO_ROOT = Path(__file__).resolve().parents[1]
REVIEW_PATH = REPO_ROOT / "alphagenome_custom/metadata/v2/audits/P19/review.json"
RECORDS_PATH = REPO_ROOT / "results/v2_p19_external_scoring/p19_external_records.tsv"
BY_SOURCE_PATH = REPO_ROOT / "results/v2_p19_external_scoring/p19_external_by_source.tsv"
DEFAULT_OUTPUT = REPO_ROOT / "results/v2_p19_manuscript_figure_20260823T060118Z"
WIDTH_MM = 183
DPI = 600
SOURCES = ("SRR18463404", "SRR18463405", "SRR18463406", "SRR3560831")
COLORS = {"SRR18463404": "#007C91", "SRR18463405": "#4C78A8", "SRR18463406": "#F58518", "SRR3560831": "#6B7280"}
LABELS = {"SRR18463404": "GSE199326 / run 04", "SRR18463405": "GSE199326 / run 05", "SRR18463406": "GSE199326 / run 06", "SRR3560831": "GSE81707 / run 31"}

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
        raise FileExistsError(f"R19 figure output is immutable: {output}")
    if not REVIEW_PATH.is_file() or json.loads(REVIEW_PATH.read_text()).get("status") != "PASS":
        raise RuntimeError("R19 figure requires a PASS review")
    records = read_tsv(RECORDS_PATH)
    by_source = read_tsv(BY_SOURCE_PATH)
    if len(records) != 60 or len(by_source) != 4:
        raise RuntimeError("R19 source data must contain 60 records and four source summaries")
    if {row["run_accession"] for row in records} != set(SOURCES):
        raise RuntimeError("R19 source set changed")
    output.mkdir(parents=True)
    source_copy = output / "p19_figure_source_data.tsv"
    with source_copy.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(records)
    fig = plt.figure(figsize=(WIDTH_MM / 25.4, 67 / 25.4))
    grid = fig.add_gridspec(1, 3, left=0.08, right=0.99, bottom=0.24, top=0.84, width_ratios=(1.3, 0.95, 0.95), wspace=0.52)

    axis = fig.add_subplot(grid[0, 0]); panel_label(axis, "a"); axis.set_title("External represented-head scores", loc="left", fontweight="bold")
    x = np.arange(1, 5)
    for index, run in enumerate(SOURCES, start=1):
        values = [float(row["primary_biological_score"]) for row in records if row["run_accession"] == run]
        axis.scatter(np.full(len(values), index), values, s=17, color=COLORS[run], edgecolor="white", linewidth=0.4, zorder=3)
        axis.hlines(np.mean(values), index - 0.22, index + 0.22, color="#202020", linewidth=1.0, zorder=4)
    setup(axis); axis.set_xlim(0.5, 4.5); axis.set_ylim(0.60, 0.71); axis.set_xticks(x, ["04", "05", "06", "31"]); axis.set_xlabel("External source run"); axis.set_ylabel("Primary score")
    axis.text(0.02, -0.37, "Points: 15 frozen P17 checkpoints\nHorizontal marks: source mean; no source-level hypothesis test.", transform=axis.transAxes, fontsize=5.3, va="top")

    axis = fig.add_subplot(grid[0, 1]); panel_label(axis, "b"); axis.set_title("Fold-level source means", loc="left", fontweight="bold")
    folds = np.arange(1, 6)
    for run in SOURCES:
        values = [np.mean([float(row["primary_biological_score"]) for row in records if row["run_accession"] == run and int(row["fold"]) == fold]) for fold in folds]
        axis.plot(folds, values, marker="o", markersize=3.5, linewidth=1.0, color=COLORS[run], label=LABELS[run])
    setup(axis); axis.set_xlim(0.8, 5.2); axis.set_ylim(0.60, 0.71); axis.set_xticks(folds); axis.set_xlabel("Predefined genomic fold"); axis.set_ylabel("Primary score")
    axis.legend(frameon=False, fontsize=5.0, loc="lower left", handlelength=1.4)
    axis.text(0.02, -0.37, "Each point averages three nested seeds.\nAll sources map to RNA_V2_G0054.", transform=axis.transAxes, fontsize=5.3, va="top")

    axis = fig.add_subplot(grid[0, 2]); panel_label(axis, "c"); axis.set_title("Metric components", loc="left", fontweight="bold")
    y = np.arange(len(SOURCES))
    gene = [float(row["gene_exon_mean"]) for row in by_source]
    bp = [float(row["pearson_128bp_mean"]) for row in by_source]
    axis.barh(y + 0.18, gene, height=0.30, color="#C7D7D9", label="Gene-exon")
    axis.barh(y - 0.18, bp, height=0.30, color="#B23A48", label="128-bp")
    setup(axis); axis.set_xlim(0.60, 0.72); axis.set_yticks(y, ["04", "05", "06", "31"]); axis.set_xlabel("Pearson correlation (log1p)"); axis.invert_yaxis(); axis.legend(frameon=False, fontsize=5.0, loc="lower right")
    axis.text(0.02, -0.37, "Bars: mean across 15 checkpoints.\nSource labels are descriptive, not independent n=4 inference.", transform=axis.transAxes, fontsize=5.3, va="top")

    paths = {name: output / f"figure_p19_external_represented_head.{name}" for name in ("pdf", "svg", "png", "tiff")}
    fig.savefig(paths["pdf"], bbox_inches="tight", facecolor="white")
    fig.savefig(paths["svg"], bbox_inches="tight", facecolor="white")
    fig.savefig(paths["png"], dpi=DPI, bbox_inches="tight", facecolor="white")
    fig.savefig(paths["tiff"], dpi=DPI, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    with Image.open(paths["png"]) as image:
        if image.width < 1000 or image.height < 400:
            raise RuntimeError("R19 figure raster export is unexpectedly small")
    qa = {"schema_version": 1, "backend": "python/matplotlib", "figure_contract": {"core_conclusion": "Frozen Model B scores remain measurable across four independent external RNA measurements for a represented RNA_V2_G0054 target identity.", "archetype": "quantitative grid", "evidence_class": "R19-verified external source benchmark; represented-head only", "statistics": "Four sources and 15 frozen checkpoints per source are displayed descriptively; no source-level population inference is claimed."}, "input_source_data": {"records": str(RECORDS_PATH.relative_to(REPO_ROOT)), "records_sha256": sha256(RECORDS_PATH), "review": str(REVIEW_PATH.relative_to(REPO_ROOT)), "review_sha256": sha256(REVIEW_PATH)}, "output": {name: str(path.relative_to(REPO_ROOT)) for name, path in paths.items()}, "exclusions": "None; all 60 registered source-by-fold-by-seed observations are represented in panel a and panel b.", "image_integrity": "No biological raster images; all marks are drawn from R19 tabular output."}
    (output / "figure_p19_qa.json").write_text(json.dumps(qa, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"status": "rendered", "output": str(output.relative_to(REPO_ROOT)), "figure": str(paths["pdf"].relative_to(REPO_ROOT))}, indent=2))


if __name__ == "__main__":
    main()
