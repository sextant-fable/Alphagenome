#!/usr/bin/env python3
"""Render a manifest-derived overview of the v2 worm RNA prediction task."""

from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
from typing import Any

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
from PIL import Image


REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MANIFEST = REPO_ROOT / "alphagenome_custom/metadata/v2/rna_seq_groups_v2_final.tsv"
DEFAULT_OUTPUT_DIR = REPO_ROOT / "results/v2_resource_overview_figure"
FIGURE_STEM = "figure_v2_resource_overview"
WIDTH_IN = 183 / 25.4
HEIGHT_IN = 112 / 25.4
RASTER_DPI = 600

COLORS = {
    "ink": "#242424",
    "muted": "#727272",
    "grid": "#dddddd",
    "teal": "#177e89",
    "orange": "#cf873d",
    "red": "#c4473f",
    "purple": "#716fa8",
    "blue": "#4477aa",
    "light_teal": "#dceff0",
    "light_grey": "#f3f3f3",
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
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
        "axes.linewidth": 0.7,
    }
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(1024 * 1024):
            digest.update(block)
    return digest.hexdigest()


def read_manifest(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    required = {
        "group_id",
        "include_formal_v2",
        "n_runs",
        "sra_study_accession",
        "development_stage",
    }
    if not rows or required - set(rows[0]):
        raise ValueError(f"Manifest is missing required columns: {path}")
    formal = [row for row in rows if row["include_formal_v2"] == "True"]
    if len(formal) != 241:
        raise ValueError(f"Expected 241 formal v2 groups, observed {len(formal)}")
    if sum(int(row["n_runs"]) for row in formal) != 482:
        raise ValueError("Expected 482 runs across formal v2 groups")
    return formal


def stage_counts(rows: list[dict[str, str]]) -> list[tuple[str, int]]:
    counts: dict[str, int] = {}
    for row in rows:
        stage = row["development_stage"] or "Not reported"
        counts[stage] = counts.get(stage, 0) + 1
    ordered = sorted(counts.items(), key=lambda item: (-item[1], item[0]))
    top = ordered[:5]
    other = sum(value for _, value in ordered[5:])
    return top + ([("Other stages", other)] if other else [])


def add_panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(
        -0.02,
        1.04,
        label,
        transform=ax.transAxes,
        fontweight="bold",
        fontsize=10,
        va="bottom",
        ha="right",
    )


def draw_box(
    ax: plt.Axes,
    xy: tuple[float, float],
    width: float,
    height: float,
    text: str,
    color: str,
    *,
    text_color: str = COLORS["ink"],
) -> None:
    patch = FancyBboxPatch(
        xy,
        width,
        height,
        boxstyle="round,pad=0.015,rounding_size=0.02",
        linewidth=0.8,
        edgecolor=color,
        facecolor="white",
        transform=ax.transAxes,
        clip_on=False,
    )
    ax.add_patch(patch)
    ax.text(
        xy[0] + width / 2,
        xy[1] + height / 2,
        text,
        transform=ax.transAxes,
        ha="center",
        va="center",
        fontsize=6.6,
        color=text_color,
        linespacing=1.15,
    )


def draw_arrow(ax: plt.Axes, start: tuple[float, float], end: tuple[float, float]) -> None:
    ax.add_patch(
        FancyArrowPatch(
            start,
            end,
            transform=ax.transAxes,
            arrowstyle="-|>",
            mutation_scale=9,
            linewidth=0.8,
            color=COLORS["muted"],
            shrinkA=2,
            shrinkB=2,
        )
    )


def render(rows: list[dict[str, str]]) -> tuple[plt.Figure, list[dict[str, Any]]]:
    figure = plt.figure(figsize=(WIDTH_IN, HEIGHT_IN))
    grid = figure.add_gridspec(
        2,
        2,
        height_ratios=(0.88, 1.12),
        width_ratios=(0.85, 1.15),
        hspace=0.54,
        wspace=0.36,
        left=0.08,
        right=0.98,
        top=0.92,
        bottom=0.12,
    )
    ax_resource = figure.add_subplot(grid[:, 0])
    ax_task = figure.add_subplot(grid[0, 1])
    ax_model = figure.add_subplot(grid[1, 1])
    source_rows: list[dict[str, Any]] = []

    stages = stage_counts(rows)
    stage_palette = [
        COLORS["teal"],
        COLORS["blue"],
        COLORS["orange"],
        COLORS["red"],
        COLORS["purple"],
        COLORS["muted"],
    ]
    total_runs = sum(int(row["n_runs"]) for row in rows)
    studies = len({row["sra_study_accession"] for row in rows})
    ax_resource.set_title("Traceable RNA-seq resource", loc="left", pad=8, fontweight="bold")
    ax_resource.text(0.00, 0.96, f"{total_runs} verified RNA-seq runs", transform=ax_resource.transAxes, fontsize=10, fontweight="bold", color=COLORS["teal"])
    ax_resource.text(0.00, 0.91, f"{studies} studies  |  241 biological RNA tracks", transform=ax_resource.transAxes, fontsize=6.8, color=COLORS["ink"])
    ax_resource.text(0.00, 0.86, "Uniform reprocessing and manifest-defined grouping", transform=ax_resource.transAxes, fontsize=6.3, color=COLORS["muted"])
    names = [name for name, _ in stages][::-1]
    values = [value for _, value in stages][::-1]
    colors = stage_palette[: len(stages)][::-1]
    bars = ax_resource.barh(range(len(names)), values, color=colors, height=0.58)
    ax_resource.set_yticks(range(len(names)), names)
    ax_resource.set_xlabel("Grouped RNA tracks")
    ax_resource.set_xlim(0, max(values) * 1.22)
    ax_resource.set_ylim(-0.6, len(names) + 1.9)
    ax_resource.grid(axis="x", color=COLORS["grid"], linewidth=0.6)
    ax_resource.set_axisbelow(True)
    ax_resource.spines[["top", "right", "left"]].set_visible(False)
    ax_resource.tick_params(axis="y", length=0)
    for bar, value in zip(bars, values):
        ax_resource.text(bar.get_width() + 1.2, bar.get_y() + bar.get_height() / 2, str(value), va="center", fontsize=6)
    ax_resource.text(0.00, -0.11, "All 241 formal-v2 groups are retained; bars aggregate developmental-stage labels from the manifest.", transform=ax_resource.transAxes, fontsize=5.7, color=COLORS["muted"], wrap=True)
    add_panel_label(ax_resource, "a")
    source_rows.extend(
        [
            {"panel": "a", "row_type": "summary", "label": "verified_runs", "value": total_runs},
            {"panel": "a", "row_type": "summary", "label": "studies", "value": studies},
            {"panel": "a", "row_type": "summary", "label": "formal_tracks", "value": len(rows)},
        ]
    )
    source_rows.extend(
        {"panel": "a", "row_type": "stage", "label": stage, "value": count}
        for stage, count in stages
    )
    source_rows.extend(
        {
            "panel": "a", "row_type": "track", "label": row["group_id"], "value": int(row["n_runs"]),
            "stage": row["development_stage"] or "Not reported", "study": row["sra_study_accession"],
        }
        for row in rows
    )

    for axis in (ax_task, ax_model):
        axis.axis("off")
        axis.set_xlim(0, 1)
        axis.set_ylim(0, 1)

    ax_task.set_title("Prediction task and blocked development", loc="left", pad=8, fontweight="bold")
    draw_box(ax_task, (0.03, 0.32), 0.22, 0.33, "131,072-bp\nDNA window", COLORS["blue"])
    draw_box(ax_task, (0.38, 0.32), 0.24, 0.33, "241 RNA tracks\n1-bp and 128-bp", COLORS["teal"])
    draw_box(ax_task, (0.75, 0.32), 0.20, 0.33, "Gene/exon and\ncoverage scores", COLORS["orange"])
    draw_arrow(ax_task, (0.255, 0.485), (0.375, 0.485))
    draw_arrow(ax_task, (0.625, 0.485), (0.745, 0.485))
    ax_task.text(0.03, 0.84, "Development: five held-out genomic folds on chromosomes I-V", fontsize=7, color=COLORS["ink"], transform=ax_task.transAxes)
    for index, chromosome in enumerate(("I", "II", "III", "IV", "V")):
        x0 = 0.03 + index * 0.185
        draw_box(ax_task, (x0, 0.05), 0.13, 0.12, f"Fold {index + 1}\n{chromosome}", COLORS["muted"])
    ax_task.text(0.03, 0.21, "Three training seeds are nested within each fold.", fontsize=6.1, color=COLORS["muted"], transform=ax_task.transAxes)
    add_panel_label(ax_task, "b")
    source_rows.extend(
        [
            {"panel": "b", "row_type": "protocol", "label": "input_bp", "value": 131072},
            {"panel": "b", "row_type": "protocol", "label": "tracks", "value": 241},
            {"panel": "b", "row_type": "protocol", "label": "development_folds", "value": 5},
            {"panel": "b", "row_type": "protocol", "label": "seeds_per_fold", "value": 3},
        ]
    )

    ax_model.set_title("Species-adapted sequence-to-RNA model", loc="left", pad=8, fontweight="bold")
    draw_box(ax_model, (0.02, 0.35), 0.20, 0.35, "One-hot\nDNA", COLORS["blue"])
    draw_box(ax_model, (0.30, 0.35), 0.22, 0.35, "Frozen\nsequence trunk", COLORS["muted"])
    draw_box(ax_model, (0.60, 0.52), 0.17, 0.20, "Worm\nembedding", COLORS["teal"])
    draw_box(ax_model, (0.60, 0.20), 0.17, 0.20, "LoRA\nadapters", COLORS["orange"])
    draw_box(ax_model, (0.84, 0.35), 0.14, 0.35, "Dual RNA head\n1 bp | 128 bp", COLORS["red"])
    draw_arrow(ax_model, (0.225, 0.525), (0.295, 0.525))
    draw_arrow(ax_model, (0.525, 0.59), (0.595, 0.62))
    draw_arrow(ax_model, (0.525, 0.46), (0.595, 0.30))
    draw_arrow(ax_model, (0.775, 0.62), (0.835, 0.56))
    draw_arrow(ax_model, (0.775, 0.30), (0.835, 0.46))
    ax_model.text(0.02, 0.05, "The fixed 241-track output head defines the supported target space; it is not an unseen-condition decoder.", fontsize=5.9, color=COLORS["muted"], transform=ax_model.transAxes, wrap=True)
    add_panel_label(ax_model, "c")
    source_rows.extend(
        [
            {"panel": "c", "row_type": "architecture", "label": "sequence_input_bp", "value": 131072},
            {"panel": "c", "row_type": "architecture", "label": "output_tracks", "value": 241},
            {"panel": "c", "row_type": "architecture", "label": "output_resolutions_bp", "value": "1;128"},
        ]
    )
    return figure, source_rows


def write_source_data(path: Path, rows: list[dict[str, Any]]) -> None:
    fields = ("panel", "row_type", "label", "value", "stage", "study")
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def save_outputs(figure: plt.Figure, output_dir: Path) -> dict[str, Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    outputs = {
        "pdf": output_dir / f"{FIGURE_STEM}.pdf",
        "svg": output_dir / f"{FIGURE_STEM}.svg",
        "png": output_dir / f"{FIGURE_STEM}.png",
        "tiff": output_dir / f"{FIGURE_STEM}.tiff",
    }
    figure.savefig(outputs["pdf"], bbox_inches="tight")
    figure.savefig(outputs["svg"], bbox_inches="tight")
    figure.savefig(outputs["png"], dpi=RASTER_DPI, bbox_inches="tight")
    figure.savefig(outputs["tiff"], dpi=RASTER_DPI, bbox_inches="tight")
    return outputs


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    rows = read_manifest(args.manifest)
    figure, source_rows = render(rows)
    outputs = save_outputs(figure, args.output_dir)
    plt.close(figure)
    source_path = args.output_dir / "source_data.tsv"
    write_source_data(source_path, source_rows)
    image = Image.open(outputs["png"])
    if image.width < 3000 or image.height < 1500:
        raise RuntimeError(f"Unexpectedly small raster output: {image.size}")
    manifest = {
        "figure": FIGURE_STEM,
        "created_at_utc": datetime.now(timezone.utc).replace(microsecond=0).isoformat(),
        "input_manifest": str(args.manifest.relative_to(REPO_ROOT)),
        "input_manifest_sha256": sha256(args.manifest),
        "formal_groups": len(rows),
        "verified_runs": sum(int(row["n_runs"]) for row in rows),
        "artifacts": {key: {"path": str(path.relative_to(REPO_ROOT)), "sha256": sha256(path), "bytes": path.stat().st_size} for key, path in outputs.items()},
        "source_data": {"path": str(source_path.relative_to(REPO_ROOT)), "sha256": sha256(source_path), "rows": len(source_rows)},
    }
    (args.output_dir / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
