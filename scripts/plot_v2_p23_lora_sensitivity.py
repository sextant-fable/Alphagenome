#!/usr/bin/env python3
"""Render the P23 LoRA rank/location sensitivity pilot as Extended Data 5."""

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
P23_SOURCE = REPO_ROOT / "results/v2_p23_lora_sensitivity_pilot/p23_lora_sensitivity_pilot.tsv"
P17_RECORDS = REPO_ROOT / "results/v2_p17_iv_controlled_matrix/p17_iv_controlled_records.tsv"
P23_REVIEW = REPO_ROOT / "alphagenome_custom/metadata/v2/audits/P23/review.json"
DEFAULT_OUTPUT = REPO_ROOT / "results/v2_p23_manuscript_figure_20260824T000000Z"

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
    "font.size": 7,
    "axes.labelsize": 7,
    "axes.titlesize": 8,
    "xtick.labelsize": 6.2,
    "ytick.labelsize": 6.2,
    "axes.linewidth": 0.7,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
})


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
    axis.grid(axis="y", color="#D9D9D9", linewidth=0.45)
    axis.set_axisbelow(True)


def panel_label(axis: plt.Axes, label: str) -> None:
    axis.text(-0.15, 1.08, label, transform=axis.transAxes, fontsize=9,
              fontweight="bold", va="top")


def p17_reference() -> dict[str, str]:
    rows = read_tsv(P17_RECORDS)
    for row in rows:
        if (
            row.get("configuration") == "B_iv_dual"
            and row.get("fold") == "1"
            and row.get("seed") == "20260714"
        ):
            return {
                "configuration": "P17 B reference",
                "fold": row["fold"],
                "seed": row["seed"],
                "physical_gpu": row.get("physical_gpu", ""),
                "primary_biological_score": row["primary_biological_score"],
                "gene_exon_coverage_pearson_log1p": row["gene_exon_coverage_pearson_log1p"],
                "per_track_pearson_128bp_log1p": row["per_track_pearson_128bp_log1p"],
                "trainable_parameters": row["trainable_parameters"],
                "evaluated_chromosomes": row["evaluated_chromosomes"],
                "locked_test_block_signal_reads": row["locked_test_block_signal_reads"],
                "source_class": "P17 same fold/seed reference",
            }
    raise RuntimeError("P17 fold-1 seed-20260714 Model B reference is missing")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    output = args.output.resolve()
    if output.exists():
        raise FileExistsError(f"P23 figure output is immutable: {output}")
    if not P23_SOURCE.is_file() or not P17_RECORDS.is_file() or not P23_REVIEW.is_file():
        raise FileNotFoundError("P23 source, P17 reference, or review is missing")
    review = json.loads(P23_REVIEW.read_text())
    if review.get("status") != "PASS":
        raise RuntimeError("P23 figure requires R23 PASS")

    p23_rows = read_tsv(P23_SOURCE)
    expected = {"B_rank4_final_block", "B_rank8_previous_block", "B_rank16_final_block"}
    observed = {row["configuration"] for row in p23_rows}
    if observed != expected or len(p23_rows) != 3:
        raise RuntimeError(f"Unexpected P23 configurations: {sorted(observed)}")
    reference = p17_reference()
    rows = p23_rows + [reference]
    labels = {
        "B_rank4_final_block": "rank 4\nfinal",
        "B_rank8_previous_block": "rank 8\nprevious",
        "B_rank16_final_block": "rank 16\nfinal",
        "P17 B reference": "P17 B\nreference",
    }
    colors = {
        "B_rank4_final_block": "#5A8F9E",
        "B_rank8_previous_block": "#B07A3C",
        "B_rank16_final_block": "#B23A48",
        "P17 B reference": "#707070",
    }
    order = ["B_rank4_final_block", "B_rank8_previous_block", "B_rank16_final_block", "P17 B reference"]
    by_config = {row["configuration"]: row for row in rows}
    x = np.arange(len(order))

    # 183 x 76 mm double-column figure size expressed in inches for matplotlib.
    fig, axes = plt.subplots(1, 3, figsize=(7.2047, 2.9921),
                             gridspec_kw={"width_ratios": [1.05, 1.15, 0.9]})
    axis = axes[0]
    panel_label(axis, "a")
    axis.set_title("Primary score", loc="left", fontweight="bold")
    values = np.array([float(by_config[key]["primary_biological_score"]) for key in order])
    axis.plot(x, values, color="#A9A9A9", linewidth=0.8, zorder=1)
    axis.scatter(x, values, s=40, c=[colors[key] for key in order], edgecolor="white",
                 linewidth=0.5, zorder=2)
    for index, value in enumerate(values):
        axis.text(index, value + 0.0022, f"{value:.3f}", ha="center", va="bottom", fontsize=5.5)
    setup(axis)
    axis.set_xticks(x, [labels[key] for key in order])
    axis.set_ylabel("Primary biological score")
    axis.set_ylim(0.605, 0.625)

    axis = axes[1]
    panel_label(axis, "b")
    axis.set_title("Endpoint components", loc="left", fontweight="bold")
    width = 0.34
    gene = np.array([float(by_config[key]["gene_exon_coverage_pearson_log1p"]) for key in order])
    bp = np.array([float(by_config[key]["per_track_pearson_128bp_log1p"]) for key in order])
    axis.bar(x - width / 2, gene, width, color="#C5D8DA", label="Gene-exon")
    axis.bar(x + width / 2, bp, width, color="#B23A48", label="Fine-resolution (128 bp)")
    setup(axis)
    axis.set_xticks(x, [labels[key] for key in order])
    axis.set_ylabel("Pearson correlation (log1p)")
    axis.set_ylim(0.56, 0.67)
    axis.legend(frameon=False, fontsize=5.5, loc="lower left")

    axis = axes[2]
    panel_label(axis, "c")
    axis.set_title("Trainable parameters", loc="left", fontweight="bold")
    params = np.array([int(by_config[key]["trainable_parameters"]) for key in order]) / 1e6
    axis.bar(x, params, color=[colors[key] for key in order], width=0.58)
    for index, value in enumerate(params):
        axis.text(index, value + 0.025, f"{value:.2f}M", ha="center", va="bottom", fontsize=5.5)
    setup(axis)
    axis.set_xticks(x, [labels[key] for key in order])
    axis.set_ylabel("Trainable parameters (millions)")
    axis.set_ylim(0, 1.65)

    fig.suptitle("LoRA rank and insertion-location sensitivity pilot", x=0.01, ha="left",
                 fontsize=9, fontweight="bold")
    fig.text(0.01, 0.015,
             "Single I-V fold/seed pilot; P17 B is the matched fold/seed reference. "
             "No inferential interval is shown.",
             fontsize=5.6, color="#555555", ha="left")
    fig.tight_layout(rect=(0, 0.10, 1, 0.93), w_pad=1.8)
    output.mkdir(parents=True)
    paths = {
        "pdf": output / "figure_ed5_lora_sensitivity.pdf",
        "svg": output / "figure_ed5_lora_sensitivity.svg",
        "png": output / "figure_ed5_lora_sensitivity.png",
        "tiff": output / "figure_ed5_lora_sensitivity.tiff",
        "source_data": output / "ed5_lora_sensitivity_source_data.tsv",
    }
    with paths["source_data"].open("w", newline="", encoding="utf-8") as handle:
        fields = list(rows[0].keys())
        handle.write("\t".join(fields) + "\n")
        for row in rows:
            handle.write("\t".join(str(row.get(field, "")) for field in fields) + "\n")
    fig.savefig(paths["pdf"], bbox_inches="tight", facecolor="white")
    fig.savefig(paths["svg"], bbox_inches="tight", facecolor="white")
    fig.savefig(paths["png"], dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(paths["tiff"], dpi=600, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    with Image.open(paths["png"]) as image:
        if image.width < 1200 or image.height < 350:
            raise RuntimeError("P23 figure raster export is unexpectedly small")
    qa = {
        "schema_version": 1,
        "backend": "python/matplotlib",
        "figure_contract": {
            "core_conclusion": "The registered P23 pilot does not show a large primary-score change across the tested LoRA ranks and insertion locations on one I-V fold/seed.",
            "archetype": "quantitative grid",
            "hero_evidence": "Primary scores for the three P23 configurations with a same-fold/seed P17 Model B reference.",
            "boundary": "This is a single fold/seed pilot, not a full five-fold or three-seed sensitivity matrix and not a biological mechanism test.",
            "statistics": "No inferential interval is shown; points and bars are descriptive for one registered fold/seed.",
        },
        "input_source_data": {
            "p23": str(P23_SOURCE.relative_to(REPO_ROOT)),
            "p23_sha256": sha256(P23_SOURCE),
            "p17_reference": str(P17_RECORDS.relative_to(REPO_ROOT)),
            "p17_reference_sha256": sha256(P17_RECORDS),
            "review": str(P23_REVIEW.relative_to(REPO_ROOT)),
            "review_sha256": sha256(P23_REVIEW),
        },
        "output": {key: str(path.relative_to(REPO_ROOT)) for key, path in paths.items()},
        "exclusions": "None; all three P23 configurations and the matched P17 reference are shown.",
        "image_integrity": "No biological raster images; all marks are rendered directly from tabular run records.",
    }
    (output / "figure_ed5_qa.json").write_text(json.dumps(qa, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"status": "rendered", "output": str(output.relative_to(REPO_ROOT))}, indent=2))


if __name__ == "__main__":
    main()
