#!/usr/bin/env python3
"""Render observed-versus-predicted P12 coverage examples from Source Data."""

from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import subprocess
import sys
from typing import Any

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image


REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = REPO_ROOT / "results/v2_p12_replicate_holdout_figure/coverage_examples"
SKILL_ROOT = Path("/home/zelinli6/.codex/skills/nature-figure")
WIDTH_MM = 183
RASTER_DPI = 600
COLORS = {"observed": "#4D4D4D", "predicted": "#B23A48", "grid": "#D8D8D8", "ink": "#2F2F2F"}

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


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise ValueError(f"Empty coverage Source Data: {path}")
    return rows


def validate_rows(rows: list[dict[str, str]]) -> dict[str, dict[str, Any]]:
    required = {"group_id", "selection_quantile", "heldout_pearson_128bp", "chromosome", "interval_start", "interval_end", "bin_index", "genomic_start", "observed_128bp", "predicted_128bp"}
    if not required.issubset(rows[0]):
        raise ValueError(f"Coverage Source Data is missing columns: {sorted(required - set(rows[0]))}")
    groups: dict[str, list[dict[str, str]]] = {}
    for row in rows:
        groups.setdefault(row["group_id"], []).append(row)
    if len(groups) != 3:
        raise ValueError(f"Expected three deterministic coverage examples, found {len(groups)}")
    by_quantile: dict[str, dict[str, Any]] = {}
    for group_id, group_rows in groups.items():
        if len(group_rows) != 1024:
            raise ValueError(f"Expected 1024 128-bp bins for {group_id}, found {len(group_rows)}")
        group_rows.sort(key=lambda row: int(row["bin_index"]))
        if [int(row["bin_index"]) for row in group_rows] != list(range(1024)):
            raise ValueError(f"Bin indices are not contiguous for {group_id}")
        observed = np.asarray([float(row["observed_128bp"]) for row in group_rows], dtype=float)
        predicted = np.asarray([float(row["predicted_128bp"]) for row in group_rows], dtype=float)
        if not np.isfinite(observed).all() or not np.isfinite(predicted).all() or (observed < 0).any() or (predicted < 0).any():
            raise ValueError(f"Non-finite or negative coverage values for {group_id}")
        corr = float(np.corrcoef(observed, predicted)[0, 1])
        quantile = group_rows[0]["selection_quantile"]
        if quantile in by_quantile:
            raise ValueError(f"Duplicate selection quantile: {quantile}")
        by_quantile[quantile] = {
            "group_id": group_id,
            "rows": group_rows,
            "observed": observed,
            "predicted": predicted,
            "correlation": corr,
            "chromosome": group_rows[0]["chromosome"],
            "interval_start": int(group_rows[0]["interval_start"]),
            "interval_end": int(group_rows[0]["interval_end"]),
            "study": group_rows[0]["sra_study_accession"],
            "stage": group_rows[0]["development_stage"],
            "tissue": group_rows[0]["tissue_or_cell_type"],
        }
    if set(by_quantile) != {"lower", "median", "upper"}:
        raise ValueError(f"Unexpected quantile labels: {sorted(by_quantile)}")
    return by_quantile


def export_figure(fig: plt.Figure, output: Path) -> dict[str, str]:
    output.mkdir(parents=True, exist_ok=True)
    fig.set_size_inches(WIDTH_MM / 25.4, 104 / 25.4)
    paths = {key: output / f"figure_p12_coverage_examples.{key}" for key in ("pdf", "svg", "png", "tiff")}
    fig.savefig(paths["pdf"], bbox_inches="tight", facecolor="white")
    fig.savefig(paths["svg"], bbox_inches="tight", facecolor="white")
    fig.savefig(paths["png"], dpi=RASTER_DPI, bbox_inches="tight", facecolor="white")
    fig.savefig(paths["tiff"], dpi=RASTER_DPI, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return {key: str(path.relative_to(REPO_ROOT)) for key, path in paths.items()}


def make_figure(examples: dict[str, dict[str, Any]], output: Path) -> dict[str, str]:
    order = ("lower", "median", "upper")
    fig, axes = plt.subplots(3, 1, figsize=(WIDTH_MM / 25.4, 104 / 25.4), sharex=True, gridspec_kw={"hspace": 0.28})
    for index, (axis, quantile) in enumerate(zip(axes, order)):
        item = examples[quantile]
        observed = item["observed"]
        predicted = item["predicted"]
        x = np.arange(1024) * 0.128
        axis.spines["top"].set_visible(False)
        axis.spines["right"].set_visible(False)
        axis.grid(axis="y", color=COLORS["grid"], linewidth=0.45, zorder=0)
        axis.set_axisbelow(True)
        axis.plot(x, np.log1p(observed), color=COLORS["observed"], linewidth=0.75, label="Observed" if index == 0 else None)
        axis.plot(x, np.log1p(predicted), color=COLORS["predicted"], linewidth=0.75, label="Predicted" if index == 0 else None)
        axis.set_ylabel("log1p coverage")
        axis.set_title(f"{quantile.capitalize()} held-out track: {item['stage']} | {item['tissue']} | {item['chromosome']}:{item['interval_start']:,}-{item['interval_end']:,} | Pearson = {item['correlation']:.3f}", loc="left", fontsize=7.2, fontweight="bold")
        axis.text(-0.055, 1.02, chr(ord("a") + index), transform=axis.transAxes, fontsize=9, fontweight="bold", va="bottom")
    axes[-1].set_xlabel("Genomic position within selected 131-kb validation window (kb)")
    axes[0].legend(frameon=False, loc="upper right", ncol=2, handlelength=2.2, columnspacing=0.8)
    return export_figure(fig, output)


def run_qa(output: Path, source_script: Path) -> dict[str, Any]:
    validator = SKILL_ROOT / "scripts/validate_figure.py"
    auditor = SKILL_ROOT / "scripts/audit_pdf_text.py"
    qa: dict[str, Any] = {"status": "passed", "backend": "python_matplotlib", "created_at": datetime.now(timezone.utc).replace(microsecond=0).isoformat(), "figures": {}}
    pdf = output / "figure_p12_coverage_examples.pdf"
    validation = subprocess.run([sys.executable, str(validator), str(source_script), "--backend", "python", "--json"], capture_output=True, text=True, check=False)
    audit = subprocess.run([sys.executable, str(auditor), str(pdf), "--min-pt", "5", "--json"], capture_output=True, text=True, check=False)
    try:
        validation_json = json.loads(validation.stdout)
    except json.JSONDecodeError:
        validation_json = {"status": "invalid", "stdout": validation.stdout, "stderr": validation.stderr}
    try:
        audit_json = json.loads(audit.stdout)
    except json.JSONDecodeError:
        audit_json = {"status": "invalid", "stdout": audit.stdout, "stderr": audit.stderr}
    raster = {}
    for suffix in ("png", "tiff"):
        path = output / f"figure_p12_coverage_examples.{suffix}"
        with Image.open(path) as image:
            raster[suffix] = {"size": list(image.size), "dpi": [float(v) if v is not None else None for v in image.info.get("dpi", (None, None))], "mode": image.mode}
    if validation.returncode or audit.returncode:
        qa["status"] = "failed"
    qa["figures"][pdf.stem] = {
        "pdf": str(pdf.relative_to(REPO_ROOT)),
        "pdf_sha256": sha256(pdf),
        "svg_sha256": sha256(output / "figure_p12_coverage_examples.svg"),
        "png_sha256": sha256(output / "figure_p12_coverage_examples.png"),
        "tiff_sha256": sha256(output / "figure_p12_coverage_examples.tiff"),
        "source_validation": validation_json,
        "pdf_text_audit": audit_json,
        "raster": raster,
    }
    return qa


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", default=str(DEFAULT_OUTPUT))
    args = parser.parse_args()
    output = Path(args.output_dir).resolve()
    source = output / "coverage_examples_source_data.tsv"
    manifest_path = output / "coverage_examples_manifest.json"
    rows = read_rows(source)
    examples = validate_rows(rows)
    paths = make_figure(examples, output)
    qa = run_qa(output, Path(__file__).resolve())
    qa_path = output / "coverage_figure_qa.json"
    qa_path.write_text(json.dumps(qa, indent=2, sort_keys=True) + "\n")
    metadata = json.loads(manifest_path.read_text())
    metadata["figure"] = paths
    metadata["figure_qa"] = str(qa_path.relative_to(REPO_ROOT))
    metadata["figure_qa_sha256"] = sha256(qa_path)
    manifest_path.write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n")
    main_manifest_path = output.parent / "manifest.json"
    if main_manifest_path.is_file():
        main_manifest = json.loads(main_manifest_path.read_text())
        main_manifest["coverage_examples"] = metadata
        main_manifest.setdefault("figures", {})["coverage"] = paths
        main_manifest.setdefault("source_data", {})["coverage"] = str(source.relative_to(REPO_ROOT))
        main_manifest_path.write_text(json.dumps(main_manifest, indent=2, sort_keys=True) + "\n")
        main_qa_path = output.parent / "figure_qa.json"
        if main_qa_path.is_file():
            main_qa = json.loads(main_qa_path.read_text())
            main_qa.setdefault("figures", {}).update(qa.get("figures", {}))
            main_qa["coverage_qa"] = str(qa_path.relative_to(REPO_ROOT))
            main_qa["manifest_sha256"] = sha256(main_manifest_path)
            main_qa["status"] = "passed" if main_qa.get("status") == "passed" and qa.get("status") == "passed" else "failed"
            main_qa_path.write_text(json.dumps(main_qa, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"figure": paths, "qa": qa, "source_data": str(source.relative_to(REPO_ROOT))}, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
