#!/usr/bin/env python3
"""Render and audit the P10 DPY-27 internal-validation figure from Source Data."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image


REPO_ROOT = Path(__file__).resolve().parents[1]
SPEC_PATH = REPO_ROOT / "alphagenome_custom/metadata/v2/p10_dpy27_internal_application_spec.json"
FIG_WIDTH_MM = 183
FIG_HEIGHT_MM = 82
RASTER_DPI = 600
FONT_SIZE_PT = 7
REQUIRED_EXPORT_SUFFIXES = (".pdf", ".svg", ".tiff", ".png")
COLORS = {
    "autosome": "#595959",
    "x": "#C43C39",
    "observed": "#007C83",
    "predicted": "#C43C39",
    "line": "#333333",
    "grid": "#D7D7D7",
}

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


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(8 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def read_json(path: Path) -> dict[str, object]:
    payload = json.loads(path.read_text())
    if not isinstance(payload, dict):
        raise ValueError(f"Expected JSON object: {path}")
    return payload


def read_source_data(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    required = {
        "panel",
        "fold",
        "seed",
        "gene_id",
        "chromosome",
        "endpoint",
        "observed",
        "predicted",
        "estimate",
        "ci_95_low",
        "ci_95_high",
    }
    if not rows or required - set(rows[0]):
        raise ValueError("Figure Source Data is empty or has an incompatible schema")
    return rows


def validate_source_data_contract(
    rows: list[dict[str, str]], spec: dict[str, object]
) -> None:
    contract = spec["figure_source_data_contract"]
    expected_panels = set(contract["panels"])
    observed_panels = {row["panel"] for row in rows}
    if observed_panels != expected_panels:
        raise ValueError(
            f"Source Data panels differ from the frozen contract: {observed_panels}"
        )
    endpoints = contract["panel_endpoints"]
    for panel in contract["panels"]:
        observed = {row["endpoint"] for row in rows if row["panel"] == panel}
        if observed != set(endpoints[panel]):
            raise ValueError(f"Source Data endpoint set changed for panel {panel}")
    expected_counts = {
        "b_fold_x_minus_autosomes": 10,
        "c_fold_agreement": 10,
        "c_fold_primary_summary": 4,
    }
    for panel, expected_count in expected_counts.items():
        actual = sum(row["panel"] == panel for row in rows)
        if actual != expected_count:
            raise ValueError(
                f"Source Data row count changed for panel {panel}: {actual}"
            )
    gene_folds = {
        int(row["fold"])
        for row in rows
        if row["panel"] == "a_gene_contrasts"
    }
    if gene_folds != {1, 2, 3, 4, 5}:
        raise ValueError("Panel a must include gene rows from all five folds")


def _numeric(value: str) -> float:
    result = float(value)
    if not np.isfinite(result):
        raise ValueError("Figure Source Data contains a non-finite plotted value")
    return result


def _summary(rows: list[dict[str, str]], endpoint: str) -> tuple[float, float, float]:
    matches = [
        row
        for row in rows
        if row["panel"] == "c_fold_primary_summary" and row["endpoint"] == endpoint
    ]
    if len(matches) != 1:
        raise ValueError(f"Expected one summary row for {endpoint}")
    row = matches[0]
    return _numeric(row["estimate"]), _numeric(row["ci_95_low"]), _numeric(row["ci_95_high"])


def _panel_label(axis: plt.Axes, label: str) -> None:
    axis.text(
        -0.16,
        1.05,
        label,
        transform=axis.transAxes,
        fontsize=9,
        fontweight="bold",
        va="top",
    )


def build_figure(rows: list[dict[str, str]]) -> plt.Figure:
    gene_rows = [row for row in rows if row["panel"] == "a_gene_contrasts"]
    fold_rows = [row for row in rows if row["panel"] == "b_fold_x_minus_autosomes"]
    if not gene_rows or len(fold_rows) != 10:
        raise ValueError("Figure Source Data lacks the expected gene or fold rows")
    fig, axes = plt.subplots(
        1,
        3,
        figsize=(FIG_WIDTH_MM / 25.4, FIG_HEIGHT_MM / 25.4),
        gridspec_kw={"width_ratios": [1.25, 0.9, 1.0], "wspace": 0.48},
    )

    ax = axes[0]
    for chromosome, color, label, zorder in (
        ("autosome", COLORS["autosome"], "Autosomes", 2),
        ("X", COLORS["x"], "Chromosome X", 3),
    ):
        selected = [
            row
            for row in gene_rows
            if (row["chromosome"] == "X") == (chromosome == "X")
        ]
        x = np.array([_numeric(row["observed"]) for row in selected])
        y = np.array([_numeric(row["predicted"]) for row in selected])
        ax.scatter(
            x,
            y,
            s=7 if chromosome == "X" else 4,
            color=color,
            alpha=0.68 if chromosome == "X" else 0.32,
            edgecolors="none",
            label=label,
            rasterized=True,
            zorder=zorder,
        )
    all_values = np.array(
        [_numeric(row[key]) for row in gene_rows for key in ("observed", "predicted")]
    )
    low, high = float(all_values.min()), float(all_values.max())
    padding = max((high - low) * 0.05, 0.02)
    ax.plot(
        [low - padding, high + padding],
        [low - padding, high + padding],
        color=COLORS["line"],
        linewidth=0.7,
        linestyle="--",
        label="Identity",
        zorder=1,
    )
    ax.set_xlim(low - padding, high + padding)
    ax.set_ylim(low - padding, high + padding)
    ax.set_xlabel("Observed gene contrast")
    ax.set_ylabel("Predicted gene contrast")
    ax.set_title("Gene-level response")
    ax.legend(frameon=False, loc="best", handletextpad=0.4)
    _panel_label(ax, "a")

    ax = axes[1]
    by_fold: dict[int, dict[str, float]] = {}
    for row in fold_rows:
        fold = int(row["fold"])
        by_fold.setdefault(fold, {})[row["endpoint"].split("_", 1)[0]] = _numeric(
            row["estimate"]
        )
    if sorted(by_fold) != [1, 2, 3, 4, 5] or any(set(value) != {"observed", "predicted"} for value in by_fold.values()):
        raise ValueError("Panel b requires observed and predicted values for folds 1-5")
    for fold, values in sorted(by_fold.items()):
        ax.plot(
            [0, 1],
            [values["observed"], values["predicted"]],
            color="#A8A8A8",
            linewidth=0.7,
            zorder=1,
        )
        ax.scatter(0, values["observed"], s=18, color=COLORS["observed"], zorder=2)
        ax.scatter(1, values["predicted"], s=18, color=COLORS["predicted"], zorder=2)
    for x_position, endpoint, color in (
        (0, "observed_median_contrast_x_minus_autosomes", COLORS["observed"]),
        (1, "predicted_median_contrast_x_minus_autosomes", COLORS["predicted"]),
    ):
        estimate, lower, upper = _summary(rows, endpoint)
        ax.errorbar(
            x_position,
            estimate,
            yerr=np.array([[estimate - lower], [upper - estimate]]),
            color=color,
            marker="D",
            markersize=4,
            markeredgecolor="white",
            markeredgewidth=0.5,
            capsize=2,
            linewidth=1.1,
            zorder=4,
        )
    ax.axhline(0, color=COLORS["line"], linewidth=0.6, linestyle=":")
    ax.set_xticks([0, 1], ["Observed", "Predicted"])
    ax.set_xlim(-0.45, 1.45)
    ax.set_ylabel("Median contrast: X - autosomes")
    ax.set_title("Chromosome-level recovery")
    _panel_label(ax, "b")

    ax = axes[2]
    endpoint_display = (
        ("x_predicted_observed_direction_concordance", "X direction\nconcordance"),
        ("predicted_observed_gene_spearman", "Gene-level\nSpearman"),
    )
    endpoint_positions = np.arange(len(endpoint_display), dtype=float)
    for endpoint_index, (endpoint, _) in enumerate(endpoint_display):
        values = np.array([_numeric(str(by_fold_row[endpoint])) for by_fold_row in _fold_endpoint_rows(rows)])
        offsets = np.linspace(-0.09, 0.09, len(values))
        ax.scatter(
            endpoint_positions[endpoint_index] + offsets,
            values,
            s=15,
            color=COLORS["autosome"],
            alpha=0.7,
            zorder=2,
        )
        estimate, lower, upper = _summary(rows, endpoint)
        ax.errorbar(
            endpoint_positions[endpoint_index],
            estimate,
            yerr=np.array([[estimate - lower], [upper - estimate]]),
            color=COLORS["predicted"],
            marker="D",
            markersize=4,
            markeredgecolor="white",
            markeredgewidth=0.5,
            capsize=2,
            linewidth=1.1,
            zorder=4,
        )
    ax.axhline(0, color=COLORS["line"], linewidth=0.6, linestyle=":")
    ax.set_xticks(endpoint_positions, [display for _, display in endpoint_display])
    ax.set_xlim(-0.5, len(endpoint_display) - 0.5)
    ax.set_ylim(-1.05, 1.05)
    ax.set_ylabel("Fold-level estimate")
    ax.set_title("Agreement endpoints")
    _panel_label(ax, "c")

    for axis in axes:
        axis.spines["top"].set_visible(False)
        axis.spines["right"].set_visible(False)
        axis.grid(axis="y", color=COLORS["grid"], linewidth=0.45, alpha=0.55)
        axis.set_axisbelow(True)
    fig.suptitle(
        "DPY-27 response recovery in held-out genomic blocks",
        fontsize=9,
        fontweight="bold",
        y=1.0,
    )
    fig.text(
        0.5,
        0.005,
        "Internal genomic-block validation; both conditions were training targets.",
        ha="center",
        va="bottom",
        fontsize=6,
        color=COLORS["autosome"],
    )
    fig.subplots_adjust(left=0.075, right=0.99, bottom=0.23, top=0.86)
    return fig


def _fold_endpoint_rows(source_rows: list[dict[str, str]]) -> list[dict[str, float]]:
    endpoints = {
        "x_predicted_observed_direction_concordance",
        "predicted_observed_gene_spearman",
    }
    result: list[dict[str, float]] = []
    for fold in range(1, 6):
        matches = [
            row
            for row in source_rows
            if row["panel"] == "c_fold_agreement" and int(row["fold"]) == fold
        ]
        values = {row["endpoint"]: _numeric(row["estimate"]) for row in matches}
        if set(values) != endpoints:
            raise ValueError(f"Panel c fold {fold} lacks the exact agreement endpoints")
        result.append(values)
    return result


def _run_json(command: list[str]) -> dict[str, object]:
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


def _image_audit(path: Path) -> dict[str, object]:
    with Image.open(path) as image:
        grayscale = np.asarray(image.convert("L"), dtype=np.uint8)
        width, height = image.size
        dpi = image.info.get("dpi", (None, None))
    nonwhite_fraction = float(np.mean(grayscale < 250))
    standard_deviation = float(grayscale.std())
    result = {
        "width_pixels": width,
        "height_pixels": height,
        "dpi": list(dpi),
        "nonwhite_fraction": nonwhite_fraction,
        "grayscale_standard_deviation": standard_deviation,
        "passed": (
            width >= 3000
            and height >= 1500
            and nonwhite_fraction >= 0.01
            and standard_deviation >= 5.0
        ),
    }
    return result


def render_and_audit(spec_path: Path, expected_spec_sha256: str) -> dict[str, object]:
    if sha256(spec_path) != expected_spec_sha256:
        raise RuntimeError("P10 spec changed before figure rendering")
    spec = read_json(spec_path)
    outputs = spec["output_paths"]
    runtime = spec["plot_runtime"]
    source_path = REPO_ROOT / outputs["figure_source_data"]
    paths = {
        "pdf": REPO_ROOT / outputs["figure_pdf"],
        "svg": REPO_ROOT / outputs["figure_svg"],
        "tiff": REPO_ROOT / outputs["figure_tiff"],
        "png": REPO_ROOT / outputs["figure_png"],
    }
    for path in paths.values():
        path.parent.mkdir(parents=True, exist_ok=True)
    validator_path = Path(runtime["validate_figure_path"])
    pdf_audit_path = Path(runtime["audit_pdf_text_path"])
    if sha256(validator_path) != runtime["validate_figure_sha256"]:
        raise RuntimeError("Figure validator changed after P10 spec freeze")
    if sha256(pdf_audit_path) != runtime["audit_pdf_text_sha256"]:
        raise RuntimeError("PDF text auditor changed after P10 spec freeze")
    static = _run_json(
        [sys.executable, str(validator_path), str(Path(__file__).resolve()), "--backend", "python", "--json"]
    )
    if not static.get("summary", {}).get("ready"):
        raise RuntimeError("Figure static preflight did not pass")
    rows = read_source_data(source_path)
    validate_source_data_contract(rows, spec)
    figure = build_figure(rows)
    try:
        figure.savefig(paths["pdf"], format="pdf", bbox_inches="tight")
        figure.savefig(paths["svg"], format="svg", bbox_inches="tight")
        figure.savefig(paths["tiff"], format="tiff", dpi=RASTER_DPI, bbox_inches="tight")
        figure.savefig(paths["png"], format="png", dpi=RASTER_DPI, bbox_inches="tight")
    finally:
        plt.close(figure)
    pdf_text = _run_json(
        [sys.executable, str(pdf_audit_path), str(paths["pdf"]), "--min-pt", "5", "--json"]
    )
    if not pdf_text.get("auditable") or pdf_text.get("below_minimum_count") != 0:
        raise RuntimeError("PDF text-size audit did not pass")
    png_audit = _image_audit(paths["png"])
    tiff_audit = _image_audit(paths["tiff"])
    if not png_audit["passed"] or not tiff_audit["passed"]:
        raise RuntimeError("Raster pixel audit did not pass")
    qa = {
        "schema_version": 1,
        "status": "passed",
        "spec_sha256": expected_spec_sha256,
        "source_data_path": str(source_path.relative_to(REPO_ROOT)),
        "source_data_sha256": sha256(source_path),
        "plot_runtime": {
            "python_executable": sys.executable,
            "python_version": sys.version.split()[0],
            "numpy_version": np.__version__,
            "matplotlib_version": mpl.__version__,
            "pillow_version": Image.__version__,
        },
        "static_preflight": static,
        "pdf_text_audit": pdf_text,
        "png_pixel_audit": png_audit,
        "tiff_pixel_audit": tiff_audit,
        "artifacts": {
            key: {
                "path": str(path.relative_to(REPO_ROOT)),
                "sha256": sha256(path),
                "size_bytes": path.stat().st_size,
            }
            for key, path in paths.items()
        },
    }
    qa_path = REPO_ROOT / outputs["figure_qa"]
    temporary = qa_path.with_suffix(qa_path.suffix + ".tmp")
    temporary.write_text(json.dumps(qa, indent=2, sort_keys=True) + "\n")
    temporary.replace(qa_path)
    return qa


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--spec", type=Path, default=SPEC_PATH)
    parser.add_argument("--spec-sha256", required=True)
    return parser


def main(argv: list[str] | None = None) -> None:
    args = build_parser().parse_args(argv)
    result = render_and_audit(args.spec.resolve(), args.spec_sha256)
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
