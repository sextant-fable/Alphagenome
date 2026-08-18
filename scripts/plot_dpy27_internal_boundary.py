#!/usr/bin/env python3
"""Render a neutral-title manuscript display from the audited P10 Source Data."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path

import matplotlib.pyplot as plt
from PIL import Image

from scripts import plot_dpy27_internal_application as p10_plot


REPO_ROOT = Path(__file__).resolve().parents[1]
SOURCE_DATA = REPO_ROOT / "runs/v2_p10_dpy27_internal_application_20260813/figure_source_data.tsv"
REVIEW = REPO_ROOT / "alphagenome_custom/metadata/v2/audits/P10/review.json"
DEFAULT_OUTPUT_DIR = REPO_ROOT / "results/v2_p10_dpy27_boundary_figure"
STEM = "figure_dpy27_internal_boundary"


def render(output_dir: Path) -> dict[str, object]:
    rows = p10_plot.read_source_data(SOURCE_DATA)
    spec = p10_plot.read_json(p10_plot.SPEC_PATH)
    p10_plot.validate_source_data_contract(rows, spec)
    review = p10_plot.read_json(REVIEW)
    if review.get("phase") != "P10" or review.get("status") != "PASS":
        raise RuntimeError("The derived boundary figure requires the formal P10 R10 PASS")
    figure = p10_plot.build_figure(rows)
    try:
        figure.axes[0].set_title("Gene-level contrast")
        figure.axes[1].set_title("Chromosome-level contrast")
        if figure._suptitle is None:
            raise RuntimeError("Expected a figure title from the audited source renderer")
        figure._suptitle.set_text("DPY-27 contrast in held-out genomic blocks")
        output_dir.mkdir(parents=True, exist_ok=True)
        paths = {
            "pdf": output_dir / f"{STEM}.pdf",
            "svg": output_dir / f"{STEM}.svg",
            "png": output_dir / f"{STEM}.png",
            "tiff": output_dir / f"{STEM}.tiff",
        }
        figure.savefig(paths["pdf"], bbox_inches="tight")
        figure.savefig(paths["svg"], bbox_inches="tight")
        figure.savefig(paths["png"], dpi=p10_plot.RASTER_DPI, bbox_inches="tight")
        figure.savefig(paths["tiff"], dpi=p10_plot.RASTER_DPI, bbox_inches="tight")
    finally:
        plt.close(figure)
    with Image.open(paths["png"]) as image:
        if image.width < 3000 or image.height < 1500:
            raise RuntimeError(f"Unexpectedly small raster output: {image.size}")
    manifest = {
        "figure": STEM,
        "created_at_utc": datetime.now(timezone.utc).replace(microsecond=0).isoformat(),
        "source_data": {
            "path": str(SOURCE_DATA.relative_to(REPO_ROOT)),
            "sha256": p10_plot.sha256(SOURCE_DATA),
            "rows": len(rows),
        },
        "formal_review": {
            "path": str(REVIEW.relative_to(REPO_ROOT)),
            "sha256": p10_plot.sha256(REVIEW),
            "status": review["status"],
        },
        "display_change": "Neutralized recovery wording only; plotted values and source data are unchanged.",
        "artifacts": {
            key: {
                "path": str(path.relative_to(REPO_ROOT)),
                "sha256": p10_plot.sha256(path),
                "bytes": path.stat().st_size,
            }
            for key, path in paths.items()
        },
    }
    (output_dir / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    print(json.dumps(render(args.output_dir), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
