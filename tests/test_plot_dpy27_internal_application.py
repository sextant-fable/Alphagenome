from __future__ import annotations

import json
from pathlib import Path
import sys
import tempfile
import unittest

try:
    import matplotlib  # noqa: F401

    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False


@unittest.skipUnless(HAS_MATPLOTLIB, "matplotlib is available only in the frozen plot environment")
class Dpy27FigureTests(unittest.TestCase):
    @staticmethod
    def source_rows() -> list[dict[str, str]]:
        rows: list[dict[str, str]] = []
        for fold in range(1, 6):
            for gene_index, chromosome in enumerate(("I", "II", "X", "X"), 1):
                observed = -0.35 + 0.16 * gene_index + 0.02 * fold
                predicted = observed * 0.84 + 0.025 * (fold - 3)
                rows.append(
                    {
                        "panel": "a_gene_contrasts",
                        "fold": str(fold),
                        "seed": "seed_mean",
                        "gene_id": f"g{fold}_{gene_index}",
                        "gene_name": f"g{fold}_{gene_index}",
                        "chromosome": chromosome,
                        "endpoint": "gene_log1p_contrast",
                        "observed": str(observed),
                        "predicted": str(predicted),
                        "estimate": "",
                        "ci_95_low": "",
                        "ci_95_high": "",
                    }
                )
            for kind, value in (
                ("observed", 0.12 + 0.01 * fold),
                ("predicted", 0.09 + 0.012 * fold),
            ):
                rows.append(
                    {
                        "panel": "b_fold_x_minus_autosomes",
                        "fold": str(fold),
                        "seed": "seed_mean",
                        "gene_id": "",
                        "gene_name": "",
                        "chromosome": "X_minus_autosomes",
                        "endpoint": f"{kind}_median_contrast_x_minus_autosomes",
                        "observed": str(value) if kind == "observed" else "",
                        "predicted": str(value) if kind == "predicted" else "",
                        "estimate": str(value),
                        "ci_95_low": "",
                        "ci_95_high": "",
                    }
                )
            for endpoint, value in (
                ("x_predicted_observed_direction_concordance", 0.68 + 0.03 * fold),
                ("predicted_observed_gene_spearman", 0.48 + 0.04 * fold),
            ):
                rows.append(
                    {
                        "panel": "c_fold_agreement",
                        "fold": str(fold),
                        "seed": "seed_mean",
                        "gene_id": "",
                        "gene_name": "",
                        "chromosome": "",
                        "endpoint": endpoint,
                        "observed": "",
                        "predicted": "",
                        "estimate": str(value),
                        "ci_95_low": "",
                        "ci_95_high": "",
                    }
                )
        summaries = (
            ("observed_median_contrast_x_minus_autosomes", 0.15, 0.13, 0.17),
            ("predicted_median_contrast_x_minus_autosomes", 0.126, 0.105, 0.145),
            ("x_predicted_observed_direction_concordance", 0.77, 0.71, 0.83),
            ("predicted_observed_gene_spearman", 0.60, 0.52, 0.68),
        )
        for endpoint, estimate, low, high in summaries:
            rows.append(
                {
                    "panel": "c_fold_primary_summary",
                    "fold": "all_folds",
                    "seed": "within_fold_seed_mean",
                    "gene_id": "",
                    "gene_name": "",
                    "chromosome": "",
                    "endpoint": endpoint,
                    "observed": "",
                    "predicted": "",
                    "estimate": str(estimate),
                    "ci_95_low": str(low),
                    "ci_95_high": str(high),
                }
            )
        return rows

    def test_four_format_export_and_raster_pdf_qa(self) -> None:
        from scripts import plot_dpy27_internal_application as plot

        spec = json.loads(plot.SPEC_PATH.read_text())
        runtime = spec["plot_runtime"]
        plot.validate_source_data_contract(self.source_rows(), spec)
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            paths = {
                suffix: root / f"synthetic.{suffix}"
                for suffix in ("pdf", "svg", "tiff", "png")
            }
            figure = plot.build_figure(self.source_rows())
            try:
                figure.savefig(paths["pdf"], format="pdf", bbox_inches="tight")
                figure.savefig(paths["svg"], format="svg", bbox_inches="tight")
                figure.savefig(paths["tiff"], format="tiff", dpi=plot.RASTER_DPI, bbox_inches="tight")
                figure.savefig(paths["png"], format="png", dpi=plot.RASTER_DPI, bbox_inches="tight")
            finally:
                plot.plt.close(figure)
            self.assertTrue(all(path.stat().st_size > 0 for path in paths.values()))
            self.assertIn("<text", paths["svg"].read_text())
            png_audit = plot._image_audit(paths["png"])
            tiff_audit = plot._image_audit(paths["tiff"])
            self.assertTrue(png_audit["passed"])
            self.assertTrue(tiff_audit["passed"])
            json.dumps({"png": png_audit, "tiff": tiff_audit})
            pdf_audit = plot._run_json(
                [
                    sys.executable,
                    runtime["audit_pdf_text_path"],
                    str(paths["pdf"]),
                    "--min-pt",
                    "5",
                    "--json",
                ]
            )
            self.assertTrue(pdf_audit["auditable"])
            self.assertEqual(pdf_audit["below_minimum_count"], 0)

    def test_ifd_rational_dpi_is_normalized_for_json(self) -> None:
        from PIL.TiffImagePlugin import IFDRational

        from scripts import plot_dpy27_internal_application as plot

        raw_dpi = (IFDRational(600, 1), IFDRational(300, 1))
        with self.assertRaisesRegex(TypeError, "IFDRational"):
            json.dumps({"dpi": list(raw_dpi)})

        normalized = plot._json_safe_dpi(raw_dpi)
        self.assertEqual(normalized, [600.0, 300.0])
        self.assertTrue(all(type(value) is float for value in normalized))
        self.assertEqual(
            json.loads(json.dumps({"dpi": normalized})),
            {"dpi": [600.0, 300.0]},
        )


if __name__ == "__main__":
    unittest.main()
