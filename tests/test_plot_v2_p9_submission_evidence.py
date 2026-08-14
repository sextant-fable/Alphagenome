from __future__ import annotations

import json
from pathlib import Path
import tempfile
import unittest

try:
    import matplotlib  # noqa: F401

    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False


@unittest.skipUnless(HAS_MATPLOTLIB, "matplotlib is available in the plot environment")
class P9SubmissionFigureTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        from scripts import plot_v2_p9_submission_evidence as plot

        cls.plot = plot
        cls.inputs = plot.load_figure_inputs()
        cls.data = plot.prepare_figure_data(cls.inputs)

    def test_real_r9_contract_and_source_data_are_complete(self) -> None:
        rows = self.plot.build_source_rows(self.data)
        self.assertEqual(self.inputs["review"]["status"], "PASS")
        self.assertEqual(len(self.data["paired_summary"]), 5)
        self.assertEqual(len(self.data["factorial_summary"]), 3)
        self.assertEqual(len(self.data["track_effects"]), 5 * 241)
        self.assertEqual(len(rows), 1373)
        self.assertEqual(
            len({row["track_id"] for row in rows if row["panel"] == "c_track_distribution"}),
            241,
        )
        self.assertTrue(
            all(
                row["inference_role"] == "descriptive_track_outcome"
                for row in rows
                if row["panel"] == "c_track_distribution"
            )
        )

    def test_required_claim_boundaries_are_encoded(self) -> None:
        head = self.data["paired_summary"]["learned_1bp_head_only"]
        self.assertLess(float(head["ci_low"]), 0.0)
        self.assertGreater(float(head["ci_high"]), 0.0)
        interaction = self.data["factorial_summary"]["worm_lora_interaction"]
        self.assertLess(float(interaction["mean_effect"]), 0.0)
        self.assertEqual(self.data["track_count"], 241)

    def test_four_format_export_preserves_fixed_geometry(self) -> None:
        figure = self.plot.build_figure(self.data)
        try:
            with tempfile.TemporaryDirectory() as directory:
                paths = self.plot.export_figure(
                    figure, Path(directory), raster_dpi=150
                )
                self.assertEqual(set(paths), {"pdf", "svg", "png", "tiff"})
                self.assertTrue(all(path.stat().st_size > 0 for path in paths.values()))
        finally:
            self.plot.plt.close(figure)

    def test_non_pass_review_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            review = dict(self.inputs["review"])
            review["status"] = "FAIL"
            path = Path(directory) / "review.json"
            path.write_text(json.dumps(review))
            with self.assertRaisesRegex(RuntimeError, "formal R9 PASS"):
                self.plot.load_figure_inputs(review_path=path)


if __name__ == "__main__":
    unittest.main()
