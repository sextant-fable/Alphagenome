from __future__ import annotations

import tempfile
from pathlib import Path
import unittest

import numpy as np

from scripts.v2_validation_metrics import (
    GeneExonAccumulator,
    PerTrackStats,
    average_ranks,
    distribution_metrics,
    load_gene_exons,
)


class V2ValidationMetricsTests(unittest.TestCase):
    def test_gtf_exons_are_merged_per_gene(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "test.gtf"
            path.write_text(
                'I\ttest\texon\t1\t10\t.\t+\t.\tgene_id "g1"; transcript_id "a";\n'
                'I\ttest\texon\t6\t15\t.\t+\t.\tgene_id "g1"; transcript_id "b";\n'
                'I\ttest\texon\t21\t25\t.\t+\t.\tgene_id "g1"; transcript_id "a";\n'
                'II\ttest\texon\t3\t4\t.\t-\t.\tgene_id "g2";\n'
            )
            genes, by_chromosome = load_gene_exons(path, {"I"})
        self.assertEqual(len(genes), 1)
        self.assertEqual(genes[0].gene_id, "g1")
        self.assertEqual(genes[0].intervals, ((0, 15), (20, 25)))
        self.assertEqual(by_chromosome, {"I": [0]})

    def test_gene_exon_accumulator_uses_unique_bases_across_windows(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "test.gtf"
            path.write_text(
                'I\ttest\texon\t3\t6\t.\t+\t.\tgene_id "g1";\n'
                'I\ttest\texon\t11\t14\t.\t+\t.\tgene_id "g2";\n'
            )
            genes, by_chromosome = load_gene_exons(path)
        accumulator = GeneExonAccumulator(genes, by_chromosome, n_tracks=2)
        first_target = np.vstack([np.arange(8), np.arange(8) * 2]).astype(float)
        accumulator.update("I", 0, first_target + 1, first_target)
        second_target = np.vstack([np.arange(8, 16), np.arange(8, 16) * 2]).astype(float)
        accumulator.update("I", 8, second_target + 1, second_target)
        self.assertEqual(accumulator.base_count.tolist(), [4, 4])
        metrics = accumulator.metrics(log1p=False)
        self.assertEqual(metrics["genes_evaluated"], 2)
        self.assertEqual(metrics["exon_bases_evaluated"], 8)
        self.assertTrue(np.isfinite(metrics["mean_per_track_pearson"]))

    def test_per_track_stats_match_direct_calculation(self) -> None:
        prediction = np.array([[1.0, 2.0, 4.0], [2.0, 4.0, 6.0]])
        target = np.array([[1.0, 3.0, 5.0], [1.0, 2.0, 3.0]])
        stats = PerTrackStats(2)
        stats.update(prediction[:, :2], target[:, :2])
        stats.update(prediction[:, 2:], target[:, 2:])
        result = stats.per_track()
        np.testing.assert_allclose(
            result["mse"], np.square(prediction - target).mean(axis=1)
        )
        np.testing.assert_allclose(result["pearson"], [0.9819805061, 1.0])

    def test_average_ranks_assigns_tie_midpoints(self) -> None:
        np.testing.assert_allclose(
            average_ranks(np.array([10.0, 5.0, 5.0, 20.0])),
            [2.0, 0.5, 0.5, 3.0],
        )

    def test_distribution_metrics_are_deterministic_and_calibrated(self) -> None:
        target = np.vstack([np.arange(200), np.arange(200)[::-1]]).astype(float)
        prediction = target * 2
        first = distribution_metrics(prediction, target, max_spearman_points=50, seed=7)
        second = distribution_metrics(prediction, target, max_spearman_points=50, seed=7)
        np.testing.assert_allclose(first["per_track_spearman_128bp"], [1.0, 1.0])
        np.testing.assert_allclose(
            first["per_track_top1_calibration_ratio_128bp"], [2.0, 2.0]
        )
        np.testing.assert_allclose(
            first["per_track_spearman_128bp"], second["per_track_spearman_128bp"]
        )


if __name__ == "__main__":
    unittest.main()
