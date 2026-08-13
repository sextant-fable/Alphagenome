from __future__ import annotations

import json
from pathlib import Path
import unittest

import numpy as np

from scripts.analyze_v2_p6b_submission import (
    DEFAULT_SPEC,
    UNKNOWN,
    clean_text,
    iter_track_values,
    load_and_validate_inputs,
    stage_stratum,
    tissue_stratum,
    treatment_stratum,
)
from scripts.v2_submission_statistics import (
    fold_primary_estimate,
    hierarchical_block_bootstrap,
    make_hierarchical_resample_plan,
    matched_differences,
    records_to_balanced_matrix,
)


class SubmissionStatisticsTests(unittest.TestCase):
    def test_balanced_matrix_preserves_block_and_repeat_order(self) -> None:
        records = [
            {"fold": fold, "seed": seed, "score": fold * 10 + seed}
            for fold in (2, 1)
            for seed in (8, 7)
        ]
        matrix = records_to_balanced_matrix(
            records,
            value_key="score",
            expected_blocks=(1, 2),
            expected_repeats=(7, 8),
        )
        np.testing.assert_allclose(matrix.values, [[17, 18], [27, 28]])
        self.assertEqual(matrix.blocks, (1, 2))
        self.assertEqual(matrix.repeats, (7, 8))

    def test_balanced_matrix_rejects_missing_and_duplicate_cells(self) -> None:
        with self.assertRaisesRegex(ValueError, "unbalanced"):
            records_to_balanced_matrix(
                [{"fold": 1, "seed": 7, "score": 1.0}],
                value_key="score",
                expected_blocks=(1, 2),
                expected_repeats=(7,),
            )
        with self.assertRaisesRegex(ValueError, "duplicate"):
            records_to_balanced_matrix(
                [
                    {"fold": 1, "seed": 7, "score": 1.0},
                    {"fold": 1, "seed": 7, "score": 2.0},
                ],
                value_key="score",
            )

    def test_matched_differences_require_identical_fold_seed_keys(self) -> None:
        candidate = [
            {"fold": 1, "seed": 7, "score": 3.0},
            {"fold": 2, "seed": 7, "score": 6.0},
        ]
        comparator = [
            {"fold": 1, "seed": 7, "score": 1.0},
            {"fold": 2, "seed": 7, "score": 4.0},
        ]
        differences = matched_differences(
            candidate, comparator, value_key="score"
        )
        self.assertEqual([row["difference"] for row in differences], [2.0, 2.0])
        with self.assertRaisesRegex(ValueError, "paired keys differ"):
            matched_differences(candidate, comparator[:1], value_key="score")

    def test_fold_primary_estimate_weights_folds_equally(self) -> None:
        values = np.asarray([[0.0, 0.0, 3.0], [10.0, 10.0, 10.0]])
        self.assertAlmostEqual(fold_primary_estimate(values), 5.5)

    def test_hierarchical_bootstrap_is_reproducible_and_fold_first(self) -> None:
        values = np.asarray(
            [
                [0.0, 1.0, 2.0],
                [10.0, 11.0, 12.0],
                [20.0, 21.0, 22.0],
                [30.0, 31.0, 32.0],
                [40.0, 41.0, 42.0],
            ]
        )
        first_plan = make_hierarchical_resample_plan(
            n_blocks=5, n_repeats=3, iterations=2000, seed=13
        )
        second_plan = make_hierarchical_resample_plan(
            n_blocks=5, n_repeats=3, iterations=2000, seed=13
        )
        np.testing.assert_array_equal(first_plan.block_indices, second_plan.block_indices)
        np.testing.assert_array_equal(first_plan.repeat_indices, second_plan.repeat_indices)
        self.assertEqual(first_plan.block_indices.shape, (2000, 5))
        self.assertEqual(first_plan.repeat_indices.shape, (2000, 5, 3))
        first = hierarchical_block_bootstrap(values, first_plan)
        second = hierarchical_block_bootstrap(values, second_plan)
        self.assertEqual(first, second)
        self.assertAlmostEqual(first.estimate, 21.0)
        self.assertLess(first.ci_lower, first.estimate)
        self.assertGreater(first.ci_upper, first.estimate)

    def test_constant_paired_effect_has_degenerate_interval(self) -> None:
        values = np.full((5, 3), 0.125)
        plan = make_hierarchical_resample_plan(
            n_blocks=5, n_repeats=3, iterations=1000, seed=17
        )
        estimate = hierarchical_block_bootstrap(values, plan)
        self.assertAlmostEqual(estimate.estimate, 0.125)
        self.assertAlmostEqual(estimate.ci_lower, 0.125)
        self.assertAlmostEqual(estimate.ci_upper, 0.125)


class SubmissionStrataTests(unittest.TestCase):
    def test_missing_values_are_retained_explicitly(self) -> None:
        self.assertEqual(clean_text(""), UNKNOWN)
        self.assertEqual(clean_text("unknown"), UNKNOWN)
        self.assertEqual(treatment_stratum("", "", ""), UNKNOWN)

    def test_treatment_rule_is_deterministic(self) -> None:
        self.assertEqual(
            treatment_stratum("DMSO", "WT", "N2"),
            "Reported treatment/selection",
        )
        self.assertEqual(
            treatment_stratum("RNA-seq", "lin-35 RNAi", "N2"),
            "Reported treatment/selection",
        )
        self.assertEqual(
            treatment_stratum("synchronized L4 stage animals", "WT", "N2"),
            "No reported treatment",
        )

    def test_stage_and_tissue_collapse_only_documented_synonyms(self) -> None:
        self.assertEqual(stage_stratum("Embryo T4"), "Embryo")
        self.assertEqual(stage_stratum("L1-L3(10H)"), "Larval time course")
        self.assertEqual(stage_stratum("Young Adult"), "Young adult")
        self.assertEqual(tissue_stratum("Whole worm"), "Whole organism")
        self.assertEqual(tissue_stratum("early embryo"), "Embryo")
        self.assertEqual(tissue_stratum("L3 seamcell"), "L3 seamcell")


@unittest.skipUnless(
    DEFAULT_SPEC.is_file()
    and Path("runs/v2_p6b_six_chromosome_corefix_20260722/formal").is_dir(),
    "local P6B validation outputs are not available",
)
class SubmissionRealInputContractTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.spec = json.loads(DEFAULT_SPEC.read_text())
        cls.loaded, cls.metadata, cls.signal_audit = load_and_validate_inputs(cls.spec)

    def test_real_matrix_has_all_runs_tracks_and_guards(self) -> None:
        self.assertEqual(len(self.loaded), 90)
        self.assertEqual(len(self.metadata), 241)
        self.assertTrue(
            all(
                item["validation"]["locked_test_block_signal_reads"] is False
                for item in self.loaded
            )
        )
        self.assertTrue(
            all(
                item["validation"]["validation_core_coverage_fraction"] == 1.0
                for item in self.loaded
            )
        )
        self.assertEqual(self.signal_audit["counts"], {"High": 80, "Low": 81, "Middle": 80})

    def test_real_validation_exposes_all_track_metric_dictionaries(self) -> None:
        values = list(
            iter_track_values(
                self.loaded[0]["validation"], set(self.metadata)
            )
        )
        self.assertEqual(len(values), 241 * 39)
        self.assertEqual(len({value[1] for value in values}), 39)
        self.assertEqual(len({value[0] for value in values}), 241)


if __name__ == "__main__":
    unittest.main()
