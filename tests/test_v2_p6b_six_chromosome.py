from __future__ import annotations

import json
import unittest

from scripts import run_v2_p6b_six_chromosome as six_chromosome


class P6BSixChromosomeTests(unittest.TestCase):
    def setUp(self) -> None:
        self.spec = json.loads(six_chromosome.SPEC_PATH.read_text())

    def test_registered_matrix_uses_all_five_block_folds(self) -> None:
        jobs = six_chromosome.formal_jobs(self.spec)
        self.assertEqual(len(jobs), 90)
        self.assertEqual(len({job["job_id"] for job in jobs}), 90)
        self.assertTrue(all(job["base_job"] is False for job in jobs))
        self.assertEqual({job["fold"] for job in jobs}, {1, 2, 3, 4, 5})
        self.assertEqual(
            {job["mean_column"] for job in jobs},
            {
                "fold_1_train_nonzero_mean",
                "fold_2_train_nonzero_mean",
                "fold_3_train_nonzero_mean",
                "fold_4_train_nonzero_mean",
                "fold_5_train_nonzero_mean",
            },
        )

    def test_ablation_jobs_do_not_reference_locked_test_blocks(self) -> None:
        jobs = six_chromosome.ablation_jobs(self.spec)
        self.assertEqual(len(jobs), 15)
        self.assertEqual(
            {job["ablation_id"] for job in jobs},
            {"no_augmentation", "no_gene_loss", "development_pool_mean"},
        )
        serialized = json.dumps(jobs)
        self.assertNotIn("test_locked", serialized)
        self.assertNotIn('"fold": 0', serialized)


if __name__ == "__main__":
    unittest.main()
