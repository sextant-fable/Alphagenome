from __future__ import annotations

import json
import unittest

from scripts import run_v2_p6b_amendment as amendment


class P6BAmendmentTests(unittest.TestCase):
    def setUp(self) -> None:
        self.spec = json.loads(amendment.SPEC_PATH.read_text())

    def test_formal_matrix_has_three_seeds_and_reuses_only_original_training(self) -> None:
        jobs = amendment.formal_jobs(self.spec)
        self.assertEqual(len(jobs), 90)
        self.assertEqual(len({job["job_id"] for job in jobs}), 90)
        self.assertEqual(sum(bool(job["base_job"]) for job in jobs), 30)
        self.assertEqual(sum(not bool(job["base_job"]) for job in jobs), 60)
        self.assertEqual({job["seed"] for job in jobs}, {20260714, 20260715, 20260716})
        self.assertEqual({job["fold"] for job in jobs}, {1, 2, 3, 4, 5})

    def test_ablation_matrix_has_all_registered_five_fold_controls(self) -> None:
        jobs = amendment.ablation_jobs(self.spec)
        self.assertEqual(len(jobs), 15)
        self.assertEqual(len({job["job_id"] for job in jobs}), 15)
        self.assertEqual(
            {job["ablation_id"] for job in jobs},
            {"no_augmentation", "no_gene_loss", "whole_I_V_mean"},
        )
        self.assertEqual({job["fold"] for job in jobs}, {1, 2, 3, 4, 5})

    def test_no_amendment_job_references_final_test_inputs(self) -> None:
        jobs = amendment.formal_jobs(self.spec) + amendment.ablation_jobs(self.spec)
        serialized = json.dumps(jobs)
        self.assertNotIn("test_locked", serialized)
        self.assertNotIn('"fold": 0', serialized)
        self.assertEqual(self.spec["chromosome_x_access"], "prohibited")


if __name__ == "__main__":
    unittest.main()
