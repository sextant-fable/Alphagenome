from __future__ import annotations

import json
import unittest

from scripts import build_v2_splits
from scripts import run_v2_p6b_six_chromosome


class SixChromosomeSplitTest(unittest.TestCase):
    def setUp(self) -> None:
        self.spec = {
            "block_labels": [
                "cv_fold_1",
                "cv_fold_2",
                "cv_fold_3",
                "cv_fold_4",
                "cv_fold_5",
                "test_locked",
            ],
            "exclusion_buffer_bp": build_v2_splits.EXCLUSION_BUFFER_BP,
            "minimum_effective_block_bp": build_v2_splits.WINDOW_SIZE
            + 2 * build_v2_splits.MAX_SHIFT_BP,
            "seed": 20260721,
        }

    def test_each_chromosome_has_six_labeled_blocks_and_safe_buffers(self) -> None:
        aligned_minimum = (
            (
                self.spec["minimum_effective_block_bp"]
                + build_v2_splits.EVALUATION_SUBWINDOW_BP
                - 1
            )
            // build_v2_splits.EVALUATION_SUBWINDOW_BP
            * build_v2_splits.EVALUATION_SUBWINDOW_BP
        )
        minimum_length = (
            len(self.spec["block_labels"])
            * aligned_minimum
            + (len(self.spec["block_labels"]) - 1)
            * self.spec["exclusion_buffer_bp"]
        )
        rows = build_v2_splits.split_blocks(
            {
                chromosome: minimum_length + index
                for index, chromosome in enumerate(build_v2_splits.CHROMOSOMES)
            },
            self.spec,
        )
        for chromosome in build_v2_splits.CHROMOSOMES:
            blocks = sorted(
                (row for row in rows if row["chromosome"] == chromosome),
                key=lambda row: row["physical_order"],
            )
            self.assertEqual(
                {row["assignment"] for row in blocks},
                set(self.spec["block_labels"]),
            )
            self.assertTrue(
                all(
                    int(row["effective_bases"])
                    >= aligned_minimum
                    for row in blocks
                )
            )
            self.assertTrue(
                all(
                    int(second["block_start"]) - int(first["block_end"])
                    == self.spec["exclusion_buffer_bp"]
                    for first, second in zip(blocks, blocks[1:])
                )
            )

    def test_train_windows_and_eval_cores_stay_inside_one_block(self) -> None:
        block_start = 5 * build_v2_splits.EXCLUSION_BUFFER_BP
        block_end = block_start + 2 * build_v2_splits.WINDOW_SIZE + 2 * build_v2_splits.MAX_SHIFT_BP
        length = block_end + build_v2_splits.WINDOW_SIZE
        block = {
            "block_id": "I_block_1",
            "block_start": block_start,
            "block_end": block_end,
        }
        train = build_v2_splits.interval_rows(
            "I", length, 1, "train", block, "six_chromosome_blocks_v1"
        )
        self.assertTrue(train)
        self.assertTrue(
            all(
                int(row["start"]) - build_v2_splits.MAX_SHIFT_BP >= block_start
                and int(row["end"]) + build_v2_splits.MAX_SHIFT_BP <= block_end
                for row in train
            )
        )
        valid = build_v2_splits.interval_rows(
            "I", length, 1, "valid", block, "six_chromosome_blocks_v1"
        )
        ordered = sorted(valid, key=lambda row: int(row["core_start"]))
        self.assertEqual(int(ordered[0]["core_start"]), block_start)
        self.assertEqual(int(ordered[-1]["core_end"]), block_end)
        self.assertTrue(
            all(
                int(first["core_end"]) == int(second["core_start"])
                for first, second in zip(ordered, ordered[1:])
            )
        )

    def test_p6b_matrix_is_preregistered_and_excludes_test_manifests(self) -> None:
        spec = json.loads(run_v2_p6b_six_chromosome.SPEC_PATH.read_text())
        formal = run_v2_p6b_six_chromosome.formal_jobs(spec)
        ablations = run_v2_p6b_six_chromosome.ablation_jobs(spec)
        self.assertEqual(len(formal), 90)
        self.assertEqual(len(ablations), 15)
        self.assertEqual(
            {(row["model"], row["loss"]) for row in formal},
            {
                (model, loss)
                for model in ("A", "B", "C")
                for loss in ("paper", "log1p_mse")
            },
        )
        self.assertEqual(
            {row["ablation_id"] for row in ablations},
            {"no_augmentation", "no_gene_loss", "development_pool_mean"},
        )
        self.assertTrue(
            all("test_locked" not in json.dumps(row) for row in formal + ablations)
        )


if __name__ == "__main__":
    unittest.main()
