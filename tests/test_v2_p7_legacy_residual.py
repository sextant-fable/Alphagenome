from __future__ import annotations

import json
import unittest

from scripts import run_v2_p7_legacy_residual as p7


class P7LegacyResidualTest(unittest.TestCase):
    def test_registered_commands_are_fold_one_and_never_request_final_test(self) -> None:
        spec = json.loads(p7.SPEC_PATH.read_text())
        base_train, residual_train, evaluate = p7.commands(spec)
        self.assertIn("D_base", base_train)
        self.assertIn("D", residual_train)
        self.assertIn("--fold", residual_train)
        self.assertEqual(residual_train[residual_train.index("--fold") + 1], "1")
        self.assertEqual(evaluate[evaluate.index("--fold") + 1], "1")
        self.assertIn("--frozen-base-checkpoint", residual_train)
        self.assertIn("--frozen-base-checkpoint", evaluate)
        self.assertEqual(spec["training"]["base_max_steps"], 2000)
        self.assertEqual(spec["training"]["residual_max_steps"], 1500)
        self.assertNotIn("--final-test", base_train)
        self.assertNotIn("--final-test", residual_train)
        self.assertNotIn("--final-test", evaluate)
        self.assertIn("fold_1_train_nonzero_mean", residual_train)


if __name__ == "__main__":
    unittest.main()
