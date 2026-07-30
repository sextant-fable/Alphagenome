from __future__ import annotations

import json
import unittest

from scripts import run_v2_p8_b_no_lora as p8


class P8BNoLoRATest(unittest.TestCase):
    def test_registered_commands_are_matched_and_never_request_final_test(self) -> None:
        spec = json.loads(p8.SPEC_PATH.read_text())
        train, evaluate = p8.commands(spec)
        self.assertIn("B_no_lora", train)
        self.assertIn("--fold", train)
        self.assertEqual(train[train.index("--fold") + 1], "1")
        self.assertEqual(evaluate[evaluate.index("--fold") + 1], "1")
        self.assertEqual(spec["training"]["max_steps"], 2000)
        self.assertIn("fold_1_train_nonzero_mean", train)
        self.assertNotIn("--final-test", train)
        self.assertNotIn("--final-test", evaluate)

    def test_embedding_only_contract_rejects_lora(self) -> None:
        self.assertTrue(
            p8.is_b_no_lora(
                {
                    "base_organism_index": 2,
                    "c_elegans_organism_embedding": True,
                    "lora_enabled": False,
                    "lora_target_modules": [],
                    "trunk_policy": "worm_embeddings_only",
                }
            )
        )
        self.assertFalse(
            p8.is_b_no_lora(
                {
                    "base_organism_index": 2,
                    "c_elegans_organism_embedding": True,
                    "lora_enabled": True,
                    "lora_target_modules": ["tower.blocks.8.mha"],
                    "trunk_policy": "worm_embeddings_and_lora_only",
                }
            )
        )


if __name__ == "__main__":
    unittest.main()
