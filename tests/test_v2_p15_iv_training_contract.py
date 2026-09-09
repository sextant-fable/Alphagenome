from __future__ import annotations

import csv
import unittest
from pathlib import Path

from scripts.prepare_v2_p15_iv_training_intervals import ALLOWED_SET, select_rows


REPO_ROOT = Path(__file__).resolve().parents[1]
INTERVAL_ROOT = REPO_ROOT / "alphagenome_custom" / "intervals" / "v2"


class P15IVTrainingContractTest(unittest.TestCase):
    def test_filter_excludes_x_from_train_and_validation_sources(self) -> None:
        for fold in range(1, 6):
            for role in ("train", "valid"):
                source = INTERVAL_ROOT / f"fold_{fold}" / f"{role}.tsv"
                with source.open(newline="", encoding="utf-8") as handle:
                    rows = list(csv.DictReader(handle, delimiter="\t"))
                self.assertTrue(any(row["chromosome"] == "X" for row in rows))
                selected = select_rows(rows, expected_role=role, source=source)
                self.assertTrue(selected)
                self.assertTrue(all(row["chromosome"] in ALLOWED_SET for row in selected))
                self.assertTrue(all(row["role"] == role for row in selected))
                self.assertTrue(all(row["role"] != "test_locked" for row in selected))


if __name__ == "__main__":
    unittest.main()
