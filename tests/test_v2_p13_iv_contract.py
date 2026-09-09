from __future__ import annotations

import csv
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
INTERVAL_ROOT = REPO_ROOT / "alphagenome_custom" / "intervals" / "v2"


class P13IVIntervalContractTest(unittest.TestCase):
    def test_source_validation_manifests_require_an_x_exclusion(self) -> None:
        for fold in range(1, 6):
            path = INTERVAL_ROOT / f"fold_{fold}" / "valid.tsv"
            with path.open(newline="", encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertTrue(any(row["chromosome"] == "X" for row in rows))
            selected = [
                row
                for row in rows
                if row["chromosome"] in {"I", "II", "III", "IV", "V"}
                and row["role"] == "valid"
            ]
            self.assertTrue(selected)
            self.assertTrue(all(row["chromosome"] != "X" for row in selected))
            self.assertTrue(all(row["role"] != "test_locked" for row in selected))


if __name__ == "__main__":
    unittest.main()
