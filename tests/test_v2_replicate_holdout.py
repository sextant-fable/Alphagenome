from __future__ import annotations

import unittest

from scripts import build_v2_replicate_holdout as holdout


class ReplicateHoldoutContractTest(unittest.TestCase):
    def test_biological_unit_prefers_replicate_then_technical_fallback(self) -> None:
        self.assertEqual(
            holdout.biological_unit(
                {"biological_replicate": "bio-1", "technical_unit_id": "tech-1"}
            ),
            ("biological_replicate", "bio-1"),
        )
        self.assertEqual(
            holdout.biological_unit(
                {"biological_replicate": "", "technical_unit_id": "tech-1"}
            ),
            ("technical_unit_fallback", "tech-1"),
        )

    def test_duplicate_classes_cannot_cross_role_boundary(self) -> None:
        rows = [
            {
                "run_accession": "A",
                "group_id": "G1",
                "role": "training",
                "duplicate_source_class": "class-1",
            },
            {
                "run_accession": "B",
                "group_id": "G2",
                "role": "heldout",
                "duplicate_source_class": "class-1",
            },
        ]
        with self.assertRaisesRegex(RuntimeError, "crosses"):
            holdout.validate_role_boundaries(rows)

    def test_prepared_assignment_has_expected_counts_when_available(self) -> None:
        path = holdout.HOLDOUT_METADATA / "holdout_assignments.tsv"
        if not path.exists():
            self.skipTest("P11 metadata has not been generated in this checkout")
        rows = holdout.read_tsv(path)
        counts = {role: sum(row["assignment_role"] == role for row in rows) for role in {
            "primary_candidate",
            "primary_pending_qc",
            "supplementary_candidate",
            "training_only",
        }}
        self.assertEqual(
            counts,
            {
                "primary_candidate": 57,
                "primary_pending_qc": 1,
                "supplementary_candidate": 22,
                "training_only": 161,
            },
        )
        self.assertEqual(len(rows), 241)


if __name__ == "__main__":
    unittest.main()
