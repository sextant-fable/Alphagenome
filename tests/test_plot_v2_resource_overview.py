"""Focused regression tests for the manifest-derived resource overview."""

from __future__ import annotations

import unittest

from scripts.plot_v2_resource_overview import DEFAULT_MANIFEST, read_manifest, stage_counts


class ResourceOverviewTests(unittest.TestCase):
    def test_formal_manifest_has_expected_inventory(self) -> None:
        rows = read_manifest(DEFAULT_MANIFEST)
        self.assertEqual(len(rows), 241)
        self.assertEqual(sum(int(row["n_runs"]) for row in rows), 482)
        self.assertEqual(len({row["sra_study_accession"] for row in rows}), 27)

    def test_stage_summary_is_complete(self) -> None:
        rows = read_manifest(DEFAULT_MANIFEST)
        counts = stage_counts(rows)
        self.assertEqual(sum(value for _, value in counts), 241)
        self.assertEqual(counts[0], ("Young Adult", 65))
        self.assertEqual(counts[1], ("L1", 61))
        self.assertIn(("Other stages", 43), counts)


if __name__ == "__main__":
    unittest.main()
