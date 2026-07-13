from __future__ import annotations

from pathlib import Path
import tempfile
import unittest

import numpy as np
import pyBigWig

from scripts import build_rna_seq_groups_v2 as p2
from scripts import finalize_v2_reprocessed_groups as finalizer


class V2GroupFinalizationTest(unittest.TestCase):
    def test_full_audit_schema_separates_reference_archive_and_extraction(self) -> None:
        row = {
            "run_accession": "SRR1",
            "source_transport_backend": "ncbi_sra",
            "source_fastq_urls": "ena.example/SRR1.fastq.gz",
            "source_fastq_md5": "a" * 32,
            "source_fastq_sha256": "b" * 64,
            "source_fastq_bytes": 123,
            "extracted_fastq_sha256": "b" * 64,
        }
        updated = finalizer.audit_schema_v2(row)
        self.assertEqual(updated["audit_schema_version"], 2)
        self.assertEqual(updated["source_reference_ena_fastq_md5"], "a" * 32)
        self.assertEqual(updated["extracted_fastq_sha256"], "b" * 64)
        self.assertNotIn("source_fastq_md5", updated)
        self.assertNotIn("source_fastq_sha256", updated)

    def test_final_hierarchy_restores_three_distinct_study_runs(self) -> None:
        samples = p2.read_tsv(finalizer.SAMPLE_MANIFEST)
        contexts = p2.read_tsv(finalizer.CONTEXT_MANIFEST)
        ontology = p2.read_tsv(finalizer.ONTOLOGY_MANIFEST)
        groups, members, _, _ = p2.build_groups(
            samples,
            contexts,
            ontology,
            exclude_duplicate_current_signals=False,
        )
        self.assertEqual(len(groups), 241)
        self.assertEqual(len(members), 482)
        restored = set(finalizer.DUPLICATE_PAIRS)
        restored_groups = {
            row["group_id"] for row in members if row["run_accession"] in restored
        }
        self.assertEqual(len(restored_groups), 1)
        restored_group = next(
            row for row in groups if row["group_id"] in restored_groups
        )
        self.assertEqual(restored_group["sra_study_accession"], "SRP310676")

    def test_nonzero_run_writer_preserves_piecewise_signal(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "test.bw"
            output = pyBigWig.open(str(path), "w")
            output.addHeader([("I", 8)])
            finalizer.write_nonzero_runs(
                output,
                "I",
                0,
                np.asarray([0, 1, 1, 0, 2, 2, 2, 0], dtype=np.float32),
            )
            output.close()
            with pyBigWig.open(str(path)) as bigwig:
                values = np.nan_to_num(bigwig.values("I", 0, 8, numpy=True))
            np.testing.assert_array_equal(values, [0, 1, 1, 0, 2, 2, 2, 0])


if __name__ == "__main__":
    unittest.main()
