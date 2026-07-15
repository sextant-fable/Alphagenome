from __future__ import annotations

from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

import numpy as np
import pyBigWig

from scripts import build_rna_seq_groups_v2 as p2
from scripts import finalize_v2_reprocessed_groups as finalizer


class V2GroupFinalizationTest(unittest.TestCase):
    def test_finalizer_imports_when_launched_from_scripts_directory(self) -> None:
        scripts_dir = Path(__file__).resolve().parents[1] / "scripts"
        result = subprocess.run(
            [sys.executable, "-c", "import finalize_v2_reprocessed_groups"],
            cwd=scripts_dir,
            capture_output=True,
            text=True,
            check=False,
        )
        self.assertEqual(result.returncode, 0, result.stderr)

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

    def test_full_audit_schema_records_ena_fallback_without_calling_it_extraction(self) -> None:
        row = {
            "run_accession": "SRR2",
            "source_transport_backend": "ena_fastq_fallback_after_sra_extraction_failure",
            "source_fastq_urls": "ena.example/SRR2_1.fastq.gz;ena.example/SRR2_2.fastq.gz",
            "source_fastq_md5": f"{'a' * 32};{'b' * 32}",
            "source_fastq_sha256": f"{'c' * 64};{'d' * 64}",
            "source_fastq_bytes": 456,
            "source_archive_md5": "e" * 32,
            "source_archive_sha256": "f" * 64,
            "downloaded_ena_fastq_sha256": f"{'c' * 64};{'d' * 64}",
            "downloaded_ena_fastq_bytes": 456,
            "ena_fastq_expected_record_count": 20,
            "sra_extraction_return_code": 3,
        }
        updated = finalizer.audit_schema_v2(row)
        self.assertEqual(updated["audit_schema_version"], 2)
        self.assertEqual(updated["downloaded_ena_fastq_bytes"], 456)
        self.assertEqual(updated["extracted_fastq_sha256"], "")
        self.assertEqual(
            updated["downloaded_ena_fastq_integrity"],
            "ENA_MD5_and_local_SHA256_plus_STAR_spot_count",
        )
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

    def test_group_aggregation_rescales_decoded_members_to_target(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            member_paths = []
            for index, values in enumerate(([1.0, 1.0], [2.0, 2.0]), start=1):
                path = root / f"member_{index}.bw"
                output = pyBigWig.open(str(path), "w")
                output.addHeader([("I", 2)])
                output.addEntries(
                    ["I", "I"], [0, 1], ends=[1, 2], values=list(values)
                )
                output.close()
                member_paths.append(path)
            members = [
                {
                    "run_accession": "SRR1",
                    "local_path": str(member_paths[0]),
                    "sha256": "a" * 64,
                    "raw_bedgraph_total": "10",
                    "decoded_total_signal": "2",
                    "decoded_rescale_to_1e8": "50000000",
                    "biological_replicate": "1",
                    "technical_unit_id": "SRX1",
                },
                {
                    "run_accession": "SRR2",
                    "local_path": str(member_paths[1]),
                    "sha256": "b" * 64,
                    "raw_bedgraph_total": "20",
                    "decoded_total_signal": "4",
                    "decoded_rescale_to_1e8": "25000000",
                    "biological_replicate": "2",
                    "technical_unit_id": "SRX2",
                },
            ]
            with (
                mock.patch.object(finalizer, "REPO_ROOT", root),
                mock.patch.object(finalizer, "GROUPED_DIR", root / "grouped"),
            ):
                audit = finalizer.aggregate_group(
                    {"group_id": "RNA_V2_TEST"}, members, {"I": 2}
                )
                output_path = root / "grouped/RNA_V2_TEST.bw"
                stats = finalizer.validate_bigwig(
                    output_path, {"I": 2}, expected_decoded_total=finalizer.TARGET_TOTAL
                )
            self.assertFalse(output_path.is_symlink())
            self.assertLess(
                stats["decoded_relative_target_error"],
                finalizer.DECODED_TOTAL_RELATIVE_TOLERANCE,
            )
            self.assertEqual(
                audit["aggregation_policy"],
                "decoded_signal_rescaled_to_1e8_then_raw_coverage_weighted_within_"
                "biological_unit_then_equal_mean_across_units",
            )


if __name__ == "__main__":
    unittest.main()
