from __future__ import annotations

import hashlib
from pathlib import Path
import tempfile
import unittest

from scripts import run_v2_reprocessing_full as full


class V2ReprocessingFullTest(unittest.TestCase):
    def test_full_source_manifest_is_exactly_the_rna_scope(self) -> None:
        sources = full.read_tsv(full.SOURCE_MANIFEST)
        samples = full.read_tsv(full.SAMPLE_MANIFEST)
        sample_by_run = full.validate_manifests(sources, samples)
        self.assertEqual(len(sources), 482)
        self.assertEqual(
            {row["run_accession"] for row in sources},
            {
                row["run_accession"]
                for row in samples
                if row["assay"] == "RNA-Seq"
            },
        )
        self.assertEqual(len(sample_by_run), 485)

    def test_full_scope_excludes_the_three_chip_runs(self) -> None:
        sources = full.read_tsv(full.SOURCE_MANIFEST)
        accessions = {row["run_accession"] for row in sources}
        self.assertTrue(
            {"SRR3535777", "SRR3535778", "SRR3535779"}.isdisjoint(
                accessions
            )
        )

    def test_complete_fastq_is_reused_without_network_call(self) -> None:
        payload = b"verified-fastq"
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "sample.fastq.gz"
            path.write_bytes(payload)

            class NoNetworkRunner:
                def run(self, command: list[str], *, stdout=None) -> None:
                    raise AssertionError(command)

            actual_sha, actual_md5 = full.download_fastq(
                NoNetworkRunner(),
                "example.invalid/sample.fastq.gz",
                hashlib.md5(payload).hexdigest(),
                len(payload),
                path,
            )
            self.assertEqual(actual_sha, hashlib.sha256(payload).hexdigest())
            self.assertEqual(actual_md5, hashlib.md5(payload).hexdigest())


if __name__ == "__main__":
    unittest.main()
