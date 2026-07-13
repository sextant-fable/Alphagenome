from __future__ import annotations

import hashlib
from pathlib import Path
import tempfile
import unittest

from scripts import run_v2_reprocessing_full as full


class V2ReprocessingFullTest(unittest.TestCase):
    def test_ncbi_sdl_selects_full_sra_not_noqual_lite(self) -> None:
        payload = {
            "result": [
                {
                    "status": 200,
                    "files": [
                        {
                            "type": "sra",
                            "name": "SRR1",
                            "size": 123,
                            "md5": "a" * 32,
                            "locations": [
                                {
                                    "service": "s3",
                                    "link": "https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR1/SRR1",
                                }
                            ],
                        },
                        {
                            "type": "sra",
                            "name": "SRR1.lite",
                            "size": 12,
                            "md5": "b" * 32,
                            "noqual": True,
                            "locations": [],
                        },
                    ],
                }
            ]
        }
        self.assertEqual(
            full.select_sra_archive(payload, "SRR1"),
            {
                "url": "https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR1/SRR1",
                "md5": "a" * 32,
                "size": 123,
            },
        )

    def test_fastq_stats_hashes_and_counts_reads(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            first = Path(directory) / "first.fastq"
            second = Path(directory) / "second.fastq"
            first.write_bytes(b"@r1\nA\n+\n!\n")
            second.write_bytes(b"@r2\nT\n+\n!\n")
            hashes, reads, total_bytes = full.fastq_stats([first, second])
            self.assertEqual(reads, 2)
            self.assertEqual(total_bytes, first.stat().st_size + second.stat().st_size)
            self.assertEqual(
                hashes,
                [
                    hashlib.sha256(first.read_bytes()).hexdigest(),
                    hashlib.sha256(second.read_bytes()).hexdigest(),
                ],
            )

    def test_expected_fastq_records_respects_layout(self) -> None:
        self.assertEqual(full.expected_fastq_records("SINGLE", 11), 11)
        self.assertEqual(full.expected_fastq_records("PAIRED", 11), 22)
        with self.assertRaises(ValueError):
            full.expected_fastq_records("UNKNOWN", 11)

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

    def test_locked_sra_manifest_matches_full_scope_when_present(self) -> None:
        if not full.SRA_MANIFEST.is_file():
            self.skipTest("SRA transport manifest is generated at P3B start")
        rows = full.read_tsv(full.SRA_MANIFEST)
        by_run = full.validate_sra_manifest(rows)
        self.assertEqual(len(by_run), 482)

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

    def test_interrupted_curl_is_restarted_as_a_resuming_process(self) -> None:
        payload = b"resume-across-processes"
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "sample.fastq.gz"
            partial = Path(str(path) + ".part")

            class InterruptedRunner:
                def __init__(self) -> None:
                    self.calls = 0

                def run(self, command: list[str], *, stdout=None) -> None:
                    self.calls += 1
                    self_command = command
                    if self_command[self_command.index("--continue-at") + 1] != "-":
                        raise AssertionError(self_command)
                    if "--retry" in self_command:
                        raise AssertionError(self_command)
                    if self.calls == 1:
                        partial.write_bytes(payload[:8])
                        raise full.subprocess.CalledProcessError(18, command)
                    partial.write_bytes(payload)

            runner = InterruptedRunner()
            actual_sha, actual_md5 = full.download_fastq(
                runner,
                "example.invalid/sample.fastq.gz",
                hashlib.md5(payload).hexdigest(),
                len(payload),
                path,
            )
            self.assertEqual(runner.calls, 2)
            self.assertEqual(path.read_bytes(), payload)
            self.assertEqual(actual_sha, hashlib.sha256(payload).hexdigest())
            self.assertEqual(actual_md5, hashlib.md5(payload).hexdigest())


if __name__ == "__main__":
    unittest.main()
