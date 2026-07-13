from __future__ import annotations

from pathlib import Path
from unittest import mock
import tempfile
import unittest

import pyBigWig

from scripts import run_v2_reprocessing_pilot as pilot


class V2ReprocessingPilotTest(unittest.TestCase):
    def test_checked_in_pilot_manifest_is_exactly_bounded(self) -> None:
        pilot_rows = pilot.read_tsv(
            pilot.REPO_ROOT / "alphagenome_custom/metadata/v2/p3_pilot_sources.tsv"
        )
        sample_rows = pilot.read_tsv(pilot.SAMPLE_MANIFEST)
        sample_by_run = pilot.validate_pilot_manifest(pilot_rows, sample_rows)
        self.assertEqual(
            set(sample_by_run).intersection(pilot.EXPECTED_PILOT_RUNS),
            pilot.EXPECTED_PILOT_RUNS,
        )
        self.assertEqual(
            sum(
                int(value)
                for row in pilot_rows
                for value in row["fastq_bytes"].split(";")
            ),
            6_221_459_671,
        )

    def test_normalize_bedgraph_targets_one_hundred_million(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory) / "source.bedGraph"
            output = Path(directory) / "normalized.bedGraph"
            source.write_text("I\t0\t10\t2\nI\t10\t20\t3\n")
            scale, raw_total = pilot.normalize_bedgraph(source, output)
            self.assertEqual(raw_total, 50.0)
            self.assertEqual(scale, 2_000_000.0)
            values = [float(line.split("\t")[3]) for line in output.read_text().splitlines()]
            self.assertEqual(values, [4_000_000.0, 6_000_000.0])

    def test_normalize_bedgraph_rejects_negative_signal(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory) / "source.bedGraph"
            output = Path(directory) / "normalized.bedGraph"
            source.write_text("I\t0\t10\t-1\n")
            with self.assertRaisesRegex(RuntimeError, "Invalid bedGraph value"):
                pilot.normalize_bedgraph(source, output)

    def test_bigwig_readback_accepts_locked_total(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "normalized.bw"
            with pyBigWig.open(str(path), "w") as bigwig:
                bigwig.addHeader([("I", 10)])
                bigwig.addEntries(["I"], [0], ends=[10], values=[10_000_000.0])
            summary = pilot.validate_output_bigwig(path, {"I": 10})
            self.assertEqual(summary["sumData"], 100_000_000.0)
            self.assertEqual(summary["relative_total_error"], 0.0)

    def test_download_resumes_existing_partial_file(self) -> None:
        payload = b"complete-fastq"
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "sample.fastq.gz"
            partial = Path(str(path) + ".part")
            partial.write_bytes(payload[:5])

            def complete_download(command: list[str]) -> None:
                self.assertEqual(command[command.index("--continue-at") + 1], "-")
                self.assertEqual(command[command.index("--retry") + 1], "20")
                partial.write_bytes(payload)

            with mock.patch.object(pilot, "run", side_effect=complete_download):
                pilot.download_fastq(
                    "example.invalid/sample.fastq.gz",
                    pilot.hashlib.md5(payload).hexdigest(),
                    len(payload),
                    path,
                )
            self.assertEqual(path.read_bytes(), payload)
            self.assertFalse(partial.exists())


if __name__ == "__main__":
    unittest.main()
