from __future__ import annotations

import csv
from pathlib import Path
import tempfile
import unittest

import numpy as np
import pyBigWig

from scripts import build_v2_splits
from scripts.v2_bigwig_dataset import V2BigWigDataset


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


class V2BigWigDatasetTest(unittest.TestCase):
    def test_eval_cores_partition_chromosome_without_overlap(self) -> None:
        length = 3 * build_v2_splits.WINDOW_SIZE + 12345
        windows = build_v2_splits.make_windows(length)
        cores = build_v2_splits.eval_cores(windows, length)
        self.assertEqual(cores[0][0], 0)
        self.assertEqual(cores[-1][1], length)
        self.assertTrue(
            all(first[1] == second[0] for first, second in zip(cores, cores[1:]))
        )
        self.assertEqual(sum(end - start for start, end in cores), length)

    def test_loader_matches_bigwig_sum_pooling_and_gene_strands(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            sequence = ("ACGT" * 64).encode()
            fasta = root / "genome.fa"
            fasta.write_bytes(b">I\n" + sequence + b"\n")
            fai = root / "genome.fa.fai"
            fai.write_text("I\t256\t3\t256\t257\n")
            gtf = root / "annotation.gtf"
            gtf.write_text(
                'I\ttest\tgene\t1\t64\t.\t+\t.\tgene_id "plus";\n'
                'I\ttest\tgene\t129\t192\t.\t-\t.\tgene_id "minus";\n'
            )
            tracks = []
            for index, values in enumerate(
                (
                    np.arange(256, dtype=np.float32),
                    np.full(256, 2.0, dtype=np.float32),
                )
            ):
                path = root / f"track_{index}.bw"
                with pyBigWig.open(str(path), "w") as bigwig:
                    bigwig.addHeader([("I", 256)])
                    bigwig.addEntries("I", 0, values=values.tolist(), span=1, step=1)
                tracks.append(
                    {"group_id": f"G{index}", "output_path": str(path)}
                )
            track_manifest = root / "tracks.tsv"
            group_manifest = root / "groups.tsv"
            interval_manifest = root / "intervals.tsv"
            write_tsv(track_manifest, tracks)
            write_tsv(group_manifest, [{"group_id": "G0"}, {"group_id": "G1"}])
            write_tsv(
                interval_manifest,
                [
                    {
                        "chromosome": "I",
                        "start": 0,
                        "end": 256,
                        "core_start": 64,
                        "core_end": 192,
                        "fold": 1,
                        "role": "valid",
                    }
                ],
            )
            dataset = V2BigWigDataset(
                interval_manifest,
                track_manifest_path=track_manifest,
                group_manifest_path=group_manifest,
                fasta_path=fasta,
                fai_path=fai,
                gtf_path=gtf,
                max_io_workers=2,
            )
            item = dataset[0]
            np.testing.assert_array_equal(
                item["target_1bp"][0].numpy(), np.arange(256, dtype=np.float32)
            )
            np.testing.assert_array_equal(
                item["target_128bp"].numpy(),
                [[sum(range(128)), sum(range(128, 256))], [256, 256]],
            )
            self.assertEqual(tuple(item["dna_sequence"].shape), (4, 256))
            self.assertEqual(tuple(item["gene_mask"].shape), (2, 256))
            self.assertEqual(int(item["gene_mask"][0].sum()), 64)
            self.assertEqual(int(item["gene_mask"][1].sum()), 64)
            self.assertEqual(int(item["core_mask"].sum()), 128)
            dataset.close()

    def test_chromosome_x_is_refused_without_final_lock(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            intervals = Path(directory) / "test.tsv"
            write_tsv(
                intervals,
                [
                    {
                        "chromosome": "X",
                        "start": 0,
                        "end": 128,
                        "core_start": 0,
                        "core_end": 128,
                        "fold": "final",
                        "role": "test_locked",
                    }
                ],
            )
            with self.assertRaisesRegex(PermissionError, "embargoed"):
                V2BigWigDataset(intervals)


if __name__ == "__main__":
    unittest.main()
