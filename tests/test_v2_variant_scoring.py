from __future__ import annotations

import csv
import gzip
import hashlib
import json
import tempfile
from pathlib import Path
import unittest
from unittest.mock import patch

import numpy as np

from scripts import v2_biological_validation as biology
from scripts import score_v2_variants as scorer
from scripts.score_v2_variants import validate_cuda_authorization, validate_variant_spec
from scripts import v2_variant_scoring as variants


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_fasta(path: Path, name: str, sequence: str, line_bases: int = 64) -> Path:
    lines = [sequence[index : index + line_bases] for index in range(0, len(sequence), line_bases)]
    path.write_text(f">{name}\n" + "\n".join(lines) + "\n")
    fai = path.with_suffix(".fa.fai")
    fai.write_text(f"{name}\t{len(sequence)}\t{len(name) + 2}\t{line_bases}\t{line_bases + 1}\n")
    return fai


class SyntheticModelBPredictor:
    def __init__(self, spec, paths) -> None:
        self.track_ids = [f"G{index:04d}" for index in range(1, 242)]

    def __call__(self, batch: np.ndarray) -> dict[int, np.ndarray]:
        gc_signal = batch[:, 1] + batch[:, 2]
        fine = np.repeat(gc_signal[:, None, :], len(self.track_ids), axis=1)
        coarse = fine.reshape(fine.shape[0], fine.shape[1], -1, 128).sum(axis=-1)
        return {1: fine.astype(np.float32), 128: coarse.astype(np.float32)}


class FailingSyntheticModelBPredictor(SyntheticModelBPredictor):
    def __call__(self, batch: np.ndarray) -> dict[int, np.ndarray]:
        raise RuntimeError("synthetic inference failure")


def write_synthetic_scoring_spec(root: Path) -> tuple[Path, dict[str, Path]]:
    sequence = "ACGT" * 128
    fasta = root / "genome.fa"
    fai = write_fasta(fasta, "I", sequence)
    vcf = root / "input.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        f"I\t129\tsnv\t{sequence[128]}\tC\t.\tPASS\t.\n"
        f"I\t133\tins\t{sequence[132]}\t{sequence[132]}GG\t.\tPASS\t.\n"
    )
    checkpoint = root / "checkpoint.pt"
    checkpoint.write_bytes(b"checkpoint-placeholder")
    model_weights = root / "weights.safetensors"
    model_weights.write_bytes(b"weights-placeholder")
    means = root / "means.tsv"
    means.write_text(
        "group_id\tfold_1_train_nonzero_mean\n"
        + "".join(f"G{index:04d}\t1.0\n" for index in range(1, 242))
    )
    tracks = root / "tracks.tsv"
    tracks.write_text(
        "group_id\toutput_path\n"
        + "".join(f"G{index:04d}\tunused_{index}.bw\n" for index in range(1, 242))
    )
    groups = root / "groups.tsv"
    groups.write_text(
        "group_id\tstrand\n"
        + "".join(f"G{index:04d}\t.\n" for index in range(1, 242))
    )
    gtf = root / "genes.gtf"
    gtf.write_text(
        'I\ttest\tgene\t101\t160\t.\t+\t.\tgene_id "g1"; gene_name "g1";\n'
        'I\ttest\texon\t111\t150\t.\t+\t.\tgene_id "g1"; gene_name "g1";\n'
    )
    model_specs = root / "model_specs.json"
    model_specs.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "n_tracks": 241,
                "resolutions": [1, 128],
                "model_b_organism_index": 2,
                "model_b_lora_rank": 8,
                "model_b_lora_alpha": 16,
                "model_b_lora_targets": [
                    "tower.blocks.8.mha",
                    "tower.blocks.8.mlp",
                ],
                "models": [
                    {
                        "model_id": "B",
                        "base_organism_index": 2,
                        "trunk_policy": "worm_embeddings_and_lora_only",
                    }
                ],
            }
        )
        + "\n"
    )
    run = root / "run.json"
    run.write_text(
        json.dumps(
            {
                "model": "B",
                "fold": 1,
                "seed": 7,
                "loss": "paper",
                "sequence_length": 128,
                "hidden_channels": 64,
                "mean_column": "fold_1_train_nonzero_mean",
                "checkpoint_sha256": sha256(checkpoint),
                "checkpoint_reload_verified": True,
            }
        )
        + "\n"
    )
    locked = {
        "checkpoint": checkpoint,
        "checkpoint_run": run,
        "fasta": fasta,
        "fai": fai,
        "gtf": gtf,
        "group_manifest": groups,
        "means": means,
        "track_manifest": tracks,
        "model_specs": model_specs,
        "model_weights": model_weights,
        "vcf": vcf,
    }
    output_dir = root / "output"
    payload = {
        "schema_version": 1,
        "contract": variants.VARIANT_SCORING_CONTRACT,
        "final_test_access": "prohibited",
        "device": "cpu",
        "model": {
            "model_id": "B",
            "loss": "paper",
            "fold": 1,
            "seed": 7,
            "sequence_length": 128,
            "hidden_channels": 64,
            "mean_column": "fold_1_train_nonzero_mean",
        },
        "inference": {
            "forward_reverse_complement_ensemble": True,
            "allele_window_policy": "variant_start_anchored_right_context_crop_or_extend",
            "track_semantics": "all_241_tracks_unstranded_no_rc_permutation",
            "maximum_alleles": 10,
            "require_pass": True,
        },
        "ism": {
            "enabled": True,
            "radius_bp": 0,
            "max_mutants_per_variant": 3,
        },
        "locked_inputs": {key: str(path) for key, path in locked.items()},
        "input_sha256": {key: sha256(path) for key, path in locked.items()},
        "output_dir": str(output_dir),
    }
    spec_path = root / "spec.json"
    spec_path.write_text(json.dumps(payload) + "\n")
    return spec_path, {**locked, "output_dir": output_dir}


class V2VariantScoringTests(unittest.TestCase):
    def test_spec_rejects_unknown_device_before_model_load(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "spec.json"
            path.write_text(
                json.dumps(
                    {
                        "schema_version": 1,
                        "contract": variants.VARIANT_SCORING_CONTRACT,
                        "final_test_access": "prohibited",
                        "device": "tpu",
                    }
                )
                + "\n"
            )
            with self.assertRaisesRegex(RuntimeError, "frozen contract"):
                validate_variant_spec(path)

    def test_vcf_parser_splits_multiallelic_records_and_reads_gzip(self) -> None:
        text = (
            "##fileformat=VCFv4.2\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "I\t3\trs1\tA\tC,G\t.\tPASS\t.\n"
        )
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "test.vcf.gz"
            with gzip.open(path, "wt") as handle:
                handle.write(text)
            records = list(variants.parse_vcf(path))
        self.assertEqual([record.alternate for record in records], ["C", "G"])
        self.assertEqual([record.alternate_index for record in records], [1, 2])

    def test_fixed_windows_validate_reference_and_handle_deletion(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            fasta = Path(directory) / "genome.fa"
            sequence = "ACGT" * 100
            fai_path = write_fasta(fasta, "I", sequence)
            fai = variants.read_fai(fai_path)
            allele = variants.VcfAllele(
                "I", 101, "del", sequence[100:103], sequence[100], ".", "PASS", ".", 3, 1
            )
            window = variants.build_allele_window(
                allele, fasta_path=fasta, fai=fai, sequence_length=64, anchor_index=20
            )
            self.assertEqual(len(window.reference_sequence), 64)
            self.assertEqual(len(window.alternate_sequence), 64)
            self.assertEqual(window.reference_sequence[20:23], allele.reference)
            self.assertEqual(window.alternate_sequence[20], allele.alternate)
            self.assertEqual(
                window.alignment_policy, "fixed_window_index_aligned_diagnostic"
            )
            self.assertFalse(variants.supports_coordinate_gene_aggregation(allele))
            self.assertEqual(
                variants.aggregate_variant_gene_exon_delta(
                    window,
                    np.zeros((1, 64)),
                    genes=[],
                    track_ids=["track"],
                ),
                [],
            )
            mismatch = variants.VcfAllele("I", 101, "bad", "AAA", "A", ".", "PASS", ".", 4, 1)
            with self.assertRaisesRegex(ValueError, "Reference allele mismatch"):
                variants.build_allele_window(
                    mismatch, fasta_path=fasta, fai=fai, sequence_length=64, anchor_index=20
                )

    def test_vcf_rejects_n_in_ref_or_alt_alleles(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            for reference, alternate, message in (
                ("N", "A", "Unsupported REF"),
                ("A", "N", "Unsupported ALT"),
            ):
                with self.subTest(reference=reference, alternate=alternate):
                    path = Path(directory) / f"{reference}_{alternate}.vcf"
                    path.write_text(
                        "##fileformat=VCFv4.2\n"
                        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                        f"I\t3\tv1\t{reference}\t{alternate}\t.\tPASS\t.\n"
                    )
                    with self.assertRaisesRegex(ValueError, message):
                        list(variants.parse_vcf(path))

    def test_forward_reverse_ensemble_restores_output_orientation(self) -> None:
        calls: list[np.ndarray] = []

        def predictor(batch: np.ndarray):
            calls.append(batch.copy())
            fine = batch[:, 0:1, :] + 2 * batch[:, 3:4, :]
            return {1: fine, 128: fine.sum(axis=-1, keepdims=True)}

        result = variants.forward_reverse_complement_ensemble(predictor, ["AAAC"])
        np.testing.assert_allclose(result[1][0, 0], [1.5, 1.5, 1.5, 0.0])
        self.assertEqual(set(result), {1, 128})
        self.assertEqual(len(calls), 2)

        with self.assertRaisesRegex(ValueError, "exactly the 1-bp and 128-bp"):
            variants.forward_reverse_complement_ensemble(
                lambda batch: {1: batch[:, :1, :]}, ["AAAC"]
            )
        with self.assertRaisesRegex(ValueError, "exactly the 1-bp and 128-bp"):
            variants.forward_reverse_complement_ensemble(
                lambda batch: {
                    1: batch[:, :1, :],
                    128: batch[:, :1, :1],
                    256: batch[:, :1, :1],
                },
                ["AAAC"],
            )

    def test_gene_and_exon_delta_aggregation(self) -> None:
        gene = biology.AnnotatedGene("g1", "gene-1", "I", "+", 2, 8, ((2, 4), (6, 8)))
        delta = np.array([[0, 0, 1, 1, 2, 2, 3, 3]], dtype=float)
        rows = variants.aggregate_gene_exon_delta(
            delta,
            chromosome="I",
            window_start=0,
            genes=[gene],
            track_ids=["track"],
            resolution=1,
        )
        self.assertEqual(rows[0]["gene_body_bases"], 6)
        self.assertEqual(rows[0]["exon_bases"], 4)
        self.assertAlmostEqual(rows[0]["exon_delta_sum"], 8.0)

        allele = variants.VcfAllele(
            "I", 4, "snv", "T", "A", ".", "PASS", ".", 3, 1
        )
        window = variants.AlleleWindow(
            allele=allele,
            start=0,
            end=8,
            anchor_index=3,
            reference_sequence="ACGTACGT",
            alternate_sequence="ACGAACGT",
            alignment_policy="coordinate_aligned_substitution",
        )
        snv_rows = variants.aggregate_variant_gene_exon_delta(
            window,
            delta,
            genes=[gene],
            track_ids=["track"],
        )
        self.assertEqual(len(snv_rows), 1)
        self.assertAlmostEqual(snv_rows[0]["exon_delta_sum"], 8.0)

    def test_ism_enumerates_three_substitutions_and_enforces_cap(self) -> None:
        mutations = list(
            variants.enumerate_ism_mutations("AN", window_start=10, positions=[0, 1])
        )
        self.assertEqual(len(mutations), 3)
        self.assertEqual({mutation.alternate for mutation in mutations}, {"C", "G", "T"})

        def predictor(batch: np.ndarray):
            fine = batch[:, :1, :]
            return {1: fine, 128: fine.sum(axis=-1, keepdims=True)}

        with self.assertRaisesRegex(RuntimeError, "max_mutants"):
            list(
                variants.score_ism(
                    predictor,
                    reference_sequence="AA",
                    window_start=0,
                    positions=[0, 1],
                    max_mutants=2,
                )
            )

    def test_run_scoring_synthetic_end_to_end(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            spec_path, paths = write_synthetic_scoring_spec(Path(directory))
            with patch.object(scorer, "ModelBPredictor", SyntheticModelBPredictor):
                audit = scorer.run_scoring(spec_path)

            self.assertEqual(audit["status"], "completed")
            self.assertEqual(audit["model_loads"], 1)
            self.assertEqual(audit["variant_count"], 2)
            self.assertEqual(audit["variant_track_rows"], 964)
            self.assertEqual(audit["variant_gene_exon_rows"], 482)
            self.assertEqual(audit["ism_track_rows"], 2892)
            self.assertEqual(audit["skipped_indel_gene_aggregation"], 1)
            for record in audit["outputs"].values():
                self.assertEqual(sha256(Path(record["path"])), record["sha256"])

            with Path(audit["outputs"]["variant_track_deltas"]["path"]).open(
                newline=""
            ) as handle:
                track_rows = list(csv.DictReader(handle, delimiter="\t"))
            with Path(audit["outputs"]["variant_gene_exon_deltas"]["path"]).open(
                newline=""
            ) as handle:
                gene_rows = list(csv.DictReader(handle, delimiter="\t"))
            with Path(audit["outputs"]["ism_track_deltas"]["path"]).open(
                newline=""
            ) as handle:
                ism_rows = list(csv.DictReader(handle, delimiter="\t"))

            insertion_key = "I:133:A>AGG:alt1"
            insertion_track_rows = [
                row for row in track_rows if row["variant_key"] == insertion_key
            ]
            self.assertTrue(insertion_track_rows)
            self.assertEqual(
                {row["result_scope"] for row in insertion_track_rows},
                {"fixed_window_index_aligned_diagnostic"},
            )
            self.assertFalse(
                any(row["variant_key"] == insertion_key for row in gene_rows)
            )

            snv_key = "I:129:A>C:alt1"
            snv_gene_rows = {
                int(row["resolution"]): row
                for row in gene_rows
                if row["variant_key"] == snv_key and row["track_id"] == "G0001"
            }
            self.assertAlmostEqual(
                float(snv_gene_rows[1]["gene_body_delta_sum"]), 1.0
            )
            self.assertAlmostEqual(
                float(snv_gene_rows[128]["gene_body_delta_sum"]), 60.0 / 128.0
            )

            snv_ism = {
                row["alternate"]: float(row["delta_sum"])
                for row in ism_rows
                if row["parent_variant_key"] == snv_key
                and row["track_id"] == "G0001"
                and row["resolution"] == "1"
            }
            self.assertEqual(set(snv_ism), {"C", "G", "T"})
            self.assertAlmostEqual(snv_ism["C"], 1.0)
            self.assertAlmostEqual(snv_ism["G"], 1.0)
            self.assertAlmostEqual(snv_ism["T"], 0.0)
            persisted = scorer.read_json(paths["output_dir"] / "scoring_audit.json")
            self.assertEqual(persisted["status"], "completed")
            self.assertEqual(persisted["model_loads"], 1)

    def test_preflight_rejects_input_hash_mismatch(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            spec_path, paths = write_synthetic_scoring_spec(Path(directory))
            paths["vcf"].write_text(paths["vcf"].read_text() + "# changed\n")
            with self.assertRaisesRegex(RuntimeError, "input hash mismatch: vcf"):
                validate_variant_spec(spec_path)

    def test_failed_scoring_audit_records_loaded_model(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            spec_path, paths = write_synthetic_scoring_spec(Path(directory))
            with patch.object(
                scorer, "ModelBPredictor", FailingSyntheticModelBPredictor
            ):
                with self.assertRaisesRegex(RuntimeError, "synthetic inference failure"):
                    scorer.run_scoring(spec_path)
            failed = scorer.read_json(paths["output_dir"] / "scoring_audit.json")
            self.assertEqual(failed["status"], "failed")
            self.assertEqual(failed["model_loads"], 1)
            self.assertEqual(failed["failures"], ["synthetic inference failure"])

    def test_full_preflight_locks_hashes_and_validates_vcf_reference_without_model(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            fasta = root / "genome.fa"
            sequence = "ACGT" * 100
            fai = write_fasta(fasta, "I", sequence)
            vcf = root / "input.vcf"
            vcf.write_text(
                "##fileformat=VCFv4.2\n"
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                f"I\t129\tv1\t{sequence[128]}\tC\t.\tPASS\t.\n"
                f"I\t133\tv2\t{sequence[132:134]}\t{sequence[132]}\t.\tPASS\t.\n"
            )
            checkpoint = root / "checkpoint.pt"
            checkpoint.write_bytes(b"checkpoint-placeholder")
            model_weights = root / "weights.safetensors"
            model_weights.write_bytes(b"weights-placeholder")
            means = root / "means.tsv"
            means.write_text(
                "group_id\tfold_1_train_nonzero_mean\n"
                + "".join(f"G{index:04d}\t1.0\n" for index in range(1, 242))
            )
            tracks = root / "tracks.tsv"
            tracks.write_text(
                "group_id\toutput_path\n"
                + "".join(f"G{index:04d}\tunused_{index}.bw\n" for index in range(1, 242))
            )
            groups = root / "groups.tsv"
            groups.write_text(
                "group_id\tstrand\n"
                + "".join(f"G{index:04d}\t.\n" for index in range(1, 242))
            )
            gtf = root / "genes.gtf"
            gtf.write_text(
                'I\ttest\tgene\t100\t140\t.\t+\t.\tgene_id "g1"; gene_name "g1";\n'
                'I\ttest\texon\t110\t120\t.\t+\t.\tgene_id "g1"; gene_name "g1";\n'
            )
            model_specs = root / "model_specs.json"
            model_specs.write_text(
                json.dumps(
                    {
                        "schema_version": 1,
                        "n_tracks": 241,
                        "resolutions": [1, 128],
                        "model_b_organism_index": 2,
                        "model_b_lora_rank": 8,
                        "model_b_lora_alpha": 16,
                        "model_b_lora_targets": [
                            "tower.blocks.8.mha",
                            "tower.blocks.8.mlp",
                        ],
                        "models": [
                            {
                                "model_id": "B",
                                "base_organism_index": 2,
                                "trunk_policy": "worm_embeddings_and_lora_only",
                            }
                        ],
                    }
                )
                + "\n"
            )
            run = root / "run.json"
            run_payload = {
                "model": "B",
                "fold": 1,
                "seed": 7,
                "loss": "paper",
                "sequence_length": 128,
                "hidden_channels": 64,
                "mean_column": "fold_1_train_nonzero_mean",
                "checkpoint_sha256": sha256(checkpoint),
                "checkpoint_reload_verified": True,
            }
            run.write_text(json.dumps(run_payload) + "\n")
            locked = {
                "checkpoint": str(checkpoint),
                "checkpoint_run": str(run),
                "fasta": str(fasta),
                "fai": str(fai),
                "gtf": str(gtf),
                "group_manifest": str(groups),
                "means": str(means),
                "track_manifest": str(tracks),
                "model_specs": str(model_specs),
                "model_weights": str(model_weights),
                "vcf": str(vcf),
            }
            spec = {
                "schema_version": 1,
                "contract": variants.VARIANT_SCORING_CONTRACT,
                "final_test_access": "prohibited",
                "device": "cpu",
                "model": {
                    "model_id": "B",
                    "loss": "paper",
                    "fold": 1,
                    "seed": 7,
                    "sequence_length": 128,
                    "hidden_channels": 64,
                    "mean_column": "fold_1_train_nonzero_mean",
                },
                "inference": {
                    "forward_reverse_complement_ensemble": True,
                    "allele_window_policy": "variant_start_anchored_right_context_crop_or_extend",
                    "track_semantics": "all_241_tracks_unstranded_no_rc_permutation",
                    "maximum_alleles": 10,
                    "require_pass": True,
                },
                "ism": {"enabled": False},
                "locked_inputs": locked,
                "input_sha256": {key: sha256(Path(value)) for key, value in locked.items()},
                "output_dir": str(root / "output"),
            }
            spec_path = root / "spec.json"
            spec_path.write_text(json.dumps(spec) + "\n")
            _, _, windows, audit = validate_variant_spec(spec_path)
        self.assertEqual(len(windows), 2)
        self.assertEqual(audit["status"], "passed")
        self.assertEqual(audit["model_loads"], 0)
        self.assertEqual(audit["final_test_access"], "prohibited")
        self.assertIsNone(audit["cuda_authorization"])
        self.assertEqual(audit["coordinate_aligned_gene_aggregation_alleles"], 1)
        self.assertEqual(audit["skipped_indel_gene_aggregation"], 1)
        self.assertTrue(audit["git_commit"])
        self.assertTrue(audit["hostname"])
        self.assertEqual(
            set(audit["implementation_sha256"]),
            {
                "scripts/score_v2_variants.py",
                "scripts/v2_variant_scoring.py",
                "scripts/v2_biological_validation.py",
                "scripts/train_v2_model.py",
                "scripts/v2_training_components.py",
            },
        )

    def test_cuda_requires_hash_bound_authorization(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            checkpoint = root / "checkpoint.pt"
            checkpoint.write_bytes(b"checkpoint")
            vcf = root / "input.vcf"
            vcf.write_text("placeholder")
            spec_path = root / "spec.json"
            spec_path.write_text("{}\n")
            spec = {
                "device": "cuda",
                "physical_gpu": 2,
                "execution_id": "candidate-1",
            }
            with self.assertRaisesRegex(RuntimeError, "cuda_authorization_path"):
                validate_cuda_authorization(
                    spec_path, spec, {"checkpoint": checkpoint, "vcf": vcf}
                )

            authorization = root / "authorization.json"
            spec["cuda_authorization_path"] = str(authorization)
            payload = {
                "schema_version": 1,
                "approved": True,
                "scope": "variant_scoring:candidate-1",
                "execution_id": "candidate-1",
                "spec_sha256": sha256(spec_path),
                "physical_gpu": 2,
                "checkpoint_sha256": sha256(checkpoint),
                "vcf_sha256": "wrong-hash",
            }
            authorization.write_text(json.dumps(payload) + "\n")
            with self.assertRaisesRegex(RuntimeError, "authorization mismatch"):
                validate_cuda_authorization(
                    spec_path, spec, {"checkpoint": checkpoint, "vcf": vcf}
                )
            payload["vcf_sha256"] = sha256(vcf)
            authorization.write_text(json.dumps(payload) + "\n")
            accepted = validate_cuda_authorization(
                spec_path, spec, {"checkpoint": checkpoint, "vcf": vcf}
            )
            self.assertEqual(accepted["scope"], "variant_scoring:candidate-1")


if __name__ == "__main__":
    unittest.main()
