from __future__ import annotations

import csv
import json
import tempfile
from pathlib import Path
import unittest

import numpy as np

from scripts import v2_biological_validation as biology
from scripts import review_v2_phase
from scripts.run_dpy27_internal_application import (
    CONTROLLER_INVOCATION,
    REGISTERED_PHASE_COMMAND,
    _figure_source_rows,
    execution_provenance,
    implementation_sha256,
    preflight_contract,
)


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


class Dpy27InternalApplicationTests(unittest.TestCase):
    def test_validation_loader_is_fold_valid_only_and_nonoverlapping(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "valid.tsv"
            rows = [
                {
                    "chromosome": "I",
                    "start": 0,
                    "end": 256,
                    "core_start": 0,
                    "core_end": 256,
                    "fold": 1,
                    "role": "valid",
                    "block_id": "I_block_1",
                }
            ]
            write_tsv(path, rows)
            _, windows = biology.load_fold_validation_subwindows(
                path, expected_fold=1, sequence_length=128
            )
            self.assertEqual([(row.start, row.end) for row in windows], [(0, 128), (128, 256)])
            rows[0]["role"] = "test_locked"
            write_tsv(path, rows)
            with self.assertRaisesRegex(RuntimeError, "Only role=valid"):
                biology.load_fold_validation_subwindows(
                    path, expected_fold=1, sequence_length=128
                )

    def test_final_test_named_path_is_rejected_before_read(self) -> None:
        with self.assertRaisesRegex(RuntimeError, "Final-test path"):
            biology.load_fold_validation_subwindows(
                "/tmp/test_locked.tsv", expected_fold=1, sequence_length=128
            )

    def test_gene_exon_accumulator_reports_only_unique_exon_bases(self) -> None:
        gene = biology.AnnotatedGene(
            "g1", "gene-1", "I", "+", 0, 8, ((1, 4), (6, 8))
        )
        accumulator = biology.Dpy27GeneExonAccumulator([gene], {"I": [0]})
        prediction = np.vstack([np.arange(8), np.arange(8) + 1]).astype(float)
        observation = prediction + 2
        accumulator.update(
            "I", 0, prediction, observation, block_id="I_block_1"
        )
        rows = accumulator.rows(fold=1, seed=20260714)
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["exon_bases_evaluated"], 5)
        self.assertTrue(rows[0]["complete_exon_coverage"])
        np.testing.assert_allclose(rows[0]["predicted_dpy27_rnai_mean"], 3.8)

    def test_spearman_handles_ties_without_scipy(self) -> None:
        self.assertAlmostEqual(
            biology.spearman_correlation([1, 2, 2, 4], [10, 20, 20, 40]), 1.0
        )

    @staticmethod
    def _gene_rows() -> list[dict[str, object]]:
        rows: list[dict[str, object]] = []
        genes = (
            ("x1", "X", 0.5, (0.0, 3.0, 6.0)),
            ("x2", "X", 1.5, (1.0, 4.0, 8.0)),
            ("a1", "I", 0.2, (0.2, 0.4, 0.6)),
            ("a2", "II", 0.7, (0.3, 0.8, 1.2)),
        )
        for fold in range(1, 6):
            for seed_index, seed in enumerate((20260714, 20260715, 20260716)):
                for gene_id, chromosome, observed_dpy, seed_predictions in genes:
                    predicted_dpy = seed_predictions[seed_index] + fold * 0.01
                    rows.append(
                        {
                            "analysis_label": "internal_genomic_block_validation",
                            "fold": fold,
                            "seed": seed,
                            "gene_id": gene_id,
                            "gene_name": gene_id,
                            "chromosome": chromosome,
                            "strand": "+",
                            "block_id": f"{chromosome}_block_{fold}",
                            "exon_bases_evaluated": 10,
                            "exon_bases_total": 10,
                            "exon_coverage_fraction": 1.0,
                            "complete_exon_coverage": True,
                            "predicted_dpy27_rnai_mean": predicted_dpy,
                            "predicted_vector_rnai_mean": 0.0,
                            "observed_dpy27_rnai_mean": observed_dpy,
                            "observed_vector_rnai_mean": 0.0,
                            "predicted_log1p_contrast_g0005_minus_g0006": np.log1p(
                                predicted_dpy
                            ),
                            "observed_log1p_contrast_g0005_minus_g0006": np.log1p(
                                observed_dpy
                            ),
                        }
                    )
        return rows

    def test_fold_endpoints_are_computed_after_gene_condition_seed_mean(self) -> None:
        rows = self._gene_rows()
        fold_genes, fold_endpoints, summaries = biology.fold_primary_summary(
            rows, bootstrap_replicates=500, bootstrap_seed=7
        )
        fold_1_x1 = next(
            row for row in fold_genes if row["fold"] == 1 and row["gene_id"] == "x1"
        )
        expected = np.log1p(np.mean([0.01, 3.01, 6.01]))
        mean_of_seed_contrasts = np.mean(np.log1p([0.01, 3.01, 6.01]))
        self.assertAlmostEqual(
            fold_1_x1["predicted_log1p_contrast_g0005_minus_g0006"], expected
        )
        self.assertNotAlmostEqual(expected, mean_of_seed_contrasts)
        self.assertEqual(len(fold_endpoints), 5)
        self.assertEqual(len(summaries), 4)
        self.assertTrue(
            all(
                row["inference_unit"]
                == "fold_after_within_fold_gene_condition_seed_mean"
                for row in summaries
            )
        )

    def test_observations_must_be_identical_across_seeds(self) -> None:
        rows = self._gene_rows()
        rows[0]["observed_dpy27_rnai_mean"] = 999.0
        with self.assertRaisesRegex(RuntimeError, "Observed signals differ"):
            biology.fold_primary_summary(rows, bootstrap_replicates=100)

    def test_fold_primary_requires_finite_endpoint_in_all_five_folds(self) -> None:
        rows = self._gene_rows()
        for row in rows:
            if row["fold"] == 1 and row["gene_id"] == "x1":
                row["predicted_dpy27_rnai_mean"] = 1.0
                row["predicted_vector_rnai_mean"] = 1.0
            if row["fold"] == 1 and row["gene_id"] == "x2":
                row["predicted_dpy27_rnai_mean"] = 1.0
                row["predicted_vector_rnai_mean"] = 1.0
        with self.assertRaisesRegex(RuntimeError, "must be finite in all five folds"):
            biology.fold_primary_summary(rows, bootstrap_replicates=100)

    def test_execution_provenance_has_stable_controller_contract(self) -> None:
        provenance = execution_provenance()
        self.assertEqual(provenance["controller_invocation"], CONTROLLER_INVOCATION)
        self.assertEqual(
            provenance["registered_phase_command"], REGISTERED_PHASE_COMMAND
        )
        for field in (
            "git_commit",
            "hostname",
            "working_directory",
            "python_executable",
        ):
            self.assertTrue(provenance[field])
        self.assertTrue(review_v2_phase._p10_execution_provenance_ok(provenance))
        provenance["registered_phase_command"] = "python bypass.py"
        self.assertFalse(review_v2_phase._p10_execution_provenance_ok(provenance))

    def test_r10_independently_recomputes_bootstrap_summary(self) -> None:
        genes = self._gene_rows()
        _, fold_rows, summaries = biology.fold_primary_summary(
            genes, bootstrap_replicates=500, bootstrap_seed=7
        )
        statistics = {
            "endpoints": list(biology.PRIMARY_ENDPOINTS),
            "bootstrap_replicates": 500,
            "bootstrap_seed": 7,
        }
        self.assertTrue(
            review_v2_phase._p10_summary_recomputation_ok(
                genes, fold_rows, summaries, statistics
            )
        )
        tampered = [dict(row) for row in summaries]
        tampered[0]["ci_95_high"] = float(tampered[0]["ci_95_high"]) + 0.01
        self.assertFalse(
            review_v2_phase._p10_summary_recomputation_ok(
                genes, fold_rows, tampered, statistics
            )
        )

    def test_r10_rebuilds_exact_four_panel_source_data(self) -> None:
        genes = self._gene_rows()
        fold_genes, fold_rows, summaries = biology.fold_primary_summary(
            genes, bootstrap_replicates=500, bootstrap_seed=7
        )
        source = _figure_source_rows(fold_genes, fold_rows, summaries)
        contract = json.loads(
            Path(
                "alphagenome_custom/metadata/v2/p10_dpy27_internal_application_spec.json"
            ).read_text()
        )["figure_source_data_contract"]
        self.assertTrue(
            review_v2_phase._p10_source_data_ok(
                source, fold_genes, fold_rows, summaries, contract
            )
        )
        missing_c = [row for row in source if row["panel"] != "c_fold_agreement"]
        self.assertFalse(
            review_v2_phase._p10_source_data_ok(
                missing_c, fold_genes, fold_rows, summaries, contract
            )
        )

    def test_r10_implementation_hash_matches_runner_contract(self) -> None:
        spec = json.loads(
            Path(
                "alphagenome_custom/metadata/v2/p10_dpy27_internal_application_spec.json"
            ).read_text()
        )
        self.assertEqual(
            review_v2_phase._p10_implementation_sha256(), implementation_sha256()
        )
        self.assertEqual(spec["implementation_sha256"], implementation_sha256())

    def test_real_p10_preflight_checks_15_runs_without_model_or_bigwig_reads(self) -> None:
        spec = json.loads(
            Path(
                "alphagenome_custom/metadata/v2/p10_dpy27_internal_application_spec.json"
            ).read_text()
        )
        self.assertEqual(
            spec["figure_source_data_contract"]["panels"],
            [
                "a_gene_contrasts",
                "b_fold_x_minus_autosomes",
                "c_fold_agreement",
                "c_fold_primary_summary",
            ],
        )
        result = preflight_contract(enforce_controller=False)
        self.assertEqual(result["status"], "passed")
        self.assertEqual(result["checkpoint_count"], 15)
        self.assertEqual(result["dataset"]["track_count"], 241)
        self.assertEqual(result["model_or_bigwig_reads"], 0)
        self.assertEqual(
            [result["dataset"]["intervals"][fold]["subwindows"] for fold in range(1, 6)],
            [83] * 5,
        )


if __name__ == "__main__":
    unittest.main()
