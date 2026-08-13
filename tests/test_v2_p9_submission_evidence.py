from __future__ import annotations

import copy
import json
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import torch

from scripts import run_v2_p9_submission_evidence as p9
from scripts import review_v2_phase
from scripts import train_v2_model
from scripts import v2_training_components as components


class P9SubmissionEvidenceTest(unittest.TestCase):
    def setUp(self) -> None:
        self.spec = json.loads(p9.SPEC_PATH.read_text())

    def test_matrix_has_90_records_and_exactly_59_new_jobs(self) -> None:
        jobs = p9.matrix_jobs(self.spec)
        self.assertEqual(len(jobs), 90)
        self.assertEqual(len({job["job_id"] for job in jobs}), 90)
        self.assertEqual(sum(not job["base_job"] for job in jobs), 59)
        self.assertEqual(sum(job["base_job"] for job in jobs), 31)
        reused_counts = {
            configuration: sum(
                job["base_job"] and job["configuration"] == configuration
                for job in jobs
            )
            for configuration in {
                "frozen_no_worm_no_lora",
                "full_worm_lora",
                "no_lora",
            }
        }
        self.assertEqual(
            reused_counts,
            {
                "frozen_no_worm_no_lora": 15,
                "full_worm_lora": 15,
                "no_lora": 1,
            },
        )
        reused = [job for job in jobs if job["base_job"]]
        self.assertEqual(len(reused), 31)
        p8_reused = [job for job in reused if job["reused_source"] == "P8"]
        self.assertEqual(len(p8_reused), 1)
        self.assertEqual(
            (p8_reused[0]["configuration"], p8_reused[0]["fold"], p8_reused[0]["seed"]),
            ("no_lora", 1, 20260714),
        )
        self.assertEqual({job["fold"] for job in jobs}, {1, 2, 3, 4, 5})
        self.assertEqual(
            {job["seed"] for job in jobs}, {20260714, 20260715, 20260716}
        )
        self.assertEqual(
            {job["configuration"] for job in jobs if job["reused_source"] == "P6B"},
            {"frozen_no_worm_no_lora", "full_worm_lora"},
        )

    def test_new_commands_are_validation_only_and_never_request_final_test(self) -> None:
        for job in p9.matrix_jobs(self.spec):
            train, evaluate = p9.commands_for_job(job, self.spec)
            if job["base_job"]:
                self.assertIsNone(train)
                self.assertIsNone(evaluate)
                continue
            serialized = " ".join([*train, *evaluate])
            self.assertNotIn("--final-test", serialized)
            self.assertNotIn("test_locked.tsv", serialized)
            self.assertIn(f"fold_{job['fold']}_train_nonzero_mean", serialized)
            self.assertEqual(
                train[train.index("--hidden-channels") + 1],
                str(job["hidden_channels"]),
            )

    def test_lora_only_uses_original_embedding_and_applies_lora(self) -> None:
        class FakePart:
            gradient_checkpointing = False

        class FakeBase(torch.nn.Module):
            def __init__(self) -> None:
                super().__init__()
                self.anchor = torch.nn.Parameter(torch.ones(1))
                self.encoder = FakePart()
                self.tower = FakePart()
                self.decoder = FakePart()

        base = FakeBase()
        with (
            mock.patch.object(train_v2_model.AlphaGenome, "from_pretrained", return_value=base),
            mock.patch.object(components, "add_c_elegans_organism_embeddings") as add_embedding,
            mock.patch.object(train_v2_model, "apply_lora") as apply_lora,
        ):
            model = train_v2_model.build_model(
                "B_lora_only", torch.ones(3), torch.device("cpu"), 64
            )
        add_embedding.assert_not_called()
        apply_lora.assert_called_once()
        self.assertEqual(model.base_organism_index, 0)
        self.assertEqual(
            train_v2_model.adaptation_metadata("B_lora_only")["trunk_policy"],
            "lora_only",
        )

    def test_one_bp_head_derives_model_space_128bp_from_unscaled_sum(self) -> None:
        torch.manual_seed(9)
        means = torch.tensor([2.0, 4.0])
        experimental = torch.rand(1, 2, 256) * 5
        model_1bp = components.scale_targets_model_space(experimental, means, 1)
        derived = components.derive_128bp_model_space_from_1bp(model_1bp, means)
        expected_experimental = experimental.reshape(1, 2, 2, 128).sum(dim=-1)
        expected = components.scale_targets_model_space(expected_experimental, means, 128)
        torch.testing.assert_close(derived, expected, rtol=1e-5, atol=1e-5)

    def test_one_bp_output_head_keeps_full_b_adaptation_and_exact_parameter_count(self) -> None:
        metadata = train_v2_model.adaptation_metadata("B_1bp_head_only")
        self.assertTrue(metadata["c_elegans_organism_embedding"])
        self.assertTrue(metadata["lora_enabled"])
        self.assertEqual(metadata["learned_head_resolutions"], [1])
        self.assertEqual(metadata["derived_output_resolutions"], [128])
        head = components.DualResolutionRnaHead(
            241, torch.ones(241), resolutions=(1,)
        )
        total = sum(parameter.numel() for parameter in head.parameters())
        total += 18_816 + 145_920
        config = next(
            row
            for row in self.spec["matrix"]["configurations"]
            if row["id"] == "learned_1bp_head_only"
        )
        self.assertEqual(total, 535_394)
        self.assertEqual(total, config["expected_trainable_parameters"])

    def test_size_matched_baseline_parameter_contract(self) -> None:
        model = train_v2_model.build_model(
            "C_size_matched", torch.ones(241), torch.device("cpu"), 152
        )
        observed = train_v2_model.trainable_parameter_count(model)
        registered = self.spec["matrix"]["parameter_match_reference"]
        self.assertEqual(observed, registered["candidate_trainable_parameters"])
        relative = abs(observed - registered["model_b_trainable_parameters"]) / registered[
            "model_b_trainable_parameters"
        ]
        self.assertLessEqual(relative, registered["maximum_relative_difference"])
        with self.assertRaisesRegex(ValueError, "hidden_channels=152"):
            train_v2_model.build_model(
                "C_size_matched", torch.ones(241), torch.device("cpu"), 144
            )

    def test_paired_interval_uses_fold_then_seed_hierarchical_resampling(self) -> None:
        pairs = [
            {"fold": fold, "seed": seed, "difference": fold + seed / 1e9}
            for fold in range(1, 6)
            for seed in (20260714, 20260715, 20260716)
        ]
        with mock.patch.object(
            p9.submission_statistics,
            "make_hierarchical_resample_plan",
            wraps=p9.submission_statistics.make_hierarchical_resample_plan,
        ) as make_plan:
            first = p9.paired_bootstrap_interval(pairs, "metric", self.spec["analysis"])
        second = p9.paired_bootstrap_interval(pairs, "metric", self.spec["analysis"])
        make_plan.assert_called_once_with(
            n_blocks=5,
            n_repeats=3,
            iterations=10_000,
            seed=20260813,
        )
        self.assertEqual(first, second)
        self.assertLessEqual(first[1], first[0])
        self.assertGreaterEqual(first[2], first[0])

    def test_factorial_effect_formulas_are_exact(self) -> None:
        values = {
            "frozen_no_worm_no_lora": 1.0,
            "lora_only": 3.0,
            "no_lora": 4.0,
            "full_worm_lora": 10.0,
        }
        records = []
        for configuration, value in values.items():
            for fold in range(1, 6):
                for seed in (20260714, 20260715, 20260716):
                    row = {metric: value for metric in p9.METRIC_NAMES}
                    row.update({"configuration": configuration, "fold": fold, "seed": seed})
                    records.append(row)
        raw, summaries = p9.factorial_effect_rows(records, self.spec)
        selected = {
            row["effect"]: row["mean_effect"]
            for row in summaries
            if row["metric"] == "primary_biological_score"
        }
        self.assertEqual(selected["worm_embedding_main"], 5.0)
        self.assertEqual(selected["lora_main"], 4.0)
        self.assertEqual(selected["worm_lora_interaction"], 4.0)
        self.assertEqual(len(raw), len(p9.METRIC_NAMES) * 15 * 3)

    def test_completed_job_recomputes_metrics_and_requires_all_artifact_hashes(self) -> None:
        job = {"base_job": False, "job_path": "job.json"}
        old = {
            "status": "completed",
            "spec_sha256": "spec",
            "checkpoint_sha256": "checkpoint",
            "run_sha256": "run",
            "validation_sha256": "validation",
            "physical_gpu": 2,
            "elapsed_seconds": 7.5,
            "primary_biological_score": -999,
        }
        fresh = {
            "checkpoint_sha256": "checkpoint",
            "run_sha256": "run",
            "validation_sha256": "validation",
            "primary_biological_score": 0.5,
        }
        with (
            mock.patch.object(p9.Path, "is_file", return_value=True),
            mock.patch.object(p9, "read_json", return_value=old),
            mock.patch.object(p9, "_record_from_validation", return_value=fresh),
        ):
            result = p9.completed_job(job, "spec")
        self.assertEqual(result["primary_biological_score"], 0.5)
        self.assertEqual(result["physical_gpu"], 2)
        changed = dict(fresh, validation_sha256="changed")
        with (
            mock.patch.object(p9.Path, "is_file", return_value=True),
            mock.patch.object(p9, "read_json", return_value=old),
            mock.patch.object(p9, "_record_from_validation", return_value=changed),
        ):
            self.assertIsNone(p9.completed_job(job, "spec"))

    def test_contract_rejects_wrong_scope_before_inputs(self) -> None:
        state = {
            "current_phase": "P9",
            "status": "RUNNING",
            "approvals": {
                "G4_gpu_experiments": {
                    "approved": True,
                    "scope": "p8_b_no_lora_fold1",
                }
            },
        }
        with mock.patch.object(p9, "read_json", return_value=state):
            with self.assertRaisesRegex(RuntimeError, "contract mismatch"):
                p9.require_contract(copy.deepcopy(self.spec))

    def test_p10_review_is_registered_and_fails_closed_on_missing_outputs(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            metadata = root / "alphagenome_custom/metadata/v2"
            metadata.mkdir(parents=True)
            (metadata / "p10_dpy27_internal_application_spec.json").write_text(
                json.dumps(
                    {
                        "output_paths": {
                            key: f"missing/{key}"
                            for key in (
                                "analysis_audit",
                                "preflight_audit",
                                "gene_contrasts",
                                "run_endpoints",
                                "fold_endpoints",
                                "fold_gene_contrasts",
                                "summary_endpoints",
                                "figure_source_data",
                                "figure_pdf",
                                "figure_svg",
                                "figure_tiff",
                                "figure_png",
                                "figure_qa",
                            )
                        }
                    }
                )
            )
            with mock.patch.object(review_v2_phase, "REPO_ROOT", root):
                report = review_v2_phase.review_p10()
        self.assertEqual(report["status"], "FAIL")
        self.assertEqual(report["checks"][0]["check_id"], "R10.01_required_outputs")


if __name__ == "__main__":
    unittest.main()
