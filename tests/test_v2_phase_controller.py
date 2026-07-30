from __future__ import annotations

import copy
import unittest
from unittest import mock

from scripts import run_v2_p3a
from scripts import v2_subprocess
from scripts import v2_phase_controller as controller


class V2PhaseControllerTest(unittest.TestCase):
    def test_phase_sequence_separates_pilot_bulk_and_final_test(self) -> None:
        self.assertEqual(controller.next_phase("P2"), "P3A")
        self.assertEqual(controller.next_phase("P3A"), "P3B")
        self.assertEqual(controller.next_phase("P5"), "P6A")
        self.assertEqual(controller.next_phase("P6B"), "P6C")
        self.assertEqual(controller.next_phase("P6C"), "P7")
        self.assertIsNone(controller.next_phase("P7"))

    def test_p3a_requires_all_three_exact_scopes(self) -> None:
        state = controller.initial_state()
        state["current_phase"] = "P3A"
        self.assertEqual(
            controller.missing_approvals(state, "P3A"),
            [
                "G1:p3_full_rna_streaming_482",
                "G2:v2_manifest_candidate_hierarchy",
                "G3:p3_full_outputs_and_p4_loader",
            ],
        )
        for gate, scope in controller.PHASE_REQUIRED_APPROVALS["P3A"]:
            state["approvals"][controller.GATE_KEYS[gate]] = {
                "approved": True,
                "scope": scope,
                "note": "test",
                "updated_at": "test",
            }
        self.assertEqual(controller.missing_approvals(state, "P3A"), [])

    def test_pilot_scope_does_not_authorize_full_streaming(self) -> None:
        state = controller.initial_state()
        state["approvals"][controller.GATE_KEYS["G1"]] = {
            "approved": True,
            "scope": "p3a_five_run_pilot",
            "note": "test",
            "updated_at": "test",
        }
        self.assertFalse(
            controller.approval_satisfies(state, "G1", "p3_full_rna_streaming_482")
        )

    def test_p3b_full_processing_and_review_are_registered(self) -> None:
        self.assertIn("P3B", controller.PHASE_COMMANDS)
        self.assertIn("P3B", controller.REVIEW_COMMANDS)
        self.assertTrue(
            "scripts.run_v2_p3b" in controller.PHASE_COMMANDS["P3B"]
        )
        self.assertTrue(
            any(value == "P3B" for value in controller.REVIEW_COMMANDS["P3B"])
        )

    def test_p4_loader_and_review_are_registered(self) -> None:
        self.assertTrue(
            "scripts.run_v2_p4" in controller.PHASE_COMMANDS["P4"]
        )
        self.assertTrue(
            any(value == "P4" for value in controller.REVIEW_COMMANDS["P4"])
        )
        state = controller.initial_state()
        state["current_phase"] = "P4"
        self.assertEqual(
            controller.missing_approvals(state, "P4"),
            ["G3:p4_six_chromosome_block_split"],
        )

    def test_revised_p4_scope_can_be_recorded_for_p6c_migration(self) -> None:
        controller.validate_approval_scope(
            "G3", "p4_six_chromosome_block_split", "P6C"
        )
        controller.validate_approval_scope(
            "G3", "p4_six_chromosome_block_split", "P4"
        )

    def test_six_chromosome_migration_refuses_without_new_g3_scope(self) -> None:
        state = controller.initial_state()
        state["current_phase"] = "P6C"
        state["status"] = "APPROVAL_REQUIRED"
        with mock.patch.object(controller, "save_state"):
            self.assertEqual(controller.prepare_six_chromosome_revision(state), 2)
        self.assertEqual(state["current_phase"], "P6C")
        self.assertEqual(state["status"], "APPROVAL_REQUIRED")
        self.assertIn("G3:p4_six_chromosome_block_split", state["history"][-1]["detail"])

    def test_six_chromosome_migration_reopens_p4_after_archival(self) -> None:
        state = controller.initial_state()
        state["current_phase"] = "P6C"
        state["status"] = "APPROVAL_REQUIRED"
        state["approvals"][controller.GATE_KEYS["G3"]] = {
            "approved": True,
            "scope": "p4_six_chromosome_block_split",
            "note": "test",
            "updated_at": "test",
        }
        with (
            mock.patch.object(controller, "save_state"),
            mock.patch.object(controller, "run_command", return_value=0) as run,
        ):
            self.assertEqual(controller.prepare_six_chromosome_revision(state), 0)
        run.assert_called_once_with(controller.SIX_CHROMOSOME_MIGRATION_COMMAND)
        self.assertEqual(state["current_phase"], "P4")
        self.assertEqual(state["status"], "PENDING")
        self.assertEqual(state["history"][-1]["status"], "REOPENED")

    def test_p5_components_and_review_are_registered(self) -> None:
        self.assertTrue(
            "scripts.run_v2_p5" in controller.PHASE_COMMANDS["P5"]
        )
        self.assertTrue(
            any(value == "P5" for value in controller.REVIEW_COMMANDS["P5"])
        )

    def test_p6a_gpu_smoke_and_review_are_registered(self) -> None:
        self.assertTrue(
            "scripts.run_v2_p6a" in controller.PHASE_COMMANDS["P6A"]
        )
        self.assertTrue(
            any(value == "P6A" for value in controller.REVIEW_COMMANDS["P6A"])
        )

    def test_p6b_formal_matrix_and_review_are_registered(self) -> None:
        self.assertTrue(
            "scripts.run_v2_p6b_six_chromosome"
            in controller.PHASE_COMMANDS["P6B"]
        )
        self.assertTrue(
            any(value == "P6B" for value in controller.REVIEW_COMMANDS["P6B"])
        )

    def test_p6c_final_test_and_review_are_registered(self) -> None:
        self.assertTrue(
            "scripts.run_v2_p6c" in controller.PHASE_COMMANDS["P6C"]
        )
        self.assertTrue(
            any(value == "P6C" for value in controller.REVIEW_COMMANDS["P6C"])
        )

    def test_p7_is_registered_as_a_scoped_development_only_follow_up(self) -> None:
        self.assertTrue(
            "scripts.run_v2_p7_legacy_residual" in controller.PHASE_COMMANDS["P7"]
        )
        self.assertTrue(
            any(value == "P7" for value in controller.REVIEW_COMMANDS["P7"])
        )
        controller.validate_approval_scope(
            "G4", "p7_single_legacy_residual_fold1", "P7"
        )
        state = controller.initial_state()
        state["current_phase"] = "P7"
        self.assertEqual(
            controller.missing_approvals(state, "P7"),
            ["G4:p7_single_legacy_residual_fold1"],
        )

    def test_final_test_scope_cannot_be_approved_early(self) -> None:
        with self.assertRaisesRegex(ValueError, "only in P6C"):
            controller.validate_approval_scope(
                "G5", "r6c_single_six_chromosome_block_test", "P3A"
            )
        controller.validate_approval_scope(
            "G5", "r6c_single_six_chromosome_block_test", "P6C"
        )

    def test_old_chr_x_scope_does_not_authorize_revised_final_test(self) -> None:
        state = controller.initial_state()
        state["current_phase"] = "P6C"
        state["approvals"][controller.GATE_KEYS["G4"]] = {
            "approved": True,
            "scope": "r6_gpu_auto_available_2_3",
            "note": "test",
            "updated_at": "test",
        }
        state["approvals"][controller.GATE_KEYS["G5"]] = {
            "approved": True,
            "scope": "r6c_single_chr_x_test",
            "note": "historical scope",
            "updated_at": "test",
        }
        self.assertEqual(
            controller.missing_approvals(state, "P6C"),
            ["G5:r6c_single_six_chromosome_block_test"],
        )

    def test_all_controller_entries_use_repository_module_execution(self) -> None:
        for commands in (controller.PHASE_COMMANDS, controller.REVIEW_COMMANDS):
            for command in commands.values():
                self.assertEqual(command[1], "-m")
                self.assertTrue(command[2].startswith("scripts."))

    def test_subprocess_helper_uses_repository_module_execution(self) -> None:
        command = v2_subprocess.module_command(
            "train_v2_model", "--model", "A"
        )
        self.assertEqual(
            command[1:],
            ["-m", "scripts.train_v2_model", "--model", "A"],
        )
        with self.assertRaises(ValueError):
            v2_subprocess.module_command("scripts.train_v2_model")

    def test_schema_v1_migration_is_fail_closed(self) -> None:
        state = {
            "schema_version": 1,
            "workflow": "c_elegans_rna_seq_v2",
            "current_phase": "P3",
            "status": "APPROVAL_REQUIRED",
            "updated_at": "test",
            "approvals": {
                key: False for key in controller.GATE_KEYS.values()
            },
            "history": [],
        }
        migrated = copy.deepcopy(state)
        self.assertTrue(controller.migrate_state(migrated))
        self.assertEqual(migrated["current_phase"], "P3A")
        self.assertEqual(migrated["schema_version"], 2)
        self.assertEqual(controller.missing_approvals(migrated, "P3A"), [
            "G1:p3_full_rna_streaming_482",
            "G2:v2_manifest_candidate_hierarchy",
            "G3:p3_full_outputs_and_p4_loader",
        ])

    def test_p3a_command_plan_is_bounded_and_cpu_only(self) -> None:
        plan = run_v2_p3a.command_plan("conda")
        pilot_command = plan["pilot_command"]
        self.assertEqual(plan["environment"], "alphagenome")
        self.assertIn("--freeze-installed", plan["install_command"])
        self.assertNotIn("create", plan["install_command"])
        self.assertTrue(
            any(value.endswith("p3_pilot_sources.tsv") for value in pilot_command)
        )
        self.assertTrue(
            any(value.endswith("rna_seq_v2_normalized_pilot") for value in pilot_command)
        )
        self.assertNotIn("CUDA_VISIBLE_DEVICES", " ".join(pilot_command))
        self.assertEqual(plan["minimum_free_bytes"], 50 * 1024**3)

    def test_environment_plan_rejects_model_package_changes(self) -> None:
        with self.assertRaisesRegex(RuntimeError, "protected model packages"):
            run_v2_p3a.reject_protected_changes(
                {"actions": {"LINK": [{"name": "python"}]}}
            )
        run_v2_p3a.reject_protected_changes(
            {"actions": {"LINK": [{"name": "htslib"}]}}
        )


if __name__ == "__main__":
    unittest.main()
