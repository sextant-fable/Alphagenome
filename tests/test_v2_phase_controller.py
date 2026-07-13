from __future__ import annotations

import copy
import unittest

from scripts import run_v2_p3a
from scripts import v2_phase_controller as controller


class V2PhaseControllerTest(unittest.TestCase):
    def test_phase_sequence_separates_pilot_bulk_and_final_test(self) -> None:
        self.assertEqual(controller.next_phase("P2"), "P3A")
        self.assertEqual(controller.next_phase("P3A"), "P3B")
        self.assertEqual(controller.next_phase("P5"), "P6A")
        self.assertEqual(controller.next_phase("P6B"), "P6C")
        self.assertIsNone(controller.next_phase("P6C"))

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
            any(value.endswith("run_v2_p3b.py") for value in controller.PHASE_COMMANDS["P3B"])
        )
        self.assertTrue(
            any(value == "P3B" for value in controller.REVIEW_COMMANDS["P3B"])
        )

    def test_p4_loader_and_review_are_registered(self) -> None:
        self.assertTrue(
            any(value.endswith("run_v2_p4.py") for value in controller.PHASE_COMMANDS["P4"])
        )
        self.assertTrue(
            any(value == "P4" for value in controller.REVIEW_COMMANDS["P4"])
        )

    def test_final_test_scope_cannot_be_approved_early(self) -> None:
        with self.assertRaisesRegex(ValueError, "only in P6C"):
            controller.validate_approval_scope(
                "G5", "r6c_single_chr_x_test", "P3A"
            )
        controller.validate_approval_scope(
            "G5", "r6c_single_chr_x_test", "P6C"
        )

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
