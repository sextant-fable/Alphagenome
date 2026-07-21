from __future__ import annotations

import hashlib
import json
from pathlib import Path
import tempfile
import unittest
from unittest import mock

from scripts import evaluate_v2_model
from scripts import run_v2_p6c
from scripts import v2_bigwig_dataset


def write_json(path: Path, payload: dict[str, object]) -> None:
    path.write_text(json.dumps(payload) + "\n")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class FinalTestGateTest(unittest.TestCase):
    def test_loader_requires_claim_then_consumed_running_state(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            state = root / "state.json"
            lock = root / "lock.json"
            claim = root / "claim.json"
            write_json(
                state,
                {
                    "current_phase": "P6C",
                    "approvals": {
                        "G5_final_test": {
                            "approved": True,
                            "scope": "r6c_single_chr_x_test",
                        }
                    },
                },
            )
            write_json(
                lock,
                {
                    "checkpoint_sha256": "checkpoint",
                    "superseded": False,
                    "test_consumed": False,
                },
            )
            write_json(
                claim,
                {
                    "status": "claimed",
                    "execution_id": "execution",
                    "checkpoint_sha256": "checkpoint",
                },
            )
            with (
                mock.patch.object(v2_bigwig_dataset, "STATE_PATH", state),
                mock.patch.object(v2_bigwig_dataset, "FINAL_LOCK_PATH", lock),
                mock.patch.object(v2_bigwig_dataset, "FINAL_CLAIM_PATH", claim),
            ):
                v2_bigwig_dataset.require_final_test_access(
                    "checkpoint", "execution"
                )
                with self.assertRaisesRegex(PermissionError, "before.*consumed"):
                    v2_bigwig_dataset.require_final_test_read_access(
                        "checkpoint", "execution"
                    )
                current = json.loads(lock.read_text())
                current.update(
                    {
                        "test_consumed": True,
                        "test_status": "running",
                        "test_execution_id": "execution",
                    }
                )
                write_json(lock, current)
                v2_bigwig_dataset.require_final_test_read_access(
                    "checkpoint", "execution"
                )
                current["test_status"] = "completed"
                write_json(lock, current)
                with self.assertRaisesRegex(PermissionError, "before.*consumed"):
                    v2_bigwig_dataset.require_final_test_read_access(
                        "checkpoint", "execution"
                    )

    def test_consume_is_bound_to_claim_and_is_one_time(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            lock = root / "final_test_lock.json"
            claim = root / "final_test_claim.json"
            write_json(
                lock,
                {
                    "checkpoint_sha256": "checkpoint",
                    "superseded": False,
                    "test_consumed": False,
                },
            )
            write_json(
                claim,
                {
                    "status": "claimed",
                    "execution_id": "execution",
                    "checkpoint_sha256": "checkpoint",
                    "physical_gpu": 2,
                },
            )
            with (
                mock.patch.object(evaluate_v2_model, "METADATA_DIR", root),
                mock.patch.object(evaluate_v2_model, "FINAL_CLAIM_PATH", claim),
            ):
                evaluate_v2_model.consume_final_test("checkpoint", "execution", 2)
                consumed = json.loads(lock.read_text())
                self.assertIs(consumed["test_consumed"], True)
                self.assertEqual(consumed["test_status"], "running")
                self.assertEqual(consumed["test_execution_id"], "execution")
                self.assertEqual(consumed["test_physical_gpu"], 2)
                with self.assertRaisesRegex(RuntimeError, "already consumed"):
                    evaluate_v2_model.consume_final_test(
                        "checkpoint", "execution", 2
                    )

    def test_locked_input_preflight_binds_registered_files(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            metadata = root / "metadata"
            intervals = root / "intervals"
            metadata.mkdir()
            intervals.mkdir()
            checkpoint = root / "checkpoint.pt"
            selection = metadata / "selection.json"
            test_intervals = intervals / "test.tsv"
            means = metadata / "means.tsv"
            split_registry = metadata / "split.json"
            means_summary = metadata / "means_summary.json"
            checkpoint.write_bytes(b"checkpoint")
            selection.write_text("{}\n")
            test_intervals.write_text("chromosome\tstart\tend\nX\t0\t128\n")
            means.write_text("group_id\tdevelopment_I_V_nonzero_mean\nG1\t1\n")
            relative_test = str(test_intervals.relative_to(root))
            write_json(
                split_registry,
                {"files": {relative_test: sha256(test_intervals)}},
            )
            write_json(
                means_summary,
                {
                    "output_sha256": sha256(means),
                    "split_registry_sha256": sha256(split_registry),
                },
            )
            lock = {
                "checkpoint_path": str(checkpoint.relative_to(root)),
                "checkpoint_sha256": sha256(checkpoint),
                "selection_path": str(selection.relative_to(root)),
                "selection_sha256": sha256(selection),
            }
            with (
                mock.patch.object(run_v2_p6c, "REPO_ROOT", root),
                mock.patch.object(run_v2_p6c, "SPLIT_REGISTRY_PATH", split_registry),
                mock.patch.object(run_v2_p6c, "MEANS_PATH", means),
                mock.patch.object(run_v2_p6c, "MEANS_SUMMARY_PATH", means_summary),
                mock.patch.object(run_v2_p6c, "TEST_INTERVALS_PATH", test_intervals),
            ):
                observed = run_v2_p6c.validate_locked_inputs(lock)
                self.assertEqual(observed["checkpoint_sha256"], sha256(checkpoint))
                self.assertEqual(
                    observed["test_intervals_sha256"], sha256(test_intervals)
                )
                checkpoint.write_bytes(b"changed")
                with self.assertRaisesRegex(RuntimeError, "checkpoint SHA-256"):
                    run_v2_p6c.validate_locked_inputs(lock)


if __name__ == "__main__":
    unittest.main()
