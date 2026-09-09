#!/usr/bin/env python3
"""Run the CPU-only P11 biological-replicate holdout data phase."""

from __future__ import annotations

from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import subprocess
import sys

from scripts import build_v2_replicate_holdout as build


REPO_ROOT = Path(__file__).resolve().parents[1]
METADATA_ROOT = REPO_ROOT / "alphagenome_custom/metadata/v2/replicate_holdout_v1"
SPEC_PATH = METADATA_ROOT / "p11_replicate_holdout_spec.json"
EXECUTION_PATH = METADATA_ROOT / "p11_replicate_holdout_execution.json"
MEANS_PATH = METADATA_ROOT / "track_nonzero_means.tsv"
MEANS_SUMMARY_PATH = METADATA_ROOT / "track_nonzero_means_summary.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(8 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def main() -> None:
    if not SPEC_PATH.is_file():
        result = build.build(type("Args", (), {"prepare_only": False})())
        if result.get("status") != "complete":
            raise RuntimeError("P11 data builder did not complete")
    elif not EXECUTION_PATH.is_file():
        result = build.build(type("Args", (), {"prepare_only": False})())
        if result.get("status") != "complete":
            raise RuntimeError("P11 data builder did not complete")
    else:
        build.verify_source_members(build.read_tsv(build.SOURCE_MEMBERS))

    if not (METADATA_ROOT / "heldout_primary_track_manifest.tsv").is_file():
        role_command = [
            sys.executable,
            "-m",
            "scripts.derive_v2_replicate_holdout_role_manifests",
        ]
        role_completed = subprocess.run(role_command, cwd=REPO_ROOT, check=False)
        if role_completed.returncode != 0:
            raise SystemExit(role_completed.returncode)
    if not MEANS_PATH.is_file() or not MEANS_SUMMARY_PATH.is_file():
        command = [
            sys.executable,
            "-m",
            "scripts.compute_v2_replicate_holdout_means",
        ]
        completed = subprocess.run(command, cwd=REPO_ROOT, check=False)
        if completed.returncode != 0:
            raise SystemExit(completed.returncode)

    execution = json.loads(EXECUTION_PATH.read_text())
    means_summary = json.loads(MEANS_SUMMARY_PATH.read_text())
    if execution.get("phase") != "P11" or execution.get("locked_test_access") != "prohibited":
        raise RuntimeError("P11 execution contract is invalid")
    if execution.get("locked_test_signal_reads") != 0:
        raise RuntimeError("P11 recorded a locked-test signal read")
    if means_summary.get("locked_test_block_signal_reads") != 0:
        raise RuntimeError("P11 means computation recorded a locked-test read")
    execution["means_path"] = str(MEANS_PATH.relative_to(REPO_ROOT))
    execution["means_sha256"] = sha256(MEANS_PATH)
    execution["means_summary_path"] = str(MEANS_SUMMARY_PATH.relative_to(REPO_ROOT))
    execution["means_summary_sha256"] = sha256(MEANS_SUMMARY_PATH)
    execution["controller_phase"] = "P11"
    execution["completed_at"] = execution.get("completed_at") or utc_now()
    EXECUTION_PATH.write_text(json.dumps(execution, indent=2, sort_keys=True) + "\n")
    print(json.dumps(execution, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
