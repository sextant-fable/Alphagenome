#!/usr/bin/env python3
"""Archive the chromosome-holdout P4-P6 evidence before split revision."""

from __future__ import annotations

from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import shutil
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[1]
METADATA_DIR = REPO_ROOT / "alphagenome_custom/metadata/v2"
INTERVAL_DIR = REPO_ROOT / "alphagenome_custom/intervals/v2"
ARCHIVE_DIR = METADATA_DIR / "chromosome_holdout_original"
ARCHIVE_INTERVAL_DIR = REPO_ROOT / "alphagenome_custom/intervals/v2_chromosome_holdout_original"
STATE_PATH = METADATA_DIR / "execution_state.json"
LOCK_PATH = METADATA_DIR / "final_test_lock.json"
MIGRATION_PATH = METADATA_DIR / "six_chromosome_revision_migration.json"
REVISION_ID = "six_chromosome_blocks_v1"


OVERWRITTEN_METADATA = (
    "split_registry_v2.json",
    "p4_loader_benchmark.json",
    "track_nonzero_means_v2.tsv",
    "track_nonzero_means_v2_summary.json",
    "model_specs_v2.json",
    "p5_component_audit.json",
    "p6a_execution.json",
    "final_test_lock.json",
)
PRESERVED_P6B = (
    "p6b_amendment_spec.json",
    "p6b_amendment_execution.json",
    "p6b_amendment_cv_results.tsv",
    "p6b_amendment_ablation_results.tsv",
    "p6b_amendment_selection.json",
    "p6b_legacy_v1_disposition.json",
)
OVERWRITTEN_AUDITS = tuple(
    METADATA_DIR / "audits" / phase / name
    for phase in ("P4", "P5", "P6A", "P6B")
    for name in ("review.json", "review.md")
)


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(8 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def atomic_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)


def copy_verified(source: Path, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists():
        if sha256(source) != sha256(destination):
            raise RuntimeError(f"Existing archive differs from source: {destination}")
        return
    temporary = destination.with_suffix(destination.suffix + ".tmp")
    shutil.copy2(source, temporary)
    if sha256(source) != sha256(temporary):
        raise RuntimeError(f"Archive copy hash mismatch: {source}")
    temporary.replace(destination)


def inventory(path: Path) -> dict[str, Any]:
    return {
        "path": str(path.relative_to(REPO_ROOT)),
        "bytes": path.stat().st_size,
        "sha256": sha256(path),
    }


def main() -> None:
    state = json.loads(STATE_PATH.read_text())
    lock = json.loads(LOCK_PATH.read_text())
    g3 = state.get("approvals", {}).get("G3_dataset_or_cache_generation", {})
    if MIGRATION_PATH.is_file():
        migration = json.loads(MIGRATION_PATH.read_text())
        if (
            migration.get("revision_id") != REVISION_ID
            or lock.get("superseded_for_revision") != REVISION_ID
            or lock.get("test_consumed") is not False
        ):
            raise RuntimeError("Existing six-chromosome migration is inconsistent")
        print(json.dumps(migration, indent=2, sort_keys=True))
        return
    g5 = state.get("approvals", {}).get("G5_final_test", {})
    if not (
        state.get("current_phase") == "P6C"
        and g3.get("approved") is True
        and g3.get("scope") == "p4_six_chromosome_block_split"
        and g5.get("approved") is not True
        and lock.get("test_consumed") is False
    ):
        raise RuntimeError(
            "Six-chromosome revision preparation requires exact G3, unconsumed P6C, "
            "and G5 unapproved"
        )
    claim_path = METADATA_DIR / "final_test_claim.json"
    if claim_path.exists():
        raise RuntimeError("A final-test claim exists; revision cannot proceed")

    required_metadata = [METADATA_DIR / name for name in OVERWRITTEN_METADATA]
    required_p6b = [METADATA_DIR / name for name in PRESERVED_P6B]
    required_intervals = sorted(path for path in INTERVAL_DIR.rglob("*") if path.is_file())
    missing = [
        str(path)
        for path in [*required_metadata, *required_p6b, *OVERWRITTEN_AUDITS]
        if not path.is_file()
    ]
    if missing or not required_intervals:
        raise RuntimeError(f"Cannot archive incomplete chromosome-holdout evidence: {missing}")

    source_lock_sha256 = sha256(LOCK_PATH)
    copy_verified(
        LOCK_PATH, ARCHIVE_DIR / "final_test_lock_pre_revision.json"
    )
    for source in required_metadata:
        copy_verified(source, ARCHIVE_DIR / source.name)
    for source in OVERWRITTEN_AUDITS:
        copy_verified(
            source,
            ARCHIVE_DIR / "audits" / source.relative_to(METADATA_DIR / "audits"),
        )
    for source in required_intervals:
        copy_verified(source, ARCHIVE_INTERVAL_DIR / source.relative_to(INTERVAL_DIR))

    archived_lock_path = ARCHIVE_DIR / "final_test_lock.json"
    archived_lock = json.loads(archived_lock_path.read_text())
    archived_lock.update(
        {
            "superseded": True,
            "superseded_at": utc_now(),
            "superseded_for_revision": REVISION_ID,
            "superseded_reason": (
                "User-required six-chromosome within-chromosome block split replaces "
                "the I-V chromosome CV and X-only final test before G5 consumption."
            ),
        }
    )
    atomic_json(archived_lock_path, archived_lock)
    atomic_json(LOCK_PATH, archived_lock)

    migration = {
        "schema_version": 1,
        "revision_id": REVISION_ID,
        "prepared_at": utc_now(),
        "source_phase": "P6C",
        "source_g5_approved": False,
        "source_test_consumed": False,
        "source_lock_sha256": source_lock_sha256,
        "source_lock_archive": inventory(
            ARCHIVE_DIR / "final_test_lock_pre_revision.json"
        ),
        "old_split_disposition": "superseded_for_six_chromosome_split",
        "archived_metadata": [inventory(ARCHIVE_DIR / path.name) for path in required_metadata],
        "archived_intervals": [
            inventory(ARCHIVE_INTERVAL_DIR / path.relative_to(INTERVAL_DIR))
            for path in required_intervals
        ],
        "archived_audits": [
            inventory(
                ARCHIVE_DIR
                / "audits"
                / path.relative_to(METADATA_DIR / "audits")
            )
            for path in OVERWRITTEN_AUDITS
        ],
        "preserved_p6b": [inventory(path) for path in required_p6b],
        "active_lock_sha256": sha256(LOCK_PATH),
    }
    atomic_json(MIGRATION_PATH, migration)
    print(json.dumps(migration, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
