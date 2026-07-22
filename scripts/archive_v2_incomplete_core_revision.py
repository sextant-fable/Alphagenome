#!/usr/bin/env python3
"""Archive the invalid six-chromosome v1 core geometry before P4 rebuild."""

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
ARCHIVE_DIR = METADATA_DIR / "six_chromosome_blocks_v1_incomplete_core"
ARCHIVE_INTERVAL_DIR = (
    REPO_ROOT / "alphagenome_custom/intervals/v2_six_chromosome_blocks_v1_incomplete_core"
)
MIGRATION_PATH = METADATA_DIR / "six_chromosome_core_coverage_revision.json"
REGISTRY_PATH = METADATA_DIR / "split_registry_v2.json"
EXECUTION_PATH = METADATA_DIR / "p6b_six_chromosome_execution.json"

CURRENT_METADATA = (
    "split_registry_v2.json",
    "p4_loader_benchmark.json",
    "track_nonzero_means_v2.tsv",
    "track_nonzero_means_v2_summary.json",
    "model_specs_v2.json",
    "p5_component_audit.json",
    "p6a_execution.json",
    "p6b_six_chromosome_execution.json",
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
    if MIGRATION_PATH.is_file():
        migration = json.loads(MIGRATION_PATH.read_text())
        if migration.get("source_revision") != "six_chromosome_blocks_v1":
            raise RuntimeError("Existing core-coverage migration is inconsistent")
        print(json.dumps(migration, indent=2, sort_keys=True))
        return
    registry = json.loads(REGISTRY_PATH.read_text())
    if registry.get("revision_id") != "six_chromosome_blocks_v1":
        raise RuntimeError("Core-coverage archive requires the v1 split registry")

    metadata_sources = [METADATA_DIR / name for name in CURRENT_METADATA]
    audit_sources = [
        METADATA_DIR / "audits" / phase / name
        for phase in ("P4", "P5", "P6A")
        for name in ("review.json", "review.md")
    ]
    interval_sources = sorted(path for path in INTERVAL_DIR.rglob("*") if path.is_file())
    required = [*metadata_sources, *audit_sources]
    missing = [str(path) for path in required if not path.is_file()]
    if missing or not interval_sources:
        raise RuntimeError(f"Cannot archive incomplete v1 evidence: {missing}")

    execution = json.loads(EXECUTION_PATH.read_text())
    if (
        execution.get("split_revision") != "six_chromosome_blocks_v1"
        or execution.get("locked_test_block_signal_reads") != 0
    ):
        raise RuntimeError("Invalid v1 P6B execution record cannot be archived safely")
    completed = execution.get("jobs", [])
    observed_coverage = [
        json.loads((REPO_ROOT / job["validation_path"]).read_text()).get(
            "validation_core_coverage_fraction"
        )
        for job in completed
    ]
    if not observed_coverage or min(float(value) for value in observed_coverage) >= 0.999:
        raise RuntimeError("Expected incomplete v1 metric-core coverage evidence")

    for source in metadata_sources:
        copy_verified(source, ARCHIVE_DIR / source.name)
    for source in audit_sources:
        copy_verified(
            source,
            ARCHIVE_DIR / "audits" / source.relative_to(METADATA_DIR / "audits"),
        )
    for source in interval_sources:
        copy_verified(source, ARCHIVE_INTERVAL_DIR / source.relative_to(INTERVAL_DIR))

    run_evidence = []
    for job in completed:
        for key in ("job_path", "validation_path"):
            path = REPO_ROOT / job[key]
            if path.is_file():
                run_evidence.append(inventory(path))
        run_path = REPO_ROOT / job["output_dir"] / "run.json"
        if run_path.is_file():
            run_evidence.append(inventory(run_path))
        checkpoint = REPO_ROOT / job["checkpoint_path"]
        if checkpoint.is_file():
            run_evidence.append(
                {
                    "path": str(checkpoint.relative_to(REPO_ROOT)),
                    "bytes": checkpoint.stat().st_size,
                    "sha256": job["checkpoint_sha256"],
                }
            )

    migration = {
        "schema_version": 1,
        "source_revision": "six_chromosome_blocks_v1",
        "target_revision": "six_chromosome_blocks_v2",
        "prepared_at": utc_now(),
        "reason": (
            "Nearest-window eval cores omitted eligible bases; the first formal job "
            "measured validation_core_coverage_fraction=0.8435322914705329."
        ),
        "source_execution_status": execution.get("status"),
        "source_completed_jobs": len(completed),
        "source_failures": execution.get("failures", []),
        "source_locked_test_block_signal_reads": 0,
        "observed_validation_core_coverage_fractions": observed_coverage,
        "archived_metadata": [inventory(ARCHIVE_DIR / path.name) for path in metadata_sources],
        "archived_audits": [
            inventory(
                ARCHIVE_DIR / "audits" / path.relative_to(METADATA_DIR / "audits")
            )
            for path in audit_sources
        ],
        "archived_intervals": [
            inventory(ARCHIVE_INTERVAL_DIR / path.relative_to(INTERVAL_DIR))
            for path in interval_sources
        ],
        "retained_run_evidence": run_evidence,
        "retained_log_root": "logs/v2_p6b_six_chromosome_20260721",
        "retained_run_root": "runs/v2_p6b_six_chromosome_20260721",
    }
    atomic_json(MIGRATION_PATH, migration)

    execution.update(
        {
            "status": "controlled_invalidated",
            "invalidated_at": utc_now(),
            "invalidation_reason": migration["reason"],
            "replacement_revision": "six_chromosome_blocks_v2",
            "completed_at": utc_now(),
        }
    )
    atomic_json(EXECUTION_PATH, execution)
    print(json.dumps(migration, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
