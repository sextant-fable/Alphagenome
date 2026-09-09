#!/usr/bin/env python3
"""Freeze the P14 I--V model-to-held-out comparison before GPU inference.

P14 reuses only the completed P12 ``B_paper`` checkpoints.  It evaluates them
on P13's immutable I--V-only interval contract so model-to-held-out scores and
training-aggregate-to-held-out agreement share a target, interval and metric
definition.  This preparation step reads metadata and checkpoint hashes only.
"""

from __future__ import annotations

import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[1]
METADATA_ROOT = REPO_ROOT / "alphagenome_custom" / "metadata" / "v2"
P12_ROOT = METADATA_ROOT / "replicate_holdout_v1"
P13_CONTRACT = REPO_ROOT / "results/v2_p13_iv_interval_contract/p13_iv_interval_contract.json"
P13_EXECUTION = METADATA_ROOT / "p13_submission_evidence_execution.json"
P13_REVIEW = METADATA_ROOT / "audits/P13/review.json"
SPEC_PATH = METADATA_ROOT / "p14_iv_frozen_model_reference_spec.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def relative(path: Path) -> str:
    return str(path.relative_to(REPO_ROOT))


def main() -> None:
    if SPEC_PATH.exists():
        raise FileExistsError(f"P14 contract already exists and is immutable: {SPEC_PATH}")
    required = (P13_CONTRACT, P13_EXECUTION, P13_REVIEW, P12_ROOT / "p12_replicate_holdout_job_records.tsv")
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        raise FileNotFoundError(f"P14 requires completed P12/P13 artifacts: {missing}")
    p13_contract = json.loads(P13_CONTRACT.read_text())
    p13_execution = json.loads(P13_EXECUTION.read_text())
    p13_review = json.loads(P13_REVIEW.read_text())
    if p13_contract.get("allowed_chromosomes") != ["I", "II", "III", "IV", "V"]:
        raise RuntimeError("P14 requires P13's I--V-only interval contract")
    if p13_execution.get("status") != "completed" or p13_review.get("status") != "PASS":
        raise RuntimeError("P14 requires completed, passing P13 execution and review")

    job_records = read_tsv(P12_ROOT / "p12_replicate_holdout_job_records.tsv")
    selected = [
        row
        for row in job_records
        if row["configuration"] == "B_paper" and row["model"] == "B" and row["loss"] == "paper"
    ]
    selected.sort(key=lambda row: (int(row["fold"]), int(row["seed"])))
    if len(selected) != 15 or {(int(row["fold"]), int(row["seed"])) for row in selected} != {
        (fold, seed) for fold in range(1, 6) for seed in (20260714, 20260715, 20260716)
    }:
        raise RuntimeError("P14 requires exactly the 15 completed P12 B_paper fold/seed checkpoints")

    jobs: list[dict[str, Any]] = []
    for row in selected:
        checkpoint = REPO_ROOT / row["checkpoint_path"]
        if not checkpoint.is_file() or sha256(checkpoint) != row["checkpoint_sha256"]:
            raise RuntimeError(f"P12 checkpoint hash mismatch for {row['job_id']}")
        interval_path = P13_CONTRACT.parent / "intervals" / f"fold_{row['fold']}" / "valid.tsv"
        if not interval_path.is_file():
            raise FileNotFoundError(interval_path)
        jobs.append(
            {
                "job_id": f"p14:iv_model_reference:B_paper:seed{row['seed']}:fold{row['fold']}",
                "source_p12_job_id": row["job_id"],
                "model": "B",
                "configuration": "B_paper",
                "training_loss": "paper",
                "fold": int(row["fold"]),
                "seed": int(row["seed"]),
                "hidden_channels": int(row["hidden_channels"]),
                "mean_column": row["mean_column"],
                "checkpoint_path": row["checkpoint_path"],
                "checkpoint_sha256": row["checkpoint_sha256"],
                "interval_path": relative(interval_path),
                "interval_sha256": sha256(interval_path),
            }
        )

    p12_spec = METADATA_ROOT / "p12_replicate_holdout_spec.json"
    payload = {
        "schema_version": 1,
        "phase": "P14",
        "status": "preregistered",
        "created_at": utc_now(),
        "purpose": "Frozen P12 Model B inference on P13 I--V-only held-out biological-unit cores for a metric-compatible replicate-agreement reference.",
        "claim_boundary": [
            "Within-collection biological-unit-held-out comparison only.",
            "No model training, hyperparameter selection, chromosome-X read or locked-test access.",
            "Not an external study/laboratory benchmark, unseen-condition result or regulatory-variant result.",
        ],
        "gpu_policy": {
            "allowed_physical_indices": [2, 3],
            "minimum_free_mib": 70000,
            "maximum_concurrent_jobs_per_gpu": 1,
        },
        "source_contract": {
            "p12_job_records": relative(P12_ROOT / "p12_replicate_holdout_job_records.tsv"),
            "p12_job_records_sha256": sha256(P12_ROOT / "p12_replicate_holdout_job_records.tsv"),
            "p12_spec": relative(p12_spec),
            "p12_spec_sha256": sha256(p12_spec),
            "p13_interval_contract": relative(P13_CONTRACT),
            "p13_interval_contract_sha256": sha256(P13_CONTRACT),
            "p13_execution": relative(P13_EXECUTION),
            "p13_execution_sha256": sha256(P13_EXECUTION),
            "p13_review": relative(P13_REVIEW),
            "p13_review_sha256": sha256(P13_REVIEW),
        },
        "evaluation": {
            "track_manifest": "alphagenome_custom/metadata/v2/replicate_holdout_v1/heldout_track_manifest.tsv",
            "group_manifest": "alphagenome_custom/metadata/v2/replicate_holdout_v1/heldout_group_manifest.tsv",
            "prediction_track_indices": "alphagenome_custom/metadata/v2/replicate_holdout_v1/heldout_track_indices.tsv",
            "means_path": "alphagenome_custom/metadata/v2/replicate_holdout_v1/track_nonzero_means.tsv",
            "primary_metric": "0.5 * gene-exon log1p Pearson + 0.5 * 128-bp log1p Pearson",
            "allowed_chromosomes": ["I", "II", "III", "IV", "V"],
            "final_test_access": "prohibited",
        },
        "jobs": jobs,
    }
    SPEC_PATH.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"status": "preregistered", "spec": relative(SPEC_PATH), "jobs": len(jobs)}, indent=2))


if __name__ == "__main__":
    main()
