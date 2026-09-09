#!/usr/bin/env python3
"""Preregister the strict I--V controlled ablation and baseline matrix.

P17 is deliberately limited to three retrained configurations under one data
contract: full dual-resolution Model B, a 128-bp-only Model B, and a
Basenji2-style from-scratch baseline. It is not an external benchmark.
"""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
METADATA_ROOT = REPO_ROOT / "alphagenome_custom" / "metadata" / "v2"
SPEC_PATH = METADATA_ROOT / "p17_iv_controlled_matrix_spec.json"
P15_ROOT = REPO_ROOT / "results/v2_p15_iv_training_interval_contract"
P16_ROOT = REPO_ROOT / "results/v2_p16_iv_training_normalization"
ALLOWED_CHROMOSOMES = ["I", "II", "III", "IV", "V"]
SEEDS = [20260714, 20260715, 20260716]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def artifact(path: Path) -> dict[str, str]:
    if not path.is_file():
        raise FileNotFoundError(path)
    return {"path": str(path.relative_to(REPO_ROOT)), "sha256": sha256(path)}


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def main() -> None:
    if SPEC_PATH.exists():
        raise FileExistsError(f"P17 preregistration already exists: {SPEC_PATH}")
    p16_execution = METADATA_ROOT / "p16_iv_training_normalization_execution.json"
    p16_review = METADATA_ROOT / "audits/P16/review.json"
    p16_summary = P16_ROOT / "track_nonzero_means_summary.json"
    p16_means = P16_ROOT / "track_nonzero_means.tsv"
    if json.loads(p16_execution.read_text()).get("status") != "completed":
        raise RuntimeError("P17 requires completed P16 normalization")
    if json.loads(p16_review.read_text()).get("status") != "PASS":
        raise RuntimeError("P17 requires the R16 pass record")
    summary = json.loads(p16_summary.read_text())
    if (
        summary.get("contract") != "strict_interval_root"
        or summary.get("allowed_chromosomes") != ALLOWED_CHROMOSOMES
        or summary.get("locked_test_block_signal_reads") != 0
    ):
        raise RuntimeError("P16 normalization does not satisfy the strict I--V contract")

    intervals = []
    for fold in range(1, 6):
        for role in ("train", "valid"):
            path = P15_ROOT / "intervals" / f"fold_{fold}" / f"{role}.tsv"
            intervals.append({"fold": fold, "role": role, **artifact(path)})

    source = {
        "p15_contract": artifact(P15_ROOT / "p15_iv_training_interval_contract.json"),
        "p15_execution": artifact(METADATA_ROOT / "p15_iv_training_contract_execution.json"),
        "p15_review": artifact(METADATA_ROOT / "audits/P15/review.json"),
        "p16_spec": artifact(METADATA_ROOT / "p16_iv_training_normalization_spec.json"),
        "p16_execution": artifact(p16_execution),
        "p16_review": artifact(p16_review),
        "p16_means": artifact(p16_means),
        "p16_summary": artifact(p16_summary),
        "training_track_manifest": artifact(METADATA_ROOT / "replicate_holdout_v1/training_track_manifest.tsv"),
        "training_group_manifest": artifact(METADATA_ROOT / "replicate_holdout_v1/training_group_manifest.tsv"),
        "heldout_track_manifest": artifact(METADATA_ROOT / "replicate_holdout_v1/heldout_track_manifest.tsv"),
        "heldout_group_manifest": artifact(METADATA_ROOT / "replicate_holdout_v1/heldout_group_manifest.tsv"),
        "heldout_track_indices": artifact(METADATA_ROOT / "replicate_holdout_v1/heldout_track_indices.tsv"),
        "intervals": intervals,
    }
    spec = {
        "schema_version": 1,
        "phase": "P17",
        "status": "preregistered",
        "created_at": utc_now(),
        "purpose": "Minimal fair I--V-only controlled matrix for dual-resolution utility and a public Basenji2-style architecture baseline.",
        "allowed_chromosomes": ALLOWED_CHROMOSOMES,
        "claim_boundary": [
            "This is within-collection biological-unit holdout evaluation, not an external, study-isolated, laboratory-isolated or unseen-condition benchmark.",
            "No chromosome-X or locked-test signal may be read.",
            "The Basenji2-style configuration is a task-matched from-scratch public-architecture baseline, not a native Basenji human/mouse checkpoint comparison.",
            "The 128-bp-only configuration derives 1-bp values by uniform expansion solely to retain a common metric interface; it does not demonstrate base-resolution utility.",
        ],
        "source": source,
        "training": {
            "means_path": "results/v2_p16_iv_training_normalization/track_nonzero_means.tsv",
            "sequence_length": 131072,
            "max_steps": 2000,
            "learning_rate": 0.0001,
            "weight_decay": 0.0001,
            "max_shift_bp": 1024,
            "reverse_complement_probability": 0.5,
            "hidden_channels": 64,
        },
        "evaluation": {
            "track_manifest": "alphagenome_custom/metadata/v2/replicate_holdout_v1/heldout_track_manifest.tsv",
            "group_manifest": "alphagenome_custom/metadata/v2/replicate_holdout_v1/heldout_group_manifest.tsv",
            "prediction_track_indices": "alphagenome_custom/metadata/v2/replicate_holdout_v1/heldout_track_indices.tsv",
            "metric": "frozen primary score = 0.5 * gene-exon Pearson(log1p) + 0.5 * 128-bp Pearson(log1p)",
        },
        "matrix": {
            "folds": [1, 2, 3, 4, 5],
            "seeds": SEEDS,
            "configurations": [
                {
                    "id": "B_iv_dual",
                    "model": "B",
                    "loss": "paper",
                    "gene_weight": 0.1,
                    "loss_1bp_weight": 1.0,
                    "loss_128bp_weight": 1.0,
                    "interpretation": "Strict I--V retrained full Model B control.",
                },
                {
                    "id": "B_128bp_only",
                    "model": "B_128bp_head_only",
                    "loss": "paper",
                    "gene_weight": 0.0,
                    "loss_1bp_weight": 0.0,
                    "loss_128bp_weight": 1.0,
                    "interpretation": "128-bp-only Model B; derived 1-bp output is metric compatibility only.",
                },
                {
                    "id": "Basenji2_style",
                    "model": "Basenji2_style",
                    "loss": "paper",
                    "gene_weight": 0.0,
                    "loss_1bp_weight": 0.0,
                    "loss_128bp_weight": 1.0,
                    "implementation_reference": "Calico Basenji repository commit 06ce5d387e20b47184d05433b3983163c5f923cd; architectural reference only.",
                    "interpretation": "Task-matched from-scratch public-architecture baseline.",
                },
            ],
            "total_jobs": 45,
        },
        "gpu_policy": {
            "allowed_physical_indices": [2, 3],
            "maximum_concurrent_jobs_per_gpu": 1,
            "minimum_free_mib": 70000,
        },
        "final_test_access": "prohibited",
    }
    SPEC_PATH.write_text(json.dumps(spec, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"status": "preregistered", "spec": str(SPEC_PATH.relative_to(REPO_ROOT))}, indent=2))


if __name__ == "__main__":
    main()
