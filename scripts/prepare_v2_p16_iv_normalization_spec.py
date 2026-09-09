#!/usr/bin/env python3
"""Preregister P16's strict I--V-only fold-normalization calculation."""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
METADATA_ROOT = REPO_ROOT / "alphagenome_custom" / "metadata" / "v2"
P15_ROOT = REPO_ROOT / "results" / "v2_p15_iv_training_interval_contract"
SPEC_PATH = METADATA_ROOT / "p16_iv_training_normalization_spec.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(REPO_ROOT.resolve()))


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    if SPEC_PATH.exists() and not args.force:
        raise FileExistsError(f"P16 specification already exists: {SPEC_PATH}")

    p15_execution = METADATA_ROOT / "p15_iv_training_contract_execution.json"
    p15_review = METADATA_ROOT / "audits" / "P15" / "review.json"
    p15_contract = P15_ROOT / "p15_iv_training_interval_contract.json"
    required = (p15_execution, p15_review, p15_contract)
    if not all(path.is_file() for path in required):
        raise FileNotFoundError("P16 requires completed P15 execution, review and contract artifacts")
    if json.loads(p15_execution.read_text()).get("status") != "completed":
        raise RuntimeError("P15 execution is not complete")
    if json.loads(p15_review.read_text()).get("status") != "PASS":
        raise RuntimeError("P15 review did not pass")
    p15_payload = json.loads(p15_contract.read_text())
    if p15_payload.get("allowed_chromosomes") != ["I", "II", "III", "IV", "V"]:
        raise RuntimeError("P15 contract does not enforce I--V-only intervals")

    intervals = []
    for fold in range(1, 6):
        path = P15_ROOT / "intervals" / f"fold_{fold}" / "train.tsv"
        if not path.is_file():
            raise FileNotFoundError(path)
        intervals.append({"fold": fold, "path": relative(path), "sha256": sha256(path)})

    result_root = REPO_ROOT / "results" / "v2_p16_iv_training_normalization"
    payload = {
        "schema_version": 1,
        "phase": "P16",
        "status": "preregistered",
        "created_at": utc_now(),
        "purpose": "Recompute 241 fold-specific normalization means from P11 training aggregates on P15's I--V-only train blocks before new ablation or public-baseline training.",
        "claim_boundary": [
            "No chromosome-X or locked-test signal may be read.",
            "No model, checkpoint, external source or GPU is used.",
            "The output is a normalization contract, not an external benchmark or a model result.",
        ],
        "allowed_chromosomes": ["I", "II", "III", "IV", "V"],
        "source": {
            "p15_execution": relative(p15_execution),
            "p15_execution_sha256": sha256(p15_execution),
            "p15_review": relative(p15_review),
            "p15_review_sha256": sha256(p15_review),
            "p15_contract": relative(p15_contract),
            "p15_contract_sha256": sha256(p15_contract),
            "training_track_manifest": "alphagenome_custom/metadata/v2/replicate_holdout_v1/training_track_manifest.tsv",
            "training_track_manifest_sha256": sha256(METADATA_ROOT / "replicate_holdout_v1" / "training_track_manifest.tsv"),
            "intervals": intervals,
        },
        "outputs": {
            "means": relative(result_root / "track_nonzero_means.tsv"),
            "summary": relative(result_root / "track_nonzero_means_summary.json"),
        },
        "implementation": {
            "mean_calculator": "scripts/compute_v2_replicate_holdout_means.py",
            "mean_calculator_sha256": sha256(REPO_ROOT / "scripts/compute_v2_replicate_holdout_means.py"),
        },
    }
    SPEC_PATH.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"status": "preregistered", "spec": relative(SPEC_PATH)}, indent=2))


if __name__ == "__main__":
    main()
