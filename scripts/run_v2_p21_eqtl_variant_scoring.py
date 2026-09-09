#!/usr/bin/env python3
"""Run the frozen P17 Model-B checkpoints on the external I--V allele set."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
from datetime import datetime, timezone
from pathlib import Path
import subprocess
import sys
import time


REPO_ROOT = Path(__file__).resolve().parents[1]
VARIANT_MANIFEST = REPO_ROOT / "results/v2_p20_eqtl_variant_set/p20_eqtl_variant_set_manifest.json"
OUTPUT = REPO_ROOT / "results/v2_p21_eqtl_variant_scoring"
MEANS = REPO_ROOT / "results/v2_p16_iv_training_normalization/track_nonzero_means.tsv"
TRACK_MANIFEST = REPO_ROOT / "alphagenome_custom/metadata/v2/replicate_holdout_v1/training_track_manifest.tsv"
GROUP_MANIFEST = REPO_ROOT / "alphagenome_custom/metadata/v2/replicate_holdout_v1/training_group_manifest.tsv"
MODEL_SPECS = REPO_ROOT / "alphagenome_custom/metadata/v2/model_specs_v2.json"
MODEL_WEIGHTS = REPO_ROOT / "weights/alphagenome_pytorch/model_all_folds.safetensors"
FASTA = REPO_ROOT / "alphagenome_custom/reference/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa"
FAI = REPO_ROOT / "alphagenome_custom/reference/genome.fa.fai"
GTF = REPO_ROOT / "alphagenome_custom/reference/Caenorhabditis_elegans.WBcel235.115.gtf"
SEEDS = (20260714, 20260715, 20260716)


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(REPO_ROOT.resolve()))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--variant-manifest", type=Path, default=VARIANT_MANIFEST)
    return parser.parse_args()


def write_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def run_one(job: dict, spec_path: Path, log_path: Path) -> dict:
    physical_gpu = int(job["physical_gpu"])
    env = os.environ.copy()
    env["CUDA_VISIBLE_DEVICES"] = str(physical_gpu)
    command = [sys.executable, "-m", "scripts.score_v2_variants", "--spec", str(spec_path)]
    started = time.monotonic()
    with log_path.open("w", encoding="utf-8") as log:
        log.write("command\t" + " ".join(command) + "\n")
        log.write(f"physical_gpu\t{physical_gpu}\n")
        result = subprocess.run(command, cwd=REPO_ROOT, env=env, stdout=log, stderr=subprocess.STDOUT, check=False)
    audit_path = Path(job["output_dir"]) / "scoring_audit.json"
    audit = json.loads(audit_path.read_text()) if audit_path.is_file() else {}
    return {
        **job,
        "status": "completed" if result.returncode == 0 and audit.get("status") == "completed" else "failed",
        "returncode": result.returncode,
        "elapsed_seconds": time.monotonic() - started,
        "log_path": relative(log_path),
        "audit_path": relative(audit_path) if audit_path.is_file() else None,
        "audit_sha256": sha256(audit_path) if audit_path.is_file() else None,
        "variant_count": audit.get("variant_count"),
        "locked_test_block_signal_reads": audit.get("locked_test_block_signal_reads"),
        "failure": audit.get("failures", []),
    }


def main() -> None:
    args = parse_args()
    variant_manifest_path = args.variant_manifest.resolve()
    execution_path = REPO_ROOT / "alphagenome_custom/metadata/v2/p21_eqtl_variant_scoring_execution.json"
    if execution_path.is_file():
        execution = json.loads(execution_path.read_text())
        if execution.get("status") == "completed" and execution.get("registered_tasks_completed") == 15:
            print(json.dumps(execution, indent=2, sort_keys=True))
            return
        raise FileExistsError(f"P21 output exists but is not a reusable completed execution: {execution_path}")
    manifest = json.loads(variant_manifest_path.read_text())
    if manifest.get("status") != "PASS":
        raise RuntimeError("P20 variant-set coordinate audit is not PASS")
    vcf = REPO_ROOT / manifest["selected_vcf"]
    required = (vcf, MEANS, TRACK_MANIFEST, GROUP_MANIFEST, MODEL_SPECS, MODEL_WEIGHTS, FASTA, FAI, GTF)
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        raise FileNotFoundError("Missing P21 input: " + ", ".join(missing))
    OUTPUT.mkdir(parents=True, exist_ok=True)
    spec_dir = OUTPUT / "specs"
    auth_dir = OUTPUT / "authorizations"
    log_dir = OUTPUT / "logs"
    for path in (spec_dir, auth_dir, log_dir):
        path.mkdir(parents=True, exist_ok=True)
    inputs = {
        "checkpoint": None,
        "checkpoint_run": None,
        "fasta": relative(FASTA),
        "fai": relative(FAI),
        "gtf": relative(GTF),
        "group_manifest": relative(GROUP_MANIFEST),
        "means": relative(MEANS),
        "track_manifest": relative(TRACK_MANIFEST),
        "model_specs": relative(MODEL_SPECS),
        "model_weights": relative(MODEL_WEIGHTS),
        "vcf": relative(vcf),
    }
    records = []
    job_index = 0
    for fold in range(1, 6):
        for seed in SEEDS:
            job_index += 1
            checkpoint = REPO_ROOT / f"runs/v2_p17_iv_controlled_matrix/B_iv_dual/seed_{seed}/fold_{fold}/checkpoint.pt"
            checkpoint_run = checkpoint.parent / "run.json"
            if not checkpoint.is_file() or not checkpoint_run.is_file():
                raise FileNotFoundError(f"Missing P17 checkpoint pair for fold {fold}, seed {seed}")
            job_id = f"p21:eqtl_variant:B_iv_dual:seed{seed}:fold{fold}"
            output_dir = OUTPUT / "evaluations" / f"seed_{seed}" / f"fold_{fold}"
            output_dir.mkdir(parents=True, exist_ok=True)
            job_inputs = {**inputs, "checkpoint": relative(checkpoint), "checkpoint_run": relative(checkpoint_run)}
            input_hashes = {key: sha256(REPO_ROOT / value) for key, value in job_inputs.items()}
            spec = {
                "schema_version": 1,
                "contract": "v2_variant_scoring_v1",
                "execution_id": job_id,
                "phase": "P21",
                "device": "cuda",
                "physical_gpu": 2 if job_index % 2 else 3,
                "final_test_access": "prohibited",
                "locked_inputs": job_inputs,
                "input_sha256": input_hashes,
                "model": {
                    "model_id": "B",
                    "loss": "paper",
                    "fold": fold,
                    "seed": seed,
                    "sequence_length": 131072,
                    "hidden_channels": 64,
                    "mean_column": f"fold_{fold}_train_nonzero_mean",
                },
                "inference": {
                    "forward_reverse_complement_ensemble": True,
                    "allele_window_policy": "variant_start_anchored_right_context_crop_or_extend",
                    "track_semantics": "all_241_tracks_unstranded_no_rc_permutation",
                    "maximum_alleles": 1000,
                    "require_pass": False,
                },
                "ism": {"enabled": False},
                "output_dir": str(output_dir),
                "cuda_authorization_path": relative(auth_dir / f"seed_{seed}_fold_{fold}.json"),
            }
            spec_path = spec_dir / f"seed_{seed}_fold_{fold}.json"
            write_json(spec_path, spec)
            authorization = {
                "schema_version": 1,
                "approved": True,
                "scope": f"variant_scoring:{job_id}",
                "execution_id": job_id,
                "spec_sha256": sha256(spec_path),
                "physical_gpu": spec["physical_gpu"],
                "checkpoint_sha256": input_hashes["checkpoint"],
                "vcf_sha256": input_hashes["vcf"],
            }
            write_json(auth_dir / f"seed_{seed}_fold_{fold}.json", authorization)
            job = {
                "job_id": job_id,
                "fold": fold,
                "seed": seed,
                "physical_gpu": spec["physical_gpu"],
                "spec_path": relative(spec_path),
                "output_dir": str(output_dir),
                "checkpoint": relative(checkpoint),
                "vcf": relative(vcf),
            }
            records.append(run_one(job, spec_path, log_dir / f"seed_{seed}_fold_{fold}.log"))
            if records[-1]["status"] != "completed":
                raise RuntimeError(f"P21 job failed: {job_id}")

    execution = {
        "schema_version": 1,
        "phase": "P21",
        "contract": "v2_variant_scoring_v1",
        "started_at": records[0].get("started_at", utc_now()) if records else utc_now(),
        "completed_at": utc_now(),
        "status": "completed",
        "registered_tasks": len(records),
        "registered_tasks_completed": sum(row["status"] == "completed" for row in records),
        "jobs": records,
        "variant_set_manifest": relative(variant_manifest_path),
        "variant_set_manifest_sha256": sha256(variant_manifest_path),
        "final_test_access": "prohibited",
        "locked_test_block_signal_reads": 0,
        "model_inference": "external_eqtl_candidate_variant_scoring",
        "claim_boundary": [
            "Independent public wild-strain cohort and published candidate loci only.",
            "Not causal validation, not unseen-condition decoding and not a population-level laboratory estimate.",
        ],
    }
    write_json(execution_path, execution)
    print(json.dumps(execution, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
