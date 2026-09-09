#!/usr/bin/env python3
"""Preregister external represented-head scoring from the completed P18 source set.

This preparation step creates one-track manifests for each checksum-verified
external RNA measurement and binds them to the original 241-head index.  It
does not load model weights or read signal values.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
METADATA_ROOT = REPO_ROOT / "alphagenome_custom/metadata/v2"
SPEC_PATH = METADATA_ROOT / "p19_external_represented_head_scoring_spec.json"
P17_EXECUTION = METADATA_ROOT / "p17_iv_controlled_matrix_execution.json"
P18_EXECUTION = METADATA_ROOT / "p18_external_continuation_execution.json"
P18_SUMMARY = REPO_ROOT / "alphagenome_custom/tracks/rna_seq_v2_external_p18_continuation/p18_reprocessing_summary.json"
TRAIN_TRACKS = REPO_ROOT / "alphagenome_custom/metadata/v2/replicate_holdout_v1/training_track_manifest.tsv"
TRAIN_GROUPS = REPO_ROOT / "alphagenome_custom/metadata/v2/replicate_holdout_v1/training_group_manifest.tsv"
MEANS = REPO_ROOT / "results/v2_p16_iv_training_normalization/track_nonzero_means.tsv"
INTERVAL_ROOT = REPO_ROOT / "results/v2_p15_iv_training_interval_contract/intervals"
EXTERNAL_ROOT = REPO_ROOT / "alphagenome_custom/tracks/rna_seq_v2_external_p18_continuation"
RESULT_ROOT = REPO_ROOT / "results/v2_p19_external_scoring"

ALLOWED_CHROMOSOMES = ["I", "II", "III", "IV", "V"]
EXPECTED_RUNS = ["SRR18463404", "SRR18463405", "SRR18463406", "SRR3560831"]
SEEDS = [20260714, 20260715, 20260716]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, str]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def repo_relative(path: Path) -> str:
    return str(path.resolve().relative_to(REPO_ROOT.resolve()))


def validate_i_to_v(interval_path: Path) -> None:
    rows = read_tsv(interval_path)
    if not rows or any(row.get("chromosome") not in ALLOWED_CHROMOSOMES for row in rows):
        raise RuntimeError(f"P19 interval contract is not I--V-only: {interval_path}")
    if any(row.get("role") == "test_locked" for row in rows):
        raise RuntimeError(f"P19 interval contract contains locked rows: {interval_path}")


def main() -> None:
    if SPEC_PATH.exists():
        raise FileExistsError(f"P19 specification is immutable: {SPEC_PATH}")
    p17 = json.loads(P17_EXECUTION.read_text())
    p18 = json.loads(P18_EXECUTION.read_text())
    p18_summary = json.loads(P18_SUMMARY.read_text())
    if p18.get("status") != "completed" or p18.get("locked_test_block_signal_reads") != 0:
        raise RuntimeError("P19 requires completed P18 with zero locked-test reads")
    if p18_summary.get("completed_accessions") != EXPECTED_RUNS:
        raise RuntimeError("P18 source set does not match the registered four-run set")

    tracks = read_tsv(TRAIN_TRACKS)
    groups = read_tsv(TRAIN_GROUPS)
    track_index = next((i for i, row in enumerate(tracks) if row["group_id"] == "RNA_V2_G0054"), None)
    if track_index is None:
        raise RuntimeError("Represented head RNA_V2_G0054 is absent from the 241-head manifest")
    group_row = next(row for row in groups if row["group_id"] == "RNA_V2_G0054")
    checkpoint_jobs = [
        row for row in p17.get("jobs", []) if row.get("configuration") == "B_iv_dual"
    ]
    if len(checkpoint_jobs) != 15:
        raise RuntimeError(f"Expected 15 P17 B_iv_dual checkpoints, observed {len(checkpoint_jobs)}")
    if not MEANS.is_file() or len(read_tsv(MEANS)) != 241:
        raise RuntimeError("P19 requires the 241-row P16 I--V means table")

    p18_transport = REPO_ROOT / p18["transport_manifest"]
    transport = {row["run_accession"]: row for row in read_tsv(p18_transport)}
    jobs: list[dict[str, object]] = []
    source_records: list[dict[str, object]] = []
    manifest_fields = list(tracks[0])
    group_fields = list(groups[0])
    index_fields = ["track_index"]
    for run in EXPECTED_RUNS:
        bw = EXTERNAL_ROOT / f"{run}.bw"
        audit = EXTERNAL_ROOT / "sample_audits" / f"{run}.json"
        if not bw.is_file() or not audit.is_file():
            raise FileNotFoundError(f"P18 source artifact missing for {run}")
        audit_payload = json.loads(audit.read_text())
        if audit_payload.get("output_sha256") != sha256(bw):
            raise RuntimeError(f"P18 output checksum mismatch for {run}")
        if audit_payload.get("coverage_chromosomes") != ALLOWED_CHROMOSOMES:
            raise RuntimeError(f"P18 output is not I--V-only for {run}")
        source_track_dir = RESULT_ROOT / "manifests" / run
        source_track = source_track_dir / "track_manifest.tsv"
        source_group = source_track_dir / "group_manifest.tsv"
        source_indices = source_track_dir / "prediction_track_indices.tsv"
        track_row = dict(tracks[track_index])
        track_row.update({"output_path": repo_relative(bw), "output_size_bytes": str(bw.stat().st_size), "output_sha256": sha256(bw), "run_accessions": run, "n_runs": "1", "n_biological_units": "1"})
        external_group = dict(group_row)
        external_group.update({"sample_ids": run, "n_members": "1", "n_runs": "1", "n_biological_units": "1", "source_n_biological_units": "1", "training_n_biological_units": "0", "holdout_unit_id": run, "replicate_holdout_role": "external_represented_head_measurement"})
        write_tsv(source_track, [track_row], manifest_fields)
        write_tsv(source_group, [external_group], group_fields)
        write_tsv(source_indices, [{"track_index": str(track_index)}], index_fields)
        source_records.append({
            "run_accession": run,
            "study_accession": transport[run]["geo_accession"],
            "sra_study_accession": transport[run]["sra_study_accession"],
            "bioproject_accession": transport[run]["bioproject_accession"],
            "represented_group_id": "RNA_V2_G0054",
            "track_index": track_index,
            "coverage_path": repo_relative(bw),
            "coverage_sha256": sha256(bw),
            "audit_path": repo_relative(audit),
            "track_manifest": repo_relative(source_track),
            "group_manifest": repo_relative(source_group),
            "prediction_track_indices": repo_relative(source_indices),
        })
        for checkpoint in checkpoint_jobs:
            fold = int(checkpoint["fold"])
            interval = INTERVAL_ROOT / f"fold_{fold}" / "valid.tsv"
            validate_i_to_v(interval)
            jobs.append({
                "job_id": f"p19:{run}:seed{checkpoint['seed']}:fold{fold}",
                "run_accession": run,
                "study_accession": transport[run]["geo_accession"],
                "source_p17_job_id": checkpoint["job_id"],
                "fold": fold,
                "seed": int(checkpoint["seed"]),
                "model": "B",
                "training_loss": "paper",
                "hidden_channels": 64,
                "mean_column": f"fold_{fold}_train_nonzero_mean",
                "checkpoint_path": checkpoint["checkpoint_path"],
                "checkpoint_sha256": checkpoint["checkpoint_sha256"],
                "interval_path": repo_relative(interval),
                "interval_sha256": sha256(interval),
                "track_manifest": repo_relative(source_track),
                "group_manifest": repo_relative(source_group),
                "prediction_track_indices": repo_relative(source_indices),
            })

    spec = {
        "schema_version": 1,
        "phase": "P19",
        "status": "preregistered",
        "created_at": utc_now(),
        "purpose": "External I--V-only scoring of four checksum-verified RNA measurements against the represented RNA_V2_G0054 output head.",
        "represented_head": {"group_id": "RNA_V2_G0054", "track_index": track_index, "context": "N2 whole-animal young-adult broad context"},
        "allowed_chromosomes": ALLOWED_CHROMOSOMES,
        "final_test_access": "prohibited",
        "model_inference": "inference_only; frozen P17 B_iv_dual checkpoints; no training",
        "claim_boundary": [
            "This is an external measurement benchmark for a represented target identity, not an unseen-target or unseen-condition decoder.",
            "The four sources span two study accessions and are not a population-level laboratory generalization estimate.",
            "No chromosome-X or locked-test signal may be read.",
            "P19 performance is not claimable until R19 verifies all source, checkpoint and interval hashes.",
        ],
        "gpu_policy": {"allowed_physical_indices": [2, 3], "maximum_concurrent_jobs_per_gpu": 1, "minimum_free_mib": 70000},
        "source": {
            "p17_execution": {"path": repo_relative(P17_EXECUTION), "sha256": sha256(P17_EXECUTION)},
            "p18_execution": {"path": repo_relative(P18_EXECUTION), "sha256": sha256(P18_EXECUTION)},
            "p18_summary": {"path": repo_relative(P18_SUMMARY), "sha256": sha256(P18_SUMMARY)},
            "p18_transport": {"path": repo_relative(p18_transport), "sha256": sha256(p18_transport)},
            "training_track_manifest": {"path": repo_relative(TRAIN_TRACKS), "sha256": sha256(TRAIN_TRACKS)},
            "training_group_manifest": {"path": repo_relative(TRAIN_GROUPS), "sha256": sha256(TRAIN_GROUPS)},
            "means": {"path": repo_relative(MEANS), "sha256": sha256(MEANS)},
            "source_records": source_records,
        },
        "evaluation": {"sequence_length": 131072, "metric": "0.5 * gene-exon Pearson(log1p) + 0.5 * 128-bp Pearson(log1p)", "jobs": len(jobs)},
        "jobs": jobs,
    }
    SPEC_PATH.write_text(json.dumps(spec, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"status": "prepared", "phase": "P19", "jobs": len(jobs), "sources": EXPECTED_RUNS, "track_index": track_index, "spec": repo_relative(SPEC_PATH)}, indent=2))


if __name__ == "__main__":
    main()
