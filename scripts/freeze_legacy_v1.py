#!/usr/bin/env python3
"""Create a reproducible, metadata-only freeze of the legacy RNA-seq11 work."""

from __future__ import annotations

import csv
from datetime import datetime, timezone
import hashlib
from importlib import metadata as importlib_metadata
import json
from pathlib import Path
import platform
import subprocess
import sys
from typing import Any, Iterable


REPO_ROOT = Path(__file__).resolve().parents[1]
OUTPUT_DIR = REPO_ROOT / "alphagenome_custom/metadata/legacy_v1"
AUDIT_DIR = REPO_ROOT / "alphagenome_custom/metadata/v2/audits/P0"

RAW_DIR = REPO_ROOT / "alphagenome_custom/tracks/bw_2026.4.20.training_input_bigwig"
GROUPED_DIR = REPO_ROOT / "alphagenome_custom/tracks/rna_seq_grouped"
WEIGHTS_PATH = REPO_ROOT / "weights/alphagenome_pytorch/model_all_folds.safetensors"
CANDIDATE_SOURCE = (
    REPO_ROOT
    / "runs/rna_seq11_all_representation_inventory_20260530/top15_models.tsv"
)
UTILITY_SOURCE = (
    REPO_ROOT
    / "runs/rna_seq11_representation_utility_summary_20260530/"
    "all_candidate_utility_summary.tsv"
)


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def sha256(path: Path, chunk_size: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def rel(path: Path) -> str:
    try:
        return str(path.relative_to(REPO_ROOT))
    except ValueError:
        return str(path)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, Any]], fields: Iterable[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(fields), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def git_output(*args: str) -> str:
    return subprocess.check_output(
        ["git", *args], cwd=REPO_ROOT, text=True
    ).strip()


def file_record(path: Path) -> dict[str, Any]:
    exists = path.is_file()
    return {
        "path": rel(path),
        "exists": str(exists),
        "size_bytes": path.stat().st_size if exists else "",
        "sha256": sha256(path) if exists else "",
    }


def raw_inventory() -> list[dict[str, Any]]:
    metadata_rows = {
        row["sample_id"]: row
        for row in read_tsv(REPO_ROOT / "alphagenome_custom/metadata/track_metadata.tsv")
    }
    rows = []
    for path in sorted(RAW_DIR.glob("*.bw")):
        sample_id = path.stem
        source = metadata_rows.get(sample_id, {})
        declared = REPO_ROOT / "alphagenome_custom" / source.get("file_path", "")
        record = file_record(path)
        record.update(
            {
                "sample_id": sample_id,
                "declared_legacy_path": rel(declared),
                "declared_legacy_path_exists": str(declared.is_file()),
                "normalization_status": source.get("normalization_status", ""),
                "data_source": source.get("data_source", ""),
                "stage": source.get("stage", ""),
                "tissue": source.get("biosample_name", ""),
            }
        )
        rows.append(record)
    return rows


def grouped_inventory() -> list[dict[str, Any]]:
    groups = read_tsv(REPO_ROOT / "alphagenome_custom/metadata/track_groups.tsv")
    metadata_rows = {
        row["group_id"]: row
        for row in read_tsv(
            REPO_ROOT / "alphagenome_custom/metadata/track_metadata_grouped.tsv"
        )
    }
    rows = []
    for group in groups:
        path = REPO_ROOT / "alphagenome_custom" / group["aggregated_file_path"]
        metadata_row = metadata_rows.get(group["group_id"], {})
        record = file_record(path)
        record.update(
            {
                "group_id": group["group_id"],
                "sample_ids": group["sample_ids"],
                "n_replicates": group["n_replicates"],
                "input_files_declared": group["input_files"],
                "normalization_plan": group["normalization_plan"],
                "normalization_status": metadata_row.get("normalization_status", ""),
                "ontology_curie": metadata_row.get("ontology_curie") or "not_available",
            }
        )
        rows.append(record)
    return rows


def reference_inventory() -> list[dict[str, Any]]:
    reference_dir = REPO_ROOT / "alphagenome_custom/reference"
    names = [
        "Caenorhabditis_elegans.WBcel235.dna.toplevel.fa",
        "Caenorhabditis_elegans.WBcel235.dna.toplevel.fa.gz",
        "genome.fa.fai",
        "Caenorhabditis_elegans.WBcel235.115.gtf",
        "Caenorhabditis_elegans.WBcel235.115.gtf.gz",
    ]
    rows = []
    for name in names:
        record = file_record(reference_dir / name)
        record["asset_type"] = (
            "FASTA_INDEX"
            if name.endswith(".fai")
            else "GTF"
            if ".gtf" in name
            else "FASTA"
        )
        rows.append(record)
    return rows


def split_inventory() -> list[dict[str, Any]]:
    intervals_dir = REPO_ROOT / "alphagenome_custom/intervals"
    rows = []
    for split, filename in (
        ("train", "train.bed"),
        ("valid", "valid.bed"),
        ("test", "test.bed"),
        ("all", "all_intervals.bed"),
    ):
        path = intervals_dir / filename
        record = file_record(path)
        chromosomes: set[str] = set()
        line_count = 0
        if path.is_file():
            with path.open() as handle:
                for line in handle:
                    if line.strip():
                        line_count += 1
                        chromosomes.add(line.split("\t", 1)[0])
        record.update(
            {
                "split": split,
                "interval_count": line_count,
                "chromosomes": ",".join(sorted(chromosomes)),
            }
        )
        rows.append(record)
    return rows


def dataset_inventory() -> list[dict[str, Any]]:
    rows = []
    for split in ("train", "valid", "test"):
        dataset_dir = REPO_ROOT / f"alphagenome_custom/datasets/rna_seq_npz_{split}"
        manifest = dataset_dir / "manifest.tsv"
        example_paths = sorted((dataset_dir / "examples").glob("*.npz"))
        rows.append(
            {
                "split": split,
                "dataset_path": rel(dataset_dir),
                "exists": str(dataset_dir.is_dir()),
                "example_count": len(example_paths),
                "total_example_bytes": sum(p.stat().st_size for p in example_paths),
                "manifest_path": rel(manifest),
                "manifest_sha256": sha256(manifest) if manifest.is_file() else "",
                "manifest_rows": len(read_tsv(manifest)) if manifest.is_file() else 0,
                "content_hash_scope": "manifest_only_large_npz_not_hashed",
            }
        )
    return rows


def candidate_inventory() -> list[dict[str, Any]]:
    rows = []
    for source in read_tsv(CANDIDATE_SOURCE):
        checkpoint = REPO_ROOT / source["checkpoint"]
        run_dir = checkpoint.parent
        config = run_dir / "config.json"
        metrics = run_dir / "metrics.tsv"
        rows.append(
            {
                "rank": source["rank"],
                "run_id": source["run_id"],
                "checkpoint": source["checkpoint"],
                "checkpoint_exists": str(checkpoint.is_file()),
                "checkpoint_size_bytes": checkpoint.stat().st_size
                if checkpoint.is_file()
                else "",
                "checkpoint_sha256": sha256(checkpoint)
                if checkpoint.is_file()
                else "",
                "config_exists": str(config.is_file()),
                "config_sha256": sha256(config) if config.is_file() else "",
                "metrics_exists": str(metrics.is_file()),
                "metrics_sha256": sha256(metrics) if metrics.is_file() else "",
                "best_step": source.get("best_step", ""),
                "full_mse": source.get("full_mse", ""),
                "pearson": source.get("pearson", ""),
            }
        )
    return rows


def selected_models() -> list[dict[str, Any]]:
    candidate_rows = read_tsv(CANDIDATE_SOURCE)
    utility_rows = read_tsv(UTILITY_SOURCE)
    mse = min(candidate_rows, key=lambda row: int(row["rank"]))
    utility = max(
        utility_rows, key=lambda row: float(row["representation_utility_score"])
    )
    selected = [
        (
            "base_weights",
            WEIGHTS_PATH,
            "Hugging Face gtca/alphagenome_pytorch; converted all-folds weights",
            "not_applicable",
        ),
        (
            "legacy_tested_adapter",
            REPO_ROOT
            / "runs/rna_seq11_adapter_formal_20260515_5000steps_lr3e-4_128bp/"
            "adapter_head_best.pt",
            "Early selected adapter with recorded chr X metrics",
            "read_once_recorded",
        ),
        (
            "validation_mse_winner",
            REPO_ROOT / mse["checkpoint"],
            mse["run_id"],
            "not_read_for_latest_selection",
        ),
        (
            "representation_utility_winner",
            REPO_ROOT / utility["checkpoint"],
            utility["run_id"],
            "not_read_for_latest_selection",
        ),
    ]
    rows = []
    for role, path, description, test_status in selected:
        record = file_record(path)
        record.update(
            {
                "role": role,
                "description": description,
                "chr_x_test_status": test_status,
            }
        )
        rows.append(record)
    return rows


def source_inventory() -> list[dict[str, Any]]:
    paths = [
        "alphagenome_custom/metadata/track_metadata.tsv",
        "alphagenome_custom/metadata/track_groups.tsv",
        "alphagenome_custom/metadata/track_metadata_grouped.tsv",
        "alphagenome_custom/metadata/group_signal_summary.tsv",
        "alphagenome_custom/metadata/grouped_bigwig_qc.tsv",
        "alphagenome_custom/metadata/interval_split_summary.tsv",
        "docs/rna_seq11_results_summary.md",
        "docs/rna_seq11_top15_extended_diagnostics_20260529.md",
        "docs/rna_seq11_all_representation_utility_20260530.md",
        "docs/experiment_log.md",
    ]
    return [file_record(REPO_ROOT / path) for path in paths]


def environment_inventory() -> list[dict[str, Any]]:
    packages = ("alphagenome-pytorch", "numpy", "pyBigWig", "torch")
    rows = [
        {"key": "snapshot_created_utc", "value": utc_now()},
        {"key": "hostname", "value": platform.node()},
        {"key": "python", "value": sys.version.replace("\n", " ")},
        {"key": "python_executable", "value": sys.executable},
        {"key": "git_branch", "value": git_output("branch", "--show-current")},
        {"key": "git_commit", "value": git_output("rev-parse", "HEAD")},
        {"key": "git_status_porcelain", "value": git_output("status", "--porcelain")},
    ]
    for package in packages:
        try:
            value = importlib_metadata.version(package)
        except importlib_metadata.PackageNotFoundError:
            value = "missing"
        rows.append({"key": f"package:{package}", "value": value})
    return rows


def write_freeze_document(summary: dict[str, Any]) -> None:
    path = REPO_ROOT / "docs/legacy_v1_freeze.md"
    path.write_text(
        f"""# Legacy RNA-seq11 v1 Freeze

Date: {summary['snapshot_created_utc']}

Host: `{summary['hostname']}`

Git commit: `{summary['git_commit']}` on `setup/agent-maintenance`

## Scope

This freeze records the legacy 2026-04-20 C. elegans RNA-seq11 data and model-development lineage without copying or modifying large assets. The authoritative machine-readable inventories are under `alphagenome_custom/metadata/legacy_v1/`.

## Frozen Evidence

- Original sample-level bigWigs: `{summary['raw_bigwig_count']}` actual files under `alphagenome_custom/tracks/bw_2026.4.20.training_input_bigwig`.
- Legacy grouped RNA-seq tracks: `{summary['grouped_bigwig_count']}`.
- Intervals: train `{summary['train_intervals']}`, valid `{summary['valid_intervals']}`, test `{summary['test_intervals']}`.
- Eligible validation-tuned checkpoint inventory: `{summary['candidate_checkpoint_count']}`.
- Raw-data before/after hash agreement during freeze: `{summary['raw_hashes_unchanged']}`.

## Interpretation Boundary

The legacy pipeline directly averaged existing bigWig scales and subsequently trained custom log1p regression heads on frozen human-index AlphaGenome PyTorch embeddings. It is an exploratory transfer-learning result, not a paper-faithful C. elegans AlphaGenome reproduction.

Chromosome V was reused for extensive checkpoint and objective selection. Chromosome X was read for the early selected 128 bp adapter and related GenomeTracks audits. The later 171-candidate representation benchmark did not read chromosome X, but chromosome X is not considered pristine for the project as a whole.

No new legacy training or dataset generation is permitted after this freeze. Corrections must be appended to the v2 execution log and must not rewrite the inventories silently.
"""
    )


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    AUDIT_DIR.mkdir(parents=True, exist_ok=True)
    started_at = utc_now()
    protected_roots = [RAW_DIR, GROUPED_DIR, REPO_ROOT / "runs"]
    protected_before = {
        rel(path): (path.stat().st_size, path.stat().st_mtime_ns)
        for root in protected_roots
        if root.exists()
        for path in root.rglob("*")
        if path.is_file()
    }

    raw_before = raw_inventory()
    grouped = grouped_inventory()
    references = reference_inventory()
    splits = split_inventory()
    datasets = dataset_inventory()
    candidates = candidate_inventory()
    selected = selected_models()
    sources = source_inventory()
    environment = environment_inventory()
    raw_after = raw_inventory()

    raw_before_hashes = {row["sample_id"]: row["sha256"] for row in raw_before}
    raw_after_hashes = {row["sample_id"]: row["sha256"] for row in raw_after}
    protected_after = {
        rel(path): (path.stat().st_size, path.stat().st_mtime_ns)
        for root in protected_roots
        if root.exists()
        for path in root.rglob("*")
        if path.is_file()
    }
    protected_changed = sorted(
        path
        for path in set(protected_before) | set(protected_after)
        if protected_before.get(path) != protected_after.get(path)
    )

    write_tsv(
        OUTPUT_DIR / "raw_bigwig_inventory.tsv",
        raw_after,
        [
            "sample_id",
            "path",
            "exists",
            "size_bytes",
            "sha256",
            "declared_legacy_path",
            "declared_legacy_path_exists",
            "normalization_status",
            "data_source",
            "stage",
            "tissue",
        ],
    )
    write_tsv(
        OUTPUT_DIR / "grouped_bigwig_inventory.tsv",
        grouped,
        [
            "group_id",
            "path",
            "exists",
            "size_bytes",
            "sha256",
            "sample_ids",
            "n_replicates",
            "input_files_declared",
            "normalization_plan",
            "normalization_status",
            "ontology_curie",
        ],
    )
    write_tsv(
        OUTPUT_DIR / "reference_inventory.tsv",
        references,
        ["asset_type", "path", "exists", "size_bytes", "sha256"],
    )
    write_tsv(
        OUTPUT_DIR / "split_inventory.tsv",
        splits,
        [
            "split",
            "path",
            "exists",
            "size_bytes",
            "sha256",
            "interval_count",
            "chromosomes",
        ],
    )
    write_tsv(
        OUTPUT_DIR / "dataset_inventory.tsv",
        datasets,
        [
            "split",
            "dataset_path",
            "exists",
            "example_count",
            "total_example_bytes",
            "manifest_path",
            "manifest_sha256",
            "manifest_rows",
            "content_hash_scope",
        ],
    )
    write_tsv(
        OUTPUT_DIR / "candidate_checkpoint_inventory.tsv",
        candidates,
        [
            "rank",
            "run_id",
            "checkpoint",
            "checkpoint_exists",
            "checkpoint_size_bytes",
            "checkpoint_sha256",
            "config_exists",
            "config_sha256",
            "metrics_exists",
            "metrics_sha256",
            "best_step",
            "full_mse",
            "pearson",
        ],
    )
    write_tsv(
        OUTPUT_DIR / "selected_model_inventory.tsv",
        selected,
        [
            "role",
            "path",
            "exists",
            "size_bytes",
            "sha256",
            "description",
            "chr_x_test_status",
        ],
    )
    write_tsv(
        OUTPUT_DIR / "source_file_inventory.tsv",
        sources,
        ["path", "exists", "size_bytes", "sha256"],
    )
    write_tsv(OUTPUT_DIR / "environment.tsv", environment, ["key", "value"])

    split_counts = {row["split"]: row["interval_count"] for row in splits}
    env = {row["key"]: row["value"] for row in environment}
    summary = {
        "schema_version": 1,
        "snapshot_created_utc": started_at,
        "hostname": platform.node(),
        "git_commit": env["git_commit"],
        "raw_bigwig_count": len(raw_after),
        "grouped_bigwig_count": len(grouped),
        "candidate_checkpoint_count": len(candidates),
        "train_intervals": split_counts.get("train", 0),
        "valid_intervals": split_counts.get("valid", 0),
        "test_intervals": split_counts.get("test", 0),
        "all_intervals": split_counts.get("all", 0),
        "raw_hashes_unchanged": raw_before_hashes == raw_after_hashes,
        "protected_files_changed_during_freeze": protected_changed,
        "candidate_source": rel(CANDIDATE_SOURCE),
        "candidate_source_sha256": sha256(CANDIDATE_SOURCE),
    }
    (OUTPUT_DIR / "freeze_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n"
    )
    write_freeze_document(summary)
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
