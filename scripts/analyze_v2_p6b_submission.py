#!/usr/bin/env python3
"""Create the locked P6B submission analysis and Nature-style figure bundle.

Data extraction and statistical analysis require only the Python standard
library and NumPy.  Matplotlib is imported lazily for figure rendering, so the
tables can be reproduced in the lean AlphaGenome environment and rendered in a
separate Python plotting environment without installing model dependencies.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
from dataclasses import dataclass
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import platform
import re
import sys
from typing import Any, Iterable, Mapping, Sequence

import numpy as np

from scripts.v2_submission_statistics import (
    hierarchical_block_bootstrap,
    make_hierarchical_resample_plan,
    matched_differences,
    records_to_balanced_matrix,
)


REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_SPEC = (
    REPO_ROOT
    / "alphagenome_custom/metadata/v2/p6b_submission_analysis_spec.json"
)
DEFAULT_OUTPUT = REPO_ROOT / "results/v2_p6b_submission_analysis"
UNKNOWN = "Missing/unknown"

RUN_ENDPOINT_EXTRACTORS: dict[str, tuple[str, ...]] = {
    "gene_exon_coverage_pearson_log1p": (
        "full_metrics",
        "primary",
        "mean_per_track_gene_exon_coverage_pearson_log1p",
    ),
    "per_track_pearson_128bp_log1p": (
        "full_metrics",
        "primary",
        "mean_per_track_pearson_128bp_log1p",
    ),
    "spearman_128bp": (
        "full_metrics",
        "distribution_128bp",
        "mean_per_track_spearman",
    ),
    "gene_body_mse_log1p_1bp": (
        "full_metrics",
        "gene_body_log1p_1bp",
        "mean_per_track_mse",
    ),
    "exon_mse_log1p_1bp": (
        "full_metrics",
        "exon_log1p_1bp",
        "mean_per_track_mse",
    ),
    "local_gradient_mse_log1p_1bp": (
        "full_metrics",
        "local_gradient_log1p_1bp",
        "mean_per_track_mse",
    ),
    "log1p_mse": ("mean_metrics", "log1p_mse"),
    "paper_loss": ("mean_metrics", "paper_loss"),
    "top1_mse_128bp": (
        "full_metrics",
        "distribution_128bp",
        "mean_per_track_top1_mse",
    ),
}

STRATIFIED_METRICS = (
    "primary.primary_biological_score",
    "gene_exon_coverage.pearson",
    "log1p_128bp_sum.pearson",
    "distribution_128bp.spearman",
)

DISPLAY_CONFIG = {
    "A/paper": "A, paper",
    "A/log1p_mse": "A, log1p MSE",
    "B/paper": "B, paper",
    "B/log1p_mse": "B, log1p MSE",
    "C/paper": "C, paper",
    "C/log1p_mse": "C, log1p MSE",
}


@dataclass
class RunningMoments:
    count: int = 0
    total: float = 0.0
    total_square: float = 0.0
    minimum: float = math.inf
    maximum: float = -math.inf

    def add(self, value: float) -> None:
        if not math.isfinite(value):
            return
        self.count += 1
        self.total += value
        self.total_square += value * value
        self.minimum = min(self.minimum, value)
        self.maximum = max(self.maximum, value)

    def summary(self) -> dict[str, float | int]:
        if self.count == 0:
            return {
                "n_finite": 0,
                "mean": math.nan,
                "sd": math.nan,
                "min": math.nan,
                "max": math.nan,
            }
        mean = self.total / self.count
        variance = (
            (self.total_square - self.count * mean * mean) / (self.count - 1)
            if self.count > 1
            else math.nan
        )
        return {
            "n_finite": self.count,
            "mean": mean,
            "sd": math.sqrt(max(0.0, variance)) if self.count > 1 else math.nan,
            "min": self.minimum,
            "max": self.maximum,
        }


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def validation_manifest_sha256(paths: Sequence[Path]) -> str:
    """Hash sorted validation hashes and repository-relative paths."""

    digest = hashlib.sha256()
    for path in sorted(paths):
        relative = path.resolve().relative_to(REPO_ROOT.resolve()).as_posix()
        digest.update(f"{sha256_file(path)}  {relative}\n".encode("utf-8"))
    return digest.hexdigest()


def require_hash(path: Path, expected: str) -> None:
    observed = sha256_file(path)
    if observed != expected:
        raise RuntimeError(
            f"input hash mismatch for {path.relative_to(REPO_ROOT)}: "
            f"expected {expected}, observed {observed}"
        )


def repo_path(relative: str) -> Path:
    candidate = (REPO_ROOT / relative).resolve()
    candidate.relative_to(REPO_ROOT.resolve())
    return candidate


def value_at(record: Mapping[str, Any], keys: Sequence[str]) -> Any:
    value: Any = record
    for key in keys:
        value = value[key]
    return value


def clean_text(value: Any) -> str:
    text = "" if value is None else str(value).strip()
    if not text or text.lower() in {".", "na", "n/a", "none", "null", "unknown"}:
        return UNKNOWN
    return text


def stage_stratum(value: Any) -> str:
    raw = clean_text(value)
    if raw == UNKNOWN:
        return UNKNOWN
    lower = raw.casefold()
    if lower.startswith("embryo") or lower.startswith("embryonic"):
        return "Embryo"
    if lower.startswith("l1-l3"):
        return "Larval time course"
    match = re.match(r"^l([1-4])(?:\b|\s|-)", lower)
    if match:
        return f"L{match.group(1)}"
    if lower == "young adult" or lower.startswith("young adult "):
        return "Young adult"
    if "dauer" in lower:
        return "Dauer"
    return raw


def tissue_stratum(value: Any) -> str:
    raw = clean_text(value)
    if raw == UNKNOWN:
        return UNKNOWN
    lower = raw.casefold()
    if "whole animal" in lower or "whole worm" in lower or "whole organism" in lower:
        return "Whole organism"
    if "embryo" in lower:
        return "Embryo"
    return raw


def treatment_stratum(condition: Any, genotype: Any, strain: Any) -> str:
    raw_condition = clean_text(condition)
    raw_genotype = clean_text(genotype)
    raw_strain = clean_text(strain)
    if raw_condition == raw_genotype == raw_strain == UNKNOWN:
        return UNKNOWN

    administrative = {"rna-seq", "synchronized l4 stage animals"}
    if (
        raw_condition != UNKNOWN
        and raw_condition.casefold() not in administrative
    ):
        return "Reported treatment/selection"
    explicit = " ".join(
        value.casefold()
        for value in (raw_genotype, raw_strain)
        if value != UNKNOWN
    )
    if "treated" in explicit or "rnai" in explicit:
        return "Reported treatment/selection"
    return "No reported treatment"


def write_delimited(
    path: Path,
    rows: Iterable[Mapping[str, Any]],
    fieldnames: Sequence[str],
    *,
    delimiter: str,
) -> int:
    path.parent.mkdir(parents=True, exist_ok=True)
    count = 0
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter=delimiter,
            extrasaction="ignore",
            lineterminator="\n",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    key: ""
                    if isinstance(row.get(key), float) and not math.isfinite(row[key])
                    else row.get(key, "")
                    for key in fieldnames
                }
            )
            count += 1
    return count


def read_delimited(path: Path, *, delimiter: str = "\t") -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def json_dump(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def flatten_scalar_metrics(validation: Mapping[str, Any]) -> dict[str, Any]:
    flattened: dict[str, Any] = {}
    for key, value in validation.get("mean_metrics", {}).items():
        if isinstance(value, (bool, int, float, str)):
            flattened[f"mean_metrics.{key}"] = value
    for section, payload in validation.get("full_metrics", {}).items():
        if not isinstance(payload, Mapping):
            continue
        for key, value in payload.items():
            if isinstance(value, (bool, int, float, str)):
                flattened[f"full_metrics.{section}.{key}"] = value
    return flattened


def load_track_metadata(
    manifest_path: Path, signal_path: Path
) -> tuple[dict[str, dict[str, Any]], dict[str, Any]]:
    manifest_rows = read_delimited(manifest_path)
    signal_rows = read_delimited(signal_path)
    signals = {row["group_id"]: float(row["development_train_nonzero_mean"]) for row in signal_rows}
    group_ids = [row["group_id"] for row in manifest_rows]
    if len(group_ids) != len(set(group_ids)):
        raise RuntimeError("duplicate group_id in track manifest")
    if set(group_ids) != set(signals):
        raise RuntimeError("track manifest and signal-mean IDs differ")

    signal_values = np.asarray([signals[group_id] for group_id in group_ids], dtype=float)
    q1, q2 = np.quantile(signal_values, [1.0 / 3.0, 2.0 / 3.0])
    metadata: dict[str, dict[str, Any]] = {}
    for row in manifest_rows:
        group_id = row["group_id"]
        signal = signals[group_id]
        if signal <= q1:
            signal_label = "Low"
        elif signal <= q2:
            signal_label = "Middle"
        else:
            signal_label = "High"
        metadata[group_id] = {
            "group_id": group_id,
            "study_raw": clean_text(row.get("sra_study_accession")),
            "stage_raw": clean_text(row.get("development_stage")),
            "tissue_raw": clean_text(row.get("tissue_or_cell_type")),
            "condition_raw": clean_text(row.get("condition")),
            "genotype_raw": clean_text(row.get("genotype")),
            "strain_raw": clean_text(row.get("strain")),
            "study_stratum": clean_text(row.get("sra_study_accession")),
            "stage_stratum": stage_stratum(row.get("development_stage")),
            "tissue_stratum": tissue_stratum(row.get("tissue_or_cell_type")),
            "treatment_stratum": treatment_stratum(
                row.get("condition"), row.get("genotype"), row.get("strain")
            ),
            "signal_strength_value": signal,
            "signal_strength_stratum": signal_label,
            "n_runs": row.get("n_runs", ""),
            "n_biological_units": row.get("n_biological_units", ""),
            "bioproject_accession": clean_text(row.get("bioproject_accession")),
            "geo_accession": clean_text(row.get("geo_accession")),
        }
    counts = Counter(row["signal_strength_stratum"] for row in metadata.values())
    signal_audit = {
        "method": "linear tertiles of development_train_nonzero_mean",
        "lower_boundary": float(q1),
        "upper_boundary": float(q2),
        "assignment": "Low <= lower; Middle > lower and <= upper; High > upper",
        "counts": dict(sorted(counts.items())),
    }
    return metadata, signal_audit


def load_and_validate_inputs(spec: Mapping[str, Any]) -> tuple[list[dict[str, Any]], dict[str, dict[str, Any]], dict[str, Any]]:
    contract = spec["input_contract"]
    for path_key, hash_key in (
        ("execution_record", "execution_record_sha256"),
        ("training_spec", "training_spec_sha256"),
        ("track_manifest", "track_manifest_sha256"),
        ("track_signal_means", "track_signal_means_sha256"),
    ):
        require_hash(repo_path(contract[path_key]), contract[hash_key])

    execution = json.loads(repo_path(contract["execution_record"]).read_text())
    jobs = execution.get("formal_jobs", [])
    if len(jobs) != int(contract["expected_jobs"]):
        raise RuntimeError(f"expected {contract['expected_jobs']} formal jobs, found {len(jobs)}")

    formal_root = repo_path(contract["validation_root"])
    prohibited = tuple(contract["prohibited_inputs"])
    paths: list[Path] = []
    loaded: list[dict[str, Any]] = []
    seen: set[tuple[str, str, int, int]] = set()
    for job in jobs:
        if job.get("status") != "completed":
            raise RuntimeError(f"formal job is not completed: {job.get('job_id')}")
        relative = str(job["validation_path"])
        if any(token in relative for token in prohibited):
            raise RuntimeError(f"prohibited analysis input: {relative}")
        path = repo_path(relative)
        path.relative_to(formal_root)
        if path.name != contract["validation_filename"] or not path.is_file():
            raise RuntimeError(f"missing or unexpected validation input: {relative}")
        validation = json.loads(path.read_text())
        key = (
            str(validation["model"]),
            str(validation["training_loss"]),
            int(validation["seed"]),
            int(validation["fold"]),
        )
        expected_key = (
            str(job["model"]),
            str(job["loss"]),
            int(job["seed"]),
            int(job["fold"]),
        )
        if key != expected_key:
            raise RuntimeError(f"execution/validation identity mismatch at {relative}")
        if key in seen:
            raise RuntimeError(f"duplicate formal validation: {key}")
        seen.add(key)
        guards = contract["required_validation_guards"]
        if validation.get("locked_test_block_signal_reads") is not guards["locked_test_block_signal_reads"]:
            raise RuntimeError(f"locked-test read guard failed for {relative}")
        if float(validation.get("validation_core_coverage_fraction", 0.0)) < float(
            guards["minimum_validation_core_coverage_fraction"]
        ):
            raise RuntimeError(f"validation coverage guard failed for {relative}")
        paths.append(path)
        loaded.append({"job": job, "validation": validation, "path": path})

    expected_keys = {
        (model, loss, seed, fold)
        for model in contract["models"]
        for loss in contract["training_losses"]
        for seed in contract["seeds"]
        for fold in contract["folds"]
    }
    if seen != expected_keys:
        raise RuntimeError(
            f"formal matrix mismatch; missing={sorted(expected_keys-seen)}, extra={sorted(seen-expected_keys)}"
        )
    manifest_digest = validation_manifest_sha256(paths)
    if manifest_digest != contract["validation_manifest_sha256"]:
        raise RuntimeError(
            "validation manifest digest mismatch: "
            f"expected {contract['validation_manifest_sha256']}, observed {manifest_digest}"
        )

    track_metadata, signal_audit = load_track_metadata(
        repo_path(contract["track_manifest"]), repo_path(contract["track_signal_means"])
    )
    if len(track_metadata) != int(contract["expected_tracks"]):
        raise RuntimeError(
            f"expected {contract['expected_tracks']} tracks, found {len(track_metadata)}"
        )
    return loaded, track_metadata, signal_audit


def make_run_row(item: Mapping[str, Any]) -> dict[str, Any]:
    job = item["job"]
    validation = item["validation"]
    gene = float(value_at(validation, RUN_ENDPOINT_EXTRACTORS["gene_exon_coverage_pearson_log1p"]))
    resolution = float(value_at(validation, RUN_ENDPOINT_EXTRACTORS["per_track_pearson_128bp_log1p"]))
    row: dict[str, Any] = {
        "job_id": job["job_id"],
        "config": f"{validation['model']}/{validation['training_loss']}",
        "model": validation["model"],
        "training_loss": validation["training_loss"],
        "fold": int(validation["fold"]),
        "seed": int(validation["seed"]),
        "validation_path": item["path"].relative_to(REPO_ROOT).as_posix(),
        "validation_sha256": sha256_file(item["path"]),
        "checkpoint_sha256": validation.get("checkpoint_sha256", ""),
        "elapsed_seconds": job.get("elapsed_seconds", ""),
        "validation_bases": validation.get("validation_bases", ""),
        "validation_subwindows": validation.get("validation_subwindows", ""),
        "validation_core_coverage_fraction": validation.get(
            "validation_core_coverage_fraction", ""
        ),
        "locked_test_block_signal_reads": validation.get(
            "locked_test_block_signal_reads", ""
        ),
        "primary_biological_score": 0.5 * gene + 0.5 * resolution,
    }
    for endpoint, keys in RUN_ENDPOINT_EXTRACTORS.items():
        row[endpoint] = float(value_at(validation, keys))
    row.update(flatten_scalar_metrics(validation))
    return row


def iter_track_values(
    validation: Mapping[str, Any], expected_group_ids: set[str]
) -> Iterable[tuple[str, str, str, str, float]]:
    """Yield group, metric ID, section, metric name and numeric value."""

    full_metrics = validation["full_metrics"]
    for section in sorted(full_metrics):
        payload = full_metrics[section]
        if not isinstance(payload, Mapping):
            continue
        for key in sorted(payload):
            values = payload[key]
            if not isinstance(values, Mapping):
                continue
            if set(values) != expected_group_ids:
                continue
            if key.startswith("per_track_"):
                metric = key[len("per_track_") :]
            elif key == "positions_per_track":
                metric = "positions"
            else:
                metric = key
            metric_id = f"{section}.{metric}"
            if metric == "positions":
                role = "evaluation_count"
            elif "target_threshold" in metric:
                role = "target_reference"
            else:
                role = "performance_metric"
            for group_id in sorted(values):
                yield group_id, metric_id, section, metric, float(values[group_id]), role

    gene = full_metrics["gene_exon_coverage"]["per_track_pearson"]
    resolution = full_metrics["log1p_128bp_sum"]["per_track_pearson"]
    if set(gene) != expected_group_ids or set(resolution) != expected_group_ids:
        raise RuntimeError("primary per-track component IDs do not match manifest")
    for group_id in sorted(expected_group_ids):
        value = 0.5 * float(gene[group_id]) + 0.5 * float(resolution[group_id])
        yield (
            group_id,
            "primary.primary_biological_score",
            "primary",
            "primary_biological_score",
            value,
            "derived_performance_metric",
        )


def quantile_summary(values: Sequence[float]) -> dict[str, Any]:
    finite = np.asarray([value for value in values if math.isfinite(value)], dtype=float)
    if finite.size == 0:
        return {
            "n_tracks": 0,
            "mean": math.nan,
            "median": math.nan,
            "q1": math.nan,
            "q3": math.nan,
            "min": math.nan,
            "max": math.nan,
        }
    q1, median, q3 = np.quantile(finite, [0.25, 0.5, 0.75])
    return {
        "n_tracks": int(finite.size),
        "mean": float(finite.mean()),
        "median": float(median),
        "q1": float(q1),
        "q3": float(q3),
        "min": float(finite.min()),
        "max": float(finite.max()),
    }


def metric_direction(metric_id: str) -> str:
    if metric_id.endswith(("pearson", "spearman", "primary_biological_score")):
        return "higher"
    if metric_id.endswith(("mse", "mae")):
        return "lower"
    if "calibration_ratio" in metric_id:
        return "closer_to_one"
    return "descriptive"


def config_order(spec: Mapping[str, Any]) -> dict[str, int]:
    order = []
    for model in spec["input_contract"]["models"]:
        for loss in spec["input_contract"]["training_losses"]:
            order.append(f"{model}/{loss}")
    return {config: index for index, config in enumerate(order)}


def analyse(spec: Mapping[str, Any], output_dir: Path) -> dict[str, Any]:
    loaded, track_metadata, signal_audit = load_and_validate_inputs(spec)
    output_dir.mkdir(parents=True, exist_ok=True)
    order = config_order(spec)
    endpoints = {row["id"]: row for row in spec["run_level_endpoints"]}
    folds = tuple(spec["input_contract"]["folds"])
    seeds = tuple(spec["input_contract"]["seeds"])
    bootstrap_spec = spec["bootstrap"]
    plan = make_hierarchical_resample_plan(
        n_blocks=len(folds),
        n_repeats=len(seeds),
        iterations=int(bootstrap_spec["iterations"]),
        seed=int(bootstrap_spec["seed"]),
    )

    run_rows = [make_run_row(item) for item in loaded]
    run_rows.sort(key=lambda row: (order[row["config"]], row["fold"], row["seed"]))
    run_fields = list(run_rows[0])
    run_count = write_delimited(
        output_dir / "run_level_metrics.tsv", run_rows, run_fields, delimiter="\t"
    )
    write_delimited(
        output_dir / "run_level_metrics.csv", run_rows, run_fields, delimiter=","
    )

    track_metadata_fields = list(next(iter(track_metadata.values())))
    track_metadata_rows = [track_metadata[key] for key in sorted(track_metadata)]
    write_delimited(
        output_dir / "track_metadata_strata.tsv",
        track_metadata_rows,
        track_metadata_fields,
        delimiter="\t",
    )

    track_long_path = output_dir / "track_metrics_long.tsv"
    track_long_fields = [
        "config",
        "model",
        "training_loss",
        "fold",
        "seed",
        "group_id",
        "metric_id",
        "section",
        "metric",
        "value",
        "finite",
        "value_role",
    ]
    expected_group_ids = set(track_metadata)
    accumulators: dict[tuple[str, str, str], RunningMoments] = defaultdict(RunningMoments)
    metric_ids_seen: set[str] = set()
    track_long_count = 0
    with track_long_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=track_long_fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for item in sorted(
            loaded,
            key=lambda value: (
                order[f"{value['validation']['model']}/{value['validation']['training_loss']}"],
                int(value["validation"]["fold"]),
                int(value["validation"]["seed"]),
            ),
        ):
            validation = item["validation"]
            config = f"{validation['model']}/{validation['training_loss']}"
            for group_id, metric_id, section, metric, value, role in iter_track_values(
                validation, expected_group_ids
            ):
                finite = math.isfinite(value)
                writer.writerow(
                    {
                        "config": config,
                        "model": validation["model"],
                        "training_loss": validation["training_loss"],
                        "fold": validation["fold"],
                        "seed": validation["seed"],
                        "group_id": group_id,
                        "metric_id": metric_id,
                        "section": section,
                        "metric": metric,
                        "value": value if finite else "",
                        "finite": finite,
                        "value_role": role,
                    }
                )
                accumulators[(config, group_id, metric_id)].add(value)
                metric_ids_seen.add(metric_id)
                track_long_count += 1

    track_summary_rows: list[dict[str, Any]] = []
    track_lookup: dict[tuple[str, str, str], float] = {}
    for (config, group_id, metric_id), moments in sorted(
        accumulators.items(), key=lambda value: (order[value[0][0]], value[0][1], value[0][2])
    ):
        summary = moments.summary()
        track_lookup[(config, group_id, metric_id)] = float(summary["mean"])
        track_summary_rows.append(
            {
                "config": config,
                "model": config.split("/")[0],
                "training_loss": config.split("/")[1],
                "group_id": group_id,
                "metric_id": metric_id,
                "direction": metric_direction(metric_id),
                "n_finite_fold_seed_values": summary["n_finite"],
                "n_expected_fold_seed_values": len(folds) * len(seeds),
                "mean": summary["mean"],
                "sd": summary["sd"],
                "min": summary["min"],
                "max": summary["max"],
            }
        )
    track_summary_fields = list(track_summary_rows[0])
    track_summary_count = write_delimited(
        output_dir / "track_configuration_summaries.tsv",
        track_summary_rows,
        track_summary_fields,
        delimiter="\t",
    )

    configuration_rows: list[dict[str, Any]] = []
    for config in sorted(order, key=order.get):
        config_records = [row for row in run_rows if row["config"] == config]
        for metric_id, metric_spec in endpoints.items():
            matrix = records_to_balanced_matrix(
                config_records,
                value_key=metric_id,
                expected_blocks=folds,
                expected_repeats=seeds,
            )
            estimate = hierarchical_block_bootstrap(
                matrix.values,
                plan,
                confidence_level=float(bootstrap_spec["confidence_level"]),
            )
            fold_means = matrix.values.mean(axis=1)
            configuration_rows.append(
                {
                    "config": config,
                    "model": config.split("/")[0],
                    "training_loss": config.split("/")[1],
                    "metric_id": metric_id,
                    "metric_role": metric_spec["role"],
                    "direction": metric_spec["direction"],
                    "estimate": estimate.estimate,
                    "ci_lower": estimate.ci_lower,
                    "ci_upper": estimate.ci_upper,
                    "confidence_level": estimate.confidence_level,
                    "n_folds": len(folds),
                    "n_seeds_per_fold": len(seeds),
                    "n_runs": len(config_records),
                    "bootstrap_iterations": estimate.iterations,
                    "bootstrap_seed": estimate.seed,
                    "fold_means": ";".join(f"{value:.12g}" for value in fold_means),
                }
            )
    configuration_fields = list(configuration_rows[0])
    configuration_count = write_delimited(
        output_dir / "configuration_bootstrap_estimates.tsv",
        configuration_rows,
        configuration_fields,
        delimiter="\t",
    )

    contrast_pair_rows: list[dict[str, Any]] = []
    contrast_summary_rows: list[dict[str, Any]] = []
    for contrast in spec["paired_contrasts"]:
        candidate_records = [
            row for row in run_rows if row["config"] == contrast["candidate"]
        ]
        comparator_records = [
            row for row in run_rows if row["config"] == contrast["comparator"]
        ]
        for metric_id, metric_spec in endpoints.items():
            pairs = matched_differences(
                candidate_records, comparator_records, value_key=metric_id
            )
            direction_sign = 1.0 if metric_spec["direction"] == "higher" else -1.0
            for pair in pairs:
                contrast_pair_rows.append(
                    {
                        "contrast_id": contrast["id"],
                        "candidate": contrast["candidate"],
                        "comparator": contrast["comparator"],
                        "contrast_family": contrast["family"],
                        "priority": contrast["priority"],
                        "metric_id": metric_id,
                        "metric_role": metric_spec["role"],
                        "direction": metric_spec["direction"],
                        "fold": pair["fold"],
                        "seed": pair["seed"],
                        "candidate_value": pair["candidate_value"],
                        "comparator_value": pair["comparator_value"],
                        "raw_difference_candidate_minus_comparator": pair["difference"],
                        "favorable_effect": direction_sign * pair["difference"],
                    }
                )
            difference_matrix = records_to_balanced_matrix(
                pairs,
                value_key="difference",
                expected_blocks=folds,
                expected_repeats=seeds,
            )
            estimate = hierarchical_block_bootstrap(
                difference_matrix.values,
                plan,
                confidence_level=float(bootstrap_spec["confidence_level"]),
            )
            if direction_sign > 0:
                favorable_lower, favorable_upper = estimate.ci_lower, estimate.ci_upper
            else:
                favorable_lower, favorable_upper = -estimate.ci_upper, -estimate.ci_lower
            contrast_summary_rows.append(
                {
                    "contrast_id": contrast["id"],
                    "candidate": contrast["candidate"],
                    "comparator": contrast["comparator"],
                    "contrast_family": contrast["family"],
                    "priority": contrast["priority"],
                    "metric_id": metric_id,
                    "metric_role": metric_spec["role"],
                    "direction": metric_spec["direction"],
                    "raw_difference_estimate": estimate.estimate,
                    "raw_difference_ci_lower": estimate.ci_lower,
                    "raw_difference_ci_upper": estimate.ci_upper,
                    "favorable_effect_estimate": direction_sign * estimate.estimate,
                    "favorable_effect_ci_lower": favorable_lower,
                    "favorable_effect_ci_upper": favorable_upper,
                    "confidence_level": estimate.confidence_level,
                    "n_folds": len(folds),
                    "n_seeds_per_fold": len(seeds),
                    "n_paired_runs": len(pairs),
                    "bootstrap_iterations": estimate.iterations,
                    "bootstrap_seed": estimate.seed,
                }
            )
    contrast_pair_fields = list(contrast_pair_rows[0])
    contrast_pair_count = write_delimited(
        output_dir / "paired_fold_seed_contrasts.tsv",
        contrast_pair_rows,
        contrast_pair_fields,
        delimiter="\t",
    )
    contrast_summary_fields = list(contrast_summary_rows[0])
    contrast_summary_count = write_delimited(
        output_dir / "paired_contrast_bootstrap_estimates.tsv",
        contrast_summary_rows,
        contrast_summary_fields,
        delimiter="\t",
    )

    stratum_fields = {
        "study": "study_stratum",
        "stage": "stage_stratum",
        "tissue": "tissue_stratum",
        "treatment": "treatment_stratum",
        "signal_strength": "signal_strength_stratum",
    }
    stratified_rows: list[dict[str, Any]] = []
    categories_by_dimension: dict[str, list[str]] = {}
    for dimension, metadata_field in stratum_fields.items():
        categories = sorted({row[metadata_field] for row in track_metadata.values()})
        categories_by_dimension[dimension] = categories
        for category in categories:
            group_ids = [
                group_id
                for group_id, metadata in track_metadata.items()
                if metadata[metadata_field] == category
            ]
            for config in sorted(order, key=order.get):
                for metric_id in STRATIFIED_METRICS:
                    values = [
                        track_lookup[(config, group_id, metric_id)] for group_id in group_ids
                    ]
                    stratified_rows.append(
                        {
                            "summary_type": "configuration",
                            "dimension": dimension,
                            "category": category,
                            "config_or_contrast": config,
                            "metric_id": metric_id,
                            **quantile_summary(values),
                            "inference": "descriptive_tracks_not_independent_n",
                        }
                    )
            for contrast_id, candidate, comparator in (
                ("B_paper_vs_A_paper", "B/paper", "A/paper"),
                ("B_paper_vs_C_paper", "B/paper", "C/paper"),
            ):
                for metric_id in STRATIFIED_METRICS:
                    values = [
                        track_lookup[(candidate, group_id, metric_id)]
                        - track_lookup[(comparator, group_id, metric_id)]
                        for group_id in group_ids
                    ]
                    stratified_rows.append(
                        {
                            "summary_type": "track_level_contrast",
                            "dimension": dimension,
                            "category": category,
                            "config_or_contrast": contrast_id,
                            "metric_id": metric_id,
                            **quantile_summary(values),
                            "inference": "descriptive_tracks_not_independent_n",
                        }
                    )
    stratified_fields = list(stratified_rows[0])
    stratified_count = write_delimited(
        output_dir / "stratified_track_summaries.tsv",
        stratified_rows,
        stratified_fields,
        delimiter="\t",
    )

    source_rows = build_source_data(
        spec=spec,
        run_rows=run_rows,
        configuration_rows=configuration_rows,
        contrast_pair_rows=contrast_pair_rows,
        contrast_summary_rows=contrast_summary_rows,
        track_metadata=track_metadata,
        track_lookup=track_lookup,
    )
    source_fields = [
        "panel",
        "data_role",
        "series",
        "subgroup",
        "metric_id",
        "fold",
        "seed",
        "group_id",
        "value",
        "ci_lower",
        "ci_upper",
        "n",
        "notes",
    ]
    source_count = write_delimited(
        output_dir / "Source_Data_Figure_1.tsv",
        source_rows,
        source_fields,
        delimiter="\t",
    )
    write_delimited(
        output_dir / "Source_Data_Figure_1.csv",
        source_rows,
        source_fields,
        delimiter=",",
    )

    primary_config = {
        row["config"]: {
            "estimate": row["estimate"],
            "ci_lower": row["ci_lower"],
            "ci_upper": row["ci_upper"],
        }
        for row in configuration_rows
        if row["metric_id"] == "primary_biological_score"
    }
    primary_contrasts = {
        row["contrast_id"]: {
            "favorable_effect_estimate": row["favorable_effect_estimate"],
            "ci_lower": row["favorable_effect_ci_lower"],
            "ci_upper": row["favorable_effect_ci_upper"],
        }
        for row in contrast_summary_rows
        if row["priority"] == "primary"
        and row["metric_id"] == "primary_biological_score"
    }
    counts = {
        "formal_validation_jsons": len(loaded),
        "run_level_rows": run_count,
        "tracks": len(track_metadata),
        "per_track_metric_ids_including_derived_primary": len(metric_ids_seen),
        "track_long_rows": track_long_count,
        "track_configuration_summary_rows": track_summary_count,
        "configuration_bootstrap_rows": configuration_count,
        "paired_fold_seed_rows": contrast_pair_count,
        "paired_contrast_bootstrap_rows": contrast_summary_count,
        "stratified_summary_rows": stratified_count,
        "source_data_rows": source_count,
    }
    manifest = {
        "analysis_id": spec["analysis_id"],
        "analysis_spec": DEFAULT_SPEC.relative_to(REPO_ROOT).as_posix(),
        "analysis_spec_sha256": sha256_file(DEFAULT_SPEC),
        "created_at": utc_now(),
        "status": "analysis_complete_figure_pending",
        "scope": "existing P6B development-validation summaries only",
        "locked_test_block_signal_reads": 0,
        "counts": counts,
        "signal_strength_tertiles": signal_audit,
        "stratum_categories": categories_by_dimension,
        "bootstrap": spec["bootstrap"],
        "software": {
            "analysis_python": platform.python_version(),
            "analysis_numpy": np.__version__,
            "analysis_executable": sys.executable,
        },
        "primary_configuration_estimates": primary_config,
        "primary_contrast_estimates": primary_contrasts,
    }
    write_methods(output_dir, spec, counts, signal_audit)
    write_legend(output_dir)
    write_readme(output_dir, counts)
    generated_files = sorted(
        path
        for path in output_dir.iterdir()
        if path.is_file()
        and path.name != "analysis_manifest.json"
        and not path.name.startswith("Figure_1_")
    )
    manifest["software"]["analysis_script_sha256"] = sha256_file(Path(__file__))
    manifest["software"]["statistics_helper_sha256"] = sha256_file(
        REPO_ROOT / "scripts/v2_submission_statistics.py"
    )
    manifest["analysis_output_files"] = {
        path.name: {"sha256": sha256_file(path), "bytes": path.stat().st_size}
        for path in generated_files
    }
    json_dump(output_dir / "analysis_manifest.json", manifest)
    return manifest


def select_panel_d_categories(
    track_metadata: Mapping[str, Mapping[str, Any]]
) -> list[tuple[str, str, str]]:
    fields = {
        "study": "study_stratum",
        "stage": "stage_stratum",
        "tissue": "tissue_stratum",
        "treatment": "treatment_stratum",
    }
    selected: list[tuple[str, str, str]] = []
    for dimension, field in fields.items():
        counts = Counter(row[field] for row in track_metadata.values())
        eligible = sorted(
            ((-count, category) for category, count in counts.items() if count >= 5)
        )[:2]
        for negative_count, category in eligible:
            selected.append((dimension, category, field))
    for category in ("Low", "Middle", "High"):
        selected.append(("signal_strength", category, "signal_strength_stratum"))
    return selected


def build_source_data(
    *,
    spec: Mapping[str, Any],
    run_rows: Sequence[Mapping[str, Any]],
    configuration_rows: Sequence[Mapping[str, Any]],
    contrast_pair_rows: Sequence[Mapping[str, Any]],
    contrast_summary_rows: Sequence[Mapping[str, Any]],
    track_metadata: Mapping[str, Mapping[str, Any]],
    track_lookup: Mapping[tuple[str, str, str], float],
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for run in run_rows:
        rows.append(
            {
                "panel": "a",
                "data_role": "run",
                "series": run["config"],
                "metric_id": "primary_biological_score",
                "fold": run["fold"],
                "seed": run["seed"],
                "value": run["primary_biological_score"],
                "notes": "algorithmic seed within development-validation fold",
            }
        )
    for config in sorted({row["config"] for row in run_rows}):
        for fold in spec["input_contract"]["folds"]:
            values = [
                float(row["primary_biological_score"])
                for row in run_rows
                if row["config"] == config and int(row["fold"]) == int(fold)
            ]
            rows.append(
                {
                    "panel": "a",
                    "data_role": "fold_mean",
                    "series": config,
                    "metric_id": "primary_biological_score",
                    "fold": fold,
                    "value": float(np.mean(values)),
                    "n": len(values),
                    "notes": "mean of three training seeds",
                }
            )
    for estimate in configuration_rows:
        if estimate["metric_id"] == "primary_biological_score":
            rows.append(
                {
                    "panel": "a",
                    "data_role": "estimate",
                    "series": estimate["config"],
                    "metric_id": estimate["metric_id"],
                    "value": estimate["estimate"],
                    "ci_lower": estimate["ci_lower"],
                    "ci_upper": estimate["ci_upper"],
                    "n": estimate["n_folds"],
                    "notes": "fold-primary hierarchical bootstrap 95% CI",
                }
            )

    panel_b_metrics = {
        "primary_biological_score",
        "gene_exon_coverage_pearson_log1p",
        "per_track_pearson_128bp_log1p",
    }
    for pair in contrast_pair_rows:
        if pair["priority"] == "primary" and pair["metric_id"] in panel_b_metrics:
            rows.append(
                {
                    "panel": "b",
                    "data_role": "paired_run_effect",
                    "series": pair["contrast_id"],
                    "subgroup": pair["metric_id"],
                    "metric_id": pair["metric_id"],
                    "fold": pair["fold"],
                    "seed": pair["seed"],
                    "value": pair["favorable_effect"],
                    "notes": "positive values favor the candidate",
                }
            )
    for contrast in spec["paired_contrasts"]:
        if contrast["priority"] != "primary":
            continue
        for metric_id in panel_b_metrics:
            pairs = [
                row
                for row in contrast_pair_rows
                if row["contrast_id"] == contrast["id"] and row["metric_id"] == metric_id
            ]
            for fold in spec["input_contract"]["folds"]:
                values = [
                    float(row["favorable_effect"])
                    for row in pairs
                    if int(row["fold"]) == int(fold)
                ]
                rows.append(
                    {
                        "panel": "b",
                        "data_role": "fold_mean_effect",
                        "series": contrast["id"],
                        "subgroup": metric_id,
                        "metric_id": metric_id,
                        "fold": fold,
                        "value": float(np.mean(values)),
                        "n": len(values),
                        "notes": "paired effect averaged across three seeds",
                    }
                )
    for estimate in contrast_summary_rows:
        if estimate["priority"] == "primary" and estimate["metric_id"] in panel_b_metrics:
            rows.append(
                {
                    "panel": "b",
                    "data_role": "estimate",
                    "series": estimate["contrast_id"],
                    "subgroup": estimate["metric_id"],
                    "metric_id": estimate["metric_id"],
                    "value": estimate["favorable_effect_estimate"],
                    "ci_lower": estimate["favorable_effect_ci_lower"],
                    "ci_upper": estimate["favorable_effect_ci_upper"],
                    "n": estimate["n_folds"],
                    "notes": "fold-primary hierarchical bootstrap 95% CI",
                }
            )

    metric_id = "primary.primary_biological_score"
    for config in ("A/paper", "B/paper", "C/paper"):
        values = []
        for group_id in sorted(track_metadata):
            value = float(track_lookup[(config, group_id, metric_id)])
            values.append(value)
            rows.append(
                {
                    "panel": "c",
                    "data_role": "track_mean",
                    "series": config,
                    "metric_id": metric_id,
                    "group_id": group_id,
                    "value": value,
                    "notes": "mean across five folds and three seeds; descriptive track",
                }
            )
        summary = quantile_summary(values)
        rows.append(
            {
                "panel": "c",
                "data_role": "median_iqr",
                "series": config,
                "metric_id": metric_id,
                "value": summary["median"],
                "ci_lower": summary["q1"],
                "ci_upper": summary["q3"],
                "n": summary["n_tracks"],
                "notes": "median and interquartile range; descriptive only",
            }
        )

    for dimension, category, metadata_field in select_panel_d_categories(track_metadata):
        differences = []
        for group_id, metadata in sorted(track_metadata.items()):
            if metadata[metadata_field] != category:
                continue
            difference = float(track_lookup[("B/paper", group_id, metric_id)]) - float(
                track_lookup[("A/paper", group_id, metric_id)]
            )
            differences.append(difference)
            rows.append(
                {
                    "panel": "d",
                    "data_role": "track_difference",
                    "series": dimension,
                    "subgroup": category,
                    "metric_id": metric_id,
                    "group_id": group_id,
                    "value": difference,
                    "notes": "B/paper minus A/paper; descriptive track",
                }
            )
        summary = quantile_summary(differences)
        rows.append(
            {
                "panel": "d",
                "data_role": "median_iqr",
                "series": dimension,
                "subgroup": category,
                "metric_id": metric_id,
                "value": summary["median"],
                "ci_lower": summary["q1"],
                "ci_upper": summary["q3"],
                "n": summary["n_tracks"],
                "notes": "median and interquartile range; descriptive only",
            }
        )
    return rows


def as_float(row: Mapping[str, str], key: str) -> float:
    return float(row[key]) if row.get(key, "") != "" else math.nan


def render_figure(output_dir: Path) -> dict[str, Any]:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.lines import Line2D
    except ImportError as error:
        raise RuntimeError(
            "Figure rendering requires matplotlib; analysis tables are already available."
        ) from error

    matplotlib.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
            "font.size": 6.5,
            "axes.labelsize": 6.5,
            "axes.titlesize": 7,
            "xtick.labelsize": 5.8,
            "ytick.labelsize": 5.8,
            "legend.fontsize": 5.5,
            "axes.linewidth": 0.65,
            "xtick.major.width": 0.55,
            "ytick.major.width": 0.55,
            "xtick.major.size": 2.5,
            "ytick.major.size": 2.5,
            "axes.spines.right": False,
            "axes.spines.top": False,
            "legend.frameon": False,
        }
    )
    rows = read_delimited(output_dir / "Source_Data_Figure_1.tsv")
    palette = {"A": "#777777", "B": "#007F73", "C": "#C65D3A"}
    metric_colors = {
        "primary_biological_score": "#007F73",
        "gene_exon_coverage_pearson_log1p": "#3C6E9E",
        "per_track_pearson_128bp_log1p": "#C65D3A",
    }
    metric_labels = {
        "primary_biological_score": "Composite",
        "gene_exon_coverage_pearson_log1p": "Gene-exon",
        "per_track_pearson_128bp_log1p": "128-bp",
    }
    contrast_labels = {
        "B_paper_vs_A_paper": "B paper - A paper",
        "B_paper_vs_C_paper": "B paper - C paper",
        "B_paper_vs_B_log1p_mse": "B paper - B log1p MSE",
    }

    fig_width_mm = 183.0
    fig_height_mm = 164.0
    width_inches = fig_width_mm / 25.4
    height_inches = fig_height_mm / 25.4
    fig = plt.figure(figsize=(width_inches, height_inches), facecolor="white")
    grid = fig.add_gridspec(
        2,
        2,
        width_ratios=(0.92, 1.08),
        height_ratios=(1.0, 1.08),
        left=0.095,
        right=0.985,
        bottom=0.09,
        top=0.955,
        wspace=0.50,
        hspace=0.38,
    )
    axes = [
        fig.add_subplot(grid[0, 0]),
        fig.add_subplot(grid[0, 1]),
        fig.add_subplot(grid[1, 0]),
        fig.add_subplot(grid[1, 1]),
    ]
    ax_a, ax_b, ax_c, ax_d = axes

    config_order_values = [
        "A/log1p_mse",
        "A/paper",
        "B/log1p_mse",
        "B/paper",
        "C/log1p_mse",
        "C/paper",
    ]
    y_positions = np.arange(len(config_order_values))[::-1]
    for y, config in zip(y_positions, config_order_values):
        model, loss = config.split("/")
        estimate = next(
            row
            for row in rows
            if row["panel"] == "a"
            and row["data_role"] == "estimate"
            and row["series"] == config
        )
        fold_rows = [
            row
            for row in rows
            if row["panel"] == "a"
            and row["data_role"] == "fold_mean"
            and row["series"] == config
        ]
        fold_values = [as_float(row, "value") for row in fold_rows]
        offsets = np.linspace(-0.11, 0.11, len(fold_values))
        marker = "o" if loss == "paper" else "s"
        ax_a.scatter(
            fold_values,
            y + offsets,
            s=9,
            marker=marker,
            facecolors=palette[model] if loss == "paper" else "white",
            edgecolors=palette[model],
            linewidths=0.65,
            alpha=0.72,
            zorder=2,
        )
        center = as_float(estimate, "value")
        lower = as_float(estimate, "ci_lower")
        upper = as_float(estimate, "ci_upper")
        ax_a.errorbar(
            center,
            y,
            xerr=np.asarray([[center - lower], [upper - center]]),
            fmt="D",
            color=palette[model],
            ecolor=palette[model],
            markersize=3.6,
            elinewidth=1.2,
            capsize=2.2,
            capthick=0.8,
            zorder=3,
        )
    ax_a.set_yticks(y_positions, [DISPLAY_CONFIG[value] for value in config_order_values])
    ax_a.set_xlabel("Composite biological score")
    ax_a.set_title("Matched development-validation performance", loc="left", pad=5)
    ax_a.grid(axis="x", color="#E6E6E6", linewidth=0.5, zorder=0)
    ax_a.tick_params(axis="y", length=0)
    handles = [
        Line2D([0], [0], marker="D", color="none", markerfacecolor="#333333", markeredgecolor="#333333", markersize=3.5, label="Fold-primary estimate (95% CI)"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor="#777777", markeredgecolor="#777777", markersize=3.5, label="Fold mean (3 seeds)"),
    ]
    ax_a.legend(handles=handles, loc="lower right", handletextpad=0.4, borderaxespad=0.2)

    primary_contrasts = [
        "B_paper_vs_A_paper",
        "B_paper_vs_C_paper",
        "B_paper_vs_B_log1p_mse",
    ]
    primary_metrics = [
        "primary_biological_score",
        "gene_exon_coverage_pearson_log1p",
        "per_track_pearson_128bp_log1p",
    ]
    entries = [(contrast, metric) for contrast in primary_contrasts for metric in primary_metrics]
    y_b = np.arange(len(entries))[::-1]
    labels_b = []
    for y, (contrast, metric) in zip(y_b, entries):
        estimate = next(
            row
            for row in rows
            if row["panel"] == "b"
            and row["data_role"] == "estimate"
            and row["series"] == contrast
            and row["metric_id"] == metric
        )
        fold_rows = [
            row
            for row in rows
            if row["panel"] == "b"
            and row["data_role"] == "fold_mean_effect"
            and row["series"] == contrast
            and row["metric_id"] == metric
        ]
        ax_b.scatter(
            [as_float(row, "value") for row in fold_rows],
            [y] * len(fold_rows),
            s=7,
            color=metric_colors[metric],
            alpha=0.42,
            linewidths=0,
            zorder=2,
        )
        center = as_float(estimate, "value")
        lower = as_float(estimate, "ci_lower")
        upper = as_float(estimate, "ci_upper")
        ax_b.errorbar(
            center,
            y,
            xerr=np.asarray([[center - lower], [upper - center]]),
            fmt="o",
            color=metric_colors[metric],
            ecolor=metric_colors[metric],
            markersize=3.2,
            elinewidth=1.1,
            capsize=2,
            zorder=3,
        )
        labels_b.append(
            f"{contrast_labels[contrast]}\n  {metric_labels[metric]}"
            if metric == primary_metrics[0]
            else f"  {metric_labels[metric]}"
        )
    ax_b.axvline(0.0, color="#555555", linewidth=0.75, linestyle="--", zorder=1)
    ax_b.set_yticks(y_b, labels_b)
    ax_b.set_xlabel("Favorable paired effect")
    ax_b.set_title("Pre-specified paired contrasts", loc="left", pad=5)
    ax_b.grid(axis="x", color="#E6E6E6", linewidth=0.5, zorder=0)
    ax_b.tick_params(axis="y", length=0)
    ax_b.text(
        0.99,
        0.02,
        "Right favors candidate",
        transform=ax_b.transAxes,
        ha="right",
        va="bottom",
        fontsize=5.5,
        color="#555555",
    )

    configs_c = ["A/paper", "B/paper", "C/paper"]
    data_c = [
        np.asarray(
            [
                as_float(row, "value")
                for row in rows
                if row["panel"] == "c"
                and row["data_role"] == "track_mean"
                and row["series"] == config
            ],
            dtype=float,
        )
        for config in configs_c
    ]
    box = ax_c.boxplot(
        data_c,
        positions=np.arange(1, 4),
        widths=0.48,
        patch_artist=True,
        showfliers=False,
        whis=1.5,
        medianprops={"color": "white", "linewidth": 1.1},
        whiskerprops={"color": "#555555", "linewidth": 0.7},
        capprops={"color": "#555555", "linewidth": 0.7},
    )
    for patch, config in zip(box["boxes"], configs_c):
        patch.set_facecolor(palette[config.split("/")[0]])
        patch.set_edgecolor(palette[config.split("/")[0]])
        patch.set_alpha(0.83)
        patch.set_linewidth(0.8)
    for index, (config, values) in enumerate(zip(configs_c, data_c), start=1):
        # A deterministic low-discrepancy offset displays every real track.
        jitter = ((np.arange(len(values)) * 0.6180339887498949) % 1.0 - 0.5) * 0.38
        ax_c.scatter(
            index + jitter,
            values,
            s=4,
            color=palette[config.split("/")[0]],
            alpha=0.22,
            linewidths=0,
            rasterized=True,
            zorder=1,
        )
    ax_c.set_xticks(np.arange(1, 4), [DISPLAY_CONFIG[value] for value in configs_c])
    ax_c.set_ylabel("Per-track composite score")
    ax_c.set_title("Track-level distribution (descriptive)", loc="left", pad=5)
    ax_c.grid(axis="y", color="#E6E6E6", linewidth=0.5, zorder=0)
    ax_c.text(
        0.02,
        0.02,
        "Each point: one of 241 tracks",
        transform=ax_c.transAxes,
        fontsize=5.5,
        color="#555555",
        va="bottom",
    )

    summary_d = [
        row
        for row in rows
        if row["panel"] == "d" and row["data_role"] == "median_iqr"
    ]
    dimension_order = {name: index for index, name in enumerate(("study", "stage", "tissue", "treatment", "signal_strength"))}
    summary_d.sort(key=lambda row: (dimension_order[row["series"]], row["subgroup"]))
    y_d = np.arange(len(summary_d))[::-1]
    label_prefix = {
        "study": "Study",
        "stage": "Stage",
        "tissue": "Tissue",
        "treatment": "Treatment",
        "signal_strength": "Signal",
    }
    compact_categories = {
        ("treatment", "No reported treatment"): "None reported",
        ("treatment", "Reported treatment/selection"): "Reported",
    }
    dimension_colors = {
        "study": "#707070",
        "stage": "#3C6E9E",
        "tissue": "#8E6C3A",
        "treatment": "#9C4F64",
        "signal_strength": "#007F73",
    }
    labels_d = []
    for y, row in zip(y_d, summary_d):
        center = as_float(row, "value")
        lower = as_float(row, "ci_lower")
        upper = as_float(row, "ci_upper")
        color = dimension_colors[row["series"]]
        ax_d.errorbar(
            center,
            y,
            xerr=np.asarray([[center - lower], [upper - center]]),
            fmt="o",
            color=color,
            ecolor=color,
            markersize=3.2,
            elinewidth=1.4,
            capsize=0,
        )
        compact_category = compact_categories.get(
            (row["series"], row["subgroup"]), row["subgroup"]
        )
        labels_d.append(
            f"{label_prefix[row['series']]} {compact_category} "
            f"({int(float(row['n']))})"
        )
    ax_d.axvline(0.0, color="#555555", linewidth=0.75, linestyle="--")
    ax_d.set_yticks(y_d, labels_d)
    ax_d.set_xlabel("B/paper - A/paper per-track score")
    ax_d.set_title("Metadata-defined strata (median and IQR)", loc="left", pad=5)
    ax_d.grid(axis="x", color="#E6E6E6", linewidth=0.5, zorder=0)
    ax_d.tick_params(axis="y", length=0)

    for label, axis in zip("abcd", axes):
        axis.text(
            -0.16,
            1.06,
            label,
            transform=axis.transAxes,
            fontsize=8,
            fontweight="bold",
            ha="left",
            va="bottom",
        )

    prefix = output_dir / "Figure_1_p6b_competitive_comparison"
    fig.savefig(prefix.with_suffix(".pdf"))
    fig.savefig(prefix.with_suffix(".svg"))
    fig.savefig(prefix.with_suffix(".png"), dpi=300)
    fig.savefig(
        prefix.with_suffix(".tiff"),
        dpi=600,
        pil_kwargs={"compression": "tiff_lzw"},
    )
    plt.close(fig)

    try:
        import PIL

        pillow_version = PIL.__version__
    except ImportError:
        pillow_version = "not_importable"
    figure_files = {
        suffix: {
            "path": prefix.with_suffix(suffix).relative_to(REPO_ROOT).as_posix(),
            "sha256": sha256_file(prefix.with_suffix(suffix)),
            "bytes": prefix.with_suffix(suffix).stat().st_size,
        }
        for suffix in (".pdf", ".svg", ".png", ".tiff")
    }
    return {
        "figure_python": platform.python_version(),
        "figure_numpy": np.__version__,
        "figure_matplotlib": matplotlib.__version__,
        "figure_pillow": pillow_version,
        "figure_executable": sys.executable,
        "figure_width_mm": 183.0,
        "figure_height_mm": 164.0,
        "figure_files": figure_files,
    }


def update_manifest_after_render(output_dir: Path, figure_audit: Mapping[str, Any]) -> None:
    manifest_path = output_dir / "analysis_manifest.json"
    manifest = json.loads(manifest_path.read_text())
    manifest["status"] = "complete"
    manifest["figure_rendered_at"] = utc_now()
    manifest["software"].update(figure_audit)
    json_dump(manifest_path, manifest)


def write_methods(
    output_dir: Path,
    spec: Mapping[str, Any],
    counts: Mapping[str, int],
    signal_audit: Mapping[str, Any],
) -> None:
    bootstrap = spec["bootstrap"]
    text = f"""# Statistical analysis methods

This analysis used the 90 completed P6B development-validation runs comprising six model/objective configurations, five blocked development folds and three training seeds per fold. It did not read the locked final-test intervals, signals or report. The analysis plan in `p6b_submission_analysis_spec.json` was frozen before this submission-analysis execution but after the training and validation results existed; it is therefore a retrospective locked analysis plan, not a prospective or public preregistration.

The primary endpoint was the arithmetic mean of per-track gene-exon coverage Pearson correlation after log1p transformation and per-track 128-bp Pearson correlation after log1p transformation. Folds were the primary resampling units. Seeds were treated as algorithmic repeats within fold, not as independent biological observations. Point estimates first averaged the three seeds within each fold and then averaged the five fold means with equal weight.

Confidence intervals used a hierarchical paired percentile bootstrap ({bootstrap['iterations']:,} iterations; NumPy PCG64 seed {bootstrap['seed']}). Each iteration sampled five folds with replacement. For every sampled fold occurrence, three seeds were independently sampled with replacement, averaged within that fold occurrence, and the five sampled fold means were then averaged. Paired comparisons were formed as candidate minus comparator at identical fold and seed before resampling. For lower-is-better endpoints, favorable effects reverse the sign so positive values consistently favor the candidate. Intervals are two-sided {int(float(bootstrap['confidence_level']) * 100)}% percentile intervals. No p values were computed. Three contrasts were designated primary; all other contrasts and diagnostic endpoints are secondary.

All {counts['tracks']} tracks were retained. The long-form table contains every 241-entry per-track dictionary from every validation JSON plus the per-track composite score ({counts['track_long_rows']:,} rows). Track-level configuration summaries average the 15 fold-by-seed values for each track. Tracks are shown descriptively and are not treated as independent inferential replicates.

Strata were derived only from manifest fields or development-training signal means. Study retained the exact SRA study accession. Stage and tissue used the deterministic synonym rules in the analysis script, with unmatched values retained verbatim. Treatment was classified as reported treatment/selection when a non-administrative condition was present or genotype/strain explicitly contained `treated` or `RNAi`; otherwise it was classified as no reported treatment, except when all three source fields were absent (`Missing/unknown`). Signal strength used tertiles of `development_train_nonzero_mean`: Low <= {signal_audit['lower_boundary']:.12g}, Middle > {signal_audit['lower_boundary']:.12g} and <= {signal_audit['upper_boundary']:.12g}, and High > {signal_audit['upper_boundary']:.12g}. Missing and unknown values were retained. All stratum summaries are descriptive.

Analysis used Python {platform.python_version()} and NumPy {np.__version__}. Figure rendering uses Python/matplotlib only; exact render versions are recorded in `analysis_manifest.json` after export.
"""
    (output_dir / "METHODS.md").write_text(text, encoding="utf-8")


def write_legend(output_dir: Path) -> None:
    text = """# Figure legend

**Figure 1 | Matched development cross-validation supports Model B with the paper objective.** **a,** Composite biological score for each architecture/objective configuration. Small points are fold means after averaging three training seeds; diamonds are equal-weight means across five folds and bars are two-sided 95% hierarchical bootstrap confidence intervals. **b,** Pre-specified paired comparisons for the composite score and its gene-exon and 128-bp components. Effects were calculated at identical fold and seed; positive values favor the named candidate. Small points are fold means and bars are fold-primary hierarchical bootstrap 95% confidence intervals. **c,** Descriptive distribution of per-track composite scores for the three paper-objective models. Each point is one of 241 RNA-seq tracks averaged across five folds and three seeds; boxes show the median and interquartile range, whiskers extend to 1.5 times the interquartile range, and all tracks are retained. **d,** Descriptive B/paper minus A/paper per-track differences within manifest-defined strata. Points show medians and bars show interquartile ranges; `n` is the number of tracks, not an independent inferential sample size. For study, stage, tissue and treatment, the two most prevalent categories containing at least five tracks were selected using metadata prevalence alone; all signal-strength tertiles are shown. Folds (n = 5) are the primary resampling units and seeds (n = 3 per fold) are algorithmic repeats. No p values were calculated. Source data are provided in `Source_Data_Figure_1.csv` and `.tsv`.
"""
    (output_dir / "FIGURE_LEGEND.md").write_text(text, encoding="utf-8")


def write_readme(output_dir: Path, counts: Mapping[str, int]) -> None:
    text = f"""# P6B submission analysis bundle

This directory is generated from the locked 90-run P6B development-validation matrix. It contains {counts['run_level_rows']} complete run rows, {counts['paired_fold_seed_rows']:,} matched fold/seed contrast rows, {counts['track_long_rows']:,} long-form track rows spanning all {counts['tracks']} tracks, descriptive metadata strata, source data and a four-panel figure bundle. No locked final-test input is read.

## Reproduce

Run extraction and statistics in the AlphaGenome environment:

```bash
/home/zelinli6/miniconda3/envs/alphagenome/bin/python -m scripts.analyze_v2_p6b_submission --skip-figures
```

Render the existing source data with the Python/matplotlib environment:

```bash
/home/zelinli6/miniconda3/envs/CellUNetr/bin/python -m scripts.analyze_v2_p6b_submission --render-only
```

The analysis stage requires Python and NumPy only. Rendering requires Python, NumPy, matplotlib and Pillow; exact versions and executables are recorded in `analysis_manifest.json`.

## Files

- `analysis_manifest.json`: input guards, exact counts, software versions, signal tertile boundaries and headline estimates.
- `run_level_metrics.tsv` / `.csv`: all 90 run identities, guards, primary/diagnostic endpoints and scalar validation metrics.
- `configuration_bootstrap_estimates.tsv`: fold-primary estimates and 95% intervals for all six configurations and ten endpoints.
- `paired_fold_seed_contrasts.tsv`: candidate/comparator values and effects at identical fold and seed.
- `paired_contrast_bootstrap_estimates.tsv`: hierarchical bootstrap intervals for all planned contrasts and endpoints.
- `track_metrics_long.tsv`: normalized long form for every 241-entry metric dictionary in every run plus the derived per-track composite.
- `track_configuration_summaries.tsv`: per-configuration, per-track averages across 15 fold/seed runs.
- `track_metadata_strata.tsv`: raw manifest fields, deterministic strata and signal-strength assignments.
- `stratified_track_summaries.tsv`: descriptive configuration and paired track summaries by study, stage, tissue, treatment and signal strength.
- `Source_Data_Figure_1.csv` / `.tsv`: plotted values and summaries for panels a-d.
- `Figure_1_p6b_competitive_comparison.*`: PDF, editable-text SVG, 600-dpi TIFF and 300-dpi PNG.
- `METHODS.md` and `FIGURE_LEGEND.md`: manuscript-ready statistical methods and legend.

Track-level distributions and strata are descriptive. The inferential unit for run-level intervals is the development-validation fold; seeds are nested algorithmic repeats.
"""
    (output_dir / "README.md").write_text(text, encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--spec", type=Path, default=DEFAULT_SPEC)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument(
        "--skip-figures",
        action="store_true",
        help="write all analysis/source-data outputs without importing matplotlib",
    )
    mode.add_argument(
        "--render-only",
        action="store_true",
        help="render from an existing Source_Data_Figure_1.tsv",
    )
    return parser.parse_args()


def main() -> None:
    global DEFAULT_SPEC
    args = parse_args()
    output_dir = args.output_dir.resolve()
    if args.render_only:
        if not (output_dir / "Source_Data_Figure_1.tsv").is_file():
            raise SystemExit(f"missing source data: {output_dir / 'Source_Data_Figure_1.tsv'}")
        audit = render_figure(output_dir)
        update_manifest_after_render(output_dir, audit)
        print(json.dumps({"status": "rendered", **audit}, indent=2, sort_keys=True))
        return

    spec_path = args.spec.resolve()
    if spec_path != DEFAULT_SPEC.resolve():
        DEFAULT_SPEC = spec_path
    spec = json.loads(spec_path.read_text())
    manifest = analyse(spec, output_dir)
    if args.skip_figures:
        print(json.dumps(manifest, indent=2, sort_keys=True))
        return
    audit = render_figure(output_dir)
    update_manifest_after_render(output_dir, audit)
    print(json.dumps({"analysis": manifest, "figure": audit}, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
