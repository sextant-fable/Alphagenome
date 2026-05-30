#!/usr/bin/env python3
"""Merge RNA-seq11 representation diagnostics into utility leaderboards."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import numpy as np


FLOAT_FORMAT = ".8g"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--inventory-tsv", required=True)
    parser.add_argument(
        "--diagnostics-dirs",
        nargs="*",
        default=[],
        help="One or more extended diagnostics output directories or shard directories.",
    )
    parser.add_argument(
        "--probe-dirs",
        nargs="*",
        default=[],
        help="One or more gene-profile probe output directories or shard directories.",
    )
    parser.add_argument(
        "--output-dir",
        default="runs/rna_seq11_representation_utility_summary_20260530",
    )
    parser.add_argument(
        "--mse-sanity-delta",
        type=float,
        default=0.02,
        help="Composite-score penalty starts when full MSE exceeds best full MSE by this amount.",
    )
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_many(paths: list[Path], filename: str) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for path in paths:
        rows.extend(read_tsv(path / filename))
    return rows


def as_float(value: str | None) -> float:
    if value is None or value == "":
        return math.nan
    try:
        return float(value)
    except ValueError:
        return math.nan


def format_float(value: float) -> str:
    if math.isnan(value):
        return "nan"
    return format(float(value), FLOAT_FORMAT)


def one(rows: list[dict[str, str]], run_id: str, **conditions: str) -> dict[str, str]:
    for row in rows:
        if row.get("run_id") != run_id:
            continue
        if all(row.get(key) == value for key, value in conditions.items()):
            return row
    return {}


def average(values: list[float]) -> float:
    finite = [value for value in values if not math.isnan(value)]
    if not finite:
        return math.nan
    return float(np.mean(finite))


def mean_rows(
    rows: list[dict[str, str]],
    run_id: str,
    value_key: str,
    *,
    track_indices: set[int] | None = None,
    **conditions: str,
) -> float:
    values: list[float] = []
    for row in rows:
        if row.get("run_id") != run_id:
            continue
        if not all(row.get(key) == value for key, value in conditions.items()):
            continue
        if track_indices is not None:
            track_index = row.get("track_index")
            if track_index in (None, ""):
                continue
            try:
                if int(track_index) not in track_indices:
                    continue
            except ValueError:
                continue
        values.append(as_float(row.get(value_key)))
    return average(values)


def rank_scores(values: dict[str, float], *, higher_is_better: bool) -> dict[str, float]:
    finite_items = [(key, value) for key, value in values.items() if not math.isnan(value)]
    if not finite_items:
        return {key: math.nan for key in values}
    finite_items.sort(key=lambda item: item[1], reverse=higher_is_better)
    n = len(finite_items)
    scores: dict[str, float] = {}
    index = 0
    while index < n:
        value = finite_items[index][1]
        end = index + 1
        while end < n and finite_items[end][1] == value:
            end += 1
        mean_rank = (index + 1 + end) / 2.0
        score = 1.0 if n == 1 else 1.0 - ((mean_rank - 1.0) / (n - 1.0))
        for key, _value in finite_items[index:end]:
            scores[key] = score
        index = end
    for key in values:
        scores.setdefault(key, math.nan)
    return scores


def averaged_rank_score(
    metrics: dict[str, dict[str, float]],
    specs: list[tuple[str, bool]],
) -> dict[str, float]:
    per_metric_scores = [
        rank_scores(metrics[name], higher_is_better=higher_is_better)
        for name, higher_is_better in specs
        if name in metrics
    ]
    run_ids = sorted({run_id for scores in per_metric_scores for run_id in scores})
    output: dict[str, float] = {}
    for run_id in run_ids:
        vals = [scores.get(run_id, math.nan) for scores in per_metric_scores]
        output[run_id] = average(vals)
    return output


def main() -> None:
    args = parse_args()
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    diagnostics_dirs = [Path(value) for value in args.diagnostics_dirs]
    probe_dirs = [Path(value) for value in args.probe_dirs]

    inventory = read_tsv(Path(args.inventory_tsv))
    run_ids = [row["run_id"] for row in inventory]
    inventory_by_run = {row["run_id"]: row for row in inventory}

    per_track = read_many(diagnostics_dirs, "per_track_metrics.tsv")
    signal = read_many(diagnostics_dirs, "signal_strata_metrics.tsv")
    region = read_many(diagnostics_dirs, "region_metrics.tsv")
    gene = read_many(diagnostics_dirs, "gene_metrics.tsv")
    resolution = read_many(diagnostics_dirs, "resolution_metrics.tsv")
    localization = read_many(diagnostics_dirs, "high_signal_localization_metrics.tsv")
    probe = read_many(probe_dirs, "gene_profile_probe_metrics.tsv")
    probe_by_run = {row["run_id"]: row for row in probe}

    metrics_for_scores: dict[str, dict[str, float]] = {
        "full_mse": {},
        "pearson": {},
        "gene_mse": {},
        "gene_spearman": {},
        "high_gt3_mse": {},
        "top5_mse": {},
        "top1_mse": {},
        "top5_overlap": {},
        "top1_overlap": {},
        "intestine_t134_gene_mse": {},
        "intestine_t134_gene_spearman": {},
        "intestine_t134_high_gt3_mse": {},
        "intestine_t134_top5_mse": {},
        "intestine_t134_top1_mse": {},
        "intestine_t134_top5_overlap": {},
        "intestine_t134_top1_overlap": {},
        "exon_mse": {},
        "gene_body_mse": {},
        "pool128_mse": {},
        "pool1024_mse": {},
        "local_gradient_pearson": {},
        "probe_top05_recall": {},
        "probe_intestine_t134_top05_recall": {},
        "probe_track_corr_spearman": {},
        "probe_family_balanced_accuracy": {},
        "probe_time_balanced_accuracy": {},
        "probe_track_balanced_accuracy": {},
    }

    summary_rows: list[dict[str, str]] = []
    intestine_t134_tracks = {1, 3, 4}
    for run_id in run_ids:
        inv = inventory_by_run[run_id]
        p = one(per_track, run_id, metric_scope="overall_full_resolution", track_index="")
        full_mse = as_float(p.get("mse")) if p else as_float(inv.get("full_mse"))
        pearson = as_float(p.get("pearson")) if p else as_float(inv.get("pearson"))
        high_gt3 = one(signal, run_id, metric_scope="signal_stratum", stratum="high_gt3", track_index="")
        top5 = one(signal, run_id, metric_scope="signal_stratum", stratum="top5_true", track_index="")
        top1 = one(signal, run_id, metric_scope="signal_stratum", stratum="top1_true", track_index="")
        exon = one(region, run_id, metric_scope="valid_core_unique_region", region="exon", track_index="")
        gene_body = one(region, run_id, metric_scope="valid_core_unique_region", region="gene_body", track_index="")
        intergenic = one(region, run_id, metric_scope="valid_core_unique_region", region="intergenic", track_index="")
        gene_overall = one(gene, run_id, metric_scope="gene_body_mean", track_index="")
        pool128 = one(resolution, run_id, metric_scope="pooled_value", pool_bp="128", track_index="")
        pool1024 = one(resolution, run_id, metric_scope="pooled_value", pool_bp="1024", track_index="")
        gradient1 = one(resolution, run_id, metric_scope="local_gradient", pool_bp="1", track_index="")
        top5_overlap = average(
            [
                as_float(row.get("precision"))
                for row in localization
                if row.get("run_id") == run_id and row.get("top_label") == "top5"
            ]
        )
        top1_overlap = average(
            [
                as_float(row.get("precision"))
                for row in localization
                if row.get("run_id") == run_id and row.get("top_label") == "top1"
            ]
        )
        intestine_t134_gene_mse = mean_rows(
            gene,
            run_id,
            "mse",
            metric_scope="gene_body_mean",
            track_indices=intestine_t134_tracks,
        )
        intestine_t134_gene_pearson = mean_rows(
            gene,
            run_id,
            "pearson",
            metric_scope="gene_body_mean",
            track_indices=intestine_t134_tracks,
        )
        intestine_t134_gene_spearman = mean_rows(
            gene,
            run_id,
            "spearman",
            metric_scope="gene_body_mean",
            track_indices=intestine_t134_tracks,
        )
        intestine_t134_high_gt3_mse = mean_rows(
            signal,
            run_id,
            "mse",
            metric_scope="signal_stratum",
            stratum="high_gt3",
            track_indices=intestine_t134_tracks,
        )
        intestine_t134_top5_mse = mean_rows(
            signal,
            run_id,
            "mse",
            metric_scope="signal_stratum",
            stratum="top5_true",
            track_indices=intestine_t134_tracks,
        )
        intestine_t134_top1_mse = mean_rows(
            signal,
            run_id,
            "mse",
            metric_scope="signal_stratum",
            stratum="top1_true",
            track_indices=intestine_t134_tracks,
        )
        intestine_t134_top5_overlap = mean_rows(
            localization,
            run_id,
            "precision",
            metric_scope="128bp_pooled_high_signal_overlap",
            top_label="top5",
            track_indices=intestine_t134_tracks,
        )
        intestine_t134_top1_overlap = mean_rows(
            localization,
            run_id,
            "precision",
            metric_scope="128bp_pooled_high_signal_overlap",
            top_label="top1",
            track_indices=intestine_t134_tracks,
        )
        probe_row = probe_by_run.get(run_id, {})
        values = {
            "full_mse": full_mse,
            "pearson": pearson,
            "gene_mse": as_float(gene_overall.get("mse")),
            "gene_pearson": as_float(gene_overall.get("pearson")),
            "gene_spearman": as_float(gene_overall.get("spearman")),
            "high_gt3_mse": as_float(high_gt3.get("mse")),
            "top5_mse": as_float(top5.get("mse")),
            "top1_mse": as_float(top1.get("mse")),
            "top5_overlap": top5_overlap,
            "top1_overlap": top1_overlap,
            "intestine_t134_gene_mse": intestine_t134_gene_mse,
            "intestine_t134_gene_pearson": intestine_t134_gene_pearson,
            "intestine_t134_gene_spearman": intestine_t134_gene_spearman,
            "intestine_t134_high_gt3_mse": intestine_t134_high_gt3_mse,
            "intestine_t134_top5_mse": intestine_t134_top5_mse,
            "intestine_t134_top1_mse": intestine_t134_top1_mse,
            "intestine_t134_top5_overlap": intestine_t134_top5_overlap,
            "intestine_t134_top1_overlap": intestine_t134_top1_overlap,
            "exon_mse": as_float(exon.get("mse")),
            "gene_body_mse": as_float(gene_body.get("mse")),
            "intergenic_mse": as_float(intergenic.get("mse")),
            "pool128_mse": as_float(pool128.get("mse")),
            "pool1024_mse": as_float(pool1024.get("mse")),
            "local_gradient_pearson": as_float(gradient1.get("pearson")),
            "probe_top05_recall": as_float(probe_row.get("valid_top05_recall")),
            "probe_intestine_t134_top05_recall": as_float(
                probe_row.get("intestine_t134_top05_recall")
            ),
            "probe_track_corr_spearman": as_float(probe_row.get("track_corr_spearman")),
            "probe_family_balanced_accuracy": as_float(
                probe_row.get("family_probe_balanced_accuracy")
            ),
            "probe_time_balanced_accuracy": as_float(
                probe_row.get("time_probe_balanced_accuracy")
            ),
            "probe_track_balanced_accuracy": as_float(
                probe_row.get("track_probe_balanced_accuracy")
            ),
        }
        for key, value in values.items():
            if key in metrics_for_scores:
                metrics_for_scores[key][run_id] = value
        summary_rows.append(
            {
                "rank": inv["rank"],
                "run_id": run_id,
                "checkpoint": inv["checkpoint"],
                "best_step": inv.get("best_step", ""),
                **{key: format_float(value) for key, value in values.items()},
            }
        )

    gene_score = averaged_rank_score(
        metrics_for_scores,
        [
            ("gene_mse", False),
            ("gene_spearman", True),
            ("intestine_t134_gene_mse", False),
            ("intestine_t134_gene_spearman", True),
            ("probe_top05_recall", True),
            ("probe_track_corr_spearman", True),
        ],
    )
    high_signal_score = averaged_rank_score(
        metrics_for_scores,
        [
            ("high_gt3_mse", False),
            ("top5_mse", False),
            ("top1_mse", False),
            ("top5_overlap", True),
            ("top1_overlap", True),
            ("intestine_t134_high_gt3_mse", False),
            ("intestine_t134_top5_mse", False),
            ("intestine_t134_top1_mse", False),
            ("intestine_t134_top5_overlap", True),
            ("intestine_t134_top1_overlap", True),
            ("probe_intestine_t134_top05_recall", True),
        ],
    )
    region_score = averaged_rank_score(
        metrics_for_scores,
        [("exon_mse", False), ("gene_body_mse", False)],
    )
    probe_score = averaged_rank_score(
        metrics_for_scores,
        [
            ("probe_family_balanced_accuracy", True),
            ("probe_time_balanced_accuracy", True),
            ("probe_track_balanced_accuracy", True),
            ("probe_track_corr_spearman", True),
        ],
    )

    best_full_mse = min(
        (as_float(row["full_mse"]) for row in summary_rows if row.get("full_mse")),
        default=math.nan,
    )
    for row in summary_rows:
        run_id = row["run_id"]
        full_mse = as_float(row.get("full_mse"))
        penalty = 0.0
        if not math.isnan(best_full_mse) and not math.isnan(full_mse):
            excess = max(0.0, full_mse - best_full_mse - args.mse_sanity_delta)
            penalty = min(0.5, excess / max(args.mse_sanity_delta, 1e-12))
        composite = average(
            [
                gene_score.get(run_id, math.nan),
                high_signal_score.get(run_id, math.nan),
                region_score.get(run_id, math.nan),
                probe_score.get(run_id, math.nan),
            ]
        )
        if not math.isnan(composite):
            composite = max(0.0, composite - penalty)
        row["gene_utility_score"] = format_float(gene_score.get(run_id, math.nan))
        row["high_signal_utility_score"] = format_float(
            high_signal_score.get(run_id, math.nan)
        )
        row["region_utility_score"] = format_float(region_score.get(run_id, math.nan))
        row["probe_utility_score"] = format_float(probe_score.get(run_id, math.nan))
        row["mse_sanity_penalty"] = format_float(penalty)
        row["representation_utility_score"] = format_float(composite)

    fieldnames = list(summary_rows[0].keys()) if summary_rows else []
    with (output_dir / "all_candidate_utility_summary.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(summary_rows)

    leaderboard_specs = [
        ("full_mse_top10.tsv", "full_mse", False),
        ("gene_utility_top10.tsv", "gene_utility_score", True),
        ("high_signal_utility_top10.tsv", "high_signal_utility_score", True),
        ("region_utility_top10.tsv", "region_utility_score", True),
        ("probe_utility_top10.tsv", "probe_utility_score", True),
        ("representation_utility_top20.tsv", "representation_utility_score", True),
    ]
    for filename, key, higher in leaderboard_specs:
        sorted_rows = sorted(
            summary_rows,
            key=lambda row: (
                math.inf
                if math.isnan(as_float(row.get(key)))
                else (-as_float(row[key]) if higher else as_float(row[key])),
                as_float(row.get("full_mse")),
                row["run_id"],
            ),
        )
        limit = 20 if "top20" in filename else 10
        with (output_dir / filename).open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()
            writer.writerows(sorted_rows[:limit])

    with (output_dir / "leaderboards.md").open("w") as handle:
        handle.write("# RNA-seq11 Representation Utility Leaderboards\n\n")
        handle.write(f"Candidates: {len(summary_rows)}\n\n")
        for filename, key, higher in leaderboard_specs:
            handle.write(f"## {filename.removesuffix('.tsv')}\n\n")
            rows = read_tsv(output_dir / filename)
            handle.write("| rank | run_id | score | full_mse | pearson |\n")
            handle.write("|---:|---|---:|---:|---:|\n")
            for index, row in enumerate(rows[:10], start=1):
                handle.write(
                    f"| {index} | `{row['run_id']}` | {row.get(key, '')} | "
                    f"{row.get('full_mse', '')} | {row.get('pearson', '')} |\n"
                )
            handle.write("\n")
    print(f"output_dir\t{output_dir}")


if __name__ == "__main__":
    main()
