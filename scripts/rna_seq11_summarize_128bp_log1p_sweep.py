#!/usr/bin/env python3
"""Summarize validation-only 128 bp log1p-MSE sweep outputs.

This parser is intentionally read-only. It consumes the launcher summary TSV and
writes small derived summaries under ignored run/log locations. It never reads
or evaluates the held-out test split.
"""

from __future__ import annotations

import argparse
import csv
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_SELECTION_METRIC = "full-mse"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--date-tag", default="20260522")
    parser.add_argument(
        "--summary",
        type=Path,
        default=None,
        help="Launcher TSV, default runs/rna_seq11_128bp_log1p_sweep_<date>/sweep_results.tsv.",
    )
    parser.add_argument(
        "--selection-metric",
        choices=["full-mse", "full-pearson", "common128-mse"],
        default=DEFAULT_SELECTION_METRIC,
        help="Validation metric for ranking rows; default is full-resolution MSE.",
    )
    parser.add_argument(
        "--ranked-output",
        type=Path,
        default=None,
        help="Optional ranked TSV path under an ignored output directory.",
    )
    parser.add_argument(
        "--markdown-output",
        type=Path,
        default=None,
        help="Optional concise Markdown summary path under docs or ignored outputs.",
    )
    return parser.parse_args()


def default_summary_path(date_tag: str) -> Path:
    return ROOT / "runs" / f"rna_seq11_128bp_log1p_sweep_{date_tag}" / "sweep_results.tsv"


def default_ranked_path(date_tag: str, metric: str) -> Path:
    return (
        ROOT
        / "runs"
        / f"rna_seq11_128bp_log1p_sweep_{date_tag}"
        / f"sweep_ranked_{metric}.tsv"
    )


def read_rows(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        raise FileNotFoundError(f"Missing sweep summary: {path}")
    with path.open("r", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def parse_float(value: str | None) -> float | None:
    if value in {None, ""}:
        return None
    try:
        return float(value)
    except ValueError:
        return None


def metric_value(row: dict[str, str], metric: str) -> float | None:
    if metric == "full-mse":
        return parse_float(row.get("full_mse"))
    if metric == "full-pearson":
        value = parse_float(row.get("full_pearson"))
        if value is None:
            return None
        return -value
    if metric == "common128-mse":
        return parse_float(row.get("common128_mse"))
    raise ValueError(f"Unsupported metric: {metric}")


def completed_ranked_rows(rows: list[dict[str, str]], metric: str) -> list[dict[str, str]]:
    candidates = [
        row
        for row in rows
        if row.get("status") == "completed" and metric_value(row, metric) is not None
    ]
    return sorted(candidates, key=lambda row: float(metric_value(row, metric) or 0.0))


def write_ranked(path: Path, rows: list[dict[str, str]], metric: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "rank",
        "selection_metric",
        "stage",
        "run_id",
        "status",
        "arch",
        "target_space",
        "loss",
        "lr",
        "weight_decay",
        "seed",
        "max_steps",
        "train_best_step",
        "train_best_loss",
        "full_mse",
        "full_mae",
        "full_pearson",
        "common128_mse",
        "common128_mae",
        "common128_pearson",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        for rank, row in enumerate(rows, start=1):
            writer.writerow(
                {
                    "rank": rank,
                    "selection_metric": metric,
                    **{field: row.get(field, "") for field in fields if field not in {"rank", "selection_metric"}},
                }
            )


def write_markdown(
    path: Path,
    *,
    summary_path: Path,
    ranked_path: Path,
    rows: list[dict[str, str]],
    ranked_rows: list[dict[str, str]],
    metric: str,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    status_counts = Counter(row.get("status", "unknown") for row in rows)
    stage_counts = Counter(row.get("stage", "unknown") for row in rows)
    lines = [
        "# 128 bp Log1p-MSE Sweep Summary",
        "",
        "- Claim status: unverified until commands, logs, and outputs are reviewed.",
        "- Split policy: validation-only summary; held-out test split is not read.",
        f"- Source TSV: `{summary_path.relative_to(ROOT)}`",
        f"- Ranked TSV: `{ranked_path.relative_to(ROOT)}`",
        f"- Selection metric: `{metric}`",
        "",
        "## Counts",
        "",
        f"- Total rows: `{len(rows)}`",
        f"- Rankable completed rows: `{len(ranked_rows)}`",
    ]
    for status, count in sorted(status_counts.items()):
        lines.append(f"- Status `{status}`: `{count}`")
    lines.extend(["", "## Stages", ""])
    for stage, count in sorted(stage_counts.items()):
        lines.append(f"- `{stage}`: `{count}`")
    lines.extend(["", "## Top Validation Candidates", ""])
    if not ranked_rows:
        lines.append("No completed rankable rows were found.")
    else:
        lines.append("| rank | run_id | full_mse | full_pearson | common128_mse |")
        lines.append("|---:|---|---:|---:|---:|")
        for rank, row in enumerate(ranked_rows[:10], start=1):
            lines.append(
                "| "
                f"{rank} | `{row.get('run_id', '')}` | "
                f"{row.get('full_mse', '')} | "
                f"{row.get('full_pearson', '')} | "
                f"{row.get('common128_mse', '')} |"
            )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    args = parse_args()
    summary_path = args.summary or default_summary_path(args.date_tag)
    if not summary_path.is_absolute():
        summary_path = ROOT / summary_path
    ranked_path = args.ranked_output or default_ranked_path(
        args.date_tag,
        args.selection_metric,
    )
    if not ranked_path.is_absolute():
        ranked_path = ROOT / ranked_path

    rows = read_rows(summary_path)
    ranked_rows = completed_ranked_rows(rows, args.selection_metric)
    write_ranked(ranked_path, ranked_rows, args.selection_metric)

    if args.markdown_output is not None:
        markdown_path = args.markdown_output
        if not markdown_path.is_absolute():
            markdown_path = ROOT / markdown_path
        write_markdown(
            markdown_path,
            summary_path=summary_path,
            ranked_path=ranked_path,
            rows=rows,
            ranked_rows=ranked_rows,
            metric=args.selection_metric,
        )

    print(f"summary\t{summary_path}")
    print(f"selection_metric\t{args.selection_metric}")
    print(f"rows\t{len(rows)}")
    print(f"rankable_completed_rows\t{len(ranked_rows)}")
    print(f"ranked_output\t{ranked_path}")
    if ranked_rows:
        best = ranked_rows[0]
        print(f"best_run_id\t{best.get('run_id', '')}")
        print(f"best_full_mse\t{best.get('full_mse', '')}")
        print(f"best_full_pearson\t{best.get('full_pearson', '')}")
        print(f"best_common128_mse\t{best.get('common128_mse', '')}")


if __name__ == "__main__":
    main()
