#!/usr/bin/env python3
"""Run validation-only RNA-seq11 adapter adaptation phases on HY-GPU.

The launcher never reads the held-out test split. It writes model outputs under
ignored runs/ paths and logs under ignored logs/ paths.
"""

from __future__ import annotations

import argparse
import csv
import os
import queue
import subprocess
import sys
import threading
import time
from dataclasses import dataclass, replace
from datetime import datetime
from pathlib import Path
from typing import Iterable


ROOT = Path(__file__).resolve().parents[1]
TRAIN_DIR = "alphagenome_custom/datasets/rna_seq_npz_train"
VALID_DIR = "alphagenome_custom/datasets/rna_seq_npz_valid"
WEIGHTS = "weights/alphagenome_pytorch/model_all_folds.safetensors"
PHASES = ("phase1", "phase2-screen", "phase2-promote", "phase3", "phase4")


@dataclass(frozen=True)
class RunConfig:
    phase: str
    stage: str
    arch: str
    seed: int
    max_steps: int
    lr: float = 1e-3
    weight_decay: float = 0.0
    loss: str = "smooth-l1"
    smooth_l1_beta: float = 1.0
    hybrid_alpha: float = 0.5
    target_space: str = "full-log1p"
    hidden_channels: int = 256
    residual_scale_init: float = 1.0
    linear_dilation: int = 2
    lr_schedule: str = "constant"
    step_decay_steps: str = ""
    step_decay_learning_rates: str = ""
    init_checkpoint: str = ""
    eval_every: int = 250
    selection_metric: str = "full-mse"
    rank_label: str = ""
    run_id_override: str = ""

    @property
    def run_id(self) -> str:
        if self.run_id_override:
            return self.run_id_override
        pieces = [
            "rna_seq11",
            self.phase.replace("-", ""),
            self.stage,
            self.arch.replace("-", ""),
            f"h{self.hidden_channels}",
        ]
        if self.arch == "residual-conv3":
            pieces.append(f"rs{self.residual_scale_init:g}")
        if self.arch == "dilated-conv3":
            pieces.append(f"d{self.linear_dilation}")
        if self.rank_label:
            pieces.append(self.rank_label)
        pieces.extend(
            [
                self.target_space.replace("-", ""),
                self.loss.replace("-", ""),
                f"b{self.smooth_l1_beta:g}",
                f"lr{self.lr:g}",
                f"seed{self.seed}",
                f"{self.max_steps}steps",
            ]
        )
        if self.lr_schedule != "constant":
            pieces.append(self.lr_schedule)
            if self.step_decay_steps:
                pieces.append(f"s{self.step_decay_steps.replace(',', '-')}")
        if self.init_checkpoint:
            pieces.append("warmstart")
        return "_".join(pieces)


@dataclass
class RunResult:
    config: RunConfig
    status: str
    returncode: int
    gpu: str
    best_step: int | None = None
    best_loss: float | None = None
    full_mse: float | None = None
    full_mae: float | None = None
    pearson: float | None = None
    common128_mse: float | None = None
    common128_mae: float | None = None
    selection_value: float | None = None


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--date-tag", default="20260524")
    parser.add_argument("--phase", choices=[*PHASES, "all"], default="phase1")
    parser.add_argument("--gpus", default="0,1,2,3")
    parser.add_argument("--selection-metric", choices=["full-mse", "loss", "pearson", "common128-mse"], default="full-mse")
    parser.add_argument("--summary", default=None)
    parser.add_argument("--dry-run", action="store_true")
    return parser.parse_args()


def now() -> str:
    return datetime.now().isoformat(timespec="seconds")


def rel(path: Path) -> str:
    return str(path.relative_to(ROOT))


def run_dir(config: RunConfig) -> Path:
    return ROOT / "runs" / config.run_id


def log_path(config: RunConfig, date_tag: str) -> Path:
    return ROOT / "logs" / f"rna_seq11_adaptation_phase_sweep_{date_tag}" / f"{config.run_id}.log"


def summary_path(date_tag: str, override: str | None) -> Path:
    if override:
        path = Path(override)
        return path if path.is_absolute() else ROOT / path
    return ROOT / "runs" / f"rna_seq11_adaptation_phase_sweep_{date_tag}" / "sweep_results.tsv"


def git_head() -> str:
    completed = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=ROOT,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        check=False,
    )
    return completed.stdout.strip() if completed.returncode == 0 else "unknown"


def train_command(config: RunConfig) -> list[str]:
    cmd = [
        sys.executable,
        "-u",
        "scripts/alphagenome_rna_seq11_finetune.py",
        "--train-dataset-dir",
        TRAIN_DIR,
        "--valid-dataset-dir",
        VALID_DIR,
        "--weights",
        WEIGHTS,
        "--output-dir",
        rel(run_dir(config)),
        "--head-type",
        "linear",
        "--embedding-resolution",
        "128",
        "--linear-head-architecture",
        config.arch,
        "--linear-hidden-channels",
        str(config.hidden_channels),
        "--residual-scale-init",
        str(config.residual_scale_init),
        "--linear-dilation",
        str(config.linear_dilation),
        "--linear-loss-type",
        config.loss,
        "--smooth-l1-beta",
        str(config.smooth_l1_beta),
        "--hybrid-loss-alpha",
        str(config.hybrid_alpha),
        "--linear-target-space",
        config.target_space,
        "--target-transform",
        "log1p",
        "--selection-metric",
        config.selection_metric,
        "--lr-schedule",
        config.lr_schedule,
        "--batch-size",
        "1",
        "--grad-accum-steps",
        "1",
        "--num-workers",
        "0",
        "--max-steps",
        str(config.max_steps),
        "--eval-every",
        str(config.eval_every),
        "--learning-rate",
        str(config.lr),
        "--weight-decay",
        str(config.weight_decay),
        "--seed",
        str(config.seed),
        "--device",
        "auto",
    ]
    if config.lr_schedule == "step":
        cmd.extend(
            [
                "--step-decay-steps",
                config.step_decay_steps,
                "--step-decay-learning-rates",
                config.step_decay_learning_rates,
            ]
        )
    if config.init_checkpoint:
        cmd.extend(["--init-adapter-checkpoint", config.init_checkpoint])
    return cmd


def append_log(path: Path, message: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8") as handle:
        handle.write(message)
        if not message.endswith("\n"):
            handle.write("\n")


def run_command(cmd: list[str], *, gpu: str, log: Path, dry_run: bool) -> int:
    if dry_run:
        print("dry_run_command\t" + " ".join(cmd), flush=True)
        return 0
    append_log(log, f"\n[{now()}] START gpu={gpu}")
    append_log(log, " ".join(cmd))
    env = os.environ.copy()
    env["CUDA_VISIBLE_DEVICES"] = gpu
    env["PYTHONUNBUFFERED"] = "1"
    with log.open("a", encoding="utf-8") as handle:
        completed = subprocess.run(
            cmd,
            cwd=ROOT,
            env=env,
            stdout=handle,
            stderr=subprocess.STDOUT,
            check=False,
        )
    append_log(log, f"[{now()}] END returncode={completed.returncode}")
    return int(completed.returncode)


def parse_float(value: str | None) -> float | None:
    if value in {None, ""}:
        return None
    try:
        return float(value)
    except ValueError:
        return None


def parse_int(value: str | None) -> int | None:
    if value in {None, ""}:
        return None
    try:
        return int(value)
    except ValueError:
        return None


def row_selection_value(row: dict[str, str], metric: str) -> float | None:
    key = {
        "loss": "loss",
        "full-mse": "full_mse",
        "pearson": "pearson",
        "common128-mse": "common128_mse",
    }[metric]
    value = parse_float(row.get(key))
    if value is None:
        return None
    return -value if metric == "pearson" else value


def parse_best_metrics(config: RunConfig) -> RunResult:
    path = run_dir(config) / "metrics.tsv"
    if not path.exists():
        return RunResult(config=config, status="missing_metrics", returncode=1, gpu="")
    best_row: dict[str, str] | None = None
    best_value: float | None = None
    with path.open("r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if row.get("split") != "valid":
                continue
            value = row_selection_value(row, config.selection_metric)
            if value is None:
                continue
            if best_value is None or value < best_value:
                best_value = value
                best_row = row
    if best_row is None:
        return RunResult(config=config, status="missing_valid_metrics", returncode=1, gpu="")
    return RunResult(
        config=config,
        status="completed",
        returncode=0,
        gpu="",
        best_step=parse_int(best_row.get("step")),
        best_loss=parse_float(best_row.get("loss")),
        full_mse=parse_float(best_row.get("full_mse")),
        full_mae=parse_float(best_row.get("full_mae")),
        pearson=parse_float(best_row.get("pearson")),
        common128_mse=parse_float(best_row.get("common128_mse")),
        common128_mae=parse_float(best_row.get("common128_mae")),
        selection_value=parse_float(best_row.get("selection_value")),
    )


def run_config(config: RunConfig, *, gpu: str, date_tag: str, dry_run: bool) -> RunResult:
    checkpoint = run_dir(config) / "adapter_head_best.pt"
    if checkpoint.exists() and (run_dir(config) / "metrics.tsv").exists():
        result = parse_best_metrics(config)
        result.gpu = gpu
        return result
    code = run_command(train_command(config), gpu=gpu, log=log_path(config, date_tag), dry_run=dry_run)
    if code != 0:
        return RunResult(config=config, status="train_failed", returncode=code, gpu=gpu)
    if dry_run:
        return RunResult(config=config, status="dry_run", returncode=0, gpu=gpu)
    result = parse_best_metrics(config)
    result.gpu = gpu
    return result


def write_summary(path: Path, results: Iterable[RunResult]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    rows = list(results)
    fields = [
        "phase",
        "stage",
        "run_id",
        "status",
        "returncode",
        "gpu",
        "arch",
        "hidden_channels",
        "residual_scale_init",
        "linear_dilation",
        "target_space",
        "loss",
        "smooth_l1_beta",
        "hybrid_alpha",
        "lr",
        "weight_decay",
        "lr_schedule",
        "step_decay_steps",
        "step_decay_learning_rates",
        "init_checkpoint",
        "seed",
        "max_steps",
        "selection_metric",
        "best_step",
        "best_loss",
        "full_mse",
        "full_mae",
        "pearson",
        "common128_mse",
        "common128_mae",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for result in rows:
            config = result.config
            writer.writerow(
                {
                    "phase": config.phase,
                    "stage": config.stage,
                    "run_id": config.run_id,
                    "status": result.status,
                    "returncode": result.returncode,
                    "gpu": result.gpu,
                    "arch": config.arch,
                    "hidden_channels": config.hidden_channels,
                    "residual_scale_init": config.residual_scale_init,
                    "linear_dilation": config.linear_dilation,
                    "target_space": config.target_space,
                    "loss": config.loss,
                    "smooth_l1_beta": config.smooth_l1_beta,
                    "hybrid_alpha": config.hybrid_alpha,
                    "lr": config.lr,
                    "weight_decay": config.weight_decay,
                    "lr_schedule": config.lr_schedule,
                    "step_decay_steps": config.step_decay_steps,
                    "step_decay_learning_rates": config.step_decay_learning_rates,
                    "init_checkpoint": config.init_checkpoint,
                    "seed": config.seed,
                    "max_steps": config.max_steps,
                    "selection_metric": config.selection_metric,
                    "best_step": result.best_step,
                    "best_loss": result.best_loss,
                    "full_mse": result.full_mse,
                    "full_mae": result.full_mae,
                    "pearson": result.pearson,
                    "common128_mse": result.common128_mse,
                    "common128_mae": result.common128_mae,
                }
            )


def read_summary(path: Path) -> list[RunResult]:
    if not path.exists():
        return []
    results: list[RunResult] = []
    with path.open("r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            try:
                config = RunConfig(
                    phase=row["phase"],
                    stage=row["stage"],
                    arch=row["arch"],
                    seed=int(row["seed"]),
                    max_steps=int(row["max_steps"]),
                    lr=float(row["lr"]),
                    weight_decay=float(row["weight_decay"]),
                    loss=row["loss"],
                    smooth_l1_beta=float(row["smooth_l1_beta"]),
                    hybrid_alpha=float(row["hybrid_alpha"]),
                    target_space=row["target_space"],
                    hidden_channels=int(row["hidden_channels"]),
                    residual_scale_init=float(row["residual_scale_init"]),
                    linear_dilation=int(row["linear_dilation"]),
                    lr_schedule=row["lr_schedule"],
                    step_decay_steps=row.get("step_decay_steps", ""),
                    step_decay_learning_rates=row.get("step_decay_learning_rates", ""),
                    init_checkpoint=row.get("init_checkpoint", ""),
                    selection_metric=row.get("selection_metric", "full-mse"),
                    run_id_override=row.get("run_id", ""),
                )
            except (KeyError, TypeError, ValueError):
                continue
            results.append(
                RunResult(
                    config=config,
                    status=row.get("status", ""),
                    returncode=int(row.get("returncode") or 0),
                    gpu=row.get("gpu", ""),
                    best_step=parse_int(row.get("best_step")),
                    best_loss=parse_float(row.get("best_loss")),
                    full_mse=parse_float(row.get("full_mse")),
                    full_mae=parse_float(row.get("full_mae")),
                    pearson=parse_float(row.get("pearson")),
                    common128_mse=parse_float(row.get("common128_mse")),
                    common128_mae=parse_float(row.get("common128_mae")),
                )
            )
    return results


def run_stage(
    configs: list[RunConfig],
    *,
    gpus: list[str],
    date_tag: str,
    summary: Path,
    all_results: list[RunResult],
    dry_run: bool,
) -> list[RunResult]:
    print(f"[{now()}] stage_start\tn={len(configs)}", flush=True)
    tasks: queue.Queue[RunConfig] = queue.Queue()
    for config in configs:
        tasks.put(config)
    lock = threading.Lock()
    stage_results: list[RunResult] = []

    def worker(gpu: str) -> None:
        while True:
            try:
                config = tasks.get_nowait()
            except queue.Empty:
                return
            print(f"[{now()}] run_start\tgpu={gpu}\t{config.run_id}", flush=True)
            result = run_config(config, gpu=gpu, date_tag=date_tag, dry_run=dry_run)
            print(
                f"[{now()}] run_end\tgpu={gpu}\t{config.run_id}\t{result.status}\t"
                f"full_mse={result.full_mse}\tpearson={result.pearson}\tcommon128={result.common128_mse}",
                flush=True,
            )
            with lock:
                stage_results.append(result)
                all_results.append(result)
                if not dry_run:
                    write_summary(summary, all_results)
            tasks.task_done()

    threads = [threading.Thread(target=worker, args=(gpu,), daemon=False) for gpu in gpus]
    for thread in threads:
        thread.start()
    for thread in threads:
        thread.join()
    print(f"[{now()}] stage_end", flush=True)
    return stage_results


def completed(results: Iterable[RunResult]) -> list[RunResult]:
    return [result for result in results if result.status == "completed" and result.full_mse is not None]


def best_results(results: Iterable[RunResult], *, n: int, metric: str) -> list[RunResult]:
    def value(result: RunResult) -> float:
        if metric == "pearson":
            return -float(result.pearson)
        if metric == "common128-mse":
            return float(result.common128_mse)
        if metric == "loss":
            return float(result.best_loss)
        return float(result.full_mse)

    candidates = [
        result for result in completed(results)
        if (
            (metric == "pearson" and result.pearson is not None)
            or (metric == "common128-mse" and result.common128_mse is not None)
            or (metric == "loss" and result.best_loss is not None)
            or metric == "full-mse"
        )
    ]
    return sorted(candidates, key=value)[:n]


def phase1_configs(selection_metric: str) -> list[RunConfig]:
    return [
        RunConfig(
            phase="phase1",
            stage="winner_stability_5000",
            arch=arch,
            seed=seed,
            max_steps=5000,
            lr=1e-3,
            loss="smooth-l1",
            smooth_l1_beta=1.0,
            selection_metric=selection_metric,
        )
        for arch in ["residual-conv3", "conv3", "conv5"]
        for seed in [20260515, 20260522, 20260523]
    ]


def stable_head_from_phase1(results: list[RunResult]) -> str:
    grouped: dict[str, list[RunResult]] = {}
    for result in completed(results):
        grouped.setdefault(result.config.arch, []).append(result)
    if not grouped:
        raise RuntimeError("Phase 1 has no completed runs")
    return min(
        grouped,
        key=lambda arch: (
            sum(float(item.full_mse) for item in grouped[arch]) / len(grouped[arch]),
            len(grouped[arch]) * -1,
        ),
    )


def phase2_screen_configs(selection_metric: str) -> list[RunConfig]:
    configs: list[RunConfig] = []
    for hidden in [128, 256, 512]:
        for scale in [0.01, 0.1, 1.0]:
            configs.append(
                RunConfig(
                    phase="phase2",
                    stage="head_screen_1000",
                    arch="residual-conv3",
                    hidden_channels=hidden,
                    residual_scale_init=scale,
                    seed=20260522,
                    max_steps=1000,
                    selection_metric=selection_metric,
                )
            )
    configs.extend(
        [
            RunConfig(phase="phase2", stage="head_screen_1000", arch="conv3x2", seed=20260522, max_steps=1000, selection_metric=selection_metric),
            RunConfig(phase="phase2", stage="head_screen_1000", arch="conv7", seed=20260522, max_steps=1000, selection_metric=selection_metric),
            RunConfig(phase="phase2", stage="head_screen_1000", arch="dilated-conv3", linear_dilation=2, seed=20260522, max_steps=1000, selection_metric=selection_metric),
            RunConfig(phase="phase2", stage="head_screen_1000", arch="dilated-conv3", linear_dilation=4, seed=20260522, max_steps=1000, selection_metric=selection_metric),
        ]
    )
    return configs


def promote_configs(results: list[RunResult], *, selection_metric: str, n: int) -> list[RunConfig]:
    configs: list[RunConfig] = []
    for rank, result in enumerate(best_results(results, n=n, metric=selection_metric), start=1):
        configs.append(
            replace(
                result.config,
                stage="head_promote_5000",
                max_steps=5000,
                eval_every=250,
                rank_label=f"p{rank}",
                run_id_override="",
            )
        )
    return configs


def phase3_configs(top_results: list[RunResult], *, selection_metric: str) -> list[RunConfig]:
    configs: list[RunConfig] = []
    top = best_results(top_results, n=2, metric=selection_metric)
    for rank, result in enumerate(top, start=1):
        base = result.config
        configs.extend(
            [
                replace(base, phase="phase3", stage="loss_sweep_5000", loss="smooth-l1", smooth_l1_beta=0.5, max_steps=5000, rank_label=f"top{rank}_b05", run_id_override=""),
                replace(base, phase="phase3", stage="loss_sweep_5000", loss="smooth-l1", smooth_l1_beta=2.0, max_steps=5000, rank_label=f"top{rank}_b2", run_id_override=""),
                replace(base, phase="phase3", stage="loss_sweep_5000", loss="hybrid", smooth_l1_beta=1.0, hybrid_alpha=0.5, max_steps=5000, rank_label=f"top{rank}_hybrid", run_id_override=""),
            ]
        )
    if top:
        winner = top[0]
        checkpoint = rel(run_dir(winner.config) / "adapter_head_best.pt")
        configs.extend(
            [
                replace(winner.config, phase="phase3", stage="mse_finetune_1000", loss="mse", lr=3e-4, max_steps=1000, init_checkpoint=checkpoint, rank_label="mseft3e4", run_id_override=""),
                replace(winner.config, phase="phase3", stage="mse_finetune_1000", loss="mse", lr=1e-4, max_steps=1000, init_checkpoint=checkpoint, rank_label="mseft1e4", run_id_override=""),
            ]
        )
    return configs


def phase4_configs(top_results: list[RunResult], *, selection_metric: str) -> list[RunConfig]:
    top = best_results(top_results, n=1, metric=selection_metric)
    if not top:
        raise RuntimeError("Phase 4 needs at least one completed prior run")
    base = top[0].config
    return [
        replace(base, phase="phase4", stage="lr_schedule_5000", lr_schedule="constant", lr=1e-3, max_steps=5000, rank_label="constant", run_id_override=""),
        replace(base, phase="phase4", stage="lr_schedule_5000", lr_schedule="cosine", lr=1e-3, max_steps=5000, rank_label="cosine", run_id_override=""),
        replace(base, phase="phase4", stage="lr_schedule_5000", lr_schedule="step", lr=1e-3, step_decay_steps="3000", step_decay_learning_rates="3e-4", max_steps=5000, rank_label="step3000", run_id_override=""),
        replace(base, phase="phase4", stage="lr_schedule_5000", lr_schedule="step", lr=1e-3, step_decay_steps="4000", step_decay_learning_rates="3e-4", max_steps=5000, rank_label="step4000", run_id_override=""),
    ]


def main() -> None:
    args = parse_args()
    gpus = [gpu.strip() for gpu in args.gpus.split(",") if gpu.strip()]
    if not gpus:
        raise ValueError("At least one GPU must be provided")
    summary = summary_path(args.date_tag, args.summary)
    all_results: list[RunResult] = read_summary(summary)

    print(f"root\t{ROOT}", flush=True)
    print(f"git_commit\t{git_head()}", flush=True)
    print(f"phase\t{args.phase}", flush=True)
    print(f"gpus\t{','.join(gpus)}", flush=True)
    print(f"selection_metric\t{args.selection_metric}", flush=True)
    print(f"summary\t{summary}", flush=True)
    print(f"dry_run\t{args.dry_run}", flush=True)

    phase1_results: list[RunResult] = [
        result for result in all_results if result.config.phase == "phase1"
    ]
    phase2_screen_results: list[RunResult] = [
        result for result in all_results
        if result.config.phase == "phase2" and result.config.stage == "head_screen_1000"
    ]
    phase2_promote_results: list[RunResult] = [
        result for result in all_results
        if result.config.phase == "phase2" and result.config.stage == "head_promote_5000"
    ]
    phase3_results: list[RunResult] = [
        result for result in all_results if result.config.phase == "phase3"
    ]

    if args.phase in {"phase1", "all"}:
        phase1_results = run_stage(
            phase1_configs(args.selection_metric),
            gpus=gpus,
            date_tag=args.date_tag,
            summary=summary,
            all_results=all_results,
            dry_run=args.dry_run,
        )
    if args.phase in {"phase2-screen", "all"}:
        if phase1_results and not args.dry_run:
            print(f"phase1_stable_head\t{stable_head_from_phase1(phase1_results)}", flush=True)
        phase2_screen_results = run_stage(
            phase2_screen_configs(args.selection_metric),
            gpus=gpus,
            date_tag=args.date_tag,
            summary=summary,
            all_results=all_results,
            dry_run=args.dry_run,
        )
    if args.phase in {"phase2-promote", "all"}:
        source = phase2_screen_results or phase1_results
        phase2_promote_results = run_stage(
            promote_configs(source, selection_metric=args.selection_metric, n=3),
            gpus=gpus,
            date_tag=args.date_tag,
            summary=summary,
            all_results=all_results,
            dry_run=args.dry_run,
        )
    if args.phase in {"phase3", "all"}:
        source = phase2_promote_results or phase2_screen_results or phase1_results
        phase3_results = run_stage(
            phase3_configs(source, selection_metric=args.selection_metric),
            gpus=gpus,
            date_tag=args.date_tag,
            summary=summary,
            all_results=all_results,
            dry_run=args.dry_run,
        )
    if args.phase in {"phase4", "all"}:
        source = phase3_results or phase2_promote_results or phase2_screen_results or phase1_results
        run_stage(
            phase4_configs(source, selection_metric=args.selection_metric),
            gpus=gpus,
            date_tag=args.date_tag,
            summary=summary,
            all_results=all_results,
            dry_run=args.dry_run,
        )

    if not args.dry_run:
        write_summary(summary, all_results)
    print(f"[{now()}] sweep_done", flush=True)


if __name__ == "__main__":
    main()
