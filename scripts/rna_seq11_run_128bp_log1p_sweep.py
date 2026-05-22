#!/usr/bin/env python3
"""Run the staged 128 bp log1p-MSE RNA-seq adapter sweep.

This launcher intentionally runs independent single-GPU jobs rather than DDP so
the training semantics stay close to the original selected 128 bp adapter run.
It never evaluates the test split.
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
BENCHMARK_VALID_MSE = 1.0609935
BENCHMARK_VALID_PEARSON = 0.58809888


@dataclass(frozen=True)
class SweepConfig:
    stage: str
    arch: str
    lr_label: str
    lr: float
    wd_label: str
    wd: float
    loss: str
    target_space: str
    seed: int
    max_steps: int
    eval_every: int = 250
    early_stopping_min_steps: int = 0
    early_stopping_patience: int | None = None
    rank_label: str = ""

    @property
    def target_label(self) -> str:
        if self.target_space == "binned128-log1p-mean":
            return "bin128"
        return "full"

    @property
    def loss_label(self) -> str:
        return self.loss.replace("-", "")

    @property
    def arch_label(self) -> str:
        return self.arch.replace("-", "")

    @property
    def run_id(self) -> str:
        pieces = [
            "rna_seq11_sweep",
            self.stage,
        ]
        if self.rank_label:
            pieces.append(self.rank_label)
        pieces.extend(
            [
                self.arch_label,
                self.target_label,
                self.loss_label,
                f"lr{self.lr_label}",
                f"wd{self.wd_label}",
                f"seed{self.seed}",
                f"{self.max_steps}steps",
            ]
        )
        return "_".join(pieces)

    @property
    def key(self) -> tuple[str, str, str, str, str]:
        return (
            self.arch,
            self.lr_label,
            self.wd_label,
            self.loss,
            self.target_space,
        )


@dataclass
class SweepResult:
    config: SweepConfig
    status: str
    returncode: int
    gpu: str
    train_best_step: int | None = None
    train_best_loss: float | None = None
    full_mse: float | None = None
    full_mae: float | None = None
    full_pearson: float | None = None
    common128_mse: float | None = None
    common128_mae: float | None = None
    common128_pearson: float | None = None


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--date-tag", default="20260522")
    parser.add_argument("--gpus", default="0,1,2,3")
    parser.add_argument(
        "--skip-final-seeds",
        action="store_true",
        help="Stop after promotion/10000-step candidates instead of final seed confirmation.",
    )
    return parser.parse_args()


def now() -> str:
    return datetime.now().isoformat(timespec="seconds")


def git_head() -> str:
    completed = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=ROOT,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        check=False,
    )
    if completed.returncode != 0:
        return f"unknown: {completed.stdout.strip()}"
    return completed.stdout.strip()


def append_log(path: Path, message: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8") as handle:
        handle.write(message)
        if not message.endswith("\n"):
            handle.write("\n")


def run_command(
    cmd: list[str],
    *,
    gpu: str,
    log_path: Path,
) -> int:
    env = os.environ.copy()
    env["CUDA_VISIBLE_DEVICES"] = gpu
    env["PYTHONUNBUFFERED"] = "1"
    append_log(log_path, f"\n[{now()}] START gpu={gpu}")
    append_log(log_path, " ".join(cmd))
    with log_path.open("a", encoding="utf-8") as handle:
        completed = subprocess.run(
            cmd,
            cwd=ROOT,
            env=env,
            stdout=handle,
            stderr=subprocess.STDOUT,
            check=False,
        )
    append_log(log_path, f"[{now()}] END returncode={completed.returncode}")
    return int(completed.returncode)


def run_dir(config: SweepConfig) -> Path:
    return ROOT / "runs" / config.run_id


def log_path(config: SweepConfig, date_tag: str) -> Path:
    return ROOT / "logs" / f"rna_seq11_128bp_log1p_sweep_{date_tag}" / f"{config.run_id}.log"


def train_command(config: SweepConfig) -> list[str]:
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
        str(run_dir(config).relative_to(ROOT)),
        "--head-type",
        "linear",
        "--embedding-resolution",
        "128",
        "--linear-head-architecture",
        config.arch,
        "--linear-hidden-channels",
        "256",
        "--linear-loss-type",
        config.loss,
        "--smooth-l1-beta",
        "1.0",
        "--linear-target-space",
        config.target_space,
        "--target-transform",
        "log1p",
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
        str(config.wd),
        "--seed",
        str(config.seed),
        "--device",
        "auto",
    ]
    if config.early_stopping_patience is not None:
        cmd.extend(
            [
                "--early-stopping-min-steps",
                str(config.early_stopping_min_steps),
                "--early-stopping-patience",
                str(config.early_stopping_patience),
            ]
        )
    return cmd


def eval_point_command(config: SweepConfig) -> list[str]:
    out_dir = run_dir(config)
    return [
        sys.executable,
        "-u",
        "scripts/alphagenome_rna_seq11_eval.py",
        "--dataset-dir",
        VALID_DIR,
        "--weights",
        WEIGHTS,
        "--checkpoint",
        str((out_dir / "adapter_head_best.pt").relative_to(ROOT)),
        "--point-metrics",
        "--metrics-output",
        str((out_dir / "valid_best_pointwise_metrics.tsv").relative_to(ROOT)),
        "--batch-size",
        "1",
        "--num-workers",
        "0",
        "--device",
        "auto",
    ]


def eval_common128_command(config: SweepConfig) -> list[str]:
    out_dir = run_dir(config)
    return [
        sys.executable,
        "-u",
        "scripts/alphagenome_rna_seq11_eval.py",
        "--dataset-dir",
        VALID_DIR,
        "--weights",
        WEIGHTS,
        "--checkpoint",
        str((out_dir / "adapter_head_best.pt").relative_to(ROOT)),
        "--common-128bp-metrics",
        "log1p-mean",
        "--metrics-output",
        str((out_dir / "valid_best_common128_log1pmean_metrics.tsv").relative_to(ROOT)),
        "--batch-size",
        "1",
        "--num-workers",
        "0",
        "--device",
        "auto",
    ]


def eval_diagnostic_command(config: SweepConfig) -> list[str]:
    out_dir = run_dir(config)
    return [
        sys.executable,
        "-u",
        "scripts/alphagenome_rna_seq11_eval.py",
        "--dataset-dir",
        VALID_DIR,
        "--weights",
        WEIGHTS,
        "--checkpoint",
        str((out_dir / "adapter_head_best.pt").relative_to(ROOT)),
        "--point-metrics",
        "--metrics-output",
        str((out_dir / "valid_best_pointwise_metrics.tsv").relative_to(ROOT)),
        "--diagnostic-output",
        str((out_dir / "valid_best_diagnostic.tsv").relative_to(ROOT)),
        "--batch-size",
        "1",
        "--num-workers",
        "0",
        "--device",
        "auto",
    ]


def parse_training_best(path: Path) -> tuple[int | None, float | None]:
    if not path.exists():
        return None, None
    best_step: int | None = None
    best_loss: float | None = None
    with path.open("r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if row.get("split") != "valid":
                continue
            try:
                step = int(row["step"])
                loss = float(row["loss"])
            except (KeyError, TypeError, ValueError):
                continue
            if best_loss is None or loss < best_loss:
                best_loss = loss
                best_step = step
    return best_step, best_loss


def parse_metric_tsv(path: Path) -> tuple[float | None, float | None, float | None]:
    if not path.exists():
        return None, None, None
    with path.open("r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if row.get("metric_scope") != "overall":
                continue
            try:
                return (
                    float(row["mse"]),
                    float(row["mae"]),
                    float(row["pearson"]),
                )
            except (KeyError, TypeError, ValueError):
                return None, None, None
    return None, None, None


def result_from_files(config: SweepConfig, *, status: str, gpu: str, returncode: int) -> SweepResult:
    out_dir = run_dir(config)
    best_step, best_loss = parse_training_best(out_dir / "metrics.tsv")
    full_mse, full_mae, full_pearson = parse_metric_tsv(
        out_dir / "valid_best_pointwise_metrics.tsv"
    )
    common_mse, common_mae, common_pearson = parse_metric_tsv(
        out_dir / "valid_best_common128_log1pmean_metrics.tsv"
    )
    return SweepResult(
        config=config,
        status=status,
        returncode=returncode,
        gpu=gpu,
        train_best_step=best_step,
        train_best_loss=best_loss,
        full_mse=full_mse,
        full_mae=full_mae,
        full_pearson=full_pearson,
        common128_mse=common_mse,
        common128_mae=common_mae,
        common128_pearson=common_pearson,
    )


def run_config(config: SweepConfig, *, gpu: str, date_tag: str) -> SweepResult:
    out_dir = run_dir(config)
    log = log_path(config, date_tag)
    checkpoint = out_dir / "adapter_head_best.pt"
    point_metrics = out_dir / "valid_best_pointwise_metrics.tsv"
    common_metrics = out_dir / "valid_best_common128_log1pmean_metrics.tsv"

    if not checkpoint.exists():
        code = run_command(train_command(config), gpu=gpu, log_path=log)
        if code != 0:
            return result_from_files(config, status="train_failed", gpu=gpu, returncode=code)
    else:
        append_log(log, f"[{now()}] SKIP train; checkpoint exists")

    if not point_metrics.exists():
        code = run_command(eval_point_command(config), gpu=gpu, log_path=log)
        if code != 0:
            return result_from_files(config, status="point_eval_failed", gpu=gpu, returncode=code)
    else:
        append_log(log, f"[{now()}] SKIP point eval; metrics exist")

    if not common_metrics.exists():
        code = run_command(eval_common128_command(config), gpu=gpu, log_path=log)
        if code != 0:
            return result_from_files(config, status="common128_eval_failed", gpu=gpu, returncode=code)
    else:
        append_log(log, f"[{now()}] SKIP common128 eval; metrics exist")

    return result_from_files(config, status="completed", gpu=gpu, returncode=0)


def write_summary(path: Path, results: Iterable[SweepResult]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    rows = list(results)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "stage",
                "run_id",
                "status",
                "returncode",
                "gpu",
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
        )
        for result in rows:
            config = result.config
            writer.writerow(
                [
                    config.stage,
                    config.run_id,
                    result.status,
                    result.returncode,
                    result.gpu,
                    config.arch,
                    config.target_space,
                    config.loss,
                    config.lr_label,
                    config.wd_label,
                    config.seed,
                    config.max_steps,
                    result.train_best_step,
                    result.train_best_loss,
                    result.full_mse,
                    result.full_mae,
                    result.full_pearson,
                    result.common128_mse,
                    result.common128_mae,
                    result.common128_pearson,
                ]
            )


def run_stage(
    name: str,
    configs: list[SweepConfig],
    *,
    gpus: list[str],
    date_tag: str,
    summary_path: Path,
    all_results: list[SweepResult],
) -> list[SweepResult]:
    print(f"[{now()}] stage_start\t{name}\tn={len(configs)}", flush=True)
    tasks: queue.Queue[SweepConfig] = queue.Queue()
    for config in configs:
        tasks.put(config)

    stage_results: list[SweepResult] = []
    lock = threading.Lock()

    def worker(gpu: str) -> None:
        while True:
            try:
                config = tasks.get_nowait()
            except queue.Empty:
                return
            print(f"[{now()}] run_start\tgpu={gpu}\t{config.run_id}", flush=True)
            result = run_config(config, gpu=gpu, date_tag=date_tag)
            print(
                f"[{now()}] run_end\tgpu={gpu}\t{config.run_id}\t"
                f"{result.status}\tfull_mse={result.full_mse}\t"
                f"pearson={result.full_pearson}\tcommon128={result.common128_mse}",
                flush=True,
            )
            with lock:
                stage_results.append(result)
                all_results.append(result)
                write_summary(summary_path, all_results)
            tasks.task_done()

    threads = [threading.Thread(target=worker, args=(gpu,), daemon=False) for gpu in gpus]
    for thread in threads:
        thread.start()
    for thread in threads:
        thread.join()

    failures = [result for result in stage_results if result.status != "completed"]
    print(
        f"[{now()}] stage_end\t{name}\tcompleted={len(stage_results) - len(failures)}\t"
        f"failed={len(failures)}",
        flush=True,
    )
    return stage_results


def completed(results: Iterable[SweepResult]) -> list[SweepResult]:
    return [
        result for result in results
        if result.status == "completed" and result.full_mse is not None
    ]


def best_by_full_mse(results: Iterable[SweepResult], *, n: int) -> list[SweepResult]:
    unique: dict[tuple[str, str, str, str, str], SweepResult] = {}
    for result in completed(results):
        key = result.config.key
        previous = unique.get(key)
        if previous is None or float(result.full_mse) < float(previous.full_mse):
            unique[key] = result
    return sorted(unique.values(), key=lambda item: float(item.full_mse))[:n]


def best_binned_by_common128(results: Iterable[SweepResult]) -> SweepResult | None:
    binned = [
        result for result in results
        if result.status == "completed"
        and result.config.target_space == "binned128-log1p-mean"
        and result.common128_mse is not None
    ]
    if not binned:
        return None
    return min(binned, key=lambda item: float(item.common128_mse))


def promotion_condition(result: SweepResult) -> bool:
    if result.full_mse is not None and float(result.full_mse) < BENCHMARK_VALID_MSE:
        return True
    if (
        result.full_pearson is not None
        and result.full_mse is not None
        and float(result.full_pearson) > BENCHMARK_VALID_PEARSON
        and float(result.full_mse) <= 1.07
    ):
        return True
    if result.train_best_step is not None and result.train_best_step >= result.config.max_steps - 1000:
        return True
    return False


def make_hyper_configs() -> list[SweepConfig]:
    lrs = [("1e-4", 1e-4), ("3e-4", 3e-4), ("1e-3", 1e-3)]
    wds = [("0", 0.0), ("1e-4", 1e-4), ("1e-3", 1e-3)]
    losses = ["mse", "smooth-l1"]
    return [
        SweepConfig(
            stage="hp1000",
            arch="conv1x1",
            lr_label=lr_label,
            lr=lr,
            wd_label=wd_label,
            wd=wd,
            loss=loss,
            target_space="full-log1p",
            seed=20260522,
            max_steps=1000,
        )
        for lr_label, lr in lrs
        for wd_label, wd in wds
        for loss in losses
    ]


def make_seed_stability_configs() -> list[SweepConfig]:
    return [
        SweepConfig(
            stage="seed5000",
            arch="conv1x1",
            lr_label="3e-4",
            lr=3e-4,
            wd_label="0",
            wd=0.0,
            loss="mse",
            target_space="full-log1p",
            seed=seed,
            max_steps=5000,
        )
        for seed in [20260522, 20260523]
    ]


def make_head_configs(best: SweepResult) -> list[SweepConfig]:
    base = best.config
    return [
        replace(
            base,
            stage="head1000",
            arch=arch,
            seed=20260522,
            max_steps=1000,
            early_stopping_min_steps=0,
            early_stopping_patience=None,
            rank_label="",
        )
        for arch in [
            "conv1x1",
            "mlp1x1",
            "conv3",
            "conv5",
            "residual-conv1x1",
            "residual-conv3",
        ]
    ]


def make_binned_configs(best: SweepResult) -> list[SweepConfig]:
    base = best.config
    return [
        replace(
            base,
            stage="bin1000",
            arch=arch,
            loss="mse",
            target_space="binned128-log1p-mean",
            seed=20260522,
            max_steps=1000,
            early_stopping_min_steps=0,
            early_stopping_patience=None,
            rank_label="",
        )
        for arch in ["conv1x1", "mlp1x1", "conv3"]
    ]


def make_promotion_configs(
    full_best: list[SweepResult],
    binned_best: SweepResult | None,
) -> list[SweepConfig]:
    configs: list[SweepConfig] = []
    for index, result in enumerate(full_best, start=1):
        configs.append(
            replace(
                result.config,
                stage="promote5000",
                seed=20260522,
                max_steps=5000,
                early_stopping_min_steps=4000,
                early_stopping_patience=8,
                rank_label=f"p{index}",
            )
        )
    if binned_best is not None:
        configs.append(
            replace(
                binned_best.config,
                stage="promote5000",
                seed=20260522,
                max_steps=5000,
                early_stopping_min_steps=4000,
                early_stopping_patience=8,
                rank_label="binp1",
            )
        )
    return configs


def make_10000_configs(promoted: list[SweepResult]) -> list[SweepConfig]:
    configs: list[SweepConfig] = []
    for result in promoted:
        if promotion_condition(result):
            configs.append(
                replace(
                    result.config,
                    stage="extend10000",
                    max_steps=10000,
                    early_stopping_min_steps=4000,
                    early_stopping_patience=8,
                )
            )
    return configs


def make_final_seed_configs(winner: SweepResult) -> list[SweepConfig]:
    config = winner.config
    return [
        replace(
            config,
            stage="finalseed",
            seed=seed,
            rank_label="winner",
        )
        for seed in [20260515, 20260522, 20260523]
    ]


def run_diagnostic(
    result: SweepResult,
    *,
    gpu: str,
    date_tag: str,
) -> None:
    out_dir = run_dir(result.config)
    diagnostic = out_dir / "valid_best_diagnostic.tsv"
    if diagnostic.exists():
        return
    run_command(eval_diagnostic_command(result.config), gpu=gpu, log_path=log_path(result.config, date_tag))


def main() -> None:
    args = parse_args()
    gpus = [gpu.strip() for gpu in args.gpus.split(",") if gpu.strip()]
    if not gpus:
        raise ValueError("At least one GPU id is required")

    sweep_root = ROOT / "runs" / f"rna_seq11_128bp_log1p_sweep_{args.date_tag}"
    sweep_root.mkdir(parents=True, exist_ok=True)
    summary_path = sweep_root / "sweep_results.tsv"
    all_results: list[SweepResult] = []

    print(f"[{now()}] sweep_start", flush=True)
    print(f"root\t{ROOT}", flush=True)
    print(f"git_commit\t{git_head()}", flush=True)
    print(f"gpus\t{','.join(gpus)}", flush=True)
    print(f"summary\t{summary_path}", flush=True)

    hyper_results = run_stage(
        "hyperparameter_screen_1000",
        make_hyper_configs(),
        gpus=gpus,
        date_tag=args.date_tag,
        summary_path=summary_path,
        all_results=all_results,
    )
    hyper_best = best_by_full_mse(hyper_results, n=1)
    if not hyper_best:
        raise RuntimeError("No completed hyperparameter run is available")
    best_hyper = hyper_best[0]
    print(
        f"[{now()}] best_hyper\t{best_hyper.config.run_id}\t"
        f"full_mse={best_hyper.full_mse}\tpearson={best_hyper.full_pearson}",
        flush=True,
    )

    seed_results = run_stage(
        "baseline_seed_stability_5000",
        make_seed_stability_configs(),
        gpus=gpus,
        date_tag=args.date_tag,
        summary_path=summary_path,
        all_results=all_results,
    )

    head_results = run_stage(
        "head_screen_1000",
        make_head_configs(best_hyper),
        gpus=gpus,
        date_tag=args.date_tag,
        summary_path=summary_path,
        all_results=all_results,
    )

    binned_results = run_stage(
        "binned_target_screen_1000",
        make_binned_configs(best_hyper),
        gpus=gpus,
        date_tag=args.date_tag,
        summary_path=summary_path,
        all_results=all_results,
    )

    full_best = best_by_full_mse([*hyper_results, *head_results], n=3)
    binned_best = best_binned_by_common128(binned_results)
    promotion_configs = make_promotion_configs(full_best, binned_best)
    promotion_results = run_stage(
        "promotion_5000",
        promotion_configs,
        gpus=gpus,
        date_tag=args.date_tag,
        summary_path=summary_path,
        all_results=all_results,
    )

    extend_configs = make_10000_configs(promotion_results)
    extend_results: list[SweepResult] = []
    if extend_configs:
        extend_results = run_stage(
            "extend_10000",
            extend_configs,
            gpus=gpus,
            date_tag=args.date_tag,
            summary_path=summary_path,
            all_results=all_results,
        )

    winner_pool = completed([*promotion_results, *extend_results])
    if not winner_pool:
        winner_pool = completed([*hyper_results, *head_results, *seed_results])
    if not winner_pool:
        raise RuntimeError("No completed run is available for winner selection")
    winner = min(winner_pool, key=lambda item: float(item.full_mse))
    print(
        f"[{now()}] winner_candidate\t{winner.config.run_id}\t"
        f"full_mse={winner.full_mse}\tpearson={winner.full_pearson}",
        flush=True,
    )

    if not args.skip_final_seeds:
        final_seed_results = run_stage(
            "final_seed_confirmation",
            make_final_seed_configs(winner),
            gpus=gpus,
            date_tag=args.date_tag,
            summary_path=summary_path,
            all_results=all_results,
        )
        completed_final = completed(final_seed_results)
        if completed_final:
            winner = min(completed_final, key=lambda item: float(item.full_mse))
            print(
                f"[{now()}] final_winner\t{winner.config.run_id}\t"
                f"full_mse={winner.full_mse}\tpearson={winner.full_pearson}",
                flush=True,
            )

    run_diagnostic(winner, gpu=gpus[0], date_tag=args.date_tag)
    print(f"[{now()}] sweep_done", flush=True)


if __name__ == "__main__":
    main()
