#!/usr/bin/env python3
"""Run a bounded, task-faithful Borzoi C. elegans adapter pilot.

The official Borzoi checkpoint contributes trunk weights only.  A new 241-head
softplus adapter is trained on the existing I--V v2 intervals.  Native 32-bp
outputs are summed into 128-bp bins for a secondary comparator endpoint.  This
pilot never reads chromosome X or the locked final-test intervals and does not
claim a 1-bp result.
"""

from __future__ import annotations

import csv
import hashlib
import json
import os
from pathlib import Path
import random
import socket
import time

import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[1]
METADATA = REPO_ROOT / "alphagenome_custom/metadata/v2"
SPEC_PATH = METADATA / "p25_borzoi_adapter_pilot_spec.json"
EXECUTION_PATH = METADATA / "p25_borzoi_adapter_pilot_execution.json"
VENDOR = Path("/home/zelinli6/vendor/borzoi")
WEIGHTS = Path("/home/zelinli6/vendor/borzoi_weights/model0_best.h5")
TRACK_MEANS = REPO_ROOT / "results/v2_p16_iv_training_normalization/track_nonzero_means.tsv"
TRAIN_INTERVALS = REPO_ROOT / "results/v2_p15_iv_training_interval_contract/intervals/fold_1/train.tsv"
VALID_INTERVALS = REPO_ROOT / "results/v2_p15_iv_training_interval_contract/intervals/fold_1/valid.tsv"
TRACK_MANIFEST = REPO_ROOT / "alphagenome_custom/metadata/v2/replicate_holdout_v1/training_track_manifest.tsv"
GROUP_MANIFEST = REPO_ROOT / "alphagenome_custom/metadata/v2/replicate_holdout_v1/training_group_manifest.tsv"
FASTA = REPO_ROOT / "alphagenome_custom/reference/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa"
FAI = REPO_ROOT / "alphagenome_custom/reference/genome.fa.fai"
GTF = REPO_ROOT / "alphagenome_custom/reference/Caenorhabditis_elegans.WBcel235.115.gtf"
OUTPUT_ROOT = REPO_ROOT / "runs/v2_p25_borzoi_adapter_pilot"
LOG_ROOT = REPO_ROOT / "logs/v2_p25_borzoi_adapter_pilot"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def utc_now() -> str:
    from datetime import datetime, timezone

    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def require_scope() -> None:
    state = json.loads((METADATA / "execution_state.json").read_text())
    approval = state.get("approvals", {}).get("G4_gpu_experiments", {})
    if not (
        state.get("current_phase") == "P25"
        and state.get("status") == "RUNNING"
        and approval.get("approved") is True
        and approval.get("scope") == "p25_borzoi_adapter_pilot"
    ):
        raise RuntimeError("P25 requires controller P25/RUNNING and G4:p25_borzoi_adapter_pilot")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_fai(path: Path) -> dict[str, tuple[int, int, int, int]]:
    result = {}
    for line in path.read_text().splitlines():
        name, length, offset, line_bases, line_width = line.split("\t")
        result[name] = (int(length), int(offset), int(line_bases), int(line_width))
    return result


def fetch_sequence(path: Path, fai: dict[str, tuple[int, int, int, int]], chrom: str, start: int, end: int) -> bytes:
    length, offset, line_bases, line_width = fai[chrom]
    if start < 0 or end > length or end <= start:
        raise ValueError(f"invalid sequence interval {chrom}:{start}-{end}")
    start_line = start // line_bases
    start_in_line = start % line_bases
    byte_start = offset + start_line * line_width + start_in_line
    line_count = ((start_in_line + end - start - 1) // line_bases) + 1
    byte_count = (end - start) + line_count * (line_width - line_bases)
    with path.open("rb") as handle:
        handle.seek(byte_start)
        sequence = handle.read(byte_count).replace(b"\n", b"").replace(b"\r", b"")[: end - start]
    if len(sequence) != end - start:
        raise RuntimeError(f"short FASTA read {chrom}:{start}-{end}")
    return sequence.upper()


def sequence_one_hot(sequence: bytes) -> np.ndarray:
    mapping = {ord("A"): 0, ord("C"): 1, ord("G"): 2, ord("T"): 3}
    result = np.zeros((len(sequence), 4), dtype=np.float32)
    for index, base in enumerate(sequence):
        channel = mapping.get(base)
        if channel is not None:
            result[index, channel] = 1.0
    return result


def load_target_item(row: dict[str, str], track_paths: list[Path], fai: dict[str, tuple[int, int, int, int]], fasta: Path) -> tuple[np.ndarray, np.ndarray]:
    import pyBigWig

    chrom = row["chromosome"]
    start = int(row["start"])
    end = int(row["end"])
    if end - start != 131072 or chrom not in {"I", "II", "III", "IV", "V"}:
        raise ValueError(f"P25 expects I-V 131072-bp windows, got {chrom}:{start}-{end}")
    sequence = sequence_one_hot(fetch_sequence(fasta, fai, chrom, start, end))
    targets = []
    handles = [pyBigWig.open(str(path)) for path in track_paths]
    try:
        for handle in handles:
            values = np.asarray(handle.values(chrom, start, end, numpy=True), dtype=np.float32)
            values = np.nan_to_num(values, nan=0.0, posinf=0.0, neginf=0.0)
            if values.shape != (131072,) or np.any(values < 0):
                raise RuntimeError(f"invalid target signal for {chrom}:{start}-{end}")
            targets.append(values.reshape(1024, 128).sum(axis=1))
    finally:
        for handle in handles:
            handle.close()
    return sequence, np.stack(targets, axis=0).astype(np.float32)


def audit_intervals(path: Path) -> None:
    rows = read_tsv(path)
    if not rows:
        raise RuntimeError(f"empty interval manifest: {path}")
    forbidden = {"X", "chromosome_X"}
    if any(row["chromosome"] in forbidden or row.get("role") == "test_locked" for row in rows):
        raise PermissionError(f"P25 interval manifest has forbidden X/locked rows: {path}")
    if any(row["chromosome"] not in {"I", "II", "III", "IV", "V"} for row in rows):
        raise PermissionError(f"P25 interval manifest is outside I-V: {path}")


def load_model_and_dataset():
    # Imports are intentionally delayed so the official Borzoi environment is
    # isolated from the PyTorch AlphaGenome environment.
    import tensorflow as tf
    from baskerville import seqnn
    os.environ.setdefault("TF_CPP_MIN_LOG_LEVEL", "2")
    tf.config.run_functions_eagerly(False)
    for device in tf.config.list_physical_devices("GPU"):
        tf.config.experimental.set_memory_growth(device, True)

    params = json.loads((VENDOR / "examples/params.json").read_text())
    model_params = params["model"]
    model_params["verbose"] = False
    model_params["seq_length"] = 131072
    model_params.pop("head_human", None)
    model_params.pop("head_mouse", None)
    # The official 5,120-bin crop is tied to the 524-kb native input. Scaling
    # the input to 131 kb requires the proportional 1,280-bin crop to retain
    # a 4,096-bin native 32-bp output (then 4-to-1 pooled to 128 bp).
    for block in model_params["trunk"]:
        if block.get("name") == "Cropping1D":
            block["cropping"] = 1280
    model_params["head_c_elegans"] = {
        "name": "final",
        "units": 241,
        "activation": "softplus",
    }
    model_params["augment_rc"] = False
    seqnn_model = seqnn.SeqNN(model_params)
    seqnn_model.restore(str(WEIGHTS), trunk=True)

    trunk_layer_ids = {id(layer) for layer in seqnn_model.model_trunk.layers}
    for layer in seqnn_model.model.layers:
        layer.trainable = id(layer) not in trunk_layer_ids
    trainable = seqnn_model.model.trainable_variables
    if not trainable:
        raise RuntimeError("Borzoi adapter has no trainable head variables")

    track_rows = read_tsv(TRACK_MANIFEST)
    track_paths = [REPO_ROOT / row["output_path"] for row in track_rows]
    if len(track_paths) != 241 or not all(path.is_file() for path in track_paths):
        raise RuntimeError("P25 requires all 241 v2 training BigWigs")
    fai = read_fai(FAI)
    train = read_tsv(TRAIN_INTERVALS)
    valid = read_tsv(VALID_INTERVALS)
    return tf, seqnn_model, train, valid, trainable, track_paths, fai


def main() -> None:
    require_scope()
    spec = json.loads(SPEC_PATH.read_text())
    audit_intervals(TRAIN_INTERVALS)
    audit_intervals(VALID_INTERVALS)
    for path in (VENDOR / "LICENSE", VENDOR / "README.md", WEIGHTS, TRACK_MEANS, FASTA, FAI, GTF):
        if not path.is_file():
            raise FileNotFoundError(path)

    seed = int(spec["pilot"]["seed"])
    random.seed(seed)
    np.random.seed(seed)
    tf, model, train, valid, trainable, track_paths, fai = load_model_and_dataset()
    means_rows = read_tsv(TRACK_MEANS)
    means = np.asarray([float(row["fold_1_train_nonzero_mean"]) for row in means_rows], dtype=np.float32)
    if means.shape != (241,) or np.any(means <= 0):
        raise RuntimeError(f"invalid P16 means shape={means.shape}")
    optimizer = tf.keras.optimizers.Adam(learning_rate=6e-5, clipnorm=0.15)
    log_dir = LOG_ROOT / "fold_1_seed_20260823"
    log_dir.mkdir(parents=True, exist_ok=True)
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    metrics_path = log_dir / "train.jsonl"
    max_steps = int(spec["pilot"]["max_steps"])
    started = time.monotonic()
    rows = []
    for step in range(max_steps):
        sequence, target_array = load_target_item(train[step % len(train)], track_paths, fai, FASTA)
        x = tf.convert_to_tensor(sequence[None, ...], dtype=tf.float32)
        target = tf.convert_to_tensor(target_array[None, ...], dtype=tf.float32)
        target = target / tf.convert_to_tensor(means[None, :, None], dtype=tf.float32)
        with tf.GradientTape() as tape:
            # Keep the pretrained trunk in inference mode while optimizing
            # only the new C. elegans head; this prevents BatchNorm/dropout
            # state changes from being mistaken for adapter learning.
            pred32 = model.model(x, training=False)
            pred32 = tf.reshape(pred32, [1, 1024, 4, 241])
            pred128 = tf.reduce_sum(pred32, axis=2)
            prediction = tf.transpose(pred128, [0, 2, 1])
            loss = tf.reduce_mean(tf.square(tf.math.log1p(prediction) - tf.math.log1p(target)))
            loss += tf.add_n(model.model.losses) if model.model.losses else 0.0
        gradients = tape.gradient(loss, trainable)
        optimizer.apply_gradients(zip(gradients, trainable))
        record = {"step": step + 1, "loss": float(loss.numpy())}
        rows.append(record)
        with metrics_path.open("a") as handle:
            handle.write(json.dumps(record, sort_keys=True) + "\n")

    # Small deterministic validation pass; this is a pilot diagnostic, not a
    # replacement for the registered P17/P19 primary metric implementation.
    preds = []
    trues = []
    for index in range(min(len(valid), 8)):
        sequence, true = load_target_item(valid[index], track_paths, fai, FASTA)
        x = tf.convert_to_tensor(sequence[None, ...], dtype=tf.float32)
        pred32 = model.model(x, training=False)
        pred128 = tf.reduce_sum(tf.reshape(pred32, [1, 1024, 4, 241]), axis=2)
        pred = pred128.numpy()[0].transpose(1, 0) * means[:, None]
        preds.append(pred)
        trues.append(true)
    pred_flat = np.concatenate(preds, axis=1)
    true_flat = np.concatenate(trues, axis=1)
    corrs = []
    for channel in range(241):
        if np.std(pred_flat[channel]) > 0 and np.std(true_flat[channel]) > 0:
            corrs.append(float(np.corrcoef(pred_flat[channel], true_flat[channel])[0, 1]))
    execution = {
        "schema_version": 1,
        "phase": "P25",
        "status": "completed_pilot",
        "adapter_status": "pilot_completed",
        "started_at": utc_now(),
        "hostname": socket.gethostname(),
        "working_directory": str(REPO_ROOT),
        "controller_invocation": "python -m scripts.v2_phase_controller run --phase P25",
        "registered_phase_command": "python -m scripts.run_v2_p25_borzoi_adapter_pilot",
        "spec_sha256": sha256(SPEC_PATH),
        "official_commit": subprocess_check_output(["git", "-C", str(VENDOR), "rev-parse", "HEAD"]),
        "weight_path": str(WEIGHTS),
        "weight_sha256": sha256(WEIGHTS),
        "trainable_variable_count": len(trainable),
        "trainable_parameter_count": int(sum(np.prod(v.shape) for v in trainable)),
        "steps": max_steps,
        "validation_intervals": min(len(valid), 8),
        "mean_track_pearson_128bp": float(np.mean(corrs)) if corrs else None,
        "validation_metric_status": "pilot_descriptive_only",
        "chromosome_x_reads": 0,
        "locked_test_block_signal_reads": 0,
        "native_checkpoint_head_used": False,
        "native_trunk_weights_used": True,
        "one_bp_result": "not_available",
        "output_dir": str(OUTPUT_ROOT),
        "log_path": str(metrics_path),
        "elapsed_seconds": time.monotonic() - started,
        "next_gate": "R25 review and, only if justified, a larger I-V adapter matrix",
    }
    EXECUTION_PATH.write_text(json.dumps(execution, indent=2, sort_keys=True) + "\n")
    print(json.dumps(execution, indent=2, sort_keys=True))


def subprocess_check_output(command: list[str]) -> str:
    import subprocess

    return subprocess.check_output(command, text=True).strip()


if __name__ == "__main__":
    main()
