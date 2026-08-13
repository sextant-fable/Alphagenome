#!/usr/bin/env python3
"""Score a frozen VCF with an audited v2 Model-B checkpoint.

The scorer consumes one checksum-locked JSON specification.  All input hashes,
VCF coordinates, and reference alleles are validated before model creation.
No split manifest is accepted or accessed.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import socket
import subprocess
import sys
import time
from typing import Any, Mapping

import numpy as np

from scripts import v2_biological_validation as biology
from scripts import v2_variant_scoring as variants


REPO_ROOT = Path(__file__).resolve().parents[1]
REQUIRED_INPUT_KEYS = {
    "checkpoint",
    "checkpoint_run",
    "fasta",
    "fai",
    "gtf",
    "group_manifest",
    "means",
    "track_manifest",
    "model_specs",
    "model_weights",
    "vcf",
}
TRACK_DELTA_FIELDS = (
    "variant_key",
    "chromosome",
    "position_1based",
    "identifier",
    "reference",
    "alternate",
    "source_line",
    "alternate_index",
    "window_start_0based",
    "window_end_0based",
    "anchor_index",
    "alignment_policy",
    "result_scope",
    "track_id",
    "resolution",
    "delta_sum",
    "delta_mean",
    "delta_l1",
    "max_abs_delta",
    "max_abs_delta_index",
)
IMPLEMENTATION_PATHS = (
    "scripts/score_v2_variants.py",
    "scripts/v2_variant_scoring.py",
    "scripts/v2_biological_validation.py",
    "scripts/train_v2_model.py",
    "scripts/v2_training_components.py",
)
EXPECTED_LORA_TARGETS = ["tower.blocks.8.mha", "tower.blocks.8.mlp"]
GENE_DELTA_FIELDS = (
    "variant_key",
    "chromosome",
    "position_1based",
    "reference",
    "alternate",
    "alignment_policy",
    "gene_id",
    "gene_name",
    "strand",
    "track_id",
    "resolution",
    "gene_body_bases",
    "exon_bases",
    "gene_body_delta_sum",
    "gene_body_delta_mean",
    "exon_delta_sum",
    "exon_delta_mean",
)
ISM_FIELDS = (
    "parent_variant_key",
    "chromosome",
    "position_1based",
    "reference",
    "alternate",
    "sequence_index",
    "track_id",
    "resolution",
    "delta_sum",
    "delta_mean",
    "delta_l1",
    "max_abs_delta",
    "max_abs_delta_index",
)


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def read_json(path: str | Path) -> dict[str, Any]:
    payload = json.loads(Path(path).read_text())
    if not isinstance(payload, dict):
        raise ValueError(f"Expected JSON object: {path}")
    return payload


def git_commit() -> str:
    return subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=REPO_ROOT,
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()


def implementation_sha256() -> dict[str, str]:
    return {
        relative: biology.sha256(REPO_ROOT / relative)
        for relative in IMPLEMENTATION_PATHS
    }


def resolve_input(path: str) -> Path:
    candidate = Path(path)
    resolved = candidate if candidate.is_absolute() else REPO_ROOT / candidate
    biology.assert_no_final_test_path(resolved)
    return resolved


def resolve_output(path: str, spec_path: Path) -> Path:
    candidate = Path(path)
    resolved = candidate if candidate.is_absolute() else spec_path.parent / candidate
    biology.assert_no_final_test_path(resolved)
    return resolved


def validate_variant_spec(
    spec_path: str | Path,
) -> tuple[dict[str, Any], dict[str, Path], list[variants.AlleleWindow], dict[str, Any]]:
    """Validate the frozen scorer contract without loading torch or a model."""

    spec_path = Path(spec_path).resolve()
    spec = read_json(spec_path)
    model = spec.get("model", {})
    inference = spec.get("inference", {})
    ism = spec.get("ism", {})
    if not (
        spec.get("schema_version") == 1
        and spec.get("contract") == variants.VARIANT_SCORING_CONTRACT
        and spec.get("final_test_access") == "prohibited"
        and spec.get("device", "cpu") in {"cpu", "cuda"}
        and model.get("model_id") == "B"
        and model.get("loss") in {"paper", "log1p_mse"}
        and int(model.get("fold", 0)) in range(1, 6)
        and int(model.get("seed", 0)) > 0
        and int(model.get("sequence_length", 0)) > 0
        and int(model.get("sequence_length", 0)) % 128 == 0
        and int(model.get("hidden_channels", 0)) > 0
        and model.get("mean_column") == f"fold_{int(model.get('fold', 0))}_train_nonzero_mean"
        and inference.get("forward_reverse_complement_ensemble") is True
        and inference.get("allele_window_policy")
        == "variant_start_anchored_right_context_crop_or_extend"
        and inference.get("track_semantics")
        == "all_241_tracks_unstranded_no_rc_permutation"
    ):
        raise RuntimeError("Variant-scoring spec violates the v1 frozen contract")
    if bool(ism.get("enabled", False)):
        if int(ism.get("radius_bp", -1)) < 0 or int(ism.get("max_mutants_per_variant", 0)) <= 0:
            raise RuntimeError("Enabled ISM requires radius_bp >= 0 and a positive mutant cap")
    locked = spec.get("locked_inputs", {})
    hashes = spec.get("input_sha256", {})
    if set(locked) != REQUIRED_INPUT_KEYS or set(hashes) != REQUIRED_INPUT_KEYS:
        raise RuntimeError(
            f"Variant spec requires exactly locked inputs {sorted(REQUIRED_INPUT_KEYS)}"
        )
    paths = {key: resolve_input(str(value)) for key, value in locked.items()}
    verified: dict[str, str] = {}
    for key in sorted(paths):
        if not paths[key].is_file():
            raise FileNotFoundError(f"Missing variant-scoring input {key}: {paths[key]}")
        actual = biology.sha256(paths[key])
        if actual != hashes[key]:
            raise RuntimeError(f"Variant-scoring input hash mismatch: {key}")
        verified[key] = actual

    run = read_json(paths["checkpoint_run"])
    expected_run = {
        "fold": int(model["fold"]),
        "seed": int(model["seed"]),
        "loss": model["loss"],
        "sequence_length": int(model["sequence_length"]),
        "hidden_channels": int(model["hidden_channels"]),
        "mean_column": model["mean_column"],
        "checkpoint_sha256": hashes["checkpoint"],
        "checkpoint_reload_verified": True,
    }
    mismatches = {
        key: (run.get(key), value)
        for key, value in expected_run.items()
        if run.get(key) != value
    }
    if run.get("model") != "B" and run.get("model_id") != "B":
        mismatches["model"] = (run.get("model", run.get("model_id")), "B")
    if mismatches:
        raise RuntimeError(f"Checkpoint run contract mismatch: {mismatches}")

    model_specs = read_json(paths["model_specs"])
    model_b_rows = [
        row for row in model_specs.get("models", []) if row.get("model_id") == "B"
    ]
    if not (
        model_specs.get("schema_version") == 1
        and model_specs.get("n_tracks") == 241
        and model_specs.get("resolutions") == [1, 128]
        and model_specs.get("model_b_organism_index") == 2
        and model_specs.get("model_b_lora_rank") == 8
        and model_specs.get("model_b_lora_alpha") == 16
        and model_specs.get("model_b_lora_targets") == EXPECTED_LORA_TARGETS
        and len(model_b_rows) == 1
        and model_b_rows[0].get("base_organism_index") == 2
        and model_b_rows[0].get("trunk_policy")
        == "worm_embeddings_and_lora_only"
    ):
        raise RuntimeError("Frozen model_specs does not describe the expected Model B")

    means = biology.read_tsv(paths["means"])
    tracks = biology.read_tsv(paths["track_manifest"])
    if (
        len(means) != 241
        or len(tracks) != 241
        or len({row["group_id"] for row in tracks}) != 241
        or [row["group_id"] for row in means] != [row["group_id"] for row in tracks]
    ):
        raise RuntimeError("Means and track manifest order differ")
    groups = {row["group_id"]: row for row in biology.read_tsv(paths["group_manifest"])}
    if set(groups) != {row["group_id"] for row in tracks} or any(
        groups[row["group_id"]].get("strand") != "." for row in tracks
    ):
        raise RuntimeError("RC no-permutation contract requires all tracks to be unstranded")
    if model["mean_column"] not in means[0]:
        raise RuntimeError(f"Unknown track-mean column: {model['mean_column']}")
    fai = variants.read_fai(paths["fai"])
    alleles = list(
        variants.parse_vcf(paths["vcf"], require_pass=bool(inference.get("require_pass", False)))
    )
    if not alleles:
        raise RuntimeError("No scoreable VCF alleles remain after filtering")
    maximum = int(inference.get("maximum_alleles", 1000))
    if maximum <= 0 or len(alleles) > maximum:
        raise RuntimeError(f"VCF allele count {len(alleles)} exceeds maximum_alleles={maximum}")
    windows = [
        variants.build_allele_window(
            allele,
            fasta_path=paths["fasta"],
            fai=fai,
            sequence_length=int(model["sequence_length"]),
            anchor_index=(
                None
                if inference.get("anchor_index") is None
                else int(inference["anchor_index"])
            ),
        )
        for allele in alleles
    ]
    output_dir = resolve_output(str(spec.get("output_dir", "variant_scoring_output")), spec_path)
    preflight = {
        "schema_version": 1,
        "contract": variants.VARIANT_SCORING_CONTRACT,
        "status": "passed",
        "checked_at": utc_now(),
        "spec_path": str(spec_path),
        "spec_sha256": biology.sha256(spec_path),
        "verified_input_sha256": verified,
        "alleles": len(alleles),
        "coordinate_aligned_gene_aggregation_alleles": sum(
            variants.supports_coordinate_gene_aggregation(window.allele)
            for window in windows
        ),
        "skipped_indel_gene_aggregation": sum(
            not variants.supports_coordinate_gene_aggregation(window.allele)
            for window in windows
        ),
        "track_count": len(tracks),
        "track_semantics": "all_241_tracks_unstranded_no_rc_permutation",
        "sequence_length": int(model["sequence_length"]),
        "output_dir": str(output_dir),
        "model_loads": 0,
        "git_commit": git_commit(),
        "hostname": socket.gethostname(),
        "working_directory": str(REPO_ROOT),
        "python_executable": sys.executable,
        "implementation_sha256": implementation_sha256(),
        "final_test_access": "prohibited",
        "locked_test_block_signal_reads": 0,
    }
    preflight["cuda_authorization"] = validate_cuda_authorization(
        spec_path, spec, paths
    )
    return spec, paths, windows, preflight


class ModelBPredictor:
    """Callable Model-B adapter returning unscaled experimental-space tracks."""

    def __init__(self, spec: Mapping[str, Any], paths: Mapping[str, Path]) -> None:
        import torch

        from scripts import train_v2_model
        from scripts import v2_training_components as components

        self.torch = torch
        self.components = components
        self.model_contract = spec["model"]
        device_name = str(spec.get("device", "cpu"))
        if device_name == "cuda":
            physical_gpu = int(spec.get("physical_gpu", -1))
            if physical_gpu not in {2, 3}:
                raise RuntimeError("CUDA variant scoring is restricted to physical GPU 2 or 3")
            if os.environ.get("CUDA_VISIBLE_DEVICES") != str(physical_gpu):
                raise RuntimeError("CUDA_VISIBLE_DEVICES must equal the frozen physical_gpu")
            if not torch.cuda.is_available() or torch.cuda.device_count() != 1:
                raise RuntimeError("Variant scorer requires exactly one visible CUDA device")
            self.device = torch.device("cuda:0")
        elif device_name == "cpu":
            self.device = torch.device("cpu")
        else:
            raise RuntimeError("device must be cpu or cuda")
        if paths["model_weights"].resolve() != train_v2_model.WEIGHTS_PATH.resolve():
            raise RuntimeError("Frozen model_weights path differs from build_model's weight source")
        means_rows = biology.read_tsv(paths["means"])
        mean_column = str(self.model_contract["mean_column"])
        self.track_ids = [row["group_id"] for row in means_rows]
        self.fold_means = torch.tensor(
            [float(row[mean_column]) for row in means_rows], dtype=torch.float32
        )
        self.model = train_v2_model.build_model(
            "B",
            self.fold_means,
            self.device,
            int(self.model_contract["hidden_channels"]),
        )
        checkpoint = torch.load(paths["checkpoint"], map_location="cpu", weights_only=False)
        expected = {
            "model_id": "B",
            "fold": int(self.model_contract["fold"]),
            "seed": int(self.model_contract["seed"]),
            "loss": self.model_contract["loss"],
            "sequence_length": int(self.model_contract["sequence_length"]),
            "hidden_channels": int(self.model_contract["hidden_channels"]),
            "mean_column": mean_column,
        }
        mismatches = {
            key: (checkpoint.get(key), value)
            for key, value in expected.items()
            if checkpoint.get(key) != value
        }
        if mismatches:
            raise RuntimeError(f"Checkpoint tensor metadata mismatch: {mismatches}")
        trainable = {
            name for name, parameter in self.model.named_parameters() if parameter.requires_grad
        }
        trainable_state = checkpoint.get("trainable_model_state", {})
        if set(trainable_state) != trainable:
            raise RuntimeError("Checkpoint trainable parameter set mismatch")
        state = self.model.state_dict()
        for name, value in trainable_state.items():
            if name not in state or state[name].shape != value.shape:
                raise RuntimeError(f"Checkpoint tensor mismatch: {name}")
            state[name] = value
        self.model.load_state_dict(state, strict=True)
        self.model.eval()

    def __call__(self, batch: np.ndarray) -> dict[int, np.ndarray]:
        torch = self.torch
        inputs = torch.from_numpy(np.ascontiguousarray(batch)).to(self.device)
        autocast = self.device.type == "cuda"
        with torch.no_grad():
            with torch.autocast(
                device_type=self.device.type,
                dtype=torch.bfloat16,
                enabled=autocast,
            ):
                model_space = self.model(inputs)
            return {
                resolution: self.components.unscale_predictions_experimental_space(
                    values.float(), self.fold_means.to(self.device), resolution
                )
                .clamp_min(0)
                .cpu()
                .numpy()
                for resolution, values in model_space.items()
            }


def _prefix(window: variants.AlleleWindow) -> dict[str, object]:
    allele = window.allele
    return {
        "variant_key": allele.key,
        "chromosome": allele.chromosome,
        "position_1based": allele.position_1based,
        "identifier": allele.identifier,
        "reference": allele.reference,
        "alternate": allele.alternate,
        "source_line": allele.source_line,
        "alternate_index": allele.alternate_index,
        "window_start_0based": window.start,
        "window_end_0based": window.end,
        "anchor_index": window.anchor_index,
        "alignment_policy": window.alignment_policy,
        "result_scope": (
            "coordinate_aligned_variant_effect"
            if variants.supports_coordinate_gene_aggregation(allele)
            else "fixed_window_index_aligned_diagnostic"
        ),
    }


def validate_cuda_authorization(
    spec_path: Path,
    spec: Mapping[str, Any],
    paths: Mapping[str, Path],
) -> dict[str, Any] | None:
    """Require an execution-specific, hash-bound authorization for CUDA."""

    if spec.get("device", "cpu") != "cuda":
        return None
    raw_path = spec.get("cuda_authorization_path")
    if not raw_path:
        raise RuntimeError("CUDA variant scoring requires cuda_authorization_path")
    authorization_path = resolve_input(str(raw_path))
    authorization = read_json(authorization_path)
    execution_id = str(spec.get("execution_id", ""))
    physical_gpu = int(spec.get("physical_gpu", -1))
    expected = {
        "schema_version": 1,
        "approved": True,
        "scope": f"variant_scoring:{execution_id}",
        "execution_id": execution_id,
        "spec_sha256": biology.sha256(spec_path),
        "physical_gpu": physical_gpu,
        "checkpoint_sha256": biology.sha256(paths["checkpoint"]),
        "vcf_sha256": biology.sha256(paths["vcf"]),
    }
    mismatches = {
        key: (authorization.get(key), value)
        for key, value in expected.items()
        if authorization.get(key) != value
    }
    if not execution_id or physical_gpu not in {2, 3} or mismatches:
        raise RuntimeError(f"CUDA variant authorization mismatch: {mismatches}")
    return {
        "path": str(authorization_path),
        "sha256": biology.sha256(authorization_path),
        "scope": authorization["scope"],
        "execution_id": execution_id,
        "physical_gpu": physical_gpu,
    }


def run_scoring(spec_path: str | Path) -> dict[str, Any]:
    started = time.monotonic()
    spec_path = Path(spec_path).resolve()
    spec, paths, windows, preflight = validate_variant_spec(spec_path)
    output_dir = Path(preflight["output_dir"])
    output_dir.mkdir(parents=True, exist_ok=True)
    preflight_path = output_dir / "preflight_audit.json"
    biology.write_json(preflight_path, preflight)
    audit_path = output_dir / "scoring_audit.json"
    audit: dict[str, Any] = {
        **preflight,
        "status": "running",
        "started_at": utc_now(),
        "completed_at": None,
        "outputs": {},
        "failures": [],
    }
    biology.write_json(audit_path, audit)
    try:
        if biology.sha256(spec_path) != preflight["spec_sha256"]:
            raise RuntimeError("Variant-scoring spec changed after preflight")
        if implementation_sha256() != preflight["implementation_sha256"]:
            raise RuntimeError("Variant-scoring implementation changed after preflight")
        for key, path in paths.items():
            if biology.sha256(path) != preflight["verified_input_sha256"][key]:
                raise RuntimeError(f"Variant-scoring input changed after preflight: {key}")
        predictor = ModelBPredictor(spec, paths)
        genes, _ = biology.load_annotated_genes(paths["gtf"])
        track_rows: list[dict[str, object]] = []
        gene_rows: list[dict[str, object]] = []
        ism_rows: list[dict[str, object]] = []
        skipped_indel_gene_aggregation = 0
        for ordinal, window in enumerate(windows, 1):
            _, _, delta = variants.score_allele_window(predictor, window)
            prefix = _prefix(window)
            for resolution, values in sorted(delta.items()):
                for row in variants.summarize_track_deltas(
                    values, predictor.track_ids, resolution=resolution
                ):
                    track_rows.append({**prefix, **row})
                if variants.supports_coordinate_gene_aggregation(window.allele):
                    for row in variants.aggregate_variant_gene_exon_delta(
                        window,
                        values,
                        genes=genes,
                        track_ids=predictor.track_ids,
                        resolution=resolution,
                    ):
                        gene_rows.append(
                            {
                                "variant_key": window.allele.key,
                                "chromosome": window.allele.chromosome,
                                "position_1based": window.allele.position_1based,
                                "reference": window.allele.reference,
                                "alternate": window.allele.alternate,
                                "alignment_policy": window.alignment_policy,
                                **row,
                            }
                        )
            if not variants.supports_coordinate_gene_aggregation(window.allele):
                skipped_indel_gene_aggregation += 1
            ism = spec.get("ism", {})
            if bool(ism.get("enabled", False)):
                radius = int(ism["radius_bp"])
                start = max(0, window.anchor_index - radius)
                stop = min(len(window.reference_sequence), window.anchor_index + radius + 1)
                for mutation, mutation_delta in variants.score_ism(
                    predictor,
                    reference_sequence=window.reference_sequence,
                    window_start=window.start,
                    positions=range(start, stop),
                    max_mutants=int(ism["max_mutants_per_variant"]),
                ):
                    for resolution, values in sorted(mutation_delta.items()):
                        for row in variants.summarize_track_deltas(
                            values, predictor.track_ids, resolution=resolution
                        ):
                            ism_rows.append(
                                {
                                    "parent_variant_key": window.allele.key,
                                    "chromosome": window.allele.chromosome,
                                    "position_1based": mutation.genomic_position_1based,
                                    "reference": mutation.reference,
                                    "alternate": mutation.alternate,
                                    "sequence_index": mutation.sequence_index,
                                    **row,
                                }
                            )
            print(f"variant_scored\t{ordinal}/{len(windows)}\t{window.allele.key}", flush=True)

        output_paths: dict[str, Path] = {
            "variant_track_deltas": output_dir / "variant_track_deltas.tsv",
            "variant_gene_exon_deltas": output_dir / "variant_gene_exon_deltas.tsv",
            "preflight_audit": preflight_path,
        }
        biology.write_tsv(output_paths["variant_track_deltas"], track_rows, TRACK_DELTA_FIELDS)
        biology.write_tsv(output_paths["variant_gene_exon_deltas"], gene_rows, GENE_DELTA_FIELDS)
        if bool(spec.get("ism", {}).get("enabled", False)):
            output_paths["ism_track_deltas"] = output_dir / "ism_track_deltas.tsv"
            biology.write_tsv(output_paths["ism_track_deltas"], ism_rows, ISM_FIELDS)
        if skipped_indel_gene_aggregation != preflight["skipped_indel_gene_aggregation"]:
            raise RuntimeError("Indel gene-aggregation skip count changed after preflight")
        audit.update(
            {
                "status": "completed",
                "completed_at": utc_now(),
                "elapsed_seconds": time.monotonic() - started,
                "model_loads": 1,
                "variant_count": len(windows),
                "variant_track_rows": len(track_rows),
                "variant_gene_exon_rows": len(gene_rows),
                "ism_track_rows": len(ism_rows),
                "skipped_indel_gene_aggregation": skipped_indel_gene_aggregation,
                "outputs": {
                    key: {
                        "path": str(path),
                        "sha256": biology.sha256(path),
                        "size_bytes": path.stat().st_size,
                    }
                    for key, path in output_paths.items()
                },
                "final_test_access": "prohibited",
                "locked_test_block_signal_reads": 0,
            }
        )
        biology.write_json(audit_path, audit)
        return audit
    except Exception as error:
        audit.update(
            {
                "status": "failed",
                "completed_at": utc_now(),
                "elapsed_seconds": time.monotonic() - started,
                "failures": [str(error)],
            }
        )
        biology.write_json(audit_path, audit)
        raise


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--spec", type=Path, required=True)
    parser.add_argument(
        "--preflight-only",
        action="store_true",
        help="Validate hashes, VCF coordinates, and REF alleles without loading the model.",
    )
    return parser


def main(argv: list[str] | None = None) -> None:
    args = build_parser().parse_args(argv)
    if args.preflight_only:
        _, _, _, preflight = validate_variant_spec(args.spec)
        print(json.dumps(preflight, indent=2, sort_keys=True))
        return
    result = run_scoring(args.spec)
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
