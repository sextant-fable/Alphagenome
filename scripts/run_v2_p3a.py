#!/usr/bin/env python3
"""Safely extend the project environment and run the P3A checkpoint."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
from typing import Any

try:
    from .doh_connect_proxy import local_doh_proxy
except ImportError:
    from doh_connect_proxy import local_doh_proxy


REPO_ROOT = Path(__file__).resolve().parents[1]
STATE_PATH = REPO_ROOT / "alphagenome_custom/metadata/v2/execution_state.json"
ENV_AUDIT_DIR = REPO_ROOT / "alphagenome_custom/metadata/v2/environment"
ENV_NAME = "alphagenome"
MIN_FREE_BYTES = 50 * 1024**3
REQUIRED_APPROVALS = {
    "G1_large_source_download_or_realignment": "p3_full_rna_streaming_482",
    "G2_v2_manifest_scientific_review": "v2_manifest_candidate_hierarchy",
    "G3_dataset_or_cache_generation": "p3_full_outputs_and_p4_loader",
}
TOOL_PACKAGES = [
    "python=3.12.13",
    "star=2.7.11b",
    "samtools=1.23",
    "bedtools=2.31.1",
    "ucsc-bedgraphtobigwig",
]
PROTECTED_PACKAGES = {
    "alphagenome-pytorch",
    "numpy",
    "pybigwig",
    "python",
    "torch",
    "triton",
}
PROTECTED_PREFIXES = ("cuda", "nvidia-")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--plan-only",
        action="store_true",
        help="Print the immutable command plan without checking approvals or writing data.",
    )
    return parser.parse_args()


def command_plan(conda: str) -> dict[str, Any]:
    install_command = [
        conda,
        "install",
        "-y",
        "-n",
        ENV_NAME,
        "--freeze-installed",
        "-c",
        "conda-forge",
        "-c",
        "bioconda",
        *TOOL_PACKAGES,
    ]
    pilot_command = [
        conda,
        "run",
        "--no-capture-output",
        "-n",
        ENV_NAME,
        "python",
        "scripts/run_v2_reprocessing_pilot.py",
        "--manifest",
        "alphagenome_custom/metadata/v2/p3_pilot_sources.tsv",
        "--threads",
        "16",
        "--work-dir",
        "shared/source_reads/v2/pilot_20260713",
        "--output-dir",
        "alphagenome_custom/tracks/rna_seq_v2_normalized_pilot",
        "--star-index",
        "shared/reference_indexes/WBcel235_STAR_2.7.11b",
    ]
    return {
        "phase": "P3A",
        "environment": ENV_NAME,
        "minimum_free_bytes": MIN_FREE_BYTES,
        "dry_run_command": [*install_command, "--dry-run", "--json"],
        "install_command": install_command,
        "pilot_command": pilot_command,
    }


def require_scoped_approvals() -> None:
    if not STATE_PATH.is_file():
        raise RuntimeError(f"Missing workflow state: {STATE_PATH}")
    state = json.loads(STATE_PATH.read_text())
    if state.get("current_phase") != "P3A":
        raise RuntimeError(
            f"P3A may run only while current_phase=P3A, got {state.get('current_phase')}"
        )
    missing = []
    for key, scope in REQUIRED_APPROVALS.items():
        record = state.get("approvals", {}).get(key, {})
        if record.get("approved") is not True or record.get("scope") != scope:
            missing.append(f"{key}:{scope}")
    if missing:
        raise RuntimeError("Missing scoped approvals: " + ",".join(missing))


def environment_exists(conda: str) -> bool:
    completed = subprocess.run(
        [conda, "env", "list", "--json"],
        cwd=REPO_ROOT,
        check=True,
        capture_output=True,
        text=True,
    )
    environments = json.loads(completed.stdout).get("envs", [])
    return any(Path(path).name == ENV_NAME for path in environments)


def proxy_environment(proxy: str) -> dict[str, str]:
    environment = os.environ.copy()
    environment["HTTP_PROXY"] = proxy
    environment["HTTPS_PROXY"] = proxy
    environment["http_proxy"] = proxy
    environment["https_proxy"] = proxy
    environment["NO_PROXY"] = "127.0.0.1,localhost"
    environment["no_proxy"] = "127.0.0.1,localhost"
    return environment


def run_json(command: list[str], *, environment: dict[str, str] | None = None) -> Any:
    print("command\t" + " ".join(command), flush=True)
    completed = subprocess.run(
        command,
        cwd=REPO_ROOT,
        env=environment,
        check=True,
        capture_output=True,
        text=True,
    )
    return json.loads(completed.stdout)


def run(command: list[str], *, environment: dict[str, str] | None = None) -> None:
    print("command\t" + " ".join(command), flush=True)
    subprocess.run(command, cwd=REPO_ROOT, env=environment, check=True)


def environment_snapshot(conda: str) -> dict[str, Any]:
    conda_packages = run_json([conda, "list", "-n", ENV_NAME, "--json"])
    code = """
import importlib.metadata as metadata
import json
import sys
names = ['alphagenome-pytorch', 'numpy', 'pyBigWig', 'torch', 'triton']
versions = {}
for name in names:
    try:
        versions[name.lower()] = metadata.version(name)
    except metadata.PackageNotFoundError:
        versions[name.lower()] = None
try:
    import torch
    cuda = torch.version.cuda
    cuda_available = torch.cuda.is_available()
except Exception:
    cuda = None
    cuda_available = False
print(json.dumps({
    'python': sys.version.split()[0],
    'distributions': versions,
    'torch_cuda': cuda,
    'cuda_available': cuda_available,
}, sort_keys=True))
"""
    runtime = run_json(
        [conda, "run", "-n", ENV_NAME, "python", "-c", code]
    )
    return {"conda_packages": conda_packages, "runtime": runtime}


def package_name(record: dict[str, Any]) -> str:
    return str(record.get("name") or record.get("dist_name", "").split("-")[0]).lower()


def reject_protected_changes(transaction: dict[str, Any]) -> None:
    linked = {
        package_name(record): record
        for record in transaction.get("actions", {}).get("LINK", [])
    }
    unlinked = {
        package_name(record): record
        for record in transaction.get("actions", {}).get("UNLINK", [])
    }
    protected_changes = []
    for action in ("LINK", "UNLINK"):
        for record in transaction.get("actions", {}).get(action, []):
            name = package_name(record)
            if (
                name == "python"
                and name in linked
                and name in unlinked
                and str(linked[name].get("version"))
                == str(unlinked[name].get("version"))
            ):
                continue
            if name in PROTECTED_PACKAGES or name.startswith(PROTECTED_PREFIXES):
                protected_changes.append(f"{action}:{name}")
    if protected_changes:
        raise RuntimeError(
            "Conda would change protected model packages: "
            + ",".join(protected_changes)
        )


def tools_available(conda: str) -> bool:
    code = """
import shutil
names = ['STAR', 'samtools', 'bedtools', 'bedGraphToBigWig']
raise SystemExit(0 if all(shutil.which(name) for name in names) else 1)
"""
    completed = subprocess.run(
        [conda, "run", "-n", ENV_NAME, "python", "-c", code],
        cwd=REPO_ROOT,
        check=False,
    )
    return completed.returncode == 0


def verify_project_runtime(conda: str) -> None:
    code = """
import alphagenome_pytorch
import numpy
import pyBigWig
import torch
assert torch.cuda.is_available()
print('alphagenome_import_ok')
"""
    run([conda, "run", "-n", ENV_NAME, "python", "-c", code])
    run(
        [
            conda,
            "run",
            "-n",
            ENV_NAME,
            "python",
            "-m",
            "unittest",
            "discover",
            "-s",
            "tests",
            "-v",
        ]
    )


def prepare_environment(conda: str, plan: dict[str, Any]) -> None:
    ENV_AUDIT_DIR.mkdir(parents=True, exist_ok=True)
    before = environment_snapshot(conda)
    (ENV_AUDIT_DIR / "p3_before.json").write_text(
        json.dumps(before, indent=2, sort_keys=True) + "\n"
    )
    if not tools_available(conda):
        with local_doh_proxy() as proxy:
            environment = proxy_environment(proxy)
            transaction = run_json(
                plan["dry_run_command"], environment=environment
            )
            reject_protected_changes(transaction)
            (ENV_AUDIT_DIR / "p3_install_plan.json").write_text(
                json.dumps(transaction, indent=2, sort_keys=True) + "\n"
            )
            run(plan["install_command"], environment=environment)
    after = environment_snapshot(conda)
    (ENV_AUDIT_DIR / "p3_after.json").write_text(
        json.dumps(after, indent=2, sort_keys=True) + "\n"
    )
    if before["runtime"] != after["runtime"]:
        raise RuntimeError("Protected alphagenome runtime changed during tool installation")
    if not tools_available(conda):
        raise RuntimeError("RNA-seq tools are still unavailable after installation")
    verify_project_runtime(conda)


def main() -> None:
    args = parse_args()
    conda = shutil.which("conda")
    if conda is None:
        raise RuntimeError("conda is required for the existing alphagenome environment")
    plan = command_plan(conda)
    if args.plan_only:
        print(json.dumps(plan, indent=2, sort_keys=True))
        return

    require_scoped_approvals()
    if not environment_exists(conda):
        raise RuntimeError(f"Required project environment does not exist: {ENV_NAME}")
    free_bytes = shutil.disk_usage(REPO_ROOT).free
    if free_bytes < MIN_FREE_BYTES:
        raise RuntimeError(
            f"P3A requires at least {MIN_FREE_BYTES} free bytes; found {free_bytes}"
        )
    prepare_environment(conda, plan)
    run(plan["pilot_command"])


if __name__ == "__main__":
    main()
