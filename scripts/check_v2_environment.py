#!/usr/bin/env python3
"""Fail fast when the Python environment cannot run the v2 workflow."""

from __future__ import annotations

import importlib
import importlib.metadata
import sys
from typing import Iterable


MINIMUM_PYTHON = (3, 12)
EXPECTED_ALPHAGENOME_VERSION = "0.3.1"


def module_version(distribution: str) -> str | None:
    try:
        return importlib.metadata.version(distribution)
    except importlib.metadata.PackageNotFoundError:
        return None


def require_import(module_name: str, attributes: Iterable[str] = ()) -> list[str]:
    try:
        module = importlib.import_module(module_name)
    except Exception as error:
        return [f"cannot import {module_name}: {error}"]
    missing = [name for name in attributes if not hasattr(module, name)]
    if missing:
        return [f"{module_name} is missing: {', '.join(missing)}"]
    return []


def main() -> int:
    errors: list[str] = []
    print(f"python: {sys.executable}")
    print(f"python_version: {sys.version.split()[0]}")
    if sys.version_info < MINIMUM_PYTHON:
        errors.append(
            "v2 requires Python 3.12+; create a new environment instead of "
            "installing into Python 3.11"
        )

    packages = {
        "alphagenome-pytorch": "alphagenome-pytorch",
        "numpy": "numpy",
        "pyBigWig": "pyBigWig",
        "pytest": "pytest",
        "torch": "torch",
    }
    for label, distribution in packages.items():
        version = module_version(distribution)
        print(f"{label}: {version or 'NOT INSTALLED'}")
        if version is None:
            errors.append(f"missing distribution: {distribution}")

    alpha_version = module_version("alphagenome-pytorch")
    if alpha_version and alpha_version != EXPECTED_ALPHAGENOME_VERSION:
        errors.append(
            "alphagenome-pytorch must be "
            f"{EXPECTED_ALPHAGENOME_VERSION}, found {alpha_version}"
        )

    errors.extend(
        require_import(
            "alphagenome_pytorch.heads",
            ("GenomeTracksHead", "predictions_scaling", "targets_scaling"),
        )
    )
    errors.extend(
        require_import("alphagenome_pytorch.losses", ("multinomial_loss",))
    )
    errors.extend(
        require_import(
            "alphagenome_pytorch.extensions.finetuning.adapters", ("apply_lora",)
        )
    )

    if errors:
        print("v2_environment: INVALID")
        for error in errors:
            print(f"- {error}")
        print("Install with: python -m pip install -r requirements.txt")
        return 1

    print("v2_environment: OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
