#!/usr/bin/env python3
"""Fail fast when the Python environment cannot run the v2 workflow."""

from __future__ import annotations

import importlib
import importlib.metadata
import sys
from typing import Iterable


EXPECTED_PYTHON = (3, 12, 13)
EXPECTED_DISTRIBUTIONS = {
    "alphagenome-pytorch": "0.3.1",
    "numpy": "2.4.4",
    "pyBigWig": "0.3.25",
    "torch": "2.11.0+cu128",
    "triton": "3.6.0",
}


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
    if sys.version_info[:3] != EXPECTED_PYTHON:
        errors.append(
            "historical v2 runtime requires Python "
            f"{'.'.join(map(str, EXPECTED_PYTHON))}"
        )

    for distribution, expected_version in EXPECTED_DISTRIBUTIONS.items():
        version = module_version(distribution)
        print(f"{distribution}: {version or 'NOT INSTALLED'}")
        if version is None:
            errors.append(f"missing distribution: {distribution}")
        elif version != expected_version:
            errors.append(
                f"{distribution} must be {expected_version}, found {version}"
            )

    pytest_version = module_version("pytest")
    print(f"pytest: {pytest_version or 'NOT INSTALLED'}")
    if pytest_version is None:
        errors.append("missing distribution: pytest")

    try:
        alpha_module = importlib.import_module("alphagenome_pytorch")
    except Exception as error:
        errors.append(f"cannot import alphagenome_pytorch: {error}")
    else:
        print(f"alphagenome_pytorch_path: {alpha_module.__file__}")

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
