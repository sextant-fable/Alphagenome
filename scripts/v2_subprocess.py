"""Build repository-local Python subprocess commands consistently."""

from __future__ import annotations

import sys


def module_command(module: str, *arguments: str) -> list[str]:
    if not module or module.startswith("scripts."):
        raise ValueError("module must be a non-empty name relative to scripts")
    return [sys.executable, "-m", f"scripts.{module}", *arguments]
