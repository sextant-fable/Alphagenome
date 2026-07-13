#!/usr/bin/env python3
"""Select idle HY-GPU devices 2/3 without touching unrelated processes."""

from __future__ import annotations

import argparse
import csv
import io
import json
import subprocess
from typing import Any


PREFERRED_GPUS = (2, 3)


def parse_gpu_csv(text: str) -> list[dict[str, Any]]:
    rows = []
    for fields in csv.reader(io.StringIO(text)):
        if not fields:
            continue
        index, uuid, name, total, used, free, utilization = (
            field.strip() for field in fields
        )
        rows.append(
            {
                "index": int(index),
                "uuid": uuid,
                "name": name,
                "memory_total_mib": int(total),
                "memory_used_mib": int(used),
                "memory_free_mib": int(free),
                "utilization_percent": int(utilization),
            }
        )
    return rows


def parse_process_csv(text: str) -> list[dict[str, Any]]:
    rows = []
    for fields in csv.reader(io.StringIO(text)):
        if not fields:
            continue
        uuid, pid, name, memory = (field.strip() for field in fields)
        rows.append(
            {
                "gpu_uuid": uuid,
                "pid": int(pid),
                "process_name": name,
                "used_memory_mib": int(memory),
            }
        )
    return rows


def snapshot() -> dict[str, Any]:
    gpu_command = [
        "nvidia-smi",
        "--query-gpu=index,uuid,name,memory.total,memory.used,memory.free,utilization.gpu",
        "--format=csv,noheader,nounits",
    ]
    process_command = [
        "nvidia-smi",
        "--query-compute-apps=gpu_uuid,pid,process_name,used_memory",
        "--format=csv,noheader,nounits",
    ]
    gpu_output = subprocess.run(
        gpu_command, check=True, capture_output=True, text=True
    ).stdout
    process_output = subprocess.run(
        process_command, check=True, capture_output=True, text=True
    ).stdout
    return {
        "gpus": parse_gpu_csv(gpu_output),
        "processes": parse_process_csv(process_output),
        "commands": [gpu_command, process_command],
    }


def select_available(
    resource_snapshot: dict[str, Any],
    *,
    count: int = 1,
    minimum_free_mib: int = 70_000,
) -> list[int]:
    if count not in {1, 2}:
        raise ValueError("count must be 1 or 2")
    processes_by_uuid = {
        row["gpu_uuid"] for row in resource_snapshot.get("processes", [])
    }
    by_index = {row["index"]: row for row in resource_snapshot["gpus"]}
    available = []
    for index in PREFERRED_GPUS:
        row = by_index.get(index)
        if row is None:
            continue
        if (
            row["uuid"] not in processes_by_uuid
            and row["memory_free_mib"] >= minimum_free_mib
            and row["utilization_percent"] <= 10
        ):
            available.append(index)
    if len(available) < count:
        raise RuntimeError(
            f"Need {count} idle GPU(s) among {PREFERRED_GPUS} with "
            f">={minimum_free_mib} MiB free; available={available}"
        )
    return available[:count]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--count", type=int, default=1, choices=[1, 2])
    parser.add_argument("--minimum-free-mib", type=int, default=70_000)
    args = parser.parse_args()
    resources = snapshot()
    resources["selected"] = select_available(
        resources, count=args.count, minimum_free_mib=args.minimum_free_mib
    )
    print(json.dumps(resources, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
