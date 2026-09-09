#!/usr/bin/env python3
"""Download the frozen public eQTL sources and write immutable source audits."""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
import subprocess
import sys


REPO_ROOT = Path(__file__).resolve().parents[1]
SPEC = REPO_ROOT / "alphagenome_custom/metadata/v2/p20_eqtl_public_source_download_spec.json"


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--spec", type=Path, default=SPEC)
    return parser.parse_args()


def download(url: str, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    tmp = destination.with_suffix(destination.suffix + ".partial")
    command = [
        "curl", "-L", "--fail", "--retry", "8", "--retry-all-errors",
        "--continue-at", "-", "--output", str(tmp), url,
    ]
    subprocess.run(command, check=True)
    tmp.replace(destination)


def main() -> None:
    args = parse_args()
    spec_path = args.spec.resolve()
    spec = json.loads(spec_path.read_text(encoding="utf-8"))
    output = REPO_ROOT / spec["output_dir"]
    if output.exists() and (output / "p20_source_manifest.json").is_file():
        raise FileExistsError(f"P20 output is immutable: {output}")
    output.mkdir(parents=True, exist_ok=True)

    filenames = {
        "expression_counts": "GSE186719_Celegans_208strains_609samples_rawCounts.tsv.gz",
        "eqtl_truth_archive": "41467_2022_31208_MOESM4_ESM.zip",
        "cendr_isotype_vcf": "WI.20200815.hard-filter.isotype.vcf.gz",
    }
    records = []
    for key, url in spec["urls"].items():
        destination = output / filenames[key]
        print(f"downloading\t{key}\t{url}", flush=True)
        download(url, destination)
        records.append({
            "key": key,
            "url": url,
            "path": str(destination.relative_to(REPO_ROOT)),
            "bytes": destination.stat().st_size,
            "sha256": sha256(destination),
        })

    registry = REPO_ROOT / spec["local_candidate_registry"]
    collision = REPO_ROOT / spec["local_collision_audit"]
    manifest = {
        "schema_version": 1,
        "contract": spec["contract"],
        "created_at": utc_now(),
        "scope": spec["scope"],
        "spec": str(spec_path.relative_to(REPO_ROOT)),
        "spec_sha256": sha256(spec_path),
        "source_records": records,
        "candidate_registry": str(registry.relative_to(REPO_ROOT)),
        "candidate_registry_sha256": sha256(registry) if registry.is_file() else None,
        "collision_audit": str(collision.relative_to(REPO_ROOT)),
        "collision_audit_sha256": sha256(collision) if collision.is_file() else None,
        "source_status": "DOWNLOADED_PENDING_COHORT_SOURCE_EQUIVALENCE_AND_COORDINATE_AUDIT",
        "prohibited_claims": spec["prohibited_claims_until_follow_up_review"],
        "next_actions": [
            "Inspect the supplementary archive contents and freeze the eQTL truth table.",
            "Map expression samples to the 208 wild strains and verify the retained cohort.",
            "Audit CeNDR release/reference coordinates against WBcel235 on I-V before variant scoring.",
        ],
    }
    manifest_path = output / "p20_source_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"status": "downloaded", "manifest": str(manifest_path), "files": records}, indent=2))


if __name__ == "__main__":
    main()
