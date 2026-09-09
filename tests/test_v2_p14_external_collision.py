"""Regression tests for the P14 metadata-only collision screen."""

from __future__ import annotations

import csv
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = REPO_ROOT / "scripts" / "audit_v2_p14_external_collision.py"
SAMPLES = REPO_ROOT / "alphagenome_custom" / "metadata" / "v2" / "rna_seq_samples_v2.tsv"
CONTEXT = REPO_ROOT / "alphagenome_custom" / "metadata" / "v2" / "rna_seq_sample_context_v2.tsv"
MEMBERS = REPO_ROOT / "alphagenome_custom" / "metadata" / "v2" / "rna_seq_group_members_v2_final.tsv"

FIELDS = [
    "run_accession",
    "experiment_accession",
    "biosample_accession",
    "sample_alias",
    "sra_study_accession",
    "bioproject_accession",
    "geo_accession",
    "fastq_md5",
    "processed_sha256",
    "source_equivalence_audit_status",
]


def read_first(path: Path) -> dict[str, str]:
    with path.open(newline="", encoding="utf-8") as handle:
        return next(csv.DictReader(handle, delimiter="\t"))


def read_for_run(path: Path, run_accession: str) -> dict[str, str]:
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["run_accession"] == run_accession:
                return row
    raise AssertionError(f"{run_accession} is absent from {path}")


class P14ExternalCollisionTest(unittest.TestCase):
    def test_existing_v2_run_is_rejected(self) -> None:
        member = read_first(MEMBERS)
        sample = read_for_run(SAMPLES, member["run_accession"])
        context = read_for_run(CONTEXT, member["run_accession"])
        self.assertEqual(sample["run_accession"], context["run_accession"])
        self.assertEqual(sample["run_accession"], member["run_accession"])
        row = {
            "run_accession": sample["run_accession"],
            "experiment_accession": sample["experiment_accession"],
            "biosample_accession": sample["biosample_accession"],
            "sample_alias": context["sample_alias"],
            "sra_study_accession": sample["sra_study_accession"],
            "bioproject_accession": sample["bioproject_accession"],
            "geo_accession": sample["geo_accession"],
            "fastq_md5": sample["fastq_md5"],
            "processed_sha256": member["sha256"],
            "source_equivalence_audit_status": "cleared_no_equivalent_source",
        }
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            registry = root / "candidate.tsv"
            with registry.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t")
                writer.writeheader()
                writer.writerow(row)
            output = root / "audit"
            completed = subprocess.run(
                [sys.executable, str(SCRIPT), "--external-registry", str(registry), "--output-dir", str(output)],
                check=True,
                capture_output=True,
                text=True,
            )
            result = json.loads(completed.stdout)
            self.assertEqual(result["counts"], {"REJECT": 1})
            with (output / "p14_collision_audit.tsv").open(newline="", encoding="utf-8") as handle:
                report = next(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(report["overall_status"], "REJECT")
            self.assertIn("run_accession", report["collision_fields"])

    def test_raw_clear_candidate_does_not_require_unavailable_processed_checksum(self) -> None:
        row = {
            "run_accession": "SRR_EXTERNAL_RAW_CLEAR",
            "experiment_accession": "SRX_EXTERNAL_RAW_CLEAR",
            "biosample_accession": "SAMN_EXTERNAL_RAW_CLEAR",
            "sample_alias": "external-raw-clear",
            "sra_study_accession": "SRP_EXTERNAL_RAW_CLEAR",
            "bioproject_accession": "PRJNA_EXTERNAL_RAW_CLEAR",
            "geo_accession": "GSE_EXTERNAL_RAW_CLEAR",
            "fastq_md5": "0123456789abcdef0123456789abcdef",
            "processed_sha256": "",
            "source_equivalence_audit_status": "cleared_no_equivalent_source",
        }
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            registry = root / "candidate.tsv"
            with registry.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t")
                writer.writeheader()
                writer.writerow(row)
            output = root / "audit"
            subprocess.run([sys.executable, str(SCRIPT), "--external-registry", str(registry), "--output-dir", str(output)], check=True, capture_output=True, text=True)
            with (output / "p14_collision_audit.tsv").open(newline="", encoding="utf-8") as handle:
                report = next(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(report["overall_status"], "CLEAR_FOR_MANUAL_GATE0_REVIEW")
            self.assertEqual(report["raw_checksum_status"], "CLEAR_MACHINE_SCREEN")
            self.assertEqual(report["processed_checksum_status"], "PENDING_POSTPROCESSING_CHECKSUM")

    def test_component_fastq_checksum_collision_is_rejected(self) -> None:
        sample = read_first(SAMPLES)
        component = sample["fastq_md5"].split(";")[0]
        row = {
            "run_accession": "SRR_EXTERNAL_COMPONENT_COLLISION",
            "experiment_accession": "SRX_EXTERNAL_COMPONENT_COLLISION",
            "biosample_accession": "SAMN_EXTERNAL_COMPONENT_COLLISION",
            "sample_alias": "external-component-collision",
            "sra_study_accession": "SRP_EXTERNAL_COMPONENT_COLLISION",
            "bioproject_accession": "PRJNA_EXTERNAL_COMPONENT_COLLISION",
            "geo_accession": "GSE_EXTERNAL_COMPONENT_COLLISION",
            "fastq_md5": f"ffffffffffffffffffffffffffffffff;{component}",
            "processed_sha256": "",
            "source_equivalence_audit_status": "cleared_no_equivalent_source",
        }
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            registry = root / "candidate.tsv"
            with registry.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t")
                writer.writeheader()
                writer.writerow(row)
            output = root / "audit"
            subprocess.run([sys.executable, str(SCRIPT), "--external-registry", str(registry), "--output-dir", str(output)], check=True, capture_output=True, text=True)
            with (output / "p14_collision_audit.tsv").open(newline="", encoding="utf-8") as handle:
                report = next(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(report["overall_status"], "REJECT")
            self.assertEqual(report["raw_checksum_status"], "REJECT_COLLISION")
            self.assertIn("fastq_md5", report["collision_fields"])
