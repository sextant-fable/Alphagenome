#!/usr/bin/env python3
"""Uniformly reprocess the P18 sources while keeping chromosome X inaccessible.

The output comprises only chromosome I--V coverage. All source FASTQs must
already have passed the P18 checksum-verified transfer. This creates source-
level coverage for a later represented-context benchmark; it does not score a
model or make an external-validation claim.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import shutil
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pyBigWig


REPO_ROOT = Path(__file__).resolve().parents[1]
P18_SCOPE = "p18_external_rna_download_and_reprocessing"
EXPECTED_RUNS = {"SRR18463404", "SRR18463405", "SRR18463406", "SRR2005820"}
IV_CHROMOSOMES = ("I", "II", "III", "IV", "V")
STAR_SORT_RAM_BYTES = 128_000_000_000


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--transport-manifest", type=Path, required=True)
    parser.add_argument("--download-audit-dir", type=Path, required=True)
    parser.add_argument("--fastq-root", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--expected-run", action="append", default=[])
    parser.add_argument("--controller-phase", default="P17")
    parser.add_argument("--approval-scope", default=P18_SCOPE)
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sha256(path: Path) -> str:
    import hashlib

    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def require_scoped_approval(phase: str, scope: str) -> None:
    from scripts import v2_phase_controller as controller

    state = controller.load_state()
    controller.validate_approval_scope("G1", scope, phase)
    if state.get("current_phase") != phase or state.get("status") != "RUNNING":
        raise RuntimeError(f"External reprocessing requires the active {phase} RUNNING period")
    if not controller.approval_satisfies(state, "G1", scope):
        raise RuntimeError(f"External reprocessing requires G1:{scope}")


def run(command: list[str], *, cwd: Path, stdout_path: Path | None = None) -> None:
    import subprocess

    print("command\t" + " ".join(command), flush=True)
    if stdout_path is None:
        subprocess.run(command, cwd=cwd, check=True)
    else:
        with stdout_path.open("w", encoding="utf-8") as handle:
            subprocess.run(command, cwd=cwd, check=True, stdout=handle)


def parse_star_log(path: Path) -> dict[str, str]:
    values: dict[str, str] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if "|" in line:
            key, value = line.split("|", 1)
            values[key.strip()] = value.strip()
    required = ("Number of input reads", "Uniquely mapped reads number", "Uniquely mapped reads %")
    if missing := [key for key in required if key not in values]:
        raise RuntimeError(f"STAR log lacks metrics {missing}: {path}")
    return values


def normalize_bedgraph(source: Path, destination: Path) -> tuple[float, float]:
    total = 0.0
    with source.open(encoding="utf-8") as handle:
        for line in handle:
            chromosome, start, end, value = line.rstrip("\n").split("\t")
            if chromosome not in IV_CHROMOSOMES:
                raise RuntimeError(f"Chromosome outside I-V appeared in external coverage: {chromosome}")
            numeric = float(value)
            if not math.isfinite(numeric) or numeric < 0:
                raise RuntimeError(f"Invalid coverage in {source}: {line.rstrip()}")
            total += (int(end) - int(start)) * numeric
    if not math.isfinite(total) or total <= 0:
        raise RuntimeError(f"Non-positive I-V coverage total: {source}")
    scale = 100_000_000.0 / total
    temporary = destination.with_suffix(destination.suffix + ".tmp")
    with source.open(encoding="utf-8") as input_handle, temporary.open("w", encoding="utf-8") as output_handle:
        for line in input_handle:
            chromosome, start, end, value = line.rstrip("\n").split("\t")
            output_handle.write(f"{chromosome}\t{start}\t{end}\t{float(value) * scale:.9g}\n")
    temporary.replace(destination)
    return scale, total


def validate_iv_bigwig(path: Path, chromosomes: dict[str, int]) -> dict[str, float]:
    with pyBigWig.open(str(path)) as bigwig:
        if not bigwig.isBigWig() or bigwig.chroms() != chromosomes:
            raise RuntimeError(f"Invalid I-V-only external bigWig: {path}")
        header = bigwig.header()
    values = {key: float(header[key]) for key in ("minVal", "maxVal", "sumData", "sumSquared")}
    if not all(math.isfinite(value) for value in values.values()) or values["minVal"] < 0 or values["maxVal"] < 0:
        raise RuntimeError(f"Invalid external bigWig statistics: {path}")
    relative_error = abs(values["sumData"] - 100_000_000.0) / 100_000_000.0
    if relative_error > 1e-5:
        raise RuntimeError(f"External bigWig I-V total is not 1e8: {path}")
    values["relative_total_error"] = relative_error
    return values


def load_download_audit(audit_dir: Path, run: str, manifest_sha: str) -> dict[str, Any]:
    path = audit_dir / f"{run}.json"
    if not path.is_file():
        raise RuntimeError(f"Missing P18 download audit for {run}")
    audit = json.loads(path.read_text(encoding="utf-8"))
    if audit.get("transport_manifest_sha256") != manifest_sha:
        raise RuntimeError(f"Download audit transport-manifest hash mismatch for {run}")
    paths = [REPO_ROOT / relative for relative in audit.get("local_fastq_paths", [])]
    if not paths or any(not path.is_file() for path in paths):
        raise RuntimeError(f"Verified FASTQ files are unavailable for {run}")
    return audit


def main() -> None:
    args = parse_args()
    if args.threads < 1 or args.threads > 16:
        raise ValueError("P18 I-V reprocessing threads must be between 1 and 16")
    transport = args.transport_manifest.resolve()
    download_audit_dir = args.download_audit_dir.resolve()
    fastq_root = args.fastq_root.resolve()
    work_dir = args.work_dir.resolve()
    output_dir = args.output_dir.resolve()
    if not transport.is_file() or not download_audit_dir.is_dir():
        raise FileNotFoundError("P18 transport manifest or completed download audit directory is missing")
    for path in (fastq_root, work_dir, output_dir):
        path.relative_to(REPO_ROOT.resolve())
    expected_runs = set(args.expected_run) if args.expected_run else EXPECTED_RUNS
    require_scoped_approval(args.controller_phase, args.approval_scope)
    rows = read_tsv(transport)
    if {row.get("run_accession", "") for row in rows} != expected_runs or len(rows) != len(expected_runs):
        raise RuntimeError("External transport manifest does not contain exactly the registered runs")
    if any(row.get("source_equivalence_audit_status") != "cleared_no_equivalent_source" for row in rows):
        raise RuntimeError("P18 source-equivalence attestation is incomplete")
    required_tools = ("STAR", "samtools", "bedtools", "bedGraphToBigWig")
    if missing := [name for name in required_tools if shutil.which(name) is None]:
        raise RuntimeError(f"P18 external reprocessing lacks tools: {missing}")
    index_dir = REPO_ROOT / "shared/reference_indexes/WBcel235_STAR_2.7.11b"
    if not (index_dir / "Genome").is_file():
        raise RuntimeError("P18 requires the pre-existing WBcel235 STAR index; it will not regenerate an index")
    reference_sizes: dict[str, int] = {}
    with (REPO_ROOT / "alphagenome_custom/reference/genome.fa.fai").open(encoding="utf-8") as handle:
        for line in handle:
            chromosome, length, *_ = line.rstrip("\n").split("\t")
            if chromosome in IV_CHROMOSOMES:
                reference_sizes[chromosome] = int(length)
    if tuple(reference_sizes) != IV_CHROMOSOMES:
        raise RuntimeError("WBcel235 I-V reference sizes are incomplete or unordered")
    work_dir.mkdir(parents=True, exist_ok=True)
    output_dir.mkdir(parents=True, exist_ok=True)
    chrom_sizes = work_dir / "WBcel235_I_to_V.chrom.sizes"
    chrom_sizes.write_text("".join(f"{chromosome}\t{length}\n" for chromosome, length in reference_sizes.items()), encoding="utf-8")
    transport_sha = sha256(transport)
    records: list[dict[str, Any]] = []
    for row in sorted(rows, key=lambda value: value["run_accession"]):
        require_scoped_approval(args.controller_phase, args.approval_scope)
        run_id = row["run_accession"]
        audit_path = output_dir / "sample_audits" / f"{run_id}.json"
        output_path = output_dir / f"{run_id}.bw"
        if audit_path.is_file() and output_path.is_file():
            existing = json.loads(audit_path.read_text(encoding="utf-8"))
            if existing.get("transport_manifest_sha256") == transport_sha and existing.get("output_sha256") == sha256(output_path):
                validate_iv_bigwig(output_path, reference_sizes)
                records.append(existing)
                continue
        download_audit = load_download_audit(download_audit_dir, run_id, transport_sha)
        fastq_paths = [REPO_ROOT / relative for relative in download_audit["local_fastq_paths"]]
        sample_dir = work_dir / run_id
        if sample_dir.exists():
            shutil.rmtree(sample_dir)
        sample_dir.mkdir(parents=True)
        started = utc_now()
        try:
            star_prefix = sample_dir / "star_"
            run([
                "STAR", "--runThreadN", str(args.threads), "--genomeDir", str(index_dir),
                "--readFilesIn", *map(str, fastq_paths), "--readFilesCommand", "zcat",
                "--outFileNamePrefix", str(star_prefix), "--outSAMtype", "BAM", "SortedByCoordinate",
                "--outFilterMultimapNmax", "1", "--outFilterMismatchNoverReadLmax", "0.04",
                "--alignSJDBoverhangMin", "1", "--limitBAMsortRAM", str(STAR_SORT_RAM_BYTES),
            ], cwd=REPO_ROOT)
            metrics = parse_star_log(sample_dir / "star_Log.final.out")
            if int(metrics["Number of input reads"].replace(",", "")) != int(row["read_count"]):
                raise RuntimeError(f"STAR input read count differs from frozen ENA metadata for {run_id}")
            star_bam = sample_dir / "star_Aligned.sortedByCoord.out.bam"
            run([
                "samtools", "index", "-@", str(max(args.threads // 2, 1)), str(star_bam),
            ], cwd=REPO_ROOT)
            iv_bam = sample_dir / "I_to_V.primary_unique.bam"
            run([
                "samtools", "view", "-@", str(max(args.threads // 2, 1)), "-b", "-F", "2308", "-q", "1", "-o", str(iv_bam),
                str(star_bam), *IV_CHROMOSOMES,
            ], cwd=REPO_ROOT)
            run(["samtools", "quickcheck", "-v", str(iv_bam)], cwd=REPO_ROOT)
            raw_bedgraph = sample_dir / "I_to_V.coverage.raw.bedGraph"
            run(["bedtools", "genomecov", "-split", "-bg", "-ibam", str(iv_bam)], cwd=REPO_ROOT, stdout_path=raw_bedgraph)
            normalized_bedgraph = sample_dir / "I_to_V.coverage.normalized.bedGraph"
            scale, raw_total = normalize_bedgraph(raw_bedgraph, normalized_bedgraph)
            temporary_bigwig = output_path.with_suffix(".bw.tmp")
            run(["bedGraphToBigWig", str(normalized_bedgraph), str(chrom_sizes), str(temporary_bigwig)], cwd=REPO_ROOT)
            temporary_bigwig.replace(output_path)
            statistics = validate_iv_bigwig(output_path, reference_sizes)
            record = {
                "schema_version": 1,
                "run_accession": run_id,
                "proposed_group_id": row["proposed_group_id"],
                "transport_manifest": str(transport.relative_to(REPO_ROOT)),
                "transport_manifest_sha256": transport_sha,
                "download_audit": str((download_audit_dir / f"{run_id}.json").relative_to(REPO_ROOT)),
                "reference": "WBcel235",
                "coverage_chromosomes": list(IV_CHROMOSOMES),
                "chromosome_x_access": "prohibited; no X coverage was emitted or read",
                "coverage_policy": "primary_unique_spliced_unstranded_I_to_V_only",
                "normalization_formula": "100000000/sum_I_to_V((end-start)*raw_coverage)",
                "raw_I_to_V_bedgraph_total": raw_total,
                "scale_to_1e8_I_to_V": scale,
                "star_input_reads": metrics["Number of input reads"],
                "star_uniquely_mapped_reads": metrics["Uniquely mapped reads number"],
                "star_uniquely_mapped_percent": metrics["Uniquely mapped reads %"],
                "output_path": str(output_path.relative_to(REPO_ROOT)),
                "output_sha256": sha256(output_path),
                "output_total_signal": statistics["sumData"],
                "output_relative_total_error": statistics["relative_total_error"],
                "started_at": started,
                "completed_at": utc_now(),
                "claim_boundary": "I-V-only external source coverage for later represented-context comparison. No model inference or external performance result is generated here.",
            }
            audit_path.parent.mkdir(parents=True, exist_ok=True)
            temporary = audit_path.with_suffix(".json.tmp")
            temporary.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n", encoding="utf-8")
            temporary.replace(audit_path)
            records.append(record)
        finally:
            if sample_dir.exists():
                shutil.rmtree(sample_dir)
        progress = {
            "schema_version": 1,
            "scope": f"{args.controller_phase} external I-V-only uniform reprocessing",
            "expected_runs": sorted(expected_runs),
            "completed_runs": len(records),
            "completed_accessions": [item["run_accession"] for item in records],
            "updated_at": utc_now(),
        }
        (output_dir / "p18_reprocessing_progress.json").write_text(json.dumps(progress, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        print(json.dumps({"run_accession": run_id, "status": "reprocessed_iv_only", "completed": len(records), "expected": len(rows)}, sort_keys=True), flush=True)
    summary = {
        "schema_version": 1,
        "scope": f"{args.controller_phase} external I-V-only uniform reprocessing",
        "transport_manifest": str(transport.relative_to(REPO_ROOT)),
        "transport_manifest_sha256": transport_sha,
        "completed_runs": len(records),
        "completed_accessions": sorted(item["run_accession"] for item in records),
        "completed_at": utc_now(),
        "next_required_gate": "Post-processing checksum/source-equivalence audit and a frozen external evaluation specification are required before model scoring.",
    }
    (output_dir / "p18_reprocessing_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"status": "complete", "runs": len(records), "output_dir": str(output_dir.relative_to(REPO_ROOT))}, sort_keys=True))


if __name__ == "__main__":
    main()
