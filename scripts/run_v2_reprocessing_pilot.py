#!/usr/bin/env python3
"""Run the approved RNA-seq v2 source-read reprocessing pilot."""

from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import shutil
import subprocess
import sys
import time
from typing import Iterable


REPO_ROOT = Path(__file__).resolve().parents[1]
SAMPLE_MANIFEST = REPO_ROOT / "alphagenome_custom/metadata/v2/rna_seq_samples_v2.tsv"
EXPECTED_PILOT_RUNS = {
    "SRR7443583",
    "SRR941632",
    "SRR10882545",
    "SRR23049895",
    "SRR36719198",
}
COMMAND_LOG: list[dict[str, object]] = []
MONITOR_PATHS: list[Path] = []
DOH_ARGS = [
    "--doh-url",
    "https://dns.google/dns-query",
    "--resolve",
    "dns.google:443:8.8.8.8",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--manifest",
        default="alphagenome_custom/metadata/v2/p3_pilot_sources.tsv",
    )
    parser.add_argument("--threads", type=int, default=16)
    parser.add_argument(
        "--work-dir", default="shared/source_reads/v2/pilot_20260713"
    )
    parser.add_argument(
        "--output-dir",
        default="alphagenome_custom/tracks/rna_seq_v2_normalized_pilot",
    )
    parser.add_argument(
        "--star-index", default="shared/reference_indexes/WBcel235_STAR_2.7.11b"
    )
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sha256(path: Path, chunk_size: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def md5(path: Path, chunk_size: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.md5(usedforsecurity=False)
    with path.open("rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def require_tools(names: Iterable[str]) -> None:
    missing = [name for name in names if shutil.which(name) is None]
    if missing:
        raise RuntimeError(f"Missing required tools: {', '.join(missing)}")


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def tree_size(path: Path) -> int:
    if not path.exists():
        return 0
    if path.is_file():
        return path.stat().st_size
    total = 0
    for item in path.rglob("*"):
        try:
            if item.is_file():
                total += item.stat().st_size
        except OSError:
            continue
    return total


def monitored_bytes() -> int:
    return sum(tree_size(path) for path in MONITOR_PATHS)


def run(command: list[str], *, stdout=None) -> None:
    print("command\t" + " ".join(command), flush=True)
    started = time.monotonic()
    peak_bytes = monitored_bytes()
    process = subprocess.Popen(command, cwd=REPO_ROOT, stdout=stdout)
    while True:
        try:
            return_code = process.wait(timeout=2)
            peak_bytes = max(peak_bytes, monitored_bytes())
            break
        except subprocess.TimeoutExpired:
            peak_bytes = max(peak_bytes, monitored_bytes())
    COMMAND_LOG.append(
        {
            "argv": command,
            "duration_seconds": round(time.monotonic() - started, 6),
            "return_code": return_code,
            "observed_managed_bytes": peak_bytes,
        }
    )
    if return_code != 0:
        raise subprocess.CalledProcessError(return_code, command)


def tool_version(command: list[str]) -> str:
    completed = subprocess.run(
        command,
        cwd=REPO_ROOT,
        check=False,
        capture_output=True,
        text=True,
    )
    output = (completed.stdout + "\n" + completed.stderr).strip().splitlines()
    return output[0].strip() if output else f"exit_code={completed.returncode}"


def parse_star_log(path: Path) -> dict[str, str]:
    values = {}
    for line in path.read_text().splitlines():
        if "|" not in line:
            continue
        key, value = line.split("|", 1)
        values[key.strip()] = value.strip()
    required = [
        "Number of input reads",
        "Uniquely mapped reads number",
        "Uniquely mapped reads %",
    ]
    missing = [key for key in required if key not in values]
    if missing:
        raise RuntimeError(f"Missing STAR metrics in {path}: {missing}")
    return values


def download_fastq(url: str, expected_md5: str, expected_bytes: int, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not path.is_file() or path.stat().st_size != expected_bytes:
        tmp_path = path.with_suffix(path.suffix + ".part")
        if not (
            tmp_path.is_file()
            and tmp_path.stat().st_size == expected_bytes
            and md5(tmp_path) == expected_md5
        ):
            run(
                [
                    "curl",
                    "-fL",
                    "--retry",
                    "20",
                    "--retry-delay",
                    "2",
                    "--retry-all-errors",
                    "--continue-at",
                    "-",
                    *DOH_ARGS,
                    f"https://{url}",
                    "-o",
                    str(tmp_path),
                ]
            )
        if tmp_path.stat().st_size != expected_bytes or md5(tmp_path) != expected_md5:
            raise RuntimeError(f"Partial FASTQ integrity failure: {tmp_path}")
        tmp_path.replace(path)
    if path.stat().st_size != expected_bytes or md5(path) != expected_md5:
        raise RuntimeError(f"FASTQ integrity failure: {path}")


def normalize_bedgraph(input_path: Path, output_path: Path) -> tuple[float, float]:
    total = 0.0
    with input_path.open() as handle:
        for line in handle:
            chromosome, start, end, value = line.rstrip("\n").split("\t")
            numeric_value = float(value)
            if not math.isfinite(numeric_value) or numeric_value < 0:
                raise RuntimeError(
                    f"Invalid bedGraph value {value!r} at {chromosome}:{start}-{end}"
                )
            total += (int(end) - int(start)) * numeric_value
    if not math.isfinite(total) or total <= 0:
        raise RuntimeError(f"Non-positive bedGraph signal total: {input_path}")
    scale = 100_000_000.0 / total
    tmp_path = output_path.with_suffix(output_path.suffix + ".tmp")
    with input_path.open() as source, tmp_path.open("w") as destination:
        for line in source:
            chromosome, start, end, value = line.rstrip("\n").split("\t")
            normalized_value = float(value) * scale
            if not math.isfinite(normalized_value) or normalized_value < 0:
                raise RuntimeError(
                    f"Invalid normalized value at {chromosome}:{start}-{end}"
                )
            destination.write(
                f"{chromosome}\t{start}\t{end}\t{normalized_value:.9g}\n"
            )
    tmp_path.replace(output_path)
    return scale, total


def ensure_star_index(index_dir: Path, threads: int) -> None:
    if (index_dir / "Genome").is_file():
        return
    index_dir.mkdir(parents=True, exist_ok=True)
    run(
        [
            "STAR",
            "--runMode",
            "genomeGenerate",
            "--runThreadN",
            str(threads),
            "--genomeDir",
            str(index_dir),
            "--genomeFastaFiles",
            "alphagenome_custom/reference/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa",
            "--sjdbGTFfile",
            "alphagenome_custom/reference/Caenorhabditis_elegans.WBcel235.115.gtf",
            "--sjdbOverhang",
            "150",
            "--genomeSAindexNbases",
            "12",
        ]
    )


def validate_pilot_manifest(
    pilot_rows: list[dict[str, str]], sample_rows: list[dict[str, str]]
) -> dict[str, dict[str, str]]:
    if len(pilot_rows) != 5:
        raise RuntimeError(f"P3A manifest must contain exactly five runs, got {len(pilot_rows)}")
    if len({row["run_accession"] for row in pilot_rows}) != 5:
        raise RuntimeError("P3A manifest run accessions are not unique")
    if {row["run_accession"] for row in pilot_rows} != EXPECTED_PILOT_RUNS:
        raise RuntimeError("P3A manifest does not match the approved five-run set")
    sample_by_run = {row["run_accession"]: row for row in sample_rows}
    missing = [
        row["run_accession"]
        for row in pilot_rows
        if row["run_accession"] not in sample_by_run
    ]
    if missing:
        raise RuntimeError(f"P3A runs missing from v2 sample manifest: {missing}")
    for row in pilot_rows:
        urls = row["fastq_ftp"].split(";")
        checksums = row["fastq_md5"].split(";")
        sizes = row["fastq_bytes"].split(";")
        expected_files = 2 if row["library_layout"] == "PAIRED" else 1
        if not (len(urls) == len(checksums) == len(sizes) == expected_files):
            raise RuntimeError(
                f"FASTQ fields do not match {row['library_layout']} for {row['run_accession']}"
            )
    return sample_by_run


def read_reference_chromosomes() -> dict[str, int]:
    chromosomes = {}
    with (REPO_ROOT / "alphagenome_custom/reference/genome.fa.fai").open() as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            chromosomes[fields[0]] = int(fields[1])
    return chromosomes


def validate_output_bigwig(path: Path, expected_chromosomes: dict[str, int]) -> dict[str, float]:
    import pyBigWig

    with pyBigWig.open(str(path)) as bigwig:
        if not bigwig.isBigWig():
            raise RuntimeError(f"Not a bigWig: {path}")
        chromosomes = bigwig.chroms()
        if chromosomes != expected_chromosomes:
            raise RuntimeError(f"Chromosome mismatch in {path}: {chromosomes}")
        header = bigwig.header()
    numeric = {
        key: float(header[key]) for key in ("minVal", "maxVal", "sumData", "sumSquared")
    }
    if not all(math.isfinite(value) for value in numeric.values()):
        raise RuntimeError(f"Non-finite bigWig header statistics: {path}")
    if numeric["minVal"] < 0 or numeric["maxVal"] < 0:
        raise RuntimeError(f"Negative bigWig signal: {path}")
    relative_error = abs(numeric["sumData"] - 100_000_000.0) / 100_000_000.0
    if relative_error > 1e-5:
        raise RuntimeError(
            f"Unexpected normalized total for {path}: {numeric['sumData']}"
        )
    numeric["relative_total_error"] = relative_error
    return numeric


def main() -> None:
    args = parse_args()
    require_tools(["STAR", "samtools", "bedtools", "bedGraphToBigWig", "curl"])
    manifest_path = REPO_ROOT / args.manifest
    work_dir = REPO_ROOT / args.work_dir
    output_dir = REPO_ROOT / args.output_dir
    index_dir = REPO_ROOT / args.star_index
    work_dir.mkdir(parents=True, exist_ok=True)
    output_dir.mkdir(parents=True, exist_ok=True)
    MONITOR_PATHS.extend([work_dir, output_dir, index_dir])
    run_started_at = utc_now()
    run_started = time.monotonic()
    pilot_rows = read_tsv(manifest_path)
    sample_rows = read_tsv(SAMPLE_MANIFEST)
    sample_by_run = validate_pilot_manifest(pilot_rows, sample_rows)
    expected_chromosomes = read_reference_chromosomes()
    raw_before = {
        row["run_accession"]: sha256(REPO_ROOT / sample_by_run[row["run_accession"]]["local_path"])
        for row in pilot_rows
    }
    raw_mismatches = [
        accession
        for accession, digest in raw_before.items()
        if digest != sample_by_run[accession]["sha256"]
    ]
    if raw_mismatches:
        raise RuntimeError(f"Provided bigWig SHA-256 mismatch before P3A: {raw_mismatches}")
    ensure_star_index(index_dir, args.threads)

    chrom_sizes = work_dir / "WBcel235.chrom.sizes"
    if not chrom_sizes.is_file():
        with (
            REPO_ROOT / "alphagenome_custom/reference/genome.fa.fai"
        ).open() as source, chrom_sizes.open("w") as destination:
            for line in source:
                fields = line.split("\t")
                destination.write(f"{fields[0]}\t{fields[1]}\n")

    report_rows = []
    tool_versions = {
        "STAR": tool_version(["STAR", "--version"]),
        "samtools": tool_version(["samtools", "--version"]),
        "bedtools": tool_version(["bedtools", "--version"]),
        "bedGraphToBigWig": tool_version(["bedGraphToBigWig"]),
    }
    for row in pilot_rows:
        sample_started = time.monotonic()
        accession = row["run_accession"]
        sample_dir = work_dir / accession
        sample_dir.mkdir(parents=True, exist_ok=True)
        urls = row["fastq_ftp"].split(";")
        checksums = row["fastq_md5"].split(";")
        sizes = [int(value) for value in row["fastq_bytes"].split(";")]
        fastq_paths = []
        for url, checksum, size in zip(urls, checksums, sizes, strict=True):
            path = sample_dir / Path(url).name
            download_fastq(url, checksum, size, path)
            fastq_paths.append(path)
        fastq_sha256 = [sha256(path) for path in fastq_paths]

        star_prefix = sample_dir / "star_"
        aligned_bam = sample_dir / "star_Aligned.sortedByCoord.out.bam"
        if not aligned_bam.is_file():
            run(
                [
                    "STAR",
                    "--runThreadN",
                    str(args.threads),
                    "--genomeDir",
                    str(index_dir),
                    "--readFilesIn",
                    *map(str, fastq_paths),
                    "--readFilesCommand",
                    "zcat",
                    "--outFileNamePrefix",
                    str(star_prefix),
                    "--outSAMtype",
                    "BAM",
                    "SortedByCoordinate",
                    "--outFilterMultimapNmax",
                    "1",
                    "--outFilterMismatchNoverReadLmax",
                    "0.04",
                    "--alignSJDBoverhangMin",
                    "1",
                    "--limitBAMsortRAM",
                    "8000000000",
                ]
            )

        star_metrics = parse_star_log(sample_dir / "star_Log.final.out")
        primary_bam = sample_dir / "primary_unique.bam"
        if not primary_bam.is_file():
            run(
                [
                    "samtools",
                    "view",
                    "-@",
                    str(max(args.threads // 2, 1)),
                    "-b",
                    "-F",
                    "2308",
                    "-q",
                    "1",
                    "-o",
                    str(primary_bam),
                    str(aligned_bam),
                ]
            )
        primary_bai = primary_bam.with_suffix(primary_bam.suffix + ".bai")
        if not primary_bai.is_file():
            run(["samtools", "index", "-@", "4", str(primary_bam)])
        run(["samtools", "quickcheck", "-v", str(primary_bam)])

        raw_bedgraph = sample_dir / "coverage.raw.bedGraph"
        with raw_bedgraph.open("w") as handle:
            run(
                ["bedtools", "genomecov", "-split", "-bg", "-ibam", str(primary_bam)],
                stdout=handle,
            )
        normalized_bedgraph = sample_dir / "coverage.normalized.bedGraph"
        scale, raw_total = normalize_bedgraph(raw_bedgraph, normalized_bedgraph)
        bigwig_path = output_dir / f"{accession}.bw"
        tmp_bigwig = bigwig_path.with_suffix(".bw.tmp")
        run(
            [
                "bedGraphToBigWig",
                str(normalized_bedgraph),
                str(chrom_sizes),
                str(tmp_bigwig),
            ]
        )
        tmp_bigwig.replace(bigwig_path)
        output_stats = validate_output_bigwig(bigwig_path, expected_chromosomes)
        raw_after = sha256(
            REPO_ROOT / sample_by_run[accession]["local_path"]
        )
        report_rows.append(
            {
                "run_accession": accession,
                "library_layout": row["library_layout"],
                "library_selection": row["library_selection"],
                "source_fastq_paths": ";".join(
                    str(path.relative_to(REPO_ROOT)) for path in fastq_paths
                ),
                "source_fastq_md5": ";".join(checksums),
                "source_fastq_sha256": ";".join(fastq_sha256),
                "source_fastq_bytes": sum(sizes),
                "provided_bigwig_path": sample_by_run[accession]["local_path"],
                "provided_bigwig_expected_sha256": sample_by_run[accession]["sha256"],
                "provided_bigwig_sha256_before": raw_before[accession],
                "provided_bigwig_sha256_after": raw_after,
                "provided_bigwig_unchanged": str(raw_before[accession] == raw_after),
                "star_input_reads": star_metrics["Number of input reads"],
                "star_uniquely_mapped_reads": star_metrics[
                    "Uniquely mapped reads number"
                ],
                "star_uniquely_mapped_percent": star_metrics[
                    "Uniquely mapped reads %"
                ],
                "primary_bam_path": str(primary_bam.relative_to(REPO_ROOT)),
                "primary_bam_size_bytes": primary_bam.stat().st_size,
                "primary_bam_sha256": sha256(primary_bam),
                "coverage_policy": "primary_unique_spliced_unstranded",
                "normalization_formula": "100000000/sum((end-start)*raw_coverage)",
                "raw_bedgraph_total": f"{raw_total:.12g}",
                "bedgraph_scale_to_1e8": f"{scale:.12g}",
                "output_path": str(bigwig_path.relative_to(REPO_ROOT)),
                "output_size_bytes": bigwig_path.stat().st_size,
                "output_sha256": sha256(bigwig_path),
                "output_strand": ".",
                "output_min": f"{output_stats['minVal']:.12g}",
                "output_max": f"{output_stats['maxVal']:.12g}",
                "output_total_signal": f"{output_stats['sumData']:.12g}",
                "output_relative_total_error": f"{output_stats['relative_total_error']:.12g}",
                "sample_elapsed_seconds": f"{time.monotonic() - sample_started:.6f}",
                "sample_work_bytes_final": tree_size(sample_dir),
            }
        )

    report_path = output_dir / "pilot_outputs.tsv"
    with report_path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(report_rows[0]),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(report_rows)
    partial_files = [
        str(path.relative_to(REPO_ROOT))
        for root in (work_dir, output_dir, index_dir)
        for pattern in ("*.part", "*.tmp")
        for path in root.rglob(pattern)
    ]
    if partial_files:
        raise RuntimeError(f"Partial P3A files remain: {partial_files}")
    commands_path = output_dir / "pilot_commands.json"
    commands_path.write_text(
        json.dumps(COMMAND_LOG, indent=2, sort_keys=True) + "\n"
    )
    raw_after = {
        accession: sha256(REPO_ROOT / sample_by_run[accession]["local_path"])
        for accession in raw_before
    }
    (output_dir / "pilot_run.json").write_text(
        json.dumps(
            {
                "schema_version": 1,
                "phase": "P3A",
                "scope": "five_run_reprocessing_pilot_only",
                "formal_v2_outputs": False,
                "manifest": args.manifest,
                "manifest_sha256": sha256(manifest_path),
                "sample_manifest_sha256": sha256(SAMPLE_MANIFEST),
                "threads": args.threads,
                "work_dir": args.work_dir,
                "output_dir": args.output_dir,
                "star_index": args.star_index,
                "outputs": len(report_rows),
                "source_fastq_bytes": sum(
                    int(value)
                    for row in pilot_rows
                    for value in row["fastq_bytes"].split(";")
                ),
                "normalization_target_total": 100_000_000,
                "coverage_policy": "primary_unique_spliced_unstranded",
                "pipeline_parameters": {
                    "star_out_filter_multimap_nmax": 1,
                    "star_out_filter_mismatch_nover_read_lmax": 0.04,
                    "star_align_sjdb_overhang_min": 1,
                    "samtools_exclude_flags": 2308,
                    "samtools_min_mapq": 1,
                    "bedtools_genomecov_split": True,
                },
                "tool_versions": tool_versions,
                "commands_path": str(commands_path.relative_to(REPO_ROOT)),
                "commands_sha256": sha256(commands_path),
                "command_count": len(COMMAND_LOG),
                "observed_managed_peak_bytes": max(
                    (int(row["observed_managed_bytes"]) for row in COMMAND_LOG),
                    default=monitored_bytes(),
                ),
                "raw_bigwig_sha256_before": raw_before,
                "raw_bigwig_sha256_after": raw_after,
                "raw_bigwigs_unchanged": raw_before == raw_after,
                "partial_files": partial_files,
                "intermediates_retained_for_r3a_diagnosis": True,
                "started_at": run_started_at,
                "completed_at": utc_now(),
                "elapsed_seconds": round(time.monotonic() - run_started, 6),
            },
            indent=2,
            sort_keys=True,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
