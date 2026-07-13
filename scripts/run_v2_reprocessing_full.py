#!/usr/bin/env python3
"""Stream all 482 verified RNA-seq runs into comparable WBcel235 bigWigs."""

from __future__ import annotations

import argparse
from concurrent.futures import Future, ThreadPoolExecutor, as_completed
import csv
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import threading
import time
from typing import Any

import pyBigWig


REPO_ROOT = Path(__file__).resolve().parents[1]
SAMPLE_MANIFEST = REPO_ROOT / "alphagenome_custom/metadata/v2/rna_seq_samples_v2.tsv"
SOURCE_MANIFEST = REPO_ROOT / "alphagenome_custom/metadata/v2/p3_full_sources.tsv"
SRA_MANIFEST = REPO_ROOT / "alphagenome_custom/metadata/v2/p3_ncbi_sra_sources.tsv"
LEDGER_PATH = REPO_ROOT / "alphagenome_custom/metadata/v2/p3_full_outputs.tsv"
SUMMARY_PATH = REPO_ROOT / "alphagenome_custom/metadata/v2/p3_full_summary.json"
DOH_ARGS = [
    "--doh-url",
    "https://dns.google/dns-query",
    "--resolve",
    "dns.google:443:8.8.8.8",
]
EXPECTED_RUN_COUNT = 482
EXPECTED_FASTQ_BYTES = 1_014_217_532_067
MIN_FREE_BYTES = 200 * 1024**3
LEDGER_LOCK = threading.Lock()
SRA_TOOLKIT_GLOB = "shared/tools/sratoolkit.*-ubuntu64/bin"


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--threads-per-sample", type=int, default=16)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument(
        "--work-dir", default="shared/source_reads/v2/full_streaming"
    )
    parser.add_argument(
        "--output-dir", default="alphagenome_custom/tracks/rna_seq_v2_normalized"
    )
    parser.add_argument(
        "--star-index", default="shared/reference_indexes/WBcel235_STAR_2.7.11b"
    )
    parser.add_argument(
        "--download-backend",
        choices=("ncbi_sra", "ena_fastq"),
        default="ncbi_sra",
    )
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def atomic_write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)


def write_tsv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise RuntimeError(f"Refusing to write an empty ledger: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(rows[0]),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def file_hashes(path: Path, chunk_size: int = 8 * 1024 * 1024) -> tuple[str, str]:
    sha_digest = hashlib.sha256()
    md5_digest = hashlib.md5(usedforsecurity=False)
    with path.open("rb") as handle:
        while chunk := handle.read(chunk_size):
            sha_digest.update(chunk)
            md5_digest.update(chunk)
    return sha_digest.hexdigest(), md5_digest.hexdigest()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(8 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def tree_size(path: Path) -> int:
    total = 0
    if not path.exists():
        return 0
    for item in path.rglob("*"):
        try:
            if item.is_file():
                total += item.stat().st_size
        except OSError:
            continue
    return total


class SampleRunner:
    def __init__(self, sample_dir: Path) -> None:
        self.sample_dir = sample_dir
        self.commands: list[dict[str, Any]] = []
        self.peak_bytes = tree_size(sample_dir)

    def run(self, command: list[str], *, stdout=None) -> None:
        print("command\t" + " ".join(command), flush=True)
        started = time.monotonic()
        process = subprocess.Popen(command, cwd=REPO_ROOT, stdout=stdout)
        while True:
            try:
                return_code = process.wait(timeout=3)
                self.peak_bytes = max(self.peak_bytes, tree_size(self.sample_dir))
                break
            except subprocess.TimeoutExpired:
                self.peak_bytes = max(self.peak_bytes, tree_size(self.sample_dir))
        self.commands.append(
            {
                "argv": command,
                "duration_seconds": round(time.monotonic() - started, 6),
                "return_code": return_code,
            }
        )
        if return_code != 0:
            raise subprocess.CalledProcessError(return_code, command)


def find_sra_toolkit() -> Path:
    matches = sorted(REPO_ROOT.glob(SRA_TOOLKIT_GLOB))
    valid = [
        path
        for path in matches
        if (path / "fasterq-dump").is_file() and (path / "vdb-validate").is_file()
    ]
    if len(valid) != 1:
        raise RuntimeError(
            f"Expected exactly one standalone SRA Toolkit at {SRA_TOOLKIT_GLOB}; "
            f"found={valid}"
        )
    return valid[0]


def require_tools(download_backend: str) -> Path | None:
    names = ["STAR", "samtools", "bedtools", "bedGraphToBigWig", "curl"]
    missing = [name for name in names if shutil.which(name) is None]
    if missing:
        raise RuntimeError("Missing required tools: " + ",".join(missing))
    return find_sra_toolkit() if download_backend == "ncbi_sra" else None


def read_reference_chromosomes() -> dict[str, int]:
    chromosomes = {}
    with (REPO_ROOT / "alphagenome_custom/reference/genome.fa.fai").open() as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            chromosomes[fields[0]] = int(fields[1])
    return chromosomes


def ensure_star_index(index_dir: Path, threads: int) -> None:
    if (index_dir / "Genome").is_file():
        return
    index_dir.mkdir(parents=True, exist_ok=True)
    subprocess.run(
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
        ],
        cwd=REPO_ROOT,
        check=True,
    )


def download_fastq(
    runner: SampleRunner,
    url: str,
    expected_md5: str,
    expected_bytes: int,
    path: Path,
) -> tuple[str, str]:
    if path.is_file() and path.stat().st_size == expected_bytes:
        actual_sha, actual_md5 = file_hashes(path)
        if actual_md5 == expected_md5:
            return actual_sha, actual_md5
    temporary = path.with_suffix(path.suffix + ".part")
    partial_complete = False
    if temporary.is_file() and temporary.stat().st_size == expected_bytes:
        _, partial_md5 = file_hashes(temporary)
        partial_complete = partial_md5 == expected_md5
    if not partial_complete:
        stagnant_attempts = 0
        previous_size = temporary.stat().st_size if temporary.exists() else 0
        for attempt in range(1, 1001):
            try:
                runner.run(
                    [
                        "curl",
                        "-fL",
                        "--connect-timeout",
                        "30",
                        "--continue-at",
                        "-",
                        *DOH_ARGS,
                        f"https://{url}",
                        "-o",
                        str(temporary),
                    ]
                )
            except subprocess.CalledProcessError:
                pass
            current_size = temporary.stat().st_size if temporary.exists() else 0
            if current_size == expected_bytes:
                break
            if current_size > expected_bytes:
                raise RuntimeError(f"Partial FASTQ exceeds expected size: {temporary}")
            stagnant_attempts = (
                stagnant_attempts + 1 if current_size <= previous_size else 0
            )
            if stagnant_attempts >= 20:
                raise RuntimeError(
                    f"FASTQ transfer made no progress for 20 attempts: {temporary}"
                )
            print(
                f"download_resume\t{path.name}\tattempt={attempt}\t"
                f"bytes={current_size}/{expected_bytes}",
                flush=True,
            )
            previous_size = current_size
        else:
            raise RuntimeError(f"FASTQ transfer attempt limit reached: {temporary}")
    if temporary.stat().st_size != expected_bytes:
        raise RuntimeError(f"FASTQ byte mismatch: {temporary}")
    actual_sha, actual_md5 = file_hashes(temporary)
    if actual_md5 != expected_md5:
        raise RuntimeError(f"FASTQ MD5 mismatch: {temporary}")
    temporary.replace(path)
    return actual_sha, actual_md5


def select_sra_archive(payload: dict[str, Any], accession: str) -> dict[str, Any]:
    results = payload.get("result", [])
    if len(results) != 1 or results[0].get("status") != 200:
        raise RuntimeError(f"NCBI SDL lookup failed for {accession}: {payload}")
    candidates = [
        item
        for item in results[0].get("files", [])
        if item.get("type") == "sra" and item.get("name") == accession
    ]
    if len(candidates) != 1:
        raise RuntimeError(f"Expected one full SRA archive for {accession}")
    candidate = candidates[0]
    links = [
        location.get("link")
        for location in candidate.get("locations", [])
        if location.get("service") == "s3"
        and str(location.get("link", "")).startswith(
            "https://sra-pub-run-odp.s3.amazonaws.com/sra/"
        )
    ]
    if len(links) != 1:
        raise RuntimeError(f"Expected one public S3 SRA location for {accession}")
    md5 = str(candidate.get("md5", ""))
    size = int(candidate.get("size", 0))
    if len(md5) != 32 or size <= 0:
        raise RuntimeError(f"Invalid SRA integrity metadata for {accession}")
    return {"url": links[0], "md5": md5, "size": size}


def locate_sra_archive(
    runner: SampleRunner, accession: str, response_path: Path
) -> dict[str, Any]:
    temporary = response_path.with_suffix(response_path.suffix + ".tmp")
    runner.run(
        [
            "curl",
            "-fL",
            "--connect-timeout",
            "30",
            *DOH_ARGS,
            (
                "https://locate.ncbi.nlm.nih.gov/sdl/2/retrieve"
                f"?acc={accession}&accept-alternate-locations=yes"
            ),
            "-o",
            str(temporary),
        ]
    )
    payload = json.loads(temporary.read_text())
    archive = select_sra_archive(payload, accession)
    temporary.replace(response_path)
    return archive


def download_sra_archive(
    runner: SampleRunner,
    archive: dict[str, Any],
    path: Path,
) -> tuple[str, str]:
    url = str(archive["url"])
    if not url.startswith("https://sra-pub-run-odp.s3.amazonaws.com/sra/"):
        raise RuntimeError(f"Unexpected SRA archive host: {url}")
    expected_md5 = str(archive["md5"])
    expected_bytes = int(archive["size"])
    return download_fastq(
        runner,
        url.removeprefix("https://"),
        expected_md5,
        expected_bytes,
        path,
    )


def fastq_stats(paths: list[Path]) -> tuple[list[str], int, int]:
    hashes = []
    total_reads = 0
    total_bytes = 0
    for path in paths:
        digest = hashlib.sha256()
        line_count = 0
        with path.open("rb") as handle:
            while chunk := handle.read(8 * 1024 * 1024):
                digest.update(chunk)
                line_count += chunk.count(b"\n")
        if line_count % 4:
            raise RuntimeError(f"FASTQ line count is not divisible by four: {path}")
        hashes.append(digest.hexdigest())
        total_reads += line_count // 4
        total_bytes += path.stat().st_size
    return hashes, total_reads, total_bytes


def expected_fastq_records(library_layout: str, manifest_read_count: int) -> int:
    if library_layout == "SINGLE":
        return manifest_read_count
    if library_layout == "PAIRED":
        return 2 * manifest_read_count
    raise ValueError(f"Unsupported library layout: {library_layout}")


def extract_sra_fastq(
    runner: SampleRunner,
    toolkit_bin: Path,
    archive_path: Path,
    accession: str,
    library_layout: str,
    expected_spots: int,
    threads: int,
) -> tuple[list[Path], list[str], int]:
    if library_layout == "PAIRED":
        paths = [
            runner.sample_dir / f"{accession}_1.fastq",
            runner.sample_dir / f"{accession}_2.fastq",
        ]
    else:
        paths = [runner.sample_dir / f"{accession}.fastq"]
    if not all(path.is_file() for path in paths):
        for path in runner.sample_dir.glob(f"{accession}*.fastq"):
            path.unlink()
        temporary = runner.sample_dir / "fasterq_tmp"
        if temporary.exists():
            shutil.rmtree(temporary)
        temporary.mkdir()
        runner.run(
            [
                str(toolkit_bin / "fasterq-dump"),
                "--threads",
                str(threads),
                "--split-files",
                "--outdir",
                str(runner.sample_dir),
                "--temp",
                str(temporary),
                str(archive_path),
            ]
        )
        shutil.rmtree(temporary)
    if not all(path.is_file() for path in paths):
        raise RuntimeError(f"Missing extracted FASTQ for {accession}: {paths}")
    hashes, read_count, total_bytes = fastq_stats(paths)
    expected_records = expected_fastq_records(library_layout, expected_spots)
    if read_count != expected_records:
        raise RuntimeError(
            f"Extracted read count mismatch for {accession}: "
            f"{read_count} != {expected_records}"
        )
    return paths, hashes, total_bytes


def normalize_bedgraph(input_path: Path, output_path: Path) -> tuple[float, float]:
    total = 0.0
    with input_path.open() as handle:
        for line in handle:
            chromosome, start, end, value = line.rstrip("\n").split("\t")
            numeric_value = float(value)
            if not math.isfinite(numeric_value) or numeric_value < 0:
                raise RuntimeError(
                    f"Invalid coverage at {chromosome}:{start}-{end}: {value}"
                )
            total += (int(end) - int(start)) * numeric_value
    if not math.isfinite(total) or total <= 0:
        raise RuntimeError(f"Non-positive coverage total: {input_path}")
    scale = 100_000_000.0 / total
    temporary = output_path.with_suffix(output_path.suffix + ".tmp")
    with input_path.open() as source, temporary.open("w") as destination:
        for line in source:
            chromosome, start, end, value = line.rstrip("\n").split("\t")
            normalized = float(value) * scale
            if not math.isfinite(normalized) or normalized < 0:
                raise RuntimeError(f"Invalid normalized value in {input_path}")
            destination.write(
                f"{chromosome}\t{start}\t{end}\t{normalized:.9g}\n"
            )
    temporary.replace(output_path)
    return scale, total


def parse_star_log(path: Path) -> dict[str, str]:
    values = {}
    for line in path.read_text().splitlines():
        if "|" in line:
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


def validate_bigwig(path: Path, chromosomes: dict[str, int]) -> dict[str, float]:
    with pyBigWig.open(str(path)) as bigwig:
        if not bigwig.isBigWig() or bigwig.chroms() != chromosomes:
            raise RuntimeError(f"Invalid bigWig chromosomes: {path}")
        header = bigwig.header()
    values = {
        key: float(header[key]) for key in ("minVal", "maxVal", "sumData", "sumSquared")
    }
    if not all(math.isfinite(value) for value in values.values()):
        raise RuntimeError(f"Non-finite bigWig statistics: {path}")
    if values["minVal"] < 0 or values["maxVal"] < 0:
        raise RuntimeError(f"Negative bigWig signal: {path}")
    relative_error = abs(values["sumData"] - 100_000_000.0) / 100_000_000.0
    if relative_error > 1e-5:
        raise RuntimeError(f"Unexpected total signal in {path}: {values['sumData']}")
    values["relative_total_error"] = relative_error
    return values


def validate_manifests(
    sources: list[dict[str, str]], samples: list[dict[str, str]]
) -> dict[str, dict[str, str]]:
    sample_by_run = {row["run_accession"]: row for row in samples}
    if len(sources) != EXPECTED_RUN_COUNT:
        raise RuntimeError(f"Expected {EXPECTED_RUN_COUNT} RNA sources, got {len(sources)}")
    if len({row["run_accession"] for row in sources}) != EXPECTED_RUN_COUNT:
        raise RuntimeError("Full source run accessions are not unique")
    expected_runs = {
        row["run_accession"] for row in samples if row["assay"] == "RNA-Seq"
    }
    if {row["run_accession"] for row in sources} != expected_runs:
        raise RuntimeError("Full source manifest does not equal the RNA-seq sample set")
    total_bytes = 0
    for row in sources:
        file_count = 2 if row["library_layout"] == "PAIRED" else 1
        fields = [row[name].split(";") for name in ("fastq_ftp", "fastq_md5", "fastq_bytes")]
        if any(len(values) != file_count for values in fields):
            raise RuntimeError(f"FASTQ list mismatch: {row['run_accession']}")
        if any(len(value) != 32 for value in fields[1]):
            raise RuntimeError(f"FASTQ MD5 invalid: {row['run_accession']}")
        total_bytes += sum(int(value) for value in fields[2])
    if total_bytes != EXPECTED_FASTQ_BYTES:
        raise RuntimeError(f"Unexpected full FASTQ bytes: {total_bytes}")
    return sample_by_run


def validate_sra_manifest(rows: list[dict[str, str]]) -> dict[str, dict[str, str]]:
    if len(rows) != EXPECTED_RUN_COUNT:
        raise RuntimeError(f"Expected {EXPECTED_RUN_COUNT} SRA sources, got {len(rows)}")
    by_run = {row["run_accession"]: row for row in rows}
    if len(by_run) != EXPECTED_RUN_COUNT:
        raise RuntimeError("SRA source run accessions are not unique")
    expected_runs = {
        row["run_accession"] for row in read_tsv(SOURCE_MANIFEST)
    }
    if set(by_run) != expected_runs:
        raise RuntimeError("SRA transport manifest does not equal the P3 RNA scope")
    for accession, row in by_run.items():
        expected_url = (
            f"https://sra-pub-run-odp.s3.amazonaws.com/sra/{accession}/{accession}"
        )
        if (
            row["sra_url"] != expected_url
            or len(row["sra_md5"]) != 32
            or int(row["sra_bytes"]) <= 0
        ):
            raise RuntimeError(f"Invalid SRA transport record: {accession}")
    return by_run


def process_sample(
    source: dict[str, str],
    sample: dict[str, str],
    work_dir: Path,
    output_dir: Path,
    index_dir: Path,
    chrom_sizes: Path,
    chromosomes: dict[str, int],
    threads: int,
    download_backend: str,
    toolkit_bin: Path | None,
    sra_source: dict[str, str] | None,
) -> dict[str, Any]:
    accession = source["run_accession"]
    sample_dir = work_dir / accession
    audit_path = output_dir / "sample_audits" / f"{accession}.json"
    output_path = output_dir / f"{accession}.bw"
    if audit_path.is_file() and output_path.is_file():
        audit = json.loads(audit_path.read_text())
        if (
            audit.get("output_sha256") == sha256(output_path)
            and audit.get("source_manifest_sha256") == sha256(SOURCE_MANIFEST)
            and (
                download_backend != "ncbi_sra"
                or audit.get("sra_manifest_sha256") == sha256(SRA_MANIFEST)
            )
        ):
            validate_bigwig(output_path, chromosomes)
            if sample_dir.exists():
                shutil.rmtree(sample_dir)
            audit["cleanup_status"] = "completed"
            atomic_write_json(audit_path, audit)
            (output_dir / "failures" / f"{accession}.json").unlink(
                missing_ok=True
            )
            return audit

    if sample_dir.exists():
        for item in sample_dir.iterdir():
            if item.is_file() and (
                item.name.endswith(".fastq.gz")
                or item.name.endswith(".fastq.gz.part")
                or item.name.endswith(".fastq")
                or item.name.endswith(".sra")
                or item.name.endswith(".sra.part")
            ):
                continue
            if item.is_dir():
                shutil.rmtree(item)
            else:
                item.unlink()
    else:
        sample_dir.mkdir(parents=True)
    runner = SampleRunner(sample_dir)
    started_at = utc_now()
    started = time.monotonic()
    try:
        provided_path = REPO_ROOT / sample["local_path"]
        provided_before = sha256(provided_path)
        if provided_before != sample["sha256"]:
            raise RuntimeError(f"Provided bigWig changed before processing: {accession}")

        urls = source["fastq_ftp"].split(";")
        expected_md5 = source["fastq_md5"].split(";")
        expected_bytes = [int(value) for value in source["fastq_bytes"].split(";")]
        source_transport: dict[str, Any]
        read_files_command: list[str]
        if download_backend == "ncbi_sra":
            if toolkit_bin is None:
                raise RuntimeError("SRA Toolkit was not configured")
            if sra_source is None:
                raise RuntimeError(f"Missing locked SRA source for {accession}")
            archive = {
                "url": sra_source["sra_url"],
                "md5": sra_source["sra_md5"],
                "size": int(sra_source["sra_bytes"]),
            }
            archive_path = sample_dir / f"{accession}.sra"
            archive_sha, archive_md5 = download_sra_archive(
                runner, archive, archive_path
            )
            runner.run([str(toolkit_bin / "vdb-validate"), str(archive_path)])
            fastq_paths, fastq_sha, extracted_bytes = extract_sra_fastq(
                runner,
                toolkit_bin,
                archive_path,
                accession,
                source["library_layout"],
                int(source["read_count"]),
                threads,
            )
            source_transport = {
                "source_transport_backend": "ncbi_sra",
                "source_archive_url": archive["url"],
                "source_archive_md5": archive_md5,
                "source_archive_sha256": archive_sha,
                "source_archive_bytes": archive["size"],
                "extracted_fastq_sha256": ";".join(fastq_sha),
                "extracted_fastq_bytes": extracted_bytes,
                "source_manifest_spot_count": int(source["read_count"]),
                "extracted_fastq_record_count": expected_fastq_records(
                    source["library_layout"], int(source["read_count"])
                ),
            }
            read_files_command = []
        else:
            fastq_paths = []
            fastq_sha = []
            for url, md5_value, byte_count in zip(
                urls, expected_md5, expected_bytes, strict=True
            ):
                path = sample_dir / Path(url).name
                sha_value, _ = download_fastq(
                    runner, url, md5_value, byte_count, path
                )
                fastq_paths.append(path)
                fastq_sha.append(sha_value)
            source_transport = {
                "source_transport_backend": "ena_fastq",
                "source_archive_url": "",
                "source_archive_md5": "",
                "source_archive_sha256": "",
                "source_archive_bytes": "",
                "extracted_fastq_sha256": ";".join(fastq_sha),
                "extracted_fastq_bytes": sum(expected_bytes),
                "source_manifest_spot_count": int(source["read_count"]),
                "extracted_fastq_record_count": expected_fastq_records(
                    source["library_layout"], int(source["read_count"])
                ),
            }
            read_files_command = ["--readFilesCommand", "zcat"]

        star_prefix = sample_dir / "star_"
        runner.run(
            [
                "STAR",
                "--runThreadN",
                str(threads),
                "--genomeDir",
                str(index_dir),
                "--readFilesIn",
                *map(str, fastq_paths),
                *read_files_command,
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
        aligned_bam = sample_dir / "star_Aligned.sortedByCoord.out.bam"
        star_metrics = parse_star_log(sample_dir / "star_Log.final.out")
        primary_bam = sample_dir / "primary_unique.bam"
        runner.run(
            [
                "samtools",
                "view",
                "-@",
                str(max(threads // 2, 1)),
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
        runner.run(["samtools", "quickcheck", "-v", str(primary_bam)])
        primary_bam_sha = sha256(primary_bam)
        primary_bam_bytes = primary_bam.stat().st_size

        raw_bedgraph = sample_dir / "coverage.raw.bedGraph"
        with raw_bedgraph.open("w") as handle:
            runner.run(
                [
                    "bedtools",
                    "genomecov",
                    "-split",
                    "-bg",
                    "-ibam",
                    str(primary_bam),
                ],
                stdout=handle,
            )
        normalized_bedgraph = sample_dir / "coverage.normalized.bedGraph"
        scale, raw_total = normalize_bedgraph(raw_bedgraph, normalized_bedgraph)
        temporary_bigwig = output_path.with_suffix(".bw.tmp")
        runner.run(
            [
                "bedGraphToBigWig",
                str(normalized_bedgraph),
                str(chrom_sizes),
                str(temporary_bigwig),
            ]
        )
        temporary_bigwig.replace(output_path)
        output_stats = validate_bigwig(output_path, chromosomes)
        provided_after = sha256(provided_path)
        if provided_after != provided_before:
            raise RuntimeError(f"Provided bigWig changed during processing: {accession}")

        audit = {
            "run_accession": accession,
            "experiment_accession": source["experiment_accession"],
            "batch": source["batch"],
            "library_layout": source["library_layout"],
            "library_selection": source["library_selection"],
            "source_fastq_urls": ";".join(urls),
            "source_fastq_md5": ";".join(expected_md5),
            "source_fastq_sha256": ";".join(fastq_sha),
            "source_fastq_bytes": sum(expected_bytes),
            **source_transport,
            "source_files_cleaned_after_verification": True,
            "primary_bam_sha256_before_cleanup": primary_bam_sha,
            "primary_bam_bytes_before_cleanup": primary_bam_bytes,
            "provided_bigwig_path": sample["local_path"],
            "provided_bigwig_sha256_before": provided_before,
            "provided_bigwig_sha256_after": provided_after,
            "star_input_reads": star_metrics["Number of input reads"],
            "star_uniquely_mapped_reads": star_metrics[
                "Uniquely mapped reads number"
            ],
            "star_uniquely_mapped_percent": star_metrics[
                "Uniquely mapped reads %"
            ],
            "coverage_policy": "primary_unique_spliced_unstranded",
            "output_strand": ".",
            "normalization_formula": "100000000/sum((end-start)*raw_coverage)",
            "raw_bedgraph_total": f"{raw_total:.12g}",
            "scale_to_1e8": f"{scale:.12g}",
            "output_path": str(output_path.relative_to(REPO_ROOT)),
            "output_size_bytes": output_path.stat().st_size,
            "output_sha256": sha256(output_path),
            "output_min": f"{output_stats['minVal']:.12g}",
            "output_max": f"{output_stats['maxVal']:.12g}",
            "output_total_signal": f"{output_stats['sumData']:.12g}",
            "output_relative_total_error": (
                f"{output_stats['relative_total_error']:.12g}"
            ),
            "source_manifest_sha256": sha256(SOURCE_MANIFEST),
            "sample_manifest_sha256": sha256(SAMPLE_MANIFEST),
            "sra_manifest_sha256": (
                sha256(SRA_MANIFEST) if download_backend == "ncbi_sra" else ""
            ),
            "commands_json": json.dumps(runner.commands, separators=(",", ":")),
            "sample_peak_work_bytes": runner.peak_bytes,
            "elapsed_seconds": f"{time.monotonic() - started:.6f}",
            "started_at": started_at,
            "completed_at": utc_now(),
            "cleanup_status": "pending",
        }
        atomic_write_json(audit_path, audit)
        shutil.rmtree(sample_dir)
        audit["cleanup_status"] = "completed"
        atomic_write_json(audit_path, audit)
        (output_dir / "failures" / f"{accession}.json").unlink(missing_ok=True)
        return audit
    except Exception as error:
        atomic_write_json(
            output_dir / "failures" / f"{accession}.json",
            {
                "run_accession": accession,
                "failed_at": utc_now(),
                "error_type": type(error).__name__,
                "error": str(error),
                "commands": runner.commands,
                "work_dir_retained": str(sample_dir.relative_to(REPO_ROOT)),
            },
        )
        raise


def refresh_ledger(output_dir: Path) -> list[dict[str, Any]]:
    rows = [
        json.loads(path.read_text())
        for path in sorted((output_dir / "sample_audits").glob("*.json"))
    ]
    if rows:
        scalar_rows = []
        for row in rows:
            scalar_rows.append(
                {
                    key: json.dumps(value, separators=(",", ":"))
                    if isinstance(value, (dict, list))
                    else value
                    for key, value in row.items()
                }
            )
        write_tsv(LEDGER_PATH, scalar_rows)
    return rows


def main() -> None:
    args = parse_args()
    if args.workers < 1 or args.threads_per_sample < 1:
        raise ValueError("workers and threads-per-sample must be positive")
    if args.workers * args.threads_per_sample > (os.cpu_count() or 1):
        raise ValueError("Requested worker threads exceed available logical CPUs")
    toolkit_bin = require_tools(args.download_backend)
    sources = read_tsv(SOURCE_MANIFEST)
    samples = read_tsv(SAMPLE_MANIFEST)
    sample_by_run = validate_manifests(sources, samples)
    sra_by_run = (
        validate_sra_manifest(read_tsv(SRA_MANIFEST))
        if args.download_backend == "ncbi_sra"
        else {}
    )
    work_dir = REPO_ROOT / args.work_dir
    output_dir = REPO_ROOT / args.output_dir
    index_dir = REPO_ROOT / args.star_index
    work_dir.mkdir(parents=True, exist_ok=True)
    output_dir.mkdir(parents=True, exist_ok=True)
    if shutil.disk_usage(REPO_ROOT).free < MIN_FREE_BYTES:
        raise RuntimeError(f"Full streaming requires at least {MIN_FREE_BYTES} free bytes")
    ensure_star_index(index_dir, args.threads_per_sample)
    chromosomes = read_reference_chromosomes()
    chrom_sizes = work_dir / "WBcel235.chrom.sizes"
    with chrom_sizes.open("w") as handle:
        for chromosome, length in chromosomes.items():
            handle.write(f"{chromosome}\t{length}\n")

    started_at = utc_now()
    started = time.monotonic()
    refresh_ledger(output_dir)
    pending = list(sources)
    pending.sort(key=lambda row: sum(map(int, row["fastq_bytes"].split(";"))))
    failures = []
    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        futures: dict[Future[dict[str, Any]], str] = {}
        source_iterator = iter(pending)

        def submit(source: dict[str, str]) -> None:
            accession = source["run_accession"]
            future = executor.submit(
                process_sample,
                source,
                sample_by_run[accession],
                work_dir,
                output_dir,
                index_dir,
                chrom_sizes,
                chromosomes,
                args.threads_per_sample,
                args.download_backend,
                toolkit_bin,
                sra_by_run.get(accession),
            )
            futures[future] = accession

        for _ in range(args.workers):
            try:
                submit(next(source_iterator))
            except StopIteration:
                break
        while futures:
            future = next(as_completed(futures))
            accession = futures.pop(future)
            try:
                future.result()
            except Exception as error:
                failures.append(f"{accession}:{type(error).__name__}:{error}")
            with LEDGER_LOCK:
                rows = refresh_ledger(output_dir)
                atomic_write_json(
                    output_dir / "progress.json",
                    {
                        "expected": EXPECTED_RUN_COUNT,
                        "completed": len(rows),
                        "failed": failures,
                        "updated_at": utc_now(),
                    },
                )
            if not failures:
                try:
                    submit(next(source_iterator))
                except StopIteration:
                    pass
    rows = refresh_ledger(output_dir)
    summary = {
        "schema_version": 1,
        "phase": "P3B",
        "scope": "full_rna_streaming_482",
        "source_manifest": str(SOURCE_MANIFEST.relative_to(REPO_ROOT)),
        "source_manifest_sha256": sha256(SOURCE_MANIFEST),
        "sample_manifest_sha256": sha256(SAMPLE_MANIFEST),
        "sra_manifest_sha256": (
            sha256(SRA_MANIFEST) if args.download_backend == "ncbi_sra" else None
        ),
        "expected_runs": EXPECTED_RUN_COUNT,
        "completed_runs": len(rows),
        "source_fastq_bytes": EXPECTED_FASTQ_BYTES,
        "download_backend": args.download_backend,
        "sra_toolkit_bin": (
            str(toolkit_bin.relative_to(REPO_ROOT)) if toolkit_bin else None
        ),
        "workers": args.workers,
        "threads_per_sample": args.threads_per_sample,
        "work_dir": args.work_dir,
        "output_dir": args.output_dir,
        "star_index": args.star_index,
        "failures": failures,
        "started_at": started_at,
        "completed_at": utc_now(),
        "elapsed_seconds": round(time.monotonic() - started, 6),
    }
    atomic_write_json(SUMMARY_PATH, summary)
    print(json.dumps(summary, indent=2, sort_keys=True))
    if failures or len(rows) != EXPECTED_RUN_COUNT:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
