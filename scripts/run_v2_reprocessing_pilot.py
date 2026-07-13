#!/usr/bin/env python3
"""Run the approved RNA-seq v2 source-read reprocessing pilot."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
from typing import Iterable


REPO_ROOT = Path(__file__).resolve().parents[1]
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


def run(command: list[str], *, stdout=None) -> None:
    print("command\t" + " ".join(command), flush=True)
    subprocess.run(command, cwd=REPO_ROOT, check=True, stdout=stdout)


def download_fastq(url: str, expected_md5: str, expected_bytes: int, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not path.is_file() or path.stat().st_size != expected_bytes:
        tmp_path = path.with_suffix(path.suffix + ".part")
        run(
            [
                "curl",
                "-fL",
                "--retry",
                "4",
                "--retry-all-errors",
                *DOH_ARGS,
                f"https://{url}",
                "-o",
                str(tmp_path),
            ]
        )
        tmp_path.replace(path)
    if path.stat().st_size != expected_bytes or md5(path) != expected_md5:
        raise RuntimeError(f"FASTQ integrity failure: {path}")


def normalize_bedgraph(input_path: Path, output_path: Path) -> float:
    total = 0.0
    with input_path.open() as handle:
        for line in handle:
            chromosome, start, end, value = line.rstrip("\n").split("\t")
            total += (int(end) - int(start)) * float(value)
    if total <= 0:
        raise RuntimeError(f"Non-positive bedGraph signal total: {input_path}")
    scale = 100_000_000.0 / total
    tmp_path = output_path.with_suffix(output_path.suffix + ".tmp")
    with input_path.open() as source, tmp_path.open("w") as destination:
        for line in source:
            chromosome, start, end, value = line.rstrip("\n").split("\t")
            destination.write(
                f"{chromosome}\t{start}\t{end}\t{float(value) * scale:.9g}\n"
            )
    tmp_path.replace(output_path)
    return scale


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


def main() -> None:
    args = parse_args()
    require_tools(["STAR", "samtools", "bedtools", "bedGraphToBigWig", "curl"])
    manifest_path = REPO_ROOT / args.manifest
    work_dir = REPO_ROOT / args.work_dir
    output_dir = REPO_ROOT / args.output_dir
    index_dir = REPO_ROOT / args.star_index
    work_dir.mkdir(parents=True, exist_ok=True)
    output_dir.mkdir(parents=True, exist_ok=True)
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
    for row in read_tsv(manifest_path):
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
            run(["samtools", "index", "-@", "4", str(primary_bam)])

        raw_bedgraph = sample_dir / "coverage.raw.bedGraph"
        with raw_bedgraph.open("w") as handle:
            run(
                ["bedtools", "genomecov", "-split", "-bg", "-ibam", str(primary_bam)],
                stdout=handle,
            )
        normalized_bedgraph = sample_dir / "coverage.normalized.bedGraph"
        scale = normalize_bedgraph(raw_bedgraph, normalized_bedgraph)
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
        report_rows.append(
            {
                "run_accession": accession,
                "output_path": str(bigwig_path.relative_to(REPO_ROOT)),
                "output_size_bytes": bigwig_path.stat().st_size,
                "output_sha256": sha256(bigwig_path),
                "bedgraph_scale_to_1e8": f"{scale:.12g}",
                "output_strand": ".",
                "source_fastq_bytes": sum(sizes),
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
    (output_dir / "pilot_run.json").write_text(
        json.dumps(
            {
                "manifest": args.manifest,
                "threads": args.threads,
                "work_dir": args.work_dir,
                "output_dir": args.output_dir,
                "star_index": args.star_index,
                "outputs": len(report_rows),
            },
            indent=2,
            sort_keys=True,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
