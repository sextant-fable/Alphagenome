#!/usr/bin/env python3
"""Freeze the small I--V published-eQTL candidate set for variant scoring.

The candidate positions come from Supplementary Data 4 of the independent
wild-strain study.  This preparation step emits only I--V records and checks
their REF alleles against the local WBcel235 FASTA.  It does not load a model
or score a variant.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from collections import defaultdict

from scripts.v2_variant_scoring import fetch_fasta, read_fai


REPO_ROOT = Path(__file__).resolve().parents[1]
SOURCE_DIR = REPO_ROOT / "shared/external_eqtl/p20_sources_v1"
SUPP_DIR = SOURCE_DIR / "supplementary/Supplementary Data 1-4"
FASTA = REPO_ROOT / "alphagenome_custom/reference/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa"
FAI = REPO_ROOT / "alphagenome_custom/reference/genome.fa.fai"
CHROM_ORDER = {chrom: index for index, chrom in enumerate(("I", "II", "III", "IV", "V"))}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=REPO_ROOT / "results/v2_p20_eqtl_variant_set")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    output = args.output_dir.resolve()
    if output.exists():
        raise FileExistsError(f"Output directory is immutable: {output}")
    counts_path = SOURCE_DIR / "GSE186719_Celegans_208strains_609samples_rawCounts.tsv.gz"
    vcf_path = SOURCE_DIR / "WI.20200815.hard-filter.isotype.vcf.gz"
    data2_path = SUPP_DIR / "Supplementary Data 2.tsv"
    data4_path = SUPP_DIR / "Supplementary Data 4.tsv"
    required = (counts_path, vcf_path, data2_path, data4_path, FASTA, FAI)
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        raise FileNotFoundError("Missing P20 source/reference files: " + ", ".join(missing))

    data2 = {row["transcript"]: row for row in read_tsv(data2_path)}
    data4_rows = read_tsv(data4_path)
    by_variant: dict[tuple[str, int], dict[str, object]] = {}
    annotations: defaultdict[tuple[str, int], list[str]] = defaultdict(list)
    for row in data4_rows:
        chromosome = row.get("eQTL_Chr", "")
        raw_position = row.get("finemapping_variant", "")
        transcript = row.get("transcript", "")
        if chromosome not in CHROM_ORDER or not raw_position.isdigit() or transcript not in data2:
            continue
        position = int(raw_position)
        key = (chromosome, position)
        joined = data2[transcript]
        annotations[key].append(
            ":".join((transcript, joined.get("GeneName", ""), joined.get("eQTL_classification", ""), row.get("finemapping_gene", "")))
        )
        if key not in by_variant:
            by_variant[key] = {
                "chromosome": chromosome,
                "position_1based": position,
                "representative_transcript": transcript,
                "representative_gene_name": joined.get("GeneName", ""),
                "representative_gene_id": joined.get("WormBaseGeneID", ""),
                "eQTL_classification": joined.get("eQTL_classification", ""),
                "published_eQTL_peak": joined.get("eQTL_peak", ""),
                "published_eQTL_logP": joined.get("logP", ""),
                "published_eQTL_var_exp": joined.get("var_exp", ""),
                "finemapping_gene": row.get("finemapping_gene", ""),
                "variant_impact": row.get("VARIANT_IMPACT", ""),
            }
    candidates = [by_variant[key] for key in sorted(by_variant, key=lambda item: (CHROM_ORDER[item[0]], item[1]))]
    for row in candidates:
        key = (str(row["chromosome"]), int(row["position_1based"]))
        row["annotation_count"] = len(annotations[key])
        row["annotations"] = ";".join(sorted(set(annotations[key])))
    if not candidates:
        raise RuntimeError("No I--V candidates survived the published eQTL truth-table filters")

    output.mkdir(parents=True)
    candidate_path = output / "published_eqtl_candidate_annotations.tsv"
    candidate_fields = list(candidates[0])
    write_tsv(candidate_path, candidates, candidate_fields)

    target_keys = {(str(row["chromosome"]), int(row["position_1based"])) for row in candidates}
    selected_vcf = output / "published_eqtl_candidates_iv.vcf"
    selected_rows = 0
    sample_count = None
    coordinate_rows: list[dict[str, object]] = []
    fai = read_fai(FAI)
    with gzip.open(vcf_path, "rt", encoding="utf-8") as source, selected_vcf.open("w", encoding="utf-8") as target:
        for line in source:
            if line.startswith("#"):
                target.write(line)
                if line.startswith("#CHROM"):
                    sample_count = max(0, len(line.rstrip("\n").split("\t")) - 9)
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                continue
            key = (fields[0], int(fields[1]))
            if key not in target_keys:
                continue
            target.write(line)
            selected_rows += 1
            reference = fields[3].upper()
            alternate = fields[4].split(",")[0].upper()
            fasta_reference = fetch_fasta(FASTA, fai, fields[0], int(fields[1]) - 1, int(fields[1]) - 1 + len(reference))
            coordinate_rows.append(
                {
                    "chromosome": fields[0],
                    "position_1based": fields[1],
                    "vcf_id": fields[2],
                    "vcf_ref": reference,
                    "vcf_alt_first": alternate,
                    "wbcel235_reference": fasta_reference,
                    "coordinate_status": "MATCH" if reference == fasta_reference else "REF_MISMATCH",
                }
            )
    if sample_count is None:
        raise RuntimeError("Selected VCF lacks a #CHROM header")
    coordinate_path = output / "published_eqtl_coordinate_audit.tsv"
    write_tsv(coordinate_path, coordinate_rows, list(coordinate_rows[0]) if coordinate_rows else ["chromosome"])
    mismatch_count = sum(row["coordinate_status"] != "MATCH" for row in coordinate_rows)
    manifest = {
        "schema_version": 1,
        "contract": "p20_eqtl_variant_set_v1",
        "created_at": utc_now(),
        "scope": "I--V candidate extraction and WBcel235 REF audit; no model inference or locked-test access",
        "source_counts_sha256": sha256(counts_path),
        "source_vcf_sha256": sha256(vcf_path),
        "data2_sha256": sha256(data2_path),
        "data4_sha256": sha256(data4_path),
        "fasta_sha256": sha256(FASTA),
        "candidate_rows": len(candidates),
        "unique_candidate_positions": len(target_keys),
        "selected_vcf_records": selected_rows,
        "vcf_sample_count": sample_count,
        "coordinate_audit_rows": len(coordinate_rows),
        "coordinate_ref_mismatches": mismatch_count,
        "candidate_annotations": str(candidate_path.relative_to(REPO_ROOT)),
        "selected_vcf": str(selected_vcf.relative_to(REPO_ROOT)),
        "coordinate_audit": str(coordinate_path.relative_to(REPO_ROOT)),
        "status": "PASS" if mismatch_count == 0 and selected_rows >= len(target_keys) else "BLOCKED_COORDINATE_OR_VCF_COVERAGE",
        "prohibited_claims": ["independent regulatory-allele validation", "causal biological mechanism", "unseen-condition prediction"],
    }
    manifest["selected_vcf_sha256"] = sha256(selected_vcf)
    manifest_path = output / "p20_eqtl_variant_set_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(manifest, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
