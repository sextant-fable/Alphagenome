#!/usr/bin/env python3
"""Build the evidence-backed 485-accession RNA-seq v2 sample manifest."""

from __future__ import annotations

from collections import defaultdict
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
import subprocess
from typing import Any, Iterable
from urllib.parse import quote
from xml.etree import ElementTree as ET
from zipfile import ZipFile


REPO_ROOT = Path(__file__).resolve().parents[1]
METADATA_DIR = REPO_ROOT / "alphagenome_custom/metadata/v2"
CACHE_DIR = REPO_ROOT / "shared/metadata_sources/v2"
SHAREPOINT_DIR = REPO_ROOT / "shared/sharepoint_2026-07-13/training_input_bigwig"

DOH_ARGS = [
    "--doh-url",
    "https://dns.google/dns-query",
    "--resolve",
    "dns.google:443:8.8.8.8",
]
ENA_FIELDS = [
    "run_accession",
    "experiment_accession",
    "sample_accession",
    "secondary_sample_accession",
    "study_accession",
    "secondary_study_accession",
    "scientific_name",
    "sample_title",
    "description",
    "experiment_title",
    "study_title",
    "age",
    "dev_stage",
    "tissue_type",
    "cell_type",
    "sex",
    "strain",
    "sub_strain",
    "host_genotype",
    "experimental_factor",
    "library_name",
    "library_strategy",
    "library_source",
    "library_selection",
    "library_layout",
    "read_strand",
    "instrument_platform",
    "instrument_model",
    "center_name",
    "read_count",
    "base_count",
    "fastq_md5",
    "fastq_bytes",
    "fastq_ftp",
    "bam_bytes",
    "bam_ftp",
    "first_public",
]


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def rel(path: Path) -> str:
    return str(path.relative_to(REPO_ROOT))


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, Any]], fields: Iterable[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(fields), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def sha256(path: Path, chunk_size: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def cell_column(reference: str) -> int:
    letters = "".join(ch for ch in reference if ch.isalpha())
    value = 0
    for letter in letters:
        value = value * 26 + ord(letter.upper()) - ord("A") + 1
    return value - 1


def read_xlsx_first_sheet(path: Path) -> list[dict[str, str]]:
    namespace = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
    ns = {"m": namespace}
    with ZipFile(path) as archive:
        shared_strings: list[str] = []
        if "xl/sharedStrings.xml" in archive.namelist():
            root = ET.fromstring(archive.read("xl/sharedStrings.xml"))
            shared_strings = [
                "".join(node.text or "" for node in item.iter(f"{{{namespace}}}t"))
                for item in root.findall("m:si", ns)
            ]
        sheet = ET.fromstring(archive.read("xl/worksheets/sheet1.xml"))
        values: list[list[str]] = []
        for row in sheet.findall(".//m:sheetData/m:row", ns):
            row_values: dict[int, str] = {}
            for cell in row.findall("m:c", ns):
                value_node = cell.find("m:v", ns)
                value = "" if value_node is None else value_node.text or ""
                if cell.get("t") == "s" and value:
                    value = shared_strings[int(value)]
                elif cell.get("t") == "inlineStr":
                    value = "".join(
                        node.text or "" for node in cell.iter(f"{{{namespace}}}t")
                    )
                row_values[cell_column(cell.get("r", "A1"))] = value.strip()
            width = max(row_values, default=-1) + 1
            values.append([row_values.get(index, "") for index in range(width)])
    if not values:
        return []
    headers = values[0]
    rows = []
    for values_row in values[1:]:
        padded = values_row + [""] * (len(headers) - len(values_row))
        row = {header: padded[index] for index, header in enumerate(headers) if header}
        if any(row.values()):
            rows.append(row)
    return rows


def curl_to_file(args: list[str], output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = output_path.with_suffix(output_path.suffix + ".tmp")
    command = [
        "curl",
        "-fsSL",
        "--retry",
        "4",
        "--retry-all-errors",
        "--connect-timeout",
        "30",
        "--max-time",
        "300",
        *DOH_ARGS,
        *args,
        "-o",
        str(tmp_path),
    ]
    subprocess.run(command, cwd=REPO_ROOT, check=True)
    if tmp_path.stat().st_size == 0:
        raise RuntimeError(f"Empty response for {' '.join(args)}")
    tmp_path.replace(output_path)


def fetch_ena_runs(accessions: list[str]) -> tuple[dict[str, dict[str, str]], list[str]]:
    rows_by_accession: dict[str, dict[str, str]] = {}
    failures: list[str] = []
    chunks_dir = CACHE_DIR / "ena_runs"
    chunks_dir.mkdir(parents=True, exist_ok=True)
    for chunk_index in range(0, len(accessions), 40):
        chunk = accessions[chunk_index : chunk_index + 40]
        output_path = chunks_dir / f"chunk_{chunk_index // 40:03d}.tsv"
        query = " OR ".join(f'run_accession="{accession}"' for accession in chunk)
        cache_valid = False
        if output_path.is_file():
            with output_path.open() as handle:
                cached_fields = handle.readline().rstrip("\n").split("\t")
            cache_valid = set(ENA_FIELDS).issubset(cached_fields)
        if not cache_valid:
            curl_to_file(
                [
                    "-G",
                    "https://www.ebi.ac.uk/ena/portal/api/search",
                    "--data-urlencode",
                    "result=read_run",
                    "--data-urlencode",
                    f"query={query}",
                    "--data-urlencode",
                    f"fields={','.join(ENA_FIELDS)}",
                    "--data-urlencode",
                    "format=tsv",
                    "--data-urlencode",
                    "limit=0",
                ],
                output_path,
            )
        for row in read_tsv(output_path):
            rows_by_accession[row["run_accession"]] = row
        failures.extend(accession for accession in chunk if accession not in rows_by_accession)
    return rows_by_accession, sorted(set(failures))


def fetch_xml_chunks(
    accessions: list[str], cache_subdir: str, entity_name: str
) -> list[Path]:
    output_paths = []
    directory = CACHE_DIR / cache_subdir
    directory.mkdir(parents=True, exist_ok=True)
    for chunk_index in range(0, len(accessions), 40):
        chunk = accessions[chunk_index : chunk_index + 40]
        output_path = directory / f"chunk_{chunk_index // 40:03d}.xml"
        if not output_path.is_file():
            joined = ",".join(quote(accession, safe="") for accession in chunk)
            curl_to_file(
                [
                    f"https://www.ebi.ac.uk/ena/browser/api/xml/{joined}?download=false"
                ],
                output_path,
            )
        try:
            ET.parse(output_path)
        except ET.ParseError as error:
            raise RuntimeError(f"Invalid {entity_name} XML: {output_path}") from error
        output_paths.append(output_path)
    return output_paths


def parse_sample_xml(paths: list[Path]) -> dict[str, dict[str, Any]]:
    samples: dict[str, dict[str, Any]] = {}
    for path in paths:
        root = ET.parse(path).getroot()
        for sample in root.iter("SAMPLE"):
            accession = sample.get("accession", "")
            attributes: dict[str, list[str]] = defaultdict(list)
            for item in sample.findall(".//SAMPLE_ATTRIBUTE"):
                tag = (item.findtext("TAG") or "").strip()
                value = (item.findtext("VALUE") or "").strip()
                if tag and value:
                    attributes[tag].append(value)
            samples[accession] = {
                "alias": sample.get("alias", ""),
                "title": (sample.findtext("TITLE") or "").strip(),
                "attributes": dict(attributes),
                "cache_path": rel(path),
            }
    return samples


def parse_study_xml(paths: list[Path]) -> dict[str, dict[str, Any]]:
    studies: dict[str, dict[str, Any]] = {}
    for path in paths:
        root = ET.parse(path).getroot()
        for study in root.iter("STUDY"):
            accession = study.get("accession", "")
            external_ids: dict[str, list[str]] = defaultdict(list)
            for identifier in study.findall(".//EXTERNAL_ID"):
                namespace = identifier.get("namespace", "")
                value = (identifier.text or "").strip()
                if namespace and value:
                    external_ids[namespace].append(value)
            xrefs: dict[str, list[str]] = defaultdict(list)
            for link in study.findall(".//XREF_LINK"):
                database = (link.findtext("DB") or "").strip()
                identifier = (link.findtext("ID") or "").strip()
                if database and identifier:
                    xrefs[database].append(identifier)
            studies[accession] = {
                "alias": study.get("alias", ""),
                "title": (study.findtext("DESCRIPTOR/STUDY_TITLE") or "").strip(),
                "abstract": (study.findtext("DESCRIPTOR/STUDY_ABSTRACT") or "").strip(),
                "external_ids": dict(external_ids),
                "xrefs": dict(xrefs),
                "cache_path": rel(path),
            }
    return studies


def first_attribute(sample: dict[str, Any], names: Iterable[str]) -> str:
    normalized = {
        re.sub(r"[^a-z0-9]+", "_", key.lower()).strip("_"): values
        for key, values in sample.get("attributes", {}).items()
    }
    for name in names:
        values = normalized.get(name)
        if values:
            return ";".join(values)
    return ""


def sum_semicolon_ints(value: str) -> int:
    total = 0
    for item in value.split(";"):
        item = item.strip()
        if item.isdigit():
            total += int(item)
    return total


def choose(*values: str) -> str:
    return next((value.strip() for value in values if value and value.strip()), "")


def evidence_row(
    sample_uid: str,
    field: str,
    value: str,
    source_type: str,
    source_locator: str,
    status: str = "source_reported",
) -> dict[str, str] | None:
    if not value:
        return None
    return {
        "sample_uid": sample_uid,
        "field": field,
        "value": value,
        "source_type": source_type,
        "source_locator": source_locator,
        "evidence_status": status,
    }


def build_local_sources() -> tuple[list[dict[str, str]], dict[str, dict[str, str]]]:
    legacy = {
        row["sample_id"]: row
        for row in read_tsv(REPO_ROOT / "alphagenome_custom/metadata/track_metadata.tsv")
    }
    batch_518 = {
        row["sample_id"]: row
        for row in read_tsv(
            REPO_ROOT
            / "alphagenome_custom/metadata/bigwig_metadata_manifest_bw_2026.5.18.tsv"
        )
    }
    batch_comparison = {
        row["sample_id"]: row
        for row in read_tsv(SHAREPOINT_DIR.parent / "batch_comparison.tsv")
    }
    local_inventory = {
        Path(row["remote_path"]).stem: row
        for row in read_tsv(SHAREPOINT_DIR.parent / "local_inventory.tsv")
        if row["remote_path"].endswith(".bw") and "/" not in row["remote_path"]
    }
    sharepoint_439 = read_xlsx_first_sheet(
        SHAREPOINT_DIR / "c_elegans_rnaseq_table.xlsx"
    )
    samples_sheet = {
        row.get("ID", ""): row
        for row in read_xlsx_first_sheet(SHAREPOINT_DIR / "Samples.xlsx")
        if row.get("ID", "").startswith("SRR")
    }

    rows: list[dict[str, str]] = []
    for accession in sorted(legacy):
        comparison = batch_comparison[accession]
        rows.append(
            {
                "run_accession": accession,
                "batch": "2026.4.20",
                "local_path": comparison["local_path"],
                "expected_size_bytes": comparison["size_bytes"],
                "expected_sha256": comparison["sha256"],
                "wbcel235_exact": comparison["wbcel235_exact"],
                "source_stage": legacy[accession].get("stage", ""),
                "source_sex": legacy[accession].get("sex", ""),
                "source_tissue": legacy[accession].get("biosample_name", ""),
                "source_condition": "",
                "source_data_source": legacy[accession].get("data_source", ""),
                "source_genotype": "",
                "source_assay": legacy[accession].get("assay", ""),
                "source_strand": legacy[accession].get("strand", ""),
                "source_project": "",
            }
        )
    for accession in sorted(batch_518):
        comparison = batch_comparison[accession]
        source = batch_518[accession]
        rows.append(
            {
                "run_accession": accession,
                "batch": "2026.5.18",
                "local_path": comparison["local_path"],
                "expected_size_bytes": comparison["size_bytes"],
                "expected_sha256": comparison["sha256"],
                "wbcel235_exact": comparison["wbcel235_exact"],
                "source_stage": source.get("stage", ""),
                "source_sex": source.get("sex", ""),
                "source_tissue": source.get("tissue", ""),
                "source_condition": source.get("condition", ""),
                "source_data_source": source.get("data_source", ""),
                "source_genotype": "N2",
                "source_assay": source.get("assay", ""),
                "source_strand": source.get("strand", ""),
                "source_project": "GSE49043",
            }
        )
    for source in sharepoint_439:
        accession = source["SRR_ID"]
        inventory = local_inventory[accession]
        rows.append(
            {
                "run_accession": accession,
                "batch": "sharepoint_2026-07-13_new439",
                "local_path": rel(SHAREPOINT_DIR / f"{accession}.bw"),
                "expected_size_bytes": inventory["local_size_bytes"],
                "expected_sha256": inventory["sha256"],
                "wbcel235_exact": inventory["wbcel235_exact"],
                "source_stage": source.get("Stage", ""),
                "source_sex": "",
                "source_tissue": "",
                "source_condition": "",
                "source_data_source": source.get("GEO Dataset", ""),
                "source_genotype": source.get("genotype", ""),
                "source_assay": source.get("Library_Type", ""),
                "source_strand": "",
                "source_project": source.get("GEO Dataset", ""),
            }
        )
    return rows, samples_sheet


def build_manifests() -> tuple[list[dict[str, Any]], list[dict[str, str]], dict[str, Any]]:
    local_rows, samples_sheet = build_local_sources()
    accessions = sorted(row["run_accession"] for row in local_rows)
    if len(accessions) != 485 or len(set(accessions)) != 485:
        raise RuntimeError(
            f"Expected 485 unique accessions, got {len(accessions)} / {len(set(accessions))}"
        )

    ena_runs, ena_failures = fetch_ena_runs(accessions)
    sample_accessions = sorted(
        {row["sample_accession"] for row in ena_runs.values() if row["sample_accession"]}
    )
    sample_xml_paths = fetch_xml_chunks(sample_accessions, "ena_samples", "sample")
    ena_samples = parse_sample_xml(sample_xml_paths)
    study_accessions = sorted(
        {
            row["secondary_study_accession"]
            for row in ena_runs.values()
            if row["secondary_study_accession"]
        }
    )
    study_xml_paths = fetch_xml_chunks(study_accessions, "ena_studies", "study")
    ena_studies = parse_study_xml(study_xml_paths)

    evidence: list[dict[str, str]] = []
    manifest: list[dict[str, Any]] = []
    raw_before: dict[str, str] = {}
    for local in local_rows:
        path = REPO_ROOT / local["local_path"]
        raw_before[local["run_accession"]] = sha256(path)

    for local in sorted(local_rows, key=lambda row: row["run_accession"]):
        accession = local["run_accession"]
        sample_uid = f"rna_seq:{accession}"
        path = REPO_ROOT / local["local_path"]
        actual_sha = raw_before[accession]
        ena = ena_runs.get(accession, {})
        ena_sample = ena_samples.get(ena.get("sample_accession", ""), {})
        study = ena_studies.get(ena.get("secondary_study_accession", ""), {})
        sheet = samples_sheet.get(accession, {})

        sample_stage = first_attribute(
            ena_sample, ("developmental_stage", "development_stage", "dev_stage", "stage", "age")
        )
        sample_tissue = first_attribute(
            ena_sample, ("tissue", "tissue_type", "cell_type", "sample_type")
        )
        sample_sex = first_attribute(ena_sample, ("sex", "host_sex"))
        sample_strain = first_attribute(
            ena_sample, ("strain", "sub_strain", "genotype", "host_genotype")
        )
        sample_replicate = first_attribute(
            ena_sample,
            ("biological_replicate", "biological_rep", "replicate", "replicate_number"),
        )
        sample_condition = first_attribute(
            ena_sample, ("condition", "treatment", "experimental_factor")
        )
        geo_ids = study.get("external_ids", {}).get("GEO", [])
        pubmed_ids = study.get("xrefs", {}).get("pubmed", [])
        fastq_bytes = sum_semicolon_ints(ena.get("fastq_bytes", ""))
        bam_bytes = sum_semicolon_ints(ena.get("bam_bytes", ""))
        read_count = int(ena["read_count"]) if ena.get("read_count", "").isdigit() else 0
        base_count = int(ena["base_count"]) if ena.get("base_count", "").isdigit() else 0
        layout = ena.get("library_layout", "")
        layout_divisor = 2 if layout == "PAIRED" else 1
        mean_read_length = (
            base_count / read_count / layout_divisor if read_count and base_count else 0.0
        )

        stage = choose(sheet.get("Stage", ""), local["source_stage"], sample_stage, ena.get("dev_stage", ""))
        tissue = choose(
            sheet.get("Tissue", ""),
            local["source_tissue"],
            ena.get("cell_type", ""),
            ena.get("tissue_type", ""),
            sample_tissue,
        )
        sex = choose(sheet.get("Sex", ""), local["source_sex"], ena.get("sex", ""), sample_sex)
        strain = choose(ena.get("strain", ""), sample_strain, local["source_genotype"])
        genotype = choose(local["source_genotype"], ena.get("sub_strain", ""), ena.get("host_genotype", ""))
        condition = choose(sheet.get("Condition", ""), local["source_condition"], sample_condition)
        assay = choose(ena.get("library_strategy", ""), local["source_assay"])
        strandedness = choose(ena.get("read_strand", ""), local["source_strand"])
        if strandedness == ".":
            strandedness = "unknown"

        flags: list[str] = []
        if local["batch"] == "2026.5.18" and sheet.get("Paper/Lab", ""):
            if sheet["Paper/Lab"] not in local["source_data_source"]:
                flags.append("data_source_local_vs_sharepoint")
        if accession in {"SRR941651", "SRR941697"}:
            flags.append("legacy_T38_collapses_distinct_sharepoint_stages")
        if local["source_project"] and ena.get("study_accession"):
            if local["source_project"].startswith("PRJ") and local["source_project"] != ena["study_accession"]:
                flags.append("provided_project_vs_ena_project")
        if accession not in ena_runs:
            flags.append("ena_run_query_failed")

        normalization_class = "C" if ena.get("fastq_ftp") or ena.get("bam_ftp") else "D"
        exclusion_reason = (
            "awaiting_uniform_reprocessing_signal_unit_unknown"
            if normalization_class == "C"
            else "signal_unit_unknown_and_no_source_reads_located"
        )
        query_status = "success" if accession in ena_runs else "not_found"
        source_locator = (
            f"ENA read_run:{accession}; sample:{ena.get('sample_accession', '')}; "
            f"study:{ena.get('secondary_study_accession', '')}"
        )

        source_values = [
            ("stage", local["source_stage"], "local_project_metadata", local["batch"]),
            ("stage", sheet.get("Stage", ""), "sharepoint_samples_xlsx", "Samples.xlsx"),
            ("stage", sample_stage, "ena_sample", ena_sample.get("cache_path", "")),
            ("tissue", local["source_tissue"], "local_project_metadata", local["batch"]),
            ("tissue", sheet.get("Tissue", ""), "sharepoint_samples_xlsx", "Samples.xlsx"),
            ("tissue", sample_tissue, "ena_sample", ena_sample.get("cache_path", "")),
            ("sex", local["source_sex"], "local_project_metadata", local["batch"]),
            ("sex", sheet.get("Sex", ""), "sharepoint_samples_xlsx", "Samples.xlsx"),
            ("sex", sample_sex, "ena_sample", ena_sample.get("cache_path", "")),
            ("strain", sample_strain, "ena_sample", ena_sample.get("cache_path", "")),
            ("genotype", local["source_genotype"], "provided_batch_table", local["batch"]),
            ("condition", condition, "resolved_source_reported", source_locator),
            ("biological_replicate", sample_replicate, "ena_sample", ena_sample.get("cache_path", "")),
            ("data_source", local["source_data_source"], "provided_batch_table", local["batch"]),
            ("data_source", sheet.get("Paper/Lab", ""), "sharepoint_samples_xlsx", "Samples.xlsx"),
            ("study_title", study.get("title", ""), "ena_study", study.get("cache_path", "")),
            ("library_strategy", ena.get("library_strategy", ""), "ena_read_run", source_locator),
            ("library_layout", layout, "ena_read_run", source_locator),
            ("strandedness", ena.get("read_strand", ""), "ena_read_run", source_locator),
        ]
        for field, value, source_type, locator in source_values:
            row = evidence_row(sample_uid, field, value, source_type, locator)
            if row:
                evidence.append(row)

        manifest.append(
            {
                "sample_uid": sample_uid,
                "batch": local["batch"],
                "run_accession": accession,
                "experiment_accession": ena.get("experiment_accession", ""),
                "biosample_accession": ena.get("sample_accession", ""),
                "sra_sample_accession": ena.get("secondary_sample_accession", ""),
                "bioproject_accession": ena.get("study_accession", ""),
                "sra_study_accession": ena.get("secondary_study_accession", ""),
                "geo_accession": ";".join(geo_ids),
                "pubmed_id": ";".join(pubmed_ids),
                "local_path": local["local_path"],
                "local_size_bytes": path.stat().st_size,
                "sha256": actual_sha,
                "sha256_matches_prior_inventory": str(actual_sha == local["expected_sha256"]),
                "wbcel235_exact": local["wbcel235_exact"],
                "organism": ena.get("scientific_name", "Caenorhabditis elegans"),
                "strain": strain,
                "genotype": genotype,
                "sex": sex,
                "tissue_or_cell_type": tissue,
                "development_stage": stage,
                "condition": condition,
                "assay": assay,
                "library_source": ena.get("library_source", ""),
                "library_selection": ena.get("library_selection", ""),
                "library_layout": layout,
                "mean_read_length_bp_derived": f"{mean_read_length:.6g}" if mean_read_length else "",
                "read_length_status": "derived_from_ENA_base_count_read_count_layout" if mean_read_length else "unknown",
                "strandedness": strandedness or "unknown",
                "biological_replicate": sample_replicate,
                "technical_replicate": "unknown",
                "study_title": study.get("title", ena.get("study_title", "")),
                "provided_data_source": local["source_data_source"],
                "sharepoint_paper_lab": sheet.get("Paper/Lab", ""),
                "provided_stage": local["source_stage"],
                "sharepoint_stage": sheet.get("Stage", ""),
                "ena_sample_stage": sample_stage,
                "provided_tissue": local["source_tissue"],
                "ena_sample_tissue": sample_tissue,
                "bigwig_generator": "unknown_from_provided_metadata",
                "bigwig_reference_genome": "WBcel235_by_chromosome_header_validation",
                "bigwig_aligner": "unknown_from_provided_metadata",
                "bigwig_generation_command": "unknown_from_provided_metadata",
                "signal_unit": "unknown",
                "normalization_class": normalization_class,
                "mathematically_convertible_without_source_reads": "False",
                "read_count": read_count,
                "base_count": base_count,
                "fastq_available": str(bool(ena.get("fastq_ftp"))),
                "fastq_md5": ena.get("fastq_md5", ""),
                "fastq_file_bytes": ena.get("fastq_bytes", ""),
                "fastq_bytes": fastq_bytes,
                "fastq_ftp": ena.get("fastq_ftp", ""),
                "bam_available": str(bool(ena.get("bam_ftp"))),
                "bam_bytes": bam_bytes,
                "bam_ftp": ena.get("bam_ftp", ""),
                "ena_query_status": query_status,
                "metadata_conflict_flags": ";".join(flags),
                "metadata_evidence_status": "source_reported_pending_harmonization",
                "include_formal_v2": "False",
                "exclusion_reason": exclusion_reason,
            }
        )

    raw_after = {
        row["run_accession"]: sha256(REPO_ROOT / row["local_path"])
        for row in local_rows
    }
    counts = defaultdict(int)
    for row in manifest:
        counts[row["normalization_class"]] += 1
    summary = {
        "schema_version": 1,
        "created_at": utc_now(),
        "sample_count": len(manifest),
        "unique_run_accessions": len({row["run_accession"] for row in manifest}),
        "batch_counts": {
            batch: sum(row["batch"] == batch for row in manifest)
            for batch in sorted({row["batch"] for row in manifest})
        },
        "normalization_class_counts": {
            key: counts.get(key, 0) for key in ("A", "B", "C", "D")
        },
        "ena_success_count": sum(row["ena_query_status"] == "success" for row in manifest),
        "ena_failure_accessions": ena_failures,
        "fastq_total_bytes": sum(int(row["fastq_bytes"]) for row in manifest),
        "bam_total_bytes_available": sum(int(row["bam_bytes"]) for row in manifest),
        "total_read_count": sum(int(row["read_count"]) for row in manifest),
        "total_base_count": sum(int(row["base_count"]) for row in manifest),
        "raw_hashes_unchanged_during_p1": raw_before == raw_after,
        "formal_v2_included_count": sum(row["include_formal_v2"] == "True" for row in manifest),
        "conflict_flagged_count": sum(bool(row["metadata_conflict_flags"]) for row in manifest),
    }
    return manifest, evidence, summary


def write_approval_request(summary: dict[str, Any]) -> None:
    gib = summary["fastq_total_bytes"] / 1024**3
    working_gib = gib * 4
    class_counts = summary["normalization_class_counts"]
    path = REPO_ROOT / "docs/v2_g1_source_reprocessing_request.md"
    path.write_text(
        f"""# G1 Source-read Reprocessing Approval Packet

Status: **P1 PLANNING ESTIMATE - CURRENT AUTHORIZATION IS RECORDED IN `execution_state.json`**

Generated: `{summary['created_at']}`

## Classification

| Class | Samples |
| --- | ---: |
| A: confirmed per-base RPM | {class_counts['A']} |
| B: known convertible unit | {class_counts['B']} |
| C: unknown bigWig unit with source reads | {class_counts['C']} |
| D: unknown unit without located source reads | {class_counts['D']} |

## Located Source Volume

- ENA run metadata found: `{summary['ena_success_count']}` / `{summary['sample_count']}`.
- Compressed FASTQ bytes reported by ENA: `{summary['fastq_total_bytes']}` ({gib:.2f} GiB).
- Total submitted sequence bases: `{summary['total_base_count']}`.
- Planning-only working-space allowance at 4x compressed FASTQ: approximately `{working_gib:.2f} GiB`.

The 4x allowance is a storage planning estimate, not a measured pipeline requirement. CPU/GPU time is intentionally not fabricated: a representative single-end/paired-end pilot must measure alignment and bigWig-generation throughput before a bulk compute estimate is approved.

## Proposed Target Locations

- Source FASTQ: `shared/source_reads/v2/fastq/`
- Alignment outputs: `alphagenome_custom/tracks/rna_seq_v2_alignment/`
- Normalized bigWigs: `alphagenome_custom/tracks/rna_seq_v2_normalized/`

All locations are ignored by Git. No download or realignment command has been run. The exact downloader, aligner, strandedness handling, multimapping policy, and per-base coverage command will be locked after the metadata/group review and a small benchmark design.
"""
    )


def main() -> None:
    METADATA_DIR.mkdir(parents=True, exist_ok=True)
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    manifest, evidence, summary = build_manifests()
    manifest_fields = list(manifest[0])
    write_tsv(METADATA_DIR / "rna_seq_samples_v2.tsv", manifest, manifest_fields)
    source_fields = [
        "run_accession",
        "experiment_accession",
        "batch",
        "library_layout",
        "library_selection",
        "read_count",
        "base_count",
        "fastq_ftp",
        "fastq_md5",
        "fastq_bytes",
    ]
    full_sources = []
    for row in manifest:
        if row["assay"] != "RNA-Seq":
            continue
        source_row = {field: row[field] for field in source_fields}
        source_row["fastq_bytes"] = row["fastq_file_bytes"]
        full_sources.append(source_row)
    write_tsv(
        METADATA_DIR / "p3_full_sources.tsv",
        full_sources,
        source_fields,
    )
    write_tsv(
        METADATA_DIR / "rna_seq_metadata_evidence_v2.tsv",
        evidence,
        [
            "sample_uid",
            "field",
            "value",
            "source_type",
            "source_locator",
            "evidence_status",
        ],
    )
    (METADATA_DIR / "p1_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n"
    )
    write_approval_request(summary)
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
