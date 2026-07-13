#!/usr/bin/env python3
"""Build conservative v2 RNA-seq groups, decisions, ontology, and replicate QC."""

from __future__ import annotations

from collections import defaultdict
import csv
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import re
from typing import Any, Iterable
from xml.etree import ElementTree as ET

import numpy as np
import pyBigWig


REPO_ROOT = Path(__file__).resolve().parents[1]
METADATA_DIR = REPO_ROOT / "alphagenome_custom/metadata/v2"
CACHE_DIR = REPO_ROOT / "shared/metadata_sources/v2"
PROFILE_CACHE_DIR = CACHE_DIR / "replicate_qc_profiles"
GTF_PATH = (
    REPO_ROOT / "alphagenome_custom/reference/Caenorhabditis_elegans.WBcel235.115.gtf"
)
WBBT_PATH = CACHE_DIR / "ontologies/wbbt.obo"
WBLS_PATH = CACHE_DIR / "ontologies/wbls.obo"

DUPLICATE_CANONICAL = {
    "SRR13964554": "SRR12483327",
    "SRR13964555": "SRR12483328",
    "SRR13964556": "SRR12483329",
}


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


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


def parse_obo(path: Path) -> dict[str, str]:
    terms: dict[str, str] = {}
    term_id = ""
    term_name = ""
    obsolete = False
    for line in path.read_text().splitlines() + ["[Term]"]:
        if line == "[Term]":
            if term_id and term_name and not obsolete:
                terms[term_id] = term_name
            term_id = ""
            term_name = ""
            obsolete = False
        elif line.startswith("id: "):
            term_id = line[4:]
        elif line.startswith("name: "):
            term_name = line[6:]
        elif line == "is_obsolete: true":
            obsolete = True
    return terms


def stage_mapping(value: str) -> tuple[str, str, str]:
    normalized = value.strip().lower()
    if normalized.startswith("embryo"):
        return "WBls:0000003", "embryo Ce", "supported_broad_stage"
    if normalized.startswith("l1-l3"):
        return "WBls:0000023", "larva Ce", "supported_broad_stage"
    if normalized == "l1":
        return "WBls:0000024", "L1 larva Ce", "exact_stage"
    if normalized == "l2":
        return "WBls:0000027", "L2 larva Ce", "exact_stage"
    if normalized == "l3":
        return "WBls:0000035", "L3 larva Ce", "exact_stage"
    if normalized == "l4" or normalized.startswith("l4("):
        return "WBls:0000038", "L4 larva Ce", "exact_or_timed_stage"
    if normalized == "young adult":
        return "WBls:0000041", "adult Ce", "supported_broad_stage"
    return "", "", "unmapped"


def anatomy_mapping(value: str) -> tuple[str, str, str]:
    normalized = value.strip().lower()
    if not normalized:
        return "", "", "missing_source_anatomy"
    if "seam" in normalized:
        return "WBbt:0005753", "seam cell", "supported_anatomy"
    if "body wall muscle" in normalized:
        return "WBbt:0006804", "body wall muscle cell", "supported_anatomy"
    if normalized.startswith("muscle") or normalized == "muscle":
        return "WBbt:0003675", "muscle cell", "supported_broad_anatomy"
    if "intestine" in normalized:
        return "WBbt:0005772", "intestine", "supported_anatomy"
    if "pharynx" in normalized:
        return "WBbt:0003681", "pharynx", "supported_anatomy"
    whole_terms = (
        "whole animal",
        "whole worm",
        "whole organism",
        "early embryo",
        "embryo",
        "embryos",
        "embryonic total rna",
        "whole animal l1 larvae",
    )
    if normalized in whole_terms:
        return (
            "WBbt:0000100",
            "C. elegans anatomical entity",
            "supported_root_for_whole_organism_sample",
        )
    return "", "", "unmapped"


def context_tuple(
    row: dict[str, str], anatomy_curie: str, stage_curie: str
) -> tuple[str, ...]:
    return (
        row["sra_study_accession"],
        row["assay"],
        row["strandedness"],
        row["tissue_or_cell_type"],
        anatomy_curie,
        row["development_stage"],
        stage_curie,
        row["sex"],
        row["strain"],
        row["genotype"],
        row["condition"],
    )


def load_ena_sample_contexts() -> dict[str, dict[str, str]]:
    contexts: dict[str, dict[str, str]] = {}
    for path in sorted((CACHE_DIR / "ena_samples").glob("*.xml")):
        root = ET.parse(path).getroot()
        for sample in root.iter("SAMPLE"):
            attributes: dict[str, list[str]] = defaultdict(list)
            for item in sample.findall(".//SAMPLE_ATTRIBUTE"):
                tag = (item.findtext("TAG") or "").strip()
                value = (item.findtext("VALUE") or "").strip()
                if tag and value:
                    normalized_tag = re.sub(
                        r"[^a-z0-9]+", "_", tag.lower()
                    ).strip("_")
                    attributes[normalized_tag].append(value)

            def first(*names: str) -> str:
                for name in names:
                    if attributes.get(name):
                        return ";".join(attributes[name])
                return ""

            alias = sample.get("alias", "")
            title = (sample.findtext("TITLE") or "").strip()
            replicate = first("biological_replicate", "replicate", "replicate_number")
            replicate_evidence = "ENA_sample_attribute" if replicate else ""
            if not replicate:
                match = re.search(
                    r"(?i)(?:biol(?:ogical)?[ _-]*)?rep(?:licate)?[ _-]*(\d+)",
                    f"{title} {alias}",
                )
                if match:
                    replicate = match.group(1)
                    replicate_evidence = "ENA_sample_title_or_alias"
            gsm = alias if alias.startswith("GSM") else ""
            contexts[sample.get("accession", "")] = {
                "sample_alias": alias,
                "sample_title": title,
                "geo_sample_accession": gsm,
                "source_name": first("source_name", "isolation_source"),
                "sample_tissue": first("tissue", "tissue_type", "cell_type"),
                "sample_stage": first(
                    "developmental_stage",
                    "developemental_stage",
                    "dev_stage",
                    "stage",
                ),
                "sample_sex": first("sex"),
                "sample_strain": first("strain", "strain_name", "strain_genotype"),
                "sample_genotype": first(
                    "genotype",
                    "genotype_variation",
                    "strain_genotype",
                    "genetic_background",
                ),
                "sample_condition": first("treatment", "condition"),
                "biological_replicate_resolved": replicate,
                "biological_replicate_evidence": replicate_evidence,
                "sample_context_cache": str(path.relative_to(REPO_ROOT)),
            }
    return contexts


def tissue_from_source_name(value: str) -> str:
    normalized = value.lower()
    if not normalized:
        return ""
    if "seamcell" in normalized or "seam cell" in normalized:
        return "seam cell"
    if "muscle" in normalized and "non-muscle" not in normalized:
        return "muscle"
    whole_markers = (
        "whole animal",
        "whole-animal",
        "whole worm",
        "whole organism",
        "entire young adult",
        "young adult",
        "larvae",
        "embryo",
    )
    if any(marker in normalized for marker in whole_markers):
        return "whole animal"
    return ""


def build_sample_context_rows(samples: list[dict[str, str]]) -> list[dict[str, str]]:
    ena_contexts = load_ena_sample_contexts()
    experiment_counts: dict[str, int] = defaultdict(int)
    for sample in samples:
        experiment_counts[sample["experiment_accession"]] += 1
    rows = []
    for sample in samples:
        ena = ena_contexts.get(sample["biosample_accession"], {})
        tissue = sample["tissue_or_cell_type"] or ena.get("sample_tissue", "")
        if not tissue:
            tissue = tissue_from_source_name(ena.get("source_name", ""))
        replicate = (
            sample["biological_replicate"]
            or ena.get("biological_replicate_resolved", "")
        )
        replicate_evidence = (
            "ENA_sample_attribute"
            if sample["biological_replicate"]
            else ena.get("biological_replicate_evidence", "")
        )
        technical_count = experiment_counts[sample["experiment_accession"]]
        rows.append(
            {
                "sample_uid": sample["sample_uid"],
                "run_accession": sample["run_accession"],
                "experiment_accession": sample["experiment_accession"],
                "biosample_accession": sample["biosample_accession"],
                "sample_alias": ena.get("sample_alias", ""),
                "sample_title": ena.get("sample_title", ""),
                "geo_sample_accession": ena.get("geo_sample_accession", ""),
                "source_name": ena.get("source_name", ""),
                "sra_study_accession": sample["sra_study_accession"],
                "bioproject_accession": sample["bioproject_accession"],
                "geo_accession": sample["geo_accession"],
                "assay": sample["assay"],
                "strandedness": sample["strandedness"],
                "tissue_or_cell_type": tissue,
                "development_stage": sample["development_stage"]
                or ena.get("sample_stage", ""),
                "sex": sample["sex"] or ena.get("sample_sex", ""),
                "strain": sample["strain"] or ena.get("sample_strain", ""),
                "genotype": ena.get("sample_genotype", "") or sample["genotype"],
                "condition": sample["condition"] or ena.get("sample_condition", ""),
                "biological_replicate": replicate,
                "biological_replicate_evidence": replicate_evidence,
                "technical_unit_id": sample["experiment_accession"],
                "technical_run_count": str(technical_count),
                "technical_replicate_status": (
                    "same_ENA_experiment_multiple_runs"
                    if technical_count > 1
                    else "single_run_for_ENA_experiment"
                ),
                "local_path": sample["local_path"],
                "sha256": sample["sha256"],
                "technical_replicate": sample["technical_replicate"],
                "context_evidence": ena.get("sample_context_cache", ""),
            }
        )
    return rows


def parse_gene_intervals() -> dict[str, list[tuple[int, int]]]:
    intervals: dict[str, list[tuple[int, int]]] = defaultdict(list)
    with GTF_PATH.open() as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 5 and fields[2] == "gene" and fields[0] in {
                "I",
                "II",
                "III",
                "IV",
                "V",
                "X",
            }:
                intervals[fields[0]].append((int(fields[3]) - 1, int(fields[4])))
    for chromosome in intervals:
        intervals[chromosome].sort()
    return intervals


def mean_bins(values: np.ndarray, bin_size: int = 1000) -> np.ndarray:
    size = values.size
    full_bins = size // bin_size
    outputs = []
    if full_bins:
        outputs.append(values[: full_bins * bin_size].reshape(full_bins, bin_size).mean(axis=1))
    if full_bins * bin_size < size:
        outputs.append(np.array([values[full_bins * bin_size :].mean()], dtype=np.float32))
    return np.concatenate(outputs).astype(np.float32, copy=False)


def signal_profile(
    row: dict[str, str], gene_intervals: dict[str, list[tuple[int, int]]]
) -> dict[str, Any]:
    accession = row["run_accession"]
    cache_path = PROFILE_CACHE_DIR / f"{accession}.npz"
    if cache_path.is_file():
        with np.load(cache_path, allow_pickle=False) as data:
            if str(data["source_sha256"]) == row["sha256"]:
                return {
                    "bin_means": data["bin_means"],
                    "gene_means": data["gene_means"],
                    "zero_fraction": float(data["zero_fraction"]),
                    "sum_signal": float(data["sum_signal"]),
                    "nonzero_mean": float(data["nonzero_mean"]),
                }

    PROFILE_CACHE_DIR.mkdir(parents=True, exist_ok=True)
    bin_parts: list[np.ndarray] = []
    gene_values: list[float] = []
    zero_count = 0
    nonzero_count = 0
    total_count = 0
    sum_signal = 0.0
    with pyBigWig.open(str(REPO_ROOT / row["local_path"])) as bigwig:
        chromosomes = bigwig.chroms()
        for chromosome in ("I", "II", "III", "IV", "V", "X"):
            length = chromosomes[chromosome]
            values = bigwig.values(chromosome, 0, length, numpy=True)
            values = np.nan_to_num(
                values, nan=0.0, posinf=0.0, neginf=0.0
            ).astype(np.float32, copy=False)
            bin_parts.append(mean_bins(values))
            zero_count += int(np.count_nonzero(values == 0))
            nonzero_count += int(np.count_nonzero(values > 0))
            total_count += values.size
            sum_signal += float(values.sum(dtype=np.float64))
            gene_values.extend(
                float(values[start:end].mean(dtype=np.float64))
                for start, end in gene_intervals.get(chromosome, [])
                if end > start
            )
    bin_means = np.concatenate(bin_parts)
    gene_means = np.asarray(gene_values, dtype=np.float32)
    zero_fraction = zero_count / max(total_count, 1)
    nonzero_mean = sum_signal / max(nonzero_count, 1)
    np.savez_compressed(
        cache_path,
        source_sha256=np.array(row["sha256"]),
        bin_means=bin_means,
        gene_means=gene_means,
        zero_fraction=np.array(zero_fraction),
        sum_signal=np.array(sum_signal),
        nonzero_mean=np.array(nonzero_mean),
    )
    return {
        "bin_means": bin_means,
        "gene_means": gene_means,
        "zero_fraction": zero_fraction,
        "sum_signal": sum_signal,
        "nonzero_mean": nonzero_mean,
    }


def safe_pearson(first: np.ndarray, second: np.ndarray) -> float:
    if first.size != second.size or first.size == 0:
        return float("nan")
    if float(first.std()) == 0.0 or float(second.std()) == 0.0:
        return float("nan")
    return float(np.corrcoef(first, second)[0, 1])


def build_groups(
    samples: list[dict[str, str]],
    sample_contexts: list[dict[str, str]],
    ontology_rows: list[dict[str, str]],
    *,
    exclude_duplicate_current_signals: bool = True,
) -> tuple[
    list[dict[str, Any]],
    list[dict[str, Any]],
    list[dict[str, Any]],
    list[dict[str, Any]],
]:
    ontology_by_key = {
        (row["ontology_type"], row["source_value"]): row for row in ontology_rows
    }
    context_by_uid = {row["sample_uid"]: row for row in sample_contexts}
    eligible: list[dict[str, str]] = []
    decisions: list[dict[str, Any]] = []
    exclusions: list[dict[str, Any]] = []
    for sample in samples:
        accession = sample["run_accession"]
        reasons: list[str] = []
        current_signal_eligible = True
        if sample["assay"] != "RNA-Seq":
            current_signal_eligible = False
            reasons.append("non_rna_seq_assay")
            exclusions.append(
                {
                    "sample_uid": sample["sample_uid"],
                    "scope": "RNA_SEQ_v2",
                    "decision": "exclude",
                    "reason": "ENA library strategy is ChIP-Seq",
                    "evidence": f"ENA run {accession}; assay={sample['assay']}",
                    "review_status": "resolved",
                }
            )
        if exclude_duplicate_current_signals and accession in DUPLICATE_CANONICAL:
            current_signal_eligible = False
            canonical = DUPLICATE_CANONICAL[accession]
            reasons.append(f"duplicate_current_bigwig_of_{canonical}")
            exclusions.append(
                {
                    "sample_uid": sample["sample_uid"],
                    "scope": "current_bigwig_signal",
                    "decision": "exclude_duplicate_signal",
                    "reason": f"byte-identical current bigWig; canonical={canonical}",
                    "evidence": f"SHA-256={sample['sha256']}",
                    "review_status": "resolved_source_reuse_requires_raw_read_review",
                }
            )
        exclusions.append(
            {
                "sample_uid": sample["sample_uid"],
                "scope": "formal_v2_until_reprocessing",
                "decision": "hold",
                "reason": sample["exclusion_reason"],
                "evidence": f"normalization_class={sample['normalization_class']}; signal_unit={sample['signal_unit']}",
                "review_status": "resolved_hold_pending_G1",
            }
        )
        decisions.append(
            {
                "sample_uid": sample["sample_uid"],
                "run_accession": accession,
                "assay": sample["assay"],
                "current_bigwig_group_eligible": str(current_signal_eligible),
                "canonical_signal_accession": DUPLICATE_CANONICAL.get(accession, accession),
                "formal_v2_eligible": "False",
                "decision_reasons": ";".join(reasons)
                if reasons
                else "awaiting_uniform_reprocessing",
            }
        )
        if current_signal_eligible and sample["assay"] == "RNA-Seq":
            eligible.append(context_by_uid[sample["sample_uid"]])

    runs_by_experiment: dict[str, list[dict[str, str]]] = defaultdict(list)
    for sample in eligible:
        runs_by_experiment[sample["experiment_accession"]].append(sample)

    biological_contexts: dict[
        tuple[str, ...], list[list[dict[str, str]]]
    ] = defaultdict(list)
    for experiment, runs in sorted(runs_by_experiment.items()):
        first = runs[0]
        anatomy = ontology_by_key[("anatomy", first["tissue_or_cell_type"])]
        stage = ontology_by_key[("life_stage", first["development_stage"])]
        base_context = context_tuple(first, anatomy["curie"], stage["curie"])
        if any(
            context_tuple(
                run,
                ontology_by_key[("anatomy", run["tissue_or_cell_type"])]["curie"],
                ontology_by_key[("life_stage", run["development_stage"])]["curie"],
            )
            != base_context
            for run in runs[1:]
        ):
            raise RuntimeError(f"Runs in {experiment} have inconsistent biological context")
        if first["biological_replicate"]:
            key = ("replicate_context", *base_context)
        else:
            key = ("singleton_biological_unit", *base_context, experiment)
        biological_contexts[key].append(runs)

    groups: list[dict[str, Any]] = []
    memberships: list[dict[str, Any]] = []
    for index, (key, units) in enumerate(sorted(biological_contexts.items()), start=1):
        group_id = f"RNA_V2_G{index:04d}"
        members = [run for unit in units for run in unit]
        first = members[0]
        anatomy = ontology_by_key[("anatomy", first["tissue_or_cell_type"])]
        stage = ontology_by_key[("life_stage", first["development_stage"])]
        ontology_reasons = []
        if not anatomy["curie"]:
            ontology_reasons.append(f"anatomy:{anatomy['mapping_status']}")
        if not stage["curie"]:
            ontology_reasons.append(f"life_stage:{stage['mapping_status']}")
        if len(units) > 1:
            evidence_types = {
                unit[0]["biological_replicate_evidence"] for unit in units
            }
            relation = (
                "source_attribute_biological_replicates"
                if evidence_types == {"ENA_sample_attribute"}
                else "source_title_or_alias_biological_replicates"
            )
        elif len(members) > 1:
            relation = "single_biological_sample_with_technical_runs"
        else:
            relation = "singleton_no_replicate_aggregation"
        groups.append(
            {
                "group_id": group_id,
                "group_status": "candidate_after_uniform_reprocessing",
                "sra_study_accession": first["sra_study_accession"],
                "bioproject_accession": first["bioproject_accession"],
                "geo_accession": first["geo_accession"],
                "assay": first["assay"],
                "strand": first["strandedness"],
                "tissue_or_cell_type": first["tissue_or_cell_type"],
                "anatomy_curie": anatomy["curie"],
                "anatomy_label": anatomy["ontology_label"],
                "anatomy_mapping_status": anatomy["mapping_status"],
                "development_stage": first["development_stage"],
                "life_stage_curie": stage["curie"],
                "life_stage_label": stage["ontology_label"],
                "life_stage_mapping_status": stage["mapping_status"],
                "sex": first["sex"],
                "strain": first["strain"],
                "genotype": first["genotype"],
                "condition": first["condition"],
                "replicate_relation_status": relation,
                "n_members": len(members),
                "n_runs": len(members),
                "n_biological_units": len(units),
                "sample_ids": ",".join(
                    sorted(member["run_accession"] for member in members)
                ),
                "ontology_exclusion_reason": ";".join(ontology_reasons),
                "include_formal_v2": "False",
                "formal_block_reason": "normalization_class_C_requires_uniform_reprocessing",
            }
        )
        for member in sorted(members, key=lambda row: row["run_accession"]):
            memberships.append(
                {
                    "group_id": group_id,
                    "sample_uid": member["sample_uid"],
                    "run_accession": member["run_accession"],
                    "local_path": member["local_path"],
                    "sha256": member["sha256"],
                    "biological_replicate": member["biological_replicate"],
                    "technical_replicate": member["technical_replicate"],
                    "biological_replicate_evidence": member[
                        "biological_replicate_evidence"
                    ],
                    "technical_unit_id": member["technical_unit_id"],
                    "technical_run_count": member["technical_run_count"],
                    "technical_replicate_status": member[
                        "technical_replicate_status"
                    ],
                    "membership_status": "candidate_after_uniform_reprocessing",
                }
            )
    return groups, memberships, decisions, exclusions


def build_ontology_rows(samples: list[dict[str, str]]) -> list[dict[str, str]]:
    wbbt = parse_obo(WBBT_PATH)
    wbls = parse_obo(WBLS_PATH)
    rows: list[dict[str, str]] = []
    for ontology_type, values, mapper, terms, path in (
        (
            "anatomy",
            sorted({row["tissue_or_cell_type"] for row in samples}),
            anatomy_mapping,
            wbbt,
            WBBT_PATH,
        ),
        (
            "life_stage",
            sorted({row["development_stage"] for row in samples}),
            stage_mapping,
            wbls,
            WBLS_PATH,
        ),
    ):
        for value in values:
            curie, label, status = mapper(value)
            if curie and terms.get(curie) != label:
                raise RuntimeError(
                    f"Ontology label mismatch for {curie}: {label!r} vs {terms.get(curie)!r}"
                )
            rows.append(
                {
                    "ontology_type": ontology_type,
                    "source_value": value,
                    "curie": curie,
                    "ontology_label": label,
                    "mapping_status": status,
                    "ontology_source": str(path.relative_to(REPO_ROOT)),
                    "ontology_sha256": sha256(path),
                }
            )
    return rows


def build_metadata_resolutions(samples: list[dict[str, str]]) -> list[dict[str, str]]:
    rows = [
        {
            "issue_id": "P2-T38",
            "scope": "SRR941651;SRR941697",
            "issue": "legacy local grouping assigned both runs to 38hr",
            "resolution": "keep as separate samples and groups: L1-L3(38H) vs L4(38H)",
            "authority": "SharePoint Samples.xlsx plus distinct ENA experiments/BioSamples",
            "status": "resolved",
        },
        {
            "issue_id": "P2-518-SOURCE",
            "scope": "2026.5.18 batch",
            "issue": "local author label differs from SharePoint author label",
            "resolution": "use GSE49043/SRP027592 study identity and PMID 24036951; preserve both labels as evidence",
            "authority": "ENA study record and PubMed identifier",
            "status": "resolved",
        },
    ]
    project_pairs = sorted(
        {
            (
                row["provided_data_source"],
                row["bioproject_accession"],
                row["sra_study_accession"],
                row["geo_accession"],
            )
            for row in samples
            if "provided_project_vs_ena_project" in row["metadata_conflict_flags"]
        }
    )
    for index, pair in enumerate(project_pairs, start=1):
        provided, bioproject, study, geo = pair
        rows.append(
            {
                "issue_id": f"P2-PROJECT-{index:02d}",
                "scope": provided,
                "issue": "provided project label is broader or differs from run-level INSDC project",
                "resolution": f"use run-level {bioproject}/{study}/{geo} for grouping; retain provided label as evidence",
                "authority": "ENA run and study records",
                "status": "resolved",
            }
        )
    for duplicate, canonical in sorted(DUPLICATE_CANONICAL.items()):
        rows.append(
            {
                "issue_id": f"P2-DUP-{duplicate}",
                "scope": f"{canonical};{duplicate}",
                "issue": "different accessions have byte-identical current bigWigs and identical read/base counts",
                "resolution": f"retain provenance for both; use {canonical} as canonical current signal and hold {duplicate}",
                "authority": "local SHA-256 plus ENA run counts",
                "status": "resolved_current_signal_raw_reads_require_G1_review",
            }
        )
    return rows


def build_replicate_qc(
    samples: list[dict[str, str]],
    groups: list[dict[str, Any]],
    memberships: list[dict[str, Any]],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    sample_by_uid = {row["sample_uid"]: row for row in samples}
    members_by_group: dict[str, list[dict[str, str]]] = defaultdict(list)
    membership_by_run = {row["run_accession"]: row for row in memberships}
    for membership in memberships:
        members_by_group[membership["group_id"]].append(
            sample_by_uid[membership["sample_uid"]]
        )
    gene_intervals = parse_gene_intervals()
    pair_rows: list[dict[str, Any]] = []
    summary_rows: list[dict[str, Any]] = []
    for group in groups:
        members = members_by_group[group["group_id"]]
        if len(members) < 2:
            continue
        profiles = {
            member["run_accession"]: signal_profile(member, gene_intervals)
            for member in members
        }
        group_pairs = []
        for first_index, first in enumerate(sorted(members, key=lambda row: row["run_accession"])):
            for second in sorted(members, key=lambda row: row["run_accession"])[first_index + 1 :]:
                first_profile = profiles[first["run_accession"]]
                second_profile = profiles[second["run_accession"]]
                bin_pearson = safe_pearson(
                    np.log1p(first_profile["bin_means"]),
                    np.log1p(second_profile["bin_means"]),
                )
                gene_pearson = safe_pearson(
                    np.log1p(first_profile["gene_means"]),
                    np.log1p(second_profile["gene_means"]),
                )
                sums = [first_profile["sum_signal"], second_profile["sum_signal"]]
                ratio = max(sums) / max(min(sums), 1e-12)
                zero_delta = abs(
                    first_profile["zero_fraction"] - second_profile["zero_fraction"]
                )
                review_flag = (
                    not math.isfinite(bin_pearson)
                    or not math.isfinite(gene_pearson)
                    or bin_pearson < 0.8
                    or gene_pearson < 0.8
                    or ratio > 2.0
                )
                row = {
                    "group_id": group["group_id"],
                    "sample_a": first["run_accession"],
                    "sample_b": second["run_accession"],
                    "pair_type": (
                        "technical_runs_same_experiment"
                        if membership_by_run[first["run_accession"]]["technical_unit_id"]
                        == membership_by_run[second["run_accession"]]["technical_unit_id"]
                        else "biological_replicate_runs"
                    ),
                    "log1p_1kb_pearson": f"{bin_pearson:.8g}",
                    "log1p_gene_body_pearson": f"{gene_pearson:.8g}",
                    "total_signal_ratio": f"{ratio:.8g}",
                    "zero_fraction_a": f"{first_profile['zero_fraction']:.8g}",
                    "zero_fraction_b": f"{second_profile['zero_fraction']:.8g}",
                    "zero_fraction_delta": f"{zero_delta:.8g}",
                    "review_flag": str(review_flag),
                    "review_action": "manual_review_no_automatic_exclusion"
                    if review_flag
                    else "none",
                }
                pair_rows.append(row)
                group_pairs.append(row)
        bin_values = [float(row["log1p_1kb_pearson"]) for row in group_pairs]
        gene_values = [float(row["log1p_gene_body_pearson"]) for row in group_pairs]
        summary_rows.append(
            {
                "group_id": group["group_id"],
                "n_members": len(members),
                "n_pairs": len(group_pairs),
                "min_log1p_1kb_pearson": f"{min(bin_values):.8g}",
                "mean_log1p_1kb_pearson": f"{np.mean(bin_values):.8g}",
                "min_log1p_gene_body_pearson": f"{min(gene_values):.8g}",
                "mean_log1p_gene_body_pearson": f"{np.mean(gene_values):.8g}",
                "flagged_pairs": sum(row["review_flag"] == "True" for row in group_pairs),
                "group_decision": "candidate_keep_pending_manual_qc_review"
                if any(row["review_flag"] == "True" for row in group_pairs)
                else "candidate_keep",
            }
        )
    return pair_rows, summary_rows


def stable_digest(payload: Any) -> str:
    return hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def main() -> None:
    samples = read_tsv(METADATA_DIR / "rna_seq_samples_v2.tsv")
    sample_contexts = build_sample_context_rows(samples)
    ontology_rows = build_ontology_rows(sample_contexts)
    groups, memberships, decisions, exclusions = build_groups(
        samples, sample_contexts, ontology_rows
    )
    second_groups, second_memberships, second_decisions, second_exclusions = build_groups(
        samples, sample_contexts, ontology_rows
    )
    deterministic_match = stable_digest(
        [groups, memberships, decisions, exclusions]
    ) == stable_digest(
        [second_groups, second_memberships, second_decisions, second_exclusions]
    )
    if not deterministic_match:
        raise RuntimeError("Group manifest construction is not deterministic")

    resolutions = build_metadata_resolutions(samples)
    replicate_pairs, replicate_summary = build_replicate_qc(
        samples, groups, memberships
    )
    scientific_review = []
    for group in groups:
        issues = []
        for field in ("sex", "strain", "genotype", "condition"):
            if not group[field]:
                issues.append(f"missing_{field}")
        if group["ontology_exclusion_reason"]:
            issues.append(group["ontology_exclusion_reason"])
        if group["n_members"] > 1:
            qc = next(
                (row for row in replicate_summary if row["group_id"] == group["group_id"]),
                None,
            )
            if qc and qc["flagged_pairs"]:
                issues.append("replicate_qc_flag")
        scientific_review.append(
            {
                "group_id": group["group_id"],
                "n_members": group["n_members"],
                "sample_ids": group["sample_ids"],
                "formal_status": "not_eligible_before_uniform_reprocessing",
                "review_issues": ";".join(issues),
                "review_outcome": "candidate_only_do_not_generate_formal_track",
            }
        )

    write_tsv(METADATA_DIR / "rna_seq_groups_v2.tsv", groups, list(groups[0]))
    write_tsv(
        METADATA_DIR / "rna_seq_sample_context_v2.tsv",
        sample_contexts,
        list(sample_contexts[0]),
    )
    write_tsv(
        METADATA_DIR / "rna_seq_group_members_v2.tsv",
        memberships,
        list(memberships[0]),
    )
    write_tsv(
        METADATA_DIR / "rna_seq_sample_decisions_v2.tsv",
        decisions,
        list(decisions[0]),
    )
    write_tsv(
        METADATA_DIR / "rna_seq_exclusions_v2.tsv",
        exclusions,
        list(exclusions[0]),
    )
    write_tsv(
        METADATA_DIR / "ontology_mapping_v2.tsv",
        ontology_rows,
        list(ontology_rows[0]),
    )
    write_tsv(
        METADATA_DIR / "metadata_resolutions_v2.tsv",
        resolutions,
        list(resolutions[0]),
    )
    pair_fields = [
        "group_id",
        "sample_a",
        "sample_b",
        "pair_type",
        "log1p_1kb_pearson",
        "log1p_gene_body_pearson",
        "total_signal_ratio",
        "zero_fraction_a",
        "zero_fraction_b",
        "zero_fraction_delta",
        "review_flag",
        "review_action",
    ]
    write_tsv(METADATA_DIR / "replicate_pair_qc_v2.tsv", replicate_pairs, pair_fields)
    summary_fields = [
        "group_id",
        "n_members",
        "n_pairs",
        "min_log1p_1kb_pearson",
        "mean_log1p_1kb_pearson",
        "min_log1p_gene_body_pearson",
        "mean_log1p_gene_body_pearson",
        "flagged_pairs",
        "group_decision",
    ]
    write_tsv(
        METADATA_DIR / "replicate_group_qc_v2.tsv", replicate_summary, summary_fields
    )
    write_tsv(
        METADATA_DIR / "scientific_group_review_v2.tsv",
        scientific_review,
        list(scientific_review[0]),
    )

    summary = {
        "schema_version": 1,
        "created_at": utc_now(),
        "input_samples": len(samples),
        "rna_seq_samples": sum(row["assay"] == "RNA-Seq" for row in samples),
        "non_rna_excluded": sum(row["assay"] != "RNA-Seq" for row in samples),
        "duplicate_current_signals_excluded": len(DUPLICATE_CANONICAL),
        "candidate_signal_members": len(memberships),
        "candidate_groups": len(groups),
        "multi_sample_groups": sum(int(row["n_members"]) > 1 for row in groups),
        "multi_biological_unit_groups": sum(
            int(row["n_biological_units"]) > 1 for row in groups
        ),
        "technical_multi_run_groups": sum(
            int(row["n_runs"]) > int(row["n_biological_units"]) for row in groups
        ),
        "singleton_groups": sum(int(row["n_members"]) == 1 for row in groups),
        "formal_groups": sum(row["include_formal_v2"] == "True" for row in groups),
        "replicate_pair_qc_rows": len(replicate_pairs),
        "replicate_groups_flagged": sum(
            int(row["flagged_pairs"]) > 0 for row in replicate_summary
        ),
        "ontology_source_sha256": {
            "WBbt": sha256(WBBT_PATH),
            "WBls": sha256(WBLS_PATH),
        },
        "deterministic_manifest_hash_match": deterministic_match,
        "sample_manifest_sha256": sha256(METADATA_DIR / "rna_seq_samples_v2.tsv"),
    }
    (METADATA_DIR / "p2_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n"
    )
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
