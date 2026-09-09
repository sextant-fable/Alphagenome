#!/usr/bin/env python3
"""Associate public wild-strain expression with scored candidate alleles."""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import math
from datetime import datetime, timezone
from pathlib import Path
from statistics import mean


REPO_ROOT = Path(__file__).resolve().parents[1]
SOURCE_ROOT = REPO_ROOT / "shared/external_eqtl/p20_sources_v1"
COUNTS = SOURCE_ROOT / "GSE186719_Celegans_208strains_609samples_rawCounts.tsv.gz"
VCF = REPO_ROOT / "results/v2_p20_eqtl_variant_set/published_eqtl_candidates_iv.vcf"
CANDIDATES = REPO_ROOT / "results/v2_p20_eqtl_variant_set/published_eqtl_candidate_annotations.tsv"
GTF = REPO_ROOT / "alphagenome_custom/reference/Caenorhabditis_elegans.WBcel235.115.gtf"
P21_EXECUTION = REPO_ROOT / "alphagenome_custom/metadata/v2/p21_eqtl_variant_scoring_execution.json"
OUTPUT = REPO_ROOT / "results/v2_p22_eqtl_effect_association"
TARGET_TRACK = "RNA_V2_G0054"
SEEDS = (20260714, 20260715, 20260716)


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def strain_name(sample: str) -> str:
    return sample.split("_GZ", 1)[0]


def parse_float(value: str) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return 0.0
    return number if math.isfinite(number) and number >= 0 else 0.0


def gtf_transcript_to_gene(target_gene_ids: set[str]) -> dict[str, set[str]]:
    mapping: dict[str, set[str]] = {}
    for line in GTF.read_text(encoding="utf-8").splitlines():
        if not line or line.startswith("#") or "\ttranscript\t" not in line:
            continue
        attributes = line.split("\t", 8)[8]
        parsed = {}
        for field in attributes.split(";"):
            field = field.strip()
            if not field or " " not in field:
                continue
            key, value = field.split(" ", 1)
            parsed[key] = value.strip().strip('"')
        gene_id = parsed.get("gene_id", "")
        transcript_id = parsed.get("transcript_id", "")
        if gene_id in target_gene_ids and transcript_id:
            mapping.setdefault(transcript_id, set()).add(gene_id)
    return mapping


def load_expression(target_gene_ids: set[str]) -> tuple[dict[str, dict[str, float]], dict[str, list[int]], int]:
    """Return CPM-log expression by gene name, plus sample indices."""

    transcript_to_gene = gtf_transcript_to_gene(target_gene_ids)

    with gzip.open(COUNTS, "rt", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        # The header names the transcript column, while data rows contain an
        # additional integer row id before the transcript identifier.
        sample_names = header[1:]
        sample_strains = [strain_name(sample) for sample in sample_names]
        strain_to_indices: dict[str, list[int]] = {}
        for index, strain in enumerate(sample_strains):
            strain_to_indices.setdefault(strain, []).append(index)
        totals = [0.0] * len(sample_names)
        raw_targets: dict[str, list[float]] = {}
        has_row_id: bool | None = None
        for row in reader:
            if not row:
                continue
            if has_row_id is None:
                has_row_id = len(row) == len(header) + 1
            transcript = row[1] if has_row_id and len(row) > 1 else row[0]
            values = [parse_float(value) for value in row[2:] if has_row_id] if has_row_id else [parse_float(value) for value in row[1:]]
            if len(values) != len(sample_names):
                continue
            for index, value in enumerate(values):
                totals[index] += value
            for gene_name in transcript_to_gene.get(transcript, set()):
                if gene_name not in raw_targets:
                    raw_targets[gene_name] = [0.0] * len(values)
                raw_targets[gene_name] = [a + b for a, b in zip(raw_targets[gene_name], values)]
    expression: dict[str, dict[str, float]] = {}
    for gene_name, values in raw_targets.items():
        by_strain: dict[str, list[float]] = {}
        for index, strain in enumerate(sample_strains):
            cpm = values[index] * 1_000_000.0 / totals[index] if totals[index] > 0 else 0.0
            by_strain.setdefault(strain, []).append(math.log1p(cpm))
        expression[gene_name] = {strain: mean(values) for strain, values in by_strain.items()}
    return expression, strain_to_indices, len(sample_names)


def load_vcf() -> tuple[dict[tuple[str, int], dict[str, object]], list[str]]:
    records: dict[tuple[str, int], dict[str, object]] = {}
    sample_names: list[str] = []
    with VCF.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#CHROM"):
                fields = line.rstrip("\n").split("\t")
                sample_names = fields[9:]
                continue
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 10:
                continue
            chrom, position = fields[0], int(fields[1])
            fmt = fields[8].split(":")
            gt_index = fmt.index("GT") if "GT" in fmt else 0
            dosage: dict[str, float] = {}
            for sample, value in zip(sample_names, fields[9:]):
                parts = value.split(":")
                gt = parts[gt_index] if gt_index < len(parts) else "."
                alleles = gt.replace("|", "/").split("/")
                if len(alleles) != 2 or any(allele not in {"0", "1", "2"} for allele in alleles):
                    continue
                dosage[sample] = float(int(alleles[0]) + int(alleles[1]))
            records[(chrom, position)] = {
                "chromosome": chrom,
                "position_1based": position,
                "reference": fields[3],
                "alternate": fields[4].split(",")[0],
                "variant_key": f"{chrom}:{position}:{fields[3]}>{fields[4].split(',')[0]}:alt1",
                "dosage": dosage,
            }
    return records, sample_names


def pearson_slope(x: list[float], y: list[float]) -> tuple[float, float]:
    if len(x) < 3:
        return float("nan"), float("nan")
    x_mean, y_mean = mean(x), mean(y)
    xx = sum((value - x_mean) ** 2 for value in x)
    yy = sum((value - y_mean) ** 2 for value in y)
    if xx <= 0 or yy <= 0:
        return float("nan"), float("nan")
    covariance = sum((a - x_mean) * (b - y_mean) for a, b in zip(x, y))
    return covariance / math.sqrt(xx * yy), covariance / xx


def load_model_effects() -> dict[tuple[int, int], dict[str, dict[str, float]]]:
    effects: dict[tuple[int, int], dict[str, dict[str, float]]] = {}
    for fold in range(1, 6):
        for seed in SEEDS:
            path = REPO_ROOT / "results/v2_p21_eqtl_variant_scoring" / "evaluations" / f"seed_{seed}" / f"fold_{fold}" / "variant_gene_exon_deltas.tsv"
            if not path.is_file():
                continue
            rows: dict[str, dict[str, float]] = {}
            with path.open(newline="", encoding="utf-8") as handle:
                for row in csv.DictReader(handle, delimiter="\t"):
                    if row.get("track_id") != TARGET_TRACK or row.get("resolution") != "1":
                        continue
                    key = f"{row.get('chromosome')}:{row.get('position_1based')}"
                    gene = row.get("gene_name", "")
                    if not gene:
                        continue
                    rows[f"{key}:{gene}"] = {
                        "model_gene_body_delta_mean": float(row.get("gene_body_delta_mean", "nan")),
                        "model_exon_delta_mean": float(row.get("exon_delta_mean", "nan")),
                        "model_gene_id": row.get("gene_id", ""),
                    }
            effects[(fold, seed)] = rows
    return effects


def main() -> None:
    if OUTPUT.exists():
        raise FileExistsError(f"P22 output is immutable: {OUTPUT}")
    candidates = read_tsv(CANDIDATES)
    target_gene_ids = {row["representative_gene_id"] for row in candidates}
    expression, _, sample_count = load_expression(target_gene_ids)
    vcf_records, vcf_samples = load_vcf()
    model_effects = load_model_effects()
    if len(model_effects) != 15:
        raise RuntimeError(f"Expected 15 completed P21 model effect files, found {len(model_effects)}")

    by_position = {(row["chromosome"], int(row["position_1based"])): row for row in candidates}
    model_rows: list[dict[str, object]] = []
    summary_rows: list[dict[str, object]] = []
    for key, candidate in sorted(by_position.items()):
        vcf = vcf_records.get(key)
        if not vcf:
            continue
        transcript = candidate["representative_transcript"]
        observed = expression.get(candidate["representative_gene_id"], {})
        dosage = vcf["dosage"]
        pairs = [(float(dosage[sample]), observed[strain_name(sample)]) for sample in vcf_samples if sample in dosage and strain_name(sample) in observed]
        x, y = [pair[0] for pair in pairs], [pair[1] for pair in pairs]
        observed_r, observed_slope = pearson_slope(x, y)
        model_values: list[float] = []
        model_exon_values: list[float] = []
        for (fold, seed), rows in sorted(model_effects.items()):
            model = rows.get(f"{key[0]}:{key[1]}:{candidate['finemapping_gene']}")
            if not model:
                continue
            body = float(model["model_gene_body_delta_mean"])
            exon = float(model["model_exon_delta_mean"])
            model_values.append(body)
            model_exon_values.append(exon)
            model_rows.append({
                "chromosome": key[0],
                "position_1based": key[1],
                "variant_key": vcf["variant_key"],
                "representative_transcript": transcript,
                "representative_gene_name": candidate["representative_gene_name"],
                "finemapping_gene": candidate["finemapping_gene"],
                "fold": fold,
                "seed": seed,
                "track_id": TARGET_TRACK,
                "resolution": 1,
                **model,
            })
        signs = [1 if value > 0 else -1 if value < 0 else 0 for value in model_values]
        sign_concordance = sum(sign * (1 if observed_slope > 0 else -1 if observed_slope < 0 else 0) > 0 for sign in signs) / len(signs) if signs and math.isfinite(observed_slope) else float("nan")
        summary_rows.append({
            "chromosome": key[0],
            "position_1based": key[1],
            "variant_key": vcf["variant_key"],
            "reference": vcf["reference"],
            "alternate": vcf["alternate"],
            "representative_transcript": transcript,
            "representative_gene_name": candidate["representative_gene_name"],
            "finemapping_gene": candidate["finemapping_gene"],
            "eQTL_classification": candidate["eQTL_classification"],
            "published_eQTL_logP": candidate["published_eQTL_logP"],
            "strain_pairs": len(pairs),
            "observed_pearson_dosage_expression": observed_r,
            "observed_signed_slope": observed_slope,
            "model_jobs_mapped": len(model_values),
            "model_gene_body_delta_mean": mean(model_values) if model_values else float("nan"),
            "model_exon_delta_mean": mean(model_exon_values) if model_exon_values else float("nan"),
            "model_observed_sign_concordance": sign_concordance,
        })
    OUTPUT.mkdir(parents=True)
    summary_fields = list(summary_rows[0]) if summary_rows else ["variant_key"]
    model_fields = list(model_rows[0]) if model_rows else ["variant_key"]
    with (OUTPUT / "eqtl_variant_effects.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=summary_fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(summary_rows)
    with (OUTPUT / "eqtl_variant_effects_by_model.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=model_fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(model_rows)
    analyzable = [row for row in summary_rows if int(row["strain_pairs"]) >= 20 and int(row["model_jobs_mapped"]) == 15 and math.isfinite(float(row["observed_signed_slope"]))]
    manifest = {
        "schema_version": 1,
        "phase": "P22",
        "contract": "p22_public_eqtl_effect_association_v1",
        "created_at": utc_now(),
        "status": "PASS" if len(analyzable) >= 20 else "BLOCKED_INSUFFICIENT_MAPPED_EFFECTS",
        "scope": "public wild-strain dosage-expression association and frozen represented-head variant effects; not causal validation",
        "target_track": TARGET_TRACK,
        "counts_sha256": sha256(COUNTS),
        "vcf_sha256": sha256(VCF),
        "candidate_annotations_sha256": sha256(CANDIDATES),
        "p21_execution_sha256": sha256(P21_EXECUTION),
        "sample_count": sample_count,
        "vcf_sample_count": len(vcf_samples),
        "candidate_count": len(candidates),
        "summary_count": len(summary_rows),
        "analyzable_count": len(analyzable),
        "mapping_rule": "Supplementary Data 4 finemapping_gene matched to WBcel235 GTF gene_name; representative transcript supplies observed expression",
        "observed_expression_transform": "sample-level log1p(CPM), then mean within strain",
        "prohibited_claims": ["causal regulatory mechanism", "unseen-condition prediction", "population-level laboratory generalization"],
        "outputs": {
            "summary": "results/v2_p22_eqtl_effect_association/eqtl_variant_effects.tsv",
            "by_model": "results/v2_p22_eqtl_effect_association/eqtl_variant_effects_by_model.tsv",
        },
    }
    (OUTPUT / "p22_eqtl_effect_association_manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(manifest, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
