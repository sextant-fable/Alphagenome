# External Benchmark Placeholder Contract

## Purpose

This document lets the manuscript be planned as a complete Nature Methods
Article without treating planned external evidence as if it had already been
observed. All fields below are write-once values taken from a frozen benchmark
manifest, review record, source-data package and figure artifact.

## Freeze Rule

No title, abstract, result sentence, figure label, table cell or cover-letter
claim may replace a placeholder until all of the following are true:

1. The external data source and isolation rule are committed to a dated
   manifest.
2. The model, comparator versions, input/output adapters, training budget,
   hyperparameter budget and metrics are frozen.
3. The completed run has a signed output manifest, source data and independent
   review/checksum report.
4. The planned endpoint and its primary statistical unit have passed a
   predeclared analysis check.
5. Any new biology claim has a verified external truth source and an exact
   result, not a locus selected after inspection.

## Required External Benchmark Fields

| Placeholder | Required value before drafting | Used in |
| --- | --- | --- |
| `EXTERNAL_BENCHMARK_DATASET_ID` | Public accession(s), version, license and raw/processed route. | Fig. 1, Table 1, Methods. |
| `EXTERNAL_SOURCE_ISOLATION_RULE` | Study/lab/submitter isolation and proof no source enters training, selection or tuning. | Fig. 1, Fig. 2, Methods. |
| `EXTERNAL_UNSEEN_CONDITION` | Important biological condition absent from development, with an objective definition. | Fig. 1, Fig. 2, Discussion. |
| `EXTERNAL_REFERENCE_AND_COORDINATE_RULE` | Genome build, annotation version, liftover if any and input coordinate convention. | Table 1, Methods. |
| `EXTERNAL_LEAKAGE_AUDIT_RESULT` | Window, receptive-field, run/track and study overlap result. | Fig. 1, ED Fig. 3, Methods. |
| `EXTERNAL_PROTOCOL_VERSION` | Frozen input window, target/bin, normalization, training and evaluation protocol SHA-256. | Table 1, Methods. |
| `EXTERNAL_PUBLIC_BASELINE_RESULTS` | Per-model metrics, total/trainable parameters, compute, adaptation notes and CIs. | Fig. 2, Table 1, Methods. |
| `EXTERNAL_MODEL_B_EFFECT_AND_CI` | Primary effect direction, estimate, 95% CI, unit and statistic. | Abstract, Fig. 2, Results. |
| `EXTERNAL_PER_STUDY_AND_CONDITION_RESULTS` | Source/condition effect tables and planned stratified intervals. | Fig. 2, Fig. 5, ED Fig. 4. |
| `EXTERNAL_VARIANT_DATASET_ID` | eQTL/allele/perturbation truth accession and independence rule. | Fig. 4, Methods. |
| `EXTERNAL_VARIANT_OR_PERTURBATION_RESULT` | Direction, ranking, matched-null and locus endpoint results with CI. | Abstract, Fig. 4, Discussion. |
| `EXTERNAL_RELEASE_MANIFEST` | Code, data, model, environment and clean reproduction identifiers. | Fig. 5, Availability. |

## Comparator Eligibility Matrix

The benchmark must explicitly classify every candidate model before results:

| Comparator | Eligible for shared main table only if | Otherwise |
| --- | --- | --- |
| Model A | Frozen trunk and same data/target/protocol adapter are fixed. | Treat as an internal transfer baseline. |
| P9 size-matched Model C | Trainable parameter match, data, target, optimizer budget and evaluation match. | Do not call an unmatched scratch model size matched. |
| Enformer | Input/output adapter preserves scientifically meaningful RNA readout and protocol differences are disclosed. | Put native-window sensitivity analysis in supplement or list as ineligible. |
| Borzoi | RNA coverage adapter, reference/target definitions and compute allowance are technically faithful. | Include only an explicitly non-comparable native protocol outside the main winner table. |
| Other public model | Same reference build, windows, target definition and no hidden information advantage. | State exact technical reason for non-inclusion. |

## External Biology Contract

The primary biological validation is a same-species cis-eQTL endpoint selected
before predictive scores are inspected. An exact regulatory-allele data set may
be an orthogonal supplementary endpoint, but cannot replace the primary test.
Predeclare:

- trait/readout and molecular condition;
- cis window and target-gene assignment rule;
- reference/alternative and haplotype policy;
- ambiguous site, indel and hyper-divergent-region handling;
- null matching variables;
- primary direction endpoint;
- primary ranking endpoint;
- blocked bootstrap/resampling unit;
- preselected illustrative loci and the rule that selected them.

## Existing Internal Evidence Boundary

- P6B/P9 can populate internal development and component-attribution displays.
- The P6C locked evaluation is a frozen within-collection performance record.
- P10 DPY-27 is not a positive application: its predicted chromosome-level
  contrast had the opposite sign to observed data. It must remain an Extended
  Data diagnostic and a Discussion boundary.
- The variant-scoring CLI has engineering and synthetic-test evidence only until
  a real external truth run fills the fields above.
