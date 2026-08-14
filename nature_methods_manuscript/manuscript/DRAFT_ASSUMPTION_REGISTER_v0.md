# Initial Draft Assumption Register

Status: owner-authorized internal scaffold, 2026-08-14.

This register allows the manuscript story to be drafted before the external
benchmark and independent cis-eQTL run are delivered. It is not evidence and
must not be uploaded with a submission. The manuscript keeps the expected
result direction in prose. The source file maps the registered fields through
`Ext*` macros, while the readable PDF uses fluent assumed-result wording so the
story is not interrupted by bracketed text. Replace those macro values
atomically from the signed benchmark manifest, source-data tables and
independent review record.

## Replacement tokens

| Token | Required replacement |
| --- | --- |
| `{{EXT_BENCHMARK_DATASET_ID}}` | Public accession/version, data route and license. |
| `{{EXT_SOURCE_ISOLATION_RULE}}` | Study/laboratory isolation rule and audit result. |
| `{{EXT_UNSEEN_CONDITION}}` | Important biological condition absent from development. |
| `{{EXT_REFERENCE_COORDINATE_RULE}}` | Reference build, annotation and coordinate rule. |
| `{{EXT_LEAKAGE_AUDIT_RESULT}}` | Window, source and study overlap audit. |
| `{{EXT_PROTOCOL_VERSION}}` | Frozen protocol identifier and checksum. |
| `{{EXT_PRIMARY_EFFECT_AND_CI}}` | Primary Model B effect, unit and 95% CI. |
| `{{EXT_PUBLIC_BASELINE_RESULTS}}` | Per-model metrics, parameters, compute and eligibility. |
| `{{EXT_STRATIFIED_RESULTS}}` | Study/condition/annotation strata and intervals. |
| `{{EXT_VARIANT_DATASET_ID}}` | Independent cis-eQTL truth source and version. |
| `{{EXT_VARIANT_RESULT}}` | Direction, ranking and matched-null results with CIs. |
| `{{EXT_RELEASE_MANIFEST}}` | Code, model, data, environment and reproduction IDs. |

## Drafting rule

The assumed direction is “full Model B improves the primary external endpoint,
and frozen sequence-difference scores recover the prespecified cis-eQTL
signal”. No numeric external result, accession, author, DOI, p-value, sample
size or repository is invented in this draft. The exact result sentence is
replaced only after the external artifact passes the five freeze checks in
`plans/08_external_benchmark_placeholder_contract.md`.

## Verified internal values already used in the draft

- The source inventory contains 482 verified RNA-seq runs grouped into 241
  biological tracks on WBcel235.
- P9 reports the registered paired comparison of size-matched Model C minus
  full Model B as -0.0920167 (95% CI -0.0997615 to -0.0849763), a LoRA main
  effect of 0.1013801 (0.0979171 to 0.1055677), a worm-embedding main effect
  of 0.0098003 (0.0089332 to 0.0108411), and a non-additive interaction of
  -0.0139765 (-0.0158702 to -0.0120685).
- The P10 DPY-27 internal diagnostic did not recover the prespecified
  chromosome-level response: predicted X-minus-autosome contrast -0.002740
  (-0.005352 to -0.000369) versus observed 0.003065 (0.001980 to 0.004237),
  with direction concordance 0.6674 and gene-level Spearman correlation
  0.0846. It is kept in Extended Data and the Discussion boundary paragraph.

## Release gate

The draft can be read for story and prose now. It cannot be labeled “final” or
used to populate a submission system until all `{{EXT_*}}` tokens are removed,
the corresponding source-data rows exist, and the manuscript-wide metrics pass
the statistics and citation audits.
