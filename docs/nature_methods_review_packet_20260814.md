# Immutable Nature Methods review packet

## Packet scope

This packet freezes the common evidence supplied to three mutually blind
pre-submission reviewers. It describes the repository and verified results as
of 2026-08-14. It is not a manuscript and contains no reviewer concerns,
recommendation, rebuttal, or synthesis.

The proposed paper is evaluated as a prospective Nature Methods Article about
an AlphaGenome-style sequence-to-RNA-track adaptation for *C. elegans*. A full
manuscript, abstract, author-defined central claim, prior-art section, and final
figure order were not supplied. Reviewers must mark those elements as not
assessable rather than infer them.

## Common journal criteria

Nature Methods states that it publishes novel methods and significant
improvements to established life-science techniques, including computational,
statistical, and machine-learning methods for biological data. Its current
scope requires strong validation, an application to an important biological
question, and performance comparisons with available approaches. The Article
description additionally names performance, reproducibility, general
applicability, and potential for discovering new biology.

- https://www.nature.com/nmeth/submission-guidelines/about/aims
- https://www.nature.com/nmeth/content
- https://www.nature.com/articles/nmeth.3686
- https://www.nature.com/articles/nmeth.3359

The common review axes are originality, scientific importance,
interdisciplinary interest, technical soundness, and readability for
nonspecialists. Reviewers may assess who would use the work and identify
technical failings, but must not claim to make the editor's final decision.

## Shared project description

The active v2 workflow adapts an AlphaGenome-style sequence-to-track model to
WBcel235 *C. elegans* RNA-seq. The audited provenance inventory contains 485
accessions, of which 482 RNA-seq runs were uniformly reprocessed into 241
biological RNA-seq groups. Targets are loaded from grouped bigWigs on demand at
1-bp and 128-bp resolution. Development uses five leave-one-chromosome-out
folds over chromosomes I through V. A separately gated one-time evaluation
covers chromosomes I through V and X. Chromosome X was exposed during legacy
work and is not pristine for the project as a whole.

The main model configurations are frozen-trunk head A, worm-embedding and LoRA
model B with a dual-resolution RNA head, and from-scratch model C. The formal
P9 extension adds LoRA-only, no-LoRA, learned-1-bp-head-only, and a
size-matched from-scratch C configuration.

## Verified evidence base

### Final locked evaluation

The selected Model B paper-loss checkpoint was evaluated once under the P6C
gate. It covered 83 aligned cores, 10,878,976 eligible bases, all six named
chromosomes, and all 241 tracks. Coverage fraction was 1.0. Mean per-track
gene-exon log1p Pearson correlation was `0.675019`; mean per-track 128-bp log1p
Pearson correlation was `0.664130`.

Evidence pointers:

- `alphagenome_custom/metadata/v2/final_test_lock.json`
- `alphagenome_custom/metadata/v2/final_test_report.json`
- lock SHA-256 `602988f99a952e23e98a205ef3e264c11b94cc22e06ad2213a109339cd9cb421`
- report SHA-256 `9c9b06c3eba15d95c30ff4e7403d97dbfe1b72cbc0e8dcd678fe98de762d3a6e`

### Matched development comparisons

The earlier 90-run P6B matrix contains A, B, and the original C under two
losses, five folds, and three seeds. A post hoc submission analysis uses folds
as the primary blocking unit and seeds within folds. For B with paper loss, the
composite score was `0.626602` (95% CI `0.611447` to `0.640932`). Its paired
gain was `0.111180` over A with paper loss, `0.139057` over the original C with
paper loss, and `0.011228` over B with log1p-MSE loss. The original C was not
parameter matched.

Evidence pointers:

- `alphagenome_custom/metadata/v2/p6b_submission_analysis_execution.json`
- `results/v2_p6b_submission_analysis/`
- `docs/v2_p6b_submission_analysis.md`

### P9 component and size-matched evidence

P9 registered 90 unique configuration-by-fold-by-seed records, including 59
new training jobs and 31 hash-verified reused records. R9 passed all nine
checks. All effects below use fold-first, seed-within-fold 10,000-replicate
bootstrap intervals and are candidate minus full Model B unless named as a
factorial main effect.

- Size-matched C primary effect `-0.092017` (95% CI `-0.099761` to
  `-0.084976`). C had 1,265,122 trainable parameters and B had 1,276,228.
- Frozen A primary effect `-0.111180` (95% CI `-0.116151` to `-0.107236`).
- No-LoRA primary effect `-0.094392` (95% CI `-0.098300` to `-0.091226`).
- LoRA-only primary effect `-0.002812` (95% CI `-0.004231` to `-0.001550`).
- Learned-1-bp-head-only primary effect `-0.001808` (95% CI `-0.004205` to
  `0.000707`).
- Worm-embedding factorial main effect `0.009800` (95% CI `0.008933` to
  `0.010841`).
- LoRA factorial main effect `0.101380` (95% CI `0.097917` to `0.105568`).
- Worm-by-LoRA statistical interaction `-0.013976` (95% CI `-0.015870` to
  `-0.012068`).

The inference unit is five genomic folds. The three seeds are algorithmic
repeats nested within folds. The 241 track outcomes are descriptive, not 241
independent replicates. The interaction is a statistical non-additivity term,
not evidence of a biological mechanism.

Evidence pointers:

- `alphagenome_custom/metadata/v2/audits/P9/review.json`
- `alphagenome_custom/metadata/v2/p9_submission_evidence_execution.json`
- `alphagenome_custom/metadata/v2/p9_submission_evidence_paired_effects.tsv`
- `alphagenome_custom/metadata/v2/p9_submission_evidence_factorial_effects.tsv`
- `alphagenome_custom/metadata/v2/p9_submission_evidence_per_track.tsv`

### P10 internal DPY-27 application

P10 evaluated 15 existing B paper-loss checkpoints on their matching
development validation cores. G0005 is DPY-27 RNAi and G0006 is vector RNAi.
Both conditions were training targets and each track had been independently
normalized to total signal of approximately 1e8. R10 passed all nine checks.

The observed chromosome-X-minus-autosome median contrast was `0.003065` (95%
CI `0.001980` to `0.004237`). The predicted contrast was `-0.002740` (95% CI
`-0.005352` to `-0.000369`). X-linked gene direction concordance was `0.6674`
(95% CI `0.6120` to `0.7176`), and predicted-observed gene-level Spearman
correlation was `0.0846` (95% CI `0.0707` to `0.0945`).

Evidence pointers:

- `alphagenome_custom/metadata/v2/audits/P10/review.json`
- `alphagenome_custom/metadata/v2/p10_dpy27_internal_application_execution.json`
- `runs/v2_p10_dpy27_internal_application_20260813/summary_endpoints.tsv`
- `runs/v2_p10_dpy27_internal_application_20260813/figure_qa.json`
- `docs/dpy27_internal_application_methods.md`

### Variant-scoring software

The repository contains a frozen-spec variant-scoring CLI with VCF and
reference checks, fixed-length ref/alt windows, forward/reverse-complement
ensembling, substitution-level gene and exon deltas, track-level indel
diagnostics, optional ISM, and checksum-based audit output. CPU synthetic tests
are available. No external eQTL dataset has been downloaded and no real VCF
plus checkpoint scoring run has been registered.

Evidence pointers:

- `scripts/v2_variant_scoring.py`
- `scripts/score_v2_variants.py`
- `tests/test_v2_variant_scoring.py`
- `docs/v2_variant_scoring_methods.md`

## Reproducibility and integrity controls

P9 and P10 ran through the controlled phase controller on host `hy8`, using
physical GPUs 2 and 3. P9 recorded commit `8063a06734833814f7ab17f8fee8bfa7286f2568`.
The successful P10 retry recorded commit
`9db53d20ad59929c8ccf8a2d34228fcb931345eb`. Both prohibit final-test access
and record zero locked-test signal reads. P10 preserves the failed first
attempt caused only by non-serializable TIFF DPI metadata, plus a clean rerun
under a corrected frozen plot hash.

## Materials not present in this packet

- No full manuscript, abstract, exact title, author-defined claim hierarchy, or
  final figure sequence.
- No benchmark held outside the 482-run provenance collection and isolated by
  study or laboratory.
- No validation on an unseen biological condition or independent species.
- No completed Enformer, Borzoi, or other technically applicable public-model
  comparison under a harmonized benchmark.
- No external cis-eQTL, regulatory-allele, reporter, CRISPR, or independent RNA
  validation result.
- No real variant-scoring execution with a frozen public VCF or haplotype set.
- No evidence packet establishing software installation on common operating
  systems, container distribution, public archival release, or independent
  user reproduction.
- No complete prior-art analysis sufficient to assess the novelty of each
  model component against current published methods.

These items are an inventory of supplied material. Reviewers must independently
decide their relevance and severity without receiving a shared concern list.
