# Main-text Writing and Supplementary-evidence Reading Report

## Scope and Reading Record

This report separates two jobs that are often confused. The main text of the
comparison corpus was read to learn argument structure, figure order, paragraph
rhythm and claim density. The Supplementary Information (SI) was downloaded and
audited to identify the controls, stratifications and source-data practices that
make the main-text claims credible. SI prose is not a stylistic model for the
AlphaGenome manuscript.

The corpus comprises six open-access papers downloaded from their official
Nature Portfolio pages on 18 August 2026. The local package contains six article
PDFs, 12 PDF supplementary or reporting-summary files, two XLSX supplementary
workbooks and one CSV table. All files passed format checks and SHA-256 capture.
The generated local files live under `literature/full_text_20260818/`, which is
intentionally ignored by Git. The official article URLs and DOIs are recorded in
the table below.

| Short name | Journal and year | DOI | Why it was read |
| --- | --- | --- | --- |
| Borzoi | *Nature Genetics* (2025) | 10.1038/s41588-024-02053-6 | The closest RNA-seq-coverage method paper. |
| AlphaGenome | *Nature* (2026) | 10.1038/s41586-025-10014-0 | Current multi-modal sequence-to-function reference. |
| Nucleotide Transformer | *Nature Methods* (2025) | 10.1038/s41592-024-02523-z | Foundation-model benchmark and efficiency narrative. |
| P-SAGE-net | *Nature Methods* (2026) | 10.1038/s41592-026-03124-8 | Correct personal-genome split and expression generalization boundaries. |
| scooby | *Nature Methods* (2025) | 10.1038/s41592-025-02854-5 | Conditional sequence-to-RNA modelling and variant-effect use. |
| CREsted | *Nature Methods* (2026) | 10.1038/s41592-026-03057-2 | Broadly applicable software plus cross-system and experimental validation. |

## What the Main Texts Actually Do

### Borzoi: RNA coverage is the object, not an aggregate afterthought

The opening establishes why base-resolution RNA coverage is a harder and more
useful target than a single gene-expression scalar. Figure 1 then does three
things at once: it introduces the architecture, places observed and predicted
coverage at a real held-out locus in front of the reader, and shows the global
test-set statistics. This is the key visual pattern for the present project.
The next figures expand the same output into distinct biological questions:
tissue-specific transcription and isoforms (Fig. 2), sequence motifs and
attribution (Fig. 3), distal regulation and enhancer-gene links (Fig. 4),
expression-QTL effects (Fig. 5), alternative polyadenylation (Fig. 6) and
splicing/polyadenylation QTL classification (Fig. 7). Each figure has a single
new endpoint; none is simply another correlation plot.

The SI shows what supports this concise main text: per-gene and per-bin
stratification by transcript length, exon number, expression variance and
tissue; attention and attribution robustness; more coverage examples; benchmark
breakdowns; CRISPR enhancer perturbation; eQTL, MPRA and QTL sensitivity
analyses. The lesson is not to copy seven main figures. It is to make a compact
main-text ladder, then put the broad diagnostic surface behind it.

**Transfer to this project.** The equivalent first evidence figure must contain
an original model-and-data schematic, genuinely held-out genomic validation,
and observed-versus-predicted 128-bp/1-bp coverage at pre-specified development
loci. The existing aggregate P6B/P9 figures cannot substitute for the locus
panel. A later figure can introduce a real *C. elegans* variant application;
the current synthetic-only variant scorer cannot.

### AlphaGenome: broad method, deliberately separated evidence tiers

The main text begins with a unified sequence model and its training regimes
(Fig. 1), immediately makes the scale and competitive context visible, then
moves to representative track fidelity plus detailed performance (Fig. 2). It
does not conflate track prediction with variant interpretation. Splicing-QTL
prediction (Fig. 3), expression-QTL prediction (Fig. 4), chromatin-QTL
prediction (Fig. 5), cross-modal mechanism at loci (Fig. 6), and explicit
resolution/sequence-length/ensembling/distillation ablations (Fig. 7) are
separate claims.

Its SI is extensive because the breadth of the main argument demands it:
multi-assay prediction examples, per-track performance strata, scorer
definitions, QTL and MPRA analyses, failure examples, data-processing tables,
baseline settings and all source tables. The work does not hide fairness
decisions in a sentence; comparator input length, outputs and scoring are made
auditable.

**Transfer to this project.** The manuscript must give competitive performance,
component attribution, external generalization and variant use separate visual
roles. Model B should be compared to frozen A, size-matched C and genuinely
adapted public methods under an explicit common protocol. The current P9 figure
is already the right evidence type for the component claim, but it should not
carry external or biological claims.

### Nucleotide Transformer: benchmark discipline is a result

The main text introduces a pretraining and fine-tuning strategy with explicit
model scales (Fig. 1), benchmarks it across diverse downstream tasks (Fig. 2),
then shows representations and genomic-element knowledge (Fig. 3), variant
prioritization (Fig. 4), and compute-efficient architecture choices (Fig. 5).
The prose repeatedly names the comparison unit, not just the winning number.

The SI contains the full task matrix, layer selection, hyperparameters,
per-task results, size comparisons and supplementary interpretation. Its main
text is therefore readable because the detail is discoverable without being
repeated.

**Transfer to this project.** State the training budget, trainable parameter
count and exact matching conditions in the first comparison legend and Methods.
The original small C must never be described as parameter-matched; only the P9
`C_size_matched` control earns that description. Per-track results belong in
Source Data and SI-style extended material, not as a claim of 241 independent
biological replicates.

### P-SAGE-net: the evaluation split defines the claim

This paper makes the unseen unit explicit in Fig. 1: personal genomes are split
by individuals and chromosomes, and the model is evaluated on the claim it
actually makes. Figure 2 then tests performance across experimental conditions,
gene-set size and variants. The article is especially useful because it makes a
boundary scientifically productive rather than burying it: across-gene and
across-person prediction are not interchangeable.

The SI adds alternative variant sets, chromosome-level analyses, capacity and
input-length sensitivity, ablations, and attribution enrichment. These analyses
do not retroactively change the definition of a held-out personal genome.

**Transfer to this project.** The manuscript must say that P6B/P9 and the
DPY-27 analysis are genomic-block validation, while a study/laboratory-held-out
benchmark is a different claim. Because the current model consumes DNA plus a
fixed output-track head, it does not predict an entirely unseen perturbation
zero-shot. The external benchmark must be designed around a valid query:
unseen study provenance for a target representation already supported by the
model, or a new condition-aware/few-shot method with its own pre-registered
protocol.

### scooby: conditioning turns sequence prediction into context-specific use

The main text begins with a cell-state-conditioned modelling problem and shows
cell-state-specific profile fidelity in Fig. 1. It escalates from cell-type gene
counts and between-cell-type effects (Fig. 2), to motif perturbation (Fig. 3),
known variant effects (Fig. 4), and cell-type resolution of bulk eQTLs (Fig. 5).
The narrative has a clean direction: prediction first, mechanism second, real
variant use last.

The SI provides the architecture, memory representation, more cell types,
ablation results, comparison to another method and additional eQTL examples.
These checks reinforce rather than replace the main story.

**Transfer to this project.** This is the most consequential architectural
comparison. Model B's public-data head is fixed over 241 grouped tracks and its
forward interface takes DNA plus a fixed organism index; it has no treatment or
cell-state input. A Nature Methods-level claim about generalizing to unseen
conditions therefore requires either a conditional decoder/track embedding or
a carefully frozen adaptation task. Until that extension exists, the manuscript
should present its verified achievement as organism-adapted multi-track RNA
coverage prediction, not universal condition generalization.

### CREsted: software is credible when it crosses systems and reaches an assay

The main text starts with a usable package and a concrete biological dataset
(Fig. 1), moves across mouse cortex, human PBMC and cancer state examples
(Figs. 2-4), compares transfer learning and a large pretrained baseline (Fig.
5), and finishes with prospective synthetic enhancer design and in vivo testing
(Fig. 6). The work earns broad-applicability language through evidence across
systems, not through the package description.

The SI supplies ten supplementary figures and notes, additional benchmarks,
model settings and source tables. The main text retains only the evidence that
moves the argument forward.

**Transfer to this project.** A strong final version should finish with one
compact real biological application. The immediate computational route is a
held-out public *C. elegans* cis-eQTL/allelic-effect benchmark. The stronger
route is a pre-registered blind selection of 8-20 regulatory candidates for
reporter, CRISPR or independent RNA quantification. The current DPY-27 panel
is an internal negative diagnostic and must not be rewritten as that endpoint.

## Shared Main-text Grammar

The six articles share a stable evidence ladder:

1. Define a specific biological prediction problem and show the data/model
   interface.
2. Demonstrate direct held-out fidelity at both a locus and global scale.
3. Establish fair comparative performance and expose the exact unseen unit.
4. Explain which methodological design choices produce the gain.
5. Use the model for a genuinely independent biological or genetic prediction.
6. Move exhaustive stratification, source-data tables, robustness checks and
   reporting detail to Extended Data, SI and Methods.

Their prose is direct. A paragraph begins with the question answered by its
figure, gives the decisive design and observation, then explains the biological
meaning. It does not narrate every preprocessing operation or pre-empt every
objection. Confidence comes from a visible evidence chain, not from adjectives.

## Visual and Analysis Gap Audit for AlphaGenome v2

The project does **not** have too little quantitative data. It has 482
uniformly processed RNA-seq runs aggregated into 241 biological groups, 90
formal development cells for the P9 matrix, five blocked folds with three
nested seeds, and a locked six-chromosome evaluation. What is currently sparse
is the *diversity of evidence roles* made visible to a reader.

| Evidence layer seen in the corpus | Current verified AlphaGenome v2 status | Suitable manuscript location | Required action |
| --- | --- | --- | --- |
| Data resource, task definition and architecture | Data and model metadata are complete; no integrated figure yet. | Main Fig. 1 | Draw an original data-to-track-to-model schematic; show the 482-to-241 hierarchy and blocked genomic split. |
| Observed-versus-predicted held-out coverage at loci | Not yet rendered for v2. | Main Fig. 1 or Fig. 2 | New development-only P11 analysis: pre-register locus selection, run frozen validation checkpoints, export 1-bp and 128-bp coverage panels with source data. |
| Global matched model comparison | Complete: P6B and P9; B exceeds A and size-matched C under paired fold/seed design. | Main Fig. 2 | Redraw as a unified visual with the exact matched protocol and a compact track-strata panel. |
| Resolution, calibration and error anatomy | 1-bp head control is complete; calibration, gene-architecture strata and coverage-shape diagnostics are not yet a coherent package. | Main Fig. 3 and Extended Data | P11 development analysis from existing checkpoints and GTF: exon boundaries, expression/length/exon-count strata, observed-predicted quantiles, residual structure. |
| Component attribution | Complete: P9 demonstrates a dominant LoRA effect, smaller worm-embedding effect, non-additivity and no decisive primary-score gain for the 1-bp-only head. | Main Fig. 4 | Use the completed P9 figure after typography harmonization; move the full 241-track distribution/source tables to extended material. |
| Internal biology diagnostic | Complete but fails its pre-specified chromosome-level recovery endpoint. | Extended Data / limitations record | Retain transparently as an internal stress test; retitle neutrally. It cannot support a positive biology figure. |
| Independent study/laboratory benchmark | Absent. No locally held-out early/new run is independent because all 482 runs entered the current v2 dataset. | Main Fig. 5 | Freeze a study/lab-held-out protocol before tuning; obtain data not represented in the 482-run manifest. |
| Public-method benchmark | Absent. Enformer/Borzoi have not been harmonized to the project task. | Main Fig. 5 | Build adapters and a fixed 131,072-bp-to-128-bp score contract; report unsupported tasks separately rather than forcing a surface comparison. |
| Variant-effect biological application | Scorer code and synthetic tests exist; no real VCF/eQTL truth evaluation. | Main Fig. 6 | Pre-register a public cis-eQTL/haplotype benchmark, then score ref/alt or personalized haplotypes with frozen weights and compare direction/ranking. |
| Experimental confirmation | Absent. | Main Fig. 6 or separate follow-up | If possible, blind-select regulatory candidates from the frozen model for reporter/CRISPR/qPCR validation. |
| Reuse and reproducibility | Controller, manifests and source-data structures are strong; clean independent install/reproduction is still absent. | Methods, Code Availability and Extended Data | Make a clean environment run and release checksum-indexed inputs, scripts, model card and source data. |

## Recommended Figure Architecture

The final Article should use six non-redundant main figures. This is an
evidence plan, not a claim that all six are already supported.

### Figure 1: A worm RNA-coverage prediction resource with a transparent held-out task

- **a** Source-to-target hierarchy: 482 runs, metadata harmonization, 241
  grouped RNA-seq tracks, staged/tissue/treatment coverage.
- **b** Model B input, adaptation modules and dual-resolution heads. Include
  trainable-parameter count and fixed-track output interface.
- **c** Leave-one-chromosome-out development protocol over I-V and the distinct
  historical one-time six-chromosome locked evaluation. Do not imply that X is
  pristine.
- **d-e** Pre-specified observed and predicted 128-bp plus 1-bp coverage
  examples from frozen development validation checkpoints.

### Figure 2: Matched model comparison establishes the adapted model's gain

- **a** Global paired fold/seed scores for A, B, B without paper loss,
  C-size-matched and other registered comparisons.
- **b** Forest plot of B versus A and B versus C-size-matched, with five
  genomic folds as the inference unit and three nested seeds.
- **c** 241-track distribution, explicitly descriptive.
- **d** Pre-defined source/stage/tissue/treatment/signal strata.

The current P6B figure contains much of this content. P9 now supplies the
parameter-matched C result and should replace the earlier non-size-matched C
claim in the final integrated display.

### Figure 3: Coverage fidelity is resolved by genomic structure and scale

- **a** Gene-level and 128-bp metrics on the locked evaluation, labeled once as
  locked evaluation rather than repeatedly treated as new test data.
- **b** 1-bp exon-boundary profiles at pre-specified development loci.
- **c** Performance across expression strength, transcript length and exon
  count, shown as strata rather than inflated replicates.
- **d** Calibration/residual diagnostics and a failure example.

The first two panels require P11 outputs. The locked values may be reported as
the one-time endpoint already preserved in P6C; no additional read is allowed.

### Figure 4: Matched component experiments identify what makes adaptation work

- **a** B versus frozen A, size-matched C, no-LoRA, LoRA-only and 1-bp-only
  configurations.
- **b** Worm embedding, LoRA and their statistical non-additivity.
- **c** Track-wise effect heterogeneity.

This is the existing P9 figure, with captions revised for Article language.
The conclusion is specific: LoRA supplies the major reproducible gain; the
worm embedding adds a smaller reproducible increment; the 1-bp-only head does
not establish necessity.

### Figure 5: A frozen external benchmark establishes portability and fair competition

- **a** Study/laboratory-held-out design with no overlapping accession,
  study or biological source in the 482-run training resource.
- **b** Aggregated external performance versus frozen B, A, size-matched C and
  technically applicable public baselines, under a common score contract.
- **c** Per-track/study distribution and prespecified failure strata.
- **d** One external observed-predicted coverage locus.

This figure is currently missing. It is not replaced by the current 20 legacy
or 26 later-added runs because they were checked and are included in v2.

### Figure 6: Frozen sequence variants predict independent regulatory effects

- **a** Ref/alt or personalized-haplotype scoring protocol and external
  cis-eQTL truth set.
- **b** Direction concordance and ranking (AUROC/AUPRC or a prespecified
  matched-null comparison) with gene/block bootstrap intervals.
- **c** Representative allele with predicted coverage change and observed
  expression effect.
- **d** ISM/attribution at the same locus, followed by blind experimental
  validation when available.

This is currently missing. A real public cis-eQTL set is the minimum
computational endpoint; reporter/CRISPR/qPCR turns it into the strongest final
figure.

## What Can Be Added Now, Without Pretending It Is External Evidence

1. **P11 development-only fidelity and diagnostics package.** Freeze the
   checkpoints, validation manifests, locus-selection rule, GTF version,
   statistics and output hashes before inference. It should generate Fig. 1d-e
   and Fig. 3 development panels, source tables and audit files. It must not
   access `test_locked.tsv` or alter P6C.
2. **Main-figure redesign from finished P6B/P9 evidence.** The current figures
   are statistically sound but look like an internal report because they lead
   with a matrix. Redraw them in the figure roles above, with a common palette,
   consistent font, legends and visual hierarchy. Do not change values.
3. **Variant scorer production hardening.** Finish VCF/reference checks,
   coordinate audit, RC ensemble, ref/alt windows, gene aggregation, ISM and
   per-variant source data. This is software readiness, not an application
   result until external truth is scored.
4. **Release readiness.** Create a clean reproduction environment, model-card
   terminology, data/track manifests, source-data index and a precise code/data
   availability plan.

## Work That Cannot Be Replaced by Better Writing or More Internal Plots

1. An accession-, study- and laboratory-isolated dataset evaluated once after
   freezing the model and analysis plan.
2. A fair public-method comparison on a technically valid common task.
3. An independent genetic truth set for variant effects, ideally followed by a
   targeted experimental assay.
4. If broad unseen-condition generalization is central to the claim, a
   condition-aware architecture or a separately pre-registered adaptation
   protocol. The current fixed output head cannot justify that claim.

## Editorial Decision for the Rewrite

The next manuscript version should be rebuilt, not polished. It must use real
figures rather than boxes, write from the verified P6B/P9/P6C evidence, and
reserve a clearly structured Results slot for the future external benchmark and
variant study. That slot must remain outside any submission-facing PDF until
there are signed result artifacts. The stronger narrative is not weakened by
this discipline: it clarifies that the method establishes a substantial
organism-adapted RNA-coverage modelling result today and has a defined path to
the evidence needed for a full Nature Methods claim.
