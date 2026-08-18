# Proposed P11: Development-Only Coverage Fidelity and Diagnostic Package

## Purpose

P11 fills the most important missing internal evidence layer: a reader must be
able to inspect real observed and predicted coverage, not only an aggregate
score. It will create the development-only locus and error-structure panels
needed for Figure 1 and Figure 3 of the final Article plan.

P11 is a proposed controlled analysis phase. It is not yet a registered
controller phase, a GPU authorization or an executed experiment.

## Scientific Question

How well does frozen full Model B reproduce observed RNA-seq coverage across
pre-specified exon structures, gene classes and signal regimes in held-out
development genomic blocks?

The question is deliberately narrower than external generalization. It uses the
five development folds only and cannot read or alter the locked final test.

## Frozen Inputs Before Execution

1. The 15 full Model B/paper checkpoints, one per development fold and seed.
2. The corresponding five `valid.tsv` manifests only; paths containing
   `test_locked` or a final-test role fail closed before any model or bigWig
   access.
3. WBcel235 FASTA, FAI, GTF, grouped-track manifest, means and the 241 grouped
   bigWigs, each checksum-registered.
4. A versioned implementation list for the runner, coverage aggregation,
   plotting and review code.
5. The P6C lock/report hashes, which must be unchanged after P11.

## Pre-Specified Displays

### Development Coverage Examples

Display three loci selected before prediction inspection:

- one multi-exon young-adult whole-organism target;
- one embryonic target with a sharp exon boundary;
- one low-to-middle signal target with a documented failure pattern.

For each fold, candidate loci are eligible only when their gene is fully
contained in a validation core, has at least five annotated exons and has
nonzero observed signal in the nominated track. Within each display class,
select the lexicographically first `(fold, track_id, gene_id)` after sorting by
the class's frozen coverage quantile rule. The runner must write the candidate
population and selection table before it writes predictions. No selection may
depend on full-model accuracy, attribution, or the comparison model.

Each locus display contains observed and predicted coverage with identical
axes, exon annotation and a 128-bp summary. It reports the checkpoint fold,
seed aggregation rule, track ID, gene ID, coordinate convention and exact
source-data row range.

### Structural Accuracy and Error Anatomy

For every completed fold/checkpoint pair, produce:

- gene-exon and 128-bp performance by predeclared observed-signal tertile;
- gene length and exon-count strata defined from the frozen GTF;
- observed-predicted calibration curves and residual quantiles;
- exon-boundary profiles at 1-bp resolution;
- a fixed failure taxonomy: low signal, short gene, high exon count and
  coverage-edge context.

Folds are the inferential units and seeds are nested algorithmic repeats. Track
and gene distributions remain descriptive unless the registered resampling plan
explicitly blocks by fold and genomic locus.

## Analysis Rules

- Do not use the six-chromosome locked final evaluation, chromosome X as new
  test evidence, or `test_locked.tsv`.
- Do not tune model weights, target transforms, selection thresholds or figure
  panel boundaries after inspecting P11 predictions.
- Never average incompatible normalized tracks to claim an absolute expression
  scale.
- Keep observed and predicted arrays, masks, aggregation tables, selector
  output, figure source data, script hashes and exact invocation in the output
  manifest.
- The P11 title and legends describe genomic-block validation. They do not use
  external-validation, unseen-condition or biological-validation language.

## Controller and Review Contract

Before any GPU inference, add an explicit P11 controller transition and require
a distinct `G4:p11_development_fidelity_diagnostics` approval. The controller
must fail closed on final-test paths, verify the current hashes of inputs and
implementation files, record physical GPU 2/3 selection, and write one atomic
shard per fold/seed. A reviewer must independently recompute source-table row
counts, selected-locus rules, fold endpoints, figure hashes and the unchanged
P6C lock/report hashes.

## Deliverables

- Figure 1d-e coverage examples plus Source Data.
- Figure 3 development-scale and error-structure panels plus Source Data.
- Extended Data coverage examples and complete strata tables.
- P11 specification, execution record, review report and reproducibility
  manifest.

## What P11 Does Not Supply

P11 improves the manuscript's empirical clarity but cannot replace a
study/laboratory-isolated benchmark, public baseline comparison, independent
variant truth set, or experimental regulatory validation.
