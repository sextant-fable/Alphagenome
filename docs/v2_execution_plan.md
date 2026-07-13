# C. elegans RNA-seq v2 Execution Plan

This plan implements a fail-closed workflow for replacing the exploratory legacy RNA-seq11 pipeline with an auditable v2 dataset and model comparison.

## State Model

Each phase follows `PENDING -> RUNNING -> REVIEWING -> PASS`. A failed review moves the phase to `FAIL`; a protected resource boundary moves it to `APPROVAL_REQUIRED`. Only `PASS` advances the controller.

The canonical state is `alphagenome_custom/metadata/v2/execution_state.json`. Phase evidence and append-only reviews are stored under `alphagenome_custom/metadata/v2/audits/`.

## Phases and Reviews

1. `P0/R0`: freeze the legacy 20-sample, 11-track, 188-window data and model lineage without modifying large assets.
2. `P1/R1`: build a 485-accession provenance manifest with evidence and A/B/C/D normalization classification.
3. `P2/R2`: resolve duplicate content, the invalid T38 grouping, ontology, replicate context, inclusion, and exclusion decisions.
4. `P3/R3`: create immutable normalized and grouped bigWigs after an explicit data-write approval.
5. `P4/R4`: implement manifest-driven BigWig loading, five autosomal validation folds, and a locked chromosome-X test embargo.
6. `P5/R5`: implement model-space two-resolution count loss, gene loss, augmentation, and genuine C. elegans organism adaptation.
7. `P6/R6`: run approved smoke tests, five-fold/three-seed model comparisons, lock one model, and evaluate chromosome X once after final-test approval.

## Mandatory Approval Boundaries

- `G1`: large FASTQ/BAM downloads or realignment.
- `G2`: scientific approval of the final sample/group manifests.
- `G3`: large normalized data or cache generation.
- `G4`: any GPU smoke test or formal experiment.
- `G5`: the single final chromosome-X evaluation.

Passing an automated review does not override these boundaries. After an approval is recorded, the bounded phase may run automatically through its next review.

## Completion Standard

Completion requires an auditable v2 provenance chain, a common and justified RNA-seq scale, deterministic on-demand loading, locked blocked cross-validation, corrected AlphaGenome-style losses, fair comparison of frozen-head, worm-embedding/LoRA, and from-scratch models, and one immutable final-test report.
