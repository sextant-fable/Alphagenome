# C. elegans RNA-seq v2 Execution Plan

This plan implements a fail-closed workflow for replacing the exploratory legacy RNA-seq11 pipeline with an auditable v2 dataset and model comparison.

## State Model

Each phase follows `PENDING -> RUNNING -> REVIEWING -> PASS`. A failed review moves the phase to `FAIL`; a protected resource boundary moves it to `APPROVAL_REQUIRED`. Only `PASS` advances the controller. `run --auto` continues through consecutive implemented phases and stops at the first failed review, missing scoped approval, unregistered phase, or completed workflow.

The canonical state is `alphagenome_custom/metadata/v2/execution_state.json`. Phase evidence and append-only reviews are stored under `alphagenome_custom/metadata/v2/audits/`.

## Phases and Reviews

1. `P0/R0`: freeze the legacy 20-sample, 11-track, 188-window data and model lineage without modifying large assets.
2. `P1/R1`: build a 485-accession provenance manifest with evidence and A/B/C/D normalization classification.
3. `P2/R2`: resolve duplicate content, the invalid T38 grouping, ontology, replicate context, inclusion, and exclusion decisions.
4. `P3A/R3A`: run the bounded five-accession source-read pilot and measure integrity, mapping, scale, I/O, runtime, disk use, and agreement with the provided signals.
5. `P3B/R3`: after a separate bulk approval informed by R3A, create immutable normalized and grouped bigWigs for the approved formal scope.
6. `P4/R4`: implement manifest-driven BigWig loading, five autosomal validation folds, and a locked chromosome-X test embargo.
7. `P5/R5`: implement model-space two-resolution count loss, gene loss, augmentation, and genuine C. elegans organism adaptation.
8. `P6A/R6A`: run the approved GPU smoke/pilot and validate environment, numerics, resource use, logging, and checkpoint reload.
9. `P6B/R6B`: run the pre-registered ablations and fair A/B/C blocked-CV comparison without reading chromosome X; select and lock one checkpoint.
10. `P6C/R6C`: after a separate final-test approval, read chromosome X once, permanently close the test entry point, and issue the final report.

## Mandatory Approval Boundaries

- `G1`: large FASTQ/BAM downloads or realignment.
- `G2`: scientific approval of the final sample/group manifests.
- `G3`: large normalized data or cache generation.
- `G4`: any GPU smoke test or formal experiment.
- `G5`: the single final chromosome-X evaluation.

Approvals are records with an exact scope, note, and timestamp rather than reusable booleans. The P3A pilot scopes do not authorize P3B bulk processing, and G4 does not authorize G5. Passing an automated review does not override these boundaries. After the exact approval is recorded, the bounded phase may run automatically through its next review.

## Completion Standard

Completion requires an auditable v2 provenance chain, a common and justified RNA-seq scale, deterministic on-demand loading, locked blocked cross-validation, corrected AlphaGenome-style losses, fair comparison of frozen-head, worm-embedding/LoRA, and from-scratch models, and one immutable final-test report.
