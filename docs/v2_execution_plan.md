# C. elegans RNA-seq v2 Execution Plan

This plan implements a fail-closed workflow for replacing the exploratory legacy RNA-seq11 pipeline with an auditable v2 dataset and model comparison.

## State Model

Each phase follows `PENDING -> RUNNING -> REVIEWING -> PASS`. A failed review moves the phase to `FAIL`; a protected resource boundary moves it to `APPROVAL_REQUIRED`. Only `PASS` advances the controller. `run --auto` continues through consecutive implemented phases and stops at the first failed review, missing scoped approval, unregistered phase, or completed workflow.

The canonical state is `alphagenome_custom/metadata/v2/execution_state.json`. Phase evidence and append-only reviews are stored under `alphagenome_custom/metadata/v2/audits/`.

## Phases and Reviews

1. `P0/R0`: freeze the legacy 20-sample, 11-track, 188-window data and model lineage without modifying large assets.
2. `P1/R1`: build a 485-accession provenance manifest with evidence and A/B/C/D normalization classification.
3. `P2/R2`: resolve duplicate content, the invalid T38 grouping, ontology, replicate context, inclusion, and exclusion decisions.
4. `P3A/R3A`: run five representative accessions as an automatic checkpoint and measure integrity, mapping, scale, I/O, runtime, disk use, and agreement with the provided signals.
5. `P3B/R3`: on R3A PASS, stream all 482 verified RNA-seq runs under the same full-processing approval, retaining normalized/grouped bigWigs while cleaning per-run FASTQ/BAM/bedGraph intermediates. Technical runs are pooled by raw coverage mass within biological units; true biological units are then averaged equally. The three reused provided signals are regenerated and retained as a separate source-study group rather than merged across studies.
6. `P4/R4`: implement manifest-driven BigWig loading and five within-chromosome blocked folds. Every train and validation fold contains disjoint blocks from `I`, `II`, `III`, `IV`, `V`, and `X`; a separate locked block from each chromosome forms the one-time final test. Valid/test effective blocks are tiled completely by context-relative, non-overlapping 131,072 bp metric cores; partial coverage is a review failure.
7. `P5/R5`: implement model-space two-resolution count loss, gene loss, augmentation, and genuine C. elegans organism adaptation.
8. `P6A/R6A`: run the approved GPU smoke/pilot and validate environment, numerics, resource use, logging, and checkpoint reload.
9. `P6B/R6B`: run the pre-registered ablations and fair A/B/C blocked-CV comparison without reading any locked test block; select and lock one checkpoint.
10. `P6C/R6C`: after a separate final-test approval, read the six locked chromosome blocks once, permanently close the test entry point, and issue the final report.

## Mandatory Approval Boundaries

- `G1`: large FASTQ/BAM downloads or realignment.
- `G2`: scientific approval of the final sample/group manifests.
- `G3`: large normalized data or cache generation.
- `G4`: any GPU smoke test or formal experiment.
- `G5`: the single final six-chromosome locked-block evaluation.

Approvals are records with an exact scope, note, and timestamp rather than reusable booleans. The user-approved `p3_full_rna_streaming_482` scope covers both P3A and P3B. The revised P4 scope is `p4_six_chromosome_block_split` and still prohibits a monolithic NPZ. The standing G4 scope permits later phases to select available GPU 2/3 but does not authorize G5. The old `r6c_single_chr_x_test` scope cannot authorize the revised workflow; only `r6c_single_six_chromosome_block_test` can unlock the final blocks after R6B.

The controller records `p4_six_chromosome_block_split` while the superseded P6C state is still active, then runs `prepare-six-chromosome-revision` to archive the holdout evidence and reopen P4. The same scoped record remains required for P4 itself; it grants neither G4 nor G5.

The original `six_chromosome_blocks_v1` evaluation geometry is invalid because its nearest-window cores covered only 84.35% of eligible validation bases in the first P6B job. `six_chromosome_blocks_v2` preserves block assignments and training data but requires one complete aligned metric core per valid/test context row and excludes every v1 P6B artifact from selection.

## Completion Standard

Completion requires an auditable v2 provenance chain, a common and justified RNA-seq scale, deterministic on-demand loading, leakage-buffered `I`/`II`/`III`/`IV`/`V`/`X` coverage in train, validation, and test, corrected AlphaGenome-style losses, fair comparison of frozen-head, worm-embedding/LoRA, and from-scratch models on that split, and one immutable six-chromosome final-test report.
