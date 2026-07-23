# C. elegans RNA-seq v2 Final Report

Status: `completed_one_time_final_test`

Execution ID: `p6c-20260723T045025Z-8f9262a0ef92`

## Data

The auditable input set contains 485 runs: {'2026.4.20': 20, '2026.5.18': 26, 'sharepoint_2026-07-13_new439': 439}. Uniform reprocessing retained 482 RNA-seq runs, excluded 3 non-RNA runs, and produced 241 grouped tracks at a target total of 100000000 (1e6 times 100-bp coverage).

The model reads these tracks on demand. No monolithic v2 NPZ was generated. Each of I, II, III, IV, V, and X contributes leakage-buffered train, validation, and locked final-test blocks. The final blocks were consumed once by the locked final execution.

## Model Selection

The formal matrix selected model B with paper loss over 5 folds, 3 seeds, and 15 jobs. Its development primary biological score was 0.626602.

## Final Test

Gene-exon coverage Pearson (log1p): 0.675019

Per-track 128-bp Pearson (log1p): 0.664130

Paper loss: 15.732128; log1p MSE: 2.943107.

## Limitations

- Chromosome X was viewed in legacy exploratory work, so its locked test block is not fully pristine; every final test block remains one-time-use within the revised v2 workflow.
- All 485 source runs were class C before reprocessing; 482 RNA-seq runs were uniformly rebuilt and three non-RNA runs were excluded.
- The augmentation and gene-loss ablations used one seed across five folds and are directional rather than high-power estimates.
- The selected B/paper configuration won the preregistered biological primary score but need not dominate every secondary metric.
- The one-time six-chromosome locked-block result must be reported regardless of outcome and cannot be used for further tuning in this study cycle.

## Reproduction

The exact command, code commit, input hashes, complete metric report, formal matrix, ablations, and failure record are stored in `alphagenome_custom/metadata/v2/v2_final_report.json`.
