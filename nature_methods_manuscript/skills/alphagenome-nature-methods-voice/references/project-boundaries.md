# Project Evidence Boundaries

## Verified Internal Evidence

- The active v2 resource contains 482 verified RNA-seq runs grouped into 241
  biological RNA tracks on WBcel235.
- The selected full Model B combines a frozen AlphaGenome-style trunk, worm
  embedding, LoRA and a dual-resolution RNA head.
- P9 is a five-fold, three-seed component matrix. The full model exceeds the
  size-matched from-scratch Model C on the registered primary score; LoRA is the
  largest measured main effect and worm embedding gives a smaller positive main
  effect.
- The P6C result is a frozen within-collection evaluation. It is not an external
  benchmark or a project-wide pristine test claim.

## Protected Negative and Limited Evidence

- P10 DPY-27 compared DPY-27 RNAi with vector RNAi on development-validation
  blocks. Its predicted chromosome-level X-minus-autosome contrast had the
  opposite sign to the observed contrast. It is a negative internal diagnostic,
  not positive biological validation.
- The learned-1-bp-head-only comparison has a primary-score interval crossing
  zero. Do not write that the dual-resolution head is necessary.
- The variant-scoring CLI has validation of engineering behavior and synthetic
  tests, not a completed external functional-variant truth evaluation.

## External Evidence Placeholders

The manuscript plan assumes a completed external benchmark and external
biological application only as a structural input. Until their signed artifacts
exist, retain the exact `EXTERNAL_*` placeholders from
`../../../plans/08_external_benchmark_placeholder_contract.md`.

## Prohibited Reframings

- `external validation` for internal genomic-block validation.
- `biological validation` or `dosage-compensation recovery` for P10.
- `parameter matched` for the original P6B Model C.
- `causal variant` for a model score without a causal experiment.
- `independent reproduction` before an actual clean third-party or clean-machine
  reproduction record exists.
