# Terminology Ledger

Lock this ledger before the first manuscript draft. Any author-approved method
name will replace the temporary labels below everywhere at once.

| Canonical term | First-use definition | Do not use |
| --- | --- | --- |
| *Caenorhabditis elegans* | Italicize full species at first use; use `C. elegans` subsequently. | `worm` in formal result claims unless context makes it unambiguous. |
| sequence-to-RNA-track model | A model that maps DNA sequence to RNA-seq coverage tracks. | `expression predictor` when coverage is the output. |
| RNA-seq coverage | The continuous genomic coverage target derived from uniformly processed RNA-seq. | `RNA expression` without defining aggregation. |
| full Model B | Frozen trunk plus worm embedding, LoRA and dual-resolution RNA head. | `our model`, `enhanced model`, `B model` without definition. |
| Model A | Frozen-trunk RNA-head baseline. | `base model` once a public baseline is discussed. |
| size-matched Model C | P9 from-scratch control matched to Model B trainable parameters. | `Model C` alone when parameter matching matters. |
| worm embedding | Learnable organism-conditioning embedding for the C. elegans adaptation. | `species token` unless the implementation changes. |
| LoRA | Low-rank adaptation in the registered target modules. | `fine-tuning` when only low-rank parameters were optimized. |
| dual-resolution RNA head | Joint 1-bp and 128-bp RNA prediction head. | `multi-scale head` before it is defined. |
| genomic-block validation | Evaluation on held-out genomic cores inside the source collection. | `external validation`. |
| locked final evaluation | One-time v2-gated final evaluation reported in the repository. | `pristine external test`. |
| study-isolated benchmark | Benchmark held out by declared study/laboratory source. | `OOD` without a precise isolation definition. |
| unseen biological condition | Condition absent from development as defined in the frozen external manifest. | `novel condition` without its identity. |
| primary composite score | Registered development endpoint, with formula and direction specified in Methods. | `accuracy` as a catch-all. |
| paired fold effect | Difference between matched configurations, summarized with folds as primary blocked units and seeds nested within fold. | `n = 241 tracks` as inferential sample size. |
| track-wise distribution | Descriptive distribution across 241 RNA tracks. | `independent per-track replication`. |
| regulatory-variant score | Frozen ref/alt or haplotype prediction difference aggregated to a stated gene-level statistic. | `causal variant score`. |
| DPY-27 internal diagnostic | P10 DPY-27 RNAi versus vector RNAi developmental-validation comparison. | `biological validation`, `dosage-compensation recovery`. |
