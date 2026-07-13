# Legacy RNA-seq11 v1 Freeze

Date: 2026-07-13T14:39:57+00:00

Host: `hy8`

Git commit: `8b67b775c88002434d952d4ba0ebb86b0fbd3642` on `setup/agent-maintenance`

## Scope

This freeze records the legacy 2026-04-20 C. elegans RNA-seq11 data and model-development lineage without copying or modifying large assets. The authoritative machine-readable inventories are under `alphagenome_custom/metadata/legacy_v1/`.

## Frozen Evidence

- Original sample-level bigWigs: `20` actual files under `alphagenome_custom/tracks/bw_2026.4.20.training_input_bigwig`.
- Legacy grouped RNA-seq tracks: `11`.
- Intervals: train `116`, valid `39`, test `33`.
- Eligible validation-tuned checkpoint inventory: `171`.
- Raw-data before/after hash agreement during freeze: `True`.

## Interpretation Boundary

The legacy pipeline directly averaged existing bigWig scales and subsequently trained custom log1p regression heads on frozen human-index AlphaGenome PyTorch embeddings. It is an exploratory transfer-learning result, not a paper-faithful C. elegans AlphaGenome reproduction.

Chromosome V was reused for extensive checkpoint and objective selection. Chromosome X was read for the early selected 128 bp adapter and related GenomeTracks audits. The later 171-candidate representation benchmark did not read chromosome X, but chromosome X is not considered pristine for the project as a whole.

No new legacy training or dataset generation is permitted after this freeze. Corrections must be appended to the v2 execution log and must not rewrite the inventories silently.
