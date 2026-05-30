# RNA-seq11 All-Candidate Representation Utility Benchmark

Date: 2026-05-30/31

Branch: `setup/agent-maintenance`

Host: `HY-GPU`

Claim status: verified for train/valid evaluation outputs listed below. No held-out test split was read.

## Scope

This benchmark reframes the RNA-seq11 checkpoint review away from a single pixel-level valid MSE criterion and toward biological signal utility. The full-resolution unweighted valid MSE is kept as a sanity reference, but it is not used as the only promotion rule.

The inventory contains 171 eligible checkpoints with `adapter_head_best.pt` and recorded valid metrics under `runs/`.

Primary outputs:

- Inventory: `runs/rna_seq11_all_representation_inventory_20260530/top15_models.tsv`
- Extended valid diagnostics: `runs/rna_seq11_all_representation_diagnostics_20260530/shard0` and `shard1`
- Downstream gene-profile probe: `runs/rna_seq11_gene_profile_probe_all_20260530/shard0` and `shard1`
- Utility summary and leaderboards: `runs/rna_seq11_representation_utility_summary_20260530`

## Policy Checks

- Split usage: train chromosomes I/II/III/IV and valid chromosome V only.
- Test split: not used.
- GPU policy: all GPU work used `CUDA_VISIBLE_DEVICES=2` or `CUDA_VISIBLE_DEVICES=3`; GPUs 0/1 were not used.
- Large artifacts: generated `runs/` and `logs/` outputs remain untracked.
- Candidate count: 171 unique run IDs in the final probe and summary tables.

## Methods

The benchmark has three evidence layers.

1. Extended valid diagnostics for every candidate:
   - full-resolution MSE/MAE/Pearson and sampled Spearman
   - high-signal strata (`high_gt3`, `top5_true`, `top1_true`)
   - 128 bp high-signal overlap
   - exon, gene-body, promoter, intron, intergenic region metrics
   - gene-body mean metrics
   - 128 bp and 1 kb pooled metrics
   - local-gradient metrics as boundary reference only

2. Downstream gene-profile probe for every candidate:
   - predicted and true gene-level 11-track profiles on train and valid splits
   - valid top expressed gene retrieval at top 1%, 5%, and 10%
   - Intestine T1/T3/T4 top expressed gene retrieval
   - cross-track correlation preservation
   - lightweight linear probes trained on train chromosomes and evaluated on chr V for family, time, and dominant-track labels

3. Rank-based utility summary:
   - `gene_utility_score`
   - `high_signal_utility_score`
   - `region_utility_score`
   - `probe_utility_score`
   - `representation_utility_score`

The composite score is a decision aid, not a biological truth claim. It averages rank-normalized evidence and applies only a mild full-MSE sanity penalty when a model exceeds the best full MSE by more than 0.02.

## Key Winners

### MSE Reference

The full-MSE reference remains:

`rna_seq11_1bpB_broad_w12_track_hard15_beta05_from_w10best1000_lr1e-5_seedrand_1500steps`

- full MSE: 0.83557710
- Pearson: 0.69463443
- gene MSE/Spearman: 0.71269831 / 0.58345970
- high_gt3/top5/top1 MSE: 4.9257477 / 5.8610266 / 10.84144
- exon/gene_body MSE: 1.6658393 / 1.1103613
- 128 bp / 1 kb pooled MSE: 0.72015391 / 0.54514746
- probe top5 recall: 0.49907236
- Intestine T1/T3/T4 probe top5 recall: 0.46841594

### Representation Utility Winner

The top composite representation winner is:

`rna_seq11_1bpB_region_w1_rank1_exongene_v2_lr1e-5_1000steps`

- representation utility score: 0.88814617
- full MSE: 0.83667865
- Pearson: 0.69444571
- gene MSE/Spearman: 0.71475527 / 0.58397658
- high_gt3/top5/top1 MSE: 4.8526769 / 5.7918645 / 10.764347
- exon/gene_body MSE: 1.6591341 / 1.1083220
- 128 bp / 1 kb pooled MSE: 0.72127108 / 0.54660438
- probe top5 recall: 0.49986748
- Intestine T1/T3/T4 probe top5 recall: 0.46938776
- cross-track correlation Spearman: 0.87777778

Interpretation: this is the best current main biological-utility candidate. It gives up only +0.00110 full MSE relative to rank1 while improving the composite of gene, high-signal, region, and probe evidence.

### Region / Exon-Gene Winner

The strongest region-specific candidate is:

`rna_seq11_1bpB_region_w1_rank1_exongene_v1_lr1e-5_1500steps`

- region utility score: 1.00000000
- full MSE: 0.83771714
- exon/gene_body MSE: 1.6549639 / 1.1079007
- high_gt3/top5/top1 MSE: 4.7828516 / 5.7118869 / 10.651411

Interpretation: v1 is the clean exon/gene-body winner. v2 is the better all-around representation winner because it has a better probe score and lower full MSE.

### High-Signal Amplitude Winner

The best pure high-signal amplitude candidate is:

`rna_seq11_1bpB_sigweight_w1_rank1_top5x3_a075_b05_lr1e-5_1500steps`

- high_gt3/top5/top1 MSE: 4.1604688 / 4.9060963 / 9.2685486
- full MSE: 0.86322585
- Pearson: 0.69495921

The similar rank3 top5x3 run is also strong:

`rna_seq11_1bpB_sigweight_w1_rank3_top5x3_a075_b05_lr1e-5_1500steps`

- high_gt3/top5/top1 MSE: 4.1873592 / 4.9452193 / 9.3397586
- full MSE: 0.86119161
- Pearson: 0.69573012
- probe top5 recall: 0.50251789

Interpretation: signal weighting really does improve amplitude and top-gene retrieval, but it overcorrects enough to hurt full MSE, exon/gene-body, and pooled MSE. These are high-signal co-winners, not main all-purpose winners.

### Intestine T1/T3/T4 Winners

For Intestine T1/T3/T4 high-signal MSE, the best candidates are intestine/high-signal weighted:

- `rna_seq11_1bpB_sigweight_w1_rank3_intestinehighx3_a075_b05_lr1e-5_1500steps`: Intestine T1/T3/T4 high_gt3 MSE 4.4585680
- `rna_seq11_1bpB_sigweight_w1_rank1_intestinehighx3_a075_b05_lr1e-5_1500steps`: Intestine T1/T3/T4 high_gt3 MSE 4.4590687

For Intestine T1/T3/T4 gene rank preservation, the best candidates are signal-weighted:

- `rna_seq11_1bpB_sigweight_w1_rank1_top5x3_a075_b05_lr1e-5_1500steps`: Intestine T1/T3/T4 gene Spearman 0.55130631
- `rna_seq11_1bpB_sigweight_w1_rank1_intestinehighx3_a075_b05_lr1e-5_1500steps`: Intestine T1/T3/T4 gene Spearman 0.55128909

Interpretation: if the biological question is specifically intestine high-expression ordering, the signal-weighted candidates deserve co-winner status despite worse full MSE.

### Downstream Probe Winner

The top aggregate probe utility candidate is:

`rna_seq11_1bpB_broad_w5_track_hard15_from_w4best2000_lr1e-5_seedrand_1500steps`

- probe utility score: 0.85808824
- full MSE: 0.83688281
- probe top5 recall: 0.49854227
- Intestine T1/T3/T4 probe top5 recall: 0.46938776
- cross-track correlation Spearman: 0.87683983
- family/time/track balanced accuracy: 0.42509607 / 0.25960569 / 0.15634862

Interpretation: this is the best current representation-probe co-winner among near-frontier full-MSE models. The probe differences are small, so it should be treated as a useful co-reference rather than a decisive replacement.

## Important Observations

High-signal amplitude is still underpredicted by all frontier models. The rank1 reference has average predicted/true threshold ratios of about 0.799 for top5 and 0.811 for top1. Signal-weighted top5x3 improves these ratios to about 0.876 and 0.869, respectively.

High-signal localization and high-signal amplitude are not the same objective. Some older models have slightly better top-bin overlap but much worse full MSE/Pearson. Signal-weighted models most clearly fix amplitude; region-weighted models give the best composite utility.

The pure best gene-body mean MSE is from an older 128 bp phase4 constant schedule run (`0.70766572`) with full MSE `0.88462103`. This supports keeping a separate gene-expression reference, but it is not the all-around representation winner.

The local-gradient metric remains a boundary reference only. It was not included in the main composite score.

## Recommended Winner Set

- MSE reference: `rna_seq11_1bpB_broad_w12_track_hard15_beta05_from_w10best1000_lr1e-5_seedrand_1500steps`
- Main representation winner: `rna_seq11_1bpB_region_w1_rank1_exongene_v2_lr1e-5_1000steps`
- Exon/gene-body winner: `rna_seq11_1bpB_region_w1_rank1_exongene_v1_lr1e-5_1500steps`
- High-signal amplitude winner: `rna_seq11_1bpB_sigweight_w1_rank1_top5x3_a075_b05_lr1e-5_1500steps`
- Intestine high-signal/rank winner: `rna_seq11_1bpB_sigweight_w1_rank1_intestinehighx3_a075_b05_lr1e-5_1500steps`
- Probe co-winner: `rna_seq11_1bpB_broad_w5_track_hard15_from_w4best2000_lr1e-5_seedrand_1500steps`

## Reproducibility Notes

Representative commands:

```bash
conda activate alphagenome
cd /home/zelinli6/Alphagenome

conda run -n alphagenome python scripts/rna_seq11_top15_extended_diagnostics.py \
  --runs-dir runs \
  --output-dir runs/rna_seq11_all_representation_inventory_20260530 \
  --top-n 999 \
  --only-write-top15

CUDA_VISIBLE_DEVICES=2 PYTHONUNBUFFERED=1 \
conda run -n alphagenome python -u scripts/rna_seq11_top15_extended_diagnostics.py \
  --runs-dir runs \
  --output-dir runs/rna_seq11_all_representation_diagnostics_20260530/shard0 \
  --top-n 999 \
  --num-shards 2 \
  --rank-shard 0 \
  --device auto \
  --batch-size 1 \
  --num-workers 0 \
  --spearman-sample-size 300000

CUDA_VISIBLE_DEVICES=3 PYTHONUNBUFFERED=1 \
conda run -n alphagenome python -u scripts/rna_seq11_top15_extended_diagnostics.py \
  --runs-dir runs \
  --output-dir runs/rna_seq11_all_representation_diagnostics_20260530/shard1 \
  --top-n 999 \
  --num-shards 2 \
  --rank-shard 1 \
  --device auto \
  --batch-size 1 \
  --num-workers 0 \
  --spearman-sample-size 300000

CUDA_VISIBLE_DEVICES=<2-or-3> PYTHONUNBUFFERED=1 \
conda run -n alphagenome python -u scripts/rna_seq11_gene_profile_probe.py \
  --candidate-tsv runs/rna_seq11_all_representation_inventory_20260530/top15_models.tsv \
  --top-n 999 \
  --output-dir runs/rna_seq11_gene_profile_probe_all_20260530/shard<0-or-1> \
  --rank-shard <0-or-1> \
  --num-shards 2 \
  --batch-size 1 \
  --num-workers 0 \
  --probe-steps 300 \
  --device auto \
  --eval-mode shared-embeddings \
  --skip-existing

conda run -n alphagenome python scripts/rna_seq11_representation_utility_summary.py \
  --inventory-tsv runs/rna_seq11_all_representation_inventory_20260530/top15_models.tsv \
  --diagnostics-dirs \
    runs/rna_seq11_all_representation_diagnostics_20260530/shard0 \
    runs/rna_seq11_all_representation_diagnostics_20260530/shard1 \
  --probe-dirs \
    runs/rna_seq11_gene_profile_probe_all_20260530/shard0 \
    runs/rna_seq11_gene_profile_probe_all_20260530/shard1 \
  --output-dir runs/rna_seq11_representation_utility_summary_20260530
```

Verification:

- `python -m py_compile` passed for the new probe and summary scripts.
- Shared-embedding probe smoke completed with `top-n 2`, two train windows, two valid windows, and two candidate outputs.
- Final extended diagnostics row counts: `per_track_metrics.tsv` shard0 1033 lines, shard1 1021 lines.
- Final probe summary row counts: shard0 87 lines, shard1 86 lines.
- Final probe per-track row counts: shard0 2839 lines, shard1 2806 lines.
- Final utility summary row count: 172 lines including header, representing 171 candidates.

