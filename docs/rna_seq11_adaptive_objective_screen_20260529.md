# RNA-seq11 Adaptive Objective Screen, 2026-05-29

## Scope

This record summarizes the validation-only adaptive objective screen requested for the C. elegans RNA-seq11 AlphaGenome-style adapter. The held-out test split was not read or evaluated. All model selection below uses unweighted valid full-resolution log1p MSE. Weighted losses were used only as training objectives.

Host: `HY-GPU`

Branch: `setup/agent-maintenance`

Code commit: `71ee5dd` (`Add RNA-seq11 adaptive loss training options`)

Environment: conda environment `alphagenome`, Python `3.12.13`, PyTorch `2.11.0+cu128`, CUDA available.

Allowed GPUs used: `CUDA_VISIBLE_DEVICES=2` and `CUDA_VISIBLE_DEVICES=3` only.

## Code Changes

Phase 0 added the following training/evaluation capabilities:

- Signal-weighted loss presets: `high-gt3-x1.5`, `high-gt3-x2`, `high-gt3-x3`, `top5-x2`, `top5-x3`, `top1-x5`, and `intestine-high-x3`.
- Region-weighted loss presets: `exon-gene-v1` and `exon-gene-v2`, using train-batch intervals and the WBcel235 GTF.
- Gene-body mean auxiliary loss with `--gene-aux-loss-weight`.
- Linear output calibration modes: `track-affine` and `high-signal-gated`.
- Checkpoint loading support for old and calibrated checkpoints in training, evaluation, and diagnostics.

Edited scripts:

- `scripts/alphagenome_rna_seq11_adapter.py`
- `scripts/alphagenome_rna_seq11_finetune.py`
- `scripts/alphagenome_rna_seq11_eval.py`
- `scripts/rna_seq11_top15_extended_diagnostics.py`

## Verification

Code checks:

```bash
conda run -n alphagenome python -m py_compile \
  scripts/alphagenome_rna_seq11_adapter.py \
  scripts/alphagenome_rna_seq11_finetune.py \
  scripts/alphagenome_rna_seq11_eval.py \
  scripts/rna_seq11_top15_extended_diagnostics.py

git diff --check
```

Smoke run:

- Output: `runs/rna_seq11_phase0_code_smoke_20260529_sig_region_gene_calib_2steps`
- Purpose: exercise signal loss, region loss, gene auxiliary loss, and calibration loading/training for 2 steps on GPU 2.
- Status: completed; treated only as code/environment validation, not as a model result.

## Main Baselines

| Model | Valid full MSE | MAE | Pearson | Notes |
|---|---:|---:|---:|---|
| `rank1_baseline` | 0.83557710 | 0.58180543 | 0.69463443 | Current full-MSE winner |
| `rank3_baseline` | 0.83564435 | 0.58202288 | 0.69462896 | Slightly better gene MSE and pooled metrics |

Rank1 checkpoint:

`runs/rna_seq11_1bpB_broad_w12_track_hard15_beta05_from_w10best1000_lr1e-5_seedrand_1500steps/adapter_head_best.pt`

Rank3 checkpoint:

`runs/rna_seq11_1bpB_broad_w10_track_hard15_beta05_from_w6beta05best1500_lr1e-5_seedrand_1500steps/adapter_head_best.pt`

## Phase Summaries

### Phase 1: Signal-Weighted Loss

Diagnostics:

- `runs/rna_seq11_phase1_sigweight_diagnostics_20260529`

Best valid full-MSE results from the screen:

| Run | Best step | Full MSE | Pearson | Common128 MSE |
|---|---:|---:|---:|---:|
| `rank1_highgt3x2_lr3e-6` | 250 | 0.84762559 | 0.69586601 | 0.83972094 |
| `rank1_top1x5_lr3e-6` | 250 | 0.84503946 | 0.69598754 | 0.83944783 |
| `fallback_rank1_highgt3x15` | 250 | 0.84062108 | 0.69566999 | 0.83594618 |
| `fallback_rank1_top5x2` | 250 | 0.84182041 | 0.69578491 | 0.83705817 |

Conclusion: signal-weighted objectives improved high-signal amplitude metrics but consistently hurt unweighted full MSE. The strongest amplitude candidate was `rank1_sig_top1x5`, with high_gt3 MSE `4.48196560`, top5 MSE `5.29603140`, and top1 MSE `9.89913290`, but full MSE worsened to `0.84503945`. These runs are secondary diagnostics, not replacement models.

### Phase 2: Region-Weighted Loss

Diagnostics:

- `runs/rna_seq11_phase2_region_diagnostics_v1v2_20260529`

| Run | Best step | Full MSE | Pearson | Exon MSE | Gene-body MSE | Intergenic MSE |
|---|---:|---:|---:|---:|---:|---:|
| `rank1_region_exongene_v2` | 500 | 0.83667865 | 0.69444571 | 1.65913410 | 1.10832200 | 0.30478636 |
| `rank1_region_exongene_v1` | 250 | 0.83771714 | 0.69460584 | 1.65496390 | 1.10790070 | 0.30843821 |

Conclusion: region weighting improved exon and gene-body MSE, but shifted error into intergenic bases and worsened full MSE. `exon-gene-v2` is the gentler region secondary candidate; `exon-gene-v1` is the strongest exon/gene-body candidate.

### Phase 3: Gene-Level Auxiliary Loss

Diagnostics:

- `runs/rna_seq11_phase3_geneaux_diagnostics_20260529`

| Run | Best step | Full MSE | Gene MSE | Gene Pearson | Gene Spearman |
|---|---:|---:|---:|---:|---:|
| `rank1_geneaux_lambda003` | 250 | 0.83582663 | 0.71335302 | 0.70930106 | 0.58352324 |
| `rank1_regionv2_geneaux_lambda003` | 250 | 0.83701929 | 0.71504905 | 0.70927970 | 0.58411273 |

Conclusion: `lambda=0.03` did not improve gene-body mean MSE. Larger `lambda=0.1` and `0.2` were not run because the low-lambda result lacked a positive gene-MSE signal.

### Phase 4: Track-Affine Calibration

Final diagnostics:

- `runs/rna_seq11_phase5_final_diagnostics_20260529`

| Run | Best step | Full MSE | Pearson | high_gt3 MSE | top5 MSE | top1 MSE |
|---|---:|---:|---:|---:|---:|---:|
| `rank1_calib_trackaffine_lr3e5` | 250 | 0.83577257 | 0.69463217 | 4.91200900 | 5.84397740 | 10.81207700 |
| `rank3_calib_trackaffine_lr3e5` | 250 | 0.83581002 | 0.69462708 | 4.91356110 | 5.84573660 | 10.81660700 |
| `rank1_calib_trackaffine_lr1e4` | 250 | 0.83607037 | 0.69462573 | 4.89313570 | 5.82073130 | 10.77249100 |

Conclusion: track-affine calibration gave small amplitude improvements but no full-MSE improvement. Because the conservative calibration did not produce a positive full-MSE or strong high-signal tradeoff, high-signal-gated calibration was not run.

## Final Candidate Diagnostics

Final diagnostics were run on:

- `rank1_baseline`
- `rank3_baseline`
- `rank1_calib_trackaffine_lr3e5`
- `rank3_calib_trackaffine_lr3e5`
- `rank1_calib_trackaffine_lr1e4`
- `rank1_region_exongene_v2`
- `rank1_region_exongene_v1`
- `rank1_sig_top1x5`

Output:

- `runs/rna_seq11_phase5_final_diagnostics_20260529`

### Overall and Weakness Metrics

| Run | Full MSE | Pearson | high_gt3 MSE | top5 MSE | top1 MSE | Exon MSE | Gene-body MSE | Gene MSE | Local-gradient Pearson |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `rank1_baseline` | 0.83557710 | 0.69463443 | 4.92574770 | 5.86102660 | 10.84144000 | 1.66583930 | 1.11036130 | 0.71269831 | 0.34511063 |
| `rank3_baseline` | 0.83564435 | 0.69462896 | 4.92478850 | 5.85962170 | 10.84044900 | 1.66586480 | 1.11048360 | 0.71264155 | 0.34460654 |
| `rank1_calib_trackaffine_lr3e5` | 0.83577257 | 0.69463217 | 4.91200900 | 5.84397740 | 10.81207700 | 1.66537020 | 1.11031290 | 0.71270465 | 0.34510730 |
| `rank3_calib_trackaffine_lr3e5` | 0.83581002 | 0.69462708 | 4.91356110 | 5.84573660 | 10.81660700 | 1.66548420 | 1.11045100 | 0.71264447 | 0.34460253 |
| `rank1_calib_trackaffine_lr1e4` | 0.83607037 | 0.69462573 | 4.89313570 | 5.82073130 | 10.77249100 | 1.66475100 | 1.11026390 | 0.71274120 | 0.34509961 |
| `rank1_region_exongene_v2` | 0.83667865 | 0.69444571 | 4.85267690 | 5.79186450 | 10.76434700 | 1.65913410 | 1.10832200 | 0.71475527 | 0.34308571 |
| `rank1_region_exongene_v1` | 0.83771714 | 0.69460584 | 4.78285160 | 5.71188690 | 10.65141100 | 1.65496390 | 1.10790070 | 0.71569129 | 0.34249607 |
| `rank1_sig_top1x5` | 0.84503945 | 0.69598754 | 4.48196560 | 5.29603140 | 9.89913290 | 1.66674760 | 1.11666410 | 0.71833957 | 0.34171745 |

### Final Rankings

Full-MSE top candidates:

1. `rank1_baseline`: `0.83557710`
2. `rank3_baseline`: `0.83564435`
3. `rank1_calib_trackaffine_lr3e5`: `0.83577257`
4. `rank3_calib_trackaffine_lr3e5`: `0.83581002`
5. `rank1_geneaux_lambda003`: `0.83582663`

Gene-MSE top candidates from final diagnostics:

1. `rank3_baseline`: `0.71264155`
2. `rank3_calib_trackaffine_lr3e5`: `0.71264447`
3. `rank1_baseline`: `0.71269831`
4. `rank1_calib_trackaffine_lr3e5`: `0.71270465`
5. `rank1_calib_trackaffine_lr1e4`: `0.71274120`

High-signal top candidates by top1 MSE from final diagnostics:

1. `rank1_sig_top1x5`: `9.89913290`
2. `rank1_region_exongene_v1`: `10.65141100`
3. `rank1_region_exongene_v2`: `10.76434700`
4. `rank1_calib_trackaffine_lr1e4`: `10.77249100`
5. `rank1_calib_trackaffine_lr3e5`: `10.81207700`

Exon/gene-body top candidates:

1. `rank1_region_exongene_v1`: exon `1.65496390`, gene-body `1.10790070`
2. `rank1_region_exongene_v2`: exon `1.65913410`, gene-body `1.10832200`
3. `rank1_calib_trackaffine_lr1e4`: exon `1.66475100`, gene-body `1.11026390`
4. `rank1_calib_trackaffine_lr3e5`: exon `1.66537020`, gene-body `1.11031290`
5. `rank3_calib_trackaffine_lr3e5`: exon `1.66548420`, gene-body `1.11045100`

## Final Decision

Main winner remains:

`runs/rna_seq11_1bpB_broad_w12_track_hard15_beta05_from_w10best1000_lr1e-5_seedrand_1500steps/adapter_head_best.pt`

Reason: it still has the lowest unweighted valid full-resolution MSE, `0.83557710`.

Co-winner / reference candidate:

`runs/rna_seq11_1bpB_broad_w10_track_hard15_beta05_from_w6beta05best1500_lr1e-5_seedrand_1500steps/adapter_head_best.pt`

Reason: it is essentially tied on full MSE and remains the best gene-MSE candidate in final diagnostics.

Secondary analysis candidates:

- Best high-signal amplitude tradeoff: `rank1_sig_top1x5`; not promoted because full MSE worsens by about `+0.00946`.
- Best exon/gene-body candidate: `rank1_region_exongene_v1`; not promoted because full MSE worsens by about `+0.00214`.
- Best new full-MSE candidate: `rank1_calib_trackaffine_lr3e5`; not promoted because full MSE worsens by about `+0.00020`.

No test-split evaluation was performed.

## Next Actions

- Do not replace the current rank1 full-MSE winner.
- Keep rank3 as a biologically useful co-reference for gene-level comparisons.
- Treat high-signal and region-weighted losses as diagnostic levers: they expose real tradeoffs but do not yet improve the selection metric.
- A future attempt should likely change model capacity or target parameterization for high-signal amplitude and local boundaries, rather than increasing loss weights further on the same plateaued residual head.
