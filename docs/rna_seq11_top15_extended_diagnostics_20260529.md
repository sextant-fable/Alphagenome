# RNA-seq11 Top-15 Extended Diagnostics - 2026-05-29

## Scope

- Split: valid only (`V` chromosome); test split was not read.
- Ranking: valid full-resolution log1p MSE ascending, Pearson descending tie-breaker.
- Selected models: top 15 existing checkpoints under `runs/`, excluding smoke and diagnostics runs.
- Output directory: `runs/rna_seq11_top15_extended_diagnostics_20260529`
- Diagnostic script: `scripts/rna_seq11_top15_extended_diagnostics.py`
- Run host: `HY-GPU`
- Run start Git HEAD: `be22349`; the diagnostic script and this report were added after the run for reproducibility.

## Output Tables

- `top15_models.tsv`: selected checkpoint list and source valid MSE/Pearson.
- `per_track_metrics.tsv`: full-resolution per-track MSE, MAE, Pearson, sampled Spearman.
- `signal_strata_metrics.tsv`: zero, low, medium, high, top 5%, and top 1% true-signal strata.
- `gene_metrics.tsv`: GTF gene-body mean expression MSE, MAE, Pearson, Spearman.
- `region_metrics.tsv`: TSS/promoter, gene body, exon, intron, and intergenic metrics.
- `high_signal_localization_metrics.tsv`: 128 bp pooled top 5% and top 1% overlap.
- `resolution_metrics.tsv`: 1 bp, 128 bp, 512 bp, 1 kb pooled metrics plus local-gradient metrics.
- `window_metrics.tsv`: per-valid-window full-resolution metrics.

## Top-15 Source Ranking

| Rank | Valid MSE | Pearson | Method summary |
|---:|---:|---:|---|
| 1 | 0.83557710 | 0.69463443 | 1bp-B residual correction, conv5 bottleneck128 hidden128, hard15, beta0.5, lr1e-5, from w10 best1000 |
| 2 | 0.83563143 | 0.69445160 | same family, hard15, beta0.5, lr3e-6, from w10 best1000 |
| 3 | 0.83564435 | 0.69462896 | same family, hard15, beta0.5, lr1e-5, from w6 beta0.5 best1500 |
| 4 | 0.83573077 | 0.69436124 | same family, hard15, beta0.5, lr1e-6, from w10 best750 |
| 5 | 0.83574712 | 0.69432326 | same family, uniform weights, beta0.5, lr3e-6 |
| 6 | 0.83582169 | 0.69423686 | same family, hard15, beta0.5, lr1e-5, from w5 best |
| 7 | 0.83583516 | 0.69431161 | same family, hard15, beta0.5, lr3e-6 |
| 8 | 0.83596947 | 0.69424414 | same family, hard15, beta0.5, lr3e-6, from w6 best750 |
| 9 | 0.83599538 | 0.69421002 | same family, hard15, beta0.5, lr1e-6, from w6 best750 |
| 10 | 0.83599816 | 0.69420599 | same family, hard15, beta0.5, lr3e-6, from w6 best1000 |
| 11 | 0.83605101 | 0.69417615 | same family, hard15, beta0.5, lr3e-6 |
| 12 | 0.83638968 | 0.69386931 | same family, hard15, hybrid alpha0.25, beta1, lr1e-5 |
| 13 | 0.83688281 | 0.69407776 | same family, hard15, beta1, lr1e-5, from w4 best2000 |
| 14 | 0.83691050 | 0.69398340 | same family, uniform weights, beta1, lr1e-5 |
| 15 | 0.83694736 | 0.69398317 | same family, hard15, beta1, lr1e-5, from w4 best1500 |

All top-15 models are the same broad method: frozen 128 bp winner prediction plus a 1 bp residual-correction head. The ranking is therefore mostly comparing continuation point, learning rate, beta, and track-weight preset rather than fundamentally different architectures.

## Main Result

The best model is rank 1:

- Full-resolution valid MSE: `0.83557710`
- MAE: `0.58180543`
- Pearson: `0.69463443`
- Common 128 bp MSE: `0.83737121`

Compared with the prior 128 bp conv5 hybrid step4000 baseline (`0.88032581` MSE, `0.67391665` Pearson), rank 1 improves valid MSE by `0.04474871` and Pearson by `0.02071778`. Every track improves versus that 128 bp baseline.

Largest per-track MSE reductions versus the 128 bp baseline:

| Track | 128 bp MSE | Rank 1 MSE | Delta |
|---|---:|---:|---:|
| Intestine T4 | 1.081604 | 1.025233 | -0.056370 |
| Muscle T4 | 0.860426 | 0.805450 | -0.054977 |
| Muscle T3 | 0.781381 | 0.732684 | -0.048696 |
| Intestine T3 | 1.054880 | 1.008271 | -0.046609 |

## Remaining Weak Spots

Per-track error is still dominated by intestine tracks:

| Track | Rank 1 MSE | Pearson | Sampled Spearman |
|---|---:|---:|---:|
| Intestine T4 | 1.025234 | 0.659890 | 0.504577 |
| Intestine T3 | 1.008271 | 0.643327 | 0.471969 |
| Intestine T1 | 0.981259 | 0.665036 | 0.490107 |
| Intestine T0 | 0.881485 | 0.659681 | 0.467845 |

The strongest muscle tracks are much lower MSE: Muscle T1 `0.688205`, Muscle T2 `0.705710`, Muscle T3 `0.732684`.

Signal-stratum diagnostics show the main failure is high-signal modeling, not just low background:

| Stratum | MSE | MAE | Pearson | Values |
|---|---:|---:|---:|---:|
| zero | 0.319119 | 0.373289 | nan | 271,275,460 |
| low 0.1-1 | 0.419998 | 0.431020 | 0.134069 | 78,969,624 |
| medium 1-3 | 1.349904 | 0.985095 | 0.282627 | 65,561,006 |
| high >3 | 4.925748 | 1.816868 | 0.351440 | 34,033,013 |
| top 5% true | 5.861027 | 2.001077 | 0.303556 | 22,623,155 |
| top 1% true | 10.841440 | 2.876432 | 0.359815 | 4,507,130 |

The zero/low regions have much smaller MSE, while top 1% signal has very large amplitude error. This means global MSE can look stable while high-expression loci remain poorly calibrated.

Region diagnostics point to the same issue:

| Region | MSE | MAE | Pearson |
|---|---:|---:|---:|
| exon | 1.665839 | 0.889688 | 0.685561 |
| gene body | 1.110361 | 0.695507 | 0.691554 |
| promoter/TSS +/-1kb | 0.790466 | 0.558717 | 0.688075 |
| intron | 0.449854 | 0.464611 | 0.414054 |
| intergenic | 0.297975 | 0.360765 | 0.367563 |

Exons and gene bodies carry most of the residual error. Intergenic MSE is low but correlation is also low, so the model is mostly helped there by low signal magnitude.

Gene-level expression is better than pointwise correlation but still limited:

- Overall gene-body mean MSE: `0.71269831`
- Gene-level Pearson: `0.70926062`
- Gene-level Spearman: `0.58345970`

Worst gene-level tracks are again Intestine T1/T4/T3.

High-signal localization at 128 bp pooled resolution:

| Metric | Precision/Recall | Jaccard | True threshold | Pred threshold |
|---|---:|---:|---:|---:|
| top 5% overlap | 0.539532 | 0.369424 | 3.504081 | 2.791328 |
| top 1% overlap | 0.362252 | 0.221189 | 5.172837 | 4.210881 |

The prediction thresholds are lower than true thresholds, especially in the top 1%, which is consistent with peak/high-expression underestimation.

Resolution-specific metrics:

| Scope | Pool bp | MSE | Pearson |
|---|---:|---:|---:|
| pooled value | 1 | 0.835577 | 0.694634 |
| pooled value | 128 | 0.720154 | 0.709227 |
| pooled value | 512 | 0.610981 | 0.715742 |
| pooled value | 1024 | 0.545147 | 0.719211 |
| local gradient | 1 | 0.005857 | 0.345111 |
| local gradient | 128 | 0.205521 | 0.650935 |
| local gradient | 512 | 0.269203 | 0.689431 |
| local gradient | 1024 | 0.296807 | 0.659439 |

Pooling improves MSE and Pearson steadily, so the model captures broad expression better than fine positional variation. The full-resolution gradient Pearson is only `0.345`, indicating weak base-scale/local-boundary behavior.

## Cross-Metric Ranking Notes

Rank 1 is the best by full MSE. Rank 3 is essentially tied and is best by several secondary metrics:

- Best gene MSE: rank 3, `0.71264155`
- Best gene Pearson: rank 3, `0.70927212`
- Best 128 bp pooled MSE: rank 3, `0.72012134`
- Best 1 kb pooled MSE: rank 3, `0.54513118`

High-signal overlap slightly favors lower-ranked beta1/uniform variants:

- Best top 5% overlap: rank 15, precision `0.53968552`
- Best top 1% overlap: rank 14, precision `0.36304917`

These differences are tiny. The top-15 spread is only about `0.00137` full MSE, so these are near-ties rather than strong evidence for a new direction.

## Interpretation

The 1 bp residual correction is genuinely useful: it improves every track over the 128 bp conv5 hybrid baseline and gives the current best valid MSE/Pearson. However, the top-15 set has converged into a narrow plateau. The remaining weakness is not broad background prediction; it is high-signal amplitude, exon/gene-body expression, and intestine-specific tracks.

The most likely next useful changes are targeted objectives, not another small learning-rate continuation:

- Add high-signal or top-quantile weighting while keeping unweighted valid MSE for selection.
- Add a gene-body auxiliary loss or gene-level calibration diagnostic during validation.
- Try exon/gene-body weighted loss using GTF masks, with careful valid-only monitoring.
- Add peak/high-expression calibration metrics such as predicted-vs-true top quantile slope.
- Keep rank 1 as the MSE winner; keep rank 3 as a co-best candidate for gene-level and pooled-resolution behavior.

Claim status: verified on valid split only.
