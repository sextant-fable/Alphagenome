# P6B submission analysis

The existing P6B development-validation matrix has now been analysed as a complete matched design: A/B/C x paper/log1p MSE x five folds x three seeds (90 runs). The analysis is locked by `alphagenome_custom/metadata/v2/p6b_submission_analysis_spec.json`, uses folds as the primary bootstrap blocks, resamples seeds within sampled folds, and does not read the locked final test.

## Main result

Model B with paper loss achieved a composite biological score of 0.626602 (fold-primary hierarchical bootstrap 95% CI 0.611447-0.640932). Its paired composite advantage was 0.111180 over A/paper (95% CI 0.107264-0.116186), 0.139057 over C/paper (0.132762-0.145395), and 0.011228 over B/log1p MSE (0.004970-0.016700).

The objective comparison has an important, interpretable trade-off. Relative to B/log1p MSE, B/paper improved gene-exon coverage Pearson by 0.047952 (95% CI 0.041355-0.054702) while reducing 128-bp Pearson by 0.025497 (95% CI -0.032343 to -0.020178). The pre-specified composite favors B/paper because its gene-level gain is larger than its 128-bp loss. This supports a biological-objective claim, not uniform dominance on every endpoint.

Track-level analyses retain all 241 tracks and are descriptive rather than inferential. Full study, stage, tissue, treatment and development-signal strata are available in `stratified_track_summaries.tsv`; missing/unknown values are retained. Signal strata contain 81 low, 80 middle and 80 high tracks.

## Reproducibility

- Analysis entry point: `scripts/analyze_v2_p6b_submission.py`
- Generic statistics API for P9 reuse: `scripts/v2_submission_statistics.py`
- Focused tests: `tests/test_v2_submission_statistics.py`
- Locked analysis plan: `alphagenome_custom/metadata/v2/p6b_submission_analysis_spec.json`
- Generated bundle: `results/v2_p6b_submission_analysis/`
- Manifest and checksums: `results/v2_p6b_submission_analysis/analysis_manifest.json`
- Methods and legend: `results/v2_p6b_submission_analysis/METHODS.md` and `FIGURE_LEGEND.md`
- Figure: PDF, SVG, 600-dpi TIFF and 300-dpi PNG under the generated bundle
- Figure Source Data: CSV and TSV under the generated bundle

The extraction/statistical path runs in the `alphagenome` environment with standard-library modules and NumPy only. Rendering used the existing `CellUNetr` Python environment with Python 3.9.16, NumPy 1.21.0, matplotlib 3.7.0 and Pillow 9.5.0; no packages were installed.

## Interpretation boundary

This analysis materially strengthens the internal A/B/C comparison and makes the current evidence submission-ready. It does not replace the requested independent study/laboratory benchmark. Because training sets overlap across cross-validation folds, the fold-primary intervals quantify variability under the locked development split and should not be described as independent-study generalization intervals.
