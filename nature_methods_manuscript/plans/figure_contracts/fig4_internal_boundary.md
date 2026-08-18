# Figure Contract: DPY-27 Internal Boundary Diagnostic

## Claim

The DPY-27 internal analysis tests a condition contrast that the present model
does not recover at the chromosome level. It defines a use-domain boundary; it
does not provide biological validation.

## Inputs and Invariance

- Source data: `runs/v2_p10_dpy27_internal_application_20260813/figure_source_data.tsv`
- Formal review: `alphagenome_custom/metadata/v2/audits/P10/review.json`
- No prediction, bigWig, checkpoint or final-test input is read.
- The derived figure changes only display titles from recovery language to
  neutral contrast language. Coordinates, points, intervals and source data
  remain identical to the R10-PASS formal figure.

## Panels

| Panel | Question | Interpretation |
| --- | --- | --- |
| a | Do observed and predicted gene contrasts agree? | Weak gene-level diagnostic, not a condition-generalization test. |
| b | Is the X-minus-autosome direction recovered? | Primary internal boundary endpoint; observed and predicted signs disagree. |
| c | What secondary agreement is present? | Descriptive fold endpoints with blocked uncertainty. |

## Output Contract

- Python/matplotlib, 183 mm x 82 mm, vector PDF/SVG and 600 dpi PNG/TIFF.
- Manifest records the formal source-data and review SHA-256 values.
- Generated under `results/v2_p10_dpy27_boundary_figure/`; it does not replace
  the formal P10 artifact.
