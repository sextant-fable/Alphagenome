# P9 Component-Attribution Figure

## Figure contract

The figure supports one bounded conclusion: under the frozen genomic-block
cross-validation protocol, full model B outperforms the size-matched
from-scratch model and major component removals; LoRA provides the largest
factorial contribution, the worm embedding provides a smaller contribution,
and the primary-score interval for the 1-bp-head-only variant includes zero.

This is an internal development comparison. It is not an independent
study/laboratory benchmark, a new final-test result, or a comparison against all
public models.

## Methods

The plotting script reads only the formal R9-PASS review and the frozen P9
paired-effect, factorial-effect and per-track TSV files. It does not load model
weights, interval manifests, BigWigs or final-test artifacts. All model effects
are paired by held-out genomic fold and seed and use the registered direction
`candidate minus B/paper`. The independent blocking unit is the held-out
genomic fold (`n = 5`); three seeds are nested algorithmic repeats within each
fold. Diamonds and horizontal intervals show fold-primary means and 95%
hierarchical percentile bootstrap intervals obtained from 10,000 resamples of
folds followed by seeds within each sampled fold.

Panel a shows the prespecified primary biological score for size-matched C,
frozen A, B without LoRA, LoRA only and the 1-bp-head-only model relative to
full B. Small circles are the five fold means. Panel b shows the registered
`2 x 2` worm-embedding main effect, LoRA main effect and worm-by-LoRA
interaction. The interaction is a statistical non-additivity contrast and does
not establish a mechanistic interaction. Panel c shows all 241 track-wise
primary-score effects after averaging seeds within fold and folds equally.
Tracks are dependent descriptive outcomes and are not treated as independent
replicates; medians and interquartile ranges summarize their distribution.

No p values or multiplicity-adjusted tests are displayed. The primary endpoint
is the principal inferential result; component diagnostics and track
distributions support interpretation. The 1-bp-head-only interval crossing zero
is not an equivalence result because no equivalence margin was prespecified.

## Ready-to-use legend

**Fig. X | Component attribution under matched genomic-block validation.**
**a,** Same-fold, same-seed paired effects on the primary biological score for
the size-matched from-scratch model C and component variants relative to full
model B. Small circles show seed-averaged fold effects; diamonds and horizontal
lines show fold-primary means and 95% hierarchical bootstrap confidence
intervals. The interval for the 1-bp-head-only model includes zero and does not
establish equivalence. **b,** Registered `2 x 2` factorial effects of the worm
embedding, LoRA and their interaction. The interaction represents statistical
non-additivity rather than a mechanistic interaction. **c,** Descriptive
distributions of track-wise paired primary-score effects across all 241 RNA-seq
tracks; central bars show medians and interquartile ranges. For a and b,
`n = 5` held-out genomic folds with three seeds nested within each fold;
confidence intervals use 10,000 fold-first, seed-within-fold bootstrap
resamples. Tracks in c are outcomes, not independent replicates. Higher effects
indicate better primary-score performance relative to B. Source data are
provided as a Source Data file.

## Reproduction and outputs

```bash
conda run -n CellUNetr python -m scripts.plot_v2_p9_submission_evidence
```

The ignored output directory is
`results/v2_p9_submission_evidence_figure/`. It contains editable PDF and SVG,
600-dpi PNG and TIFF, Figure Source Data, a machine-readable QA report and an
artifact manifest. A visual-review pass is recorded only after inspecting each
panel and the complete raster at final physical size.
