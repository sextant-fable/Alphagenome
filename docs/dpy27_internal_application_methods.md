# DPY-27 internal genomic-block application

## Scope and interpretation

This application tests whether the 15 existing Model B / paper-loss blocked-CV
checkpoints recover the DPY-27 RNAi versus vector RNAi expression contrast on
their matching validation cores. It is labelled **internal genomic-block
validation** throughout. `RNA_V2_G0005` and `RNA_V2_G0006` were among the 241
training targets, so this analysis is not independent-study, unseen-condition,
causal, or chromosome-X final-test validation.

Both grouped tracks were independently rescaled to a decoded total signal of
1e8. The analysis therefore supports relative contrasts within this
total-signal-normalized space. It cannot establish an absolute change in total
RNA output without an external spike-in or equivalent calibration.

## Frozen inputs and execution

The complete contract is frozen in
`alphagenome_custom/metadata/v2/p10_dpy27_internal_application_spec.json`.
It locks all 15 checkpoint and run-record hashes, five validation manifests,
the two grouped BigWigs, 241-track order and fold means, WBcel235 FASTA/FAI/GTF,
AlphaGenome weights, and the completed final-test lock/report hashes. It does
not name or resolve a final-test interval file.

P10 is launched only by the controller after exact G4 scope
`p10_dpy27_internal_application` is active:

```bash
/home/zelinli6/miniconda3/envs/alphagenome/bin/python \
  -m scripts.v2_phase_controller run --phase P10
```

The registered no-argument phase command is:

```bash
/home/zelinli6/miniconda3/envs/alphagenome/bin/python \
  -m scripts.run_dpy27_internal_application
```

Before GPU selection, the runner verifies every frozen hash and confirms 15
Model B / paper-loss jobs, 241 ordered tracks, and 83 non-overlapping 131,072-bp
validation-core subwindows per fold. Physical GPUs are restricted to 2 and 3,
with one serial worker queue per GPU. Each fold/seed shard is written
atomically and is reused only when its contract and output hashes match.

## Gene-level quantities and endpoints

Predicted and observed 1-bp signals are accumulated over the union of annotated
exons for each WBcel235 gene. A gene enters the primary analysis only when all
of its unique exon bases lie in the evaluated validation cores. For each
condition, the exon mean is computed before forming

`log1p(G0005 exon mean) - log1p(G0006 exon mean)`.

Within each fold and gene, the G0005 and G0006 predicted exon means are first
averaged across the three seeds. The contrast is then recomputed from those
seed-mean condition signals. Observed condition means must be numerically
identical across seeds. This ordering is important because medians, sign
agreement, log contrasts, and Spearman correlation are nonlinear.

The four prespecified fold-level endpoints are:

1. Predicted median gene contrast on chromosome X minus the autosomal median.
2. Observed median gene contrast on chromosome X minus the autosomal median.
3. Predicted-observed sign concordance among X-linked genes with nonzero contrasts.
4. Predicted-observed Spearman correlation across all complete-coverage genes.

Each endpoint is computed once per fold. The reported estimate is the mean of
the five fold estimates; its 95% interval is a percentile block-bootstrap
interval from 10,000 resamples of the five folds. Individual checkpoint
endpoints are retained as descriptive diagnostics and are not inferential
replicates.

## Output contract

The run root is `runs/v2_p10_dpy27_internal_application_20260813/`:

- `gene_contrasts.tsv`: all 15 checkpoint-by-gene rows.
- `run_endpoints.tsv`: descriptive per-checkpoint endpoints.
- `fold_gene_contrasts.tsv`: condition signals and recomputed contrasts after within-gene three-seed averaging.
- `fold_endpoints.tsv`: the five primary independent block estimates.
- `summary_endpoints.tsv`: fold-primary estimates and bootstrap intervals.
- `figure_source_data.tsv`: exact plotted values, including every gene and fold.
- `figure_dpy27_internal_application.{pdf,svg,tiff,png}`: publication and preview exports.
- `figure_qa.json`: static-source, PDF-font, and raster-pixel QA.
- `shards/seed_*/fold_*`: restartable per-checkpoint tables and audits.

The terminal execution record is
`alphagenome_custom/metadata/v2/p10_dpy27_internal_application_execution.json`.
It reaches `status: completed` only after all 15 shards, aggregate tables,
four figure exports, output hashes, and figure QA pass.
The record also preserves `git_commit`, `hostname`, `working_directory`,
`python_executable`, `controller_invocation`, and `registered_phase_command`
as stable top-level provenance fields.

## Figure legend

**Figure: Model B recovers the DPY-27 response pattern in held-out genomic
blocks.** **a,** Predicted versus observed gene-level log1p contrasts for
DPY-27 RNAi relative to vector RNAi after averaging the two condition signals
across three seeds within each fold and gene. Chromosome-X genes are shown
separately from autosomal genes; the dashed line is identity. **b,** Observed
and predicted chromosome-X-minus-autosome median contrast for each of five
validation folds (connected points). Diamonds show the across-fold mean and
95% fold-block-bootstrap interval. **c,** Fold-level X-linked direction
concordance and predicted-observed gene-level Spearman correlation. Grey points
are the five folds; diamonds and intervals show fold-primary estimates and 95%
bootstrap intervals. Only genes with complete unique-exon coverage are shown.
Both RNA-seq tracks are independently total-signal normalized; the figure is an
internal genomic-block validation, not an independent-condition or causal test.

This methods document does not assert that P10 has run. Execution status and
biological results are established only by the checksum-verified P10 execution
record, terminal R10 review, and generated Source Data artifacts.
