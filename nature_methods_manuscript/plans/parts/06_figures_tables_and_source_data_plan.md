# Part 6: Figures, Tables and Source Data Plan

## Design System

No plotting backend has been selected because this is a planning-only stage.
Before implementation, invoke `nature-figure`, resolve the persistent Python/R
preference, write a separate figure contract, then draw every visual with that
one backend. All figures are original project figures; published figures are
structural references only.

Use white backgrounds, 5-7 pt sans-serif labels at final size, editable vector
text and a restrained, stable palette:

- Full Model B: focal blue-green.
- Model A and size-matched Model C: charcoal and neutral grey.
- Public external baselines: a muted contrasting color family.
- Observed signal: black/grey; predicted signal: focal model color.
- Biological alternative/reference allele: two colorblind-safe accents.

## Main Display Budget

Five figures plus one compact table are planned. This uses the Nature Methods
maximum of six display items without forcing every available analysis into the
main text.

### Figure 1: A Traceable C. elegans Resource and Species-Adapted RNA Model

**Core conclusion:** The resource and full model map genomic sequence to
multi-scale worm RNA coverage under a frozen evaluation design.

**Archetype:** Schematic-led composite.

| Panel | Question answered | Planned visual | Data / source-data requirement |
| --- | --- | --- | --- |
| a | What data enter the task? | Provenance flow from RNA-seq runs to biological tracks, grouped by study/stage/tissue/condition at high level. | Counts and group manifest. |
| b | What does the model receive and return? | 131,072-bp DNA window -> frozen trunk + worm embedding + LoRA -> 1-bp and 128-bp RNA outputs. | Model dimensions and parameter table. |
| c | What is externally held out? | Study/lab and unseen-condition split diagram with a hard visual boundary. | External split manifest and leakage audit. |
| d | What does a prediction look like? | Preselected external locus: observed versus predicted coverage, gene/exon annotation and 128-bp summary. | Locus coordinates, track IDs, prediction and observation arrays. |

### Figure 2: The Full Model Wins a Harmonized Study-Isolated Benchmark

**Core conclusion:** Under a single frozen protocol, the full model improves
external RNA prediction over eligible public and from-scratch alternatives.

**Archetype:** Asymmetric quantitative grid; panel a is the hero.

| Panel | Question answered | Planned visual | Statistics |
| --- | --- | --- | --- |
| a | What is the primary comparative effect? | Forest/paired-effect plot of full Model B versus every eligible comparator. | Source- or fold-blocked 95% CI; exact unit shown. |
| b | Is the advantage stable across source and condition? | Faceted study x condition effect plot with a common x-axis. | Prespecified stratified estimates and CI. |
| c | Does the model improve both gene and coverage tasks? | Paired dot/interval comparison for gene/exon and 128-bp endpoints. | Same blocked uncertainty definition as a. |
| d | What is the accuracy-cost position? | Parameter/compute versus primary-score scatter with direct labels. | Total/trainable parameters, wall time, memory and final score. |

### Table 1: Frozen External Benchmark and Comparator Contract

**Core conclusion:** The comparison is technically interpretable and fair.

Rows: external datasets; columns: source isolation, unseen condition, reference
build, input context, target/bin, trainable data, training budget, total and
trainable parameters, tuning rule, metrics, eligibility status and reason for
any exclusion. Keep it compact; move all hyperparameters to Supplementary
Table 2.

### Figure 3: Controlled Attribution Identifies the Measured Value of Adaptation

**Core conclusion:** LoRA supplies the dominant measured gain, worm embedding
adds a smaller gain, and the current primary endpoint does not establish a
necessity claim for the learned 1-bp head.

**Archetype:** Quantitative grid, adapted from the audited P9 figure only after
terminology and style harmonization.

| Panel | Question answered | Planned visual | Statistics |
| --- | --- | --- | --- |
| a | How do complete and removal/replacement models compare? | Five-effect paired forest, full Model B as reference. | Five fold means, seeds nested within fold, 10,000 hierarchical bootstrap. |
| b | What are the worm-embedding, LoRA and interaction effects? | 2 x 2 factorial effect panel. | Paired factorial estimates and CI; interaction labeled non-additive. |
| c | Is the result broad across RNA tracks? | Track-wise paired-effect distribution with individual track points and density. | Descriptive only; do not label tracks as independent n. |

### Figure 4: Frozen Sequence-Difference Scores Recover an Independent Regulatory Endpoint

**Core conclusion:** A ref/alt or haplotype scoring procedure ranks and signs a
prespecified independent same-species cis-eQTL signal.

**Archetype:** Asymmetric mixed-modality figure; panel b is the biological hero.

| Panel | Question answered | Planned visual | Statistics |
| --- | --- | --- | --- |
| a | What was frozen before labels were inspected? | Variant scoring workflow: VCF/reference check -> paired inference -> gene score -> null match -> endpoint. | Manifest and rule IDs. |
| b | Are effects directionally correct? | Signed observed versus predicted effect scatter with quadrants and a neutral diagonal. | Direction concordance and blocked CI. |
| c | Can scores distinguish credible cis variants from matched nulls? | PR and ROC curves with per-block uncertainty. | AUROC/AUPRC and CI. |
| d | What does one prespecified regulatory locus show? | Coverage, ref/alt difference and ISM/motif annotation. | Exact locus and raw prediction arrays. |

### Figure 5: External Stratification Defines the Model's Useful Domain

**Core conclusion:** The external model advantage has a measurable domain of
reliability rather than an unqualified global average.

**Archetype:** Quantitative grid.

| Panel | Question answered | Planned visual | Statistics |
| --- | --- | --- | --- |
| a | Which source/condition strata retain the gain? | Ordered forest of predefined strata. | Same primary blocked unit and CI. |
| b | How does performance vary with signal/annotation features? | Binned trend with raw source-level points. | Predefined bins and no post hoc thresholding. |
| c | Are predicted scores calibrated enough for prioritization? | Calibration or risk-coverage curve. | Frozen calibration measure. |
| d | What failures matter to users? | Three representative, preselected external failure/success examples with identical axes. | Selection rule and raw arrays. |

## Extended Data Plan

| Item | Purpose |
| --- | --- |
| ED Fig. 1 | Complete P6B A/B/C x loss matrix, fold/seed paired effects, track distribution and strata. |
| ED Fig. 2 | Architecture dimensions, LoRA target modules and total/trainable parameter accounting. |
| ED Fig. 3 | External benchmark leakage, reference/coordinate and window-context overlap audit. |
| ED Fig. 4 | Full external performance strata and all predeclared metrics. |
| ED Fig. 5 | Additional external locus coverage examples and error modes. |
| ED Fig. 6 | DPY-27 internal diagnostic, including the observed/predicted sign disagreement and its normalized-scale boundary. |
| ED Fig. 7 | Variant score sensitivity, null matching and indel/haplotype handling. |
| ED Fig. 8 | Clean-room reproduction, environment and release verification. |

## Source Data Standard

Each main figure receives a tabular source-data package with:

- one row per raw analytical unit;
- immutable input and output manifest identifiers;
- plotted values before rounding;
- group labels, exclusions and exact selection rules;
- a README defining units, transformations and statistic;
- hash of the final rendered figure and script revision.

No source-data table is allowed to hide a selected subset without the complete
population table and a deterministic selection rule.
