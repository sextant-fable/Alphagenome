# Reviewer 1

- **Packet SHA-256** `ed44c9d585b9071fba54cd86fded832f00b51c5d221b36bd25b42f480b5a65fa`
- **Assessment boundary** This review is based only on the immutable packet dated 2026-08-14. The packet is a repository evidence summary rather than a manuscript. No abstract, author-defined central claim, prior-art analysis, complete Methods, figures, code release, or supplementary material was supplied. Novelty and manuscript readability are therefore not fully assessable. This report evaluates technical readiness for the prospective Nature Methods framing described in the packet and does not infer evidence beyond it.
- **Emphasis brief** Technical soundness and technical failings, with all five Nature-style axes also assessed.
- **Overall assessment** The project shows unusually careful internal governance. Its audited provenance, controller-gated workflow, immutable test records, explicit inference unit, matched parameter control, and component ablations provide a credible internal comparison of the proposed Model B against the supplied A and C configurations. The present packet does not yet establish a complete methods-paper case. The locked evaluation remains within the collected resource, comparison to available public approaches is absent, the DPY-27 application is internal and gives weak or directionally discordant results, and independent software reproduction has not been demonstrated. The precision of several inferential claims is also difficult to judge with five chromosome folds. These limitations do not erase the internal results, but they constrain them to a development study rather than evidence of broad, externally validated methodological utility.
- **Who would be interested in the results, and why** Researchers in regulatory genomics, sequence-to-function modelling, model adaptation across species, and *C. elegans* functional genomics would be interested in the controlled transfer of a large sequence model to 241 worm RNA-seq tracks. Method developers would also value the matched ablations because they distinguish the contributions of LoRA, the worm embedding, and the learned 1-bp head within the supplied design.
- **Major strengths** The 485-accession provenance audit and uniform reduction of 482 RNA-seq runs to 241 biological groups provide a traceable data foundation. The one-time locked evaluation records complete coverage, fixed hashes, and zero later test-signal reads. P9 improves substantially on the original unmatched comparison by adding a size-matched C model and targeted component ablations across five folds and three seeds. The packet correctly treats folds as the inference unit, seeds as nested algorithmic repeats, tracks as descriptive outcomes, and the factorial interaction as statistical non-additivity rather than mechanism. Preservation of the failed P10 rendering attempt and its corrected rerun further supports procedural integrity.
- **Major Concerns**

## R1-M1 [experimental-design]

**Severity** Major  
**Blocking** Yes  
**Axis** Technical soundness and scientific importance  
**Claim pointer** The prospective methods case uses the locked performance across 83 aligned cores, six chromosomes, and 241 tracks as validation of the adapted Model B.  
**Evidence pointer** Packet sections `Shared project description`, `Final locked evaluation`, and `Materials not present in this packet`  
**Concern** The supplied validation remains within the 482-run provenance collection. No study-held-out or laboratory-held-out benchmark, unseen biological condition, or independent species is present. Chromosome X was exposed during legacy work, so it is explicitly not pristine for the project as a whole. The one-time gate and hash controls establish disciplined access within v2, but they do not establish distributional independence or external generalization. The packet also does not provide the genomic split definitions needed to verify that every test core and its receptive-field context were excluded from training and model selection.  
**Why it matters** General applicability is part of the stated methods-paper standard. Performance on protected cores from the same resource can support internal interpolation, but it cannot by itself support broad transfer, robustness across studies, or an unbiased project-level chromosome-X claim.  
**Resolution test** Add a prospectively frozen benchmark isolated by study or laboratory and, preferably, by an unseen biological condition. Report source-level independence, window and receptive-field overlap audits, all preprocessing decisions, and performance stratified by study and chromosome. Treat chromosome X as historically exposed or exclude it from claims of pristine testing. If independent validation cannot be added, narrow the central claim to internal development performance.

## R1-M2 [novelty-significance]

**Severity** Major  
**Blocking** Yes  
**Axis** Originality, scientific importance, and technical soundness  
**Claim pointer** The P6B and P9 comparisons are used to support a significant methodological advance for sequence-to-RNA-track prediction.  
**Evidence pointer** Packet sections `Matched development comparisons`, `P9 component and size-matched evidence`, and `Materials not present in this packet`  
**Concern** The evidence compares only repository-defined A, B, C, and component variants. No Enformer, Borzoi, other technically applicable public model, or strong conventional baseline has been evaluated under a harmonized benchmark. The size-matched C comparison controls trainable parameter count, but the packet does not establish matching of total parameters, pretraining information, optimization budget, or inference cost. Internal ablations identify which supplied components matter, but they do not establish performance relative to available approaches. The absence of a complete prior-art analysis also prevents an independent novelty assessment.  
**Why it matters** A methods paper needs both a credible advance over available approaches and a clear account of what is new. Without external or public-method comparators, the reported gains show superiority only within the authors' design space.  
**Resolution test** Run technically applicable public and conventional baselines on the same frozen inputs, splits, targets, and metrics. Report tuning allowances, pretraining access, total and trainable parameters, compute, and uncertainty. Give a component-level prior-art analysis that distinguishes the worm embedding, LoRA adaptation, and dual-resolution head from published alternatives. Where a comparator cannot be run, provide a technically specific exclusion rationale and moderate the comparative claim.

## R1-M3 [statistical-rigor]

**Severity** Major  
**Blocking** No  
**Axis** Technical soundness  
**Claim pointer** The P9 fold-first bootstrap intervals quantify component effects, while the P6C means summarize final performance across tracks.  
**Evidence pointer** Packet sections `Final locked evaluation` and `P9 component and size-matched evidence`  
**Concern** Treating seeds as nested repeats and folds as the primary unit is appropriate, but the inferential sample remains five chromosome folds. Ten thousand bootstrap draws do not increase the number of independent genomic units. The packet gives no fold-level effects, leave-one-fold-out sensitivity, or uncertainty for the two locked final-test means. This matters especially for small effects such as the LoRA-only and learned-1-bp-head-only comparisons. The formula and weighting for the composite score are also absent from the packet.  
**Why it matters** Narrow intervals generated from five clusters may convey more precision than the design can support, and aggregate means can conceal chromosome-specific or track-specific failure. The direction of the large internal effects may be credible while their precision and generality remain uncertain.  
**Resolution test** Show every fold-by-seed result and the per-track distributions, define the composite score and aggregation hierarchy, and provide leave-one-chromosome-out sensitivity. Use an interval or randomization analysis whose assumptions are explicit for five clustered units, and qualify precision accordingly. For P6C, report chromosome-stratified and track-stratified uncertainty together with sensitivity to low-signal or low-variance tracks.

## R1-M4 [experimental-design]

**Severity** Major  
**Blocking** Yes  
**Axis** Technical soundness, scientific importance, and interdisciplinary interest  
**Claim pointer** P10 is presented as the internal DPY-27 biological application for the proposed methods paper.  
**Evidence pointer** Packet section `P10 internal DPY-27 application` and the independent-validation items in `Materials not present in this packet`  
**Concern** P10 reuses checkpoints on matching development validation cores, and both DPY-27 RNAi and vector RNAi were training targets. The predicted chromosome-X-minus-autosome contrast has the opposite sign to the observed contrast. Predicted-observed gene-level Spearman correlation is `0.0846`, despite X-linked direction concordance of `0.6674`. Each condition was independently normalized to total signal of approximately `1e8`, which makes the contrast compositional and leaves its relation to global expression changes unresolved. No orthogonal RNA measurement, perturbation assay, reporter, or independent condition validates the prediction. P10 therefore functions as an informative internal diagnostic, not yet as a successful demonstration of biological discovery or out-of-condition utility.  
**Why it matters** The stated Nature Methods scope expects application to an important biological question and potential for discovering new biology. The current application does not establish that capability and includes a directionally discordant headline contrast.  
**Resolution test** Evaluate a prospectively held-out biological condition or independent study with prespecified endpoints. Include normalization sensitivity on an appropriate non-compositional scale where possible, simple biological and sequence-only baselines, and orthogonal validation of the central predictions. If these tests are unavailable, present P10 as a limitation or failure analysis and provide a different validated biological application.

## R1-M5 [reproducibility]

**Severity** Major  
**Blocking** Yes  
**Axis** Technical soundness and interdisciplinary interest  
**Claim pointer** The controlled repository workflow and frozen records are intended to support reproducibility and reuse of the method.  
**Evidence pointer** Packet sections `Reproducibility and integrity controls` and `Materials not present in this packet`  
**Concern** Commit identifiers, hashes, phase controls, and preserved failures provide strong internal auditability. They are not yet evidence that an independent user can install the software, obtain compatible data and weights, recreate preprocessing, train or load the selected model, and reproduce the reported metrics. The packet contains no public archival release, common-operating-system installation test, container, or independent reproduction record. A full manuscript and complete Methods were also not supplied, so reporting completeness cannot be assessed.  
**Why it matters** Reproducibility and practical accessibility are core parts of the stated methods-paper case. A repository that is reproducible only inside its originating environment has limited community impact and cannot yet support claims of reusable methodology.  
**Resolution test** Produce a versioned public archive with immutable code, environment and container specifications, data and weight access instructions, licenses, hardware requirements, deterministic evaluation commands, expected outputs, and checksums. Demonstrate installation and end-to-end metric reproduction in a clean environment by an independent user or automated runner, and report any irreducible nondeterminism.

## R1-M6 [data-resource-quality]

**Severity** Major  
**Blocking** No  
**Axis** Technical soundness  
**Claim pointer** The variant-scoring CLI is part of the method's functional capability.  
**Evidence pointer** Packet section `Variant-scoring software` and the variant-validation items in `Materials not present in this packet`  
**Concern** The CLI has useful input checks, strand ensembling, substitution deltas, indel diagnostics, optional ISM, and synthetic CPU tests. No real VCF plus checkpoint run has been registered, and no external eQTL, regulatory-allele, reporter, CRISPR, or independent RNA result validates its scores. Synthetic tests establish software behaviour on constructed cases, not biological calibration, allele ranking, or indel validity.  
**Why it matters** If variant scoring is presented as an empirical capability rather than prospective software, the current evidence does not establish accuracy or utility on real variants. This gap is secondary if variant scoring is removed from the main claim.  
**Resolution test** Run the frozen CLI on a predefined public VCF or haplotype set with the selected checkpoint. Validate score direction and ranking against an independent molecular-QTL or perturbation benchmark, include appropriate null and baseline methods, and report calibration by variant class. Otherwise label the feature as unvalidated software and exclude performance or biological-utility claims.

- **Minor Comments**

## R1-m1 [writing-clarity]

**Severity** Minor  
**Axis** Readability for nonspecialists  
**Affected element** Definitions and aggregation of the reported performance quantities  
**Evidence pointer** Packet sections `Final locked evaluation` and `Matched development comparisons`  
**Issue** Terms including `aligned core`, `eligible base`, `coverage fraction`, `composite score`, and `gene-exon log1p Pearson correlation` are not defined in the packet. A reader cannot tell how regions and tracks are weighted or which exclusions precede averaging.  
**Required correction** Define each estimand, transformation, inclusion rule, denominator, weighting choice, and aggregation order in the manuscript. Include a compact schematic linking 1-bp and 128-bp targets to the two evaluation metrics.

## R1-m2 [writing-clarity]

**Severity** Minor  
**Axis** Readability for nonspecialists  
**Affected element** Relationship among development folds, the final-test cores, and chromosome X  
**Evidence pointer** Packet sections `Shared project description` and `Final locked evaluation`  
**Issue** The text states that development uses leave-one-chromosome-out folds on chromosomes I through V, while the one-time evaluation covers chromosomes I through V and X. Without a split diagram, the distinction between held-out fold regions, final-test cores, and historical exposure of chromosome X is difficult to follow.  
**Required correction** Add a single split schematic and timeline that names the chromosomes, genomic units, access gates, model-selection steps, and historical chromosome-X exposure. Use distinct terms for fold validation and final locked evaluation throughout.

- **Technical failings that need to be addressed before the case is established** `R1-M1`, `R1-M2`, `R1-M4`, and `R1-M5` are blocking under the prospective Nature Methods framing. They concern independent validation, comparison with available approaches, a biologically convincing application, and independent reproducibility. `R1-M3` and `R1-M6` remain major limitations that require resolution or explicit narrowing of the relevant claims.
- **Assessment against Nature-style criteria**

  **Originality** Not assessable from the provided material. The packet does not contain a complete prior-art analysis or an author-defined novelty claim. The matched component evidence is useful, but component attribution is not equivalent to novelty relative to published methods.

  **Scientific importance** Potentially substantial for cross-species regulatory modelling and worm genomics. The current evidence establishes internal improvements, not yet outstanding or broadly generalizable importance. The weak internal DPY-27 result limits the biological-impact case.

  **Interdisciplinary readership** The work could interest computational biologists, functional genomics researchers, and model-transfer researchers. Wider interest depends on showing that the method generalizes beyond the originating collection and produces a validated biological insight or broadly useful predictive capability.

  **Technical soundness** Internal governance and matched ablations are strong. The case remains incomplete because external validation, public-method comparisons, a convincing biological application, robust uncertainty characterization, and independent reproduction are absent or insufficient in the supplied evidence.

  **Readability for nonspecialists** Not assessable for the manuscript because no manuscript, abstract, or figures were supplied. The packet is organized and appropriately caveated, but central split and metric terminology would need a schematic and explicit definitions for nonspecialists.

- **Recommendation posture** Currently not established from the provided evidence. I would support reassessment after the blocking validation, comparator, biological-application, and reproducibility requirements are met, with the statistical and variant-scoring claims either strengthened or narrowed. This is a technical readiness judgment within the supplied boundary, not an editorial decision.
