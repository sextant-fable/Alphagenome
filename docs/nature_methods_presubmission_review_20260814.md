# Nature Methods pre-submission review

## Review setup

- **Input scope** Prospective Nature Methods Article on an AlphaGenome-style sequence-to-RNA-track adaptation for *C. elegans*, assessed only from the immutable repository evidence packet dated 2026-08-14.
- **Assessment boundary** This is a bounded pre-submission review of a repository evidence summary, not a full manuscript review. No full manuscript, abstract, exact title, author-defined claim hierarchy, complete prior-art analysis, or final figure sequence was supplied. Manuscript-level novelty positioning, prose readability, figure accessibility, and abstract-to-results consistency are therefore not fully assessable.
- **Shared manuscript claim summary** The prospective paper would present a controlled AlphaGenome-style adaptation to 241 *C. elegans* RNA-seq tracks, with chromosome-blocked development, a separately gated final evaluation, matched model comparisons and component ablations, an internal DPY-27 application, and a variant-scoring software extension.
- **Visible evidence base** The immutable packet includes the locked P6C evaluation, P6B and P9 development comparisons, the P10 internal DPY-27 application, synthetic tests for variant-scoring software, and selected reproducibility and integrity records.
- **Missing materials affecting confidence** The packet contains no study- or laboratory-isolated external benchmark, unseen-condition or independent-species validation, harmonized public-model comparison, external molecular or variant validation, real variant-scoring run, public archival software release, independent user reproduction, complete prior-art analysis, or full manuscript materials.
- **Isolation and freeze status** The three reports were produced in separate mutually blind contexts from the same immutable packet and were frozen before this synthesis pass. The synthesis was generated only after all three frozen files were available and was not shown to the reviewers.
- **Reproduction note** Each frozen report is reproduced fully and faithfully below. Only Markdown heading levels were adjusted to keep the report under its reviewer section. No reviewer prose, concern ID, severity, blocking flag, evidence pointer, or resolution test was rewritten.
- **Immutable packet SHA-256** `ed44c9d585b9071fba54cd86fded832f00b51c5d221b36bd25b42f480b5a65fa`
- **Reviewer 1 report SHA-256** `f36bcd943a34d8a8d1bfd4ce113c975ba4db2f7e0b3b9f83467745ba42502840`
- **Reviewer 2 report SHA-256** `2213544e20ad7df5b71c019e07e5b14870997cd7c373f3a7976a1d98cf9b14da`
- **Reviewer 3 report SHA-256** `a0a34a22fa639e2cd0088cbaf0a239768ac9d9f693db7ae990e3a6d46a571e69`

## Reviewer 1

- **Packet SHA-256** `ed44c9d585b9071fba54cd86fded832f00b51c5d221b36bd25b42f480b5a65fa`
- **Assessment boundary** This review is based only on the immutable packet dated 2026-08-14. The packet is a repository evidence summary rather than a manuscript. No abstract, author-defined central claim, prior-art analysis, complete Methods, figures, code release, or supplementary material was supplied. Novelty and manuscript readability are therefore not fully assessable. This report evaluates technical readiness for the prospective Nature Methods framing described in the packet and does not infer evidence beyond it.
- **Emphasis brief** Technical soundness and technical failings, with all five Nature-style axes also assessed.
- **Overall assessment** The project shows unusually careful internal governance. Its audited provenance, controller-gated workflow, immutable test records, explicit inference unit, matched parameter control, and component ablations provide a credible internal comparison of the proposed Model B against the supplied A and C configurations. The present packet does not yet establish a complete methods-paper case. The locked evaluation remains within the collected resource, comparison to available public approaches is absent, the DPY-27 application is internal and gives weak or directionally discordant results, and independent software reproduction has not been demonstrated. The precision of several inferential claims is also difficult to judge with five chromosome folds. These limitations do not erase the internal results, but they constrain them to a development study rather than evidence of broad, externally validated methodological utility.
- **Who would be interested in the results, and why** Researchers in regulatory genomics, sequence-to-function modelling, model adaptation across species, and *C. elegans* functional genomics would be interested in the controlled transfer of a large sequence model to 241 worm RNA-seq tracks. Method developers would also value the matched ablations because they distinguish the contributions of LoRA, the worm embedding, and the learned 1-bp head within the supplied design.
- **Major strengths** The 485-accession provenance audit and uniform reduction of 482 RNA-seq runs to 241 biological groups provide a traceable data foundation. The one-time locked evaluation records complete coverage, fixed hashes, and zero later test-signal reads. P9 improves substantially on the original unmatched comparison by adding a size-matched C model and targeted component ablations across five folds and three seeds. The packet correctly treats folds as the inference unit, seeds as nested algorithmic repeats, tracks as descriptive outcomes, and the factorial interaction as statistical non-additivity rather than mechanism. Preservation of the failed P10 rendering attempt and its corrected rerun further supports procedural integrity.
- **Major Concerns**

### R1-M1 [experimental-design]

**Severity** Major  
**Blocking** Yes  
**Axis** Technical soundness and scientific importance  
**Claim pointer** The prospective methods case uses the locked performance across 83 aligned cores, six chromosomes, and 241 tracks as validation of the adapted Model B.  
**Evidence pointer** Packet sections `Shared project description`, `Final locked evaluation`, and `Materials not present in this packet`  
**Concern** The supplied validation remains within the 482-run provenance collection. No study-held-out or laboratory-held-out benchmark, unseen biological condition, or independent species is present. Chromosome X was exposed during legacy work, so it is explicitly not pristine for the project as a whole. The one-time gate and hash controls establish disciplined access within v2, but they do not establish distributional independence or external generalization. The packet also does not provide the genomic split definitions needed to verify that every test core and its receptive-field context were excluded from training and model selection.  
**Why it matters** General applicability is part of the stated methods-paper standard. Performance on protected cores from the same resource can support internal interpolation, but it cannot by itself support broad transfer, robustness across studies, or an unbiased project-level chromosome-X claim.  
**Resolution test** Add a prospectively frozen benchmark isolated by study or laboratory and, preferably, by an unseen biological condition. Report source-level independence, window and receptive-field overlap audits, all preprocessing decisions, and performance stratified by study and chromosome. Treat chromosome X as historically exposed or exclude it from claims of pristine testing. If independent validation cannot be added, narrow the central claim to internal development performance.

### R1-M2 [novelty-significance]

**Severity** Major  
**Blocking** Yes  
**Axis** Originality, scientific importance, and technical soundness  
**Claim pointer** The P6B and P9 comparisons are used to support a significant methodological advance for sequence-to-RNA-track prediction.  
**Evidence pointer** Packet sections `Matched development comparisons`, `P9 component and size-matched evidence`, and `Materials not present in this packet`  
**Concern** The evidence compares only repository-defined A, B, C, and component variants. No Enformer, Borzoi, other technically applicable public model, or strong conventional baseline has been evaluated under a harmonized benchmark. The size-matched C comparison controls trainable parameter count, but the packet does not establish matching of total parameters, pretraining information, optimization budget, or inference cost. Internal ablations identify which supplied components matter, but they do not establish performance relative to available approaches. The absence of a complete prior-art analysis also prevents an independent novelty assessment.  
**Why it matters** A methods paper needs both a credible advance over available approaches and a clear account of what is new. Without external or public-method comparators, the reported gains show superiority only within the authors' design space.  
**Resolution test** Run technically applicable public and conventional baselines on the same frozen inputs, splits, targets, and metrics. Report tuning allowances, pretraining access, total and trainable parameters, compute, and uncertainty. Give a component-level prior-art analysis that distinguishes the worm embedding, LoRA adaptation, and dual-resolution head from published alternatives. Where a comparator cannot be run, provide a technically specific exclusion rationale and moderate the comparative claim.

### R1-M3 [statistical-rigor]

**Severity** Major  
**Blocking** No  
**Axis** Technical soundness  
**Claim pointer** The P9 fold-first bootstrap intervals quantify component effects, while the P6C means summarize final performance across tracks.  
**Evidence pointer** Packet sections `Final locked evaluation` and `P9 component and size-matched evidence`  
**Concern** Treating seeds as nested repeats and folds as the primary unit is appropriate, but the inferential sample remains five chromosome folds. Ten thousand bootstrap draws do not increase the number of independent genomic units. The packet gives no fold-level effects, leave-one-fold-out sensitivity, or uncertainty for the two locked final-test means. This matters especially for small effects such as the LoRA-only and learned-1-bp-head-only comparisons. The formula and weighting for the composite score are also absent from the packet.  
**Why it matters** Narrow intervals generated from five clusters may convey more precision than the design can support, and aggregate means can conceal chromosome-specific or track-specific failure. The direction of the large internal effects may be credible while their precision and generality remain uncertain.  
**Resolution test** Show every fold-by-seed result and the per-track distributions, define the composite score and aggregation hierarchy, and provide leave-one-chromosome-out sensitivity. Use an interval or randomization analysis whose assumptions are explicit for five clustered units, and qualify precision accordingly. For P6C, report chromosome-stratified and track-stratified uncertainty together with sensitivity to low-signal or low-variance tracks.

### R1-M4 [experimental-design]

**Severity** Major  
**Blocking** Yes  
**Axis** Technical soundness, scientific importance, and interdisciplinary interest  
**Claim pointer** P10 is presented as the internal DPY-27 biological application for the proposed methods paper.  
**Evidence pointer** Packet section `P10 internal DPY-27 application` and the independent-validation items in `Materials not present in this packet`  
**Concern** P10 reuses checkpoints on matching development validation cores, and both DPY-27 RNAi and vector RNAi were training targets. The predicted chromosome-X-minus-autosome contrast has the opposite sign to the observed contrast. Predicted-observed gene-level Spearman correlation is `0.0846`, despite X-linked direction concordance of `0.6674`. Each condition was independently normalized to total signal of approximately `1e8`, which makes the contrast compositional and leaves its relation to global expression changes unresolved. No orthogonal RNA measurement, perturbation assay, reporter, or independent condition validates the prediction. P10 therefore functions as an informative internal diagnostic, not yet as a successful demonstration of biological discovery or out-of-condition utility.  
**Why it matters** The stated Nature Methods scope expects application to an important biological question and potential for discovering new biology. The current application does not establish that capability and includes a directionally discordant headline contrast.  
**Resolution test** Evaluate a prospectively held-out biological condition or independent study with prespecified endpoints. Include normalization sensitivity on an appropriate non-compositional scale where possible, simple biological and sequence-only baselines, and orthogonal validation of the central predictions. If these tests are unavailable, present P10 as a limitation or failure analysis and provide a different validated biological application.

### R1-M5 [reproducibility]

**Severity** Major  
**Blocking** Yes  
**Axis** Technical soundness and interdisciplinary interest  
**Claim pointer** The controlled repository workflow and frozen records are intended to support reproducibility and reuse of the method.  
**Evidence pointer** Packet sections `Reproducibility and integrity controls` and `Materials not present in this packet`  
**Concern** Commit identifiers, hashes, phase controls, and preserved failures provide strong internal auditability. They are not yet evidence that an independent user can install the software, obtain compatible data and weights, recreate preprocessing, train or load the selected model, and reproduce the reported metrics. The packet contains no public archival release, common-operating-system installation test, container, or independent reproduction record. A full manuscript and complete Methods were also not supplied, so reporting completeness cannot be assessed.  
**Why it matters** Reproducibility and practical accessibility are core parts of the stated methods-paper case. A repository that is reproducible only inside its originating environment has limited community impact and cannot yet support claims of reusable methodology.  
**Resolution test** Produce a versioned public archive with immutable code, environment and container specifications, data and weight access instructions, licenses, hardware requirements, deterministic evaluation commands, expected outputs, and checksums. Demonstrate installation and end-to-end metric reproduction in a clean environment by an independent user or automated runner, and report any irreducible nondeterminism.

### R1-M6 [data-resource-quality]

**Severity** Major  
**Blocking** No  
**Axis** Technical soundness  
**Claim pointer** The variant-scoring CLI is part of the method's functional capability.  
**Evidence pointer** Packet section `Variant-scoring software` and the variant-validation items in `Materials not present in this packet`  
**Concern** The CLI has useful input checks, strand ensembling, substitution deltas, indel diagnostics, optional ISM, and synthetic CPU tests. No real VCF plus checkpoint run has been registered, and no external eQTL, regulatory-allele, reporter, CRISPR, or independent RNA result validates its scores. Synthetic tests establish software behaviour on constructed cases, not biological calibration, allele ranking, or indel validity.  
**Why it matters** If variant scoring is presented as an empirical capability rather than prospective software, the current evidence does not establish accuracy or utility on real variants. This gap is secondary if variant scoring is removed from the main claim.  
**Resolution test** Run the frozen CLI on a predefined public VCF or haplotype set with the selected checkpoint. Validate score direction and ranking against an independent molecular-QTL or perturbation benchmark, include appropriate null and baseline methods, and report calibration by variant class. Otherwise label the feature as unvalidated software and exclude performance or biological-utility claims.

- **Minor Comments**

### R1-m1 [writing-clarity]

**Severity** Minor  
**Axis** Readability for nonspecialists  
**Affected element** Definitions and aggregation of the reported performance quantities  
**Evidence pointer** Packet sections `Final locked evaluation` and `Matched development comparisons`  
**Issue** Terms including `aligned core`, `eligible base`, `coverage fraction`, `composite score`, and `gene-exon log1p Pearson correlation` are not defined in the packet. A reader cannot tell how regions and tracks are weighted or which exclusions precede averaging.  
**Required correction** Define each estimand, transformation, inclusion rule, denominator, weighting choice, and aggregation order in the manuscript. Include a compact schematic linking 1-bp and 128-bp targets to the two evaluation metrics.

### R1-m2 [writing-clarity]

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

## Reviewer 2

**Review status** Frozen individual report. No other reviewer report, concern ledger, or synthesis was consulted.

**Packet SHA-256** `ed44c9d585b9071fba54cd86fded832f00b51c5d221b36bd25b42f480b5a65fa`

**Assessment boundary** This assessment uses only `docs/nature_methods_review_packet_20260814.md`, which summarizes repository evidence verified as of 2026-08-14. The packet is not a manuscript. No abstract, full text, exact title, author-defined claim hierarchy, complete prior-art analysis, figures, or final figure order was supplied. Manuscript-level novelty positioning and nonspecialist readability are therefore not fully assessable. My preassigned emphasis is originality and scientific importance, while all five Nature-style axes are considered.

### Overall assessment

The packet supports a carefully governed and internally informative adaptation study. Model B is evaluated through chromosome-blocked development folds, matched seeds, a size-matched from-scratch comparator, component ablations, and a separately gated final evaluation. The provenance, phase controls, hashes, and explicit distinction between folds, seeds, and tracks are notable strengths.

The present evidence does not yet establish the central case expected for a Nature Methods Article. The internal comparisons show that the full Model B performs better than the supplied A and C configurations, and the P9 factorial analysis usefully localizes much of that advantage to LoRA. They do not establish field-level originality or superiority over available approaches. No complete prior-art analysis or harmonized public-model comparison is supplied. General applicability is also untested outside the same provenance collection, and the only biological application is an internal evaluation on training targets with a sign-discordant chromosome contrast and weak gene-level rank correlation. These limitations constrain both scientific importance and the claim that the method can discover new biology.

### Who would be interested in the results, and why

The work should interest researchers adapting large sequence-to-track models to organisms with smaller data resources, investigators in *C. elegans* functional genomics, and developers of audited computational genomics workflows. The fold-first ablation design and controlled final-test access may also interest groups concerned with leakage and reproducibility in genomic machine learning. Interest beyond these communities will depend on showing that the adaptation provides a distinctive methodological advance, transfers across studies or biological conditions, and enables a biological conclusion that existing approaches could not support.

### Major strengths

- The source inventory and target construction are unusually explicit. The packet records 482 verified RNA-seq runs, 241 biological groups, uniform reprocessing, and on-demand 1-bp and 128-bp loading.
- The P9 design improves interpretability over a simple model leaderboard. It includes a size-matched C model, component configurations, five chromosome folds, three nested algorithmic seeds, and fold-first uncertainty intervals.
- The packet does not treat 241 tracks as independent replicates and does not interpret the worm-by-LoRA interaction as a biological mechanism.
- The P6C evaluation records complete coverage over 83 aligned cores, 10,878,976 eligible bases, six named chromosomes, and all 241 tracks, with immutable lock and report hashes.
- The execution controls preserve commit identifiers, GPU scope, zero locked-test reads during P9 and P10, and the failed first P10 attempt. This is strong internal auditability.
- Limitations are stated candidly. The packet identifies prior chromosome-X exposure, the internal status of P10, the absence of external validation, and the lack of a real variant-scoring run.

### Major Concerns

#### R2-M1 [novelty-significance]

**Severity** Major

**Blocking** Yes

**Claim pointer** The proposed AlphaGenome-style adaptation, including worm embedding, LoRA, and a dual-resolution RNA head, is positioned as a novel and significant methods advance.

**Evidence pointer** Packet sections `Shared project description`, `Matched development comparisons`, `P9 component and size-matched evidence`, and `Materials not present in this packet`.

**Concern** The packet provides only within-project comparisons among A, B, C, and component variants. It supplies neither a complete prior-art analysis nor a harmonized comparison with technically applicable public models. The P9 analysis is valuable for attribution, but attribution within one implementation is not evidence of originality relative to the field. It also shows a large LoRA factorial effect of `0.101380`, a much smaller worm-embedding effect of `0.009800`, and a learned-1-bp-head-only contrast whose interval includes zero. Without an explicit account of which components are established practice, which are newly introduced here, and which empirical gains are distinctive, the claimed methods advance cannot be judged.

**Why it matters** Originality and significant improvement over established techniques are central to the proposed article type. The current evidence can support an internally optimized worm model, but it cannot establish that the method changes what the broader community can do.

**Resolution test** Provide a claim-by-claim prior-art analysis and a harmonized benchmark against technically applicable available approaches under the same inputs, splits, targets, metrics, and reporting rules. Separate novelty of the architecture, adaptation strategy, loss, data resource, and governance workflow. Report compute and trainable-parameter differences. If a fair implementation of a comparator is impossible, document the incompatibility and narrow the originality claim rather than treating the comparator as implicitly inferior.

#### R2-M2 [experimental-design]

**Severity** Major

**Blocking** Yes

**Claim pointer** The reported performance supports general applicability beyond the data used to develop and select Model B.

**Evidence pointer** Packet sections `Shared project description`, `Final locked evaluation`, `Matched development comparisons`, and `Materials not present in this packet`.

**Concern** Every reported training, development, and final result remains within the 482-run provenance collection. The packet explicitly provides no benchmark isolated by study or laboratory, no unseen biological condition, and no independent species. Leave-one-chromosome-out folds test transfer across genomic partitions, not robustness to laboratory, study, assay, or biological-condition shift. The gated P6C evaluation is procedurally valuable, but it reports only the selected Model B and is not an external comparative benchmark. Chromosome X was also exposed during legacy work, so the final block is not pristine for the project as a whole even though the v2 gate is documented.

**Why it matters** A method intended for broad reuse must demonstrate that its advantage is not specific to one provenance collection, normalization regime, or genomic partition. Without independent validation, the narrow confidence intervals quantify variation across five chromosome folds rather than generalization across studies or biological settings.

**Resolution test** Evaluate the frozen method on a genuinely external, prespecified benchmark held out by study or laboratory, or on an unseen biological condition with a defensible distribution-shift design. Compare all relevant methods under identical processing and report performance by study, condition, and track class. Preserve the existing P6C lock without reopening it, and describe chromosome X as v2-gated but not project-wide pristine.

#### R2-M3 [claim-moderation]

**Severity** Major

**Blocking** Yes

**Claim pointer** The P10 DPY-27 analysis supplies a strong application to an important biological question and supports potential discovery of new biology.

**Evidence pointer** Packet section `P10 internal DPY-27 application` and the corresponding limitations in `Materials not present in this packet`.

**Concern** P10 is explicitly an internal application. Both DPY-27 RNAi and vector RNAi were training targets, and the checkpoints were evaluated on their matching development validation cores. The observed chromosome-X-minus-autosome contrast is positive at `0.003065`, whereas the predicted contrast is negative at `-0.002740`. The predicted-observed gene-level Spearman correlation is `0.0846`, despite X-linked direction concordance of `0.6674`. Independent per-track normalization further limits interpretation of absolute condition effects. No orthogonal regulatory, perturbational, or independent RNA validation is supplied. The analysis therefore demonstrates partial directional information under an internal setting, not a convincing biological discovery or unseen-condition prediction.

**Why it matters** The common criteria require strong validation and an application to an important biological question. The current application does not establish that the method captures the central chromosome-level response or yields reliable gene-level biological insight.

**Resolution test** Add a prespecified application on an unseen condition or independent dataset, compare against relevant baselines, and validate the key biological conclusion with orthogonal evidence. Define an effect-size and direction criterion before evaluation. If no stronger application is available, present P10 as an internal diagnostic, state the sign discordance and low rank correlation prominently, and remove claims of biological discovery.

#### R2-M4 [reproducibility]

**Severity** Major

**Blocking** No

**Claim pointer** The controlled repository workflow is reproducible and generally usable by independent researchers.

**Evidence pointer** Packet sections `Reproducibility and integrity controls` and `Materials not present in this packet`.

**Concern** The recorded commits, hashes, phase controller, and access checks establish strong internal provenance. They do not yet establish independent reproducibility. The packet provides no public archival release, common-platform installation evidence, container or equivalent environment distribution, or independent user reproduction. Large inputs and checkpoints are referenced only as repository artifacts in the packet.

**Why it matters** Reproducibility and practical reuse are part of the methods contribution. Internal auditability alone does not show that another group can install the software, obtain the required assets, and reproduce a representative result.

**Resolution test** Release a versioned, archived software package with an environment lock or container, installation checks on declared operating systems, data and checkpoint access instructions, checksums, and a small end-to-end example. Record at least one clean reproduction outside the development environment, including expected outputs and tolerances.

#### R2-M5 [claim-moderation]

**Severity** Major

**Blocking** No

**Claim pointer** The variant-scoring workflow is a validated downstream capability of the proposed method.

**Evidence pointer** Packet sections `Variant-scoring software` and `Materials not present in this packet`.

**Concern** The available evidence establishes a frozen software specification and CPU synthetic tests. It includes no registered scoring run on a real VCF and checkpoint, no frozen public allele set, and no external eQTL, reporter, CRISPR, or regulatory-allele validation. Software correctness on synthetic cases does not establish empirical validity or biological utility.

**Why it matters** Variant-effect interpretation would materially broaden the paper's significance, but an unexecuted and empirically unvalidated capability cannot support that broader claim.

**Resolution test** Run the frozen workflow on a prespecified real variant set with allele and reference QC, report complete execution and failure accounting, and evaluate predictions against independent functional or genetic evidence. Otherwise describe the CLI as unvalidated software infrastructure and exclude variant-scoring performance or utility from the central claims.

### Minor Comments

#### R2-m1 [writing-clarity]

**Severity** Minor

**Affected element** The relationship between the P6B composite score, the P9 component effects, and the two P6C final metrics.

**Evidence pointer** Packet sections `Final locked evaluation`, `Matched development comparisons`, and `P9 component and size-matched evidence`.

**Issue** The packet does not define the composite score or explain how it relates to gene-exon and 128-bp log1p Pearson correlations. A reader cannot reconstruct the selection objective or determine whether development and final reporting emphasize the same estimand.

**Required correction** Define the composite formula, aggregation order, weighting, direction, and treatment of missing values. Explain how it maps to the two final metrics and keep metric names and precision consistent across the abstract, text, figures, and tables.

#### R2-m2 [claim-moderation]

**Severity** Minor

**Affected element** Terminology for the P6C evaluation boundary.

**Evidence pointer** Packet sections `Shared project description` and `Final locked evaluation`.

**Issue** The evidence supports a one-time v2-gated evaluation, but chromosome X was exposed in legacy work. Calling the result an untouched or project-wide pristine test would exceed the supplied evidence.

**Required correction** Use wording that distinguishes the v2 lock from the history of the full project. State the prior chromosome-X exposure wherever the independence of the final evaluation is summarized.

### Technical failings that need to be addressed before the case is established

The blocking failings are R2-M1, R2-M2, and R2-M3. The current packet does not establish field-level originality, independent generalization, or a convincing biological application.

### Assessment against Nature-style criteria

**Originality** Not established from the supplied material. The component matrix supports internal attribution, especially a strong LoRA contribution, but no complete prior-art analysis or harmonized available-method comparison is present.

**Scientific importance** Potentially meaningful for organism-specific sequence modeling and *C. elegans* genomics, but outstanding importance is not yet demonstrated. The evidence remains internal to one provenance collection, and the biological application shows weak gene-level agreement and a reversed chromosome-level contrast.

**Interdisciplinary readership** The immediate audience is clear across computational genomics, transfer learning, and worm functional genomics. A conclusion of wider life-science interest depends on external transfer and a biologically consequential application.

**Technical soundness** Internal design and audit controls are strong. Fold-first inference, nested seeds, size matching, immutable locks, and explicit negative evidence improve confidence in the reported internal results. The major technical limitations are the absence of state-of-the-art external comparisons, study-level or condition-level validation, independent biological corroboration, and independent software reproduction.

**Readability for nonspecialists** Not assessable from the supplied material because no manuscript, abstract, figures, or final claim hierarchy was provided. The review packet itself is structured and candid, but it cannot substitute for assessing the paper. The eventual manuscript will need to define the model configurations, composite score, validation hierarchy, and biological interpretation without assuming familiarity with the repository phase system.

### Recommendation posture

Currently not established from the provided evidence as a strong Nature Methods case. I would regard the project as promising and internally disciplined, but substantial new evidence is needed for originality, comparative performance, general applicability, and biological utility. This is a reviewer assessment posture, not an editorial decision.

## Reviewer 3

### Frozen report review setup

- **Input scope** The immutable common review packet describes a prospective Nature Methods Article on an AlphaGenome-style sequence-to-RNA-track adaptation for *C. elegans*. It reports repository state and verified results as of 2026-08-14.
- **Packet integrity** SHA-256 `ed44c9d585b9071fba54cd86fded832f00b51c5d221b36bd25b42f480b5a65fa`.
- **Assessment boundary** This is a bounded assessment of the packet, not a full manuscript review. No manuscript, abstract, exact title, author-defined claim hierarchy, complete prior-art analysis, or final figure sequence was supplied. Direct assessment of prose readability, figure accessibility, abstract-to-results consistency, and component-level novelty is therefore limited. No evidence outside the packet was considered.
- **Visible evidence base** The packet includes the locked P6C evaluation, the P6B and P9 internal comparisons and component analyses, the P10 internal DPY-27 application, synthetic tests for variant-scoring software, and selected reproducibility and integrity records.
- **Missing materials affecting confidence** The packet contains no study- or laboratory-isolated external benchmark, unseen-condition or independent-species validation, harmonized comparison with technically applicable public models, external molecular or variant validation, real variant-scoring run, public archival software release, or independent user reproduction.
- **Report status** Frozen after independent review and groundedness checks. This report contains no information from another reviewer and no cross-review synthesis.

- **Overall assessment** Reviewer 3 places greatest weight on interdisciplinary readership and readability for nonspecialists while assessing all five common axes. The project shows unusually careful internal auditability, a disciplined chromosome-fold design, and informative matched component comparisons. The size-matched from-scratch control and the explicit separation of genomic folds, algorithmic seeds, and descriptive track outcomes are particular strengths. The current evidence nevertheless does not establish the central case expected for a broadly influential methods paper. The only supplied biological application uses training-target conditions and fails to recover the sign of the reported chromosome-level contrast, while no harmonized comparison with available public approaches establishes the claimed methods advance. The work could interest communities beyond *C. elegans* if it demonstrates reliable transfer of a sequence model to a distinct organismal RNA-seq setting and enables new biological inference. The packet currently demonstrates an internally controlled engineering and evaluation framework more clearly than that broader conclusion.
- **Who would be interested in the results, and why** Immediate interest would come from researchers developing sequence-to-track models, studying *C. elegans* functional genomics, adapting large genomic models across organisms, or building reproducible computational biology workflows. A wider methods readership could care if the work shows that model adaptation delivers reliable gains over available approaches, generalizes beyond the development collection, and supports a biological conclusion that could not be reached otherwise. Those wider-interest conditions are not yet demonstrated in the supplied evidence.
- **Major strengths** The provenance inventory and uniform grouping provide a clearly counted data basis. The one-time evaluation records full stated coverage across 83 aligned cores, six chromosomes, and 241 tracks. P9 includes a size-matched from-scratch model and factorial component analysis, with folds treated as the inference unit and seeds nested within folds. The packet also discloses important limitations rather than obscuring them. These include prior exposure of chromosome X, the internal status of P10, the descriptive status of track outcomes, the non-mechanistic meaning of the interaction term, and the absence of real variant-scoring validation. Commit hashes, locked-test access controls, and preservation of a failed P10 attempt strengthen the integrity trail.

### Major Concerns

#### R3-M1 [experimental-design]

**Severity** Major

**Blocking** Yes

**Claim pointer** The prospective Article would use P10 to establish an application to an important biological question and support the method's capacity to recover biologically meaningful regulation.

**Evidence pointer** Packet section `P10 internal DPY-27 application`, including the condition status and four reported endpoints. Packet section `Materials not present in this packet`, including the absent unseen-condition and independent biological validations.

**Concern** P10 is an internal evaluation on matching development validation cores. Both DPY-27 RNAi and vector RNAi were training targets. More importantly, the observed chromosome-X-minus-autosome contrast is positive at `0.003065`, whereas the predicted contrast is negative at `-0.002740`. The predicted-observed gene-level Spearman correlation is only `0.0846`. Direction concordance of `0.6674` does not resolve the contradictory chromosome-level sign or the weak gene-level ranking. On the evidence shown, P10 is an informative failure analysis or internal diagnostic, not strong validation of biological fidelity or discovery.

**Why it matters** The supplied journal criteria require strong validation and an application to an important biological question. The only visible application does not establish that the method reproduces the principal reported biological pattern or generalizes to a new condition. Presenting it as affirmative biological validation would overstate the evidence and weaken both technical soundness and scientific importance.

**Resolution test** Demonstrate a prespecified application on data isolated from training by biological condition and preferably by study or laboratory. Show recovery of the direction, magnitude, and gene-level structure of biologically relevant endpoints, with orthogonal RNA, reporter, perturbational, or equivalent independent validation appropriate to the claim. Otherwise, describe P10 explicitly as an internal limitation analysis and remove or narrow claims of biological discovery.

#### R3-M2 [novelty-significance]

**Severity** Major

**Blocking** Yes

**Claim pointer** The AlphaGenome-style adaptation is a novel method or a significant improvement over established sequence-to-track approaches for biological data.

**Evidence pointer** Packet sections `Matched development comparisons` and `P9 component and size-matched evidence`. Packet section `Materials not present in this packet`, including the absent harmonized public-model comparison and complete prior-art analysis.

**Concern** The supplied comparisons establish differences among repository configurations A, B, and C and isolate several components of Model B. They do not establish performance relative to technically applicable available approaches. The packet explicitly reports no completed Enformer, Borzoi, or other public-model comparison under a harmonized benchmark, and no prior-art analysis sufficient to determine which components are novel. The size-matched C comparison addresses parameter count within the project but cannot substitute for comparison with existing methods or for novelty positioning.

**Why it matters** Internal ablation can explain where an in-repository gain originates, but it cannot show that the method advances the state of the art or delivers the significant improvement expected of a methods Article. Without this distinction, originality and broad scientific importance remain unassessable.

**Resolution test** Define the claimed advance component by component, provide a complete and accurate prior-art comparison, and evaluate technically applicable public methods and simple baselines under the same input data, splits, target definitions, metrics, and leakage controls. Report accuracy together with trainable parameters, compute, and adaptation requirements. If a public model cannot be evaluated fairly, document the incompatibility and narrow the comparative claim.

#### R3-M3 [claim-moderation]

**Severity** Major

**Blocking** No

**Claim pointer** The method has general applicability beyond the particular genomic regions and 241-track collection used for development.

**Evidence pointer** Packet sections `Shared project description` and `Final locked evaluation`. Packet section `Materials not present in this packet`, including the absent study- or laboratory-isolated benchmark, unseen biological condition, and independent species validation.

**Concern** The chromosome-fold design tests performance across held-out genomic regions within the same provenance collection and fixed track set. The packet provides no evidence for transfer to an unseen study, laboratory, biological condition, assay collection, or species. Chromosome X also was exposed during legacy work and is not pristine for the project as a whole. The reported final evaluation therefore supports a bounded within-collection genomic generalization claim, not a broad claim of experimental or cross-organism applicability.

**Why it matters** Readers outside the immediate project need to know what kind of transfer the method actually achieves. Conflating genomic-region holdout with transfer across biological or technical distributions would inflate the method's expected usefulness and obscure the conditions under which another laboratory could rely on it.

**Resolution test** Add an evaluation isolated by study or laboratory and by a biologically meaningful unseen condition, with an explicit leakage audit and prespecified endpoints. Cross-species evidence would strengthen a broader adaptation claim but is not required if the Article consistently limits generality to the demonstrated *C. elegans* setting. Qualify chromosome X as project-exposed wherever final-test independence is described.

#### R3-M4 [statistical-rigor]

**Severity** Major

**Blocking** No

**Claim pointer** The fold-first bootstrap intervals provide reliable uncertainty for the P6B and P9 comparative and component effects.

**Evidence pointer** Packet section `Matched development comparisons`, including the post hoc submission analysis. Packet section `P9 component and size-matched evidence`, including the five-fold inference unit and 10,000-replicate nested bootstrap.

**Concern** The packet correctly identifies only five genomic folds as independent inference units, but it reports aggregate bootstrap intervals without the five fold-level paired effects or sensitivity to any single chromosome. Increasing bootstrap replicates does not increase the number of independent folds. The adequacy and stability of the intervals are therefore difficult to assess, particularly for effects near zero such as LoRA-only and learned-1-bp-head-only.

**Why it matters** Component attribution is central to explaining why Model B works. If small estimated effects or their intervals are driven by one chromosome, the architectural interpretation and its portability may change even when the overall Model B advantage remains.

**Resolution test** Show all fold-level paired effects and seed dispersion, report leave-one-fold-out sensitivity, and state the inferential population supported by five named chromosomes. Provide a justification or sensitivity analysis for interval construction with five independent units. Moderate claims for small component effects when their direction is not stable across folds.

#### R3-M5 [reproducibility]

**Severity** Major

**Blocking** No

**Claim pointer** The method is reproducible and usable by communities beyond the originating repository.

**Evidence pointer** Packet section `Reproducibility and integrity controls`. Packet section `Materials not present in this packet`, including the absent common-system installation evidence, container, public archival release, and independent user reproduction.

**Concern** Controller records, commit hashes, output hashes, and test-lock checks provide strong internal provenance. They do not yet establish that another group can obtain the software, install it, reconstruct the environment, run a small example, or reproduce a principal result. No versioned archival release, containerized environment, common-system installation check, or independent reproduction is supplied.

**Why it matters** Practical adoption is part of the significance case for a technical paper. Internal auditability protects integrity, while external reproducibility determines whether the method can affect a community of researchers.

**Resolution test** Release a versioned and archived code package with environment or container specifications, model and data access instructions, a tractable end-to-end example, expected outputs and checksums, and documented resource requirements. Include an independent clean-environment reproduction of at least one principal evaluation or application result.

#### R3-M6 [claim-moderation]

**Severity** Major

**Blocking** No

**Claim pointer** The variant-scoring functionality can support scientific inference about regulatory variants.

**Evidence pointer** Packet section `Variant-scoring software`. Packet section `Materials not present in this packet`, including the absent external variant evidence and real frozen-VCF execution.

**Concern** The variant-scoring implementation has synthetic CPU tests but no registered real VCF plus checkpoint run and no evaluation against cis-eQTLs, regulatory alleles, reporters, CRISPR results, or another independent molecular readout. This establishes software behavior under synthetic tests, not empirical validity of variant-effect scores.

**Why it matters** Variant interpretation would substantially expand the paper's biological and interdisciplinary appeal, but an unvalidated capability can mislead readers about the present method's scientific reach. Its inclusion as a headline result would outpace the evidence.

**Resolution test** Evaluate a frozen public variant or haplotype set with prespecified metrics, matched baselines, and independent measured effects. Until then, identify variant scoring as an unvalidated software extension or omit it from the principal claims.

### Minor Comments

#### R3-m1 [writing-clarity]

**Severity** Minor

**Affected element** The naming and organization of models, losses, phases, gates, and evaluation units.

**Evidence pointer** Packet sections `Shared project description`, `Matched development comparisons`, `P9 component and size-matched evidence`, and `P10 internal DPY-27 application`.

**Issue** Labels such as P6B, P6C, P9, P10, R9, R10, A, B, C, G0005, G0006, `paper loss`, `aligned cores`, and `composite score` are efficient audit identifiers but do not communicate the scientific design to a nonspecialist. Their density obscures which intervention, comparator, target, and outcome belong to each result.

**Required correction** Lead each section and figure with descriptive scientific names. Define each metric and unit at first use, state what a higher value means, and retain phase or gate identifiers only as secondary reproducibility labels.

#### R3-m2 [writing-clarity]

**Severity** Minor

**Affected element** Explanation of the dual-resolution prediction task and headline correlations.

**Evidence pointer** Packet sections `Shared project description` and `Final locked evaluation`.

**Issue** The packet reports 1-bp and 128-bp targets and two mean per-track correlations without explaining their distinct biological purpose, the construction of gene-exon versus 128-bp outcomes, or what performance at the reported magnitude enables in practice.

**Required correction** Add a concise conceptual explanation of both output resolutions, their intended biological uses, the denominators and averaging scheme for each correlation, and an interpretable baseline or reference point. A nonspecialist should be able to connect each metric to a concrete scientific task.

#### R3-m3 [claim-moderation]

**Severity** Minor

**Affected element** Description of the locked final evaluation and chromosome X.

**Evidence pointer** Packet sections `Shared project description` and `Final locked evaluation`.

**Issue** The active workflow uses a separately gated one-time evaluation, but chromosome X was exposed during legacy work and is not pristine for the project as a whole. An unqualified phrase such as `unseen test chromosome` or `pristine held-out test` would therefore be inaccurate.

**Required correction** Distinguish the active v2 access lock from project-wide historical exposure every time test independence is summarized. State precisely what information was unavailable to v2 model selection and what had been exposed previously.

#### R3-m4 [figures-and-tables]

**Severity** Minor

**Affected element** The final figure sequence and visual explanation of the method.

**Evidence pointer** Packet section `Packet scope`, which states that no final figure order was supplied. Packet sections `Shared project description` and `Verified evidence base`, which contain the elements that require visual linkage.

**Issue** Figure accessibility is not assessable from the packet. The prospective paper nevertheless needs a visual route from sequence inputs and dual-resolution targets through chromosome-based development, locked evaluation, component tests, and biological application. Without that route, the phase-oriented evidence structure is likely to remain inaccessible to readers outside the immediate modeling workflow.

**Required correction** Include a simple conceptual schematic and a figure sequence organized around scientific questions rather than project phases. Show paired fold-level comparisons and the agreement or disagreement between observed and predicted biological endpoints, while defining uncertainty and replicate units in each legend.

- **Technical failings that need to be addressed before the case is established** R3-M1 and R3-M2. The present packet lacks a convincing independent biological application and a harmonized comparison that establishes novelty and performance relative to available approaches.
- **Assessment against Nature-style criteria** **Originality** Not assessable from the supplied material because component-level prior art and external method comparisons are absent. The internal component analysis is informative but does not establish novelty. **Scientific importance** Potentially substantial if the adaptation enables reliable RNA-track prediction and new biology in *C. elegans*. The supplied P10 result does not yet demonstrate that impact. **Interdisciplinary readership** Immediate interest is credible across computational genomics, model-organism biology, and sequence-model adaptation. Broader interest remains conditional on external performance, transferable validation, and a biologically informative application. **Technical soundness** Internal controls, audit records, matched analyses, and explicit inference units are strengths. R3-M1 and R3-M2 prevent the present evidence from establishing the central methods case, while R3-M3 through R3-M6 limit generality, inferential confidence, reuse, and application scope. **Readability for nonspecialists** Direct assessment is not possible without the manuscript, abstract, and figures. The packet's reliance on internal phase codes and unexplained metrics indicates that a substantial conceptual rewrite will be needed if this vocabulary carries into the Article.
- **Recommendation posture** Promising but the broad-interest and technical case is currently not established from the provided evidence. I would be supportive of further consideration after the blocking biological-validation and comparative-benchmark issues are resolved and the central claim is made explicit. Final journal fit and prose quality remain outside this bounded assessment.

## Cross-review synthesis (post-review; not shown to reviewers)

本节只对三份已冻结报告进行事后归并。只有至少两份报告独立提出同一底层问题时，才标记为共识。单份报告中的实质性要求保留为侧重点差异。这里不代表编辑决定，也不判断最终期刊去向。

### Consensus strengths

- 三份报告一致认可内部审计与研究治理。共同依据包括可追溯的 485 个 accession 清单、482 个 RNA-seq run 到 241 个 biological group 的统一处理、P6C 的一次性锁定与哈希记录、P9 和 P10 的零 locked-test signal read，以及失败 P10 尝试的保留。对应各报告的 `Major strengths`。
- 三份报告一致认可开发比较的设计透明度。共同强调五个 chromosome fold 是推断单位，三个 seed 是 fold 内的算法重复，241 个 track 只作描述性结果，size-matched C 和组件消融提高了内部归因能力，worm-by-LoRA interaction 只被解释为统计非加和项。对应各报告的 `Major strengths`。
- 三份报告一致认为现有证据能够支持一个审计严格、内部比较信息充分的开发研究。它们也一致肯定证据包主动披露 chromosome X 的历史暴露、P10 的内部性质和 variant scoring 尚未实证验证等限制。该优势不能替代外部有效性、方法新颖性或生物学效用证据。

### Consensus blocking concerns

1. **缺少足以建立领域级新颖性和相对性能的先验技术定位与同台比较**。来源为 `R1-M2`、`R2-M1` 和 `R3-M2`，三者均为 `Blocking Yes`。证据指针为 packet 的 `Shared project description`、`Matched development comparisons`、`P9 component and size-matched evidence` 和 `Materials not present in this packet`。A、B、C 及组件变体之间的比较能够回答项目内部归因问题，但不能单独证明相对于可用方法的原创性或显著改进。当前还没有完整 prior-art analysis，也没有在统一输入、split、target 和 metric 下完成 technically applicable public model 或强 conventional baseline 的比较。
2. **缺少对开发来源集合之外泛化能力的独立验证**。来源为 `R1-M1`、`R2-M2` 和 `R3-M3`。`R1-M1` 与 `R2-M2` 为 `Blocking Yes`，`R3-M3` 为 `Blocking No`，因此底层问题有三方共识，而阻断等级存在一份较宽松判断。证据指针为 packet 的 `Shared project description`、`Final locked evaluation`、`Matched development comparisons` 和 `Materials not present in this packet`。现有结果均位于同一 482-run provenance collection。leave-one-chromosome-out 设计检验的是集合内 genomic partition 转移，不等同于跨 study、laboratory、condition、assay collection 或 species 的分布外泛化。chromosome X 在 legacy 工作中曾被暴露，也不能被表述为 project-wide pristine。
3. **P10 目前不能作为肯定性的生物学应用或新生物学发现证据**。来源为 `R1-M4`、`R2-M3` 和 `R3-M1`，三者均为 `Blocking Yes`。证据指针为 packet 的 `P10 internal DPY-27 application` 和 `Materials not present in this packet`。DPY-27 RNAi 与 vector RNAi 都是训练 target，评估使用 matching development validation cores。observed chromosome-X-minus-autosome contrast 为 `0.003065`，predicted contrast 为 `-0.002740`，方向相反。gene-level Spearman correlation 为 `0.0846`。`0.6674` 的 X-linked direction concordance 不能消除上述矛盾，也没有独立或正交验证支持生物学发现主张。按当前证据，P10 更适合被界定为内部诊断或限制分析。

### Other consensus major concerns

1. **五个独立 chromosome fold 对区间精度和组件归因稳定性的约束**。来源为 `R1-M3` 和 `R3-M4`，两者均为 `Blocking No`。证据指针为 packet 的 `Final locked evaluation`、`Matched development comparisons` 和 `P9 component and size-matched evidence`。10,000 次 bootstrap 不会增加独立 genomic unit 数量。packet 未展示五个 fold 的 paired effect、leave-one-fold-out sensitivity，也没有给出两个 locked final-test mean 的不确定性。小效应尤其需要谨慎解释。
2. **内部可审计性尚未转化为独立可复现性**。来源为 `R1-M5`、`R2-M4` 和 `R3-M5`。三者均为 Major，但 `R1-M5` 是 `Blocking Yes`，`R2-M4` 与 `R3-M5` 是 `Blocking No`。证据指针为 packet 的 `Reproducibility and integrity controls` 和 `Materials not present in this packet`。commit、hash、controller 和 access check 证明了强内部 provenance，但 packet 未提供公共归档版本、可移植环境、常见系统安装验证、外部用户复现或干净环境端到端复现。
3. **variant-scoring 只建立了软件行为，尚未建立经验有效性**。来源为 `R1-M6`、`R2-M5` 和 `R3-M6`，三者均为 `Blocking No`。证据指针为 packet 的 `Variant-scoring software` 和 `Materials not present in this packet`。现有 synthetic CPU tests 不能替代 real VCF plus checkpoint execution，也不能支持 allele ranking、calibration 或 biological utility 主张。该功能若不补充实证验证，应明确作为 unvalidated software extension，并移出核心性能或生物学结论。

### Where emphasis differs across reviewers

- Reviewer 1 的技术有效性侧重点更严格。`R1-M1` 额外要求核验 test core 及其 receptive-field context 与训练和选择过程之间的重叠，`R1-M4` 还明确指出每个条件独立归一化到约 `1e8` 会使 P10 contrast 具有 compositional 限制。这些细节没有被其他报告以同等明确程度单列，但应保留为技术审计重点。
- Reviewer 2 更强调方法贡献必须逐组件定位。`R2-M1` 用 P9 数值区分 LoRA、worm embedding 和 learned-1-bp-head-only 的贡献，并要求分别界定 architecture、adaptation strategy、loss、data resource 和 governance workflow 的新颖性。它还特别强调 P6C 只报告选定 Model B，因此本身不是外部 comparative benchmark。
- Reviewer 3 更强调跨学科可读性和主张边界。`R3-M3` 明确指出，如果文章始终将 generality 限定在已证明的 *C. elegans* setting，cross-species evidence 并非该窄主张的必要条件。`R3-m1` 至 `R3-m4` 则集中指出 phase code、metric、dual-resolution task 和 figure sequence 对非专业读者构成的理解障碍。
- 阻断等级并非所有问题都一致。外部泛化被 Reviewer 1 和 Reviewer 2 判为 blocking，而 Reviewer 3 判为 nonblocking。独立可复现性只被 Reviewer 1 判为 blocking。综合结果保留这些差异，不把多数意见改写成三方一致的阻断判断。

### Minor revision checklist

1. 定义 `aligned core`、`eligible base`、`coverage fraction`、`composite score`、`gene-exon log1p Pearson correlation`、1-bp 和 128-bp 输出。写清 transformation、denominator、inclusion rule、weighting、aggregation order、missing value 处理及 score direction。来源为 `R1-m1`、`R2-m1`、`R3-m2`。证据指针为 `Shared project description`、`Final locked evaluation`、`Matched development comparisons` 和 `P9 component and size-matched evidence`。
2. 用一个 split schematic 和 timeline 区分 fold validation、final-test core、model selection、access gate 与 chromosome X 的 legacy exposure。所有 test-independence 表述都必须区分 v2 one-time lock 与 project-wide historical exposure。来源为 `R1-m2`、`R2-m2`、`R3-m3`。证据指针为 `Shared project description` 和 `Final locked evaluation`。
3. 科学叙事优先使用描述性名称解释 intervention、comparator、target 和 outcome。P6B、P6C、P9、P10、R9、R10、A、B、C、G0005 和 G0006 只作为次级审计标签。来源为 `R3-m1`。证据指针为 `Shared project description`、`Matched development comparisons`、`P9 component and size-matched evidence` 和 `P10 internal DPY-27 application`。
4. 建立以科学问题而不是项目 phase 为主线的 conceptual schematic 和 final figure sequence。图中需要呈现 dual-resolution task、paired fold-level comparison、uncertainty unit，以及 P10 observed 与 predicted endpoint 的一致和不一致部分。来源为 `R3-m4`。证据指针为 `Packet scope`、`Shared project description` 和 `Verified evidence base`。

### Broad-interest / significance readout

三份报告都识别出明确的直接受众，包括 sequence-to-track model 开发者、跨物种或低资源物种模型适配研究者、*C. elegans* functional genomics 社区，以及重视 genomic machine-learning auditability 的研究者。现有材料已经支持这些群体对内部方法设计和审计框架的兴趣。

更广泛的方法学和生命科学意义仍然是有条件的。它依赖于对 available approaches 的公平比较、开发来源集合之外的验证、一个能够正确恢复关键生物学 endpoint 并有独立支持的应用，以及外部用户可执行的复现路径。由于没有 manuscript、abstract、完整 prior-art analysis 和 figures，field-level originality、非专业读者可读性和最终 broad-readership presentation 均不能从当前 packet 直接判定。这里不能推出确定的编辑结论。

### Most important issues to resolve before a strong Nature Methods case is established

1. 用完整 prior-art analysis 和统一 benchmark 界定相对于 available public and conventional approaches 的具体方法增量。对应 `R1-M2`、`R2-M1`、`R3-M2`。
2. 增加按 study 或 laboratory 隔离、最好同时具有 unseen biological condition 的预先冻结验证，并完成 source、window 和 receptive-field overlap audit。否则把主张严格收窄为同一 provenance collection 内的 genomic generalization。对应 `R1-M1`、`R2-M2`、`R3-M3`。
3. 将 P10 如实定位为内部诊断或限制分析，除非能够提供训练条件之外、endpoint 预先规定并有独立或正交证据支持的生物学应用。对应 `R1-M4`、`R2-M3`、`R3-M1`。
4. 把内部 provenance 转化为外部可复现性。需要可归档版本、环境或 container、data 和 checkpoint access、预期输出与 checksum，以及干净环境或独立用户的代表性复现。对应 `R1-M5`、`R2-M4`、`R3-M5`。其 blocking 等级在报告之间不同。
5. 展示 fold-by-seed 结果、fold-level paired effect、leave-one-fold-out sensitivity、composite score 定义和适合五个独立 fold 的不确定性说明，并对小组件效应限制精度主张。对应 `R1-M3`、`R3-M4`。
6. 对 variant-scoring 完成预先定义的 real-data execution 和独立经验验证，或将其明确限制为尚未验证的软件扩展并移出 headline claim。对应 `R1-M6`、`R2-M5`、`R3-M6`。

## Risk / unsupported claims

- **Field-level novelty or state-of-the-art superiority** 当前不受支持。packet 只有项目内部模型和组件比较，没有完整 prior-art analysis 或 harmonized public-model benchmark。来源为 `R1-M2`、`R2-M1`、`R3-M2`。
- **Broad general applicability** 当前不受支持。已验证的是同一 provenance collection 和固定 241-track 集合中的 genomic partition transfer，不是跨 study、laboratory、condition、assay collection 或 species 的泛化。来源为 `R1-M1`、`R2-M2`、`R3-M3`。
- **P10 as biological discovery or affirmative biological validation** 当前不受支持。training-target conditions、相反符号的 chromosome-level contrast、`0.0846` 的 gene-level Spearman correlation 和缺少独立验证共同限制该主张。来源为 `R1-M4`、`R2-M3`、`R3-M1`。
- **Project-wide pristine chromosome-X test** 不受支持。只能陈述 v2 的 one-time gate 与锁定记录，同时必须披露 chromosome X 的 legacy exposure。来源为 `R1-M1`、`R2-M2`、`R3-M3`、`R1-m2`、`R2-m2`、`R3-m3`。
- **Independent software reproducibility and general usability** 尚未建立。内部 controller、commit 和 hash 记录不能代替公共归档、可移植安装和独立复现。来源为 `R1-M5`、`R2-M4`、`R3-M5`。
- **Empirically validated variant-effect prediction** 不受支持。只有 synthetic CPU tests，没有 registered real VCF plus checkpoint run 或外部 molecular validation。来源为 `R1-M6`、`R2-M5`、`R3-M6`。
- **High-precision inference from five folds** 需要限制。bootstrap replicate 数量不会增加独立 fold 数，且缺少 fold-level effect 和 sensitivity display。来源为 `R1-M3`、`R3-M4`。
- **Verified absence of genomic or receptive-field leakage** 从当前 packet 无法完整评估。`R1-M1` 指出没有提供足以核验 test core 与 receptive-field context 排除情况的 split definitions。这是未能评估的边界，不是已证明存在 leakage。
- **Manuscript-level readability, figure accessibility, claim hierarchy, and abstract-to-results consistency** 无法从当前材料评估，因为 manuscript、abstract、final figures 和 final claim hierarchy 均未提供。来源为三份报告的 assessment boundary 和 `R3-m4`。
- 本综合只基于 packet 与三份冻结报告，不引入外部文献、额外项目证据、未提供的实验或编辑判断。
