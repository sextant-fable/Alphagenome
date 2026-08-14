# Reviewer 2

**Review status** Frozen individual report. No other reviewer report, concern ledger, or synthesis was consulted.

**Packet SHA-256** `ed44c9d585b9071fba54cd86fded832f00b51c5d221b36bd25b42f480b5a65fa`

**Assessment boundary** This assessment uses only `docs/nature_methods_review_packet_20260814.md`, which summarizes repository evidence verified as of 2026-08-14. The packet is not a manuscript. No abstract, full text, exact title, author-defined claim hierarchy, complete prior-art analysis, figures, or final figure order was supplied. Manuscript-level novelty positioning and nonspecialist readability are therefore not fully assessable. My preassigned emphasis is originality and scientific importance, while all five Nature-style axes are considered.

## Overall assessment

The packet supports a carefully governed and internally informative adaptation study. Model B is evaluated through chromosome-blocked development folds, matched seeds, a size-matched from-scratch comparator, component ablations, and a separately gated final evaluation. The provenance, phase controls, hashes, and explicit distinction between folds, seeds, and tracks are notable strengths.

The present evidence does not yet establish the central case expected for a Nature Methods Article. The internal comparisons show that the full Model B performs better than the supplied A and C configurations, and the P9 factorial analysis usefully localizes much of that advantage to LoRA. They do not establish field-level originality or superiority over available approaches. No complete prior-art analysis or harmonized public-model comparison is supplied. General applicability is also untested outside the same provenance collection, and the only biological application is an internal evaluation on training targets with a sign-discordant chromosome contrast and weak gene-level rank correlation. These limitations constrain both scientific importance and the claim that the method can discover new biology.

## Who would be interested in the results, and why

The work should interest researchers adapting large sequence-to-track models to organisms with smaller data resources, investigators in *C. elegans* functional genomics, and developers of audited computational genomics workflows. The fold-first ablation design and controlled final-test access may also interest groups concerned with leakage and reproducibility in genomic machine learning. Interest beyond these communities will depend on showing that the adaptation provides a distinctive methodological advance, transfers across studies or biological conditions, and enables a biological conclusion that existing approaches could not support.

## Major strengths

- The source inventory and target construction are unusually explicit. The packet records 482 verified RNA-seq runs, 241 biological groups, uniform reprocessing, and on-demand 1-bp and 128-bp loading.
- The P9 design improves interpretability over a simple model leaderboard. It includes a size-matched C model, component configurations, five chromosome folds, three nested algorithmic seeds, and fold-first uncertainty intervals.
- The packet does not treat 241 tracks as independent replicates and does not interpret the worm-by-LoRA interaction as a biological mechanism.
- The P6C evaluation records complete coverage over 83 aligned cores, 10,878,976 eligible bases, six named chromosomes, and all 241 tracks, with immutable lock and report hashes.
- The execution controls preserve commit identifiers, GPU scope, zero locked-test reads during P9 and P10, and the failed first P10 attempt. This is strong internal auditability.
- Limitations are stated candidly. The packet identifies prior chromosome-X exposure, the internal status of P10, the absence of external validation, and the lack of a real variant-scoring run.

## Major Concerns

### R2-M1 [novelty-significance]

**Severity** Major

**Blocking** Yes

**Claim pointer** The proposed AlphaGenome-style adaptation, including worm embedding, LoRA, and a dual-resolution RNA head, is positioned as a novel and significant methods advance.

**Evidence pointer** Packet sections `Shared project description`, `Matched development comparisons`, `P9 component and size-matched evidence`, and `Materials not present in this packet`.

**Concern** The packet provides only within-project comparisons among A, B, C, and component variants. It supplies neither a complete prior-art analysis nor a harmonized comparison with technically applicable public models. The P9 analysis is valuable for attribution, but attribution within one implementation is not evidence of originality relative to the field. It also shows a large LoRA factorial effect of `0.101380`, a much smaller worm-embedding effect of `0.009800`, and a learned-1-bp-head-only contrast whose interval includes zero. Without an explicit account of which components are established practice, which are newly introduced here, and which empirical gains are distinctive, the claimed methods advance cannot be judged.

**Why it matters** Originality and significant improvement over established techniques are central to the proposed article type. The current evidence can support an internally optimized worm model, but it cannot establish that the method changes what the broader community can do.

**Resolution test** Provide a claim-by-claim prior-art analysis and a harmonized benchmark against technically applicable available approaches under the same inputs, splits, targets, metrics, and reporting rules. Separate novelty of the architecture, adaptation strategy, loss, data resource, and governance workflow. Report compute and trainable-parameter differences. If a fair implementation of a comparator is impossible, document the incompatibility and narrow the originality claim rather than treating the comparator as implicitly inferior.

### R2-M2 [experimental-design]

**Severity** Major

**Blocking** Yes

**Claim pointer** The reported performance supports general applicability beyond the data used to develop and select Model B.

**Evidence pointer** Packet sections `Shared project description`, `Final locked evaluation`, `Matched development comparisons`, and `Materials not present in this packet`.

**Concern** Every reported training, development, and final result remains within the 482-run provenance collection. The packet explicitly provides no benchmark isolated by study or laboratory, no unseen biological condition, and no independent species. Leave-one-chromosome-out folds test transfer across genomic partitions, not robustness to laboratory, study, assay, or biological-condition shift. The gated P6C evaluation is procedurally valuable, but it reports only the selected Model B and is not an external comparative benchmark. Chromosome X was also exposed during legacy work, so the final block is not pristine for the project as a whole even though the v2 gate is documented.

**Why it matters** A method intended for broad reuse must demonstrate that its advantage is not specific to one provenance collection, normalization regime, or genomic partition. Without independent validation, the narrow confidence intervals quantify variation across five chromosome folds rather than generalization across studies or biological settings.

**Resolution test** Evaluate the frozen method on a genuinely external, prespecified benchmark held out by study or laboratory, or on an unseen biological condition with a defensible distribution-shift design. Compare all relevant methods under identical processing and report performance by study, condition, and track class. Preserve the existing P6C lock without reopening it, and describe chromosome X as v2-gated but not project-wide pristine.

### R2-M3 [claim-moderation]

**Severity** Major

**Blocking** Yes

**Claim pointer** The P10 DPY-27 analysis supplies a strong application to an important biological question and supports potential discovery of new biology.

**Evidence pointer** Packet section `P10 internal DPY-27 application` and the corresponding limitations in `Materials not present in this packet`.

**Concern** P10 is explicitly an internal application. Both DPY-27 RNAi and vector RNAi were training targets, and the checkpoints were evaluated on their matching development validation cores. The observed chromosome-X-minus-autosome contrast is positive at `0.003065`, whereas the predicted contrast is negative at `-0.002740`. The predicted-observed gene-level Spearman correlation is `0.0846`, despite X-linked direction concordance of `0.6674`. Independent per-track normalization further limits interpretation of absolute condition effects. No orthogonal regulatory, perturbational, or independent RNA validation is supplied. The analysis therefore demonstrates partial directional information under an internal setting, not a convincing biological discovery or unseen-condition prediction.

**Why it matters** The common criteria require strong validation and an application to an important biological question. The current application does not establish that the method captures the central chromosome-level response or yields reliable gene-level biological insight.

**Resolution test** Add a prespecified application on an unseen condition or independent dataset, compare against relevant baselines, and validate the key biological conclusion with orthogonal evidence. Define an effect-size and direction criterion before evaluation. If no stronger application is available, present P10 as an internal diagnostic, state the sign discordance and low rank correlation prominently, and remove claims of biological discovery.

### R2-M4 [reproducibility]

**Severity** Major

**Blocking** No

**Claim pointer** The controlled repository workflow is reproducible and generally usable by independent researchers.

**Evidence pointer** Packet sections `Reproducibility and integrity controls` and `Materials not present in this packet`.

**Concern** The recorded commits, hashes, phase controller, and access checks establish strong internal provenance. They do not yet establish independent reproducibility. The packet provides no public archival release, common-platform installation evidence, container or equivalent environment distribution, or independent user reproduction. Large inputs and checkpoints are referenced only as repository artifacts in the packet.

**Why it matters** Reproducibility and practical reuse are part of the methods contribution. Internal auditability alone does not show that another group can install the software, obtain the required assets, and reproduce a representative result.

**Resolution test** Release a versioned, archived software package with an environment lock or container, installation checks on declared operating systems, data and checkpoint access instructions, checksums, and a small end-to-end example. Record at least one clean reproduction outside the development environment, including expected outputs and tolerances.

### R2-M5 [claim-moderation]

**Severity** Major

**Blocking** No

**Claim pointer** The variant-scoring workflow is a validated downstream capability of the proposed method.

**Evidence pointer** Packet sections `Variant-scoring software` and `Materials not present in this packet`.

**Concern** The available evidence establishes a frozen software specification and CPU synthetic tests. It includes no registered scoring run on a real VCF and checkpoint, no frozen public allele set, and no external eQTL, reporter, CRISPR, or regulatory-allele validation. Software correctness on synthetic cases does not establish empirical validity or biological utility.

**Why it matters** Variant-effect interpretation would materially broaden the paper's significance, but an unexecuted and empirically unvalidated capability cannot support that broader claim.

**Resolution test** Run the frozen workflow on a prespecified real variant set with allele and reference QC, report complete execution and failure accounting, and evaluate predictions against independent functional or genetic evidence. Otherwise describe the CLI as unvalidated software infrastructure and exclude variant-scoring performance or utility from the central claims.

## Minor Comments

### R2-m1 [writing-clarity]

**Severity** Minor

**Affected element** The relationship between the P6B composite score, the P9 component effects, and the two P6C final metrics.

**Evidence pointer** Packet sections `Final locked evaluation`, `Matched development comparisons`, and `P9 component and size-matched evidence`.

**Issue** The packet does not define the composite score or explain how it relates to gene-exon and 128-bp log1p Pearson correlations. A reader cannot reconstruct the selection objective or determine whether development and final reporting emphasize the same estimand.

**Required correction** Define the composite formula, aggregation order, weighting, direction, and treatment of missing values. Explain how it maps to the two final metrics and keep metric names and precision consistent across the abstract, text, figures, and tables.

### R2-m2 [claim-moderation]

**Severity** Minor

**Affected element** Terminology for the P6C evaluation boundary.

**Evidence pointer** Packet sections `Shared project description` and `Final locked evaluation`.

**Issue** The evidence supports a one-time v2-gated evaluation, but chromosome X was exposed in legacy work. Calling the result an untouched or project-wide pristine test would exceed the supplied evidence.

**Required correction** Use wording that distinguishes the v2 lock from the history of the full project. State the prior chromosome-X exposure wherever the independence of the final evaluation is summarized.

## Technical failings that need to be addressed before the case is established

The blocking failings are R2-M1, R2-M2, and R2-M3. The current packet does not establish field-level originality, independent generalization, or a convincing biological application.

## Assessment against Nature-style criteria

**Originality** Not established from the supplied material. The component matrix supports internal attribution, especially a strong LoRA contribution, but no complete prior-art analysis or harmonized available-method comparison is present.

**Scientific importance** Potentially meaningful for organism-specific sequence modeling and *C. elegans* genomics, but outstanding importance is not yet demonstrated. The evidence remains internal to one provenance collection, and the biological application shows weak gene-level agreement and a reversed chromosome-level contrast.

**Interdisciplinary readership** The immediate audience is clear across computational genomics, transfer learning, and worm functional genomics. A conclusion of wider life-science interest depends on external transfer and a biologically consequential application.

**Technical soundness** Internal design and audit controls are strong. Fold-first inference, nested seeds, size matching, immutable locks, and explicit negative evidence improve confidence in the reported internal results. The major technical limitations are the absence of state-of-the-art external comparisons, study-level or condition-level validation, independent biological corroboration, and independent software reproduction.

**Readability for nonspecialists** Not assessable from the supplied material because no manuscript, abstract, figures, or final claim hierarchy was provided. The review packet itself is structured and candid, but it cannot substitute for assessing the paper. The eventual manuscript will need to define the model configurations, composite score, validation hierarchy, and biological interpretation without assuming familiarity with the repository phase system.

## Recommendation posture

Currently not established from the provided evidence as a strong Nature Methods case. I would regard the project as promising and internally disciplined, but substantial new evidence is needed for originality, comparative performance, general applicability, and biological utility. This is a reviewer assessment posture, not an editorial decision.
