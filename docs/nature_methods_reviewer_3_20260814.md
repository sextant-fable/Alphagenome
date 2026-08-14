# Review setup

- **Input scope** The immutable common review packet describes a prospective Nature Methods Article on an AlphaGenome-style sequence-to-RNA-track adaptation for *C. elegans*. It reports repository state and verified results as of 2026-08-14.
- **Packet integrity** SHA-256 `ed44c9d585b9071fba54cd86fded832f00b51c5d221b36bd25b42f480b5a65fa`.
- **Assessment boundary** This is a bounded assessment of the packet, not a full manuscript review. No manuscript, abstract, exact title, author-defined claim hierarchy, complete prior-art analysis, or final figure sequence was supplied. Direct assessment of prose readability, figure accessibility, abstract-to-results consistency, and component-level novelty is therefore limited. No evidence outside the packet was considered.
- **Visible evidence base** The packet includes the locked P6C evaluation, the P6B and P9 internal comparisons and component analyses, the P10 internal DPY-27 application, synthetic tests for variant-scoring software, and selected reproducibility and integrity records.
- **Missing materials affecting confidence** The packet contains no study- or laboratory-isolated external benchmark, unseen-condition or independent-species validation, harmonized comparison with technically applicable public models, external molecular or variant validation, real variant-scoring run, public archival software release, or independent user reproduction.
- **Report status** Frozen after independent review and groundedness checks. This report contains no information from another reviewer and no cross-review synthesis.

# Reviewer 3

- **Overall assessment** Reviewer 3 places greatest weight on interdisciplinary readership and readability for nonspecialists while assessing all five common axes. The project shows unusually careful internal auditability, a disciplined chromosome-fold design, and informative matched component comparisons. The size-matched from-scratch control and the explicit separation of genomic folds, algorithmic seeds, and descriptive track outcomes are particular strengths. The current evidence nevertheless does not establish the central case expected for a broadly influential methods paper. The only supplied biological application uses training-target conditions and fails to recover the sign of the reported chromosome-level contrast, while no harmonized comparison with available public approaches establishes the claimed methods advance. The work could interest communities beyond *C. elegans* if it demonstrates reliable transfer of a sequence model to a distinct organismal RNA-seq setting and enables new biological inference. The packet currently demonstrates an internally controlled engineering and evaluation framework more clearly than that broader conclusion.
- **Who would be interested in the results, and why** Immediate interest would come from researchers developing sequence-to-track models, studying *C. elegans* functional genomics, adapting large genomic models across organisms, or building reproducible computational biology workflows. A wider methods readership could care if the work shows that model adaptation delivers reliable gains over available approaches, generalizes beyond the development collection, and supports a biological conclusion that could not be reached otherwise. Those wider-interest conditions are not yet demonstrated in the supplied evidence.
- **Major strengths** The provenance inventory and uniform grouping provide a clearly counted data basis. The one-time evaluation records full stated coverage across 83 aligned cores, six chromosomes, and 241 tracks. P9 includes a size-matched from-scratch model and factorial component analysis, with folds treated as the inference unit and seeds nested within folds. The packet also discloses important limitations rather than obscuring them. These include prior exposure of chromosome X, the internal status of P10, the descriptive status of track outcomes, the non-mechanistic meaning of the interaction term, and the absence of real variant-scoring validation. Commit hashes, locked-test access controls, and preservation of a failed P10 attempt strengthen the integrity trail.

## Major Concerns

### R3-M1 [experimental-design]

**Severity** Major

**Blocking** Yes

**Claim pointer** The prospective Article would use P10 to establish an application to an important biological question and support the method's capacity to recover biologically meaningful regulation.

**Evidence pointer** Packet section `P10 internal DPY-27 application`, including the condition status and four reported endpoints. Packet section `Materials not present in this packet`, including the absent unseen-condition and independent biological validations.

**Concern** P10 is an internal evaluation on matching development validation cores. Both DPY-27 RNAi and vector RNAi were training targets. More importantly, the observed chromosome-X-minus-autosome contrast is positive at `0.003065`, whereas the predicted contrast is negative at `-0.002740`. The predicted-observed gene-level Spearman correlation is only `0.0846`. Direction concordance of `0.6674` does not resolve the contradictory chromosome-level sign or the weak gene-level ranking. On the evidence shown, P10 is an informative failure analysis or internal diagnostic, not strong validation of biological fidelity or discovery.

**Why it matters** The supplied journal criteria require strong validation and an application to an important biological question. The only visible application does not establish that the method reproduces the principal reported biological pattern or generalizes to a new condition. Presenting it as affirmative biological validation would overstate the evidence and weaken both technical soundness and scientific importance.

**Resolution test** Demonstrate a prespecified application on data isolated from training by biological condition and preferably by study or laboratory. Show recovery of the direction, magnitude, and gene-level structure of biologically relevant endpoints, with orthogonal RNA, reporter, perturbational, or equivalent independent validation appropriate to the claim. Otherwise, describe P10 explicitly as an internal limitation analysis and remove or narrow claims of biological discovery.

### R3-M2 [novelty-significance]

**Severity** Major

**Blocking** Yes

**Claim pointer** The AlphaGenome-style adaptation is a novel method or a significant improvement over established sequence-to-track approaches for biological data.

**Evidence pointer** Packet sections `Matched development comparisons` and `P9 component and size-matched evidence`. Packet section `Materials not present in this packet`, including the absent harmonized public-model comparison and complete prior-art analysis.

**Concern** The supplied comparisons establish differences among repository configurations A, B, and C and isolate several components of Model B. They do not establish performance relative to technically applicable available approaches. The packet explicitly reports no completed Enformer, Borzoi, or other public-model comparison under a harmonized benchmark, and no prior-art analysis sufficient to determine which components are novel. The size-matched C comparison addresses parameter count within the project but cannot substitute for comparison with existing methods or for novelty positioning.

**Why it matters** Internal ablation can explain where an in-repository gain originates, but it cannot show that the method advances the state of the art or delivers the significant improvement expected of a methods Article. Without this distinction, originality and broad scientific importance remain unassessable.

**Resolution test** Define the claimed advance component by component, provide a complete and accurate prior-art comparison, and evaluate technically applicable public methods and simple baselines under the same input data, splits, target definitions, metrics, and leakage controls. Report accuracy together with trainable parameters, compute, and adaptation requirements. If a public model cannot be evaluated fairly, document the incompatibility and narrow the comparative claim.

### R3-M3 [claim-moderation]

**Severity** Major

**Blocking** No

**Claim pointer** The method has general applicability beyond the particular genomic regions and 241-track collection used for development.

**Evidence pointer** Packet sections `Shared project description` and `Final locked evaluation`. Packet section `Materials not present in this packet`, including the absent study- or laboratory-isolated benchmark, unseen biological condition, and independent species validation.

**Concern** The chromosome-fold design tests performance across held-out genomic regions within the same provenance collection and fixed track set. The packet provides no evidence for transfer to an unseen study, laboratory, biological condition, assay collection, or species. Chromosome X also was exposed during legacy work and is not pristine for the project as a whole. The reported final evaluation therefore supports a bounded within-collection genomic generalization claim, not a broad claim of experimental or cross-organism applicability.

**Why it matters** Readers outside the immediate project need to know what kind of transfer the method actually achieves. Conflating genomic-region holdout with transfer across biological or technical distributions would inflate the method's expected usefulness and obscure the conditions under which another laboratory could rely on it.

**Resolution test** Add an evaluation isolated by study or laboratory and by a biologically meaningful unseen condition, with an explicit leakage audit and prespecified endpoints. Cross-species evidence would strengthen a broader adaptation claim but is not required if the Article consistently limits generality to the demonstrated *C. elegans* setting. Qualify chromosome X as project-exposed wherever final-test independence is described.

### R3-M4 [statistical-rigor]

**Severity** Major

**Blocking** No

**Claim pointer** The fold-first bootstrap intervals provide reliable uncertainty for the P6B and P9 comparative and component effects.

**Evidence pointer** Packet section `Matched development comparisons`, including the post hoc submission analysis. Packet section `P9 component and size-matched evidence`, including the five-fold inference unit and 10,000-replicate nested bootstrap.

**Concern** The packet correctly identifies only five genomic folds as independent inference units, but it reports aggregate bootstrap intervals without the five fold-level paired effects or sensitivity to any single chromosome. Increasing bootstrap replicates does not increase the number of independent folds. The adequacy and stability of the intervals are therefore difficult to assess, particularly for effects near zero such as LoRA-only and learned-1-bp-head-only.

**Why it matters** Component attribution is central to explaining why Model B works. If small estimated effects or their intervals are driven by one chromosome, the architectural interpretation and its portability may change even when the overall Model B advantage remains.

**Resolution test** Show all fold-level paired effects and seed dispersion, report leave-one-fold-out sensitivity, and state the inferential population supported by five named chromosomes. Provide a justification or sensitivity analysis for interval construction with five independent units. Moderate claims for small component effects when their direction is not stable across folds.

### R3-M5 [reproducibility]

**Severity** Major

**Blocking** No

**Claim pointer** The method is reproducible and usable by communities beyond the originating repository.

**Evidence pointer** Packet section `Reproducibility and integrity controls`. Packet section `Materials not present in this packet`, including the absent common-system installation evidence, container, public archival release, and independent user reproduction.

**Concern** Controller records, commit hashes, output hashes, and test-lock checks provide strong internal provenance. They do not yet establish that another group can obtain the software, install it, reconstruct the environment, run a small example, or reproduce a principal result. No versioned archival release, containerized environment, common-system installation check, or independent reproduction is supplied.

**Why it matters** Practical adoption is part of the significance case for a technical paper. Internal auditability protects integrity, while external reproducibility determines whether the method can affect a community of researchers.

**Resolution test** Release a versioned and archived code package with environment or container specifications, model and data access instructions, a tractable end-to-end example, expected outputs and checksums, and documented resource requirements. Include an independent clean-environment reproduction of at least one principal evaluation or application result.

### R3-M6 [claim-moderation]

**Severity** Major

**Blocking** No

**Claim pointer** The variant-scoring functionality can support scientific inference about regulatory variants.

**Evidence pointer** Packet section `Variant-scoring software`. Packet section `Materials not present in this packet`, including the absent external variant evidence and real frozen-VCF execution.

**Concern** The variant-scoring implementation has synthetic CPU tests but no registered real VCF plus checkpoint run and no evaluation against cis-eQTLs, regulatory alleles, reporters, CRISPR results, or another independent molecular readout. This establishes software behavior under synthetic tests, not empirical validity of variant-effect scores.

**Why it matters** Variant interpretation would substantially expand the paper's biological and interdisciplinary appeal, but an unvalidated capability can mislead readers about the present method's scientific reach. Its inclusion as a headline result would outpace the evidence.

**Resolution test** Evaluate a frozen public variant or haplotype set with prespecified metrics, matched baselines, and independent measured effects. Until then, identify variant scoring as an unvalidated software extension or omit it from the principal claims.

## Minor Comments

### R3-m1 [writing-clarity]

**Severity** Minor

**Affected element** The naming and organization of models, losses, phases, gates, and evaluation units.

**Evidence pointer** Packet sections `Shared project description`, `Matched development comparisons`, `P9 component and size-matched evidence`, and `P10 internal DPY-27 application`.

**Issue** Labels such as P6B, P6C, P9, P10, R9, R10, A, B, C, G0005, G0006, `paper loss`, `aligned cores`, and `composite score` are efficient audit identifiers but do not communicate the scientific design to a nonspecialist. Their density obscures which intervention, comparator, target, and outcome belong to each result.

**Required correction** Lead each section and figure with descriptive scientific names. Define each metric and unit at first use, state what a higher value means, and retain phase or gate identifiers only as secondary reproducibility labels.

### R3-m2 [writing-clarity]

**Severity** Minor

**Affected element** Explanation of the dual-resolution prediction task and headline correlations.

**Evidence pointer** Packet sections `Shared project description` and `Final locked evaluation`.

**Issue** The packet reports 1-bp and 128-bp targets and two mean per-track correlations without explaining their distinct biological purpose, the construction of gene-exon versus 128-bp outcomes, or what performance at the reported magnitude enables in practice.

**Required correction** Add a concise conceptual explanation of both output resolutions, their intended biological uses, the denominators and averaging scheme for each correlation, and an interpretable baseline or reference point. A nonspecialist should be able to connect each metric to a concrete scientific task.

### R3-m3 [claim-moderation]

**Severity** Minor

**Affected element** Description of the locked final evaluation and chromosome X.

**Evidence pointer** Packet sections `Shared project description` and `Final locked evaluation`.

**Issue** The active workflow uses a separately gated one-time evaluation, but chromosome X was exposed during legacy work and is not pristine for the project as a whole. An unqualified phrase such as `unseen test chromosome` or `pristine held-out test` would therefore be inaccurate.

**Required correction** Distinguish the active v2 access lock from project-wide historical exposure every time test independence is summarized. State precisely what information was unavailable to v2 model selection and what had been exposed previously.

### R3-m4 [figures-and-tables]

**Severity** Minor

**Affected element** The final figure sequence and visual explanation of the method.

**Evidence pointer** Packet section `Packet scope`, which states that no final figure order was supplied. Packet sections `Shared project description` and `Verified evidence base`, which contain the elements that require visual linkage.

**Issue** Figure accessibility is not assessable from the packet. The prospective paper nevertheless needs a visual route from sequence inputs and dual-resolution targets through chromosome-based development, locked evaluation, component tests, and biological application. Without that route, the phase-oriented evidence structure is likely to remain inaccessible to readers outside the immediate modeling workflow.

**Required correction** Include a simple conceptual schematic and a figure sequence organized around scientific questions rather than project phases. Show paired fold-level comparisons and the agreement or disagreement between observed and predicted biological endpoints, while defining uncertainty and replicate units in each legend.

- **Technical failings that need to be addressed before the case is established** R3-M1 and R3-M2. The present packet lacks a convincing independent biological application and a harmonized comparison that establishes novelty and performance relative to available approaches.
- **Assessment against Nature-style criteria** **Originality** Not assessable from the supplied material because component-level prior art and external method comparisons are absent. The internal component analysis is informative but does not establish novelty. **Scientific importance** Potentially substantial if the adaptation enables reliable RNA-track prediction and new biology in *C. elegans*. The supplied P10 result does not yet demonstrate that impact. **Interdisciplinary readership** Immediate interest is credible across computational genomics, model-organism biology, and sequence-model adaptation. Broader interest remains conditional on external performance, transferable validation, and a biologically informative application. **Technical soundness** Internal controls, audit records, matched analyses, and explicit inference units are strengths. R3-M1 and R3-M2 prevent the present evidence from establishing the central methods case, while R3-M3 through R3-M6 limit generality, inferential confidence, reuse, and application scope. **Readability for nonspecialists** Direct assessment is not possible without the manuscript, abstract, and figures. The packet's reliance on internal phase codes and unexplained metrics indicates that a substantial conceptual rewrite will be needed if this vocabulary carries into the Article.
- **Recommendation posture** Promising but the broad-interest and technical case is currently not established from the provided evidence. I would be supportive of further consideration after the blocking biological-validation and comparative-benchmark issues are resolved and the central claim is made explicit. Final journal fit and prose quality remain outside this bounded assessment.
