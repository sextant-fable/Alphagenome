# Nature Methods manuscript strategy

## Draft

### One-sentence argument supported now

In *C. elegans* sequence-to-RNA-track modeling, worm-specific embedding and
LoRA adaptation improved within-collection genomic-block prediction over a
frozen-trunk baseline and a size-matched from-scratch model, with LoRA
accounting for most of the measured gain; external generalization and
biological utility remain unestablished.

This is the strongest defensible argument from the completed evidence. It is
not yet a Nature Methods Article argument because it does not contain an
independent benchmark, a comparison with available public approaches, or a
successful biological application.

### Target argument after the required evidence

In *C. elegans* sequence-to-function modeling, we show that a species-adapted
sequence-to-track model improves RNA prediction over harmonized public and
from-scratch baselines, generalizes to an independently held-out study or
condition, and ranks regulatory variant effects validated by external
expression genetics or perturbation data.

The final clause is conditional. It may enter the title, abstract, or cover
letter only after the external benchmark and biological validation pass their
frozen primary endpoints.

## Terminology ledger

| Canonical term | First-use definition | Decision |
|---|---|---|
| Model A | frozen-trunk RNA head baseline | Use `Model A` only after this definition. |
| Model B | worm-embedding, LoRA-adapted, dual-resolution RNA model | Use as the current full model; do not coin a method name without author approval. |
| size-matched Model C | from-scratch baseline matched to Model B trainable parameter count | Do not call the original P6B C configuration parameter matched. |
| genomic-block validation | evaluation on held-out genomic cores within the 482-run provenance collection | Do not shorten to external validation. |
| locked final evaluation | the one-time v2-gated P6C evaluation | Do not call it project-wide pristine because chromosome X was exposed in legacy work. |
| independent benchmark | a study- or laboratory-isolated benchmark outside the 482-run collection | Reserve this term for future external data. |
| primary biological score | the registered composite development endpoint | Define its formula and direction at first use. |
| DPY-27 internal diagnostic | the P10 comparison of DPY-27 RNAi with vector RNAi on training-target conditions | Do not call it biological validation or response recovery. |
| variant-scoring CLI | frozen-spec ref/alt and ISM software | Call it validated only after a real frozen VCF plus external truth run. |

## Section outline

1. **Introduction.** Define the need for organism-specific sequence-to-track
   modeling, identify the data and transfer bottleneck, position the proposed
   adaptation against verified prior work, and state a claim that matches the
   completed external evidence.
2. **Method overview.** Present one scientific schematic from WBcel235 sequence
   and 241 RNA targets through worm embedding, LoRA adaptation, dual-resolution
   outputs, frozen splits, and evaluation. Keep controller phase identifiers in
   Methods and reproducibility notes, not in the main scientific narrative.
3. **Audited resource and evaluation design.** Establish 482 RNA-seq runs, 241
   groups, on-demand 1-bp and 128-bp targets, leakage controls, fold units, and
   the historical chromosome-X boundary.
4. **Fair competitive benchmark.** Lead with the future independent benchmark
   and harmonized public baselines. Report accuracy, uncertainty, total and
   trainable parameters, compute, preprocessing, tuning allowances, and failure
   cases. The internal A/B/C comparison is supporting evidence, not a substitute.
5. **Component attribution.** Use the completed P9 matrix to show that LoRA is
   the dominant measured contribution, worm embedding gives a smaller stable
   gain, and the learned-1-bp-head-only primary interval includes zero. Keep the
   interaction explicitly statistical.
6. **Independent biological application.** Use a prespecified external cis-eQTL
   or regulatory-allele application, not P10, as the main biological test.
   Report direction, ranking, null discrimination, calibration, and blocked
   uncertainty. Add reporter, CRISPR, or independent RNA validation if available.
7. **Reproducibility and practical use.** Provide an archived release,
   environment or container, checkpoint and data access, a small end-to-end
   example, expected hashes, resource requirements, and an independent clean
   reproduction.
8. **Discussion.** Interpret the demonstrated transfer regime, explain what
   LoRA contributes empirically, state where generalization stops, disclose the
   negative DPY-27 diagnostic, and separate future variant capability from
   validated utility.

## Proposed figure sequence

| Display | Scientific job | Current status |
|---|---|---|
| Figure 1 | Data, model, split, and evaluation schematic | Needed; organize by scientific workflow rather than P phases. |
| Figure 2 | Independent benchmark and harmonized public-model comparison | Blocking evidence not yet generated. |
| Figure 3 | Component attribution across five folds | Ready in `results/v2_p9_submission_evidence_figure/`. |
| Figure 4 | Independent regulatory-variant or perturbation application | Blocking evidence not yet generated. |
| Extended Data 1 | Complete internal P6B comparison, strata, and track distributions | Ready in `results/v2_p6b_submission_analysis/`. |
| Extended Data 2 | DPY-27 internal failure analysis | Ready as a diagnostic, not as affirmative validation. |
| Extended Data 3 | Reproducibility, hash, leakage, and sensitivity audits | Partly ready; independent clean reproduction is missing. |

## Title strategy

### Defensible now

`LoRA adaptation of a sequence-to-track model for C. elegans RNA prediction`

This title is accurate but may not support the breadth expected for Nature
Methods without the missing comparisons and application.

### Conditional after successful external validation

1. `Species-adapted sequence-to-track modeling predicts C. elegans RNA programs`
2. `A species-adapted sequence model links C. elegans genomes to RNA regulation`
3. `Sequence-to-track adaptation enables regulatory variant scoring in C. elegans`

Titles 2 and 3 require external biological evidence. They must not be used if
the cis-eQTL or perturbation endpoint fails.

## Abstract architecture

Use `challenge -> contribution -> decisive comparison -> independent
application -> boundary`.

1. One sentence on the limited availability of sequence-to-track models for
   compact, deeply annotated model organisms.
2. One sentence defining the transfer and validation gap in existing practice,
   grounded in verified prior art.
3. One sentence introducing the worm embedding, LoRA adaptation, and
   dual-resolution RNA prediction workflow.
4. One sentence reporting the harmonized external comparison and uncertainty.
5. One sentence reporting the independent biological endpoint.
6. One bounded implication sentence naming the demonstrated species, assay,
   data regime, and released software.

Do not lead the abstract with P phase codes, the 0.675/0.664 internal final
correlations, or the negative DPY-27 result. Those quantities establish trust
and boundary, but neither supplies the missing field-level comparison or
biological application.

## Assumptions or missing inputs

- The paper is treated as a prospective Nature Methods Article rather than an
  Analysis or Resource.
- No author-approved method name or exact novelty claim has been supplied.
- Enformer and Borzoi are planned comparators, but their final eligibility and
  harmonized protocols remain to be frozen.
- The independent benchmark and external biological application remain
  blocking experiments.
- A full manuscript, author list, affiliations, code/data release locations,
  contribution statement, and competing-interest statement are not available.

## Claim-evidence map

| Claim | Evidence | Status |
|---|---|---|
| Model B exceeds frozen A within the registered development design. | Paired primary effect `+0.111180`, five-fold hierarchical CI excludes zero. | supported |
| Model B exceeds a trainable-parameter-matched from-scratch model. | C-minus-B `-0.092017`; trainable parameter difference 0.87%. | supported within the internal design |
| LoRA is the dominant measured component contribution. | LoRA main effect `+0.101380`; no-LoRA-minus-B `-0.094392`. | supported within the 2 x 2 P9 design |
| Worm embedding adds a smaller gain. | Main effect `+0.009800`; LoRA-only-minus-B `-0.002812`. | supported within the P9 design |
| The dual-resolution head is necessary. | Learned-1-bp-head-only effect `-0.001808`, CI includes zero. | not supported as a primary-score necessity claim |
| The method recovers DPY-27 dosage compensation. | Predicted X-autosome contrast has the opposite sign from observed. | contradicted by the prespecified P10 primary endpoint |
| The method generalizes across studies or conditions. | No study-isolated or unseen-condition benchmark. | needs evidence |
| The method outperforms available public approaches. | No completed harmonized public-model comparison. | needs evidence |
| Variant scores identify functional regulatory alleles. | Software and synthetic tests only. | needs evidence |
| The workflow is independently reproducible. | Strong internal audit trail, no clean external reproduction. | needs evidence |

## Why this structure

- It answers the methods reader's questions in order: relevance, novelty,
  trust, reuse, and meaning.
- It places fair comparison and independent application before internal
  ablation, so the paper does not mistake engineering attribution for field
  superiority.
- It uses the negative DPY-27 result to define a real failure boundary instead
  of weakening the main biological story with an unsupported positive claim.
- It keeps every ambitious sentence conditional on the experiment that would
  make it true.

## To redirect this strategy

Identify the claim, figure job, or section order that should change. Revise
that element only and keep the verified terminology and evidence boundaries
fixed.
