# Part 3: Results Plan

## Results Architecture

Target 1,850-2,050 words across five topical subsections. Each subsection
opens with the question being tested, then reports the corresponding figure
before interpretation. The text must use past tense for observations and keep
the external benchmark as the primary performance evidence.

| Subsection | Question | Main display | Target words | Evidence status |
| --- | --- | --- | --- | --- |
| 3.1 | What resource and prediction task were created? | Fig. 1 | 300-350 | Internal resource/model specification verified. |
| 3.2 | Does the full model outperform fair alternatives on genuinely unseen data? | Fig. 2 and Table 1 | 520-600 | External benchmark placeholders until freeze. |
| 3.3 | Which design choices account for the gain? | Fig. 3 | 360-430 | P9 verified. |
| 3.4 | Can frozen sequence differences prioritize an independent regulatory endpoint? | Fig. 4 | 420-500 | External biology placeholders until freeze. |
| 3.5 | In which external contexts is the model reliable enough to use? | Fig. 5 | 260-320 | External robustness placeholders until freeze. |

## 3.1 Build a Species-Adapted RNA Prediction Resource and Model

**Question:** Can a heterogeneous public worm RNA-seq collection be converted
into a traceable multi-track sequence prediction task without sacrificing
resolution or provenance?

**Paragraph jobs:**

1. State the finalized resource scale: 482 verified RNA-seq runs, uniformly
   grouped into 241 biological RNA tracks, all on WBcel235. Explain the
   biological axes represented only with figures/tables, not a long list.
2. Define the task: a 131,072-bp sequence window produces 1-bp and 128-bp RNA
   coverage targets, with gene/exon aggregation as a biologically interpretable
   readout.
3. Define the technical move: the frozen trunk gains a worm embedding, final
   tower LoRA and dual-resolution RNA head. Explain each module by the question
   it answers, not by implementation trivia.
4. State that training, validation and the external protocol were frozen before
   result inspection. Keep detailed split logic in Fig. 1 and Online Methods.

**Figure logic:** Fig. 1 must contain one held-out external coverage example
when available. A schematic without a real prediction is insufficient.

## 3.2 Establish Performance on a Study-Isolated Benchmark

**Question:** Does the full model retain an advantage when all model choices are
fixed and the test set comes from an unseen study/laboratory source including
an important unseen biological condition?

**Paragraph jobs:**

1. Lead with the isolation definition: `EXTERNAL_SOURCE_ISOLATION_RULE`,
   `EXTERNAL_UNSEEN_CONDITION` and `EXTERNAL_LEAKAGE_AUDIT_RESULT`. Name what
   was excluded from training, tuning and selection.
2. State the shared protocol: reference build, input window, target definition,
   trainable-data budget, optimization budget and test set were matched across
   eligible models. Refer readers to Table 1 for eligibility and non-adaptable
   methods.
3. Report the primary full-model comparison as
   `EXTERNAL_MODEL_B_EFFECT_AND_CI`; then report all eligible public baselines
   as `EXTERNAL_PUBLIC_BASELINE_RESULTS`. Use a comparison verb only after the
   direction and interval are populated.
4. Report study-wise and condition-wise consistency. Do not replace a
   source-level result with 241 track points.
5. Close with one sentence explaining the practical performance object: gene
   abundance, binned coverage or both, in the exact order prespecified.

**Mandatory boundaries:**

- Enformer and Borzoi appear only if their adaptation is technically faithful
  under the shared protocol. List ineligible models with a precise reason, not
  an empty comparison row.
- The internal six-chromosome locked evaluation may appear as a supporting
  reference, not as external evidence.
- Parameters and compute are reported beside accuracy; an advantage must not be
  hidden behind a materially larger training budget.

## 3.3 Attribute the Full-Model Gain to the Adaptation Components

**Question:** Which adaptation components produce the measured advantage once
the external result establishes that the method is useful?

**Paragraph jobs:**

1. State that this is a controlled internal component-attribution experiment,
   not an independent generalization test.
2. Report the full-model advantage over size-matched Model C: primary paired
   effect `+0.092017` in favor of full Model B, with the registered five-fold
   interval from P9.
3. Report LoRA as the dominant main effect (`+0.101380`) and worm embedding as
   a smaller positive main effect (`+0.009800`), using paired fold/seed wording.
4. Report the interaction strictly as statistical non-additivity; do not assign
   a biological mechanism to it.
5. State that the learned-1-bp-head-only primary interval crossed zero, so the
   data do not support a necessity claim for that head on the registered primary
   score.

**Source evidence:** P9 `paired_effects.tsv`, `factorial_effects.tsv`,
`per_track.tsv` and the audited component-attribution figure bundle.

## 3.4 Validate Regulatory Variant Prioritization Outside the Training Collection

**Question:** Do frozen ref/alt or haplotype predictions recover a prespecified
independent regulatory signal rather than merely changing a model output?

**Paragraph jobs:**

1. Define the external truth source and its independence:
   `EXTERNAL_VARIANT_DATASET_ID`, `EXTERNAL_VARIANT_ISOLATION_RULE` and the
   frozen annotation/build alignment.
2. State the score before results: ref/alt window prediction, reverse-complement
   ensemble, gene/exon aggregation and optional ISM, all measured without
   retuning to labels.
3. Report direction concordance and its blocked confidence interval.
4. Report ranking/discrimination against a matched-null set using AUROC/AUPRC
   and its confidence interval. Explain the matching variables.
5. Present one preselected locus-level coverage example and one interpretable
   sequence attribution. The example illustrates the quantitative test; it does
   not replace it.

**Primary endpoint:** a same-species cis-eQTL set. Do not use trans-eQTLs as
the primary proof of a local sequence-difference claim. An exact
regulatory-allele perturbation panel can enter as an orthogonal supplementary
analysis only if it is independently frozen.

## 3.5 Define the External Use Domain and Release-Ready Behavior

**Question:** In which study, condition, signal and genomic-context strata is
the model's external advantage stable, and what should users expect from it?

**Paragraph jobs:**

1. Report source- and condition-level effects rather than only one pooled score.
2. Show performance/calibration across prespecified signal and annotation
   strata, including the difficult strata that do not retain the global gain.
3. Report inference cost, parameter counts and output format for a practical
   scoring run.
4. Point to the release artifact: model weights, track metadata, exact software
   environment, test fixtures and a minimum reproducible prediction command.

This subsection is not a software screenshot. Its evidence must establish a
usable domain of operation and visible failure modes.

## Results-Level Prohibitions

- Do not turn a favorable global average into a claim about every study,
  condition, chromosome or track.
- Do not call the 241 tracks independent biological replicates.
- Do not discuss why the method works in the external benchmark beyond the
  P9-supported component statements; save broader interpretation for Discussion.
- Do not mention the DPY-27 diagnostic here as a successful result. It belongs
  in Extended Data and the boundary paragraph of Discussion.
