---
name: alphagenome-nature-methods-voice
description: Preserve evidence-bound, natural English prose for the AlphaGenome C. elegans Nature Methods manuscript. Use when drafting, revising, translating, polishing or auditing its title, abstract, introduction, Results, Discussion, Online Methods, figure legends, Data/Code Availability statements, cover letter or reviewer response; also use when a draft sounds templated, over-promotional or AI-generated while scientific claims, terminology, numerical values and evidence boundaries must remain exact.
---

# AlphaGenome Nature Methods Voice

Use this project-local skill after the scientific evidence is assembled. Do not
use it to invent a method name, benchmark result, citation, biological
interpretation or missing experiment.

## Start With the Evidence Contract

1. Read `../../plans/00_storyline_and_claim_contract.md` and
   `../../plans/TERMINOLOGY_LEDGER.md`.
2. Read the plan for the requested manuscript part under `../../plans/parts/`.
3. Load `references/project-boundaries.md` before changing scientific claims.
4. If the text contains an `EXTERNAL_*` placeholder, preserve it and report a
   drafting blocker. Do not replace it with generic positive wording.
5. If a result is unsupported, mark the sentence `[EVIDENCE NEEDED]` rather
   than softening it into an apparently factual sentence.

## Choose the Task

| Task | Do this |
| --- | --- |
| Draft from approved evidence | Build one paragraph job at a time. State claim, evidence, condition and implication in that order. |
| Revise a human draft | Keep the original when it is already direct, precise and natural. Make only changes that improve clarity, logic or scientific fit. |
| Translate Chinese notes | Split each note into claim, evidence, condition and implication before writing English. Translate intent, not Chinese syntax. |
| De-AI a section | Read `references/voice-rubric.md`, remove empty emphasis and mechanical connectors, then run the local audit script. |
| Audit before sign-off | Run `scripts/audit_voice.py <file> --strict`; resolve blockers manually and treat style flags as human review prompts. |

## Revision Workflow

1. **Lock facts.** Preserve all numbers, units, p values, confidence intervals,
   figure/table references, method labels and citations unless the owner supplies
   corrected evidence.
2. **State one message.** Give each paragraph exactly one job: context, gap,
   approach, result, comparison, interpretation, implication or boundary.
3. **Use concrete syntax.** Prefer named actors, measured objects and observed
   outcomes. Replace vague praise with the comparator, task or effect that
   supports the claim.
4. **Repair flow.** Let causal, contrastive or chronological relations connect
   sentences. Delete formulaic signposts rather than replacing them mechanically.
5. **Preserve a human rhythm.** Mix short factual sentences with longer
   explanatory ones when the logic benefits. Do not make every sentence follow
   the same `We verb X, which verb Y` pattern.
6. **Keep scope visible.** Put the supported boundary in the sentence that makes
   the claim. Do not append generic caveat paragraphs to create a cautious tone.
7. **Report material edits.** Return revised text, a compact modification log,
   and an evidence/terminology issue list only when it contains real issues.

## AlphaGenome-Specific Rules

- Use `full Model B`, `Model A`, and `size-matched Model C` exactly as defined
  in the ledger. Do not call all three `baselines` without their role.
- Use `study-isolated benchmark` only after its source-isolation rule is stated.
  Use `genomic-block validation` for within-collection evidence.
- State five genomic folds as the primary blocked unit and three seeds as nested
  algorithmic repeats. Treat the 241 track distribution as descriptive.
- Do not claim that the 1-bp head is necessary.
- Describe P10 only as the `DPY-27 internal diagnostic`; it did not recover the
  prespecified chromosome-level response and cannot be promoted to biological
  validation.
- Attribute a variant effect to a tested molecular endpoint, not to causal
  disease biology, unless a causal experiment is present.

## Style Guardrails

Use `references/voice-rubric.md` for the detailed rubric. In brief:

- Keep plain technical words when they are more exact than impressive words.
- Remove empty intensifiers, vague attribution and generic triumphal framing.
- Prefer direct transitions based on evidence rather than `First and foremost`,
  `It is worth noting`, `Notably`, `In conclusion` or a chain of numbered
  rhetorical steps.
- Do not force prose into bullet lists in the manuscript body. Lists remain
  appropriate in Methods, checklists, figure legends and supplementary material
  when they improve reproducibility.
- Do not substitute synonyms merely to create variety; terminology consistency
  wins over stylistic variation.
- Retain a strong sentence unchanged when the audit finds no concrete problem.

## Output Contract

For a requested revision, return:

1. Revised English text.
2. A concise modification log only for material changes.
3. `EVIDENCE NEEDED` or `TERM CONFLICT` items, if any.
4. A statement that no numerical or factual content was changed, or a precise
   list of owner-authorized factual changes.

For a draft request, return the draft plus a claim-evidence map. Do not write a
full section if its external benchmark fields remain unfilled; provide a
scaffold with preserved placeholders instead.

## Resources

- `references/project-boundaries.md`: exact claim boundaries and protected
  internal results.
- `references/voice-rubric.md`: high-threshold natural-prose rules and source
  inspiration statement.
- `references/revision-output.md`: output patterns for prose, Methods and
  figure legends.
- `scripts/audit_voice.py`: deterministic text audit; it flags, but never
  rewrites, potential issues.
