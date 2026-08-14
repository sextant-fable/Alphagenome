# Part 4: Discussion Plan

## Target

Five paragraphs, approximately 550-650 words, without subheadings. This is a
confident interpretation of evidence already shown, not a defensive recitation
of caveats and not a second Results section.

## Paragraph Map

| Paragraph | Job | Content |
| --- | --- | --- |
| 1 | Central advance | Restate the method, its study-isolated performance and its independently tested regulatory use in one compact argument. |
| 2 | What changed technically | Interpret the P9 finding that LoRA is the major measured contributor and worm conditioning gives a smaller stable gain. Explain this as an empirical transfer result, not a claimed molecular mechanism. |
| 3 | Relation to prior work | Position the work relative to long-context sequence models and RNA-coverage modeling: the contribution is a rigorously evaluated species adaptation and workflow, not a claim to replace all existing models. |
| 4 | Scope and visible boundary | Define the supported data regime, external source rule and biological endpoint. State the internal DPY-27 diagnostic accurately: it did not recover the prespecified chromosome-level response and therefore is presented as a negative diagnostic rather than validation. |
| 5 | Reuse and next questions | Close with the released resource, reproducible scoring interface and the concrete experiments now enabled, such as extending to perturbation collections or another model organism. |

## Confidence Calibration

Use `show` or `demonstrate` only for results whose frozen benchmark or P9
analysis directly supports the statement. Use `supports`, `is consistent with`
or `suggests` for architectural interpretation. Do not introduce generic
hedges such as `may potentially be useful`; say what the released model can do
under the stated protocol.

## Specific Interpretive Moves

- The full model's gain is a result of controlled adaptation, not proof that a
  worm embedding encodes a particular biological process.
- The external benchmark establishes an operational generalization regime. It
  does not establish cross-species transfer unless a cross-species study was
  actually included.
- The variant result establishes utility for the exact local regulatory endpoint
  tested. It does not prove disease causality or solve trans-regulation.
- The DPY-27 result is a useful boundary because the model failed a prespecified
  internal condition-contrast endpoint. State it once, concretely, without
  letting it displace the independent biological evidence.

## Do Not Do

- Do not repeat all Fig. 2 numbers.
- Do not promise genome-wide causal inference, clinical translation or universal
  regulatory-code decoding.
- Do not claim a field first unless a completed prior-art audit establishes a
  scoped first claim.
- Do not introduce a new biological locus, external dataset or analysis.
