# Nature Methods Manuscript Workspace Rules

## Scope

This directory manages the prospective *Nature Methods* Article for the
AlphaGenome-style *C. elegans* sequence-to-RNA-track project. It is a separate
manuscript workspace, not a new model-development area. Do not run training,
change v2 phase state, read the locked final-test block, regenerate biology
results, or alter source data from this directory.

The owner has authorized an initial manuscript draft. Prose may now be created
under `manuscript/`, but this is an internal draft and not a submission-ready
record. Keep concrete section plans in `plans/`, never in this file. The
historical assumed-results scaffold is archived: a readable or compiled
manuscript may contain only verified results and real figures.

## Evidence Contract

- Treat the frozen external benchmark and independent biological application as
  planned manuscript inputs, not as verified repository results, until their
  signed manifests and source data are supplied.
- Use `plans/08_external_benchmark_placeholder_contract.md` only to plan the
  missing external benchmark and biological application. Do not place its
  placeholders, assumed effects, guessed dataset names or generic positive
  substitute language in a compiled manuscript or figure legend.
- Distinguish `external benchmark`, `internal genomic-block validation`, and
  `locked final evaluation`. They are not interchangeable.
- Preserve the P10 DPY-27 result as a negative internal diagnostic. It cannot
  be called a successful biological validation, dosage-compensation recovery or
  independent generalization.
- Use the P9 terminology and inference units exactly: five genomic folds are
  the primary blocked units; three seeds are nested algorithmic repeats; 241
  tracks are descriptive outcomes rather than independent biological repeats.
- Do not claim that the learned 1-bp head is necessary: its primary-score
  interval includes zero in P9.
- Do not call the original P6B Model C parameter matched. Only P9
  `C_size_matched` has that status.

## Local Layout

- `literature/`: verified prior-art and visual-structure audit.
- `plans/`: master narrative, section plans, figures, evidence registry and
  submission gates.
- `template/`: official template provenance. `template/vendor/` is an ignored,
  unmodified third-party artifact.
- `skills/alphagenome-nature-methods-voice/`: project-local writing skill.
- `manuscript/`: owner-approved initial draft, assumption register and later
  template-ready source. The assumption register is not part of the submitted
  manuscript.

Keep planning files in English. Use stable, descriptive filenames. Do not put
raw data, model weights, generated figures or external benchmark downloads in
Git. Store their path, checksum, command and immutable manifest in the future
release record instead.

## Nature Skills Routing

Use the smallest relevant skill set and read each selected `SKILL.md` before
acting. The following recipes are mandatory for this workspace.

| Need | Required skill sequence | Required output control |
| --- | --- | --- |
| New literature, citation or prior-art check | `nature-academic-search` -> `nature-citation` when citations enter prose | Start with `multi-source-search`; use primary papers and official journal pages; record DOI, URL and verification date. |
| Lawful paper or template download | `nature-downloader` | Obtain only official, open-access or user-authorized material; record source, format, SHA-256 and SI choice/status. Template downloads have `SI: not applicable`. |
| Manuscript outline or prose | `nature-writing` -> `$alphagenome-nature-methods-voice` | Build the one-sentence argument and terminology ledger first; draft only from frozen evidence; retain explicit boundaries. |
| Figure planning or plotting | `nature-figure` -> `nature-statistics` | Create a figure contract before code. Before the first actual plot, resolve the user's persistent Python/R preference as required by `nature-figure`. |
| Data/code availability and FAIR release | `nature-data` | Keep dataset accessions, checksums, repositories, code release, license and model-weight conditions explicit. |
| Statistical reporting | `nature-statistics` | Verify units, pairing, bootstrap design, intervals, multiplicity and figure legends against the frozen analysis. |
| Pre-submission adversarial audit | `nature-reviewer` -> `nature-response` only after a decision letter | Freeze reviewer contexts before synthesis; do not turn reviewer speculation into manuscript facts. |
| English polish after scientific sign-off | `nature-polishing` -> `$alphagenome-nature-methods-voice` | Preserve claims, citations, equations, symbols, defined terminology and honest scope. |

## Project-Local Voice Skill

For any manuscript wording, figure legend, title, abstract, cover letter or
editorial response in this project, invoke the local skill as
`$alphagenome-nature-methods-voice` and read:

`skills/alphagenome-nature-methods-voice/SKILL.md`

The skill incorporates the useful high-threshold editing principles from the
public `awesome-ai-research-writing` collection with Nature-style
claim-evidence discipline. It does not merely swap words to make prose sound
different. Its default is to retain clear human-authored sentences, remove
empty promotion and mechanical transitions, and flag unsupported claims rather
than smoothing them over.

Use it after evidence assembly, never to manufacture evidence. Run its local
audit script before a section is marked ready, then resolve each substantive
finding manually.

## Figure and Table Rules

- Plan figures by claim, not by source file. One figure must have one central
  conclusion and each panel must answer a distinct question.
- Keep a consistent method palette: full Model B as the focal color, Model A and
  from-scratch controls as neutral comparators, and clearly labeled external
  references as a separate restrained family.
- Plot paired effects for matched fold/seed comparisons. Report the five-fold,
  seed-nested bootstrap interval; do not visually inflate n with 241 tracks.
- Every quantitative display needs a source-data table, statistic/interval
  definition, sample unit, exclusion rule and versioned input manifest.
- Reuse no copyrighted figure or layout. Learn structural roles from the
  literature audit and draw original figures from project data.

## Draft and Release Gates

1. Before drafting: complete the terminology ledger and claim-evidence map.
   Planning documents may identify the evidence needed for external claims, but
   prose and figures must omit those claims until frozen artifacts exist.
2. Before a figure is rendered: complete its figure contract and select the
   plotting backend.
3. Before internal review: run the local voice audit, terminology audit,
   citation verification, statistics audit and figure/source-data audit.
4. Before submission: complete data/code availability, Reporting Summary,
   Software Submission Checklist, ethics/AI-use disclosure and clean-room
   reproduction checks.

## Git Hygiene

Check `git status --short --branch` before and after changes. Stage only files
inside `nature_methods_manuscript/` that were created for this workspace. Never
stage user-owned changes elsewhere, `template/vendor/`, raw data, checkpoints,
or rendered manuscript outputs. Use an English commit message. Inspect
ahead/behind and outgoing commits before pushing.
