# Part 7: Supplementary Information, Release and Submission Plan

## Supplementary Package

| Deliverable | Required content |
| --- | --- |
| Supplementary Methods | Expanded processing, model, baseline, statistics, variant and reproducibility details too granular for Online Methods. |
| Supplementary Table 1 | Complete 241-track metadata and biological grouping manifest. |
| Supplementary Table 2 | All model/baseline hyperparameters, parameter counts, compute and tuning permissions. |
| Supplementary Table 3 | External dataset accession, source-isolation, condition and leakage-audit inventory. |
| Supplementary Table 4 | Full per-study, per-condition and per-track benchmark results. |
| Supplementary Table 5 | Variant truth, matched-null criteria, coordinate checks and endpoint outcomes. |
| Supplementary Table 6 | Software release, environment lock, checkpoint and reproduction hashes. |
| Supplementary Figures | Expanded loci, calibration, sensitivity, failed modes and all visual evidence not essential to the main claim ladder. |
| Source Data | One workbook or structured archive per figure/table, following the Part 6 contract. |

## Release Architecture

Before submission, publish or archive the following with stable identifiers:

1. Source code with an open license, tagged commit and reproducible environment
   definition.
2. Data manifest with raw accession IDs, grouping logic, reference files,
   checksums and a lawful access route for each item.
3. Frozen split and benchmark manifests, including a machine-readable
   source-isolation and receptive-field overlap audit.
4. Checkpoint/weight release terms, model card, inference command and expected
   input/output schemas.
5. Variant scoring specification, examples, test VCF, test reference region and
   expected checksums.
6. A clean-machine run log that recreates one benchmark prediction, one figure
   source-data table and a source-data hash without access to a private
   development tree.

## Submission Components

- Main manuscript with figures and Extended Data in the journal-requested form.
- Cover letter focused on the field-level method problem, external benchmark,
  biological application and released reuse path.
- Nature Portfolio Reporting Summary.
- Software Submission Checklist for code central to the main claims.
- Data Availability and Code Availability statements verified with
  `nature-data`.
- Ethics declarations, competing interests, funding and author contributions.
- AI-use statement describing assistive use, human validation and no AI
  authorship, as required by the journal policy at submission time.
- Suggested/excluded reviewer list only after conflict and expertise checks.

## Pre-Submission Review Sequence

1. Use `nature-statistics` to audit every effect, n, interval and caption.
2. Use `nature-citation` plus primary-source checks for every literature claim.
3. Use `nature-data` for availability and repository metadata.
4. Use `$alphagenome-nature-methods-voice` on each approved section; preserve
   terminology and all numerical evidence.
5. Run independent `nature-reviewer` simulations in frozen contexts, synthesize
   only after all reports are frozen, and resolve evidence gaps rather than
   arguing around them.
6. Run figure visual and source-data QA at final physical size.

## Nature Methods Formatting Check

Recheck the current official pages immediately before submission. The current
plan uses a 150-word abstract, approximately 3,000 main-text words and at most
six main display items, but current journal instructions override this plan.
