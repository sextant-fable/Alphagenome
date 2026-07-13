# C. elegans RNA-seq v2 Execution Log

Entries are append-only. Failed phases and corrected conclusions remain in the record.

## 2026-07-13 - Workflow Initialization

- Purpose: initialize the fail-closed v2 workflow and begin the metadata-only legacy freeze.
- Host: `hy8` (`HY-GPU`).
- Starting Git commit: `8b67b775c88002434d952d4ba0ebb86b0fbd3642` on `setup/agent-maintenance`.
- Scope: controller, phase reviews, metadata inventories, and documentation only.
- Excluded actions: no raw-data modification, dataset generation, large download, GPU work, training, or test evaluation.

## 2026-07-13 - P0/R0 Legacy Freeze

- Result: PASS. Frozen 20 original bigWigs, 11 grouped tracks, 188 intervals, 171 eligible checkpoints, selected model weights, references, environment, and source reports.
- Integrity: raw hashes agreed before and after; no protected raw, grouped, or run file changed.
- Transition: advanced automatically to P1.

## 2026-07-13 - P1/R1 Provenance Audit

- Result: PASS. All 485 expected accessions and local bigWigs were matched to ENA run, experiment, BioSample, BioProject, and study records.
- Source volume: 1,020,357,970,574 compressed FASTQ bytes reported by ENA.
- Normalization decision: A=0, B=0, C=485, D=0. Current bigWig units remain unknown; no sample was promoted to formal v2.
- Transition: advanced automatically to P2.

## 2026-07-13 - P2/R2 Manifest Review and Correction

- Initial structural review: PASS, then retained as an archived review.
- Correction: the pre-G2 scientific audit identified 45 additional runs sharing ENA experiments and 267 source-title/alias replicate labels. P2 was formally reopened rather than treating the initial review as final.
- Corrected result: PASS. The hierarchy now separates technical runs from biological units, yielding 240 candidate groups from 479 current RNA-seq signals.
- Exclusions: three ChIP-seq accessions and three secondary byte-identical current signals.
- Replicate QC: 645 pairs evaluated; 108 review-only flags in 17 groups; no automatic removal.
- Transition: advanced to P3, which remains behind G1/G2/G3 approvals.

## 2026-07-13 - Scoped Controller and P3A/R3A Preflight

- State migration: schema v1 `P3` became bounded `P3A`; all ungranted Boolean gates became ungranted scoped approval records.
- Phase separation: P3A pilot and P3B bulk processing are distinct; P6A pilot, P6B formal validation, and P6C one-time chromosome-X evaluation are distinct.
- Automation: `run --auto` continues only across PASS phases and stops fail-closed at missing scopes, failed reviews, or unregistered phase implementations.
- Defense in depth: the P3A wrapper independently verifies current phase and exact G1/G2/G3 scopes before environment installation, download, alignment, or output writes.
- Preflight result: Python compilation, command-plan inspection, missing-approval refusal, scope-mismatch refusal, and approval-boundary auto-stop passed without installing tools or writing large data.
- Current boundary: `P3A / APPROVAL_REQUIRED`; no P3A data action has run.
