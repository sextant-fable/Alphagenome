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

## 2026-07-13 - Full P3 and Standing GPU Authorization

- User direction: do not stop at a five-run endpoint; process the complete verified RNA-seq scope and automatically use available project GPUs in later phases.
- G1: approved `p3_full_rna_streaming_482` for 482 RNA-seq runs and 944.56 GiB compressed FASTQ in the existing `alphagenome` environment.
- G2: approved the audited 240-group candidate hierarchy as the reprocessing starting point; unknown-unit provided bigWigs remain non-formal.
- G3: approved full normalized/grouped outputs and the subsequent manifest-driven P4 loader/split; the 347 GB monolithic NPZ remains prohibited.
- G4: approved automatic available-device selection under the repository GPU 2/3 policy. No GPU is used in P3.
- G5: remains unapproved until a single checkpoint is locked after R6B.
- Storage design: use four 16-thread streaming workers and remove verified per-run FASTQ/BAM/bedGraph intermediates instead of retaining the approximately 3.8 TiB planning envelope.
- Metadata correction: added ENA per-file FASTQ MD5 and byte lists; R1 now passes 13 checks for an exact 482-run source manifest.

## 2026-07-14 - P3A First Attempt and Resume Correction

- Environment result: the existing `alphagenome` environment retained Python `3.12.13`, AlphaGenome importability, PyTorch `2.11.0+cu128`, and CUDA availability after installing STAR `2.7.11b`, samtools `1.23.1`, bedtools `2.31.1`, and UCSC bedGraphToBigWig `482`.
- Index result: the WBcel235 STAR index completed at `shared/reference_indexes/WBcel235_STAR_2.7.11b`.
- Attempt result: FAIL before alignment because ENA closed the first HTTPS FASTQ transfer repeatedly and the pilot curl command discarded prior partial bytes on retry.
- Correction: both pilot and full download paths now use `--continue-at -`, 20 transfer retries, and a two-second retry delay; a regression test confirms continuation from an existing `.part` file.
- Data integrity: no normalized output was produced, raw bigWigs were not modified, and the retained partial FASTQ is validated against the expected byte count and MD5 before promotion.

## 2026-07-14 - P3B Full-Processing Implementation

- Full scope: the controller now registers P3B execution and R3 review for all 482 RNA-seq runs; the three ChIP-seq runs remain excluded.
- Recovery: both partial and fully downloaded FASTQs survive a failed sample attempt and are checksum-validated before reuse.
- Final hierarchy: uniform reprocessing restores the three held source-reuse accessions as a distinct SRP310676 group, producing 482 one-to-one run memberships in 241 biological groups. They are not averaged with the canonical SRP278203 group because the source studies differ.
- Aggregation: technical runs within a biological unit are weighted by pre-normalization coverage mass, then biological units are averaged equally. Singleton group paths are audited relative symlinks to avoid duplicate large files.
- Post-reprocessing review: R3 verifies all normalized and grouped bigWigs, raw-data immutability, exact group context, ontology CURIEs, replicate QC, source-reuse resolution, cleanup, manifest hashes, and the long-run execution record.
- Status: implementation and 18-test suite pass; P3B has not started and remains conditional on R3A PASS.
