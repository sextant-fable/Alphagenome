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

## 2026-07-14 - P3B NCBI SRA Transport Validation

- Reason: repeated ENA HTTPS range transfers retained bytes correctly but were too slow for the 944.56 GiB full scope. A 50 MiB read-only benchmark delivered the official NCBI SRA object at 7.84 MB/s while the corresponding ENA FASTQ endpoint delivered 2,359,691 bytes in 60 seconds.
- Tooling: installed the standalone NCBI SRA Toolkit `3.4.1` under ignored `shared/tools/`; the 89,143,120-byte official toolkit archive has SHA-256 `b950362c054765a4184af41947f022f040e94e964862017c0ecb0b0273db3596`. The project Python environment was not replaced.
- End-to-end check: downloaded the official full-quality `SRR941632` SRA object, matched SDL MD5 `57120e81e8e8b183e4fb78014d1237d9`, passed `vdb-validate`, and extracted 20,763,724 reads, exactly matching the ENA manifest. The independent validation copy was then removed.
- Full-run policy: P3B uses `ncbi_sra` transport by default and retains `ena_fastq` as an explicit fallback. Every run must match the SDL archive MD5, pass archive validation, and reproduce the manifest read count before STAR. Transport provenance and extracted FASTQ SHA-256 values are recorded per run.
- Count semantics: ENA `read_count` is treated as the SRA spot count. A single-end archive must extract one FASTQ record per spot and a paired-end archive must extract two; both the source spot count and extracted FASTQ record count are retained in each audit.
- Locked source inventory: all 482 RNA-seq accessions resolved to one public full-quality SRA object. The tracked transport manifest contains 722,781,740,735 bytes (673.14 GiB) and has SHA-256 `c04c455ebdc4b7b0daf3382533b09bbd519aa43cadb7fd77109e99aa178c71f6`.

## 2026-07-14 - Transfer Retry Semantics Correction

- Observation: curl continuation worked when P3A was restarted, but curl's internal retry logic truncated bytes received by the current curl process after a remote exit 18.
- Controlled action: the project P3A process was interrupted after one completed sample; the current `.part` and all completed intermediates were retained.
- Correction: each curl process now performs one resumable transfer attempt. Python restarts curl after a nonzero exit, so every new process resumes the retained `.part`; 200 attempts are allowed, and ten consecutive no-progress attempts fail closed.
- Regression coverage: tests verify both complete-file reuse and an exit-18 transfer that succeeds on the next process without `--retry`.
- Long-run hardening: the per-file attempt budget is 1,000 with a fail-closed threshold of 20 consecutive no-progress attempts. P3A demonstrated cumulative byte growth across repeated ENA exit-18 responses and completed the 1,772,340,596-byte SRR10882545 FASTQ without treating any equal-sized file as reusable.

## 2026-07-14 - P4 Dynamic Loader Implementation

- Split design: five leave-one-chromosome-out folds across chromosomes I-V, with chromosome X in a separately locked manifest and non-overlapping nearest-window-center evaluation cores.
- Loader interface: lazy per-process bigWig handles return 1 bp signal, 128 bp sum-pooled signal, track mask/strand, plus/minus gene masks, DNA, and evaluation core mask.
- Test embargo: chromosome X requires P6C state, scoped G5 approval, a locked checkpoint SHA-256, and an unused one-time test record.
- Benchmark design: one full 241-track CPU window plus direct pyBigWig, 128 bp pooling, single/multi-worker determinism, file descriptor, RAM, and I/O timing checks. No v2 NPZ is generated.
- Status: implementation tests pass; real R4 benchmark remains conditional on R3 PASS.

## 2026-07-14 - P5 Loss, Augmentation, and Model Contracts

- Model-space target: nonzero-mean normalization, RNA `x^0.75` compression, and soft clipping use the local AlphaGenome PyTorch formula; Poisson/multinomial receives the positive model-space prediction and scaled target directly.
- Resolution loss: both 1 bp and 128 bp sum-pooled targets use eight genomic segments, positional weight 5.0, and a strand-aware gene-body cross-track KL term weighted 0.1.
- Controlled baseline: the same model-space outputs can instead be unscaled and compared with log1p-MSE under identical data and seeds.
- Augmentation: deterministic per-epoch shifts are sampled within +/-1024 bp and chromosome bounds; reverse complement probability is 0.5, with DNA, 1/128 bp targets, masks, gene strands, track strands, and strand-paired channels transformed together.
- Model A: frozen AlphaGenome trunk with a new 241-track, two-resolution RNA head.
- Model B: four actual AlphaGenome organism embeddings are expanded from two to three rows; row 2 is initialized from the human/mouse mean and is the only embedding row receiving gradients. Final transformer modules receive LoRA. This is not an out-of-range `organism_index=2` shortcut.
- Model C: a residual dilated sequence baseline produces matched positive 1 bp and 128 bp outputs from scratch.
- Leakage control: per-track nonzero means are computed separately from each fold's four training chromosomes. An I-V development mean is retained only for the registered comparison; chromosome X is not read.
- Status: synthetic formula, gradient, shift, reverse-complement, strand, pooling, embedding, and baseline tests pass. Actual checkpoint/LoRA CPU inspection remains conditional on P4 PASS.

## 2026-07-14 - P6A Automatic GPU Smoke Implementation

- Resource snapshot: `hy8` has four A100 80 GB GPUs. GPU 0 had two unrelated processes using approximately 80 GB; GPUs 1, 2, and 3 were idle at inspection time. Project policy continues to exclude GPU 0/1 and prefer GPU 2, then GPU 3.
- Selection rule: query GPU UUIDs, memory, utilization, and active compute processes immediately before the phase; select only GPU 2/3 with at least 70,000 MiB free, utilization at most 10%, and no compute process. Never cancel or modify another process.
- Pilot scope: run A, B, and C sequentially on the selected device, with the same fold-1 data, seed `20260714`, 241 tracks, scaled paper loss, 131,072 bp diagnostic crop, and two optimizer steps.
- Integrity: each smoke records finite loss and gradient norms, peak CUDA memory, augmentation coordinates, exact interval/mean hashes, compact trainable-only checkpoint hash, reload verification, and an ignored log.
- Claim scope: P6A is environment and numerical validation only, not a model result.
- Status: implementation and 36-test suite pass; no v2 GPU job has run because P3-P5 reviews have not yet completed.

## 2026-07-14 - P6B Blocked-CV Matrix Registration

- Pre-registration: lock 30 runs from models A/B/C, training losses paper/log1p-MSE, and chromosome folds 1-5. Every run uses seed `20260714`, 2,000 optimizer steps, 131,072 bp subwindows, learning rate `1e-4`, the same shift/reverse-complement policy, and fold-training-only target means.
- Validation: evaluate both loss spaces for every checkpoint on all complete, non-overlapping 131,072 bp subwindows contained within each held-out chromosome's evaluation cores. Also record 128 bp per-track log-signal Pearson; chromosome X remains inaccessible.
- Selection: minimize five-fold mean paper loss, then five-fold mean log1p-MSE, then maximize five-fold mean per-track Pearson. This compares the two training objectives on identical data, seeds, steps, and evaluation metrics.
- Final development checkpoint: retrain only the selected model/loss configuration for 2,500 steps on chromosomes I-V with the development nonzero means, then bind its checkpoint SHA-256 into `final_test_lock.json` with `test_consumed=false`.
- GPU policy: wait for at least one idle allowed physical GPU and use up to GPUs 2 and 3 concurrently, with one job per GPU and no cancellation or oversubscription. Completed CV jobs are hash-validated and reusable after interruption.
- Controller: P6B/R6B is registered. R6B requires all 30 jobs, six five-fold aggregates, the exact locked selection rule, development retraining, checkpoint integrity, and continued chromosome-X embargo.
- Status: implementation and 45-test suite pass; no P6B GPU job has run because P3-P6A reviews have not yet completed.

## 2026-07-14 - P6C One-Time Final-Test Implementation

- Boundary: P6C requires both the standing GPU scope and a separate `r6c_single_chr_x_test` G5 approval that can be granted only after P6B locks one development checkpoint. G5 remains unapproved.
- Concurrency: create one exclusive final-test claim bound to the checkpoint SHA-256 and an idle physical GPU 2/3. A pre-existing claim or consumed lock refuses another entry.
- Consumption semantics: the loader first verifies P6C state, G5 scope, checkpoint hash, and `test_consumed=false`. The evaluator then atomically sets `test_consumed=true` before the first chromosome-X BigWig read. A later evaluation failure remains consumed and is recorded as `failed_after_consumption`.
- Evaluation: use all complete core-only 131,072 bp chromosome-X subwindows and report paper loss, log1p-MSE, and 128 bp per-track Pearson for the locked checkpoint only.
- Review: R6C verifies the single claim, GPU policy, exact locked checkpoint, chromosome-X interval hash, finite metrics, permanent consumption, final report hash, G5 scope, and disclosure that chr X had prior legacy-project exposure.
- Status: implementation and 46-test suite pass; no chromosome-X data was read by this implementation work.
