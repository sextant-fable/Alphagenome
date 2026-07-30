# C. elegans RNA-seq v2 Execution Log

Entries are append-only. Failed phases and corrected conclusions remain in the record.

## 2026-07-30 - P8 B-noLoRA Fold-1 Ablation Result

- Execution: P8 completed on HY-GPU physical GPU 3 from `2026-07-30T11:54:24+00:00` to `2026-07-30T12:31:54+00:00` under `G4:p8_b_no_lora_fold1`, using execution commit `923f58190727f79ae3c8d2f3ce6a36be4b6638e2`.
- Controlled difference: B-noLoRA preserved B's newly added trainable third-row C. elegans organism embedding and 241-track dual 1 bp/128 bp RNA head. It removed all LoRA modules (`lora_enabled=false`, empty target list) and kept original trunk policy `worm_embeddings_only`. The candidate trained `1,130,308` parameters, `145,920` fewer than matched B/paper (`1,276,228`).
- Matched development contract: fold `1`, seed `20260714`, paper loss, 2,000 steps, 131,072 bp contexts, learning rate `1e-4`, gene weight `0.1`, standard shift/reverse-complement augmentation, hidden channels `64`, and the fold-1 train-only mean. All 83 registered development cores were evaluated, covering `10,878,976` bases at coverage `1.0`.
- Result: B-noLoRA obtained primary biological score `0.509567111366191`, gene-exon coverage Pearson `0.443595944040327`, and mean 128 bp per-track Pearson `0.575538278692054`. Exact matched B/paper was `0.602161323942477`, `0.566838838705920`, and `0.637483809179034`; candidate-minus-baseline deltas were `-0.092594212576286`, `-0.123242894665593`, and `-0.061945530486979`.
- Interpretation: in this one controlled fold/seed comparison, removing LoRA substantially harms every primary biological metric, so the current evidence does not support dropping LoRA. It remains a single-fold development ablation, not a replacement for the P6B selection or a general cross-fold estimate.
- Test preservation and review: `locked_test_block_signal_reads=0`, final-test access was prohibited, and R8 passed all seven checks, including unchanged P6C lock/report hashes. The machine-readable records are `p8_b_no_lora_fold1_execution.json`, `p8_b_no_lora_fold1_comparison.json`, and `audits/P8/review.json`.

## 2026-07-30 - P8 B-noLoRA Fold-1 Ablation Registration

- User scope: one matched fold-1 B-noLoRA run retaining the third C. elegans organism embedding and dual RNA head while removing LoRA. This is development-only, does not retune the formal selection, and may not access the consumed final-test blocks.
- Locked contract: fold `1`, seed `20260714`, paper loss, 2,000 steps, 131,072 bp contexts, standard shift/reverse-complement augmentation, and fold-1 train-only means. The exact B/paper fold-1 baseline checkpoint/validation artifact, manifests, means, split registry, and completed P6C lock/report are SHA-256 bound in `p8_b_no_lora_fold1_spec.json`.
- Architecture audit: `B_no_lora` expands and trains only the new third-row C. elegans organism embeddings, retains the 241-track direct dual-resolution RNA head, and has no LoRA modules. The trained run record and R8 review require `lora_enabled=false`, an empty LoRA target list, and `trunk_policy=worm_embeddings_only`.
- Controller boundary: P8 can start only after P7 is complete and requires the distinct `G4:p8_b_no_lora_fold1` scope. Its commands reject `--final-test`; R8 verifies that the P6C lock/report hashes remain unchanged.
- Status: registered and pending execution; no P8 signal read or model result yet.

## 2026-07-30 - P7 Legacy 1bp-B Architecture Port Result

- Execution: P7 completed on HY-GPU physical GPU 2 from `2026-07-30T07:57:40+00:00` to `2026-07-30T08:33:08+00:00` under `G4:p7_single_legacy_residual_fold1`. The implementation used commit `b2501f1`; the controller/runner execution commit was `f5e375746b249bd7df25f60e51bca0b1dbef8111`.
- Data and topology: all 241 current uniformly reprocessed RNA-seq tracks were used. A new 241-track 128 bp Conv5 base head was trained for 2,000 steps, frozen (SHA-256 `7b8300b2230fed406ab930e52eed37a5a3eacb07ec72e1f64548c26552f327f2`), then its frozen 1 bp interpolation received the legacy 128-channel-bottleneck Conv5 residual trained for 1,500 steps. Historical 11-track weights were not dimension-compatible with the new targets.
- Development result: across all 83 complete fold-1 development cores (10,878,976 bases; coverage `1.0`), D/paper seed `20260714` achieved primary biological score `0.602707560944138`, gene-exon coverage Pearson `0.564589786166844`, and mean 128 bp per-track Pearson `0.640825335721432`.
- Matched B/paper comparator: `0.602161323942477`, `0.566838838705920`, and `0.637483809179034`, respectively. D minus B was `+0.000546237001661` for the composite, `-0.002249052539076` for gene-exon Pearson, and `+0.003341526542398` for 128 bp Pearson.
- Interpretation: this one-fold/one-seed architecture check is effectively mixed/near-neutral, not a new selection. It cannot supersede P6B's formal selection or P6C's already-consumed final evaluation.
- Test preservation: `locked_test_block_signal_reads=0`, final-test access was prohibited, and R7 verified unchanged P6C lock/report hashes. P7 did evaluate the six-chromosome development split, including its registered development X blocks, but not locked final-test blocks.
- Review: R7 passed all required-output, locked-input, training-contract, coverage, baseline-comparison, final-test-preservation, and controller-scope checks. Machine-readable records are `p7_legacy_residual_fold1_execution.json`, `p7_legacy_residual_fold1_comparison.json`, and `audits/P7/review.json`.

## 2026-07-30 - P7 Legacy 1bp-B Architecture Port Registration

- User scope: one fold-1 run of the pre-new-data 1bp-B structure on all 241 current v2 tracks, with a result reported against the matched existing B/paper fold-1 baseline.
- Architecture fidelity: the 128 bp Conv5 base head is trained on the new target space then frozen; the 1 bp head uses the legacy 128-channel bottleneck plus Conv5 residual correction. The old 11-track checkpoint cannot be loaded because its output dimensions do not match 241 v2 tracks.
- Locked development contract: fold `1`, seed `20260714`, paper loss, 131,072 bp windows, fold-1 train-only means, 2,000 base-head steps followed by 1,500 residual steps. The B/paper baseline checkpoint, baseline validation report, manifests, means, split registry, and completed P6C lock/report are SHA-256 bound in `p7_legacy_residual_fold1_spec.json`.
- Final-test preservation: P7 prohibits `--final-test`, uses only `fold_1/valid.tsv`, and R7 verifies that the already-consumed P6C lock/report hashes are unchanged. It may not alter the formal selection or final-test conclusion.
- Controller boundary: the new P7 entry point is allowed only from completed P6C and requires the distinct `G4:p7_single_legacy_residual_fold1` approval. The registered implementation commit is `b2501f1`.
- Status: registered and pending execution; no P7 signal read or model result yet.

## 2026-07-22 - P6B External Interruption Recovery

- Observation: the corrected six-chromosome P6B controller stopped making progress after `2026-07-22T04:51:55+00:00`. At audit time, no AlphaGenome process remained on physical GPU 2 or 3, while the live execution record was stale at `running`, with eight completed jobs, zero recorded failures, and zero locked-test reads.
- Preserved evidence: the execution JSON, both incomplete B/fold-2 logs, and the completed-training B/paper partial output were copied to timestamped interruption directories with SHA-256 values in `alphagenome_custom/metadata/v2/p6b_external_interrupt_20260722T045155Z/audit.json`.
- Disposition: eight complete jobs remain eligible for the existing hash-verified reuse path. The B/paper job stopped at validation 21/83 and the B/log1p-MSE job at training step 1345/2000; neither may count as a result, and both restart from the beginning under the identical locked specification.
- Gate status: G4 remains the same scoped GPU 2/3 approval. G5 remains unapproved, and locked-test signal reads remain zero.
- Recovery start: the first controller retry failed before GPU selection because managed-sandbox `nvidia-smi` returned exit 9. After an append-only reopen, the identical host-permission command selected physical GPUs 2/3, hash-verified and reused all eight complete jobs, and restarted both incomplete jobs from step 0.
- Durable execution: both restarted B/fold-2 jobs completed with full six-chromosome coverage, bringing the reusable total to ten. The controller was then stopped deliberately and migrated to detached tmux session `alphagenome_v2_p6b`; C/paper and C/log1p-MSE partial attempts at steps 596 and 81 were archived and restart from step 0. The machine-readable migration record is `p6b_tmux_migration_20260722T064354Z/audit.json`.

## 2026-07-22 - P6B Core-Coverage Failure and v2 Correction

- Controlled stop: the first completed `six_chromosome_blocks_v1` formal job measured `validation_core_coverage_fraction=0.8435322914705329`. The controller was interrupted before the matrix could produce a comparative result; one completed job and two partial jobs are retained as invalid/incomplete evidence. Locked-test signal reads remained zero.
- Root cause: effective blocks were multiples of `131,072 bp`, but nearest-window core boundaries were not. The evaluator accepted only complete subwindows inside each irregular core, silently omitting 13 of the 83 intended metric cores in fold 1.
- Corrected split: `six_chromosome_blocks_v2` reuses the exact chromosome block assignments, train windows, seed, buffers, data, models, losses, folds, seeds, and budgets. Only valid/test manifest geometry changes: every effective block is tiled by complete, contiguous `131,072 bp` cores, and each core receives a containing `1,048,576 bp` context row with a 128-bp-aligned relative crop offset.
- Dry-run evidence: every fold and the locked test contain 83 metric rows/subwindows over `10,878,976 bp`, for coverage `1.0`; all contain `I`, `II`, `III`, `IV`, `V`, and `X`. The v2 P6B specification rejects every v1 job and uses independent `corefix_20260722` run/log paths.
- Review amendment: R4 now verifies exact core length, count, continuity, context containment, and relative 128-bp alignment. R6B retains the >=0.999 coverage gate and separately verifies that the invalid v1 attempt was archived and excluded.

## 2026-07-21 - Six-Chromosome P4-P6A Execution

- Migration: the original chromosome-holdout intervals, P4/P5/P6A/P6B reviews, overwrite-prone metadata, and unconsumed final-test lock were copied to dedicated archive directories with SHA-256 evidence. The old lock is superseded, not deleted; G5 remained unapproved.
- R4 PASS: every train and validation fold contains `I`, `II`, `III`, `IV`, `V`, and `X`; the locked final test contains one block assignment per chromosome. Each fold has 64 train windows and 16 validation windows. The final test has 16 windows, and adjacent effective blocks have a `1,050,624 bp` exclusion buffer. Default loader construction refused the test manifest; signal reads were zero.
- R5 PASS: 241 nonzero means were recomputed from 30 development CV blocks and 24 train blocks per fold, with zero locked-test reads. The first implementation reread the same blocks for each fold; commit `561f83c` changed this to one per-block sum/count scan followed by identical fold aggregation, covered by 79 passing tests.
- R6A PASS: A, B, and C each completed the registered two-step smoke on physical GPU 2 with finite losses/gradients and checkpoint/log verification. GPU 3 was occupied by an unrelated CellUNetr process and was not touched. Smoke evidence is environmental only, not a model result.
- Operational evidence: the first P4 benchmark was interrupted after the managed sandbox denied DataLoader socket creation; an append-only reopen records the cause, and the identical host-permission controller retry passed. No monolithic NPZ was created.

## 2026-07-21 - Six-Chromosome Split Revision Preregistration

- User requirement: every formal train, validation, and test partition must contain WBcel235 blocks from `I`, `II`, `III`, `IV`, `V`, and sex chromosome `X`; `X` is not chromosome 10.
- Disposition: preserve the chromosome-holdout P4-P6B evidence, but supersede its I-V leave-one-chromosome-out promotion decision, development checkpoint, and X-only final-test lock before any G5 consumption. P1-P3 provenance, 482 uniform RNA-seq reprocessings, and 241 grouped tracks remain reusable if their hashes match.
- Preregistered split: on each chromosome, deterministically permute five CV block labels and one locked-test label with seed `20260721`. Each fold trains on four CV blocks and validates on the fifth; all six chromosomes therefore occur in both roles. Development training uses all five CV blocks, and final test uses one locked block per chromosome.
- Leakage boundary: adjacent effective blocks are separated by `1,050,624 bp`, equal to the `1,048,576 bp` context window plus two `1,024 bp` augmentation margins. Effective block sizes are multiples of `131,072 bp`, so valid/test core coverage is complete without overlapping metric subwindows. Training shifts are constrained to remain inside their registered block.
- Expected geometry: 36 effective blocks total; each fold has 64 one-megabase training windows and 16 validation windows, with `10,878,976` eligible validation bases. The final test has 16 windows across six chromosomes and the same eligible-base total.
- Gates: P4 uses the revised `G3:p4_six_chromosome_block_split` scope. The standing GPU 2/3 G4 remains applicable to P6A/P6B. The old X-only G5 is invalid; the revised final scope is `G5:r6c_single_six_chromosome_block_test` and remains unapproved until revised R6B PASS.
- Preregistered P6B: A/B/C x paper/log1p-MSE x five folds x three seeds (90 formal jobs), plus 15 five-fold single-seed augmentation, gene-loss, and development-pool-mean ablations. Training budgets and biological promotion metrics remain identical to the amended chromosome-holdout comparison.

## 2026-07-21 - P6C Pre-G5 Gate Audit

- Purpose: independently verify the one-time chromosome-X entry point after the amended R6B lock and before requesting G5. No chromosome-X BigWig was opened.
- Findings corrected: R6C still referenced the superseded `p6b_selection.json` instead of the selection path bound into the amended lock; the evaluator did not require an exclusive claim or unique execution ID; and direct loader use could read chromosome X after G5 without first consuming the one-time lock.
- Gate hardening: P6C now validates the locked checkpoint, amended selection, preregistered test intervals, development means, and split-registry SHA-256 values before creating a claim. The claim is created with `O_EXCL`, records one execution ID and physical GPU 2/3, and is serialized during consumption with a file lock.
- Read boundary: a chromosome-X dataset requires the matching claim at construction and requires the same execution ID in `test_consumed=true / test_status=running` state before every signal read. A concurrent or repeated consumer is refused before signal access.
- Failure semantics: a failure after first signal access remains permanently consumed and becomes `failed_after_consumption`; a successful run binds the metric report, final synthesis, human-readable report, checkpoint, claim, GPU, and execution ID by SHA-256.
- Final synthesis: successful P6C execution will automatically summarize the 485-source lineage, 482 uniformly reprocessed RNA-seq runs, 241 grouped tracks, blocked split, six A/B/C-loss configurations, three ablations, registered failures, selected checkpoint, final metrics, limitations, and exact reproduction command.
- Verification: 68 v2 unit tests passed; all scripts and v2 tests compiled; `git diff --check` passed. A metadata-only summary dry-run returned 485 samples, 241 tracks, six matrix rows, three ablations, and zero registered P6B failures. Locked hashes matched checkpoint `b7bdd540...25b`, amended selection `ff434c21...6cae`, test intervals `5f398d18...016`, means `3734bdf7...beb`, and split registry `b6abf1d8...875`.
- Current boundary: G5 is still unapproved, `test_consumed=false`, no final-test claim exists, and the controller remains `P6C / APPROVAL_REQUIRED`.

## 2026-07-21 - P6B Amendment Completion and R6B PASS

- Purpose: complete the pre-registered P6B amendment required by the independent pre-G5 audit, without reading chromosome X.
- Command: `/home/zelinli6/miniconda3/envs/alphagenome/bin/python -m scripts.v2_phase_controller run --phase P6B --auto`.
- Execution window: `2026-07-20T09:26:03+00:00` through `2026-07-21T05:10:10+00:00` on HY-GPU, using physical GPUs 2 and 3 without using GPUs 0 or 1. The selected development retrain used physical GPU 2.
- Code version: the run started from `3769d54` (`Reopen P6B with complete validation matrix`). Later commits through `79463c4` changed only `README.md` and `AGENTS.md`, not execution code.
- Inputs and split boundary: 241 uniformly reprocessed grouped RNA-seq tracks on WBcel235; blocked development folds I-V; chromosome X remained embargoed.
- Formal matrix: 90/90 jobs completed for A/B/C x paper/log1p-MSE x five folds x three seeds (`20260714`, `20260715`, `20260716`). The original 30 jobs were retained and backfilled with the amended metrics; 60 additional seed jobs completed. Registered amendment failures: none.
- Ablations: 15/15 jobs completed for no augmentation, no gene loss, and the whole-I-V mean comparator, each across five folds at seed `20260714`. These are single-seed directional checks, not significance tests.
- Selection: the pre-registered biological score selected B/paper with mean score `0.6322474761`, gene-exon coverage Pearson `0.6042942527`, and 128-bp per-track Pearson `0.6602006995`. B/log1p-MSE scored `0.6190121491`, `0.5539517632`, and `0.6840725349`, respectively.
- Paired interpretation: B/paper beat B/log1p-MSE on the primary score in 13/15 seed-fold pairs and on gene-exon Pearson in 15/15, but lost on 128-bp Pearson and log1p-MSE in 15/15. It is the winner under the locked composite objective, not across every metric.
- Stability: B/paper seed-level primary means were `0.632253500`, `0.629022262`, and `0.635466666` (sample SD `0.003222206`); fold-level means had sample SD `0.029768863`, so fold variation dominated seed variation.
- Ablation interpretation: relative to the exact same-seed formal baselines, disabling augmentation changed the mean primary score by `+0.001402370`, disabling gene loss by `+0.000336154`, and the whole-I-V mean changed A/paper by `-0.000343856`. Effects are small and single-seed; they do not establish that augmentation or gene loss should be removed.
- Final development retrain: B/paper, seed `20260717`, 2,500 steps on development chromosomes I-V. Checkpoint: `runs/v2_p6b_amendment_20260720/development_selected/checkpoint.pt`; SHA-256: `b7bdd54066b102026edd29dbe6dc23bd5b59566750b7f7ba94728821ec31525b`.
- Auditable outputs: `alphagenome_custom/metadata/v2/p6b_amendment_execution.json`, `p6b_amendment_cv_results.tsv`, `p6b_amendment_ablation_results.tsv`, `p6b_amendment_selection.json`, and `audits/P6B/review.json`.
- Review result: R6B-amended PASS on all ten checks. The prior P6B lock remains preserved as superseded; legacy v1 is not treated as a formal label-scale ablation.
- Test boundary: `chromosome_x_reads=0`, `test_consumed=false`, and prior exploratory chromosome-X exposure is disclosed in the lock. The controller stopped at `P6C / APPROVAL_REQUIRED`; only explicit `G5:r6c_single_chr_x_test` approval may consume the one-time final test.

## 2026-07-20 - P6B Pre-G5 Independent Audit and Reopen

- Original execution result retained: all 30 registered single-seed A/B/C x paper/log1p-MSE x five-fold jobs completed without a run failure, and the original implementation selected B/paper by minimum mean paper loss.
- Audit correction: that R6B review tested only the reduced local specification, not the complete locked research objective. It omitted formal multi-seed comparison, augmentation and gene-loss ablations, a defensible legacy-vs-v2 disposition, gene-level exon coverage Pearson, Spearman, high-signal calibration/error, gene-body/exon metrics, and local-gradient metrics. Its paper-loss-first promotion rule also conflicted with the objective's primary biological metrics.
- State action: P6B was formally reopened. The original result, selection, checkpoint, and review remain preserved as provisional evidence; the old final-test lock is marked superseded.
- Test embargo: chromosome X was not read. G5 remains unapproved, and both the loader and P6C reject a superseded lock.

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
- Audit schema lock: before group finalization, all 482 sample audits are migrated to schema v2. ENA FASTQ URL/MD5/bytes are explicitly labeled as source-reference metadata; NCBI SRA URL/MD5/SHA-256/bytes identify the downloaded archive; extracted FASTQ SHA-256/bytes/record count identify the STAR input. Ambiguous legacy `source_fastq_*` checksum fields are removed before R3.

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
