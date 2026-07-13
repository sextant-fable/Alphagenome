# Experiment Log

Use this file for concise, auditable records of project runs. Keep smoke tests, sanity checks, preprocessing, training, and evaluation clearly separated.

Do not delete failed runs. Append corrections or follow-up notes instead.

## Template

````md
## YYYY-MM-DD - Short Run Title

- Run type: smoke test | sanity check | preprocessing | training | evaluation | analysis
- Purpose:
- Git commit:
- Branch:
- Host:
- Slurm job ID:
- Slurm request:
- Environment:
- Command:

```bash

```

- Input data:
- Output path:
- Result summary:
- Verification:
- Failures or warnings:
- Next actions:
- Claim status: verified | unverified | inferred | abstract-only | secondary-source-only
````

## Runs

## 2026-07-14 - RNA-seq v2 P3A Download Failure

- Run type: preprocessing
- Purpose: Run the five-sample R3A technical checkpoint before automatically starting the authorized 482-run RNA-seq reprocessing scope.
- Git commit: `38271ca`
- Branch: `setup/agent-maintenance`
- Host: `HY-GPU` (`hy8`)
- Slurm job ID: Not applicable; HY-GPU is a non-Slurm server.
- Slurm request: Not applicable.
- Environment: Existing Conda environment `alphagenome`; STAR `2.7.11b`; samtools `1.23.1`; bedtools `2.31.1`; UCSC bedGraphToBigWig `482`; 16 CPU threads; no GPU used.
- Command:

```bash
/home/zelinli6/miniconda3/envs/alphagenome/bin/python \
  scripts/v2_phase_controller.py run --phase P3A
```

- Input data: `alphagenome_custom/metadata/v2/p3_pilot_sources.tsv`; five ENA RNA-seq runs totaling 6,221,459,671 compressed FASTQ bytes; WBcel235 reference and Ensembl release 115 GTF.
- Output path: `shared/source_reads/v2/pilot_20260713/`; intended normalized output directory `alphagenome_custom/tracks/rna_seq_v2_normalized_pilot/`.
- Result summary: Failed during the first FASTQ download before alignment. The STAR WBcel235 index completed successfully. No normalized bigWig was produced.
- Verification: The controller recorded `P3A / FAIL` at `2026-07-13T16:45:56+00:00`. curl exit code `18` showed repeated premature remote connection closure. The downloader discarded partial content on each internal retry because resume mode was absent.
- Failures or warnings: This is a transfer-resilience defect, not an RNA-seq, STAR, environment, CPU, or GPU failure. A 175,746-byte `.part` file remains for controlled resume.
- Next actions: Add and test curl continuation with a larger retry budget, reopen P3A with the failure reason retained, and resume the same source file.
- Claim status: verified

## 2026-07-14 - RNA-seq v2 P3A Controlled Transfer Restart

- Run type: preprocessing
- Purpose: Resume P3A after adding curl continuation, then verify that continuation also survives curl's own retry behavior before applying it to the 482-run scope.
- Git commit: `ec40470` at run start; P3B implementation commit `1095393` was created while the process ran.
- Branch: `setup/agent-maintenance`
- Host: `HY-GPU` (`hy8`)
- Environment: Existing Conda environment `alphagenome`; 16 CPU threads; no GPU used.
- Command: `/home/zelinli6/miniconda3/envs/alphagenome/bin/python scripts/v2_phase_controller.py run --phase P3A`
- Input data: Five-run P3A source manifest; retained partial ENA FASTQ from the first attempt.
- Output path: `shared/source_reads/v2/pilot_20260713/` and `alphagenome_custom/tracks/rna_seq_v2_normalized_pilot/`.
- Result summary: Controlled interruption after one of five samples completed alignment, BAM validation, coverage generation, and normalized bigWig output. The second FASTQ was incomplete; no R3A review ran.
- Verification: `SRR7443583.bw` was generated and its intermediates remain available for deterministic reuse. The second source remained a `.part` file. No raw provided bigWig was modified.
- Failures or warnings: curl `--continue-at -` resumed across separate invocations, but curl's internal `--retry` discarded 473,239,552 bytes and later 582,287,746 bytes from its current invocation. The process was interrupted intentionally because this behavior cannot safely handle the largest 11.77 GiB run.
- Next actions: Remove curl internal retry, restart curl from Python after each nonzero exit, verify byte growth across process boundaries, and resume P3A from the retained partial file.
- Claim status: verified

## 2026-07-14 - RNA-seq v2 P3A Controlled SRA Transport Switch

- Run type: preprocessing
- Purpose: Stop spending unbounded time on the fifth pilot accession's unstable ENA endpoint and resume the same biological run from a locked, full-quality NCBI SRA object.
- Git commit: `78c681b` for the fallback implementation and review contract.
- Branch: `setup/agent-maintenance`
- Host: `HY-GPU` (`hy8`); 16 CPU threads; no GPU used.
- Command: controlled interrupt of the project-owned P3A controller, followed by `python scripts/v2_phase_controller.py reopen --phase P3A ...` and the registered P3A command.
- Input data: Four completed ENA FASTQ pilot accessions; `SRR36719198` locked NCBI SRA URL, size, and MD5 from `p3_ncbi_sra_sources.tsv`; WBcel235 reference.
- Output path: unchanged pilot work and output directories under `shared/source_reads/v2/pilot_20260713/` and `alphagenome_custom/tracks/rna_seq_v2_normalized_pilot/`.
- Result summary: The interrupted ENA attempt had retained 140,038,604 bytes for `SRR36719198` but continued to close after short ranges. The partial file is an incomplete transport artifact, not a valid source input. Four normalized pilot bigWigs were already complete and retained.
- Verification: The alternate full-quality SRA object was pre-resolved for the same run accession, and an independent `SRR941632` end-to-end transport validation matched the official archive MD5, passed `vdb-validate`, and reproduced the exact manifest spot/read count.
- Failures or warnings: P3A will use mixed transport packaging for this technical checkpoint: four ENA FASTQs and one NCBI SRA. R3A explicitly permits only `SRR36719198` as the SRA fallback and independently verifies archive and extracted FASTQ hashes. This does not make the provided unknown-unit bigWigs formal labels.
- Next actions: Reopen P3A, complete the fifth alignment/coverage, run R3A, then start the locked 482-run P3B SRA workflow on PASS.
- Claim status: verified

## 2026-07-14 - RNA-seq v2 P3A/R3A Technical Checkpoint PASS

- Run type: preprocessing and technical review
- Purpose: Complete the five-accession checkpoint before the authorized full 482-run processing phase.
- Git commit: `78c681b` for executed pilot/review code; execution record committed immediately after PASS.
- Branch: `setup/agent-maintenance`
- Host: `HY-GPU` (`hy8`); 16 CPU threads; no GPU used.
- Input data: four checksum-verified ENA FASTQ accessions and one locked full-quality NCBI SRA fallback (`SRR36719198`); WBcel235 reference and Ensembl release 115 GTF.
- Output path: `alphagenome_custom/tracks/rna_seq_v2_normalized_pilot/` and retained diagnostics under `shared/source_reads/v2/pilot_20260713/`.
- Result summary: PASS. Five normalized bigWigs were produced at total signal 100,000,000 under the primary unique, spliced, unstranded coverage policy.
- Verification: R3A passed all 13 checks, including source/archive/FASTQ integrity, BAM hashes, raw provided-bigWig immutability, exact five-run membership, normalization, output read-back, path containment, no partial files, resource measurement, and five diagnostic provided-signal comparisons. The final resumed pass took 483.58 seconds and observed 25,792,977,581 managed bytes.
- Failures or warnings: Provided bigWigs retain unknown units and comparisons against them are diagnostic only. P3A outputs are not formal v2 labels.
- Next actions: Start P3B for all 482 locked RNA-seq accessions using the validated NCBI SRA transport and per-sample streaming cleanup.
- Claim status: verified

## 2026-07-14 - RNA-seq v2 P3B Full 482-Run Reprocessing

- Run type: preprocessing
- Purpose: Uniformly regenerate formal comparable coverage for all 482 verified RNA-seq accessions, then build the reviewed 241-track biological hierarchy.
- Git commit: `6630e111609eb68f9876cd74767d1a731c6e5a4c` at run start.
- Branch: `setup/agent-maintenance`
- Host: `HY-GPU` (`hy8`); non-Slurm CPU preprocessing; no GPU requested.
- Resources: four concurrent samples, 16 threads per sample; 112 logical CPUs, approximately 945 GiB available RAM and 2.0 TiB available disk at launch.
- Command: `/home/zelinli6/miniconda3/envs/alphagenome/bin/python scripts/v2_phase_controller.py run --phase P3B`.
- Input data: 482 locked NCBI full-quality SRA objects totaling 722,781,740,735 bytes (673.14 GiB), manifest SHA-256 `c04c455ebdc4b7b0daf3382533b09bbd519aa43cadb7fd77109e99aa178c71f6`; WBcel235 reference; 482 immutable provided bigWigs used only for identity/immutability evidence.
- Output path: normalized run tracks under `alphagenome_custom/tracks/rna_seq_v2_normalized/`; grouped tracks under `alphagenome_custom/tracks/rna_seq_v2_grouped/`; streaming work under `shared/source_reads/v2/full_streaming/`; log under `logs/v2_p3b_20260714/p3b_full.log`.
- Result summary: Running.
- Verification: Each sample must match archive MD5, pass `vdb-validate`, reproduce the layout-aware FASTQ record count, pass STAR/BAM/bigWig checks, and reach total signal 100,000,000 before its SRA/FASTQ/BAM/bedGraph intermediates are deleted. R3 must pass before P4.
- Failures or warnings: This is a long CPU/network run. A failed sample stops admission of new samples and retains its work directory for diagnosis; completed normalized outputs are resumable. Existing unknown-unit bigWigs are not reused as labels.
- Next actions: Monitor progress, append the final result without removing failures, run R3, and advance to P4 only on PASS.
- Claim status: unverified

## 2026-07-14 - RNA-seq v2 P3B SRA Extraction Failure and Audited Fallback

- Run type: preprocessing failure analysis
- Purpose: Record the first fail-closed P3B attempt and correct the transport path for aligned SRA archives that require unavailable external reference objects.
- Git commit: `6630e111609eb68f9876cd74767d1a731c6e5a4c` at failed-run start; fallback implementation committed before restart.
- Branch: `setup/agent-maintenance`
- Host: `HY-GPU` (`hy8`); non-Slurm CPU preprocessing; no GPU requested.
- Command: `/home/zelinli6/miniconda3/envs/alphagenome/bin/python scripts/v2_phase_controller.py run --phase P3B`.
- Input data: The same locked 482-run SRA and ENA source manifests used by P3B.
- Output path: Existing normalized outputs and audits under `alphagenome_custom/tracks/rna_seq_v2_normalized/`; retained failed work under `shared/source_reads/v2/full_streaming/`; log under `logs/v2_p3b_20260714/p3b_full.log`.
- Result summary: Failed closed after 73/482 runs completed. `SRR7443600` and `SRR7443602` returned `fasterq-dump` exit code 3 because their aligned archives referenced unavailable `BX28460*.4/.5` RefSeq objects. No completed output was invalidated.
- Verification: Both SRA archives matched locked NCBI SDL MD5 values and passed `vdb-validate`; the failure occurred during read reconstruction. Completed outputs retained their output hashes, immutable input-bigWig checks, STAR metrics, normalization totals, and cleanup status.
- Failures or warnings: An aligned SRA can be archive-valid yet not be self-contained for local FASTQ reconstruction. The correction falls back only after an observed SRA extraction failure to the same run accession's locked ENA FASTQ files, with exact ENA MD5/local SHA-256 and STAR spot-count verification. Transport fallback does not relax biological grouping or output QC.
- Next actions: Reopen P3B, reuse the 73 completed audited outputs, verify the first ENA fallback end to end, and continue the remaining full scope.
- Claim status: verified

## 2026-07-08 - Latest RNA-seq11 Valid/Test Prediction bigWig Export

- Run type: evaluation and generated-artifact export
- Purpose: Export the latest full-MSE validation leader's valid and test predictions as per-track bigWig files for biological genome-browser review.
- Git commit: `ac63bd8`
- Branch: `setup/agent-maintenance`
- Host: `HY-GPU` (`hy8`)
- Slurm job ID: Not applicable; HY-GPU is a non-Slurm server.
- Slurm request: Not applicable.
- Environment: Conda environment `alphagenome`; Python `/home/zelinli6/miniconda3/envs/alphagenome/bin/python`; PyTorch `2.11.0+cu128`; export restricted to `CUDA_VISIBLE_DEVICES=2`.
- Command:

```bash
mkdir -p logs/rna_seq11_bigwig_exports_20260708

CUDA_VISIBLE_DEVICES=2 conda run -n alphagenome python -u \
  scripts/alphagenome_rna_seq11_export_bigwig.py \
  --checkpoint runs/rna_seq11_1bpB_broad_w12_track_hard15_beta05_from_w10best1000_lr1e-5_seedrand_1500steps/adapter_head_best.pt \
  --weights weights/alphagenome_pytorch/model_all_folds.safetensors \
  --splits valid,test \
  --output-dir runs/rna_seq11_1bpB_broad_w12_track_hard15_beta05_from_w10best1000_lr1e-5_seedrand_1500steps/bigwig_exports/valid_test_20260708_expm1 \
  --prediction-transform expm1 \
  --clip-min 0 \
  --device cuda \
  > logs/rna_seq11_bigwig_exports_20260708/latest_valid_test_export.log 2>&1
```

- Input data: Valid NPZ dataset `alphagenome_custom/datasets/rna_seq_npz_valid` (`39` examples on chromosome `V`); test NPZ dataset `alphagenome_custom/datasets/rna_seq_npz_test` (`33` examples on chromosome `X`); base weights `weights/alphagenome_pytorch/model_all_folds.safetensors`; selected checkpoint `runs/rna_seq11_1bpB_broad_w12_track_hard15_beta05_from_w10best1000_lr1e-5_seedrand_1500steps/adapter_head_best.pt`, ranked first in `runs/rna_seq11_top15_extended_diagnostics_20260529/top15_models.tsv` with valid full MSE `0.83557710`.
- Output path: `runs/rna_seq11_1bpB_broad_w12_track_hard15_beta05_from_w10best1000_lr1e-5_seedrand_1500steps/bigwig_exports/valid_test_20260708_expm1/`; log at `logs/rna_seq11_bigwig_exports_20260708/latest_valid_test_export.log`. These generated outputs are ignored by Git.
- Result summary: Completed successfully. The export wrote 22 prediction bigWigs: 11 valid files under `valid/` and 11 test files under `test/`, plus `bigwig_manifest.tsv` for each split and `export_summary.json`. Predictions were merged by mean over overlapping 1,048,576 bp windows, converted from model log1p space back to raw signal scale with `expm1`, and clipped at zero.
- Verification: `python -m py_compile scripts/alphagenome_rna_seq11_export_bigwig.py` passed. `pyBigWig` read-back opened all 22 files as valid bigWigs. Valid files cover chromosome `V` length `20,924,180`; test files cover chromosome `X` length `17,718,942`. Manifest line counts are `12` for valid and `12` for test. No `.tmp` files remained. Post-run `nvidia-smi` showed no active project GPU process on GPU 2.
- Failures or warnings: The held-out test split was read and exported at user request for visualization. These test bigWigs must not be used for further model selection or hyperparameter tuning. The exported values are model predictions in raw-signal scale, not observed RNA-seq signal.
- Next actions: Share `valid/`, `test/`, `bigwig_manifest.tsv`, and `export_summary.json` with collaborators along with the warning that test output is for inspection/reporting only.
- Claim status: verified

## 2026-05-30 - RNA-seq11 All-Candidate Representation Utility Benchmark

- Run type: evaluation and analysis
- Purpose: Evaluate all eligible trained RNA-seq11 checkpoints by biological signal utility rather than only full-resolution valid MSE. The run added all-candidate extended valid diagnostics, train/valid gene-profile downstream probes, and composite utility leaderboards. The held-out test split was not read.
- Git commit: `13d67e6` at run start, with new probe and summary scripts added in the working tree for this run.
- Branch: `setup/agent-maintenance`
- Host: `HY-GPU`
- Slurm job ID: Not applicable; HY-GPU is a non-Slurm server.
- Slurm request: Not applicable.
- Environment: Conda environment `alphagenome`; Python `3.12.13`; PyTorch `2.11.0+cu128`; GPU jobs restricted to `CUDA_VISIBLE_DEVICES=2` and `CUDA_VISIBLE_DEVICES=3`.
- Command:

```bash
conda activate alphagenome
cd /home/zelinli6/Alphagenome

conda run -n alphagenome python scripts/rna_seq11_top15_extended_diagnostics.py \
  --runs-dir runs \
  --output-dir runs/rna_seq11_all_representation_inventory_20260530 \
  --top-n 999 \
  --only-write-top15

CUDA_VISIBLE_DEVICES=2 PYTHONUNBUFFERED=1 \
conda run -n alphagenome python -u scripts/rna_seq11_top15_extended_diagnostics.py \
  --runs-dir runs \
  --output-dir runs/rna_seq11_all_representation_diagnostics_20260530/shard0 \
  --top-n 999 --num-shards 2 --rank-shard 0 \
  --device auto --batch-size 1 --num-workers 0 \
  --spearman-sample-size 300000

CUDA_VISIBLE_DEVICES=3 PYTHONUNBUFFERED=1 \
conda run -n alphagenome python -u scripts/rna_seq11_top15_extended_diagnostics.py \
  --runs-dir runs \
  --output-dir runs/rna_seq11_all_representation_diagnostics_20260530/shard1 \
  --top-n 999 --num-shards 2 --rank-shard 1 \
  --device auto --batch-size 1 --num-workers 0 \
  --spearman-sample-size 300000

# Shared-embedding downstream probe after smoke validation.
CUDA_VISIBLE_DEVICES=<2-or-3> PYTHONUNBUFFERED=1 \
conda run -n alphagenome python -u scripts/rna_seq11_gene_profile_probe.py \
  --candidate-tsv runs/rna_seq11_all_representation_inventory_20260530/top15_models.tsv \
  --top-n 999 \
  --output-dir runs/rna_seq11_gene_profile_probe_all_20260530/shard<0-or-1> \
  --rank-shard <0-or-1> --num-shards 2 \
  --batch-size 1 --num-workers 0 \
  --probe-steps 300 --device auto \
  --eval-mode shared-embeddings --skip-existing

conda run -n alphagenome python scripts/rna_seq11_representation_utility_summary.py \
  --inventory-tsv runs/rna_seq11_all_representation_inventory_20260530/top15_models.tsv \
  --diagnostics-dirs \
    runs/rna_seq11_all_representation_diagnostics_20260530/shard0 \
    runs/rna_seq11_all_representation_diagnostics_20260530/shard1 \
  --probe-dirs \
    runs/rna_seq11_gene_profile_probe_all_20260530/shard0 \
    runs/rna_seq11_gene_profile_probe_all_20260530/shard1 \
  --output-dir runs/rna_seq11_representation_utility_summary_20260530
```

- Input data: Existing checkpoints under `runs/`; train NPZ data from `alphagenome_custom/datasets/rna_seq_npz_train`; valid NPZ data from `alphagenome_custom/datasets/rna_seq_npz_valid`; GTF annotations from `alphagenome_custom/reference/Caenorhabditis_elegans.WBcel235.115.gtf`; base weights from `weights/alphagenome_pytorch/model_all_folds.safetensors`.
- Output path: `runs/rna_seq11_all_representation_inventory_20260530`, `runs/rna_seq11_all_representation_diagnostics_20260530`, `runs/rna_seq11_gene_profile_probe_all_20260530`, `runs/rna_seq11_representation_utility_summary_20260530`; logs under `logs/rna_seq11_all_representation_20260530` and `logs/rna_seq11_gene_profile_probe_all_20260530`; report at `docs/rna_seq11_all_representation_utility_20260530.md`.
- Result summary: 171 unique eligible candidates were evaluated. The full-MSE reference remained `rna_seq11_1bpB_broad_w12_track_hard15_beta05_from_w10best1000_lr1e-5_seedrand_1500steps` with full MSE `0.83557710`. The top composite representation-utility candidate was `rna_seq11_1bpB_region_w1_rank1_exongene_v2_lr1e-5_1000steps` with representation utility score `0.88814617` and full MSE `0.83667865`. The exon/gene-body winner was `rna_seq11_1bpB_region_w1_rank1_exongene_v1_lr1e-5_1500steps`; the high-signal amplitude winner was `rna_seq11_1bpB_sigweight_w1_rank1_top5x3_a075_b05_lr1e-5_1500steps`; the near-frontier probe co-winner was `rna_seq11_1bpB_broad_w5_track_hard15_from_w4best2000_lr1e-5_seedrand_1500steps`.
- Verification: `py_compile` passed for `scripts/rna_seq11_gene_profile_probe.py` and `scripts/rna_seq11_representation_utility_summary.py`. Shared-embedding probe smoke completed for two models and small train/valid subsets. Extended diagnostics row counts were `per_track_metrics.tsv` shard0 `1033` and shard1 `1021`. Final probe metric row counts were shard0 `87` and shard1 `86`; per-track probe row counts were shard0 `2839` and shard1 `2806`. The final utility summary has `172` lines including header and `171` unique run IDs.
- Failures or warnings: An initial per-model probe implementation was stopped after writing a few rows because it repeated the frozen AlphaGenome encode for every checkpoint and would have taken many hours. A shared-embedding implementation was added, smoke-tested, and used with `--skip-existing` to complete the remaining candidates. No duplicate run IDs were present in the final probe tables. No test metrics were computed.
- Next actions: Use the representation winner for biological-utility comparisons, retain the full-MSE rank1 checkpoint as the pixel-error reference, and treat signal-weighted runs as high-expression/intestine co-winners rather than all-purpose replacements.
- Claim status: verified

## 2026-05-29 - RNA-seq11 Adaptive Objective and Calibration Screen

- Run type: training, evaluation, and analysis
- Purpose: Implement and test validation-only adaptive training objectives targeting high-signal amplitude, exon/gene-body error, intestine tracks, gene-level expression, and local 1 bp boundary/gradient behavior. Model selection remained unweighted valid full-resolution MSE; the held-out test split was not read.
- Git commit: `71ee5dd`
- Branch: `setup/agent-maintenance`
- Host: `HY-GPU`
- Slurm job ID: Not applicable; HY-GPU is a non-Slurm server.
- Slurm request: Not applicable.
- Environment: Conda environment `alphagenome`; Python `3.12.13`; PyTorch `2.11.0+cu128`; jobs restricted to `CUDA_VISIBLE_DEVICES=2` and `CUDA_VISIBLE_DEVICES=3`.
- Command:

```bash
conda activate alphagenome
cd /home/zelinli6/Alphagenome

python -m py_compile \
  scripts/alphagenome_rna_seq11_adapter.py \
  scripts/alphagenome_rna_seq11_finetune.py \
  scripts/alphagenome_rna_seq11_eval.py \
  scripts/rna_seq11_top15_extended_diagnostics.py

# Representative training command shape. Individual runs varied init checkpoint,
# output-dir, weighted objective flags, calibration flags, learning rate, and steps.
CUDA_VISIBLE_DEVICES=<2-or-3> PYTHONUNBUFFERED=1 \
conda run -n alphagenome python -u scripts/alphagenome_rna_seq11_finetune.py \
  --train-dataset-dir alphagenome_custom/datasets/rna_seq_npz_train \
  --valid-dataset-dir alphagenome_custom/datasets/rna_seq_npz_valid \
  --weights weights/alphagenome_pytorch/model_all_folds.safetensors \
  --residual-base-checkpoint runs/rna_seq11_phase4_lr_schedule_5000_conv5_h256_step4000_fulllog1p_hybrid_b1_lr0.001_seed20260522_5000steps_step_s4000/adapter_head_best.pt \
  --head-type linear \
  --embedding-resolution 1 \
  --linear-head-architecture conv5 \
  --linear-hidden-channels 128 \
  --linear-input-bottleneck-channels 128 \
  --linear-loss-type hybrid \
  --hybrid-loss-alpha 0.75 \
  --smooth-l1-beta 0.5 \
  --linear-target-space full-log1p \
  --target-transform log1p \
  --track-loss-weight-preset intestine-hard-boost-1.5 \
  --selection-metric full-mse \
  --batch-size 1 \
  --grad-accum-steps 1 \
  --num-workers 0 \
  --device auto

CUDA_VISIBLE_DEVICES=2 PYTHONUNBUFFERED=1 \
conda run -n alphagenome python -u scripts/rna_seq11_top15_extended_diagnostics.py \
  --runs-dir runs/rna_seq11_phase5_final_diag_candidates_20260529 \
  --output-dir runs/rna_seq11_phase5_final_diagnostics_20260529 \
  --top-n 8 \
  --device auto \
  --batch-size 1 \
  --num-workers 0
```

- Input data: `alphagenome_custom/datasets/rna_seq_npz_train` and `alphagenome_custom/datasets/rna_seq_npz_valid`; GTF annotations from `alphagenome_custom/reference/Caenorhabditis_elegans.WBcel235.115.gtf`; base weights from `weights/alphagenome_pytorch/model_all_folds.safetensors`; frozen 128 bp residual base checkpoint from `runs/rna_seq11_phase4_lr_schedule_5000_conv5_h256_step4000_fulllog1p_hybrid_b1_lr0.001_seed20260522_5000steps_step_s4000/adapter_head_best.pt`.
- Output path: Run outputs under `runs/rna_seq11_1bpB_sigweight_*`, `runs/rna_seq11_1bpB_region_*`, `runs/rna_seq11_1bpB_geneaux_*`, `runs/rna_seq11_1bpB_calib_*`; diagnostics under `runs/rna_seq11_phase1_sigweight_diagnostics_20260529`, `runs/rna_seq11_phase2_region_diagnostics_v1v2_20260529`, `runs/rna_seq11_phase3_geneaux_diagnostics_20260529`, and `runs/rna_seq11_phase5_final_diagnostics_20260529`; detailed run record at `docs/rna_seq11_adaptive_objective_screen_20260529.md`.
- Result summary: No new candidate improved unweighted valid full MSE over the current rank1 baseline `0.83557710`. Best new full-MSE candidate was `rank1_calib_trackaffine_lr3e5` at `0.83577257`. Signal weighting improved high-signal amplitude metrics but hurt full MSE. Region weighting improved exon/gene-body MSE but hurt full MSE and intergenic error. Gene auxiliary loss at `lambda=0.03` did not improve gene-MSE. Track-affine calibration slightly improved high-signal amplitude but did not improve full MSE, so high-signal-gated calibration was not run.
- Verification: `py_compile` passed for the edited scripts. Phase0 2-step smoke run completed at `runs/rna_seq11_phase0_code_smoke_20260529_sig_region_gene_calib_2steps`. Final diagnostics completed for 8 candidates with row counts `per_track=97`, `gene=97`, `signal_strata=673`, `region=481`, and `resolution=769`.
- Failures or warnings: Git push of commit `71ee5dd` was not attempted successfully in this session because the checkout lacks GitHub HTTPS credentials; this is not a code failure. Several runs were deliberately stopped after intermediate validation showed plateau or full-MSE degradation. No test metrics were computed.
- Next actions: Keep the current rank1 checkpoint as the main validation winner and rank3 as a gene-level co-reference. Future work should change model capacity or target parameterization for high-signal amplitude/local-boundary behavior rather than simply increasing loss weights on the same residual head.
- Claim status: verified

## 2026-05-29 - RNA-seq11 Top-15 Extended Valid Diagnostics

- Run type: evaluation and analysis
- Purpose: Compute extended valid-only diagnostics for the top 15 existing RNA-seq11 checkpoints ranked by valid full-resolution log1p MSE with Pearson as tie-breaker. Diagnostics include per-track MSE/MAE/Pearson/sampled Spearman, signal strata, gene-level expression, TSS/promoter/gene-body/exon/intron/intergenic regions, high-signal top-bin overlap, and resolution-specific pooled/gradient metrics.
- Git commit: `be22349` at run start, with `scripts/rna_seq11_top15_extended_diagnostics.py` added in the working tree for this run.
- Branch: `setup/agent-maintenance`
- Host: `HY-GPU`
- Slurm job ID: Not applicable; HY-GPU is a non-Slurm server.
- Slurm request: Not applicable.
- Environment: Conda environment `alphagenome`; GPU work restricted to `CUDA_VISIBLE_DEVICES=2` and `CUDA_VISIBLE_DEVICES=3` for two evaluation shards. GPU `1` had an unrelated VLLM process; GPU `0` was not used.
- Command:

```bash
conda activate alphagenome
cd /home/zelinli6/Alphagenome

python -m py_compile scripts/rna_seq11_top15_extended_diagnostics.py

python scripts/rna_seq11_top15_extended_diagnostics.py \
  --only-write-top15 \
  --output-dir runs/rna_seq11_top15_extended_diagnostics_20260529

CUDA_VISIBLE_DEVICES=2 python scripts/rna_seq11_top15_extended_diagnostics.py \
  --top-n 15 \
  --num-shards 2 \
  --rank-shard 0 \
  --output-dir runs/rna_seq11_top15_extended_diagnostics_20260529/shard0 \
  --device auto \
  --spearman-sample-size 300000

CUDA_VISIBLE_DEVICES=3 python scripts/rna_seq11_top15_extended_diagnostics.py \
  --top-n 15 \
  --num-shards 2 \
  --rank-shard 1 \
  --output-dir runs/rna_seq11_top15_extended_diagnostics_20260529/shard1 \
  --device auto \
  --spearman-sample-size 300000
```

- Input data: `alphagenome_custom/datasets/rna_seq_npz_valid` only; GTF annotations from `alphagenome_custom/reference/Caenorhabditis_elegans.WBcel235.115.gtf`; checkpoints selected from existing `runs/` outputs. The test split was not read.
- Output path: `runs/rna_seq11_top15_extended_diagnostics_20260529`; summary report at `docs/rna_seq11_top15_extended_diagnostics_20260529.md`.
- Result summary: Rank 1 is `runs/rna_seq11_1bpB_broad_w12_track_hard15_beta05_from_w10best1000_lr1e-5_seedrand_1500steps/adapter_head_best.pt`, valid MSE `0.83557710`, MAE `0.58180543`, Pearson `0.69463443`. It improves over the prior 128 bp conv5 hybrid step4000 baseline by `-0.04474871` MSE and `+0.02071778` Pearson. Extended diagnostics show remaining weakness in intestine tracks, high-signal/top-quantile regions, exon/gene-body regions, and local high-resolution gradient behavior.
- Verification: `py_compile` passed. Smoke diagnostics with `--top-n 1 --max-examples 1` completed before full evaluation. The merged tables cover ranks `1-15`; row counts were `per_track=180`, `signal_strata=1260`, `region=900`, `resolution=1440`, `window=585`, `gene=180`, `high_signal_localization=360`, and `model_config=15`.
- Failures or warnings: Full-resolution Spearman is sampled at about 300,000 positions per track rather than exact over all bases. Region and gene-level metrics use non-overlapping valid-window cores to avoid double-counting overlapping valid windows. No test metrics were computed.
- Next actions: Treat rank 1 as the MSE winner and rank 3 as a co-best gene/pooled-resolution candidate. Next exploration should target high-signal, exon/gene-body, and gene-level calibration objectives rather than more tiny low-learning-rate continuations.
- Claim status: verified

## 2026-05-28 - RNA-seq11 1bp-B Head, Fusion, and 1bp-C Sanity Sweep

- Run type: training and evaluation
- Purpose: Continue validation-only exploration on HY-GPU GPUs `2,3` after GPUs `0,1` were reserved for other work. The run tested whether expanding the current 1bp-B residual-correction head, adding base-prediction fusion, or using direct full-resolution 1bp-C depthwise heads could improve over the current validation leader without touching the held-out test split.
- Git commit: `0b2b0ce`
- Branch: `setup/agent-maintenance`
- Host: `HY-GPU`
- Slurm job ID: Not applicable; HY-GPU is a non-Slurm server
- Slurm request: Not applicable
- Environment: Conda environment `alphagenome`; independent single-GPU jobs with `CUDA_VISIBLE_DEVICES=2` or `CUDA_VISIBLE_DEVICES=3`; no DDP; no test split used.
- Command:

```bash
# Representative command shape. Individual runs changed output-dir,
# architecture, dilation/kernel, fusion flag, and learning rate.
CUDA_VISIBLE_DEVICES=<2-or-3> PYTHONUNBUFFERED=1 \
  /home/zelinli6/miniconda3/envs/alphagenome/bin/python -u \
  scripts/alphagenome_rna_seq11_finetune.py \
  --train-dataset-dir alphagenome_custom/datasets/rna_seq_npz_train \
  --valid-dataset-dir alphagenome_custom/datasets/rna_seq_npz_valid \
  --weights weights/alphagenome_pytorch/model_all_folds.safetensors \
  --residual-base-checkpoint runs/rna_seq11_phase4_lr_schedule_5000_conv5_h256_step4000_fulllog1p_hybrid_b1_lr0.001_seed20260522_5000steps_step_s4000/adapter_head_best.pt \
  --init-adapter-checkpoint runs/rna_seq11_1bpB_broad_w12_track_hard15_beta05_from_w10best1000_lr1e-5_seedrand_1500steps/adapter_head_best.pt \
  --init-adapter-checkpoint-mode matching-shapes \
  --head-type linear \
  --embedding-resolution 1 \
  --linear-loss-type hybrid \
  --hybrid-loss-alpha 0.5 \
  --smooth-l1-beta 0.5 \
  --track-loss-weight-preset intestine-hard-boost-1.5 \
  --linear-target-space full-log1p \
  --target-transform log1p \
  --batch-size 1 \
  --grad-accum-steps 1 \
  --num-workers 0 \
  --eval-every 250 \
  --learning-rate 3e-5 \
  --lr-schedule constant \
  --selection-metric full-mse \
  --residual-correction-l2 0.001 \
  --residual-correction-scale-init 0.01 \
  --seed -1 \
  --device auto
```

- Input data: `alphagenome_custom/datasets/rna_seq_npz_train`, `alphagenome_custom/datasets/rna_seq_npz_valid`, base weights `weights/alphagenome_pytorch/model_all_folds.safetensors`, frozen 128 bp base checkpoint `runs/rna_seq11_phase4_lr_schedule_5000_conv5_h256_step4000_fulllog1p_hybrid_b1_lr0.001_seed20260522_5000steps_step_s4000/adapter_head_best.pt`, current 1bp-B leader checkpoint `runs/rna_seq11_1bpB_broad_w12_track_hard15_beta05_from_w10best1000_lr1e-5_seedrand_1500steps/adapter_head_best.pt`
- Output path: `runs/rna_seq11_1bpB_head_sweep_20260528_*`, `runs/rna_seq11_1bpB_fusion_20260528_*`, `runs/rna_seq11_1bpC_sanity_20260528_*`; logs under `logs/rna_seq11_next_20260528/`
- Result summary: No run improved over the current validation leader `0.83557710` full MSE. Best new result was base-prediction fusion conv5 at `0.84498210` full MSE. Match-initialized head/fusion variants clustered around `0.845-0.852`; direct 1bp-C full-head sanity runs were much worse (`1.15-1.42` full MSE at 250 steps).
- Verification:
  - `fusion basepred conv5 b128 h128`: valid step 250 full MSE `0.84498210`, MAE `0.58407372`, Pearson `0.69025205`, common128 MSE `0.84939837`
  - `conv7 matchinit`: valid step 250 full MSE `0.84499249`, MAE `0.59198224`, Pearson `0.69008607`, common128 MSE `0.84404808`
  - `dilated-conv3 d4 matchinit`: valid step 250 full MSE `0.84734383`, Pearson `0.68874216`
  - `dilated-conv3 d2 matchinit`: valid step 250 full MSE `0.84758516`, Pearson `0.68842723`
  - `fusion basepred conv7`: valid step 250 full MSE `0.85084183`, Pearson `0.68730186`
  - `conv3x2 matchinit`: valid step 250 full MSE `0.85233903`, Pearson `0.68717763`
  - `conv5 b128 h256 fresh beta1 lr1e-4`: best observed valid full MSE `0.86154250` at step 500
  - `conv5 b256 h128 fresh beta1 lr1e-4`: best observed valid full MSE `0.86482886` at step 500
  - `1bp-C depthwise k31 lr3e-4`: valid step 250 full MSE `1.15021460`, Pearson `0.55297077`
  - `1bp-C depthwise k15 lr3e-4`: valid step 250 full MSE `1.27282590`, Pearson `0.53974875`
  - `1bp-C depthwise k15 lr1e-4`: valid step 250 full MSE `1.41703420`, Pearson `0.39250081`
- Failures or warnings: Two initial capacity runs were identified as weak because changed head shapes prevented strict warm-start and they effectively trained fresh correction heads at low learning rate; they were stopped after step-500 validation. Subsequent architecture/fusion runs used `--init-adapter-checkpoint-mode matching-shapes`. Weak runs were stopped after enough validation evidence to avoid unnecessary GPU use.
- Next actions: Keep the current 1bp-B leader as best validation checkpoint. Do not continue direct 1bp-C from scratch. If exploring further, prioritize diagnostics-guided objective changes around the existing leader or carefully designed residual/fusion methods that preserve the current conv5 correction behavior.
- Claim status: verified

## 2026-05-27 - RNA-seq11 1bp-B Broad Beta-0.5 Continuation

- Run type: training and evaluation
- Purpose: Continue validation-only 1 bp residual-correction exploration without using the held-out test split. The main question was whether the best 1bp-B residual-correction model could be improved by hard-intestine track-weighted training with `SmoothL1 beta=0.5`, staged continuation, and low-learning-rate or uniform-loss polish runs.
- Git commit: `13ef2c6`
- Branch: `setup/agent-maintenance`
- Host: `HY-GPU`
- Slurm job ID: Not applicable; HY-GPU is a non-Slurm server
- Slurm request: Not applicable
- Environment: Conda environment `alphagenome`; user initially allowed GPUs `0,1,2,3`, then requested GPUs `0,1` for urgent work, after which only GPUs `2,3` were used for AlphaGenome jobs. Jobs were independent single-GPU runs, not DDP.
- Command:

```bash
# Representative best continuation command. Follow-up runs changed only the
# output/init checkpoint, learning rate, or track-loss-weight preset.
CUDA_VISIBLE_DEVICES=<gpu> PYTHONUNBUFFERED=1 \
  /home/zelinli6/miniconda3/envs/alphagenome/bin/python -u \
  scripts/alphagenome_rna_seq11_finetune.py \
  --train-dataset-dir alphagenome_custom/datasets/rna_seq_npz_train \
  --valid-dataset-dir alphagenome_custom/datasets/rna_seq_npz_valid \
  --weights weights/alphagenome_pytorch/model_all_folds.safetensors \
  --output-dir runs/rna_seq11_1bpB_broad_w12_track_hard15_beta05_from_w10best1000_lr1e-5_seedrand_1500steps \
  --init-adapter-checkpoint runs/rna_seq11_1bpB_broad_w12_track_hard15_beta05_from_w10best1000_lr1e-5_seedrand_1500steps/init_beta05_w10best1000.pt \
  --residual-base-checkpoint runs/rna_seq11_phase4_lr_schedule_5000_conv5_h256_step4000_fulllog1p_hybrid_b1_lr0.001_seed20260522_5000steps_step_s4000/adapter_head_best.pt \
  --head-type linear \
  --embedding-resolution 1 \
  --linear-head-architecture conv5 \
  --linear-hidden-channels 128 \
  --linear-input-bottleneck-channels 128 \
  --linear-loss-type hybrid \
  --hybrid-loss-alpha 0.5 \
  --smooth-l1-beta 0.5 \
  --track-loss-weight-preset intestine-hard-boost-1.5 \
  --linear-target-space full-log1p \
  --target-transform log1p \
  --batch-size 1 \
  --grad-accum-steps 1 \
  --num-workers 0 \
  --max-steps 1500 \
  --eval-every 250 \
  --learning-rate 1e-5 \
  --lr-schedule constant \
  --selection-metric full-mse \
  --residual-correction-l2 0.001 \
  --early-stopping-min-steps 750 \
  --early-stopping-patience 4 \
  --seed -1 \
  --device auto

# Valid-only diagnostics for the final selected checkpoint.
CUDA_VISIBLE_DEVICES=2 PYTHONUNBUFFERED=1 \
  /home/zelinli6/miniconda3/envs/alphagenome/bin/python -u \
  scripts/alphagenome_rna_seq11_eval.py \
  --dataset-dir alphagenome_custom/datasets/rna_seq_npz_valid \
  --weights weights/alphagenome_pytorch/model_all_folds.safetensors \
  --checkpoint runs/rna_seq11_1bpB_broad_w12_track_hard15_beta05_from_w10best1000_lr1e-5_seedrand_1500steps/adapter_head_best.pt \
  --diagnostic-output runs/rna_seq11_1bpB_broad_diagnostics_20260527/w12_beta05_step750_best_diagnostics.tsv \
  --point-metrics \
  --batch-size 1 \
  --num-workers 0 \
  --device auto
```

- Input data: `alphagenome_custom/datasets/rna_seq_npz_train` and `alphagenome_custom/datasets/rna_seq_npz_valid`; base model weights `weights/alphagenome_pytorch/model_all_folds.safetensors`; frozen 128 bp residual base checkpoint `runs/rna_seq11_phase4_lr_schedule_5000_conv5_h256_step4000_fulllog1p_hybrid_b1_lr0.001_seed20260522_5000steps_step_s4000/adapter_head_best.pt`. The held-out test split was not used.
- Output path: Ignored run directories under `runs/rna_seq11_1bpB_broad_w*_...`; logs under `logs/rna_seq11_1bpB_broad_20260526/`; diagnostics under `runs/rna_seq11_1bpB_broad_diagnostics_20260527/`. These generated outputs must not be committed.
- Result summary: The best validation-only checkpoint from this round is `runs/rna_seq11_1bpB_broad_w12_track_hard15_beta05_from_w10best1000_lr1e-5_seedrand_1500steps/adapter_head_best.pt`, selected by valid full MSE at step 750. Metrics: full MSE `0.83557710`, MAE `0.58180543`, Pearson `0.69463443`, common128 MSE `0.83737121`, common128 MAE `0.57279302`.
- Additional results: The strongest direction was continuing hard-intestine weighted `hybrid alpha=0.5` with `SmoothL1 beta=0.5` at `lr=1e-5`. Low-learning-rate `3e-6` and `1e-6` polish runs from the best checkpoints improved or preserved MAE in some checkpoints but did not beat the best full MSE. Uniform-loss polish also did not beat the hard-intestine weighted objective by full MSE.
- Comparison: Versus the previous 1bp-B best (`0.83941437` full MSE, MAE `0.59241352`, Pearson `0.69303060`, common128 MSE `0.83968694`), this round improved full MSE by `0.00383727`, MAE by `0.01060809`, Pearson by `0.00160383`, and common128 MSE by `0.00231573`. Versus the 128 bp validation baseline (`0.88032581` full MSE, MAE `0.58482900`, Pearson `0.67391665`, common128 MSE `0.86829218`), it improved full MSE by `0.04474871`, MAE by `0.00302357`, Pearson by `0.02071778`, and common128 MSE by `0.03092097`.
- Diagnostics: Valid-only diagnostics were written to `runs/rna_seq11_1bpB_broad_diagnostics_20260527/w12_beta05_step750_best_diagnostics.tsv`. The hardest tracks by MSE remain intestine T4 (`1.0252335`), intestine T3 (`1.0082711`), and intestine T1 (`0.98125877`). Muscle tracks have the highest Pearson values, around `0.7256` to `0.7325`.
- Verification: The final process check showed no remaining `alphagenome_rna_seq11_finetune.py` or `alphagenome_rna_seq11_eval.py` jobs. GPU indices `0` and `1` were not used by AlphaGenome after the user requested them for urgent work.
- Failures or warnings: These are validation-only model-selection results, not held-out test results. Several dominated runs were intentionally stopped to conserve GPU time. Generated run directories, logs, diagnostics, checkpoints, and weights remain ignored artifacts and must not be committed.
- Next actions: Treat the `w12` step-750 checkpoint as the current validation leader. Recommended follow-up is a concise confirmation/diagnostic pass around same-learning-rate beta-0.5 continuation and per-track error analysis before any held-out test evaluation.
- Claim status: verified

## 2026-05-26 - RNA-seq11 1bp-B Upper-Bound Residual L2 Sweep

- Run type: training and evaluation
- Purpose: Continue validation-only 1bp-B residual-correction exploration from the best 1bp-B warm-start checkpoint, using the 128 bp `conv5` hybrid step-4000 checkpoint as a frozen base prediction and training only a 1 bp residual correction head. The held-out test split was not used.
- Git commit: `222cd09`
- Branch: `setup/agent-maintenance`
- Host: `HY-GPU`
- Slurm job ID: Not applicable; HY-GPU is a non-Slurm server
- Slurm request: Not applicable
- Environment: Conda environment `alphagenome`; user explicitly allowed GPUs `0,1,2,3`; independent single-GPU jobs were launched and monitored by the main controller, not DDP.
- Command:

```bash
# Main winning continuation, launched on one GPU.
CUDA_VISIBLE_DEVICES=2 PYTHONUNBUFFERED=1 \
  /home/zelinli6/miniconda3/envs/alphagenome/bin/python -u \
  scripts/alphagenome_rna_seq11_finetune.py \
  --train-dataset-dir alphagenome_custom/datasets/rna_seq_npz_train \
  --valid-dataset-dir alphagenome_custom/datasets/rna_seq_npz_valid \
  --weights weights/alphagenome_pytorch/model_all_folds.safetensors \
  --output-dir runs/rna_seq11_1bpB_upper_w2_resl2_1e-3_from_w1a05lr1e4_seedrand_2000steps \
  --init-adapter-checkpoint runs/rna_seq11_1bpB_upper_w1_hybrid_a05_lr1e-4_seedrand_4000steps/adapter_head_best.pt \
  --residual-base-checkpoint runs/rna_seq11_phase4_lr_schedule_5000_conv5_h256_step4000_fulllog1p_hybrid_b1_lr0.001_seed20260522_5000steps_step_s4000/adapter_head_best.pt \
  --head-type linear \
  --embedding-resolution 1 \
  --linear-head-architecture conv5 \
  --linear-hidden-channels 128 \
  --linear-input-bottleneck-channels 128 \
  --linear-loss-type hybrid \
  --hybrid-loss-alpha 0.5 \
  --smooth-l1-beta 1.0 \
  --linear-target-space full-log1p \
  --target-transform log1p \
  --batch-size 1 \
  --grad-accum-steps 1 \
  --num-workers 0 \
  --max-steps 2000 \
  --eval-every 250 \
  --learning-rate 1e-4 \
  --lr-schedule constant \
  --selection-metric full-mse \
  --residual-correction-l2 0.001 \
  --seed -1 \
  --device auto

# Follow-up variants used the same train/valid data, base checkpoint, head, target,
# selection metric, and random seed sentinel, changing only one or two fields:
# --residual-correction-l2 0.0005, 0.0015, or 0.003
# --residual-base-gate sigmoid --residual-base-gate-floor 0.5, 0.75
# --init-adapter-checkpoint <winning 1bp-B best> with --learning-rate 1e-5 or 3e-5

# Valid-only diagnostics for final 1bp-B best and the 128 bp baseline.
CUDA_VISIBLE_DEVICES=<free_gpu> PYTHONUNBUFFERED=1 \
  /home/zelinli6/miniconda3/envs/alphagenome/bin/python -u \
  scripts/alphagenome_rna_seq11_eval.py \
  --dataset-dir alphagenome_custom/datasets/rna_seq_npz_valid \
  --weights weights/alphagenome_pytorch/model_all_folds.safetensors \
  --checkpoint <checkpoint> \
  --diagnostic-output runs/rna_seq11_1bpB_upper_diagnostics_20260526/<name>_diagnostics.tsv \
  --point-metrics \
  --batch-size 1 \
  --num-workers 0 \
  --device auto
```

- Input data: `alphagenome_custom/datasets/rna_seq_npz_train` and `alphagenome_custom/datasets/rna_seq_npz_valid`; base weights `weights/alphagenome_pytorch/model_all_folds.safetensors`; frozen 128 bp base checkpoint `runs/rna_seq11_phase4_lr_schedule_5000_conv5_h256_step4000_fulllog1p_hybrid_b1_lr0.001_seed20260522_5000steps_step_s4000/adapter_head_best.pt`; 1bp-B warm-start checkpoint `runs/rna_seq11_1bpB_upper_w1_hybrid_a05_lr1e-4_seedrand_4000steps/adapter_head_best.pt`. The held-out test split was not used.
- Output path: Ignored run directories under `runs/rna_seq11_1bpB_upper_w*_...`; logs under `logs/rna_seq11_1bpB_upper_20260526/`; diagnostics under `runs/rna_seq11_1bpB_upper_diagnostics_20260526/`. These generated outputs must not be committed.
- Result summary: The best validation-only checkpoint is `runs/rna_seq11_1bpB_upper_w2_resl2_1e-3_from_w1a05lr1e4_seedrand_2000steps/adapter_head_best.pt`, selected by valid full MSE at step 2000. Metrics: full MSE `0.83941437`, MAE `0.59241352`, Pearson `0.69303060`, common128 MSE `0.83968694`, common128 MAE `0.58274677`. The step-1500 checkpoint was nearly tied on full MSE (`0.83941991`) and had better MAE/Pearson/common128 (`0.59154043`, `0.69311764`, `0.83934808`), but the stored best follows the configured `full-mse` selection metric.
- Additional results: Residual L2 `1e-3` was the strongest variant. `5e-4`, `1.5e-3`, and `3e-3` did not beat it. Base-gated residual variants improved MAE in some early checkpoints but worsened full MSE/common128. Low-learning-rate fine-tunes from the step-1500 checkpoint (`1e-5`, `3e-5`) did not beat the main constant-`1e-4` continuation by full MSE.
- Comparison: Versus the 128 bp validation baseline (`0.88032581` full MSE, Pearson `0.67391665`, common128 MSE `0.86829218`), the final 1bp-B best improved full MSE by `0.04091144`, Pearson by `0.01911395`, and common128 MSE by `0.02860524`, while MAE worsened from `0.58482900` to `0.59241352`. Versus the previous 1bp-B leader (`0.84327834` full MSE, MAE `0.59956417`, Pearson `0.69191858`, common128 MSE `0.84077399`), it improved all four headline metrics.
- Diagnostics: Final 1bp-B diagnostics were written to `runs/rna_seq11_1bpB_upper_diagnostics_20260526/w2_resl2_1e-3_final_best_diagnostics.tsv`; 128 bp baseline diagnostics were written to `runs/rna_seq11_1bpB_upper_diagnostics_20260526/baseline128_conv5_hybrid_step4000_diagnostics.tsv`. The final 1bp-B best improved per-track MSE and Pearson over the 128 bp baseline for all 11 tracks. Remaining high-MSE tracks are intestine T4 (`1.0291699`), intestine T3 (`1.0110442`), and intestine T1 (`0.98541102`). Muscle tracks remain stronger, with Pearson about `0.724` to `0.731`.
- Verification: `python -m py_compile scripts/*.py` passed in the `alphagenome` conda environment. `nvidia-smi` and process checks showed no remaining `alphagenome_rna_seq11_finetune.py` or `alphagenome_rna_seq11_eval.py` jobs after cleanup.
- Failures or warnings: These are validation-only model-selection results, not held-out test results. Several dominated jobs were intentionally terminated after enough validation evidence to conserve GPU time. The push of commit `222cd09` still requires user-side GitHub credentials.
- Next actions: Treat `runs/rna_seq11_1bpB_upper_w2_resl2_1e-3_from_w1a05lr1e4_seedrand_2000steps/adapter_head_best.pt` as the current validation leader. For the next round, consider diagnostics-guided work on intestine tracks and a small confirmation run around the step-1500/2000 tradeoff before any held-out test evaluation.
- Claim status: verified

## 2026-05-26 - 1bp-B Residual Correction Gate and Extension

- Run type: training and evaluation
- Purpose: Test whether 1 bp AlphaGenome embeddings add validation-only signal beyond the current 128 bp winner by freezing the 128 bp `conv5` hybrid step-4000 checkpoint as a base prediction and training only a small 1 bp residual correction head. The held-out test split was not used.
- Git commit: `fae06bf`
- Branch: `setup/agent-maintenance`
- Host: `HY-GPU`
- Slurm job ID: Not applicable; HY-GPU is a non-Slurm server
- Slurm request: Not applicable
- Environment: Conda environment `alphagenome`; user explicitly allowed `CUDA_VISIBLE_DEVICES=0,1,2,3`; four independent single-GPU runs were launched in parallel, not DDP.
- Command:

```bash
# Initial 1000-step residual gate; repeated across GPUs/configs.
CUDA_VISIBLE_DEVICES=<gpu> PYTHONUNBUFFERED=1 \
  /home/zelinli6/miniconda3/envs/alphagenome/bin/python -u \
  scripts/alphagenome_rna_seq11_finetune.py \
  --train-dataset-dir alphagenome_custom/datasets/rna_seq_npz_train \
  --valid-dataset-dir alphagenome_custom/datasets/rna_seq_npz_valid \
  --weights weights/alphagenome_pytorch/model_all_folds.safetensors \
  --output-dir runs/rna_seq11_1bpB_residual_20260526_<CONFIG> \
  --head-type linear \
  --embedding-resolution 1 \
  --linear-input-bottleneck-channels <64_or_128> \
  --linear-head-architecture <conv3_or_conv5> \
  --linear-hidden-channels <64_or_128> \
  --linear-target-space full-log1p \
  --target-transform log1p \
  --linear-loss-type hybrid \
  --hybrid-loss-alpha 0.5 \
  --smooth-l1-beta 1.0 \
  --residual-base-checkpoint runs/rna_seq11_phase4_lr_schedule_5000_conv5_h256_step4000_fulllog1p_hybrid_b1_lr0.001_seed20260522_5000steps_step_s4000/adapter_head_best.pt \
  --residual-correction-scale-init <0.01_or_0.1> \
  --selection-metric full-mse \
  --batch-size 1 \
  --grad-accum-steps 1 \
  --num-workers 0 \
  --max-steps 1000 \
  --eval-every 250 \
  --learning-rate 1e-3 \
  --seed 20260526 \
  --device auto

# Warm-start extension; repeated for top configs/losses.
CUDA_VISIBLE_DEVICES=<gpu> PYTHONUNBUFFERED=1 \
  /home/zelinli6/miniconda3/envs/alphagenome/bin/python -u \
  scripts/alphagenome_rna_seq11_finetune.py \
  --init-adapter-checkpoint <1bpB_best_checkpoint> \
  --learning-rate <3e-4_or_1e-4> \
  --linear-loss-type <hybrid_or_mse> \
  ...same dataset/base/residual arguments...
```

- Input data: `alphagenome_custom/datasets/rna_seq_npz_train` and `alphagenome_custom/datasets/rna_seq_npz_valid`; base weights `weights/alphagenome_pytorch/model_all_folds.safetensors`; frozen 128 bp base checkpoint `runs/rna_seq11_phase4_lr_schedule_5000_conv5_h256_step4000_fulllog1p_hybrid_b1_lr0.001_seed20260522_5000steps_step_s4000/adapter_head_best.pt`.
- Output path: ignored run directories under `runs/rna_seq11_1bpB_residual_20260526_*` and `runs/rna_seq11_1bpB_residual_extend_20260526_*`; logs under `logs/rna_seq11_1bpB_residual_20260526/` and `logs/rna_seq11_1bpB_residual_extend_20260526/`; diagnostics under `runs/rna_seq11_1bpB_best_diagnostics_20260526/`.
- Result summary: The 1 bp residual correction gate produced a clear validation-only improvement. The best initial 1000-step run was `conv5`, bottleneck/hidden `64`, residual scale init `0.01`, hybrid loss, best step `1000`, valid full MSE `0.84735946`, MAE `0.59909206`, Pearson `0.68895678`, common128 MSE `0.84613023`. Warm-start extension improved further: `conv5`, bottleneck/hidden `128`, residual scale init `0.01`, hybrid loss, lr `3e-4`, initialized from its initial 1000-step checkpoint, reached best step `1000` of the extension with valid full MSE `0.84327834`, MAE `0.59956417`, Pearson `0.69191858`, common128 MSE `0.84077399`. This is the current validation full-MSE leader. The best common128-only result among extension runs was `0.84009590` from `conv5` bottleneck/hidden `64`, hybrid lr `3e-4`, extension step `1000`, but its best full-MSE checkpoint was `0.84387702`.
- Verification: All reported training rows used train/valid only and `selection_metric=full-mse` for 1bp-B. Smoke tests before the gate verified `prediction_shape=1x11x1048576` for residual correction, frozen base with `base_has_grad=False`, and checkpoint reload for residual models. Final diagnostics for the selected 1bp-B leader reproduced valid full MSE `0.84327834`, MAE `0.59956417`, Pearson `0.69191858`, and common128 MSE `0.84077399`.
- Failures or warnings: These are validation-only model-selection results, not held-out test results. The selected 1bp-B leader improves full MSE and Pearson but has higher MAE than the previous 128 bp full-MSE leader (`0.59956417` vs `0.58482900`) because zero/low-signal strata worsened while medium/high-signal strata improved strongly. Generated checkpoints, logs, diagnostics, and run outputs remain ignored and must not be committed.
- Next actions: Treat the 1bp-B residual checkpoint as the new validation-selected candidate, but do not evaluate the held-out test split until model-selection criteria are frozen. If continuing train/valid exploration, compare the full-MSE leader against the common128 leader and inspect high-signal/zero-signal error tradeoffs before any fusion work.
- Claim status: verified

## 2026-05-26 - 1bp-A Pooled-to-128bp Gate

- Run type: training and evaluation
- Purpose: Test whether `embeddings_1bp`, bottlenecked and pooled to 128 bp bins, can match or beat direct 128 bp embeddings on the same `binned128-log1p-mean` RNA-seq target. This was a validation-only gate before considering larger 1 bp training.
- Git commit: `fae06bf`
- Branch: `setup/agent-maintenance`
- Host: `HY-GPU`
- Slurm job ID: Not applicable; HY-GPU is a non-Slurm server
- Slurm request: Not applicable
- Environment: Conda environment `alphagenome`; user explicitly allowed `CUDA_VISIBLE_DEVICES=0,1,2,3`; four independent single-GPU runs were launched in parallel.
- Command:

```bash
CUDA_VISIBLE_DEVICES=<gpu> PYTHONUNBUFFERED=1 \
  /home/zelinli6/miniconda3/envs/alphagenome/bin/python -u \
  scripts/alphagenome_rna_seq11_finetune.py \
  --train-dataset-dir alphagenome_custom/datasets/rna_seq_npz_train \
  --valid-dataset-dir alphagenome_custom/datasets/rna_seq_npz_valid \
  --weights weights/alphagenome_pytorch/model_all_folds.safetensors \
  --output-dir runs/rna_seq11_1bpA_gate_20260526_<CONFIG> \
  --head-type linear \
  --embedding-resolution 1 \
  --linear-input-bottleneck-channels <64_or_128> \
  --linear-head-architecture <conv3_or_conv5> \
  --linear-hidden-channels <64_or_128> \
  --linear-target-space binned128-log1p-mean \
  --target-transform log1p \
  --linear-loss-type hybrid \
  --hybrid-loss-alpha 0.5 \
  --smooth-l1-beta 1.0 \
  --selection-metric common128-mse \
  --batch-size 1 \
  --grad-accum-steps 1 \
  --num-workers 0 \
  --max-steps 1000 \
  --eval-every 250 \
  --learning-rate 1e-3 \
  --seed 20260526 \
  --device auto
```

- Input data: `alphagenome_custom/datasets/rna_seq_npz_train` and `alphagenome_custom/datasets/rna_seq_npz_valid`; base weights `weights/alphagenome_pytorch/model_all_folds.safetensors`.
- Output path: ignored run directories under `runs/rna_seq11_1bpA_gate_20260526_*`; logs under `logs/rna_seq11_1bp_gate_20260526/`.
- Result summary: The 1bp-A gate did not pass. Best common128 MSE was `0.99800044` from `conv5`, bottleneck/hidden `64`, hybrid loss, best step `1000`; other best common128 MSEs were `1.00116930` (`conv5` bottleneck/hidden `128`), `1.01619440` (`conv3` bottleneck/hidden `64`), and `1.02019590` (`conv3` bottleneck/hidden `128`). These are far worse than the direct 128 bp common128 reference around `0.853-0.855`.
- Verification: All runs completed with status code `0`, used `selection_metric=common128-mse`, and only read train/valid data. Prediction shape was `1x11x8192`, confirming pooling to 128 bp bins.
- Failures or warnings: Metrics are validation-only and do not use test. Because 1bp-A did not approach the 128 bp common128 reference, it was not extended to 2000 steps.
- Next actions: Do not promote 1bp-A as a standalone replacement for 128 bp embeddings. Use 1bp only as residual/correction signal unless future train/valid evidence changes this.
- Claim status: verified

## 2026-05-26 - Current 128bp Best Validation Diagnostics

- Run type: evaluation and analysis
- Purpose: Diagnose the current 128 bp full-MSE validation leader and the 128 bp common128 comparator before launching 1 bp experiments. The held-out test split was not used.
- Git commit: `fae06bf`
- Branch: `setup/agent-maintenance`
- Host: `HY-GPU`
- Slurm job ID: Not applicable; HY-GPU is a non-Slurm server
- Slurm request: Not applicable
- Environment: Conda environment `alphagenome`; user explicitly allowed `CUDA_VISIBLE_DEVICES=0,1,2,3`.
- Command:

```bash
CUDA_VISIBLE_DEVICES=<gpu> PYTHONUNBUFFERED=1 \
  /home/zelinli6/miniconda3/envs/alphagenome/bin/python -u \
  scripts/alphagenome_rna_seq11_eval.py \
  --dataset-dir alphagenome_custom/datasets/rna_seq_npz_valid \
  --weights weights/alphagenome_pytorch/model_all_folds.safetensors \
  --checkpoint <128bp_checkpoint> \
  --point-metrics \
  --spearman-sample-size 200000 \
  --spearman-seed 20260526 \
  --metrics-output <point_metrics.tsv> \
  --diagnostic-output <diagnostic.tsv> \
  --device auto

CUDA_VISIBLE_DEVICES=<gpu> PYTHONUNBUFFERED=1 \
  /home/zelinli6/miniconda3/envs/alphagenome/bin/python -u \
  scripts/alphagenome_rna_seq11_eval.py \
  --dataset-dir alphagenome_custom/datasets/rna_seq_npz_valid \
  --weights weights/alphagenome_pytorch/model_all_folds.safetensors \
  --checkpoint <128bp_checkpoint> \
  --common-128bp-metrics log1p-mean \
  --metrics-output <common128.tsv> \
  --device auto
```

- Input data: `alphagenome_custom/datasets/rna_seq_npz_valid`; base weights `weights/alphagenome_pytorch/model_all_folds.safetensors`; checkpoints `runs/rna_seq11_phase4_lr_schedule_5000_conv5_h256_step4000_fulllog1p_hybrid_b1_lr0.001_seed20260522_5000steps_step_s4000/adapter_head_best.pt` and `runs/rna_seq11_phase4_lr_schedule_5000_conv5_h256_constant_fulllog1p_hybrid_b1_lr0.001_seed20260522_5000steps/adapter_head_best.pt`.
- Output path: diagnostics under `runs/rna_seq11_best_diagnostics_20260526/`; logs under `logs/rna_seq11_best_diagnostics_20260526/`.
- Result summary: The 128 bp full-MSE leader reproduced valid full MSE `0.88032581`, MAE `0.58482900`, Pearson `0.67391665`, common128 MSE `0.86829218`. Its largest stratum error was high signal `log1p>3`, MSE `5.8507643`; worst tracks by MSE were `RNA_SEQ_005`, `RNA_SEQ_004`, and `RNA_SEQ_002`. The 128 bp constant-hybrid comparator had slightly worse full MSE `0.88462103` but better common128 MSE `0.85320788` and lower high-signal stratum MSE `5.2333010`, while worsening zero/low-signal strata.
- Verification: Diagnostics wrote full point metrics, per-track metrics, per-window rows, strata rows, and common128 metrics for validation only.
- Failures or warnings: Diagnostics are interpretive validation analyses, not held-out test results. They show a tradeoff: full-MSE selection favored lower zero/low-signal error, while common128/high-signal behavior favored the constant-hybrid comparator.
- Next actions: Use these diagnostics as the baseline for evaluating 1 bp residual correction; do not use test until model selection is frozen.
- Claim status: verified

## 2026-05-26 - 1bp Residual Adapter Code Support

- Run type: engineering smoke test and validation
- Purpose: Add and verify code paths needed for 1bp-A bottleneck pooling and 1bp-B residual correction while preserving old checkpoint loading behavior.
- Git commit: `fae06bf`
- Branch: `setup/agent-maintenance`
- Host: `HY-GPU`
- Slurm job ID: Not applicable; HY-GPU is a non-Slurm server
- Slurm request: Not applicable
- Environment: Conda environment `alphagenome`; user explicitly allowed `CUDA_VISIBLE_DEVICES=0,1,2,3`.
- Command:

```bash
/home/zelinli6/miniconda3/envs/alphagenome/bin/python -m py_compile scripts/*.py
/home/zelinli6/miniconda3/envs/alphagenome/bin/python scripts/alphagenome_rna_seq11_finetune.py --help
/home/zelinli6/miniconda3/envs/alphagenome/bin/python scripts/alphagenome_rna_seq11_eval.py --help

# Two one-step smoke tests used max one train example and one valid example:
# 1bp-A bottleneck+pool, no checkpoint;
# 1bp-B frozen 128bp base plus 1bp residual correction, no checkpoint.
# A second pair saved tiny checkpoints and reloaded them with eval.py.
```

- Input data: First train/valid NPZ examples for smoke tests; base weights `weights/alphagenome_pytorch/model_all_folds.safetensors`; frozen 128 bp base checkpoint for 1bp-B smoke.
- Output path: ignored smoke directories under `runs/rna_seq11_1bpA_bottleneck_*` and `runs/rna_seq11_1bpB_residual_*`; logs under `logs/rna_seq11_best_diagnostics_20260526/`.
- Result summary: Added `--linear-input-bottleneck-channels`, `--residual-base-checkpoint`, and `--residual-correction-scale-init`; checkpoint save/reload now records optional bottleneck and frozen residual-base metadata. Smoke tests confirmed 1bp-A prediction shape `1x11x8192` and 1bp-B prediction shape `1x11x1048576`, with `base_has_grad=False`.
- Verification: `python -m py_compile scripts/*.py` passed; `finetune.py --help` and `eval.py --help` exposed the new options; save/reload smoke tests completed successfully.
- Failures or warnings: `git push` failed because HTTPS credentials are unavailable in this environment: `fatal: could not read Username for 'https://github.com': No such device or address`.
- Next actions: Use the new code path for validation-only 1bp-A and 1bp-B pilots, keeping test split untouched.
- Claim status: verified

## 2026-05-26 - 1bp-A Binned128 Conv3 Smoke

- Run type: smoke test
- Purpose: Verify the prepared 1bp-A path on one train example and one validation example: `embeddings_1bp` with 1536 channels, average pooling to 128 bp bins before the linear conv head, `conv3` hidden `64`, and `binned128-log1p-mean` target. This was a pipeline sanity check only, not a model-comparison run.
- Git commit: `0a760f8c3acdb3e0828cf1a01ea574409b49518e`
- Branch: `setup/agent-maintenance`
- Host: `HY-GPU`
- Slurm job ID: Not applicable; HY-GPU is a non-Slurm server
- Slurm request: Not applicable
- Environment: Conda environment `alphagenome`; Python `3.12.13`; PyTorch `2.11.0+cu128`; `CUDA_VISIBLE_DEVICES=2`; user explicitly allowed GPUs `0,1,2,3`
- Command:

```bash
CUDA_VISIBLE_DEVICES=2 PYTHONUNBUFFERED=1 /home/zelinli6/miniconda3/envs/alphagenome/bin/python -u \
  scripts/alphagenome_rna_seq11_finetune.py \
  --train-dataset-dir alphagenome_custom/datasets/rna_seq_npz_train \
  --valid-dataset-dir alphagenome_custom/datasets/rna_seq_npz_valid \
  --weights weights/alphagenome_pytorch/model_all_folds.safetensors \
  --output-dir runs/rna_seq11_1bpA_binned128_conv3_h64_smoke_20260526_1step \
  --head-type linear \
  --embedding-resolution 1 \
  --linear-head-architecture conv3 \
  --linear-hidden-channels 64 \
  --linear-target-space binned128-log1p-mean \
  --target-transform log1p \
  --linear-loss-type mse \
  --selection-metric full-mse \
  --batch-size 1 \
  --grad-accum-steps 1 \
  --num-workers 0 \
  --max-train-examples 1 \
  --max-valid-examples 1 \
  --max-steps 1 \
  --eval-every 1 \
  --learning-rate 1e-4 \
  --seed 20260526 \
  --device auto \
  --no-save-checkpoint
```

- Input data: First train example from `alphagenome_custom/datasets/rna_seq_npz_train` and first validation example from `alphagenome_custom/datasets/rna_seq_npz_valid`; base weights `weights/alphagenome_pytorch/model_all_folds.safetensors`. The held-out test split was not used.
- Output path: `runs/rna_seq11_1bpA_binned128_conv3_h64_smoke_20260526_1step/` and log `logs/rna_seq11_1bpA_binned128_conv3_h64_smoke_20260526_1step.log`; no checkpoint was saved.
- Result summary: Completed successfully. Observed `embedding_resolution=1`, `linear_target_space=binned128-log1p-mean`, `linear_head_architecture=conv3`, `linear_hidden_channels=64`, `trainable_parameters=111435`, `prediction_shape=1x11x8192`, `base_has_grad=False`, one-step train loss `3.09554`, and one-example validation full MSE `1.96545`.
- Verification: This confirms the 1bp-A binned path can execute with real AlphaGenome 1 bp embeddings and produce 128 bp-bin predictions. It does not test whether 1bp-A can beat the 128 bp binned baseline.
- Failures or warnings: Smoke metrics are environment and wiring validation only, not model results. Generated logs and run outputs remain ignored and must not be committed.
- Next actions: If continuing 1bp-A, run a 1000-step train/valid screen on GPU only after deciding the exact hidden size and loss; do not touch the held-out test split.
- Claim status: verified smoke-only

## 2026-05-26 - Full-MSE Selection Phase 4 LR Schedule Sweep

- Run type: training
- Purpose: Run Phase 4 validation-only LR schedule comparison on the Phase 3 `conv5` hybrid-loss winner: constant `1e-3`, cosine decay, step decay at `3000` to `3e-4`, and step decay at `4000` to `3e-4`, selecting all checkpoints by validation full MSE and not reading the held-out test split.
- Git commit: `da8f97a77460271a0522ec7020be2446ee8ad244`
- Branch: `setup/agent-maintenance`
- Host: `HY-GPU`
- Slurm job ID: Not applicable; HY-GPU is a non-Slurm server
- Slurm request: Not applicable
- Environment: Conda environment `alphagenome`; Python `3.12.13`; PyTorch `2.11.0+cu128`; user explicitly allowed `CUDA_VISIBLE_DEVICES=0,1,2,3`
- Command:

```bash
PYTHONUNBUFFERED=1 /home/zelinli6/miniconda3/envs/alphagenome/bin/python -u \
  scripts/rna_seq11_run_adaptation_phase_sweep.py \
  --date-tag 20260524 \
  --phase phase4 \
  --gpus 0,1,2,3 \
  --selection-metric full-mse
```

- Input data: `alphagenome_custom/datasets/rna_seq_npz_train` and `alphagenome_custom/datasets/rna_seq_npz_valid`; base weights `weights/alphagenome_pytorch/model_all_folds.safetensors`
- Output path: Summary TSV `runs/rna_seq11_adaptation_phase_sweep_20260524/sweep_results.tsv`; logs under `logs/rna_seq11_adaptation_phase_sweep_20260524/`; checkpoints under ignored `runs/rna_seq11_phase4_*`
- Result summary: Completed 4/4 validation-only Phase 4 runs. New best validation run was `conv5` hybrid loss with step LR decay at step `4000` from `1e-3` to `3e-4`, best step `4500`, valid full MSE `0.88032581`, MAE `0.58482900`, Pearson `0.67391665`, common128 MSE `0.86829218`, checkpoint `runs/rna_seq11_phase4_lr_schedule_5000_conv5_h256_step4000_fulllog1p_hybrid_b1_lr0.001_seed20260522_5000steps_step_s4000/adapter_head_best.pt`. The next Phase 4 runs were step decay at `3000` full MSE `0.88392895`, constant full MSE `0.88462103`, and cosine full MSE `0.88622749`.
- Verification: Driver reported `stage_end` and `sweep_done` at `2026-05-26T00:02:01`; all Phase 4 rows have status `completed`, return code `0`, and `selection_metric=full-mse`. Held-out test split was not used.
- Failures or warnings: Metrics are validation-only model-selection results, not held-out test results. Generated logs, checkpoints, and run outputs remain ignored and must not be committed. The best Phase 4 run improves validation full MSE over the Phase 3 leader `0.88446225`, but its common128 MSE `0.86829218` is higher than the Phase 3 leader common128 MSE `0.85788156`.
- Next actions: Treat `conv5` hybrid step@4000 as the current validation-full-MSE leader. Before any test evaluation, decide whether to prioritize full MSE or common128 MSE for model selection; for 1bp work, start with small train/valid sanity only.
- Claim status: verified validation-only

## 2026-05-25 - Full-MSE Selection Phase 3 Loss Sweep and MSE Fine-Tune

- Run type: training
- Purpose: Run Phase 3 validation-only loss variants on the top distinct 5000-step heads (`conv5` and `conv3`) plus two low-LR MSE fine-tunes from the Phase 1 `conv5` winner checkpoint, selecting all checkpoints by validation full MSE and not reading the held-out test split.
- Git commit: `8a661b11a79d33cc54fa1302768d564b6aef2f45`
- Branch: `setup/agent-maintenance`
- Host: `HY-GPU`
- Slurm job ID: Not applicable; HY-GPU is a non-Slurm server
- Slurm request: Not applicable
- Environment: Conda environment `alphagenome`; Python `3.12.13`; PyTorch `2.11.0+cu128`; user explicitly allowed `CUDA_VISIBLE_DEVICES=0,1,2,3`
- Command:

```bash
PYTHONUNBUFFERED=1 /home/zelinli6/miniconda3/envs/alphagenome/bin/python -u \
  scripts/rna_seq11_run_adaptation_phase_sweep.py \
  --date-tag 20260524 \
  --phase phase3 \
  --gpus 0,1,2,3 \
  --selection-metric full-mse
```

- Input data: `alphagenome_custom/datasets/rna_seq_npz_train` and `alphagenome_custom/datasets/rna_seq_npz_valid`; base weights `weights/alphagenome_pytorch/model_all_folds.safetensors`; MSE fine-tunes initialized from `runs/rna_seq11_phase1_winner_stability_5000_conv5_h256_fulllog1p_smoothl1_b1_lr0.001_seed20260522_5000steps/adapter_head_best.pt`
- Output path: Summary TSV `runs/rna_seq11_adaptation_phase_sweep_20260524/sweep_results.tsv`; logs under `logs/rna_seq11_adaptation_phase_sweep_20260524/`; checkpoints under ignored `runs/rna_seq11_phase3_*`
- Result summary: Completed 8/8 validation-only Phase 3 runs. New best validation run was `conv5` hybrid loss (`0.5*MSE + 0.5*SmoothL1`, beta `1.0`), best step `4500`, valid full MSE `0.88446225`, MAE `0.60360136`, Pearson `0.67363609`, common128 MSE `0.85788156`, checkpoint `runs/rna_seq11_phase3_loss_sweep_5000_conv5_h256_top1_hybrid_fulllog1p_hybrid_b1_lr0.001_seed20260522_5000steps/adapter_head_best.pt`. This improves over the Phase 1 `conv5` validation leader full MSE `0.88662339` but remains validation-only.
- Verification: Driver reported `stage_end` and `sweep_done` at `2026-05-25T18:56:43`; all Phase 3 rows have status `completed`, return code `0`, and `selection_metric=full-mse`. Held-out test split was not used.
- Failures or warnings: Metrics are validation-only model-selection results, not held-out test results. Generated logs, checkpoints, and run outputs remain ignored and must not be committed. Low-LR MSE fine-tunes did not beat the hybrid-loss winner; best MSE fine-tune full MSE was `0.88974877` at lr `1e-4`.
- Next actions: Run Phase 4 LR schedule comparison on the new `conv5` hybrid configuration: constant `1e-3`, cosine decay, step decay at `3000`, and step decay at `4000`.
- Claim status: verified validation-only

## 2026-05-25 - Full-MSE Selection Phase 2 Promote

- Run type: training
- Purpose: Promote the top 3 Phase 2 1000-step head-screen configs to 5000-step validation-only runs, still selecting checkpoints by validation full MSE and not reading the held-out test split.
- Git commit: `57b4fc26c9d42e3edcb92937a8e839d1cddd141b`
- Branch: `setup/agent-maintenance`
- Host: `HY-GPU`
- Slurm job ID: Not applicable; HY-GPU is a non-Slurm server
- Slurm request: Not applicable
- Environment: Conda environment `alphagenome`; Python `3.12.13`; PyTorch `2.11.0+cu128`; user explicitly allowed `CUDA_VISIBLE_DEVICES=0,1,2,3`
- Command:

```bash
PYTHONUNBUFFERED=1 /home/zelinli6/miniconda3/envs/alphagenome/bin/python -u \
  scripts/rna_seq11_run_adaptation_phase_sweep.py \
  --date-tag 20260524 \
  --phase phase2-promote \
  --gpus 0,1,2,3 \
  --selection-metric full-mse
```

- Input data: `alphagenome_custom/datasets/rna_seq_npz_train` and `alphagenome_custom/datasets/rna_seq_npz_valid`; base weights `weights/alphagenome_pytorch/model_all_folds.safetensors`
- Output path: Summary TSV `runs/rna_seq11_adaptation_phase_sweep_20260524/sweep_results.tsv`; logs under `logs/rna_seq11_adaptation_phase_sweep_20260524/`; checkpoints under ignored `runs/rna_seq11_phase2_head_promote_5000_*`
- Result summary: Completed 3/3 validation-only 5000-step promote runs. Best promoted config by full MSE was `dilated-conv3` hidden `256`, dilation `2`, best step `4750`, valid full MSE `0.90497149`, MAE `0.60353299`, Pearson `0.66171282`, common128 MSE `0.89870273`. The other promoted configs were `residual-conv3` hidden `512`, residual scale `0.1`, best step `3750`, full MSE `0.90897760`, MAE `0.60479374`, Pearson `0.66221824`, common128 MSE `0.89215360`; and `residual-conv3` hidden `512`, residual scale `1.0`, best step `4250`, full MSE `0.91719273`, MAE `0.59587331`, Pearson `0.65961935`, common128 MSE `0.90538565`.
- Verification: Driver reported `stage_end` and `sweep_done` at `2026-05-25T10:03:54`; all Phase 2 promote rows have status `completed`, return code `0`, and `selection_metric=full-mse`. Held-out test split was not used.
- Failures or warnings: Metrics are validation-only model-selection results, not held-out test results. Generated logs, checkpoints, and run outputs remain ignored and must not be committed. None of the Phase 2 promote runs beat the Phase 1 `conv5` seed `20260522` validation leader with full MSE `0.88662339`.
- Next actions: Keep Phase 1 `conv5` seed `20260522` as the current validation leader. Run Phase 3 loss/LR fine-tune against that `conv5` winner rather than automatically using the lower-performing Phase 2 promoted heads.
- Claim status: verified validation-only

## 2026-05-25 - Full-MSE Selection Phase 2 Head Screen

- Run type: training
- Purpose: Run the 13-job Phase 2 1000-step head screen on train/valid only, selecting best checkpoints by validation full MSE for residual hidden/scale variants, `conv3x2`, `conv7`, and dilated `conv3` heads.
- Git commit: `84685050e50e0e19a122a9a09858905c8dcc0942`
- Branch: `setup/agent-maintenance`
- Host: `HY-GPU`
- Slurm job ID: Not applicable; HY-GPU is a non-Slurm server
- Slurm request: Not applicable
- Environment: Conda environment `alphagenome`; Python `3.12.13`; PyTorch `2.11.0+cu128`; user explicitly allowed `CUDA_VISIBLE_DEVICES=0,1,2,3`
- Command:

```bash
PYTHONUNBUFFERED=1 /home/zelinli6/miniconda3/envs/alphagenome/bin/python -u \
  scripts/rna_seq11_run_adaptation_phase_sweep.py \
  --date-tag 20260524 \
  --phase phase2-screen \
  --gpus 0,1,2,3 \
  --selection-metric full-mse
```

- Input data: `alphagenome_custom/datasets/rna_seq_npz_train` and `alphagenome_custom/datasets/rna_seq_npz_valid`; base weights `weights/alphagenome_pytorch/model_all_folds.safetensors`
- Output path: Summary TSV `runs/rna_seq11_adaptation_phase_sweep_20260524/sweep_results.tsv`; logs under `logs/rna_seq11_adaptation_phase_sweep_20260524/`; checkpoints under ignored `runs/rna_seq11_phase2_head_screen_1000_*`
- Result summary: Completed 13/13 validation-only 1000-step screens. Top full-MSE configs were `residual-conv3` hidden `512`, residual scale `0.1`, full MSE `0.95572846`; `dilated-conv3` dilation `2`, full MSE `0.95791135`; and `residual-conv3` hidden `512`, residual scale `1.0`, full MSE `0.95821054`. For comparison, Phase 1 `conv5` seed `20260522` had full MSE `0.959042` at 1000 steps and later reached `0.88662339` at 4500 steps.
- Verification: Driver reported `stage_end` and `sweep_done` at `2026-05-25T05:46:50`; all Phase 2 screen rows have status `completed`, return code `0`, and `selection_metric=full-mse`. Held-out test split was not used.
- Failures or warnings: Metrics are validation-only model-selection results, not held-out test results. Generated logs, checkpoints, and run outputs remain ignored and must not be committed. The branch still could not be pushed from this environment because HTTPS Git credentials are unavailable.
- Next actions: Promote the top 2-3 Phase 2 screen configs to 5000-step runs while retaining Phase 1 `conv5` as the current validation leader.
- Claim status: verified validation-only

## 2026-05-25 - Full-MSE Selection Phase 1 Winner Stability

- Run type: training
- Purpose: Run the 9-job Phase 1 winner-stability screen on train/valid only, selecting best checkpoints by validation full MSE for `residual-conv3`, `conv3`, and `conv5` heads across seeds `20260515`, `20260522`, and `20260523`.
- Git commit: `58c46c6bc81f5ec98848c54c82a8e49891b8b784`
- Branch: `setup/agent-maintenance`
- Host: `HY-GPU`
- Slurm job ID: Not applicable; HY-GPU is a non-Slurm server
- Slurm request: Not applicable
- Environment: Conda environment `alphagenome`; Python `3.12.13`; PyTorch `2.11.0+cu128`; user explicitly allowed `CUDA_VISIBLE_DEVICES=0,1,2,3`
- Command:

```bash
conda run -n alphagenome python -u \
  scripts/rna_seq11_run_adaptation_phase_sweep.py \
  --date-tag 20260524 \
  --phase phase1 \
  --gpus 0,1,2,3 \
  --selection-metric full-mse
```

- Input data: `alphagenome_custom/datasets/rna_seq_npz_train` and `alphagenome_custom/datasets/rna_seq_npz_valid`; base weights `weights/alphagenome_pytorch/model_all_folds.safetensors`
- Output path: Summary TSV `runs/rna_seq11_adaptation_phase_sweep_20260524/sweep_results.tsv`; logs under `logs/rna_seq11_adaptation_phase_sweep_20260524/`; checkpoints under ignored `runs/rna_seq11_phase1_*`
- Result summary: Completed 9/9 validation-only runs. Best single run was `conv5` seed `20260522`, best step `4500`, valid full MSE `0.88662339`, MAE `0.57668143`, Pearson `0.67223756`, common128 MSE `0.87388091`. Head-level mean full MSE was `conv5` `0.89759223`, `conv3` `0.90369701`, and `residual-conv3` `0.91272342`. Head-level mean Pearson was `conv5` `0.66943618`, `conv3` `0.66358423`, and `residual-conv3` `0.66031831`.
- Verification: Driver reported `stage_end` and `sweep_done` at `2026-05-25T02:38:59`; all runs have status `completed`, return code `0`, and `selection_metric=full-mse`. Held-out test split was not used.
- Failures or warnings: Metrics are validation-only model-selection results, not held-out test results. Generated logs, checkpoints, and run outputs remain ignored and must not be committed. The branch still could not be pushed from this environment because HTTPS Git credentials are unavailable.
- Next actions: Use `conv5` as the Phase 2 primary head family and include `conv3` as a close stability comparator; deprioritize default `residual-conv3` unless residual-scale variants recover it.
- Claim status: verified validation-only

## 2026-05-24 - Full-MSE Selection Phase 0 Smoke and Phase 1 Launch

- Run type: smoke test and training
- Purpose: Validate Phase 0 code changes for full-MSE checkpoint selection, richer linear heads, SmoothL1/hybrid loss controls, residual scale init, dilated heads, and staged validation-only launcher; then launch Phase 1 winner-stability training without reading the held-out test split.
- Git commit: `828e2290a9450aecb33e5b338e9c17e1a24c5a21`
- Branch: `setup/agent-maintenance`
- Host: `HY-GPU`
- Slurm job ID: Not applicable; HY-GPU is a non-Slurm server
- Slurm request: Not applicable
- Environment: Conda environment `alphagenome`; Python `3.12.13`; PyTorch `2.11.0+cu128`; PyTorch CUDA `12.8`; driver CUDA `13.2`; user explicitly allowed `CUDA_VISIBLE_DEVICES=0,1,2,3`
- Command:

```bash
conda run -n alphagenome python -m py_compile scripts/*.py

CUDA_VISIBLE_DEVICES=2 conda run -n alphagenome python -u \
  scripts/alphagenome_rna_seq11_finetune.py \
  --train-dataset-dir alphagenome_custom/datasets/rna_seq_npz_train \
  --valid-dataset-dir alphagenome_custom/datasets/rna_seq_npz_valid \
  --weights weights/alphagenome_pytorch/model_all_folds.safetensors \
  --output-dir runs/rna_seq11_phase0_smoke_20260524_residualconv3_fullmse \
  --head-type linear \
  --embedding-resolution 128 \
  --linear-head-architecture residual-conv3 \
  --linear-hidden-channels 256 \
  --residual-scale-init 0.1 \
  --linear-loss-type smooth-l1 \
  --smooth-l1-beta 1.0 \
  --linear-target-space full-log1p \
  --target-transform log1p \
  --selection-metric full-mse \
  --batch-size 1 \
  --grad-accum-steps 1 \
  --num-workers 0 \
  --max-train-examples 1 \
  --max-valid-examples 1 \
  --max-steps 1 \
  --eval-every 1 \
  --learning-rate 1e-3 \
  --weight-decay 0 \
  --seed 20260524 \
  --device auto

conda run -n alphagenome python -u \
  scripts/rna_seq11_run_adaptation_phase_sweep.py \
  --date-tag 20260524 \
  --phase phase1 \
  --gpus 0,1,2,3 \
  --selection-metric full-mse
```

- Input data: `alphagenome_custom/datasets/rna_seq_npz_train` and `alphagenome_custom/datasets/rna_seq_npz_valid`; base weights `weights/alphagenome_pytorch/model_all_folds.safetensors`
- Output path: Phase 0 smoke output `runs/rna_seq11_phase0_smoke_20260524_residualconv3_fullmse`; Phase 1 summary `runs/rna_seq11_adaptation_phase_sweep_20260524/sweep_results.tsv`; Phase 1 logs `logs/rna_seq11_adaptation_phase_sweep_20260524/`
- Result summary: Phase 0 smoke completed successfully. Phase 1 launch pending in this entry.
- Verification: Pre-launch checks confirmed host `HY-GPU`, repository path `/home/zelinli6/Alphagenome`, conda environment `alphagenome`, PyTorch CUDA availability, and idle GPUs 0-3. Smoke output included `valid_loss`, `valid_full_mse`, `valid_full_mae`, `valid_pearson`, `valid_common128_mse`, `valid_common128_mae`, `selection_metric=full-mse`, and saved `adapter_head_best.pt`.
- Failures or warnings: Smoke metrics are environment validation only, not model results. Generated logs, checkpoints, and run outputs remain ignored and must not be committed.
- Next actions: Launch and monitor Phase 1; summarize validation-only metrics by head and seed before Phase 2.
- Claim status: verified for Phase 0 smoke; unverified for Phase 1 pending training

## 2026-05-22 - 128 bp Log1p-MSE Staged Sweep Launch

- Run type: training and evaluation
- Purpose: Launch the agreed validation-only sweep around the strongest baseline, `128bp embedding + simple linear head + log1p MSE`, using independent single-GPU jobs on four GPUs rather than DDP. The launcher runs hyperparameter screening, baseline seed stability, head screening, binned-target screening, promotion runs, optional 10000-step extensions, final seed confirmation, and a final valid diagnostic. The held-out test split is not used.
- Git commit: Recorded by the launcher at runtime in `logs/rna_seq11_128bp_log1p_sweep_20260522/driver.log` as `git_commit`.
- Branch: `setup/agent-maintenance`
- Host: `HY-GPU`
- Slurm job ID: Not applicable; HY-GPU is a non-Slurm server
- Slurm request: Not applicable
- Environment: Conda environment `alphagenome`; user explicitly allowed `CUDA_VISIBLE_DEVICES=0,1,2,3`
- Command:

```bash
setsid /home/zelinli6/miniconda3/envs/alphagenome/bin/python -u \
  scripts/rna_seq11_run_128bp_log1p_sweep.py \
  --date-tag 20260522 \
  --gpus 0,1,2,3 \
  > logs/rna_seq11_128bp_log1p_sweep_20260522/driver.log 2>&1 < /dev/null &
```

- Input data: `alphagenome_custom/datasets/rna_seq_npz_train` and `alphagenome_custom/datasets/rna_seq_npz_valid`
- Output path: ignored run directories under `runs/rna_seq11_sweep_*`; summary TSV `runs/rna_seq11_128bp_log1p_sweep_20260522/sweep_results.tsv`; logs under `logs/rna_seq11_128bp_log1p_sweep_20260522/`
- Result summary: Launched; results are pending.
- Verification: Pre-launch checks confirmed host `HY-GPU`, repository path `/home/zelinli6/Alphagenome`, commit `b211afd2c114b11fb21b8d9c7e9c1ae8938921b3`, conda environment `alphagenome`, PyTorch CUDA availability, and idle GPUs 0-3.
- Failures or warnings: The branch could not be pushed before launch because HTTPS credentials are unavailable in this environment. A first `nohup conda run ... &` launch attempt exited without starting the driver in the Codex tool environment, so the active launch uses `setsid` and the conda environment's direct Python executable. Generated logs, checkpoints, and run outputs remain ignored and must not be committed.
- Next actions: Monitor `logs/rna_seq11_128bp_log1p_sweep_20260522/driver.log` and summarize valid-only results after the launcher completes the staged sweep.
- Claim status: verified

## 2026-05-22 - 128 bp Log1p-MSE Sweep Tooling and Smoke Validation

- Run type: engineering smoke test and validation
- Purpose: Implement the next 128 bp frozen-trunk log1p-MSE sweep tooling before launching hour-scale screening runs. This adds richer linear head variants, SmoothL1 support, 128 bp binned target training, early stopping, checkpoint metadata, and validation diagnostics while keeping test split untouched.
- Git commit: `a7ec8cbf8d9b954d567ea50aafe0c83ff66768c3`
- Branch: `setup/agent-maintenance`
- Host: `HY-GPU`
- Slurm job ID: Not applicable; HY-GPU is a non-Slurm server
- Slurm request: Not applicable
- Environment: Conda environment `alphagenome`; A100 80GB GPUs; user explicitly allowed four-GPU work with `CUDA_VISIBLE_DEVICES=0,1,2,3`
- Command:

```bash
conda run -n alphagenome python -m py_compile \
  scripts/alphagenome_rna_seq11_adapter.py \
  scripts/alphagenome_rna_seq11_finetune.py \
  scripts/alphagenome_rna_seq11_eval.py \
  scripts/alphagenome_rna_seq11_scale_diagnostic.py

# Each smoke used one train example, one valid example, max_steps=1, and one GPU.
# Architectures checked:
# conv1x1, mlp1x1, conv3, conv5, residual-conv1x1, residual-conv3.
# Additional checks:
# linear-loss-type smooth-l1 and linear-target-space binned128-log1p-mean.

CUDA_VISIBLE_DEVICES=2 conda run -n alphagenome python -u \
  scripts/alphagenome_rna_seq11_eval.py \
  --dataset-dir alphagenome_custom/datasets/rna_seq_npz_valid \
  --weights weights/alphagenome_pytorch/model_all_folds.safetensors \
  --checkpoint runs/rna_seq11_adapter_formal_20260515_5000steps_lr3e-4_128bp/adapter_head_best.pt \
  --point-metrics \
  --spearman-sample-size 200000 \
  --spearman-seed 20260515 \
  --metrics-output runs/rna_seq11_sweep_legacy_valid_repro_20260522/legacy_valid_point_metrics.tsv \
  --diagnostic-output runs/rna_seq11_sweep_legacy_valid_repro_20260522/legacy_valid_diagnostic.tsv \
  --device auto
```

- Input data: `alphagenome_custom/datasets/rna_seq_npz_train` and `alphagenome_custom/datasets/rna_seq_npz_valid`; base weights `weights/alphagenome_pytorch/model_all_folds.safetensors`; legacy selected checkpoint `runs/rna_seq11_adapter_formal_20260515_5000steps_lr3e-4_128bp/adapter_head_best.pt`
- Output paths: ignored smoke directories under `runs/rna_seq11_sweep_smoke_*_20260522_1step`; ignored logs under `logs/`; legacy validation reproduction at `runs/rna_seq11_sweep_legacy_valid_repro_20260522`
- Result summary: Code compilation passed. All new linear head architectures completed 1-step train, 1-example valid, checkpoint save, and checkpoint reload. SmoothL1 and binned128 target-space smoke checks completed. The binned checkpoint common128 and scale-diagnostic paths were fixed and verified after catching an initial shape mismatch during smoke. Combined diagnostic output wrote the expected scopes.
- Verification:
  - 1-step prediction shapes were `1x11x1048576` for full-log1p heads and `1x11x8192` for the binned128 target-space smoke.
  - Trainable/head parameters were `33803` for `conv1x1`, `789515` for `mlp1x1`, `986379` for `conv3`, `1117451` for `conv5`, `33814` for `residual-conv1x1`, and `986390` for `residual-conv3`.
  - Legacy selected checkpoint full valid reproduction matched the recorded benchmark: MSE `1.0609935`, MAE `0.71866247`, Pearson `0.58809888`.
  - Full valid diagnostic output contained `11` track rows, `39` window rows, and `4` stratum rows.
  - Binned128 smoke common128 reload reported MSE `1.9483177` and Pearson `0.070606199` on one validation example; this is an environment/checkpoint reload check, not a model result.
- Failures or warnings: The first binned128 common128/scale diagnostic attempt exposed a shape mismatch caused by upsampling binned predictions before common128 evaluation; the code was corrected and the binned smoke was rerun successfully. These smoke-test metrics are environment/tooling validation only and should not be used for model selection. The held-out test split was not used. `git push` failed with `fatal: could not read Username for 'https://github.com': No such device or address`.
- Next actions: Launch the staged validation-only sweep from this tooling: seed-stability runs, 18-run short hyperparameter screen, head screen, binned target screen, then promotion runs only for valid-selected candidates.
- Claim status: verified

## 2026-05-21 - GenomeTracks 128 bp Warm-Start Count and Hybrid Loss Sweep

- Run type: analysis, training, and evaluation
- Purpose: Use validation-only diagnostics to decide whether the 128 bp GenomeTracksHead should be extended or promoted to 1 bp + 128 bp multi-resolution training after the from-scratch Poisson-multinomial objective underperformed.
- Git commit: `784789849e75cffb59af5f769ec1c92a311e5803`
- Branch: `setup/agent-maintenance`
- Host: `HY-GPU`
- Slurm job ID: Not applicable; HY-GPU is a non-Slurm server
- Slurm request: Not applicable
- Environment: Conda environment `alphagenome`; PyTorch `2.11.0+cu128`; A100 80GB GPUs; user explicitly allowed `CUDA_VISIBLE_DEVICES=0,1,2,3`
- Command:

```bash
# Scale/common-metric diagnostics for existing checkpoints.
CUDA_VISIBLE_DEVICES=0 conda run -n alphagenome python scripts/alphagenome_rna_seq11_scale_diagnostic.py \
  --dataset-dir alphagenome_custom/datasets/rna_seq_npz_valid \
  --checkpoint <checkpoint> \
  --weights weights/alphagenome_pytorch/model_all_folds.safetensors \
  --metrics-output <run_dir>/valid_scale_diagnostic.tsv \
  --batch-size 1 \
  --num-workers 2 \
  --device cuda

# Warm-start Poisson-multinomial run.
RUN_ID=rna_seq11_genometracks_128bp_pm_warmstart_pos1_20260521_500steps_lr3e-5_4gpu_gacc4
CUDA_VISIBLE_DEVICES=0,1,2,3 conda run -n alphagenome torchrun --standalone --nproc_per_node=4 \
  scripts/alphagenome_rna_seq11_finetune.py \
  --train-dataset-dir alphagenome_custom/datasets/rna_seq_npz_train \
  --valid-dataset-dir alphagenome_custom/datasets/rna_seq_npz_valid \
  --weights weights/alphagenome_pytorch/model_all_folds.safetensors \
  --output-dir runs/${RUN_ID} \
  --init-adapter-checkpoint runs/rna_seq11_genometracks_128bp_trainnonzero_20260520_5000steps_lr3e-4_4gpu_gacc4/adapter_head_best.pt \
  --embedding-resolution 128 \
  --head-type genome-tracks \
  --head-resolutions 128 \
  --target-transform none \
  --track-means-source train-nonzero \
  --loss-type poisson-multinomial \
  --multinomial-num-segments 8 \
  --positional-weight 1.0 \
  --count-weight 1.0 \
  --batch-size 1 \
  --grad-accum-steps 4 \
  --max-steps 500 \
  --eval-every 100 \
  --learning-rate 3e-5 \
  --warmup-steps 50 \
  --lr-schedule cosine \
  --grad-clip-norm 1.0 \
  --num-workers 2

# Warm-start hybrid sweep; repeated with poisson-weight 1e-5, 3e-5, and 1e-4.
RUN_ID=rna_seq11_genometracks_128bp_hybrid_pw<WEIGHT>_20260521_500steps_lr3e-5_4gpu_gacc4
CUDA_VISIBLE_DEVICES=0,1,2,3 conda run -n alphagenome torchrun --standalone --nproc_per_node=4 \
  scripts/alphagenome_rna_seq11_finetune.py \
  --train-dataset-dir alphagenome_custom/datasets/rna_seq_npz_train \
  --valid-dataset-dir alphagenome_custom/datasets/rna_seq_npz_valid \
  --weights weights/alphagenome_pytorch/model_all_folds.safetensors \
  --output-dir runs/${RUN_ID} \
  --init-adapter-checkpoint runs/rna_seq11_genometracks_128bp_trainnonzero_20260520_5000steps_lr3e-4_4gpu_gacc4/adapter_head_best.pt \
  --embedding-resolution 128 \
  --head-type genome-tracks \
  --head-resolutions 128 \
  --target-transform none \
  --track-means-source train-nonzero \
  --loss-type hybrid-mse-poisson \
  --multinomial-num-segments 8 \
  --positional-weight 1.0 \
  --count-weight 1.0 \
  --mse-weight 1.0 \
  --poisson-weight <WEIGHT> \
  --batch-size 1 \
  --grad-accum-steps 4 \
  --max-steps 500 \
  --eval-every 100 \
  --learning-rate 3e-5 \
  --warmup-steps 50 \
  --lr-schedule cosine \
  --grad-clip-norm 1.0 \
  --num-workers 2
```

- Input data: `alphagenome_custom/datasets/rna_seq_npz_train` and `alphagenome_custom/datasets/rna_seq_npz_valid`
- Output paths:
  - `runs/rna_seq11_genometracks_128bp_pm_warmstart_pos1_20260521_500steps_lr3e-5_4gpu_gacc4`
  - `runs/rna_seq11_genometracks_128bp_hybrid_pw1e-5_20260521_500steps_lr3e-5_4gpu_gacc4`
  - `runs/rna_seq11_genometracks_128bp_hybrid_pw3e-5_20260521_500steps_lr3e-5_4gpu_gacc4`
  - `runs/rna_seq11_genometracks_128bp_hybrid_pw1e-4_20260521_500steps_lr3e-5_4gpu_gacc4`
- Result summary: The diagnostics showed that improving raw signal mean scale did not improve the common 128 bp `log1p(mean raw)` metric. The legacy 128 bp linear adapter still had the best valid common metric. The 128 bp GenomeTracks MSE checkpoint remained the best GenomeTracks-family checkpoint in this round, and every warm-start count or hybrid loss candidate was worse than it.
- Verification:
  - Scale diagnostics before warm-starting: legacy 128 bp linear valid prediction/target mean ratio `0.10651676`, common128 MSE `1.0298867`, Pearson `0.62135428`; 128 bp GenomeTracks MSE ratio `0.53258544`, common128 MSE `2.133715`, Pearson `0.53935375`; from-scratch Poisson-multinomial ratio `0.65168605`, common128 MSE `3.0527705`, Pearson `0.49153418`.
  - Warm-start Poisson-multinomial valid losses were step 100 `18511.789`, 200 `18495.892`, 300 `18489.356`, 400 `18485.56`, and 500 `18485.124`. Its scale ratio was `0.78547391`, common128 MSE `2.9095493`, MAE `1.3788249`, Pearson `0.52790218`.
  - Hybrid `poisson_weight=1e-5` valid losses were step 100 `1.688949`, 200 `1.6883764`, 300 `1.6878839`, 400 `1.6877085`, and 500 `1.6876779`. Its scale ratio was `0.55806949`, common128 MSE `2.2140787`, MAE `1.1680553`, Pearson `0.53737598`.
  - Hybrid `poisson_weight=3e-5` valid losses were step 100 `2.0615844`, 200 `2.0610508`, 300 `2.0605806`, 400 `2.0604339`, and 500 `2.0604059`. Its scale ratio was `0.5863238`, common128 MSE `2.3104292`, MAE `1.2004244`, Pearson `0.5351673`.
  - Hybrid `poisson_weight=1e-4` valid losses were step 100 `3.3626393`, 200 `3.3621139`, 300 `3.3615108`, 400 `3.3613712`, and 500 `3.3613398`. Its scale ratio was `0.64649976`, common128 MSE `2.504247`, MAE `1.2616879`, Pearson `0.53141917`.
- Failures or warnings: The held-out test split was not used in this sweep. Objective-space losses are not directly comparable across different `poisson_weight` settings, so model selection used only the common valid metric and diagnostics. All warm-started count/hybrid candidates missed the pre-specified promotion target of improving on the 128 bp GenomeTracks MSE common metric.
- Next actions: Do not extend these count/hybrid 128 bp candidates and do not start 1 bp + 128 bp multi-resolution training from them. Keep the legacy 128 bp linear adapter as the selected model; revisit target scaling/loss design only after the expanded data are available or after a separate validation-only normalization study.
- Claim status: verified

## 2026-05-21 - GenomeTracks 128 bp MSE Run Stopped at Step 3000 and Common-Metric Audit

- Run type: training and evaluation
- Purpose: Stop the 4-GPU 128 bp-only GenomeTracksHead MSE run at the user's requested 3000-step cutoff, then compare its best checkpoint with the legacy 128 bp linear adapter best checkpoint using a common 128 bp binned `log1p(mean raw signal)` metric.
- Git commit: Training was launched from `65ae17776f7ffccfdfd744840104d9038352f279`; common-metric evaluation used `8c094c57c1c59122a356348d7f64be1d4f541d86`.
- Branch: `setup/agent-maintenance`
- Host: `HY-GPU`
- Slurm job ID: Not applicable; HY-GPU is a non-Slurm server
- Slurm request: Not applicable
- Environment: Conda environment `alphagenome`; PyTorch `2.11.0+cu128`; A100 80GB GPUs; user explicitly allowed `CUDA_VISIBLE_DEVICES=0,1,2,3`
- Command:

```bash
RUN_ID=rna_seq11_genometracks_128bp_trainnonzero_20260520_5000steps_lr3e-4_4gpu_gacc4

# Training command was launched on 2026-05-20 and stopped after valid step 3000.
CUDA_VISIBLE_DEVICES=0,1,2,3 conda run -n alphagenome torchrun --standalone --nproc_per_node=4 \
  scripts/alphagenome_rna_seq11_finetune.py \
  --train-dataset-dir alphagenome_custom/datasets/rna_seq_npz_train \
  --valid-dataset-dir alphagenome_custom/datasets/rna_seq_npz_valid \
  --weights weights/alphagenome_pytorch/model_all_folds.safetensors \
  --output-dir runs/${RUN_ID} \
  --embedding-resolution 128 \
  --head-type genome-tracks \
  --head-resolutions 128 \
  --target-transform none \
  --track-means-source train-nonzero \
  --batch-size 1 \
  --grad-accum-steps 4 \
  --max-steps 5000 \
  --eval-every 250 \
  --learning-rate 3e-4 \
  --warmup-steps 500 \
  --lr-schedule cosine \
  --grad-clip-norm 1.0 \
  --num-workers 2

CUDA_VISIBLE_DEVICES=0 conda run -n alphagenome python scripts/alphagenome_rna_seq11_eval.py \
  --dataset-dir alphagenome_custom/datasets/rna_seq_npz_valid \
  --checkpoint runs/${RUN_ID}/adapter_head_best.pt \
  --common-128bp-metrics log1p-mean \
  --metrics-output runs/${RUN_ID}/valid_best_common128_log1pmean_metrics.tsv \
  --device cuda

CUDA_VISIBLE_DEVICES=1 conda run -n alphagenome python scripts/alphagenome_rna_seq11_eval.py \
  --dataset-dir alphagenome_custom/datasets/rna_seq_npz_test \
  --checkpoint runs/${RUN_ID}/adapter_head_best.pt \
  --common-128bp-metrics log1p-mean \
  --metrics-output runs/${RUN_ID}/test_best_common128_log1pmean_metrics.tsv \
  --device cuda
```

- Input data: `alphagenome_custom/datasets/rna_seq_npz_train`, `alphagenome_custom/datasets/rna_seq_npz_valid`, and audit-only `alphagenome_custom/datasets/rna_seq_npz_test`
- Output path: `runs/rna_seq11_genometracks_128bp_trainnonzero_20260520_5000steps_lr3e-4_4gpu_gacc4`
- Result summary: Training was stopped after the step-3000 validation. GenomeTracks model-space valid MSE improved monotonically from `1.808731` at step 250 to `1.5025958` at step 3000. On the common 128 bp binned `log1p(mean raw)` metric, however, this checkpoint was worse than the legacy 128 bp linear best checkpoint.
- Verification: GenomeTracks MSE common metrics were valid MSE `2.133715`, MAE `1.1418334`, Pearson `0.53935375`; test MSE `1.8786616`, MAE `1.0444137`, Pearson `0.53915857`. Legacy 128 bp linear best common metrics were valid MSE `1.0298867`, MAE `0.68138634`, Pearson `0.62135428`; test MSE `1.1380058`, MAE `0.73726858`, Pearson `0.59136682`.
- Failures or warnings: The held-out test split was used only for this audit comparison, not for model selection. GenomeTracks model-space MSE is not directly comparable to legacy log1p MSE; the common 128 bp metric was added for cross-head comparison.
- Next actions: Do not promote the MSE-trained 128 bp GenomeTracks checkpoint to 1 bp + 128 bp. First test the AlphaGenome-style count/position loss on 128 bp only.
- Claim status: verified

## 2026-05-21 - GenomeTracks 128 bp Poisson-Multinomial Loss 500-Step Sanity

- Run type: smoke test, training, and evaluation
- Purpose: Evaluate whether an AlphaGenome-style RNA-seq loss, implemented as Poisson total-count plus positional multinomial count loss over 8 sequence segments, is a better 128 bp-only GenomeTracksHead objective before attempting 1 bp + 128 bp multi-resolution training.
- Git commit: `8c094c57c1c59122a356348d7f64be1d4f541d86`
- Branch: `setup/agent-maintenance`
- Host: `HY-GPU`
- Slurm job ID: Not applicable; HY-GPU is a non-Slurm server
- Slurm request: Not applicable
- Environment: Conda environment `alphagenome`; PyTorch `2.11.0+cu128`; A100 80GB GPUs; user explicitly allowed `CUDA_VISIBLE_DEVICES=0,1,2,3`
- Command:

```bash
RUN_ID=rna_seq11_genometracks_128bp_poissonmulti_20260521_500steps_lr3e-4_4gpu_gacc4_run2

CUDA_VISIBLE_DEVICES=0,1,2,3 conda run -n alphagenome torchrun --standalone --nproc_per_node=4 \
  scripts/alphagenome_rna_seq11_finetune.py \
  --train-dataset-dir alphagenome_custom/datasets/rna_seq_npz_train \
  --valid-dataset-dir alphagenome_custom/datasets/rna_seq_npz_valid \
  --weights weights/alphagenome_pytorch/model_all_folds.safetensors \
  --output-dir runs/${RUN_ID} \
  --embedding-resolution 128 \
  --head-type genome-tracks \
  --head-resolutions 128 \
  --target-transform none \
  --track-means-source train-nonzero \
  --loss-type poisson-multinomial \
  --multinomial-num-segments 8 \
  --positional-weight 5.0 \
  --count-weight 1.0 \
  --batch-size 1 \
  --grad-accum-steps 4 \
  --max-steps 500 \
  --eval-every 100 \
  --learning-rate 3e-4 \
  --warmup-steps 50 \
  --lr-schedule cosine \
  --grad-clip-norm 1.0 \
  --num-workers 2

CUDA_VISIBLE_DEVICES=0 conda run -n alphagenome python scripts/alphagenome_rna_seq11_eval.py \
  --dataset-dir alphagenome_custom/datasets/rna_seq_npz_valid \
  --checkpoint runs/${RUN_ID}/adapter_head_best.pt \
  --common-128bp-metrics log1p-mean \
  --metrics-output runs/${RUN_ID}/valid_best_common128_log1pmean_metrics.tsv \
  --device cuda

CUDA_VISIBLE_DEVICES=1 conda run -n alphagenome python scripts/alphagenome_rna_seq11_eval.py \
  --dataset-dir alphagenome_custom/datasets/rna_seq_npz_test \
  --checkpoint runs/${RUN_ID}/adapter_head_best.pt \
  --common-128bp-metrics log1p-mean \
  --metrics-output runs/${RUN_ID}/test_best_common128_log1pmean_metrics.tsv \
  --device cuda
```

- Input data: `alphagenome_custom/datasets/rna_seq_npz_train`, `alphagenome_custom/datasets/rna_seq_npz_valid`, and audit-only `alphagenome_custom/datasets/rna_seq_npz_test`
- Output path: `runs/rna_seq11_genometracks_128bp_poissonmulti_20260521_500steps_lr3e-4_4gpu_gacc4_run2`
- Result summary: The count/position objective ran successfully and valid loss decreased from `91680.721` at step 100 to `85025.399` at step 500, but the common 128 bp `log1p(mean raw)` metrics were substantially worse than both the legacy 128 bp linear best checkpoint and the MSE-trained GenomeTracks checkpoint.
- Verification: Poisson-multinomial valid losses were step 100 `91680.721`, step 200 `86910.675`, step 300 `85539.633`, step 400 `85080.233`, and step 500 `85025.399`. Audit test Poisson-multinomial loss was `54648.254`. Common 128 bp `log1p(mean raw)` metrics were valid MSE `3.0527705`, MAE `1.4647581`, Pearson `0.49153418`; test MSE `2.5278237`, MAE `1.2957843`, Pearson `0.5309488`.
- Failures or warnings: A first background `nohup conda run` launch at `runs/rna_seq11_genometracks_128bp_poissonmulti_20260521_500steps_lr3e-4_4gpu_gacc4` exited after writing only `config.json`; the log was empty and no metrics were written. A 1-step smoke test and a 2-step full-data debug run completed successfully before the foreground 500-step run. Test metrics are audit-only and should not guide model selection.
- Next actions: Do not start 1 bp + 128 bp multi-resolution training from this objective yet. Investigate loss weighting/scaling, warm-starting from the MSE-trained GenomeTracks checkpoint, or a hybrid objective before spending 1 bp GPU time.
- Claim status: verified

## 2026-05-19 - GenomeTracks RNA-seq 128 bp Head One-Step Sanity on HY-GPU GPU 2

- Run type: smoke test
- Purpose: Verify that the new AlphaGenome-style RNA-seq `GenomeTracksHead` path can run a 128 bp-only forward/backward pass on the existing C. elegans train/valid NPZ inputs using raw targets scaled into model space, while keeping the frozen trunk unchanged and keeping the held-out test split untouched.
- Git commit: `fcc3ea7c290ffdcba4005cd12b133d094d8fb23b`; working tree contained uncommitted GenomeTracksHead implementation updates.
- Branch: `setup/agent-maintenance`
- Host: `HY-GPU`
- Slurm job ID: Not applicable; HY-GPU is a non-Slurm server
- Slurm request: Not applicable
- Environment: Conda environment `alphagenome`; Python `3.12.13`; PyTorch `2.11.0+cu128`; PyTorch CUDA `12.8`; driver CUDA `13.2`; `alphagenome_pytorch 0.3.1`; `CUDA_VISIBLE_DEVICES=2`
- Command:

```bash
cd /home/zelinli6/Alphagenome

RUN_DIR=runs/rna_seq11_genometracks_128bp_head_sanity_20260519_1step
mkdir -p "$RUN_DIR"

CUDA_VISIBLE_DEVICES=2 conda run -n alphagenome \
  python -u scripts/alphagenome_rna_seq11_finetune.py \
  --train-dataset-dir alphagenome_custom/datasets/rna_seq_npz_train \
  --valid-dataset-dir alphagenome_custom/datasets/rna_seq_npz_valid \
  --weights weights/alphagenome_pytorch/model_all_folds.safetensors \
  --output-dir "$RUN_DIR" \
  --batch-size 1 \
  --max-steps 1 \
  --eval-every 1 \
  --max-train-examples 1 \
  --max-valid-examples 1 \
  --learning-rate 1e-4 \
  --embedding-resolution 128 \
  --head-type genome-tracks \
  --head-resolutions 128 \
  --target-transform none \
  --track-means-source grouped-qc \
  --seed 20260519 \
  --device auto \
  --no-save-checkpoint 2>&1 | tee "$RUN_DIR/train.log"
```

- Input data: First train example from `alphagenome_custom/datasets/rna_seq_npz_train`; first validation example from `alphagenome_custom/datasets/rna_seq_npz_valid`; base weights `weights/alphagenome_pytorch/model_all_folds.safetensors`; grouped track means from `alphagenome_custom/metadata/grouped_bigwig_qc.tsv`. The held-out test split was not used.
- Output path: `runs/rna_seq11_genometracks_128bp_head_sanity_20260519_1step/`, containing `config.json`, `metrics.tsv`, and `train.log`. Checkpoint saving was disabled. This path is ignored by Git.
- Result summary: Completed successfully. The 128 bp GenomeTracks-style RNA-seq head produced model-space predictions with shape `128:1x11x8192`, used raw targets with `--target-transform none`, and reported model-space train loss `2.0833006` and one-example validation loss `2.7360559`.
- Verification: Reported `base_has_grad=False`, `head_weight_grad_norm=2.45218`, `head_parameters=33814`, `head_trainable_parameters=33814`, `base_non_lora_trainable_parameters=0`, and peak CUDA allocation `17585.7` MB. Post-run `nvidia-smi` showed no active GPU compute process.
- Failures or warnings: This is a one-step engineering sanity check on one train example and one validation example, not a validation-set model result. The loss is in GenomeTracksHead model-scaled space and is not directly comparable to the earlier log1p MSE values. The held-out test split was not used.
- Next actions: If this path is pursued, run a validation-only 128 bp pilot with full validation, then test the 1 bp + 128 bp multi-resolution head only after the 128 bp-only path is stable.
- Claim status: verified

## 2026-05-16 - AlphaGenome 128 bp Last-Block LoRA Scheme C 100-Step Pilots on HY-GPU GPU 2

- Run type: training
- Purpose: Test scheme C for lightweight trunk adaptation: initialize from the selected 5000-step 128 bp RNA-seq adapter head, then compare LoRA-only continuation with the head frozen against LoRA-plus-head continuation, while keeping original non-LoRA AlphaGenome trunk weights frozen and keeping the held-out test split untouched.
- Git commit: `b592b0426b158e4a13e8afd8f2c96eae86692423`; working tree contained uncommitted LoRA training-script updates and documentation updates.
- Branch: `setup/agent-maintenance`
- Host: `HY-GPU`
- Slurm job ID: Not applicable; HY-GPU is a non-Slurm server
- Slurm request: Not applicable
- Environment: Conda environment `alphagenome`; Python `3.12.13`; PyTorch `2.11.0+cu128`; PyTorch CUDA `12.8`; driver CUDA `13.2`; `alphagenome_pytorch 0.3.1`; `CUDA_VISIBLE_DEVICES=2`
- Command:

```bash
cd /home/zelinli6/Alphagenome

RUN_DIR=runs/rna_seq11_adapter_128bp_lora_phaseA_20260516_100steps_lr5e-5_freezehead
mkdir -p "$RUN_DIR"

CUDA_VISIBLE_DEVICES=2 conda run -n alphagenome \
  python -u scripts/alphagenome_rna_seq11_finetune.py \
  --train-dataset-dir alphagenome_custom/datasets/rna_seq_npz_train \
  --valid-dataset-dir alphagenome_custom/datasets/rna_seq_npz_valid \
  --weights weights/alphagenome_pytorch/model_all_folds.safetensors \
  --output-dir "$RUN_DIR" \
  --init-adapter-checkpoint runs/rna_seq11_adapter_formal_20260515_5000steps_lr3e-4_128bp/adapter_head_best.pt \
  --batch-size 1 \
  --max-steps 100 \
  --eval-every 20 \
  --learning-rate 5e-5 \
  --embedding-resolution 128 \
  --seed 20260516 \
  --grad-accum-steps 1 \
  --grad-clip-norm 1.0 \
  --lr-schedule cosine \
  --warmup-steps 20 \
  --enable-last-block-lora \
  --lora-rank 4 \
  --lora-alpha 8 \
  --freeze-head \
  --device auto 2>&1 | tee "$RUN_DIR/train.log"

RUN_DIR=runs/rna_seq11_adapter_128bp_lora_phaseB_20260516_100steps_lr5e-5_lora_head
mkdir -p "$RUN_DIR"

CUDA_VISIBLE_DEVICES=2 conda run -n alphagenome \
  python -u scripts/alphagenome_rna_seq11_finetune.py \
  --train-dataset-dir alphagenome_custom/datasets/rna_seq_npz_train \
  --valid-dataset-dir alphagenome_custom/datasets/rna_seq_npz_valid \
  --weights weights/alphagenome_pytorch/model_all_folds.safetensors \
  --output-dir "$RUN_DIR" \
  --init-adapter-checkpoint runs/rna_seq11_adapter_formal_20260515_5000steps_lr3e-4_128bp/adapter_head_best.pt \
  --batch-size 1 \
  --max-steps 100 \
  --eval-every 20 \
  --learning-rate 5e-5 \
  --embedding-resolution 128 \
  --seed 20260516 \
  --grad-accum-steps 1 \
  --grad-clip-norm 1.0 \
  --lr-schedule cosine \
  --warmup-steps 20 \
  --enable-last-block-lora \
  --lora-rank 4 \
  --lora-alpha 8 \
  --device auto 2>&1 | tee "$RUN_DIR/train.log"
```

- Input data: `alphagenome_custom/datasets/rna_seq_npz_train` (`116` examples); `alphagenome_custom/datasets/rna_seq_npz_valid` (`39` examples); base weights `weights/alphagenome_pytorch/model_all_folds.safetensors`; initialization checkpoint `runs/rna_seq11_adapter_formal_20260515_5000steps_lr3e-4_128bp/adapter_head_best.pt`. The held-out test split was not used.
- Output path: `runs/rna_seq11_adapter_128bp_lora_phaseA_20260516_100steps_lr5e-5_freezehead/` and `runs/rna_seq11_adapter_128bp_lora_phaseB_20260516_100steps_lr5e-5_lora_head/`, each containing `config.json`, `metrics.tsv`, `adapter_head.pt`, and `adapter_head_best.pt`. These paths are ignored by Git.
- Result summary: Both scheme C pilots completed successfully and confirmed gradients only through the intended trainable modules. Phase A trained only the final-block LoRA parameters with the previously selected adapter head frozen. Phase B trained final-block LoRA parameters plus the adapter head. Neither pilot improved on the selected frozen-trunk 128 bp adapter validation MSE of `1.0609935`.
- Verification: Phase A reported `head_trainable_parameters=0`, `lora_trainable_parameters=72960`, and `optimizer_trainable_parameters=72960`; full-validation MSEs were `1.0613562` at step 20, `1.0625891` at step 40, `1.0640197` at step 60, `1.0650645` at step 80, and `1.0653593` at step 100. Phase B reported `head_trainable_parameters=33803`, `lora_trainable_parameters=72960`, and `optimizer_trainable_parameters=106763`; full-validation MSEs were `1.0698467` at step 20, `1.067338` at step 40, `1.0680136` at step 60, `1.0715059` at step 80, and `1.0727408` at step 100. Both runs used `embedding_resolution=128`, LoRA rank `4`, alpha `8`, learning rate `5e-5`, 20-step warmup, cosine decay, and gradient clipping at norm `1.0`. Peak CUDA allocation was about `17820.6` to `17820.9` MB.
- Failures or warnings: These are validation-set continuation pilots, not held-out test results. Phase A was numerically very close to the original selected checkpoint but did not beat it; Phase B degraded more clearly. The held-out test split was not used.
- Next actions: Keep the selected model unchanged unless a future validation-only run beats `1.0609935`. If continuing LoRA exploration, try a lower continuation learning rate such as `1e-5` to `2e-5`, or a more conservative schedule, and select only by validation performance.
- Claim status: verified

## 2026-05-16 - AlphaGenome 128 bp Last-Block LoRA Adapter Feasibility Smoke on HY-GPU GPU 2

- Run type: smoke test
- Purpose: Verify that a parameter-efficient LoRA adaptation path can train the custom C. elegans 11-track RNA-seq adapter while keeping original AlphaGenome trunk weights frozen.
- Git commit: `b592b0426b158e4a13e8afd8f2c96eae86692423`; working tree contained uncommitted DDP/scheduler updates from the previous 1 bp pilot plus the new LoRA smoke script and adapter gradient-control update.
- Branch: `setup/agent-maintenance`
- Host: `HY-GPU`
- Slurm job ID: Not applicable; HY-GPU is a non-Slurm server
- Slurm request: Not applicable
- Environment: Conda environment `alphagenome`; Python `3.12.13`; PyTorch `2.11.0+cu128`; PyTorch CUDA `12.8`; driver CUDA `13.2`; `alphagenome_pytorch 0.3.1`; `CUDA_VISIBLE_DEVICES=2`
- Command:

```bash
cd /home/zelinli6/Alphagenome

RUN_DIR=runs/rna_seq11_adapter_128bp_lora_lastblock_smoke_20260516_1step_v2
mkdir -p "$RUN_DIR"

CUDA_VISIBLE_DEVICES=2 conda run -n alphagenome \
  python -u scripts/alphagenome_rna_seq11_lora_smoke.py \
  --dataset-dir alphagenome_custom/datasets/rna_seq_npz_train \
  --valid-dataset-dir alphagenome_custom/datasets/rna_seq_npz_valid \
  --weights weights/alphagenome_pytorch/model_all_folds.safetensors \
  --batch-size 1 \
  --max-steps 1 \
  --max-train-examples 1 \
  --max-valid-examples 1 \
  --max-valid-batches 1 \
  --learning-rate 1e-4 \
  --lora-rank 4 \
  --lora-alpha 8 \
  --device auto 2>&1 | tee "$RUN_DIR/train.log"
```

- Input data: First train example from `alphagenome_custom/datasets/rna_seq_npz_train`; first validation example from `alphagenome_custom/datasets/rna_seq_npz_valid`; base weights `weights/alphagenome_pytorch/model_all_folds.safetensors`. The held-out test split was not used.
- Output path: `runs/rna_seq11_adapter_128bp_lora_lastblock_smoke_20260516_1step_v2/train.log`. This path is ignored by Git.
- Result summary: Completed successfully. LoRA was applied only to the final transformer block sequence modules `tower.blocks.8.mha` and `tower.blocks.8.mlp`, while original non-LoRA AlphaGenome trunk parameters remained frozen.
- Verification: Observed `lora_modules_applied=6`: `tower.blocks.8.mha.q_proj`, `tower.blocks.8.mha.k_proj`, `tower.blocks.8.mha.v_proj`, `tower.blocks.8.mha.linear_embedding`, `tower.blocks.8.mlp.fc1`, and `tower.blocks.8.mlp.fc2`. Observed `lora_trainable_parameters=72960`, `head_parameters=33803`, `base_non_lora_trainable_parameters=0`, and `total_trainable_parameters=106763`. One training step produced `prediction_shape=1x11x1048576`, train loss `2.83829`, `base_non_lora_has_grad=False`, `lora_has_grad=True`, `head_has_grad=True`, LoRA grad norm `0.0144359`, head grad norm `5.5984`, peak CUDA allocation `17767.9` MB, and one-example valid loss `1.91089`. Post-run check showed no active GPU compute process.
- Failures or warnings: This is a feasibility smoke test on one train example and one validation example, not a model-performance result. It does not support model selection. The full held-out test split was not used.
- Next actions: If continuing, run a validation-only 128 bp LoRA pilot, for example 100 to 500 steps with full validation, lower learning rate, warmup/cosine, and best-validation checkpointing. Continue keeping test untouched.
- Claim status: verified

## 2026-05-16 - AlphaGenome 128 bp Last-Block LoRA Adapter Smoke Failure on HY-GPU GPU 2

- Run type: smoke test
- Purpose: First attempt at the 128 bp last-block LoRA feasibility smoke.
- Git commit: `b592b0426b158e4a13e8afd8f2c96eae86692423`; working tree contained uncommitted LoRA smoke edits.
- Branch: `setup/agent-maintenance`
- Host: `HY-GPU`
- Slurm job ID: Not applicable; HY-GPU is a non-Slurm server
- Slurm request: Not applicable
- Environment: Conda environment `alphagenome`; `CUDA_VISIBLE_DEVICES=2`
- Command:

```bash
cd /home/zelinli6/Alphagenome

RUN_DIR=runs/rna_seq11_adapter_128bp_lora_lastblock_smoke_20260516_1step
mkdir -p "$RUN_DIR"

CUDA_VISIBLE_DEVICES=2 conda run -n alphagenome \
  python -u scripts/alphagenome_rna_seq11_lora_smoke.py \
  --dataset-dir alphagenome_custom/datasets/rna_seq_npz_train \
  --valid-dataset-dir alphagenome_custom/datasets/rna_seq_npz_valid \
  --weights weights/alphagenome_pytorch/model_all_folds.safetensors \
  --batch-size 1 \
  --max-steps 1 \
  --max-train-examples 1 \
  --max-valid-examples 1 \
  --max-valid-batches 1 \
  --learning-rate 1e-4 \
  --lora-rank 4 \
  --lora-alpha 8 \
  --device auto 2>&1 | tee "$RUN_DIR/train.log"
```

- Input data: First train example and first validation example from the existing NPZ datasets; base weights `weights/alphagenome_pytorch/model_all_folds.safetensors`. The held-out test split was not used.
- Output path: `runs/rna_seq11_adapter_128bp_lora_lastblock_smoke_20260516_1step/train.log`. This path is ignored by Git.
- Result summary: Failed before completing the first forward pass.
- Verification: Module targeting worked before the failure: the script reported `lora_modules_applied=6`, `lora_trainable_parameters=72960`, `head_parameters=33803`, `base_non_lora_trainable_parameters=0`, and `total_trainable_parameters=106763`.
- Failures or warnings: Runtime error reported CPU/CUDA tensor mismatch because LoRA modules were inserted after the base model had already been moved to CUDA, leaving newly created LoRA layers on CPU. The script was fixed by moving the model to the selected device after applying LoRA, and the follow-up run succeeded.
- Next actions: See the successful `2026-05-16 - AlphaGenome 128 bp Last-Block LoRA Adapter Feasibility Smoke on HY-GPU GPU 2` entry.
- Claim status: verified

## 2026-05-15 - Frozen AlphaGenome 1 bp Adapter 500-Step DDP Pilot on HY-GPU GPUs 2 and 3

- Run type: training
- Purpose: Test a more stable 1 bp embedding adapter schedule using frozen AlphaGenome trunk, two-GPU DDP, lower learning rate, warmup, cosine decay, gradient clipping, and gradient accumulation, while keeping the held-out test split untouched.
- Git commit: `b592b0426b158e4a13e8afd8f2c96eae86692423`; working tree contained uncommitted training-script DDP/scheduler updates plus experiment-log updates.
- Branch: `setup/agent-maintenance`
- Host: `HY-GPU`
- Slurm job ID: Not applicable; HY-GPU is a non-Slurm server
- Slurm request: Not applicable
- Environment: Conda environment `alphagenome`; Python `3.12.13`; PyTorch `2.11.0+cu128`; PyTorch CUDA `12.8`; driver CUDA `13.2`; `alphagenome_pytorch 0.3.1`; `CUDA_VISIBLE_DEVICES=2,3`; `torch.distributed.run --nproc_per_node=2`
- Command:

```bash
cd /home/zelinli6/Alphagenome

RUN_DIR=runs/rna_seq11_adapter_1bp_ddp_pilot_20260515_500steps_lr1e-4_warmup50_accum4_clip1
mkdir -p "$RUN_DIR"

CUDA_VISIBLE_DEVICES=2,3 conda run -n alphagenome \
  python -m torch.distributed.run --standalone --nproc_per_node=2 \
  scripts/alphagenome_rna_seq11_finetune.py \
  --train-dataset-dir alphagenome_custom/datasets/rna_seq_npz_train \
  --valid-dataset-dir alphagenome_custom/datasets/rna_seq_npz_valid \
  --weights weights/alphagenome_pytorch/model_all_folds.safetensors \
  --output-dir "$RUN_DIR" \
  --batch-size 1 \
  --max-steps 500 \
  --eval-every 50 \
  --learning-rate 1e-4 \
  --embedding-resolution 1 \
  --seed 20260515 \
  --grad-accum-steps 4 \
  --grad-clip-norm 1.0 \
  --lr-schedule cosine \
  --warmup-steps 50 \
  --device auto 2>&1 | tee "$RUN_DIR/train.log"
```

- Input data: `alphagenome_custom/datasets/rna_seq_npz_train` (`116` examples); `alphagenome_custom/datasets/rna_seq_npz_valid` (`39` examples); base weights `weights/alphagenome_pytorch/model_all_folds.safetensors`. The held-out test split was not used.
- Output path: `runs/rna_seq11_adapter_1bp_ddp_pilot_20260515_500steps_lr1e-4_warmup50_accum4_clip1/` containing `config.json`, `metrics.tsv`, `train.log`, `adapter_head.pt`, and `adapter_head_best.pt`. This path is ignored by Git.
- Result summary: Completed successfully with two DDP ranks. The run used frozen trunk plus the 1 bp adapter head (`16907` trainable parameters), per-process batch size `1`, gradient accumulation `4`, and global effective batch size `8`. Validation MSE improved throughout the run and the best checkpoint was the final step.
- Verification: Observed `distributed=True`, `world_size=2`, `global_effective_batch_size=8`, repeated `prediction_shape=1x11x1048576`, `base_has_grad=False`, and peak torch CUDA allocation about `36039.0` MB on rank 0. Full-validation MSEs were `1.6729847` at step 50, `1.6765932` at step 100, `1.5817027` at step 150, `1.5381760` at step 200, `1.5135129` at step 250, `1.4926636` at step 300, `1.4850671` at step 350, `1.4796252` at step 400, `1.4776338` at step 450, and `1.4771859` at step 500. Best validation MSE was `1.4771859` at step 500, saved to `adapter_head_best.pt`. Post-run check showed no active GPU compute process.
- Failures or warnings: This is a validation-set pilot, not held-out test performance. The run prints a PyTorch warning that `OMP_NUM_THREADS` is set to `1` by `torch.distributed.run`; no training failure resulted. The final cosine learning rate reached `0`, so continuing this exact run without changing the schedule would not further update weights.
- Next actions: Evaluate the selected 1 bp checkpoint on the validation split with per-track MSE, MAE, Pearson, and sampled Spearman before deciding whether to run a longer 1 bp schedule or return to the stronger 128 bp adapter path. Do not use the held-out test split for this decision.
- Claim status: verified

## 2026-05-15 - Frozen AlphaGenome 1 bp Adapter DDP One-Step Smoke on HY-GPU GPUs 2 and 3

- Run type: smoke test
- Purpose: Verify the new two-process DDP training path before launching the longer 1 bp adapter pilot.
- Git commit: `b592b0426b158e4a13e8afd8f2c96eae86692423`; working tree contained uncommitted DDP/scheduler updates to the training script.
- Branch: `setup/agent-maintenance`
- Host: `HY-GPU`
- Slurm job ID: Not applicable; HY-GPU is a non-Slurm server
- Slurm request: Not applicable
- Environment: Conda environment `alphagenome`; `CUDA_VISIBLE_DEVICES=2,3`; `torch.distributed.run --nproc_per_node=2`
- Command:

```bash
cd /home/zelinli6/Alphagenome

RUN_DIR=runs/rna_seq11_adapter_1bp_ddp_smoke_20260515_1step
mkdir -p "$RUN_DIR"

CUDA_VISIBLE_DEVICES=2,3 conda run -n alphagenome \
  python -m torch.distributed.run --standalone --nproc_per_node=2 \
  scripts/alphagenome_rna_seq11_finetune.py \
  --train-dataset-dir alphagenome_custom/datasets/rna_seq_npz_train \
  --valid-dataset-dir alphagenome_custom/datasets/rna_seq_npz_valid \
  --weights weights/alphagenome_pytorch/model_all_folds.safetensors \
  --output-dir "$RUN_DIR" \
  --batch-size 1 \
  --max-steps 1 \
  --eval-every 1 \
  --max-train-examples 2 \
  --max-valid-examples 1 \
  --learning-rate 1e-4 \
  --embedding-resolution 1 \
  --seed 20260515 \
  --grad-accum-steps 4 \
  --grad-clip-norm 1.0 \
  --lr-schedule cosine \
  --warmup-steps 50 \
  --device auto \
  --no-save-checkpoint 2>&1 | tee "$RUN_DIR/train.log"
```

- Input data: First two train examples and first validation example from the existing NPZ datasets; base weights `weights/alphagenome_pytorch/model_all_folds.safetensors`. The held-out test split was not used.
- Output path: `runs/rna_seq11_adapter_1bp_ddp_smoke_20260515_1step/` containing `config.json`, `metrics.tsv`, and `train.log`. This path is ignored by Git.
- Result summary: Completed successfully. The DDP path reported `distributed=True`, `world_size=2`, per-process batch size `1`, and global effective batch size `8`.
- Verification: Observed `embedding_resolution=1`, `trainable_parameters=16907`, `prediction_shape=1x11x1048576`, `base_has_grad=False`, one-step train loss `2.55981`, one-example valid loss `1.89307`, and peak rank-0 CUDA allocation about `35934.9` MB.
- Failures or warnings: This is an environment and DDP plumbing smoke test, not a model result. PyTorch emitted the expected `OMP_NUM_THREADS` notice for `torch.distributed.run`.
- Next actions: Launch the planned 500-step validation-only DDP pilot with all train and validation NPZ examples.
- Claim status: verified

## 2026-05-15 - Frozen AlphaGenome 1 bp Embedding Adapter 100-Step Pilot on HY-GPU GPU 2

- Run type: training
- Purpose: Run a short validation-set pilot for the 1 bp embedding adapter path after the one-step feasibility smoke test, while keeping the held-out test split untouched for this new experiment round.
- Git commit: `b592b0426b158e4a13e8afd8f2c96eae86692423`; working tree contained experiment-log updates for the 1 bp smoke and this pilot.
- Branch: `setup/agent-maintenance`
- Host: `HY-GPU`
- Slurm job ID: Not applicable; HY-GPU is a non-Slurm server
- Slurm request: Not applicable
- Environment: Conda environment `alphagenome`; Python `3.12.13`; PyTorch `2.11.0+cu128`; PyTorch CUDA `12.8`; driver CUDA `13.2`; `alphagenome_pytorch 0.3.1`; `CUDA_VISIBLE_DEVICES=2`
- Command:

```bash
cd /home/zelinli6/Alphagenome

mkdir -p runs/rna_seq11_adapter_1bp_pilot_20260515_100steps_lr3e-4

CUDA_VISIBLE_DEVICES=2 conda run -n alphagenome \
  python -u scripts/alphagenome_rna_seq11_finetune.py \
  --train-dataset-dir alphagenome_custom/datasets/rna_seq_npz_train \
  --valid-dataset-dir alphagenome_custom/datasets/rna_seq_npz_valid \
  --weights weights/alphagenome_pytorch/model_all_folds.safetensors \
  --output-dir runs/rna_seq11_adapter_1bp_pilot_20260515_100steps_lr3e-4 \
  --batch-size 1 \
  --max-steps 100 \
  --eval-every 20 \
  --learning-rate 3e-4 \
  --embedding-resolution 1 \
  --seed 20260515 \
  --device auto \
  > runs/rna_seq11_adapter_1bp_pilot_20260515_100steps_lr3e-4/train.log 2>&1
```

- Input data: `alphagenome_custom/datasets/rna_seq_npz_train` (`116` examples); `alphagenome_custom/datasets/rna_seq_npz_valid` (`39` examples); base weights `weights/alphagenome_pytorch/model_all_folds.safetensors`. The held-out test split was not used.
- Output path: `runs/rna_seq11_adapter_1bp_pilot_20260515_100steps_lr3e-4/` containing `config.json` (`936` bytes), `metrics.tsv` (`4.1K`), `train.log` (`22K`), `adapter_head.pt` (`69K`), and `adapter_head_best.pt` (`69K`). This path is ignored by Git.
- Result summary: Completed successfully on a single allowed GPU (`CUDA_VISIBLE_DEVICES=2`). The 1 bp embedding adapter completed 100 training steps and full-validation passes every 20 steps.
- Verification: Observed `embedding_resolution=1`, `trainable_parameters=16907`, `head_parameters=16907`, repeated `prediction_shape=1x11x1048576`, and `base_has_grad=False`. Full-validation MSEs were `1.6738151` at step 20, `1.7223911` at step 40, `1.6576294` at step 60, `1.5281759` at step 80, and final-step `1.5943006` at step 100. Best validation MSE was `1.5281759` at step 80, saved to `adapter_head_best.pt`. Peak CUDA memory allocation was about `36038.9` MB, while `nvidia-smi` process memory was about `60911` MiB during the run. Post-run check showed no active GPU compute process.
- Failures or warnings: This is a short validation-set pilot, not held-out test performance. The 1 bp adapter path was computationally much slower than the 128 bp adapter path and did not outperform the prior 128 bp adapter pilots at comparable early steps. The held-out test split was not used.
- Next actions: Treat the current 1 bp result as feasible but not yet better. If continuing, tune only on validation data, for example a lower learning rate or longer run, and do not reuse the held-out test for selection.
- Claim status: verified

## 2026-05-15 - Frozen AlphaGenome 1 bp Embedding Adapter Feasibility Smoke on HY-GPU GPU 2

- Run type: smoke test
- Purpose: Check whether the frozen AlphaGenome 1 bp embedding path can run a minimal C. elegans 11-track RNA-seq adapter forward/backward step and a one-example validation pass on HY-GPU before attempting any longer 1 bp-resolution adapter experiment.
- Git commit: `b592b0426b158e4a13e8afd8f2c96eae86692423`; working tree contained this new experiment-log update after the pushed adapter workflow commit.
- Branch: `setup/agent-maintenance`
- Host: `HY-GPU`
- Slurm job ID: Not applicable; HY-GPU is a non-Slurm server
- Slurm request: Not applicable
- Environment: Conda environment `alphagenome`; Python `3.12.13`; PyTorch `2.11.0+cu128`; PyTorch CUDA `12.8`; driver CUDA `13.2`; `alphagenome_pytorch 0.3.1`; `CUDA_VISIBLE_DEVICES=2`
- Command:

```bash
cd /home/zelinli6/Alphagenome

mkdir -p runs/rna_seq11_adapter_1bp_smoke_20260515_1step

CUDA_VISIBLE_DEVICES=2 conda run -n alphagenome \
  python -u scripts/alphagenome_rna_seq11_finetune.py \
  --train-dataset-dir alphagenome_custom/datasets/rna_seq_npz_train \
  --valid-dataset-dir alphagenome_custom/datasets/rna_seq_npz_valid \
  --weights weights/alphagenome_pytorch/model_all_folds.safetensors \
  --output-dir runs/rna_seq11_adapter_1bp_smoke_20260515_1step \
  --batch-size 1 \
  --max-steps 1 \
  --eval-every 1 \
  --max-train-examples 1 \
  --max-valid-examples 1 \
  --learning-rate 3e-4 \
  --embedding-resolution 1 \
  --seed 20260515 \
  --device auto \
  > runs/rna_seq11_adapter_1bp_smoke_20260515_1step/train.log 2>&1
```

- Input data: First example from `alphagenome_custom/datasets/rna_seq_npz_train`; first example from `alphagenome_custom/datasets/rna_seq_npz_valid`; base weights `weights/alphagenome_pytorch/model_all_folds.safetensors`. The held-out test split was not used.
- Output path: `runs/rna_seq11_adapter_1bp_smoke_20260515_1step/` containing `config.json` (`914` bytes), `metrics.tsv` (`152` bytes), `train.log` (`1.4K`), `adapter_head.pt` (`69K`), and `adapter_head_best.pt` (`69K`). This path is ignored by Git.
- Result summary: Completed successfully on a single allowed GPU (`CUDA_VISIBLE_DEVICES=2`). The 1 bp embedding adapter produced full-length `1x11x1048576` predictions, completed one backward/optimizer step, and evaluated one validation example.
- Verification: Observed `embedding_resolution=1`, `base_parameters=450452613`, `trainable_parameters=16907`, `head_parameters=16907`, `base_has_grad=False`, `dna_sequence_shape=1x4x1048576`, `rna_seq_shape=1x11x1048576`, `prediction_shape=1x11x1048576`, one-step train loss `2.734642`, one-example valid loss `1.8309213`, and peak CUDA memory allocation about `35994.9` MB. Post-run check showed no active GPU compute process.
- Failures or warnings: This is only a feasibility smoke test on one train example and one validation example, not a model result. Exposing GPUs 2 and 3 would not automatically split the model with the current script; this run showed a single A100 80GB GPU is sufficient for the minimal 1 bp adapter path.
- Next actions: If continuing toward a closer-to-original 1 bp-resolution adapter experiment, run a validation-only pilot schedule first, for example 20 to 100 steps with full validation, and keep the held-out test untouched for this new experiment round.
- Claim status: verified

## 2026-05-15 - Frozen AlphaGenome Adapter Best Checkpoint Held-Out Test Evaluation on HY-GPU GPU 2

- Run type: evaluation
- Purpose: Perform the first held-out test evaluation using the checkpoint selected by the fixed validation-set rule after the 5000-step frozen AlphaGenome adapter-only run.
- Git commit: `f36a2c816fbe780918a1e67cdd2bcc1c1646a85e`; working tree contained uncommitted adapter scripts, baseline scripts, and experiment-log updates from the current development session.
- Branch: `setup/agent-maintenance`
- Host: `HY-GPU`
- Slurm job ID: Not applicable; HY-GPU is a non-Slurm server
- Slurm request: Not applicable
- Environment: Conda environment `alphagenome`; `CUDA_VISIBLE_DEVICES=2`; evaluation ran on an NVIDIA A100 80GB PCIe GPU
- Command:

```bash
cd /home/zelinli6/Alphagenome

CUDA_VISIBLE_DEVICES=2 conda run -n alphagenome python -u scripts/alphagenome_rna_seq11_eval.py \
  --dataset-dir alphagenome_custom/datasets/rna_seq_npz_test \
  --weights weights/alphagenome_pytorch/model_all_folds.safetensors \
  --checkpoint runs/rna_seq11_adapter_formal_20260515_5000steps_lr3e-4_128bp/adapter_head_best.pt \
  --batch-size 1 \
  --device auto \
  --per-track \
  --point-metrics \
  --spearman-sample-size 200000 \
  --spearman-seed 20260515 \
  --metrics-output runs/rna_seq11_adapter_formal_20260515_5000steps_lr3e-4_128bp/test_best_pointwise_metrics.tsv \
  > runs/rna_seq11_adapter_formal_20260515_5000steps_lr3e-4_128bp/test_best_eval.log 2>&1
```

- Input data: `alphagenome_custom/datasets/rna_seq_npz_test` (`33` examples); base weights `weights/alphagenome_pytorch/model_all_folds.safetensors`; selected checkpoint `runs/rna_seq11_adapter_formal_20260515_5000steps_lr3e-4_128bp/adapter_head_best.pt`.
- Output path: `runs/rna_seq11_adapter_formal_20260515_5000steps_lr3e-4_128bp/test_best_eval.log` (`2.2K`) and `runs/rna_seq11_adapter_formal_20260515_5000steps_lr3e-4_128bp/test_best_pointwise_metrics.tsv` (`1.5K`). These paths are ignored by Git.
- Result summary: Completed successfully. This is the first held-out test result for the current frozen-trunk adapter-only pipeline after fixing the validation-based selection rule.
- Verification: Overall held-out test metrics were `mse=1.1711332`, `mae=0.76486897`, and `pearson=0.55539986` over `33` test examples. Per-track MSE ranged from `1.0176844` (`RNA_SEQ_003`) to `1.3059002` (`RNA_SEQ_005`). Sampled per-track Spearman used `200000` positions per track with seed `20260515` and ranged from `0.39898333` (`RNA_SEQ_001`) to `0.52200972` (`RNA_SEQ_009`). Peak CUDA memory allocation was `17842.8` MB.
- Failures or warnings: This test split has now been used once for the selected model and should not be used for further model selection or tuning. Spearman is sampled, not exact over all test base positions.
- Next actions: Summarize validation/test/baseline comparison and freeze this result as the current benchmark before starting any next-round model changes.
- Claim status: verified

## 2026-05-15 - Frozen AlphaGenome RNA-seq 11-Track Adapter 5000-Step Validation Run on HY-GPU GPU 2

- Run type: training
- Purpose: Run the first longer frozen-trunk AlphaGenome adapter-only fine-tuning experiment for C. elegans 11-track RNA-seq, after the stronger learned tiny Conv1d baseline.
- Git commit: `f36a2c816fbe780918a1e67cdd2bcc1c1646a85e`; working tree contained uncommitted adapter scripts, baseline scripts, and experiment-log updates from the current development session.
- Branch: `setup/agent-maintenance`
- Host: `HY-GPU`
- Slurm job ID: Not applicable; HY-GPU is a non-Slurm server
- Slurm request: Not applicable
- Environment: Conda environment `alphagenome`; PyTorch `2.11.0+cu128`; PyTorch CUDA `12.8`; driver CUDA `13.2`; `alphagenome_pytorch 0.3.1`; `CUDA_VISIBLE_DEVICES=2`
- Command:

```bash
cd /home/zelinli6/Alphagenome

mkdir -p runs/rna_seq11_adapter_formal_20260515_5000steps_lr3e-4_128bp

CUDA_VISIBLE_DEVICES=2 conda run -n alphagenome \
  python -u scripts/alphagenome_rna_seq11_finetune.py \
  --train-dataset-dir alphagenome_custom/datasets/rna_seq_npz_train \
  --valid-dataset-dir alphagenome_custom/datasets/rna_seq_npz_valid \
  --weights weights/alphagenome_pytorch/model_all_folds.safetensors \
  --output-dir runs/rna_seq11_adapter_formal_20260515_5000steps_lr3e-4_128bp \
  --batch-size 1 \
  --max-steps 5000 \
  --eval-every 250 \
  --learning-rate 3e-4 \
  --embedding-resolution 128 \
  --seed 20260515 \
  --device auto \
  > runs/rna_seq11_adapter_formal_20260515_5000steps_lr3e-4_128bp/train.log 2>&1
```

- Input data: `alphagenome_custom/datasets/rna_seq_npz_train` (`116` examples); `alphagenome_custom/datasets/rna_seq_npz_valid` (`39` examples); base weights `weights/alphagenome_pytorch/model_all_folds.safetensors`. The held-out test split was not used.
- Output path: `runs/rna_seq11_adapter_formal_20260515_5000steps_lr3e-4_128bp/` containing `config.json` (`944` bytes), `metrics.tsv` (`201K`), `train.log` (`1003K`), final checkpoint `adapter_head.pt` (`135K`), and best-validation checkpoint `adapter_head_best.pt` (`135K`). This path is ignored by Git.
- Result summary: Completed successfully. The AlphaGenome base model stayed frozen (`base_has_grad=False` in training logs), and the 11-track adapter/head trained `33803` parameters while reusing `450452613` frozen base parameters. The run evaluated the full validation split every 250 steps.
- Verification: Validation MSEs included `1.2141448` at step 500, `1.1391212` at step 1000, `1.1319546` at step 1500, `1.1181387` at step 2000, `1.0989698` at step 2500, `1.0935162` at step 3000, `1.0837174` at step 3250, `1.0615692` at step 4000, `1.0609935` at step 4250, `1.0840007` at step 4500, `1.0764805` at step 4750, and final-step `1.1053824` at step 5000. Best validation MSE was `1.0609935` at step 4250, saved to `adapter_head_best.pt`. Final-step checkpoint was saved to `adapter_head.pt`.
- Failures or warnings: This is validation-set development, not held-out test performance. The selection rule is best full-validation MSE; the final step was worse than the best checkpoint. The base model was not unfrozen, so this run only tests whether frozen AlphaGenome representations plus a small adapter/head help on the custom C. elegans RNA-seq tracks.
- Next actions: Use `adapter_head_best.pt` for reload sanity evaluation and comparison against the learned tiny Conv1d baseline before any held-out test evaluation.
- Claim status: verified

## 2026-05-15 - Frozen AlphaGenome Adapter 5000-Step Best Checkpoint Reload Evaluation on HY-GPU GPU 2

- Run type: evaluation
- Purpose: Reload the best checkpoint from the 5000-step frozen AlphaGenome adapter-only run and compute full-validation overall and per-track metrics.
- Git commit: `f36a2c816fbe780918a1e67cdd2bcc1c1646a85e`; working tree contained uncommitted adapter scripts, baseline scripts, and experiment-log updates from the current development session.
- Branch: `setup/agent-maintenance`
- Host: `HY-GPU`
- Slurm job ID: Not applicable; HY-GPU is a non-Slurm server
- Slurm request: Not applicable
- Environment: Conda environment `alphagenome`; PyTorch `2.11.0+cu128`; PyTorch CUDA `12.8`; driver CUDA `13.2`; `alphagenome_pytorch 0.3.1`; `CUDA_VISIBLE_DEVICES=2`
- Command:

```bash
cd /home/zelinli6/Alphagenome

CUDA_VISIBLE_DEVICES=2 conda run -n alphagenome \
  python -u scripts/alphagenome_rna_seq11_eval.py \
  --dataset-dir alphagenome_custom/datasets/rna_seq_npz_valid \
  --weights weights/alphagenome_pytorch/model_all_folds.safetensors \
  --checkpoint runs/rna_seq11_adapter_formal_20260515_5000steps_lr3e-4_128bp/adapter_head_best.pt \
  --batch-size 1 \
  --device auto \
  --per-track \
  --point-metrics \
  --spearman-sample-size 200000 \
  --spearman-seed 20260515 \
  --metrics-output runs/rna_seq11_adapter_formal_20260515_5000steps_lr3e-4_128bp/valid_best_pointwise_metrics.tsv \
  > runs/rna_seq11_adapter_formal_20260515_5000steps_lr3e-4_128bp/reload_best_eval.log 2>&1
```

- Input data: `alphagenome_custom/datasets/rna_seq_npz_valid` (`39` examples); base weights `weights/alphagenome_pytorch/model_all_folds.safetensors`; checkpoint `runs/rna_seq11_adapter_formal_20260515_5000steps_lr3e-4_128bp/adapter_head_best.pt`. The held-out test split was not used.
- Output path: `runs/rna_seq11_adapter_formal_20260515_5000steps_lr3e-4_128bp/reload_best_eval.log` and `runs/rna_seq11_adapter_formal_20260515_5000steps_lr3e-4_128bp/valid_best_pointwise_metrics.tsv`. These paths are ignored by Git.
- Result summary: Completed successfully. Reloaded checkpoint reproduced the best full-validation MSE from training.
- Verification: Overall full-validation metrics were `mse=1.0609935`, `mae=0.71866247`, and `pearson=0.58809888`. Sampled per-track Spearman used `200000` positions per track with seed `20260515` and ranged from `0.4215102` (`RNA_SEQ_001`) to `0.54338671` (`RNA_SEQ_010`). Per-track MSE ranged from `0.92602274` (`RNA_SEQ_007`) to `1.2535293` (`RNA_SEQ_002`). This best adapter checkpoint was stronger than the stronger tiny Conv1d baseline best validation MSE (`1.0609935` vs `1.553865`) and stronger than its final pointwise Pearson (`0.58809888` vs `0.27318148`).
- Failures or warnings: Spearman is sampled, not exact over all validation base positions. This remains validation-set development, not held-out test performance.
- Next actions: Keep test untouched until baseline, metrics, and selection rule are finalized. Candidate selection rule after this run is best full-validation MSE, currently `adapter_head_best.pt` from step 4250.
- Claim status: verified

## 2026-05-15 - Stronger Tiny Conv1d Learned Baseline 2000-Step Validation Run on HY-GPU GPU 2

- Run type: training
- Purpose: Train a stronger learned sequence-only baseline from scratch on the same C. elegans 11-track RNA-seq train/valid NPZ split before the longer adapter-only run.
- Git commit: `f36a2c816fbe780918a1e67cdd2bcc1c1646a85e`; working tree contained uncommitted adapter scripts, baseline scripts, and experiment-log updates from the current development session.
- Branch: `setup/agent-maintenance`
- Host: `HY-GPU`
- Slurm job ID: Not applicable; HY-GPU is a non-Slurm server
- Slurm request: Not applicable
- Environment: Conda environment `alphagenome`; PyTorch `2.11.0+cu128`; PyTorch CUDA `12.8`; driver CUDA `13.2`; `CUDA_VISIBLE_DEVICES=2`
- Command:

```bash
cd /home/zelinli6/Alphagenome

mkdir -p runs/rna_seq11_tiny_conv_baseline_20260515_2000steps_lr1e-3_h64

CUDA_VISIBLE_DEVICES=2 conda run -n alphagenome \
  python -u scripts/rna_seq11_tiny_conv_baseline.py \
  --train-dataset-dir alphagenome_custom/datasets/rna_seq_npz_train \
  --valid-dataset-dir alphagenome_custom/datasets/rna_seq_npz_valid \
  --output-dir runs/rna_seq11_tiny_conv_baseline_20260515_2000steps_lr1e-3_h64 \
  --batch-size 1 \
  --max-steps 2000 \
  --eval-every 200 \
  --hidden-channels 64 \
  --learning-rate 1e-3 \
  --seed 20260515 \
  --spearman-sample-size 200000 \
  --spearman-seed 20260515 \
  --device auto \
  > runs/rna_seq11_tiny_conv_baseline_20260515_2000steps_lr1e-3_h64/train.log 2>&1
```

- Input data: `alphagenome_custom/datasets/rna_seq_npz_train` (`116` examples); `alphagenome_custom/datasets/rna_seq_npz_valid` (`39` examples). The held-out test split was not used.
- Output path: `runs/rna_seq11_tiny_conv_baseline_20260515_2000steps_lr1e-3_h64/` containing `config.json` (`941` bytes), `metrics.tsv` (`76K`), `train.log` (`61K`), checkpoint `tiny_conv.pt` (`262K`), `valid_pointwise_metrics.tsv` (`1.5K`), and comparison file `valid_adapter500_vs_tiny_h64_final.tsv`. This path is ignored by Git.
- Result summary: Completed successfully. The stronger tiny Conv1d baseline trained all `66123` model parameters from scratch and evaluated the full validation split every 200 steps.
- Verification: Observed validation MSEs `1.5847078` at step 200, `1.553865` at step 400, `1.6888396` at step 600, `1.8387403` at step 800, `1.6478015` at step 1000, `1.6544384` at step 1200, `1.6718817` at step 1400, `1.6301199` at step 1600, `1.7003493` at step 1800, and `1.6456854` at step 2000. Final pointwise metrics were `mse=1.6456854`, `mae=1.0224008`, and `pearson=0.27318148`; sampled per-track Spearman ranged from `0.1631497` (`RNA_SEQ_002`) to `0.23621425` (`RNA_SEQ_009`). Best validation MSE was `1.553865` at step 400. Compared with the 500-step frozen-AlphaGenome adapter, the adapter still had lower final overall validation MSE (`1.2418628` vs `1.6456854`) and higher Pearson (`0.52766246` vs `0.27318148`).
- Failures or warnings: This is a validation-set development baseline, not held-out test performance. The stronger tiny Conv1d baseline is still a simple local Conv1d architecture and has not been fully hyperparameter tuned. Spearman is sampled, not exact over all validation base positions.
- Next actions: Run the 5000-step frozen-AlphaGenome adapter-only training with best-valid checkpoint saving enabled, then compare against this baseline before touching the held-out test split.
- Claim status: verified

## 2026-05-15 - Tiny Conv1d Learned Baseline 500-Step Validation Run on HY-GPU GPU 2

- Run type: training
- Purpose: Train a learned sequence-only baseline from scratch on the same C. elegans 11-track RNA-seq train/valid NPZ split, so the 500-step frozen-AlphaGenome adapter can be compared against more than a train-track mean predictor.
- Git commit: `f36a2c816fbe780918a1e67cdd2bcc1c1646a85e`; working tree contained uncommitted adapter scripts, baseline scripts, and experiment-log updates from the current development session.
- Branch: `setup/agent-maintenance`
- Host: `HY-GPU`
- Slurm job ID: Not applicable; HY-GPU is a non-Slurm server
- Slurm request: Not applicable
- Environment: Conda environment `alphagenome`; PyTorch `2.11.0+cu128`; PyTorch CUDA `12.8`; driver CUDA `13.2`; `CUDA_VISIBLE_DEVICES=2`
- Command:

```bash
cd /home/zelinli6/Alphagenome

mkdir -p runs/rna_seq11_tiny_conv_baseline_20260515_500steps_lr3e-4_h16

CUDA_VISIBLE_DEVICES=2 conda run -n alphagenome \
  python scripts/rna_seq11_tiny_conv_baseline.py \
  --train-dataset-dir alphagenome_custom/datasets/rna_seq_npz_train \
  --valid-dataset-dir alphagenome_custom/datasets/rna_seq_npz_valid \
  --output-dir runs/rna_seq11_tiny_conv_baseline_20260515_500steps_lr3e-4_h16 \
  --batch-size 1 \
  --max-steps 500 \
  --eval-every 100 \
  --hidden-channels 16 \
  --learning-rate 3e-4 \
  --seed 20260515 \
  --spearman-sample-size 200000 \
  --spearman-seed 20260515 \
  --device auto \
  2>&1 | tee runs/rna_seq11_tiny_conv_baseline_20260515_500steps_lr3e-4_h16/train.log
```

- Input data: `alphagenome_custom/datasets/rna_seq_npz_train` (`116` examples); `alphagenome_custom/datasets/rna_seq_npz_valid` (`39` examples). The held-out test split was not used.
- Output path: `runs/rna_seq11_tiny_conv_baseline_20260515_500steps_lr3e-4_h16/` containing `config.json` (`940` bytes), `metrics.tsv` (`19K`), `train.log` (`16K`), checkpoint `tiny_conv.pt` (`23K`), `valid_pointwise_metrics.tsv` (`1.5K`), and comparison file `valid_adapter_vs_tiny_conv.tsv` (`2.2K`). This path is ignored by Git.
- Result summary: Completed successfully. The tiny Conv1d baseline trained all `5019` model parameters from scratch and evaluated the full validation split every 100 steps.
- Verification: Observed validation MSEs `1.7203917` at step 100, `1.7174804` at step 200, `1.719583` at step 300, `1.6157231` at step 400, and `1.6399873` at step 500. Final pointwise metrics were `mse=1.6399873`, `mae=1.0042184`, and `pearson=0.19711665`; sampled per-track Spearman ranged from `0.12864133` (`RNA_SEQ_002`) to `0.19100426` (`RNA_SEQ_009`). Compared with the 500-step frozen-AlphaGenome adapter, the adapter had lower overall validation MSE (`1.2418628` vs `1.6399873`), lower MAE (`0.84385942` vs `1.0042184`), and higher Pearson (`0.52766246` vs `0.19711665`). The adapter's relative MSE improvement over this tiny Conv1d baseline was `24.3%` overall, with per-track relative improvements from `20.0%` to `28.9%`. The tiny Conv1d baseline was still better than the train-track mean baseline MSE (`1.7282558`), by about `5.1%` relative MSE.
- Failures or warnings: This is a validation-set development baseline, not held-out test performance. The tiny Conv1d baseline is intentionally small and matched to the adapter pilot schedule; it has not been hyperparameter tuned. Spearman is sampled, not exact over all validation base positions.
- Next actions: Decide whether to tune the learned baseline further, for example `lr=1e-3` or more steps, before freezing the selection rule and using the held-out test split.
- Claim status: verified

## 2026-05-15 - Adapter 500-Step Validation Pointwise Metrics on HY-GPU GPU 2

- Run type: evaluation
- Purpose: Add pointwise validation metrics beyond MSE for the 500-step adapter-only checkpoint, including exact MAE and Pearson correlation plus deterministic sampled Spearman correlation per track.
- Git commit: `f36a2c816fbe780918a1e67cdd2bcc1c1646a85e`; working tree contained uncommitted adapter scripts, baseline script, and experiment-log updates from the current development session.
- Branch: `setup/agent-maintenance`
- Host: `HY-GPU`
- Slurm job ID: Not applicable; HY-GPU is a non-Slurm server
- Slurm request: Not applicable
- Environment: Conda environment `alphagenome`; `alphagenome-pytorch` `0.3.1`; PyTorch `2.11.0+cu128`; PyTorch CUDA `12.8`; driver CUDA `13.2`; `CUDA_VISIBLE_DEVICES=2`
- Command:

```bash
conda activate alphagenome
cd /home/zelinli6/Alphagenome

CUDA_VISIBLE_DEVICES=2 python scripts/alphagenome_rna_seq11_eval.py \
  --dataset-dir alphagenome_custom/datasets/rna_seq_npz_valid \
  --weights weights/alphagenome_pytorch/model_all_folds.safetensors \
  --checkpoint runs/rna_seq11_adapter_formal_20260515_500steps_lr3e-4/adapter_head.pt \
  --batch-size 1 \
  --device cuda \
  --point-metrics \
  --spearman-sample-size 200000 \
  --spearman-seed 20260515 \
  --metrics-output runs/rna_seq11_adapter_formal_20260515_500steps_lr3e-4/valid_pointwise_metrics.tsv
```

- Input data: `alphagenome_custom/datasets/rna_seq_npz_valid` (`39` examples); `weights/alphagenome_pytorch/model_all_folds.safetensors`; `runs/rna_seq11_adapter_formal_20260515_500steps_lr3e-4/adapter_head.pt`
- Output path: `runs/rna_seq11_adapter_formal_20260515_500steps_lr3e-4/valid_pointwise_metrics.tsv` (`1457` bytes). This path is ignored by Git.
- Result summary: Completed successfully. Exact pointwise MSE, MAE, and Pearson were computed on the full validation split. Spearman was computed on deterministic sampled positions, with `200463` sampled positions per track.
- Verification: Observed overall `mse=1.2418628`, `mae=0.84385942`, and `pearson=0.52766246`. Per-track Pearson ranged from `0.4820762` (`RNA_SEQ_001`) to `0.57240715` (`RNA_SEQ_010`). Sampled per-track Spearman ranged from `0.38375246` (`RNA_SEQ_001`) to `0.48221006` (`RNA_SEQ_009`). Peak `cuda_max_memory_allocated_mb=17842.8`. Post-run check showed no active GPU compute process.
- Failures or warnings: Spearman is sampled, not exact over all valid base positions. These are validation-set development metrics, not held-out test performance or biological conclusions.
- Next actions: Use these validation metrics to decide whether to run a learned tiny Conv1d baseline before touching the held-out test split.
- Claim status: verified

## 2026-05-15 - Train-Track Mean Baseline on Validation Split

- Run type: evaluation
- Purpose: Compute a simple baseline for interpreting the 500-step adapter's validation MSE. The baseline estimates one global `log1p` mean per RNA-seq track from the full train split and predicts that constant value at every validation position.
- Git commit: `f36a2c816fbe780918a1e67cdd2bcc1c1646a85e`; working tree contained uncommitted adapter scripts, baseline script, and experiment-log updates from the current development session.
- Branch: `setup/agent-maintenance`
- Host: `HY-GPU`
- Slurm job ID: Not applicable; HY-GPU is a non-Slurm server
- Slurm request: Not applicable
- Environment: Conda environment `alphagenome`; CPU/IO evaluation; no GPU work was used for the baseline.
- Command:

```bash
conda activate alphagenome
cd /home/zelinli6/Alphagenome

conda run -n alphagenome python scripts/rna_seq11_mean_baseline.py \
  --train-dataset-dir alphagenome_custom/datasets/rna_seq_npz_train \
  --eval-dataset-dir alphagenome_custom/datasets/rna_seq_npz_valid \
  --target-transform log1p \
  --metrics-output runs/rna_seq11_adapter_formal_20260515_500steps_lr3e-4/valid_track_mean_baseline_metrics.tsv
```

- Input data: `alphagenome_custom/datasets/rna_seq_npz_train` (`116` examples); `alphagenome_custom/datasets/rna_seq_npz_valid` (`39` examples)
- Output path: `runs/rna_seq11_adapter_formal_20260515_500steps_lr3e-4/valid_track_mean_baseline_metrics.tsv` (`1102` bytes), plus comparison file `runs/rna_seq11_adapter_formal_20260515_500steps_lr3e-4/valid_adapter_vs_mean_baseline.tsv` (`1293` bytes). These paths are ignored by Git.
- Result summary: Completed successfully. The train-track mean baseline's overall full-validation pointwise masked MSE was worse than the 500-step adapter checkpoint.
- Verification: Observed train-track mean baseline overall valid MSE `1.7282558`; 500-step adapter overall valid MSE `1.2418628`; absolute improvement `0.486393`; relative improvement `0.28143577` (`28.1%`). Per-track relative improvements ranged from `23.2%` (`RNA_SEQ_001`) to `33.7%` (`RNA_SEQ_009`).
- Failures or warnings: This is a simple pointwise baseline using `log1p` targets. It does not test higher-level biological validity, peak/profile correlation, or held-out test performance.
- Next actions: Add a stronger learned baseline, such as the existing tiny Conv1d model with the same train/valid split and per-track metrics, before using the held-out test split.
- Claim status: verified

## 2026-05-15 - Adapter 500-Step Full-Validation Per-Track Metrics on HY-GPU GPU 2

- Run type: evaluation
- Purpose: Evaluate the 500-step adapter-only checkpoint on the full validation split and report both overall and per-track masked MSE for the 11 C. elegans RNA-seq tracks.
- Git commit: `f36a2c816fbe780918a1e67cdd2bcc1c1646a85e`; working tree contained uncommitted adapter scripts and experiment-log updates from the current development session.
- Branch: `setup/agent-maintenance`
- Host: `HY-GPU`
- Slurm job ID: Not applicable; HY-GPU is a non-Slurm server
- Slurm request: Not applicable
- Environment: Conda environment `alphagenome`; `alphagenome-pytorch` `0.3.1`; PyTorch `2.11.0+cu128`; PyTorch CUDA `12.8`; driver CUDA `13.2`; `CUDA_VISIBLE_DEVICES=2`
- Command:

```bash
conda activate alphagenome
cd /home/zelinli6/Alphagenome

hostname
pwd
git rev-parse HEAD
git status --short --branch
nvidia-smi

CUDA_VISIBLE_DEVICES=2 python scripts/alphagenome_rna_seq11_eval.py \
  --dataset-dir alphagenome_custom/datasets/rna_seq_npz_valid \
  --weights weights/alphagenome_pytorch/model_all_folds.safetensors \
  --checkpoint runs/rna_seq11_adapter_formal_20260515_500steps_lr3e-4/adapter_head.pt \
  --batch-size 1 \
  --device cuda \
  --per-track \
  --metrics-output runs/rna_seq11_adapter_formal_20260515_500steps_lr3e-4/valid_per_track_metrics.tsv
```

- Input data: `alphagenome_custom/datasets/rna_seq_npz_valid` (`39` examples); `weights/alphagenome_pytorch/model_all_folds.safetensors`; `runs/rna_seq11_adapter_formal_20260515_500steps_lr3e-4/adapter_head.pt`
- Output path: `runs/rna_seq11_adapter_formal_20260515_500steps_lr3e-4/valid_per_track_metrics.tsv` (`953` bytes). This path is ignored by Git.
- Result summary: Completed successfully. The overall full-validation masked MSE matched the 500-step checkpoint reload value, and per-track masked MSE was reported for all 11 RNA-seq tracks.
- Verification: Observed `valid_batches=39`, `valid_examples=39`, `valid_loss=1.2418628`, `cuda_max_memory_allocated_mb=17604.8`, and per-track losses: `RNA_SEQ_001=1.27215`, `RNA_SEQ_002=1.4491605`, `RNA_SEQ_003=1.1401407`, `RNA_SEQ_004=1.3657115`, `RNA_SEQ_005=1.3945642`, `RNA_SEQ_006=1.3448502`, `RNA_SEQ_007=1.1034249`, `RNA_SEQ_008=1.1089619`, `RNA_SEQ_009=1.1232189`, `RNA_SEQ_010=1.2121915`, and `RNA_SEQ_011=1.1461164`. Post-run check showed no active GPU compute process.
- Failures or warnings: This is validation-set evaluation for model development, not held-out test performance. Do not treat these values as final model performance or biological conclusions.
- Next actions: Add a simple baseline comparison on the same valid split before touching the held-out test split.
- Claim status: verified

## 2026-05-15 - Adapter 100-Step Checkpoint Reload Full-Validation Sanity Check on HY-GPU GPU 2

- Run type: sanity check
- Purpose: Reload the 100-step adapter-only checkpoint and verify that full-validation masked MSE is reproducible before starting a longer adapter-only training run.
- Git commit: `f36a2c816fbe780918a1e67cdd2bcc1c1646a85e`; working tree contained uncommitted adapter scripts and experiment-log updates from the current development session.
- Branch: `setup/agent-maintenance`
- Host: `HY-GPU`
- Slurm job ID: Not applicable; HY-GPU is a non-Slurm server
- Slurm request: Not applicable
- Environment: Conda environment `alphagenome`; `alphagenome-pytorch` `0.3.1`; PyTorch `2.11.0+cu128`; PyTorch CUDA `12.8`; driver CUDA `13.2`; `CUDA_VISIBLE_DEVICES=2`
- Command:

```bash
conda activate alphagenome
cd /home/zelinli6/Alphagenome

hostname
pwd
git rev-parse HEAD
git status --short --branch
nvidia-smi

CUDA_VISIBLE_DEVICES=2 python scripts/alphagenome_rna_seq11_eval.py \
  --dataset-dir alphagenome_custom/datasets/rna_seq_npz_valid \
  --weights weights/alphagenome_pytorch/model_all_folds.safetensors \
  --checkpoint runs/rna_seq11_adapter_pilot_20260514_100steps_lr3e-4/adapter_head.pt \
  --batch-size 1 \
  --device cuda
```

- Input data: `alphagenome_custom/datasets/rna_seq_npz_valid` (`39` examples); `weights/alphagenome_pytorch/model_all_folds.safetensors`; `runs/rna_seq11_adapter_pilot_20260514_100steps_lr3e-4/adapter_head.pt`
- Output path: Interactive terminal output only; no new checkpoint, metrics file, or prediction file was written.
- Result summary: Completed successfully. The 100-step adapter checkpoint reloaded and reproduced its recorded full-validation loss.
- Verification: Observed `valid_batches=39`, `valid_examples=39`, `valid_loss=1.4401001`, `cuda_max_memory_allocated_mb=17560.9`, `base_parameters=450452613`, `trainable_parameters=33803`, and `head_parameters=33803`.
- Failures or warnings: This is a checkpoint reload sanity check, not a biological or model-performance claim. The adapter still uses `organism_index=0` as a frozen AlphaGenome trunk selector.
- Next actions: Start the 500-step adapter-only training run.
- Claim status: verified

## 2026-05-15 - Frozen AlphaGenome 11-Track RNA-seq Adapter 500-Step Formal Pilot on HY-GPU GPU 2

- Run type: training
- Purpose: Run the first longer adapter-only training pilot for the C. elegans 11-track RNA-seq head, keeping the AlphaGenome trunk frozen and using 128 bp embeddings.
- Git commit: `f36a2c816fbe780918a1e67cdd2bcc1c1646a85e`; working tree contained uncommitted adapter scripts and experiment-log updates from the current development session.
- Branch: `setup/agent-maintenance`
- Host: `HY-GPU`
- Slurm job ID: Not applicable; HY-GPU is a non-Slurm server
- Slurm request: Not applicable
- Environment: Conda environment `alphagenome`; `alphagenome-pytorch` `0.3.1`; PyTorch `2.11.0+cu128`; PyTorch CUDA `12.8`; driver CUDA `13.2`; `CUDA_VISIBLE_DEVICES=2`
- Command:

```bash
conda activate alphagenome
cd /home/zelinli6/Alphagenome

hostname
pwd
git rev-parse HEAD
git status --short --branch
nvidia-smi

set -o pipefail
mkdir -p runs/rna_seq11_adapter_formal_20260515_500steps_lr3e-4
CUDA_VISIBLE_DEVICES=2 python scripts/alphagenome_rna_seq11_finetune.py \
  --train-dataset-dir alphagenome_custom/datasets/rna_seq_npz_train \
  --valid-dataset-dir alphagenome_custom/datasets/rna_seq_npz_valid \
  --weights weights/alphagenome_pytorch/model_all_folds.safetensors \
  --batch-size 1 \
  --max-steps 500 \
  --eval-every 100 \
  --learning-rate 3e-4 \
  --embedding-resolution 128 \
  --device cuda \
  --output-dir runs/rna_seq11_adapter_formal_20260515_500steps_lr3e-4 \
  2>&1 | tee runs/rna_seq11_adapter_formal_20260515_500steps_lr3e-4/train.log
```

- Input data: `alphagenome_custom/datasets/rna_seq_npz_train` (`116` examples); `alphagenome_custom/datasets/rna_seq_npz_valid` (`39` examples); `weights/alphagenome_pytorch/model_all_folds.safetensors`
- Output path: `runs/rna_seq11_adapter_formal_20260515_500steps_lr3e-4/` containing `config.json` (`936` bytes), `metrics.tsv` (`20132` bytes), `train.log` (`102877` bytes), and adapter-only checkpoint `adapter_head.pt` (`137357` bytes). This path is ignored by Git.
- Result summary: Completed successfully. The run trained only the 11-track adapter head for 500 optimizer steps and evaluated the full validation split every 100 steps. The full AlphaGenome trunk remained frozen.
- Verification: Observed `n_train_examples=116`, `n_valid_examples=39`, `n_tracks=11`, `base_parameters=450452613`, `trainable_parameters=33803`, `head_parameters=33803`, repeated `prediction_shape=1x11x1048576`, `base_has_grad=False` at every train step, full-validation losses `1.4401001` at step 100, `1.3572929` at step 200, `1.2502746` at step 300, `1.261003` at step 400, and `1.2418628` at step 500, and peak `cuda_max_memory_allocated_mb=17665.2`.
- Failures or warnings: This remains a provisional adapter-only training run, not a biological or paper-level model-performance result. The validation loss improved through step 300, rose slightly at step 400, and ended slightly lower at step 500. The adapter still uses `organism_index=0` as a frozen AlphaGenome trunk selector, not as a verified C. elegans organism embedding.
- Next actions: Reload the 500-step checkpoint to verify checkpoint reproducibility, then decide whether to run a held-out test evaluation once model selection criteria are fixed.
- Claim status: verified

## 2026-05-15 - Adapter 500-Step Checkpoint Reload Full-Validation Sanity Check on HY-GPU GPU 2

- Run type: sanity check
- Purpose: Reload the 500-step adapter-only checkpoint and verify that full-validation masked MSE matches the training run's final validation result.
- Git commit: `f36a2c816fbe780918a1e67cdd2bcc1c1646a85e`; working tree contained uncommitted adapter scripts and experiment-log updates from the current development session.
- Branch: `setup/agent-maintenance`
- Host: `HY-GPU`
- Slurm job ID: Not applicable; HY-GPU is a non-Slurm server
- Slurm request: Not applicable
- Environment: Conda environment `alphagenome`; `alphagenome-pytorch` `0.3.1`; PyTorch `2.11.0+cu128`; PyTorch CUDA `12.8`; driver CUDA `13.2`; `CUDA_VISIBLE_DEVICES=2`
- Command:

```bash
conda activate alphagenome
cd /home/zelinli6/Alphagenome

CUDA_VISIBLE_DEVICES=2 python scripts/alphagenome_rna_seq11_eval.py \
  --dataset-dir alphagenome_custom/datasets/rna_seq_npz_valid \
  --weights weights/alphagenome_pytorch/model_all_folds.safetensors \
  --checkpoint runs/rna_seq11_adapter_formal_20260515_500steps_lr3e-4/adapter_head.pt \
  --batch-size 1 \
  --device cuda
```

- Input data: `alphagenome_custom/datasets/rna_seq_npz_valid` (`39` examples); `weights/alphagenome_pytorch/model_all_folds.safetensors`; `runs/rna_seq11_adapter_formal_20260515_500steps_lr3e-4/adapter_head.pt`
- Output path: Interactive terminal output only; no new checkpoint, metrics file, or prediction file was written.
- Result summary: Completed successfully. The 500-step adapter checkpoint reloaded and reproduced the final full-validation loss from the training run.
- Verification: Observed `valid_batches=39`, `valid_examples=39`, `valid_loss=1.2418628`, `cuda_max_memory_allocated_mb=17560.9`, `base_parameters=450452613`, `trainable_parameters=33803`, and `head_parameters=33803`. This matches the training run's step 500 full-validation loss `1.2418628`.
- Failures or warnings: This is a checkpoint reload sanity check, not a biological or model-performance claim.
- Next actions: Fix model-selection criteria before using the held-out test split; consider a tiny 1 bp embedding feasibility check only after documenting the 128 bp adapter baseline.
- Claim status: verified

## 2026-05-14 - Frozen AlphaGenome 11-Track RNA-seq Adapter 100-Step Fine-Tune Pilot on HY-GPU GPU 2

- Run type: training
- Purpose: Run a more stable adapter-only pilot for the C. elegans 11-track RNA-seq head using a lower learning rate than the first 20-step pilot, while keeping the AlphaGenome trunk frozen and using 128 bp embeddings.
- Git commit: `f36a2c816fbe780918a1e67cdd2bcc1c1646a85e`; working tree contained uncommitted adapter scripts and experiment-log updates from the current development session.
- Branch: `setup/agent-maintenance`
- Host: `HY-GPU`
- Slurm job ID: Not applicable; HY-GPU is a non-Slurm server
- Slurm request: Not applicable
- Environment: Conda environment `alphagenome`; `alphagenome-pytorch` `0.3.1`; PyTorch `2.11.0+cu128`; PyTorch CUDA `12.8`; driver CUDA `13.2`; `CUDA_VISIBLE_DEVICES=2`
- Command:

```bash
conda activate alphagenome
cd /home/zelinli6/Alphagenome

hostname
pwd
git rev-parse HEAD
git status --short --branch
nvidia-smi

CUDA_VISIBLE_DEVICES=2 python scripts/alphagenome_rna_seq11_finetune.py \
  --train-dataset-dir alphagenome_custom/datasets/rna_seq_npz_train \
  --valid-dataset-dir alphagenome_custom/datasets/rna_seq_npz_valid \
  --weights weights/alphagenome_pytorch/model_all_folds.safetensors \
  --batch-size 1 \
  --max-steps 100 \
  --eval-every 20 \
  --learning-rate 3e-4 \
  --embedding-resolution 128 \
  --device cuda \
  --output-dir runs/rna_seq11_adapter_pilot_20260514_100steps_lr3e-4
```

- Input data: `alphagenome_custom/datasets/rna_seq_npz_train` (`116` examples); `alphagenome_custom/datasets/rna_seq_npz_valid` (`39` examples); `weights/alphagenome_pytorch/model_all_folds.safetensors`
- Output path: `runs/rna_seq11_adapter_pilot_20260514_100steps_lr3e-4/` containing `config.json` (`934` bytes), `metrics.tsv` (`4158` bytes), and adapter-only checkpoint `adapter_head.pt` (`137357` bytes). This path is ignored by Git.
- Result summary: Completed successfully. The run trained only the 11-track adapter head for 100 optimizer steps and evaluated the full validation split every 20 steps. The full AlphaGenome trunk remained frozen.
- Verification: Observed `n_train_examples=116`, `n_valid_examples=39`, `n_tracks=11`, `base_parameters=450452613`, `trainable_parameters=33803`, `head_parameters=33803`, repeated `prediction_shape=1x11x1048576`, `base_has_grad=False` at every train step, full-validation losses `1.8348786` at step 20, `1.6128415` at step 40, `1.5110018` at step 60, `1.4737058` at step 80, and `1.4401001` at step 100, and peak `cuda_max_memory_allocated_mb=17665.2`. Post-run check showed no active GPU compute process.
- Failures or warnings: This is a short pilot/provisional training run, not a model result or biological conclusion. The adapter still uses `organism_index=0` as a frozen AlphaGenome trunk selector, not as a verified C. elegans organism embedding.
- Next actions: Reload this 100-step adapter checkpoint with `scripts/alphagenome_rna_seq11_eval.py` to verify checkpoint reproducibility, then compare against a tiny 1 bp embedding feasibility check before deciding whether 128 bp resolution is sufficient.
- Claim status: verified

## 2026-05-14 - Adapter Checkpoint Reload Full-Validation Sanity Check on HY-GPU GPU 2

- Run type: sanity check
- Purpose: Verify that the adapter-only checkpoint from the 20-step pilot can be reloaded through the new eval entrypoint and reproduce the recorded full-validation masked MSE.
- Git commit: `f36a2c816fbe780918a1e67cdd2bcc1c1646a85e`
- Branch: `setup/agent-maintenance`
- Host: `HY-GPU`
- Slurm job ID: Not applicable; HY-GPU is a non-Slurm server
- Slurm request: Not applicable
- Environment: Conda environment `alphagenome`; `alphagenome-pytorch` `0.3.1`; PyTorch `2.11.0+cu128`; PyTorch CUDA `12.8`; driver CUDA `13.2`; `CUDA_VISIBLE_DEVICES=2`
- Command:

```bash
conda activate alphagenome
cd /home/zelinli6/Alphagenome

hostname
pwd
git rev-parse HEAD
git status --short --branch
nvidia-smi

CUDA_VISIBLE_DEVICES=2 python scripts/alphagenome_rna_seq11_eval.py \
  --dataset-dir alphagenome_custom/datasets/rna_seq_npz_valid \
  --weights weights/alphagenome_pytorch/model_all_folds.safetensors \
  --checkpoint runs/rna_seq11_adapter_pilot_20260514_20steps/adapter_head.pt \
  --batch-size 1 \
  --device cuda
```

- Input data: `alphagenome_custom/datasets/rna_seq_npz_valid` (`39` examples); `weights/alphagenome_pytorch/model_all_folds.safetensors`; `runs/rna_seq11_adapter_pilot_20260514_20steps/adapter_head.pt`
- Output path: Interactive terminal output only; no new checkpoint, metrics file, or prediction file was written.
- Result summary: Completed successfully. The adapter checkpoint reloaded and reproduced the 20-step pilot's full-validation loss.
- Verification: Observed `n_examples=39`, `n_tracks=11`, `embedding_resolution=128`, `organism_index=0`, `base_parameters=450452613`, `trainable_parameters=33803`, `head_parameters=33803`, `valid_batches=39`, `valid_examples=39`, `valid_loss=1.4536235`, and `cuda_max_memory_allocated_mb=17560.9`. This matches the prior pilot's step 20 full-validation loss `1.4536235`. Post-run check showed no active GPU compute process.
- Failures or warnings: This is a checkpoint reload sanity check, not a biological or model-performance claim. The adapter still uses `organism_index=0` as a frozen AlphaGenome trunk selector, not as a verified C. elegans organism embedding.
- Next actions: Run a longer adapter-only pilot with lower learning rate, or run a tiny 1 bp embedding feasibility check before deciding whether 128 bp resolution is too limiting.
- Claim status: verified

## 2026-05-14 - Frozen AlphaGenome 11-Track RNA-seq Adapter 20-Step Fine-Tune Pilot on HY-GPU GPU 2

- Run type: training
- Purpose: Run the first recorded short fine-tuning pilot for a C. elegans 11-track RNA-seq adapter head on a frozen `alphagenome-pytorch` trunk, with config, metrics, and adapter-only checkpoint output.
- Git commit: `f36a2c816fbe780918a1e67cdd2bcc1c1646a85e`
- Branch: `setup/agent-maintenance`
- Host: `HY-GPU`
- Slurm job ID: Not applicable; HY-GPU is a non-Slurm server
- Slurm request: Not applicable
- Environment: Conda environment `alphagenome`; `alphagenome-pytorch` `0.3.1`; PyTorch `2.11.0+cu128`; PyTorch CUDA `12.8`; driver CUDA `13.2`; `CUDA_VISIBLE_DEVICES=2`
- Command:

```bash
conda activate alphagenome
cd /home/zelinli6/Alphagenome

hostname
pwd
git rev-parse HEAD
git status --short --branch
nvidia-smi

CUDA_VISIBLE_DEVICES=2 python scripts/alphagenome_rna_seq11_finetune.py \
  --train-dataset-dir alphagenome_custom/datasets/rna_seq_npz_train \
  --valid-dataset-dir alphagenome_custom/datasets/rna_seq_npz_valid \
  --weights weights/alphagenome_pytorch/model_all_folds.safetensors \
  --batch-size 1 \
  --max-steps 20 \
  --eval-every 5 \
  --embedding-resolution 128 \
  --device cuda \
  --output-dir runs/rna_seq11_adapter_pilot_20260514_20steps
```

- Input data: `alphagenome_custom/datasets/rna_seq_npz_train` (`116` examples); `alphagenome_custom/datasets/rna_seq_npz_valid` (`39` examples); `weights/alphagenome_pytorch/model_all_folds.safetensors`
- Output path: `runs/rna_seq11_adapter_pilot_20260514_20steps/` containing `config.json` (`923` bytes), `metrics.tsv` (`982` bytes), and adapter-only checkpoint `adapter_head.pt` (`137357` bytes). This path is ignored by Git.
- Result summary: Completed successfully. The run trained only the 11-track adapter head for 20 optimizer steps and evaluated the full validation split every 5 steps. The full AlphaGenome trunk remained frozen.
- Verification: Observed `n_train_examples=116`, `n_valid_examples=39`, `n_tracks=11`, `base_parameters=450452613`, `trainable_parameters=33803`, `head_parameters=33803`, repeated `prediction_shape=1x11x1048576`, `base_has_grad=False` at every train step, full-validation losses `1.81721` at step 5, `2.431641` at step 10, `1.7478276` at step 15, and `1.4536235` at step 20, and peak `cuda_max_memory_allocated_mb=17665.2`. Post-run check showed no active GPU compute process.
- Failures or warnings: This is a short pilot/provisional training run, not a model result or biological conclusion. The adapter still uses `organism_index=0` as a frozen AlphaGenome trunk selector, not as a verified C. elegans organism embedding.
- Next actions: Inspect `metrics.tsv`, then decide whether to run a longer adapter-only pilot with a lower learning rate and less frequent validation, or test the 1 bp embedding path on a small step count for resolution comparison.
- Claim status: verified

## 2026-05-14 - Frozen AlphaGenome 11-Track RNA-seq Adapter Pilot Train-Valid Smoke Test on HY-GPU GPU 2

- Run type: smoke test
- Purpose: Verify a minimal train/validation loop for a frozen `alphagenome-pytorch` trunk with a trainable C. elegans 11-track RNA-seq adapter head using existing NPZ dataloaders.
- Git commit: `f36a2c816fbe780918a1e67cdd2bcc1c1646a85e`
- Branch: `setup/agent-maintenance`
- Host: `HY-GPU`
- Slurm job ID: Not applicable; HY-GPU is a non-Slurm server
- Slurm request: Not applicable
- Environment: Conda environment `alphagenome`; `alphagenome-pytorch` `0.3.1`; PyTorch `2.11.0+cu128`; PyTorch CUDA `12.8`; driver CUDA `13.2`; `CUDA_VISIBLE_DEVICES=2`
- Command:

```bash
conda activate alphagenome
cd /home/zelinli6/Alphagenome

hostname
pwd
git rev-parse HEAD
git status --short --branch
nvidia-smi

CUDA_VISIBLE_DEVICES=2 python scripts/alphagenome_rna_seq11_head_smoke.py \
  --dataset-dir alphagenome_custom/datasets/rna_seq_npz_pilot_train \
  --valid-dataset-dir alphagenome_custom/datasets/rna_seq_npz_valid \
  --weights weights/alphagenome_pytorch/model_all_folds.safetensors \
  --batch-size 1 \
  --max-steps 5 \
  --max-valid-examples 2 \
  --max-valid-batches 2 \
  --embedding-resolution 128 \
  --device cuda
```

- Input data: `alphagenome_custom/datasets/rna_seq_npz_pilot_train`; first 2 examples from `alphagenome_custom/datasets/rna_seq_npz_valid`; `weights/alphagenome_pytorch/model_all_folds.safetensors`
- Output path: Interactive terminal output; no checkpoint, prediction file, or log file was written.
- Result summary: Completed successfully. The run trained only the 11-track adapter head for 5 smoke-test steps and evaluated masked MSE on 2 validation examples.
- Verification: Observed `n_examples=4`, `n_valid_examples=2`, `n_tracks=11`, `base_parameters=450452613`, `trainable_parameters=33803`, `head_parameters=33803`, repeated `prediction_shape=1x11x1048576`, `base_has_grad=False` at every train step, train losses `2.86234`, `2.83171`, `1.74204`, `2.40243`, and `1.86177`, `valid_batches=2`, `valid_loss=1.61174`, and peak `cuda_max_memory_allocated_mb=17561.2`.
- Failures or warnings: This is a smoke/provisional engineering run only. The adapter still uses `organism_index=0` as a frozen trunk selector, not as a verified C. elegans organism embedding. The validation set was restricted to 2 examples.
- Next actions: If continuing, promote this smoke script into a small fine-tuning script with explicit log output, optional checkpointing to an ignored path, full validation pass, and clearer experiment configuration capture.
- Claim status: verified

## 2026-05-14 - Frozen AlphaGenome 11-Track RNA-seq Adapter Smoke Test on HY-GPU GPU 2

- Run type: smoke test
- Purpose: Verify that a frozen `alphagenome-pytorch` trunk can feed a minimal C. elegans 11-track RNA-seq adapter head using the existing pilot NPZ dataloader, without modifying the full AlphaGenome output heads or saving checkpoints.
- Git commit: `f36a2c816fbe780918a1e67cdd2bcc1c1646a85e`
- Branch: `setup/agent-maintenance`
- Host: `HY-GPU`
- Slurm job ID: Not applicable; HY-GPU is a non-Slurm server
- Slurm request: Not applicable
- Environment: Conda environment `alphagenome`; `alphagenome-pytorch` `0.3.1`; PyTorch `2.11.0+cu128`; PyTorch CUDA `12.8`; driver CUDA `13.2`; `CUDA_VISIBLE_DEVICES=2`
- Command:

```bash
conda activate alphagenome
cd /home/zelinli6/Alphagenome

hostname
pwd
git rev-parse HEAD
git status --short --branch
nvidia-smi

CUDA_VISIBLE_DEVICES=2 python scripts/alphagenome_rna_seq11_head_smoke.py \
  --dataset-dir alphagenome_custom/datasets/rna_seq_npz_pilot_train \
  --weights weights/alphagenome_pytorch/model_all_folds.safetensors \
  --batch-size 1 \
  --max-steps 1 \
  --embedding-resolution 128 \
  --device cuda
```

- Input data: `alphagenome_custom/datasets/rna_seq_npz_pilot_train`; `weights/alphagenome_pytorch/model_all_folds.safetensors`
- Output path: Interactive terminal output; no checkpoint, prediction file, or log file was written.
- Result summary: Completed successfully. The script loaded the local converted AlphaGenome checkpoint, froze the base model, trained only a new 11-track Conv1d adapter head for one optimizer step, and upsampled 128 bp embeddings to the full NPZ target length.
- Verification: Observed `n_examples=4`, `n_tracks=11`, `base_parameters=450452613`, `trainable_parameters=33803`, `head_parameters=33803`, `dna_sequence_shape=1x4x1048576`, `rna_seq_shape=1x11x1048576`, `prediction_shape=1x11x1048576`, `loss=2.87095`, `base_has_grad=False`, `head_weight_grad_norm=5.65389`, and `cuda_max_memory_allocated_mb=17516.8`.
- Failures or warnings: The adapter uses `organism_index=0` only as a frozen AlphaGenome trunk selector; this is not a verified C. elegans organism embedding. This was a smoke test, not a biological result.
- Next actions: Review the adapter script, then decide whether to add validation-loop support and run a short train/valid pilot with recorded logs before any longer fine-tuning.
- Claim status: verified

## 2026-05-14 - AlphaGenome Converted Weight Load Sanity Check on HY-GPU

- Run type: sanity check
- Purpose: Verify that the converted `alphagenome-pytorch` all-folds safetensors checkpoint can be loaded on HY-GPU GPU 2 without running training.
- Git commit: Not recorded in terminal output
- Branch: `setup/agent-maintenance`
- Host: `HY-GPU`
- Slurm job ID: Not applicable; HY-GPU is a non-Slurm server
- Slurm request: Not applicable
- Environment: Conda environment `alphagenome`; `alphagenome-pytorch` `0.3.1`; PyTorch `2.11.0+cu128`; CUDA `12.8`; `CUDA_VISIBLE_DEVICES=2`
- Command:

```bash
conda activate alphagenome
cd /home/zelinli6/Alphagenome

python -m pip install -U "huggingface_hub[cli]"

mkdir -p weights/alphagenome_pytorch
hf download gtca/alphagenome_pytorch model_all_folds.safetensors \
  --local-dir weights/alphagenome_pytorch

CUDA_VISIBLE_DEVICES=2 python - <<'PY'
import torch
from alphagenome_pytorch import AlphaGenome

path = "weights/alphagenome_pytorch/model_all_folds.safetensors"

print("cuda_available", torch.cuda.is_available())
print("device", torch.cuda.get_device_name(0) if torch.cuda.is_available() else "no cuda")
print("loading", path)

model = AlphaGenome.from_pretrained(path, device="cuda")
model.eval()

n_params = sum(p.numel() for p in model.parameters())
print("loaded", type(model).__name__)
print("n_params", n_params)
print("device_first_param", next(model.parameters()).device)
PY
```

- Input data: `weights/alphagenome_pytorch/model_all_folds.safetensors`
- Output path: Interactive terminal output; no log file was captured for this run. The checkpoint is stored in the ignored `weights/` directory.
- Result summary: Completed successfully. The converted all-folds AlphaGenome checkpoint loaded on CUDA without running training.
- Verification: Observed `cuda_available=True`, `device=NVIDIA A100 80GB PCIe`, `loaded=AlphaGenome`, `n_params=450452613`, and `device_first_param=cuda:0`.
- Failures or warnings: None reported in the provided terminal output.
- Next actions: Inspect model methods and output heads, then design a minimal C. elegans 11-track RNA-seq head adaptation plan before making code changes.
- Claim status: verified

## 2026-05-14 - alphagenome-pytorch Import Sanity Check on HY-GPU

- Run type: sanity check
- Purpose: Verify that `alphagenome-pytorch` installs and imports in the HY-GPU `alphagenome` conda environment without downloading model weights or running training.
- Git commit: Not recorded in terminal output
- Branch: `setup/agent-maintenance`
- Host: `HY-GPU`
- Slurm job ID: Not applicable; HY-GPU is a non-Slurm server
- Slurm request: Not applicable
- Environment: Conda environment `alphagenome`; `alphagenome-pytorch` `0.3.1`; PyTorch `2.11.0+cu128`; CUDA `12.8`; `CUDA_VISIBLE_DEVICES=2`
- Command:

```bash
conda activate alphagenome
cd /home/zelinli6/Alphagenome

CUDA_VISIBLE_DEVICES=2 python - <<'PY'
import importlib.metadata as md
import torch

print("python_package_alphagenome_pytorch", md.version("alphagenome-pytorch"))
print("torch", torch.__version__)
print("cuda_available", torch.cuda.is_available())
print("torch_cuda", torch.version.cuda)
print("device_count", torch.cuda.device_count())
print("device", torch.cuda.get_device_name(0) if torch.cuda.is_available() else "no cuda")

import alphagenome_pytorch
print("import_alphagenome_pytorch", "ok")

from alphagenome_pytorch import AlphaGenome
print("import_AlphaGenome", "ok")
PY
```

- Input data: None
- Output path: Interactive terminal output; no log file was captured for this run
- Result summary: Completed successfully. `alphagenome-pytorch` imported and `AlphaGenome` imported successfully. No model weights were downloaded or loaded.
- Verification: Observed `python_package_alphagenome_pytorch=0.3.1`, `torch=2.11.0+cu128`, `cuda_available=True`, `torch_cuda=12.8`, `device_count=1`, `device=NVIDIA A100 80GB PCIe`, `import_alphagenome_pytorch=ok`, and `import_AlphaGenome=ok`.
- Failures or warnings: None reported in the provided terminal output.
- Next actions: Inspect official `alphagenome-pytorch` demo and weight-loading API before downloading any safetensors weights.
- Claim status: verified

## 2026-05-14 - Full Train NPZ PyTorch CUDA Smoke Test on HY-GPU GPU 2

- Run type: smoke test
- Purpose: Verify that the full train NPZ dataset can be read on the HY-GPU non-Slurm server and that the PyTorch dataloader, tiny Conv1d model, forward pass, backward pass, and optimizer step work for two CUDA steps.
- Git commit: Not recorded in terminal output
- Branch: `setup/agent-maintenance`
- Host: `HY-GPU`
- Slurm job ID: Not applicable; HY-GPU is a non-Slurm server
- Slurm request: Not applicable
- Environment: Conda environment `alphagenome`; CUDA device restricted with `CUDA_VISIBLE_DEVICES=2`; PyTorch version previously observed as `2.11.0+cu128` with CUDA `12.8`
- Command:

```bash
conda activate alphagenome
cd /home/zelinli6/Alphagenome

CUDA_VISIBLE_DEVICES=2 python scripts/torch_smoke_train.py \
  --dataset-dir alphagenome_custom/datasets/rna_seq_npz_train \
  --batch-size 1 \
  --max-steps 2 \
  --hidden-channels 16 \
  --device auto
```

- Input data: `alphagenome_custom/datasets/rna_seq_npz_train`
- Output path: Interactive terminal output; no log file was captured for this run
- Result summary: Completed successfully. The full train NPZ dataset loaded and the tiny PyTorch model completed two CUDA forward/backward/optimizer steps using the allowed HY-GPU GPU 2 policy.
- Verification: Observed `n_examples=116`, `n_tracks=11`, `device=cuda`, `dna_sequence_shape=1x4x1048576`, `rna_seq_shape=1x11x1048576`, `prediction_shape=1x11x1048576`, `step=1`, `step=2`, and loss values `4.7324` and `2.84982`.
- Failures or warnings: None reported in the provided terminal output.
- Next actions: Verify that valid and test NPZ datasets are present on HY-GPU, then prepare the next stage for `alphagenome-pytorch` installation and official model sanity checks.
- Claim status: verified

## 2026-05-14 - Pilot NPZ PyTorch CUDA Smoke Test on HY-GPU GPU 2

- Run type: smoke test
- Purpose: Verify that the HY-GPU non-Slurm server can run the pilot NPZ PyTorch smoke test on the allowed GPU 2 device policy.
- Git commit: Not recorded in terminal output
- Branch: `setup/agent-maintenance`
- Host: `HY-GPU`
- Slurm job ID: Not applicable; HY-GPU is a non-Slurm server
- Slurm request: Not applicable
- Environment: Conda environment `alphagenome`; CUDA device restricted with `CUDA_VISIBLE_DEVICES=2`; PyTorch version not recorded in terminal output
- Command:

```bash
conda activate alphagenome
cd /home/zelinli6/Alphagenome

CUDA_VISIBLE_DEVICES=2 python scripts/torch_smoke_train.py \
  --dataset-dir alphagenome_custom/datasets/rna_seq_npz_pilot_train \
  --batch-size 1 \
  --max-steps 2 \
  --hidden-channels 16 \
  --device auto
```

- Input data: `alphagenome_custom/datasets/rna_seq_npz_pilot_train`
- Output path: Interactive terminal output; no log file was captured for this run
- Result summary: Completed successfully. The pilot NPZ dataset loaded and the tiny PyTorch model completed two CUDA forward/backward/optimizer steps using the allowed HY-GPU GPU 2 policy.
- Verification: Observed `n_examples=4`, `n_tracks=11`, `device=cuda`, `dna_sequence_shape=1x4x1048576`, `rna_seq_shape=1x11x1048576`, `prediction_shape=1x11x1048576`, `step=1`, `step=2`, and loss values `2.35111` and `3.45807`.
- Failures or warnings: None reported in the provided terminal output.
- Next actions: Record the exact PyTorch/CUDA package versions on HY-GPU, then transfer the full NPZ train/valid/test datasets if continuing on HY-GPU.
- Claim status: verified

## 2026-05-14 - Pilot NPZ PyTorch CUDA Smoke Test on HY-GPU

- Run type: smoke test
- Purpose: Verify that the HY-GPU non-Slurm server can read the existing pilot NPZ dataset and run the PyTorch dataloader, tiny Conv1d model, forward pass, backward pass, and optimizer step for two steps on CUDA.
- Git commit: Not recorded in terminal output
- Branch: `setup/agent-maintenance`
- Host: `HY-GPU`
- Slurm job ID: Not applicable; HY-GPU is a non-Slurm server
- Slurm request: Not applicable
- Environment: Conda environment `alphagenome`; CUDA device restricted with `CUDA_VISIBLE_DEVICES=0`; PyTorch version not recorded in terminal output
- Command:

```bash
conda activate alphagenome
cd /home/zelinli6/Alphagenome

CUDA_VISIBLE_DEVICES=0 python scripts/torch_smoke_train.py \
  --dataset-dir alphagenome_custom/datasets/rna_seq_npz_pilot_train \
  --batch-size 1 \
  --max-steps 2 \
  --hidden-channels 16 \
  --device auto
```

- Input data: `alphagenome_custom/datasets/rna_seq_npz_pilot_train`
- Output path: Interactive terminal output; no log file was captured for this run
- Result summary: Completed successfully. The pilot NPZ dataset loaded and the tiny PyTorch model completed two CUDA forward/backward/optimizer steps.
- Verification: Observed `n_examples=4`, `n_tracks=11`, `device=cuda`, `dna_sequence_shape=1x4x1048576`, `rna_seq_shape=1x11x1048576`, `prediction_shape=1x11x1048576`, `step=1`, `step=2`, and loss values `3.68784` and `2.89785`.
- Failures or warnings: None reported in the provided terminal output.
- Next actions: Record the exact PyTorch/CUDA package versions on HY-GPU, add a non-Slurm HY-GPU policy to `AGENTS.md`, then transfer the full NPZ train/valid/test datasets if continuing on HY-GPU.
- Claim status: verified

## 2026-05-14 - Pilot NPZ PyTorch Smoke Test on GPU1

- Run type: smoke test
- Purpose: Verify that a compute node can read the existing pilot NPZ dataset and run the PyTorch dataloader, tiny Conv1d model, forward pass, backward pass, and optimizer step for two steps.
- Git commit: `f0083bc`
- Branch: `setup/agent-maintenance`
- Host: `gpu1`
- Slurm job ID: `57458`
- Slurm request: `gpu1`, 1 GPU, 8 CPUs, 2 smoke-test steps
- Environment: Python `3.13.9`; PyTorch `2.11.0+cu130`; node driver reports CUDA `12.2`
- Command:

```bash
mkdir -p logs/slurm
sbatch -p gpu1 --job-name=ag_torch_smoke_gpu1 --gres=gpu:1 --cpus-per-task=8 scripts/slurm_smoke_torch.sh
```

- Input data: `alphagenome_custom/datasets/rna_seq_npz_pilot_train`
- Output path: `logs/slurm/ag_torch_smoke_gpu1-57458.out` and `logs/slurm/ag_torch_smoke_gpu1-57458.err`
- Result summary: Completed with exit code `0:0`. The pilot NPZ dataset loaded successfully and the tiny PyTorch model completed two forward/backward/optimizer steps. The run used CPU because CUDA was not available to PyTorch in this environment.
- Verification: Observed `n_examples=4`, `n_tracks=11`, `dna_sequence_shape=1x4x1048576`, `rna_seq_shape=1x11x1048576`, `prediction_shape=1x11x1048576`, `step=1`, `step=2`, and loss values `3.47655` and `2.51749`.
- Failures or warnings: PyTorch reported `cuda_available=False`. The error log warned that the NVIDIA driver was too old for the installed PyTorch CUDA build. `nvidia-smi` showed driver `535.183.01` and CUDA `12.2`, while PyTorch was built for CUDA `13.0`.
- Next actions: Wait for the preferred `gpu2` A100 80GB smoke test or prepare a PyTorch environment compatible with the cluster driver before GPU fine-tuning.
- Claim status: verified

## 2026-05-14 - Pilot NPZ PyTorch GPU Smoke Test (Pending)

- Run type: smoke test
- Purpose: Verify that a GPU compute node can read the existing pilot NPZ dataset and run the PyTorch dataloader, tiny Conv1d model, forward pass, backward pass, and optimizer step for two steps.
- Git commit: `f0083bc`
- Branch: `setup/agent-maintenance`
- Host: Pending
- Slurm job ID: `57454`
- Slurm request: `gpu2`, 1 GPU, 16 CPUs, 2 smoke-test steps
- Environment: Pending
- Command:

```bash
mkdir -p logs/slurm
sbatch scripts/slurm_smoke_torch.sh
```

- Input data: `alphagenome_custom/datasets/rna_seq_npz_pilot_train`
- Output path: `logs/slurm/ag_torch_smoke-<job_id>.out` and `logs/slurm/ag_torch_smoke-<job_id>.err`
- Result summary: Pending
- Verification: Pending
- Failures or warnings: Pending
- Next actions: Submit only after explicit user approval; inspect Slurm logs and record tensor shapes, device, loss lines, and any errors.
- Claim status: unverified
