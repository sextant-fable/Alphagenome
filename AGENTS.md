# AGENTS.md

This is the canonical project guidance file for Codex and future agents working on this repository.

## Project Scope

This is a research repository for adapting AlphaGenome-style sequence-to-track modeling to C. elegans RNA-seq data.

The current repository evidence shows:

- 20 original sample-level C. elegans bigWig RNA-seq files are kept locally and ignored by Git.
- The processed project metadata groups those samples into 11 grouped `RNA_SEQ` tracks.
- The reference genome is C. elegans WBcel235.
- Intervals are 1,048,576 bp windows with 524,288 bp stride.
- The chromosome-level split is train: `I, II, III, IV`, valid: `V`, test: `X`; `MtDNA` is excluded.
- Existing NPZ datasets are local generated assets, not Git-tracked source files.

## Research Integrity Rules

- Do not fabricate citations, DOIs, paper titles, author names, years, datasets, benchmark numbers, p-values, metrics, or experimental results.
- Mark claims as verified, unverified, inferred, abstract-only, or secondary-source-only when evidence is incomplete.
- Treat smoke-test output as environment validation, not a model result.
- Treat exploratory runs as provisional until the command, commit, data inputs, environment, and outputs are recorded.
- Do not overstate biological or model conclusions from preprocessing metadata alone.

## Repository Layout

- `scripts/`: preprocessing, QC, NPZ reading, and PyTorch smoke-test utilities.
- `alphagenome_custom/metadata/`: tracked metadata, QC summaries, processing reports, and read-back summaries.
- `alphagenome_custom/intervals/`: tracked BED files for train, valid, test, and all intervals.
- `remote_inventory/`: small tracked inventory/manifests for source files.
- `docs/`: project maintenance, reproducibility, claims, and experiment records.
- Local or ignored large assets include raw bigWigs, grouped bigWigs, NPZ datasets, references, model weights, logs, and generated run outputs.

## Data Policy

- Do not modify raw data in place.
- Do not commit raw bigWigs, grouped bigWigs, NPZ datasets, reference FASTA/GTF files, model weights, checkpoints, logs, or generated experiment outputs.
- Metadata, small manifests, split BED files, and reproducibility notes may be tracked when they are useful for auditability.
- If a generated artifact is needed for reproducibility but too large for Git, record its path, checksum if available, creation command, and expected shape/counts.
- Do not regenerate datasets unless the user explicitly approves the command and target output location.

## HPC and Slurm Policy

- `login1` is a login node. Do not run training, long preprocessing, large bigWig operations, or GPU workloads on it.
- Use Slurm for GPU work. Prefer `sbatch` for production runs and short `srun` commands for controlled smoke tests.
- Before running GPU work, check `hostname`, `sinfo`, and, inside compute jobs, `nvidia-smi`.
- Never cancel or modify Slurm jobs that do not belong to the current AlphaGenome project. Before cancelling any job, verify its `WorkDir` and `Command`, confirm that it belongs to this repository, and get explicit user approval.
- Prefer finding an available GPU partition or waiting in queue over cancelling jobs.
- Prefer A100 80GB GPUs on `gpu2` for AlphaGenome smoke tests and fine-tuning when available. If `gpu2` is blocked by queue or QoS limits, check available partitions and ask before changing the target partition.
- Record Slurm job ID, partition, CPU/GPU request, command, Git commit, input data paths, output path, and result summary.
- Recommended starting point for short GPU checks is a single GPU with `--cpus-per-task=16`.
- Do not submit long jobs, download large weights, or install large dependencies without explicit user approval.

## HY-GPU Non-Slurm Policy

- `HY-GPU` is a separate non-Slurm GPU server. `squeue`, `sbatch`, and `srun` are not expected to exist there.
- Before running any GPU command on `HY-GPU`, check `hostname`, `nvidia-smi`, current GPU processes, current directory, Git commit, and active Python environment.
- Use only GPU indices `2` and `3` on `HY-GPU` unless the user explicitly approves otherwise. Do not use GPU `0` or GPU `1`.
- For single-GPU smoke tests on `HY-GPU`, prefer `CUDA_VISIBLE_DEVICES=2`.
- For multi-GPU experiments on `HY-GPU`, use only `CUDA_VISIBLE_DEVICES=2,3` unless the user explicitly changes this policy.
- Short smoke tests may run directly on `HY-GPU` after user approval. Long runs must write logs to an ignored path such as `logs/` or `runs/`, record the command in `docs/experiment_log.md`, and avoid interrupting other users' processes.
- Do not kill or modify unrelated processes on `HY-GPU`. If GPUs are busy, report the process list and ask the user before taking action.

## Coding Rules

- Keep changes minimal, reviewable, and specific to the requested task.
- Prefer existing scripts and repository patterns over new abstractions.
- Do not rewrite preprocessing or training scripts unless the user approves that phase.
- Avoid heavy dependencies unless there is a clear, concrete need.
- Do not add MCP configuration unless there is a clear, concrete reason.
- Use English for repository files, paths, commands, comments, commit messages, and documentation unless the user explicitly requests otherwise.

## Experiment Rules

- Every real run should have an entry in `docs/experiment_log.md` or a linked run record.
- Record date, purpose, command, Git commit, Slurm job ID, input data, output path, result summary, failures, and next actions.
- Separate smoke tests, sanity checks, preprocessing runs, training runs, and evaluation runs.
- Keep failed runs in the record when they explain environment or data issues.
- Do not edit past records to make results look cleaner; append corrections instead.

## Verification Commands

Safe repository checks:

```bash
git status --short --branch
git ls-files | sort
python -m py_compile scripts/*.py
```

Data presence checks:

```bash
find alphagenome_custom/datasets/rna_seq_npz_train/examples -name '*.npz' | wc -l
find alphagenome_custom/datasets/rna_seq_npz_valid/examples -name '*.npz' | wc -l
find alphagenome_custom/datasets/rna_seq_npz_test/examples -name '*.npz' | wc -l
find alphagenome_custom/tracks/rna_seq_grouped -name '*.bw' | wc -l
wc -l alphagenome_custom/metadata/track_metadata_grouped.tsv
wc -l alphagenome_custom/intervals/train.bed alphagenome_custom/intervals/valid.bed alphagenome_custom/intervals/test.bed
```

HPC checks:

```bash
hostname
sinfo
squeue -u "$USER"
```

Approved GPU smoke-test pattern:

Run this only after explicit user approval. It must be submitted through Slurm and must not be run directly on `login1`.

```bash
srun -p gpu2 --gres=gpu:1 --cpus-per-task=16 \
  python scripts/torch_smoke_train.py \
  --dataset-dir alphagenome_custom/datasets/rna_seq_npz_pilot_train \
  --batch-size 1 \
  --max-steps 2 \
  --hidden-channels 16 \
  --device auto
```

## Git and Commit Rules

- Check `git status --short --branch` before and after edits.
- Do not use `git add .` in this repository because large local assets are present.
- Stage files explicitly.
- Keep generated data, model weights, checkpoints, logs, and local settings out of Git.
- Use concise English commit messages that describe the maintenance or research change.
- Prefer small commits that separate documentation, scripts, metadata updates, and experiment records.
