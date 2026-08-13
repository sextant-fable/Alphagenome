# AGENTS.md

This is the canonical project guidance file for Codex and future agents working on this repository.

## Project Scope

This is a research repository for adapting AlphaGenome-style sequence-to-track modeling to C. elegans RNA-seq data.

The repository contains two separate research lineages:

- Legacy v1 is frozen. It uses 20 provided sample-level bigWigs grouped into 11 `RNA_SEQ` tracks, 1,048,576 bp windows with 524,288 bp stride, and local generated NPZ datasets. It is exploratory evidence, not a paper-faithful AlphaGenome reproduction. Do not generate new legacy datasets or training runs.
- v2 is the active workflow. It replaces the legacy route with audited source reprocessing, manifest-driven bigWig loading, blocked cross-validation, and a separately gated final evaluation.

The active v2 evidence shows:

- The reference genome is C. elegans WBcel235.
- The provenance inventory contains 485 accessions: 482 verified RNA-seq runs and three excluded ChIP-seq runs.
- Uniform reprocessing produces 241 biological RNA-seq groups.
- The loader reads bigWigs on demand at 1 bp and 128 bp resolution; v2 must not build a monolithic NPZ dataset.
- Development uses five leave-one-chromosome-out folds across chromosomes `I` through `V`.
- Chromosome `X` is behind the one-time G5 final-test gate. It was exposed during legacy work and is not pristine for the project as a whole.
- Model comparison covers frozen-trunk head A, worm-embedding/LoRA model B, and matched from-scratch model C.

As recorded on 2026-07-20, P6B was reopened because the original validation matrix was narrower than the locked research objective. The original selection is provisional and its final-test lock is superseded. Always read the live state with `python -m scripts.v2_phase_controller status` and the append-only `docs/v2_execution_log.md` before describing current progress.

## Research Integrity Rules

- Do not fabricate citations, DOIs, paper titles, author names, years, datasets, benchmark numbers, p-values, metrics, or experimental results.
- Mark claims as verified, unverified, inferred, abstract-only, or secondary-source-only when evidence is incomplete.
- Treat smoke-test output as environment validation, not a model result.
- Treat exploratory runs as provisional until the command, commit, data inputs, environment, and outputs are recorded.
- Do not overstate biological or model conclusions from preprocessing metadata alone.
- Do not present a provisional, failed, reopened, or superseded phase result as final.
- Do not read chromosome `X`, approve G5, or weaken the test lock unless the user explicitly authorizes the exact one-time v2 final-test scope.

## Repository Layout

- `README.md`: project entry point, current lineage summary, and navigation.
- `scripts/`: legacy utilities plus the v2 controller, reprocessing, loading, training, evaluation, and review modules.
- `tests/`: v2 unit and regression tests.
- `alphagenome_custom/metadata/`: tracked legacy metadata and read-back summaries.
- `alphagenome_custom/metadata/v2/`: v2 manifests, state, specifications, summaries, and phase audits.
- `alphagenome_custom/intervals/`: tracked legacy BED files and v2 fold/test split definitions.
- `remote_inventory/`: small tracked source-file inventories and manifests.
- `docs/`: execution plans, append-only logs, reproducibility notes, claims policy, and research reports.
- `SERVER_README.md`: EEHPC/Slurm reference for the archival host, not instructions for HY-GPU.
- Local or ignored large assets include raw and grouped bigWigs, NPZ datasets, references, model weights, checkpoints, logs, and generated run outputs.

## Data Policy

- Do not modify raw data in place.
- Do not commit raw bigWigs, grouped bigWigs, NPZ datasets, reference FASTA/GTF files, model weights, checkpoints, logs, or generated experiment outputs.
- Metadata, small manifests, split BED files, and reproducibility notes may be tracked when they are useful for auditability.
- If a generated artifact is needed for reproducibility but too large for Git, record its path, checksum if available, creation command, and expected shape/counts.
- Do not regenerate datasets unless the user explicitly approves the command and target output location.
- Do not replace the v2 dynamic loader with a monolithic dataset export.
- Do not infer that an ignored output is valid merely because it exists; verify it against its audit record and expected hash or source signature.

## Controlled v2 Workflow

- The phase sequence is `P0 -> P1 -> P2 -> P3A -> P3B -> P4 -> P5 -> P6A -> P6B -> P6C -> P7 -> P8 -> P9 -> P10`.
- Formal phase work must go through `scripts/v2_phase_controller.py`; do not bypass the controller by invoking training or evaluation entry points directly.
- G1 protects large downloads and realignment, G2 protects scientific manifest acceptance, G3 protects large generated data, G4 protects GPU experiments, and G5 protects the single chromosome-X evaluation.
- Scoped approvals are records, not reusable booleans. A prior approval does not authorize a different scope.
- The canonical state is `alphagenome_custom/metadata/v2/execution_state.json`; do not hand-edit it while a controller phase is running.
- Phase reviews and corrections are append-only. Preserve failed, original, superseded, and reopened evidence.
- P6C may consume chromosome `X` only once, after P6B locks one non-superseded checkpoint and the user separately approves G5.
- P7 is a one-off, user-authorized post-completion port of the legacy 1bp-B architecture to all 241 v2 tracks. It is development-only, must retain the P6C lock/report hashes, must not read the final-test block, and cannot replace the P6B selection or P6C result.
- P8 is a one-off, user-authorized development-only B-noLoRA ablation on fold 1 and seed `20260714`. It retains the C. elegans embedding and dual-resolution RNA head, removes only LoRA, must retain the P6C lock/report hashes, must not read the final-test block, and cannot replace the P6B selection or P6C result.
- P9 is the user-authorized development-only submission-evidence matrix. It registers 90 component records, reusing 30 P6B A/B records and one P8 no-LoRA record while training 59 new jobs. It requires `G4:p9_submission_evidence_matrix`, uses only fold validation, and preserves P6B/P6C evidence and final-test lock/report hashes.
- P10 is the user-authorized DPY-27 internal biological application. It starts automatically in `PENDING` only after P9 passes, then requires `G4:p10_dpy27_internal_application`. It reuses the 15 B/paper CV checkpoints on their matching validation cores, is explicitly not an unseen-condition benchmark, and cannot pass review until all source-data and figure-QA artifacts are complete.

## EEHPC Archival Host Policy

- The EEHPC/login1 checkout is archival and source-data only unless the user explicitly reactivates work there.
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

- `HY-GPU` (`hy8`) is the active development host and is a non-Slurm GPU server. `squeue`, `sbatch`, and `srun` are not expected to exist there.
- Before running any GPU command on `HY-GPU`, check `hostname`, `nvidia-smi`, current GPU processes, current directory, Git commit, and active Python environment.
- Default multi-GPU work on `HY-GPU` should use GPU indices `2` and `3`.
- Use all four GPUs on `HY-GPU` only when the user explicitly requests four-GPU work; in that case use `CUDA_VISIBLE_DEVICES=0,1,2,3` after checking current GPU processes.
- Do not use GPU `0` or GPU `1` unless the user explicitly requests four-GPU work.
- For single-GPU smoke tests on `HY-GPU`, prefer `CUDA_VISIBLE_DEVICES=2`.
- For multi-GPU experiments on `HY-GPU`, use `CUDA_VISIBLE_DEVICES=2,3` by default. Use `CUDA_VISIBLE_DEVICES=0,1,2,3` only for explicit four-GPU requests.
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
- Record date, purpose, command, Git commit, host, controller phase, approval scope, Slurm job ID or physical GPU index as applicable, input data, output path, result summary, failures, and next actions.
- Separate smoke tests, sanity checks, preprocessing runs, training runs, and evaluation runs.
- Keep failed runs in the record when they explain environment or data issues.
- Do not edit past records to make results look cleaner; append corrections instead.
- Treat the original P6B matrix and selection as provisional unless the reopened P6B amendment passes review.

## Verification Commands

Safe repository checks:

```bash
git status --short --branch
git ls-files | sort
python -m scripts.v2_phase_controller status
python -m py_compile scripts/*.py
python -m pytest -q
```

Legacy v1 data presence checks, only when inspecting the frozen lineage:

```bash
find alphagenome_custom/datasets/rna_seq_npz_train/examples -name '*.npz' | wc -l
find alphagenome_custom/datasets/rna_seq_npz_valid/examples -name '*.npz' | wc -l
find alphagenome_custom/datasets/rna_seq_npz_test/examples -name '*.npz' | wc -l
find alphagenome_custom/tracks/rna_seq_grouped -name '*.bw' | wc -l
wc -l alphagenome_custom/metadata/track_metadata_grouped.tsv
wc -l alphagenome_custom/intervals/train.bed alphagenome_custom/intervals/valid.bed alphagenome_custom/intervals/test.bed
```

HY-GPU checks:

```bash
hostname
nvidia-smi
```

EEHPC checks, only on the archival cluster:

```bash
hostname
sinfo
squeue -u "$USER"
```

Do not reuse the legacy NPZ smoke command as a formal v2 experiment. Approved v2 GPU work must follow the controller, G4 scope, live GPU availability check, and v2 run record.

## Git and Commit Rules

- Check `git status --short --branch` before and after edits.
- Do not use `git add .` in this repository because large local assets are present.
- Stage files explicitly.
- Keep generated data, model weights, checkpoints, logs, and local settings out of Git.
- Use concise English commit messages that describe the maintenance or research change.
- Prefer small commits that separate documentation, scripts, metadata updates, and experiment records.
- Never amend, reset, force-push, or otherwise rewrite published research history unless the user explicitly requests it and the consequences have been reviewed.

## Git Coordination

- The active development host is `HY-GPU` at `/home/zelinli6/Alphagenome`.
- The old EEHPC/login1 checkout at `/home/zelinli6/Alphagenome` is now archival/source-data only.
- Both checkouts use the same GitHub remote: `https://github.com/sextant-fable/Alphagenome.git`.
- Work on branch `setup/agent-maintenance`.
- On `HY-GPU`, after approved file edits:
  - run `git status --short --branch`
  - run `git add <explicit files>`
  - run `git commit -m "English commit message"`
  - run `git push`
- Before pushing, inspect the ahead/behind count and the commits that would be published. Do not assume a small current edit means the outgoing history is small.
- If a push would publish accumulated commits outside the current task, report the count and require explicit user confirmation before pushing them together.
- On `login1`, do not make new edits unless explicitly requested. Use only `git pull` to sync documentation/code from GitHub.
- Do not edit the same files on both hosts at the same time.
- Do not use `git add .`.
- Do not commit ignored large assets such as datasets, bigWigs, references, weights, logs, or runs.
