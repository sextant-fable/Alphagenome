# Reproducibility Notes

This project is a research workspace for AlphaGenome-style modeling on C. elegans RNA-seq data. Git tracks code, metadata, split definitions, and project notes. Git does not track the full data or model state.

## Environment Assumptions

- Work is performed on the EEHPC cluster described in `SERVER_README.md`.
- `login1` is a login node and must not be used for training or long preprocessing.
- GPU work should use Slurm.
- Python and CUDA availability may differ between login and compute nodes.
- Current lightweight PyTorch requirements are listed in `requirements-torch.txt`.
- Regenerating bigWig-derived datasets requires additional packages such as `pyBigWig`.

## Data Paths

Tracked small project assets:

- `scripts/`
- `alphagenome_custom/metadata/`
- `alphagenome_custom/intervals/`
- `remote_inventory/`
- `docs/`

Large or local assets that should not be committed directly to Git:

- `training_input_bigwig/`
- `alphagenome_custom/datasets/`
- `alphagenome_custom/tracks/`
- `alphagenome_custom/reference/`
- `shared/`
- `shared.zip`
- spreadsheet source files such as `*.xlsx`
- model weights and checkpoints such as `*.safetensors`, `*.pt`, `*.pth`, and `*.ckpt`
- logs, run outputs, and temporary results

## HY-GPU Migration State

`HY-GPU` is a separate non-Slurm GPU server with local A100 80GB GPUs. It is now the preferred execution host for CUDA smoke tests and later PyTorch fine-tuning work, while the original EEHPC server remains the source of some archived large assets.

Assets that should be present on `HY-GPU` for the current PyTorch route:

- Git-tracked code, metadata, interval files, and documentation.
- `alphagenome_custom/datasets/rna_seq_npz_pilot_train`
- `alphagenome_custom/datasets/rna_seq_npz_train`
- `alphagenome_custom/datasets/rna_seq_npz_valid`
- `alphagenome_custom/datasets/rna_seq_npz_test`
- `weights/alphagenome_pytorch/model_all_folds.safetensors`

Assets recommended to transfer to `HY-GPU` if direct BigWig/FASTA inspection or loader work is needed:

- `alphagenome_custom/reference/`
- `alphagenome_custom/tracks/rna_seq_grouped/`

Assets that can remain on the original EEHPC server unless preprocessing must be reproduced from raw inputs:

- `training_input_bigwig/`
- `shared/`
- `shared.zip`
- original spreadsheet files such as `Samples.xlsx` and `remote_inventory/Samples.drive.xlsx`

Do not copy large assets just to make the two servers identical. Transfer only the assets needed for the next approved task, and record the transfer in `docs/experiment_log.md` when it affects reproducibility.

## Current Data Contract

The tracked metadata documents the current processed data contract:

- 20 original sample-level C. elegans RNA-seq bigWigs are grouped into 11 `RNA_SEQ` tracks.
- Reference genome: C. elegans WBcel235.
- Window size: 1,048,576 bp.
- Stride: 524,288 bp.
- Split: train chromosomes `I, II, III, IV`; valid chromosome `V`; test chromosome `X`; `MtDNA` excluded.
- Logical training fields: `dna_sequence`, `rna_seq`, `rna_seq_mask`, `rna_seq_strand`, and interval coordinates.

Git alone cannot reproduce the full training state because large data files, generated NPZ datasets, reference files, and future model weights are intentionally ignored.

## Slurm Conventions

- Use `sbatch` for production runs.
- Use short `srun` commands only for controlled debugging or smoke tests.
- Record every Slurm run in `docs/experiment_log.md`.
- Record partition, GPU request, CPU request, job ID, command, Git commit, input paths, and output path.
- Start with small smoke tests before long jobs.
- Do not run GPU training directly on `login1`.

Example approved GPU smoke-test pattern:

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

## Random Seeds

When adding training or evaluation scripts, expose and record random seeds for:

- Python random state
- NumPy
- PyTorch CPU
- PyTorch CUDA
- data loader shuffling

Record the seed values in the experiment log. If a run is intentionally non-deterministic, state that explicitly.

## Input and Output Conventions

- Inputs should be read from tracked metadata plus ignored large data directories.
- New generated outputs should go to ignored directories such as `outputs/`, `runs/`, `logs/`, or `checkpoints/`.
- Do not overwrite raw data or canonical generated datasets without explicit approval.
- Prefer new timestamped output directories for experimental runs.
- For large generated artifacts, record paths and checksums when practical.

## Verification

Safe checks:

```bash
git status --short --branch
python -m py_compile scripts/*.py
hostname
sinfo
```

Data count checks:

```bash
find alphagenome_custom/datasets/rna_seq_npz_train/examples -name '*.npz' | wc -l
find alphagenome_custom/datasets/rna_seq_npz_valid/examples -name '*.npz' | wc -l
find alphagenome_custom/datasets/rna_seq_npz_test/examples -name '*.npz' | wc -l
find alphagenome_custom/tracks/rna_seq_grouped -name '*.bw' | wc -l
wc -l alphagenome_custom/metadata/track_metadata_grouped.tsv
wc -l alphagenome_custom/intervals/train.bed alphagenome_custom/intervals/valid.bed alphagenome_custom/intervals/test.bed
```
