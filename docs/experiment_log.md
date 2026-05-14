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
