#!/bin/bash
# Minimal Slurm smoke test for the custom C. elegans RNA-seq NPZ PyTorch loader.
#
# Before submitting, create the ignored log directory from the repository root:
#   mkdir -p logs/slurm

#SBATCH --job-name=ag_torch_smoke
#SBATCH --partition=gpu2
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=16
#SBATCH --output=logs/slurm/%x-%j.out
#SBATCH --error=logs/slurm/%x-%j.err

set -euo pipefail

cd "${SLURM_SUBMIT_DIR:-$PWD}"

echo "date	$(date -Is)"
echo "hostname	$(hostname)"
echo "current_directory	$(pwd)"
echo "git_commit	$(git rev-parse --short HEAD 2>/dev/null || echo unavailable)"

python - <<'PY'
import sys

print(f"python_executable\t{sys.executable}")
print(f"python_version\t{sys.version.split()[0]}")

try:
    import torch
except Exception as exc:
    print(f"torch_importable\tFalse")
    print(f"torch_import_error\t{type(exc).__name__}: {exc}")
else:
    print(f"torch_importable\tTrue")
    print(f"torch_version\t{torch.__version__}")
    print(f"cuda_available\t{torch.cuda.is_available()}")
    print(f"torch_cuda_version\t{torch.version.cuda}")
    if torch.cuda.is_available():
        print(f"cuda_device_count\t{torch.cuda.device_count()}")
        print(f"cuda_device_name\t{torch.cuda.get_device_name(0)}")
PY

if command -v nvidia-smi >/dev/null 2>&1; then
  nvidia-smi
else
  echo "nvidia_smi	unavailable"
fi

cmd=(
  python scripts/torch_smoke_train.py
  --dataset-dir alphagenome_custom/datasets/rna_seq_npz_pilot_train
  --batch-size 1
  --max-steps 2
  --hidden-channels 16
  --device auto
)

printf 'command'
printf '\t%q' "${cmd[@]}"
printf '\n'

"${cmd[@]}"
