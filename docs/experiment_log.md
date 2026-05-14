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

## 2026-05-14 - Pilot NPZ PyTorch GPU Smoke Test (Pending)

- Run type: smoke test
- Purpose: Verify that a GPU compute node can read the existing pilot NPZ dataset and run the PyTorch dataloader, tiny Conv1d model, forward pass, backward pass, and optimizer step for two steps.
- Git commit: TODO
- Branch: `setup/agent-maintenance`
- Host: Pending
- Slurm job ID: Pending
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
