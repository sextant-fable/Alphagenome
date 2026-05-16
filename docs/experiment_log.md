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
