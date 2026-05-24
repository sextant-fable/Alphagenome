# RNA-seq11 128 bp Adapter Phase Plan

Claim status: unverified planning record. This document records the validation
workflow and launcher shape only; it does not report completed training results.

## Scope

- Split policy: train/valid only; the held-out test split is not used.
- Outputs: model artifacts under ignored `runs/` paths and logs under ignored
  `logs/` paths.
- Default selection: validation `full-mse`, so `adapter_head_best.pt` is no
  longer selected by SmoothL1 loss when `--selection-metric full-mse` is used.
- HY-GPU execution: launcher can use GPU IDs `0,1,2,3` when explicitly allowed
  by the user.

## Phase 1

Winner stability, 9 runs:

- Heads: `residual-conv3`, `conv3`, `conv5`
- Seeds: `20260515`, `20260522`, `20260523`
- Target: `full-log1p`
- Loss: `smooth-l1`, beta `1.0`
- LR: `1e-3`
- Steps: `5000`
- Selection: validation `full-mse`

## Phase 2

Head screen, then promote the top 2-3 by validation `full-mse`:

- `residual-conv3` hidden channels `128`, `256`, `512`
- residual scale init `0.01`, `0.1`, `1.0`
- `conv3x2`
- `conv7`
- `dilated-conv3` dilation `2`, `4`

## Phase 3

Loss sweep on the top 1-2 heads:

- `smooth-l1` beta `0.5`
- `smooth-l1` beta `2.0`
- `hybrid`: `0.5*MSE + 0.5*SmoothL1`
- MSE fine-tune from the winner checkpoint for 1000 steps at `3e-4` and `1e-4`

## Phase 4

LR schedule sweep on the strongest configuration:

- constant `1e-3`
- cosine decay
- step decay at step `3000`: `1e-3 -> 3e-4`
- step decay at step `4000`: `1e-3 -> 3e-4`

## Launcher

Dry run:

```bash
python scripts/rna_seq11_run_adaptation_phase_sweep.py \
  --phase phase1 \
  --gpus 0,1,2,3 \
  --selection-metric full-mse \
  --dry-run
```

Real Phase 1:

```bash
python scripts/rna_seq11_run_adaptation_phase_sweep.py \
  --phase phase1 \
  --gpus 0,1,2,3 \
  --selection-metric full-mse
```
