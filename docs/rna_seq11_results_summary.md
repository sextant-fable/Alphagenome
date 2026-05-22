# C. elegans RNA-seq 11-Track AlphaGenome Adapter Summary

Date: 2026-05-22

This document summarizes the current frozen AlphaGenome adapter experiment for
custom C. elegans 11-track RNA-seq prediction. Full command records and run
details are in `docs/experiment_log.md`.

For a collaborator-facing synthesis with experimental design rationale,
interpretation, and next-step recommendations, see
`docs/rna_seq11_collaborator_update_20260518.md`.

## Status

- Current selected model: frozen AlphaGenome trunk plus 11-track RNA-seq adapter.
- Selection rule: choose the checkpoint with lowest full-validation MSE.
- Selected checkpoint: `runs/rna_seq11_adapter_formal_20260515_5000steps_lr3e-4_128bp/adapter_head_best.pt`
- Selected step: 4250
- Best validation MSE: `1.0609935`
- Held-out test has now been evaluated once for this selected model.
- Do not use the test split for additional model selection or tuning.
- A follow-up 128 bp last-block LoRA scheme C feasibility round did not improve
  validation MSE, so the selected model remains unchanged.
- Follow-up 128 bp GenomeTracksHead experiments with model-space MSE and with
  AlphaGenome-style Poisson + multinomial count loss did not beat the selected
  legacy 128 bp linear adapter under the common 128 bp binned `log1p(mean raw)`
  audit metric, so the selected model remains unchanged.
- Validation-only warm-start count and hybrid loss sweeps also did not improve
  the best 128 bp GenomeTracksHead checkpoint, so 1 bp + 128 bp multi-resolution
  GenomeTracksHead training is not recommended from these objectives yet.
- The next 128 bp log1p-MSE sweep tooling is implemented and smoke-validated,
  but the long validation-only sweep has not yet been launched. This keeps the
  selected model unchanged.

## Data

- Species/reference: C. elegans WBcel235.
- Input intervals: 1,048,576 bp windows with 524,288 bp stride.
- Terminology note: 1,048,576 bp is the sequence window length, not the model
  resolution. In the current adapter run, the frozen AlphaGenome representation
  is taken at 128 bp resolution and then interpolated back to the 1,048,576
  output positions for comparison with the target arrays.
- Chromosome split:
  - train: chromosomes `I`, `II`, `III`, `IV`
  - valid: chromosome `V`
  - test: chromosome `X`
  - `MtDNA` excluded
- NPZ examples:
  - train: `116`
  - valid: `39`
  - test: `33`
- Targets:
  - 20 original local sample-level RNA-seq bigWig files were grouped into 11
    `RNA_SEQ` tracks.
  - Targets are loaded as `rna_seq` arrays and transformed with `log1p` after
    clipping negative values to zero.
  - Tensor shapes during training/evaluation:
    - DNA input: `[B, 4, 1048576]`
    - RNA-seq target: `[B, 11, 1048576]`
    - RNA-seq mask: `[B, 11, 1]`

## Model

The adapter model is implemented in `scripts/alphagenome_rna_seq11_adapter.py`.

Construction:

1. Load the converted PyTorch AlphaGenome weights with
   `AlphaGenome.from_pretrained(weights/alphagenome_pytorch/model_all_folds.safetensors)`.
2. Freeze all AlphaGenome base-model parameters with `requires_grad=False`.
3. Run the frozen model in `torch.no_grad()` via `base_model.encode(...)`.
4. Use 128 bp embeddings: `embeddings_128bp`.
5. Apply a trainable 1x1 Conv1d head from 3072 embedding channels to 11 RNA-seq
   tracks.
6. Linearly interpolate the 128 bp-resolution adapter output back to the full
   1,048,576 bp target length.

So the evaluated prediction arrays have one value per base-position slot across
the 1,048,576 bp window, but the adapter's learned signal is bottlenecked by the
128 bp AlphaGenome embedding grid.

Parameter counts for the selected setup:

- frozen AlphaGenome base parameters: `450452613`
- trainable adapter/head parameters: `33803`
- total module parameters: `450486416`

The trunk is not fine-tuned in the current result. This is adaptation from
pretrained AlphaGenome-style weights to custom C. elegans grouped RNA-seq tracks
using only the small adapter/head.

Additional 1 bp pilot variant:

- A follow-up validation-only pilot used `embeddings_1bp` instead of
  `embeddings_128bp`.
- The 1 bp adapter head is a trainable 1x1 Conv1d from 1536 embedding channels
  to 11 RNA-seq tracks and has `16907` trainable parameters.
- This variant was run with the same frozen trunk, but with two-process DDP on
  HY-GPU GPUs 2 and 3, per-process batch size `1`, gradient accumulation `4`,
  global effective batch size `8`, learning rate `1e-4`, 50-step warmup, cosine
  decay, and gradient clipping at norm `1.0`.
- The 1 bp pilot is not the selected model; the held-out test split was not used
  for it.

LoRA feasibility check:

- A 128 bp last-block LoRA smoke test was added after the frozen-adapter pilots.
- LoRA was applied only to `tower.blocks.8.mha` and `tower.blocks.8.mlp`, covering
  six Linear layers in the final transformer block.
- The smoke used rank `4`, alpha `8`, and added `72960` trainable LoRA parameters
  on top of the `33803` trainable RNA-seq adapter-head parameters.
- The run confirmed `base_non_lora_trainable_parameters=0`,
  `base_non_lora_has_grad=False`, `lora_has_grad=True`, and
  `head_has_grad=True` for a one-step 128 bp feasibility pass.
- This is not a validation or test result. It only establishes that a
  lightweight trunk-adaptation path can run.

LoRA scheme C continuation pilots:

- Both pilots initialized the 128 bp RNA-seq adapter head from
  `runs/rna_seq11_adapter_formal_20260515_5000steps_lr3e-4_128bp/adapter_head_best.pt`.
- Phase A froze that initialized adapter head and trained only final-block LoRA
  parameters. Original non-LoRA AlphaGenome trunk weights stayed frozen.
- Phase B trained both the final-block LoRA parameters and the initialized
  adapter head. Original non-LoRA AlphaGenome trunk weights stayed frozen.
- Both pilots used LoRA rank `4`, alpha `8`, learning rate `5e-5`, 20-step
  warmup, cosine decay, gradient clipping at norm `1.0`, batch size `1`, and
  100 total steps with full validation every 20 steps.
- Neither pilot outperformed the selected frozen-trunk 128 bp adapter checkpoint.

128 bp log1p-MSE sweep tooling:

- The legacy `conv1x1` head remains the default and preserves backward
  compatibility with the selected checkpoint.
- Additional validation-only sweep heads are now available:
  `mlp1x1`, `conv3`, `conv5`, `residual-conv1x1`, and `residual-conv3`.
- Linear-head training now supports masked MSE or SmoothL1, full-resolution
  `log1p` targets, and 128 bp binned `log1p(mean raw)` targets.
- Checkpoints record the linear head architecture, hidden channels, loss type,
  target space, and SmoothL1 beta so evaluation can reconstruct the correct
  head automatically.

## Loss And Metrics

Training loss:

- masked MSE over valid RNA-seq tracks and all base positions
- denominator: number of valid masked tracks times sequence length

Reported evaluation metrics:

- MSE: exact pointwise masked mean squared error
- MAE: exact pointwise masked mean absolute error
- Pearson: exact pointwise Pearson correlation
- Spearman: sampled per-track Spearman with deterministic sampling

For Spearman runs reported here:

- `spearman_sample_size=200000`
- `spearman_seed=20260515`

New validation diagnostics for the 128 bp log1p-MSE sweep:

- per-track MSE, MAE, Pearson, and sampled Spearman
- per-window MSE, MAE, and Pearson
- stratum metrics for `zero`, `0<log1p<=1`, `1<log1p<=3`, and `log1p>3`
- common 128 bp metrics plus upsampled full-resolution metrics for binned-target
  checkpoints

## Environment

HY-GPU non-Slurm server:

- host: `HY-GPU`
- primary selected-model GPU: `CUDA_VISIBLE_DEVICES=2`
- later 128 bp GenomeTracks follow-up runs used `CUDA_VISIBLE_DEVICES=0,1,2,3`
  after explicit user approval
- GPU: NVIDIA A100 80GB PCIe
- conda environment: `alphagenome`
- Python: `3.12.13`
- PyTorch: `2.11.0+cu128`
- PyTorch CUDA: `12.8`
- NVIDIA driver CUDA reported by `nvidia-smi`: `13.2`
- `alphagenome_pytorch`: `0.3.1`
- Git commit recorded for these runs: `f36a2c816fbe780918a1e67cdd2bcc1c1646a85e`

## Training Configuration

Formal adapter-only training command:

```bash
CUDA_VISIBLE_DEVICES=2 conda run -n alphagenome python -u scripts/alphagenome_rna_seq11_finetune.py \
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
  --device auto
```

Key settings:

- optimizer: AdamW
- learning rate: `3e-4`
- weight decay: `0.0`
- batch size: `1`
- max steps: `5000`
- validation interval: every `250` steps
- target transform: `log1p`
- embedding resolution: `128`
- organism index: `0`
- num workers: `0`
- checkpoint selection: best full-validation MSE

During training, logs repeatedly confirmed:

- `base_has_grad=False`
- prediction shape: `1x11x1048576`
- CUDA peak memory around `17665.2` MB during training

LoRA scheme C continuation commands used the same train/valid NPZ splits and
the selected 5000-step adapter checkpoint as initialization:

```bash
CUDA_VISIBLE_DEVICES=2 conda run -n alphagenome python -u scripts/alphagenome_rna_seq11_finetune.py \
  --train-dataset-dir alphagenome_custom/datasets/rna_seq_npz_train \
  --valid-dataset-dir alphagenome_custom/datasets/rna_seq_npz_valid \
  --weights weights/alphagenome_pytorch/model_all_folds.safetensors \
  --output-dir runs/rna_seq11_adapter_128bp_lora_phaseA_20260516_100steps_lr5e-5_freezehead \
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
  --device auto

CUDA_VISIBLE_DEVICES=2 conda run -n alphagenome python -u scripts/alphagenome_rna_seq11_finetune.py \
  --train-dataset-dir alphagenome_custom/datasets/rna_seq_npz_train \
  --valid-dataset-dir alphagenome_custom/datasets/rna_seq_npz_valid \
  --weights weights/alphagenome_pytorch/model_all_folds.safetensors \
  --output-dir runs/rna_seq11_adapter_128bp_lora_phaseB_20260516_100steps_lr5e-5_lora_head \
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
  --device auto
```

## Baselines

### Train-track mean baseline

Implemented in `scripts/rna_seq11_mean_baseline.py`.

- Estimate one global `log1p` mean per track from train.
- Predict the train mean for every valid base position.
- Valid MSE: `1.7282558`

### Tiny Conv1d learned baseline

Implemented through `TinyConvRnaSeqModel` in `scripts/torch_rna_seq_dataset.py`.

Architecture:

```text
Conv1d(4 -> hidden_channels, kernel_size=15, padding=7)
GELU
Conv1d(hidden_channels -> hidden_channels, kernel_size=15, padding=7)
GELU
Conv1d(hidden_channels -> 11, kernel_size=1)
```

Small 500-step baseline:

- hidden channels: `16`
- trainable parameters: `5019`
- learning rate: `3e-4`
- max steps: `500`
- valid MSE at step 500: `1.6399873`
- valid MAE: `1.0042184`
- valid Pearson: `0.19711665`

Stronger 2000-step baseline:

- hidden channels: `64`
- trainable parameters: `66123`
- learning rate: `1e-3`
- max steps: `2000`
- eval every: `200`
- best valid MSE: `1.553865` at step 400
- final valid MSE: `1.6456854`
- final valid MAE: `1.0224008`
- final valid Pearson: `0.27318148`

## Main Results

Validation-set development results:

| model | split | selection | MSE | MAE | Pearson |
|---|---:|---|---:|---:|---:|
| train-track mean baseline | valid | fixed baseline | `1.7282558` | not computed | not computed |
| tiny Conv1d h16 | valid | final step 500 | `1.6399873` | `1.0042184` | `0.19711665` |
| tiny Conv1d h64 | valid | best step 400 | `1.553865` | not computed at best step | not computed at best step |
| tiny Conv1d h64 | valid | final step 2000 | `1.6456854` | `1.0224008` | `0.27318148` |
| frozen AlphaGenome 1 bp adapter | valid | best/final step 500 DDP pilot | `1.4771859` | not computed | not computed |
| frozen AlphaGenome adapter | valid | step 500 pilot | `1.2418628` | `0.84385942` | `0.52766246` |
| frozen AlphaGenome adapter | valid | best step 4250 | `1.0609935` | `0.71866247` | `0.58809888` |
| 128 bp last-block LoRA, initialized from selected adapter, head frozen | valid | best step 20 | `1.0613562` | not computed | not computed |
| 128 bp last-block LoRA, initialized from selected adapter, LoRA+head | valid | best step 40 | `1.067338` | not computed | not computed |

Tooling validation result:

| check | split | result |
|---|---:|---|
| selected checkpoint reload after sweep-tooling changes | valid | reproduced MSE `1.0609935`, MAE `0.71866247`, Pearson `0.58809888` |
| diagnostic TSV for selected checkpoint | valid | wrote `11` track rows, `39` window rows, and `4` stratum rows |
| new head/loss/target smoke tests | train/valid mini-subsets | 1-step train, 1-example valid, checkpoint save, and reload passed |

Held-out test result for the selected checkpoint:

| model | split | checkpoint | MSE | MAE | Pearson |
|---|---:|---|---:|---:|---:|
| frozen AlphaGenome adapter | test | best valid checkpoint, step 4250 | `1.1711332` | `0.76486897` | `0.55539986` |

GenomeTracksHead follow-up audit metrics:

These use a common 128 bp binned `log1p(mean raw signal)` metric so the
GenomeTracks checkpoints and legacy linear checkpoint can be compared on the
same transformed target space. These audit metrics are separate from each
run's training loss.

| model | split | checkpoint/objective | MSE | MAE | Pearson |
|---|---:|---|---:|---:|---:|
| legacy 128 bp linear adapter | valid | best step 4250, log1p MSE objective | `1.0298867` | `0.68138634` | `0.62135428` |
| legacy 128 bp linear adapter | test | best step 4250, log1p MSE objective | `1.1380058` | `0.73726858` | `0.59136682` |
| 128 bp GenomeTracksHead | valid | best step 3000, model-space MSE objective | `2.133715` | `1.1418334` | `0.53935375` |
| 128 bp GenomeTracksHead | test | best step 3000, model-space MSE objective | `1.8786616` | `1.0444137` | `0.53915857` |
| 128 bp GenomeTracksHead | valid | best step 500, Poisson + multinomial objective | `3.0527705` | `1.4647581` | `0.49153418` |
| 128 bp GenomeTracksHead | test | best step 500, Poisson + multinomial objective | `2.5278237` | `1.2957843` | `0.5309488` |
| 128 bp GenomeTracksHead | valid | warm-start from MSE, Poisson + multinomial, positional weight 1.0 | `2.9095493` | `1.3788249` | `0.52790218` |
| 128 bp GenomeTracksHead | valid | warm-start from MSE, hybrid MSE + `1e-5` Poisson-multinomial | `2.2140787` | `1.1680553` | `0.53737598` |
| 128 bp GenomeTracksHead | valid | warm-start from MSE, hybrid MSE + `3e-5` Poisson-multinomial | `2.3104292` | `1.2004244` | `0.5351673` |
| 128 bp GenomeTracksHead | valid | warm-start from MSE, hybrid MSE + `1e-4` Poisson-multinomial | `2.504247` | `1.2616879` | `0.53141917` |

GenomeTracksHead objective-specific validation losses:

| run | objective-space validation curve |
|---|---|
| 128 bp GenomeTracksHead, model-space MSE | step 250 `1.808731`, 500 `1.7033822`, 1000 `1.6219809`, 1500 `1.5778429`, 2000 `1.5444947`, 2500 `1.5202813`, 3000 `1.5025958` |
| 128 bp GenomeTracksHead, Poisson + multinomial | step 100 `91680.721`, 200 `86910.675`, 300 `85539.633`, 400 `85080.233`, 500 `85025.399` |
| 128 bp GenomeTracksHead, warm-start Poisson + multinomial, positional weight 1.0 | step 100 `18511.789`, 200 `18495.892`, 300 `18489.356`, 400 `18485.56`, 500 `18485.124` |
| 128 bp GenomeTracksHead, warm-start hybrid, Poisson weight `1e-5` | step 100 `1.688949`, 200 `1.6883764`, 300 `1.6878839`, 400 `1.6877085`, 500 `1.6876779` |
| 128 bp GenomeTracksHead, warm-start hybrid, Poisson weight `3e-5` | step 100 `2.0615844`, 200 `2.0610508`, 300 `2.0605806`, 400 `2.0604339`, 500 `2.0604059` |
| 128 bp GenomeTracksHead, warm-start hybrid, Poisson weight `1e-4` | step 100 `3.3626393`, 200 `3.3621139`, 300 `3.3615108`, 400 `3.3613712`, 500 `3.3613398` |

Relative validation MSE comparison:

- adapter best vs train-track mean baseline:
  - MSE reduction: `0.6672623`
  - relative reduction: about `38.6%`
- adapter best vs stronger tiny Conv1d best:
  - MSE reduction: `0.4928715`
  - relative reduction: about `31.7%`
- adapter best vs 500-step adapter pilot:
  - MSE reduction: `0.1808693`
  - relative reduction: about `14.6%`

## Formal Adapter Validation Curve

Selected validation points from the 5000-step adapter-only run:

| step | valid MSE |
|---:|---:|
| 500 | `1.2141448` |
| 1000 | `1.1391212` |
| 1500 | `1.1319546` |
| 2000 | `1.1181387` |
| 2500 | `1.0989698` |
| 3000 | `1.0935162` |
| 3250 | `1.0837174` |
| 3500 | `1.1061893` |
| 3750 | `1.0840694` |
| 4000 | `1.0615692` |
| 4250 | `1.0609935` |
| 4500 | `1.0840007` |
| 4750 | `1.0764805` |
| 5000 | `1.1053824` |

The final checkpoint was worse than the best validation checkpoint, so the
selected model is `adapter_head_best.pt`, not `adapter_head.pt`.

## Held-Out Test Per-Track Metrics

Checkpoint: `runs/rna_seq11_adapter_formal_20260515_5000steps_lr3e-4_128bp/adapter_head_best.pt`

| track | name | MSE | MAE | Pearson | sampled Spearman |
|---|---|---:|---:|---:|---:|
| RNA_SEQ_001 | Intestine (end1), Embryo T0 | `1.1087508` | `0.72696165` | `0.51305557` | `0.39898333` |
| RNA_SEQ_002 | Intestine (end1), Embryo T1 | `1.2813301` | `0.80888759` | `0.51497802` | `0.4206953` |
| RNA_SEQ_003 | Intestine (end1), Embryo T2 | `1.0176844` | `0.69805531` | `0.53053846` | `0.42976539` |
| RNA_SEQ_004 | Intestine (end1), Embryo T3 | `1.2760303` | `0.8365687` | `0.53190608` | `0.42524173` |
| RNA_SEQ_005 | Intestine (end1), Embryo T4 | `1.3059002` | `0.85648933` | `0.55922813` | `0.46137527` |
| RNA_SEQ_006 | Muscle (hlh-1), Embryo T0 | `1.2404023` | `0.76799258` | `0.53393285` | `0.44893427` |
| RNA_SEQ_007 | Muscle (hlh-1), Embryo T1 | `1.062132` | `0.69213948` | `0.55428168` | `0.47665318` |
| RNA_SEQ_008 | Muscle (hlh-1), Embryo T2 | `1.1175543` | `0.71924756` | `0.57475091` | `0.49999499` |
| RNA_SEQ_009 | Muscle (hlh-1), Embryo T3 | `1.1340388` | `0.73409846` | `0.59417222` | `0.52200972` |
| RNA_SEQ_010 | Muscle (hlh-1), Embryo T4 | `1.2012671` | `0.78986845` | `0.59174095` | `0.50791787` |
| RNA_SEQ_011 | Pharynx (pha-4), Embryo T4 | `1.1373745` | `0.78324954` | `0.54218363` | `0.46518732` |

## Output Artifacts

Ignored run directories:

- formal adapter run:
  `runs/rna_seq11_adapter_formal_20260515_5000steps_lr3e-4_128bp/`
- stronger tiny Conv1d baseline:
  `runs/rna_seq11_tiny_conv_baseline_20260515_2000steps_lr1e-3_h64/`

Key ignored artifacts:

- selected checkpoint:
  `runs/rna_seq11_adapter_formal_20260515_5000steps_lr3e-4_128bp/adapter_head_best.pt`
- final checkpoint:
  `runs/rna_seq11_adapter_formal_20260515_5000steps_lr3e-4_128bp/adapter_head.pt`
- validation metrics:
  `runs/rna_seq11_adapter_formal_20260515_5000steps_lr3e-4_128bp/valid_best_pointwise_metrics.tsv`
- held-out test metrics:
  `runs/rna_seq11_adapter_formal_20260515_5000steps_lr3e-4_128bp/test_best_pointwise_metrics.tsv`

These files are generated artifacts and should not be committed to Git.

## Interpretation And Limits

Verified from current runs:

- Frozen AlphaGenome embeddings plus a small 11-track adapter outperform both
  the train-track mean baseline and the learned tiny Conv1d baseline on the
  validation split.
- The selected adapter checkpoint also gives a held-out test MSE of `1.1711332`
  and test Pearson of `0.55539986`.
- The held-out test MSE is worse than the selected validation MSE
  (`1.1711332` vs `1.0609935`), so there is a generalization gap.

Limits:

- This does not reproduce original AlphaGenome training or benchmark protocol.
- This is a species- and dataset-specific adapter experiment on custom
  C. elegans RNA-seq tracks.
- The trunk was frozen; only the small adapter/head was trained.
- Spearman values are sampled, not exact over every base.
- The held-out test split has now been used once and should not guide further
  model selection.

## Recommended Next Step

For this experiment round, keep the legacy 128 bp linear adapter as the current
benchmark. Do not start 1 bp + 128 bp GenomeTracksHead training from the current
Poisson + multinomial or hybrid objectives. Launch the staged 128 bp
log1p-MSE sweep next: seed stability, short hyperparameter screening, head
screening, binned-target screening, and promotion runs chosen only by validation
metrics. Keep the held-out test split frozen.
