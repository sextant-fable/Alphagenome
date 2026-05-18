# C. elegans RNA-seq 11-Track AlphaGenome Adapter Update

Date: 2026-05-18

This is a collaborator-facing summary of the current C. elegans RNA-seq
11-track adaptation experiments. The audit trail with exact commands is in
`docs/experiment_log.md`; the shorter running results summary is in
`docs/rna_seq11_results_summary.md`.

## Executive Summary

- The current selected model is the 128 bp frozen AlphaGenome-trunk adapter:
  `runs/rna_seq11_adapter_formal_20260515_5000steps_lr3e-4_128bp/adapter_head_best.pt`.
- The selected checkpoint is step `4250`, chosen by the lowest full-validation
  MSE.
- Best validation result: MSE `1.0609935`, MAE `0.71866247`, Pearson
  `0.58809888`.
- Held-out test was evaluated once after model selection: MSE `1.1711332`,
  MAE `0.76486897`, Pearson `0.55539986`.
- Test should now stay closed for this experiment family. Further model tuning
  should use train and validation only.
- A stronger tiny Conv1d baseline, 1 bp adapter pilots, and last-block LoRA
  pilots did not beat the selected 128 bp frozen-trunk adapter.

## Data And Split

- Species/reference: C. elegans WBcel235.
- Source signal: 20 local sample-level C. elegans RNA-seq bigWig files.
- Grouped targets: 11 grouped `RNA_SEQ` tracks.
- Windowing: 1,048,576 bp input windows with 524,288 bp stride.
- Chromosome split:
  - train: `I`, `II`, `III`, `IV`
  - valid: `V`
  - test: `X`
  - `MtDNA` excluded
- NPZ example counts:
  - train: `116`
  - valid: `39`
  - test: `33`
- Tensor shapes:
  - DNA input: `[B, 4, 1048576]`
  - RNA-seq target: `[B, 11, 1048576]`
  - RNA-seq mask: `[B, 11, 1]`
- Target transform: negative values clipped to zero, then `log1p`.

## Evaluation Discipline

The workflow used validation-first model selection:

1. Train only on the train split.
2. Select checkpoints and compare model families only on the validation split.
3. Use the test split once after the model and selection rule are fixed.

The test split has already been used once for the selected 128 bp adapter
checkpoint. New experiments should be treated as a new validation-only round
until a new selection rule and model family are finalized.

## Model Design

### Selected Adapter Model

The selected model is intentionally small:

1. Load converted PyTorch AlphaGenome weights from
   `weights/alphagenome_pytorch/model_all_folds.safetensors`.
2. Freeze the full AlphaGenome base model.
3. Run `base_model.encode(...)` under `torch.no_grad()`.
4. Take `embeddings_128bp`.
5. Train a 1x1 Conv1d head from 3072 embedding channels to 11 RNA-seq tracks.
6. Linearly interpolate the 128 bp grid back to 1,048,576 output positions.

This means the output tensor has one predicted value per target base-position
slot, but the learned adapter signal is bottlenecked by the 128 bp embedding
grid. The sequence window length is 1,048,576 bp; the selected adapter
representation resolution is 128 bp.

Selected setup parameter counts:

- Frozen AlphaGenome base parameters: `450452613`
- Trainable adapter-head parameters: `33803`
- Total module parameters: `450486416`

### Why This Design Was Chosen

The first goal was to answer a practical question: do pretrained
AlphaGenome-style sequence embeddings help on these custom C. elegans RNA-seq
tracks, even with a very small downstream dataset?

The adapter-only route is conservative. It avoids immediately fine-tuning
hundreds of millions of trunk parameters on only 116 train windows, and it gives
a clean comparison against sequence-only baselines. Once that worked, we tested
two natural extensions:

- 1 bp embeddings, to see whether finer model representation resolution helps.
- Last-block LoRA, to test whether a small amount of trunk adaptation improves
  over the already trained adapter head.

So far, neither extension improved validation performance.

## Experiment Matrix

Validation results are the main comparison surface. MSE and MAE are pointwise
masked metrics over valid RNA-seq tracks and positions. Pearson is exact
pointwise Pearson. Spearman values, where used, are sampled per track with
`spearman_sample_size=200000` and `spearman_seed=20260515`.

| Experiment | Key setup | Selection | Valid MSE | Valid MAE | Valid Pearson | Notes |
|---|---|---:|---:|---:|---:|---|
| Train-track mean baseline | One global train mean per track | fixed | `1.7282558` | not computed | not computed | Sanity baseline |
| Tiny Conv1d h16 | Scratch Conv1d, h=16, 500 steps, lr `3e-4` | final step 500 for point metrics | `1.6399873` | `1.0042184` | `0.19711665` | Best MSE was `1.6157231` at step 400 |
| Tiny Conv1d h64 | Scratch Conv1d, h=64, 2000 steps, lr `1e-3` | best step 400 | `1.553865` | not computed at best | not computed at best | Final step Pearson was `0.27318148` |
| 128 bp adapter pilot | Frozen trunk, 20 steps, lr `1e-3` | step 20 | `1.4536235` | not computed | not computed | Early feasibility |
| 128 bp adapter pilot | Frozen trunk, 100 steps, lr `3e-4` | step 100 | `1.4401001` | not computed | not computed | Improved with longer pilot |
| 128 bp adapter formal | Frozen trunk, 500 steps, lr `3e-4` | step 500 | `1.2418628` | `0.84385942` | `0.52766246` | First strong adapter result |
| 128 bp adapter formal | Frozen trunk, 5000 steps, lr `3e-4` | best step 4250 | `1.0609935` | `0.71866247` | `0.58809888` | Current selected model |
| 1 bp adapter pilot | Frozen trunk, 100 steps, lr `3e-4` | best step 80 | `1.5281759` | not computed | not computed | Much more memory and slower |
| 1 bp adapter DDP pilot | Frozen trunk, 500 steps, lr `1e-4`, warmup/cosine, grad accum 4, 2 GPUs | step 500 | `1.4771859` | not computed | not computed | Better than tiny Conv, worse than 128 bp adapter |
| 128 bp LoRA phase A | Initialize from selected head, freeze head, train final-block LoRA only, lr `5e-5` | best step 20 | `1.0613562` | not computed | not computed | Very close, did not beat selected model |
| 128 bp LoRA phase B | Initialize from selected head, train LoRA plus head, lr `5e-5` | best step 40 | `1.067338` | not computed | not computed | Degraded more clearly |

Held-out test result for the selected checkpoint:

| Model | Split | Checkpoint | MSE | MAE | Pearson |
|---|---|---|---:|---:|---:|
| 128 bp frozen-trunk adapter | test | best valid checkpoint, step 4250 | `1.1711332` | `0.76486897` | `0.55539986` |

## Selected Model Per-Track Metrics

Validation, selected step 4250:

| Track | Short label | MSE | MAE | Pearson | Sampled Spearman |
|---:|---|---:|---:|---:|---:|
| 0 | Intestine end1 Embryo T0 | `1.1183739` | `0.74297227` | `0.54028078` | `0.4215102` |
| 1 | Intestine end1 Embryo T1 | `1.2535293` | `0.80409167` | `0.54484031` | `0.43816201` |
| 2 | Intestine end1 Embryo T2 | `0.99361222` | `0.68495703` | `0.55136293` | `0.4345624` |
| 3 | Intestine end1 Embryo T3 | `1.2119518` | `0.79896007` | `0.54676353` | `0.43417157` |
| 4 | Intestine end1 Embryo T4 | `1.2211032` | `0.8137499` | `0.57354249` | `0.46537366` |
| 5 | Muscle hlh-1 Embryo T0 | `1.1029962` | `0.72716302` | `0.60321` | `0.49778837` |
| 6 | Muscle hlh-1 Embryo T1 | `0.92602274` | `0.64502917` | `0.60869882` | `0.50940104` |
| 7 | Muscle hlh-1 Embryo T2 | `0.92952961` | `0.64687189` | `0.61981902` | `0.51976099` |
| 8 | Muscle hlh-1 Embryo T3 | `0.94064081` | `0.65761898` | `0.63397316` | `0.5379481` |
| 9 | Muscle hlh-1 Embryo T4 | `0.99617012` | `0.69706668` | `0.64370303` | `0.54338671` |
| 10 | Pharynx pha-4 Embryo T4 | `0.97699841` | `0.68680646` | `0.58150481` | `0.48698593` |

Test, selected step 4250:

| Track | Short label | MSE | MAE | Pearson | Sampled Spearman |
|---:|---|---:|---:|---:|---:|
| 0 | Intestine end1 Embryo T0 | `1.1087508` | `0.72696165` | `0.51305557` | `0.39898333` |
| 1 | Intestine end1 Embryo T1 | `1.2813301` | `0.80888759` | `0.51497802` | `0.4206953` |
| 2 | Intestine end1 Embryo T2 | `1.0176844` | `0.69805531` | `0.53053846` | `0.42976539` |
| 3 | Intestine end1 Embryo T3 | `1.2760303` | `0.8365687` | `0.53190608` | `0.42524173` |
| 4 | Intestine end1 Embryo T4 | `1.3059002` | `0.85648933` | `0.55922813` | `0.46137527` |
| 5 | Muscle hlh-1 Embryo T0 | `1.2404023` | `0.76799258` | `0.53393285` | `0.44893427` |
| 6 | Muscle hlh-1 Embryo T1 | `1.062132` | `0.69213948` | `0.55428168` | `0.47665318` |
| 7 | Muscle hlh-1 Embryo T2 | `1.1175543` | `0.71924756` | `0.57475091` | `0.49999499` |
| 8 | Muscle hlh-1 Embryo T3 | `1.1340388` | `0.73409846` | `0.59417222` | `0.52200972` |
| 9 | Muscle hlh-1 Embryo T4 | `1.2012671` | `0.78986845` | `0.59174095` | `0.50791787` |
| 10 | Pharynx pha-4 Embryo T4 | `1.1373745` | `0.78324954` | `0.54218363` | `0.46518732` |

## Interpretation

The current results support a cautious but positive conclusion:

- The pretrained AlphaGenome-style trunk is useful as a frozen feature extractor
  for this small C. elegans RNA-seq adaptation problem.
- The 128 bp adapter strongly outperforms both the train-mean baseline and the
  learned tiny Conv1d sequence-only baseline on validation MSE and Pearson.
- The improvement over the stronger tiny Conv1d h64 baseline is substantial:
  valid MSE improves from `1.553865` to `1.0609935`.
- The test result is worse than validation, but not catastrophically so:
  MSE `1.1711332` and Pearson `0.55539986`. This is consistent with chromosome
  holdout being harder than random-window validation would be.
- Muscle tracks are generally easier for the selected model than intestine
  tracks in this split. On validation, muscle tracks have Pearson around
  `0.60` to `0.64`, while intestine tracks are around `0.54` to `0.57`.

The results do not justify claiming that we reproduced the original
AlphaGenome paper's training recipe or benchmark. This is a downstream
adaptation experiment with different species-specific targets, grouped tracks,
small data, a custom chromosome split, and a frozen-trunk adapter objective.

## Why 1 bp Did Not Win Yet

The 1 bp embedding path is feasible but did not improve validation performance.
Likely reasons:

- The 1 bp adapter has more detailed spatial input but fewer head input channels
  in this implementation: 1536 channels versus 3072 for 128 bp embeddings.
- It is much more memory and compute intensive, so the run budget was shorter
  and needed DDP, gradient accumulation, warmup, cosine decay, and clipping.
- The dataset is small. Finer resolution gives the optimizer more high-frequency
  target variation to chase, but not necessarily more reliable supervision.
- Our current loss is dense pointwise MSE over log1p signal. It may reward smooth
  broad trends more than exact base-level sharpness, especially after grouping
  sample-level tracks.

The practical conclusion is not that 1 bp is useless. It is that, under the
current budget and objective, the 128 bp representation is the better operating
point.

## Why LoRA Did Not Improve Yet

LoRA was tested as a minimal trunk-adaptation path. The original non-LoRA trunk
weights remained frozen; LoRA added trainable low-rank deltas in the final
transformer block only.

The scheme C pilots loaded the already selected 5000-step 128 bp adapter head
before LoRA training:

- Phase A: head frozen, train LoRA only.
- Phase B: train LoRA plus head.

Phase A was very close to the selected checkpoint but did not beat it. Phase B
degraded more. The simplest interpretation is that the selected adapter head is
already near a good validation optimum for the frozen features, and `5e-5`
continuation updates are still large enough to disturb it. If LoRA is revisited,
the next validation-only trial should use a lower learning rate such as `1e-5`
to `2e-5`, or a shorter early-stopped schedule.

## Differences From The Original AlphaGenome-Style Setup

These experiments should be described as adaptation from pretrained
AlphaGenome-style weights, not as original-model reproduction.

Important differences:

- Target space: our target is 11 grouped C. elegans RNA-seq tracks, not the full
  original multi-assay target space.
- Data scale: only 116 train windows are available after chromosome split.
- Training scope: the selected model freezes the 450M-parameter trunk and trains
  only a 33,803-parameter RNA-seq adapter head.
- Objective: dense masked MSE on log1p-transformed grouped RNA-seq signal.
- Resolution: selected adapter uses 128 bp embeddings and interpolates to the
  full 1,048,576 target length for evaluation.
- Evaluation: chromosome holdout `V` for validation and `X` for test, rather
  than any paper benchmark split.

Therefore, the strongest claim is: pretrained AlphaGenome-style sequence
representations appear useful for our C. elegans 11-track RNA-seq prediction
task under this adapter setup.

## Recommended Next Steps

1. Freeze the current test result as the reported held-out result for this
   experiment family.
2. Prepare collaborator plots from existing artifacts:
   - validation MSE by step for 5000-step 128 bp adapter
   - baseline versus adapter bar chart
   - selected model per-track valid/test Pearson and MSE
3. If continuing modeling, use validation only and prioritize:
   - lower-LR LoRA-only continuation from selected head
   - adapter regularization or early stopping rather than longer unregularized
     training
   - alternative losses or metrics better aligned to RNA-seq peak/broad-signal
     behavior
4. Keep a strict rule: do not evaluate the test split again until a genuinely
   new model-selection protocol is frozen.

## Key Run Paths

- Selected checkpoint:
  `runs/rna_seq11_adapter_formal_20260515_5000steps_lr3e-4_128bp/adapter_head_best.pt`
- Selected valid metrics:
  `runs/rna_seq11_adapter_formal_20260515_5000steps_lr3e-4_128bp/valid_best_pointwise_metrics.tsv`
- One-time test metrics:
  `runs/rna_seq11_adapter_formal_20260515_5000steps_lr3e-4_128bp/test_best_pointwise_metrics.tsv`
- Mean baseline:
  `runs/rna_seq11_adapter_formal_20260515_500steps_lr3e-4/valid_track_mean_baseline_metrics.tsv`
- Tiny Conv baselines:
  `runs/rna_seq11_tiny_conv_baseline_20260515_500steps_lr3e-4_h16/`
  and `runs/rna_seq11_tiny_conv_baseline_20260515_2000steps_lr1e-3_h64/`
- 1 bp DDP pilot:
  `runs/rna_seq11_adapter_1bp_ddp_pilot_20260515_500steps_lr1e-4_warmup50_accum4_clip1/`
- LoRA scheme C pilots:
  `runs/rna_seq11_adapter_128bp_lora_phaseA_20260516_100steps_lr5e-5_freezehead/`
  and `runs/rna_seq11_adapter_128bp_lora_phaseB_20260516_100steps_lr5e-5_lora_head/`

