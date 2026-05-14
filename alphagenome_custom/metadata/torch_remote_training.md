# PyTorch Remote Training Notes

The local preprocessing outputs are framework-neutral NPZ files. They can be
trained with PyTorch on a remote server without installing TensorFlow locally.

## Data To Transfer

Minimum files needed for PyTorch training:

```text
requirements-torch.txt
scripts/torch_rna_seq_dataset.py
scripts/torch_smoke_train.py
alphagenome_custom/datasets/rna_seq_npz_train
alphagenome_custom/datasets/rna_seq_npz_valid
alphagenome_custom/datasets/rna_seq_npz_test
alphagenome_custom/metadata/track_metadata_grouped.tsv
alphagenome_custom/metadata/training_inputs.tsv
alphagenome_custom/metadata/processing_report.md
```

Optional small smoke-test data:

```text
alphagenome_custom/datasets/rna_seq_npz_pilot_train
```

The full train/valid/test NPZ dataset is about 8.7G. The pilot dataset is about
192M and is useful for checking the server environment before transferring or
training on the full set.

## Environment

Example on a remote server:

```bash
python -m venv ag_torch_venv
source ag_torch_venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements-torch.txt
```

For CUDA servers, install the PyTorch wheel matching the server CUDA version if
the generic install does not select a GPU build.

## Smoke Test

Run this first after copying the project files to the server:

```bash
python scripts/torch_smoke_train.py \
  --dataset-dir alphagenome_custom/datasets/rna_seq_npz_pilot_train \
  --batch-size 1 \
  --max-steps 2 \
  --hidden-channels 16 \
  --device auto
```

Expected tensor shapes:

```text
dna_sequence_shape      1x4x1048576
rna_seq_shape           1x11x1048576
prediction_shape        1x11x1048576
```

Then test the full training split:

```bash
python scripts/torch_smoke_train.py \
  --dataset-dir alphagenome_custom/datasets/rna_seq_npz_train \
  --batch-size 1 \
  --max-steps 2 \
  --hidden-channels 16 \
  --device auto
```

## Dataset Interface

Use `AlphaGenomeRnaSeqNpzDataset` from `scripts/torch_rna_seq_dataset.py`.
It returns channel-first tensors:

```text
dna_sequence: [4, 1048576]
rna_seq: [11, 1048576]
rna_seq_mask: [11, 1]
rna_seq_strand: [11]
```

The default target transform is `log1p`, which is appropriate for initial
training stability because the raw bigWig signal has a high dynamic range. Use
`target_transform="none"` if the model/loss should consume raw signal values.
