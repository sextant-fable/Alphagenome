# AlphaGenome Custom Data Processing Report

## Source

- Google Drive folder: `1RVBD8jusGijMsAKRUU9OZJh0bxgS154j`
- Excel metadata file: `remote_inventory/Samples.drive.xlsx`
- File inventory: `remote_inventory/gdrive_files.tsv`
- BigWig manifest: `remote_inventory/bigwig_metadata_manifest.tsv`

## Current Inventory

- Total Google Drive files found: 21
- BigWig files found: 20
- Excel rows matched by bigWig filename stem to `ID`: 20/20
- All matched tracks are `mRNA` in `Sheet1`, mapped to AlphaGenome `RNA_SEQ`.

## Biological Context Groups

AlphaGenome aggregates experiments into biological contexts before training.
For this dataset, grouping by output type, assay, strand, tissue, stage, sex,
condition, and data source yields 11 RNA-seq groups.

Group definitions are in:

- `alphagenome_custom/metadata/track_groups.tsv`

## Track Metadata

Per-file track metadata is in:

- `alphagenome_custom/metadata/track_metadata.tsv`

The current table intentionally leaves `ontology_curie` empty because a stable
C. elegans ontology mapping has not yet been assigned. This should be filled
before building model metadata/protos.

## BigWig QC Status

All 20 input bigWig files have been downloaded and validated through the
`alphagenome_custom/tracks/rna_seq` symlink to `training_input_bigwig`.

QC table:

- `alphagenome_custom/metadata/bigwig_qc.tsv`

Shared chromosome set:

- `I`, `II`, `III`, `IV`, `V`, `X`, `MtDNA`

All files cover `100286401` bases and open as valid bigWig files.

## Signal Scale Summary

Per-sample signal summaries are in:

- `alphagenome_custom/metadata/bigwig_signal_summary.tsv`

Per-biological-context replicate summaries are in:

- `alphagenome_custom/metadata/group_signal_summary.tsv`

Some replicate groups show non-trivial global signal differences, especially
muscle contexts. These differences are retained in the QC tables. The current
processing follows AlphaGenome's RNA-seq approach at the grouping level by
averaging replicate tracks within each biological context on their existing
scale.

## Grouped RNA-seq Tracks

Replicates have been averaged into 11 grouped RNA-seq tracks:

- `alphagenome_custom/tracks/rna_seq_grouped/RNA_SEQ_001.bw`
- `alphagenome_custom/tracks/rna_seq_grouped/RNA_SEQ_002.bw`
- `alphagenome_custom/tracks/rna_seq_grouped/RNA_SEQ_003.bw`
- `alphagenome_custom/tracks/rna_seq_grouped/RNA_SEQ_004.bw`
- `alphagenome_custom/tracks/rna_seq_grouped/RNA_SEQ_005.bw`
- `alphagenome_custom/tracks/rna_seq_grouped/RNA_SEQ_006.bw`
- `alphagenome_custom/tracks/rna_seq_grouped/RNA_SEQ_007.bw`
- `alphagenome_custom/tracks/rna_seq_grouped/RNA_SEQ_008.bw`
- `alphagenome_custom/tracks/rna_seq_grouped/RNA_SEQ_009.bw`
- `alphagenome_custom/tracks/rna_seq_grouped/RNA_SEQ_010.bw`
- `alphagenome_custom/tracks/rna_seq_grouped/RNA_SEQ_011.bw`

Grouped track metadata:

- `alphagenome_custom/metadata/track_metadata_grouped.tsv`

Grouped track QC:

- `alphagenome_custom/metadata/grouped_bigwig_qc.tsv`

## Reference Genome And Intervals

Reference files have been downloaded from Ensembl current release paths for
`Caenorhabditis_elegans.WBcel235`:

- `alphagenome_custom/reference/genome.fa`
- `alphagenome_custom/reference/genome.fa.fai`
- `alphagenome_custom/reference/annotation.gtf`

The FASTA chromosome lengths match the bigWig chromosome lengths exactly:

- `alphagenome_custom/metadata/reference_bigwig_validation.tsv`

Training intervals are 1,048,576 bp windows with 524,288 bp stride. The split
is chromosome-level to avoid leakage across 1 Mb windows:

- train: `I`, `II`, `III`, `IV` (`116` windows)
- valid: `V` (`39` windows)
- test: `X` (`33` windows)
- excluded: `MtDNA`

Interval files:

- `alphagenome_custom/intervals/train.bed`
- `alphagenome_custom/intervals/valid.bed`
- `alphagenome_custom/intervals/test.bed`

Training input manifest:

- `alphagenome_custom/metadata/training_inputs.tsv`

## AlphaGenome-like NPZ Dataset

TensorFlow is not installed in the current local environment, so the first
dataset封装 uses a lightweight per-interval NPZ format rather than TFRecord.
The field names and shapes are aligned with AlphaGenome's RNA_SEQ bundle:

- `dna_sequence`: `[1048576, 4]`, stored as `uint8` one-hot and loaded as `float32`
- `rna_seq`: `[1048576, 11]`, stored as `float32`
- `rna_seq_mask`: `[1, 11]`, bool
- `rna_seq_strand`: `[1, 11]`, int32, currently all `0` for unstranded
- `interval_chromosome`, `interval_start`, `interval_end`

`float32` is used for RNA-seq targets because full-dataset conversion found
signals above the `float16` finite range. The builder now errors if `float16`
would overflow instead of silently writing invalid targets.

Pilot dataset:

- `alphagenome_custom/datasets/rna_seq_npz_pilot_train`

Read-back validation:

- `alphagenome_custom/metadata/npz_pilot_readback.tsv`

Full datasets:

- train: `alphagenome_custom/datasets/rna_seq_npz_train` (`116` examples, `5.4G`)
- valid: `alphagenome_custom/datasets/rna_seq_npz_valid` (`39` examples, `1.8G`)
- test: `alphagenome_custom/datasets/rna_seq_npz_test` (`33` examples, `1.5G`)

Full read-back validation:

- `alphagenome_custom/metadata/npz_train_readback.tsv`
- `alphagenome_custom/metadata/npz_valid_readback.tsv`
- `alphagenome_custom/metadata/npz_test_readback.tsv`

Read-back result:

- batch size: `2`
- `dna_sequence`: `2 x 1048576 x 4`
- `rna_seq`: `2 x 1048576 x 11`
- `rna_seq_mask`: `2 x 1 x 11`
- `rna_seq_strand`: `2 x 1 x 11`
- non-finite `rna_seq` values: `0` in train, valid, and test

Full-split signal ranges after read-back:

- train: min `0`, max `506671`, mean `26.7011`
- valid: min `0`, max `179308`, mean `22.5533`
- test: min `0`, max `27789`, mean `14.6807`

Commands:

```bash
/tmp/ag_gdrive_venv/bin/python scripts/build_rna_seq_npz_dataset.py \
  --split train \
  --max-intervals 4 \
  --output-dir alphagenome_custom/datasets/rna_seq_npz_pilot_train \
  --target-dtype float32 \
  --force

/tmp/ag_gdrive_venv/bin/python scripts/read_rna_seq_npz_dataset.py \
  alphagenome_custom/datasets/rna_seq_npz_pilot_train \
  --batch-size 2 \
  --scan-all \
  --summary-output alphagenome_custom/metadata/npz_pilot_readback.tsv
```

Full dataset conversion commands:

```bash
/tmp/ag_gdrive_venv/bin/python scripts/build_rna_seq_npz_dataset.py \
  --split train \
  --output-dir alphagenome_custom/datasets/rna_seq_npz_train \
  --target-dtype float32 \
  --force

/tmp/ag_gdrive_venv/bin/python scripts/build_rna_seq_npz_dataset.py \
  --split valid \
  --output-dir alphagenome_custom/datasets/rna_seq_npz_valid \
  --target-dtype float32 \
  --force

/tmp/ag_gdrive_venv/bin/python scripts/build_rna_seq_npz_dataset.py \
  --split test \
  --output-dir alphagenome_custom/datasets/rna_seq_npz_test \
  --target-dtype float32 \
  --force

/tmp/ag_gdrive_venv/bin/python scripts/read_rna_seq_npz_dataset.py \
  alphagenome_custom/datasets/rna_seq_npz_train \
  --batch-size 2 \
  --scan-all \
  --summary-output alphagenome_custom/metadata/npz_train_readback.tsv

/tmp/ag_gdrive_venv/bin/python scripts/read_rna_seq_npz_dataset.py \
  alphagenome_custom/datasets/rna_seq_npz_valid \
  --batch-size 2 \
  --scan-all \
  --summary-output alphagenome_custom/metadata/npz_valid_readback.tsv

/tmp/ag_gdrive_venv/bin/python scripts/read_rna_seq_npz_dataset.py \
  alphagenome_custom/datasets/rna_seq_npz_test \
  --batch-size 2 \
  --scan-all \
  --summary-output alphagenome_custom/metadata/npz_test_readback.tsv
```

## PyTorch Training Entry Point

Although the official AlphaGenome research code uses JAX for the model and
TFRecord readers for released training data, these custom C. elegans inputs are
framework-neutral NPZ files. PyTorch training can run on a remote server without
installing TensorFlow locally.

PyTorch helper scripts:

- `scripts/torch_rna_seq_dataset.py`
- `scripts/torch_smoke_train.py`

Remote training notes:

- `alphagenome_custom/metadata/torch_remote_training.md`

The PyTorch dataset returns channel-first tensors:

- `dna_sequence`: `[4, 1048576]`
- `rna_seq`: `[11, 1048576]`
- `rna_seq_mask`: `[11, 1]`
- `rna_seq_strand`: `[11]`

The Google Drive download speed is slow and may require resumable download.
Use:

```bash
/tmp/ag_gdrive_venv/bin/python scripts/download_gdrive_bigwigs.py
```

Then run:

```bash
/tmp/ag_gdrive_venv/bin/python scripts/qc_bigwigs.py
```

To regenerate the grouped tracks, use:

```bash
/tmp/ag_gdrive_venv/bin/python scripts/summarize_bigwig_signals.py
/tmp/ag_gdrive_venv/bin/python scripts/average_bigwig_groups.py
```

## Remaining Requirements

The custom C. elegans RNA_SEQ NPZ inputs are now ready for a dataloader or
training adapter. Remaining decisions before model training/fine-tuning:

- Fill stable C. elegans ontology IDs in `ontology_curie` if downstream model
  metadata/protos require ontology labels.
- Decide whether to keep the current bigWig scale as-is or add a consistent
  normalization/log transform before training. The current files preserve the
  supplied signal scale and average replicates within biological context.
- Choose the model integration path: train a custom reader from these NPZ files,
  or install TensorFlow/protobuf and emit TFRecord/GZIP examples closer to the
  original AlphaGenome training pipeline.
- Add a species/organism adapter if reusing AlphaGenome model code that assumes
  human/mouse organism indices.
