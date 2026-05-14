# Collaborator Handoff: Processed C. elegans AlphaGenome-like Data

## Recommended Packages

There are three practical sharing levels. Choose based on what the collaborating
school needs to do.

## A. Review-Only Package

Use this if they only need to review the processing logic, metadata, QC, and data
format.

Files/directories:

```text
Samples.xlsx
remote_inventory/
alphagenome_custom/metadata/
alphagenome_custom/intervals/
scripts/
requirements-torch.txt
```

Approximate size:

```text
< 1 MB, excluding scripts cache files
```

Important files to read first:

```text
alphagenome_custom/metadata/teacher_report_data_format.md
alphagenome_custom/metadata/processing_report.md
alphagenome_custom/metadata/training_inputs.tsv
alphagenome_custom/metadata/track_metadata.tsv
alphagenome_custom/metadata/track_groups.tsv
alphagenome_custom/metadata/track_metadata_grouped.tsv
```

This package is enough for them to understand:

- What the original 20 bigWig files were.
- How Excel metadata was mapped to sample IDs.
- How 20 sample-level tracks were grouped into 11 RNA_SEQ tracks.
- How train/valid/test intervals were generated.
- What final model input shapes are.
- What QC checks passed.

## B. Processed Signal Inspection Package

Use this if they need to inspect the processed tracks in IGV/UCSC or run their
own interval extraction, but do not need the final NPZ training examples yet.

Files/directories:

```text
Samples.xlsx
remote_inventory/
alphagenome_custom/metadata/
alphagenome_custom/intervals/
alphagenome_custom/reference/
alphagenome_custom/tracks/rna_seq_grouped/
scripts/
requirements-torch.txt
```

Approximate size:

```text
~1.2 GB
```

Key processed signal files:

```text
alphagenome_custom/tracks/rna_seq_grouped/RNA_SEQ_001.bw
alphagenome_custom/tracks/rna_seq_grouped/RNA_SEQ_002.bw
alphagenome_custom/tracks/rna_seq_grouped/RNA_SEQ_003.bw
alphagenome_custom/tracks/rna_seq_grouped/RNA_SEQ_004.bw
alphagenome_custom/tracks/rna_seq_grouped/RNA_SEQ_005.bw
alphagenome_custom/tracks/rna_seq_grouped/RNA_SEQ_006.bw
alphagenome_custom/tracks/rna_seq_grouped/RNA_SEQ_007.bw
alphagenome_custom/tracks/rna_seq_grouped/RNA_SEQ_008.bw
alphagenome_custom/tracks/rna_seq_grouped/RNA_SEQ_009.bw
alphagenome_custom/tracks/rna_seq_grouped/RNA_SEQ_010.bw
alphagenome_custom/tracks/rna_seq_grouped/RNA_SEQ_011.bw
```

Track annotation:

```text
alphagenome_custom/metadata/track_metadata_grouped.tsv
```

This is the best package for biological review because grouped bigWig tracks are
standard genome-browser files and are easier to inspect than NPZ training
examples.

## C. Full Training Package

Use this if they need to directly train/fine-tune a model.

Files/directories:

```text
Samples.xlsx
remote_inventory/
alphagenome_custom/metadata/
alphagenome_custom/intervals/
alphagenome_custom/reference/
alphagenome_custom/tracks/rna_seq_grouped/
alphagenome_custom/datasets/rna_seq_npz_train/
alphagenome_custom/datasets/rna_seq_npz_valid/
alphagenome_custom/datasets/rna_seq_npz_test/
alphagenome_custom/datasets/rna_seq_npz_pilot_train/
scripts/
requirements-torch.txt
```

Approximate size:

```text
~10 GB
```

Final model-ready datasets:

```text
alphagenome_custom/datasets/rna_seq_npz_train   # 116 examples, ~5.4G
alphagenome_custom/datasets/rna_seq_npz_valid   # 39 examples, ~1.8G
alphagenome_custom/datasets/rna_seq_npz_test    # 33 examples, ~1.5G
```

Small smoke-test dataset:

```text
alphagenome_custom/datasets/rna_seq_npz_pilot_train  # 4 examples, ~192M
```

Each NPZ training example contains:

```text
dna_sequence: [1048576, 4], uint8 one-hot on disk
rna_seq: [1048576, 11], float32
rna_seq_mask: [1, 11], bool
rna_seq_strand: [1, 11], int32
interval_chromosome
interval_start
interval_end
```

PyTorch dataloader entry point:

```text
scripts/torch_rna_seq_dataset.py
scripts/torch_smoke_train.py
alphagenome_custom/metadata/torch_remote_training.md
```

## Do They Need The Original Raw BigWigs?

Usually no, if the purpose is to review the processed data or train from the
processed dataset.

Send the original raw `training_input_bigwig/` only if they need to audit the
replicate averaging step from scratch. That directory contains the 20 original
sample-level bigWig files. The grouped bigWig tracks already preserve the
processed training signal.

## Suggested First Share

Recommended first upload:

```text
Review-Only Package + Processed Signal Inspection Package
```

That gives them:

- The report and QC tables.
- The grouped 11 RNA_SEQ bigWigs.
- The reference genome and intervals.
- The scripts used to generate the outputs.

Then send the Full Training Package only if they confirm they need to run model
training.

## Suggested Message To Collaborators

```text
I am sharing the processed C. elegans RNA-seq data prepared in an
AlphaGenome-like sequence-to-track format.

The original data consisted of 20 sample-level bigWig files plus an Excel
metadata table. I mapped each bigWig to its metadata by SRR ID, standardized the
metadata into RNA_SEQ tracks, grouped replicates by tissue/stage/source, and
averaged replicates into 11 grouped RNA_SEQ bigWig tracks.

I also prepared 1 Mb genomic intervals on the WBcel235 reference genome and
converted each interval into model-ready examples:

dna_sequence: [1048576, 4]
rna_seq: [1048576, 11]
rna_seq_mask: [1, 11]
rna_seq_strand: [1, 11]

The files under alphagenome_custom/metadata document the processing steps, QC,
track metadata, interval splits, and final training input format.
```
