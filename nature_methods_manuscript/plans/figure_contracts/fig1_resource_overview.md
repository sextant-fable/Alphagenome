# Figure Contract: Resource, Task and Adapted Model Overview

## Claim

The v2 resource is a traceable, multi-track *C. elegans* RNA-seq prediction
task with a defined genomic-block development protocol and a specific
organism-adapted model interface.

This figure establishes the object being evaluated. It does not make a
performance, external-generalization or biological-utility claim.

## Inputs

- `alphagenome_custom/metadata/v2/rna_seq_groups_v2_final.tsv`
- Frozen architectural constants verified in
  `scripts/v2_training_components.py`
- Development-fold definition from the v2 protocol: five genomic folds, three
  nested seeds, chromosomes I-V for development.

## Panels

| Panel | Question | Display | Statistical role |
| --- | --- | --- | --- |
| a | What biological resource enters the task? | 482 runs, 27 studies and 241 grouped tracks with a manifest-derived developmental-stage distribution. | Descriptive inventory. |
| b | What is held out during model development? | 131,072-bp DNA window, target scales, five genomic folds and nested seeds. | Protocol diagram, not an outcome. |
| c | What changed in the model? | Frozen trunk, worm embedding, LoRA and dual-resolution 241-track head. | Architecture schematic, not a mechanism claim. |

## Selection and Exclusions

All formal-v2 groups in the manifest are represented. Stage bars use the five
largest manifest categories and one `Other stages` aggregate defined before
plotting. The figure deliberately contains no locus example: that panel is
reserved for the future P11 pre-registered development-only coverage analysis,
which must use frozen checkpoints and validation manifests.

## Output Contract

- Python/matplotlib backend.
- 183 mm x 112 mm figure with editable PDF/SVG text and 600 dpi PNG/TIFF.
- Source data reports every formal group plus derived summary rows.
- Output is generated under `results/v2_resource_overview_figure/` and is not
  a Git-tracked data artifact.
