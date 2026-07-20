# AlphaGenome for *C. elegans* RNA-seq

This repository adapts AlphaGenome-style sequence-to-track modeling to
*Caenorhabditis elegans* RNA-seq. It contains preprocessing and quality-control
utilities, tracked data manifests, chromosome split definitions, PyTorch model
code, and an auditable phase controller for the current v2 workflow.

Large biological data, reference files, model weights, checkpoints, and run
outputs are deliberately kept outside Git. A clone therefore contains the
code and research record, but not a runnable copy of every dataset or model.

## Project status

The repository contains two distinct research lineages:

- **Legacy v1 is frozen.** It groups 20 provided sample-level bigWigs into 11
  RNA-seq tracks and uses 1,048,576 bp windows. It is retained as exploratory
  evidence and must not be presented as a paper-faithful AlphaGenome
  reproduction. See [the v1 freeze record](docs/legacy_v1_freeze.md).
- **v2 is the active workflow.** It provides provenance and review gates from
  source reprocessing through blocked cross-validation and a separately
  approved chromosome-X evaluation. See the
  [execution plan](docs/v2_execution_plan.md) and
  [append-only execution log](docs/v2_execution_log.md).

As recorded on 2026-07-20, v2 phase P6B was reopened for a fuller validation
matrix. Chromosome X remains behind the separate G5 approval boundary. Check
the live machine-readable state rather than inferring progress from output
files:

```bash
python -m scripts.v2_phase_controller status
```

Smoke-test output validates the environment only. Scientific or model claims
must be tied to recorded inputs, commands, commits, and review artifacts.

## Data contracts

### Legacy v1

- Reference: WBcel235.
- Tracks: 20 provided bigWigs grouped into 11 `RNA_SEQ` contexts.
- Windows: 1,048,576 bp with a 524,288 bp stride.
- Split: chromosomes I-IV train, V validation, X test; MtDNA excluded.
- Dataset format: per-window NPZ with `dna_sequence`, `rna_seq`,
  `rna_seq_mask`, `rna_seq_strand`, and interval coordinates.

The tracked processing details and read-back checks are in
[`alphagenome_custom/metadata/processing_report.md`](alphagenome_custom/metadata/processing_report.md).

### Active v2

- Provenance inventory: 485 accessions, including 482 verified RNA-seq runs
  and three excluded ChIP-seq runs.
- Final data hierarchy: normalized run-level signals aggregated into 241
  biological groups after uniform reprocessing.
- Loader: manifest-driven bigWig access at 1 bp and 128 bp resolution; v2 does
  not build a monolithic NPZ dataset.
- Development splits: five leave-one-chromosome-out folds across chromosomes
  I-V.
- Final test: chromosome X is locked until one selected checkpoint and a
  separate, one-time G5 approval are recorded.
- Model comparison: frozen-trunk head (A), worm embedding plus LoRA (B), and a
  matched from-scratch sequence model (C).

The canonical v2 state is
[`alphagenome_custom/metadata/v2/execution_state.json`](alphagenome_custom/metadata/v2/execution_state.json).
Reviews and supporting evidence live under
`alphagenome_custom/metadata/v2/audits/`.

## Repository layout

```text
Alphagenome/
|-- scripts/                         preprocessing, workflow, training, evaluation
|-- tests/                           v2 unit and regression tests
|-- alphagenome_custom/
|   |-- intervals/                   tracked legacy and v2 split definitions
|   `-- metadata/                    manifests, QC, state, and audit records
|-- remote_inventory/                small source-file inventories
|-- docs/                            plans, logs, results, and reproducibility notes
|-- requirements-torch.txt           minimal legacy PyTorch dependencies
|-- SERVER_README.md                 EEHPC/Slurm reference for the archival host
`-- AGENTS.md                        repository operating and research-integrity rules
```

Common local-only paths include `alphagenome_custom/tracks/`,
`alphagenome_custom/reference/`, `alphagenome_custom/datasets/`, `weights/`,
`checkpoints/`, `logs/`, `runs/`, and `shared/`. They are ignored by Git and
must not be staged as substitutes for manifests or run records.

## Getting started

Run commands from the repository root. The lightweight legacy utilities need
Python 3, NumPy, and PyTorch:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements-torch.txt
```

The full v2 workflow also relies on the prepared project environment and local
assets such as bigWigs, WBcel235 FASTA/GTF files, AlphaGenome weights, and
bioinformatics tools. Do not regenerate or download those assets merely to
make a checkout look complete. Confirm the required inputs for the intended
phase in the execution plan and its audit records first.

Useful read-only checks:

```bash
git status --short --branch
python -m scripts.v2_phase_controller status
python -m py_compile scripts/*.py
```

When the development environment includes pytest, run the v2 regression suite
with:

```bash
python -m pytest -q
```

## Controlled workflow

The v2 controller advances phases only after their implementation and review
both pass:

```text
P0 -> P1 -> P2 -> P3A -> P3B -> P4 -> P5 -> P6A -> P6B -> P6C
```

The phases cover the legacy freeze, provenance, manifest review, pilot and
full reprocessing, dynamic loading, model/loss contracts, GPU smoke tests,
blocked cross-validation, and the one-time final test. Approval gates protect
large downloads or realignment (G1), scientific manifest acceptance (G2),
large generated data (G3), GPU experiments (G4), and chromosome-X access
(G5).

Inspect controller options with:

```bash
python -m scripts.v2_phase_controller --help
```

Do not bypass the controller for a formal run. Every real preprocessing,
training, or evaluation run must be recorded in
[`docs/experiment_log.md`](docs/experiment_log.md) or a linked phase-specific
record. Generated logs and checkpoints remain ignored; tracked summaries store
their provenance and hashes.

## Execution environments

The active checkout is on HY-GPU, a non-Slurm server. Before GPU work, inspect
`nvidia-smi`, use physical GPU 2 by default for a single-GPU job, and use GPUs
2 and 3 for approved two-GPU work. Never interrupt unrelated processes or
oversubscribe GPU memory.

The original EEHPC checkout is archival/source-data only. Its Slurm and storage
reference is retained in [SERVER_README.md](SERVER_README.md); do not apply
those commands to HY-GPU. More detail on which assets belong on each host is in
[the reproducibility notes](docs/reproducibility.md).

## Research and Git policy

- Do not modify raw data in place or commit large/generated biological assets.
- Do not fabricate citations, metrics, biological conclusions, or experimental
  results; follow [the claims policy](docs/claims_and_citations.md).
- Keep failed and superseded runs in the append-only record when they explain
  the research history.
- Stage files explicitly. In this repository, do not use `git add .`.
- Treat chromosome X as previously exposed in the legacy project and still
  protected by the v2 one-time final-test gate.

For repository-specific operational rules, read [AGENTS.md](AGENTS.md) before
changing code, data, workflow state, or experiment records.
