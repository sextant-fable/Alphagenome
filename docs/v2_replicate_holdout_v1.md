# V2 Replicate-Holdout Validation

P11 derives an isolated label layer from the audited 482-member final v2
manifest. Biological units are taken from `biological_replicate`, with
`technical_unit_id` as the explicit fallback. The assignment is frozen by
`SELECTION_SEED=20260819` before any derived BigWig is generated.

Groups with three or more member-level units contribute one held-out unit to
the primary candidate pool. Groups with two units contribute supplementary
held-out labels. Singleton groups remain training-only. Exact duplicate-source
equivalence classes are constrained to one side of the boundary. Training
labels exclude the selected unit; held-out labels are evaluation-only.

P12 trains all 165 fresh jobs against the training-label manifest:

- 90 formal A/B/C jobs (two losses, five folds, three seeds);
- 15 one-seed P6B augmentation/mean ablations; and
- 60 three-seed P9 component jobs.

Each checkpoint is evaluated on its fold validation blocks against the
reaggregated training labels and on the same blocks against primary and
supplementary held-out biological-unit labels. The held-out evaluator uses the
full 241-head model but slices only the registered output heads corresponding
to the held-out manifest. Fold 0 and `test_locked.tsv` are prohibited.

The principal inferential unit is the genomic fold, with seeds nested within
fold. Track-level values are retained for descriptive distributions and are
not treated as independent biological replicates. Paired effects use a
fold-first, seed-within-fold hierarchical percentile bootstrap with 10,000
resamples.

## Execution

The controller sequence is P10 -> P11 -> P12. P11 requires a new scoped G3
approval (`p11_replicate_holdout_data`); P12 requires a new scoped G4 approval
(`p12_replicate_holdout_matrix`). The P6C lock/report hashes are checked before
both phases and no phase may read the consumed final-test block.

CPU-only preparation:

```bash
python -m scripts.v2_phase_controller start-p11-replicate-holdout
python -m scripts.v2_phase_controller approve --gate G3 --scope p11_replicate_holdout_data --note "Approved replicate-holdout label generation"
python -m scripts.v2_phase_controller run --phase P11
```

After R11 passes, the fresh GPU matrix is explicitly approved and run:

```bash
python -m scripts.v2_phase_controller approve --gate G4 --scope p12_replicate_holdout_matrix --note "Approved fresh 165-job replicate-holdout matrix on physical GPUs 2 and 3"
python -m scripts.v2_phase_controller run --phase P12
```

The large BigWigs, checkpoints, logs, and generated result bundles remain
ignored local artifacts. Their paths, hashes, manifests, and review reports
are recorded under `alphagenome_custom/metadata/v2/replicate_holdout_v1/`.
