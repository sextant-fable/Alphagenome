# Audited v2 variant scoring

## Contract

`scripts/score_v2_variants.py` accepts one frozen JSON spec and writes an audit
alongside tabular results. The spec must declare contract
`v2_variant_scoring_v1`, prohibit final-test access, lock the checkpoint,
checkpoint `run.json`, VCF, WBcel235 FASTA/FAI/GTF, track order, group metadata,
fold means, model specs, and base weights by SHA-256, and identify one Model B
fold/seed/loss configuration.

The VCF parser uses only the Python standard library, supports plain or gzip
VCF, splits multiallelic records, and requires concrete A/C/G/T REF and ALT
alleles. Sequence context may contain `N`, encoded as all-zero one-hot columns.
POS is interpreted as one-based VCF coordinates and REF must exactly match the
FASTA before a model is loaded.

Every variant uses equal-length reference and alternate windows with the
variant start at a fixed anchor. Insertions crop right context; deletions extend
right context. Equal-length substitutions use coordinate-aligned gene/exon
aggregation. Indels are intentionally limited to track-level output labelled
`fixed_window_index_aligned_diagnostic`; their shifted downstream indices are
not mapped onto reference gene/exon coordinates, and the audit reports
`skipped_indel_gene_aggregation`.

Forward and reverse-complement predictions are averaged after restoring the RC
position axis. The spec must assert
`all_241_tracks_unstranded_no_rc_permutation`; preflight verifies `strand=.` for
every group, so no track-channel permutation is applied. Outputs include
per-track delta sum, mean, L1 magnitude and maximum absolute delta, plus
gene-body and unique-exon summaries for coordinate-aligned substitutions.
Optional ISM enumerates all three substitutions at each requested canonical
base and enforces a frozen mutant cap.

## Minimal spec shape

```json
{
  "schema_version": 1,
  "contract": "v2_variant_scoring_v1",
  "final_test_access": "prohibited",
  "device": "cpu",
  "model": {
    "model_id": "B",
    "loss": "paper",
    "fold": 1,
    "seed": 20260714,
    "sequence_length": 131072,
    "hidden_channels": 64,
    "mean_column": "fold_1_train_nonzero_mean"
  },
  "inference": {
    "forward_reverse_complement_ensemble": true,
    "allele_window_policy": "variant_start_anchored_right_context_crop_or_extend",
    "track_semantics": "all_241_tracks_unstranded_no_rc_permutation",
    "require_pass": true,
    "maximum_alleles": 1000
  },
  "ism": {"enabled": false},
  "locked_inputs": {"...": "repository or absolute paths"},
  "input_sha256": {"...": "matching SHA-256 values"},
  "output_dir": "variant_scoring_output"
}
```

The exact required input keys are `checkpoint`, `checkpoint_run`, `fasta`,
`fai`, `gtf`, `group_manifest`, `means`, `track_manifest`, `model_specs`,
`model_weights`, and `vcf`.

CPU preflight validates hashes, VCF coordinates, REF alleles, track semantics,
and checkpoint metadata without loading the model:

```bash
/home/zelinli6/miniconda3/envs/alphagenome/bin/python \
  -m scripts.score_v2_variants --spec path/to/spec.json --preflight-only
```

CUDA execution additionally requires `physical_gpu` 2 or 3, a nonempty
`execution_id`, and `cuda_authorization_path`. The authorization JSON must bind
`approved: true`, scope `variant_scoring:<execution_id>`, the exact spec hash,
physical GPU, checkpoint hash, and VCF hash. This prevents implementation-only
work from silently authorizing a new GPU analysis.

No external eQTL or variant dataset is downloaded by this tool.
