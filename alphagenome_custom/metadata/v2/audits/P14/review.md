# R14 frozen Model B I--V reference Review

Status: **PASS**

Reviewed at: `2026-08-22T18:27:56+00:00`

| Check | Status | Evidence |
| --- | --- | --- |
| R14.01_required_outputs | PASS | missing=[] |
| R14.02_preregistered_i_to_v_contract | PASS | {'allowed_chromosomes': ['I', 'II', 'III', 'IV', 'V'], 'final_test_access': 'prohibited', 'group_manifest': 'alphagenome_custom/metadata/v2/replicate_holdout_v1/heldout_group_manifest.tsv', 'means_path': 'alphagenome_custom/metadata/v2/replicate_holdout_v1/track_nonzero_means.tsv', 'prediction_track_indices': 'alphagenome_custom/metadata/v2/replicate_holdout_v1/heldout_track_indices.tsv', 'primary_metric': '0.5 * gene-exon log1p Pearson + 0.5 * 128-bp log1p Pearson', 'track_manifest': 'alphagenome_custom/metadata/v2/replicate_holdout_v1/heldout_track_manifest.tsv'} |
| R14.03_execution_contract | PASS | completed=15 gpus=[2, 3] |
| R14.04_checkpoint_record_coverage | PASS | records=15 |
| R14.05_fold_and_reference_summary | PASS | folds=5 summary=3 paired=5 |
