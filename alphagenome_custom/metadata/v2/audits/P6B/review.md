# R6B-amended Review

Status: **PASS**

Reviewed at: `2026-07-21T05:10:14+00:00`

| Check | Status | Evidence |
| --- | --- | --- |
| R6B.A1_required_outputs | PASS | missing=[] |
| R6B.A2_preregistered_amendment | PASS | spec_sha256=c0bd7e436f7bc5c8be33e2e7536b36c1ec974338d46a34b7cbb6fe91da60ea8e |
| R6B.A3_complete_three_seed_matrix | PASS | formal=90/90 errors=[] |
| R6B.A4_complete_ablations | PASS | ablations=15/15 configs=['no_augmentation', 'no_gene_loss', 'whole_I_V_mean'] |
| R6B.A5_full_metric_aggregation | PASS | configs=[('A', 'log1p_mse'), ('A', 'paper'), ('B', 'log1p_mse'), ('B', 'paper'), ('C', 'log1p_mse'), ('C', 'paper')] |
| R6B.A6_biological_primary_selection | PASS | selected=B/paper |
| R6B.A7_legacy_disposition | PASS | not_formally_comparable_as_a_label_scale_ablation |
| R6B.A8_development_and_final_lock | PASS | checkpoint=runs/v2_p6b_amendment_20260720/development_selected/checkpoint.pt |
| R6B.A9_old_lock_preserved | PASS | old_checkpoint=runs/v2_p6b_20260714/development_selected/checkpoint.pt |
| R6B.A10_test_embargo | PASS | chromosome_x_reads=0; G5 unapproved |
