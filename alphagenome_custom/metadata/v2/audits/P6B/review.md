# R6B-six-chromosome Review

Status: **PASS**

Reviewed at: `2026-07-23T04:33:10+00:00`

| Check | Status | Evidence |
| --- | --- | --- |
| R6B.S1_required_outputs | PASS | missing=[] |
| R6B.S2_preregistered_revision | PASS | spec_sha256=7f50967e1e53427b5b0ac6630c1b19c9b7344a79e67b60cfb36cb0c1a9289a05 |
| R6B.S3_complete_matrix | PASS | formal=90/90 ablations=15/15 errors=[] |
| R6B.S4_gpu_policy | PASS | selected=[2, 3] |
| R6B.S5_metrics_and_aggregation | PASS | configs=[('A', 'log1p_mse'), ('A', 'paper'), ('B', 'log1p_mse'), ('B', 'paper'), ('C', 'log1p_mse'), ('C', 'paper')] aggregation_errors=[] |
| R6B.S6_selection_and_development | PASS | selected=B/paper |
| R6B.S7_final_lock | PASS | checkpoint=runs/v2_p6b_six_chromosome_corefix_20260722/development_selected/checkpoint.pt |
| R6B.S8_test_embargo | PASS | locked_reads=0; G5 unapproved |
| R6B.S9_old_evidence_preserved | PASS | old_checkpoint=runs/v2_p6b_amendment_20260720/development_selected/checkpoint.pt |
| R6B.S10_incomplete_core_attempt_excluded | PASS | source_status=running coverage=[0.8435322914705329] |
