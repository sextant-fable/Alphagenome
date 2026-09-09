# R13 I--V replicate-agreement reference Review

Status: **PASS**

Reviewed at: `2026-08-22T18:27:56+00:00`

| Check | Status | Evidence |
| --- | --- | --- |
| R13.01_required_outputs | PASS | missing=[] |
| R13.02_execution_contract | PASS | CPU-only BigWig reads; no GPU, checkpoint, training, chromosome-X or locked-test access |
| R13.03_i_to_v_boundary | PASS | folds=5 |
| R13.04_group_coverage | PASS | per_group_rows=285 |
| R13.05_summary_coverage | PASS | by_fold=5 summary=3 |
| R13.06_analysis_audit | PASS | ['No checkpoint was loaded.', 'No GPU was used.', 'No training or data generation was performed.', 'No locked final-test intervals or signals were opened.'] |
