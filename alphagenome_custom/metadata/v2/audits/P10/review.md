# R10-dpy27-internal-application Review

Status: **PASS**

Reviewed at: `2026-08-14T01:51:12+00:00`

| Check | Status | Evidence |
| --- | --- | --- |
| R10.01_required_outputs | PASS | all frozen outputs present |
| R10.02_frozen_contract | PASS | input_errors=[] checkpoint_errors=[] |
| R10.03_preflight | PASS | checkpoints=15 tracks=241 |
| R10.04_execution_matrix | PASS | status=completed workers=15 gpus=[2, 3] |
| R10.05_restart_shards | PASS | shards=15/15 errors=[] |
| R10.06_fold_primary_source_data | PASS | runs=15 folds=5 summaries=4 artifact_errors=[] |
| R10.07_figure_qa | PASS | qa=passed artifact_errors=[] |
| R10.08_final_test_preserved | PASS | P6C lock and report hashes unchanged |
| R10.09_controller_scope | PASS | phase=P10 status=REVIEWING |
