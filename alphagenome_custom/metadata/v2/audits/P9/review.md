# R9-submission-evidence-matrix Review

Status: **PASS**

Reviewed at: `2026-08-14T01:17:41+00:00`

| Check | Status | Evidence |
| --- | --- | --- |
| R9.01_required_outputs | PASS | missing=[] |
| R9.02_locked_inputs | PASS | input_errors=[] dataset_errors=[] reuse_errors=[] implementation_errors=[] |
| R9.03_matrix_contract | PASS | configs={'frozen_no_worm_no_lora': 15, 'full_worm_lora': 15, 'learned_1bp_head_only': 15, 'lora_only': 15, 'no_lora': 15, 'size_matched_scratch': 15} reused={'P6B': 30, 'P8': 1} |
| R9.04_artifact_hashes | PASS | path_errors=[] validation_errors=[] |
| R9.05_parameter_contract | PASS | errors=[] |
| R9.06_paired_factorial_statistics | PASS | raw_paired=990 raw_factorial=495 summary=33 recompute_errors=[] |
| R9.07_per_track_source_data | PASS | rows=65070 tracks=241 |
| R9.08_final_test_preserved | PASS | P6C lock and report hashes unchanged |
| R9.09_controller_scope | PASS | phase=P9 |
