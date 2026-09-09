# R15 I--V-only training interval contract Review

Status: **PASS**

Reviewed at: `2026-08-22T18:28:06+00:00`

| Check | Status | Evidence |
| --- | --- | --- |
| R15.01_required_outputs | PASS | missing=[] |
| R15.02_execution_contract | PASS | metadata-only interval filtering; no GPU, DNA, BigWig, checkpoint, training, chromosome-X or locked-test access |
| R15.03_contract_structure | PASS | folds=5 |
| R15.04_no_x_or_locked_rows | PASS | rows=600 |
| R15.05_source_x_exclusion | PASS | source manifests contain X; output contract excludes it |
