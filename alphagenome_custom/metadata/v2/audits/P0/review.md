# R0 Review

Status: **PASS**

Reviewed at: `2026-07-13T14:40:02+00:00`

| Check | Status | Evidence |
| --- | --- | --- |
| R0.01_required_outputs | PASS | missing=[] |
| R0.02_raw_bigwigs | PASS | count=20 |
| R0.03_grouped_bigwigs | PASS | count=11 |
| R0.04_normalization_and_ontology_recorded | PASS | normalization status and ontology columns are present |
| R0.05_reference_assets | PASS | types=['FASTA', 'FASTA_INDEX', 'GTF'] |
| R0.06_split_counts | PASS | counts={'train': 116, 'valid': 39, 'test': 33, 'all': 188} |
| R0.07_dataset_counts | PASS | counts={'train': 116, 'valid': 39, 'test': 33}; NPZ content not hashed by policy |
| R0.08_candidate_inventory | PASS | count=171 unique=171 |
| R0.09_selected_models | PASS | roles=['base_weights', 'legacy_tested_adapter', 'representation_utility_winner', 'validation_mse_winner'] |
| R0.10_test_exposure | PASS | legacy test exposure and latest validation-only status are explicit |
| R0.11_raw_immutability | PASS | raw hashes agree before/after and protected paths were unchanged |
| R0.12_source_evidence | PASS | source_files=10 |
| R0.13_freeze_document | PASS | docs/legacy_v1_freeze.md |
