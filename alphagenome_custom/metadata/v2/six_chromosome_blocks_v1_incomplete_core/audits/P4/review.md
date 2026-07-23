# R4 Review

Status: **PASS**

Reviewed at: `2026-07-21T11:51:24+00:00`

| Check | Status | Evidence |
| --- | --- | --- |
| R4.01_required_outputs | PASS | missing=[] |
| R4.02_five_fold_split | PASS | errors=[] chromosomes=['I', 'II', 'III', 'IV', 'V', 'X'] |
| R4.03_block_buffers_and_no_leakage | PASS | block_errors=[] core_errors=[] test_valid=True |
| R4.04_split_lock_integrity | PASS | hash_errors=[] |
| R4.05_loader_shapes_and_dtypes | PASS | shapes={'core_mask': [1048576], 'dna_sequence': [4, 1048576], 'gene_mask': [2, 1048576], 'target_128bp': [241, 8192], 'target_1bp': [241, 1048576], 'track_mask': [241, 1], 'track_strand': [241]} dtypes={'core_mask': 'torch.bool', 'dna_sequence': 'torch.float32', 'gene_mask': 'torch.bool', 'target_128bp': 'torch.float32', 'target_1bp': 'torch.float32', 'track_mask': 'torch.bool', 'track_strand': 'torch.int8'} |
| R4.06_direct_bigwig_and_pooling | PASS | errors=[] |
| R4.07_worker_determinism_and_handles | PASS | worker digests agree and parent file descriptors returned near baseline |
| R4.08_test_embargo | PASS | default_loader_refused=True test_windows=16 chromosomes=['I', 'II', 'III', 'IV', 'V', 'X'] |
| R4.09_io_benchmark | PASS | window_seconds=2.765443637035787 target_bytes=1018724352 |
| R4.10_no_monolithic_npz | PASS | forbidden_paths=[] |
| R4.11_data_manifest_provenance | PASS | benchmark digests bind loader results to splits and P3 tracks |
| R4.12_superseded_holdout_preserved | PASS | archived_intervals=12 |
