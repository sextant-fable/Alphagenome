# R5 Review

Status: **PASS**

Reviewed at: `2026-07-22T02:45:18+00:00`

| Check | Status | Evidence |
| --- | --- | --- |
| R5.01_required_outputs | PASS | missing=[] |
| R5.02_leakage_safe_track_means | PASS | tracks=241 test_reads=0 |
| R5.03_model_space_scaling_and_dual_loss | PASS | head_shapes={'1': [1, 241, 1024], '128': [1, 241, 8]} loss=27.337594985961914 |
| R5.04_gene_loss_and_augmentation_lock | PASS | gene_weight=0.1 augmentation={'max_shift_bp': 1024, 'reverse_complement_probability': 0.5} |
| R5.05_three_model_contract | PASS | models={'A': {'base_organism_index': 0, 'description': 'frozen AlphaGenome trunk with new dual RNA head', 'model_id': 'A', 'trunk_policy': 'frozen'}, 'B': {'base_organism_index': 2, 'description': 'new C. elegans embedding plus final-tower LoRA and dual RNA head', 'model_id': 'B', 'trunk_policy': 'worm_embeddings_and_lora_only'}, 'C': {'base_organism_index': None, 'description': 'residual sequence baseline trained from scratch', 'model_id': 'C', 'trunk_policy': 'from_scratch'}} |
| R5.06_real_worm_embedding | PASS | paths=['organism_embed', 'embedder_128bp.organism_embed', 'embedder_1bp.organism_embed', 'embedder_pair.organism_embed'] index=2 |
| R5.07_checkpoint_freeze_and_lora_scope | PASS | lora_modules=6 unexpected=[] |
| R5.08_cpu_forward_backward | PASS | head={'1': [1, 241, 1024], '128': [1, 241, 8]} baseline={'1': [1, 16, 1024], '128': [1, 16, 8]} |
| R5.09_unit_and_golden_tests | PASS | unit_tests=80 return=0 |
| R5.10_test_embargo | PASS | x=registered_train_valid_blocks_allowed test=prohibited_until_locked_G5_P6C |
| R5.11_manifest_hashes | PASS | P5 audit hashes bind means and model specification |
