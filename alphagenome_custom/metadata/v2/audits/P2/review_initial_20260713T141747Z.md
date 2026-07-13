# R2 Review

Status: **PASS**

Reviewed at: `2026-07-13T14:17:47+00:00`

| Check | Status | Evidence |
| --- | --- | --- |
| R2.01_required_outputs | PASS | missing=[] |
| R2.02_scope_and_decisions | PASS | samples=485 decisions=485 |
| R2.03_wrong_assay_excluded | PASS | chip_ids=['rna_seq:SRR3535777', 'rna_seq:SRR3535778', 'rna_seq:SRR3535779'] |
| R2.04_duplicate_signals | PASS | three secondary byte-identical signals excluded while provenance remains |
| R2.05_single_membership | PASS | members=479 max_memberships=1 |
| R2.06_group_member_counts | PASS | groups=458 |
| R2.07_context_consistency | PASS | all members equal their group across the locked grouping fields |
| R2.08_t38_split | PASS | groups={'SRR941651': 'RNA_V2_G0032', 'SRR941697': 'RNA_V2_G0036'} stages={'SRR941651': 'L1-L3(38H)', 'SRR941697': 'L4(38H)'} |
| R2.09_ontology | PASS | mapping_rows=53 ontology_ids=13 |
| R2.10_replicate_relationship | PASS | multi_sample_groups=13 |
| R2.11_replicate_qc | PASS | pairs=33 expected=33 groups=13 |
| R2.12_no_automatic_qc_exclusion | PASS | flagged_pairs=2 |
| R2.13_conflict_resolutions | PASS | resolution_rows=8 statuses=['resolved', 'resolved_current_signal_raw_reads_require_G1_review'] |
| R2.14_formal_fail_closed | PASS | formal_groups=0 candidate_reviews=458 |
| R2.15_determinism | PASS | group construction repeated in memory with identical stable digest |
