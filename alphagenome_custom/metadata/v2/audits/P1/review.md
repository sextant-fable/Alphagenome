# R1 Review

Status: **PASS**

Reviewed at: `2026-07-13T16:34:58+00:00`

| Check | Status | Evidence |
| --- | --- | --- |
| R1.01_required_outputs | PASS | missing=[] |
| R1.02_sample_scope | PASS | rows=485 unique_runs=485 |
| R1.03_batch_counts | PASS | counts={'sharepoint_2026-07-13_new439': 439, '2026.4.20': 20, '2026.5.18': 26} |
| R1.04_local_integrity | PASS | all local paths match prior SHA-256 and WBcel235 evidence; P1 before/after hashes agree |
| R1.05_accession_provenance | PASS | ena_success=485 failures=0 |
| R1.06_query_failures_explicit | PASS | statuses=['success'] |
| R1.07_normalization_classification | PASS | classes={'A': 0, 'B': 0, 'C': 485, 'D': 0} |
| R1.08_no_unjustified_formal_inclusion | PASS | formal_included=0 |
| R1.09_evidence_coverage | PASS | evidence_rows=4265 samples=485 |
| R1.10_bigwig_provenance_explicit | PASS | unknown provenance is explicit rather than inferred from bigWig statistics |
| R1.11_conflicts_preserved | PASS | conflict_flagged=61 |
| R1.12_g1_packet | PASS | fastq_bytes=1020357970574 approval_not_granted |
| R1.13_full_rna_source_manifest | PASS | 482 RNA-seq runs have per-file URL, MD5, byte count, and layout |
