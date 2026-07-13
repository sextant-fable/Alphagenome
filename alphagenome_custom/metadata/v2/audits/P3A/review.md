# R3A Review

Status: **PASS**

Reviewed at: `2026-07-13T18:17:11+00:00`

| Check | Status | Evidence |
| --- | --- | --- |
| R3A.01_required_outputs | PASS | missing=[] |
| R3A.02_bounded_scope | PASS | runs=['SRR10882545', 'SRR23049895', 'SRR36719198', 'SRR7443583', 'SRR941632'] formal=False |
| R3A.03_manifest_integrity | PASS | pilot, sample, and locked SRA manifest SHA-256 values match the executed record |
| R3A.04_output_membership | PASS | outputs=['SRR10882545', 'SRR23049895', 'SRR36719198', 'SRR7443583', 'SRR941632'] |
| R3A.05_fastq_integrity | PASS | errors=[] |
| R3A.06_raw_bigwig_immutability | PASS | errors=[] |
| R3A.07_tool_and_command_provenance | PASS | commands=21 versions={'STAR': '2.7.11b', 'bedGraphToBigWig': 'bedGraphToBigWig v 2.10 - Convert a bedGraph file to bigWig format (bbi version: 4).', 'bedtools': 'bedtools v2.31.1', 'fasterq-dump': '/home/zelinli6/Alphagenome/shared/tools/sratoolkit.3.4.1-ubuntu64/bin/fasterq-dump : 3.4.1', 'samtools': 'samtools 1.23.1', 'vdb-validate': '/home/zelinli6/Alphagenome/shared/tools/sratoolkit.3.4.1-ubuntu64/bin/vdb-validate : 3.4.1'} |
| R3A.08_mapping_metrics | PASS | mapping_errors=[] bam_errors=[] |
| R3A.09_normalization_formula | PASS | explicit 1e6 x 100 bp total-signal scale |
| R3A.10_bigwig_integrity | PASS | errors=[] |
| R3A.11_path_and_partial_file_safety | PASS | path_errors=[] partial_files=[] |
| R3A.12_resource_measurement | PASS | elapsed=483.578713 observed_bytes=25792977581 |
| R3A.13_provided_signal_comparison | PASS | comparisons=5 errors=[] |
