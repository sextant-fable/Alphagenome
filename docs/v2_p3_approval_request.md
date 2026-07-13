# P3 Full Streaming Authorization

Status: **G1/G2/G3 APPROVED FOR ALL 482 RNA-SEQ RUNS**

Authorization recorded: `2026-07-13T16:35:12+00:00`.

## Scientific Scope

- The 485-accession inventory contains 482 RNA-seq runs and 3 ChIP-seq runs.
- The three ChIP-seq accessions (`SRR3535777`, `SRR3535778`, and `SRR3535779`) remain excluded from RNA-seq processing.
- The audited 240-group hierarchy is accepted as the starting hierarchy. Existing unknown-unit bigWigs remain non-formal until uniform source-read processing is complete.
- All 482 RNA-seq source runs are processed, including the three secondary accessions whose provided bigWigs duplicate another file; their uniformly regenerated signals are reviewed before final membership is decided.
- All 482 libraries have unknown strandedness, so the locked comparable signal layer is explicitly unstranded (`.`).

## Source and Resource Scope

- Source manifest: `alphagenome_custom/metadata/v2/p3_full_sources.tsv`.
- Exact compressed FASTQ volume: `1,014,217,532,067` bytes (`944.56 GiB`).
- Layout: 290 single-end and 192 paired-end runs.
- Largest compressed run: `SRR14701368`, `12,637,552,636` bytes (`11.77 GiB`).
- Host capacity at authorization: 112 logical CPUs, approximately 1 TiB RAM, and 2.1 TiB free disk.
- P3B runs four concurrent samples with 16 CPU threads each after P3A measurements pass.
- Per-sample FASTQ, BAM, and bedGraph intermediates are checksum-verified, audited, and removed after atomic bigWig publication. Failed-sample work is retained for diagnosis and the phase pauses.

## Existing Project Environment

The existing `alphagenome` environment is reused. Installation is guarded by:

1. pre-install environment and runtime snapshots;
2. `conda --freeze-installed --dry-run` through a localhost-only DoH proxy;
3. rejection of changes to AlphaGenome, NumPy, PyTorch, CUDA, pyBigWig, or Triton;
4. allowance for a Python build change only when the Python version remains exactly `3.12.13`;
5. post-install imports, CUDA availability check, and repository unit tests.

Locked tools:

- STAR `2.7.11b`;
- samtools `1.23.1` (samtools 1.24 conflicts with STAR's required `htslib <1.24`);
- bedtools `2.31.1`;
- UCSC `bedGraphToBigWig` package `482`.

## Processing Contract

- Reference: WBcel235 FASTA and release-115 GTF already frozen in legacy_v1.
- Alignment: STAR primary unique alignments, mismatch/read-length ratio at most `0.04`.
- BAM filter: exclude flags `2308`, minimum MAPQ `1`.
- Coverage: `bedtools genomecov -split -bg`.
- Scale: `100000000 / sum((end-start) * raw_coverage)`, equivalent to the locked `1e6 x 100 bp` total-signal convention.
- Outputs: `alphagenome_custom/tracks/rna_seq_v2_normalized/` and `alphagenome_custom/tracks/rna_seq_v2_grouped/`.
- Source and generated large files remain ignored by Git; small ledgers, commands, checksums, reviews, and summaries are tracked.

## Later GPU Scope

G4 is approved as `r6_gpu_auto_available_2_3`. At GPU phases the controller must recheck `nvidia-smi`, use GPU 2 by default for one-GPU work and GPU 2/3 for two-GPU work, and never interrupt unrelated processes. P3 remains CPU/I/O work.

G5 is not approved. Chromosome X remains locked until one checkpoint and its SHA-256 are frozen after R6B.
