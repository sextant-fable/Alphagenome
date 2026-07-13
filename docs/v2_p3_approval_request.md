# P3 Pilot Approval Request

Status: **G1/G2/G3 APPROVAL REQUIRED - NOT EXECUTED**

## G2 Scientific Manifest Review

The corrected R2 hierarchy contains 479 current RNA-seq signal members in 240 candidate groups:

- 82 groups contain multiple source-reported biological units.
- 15 groups contain multiple runs from the same ENA experiment and preserve that technical-run relationship.
- 158 groups are single biological units.
- 3 ChIP-seq accessions are excluded from RNA-seq.
- 3 secondary byte-identical bigWigs are held while their source reuse is preserved in provenance.
- `SRR941651` and `SRR941697` are separate groups and stages.
- 108 of 645 current-signal pairs trigger a review-only QC flag, affecting 17 groups. No sample was automatically removed for this reason.

All 240 groups remain non-formal because their current bigWig unit is unknown. G2 approval accepts this hierarchy for reprocessing and pilot evaluation; it does not promote the current bigWigs to formal labels.

## G1 Bounded Download and Tooling

The pilot manifest is `alphagenome_custom/metadata/v2/p3_pilot_sources.tsv` and contains five runs covering legacy/new, single/paired, cDNA, inverse-rRNA, random-selection, and PolyA libraries.

- Exact compressed FASTQ volume: `6,221,459,671` bytes (`5.79 GiB`).
- FASTQ integrity: ENA-provided per-file MD5.
- Threads: `16` CPU threads.
- No GPU is used.
- No existing raw file is modified.

Required isolated environment command:

```bash
conda create -y -n alphagenome-rnaseq-v2 \
  -c conda-forge -c bioconda \
  python=3.12 star=2.7.11b samtools=1.24 bedtools=2.31.1 \
  ucsc-bedgraphtobigwig pybigwig=0.3.25 numpy
```

The environment is not installed until G1 approval.

## G3 Pilot Data Write

Approved command after the environment is available:

```bash
conda run -n alphagenome-rnaseq-v2 \
  python scripts/run_v2_reprocessing_pilot.py \
  --manifest alphagenome_custom/metadata/v2/p3_pilot_sources.tsv \
  --threads 16 \
  --work-dir shared/source_reads/v2/pilot_20260713 \
  --output-dir alphagenome_custom/tracks/rna_seq_v2_normalized_pilot \
  --star-index shared/reference_indexes/WBcel235_STAR_2.7.11b
```

The pilot uses STAR unique primary alignments, spliced coverage, and an explicit final scaling to total per-base signal `100,000,000`. Output tracks are deliberately unstranded (`.`), so strand behavior is defined even when library strandedness is absent.

Expected temporary use is well below the available `2.1 TiB`, but the bulk 485-run workflow is not approved by this packet. Pilot R3a must measure runtime, disk peaks, mapping rate, output signal, and agreement with the provided bigWigs before a streaming bulk command and cleanup policy are proposed.
