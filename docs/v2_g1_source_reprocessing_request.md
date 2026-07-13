# G1 Source-read Reprocessing Approval Packet

Status: **APPROVAL REQUIRED - NOT EXECUTED**

Generated: `2026-07-13T14:07:48+00:00`

## Classification

| Class | Samples |
| --- | ---: |
| A: confirmed per-base RPM | 0 |
| B: known convertible unit | 0 |
| C: unknown bigWig unit with source reads | 485 |
| D: unknown unit without located source reads | 0 |

## Located Source Volume

- ENA run metadata found: `485` / `485`.
- Compressed FASTQ bytes reported by ENA: `1020357970574` (950.28 GiB).
- Total submitted sequence bases: `1950829520132`.
- Planning-only working-space allowance at 4x compressed FASTQ: approximately `3801.13 GiB`.

The 4x allowance is a storage planning estimate, not a measured pipeline requirement. CPU/GPU time is intentionally not fabricated: a representative single-end/paired-end pilot must measure alignment and bigWig-generation throughput before a bulk compute estimate is approved.

## Proposed Target Locations

- Source FASTQ: `shared/source_reads/v2/fastq/`
- Alignment outputs: `alphagenome_custom/tracks/rna_seq_v2_alignment/`
- Normalized bigWigs: `alphagenome_custom/tracks/rna_seq_v2_normalized/`

All locations are ignored by Git. No download or realignment command has been run. The exact downloader, aligner, strandedness handling, multimapping policy, and per-base coverage command will be locked after the metadata/group review and a small benchmark design.
