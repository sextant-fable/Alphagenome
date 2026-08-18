# Supplementary Information Download Provenance

## Scope

Six open-access comparison articles and their public publisher-hosted
supplementary material were retrieved on 18 August 2026. The article PDFs are
held in the ignored local reading cache. This record is tracked so that the
retrieval can be reproduced without committing publisher files.

The initial batch downloader retrieved all six article PDFs but its optional SI
browser route could not reach a local CDP proxy. The listed SI files were then
retrieved directly from the official Springer Nature media links exposed by each
article page. Every PDF was accepted by `file` and `pdftotext`; XLSX files
passed `unzip -t`; the CSV was parsed as text.

## Article Corpus

| Article | DOI | Official article page |
| --- | --- | --- |
| Borzoi | 10.1038/s41588-024-02053-6 | https://www.nature.com/articles/s41588-024-02053-6 |
| AlphaGenome | 10.1038/s41586-025-10014-0 | https://www.nature.com/articles/s41586-025-10014-0 |
| Nucleotide Transformer | 10.1038/s41592-024-02523-z | https://www.nature.com/articles/s41592-024-02523-z |
| P-SAGE-net | 10.1038/s41592-026-03124-8 | https://www.nature.com/articles/s41592-026-03124-8 |
| scooby | 10.1038/s41592-025-02854-5 | https://www.nature.com/articles/s41592-025-02854-5 |
| CREsted | 10.1038/s41592-026-03057-2 | https://www.nature.com/articles/s41592-026-03057-2 |

## Supporting Files

| Article | Local cache name | Material | SHA-256 |
| --- | --- | --- | --- |
| AlphaGenome | `alphagenome_supplementary_information.pdf` | Supplementary Information | `86b2e6a07543e3c201e2e157a3e9c6eb5235ff13eb8a5224f6fcb6fe1808d4c0` |
| AlphaGenome | `alphagenome_reporting_summary.pdf` | Reporting Summary | `02d80c781af2b9c16ff6abcd9ad660aa967113302d91b3c7a731f74a829e6b5c` |
| AlphaGenome | `alphagenome_supplementary_tables.xlsx` | Supplementary Tables | `833cb78b6ae6fe39415cfff296ac00c48d800326139f13eb307531a1cc133154` |
| Borzoi | `borzoi_supplementary_information.pdf` | Supplementary Information | `70055692d7f0694dd9e1fd303acd45267bb90a6bb86f64be5a9126700068cb8b` |
| Borzoi | `borzoi_reporting_summary.pdf` | Reporting Summary | `32441847788dd6766b0c0688e4ef007a922ab5e20c3b134d2f7878ea37ee4cf5` |
| CREsted | `crested_supplementary_information.pdf` | Supplementary Information | `e859581fa70382e3b5cb114b4d8612a24a829b1bdd1e298ad9869b960bdc5449` |
| CREsted | `crested_reporting_summary.pdf` | Reporting Summary | `bf7e8cbe6f7f3ccf51ffb55fc2660198938bee15cd221026b70edd102eaa9b09` |
| CREsted | `crested_supplementary_tables.xlsx` | Supplementary Tables | `2c4b73eecd92fd87068b7d4a184004c42100bb6a1ab49bc9b64fc6f43ce41382` |
| Nucleotide Transformer | `nucleotide_transformer_supplementary_information.pdf` | Supplementary Information | `ef0f357a1289aad1910059d455579e242a9a82ab3e32ec38d8ea2584effc436e` |
| Nucleotide Transformer | `nucleotide_transformer_reporting_summary.pdf` | Reporting Summary | `44d38b23ac3ab60a168acafeab244532bc54cea29caeeda44fce39b72db98849` |
| P-SAGE-net | `sage_net_supplementary_information.pdf` | Supplementary Information | `8eb95025ba67da494882ac9dfd901212e7247c3101bb19f0defd6ab933382adc` |
| P-SAGE-net | `sage_net_reporting_summary.pdf` | Reporting Summary | `c72c35513811739fbb7c012da18ce3a0db372f0e3a65dcaab2f9f5ad2cd4cece` |
| scooby | `scooby_supplementary_information.pdf` | Supplementary Information | `858577c5eba939d4f817c22a89d74e36cbb1e2bfd69b0c7ad3e925ccb95d45b7` |
| scooby | `scooby_reporting_summary.pdf` | Reporting Summary | `73f0a6bf5111d93375f47681dc508f6e9e008b101117f9f209c7a3ec1195d910` |
| scooby | `scooby_supplementary_table_1.csv` | Supplementary Table 1 | `45fb95f0d0a3fa0521d1b201d763e770866614737318caa6446de0ec09374690` |

## How the SI Was Used

SI was used to enumerate supplementary controls, breakdowns, comparison
conditions, source-data practices and reporting detail. The manuscript's
argument, paragraph form and main-figure sequence are derived from the article
main texts, as recorded in `main_text_and_si_comparative_reading_20260818.md`.
No publisher figure, wording or numeric result is reused in this project.
