# Part 5: Online Methods Plan

## Purpose

Make the method independently reconstructible. Write this part first once
drafting begins, then use it to constrain the main text. The final manuscript
should have no ambiguity about data provenance, split isolation, baseline
fairness, statistic, or software release.

## Required Topical Subheadings

| Method subsection | Exact content required |
| --- | --- |
| Reporting and software statement | Study type, code repository/commit, license, environment, hardware, LLM-use disclosure and human verification record. |
| RNA-seq source inventory and processing | 485 accession inventory, inclusion/exclusion rule, 482 RNA-seq runs, reference build, alignment/coverage processing and all quality controls. |
| Track grouping and normalization | Biological grouping rule, 241 track manifest, signal unit, within-track normalization, target transform and no monolithic NPZ policy. |
| Genomic windows and leakage controls | Window length, core definition, reverse-complement/shift policy, overlap/receptive-field controls, chromosome roles and access-lock boundary. |
| External benchmark design | Dataset identity, laboratory/study isolation, unseen condition, raw or processed data route, pre-registration date, no-tuning rule and freeze manifest. |
| Model architecture | Frozen trunk, organism embedding, LoRA target modules/rank/alpha, RNA head, output resolutions, parameter counts and initialization. |
| Baseline adaptation and eligibility | Model A, P9 size-matched Model C, Enformer, Borzoi and any excluded public method; exact adaptation recipe, allowed tuning, compute and ineligibility reasons. |
| Training and model selection | Training data, optimizer, steps, loss, augmentation, folds, seed nesting, selection rule and external model freeze. |
| Evaluation metrics and statistics | Gene/exon and 128-bp metrics, primary composite, paired effect definition, fold-first/seed-nested 10,000 bootstrap, strata and all exclusion/missingness rules. |
| P9 component attribution | 2 x 2 worm embedding x LoRA design, one-head test, size-matched control, interaction definition and descriptive treatment of tracks. |
| Variant scoring and interpretation | VCF/reference validation, coordinate convention, ref/alt or haplotype window, RC ensemble, gene score, ISM, null matching and endpoint calculation. |
| Biological validation analysis | External eQTL or perturbation truth source, lead/credible-set policy, hyper-divergent/ambiguous-site handling, direction and ranking statistics. |
| Internal DPY-27 diagnostic | Track identities, contrast definition, seed/fold aggregation, observed/predicted direction, normalized-scale limitation and why it is not external validation. |
| Reproducibility and availability | Checksums, controller/audit artifacts, test suite, clean installation, checkpoint release, data availability and source-data files. |

## Reproducibility Details That Must Not Be Omitted

- Genome annotation and reference versions, plus coordinate system.
- Exact train/validation/test split files and a machine-readable overlap audit.
- Every model's trainable and total parameter counts; hardware, wall time and
  memory for each formal comparator.
- All data transformation and back-transformation steps, including whether
  values can be compared in absolute scale across tracks.
- Random seeds and the distinction between seed repetition and a biological
  replicate.
- A complete data/model/code manifest with SHA-256 hashes and a read-only
  release record.

## Main-Text Discipline

Keep controller phase identifiers, failed operational attempts and exhaustive
hash tables in Online Methods, Extended Data or the release archive. The main
text needs the scientific design, not internal workflow labels.
