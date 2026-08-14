# Official Template Provenance

## Decision

*Nature Methods* does not supply a separate journal-specific initial-submission
template. Its official author guidance accepts Word or TeX/LaTeX and directs
LaTeX authors to the Springer Nature journal article template. The local
template is therefore a content-first authoring skeleton, not a claim that
Nature Methods requires a particular visual layout at first submission.

## Download Record

| Field | Value |
| --- | --- |
| Source organization | Springer Nature |
| Official landing page | https://www.springernature.com/gp/authors/campaigns/latex-author-support |
| Official download URL | https://cms-resources.apps.public.k8s.springernature.io/springer-cms/rest/v1/content/18782940/data/v12 |
| Package stated version | December 2024 journal article template package (v3.1) |
| Downloaded UTC date | 2026-08-14 |
| Local artifact | `template/vendor/springer_nature_latex_template_2024-12.zip` |
| SHA-256 | `812e76dcaa9c28dc1bff1fb6065d51729b67d4ea140552a05088317414a3ecae` |
| Integrity check | `unzip -t` passed for all 18 entries |
| Supporting Information | Not applicable; this is a publisher template, not an article download. |

The original ZIP was expanded without modification under `template/vendor/`.
That directory is ignored by Git because the package is publisher-provided
third-party material. This repository tracks only the provenance record.

## Contents Verified

- `sn-article-template/sn-article.tex`
- `sn-article-template/sn-jnl.cls`
- `sn-article-template/bst/sn-nature.bst`
- `sn-article-template/user-manual.pdf`

## Nature Methods Constraints to Apply Above the Template

- Article abstract: at most 150 words.
- Main text: target 3,000 words, with editorial discretion up to 5,000.
- Main display items: at most six figures and/or tables.
- Main architecture: unheaded Introduction, Results, Discussion and Online
  Methods.
- Use numbered references and check current journal-level instructions before
  submission.
- At final formatting, use editable figure text, at least 300 dpi rasters where
  applicable, no wider than 180 mm, and 5-7 pt sans-serif standard labels.

Official sources checked 2026-08-14:

- https://www.nature.com/nmeth/content
- https://www.nature.com/nmeth/submission-guidelines/aip-and-formatting
- https://www.springernature.com/gp/authors/campaigns/latex-author-support
