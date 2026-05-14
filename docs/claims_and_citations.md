# Claims and Citations Policy

This project prioritizes conservative research claims and traceable evidence.

## Prohibited Fabrication

Do not fabricate:

- citations
- DOIs
- paper titles
- author names
- publication years
- datasets
- benchmark scores
- p-values
- model metrics
- experimental results
- biological conclusions

If evidence is missing, say so directly.

## Claim Status Labels

Use one of these labels when discussing project evidence:

- `verified`: checked against local files, command output, or a primary source.
- `unverified`: plausible but not yet checked.
- `inferred`: reasoned from available evidence, but not directly observed.
- `abstract-only`: based only on an abstract.
- `secondary-source-only`: based only on a non-primary source.

Claims about model performance, biological interpretation, or dataset provenance should include a status label unless the evidence is obvious and local.

## Citation Rules

- Prefer primary sources: official documentation, original papers, repositories, dataset pages, or release notes.
- Do not invent missing citation fields.
- If only partial citation information is known, record only the known fields.
- If a source has not been read, do not imply that it has been read.
- For papers, distinguish between full-text reading, abstract-only reading, and secondary-source summaries.
- For software, record repository URL, version, commit, tag, or package version when available.

## Result Reporting Rules

- Do not report a result unless it came from a recorded command, file, or run log.
- Smoke tests only validate that code and environment paths work; they are not training results.
- Failed runs should be reported as failures, not omitted.
- If a metric is computed on pilot data, label it as pilot-only.
- If a metric is computed on validation or test data, record the exact split, command, commit, and input data paths.

## Uncertainty Language

Use direct uncertainty statements:

- "This has not been verified yet."
- "This is inferred from the metadata."
- "This was checked against local file counts."
- "This is based on a secondary source only."
- "No result is available yet."

Avoid stronger language when evidence is incomplete.
