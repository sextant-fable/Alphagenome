# Draft v0 Claim--Evidence Map

This map is the Nature-writing claim contract for the owner-authorized initial
draft. `Verified` means supported by the repository evidence already audited;
`assumed placeholder` means the prose direction is intentionally drafted but
the exact external artifact is still required.

| Manuscript claim | Evidence in draft | Status | Replacement/audit action |
| --- | --- | --- | --- |
| The resource contains 482 RNA-seq runs grouped into 241 biological tracks on WBcel235. | Source inventory, grouped-track manifest and v2 loader. | Verified | Recheck manifest checksum in release record. |
| The task emits 1-bp and 128-bp RNA coverage from 131,072-bp sequence. | Model/data specification and P9/P10 contracts. | Verified | Keep coordinate and target definitions identical in Methods. |
| Full Model B wins a source-isolated benchmark including an unseen condition. | `EXT_BENCHMARK_DATASET_ID`, isolation, condition, leakage and primary-effect tokens. | Assumed placeholder | Replace from signed external benchmark manifest and review. |
| Public comparators are technically eligible under one protocol. | Table 1 eligibility rows and `EXT_PUBLIC_BASELINE_RESULTS`. | Assumed placeholder | Attach adapter, parameter, compute and eligibility audit for each model. |
| LoRA is the dominant measured component effect. | P9 paired/factorial results. | Verified | Preserve five-fold/three-seed unit and 10,000 bootstrap definition. |
| Worm embedding has a smaller positive main effect. | P9 factorial result. | Verified | Do not call it mechanistic or universal. |
| The 1-bp head is necessary. | Not claimed; interval crosses zero. | Explicitly not supported | Keep this boundary in Results and Discussion. |
| Frozen sequence differences recover independent cis-eQTL signal. | `EXT_VARIANT_DATASET_ID` and `EXT_VARIANT_RESULT`. | Assumed placeholder | Replace after variant truth source, nulls and endpoint review. |
| DPY-27 is a successful biological validation. | P10 internal diagnostic has opposite chromosome-level sign. | Rejected | Keep only as Extended Data boundary and Discussion limitation. |
| The release is reproducible and reusable. | `EXT_RELEASE_MANIFEST`, source-data and clean-room clauses. | Assumed placeholder | Populate repository, code, model and environment identifiers. |

## Prose audit rule

The main text is deliberately direct and confident about the registered
scientific argument. Missing external facts are not softened into generic
claims; they remain visible tokens and are tracked here. The draft is not
submission-ready while any token remains.
