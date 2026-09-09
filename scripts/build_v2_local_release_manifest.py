"""Build a small local release manifest for the completed manuscript evidence.

The manifest records hashes and paths only; it does not copy large checkpoints or
publish data.  It is intentionally a local reproducibility aid until a public
archive and redistribution permissions are available.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "nature_methods_manuscript/reproducibility/LOCAL_RELEASE_MANIFEST_v1.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def artifact(path: str, required: bool = True) -> dict[str, object]:
    candidate = ROOT / path
    return {
        "path": path,
        "exists": candidate.exists(),
        "required": required,
        "sha256": sha256(candidate) if candidate.is_file() else None,
        "bytes": candidate.stat().st_size if candidate.is_file() else None,
    }


def main() -> None:
    paths = [
        "alphagenome_custom/metadata/v2/p17_iv_controlled_matrix_execution.json",
        "alphagenome_custom/metadata/v2/audits/P17/review.json",
        "alphagenome_custom/metadata/v2/p19_external_represented_head_scoring_execution.json",
        "alphagenome_custom/metadata/v2/audits/P19/review.json",
        "results/v2_p20_eqtl_variant_set/p20_eqtl_variant_set_manifest.json",
        "alphagenome_custom/metadata/v2/audits/P20/review.json",
        "alphagenome_custom/metadata/v2/p21_eqtl_variant_scoring_execution.json",
        "alphagenome_custom/metadata/v2/audits/P21/review.json",
        "results/v2_p17_manuscript_figure_20260823T032000Z_r2/figure_p17_qa.json",
        "results/v2_p19_manuscript_figure_20260823T060118Z/figure_p19_qa.json",
        "results/v2_p22_manuscript_figure_20260823T0710Z/p22_eqtl_analysis_manifest.json",
        "results/v2_p22_manuscript_figure_20260823T0710Z/p22_eqtl_figure_source_data.tsv",
        "results/v2_p22_manuscript_figure_20260823T0710Z/figure_p22_eqtl_effects.pdf",
        "results/v2_p22_manuscript_figure_20260823T0710Z/figure_p22_eqtl_effects.svg",
        "results/v2_p22_manuscript_figure_20260823T0710Z/figure_p22_eqtl_effects.tiff",
        "results/v2_p22_manuscript_figure_20260823T0710Z/figure_p22_qa.json",
        "results/v2_ed6_calibration_strata/figure_ed6_manifest.json",
        "results/v2_ed6_calibration_strata/figure_ed6_qa.json",
        "results/v2_ed6_calibration_strata/ed6_group_source_data.tsv",
        "results/v2_ed6_calibration_strata/ed6_strata_source_data.tsv",
        "nature_methods_manuscript/extended_data/ED7_profiles/ed7_source_data.tsv",
        "nature_methods_manuscript/extended_data/ED8_dpy27/ed8_source_data.tsv",
        "nature_methods_manuscript/manuscript/alphagenome_nature_methods_draft_v1.tex",
        "nature_methods_manuscript/manuscript/alphagenome_nature_methods_draft_v1.pdf",
        "nature_methods_manuscript/reproducibility/CLEAN_ROOM_REPRODUCTION_v1.md",
        "nature_methods_manuscript/reproducibility/EXTENDED_DATA_INDEX_v1.md",
        "nature_methods_manuscript/reproducibility/SUPPLEMENTARY_TABLE_INDEX_v1.md",
        "nature_methods_manuscript/reproducibility/PUBLIC_RELEASE_STAGING_v1.md",
        "alphagenome_custom/metadata/v2/p23_lora_sensitivity_spec.json",
        "alphagenome_custom/metadata/v2/p24_borzoi_adapter_spec.json",
        "alphagenome_custom/metadata/v2/p25_borzoi_adapter_pilot_spec.json",
        "alphagenome_custom/metadata/v2/p23_lora_sensitivity_execution.json",
        "alphagenome_custom/metadata/v2/p24_borzoi_adapter_execution.json",
        "alphagenome_custom/metadata/v2/p25_borzoi_adapter_pilot_execution.json",
    ]
    checkpoints = sorted(
        ROOT.glob("runs/v2_p17_iv_controlled_matrix/B_iv_dual/seed_*/fold_*/checkpoint.pt")
    )
    payload = {
        "schema_version": 1,
        "status": "local_manifest_only",
        "scope": "P17/P19/P20/P21/P22 manuscript evidence; no locked-test or chromosome-X outputs",
        "artifacts": [artifact(path) for path in paths],
        "p17_b_iv_dual_checkpoints": [
            artifact(str(path.relative_to(ROOT))) for path in checkpoints
        ],
        "release_blockers": [
            "public immutable archive and DOI not created",
            "fresh-machine clean-room run not completed",
            "checkpoint redistribution permissions not recorded",
        ],
    }
    OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"output": str(OUTPUT.relative_to(ROOT)), "checkpoint_count": len(checkpoints)}, indent=2))


if __name__ == "__main__":
    main()
