#!/usr/bin/env python3
"""Build an isolated biological-replicate holdout dataset from v2 run tracks.

The formal v2 grouped tracks are immutable.  This module derives a second set of
training and held-out labels from the audited single-run BigWigs and records the
member-level role assignment before any model code is allowed to read them.
"""

from __future__ import annotations

from collections import defaultdict
import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
from typing import Any

from scripts import finalize_v2_reprocessed_groups as finalized


REPO_ROOT = Path(__file__).resolve().parents[1]
METADATA_ROOT = REPO_ROOT / "alphagenome_custom/metadata/v2"
HOLDOUT_METADATA = METADATA_ROOT / "replicate_holdout_v1"
OUTPUT_ROOT = REPO_ROOT / "alphagenome_custom/tracks/rna_seq_v2_replicate_holdout_v1"
SOURCE_GROUPS = METADATA_ROOT / "rna_seq_groups_v2_final.tsv"
SOURCE_MEMBERS = METADATA_ROOT / "rna_seq_group_members_v2_final.tsv"
SOURCE_OUTPUTS = METADATA_ROOT / "p3_group_outputs.tsv"
SOURCE_GROUP_QC = METADATA_ROOT / "replicate_group_qc_v2_post_reprocessing.tsv"
SOURCE_DUPLICATES = METADATA_ROOT / "duplicate_source_reuse_review_v2.tsv"
FINAL_LOCK = METADATA_ROOT / "final_test_lock.json"
FINAL_REPORT = METADATA_ROOT / "final_test_report.json"
SELECTION_SEED = 20260819
EXPECTED_RUNS = 482
EXPECTED_GROUPS = 241
EXPECTED_PRIMARY_CANDIDATES = 58
EXPECTED_SUPPLEMENTARY = 22
EXPECTED_SINGLETONS = 161


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise RuntimeError(f"Refusing to write empty table: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    fields: list[str] = list(rows[0])
    fields.extend(sorted({key for row in rows for key in row}.difference(fields)))
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def atomic_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(8 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def biological_unit(member: dict[str, str]) -> tuple[str, str]:
    if member["biological_replicate"]:
        return "biological_replicate", member["biological_replicate"]
    if member["technical_unit_id"]:
        return "technical_unit_fallback", member["technical_unit_id"]
    raise RuntimeError(f"Member has neither biological nor technical unit: {member}")


def source_equivalence_classes() -> dict[str, str]:
    rows = read_tsv(SOURCE_DUPLICATES)
    classes: dict[str, str] = {}
    for index, row in enumerate(rows, start=1):
        class_id = f"duplicate_source_class_{index:03d}"
        for key in ("canonical_run", "secondary_run"):
            accession = row[key]
            if accession in classes and classes[accession] != class_id:
                raise RuntimeError(f"Run belongs to multiple duplicate classes: {accession}")
            classes[accession] = class_id
    return classes


def verify_source_members(members: list[dict[str, str]]) -> dict[str, Any]:
    if len(members) != EXPECTED_RUNS:
        raise RuntimeError(f"Expected {EXPECTED_RUNS} source members, got {len(members)}")
    accessions = [row["run_accession"] for row in members]
    if len(set(accessions)) != EXPECTED_RUNS:
        raise RuntimeError("Source member manifest contains duplicate run accessions")
    if any(row["membership_status"] != "formal_v2_uniform_reprocessed" for row in members):
        raise RuntimeError("Source member manifest contains non-formal rows")
    checked = 0
    mismatches = []
    for row in members:
        path = REPO_ROOT / row["local_path"]
        if not path.is_file():
            raise FileNotFoundError(path)
        observed = sha256(path)
        checked += 1
        if observed != row["sha256"]:
            mismatches.append({"run_accession": row["run_accession"], "expected": row["sha256"], "observed": observed})
    if mismatches:
        raise RuntimeError(f"Source BigWig SHA-256 mismatch: {mismatches[:3]}")
    return {"member_count": len(members), "sha256_checked": checked, "sha256_mismatches": 0}


def group_units(members: list[dict[str, str]]) -> tuple[dict[str, list[dict[str, str]]], dict[str, dict[str, list[dict[str, str]]]]]:
    by_group: dict[str, list[dict[str, str]]] = defaultdict(list)
    by_unit: dict[str, dict[str, list[dict[str, str]]]] = defaultdict(lambda: defaultdict(list))
    for member in members:
        by_group[member["group_id"]].append(member)
        unit_kind, unit_id = biological_unit(member)
        member = dict(member)
        member["derived_unit_kind"] = unit_kind
        member["derived_unit_id"] = unit_id
        by_unit[member["group_id"]][unit_id].append(member)
    return dict(by_group), {group: dict(units) for group, units in by_unit.items()}


def duplicate_forced_units(
    members_by_group: dict[str, list[dict[str, str]]],
    duplicate_classes: dict[str, str],
    units_by_group: dict[str, dict[str, list[dict[str, str]]]],
) -> dict[str, str]:
    class_members: dict[str, list[tuple[str, str]]] = defaultdict(list)
    group_classes: dict[str, set[str]] = defaultdict(set)
    for group_id, group_members in members_by_group.items():
        for member in group_members:
            class_id = duplicate_classes.get(member["run_accession"])
            if class_id is None:
                continue
            _, unit_id = biological_unit(member)
            class_members[class_id].append((group_id, unit_id))
            group_classes[group_id].add(class_id)
    if not class_members:
        return {}

    graph: dict[str, set[str]] = defaultdict(set)
    for members_in_class in class_members.values():
        groups = {group_id for group_id, _ in members_in_class}
        for group_id in groups:
            graph[group_id].update(groups - {group_id})
    forced: dict[str, str] = {}
    unseen = set(graph)
    while unseen:
        root = min(unseen)
        component = {root}
        frontier = [root]
        while frontier:
            group_id = frontier.pop()
            for neighbour in graph[group_id]:
                if neighbour not in component:
                    component.add(neighbour)
                    frontier.append(neighbour)
        unseen.difference_update(component)
        anchor = min(component)
        candidate_units = [
            unit_id
            for unit_id in units_by_group[anchor]
            if any(
                class_id in group_classes[anchor]
                and any(group_id == anchor and member_unit == unit_id for group_id, member_unit in class_members[class_id])
                for class_id in group_classes[anchor]
            )
        ]
        if not candidate_units:
            raise RuntimeError(f"Duplicate component has no candidate unit in {anchor}")
        _, anchor_unit = min(
            (
                hashlib.sha256(f"{SELECTION_SEED}:{anchor}:{unit_id}".encode()).hexdigest(),
                unit_id,
            )
            for unit_id in candidate_units
        )
        selected_classes = {
            class_id
            for class_id in group_classes[anchor]
            if any(group_id == anchor and member_unit == anchor_unit for group_id, member_unit in class_members[class_id])
        }
        for class_id in selected_classes:
            for group_id, unit_id in class_members[class_id]:
                previous = forced.get(group_id)
                if previous is not None and previous != unit_id:
                    raise RuntimeError(f"Duplicate source classes force multiple units in {group_id}")
                forced[group_id] = unit_id
    for group_id, unit_id in forced.items():
        if len(units_by_group[group_id]) < 2:
            raise RuntimeError(f"Duplicate source class cannot be held out from singleton {group_id}")
        if unit_id not in units_by_group[group_id]:
            raise RuntimeError(f"Forced duplicate unit is missing from {group_id}: {unit_id}")
    return forced


def select_units(
    groups: list[dict[str, str]],
    units_by_group: dict[str, dict[str, list[dict[str, str]]]],
    forced_units: dict[str, str],
    group_qc: dict[str, dict[str, str]],
) -> tuple[list[dict[str, Any]], dict[str, str]]:
    assignments: list[dict[str, Any]] = []
    selected: dict[str, str] = {}
    for group in sorted(groups, key=lambda row: row["group_id"]):
        group_id = group["group_id"]
        units = units_by_group[group_id]
        count = len(units)
        if count == 1:
            role = "training_only"
            selected_unit = ""
            reason = "singleton_no_replicate_holdout"
        else:
            role = "primary_candidate" if count >= 3 else "supplementary_candidate"
            if group_id in forced_units:
                selected_unit = forced_units[group_id]
                reason = "duplicate_source_equivalence_class_forced"
            else:
                ranked = sorted(
                    (
                        hashlib.sha256(
                            f"{SELECTION_SEED}:{group_id}:{unit_id}".encode()
                        ).hexdigest(),
                        unit_id,
                    )
                    for unit_id in units
                )
                _, selected_unit = ranked[0]
                reason = "minimum_sha256_ranked_unit"
            selected[group_id] = selected_unit
        qc = group_qc.get(group_id, {})
        qc_decision = qc.get("group_decision", "missing")
        if role == "primary_candidate" and qc_decision != "candidate_keep":
            role = "primary_pending_qc"
        assignments.append(
            {
                "group_id": group_id,
                "source_n_biological_units": count,
                "assignment_role": role,
                "selected_unit_id": selected_unit,
                "selection_reason": reason,
                "group_qc_decision": qc_decision,
                "group_qc_flagged_pairs": qc.get("flagged_pairs", ""),
                "source_declared_n_biological_units": group["n_biological_units"],
            }
        )
    counts = defaultdict(int)
    for row in assignments:
        counts[row["assignment_role"]] += 1
    expected = {
        "primary_candidate": EXPECTED_PRIMARY_CANDIDATES - 1,
        "primary_pending_qc": 1,
        "supplementary_candidate": EXPECTED_SUPPLEMENTARY,
        "training_only": EXPECTED_SINGLETONS,
    }
    if dict(counts) != expected:
        raise RuntimeError(f"Unexpected assignment counts: observed={dict(counts)} expected={expected}")
    return assignments, selected


def member_role_rows(
    members: list[dict[str, str]],
    assignments: dict[str, dict[str, Any]],
    duplicate_classes: dict[str, str],
) -> list[dict[str, Any]]:
    rows = []
    for member in sorted(members, key=lambda row: (row["group_id"], row["run_accession"])):
        _, unit_id = biological_unit(member)
        assignment = assignments[member["group_id"]]
        heldout = bool(assignment["selected_unit_id"] and unit_id == assignment["selected_unit_id"])
        rows.append(
            {
                "group_id": member["group_id"],
                "run_accession": member["run_accession"],
                "local_path": member["local_path"],
                "source_sha256": member["sha256"],
                "derived_unit_id": unit_id,
                "derived_unit_kind": biological_unit(member)[0],
                "role": "heldout" if heldout else "training",
                "assignment_role": assignment["assignment_role"],
                "duplicate_source_class": duplicate_classes.get(member["run_accession"], ""),
            }
        )
    return rows


def validate_role_boundaries(role_rows: list[dict[str, Any]]) -> None:
    by_run = {row["run_accession"]: row["role"] for row in role_rows}
    if len(by_run) != len(role_rows):
        raise RuntimeError("A source run has more than one role")
    classes: dict[str, set[str]] = defaultdict(set)
    for row in role_rows:
        if row["duplicate_source_class"]:
            classes[row["duplicate_source_class"]].add(row["role"])
    crossed = {key: sorted(value) for key, value in classes.items() if len(value) > 1}
    if crossed:
        raise RuntimeError(f"Duplicate source class crosses training/heldout boundary: {crossed}")
    heldout_units = defaultdict(set)
    for row in role_rows:
        if row["role"] == "heldout":
            heldout_units[row["group_id"]].add(row["derived_unit_id"])
    if any(len(units) != 1 for units in heldout_units.values()):
        raise RuntimeError("A group has zero or multiple heldout biological units")


def manifest_row(group_id: str, audit: dict[str, Any]) -> dict[str, Any]:
    row = dict(audit)
    row["group_id"] = group_id
    return row


def build(args: argparse.Namespace) -> dict[str, Any]:
    groups = read_tsv(SOURCE_GROUPS)
    members = read_tsv(SOURCE_MEMBERS)
    outputs = read_tsv(SOURCE_OUTPUTS)
    group_qc = {row["group_id"]: row for row in read_tsv(SOURCE_GROUP_QC)}
    duplicate_classes = source_equivalence_classes()
    source_audit = verify_source_members(members)
    members_by_group, units_by_group = group_units(members)
    if len(groups) != EXPECTED_GROUPS or len(members_by_group) != EXPECTED_GROUPS:
        raise RuntimeError("Source group/member manifest does not contain 241 groups")
    output_ids = [row["group_id"] for row in outputs]
    group_ids = [row["group_id"] for row in sorted(groups, key=lambda row: row["group_id"])]
    if output_ids != group_ids:
        raise RuntimeError("Source track and group manifests have different group order")
    forced_units = duplicate_forced_units(members_by_group, duplicate_classes, units_by_group)
    assignment_rows, selected_units = select_units(sorted(groups, key=lambda row: row["group_id"]), units_by_group, forced_units, group_qc)
    assignments = {row["group_id"]: row for row in assignment_rows}
    role_rows = member_role_rows(members, assignments, duplicate_classes)
    validate_role_boundaries(role_rows)

    assignment_path = HOLDOUT_METADATA / "holdout_assignments.tsv"
    role_path = HOLDOUT_METADATA / "member_roles.tsv"
    training_manifest = HOLDOUT_METADATA / "training_track_manifest.tsv"
    heldout_manifest = HOLDOUT_METADATA / "heldout_track_manifest.tsv"
    training_groups = HOLDOUT_METADATA / "training_group_manifest.tsv"
    heldout_groups = HOLDOUT_METADATA / "heldout_group_manifest.tsv"
    index_path = HOLDOUT_METADATA / "heldout_track_indices.tsv"
    assignment_by_group = {row["group_id"]: row for row in assignment_rows}
    training_group_rows = []
    heldout_group_rows = []
    for group in sorted(groups, key=lambda row: row["group_id"]):
        updated = dict(group)
        updated["source_n_biological_units"] = len(units_by_group[group["group_id"]])
        selected_unit = selected_units.get(group["group_id"], "")
        updated["training_n_biological_units"] = len(units_by_group[group["group_id"]]) - bool(selected_unit)
        updated["holdout_unit_id"] = selected_unit
        updated["replicate_holdout_role"] = assignment_by_group[group["group_id"]]["assignment_role"]
        training_group_rows.append(updated)
        if selected_unit:
            heldout_group_rows.append(updated)
    index_rows = []
    for index, group in enumerate(sorted(groups, key=lambda row: row["group_id"])):
        assignment = assignment_by_group[group["group_id"]]
        if assignment["selected_unit_id"]:
            index_rows.append(
                {
                    "track_index": index,
                    "group_id": group["group_id"],
                    "assignment_role": assignment["assignment_role"],
                    "evaluation_primary": str(assignment["assignment_role"] == "primary_candidate"),
                }
            )
    write_tsv(assignment_path, assignment_rows)
    write_tsv(role_path, role_rows)
    write_tsv(training_groups, training_group_rows)
    write_tsv(heldout_groups, heldout_group_rows)
    write_tsv(index_path, index_rows)

    spec_path = HOLDOUT_METADATA / "p11_replicate_holdout_spec.json"
    spec = {
        "schema_version": 1,
        "phase": "P11",
        "contract": "replicate_holdout_v1",
        "locked_before_generation": True,
        "selection_seed": SELECTION_SEED,
        "source_paths": {
            "groups": str(SOURCE_GROUPS.relative_to(REPO_ROOT)),
            "members": str(SOURCE_MEMBERS.relative_to(REPO_ROOT)),
            "outputs": str(SOURCE_OUTPUTS.relative_to(REPO_ROOT)),
            "group_qc": str(SOURCE_GROUP_QC.relative_to(REPO_ROOT)),
            "duplicates": str(SOURCE_DUPLICATES.relative_to(REPO_ROOT)),
        },
        "source_sha256": {
            str(path.relative_to(REPO_ROOT)): sha256(path)
            for path in (SOURCE_GROUPS, SOURCE_MEMBERS, SOURCE_OUTPUTS, SOURCE_GROUP_QC, SOURCE_DUPLICATES)
        },
        "assignment_sha256": sha256(assignment_path),
        "member_roles_sha256": sha256(role_path),
        "counts": {
            "source_runs": len(members),
            "source_groups": len(groups),
            "primary_candidate_tracks": sum(row["assignment_role"] == "primary_candidate" for row in assignment_rows),
            "primary_pending_qc_tracks": sum(row["assignment_role"] == "primary_pending_qc" for row in assignment_rows),
            "supplementary_tracks": sum(row["assignment_role"] == "supplementary_candidate" for row in assignment_rows),
            "training_only_tracks": sum(row["assignment_role"] == "training_only" for row in assignment_rows),
        },
        "duplicate_source_boundary": "every duplicate source equivalence class must remain entirely training or entirely heldout",
        "aggregation_policy": "decoded_signal_rescaled_to_1e8_then_raw_coverage_weighted_within_biological_unit_then_equal_mean_across_units",
        "training_labels_only": "selected heldout biological units are excluded from every training aggregate",
        "heldout_labels_only": "heldout aggregates are evaluation-only and excluded from model selection",
        "locked_test_access": "prohibited",
        "p6c_lock_sha256": sha256(FINAL_LOCK),
        "p6c_report_sha256": sha256(FINAL_REPORT),
        "output_root": str(OUTPUT_ROOT.relative_to(REPO_ROOT)),
    }
    if spec_path.exists():
        existing_spec = json.loads(spec_path.read_text())
        if existing_spec != spec:
            raise RuntimeError("Existing P11 spec does not match the current source or assignment contract")
    else:
        atomic_json(spec_path, spec)

    if args.prepare_only:
        return {"status": "prepared", "spec_path": str(spec_path.relative_to(REPO_ROOT)), "spec_sha256": sha256(spec_path), **spec["counts"]}

    training_dir = OUTPUT_ROOT / "training"
    heldout_dir = OUTPUT_ROOT / "heldout"
    training_audits = []
    heldout_audits = []
    expected_chromosomes = finalized.chromosomes()
    for index, group in enumerate(sorted(groups, key=lambda row: row["group_id"]), start=1):
        group_id = group["group_id"]
        selected_unit = selected_units.get(group_id, "")
        training_members = [row for row in members_by_group[group_id] if biological_unit(row)[1] != selected_unit]
        if not training_members:
            raise RuntimeError(f"No training members remain for {group_id}")
        print(f"training_group\t{index}/{len(groups)}\t{group_id}", flush=True)
        training_audits.append(
            manifest_row(
                group_id,
                finalized.aggregate_group(
                    group, training_members, expected_chromosomes, output_dir=training_dir
                ),
            )
        )
        if selected_unit:
            heldout_members = [row for row in members_by_group[group_id] if biological_unit(row)[1] == selected_unit]
            print(f"heldout_group\t{index}/{len(groups)}\t{group_id}", flush=True)
            heldout_audits.append(
                manifest_row(group_id, finalized.aggregate_group(group, heldout_members, expected_chromosomes, output_dir=heldout_dir))
            )

    write_tsv(training_manifest, training_audits)
    write_tsv(heldout_manifest, heldout_audits)

    lock_payload = json.loads(FINAL_LOCK.read_text())
    report_payload = json.loads(FINAL_REPORT.read_text())
    metadata_files = [assignment_path, role_path, training_manifest, heldout_manifest, training_groups, heldout_groups, index_path]
    execution = {
        "schema_version": 1,
        "phase": "P11",
        "contract": "replicate_holdout_v1",
        "status": "complete",
        "created_at": utc_now(),
        "selection_seed": SELECTION_SEED,
        "source_groups": str(SOURCE_GROUPS.relative_to(REPO_ROOT)),
        "source_members": str(SOURCE_MEMBERS.relative_to(REPO_ROOT)),
        "source_outputs": str(SOURCE_OUTPUTS.relative_to(REPO_ROOT)),
        "source_audit": source_audit,
        "source_group_declared_vs_member_derived_mismatches": [
            row["group_id"] for row in assignment_rows
            if row["source_declared_n_biological_units"] != str(row["source_n_biological_units"])
        ],
        "counts": {
            "source_runs": len(members),
            "source_groups": len(groups),
            "primary_candidate_tracks": sum(row["assignment_role"] == "primary_candidate" for row in assignment_rows),
            "primary_pending_qc_tracks": sum(row["assignment_role"] == "primary_pending_qc" for row in assignment_rows),
            "supplementary_tracks": sum(row["assignment_role"] == "supplementary_candidate" for row in assignment_rows),
            "training_only_tracks": sum(row["assignment_role"] == "training_only" for row in assignment_rows),
            "heldout_groups": len(heldout_audits),
        },
        "duplicate_source_classes": duplicate_classes,
        "aggregation_policy": "decoded_signal_rescaled_to_1e8_then_raw_coverage_weighted_within_biological_unit_then_equal_mean_across_units",
        "training_labels_only": "training aggregates exclude the selected heldout biological unit for every group with >=2 derived units",
        "heldout_labels_only": "heldout aggregates are evaluation-only and are never passed to training or model selection",
        "locked_test_access": "prohibited",
        "locked_test_signal_reads": 0,
        "p6c_lock_sha256": sha256(FINAL_LOCK),
        "p6c_report_sha256": sha256(FINAL_REPORT),
        "p6c_lock_test_consumed": lock_payload.get("test_consumed"),
        "p6c_report_schema_version": report_payload.get("schema_version"),
        "metadata_files": {str(path.relative_to(REPO_ROOT)): sha256(path) for path in metadata_files},
        "training_track_manifest_sha256": sha256(training_manifest),
        "heldout_track_manifest_sha256": sha256(heldout_manifest),
        "spec_sha256": sha256(spec_path),
        "training_group_manifest_sha256": sha256(training_groups),
        "heldout_group_manifest_sha256": sha256(heldout_groups),
        "output_roots": {
            "training": str(training_dir.relative_to(REPO_ROOT)),
            "heldout": str(heldout_dir.relative_to(REPO_ROOT)),
        },
    }
    atomic_json(HOLDOUT_METADATA / "p11_replicate_holdout_execution.json", execution)
    return execution


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--preflight-only", action="store_true")
    parser.add_argument("--prepare-only", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.preflight_only:
        members = read_tsv(SOURCE_MEMBERS)
        print(json.dumps(verify_source_members(members), indent=2, sort_keys=True))
        return
    print(json.dumps(build(args), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
