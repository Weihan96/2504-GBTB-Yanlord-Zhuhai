#!/usr/bin/env python3
"""Render the 15-step IFC normalization exit state from current evidence.

This command is read-only. The IFC, normalization standard and decision CSV
remain canonical; JSON reports under build/ are accepted only when their source
SHA-256 matches the current IFC. The output is a derived operator view, not a
second project-management truth source.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import ifcopenshell

from direction_noise_candidate import direction_inventory
from geometry_alignment_audit import sha256


RESOLVED_DECISION_STATUSES = {"confirmed", "implemented", "rejected", "delegated"}
NORMALIZATION_DECISION_PREFIXES = (
    "COORD-",
    "FLOOR-",
    "SLAB-SEM-",
    "COVERING-",
    "FURNITURE-ANCHOR-",
    "SANITARY-ANCHOR-",
)
SOCKET_CONTROLLED_EXCEPTION_IDS = {
    "0laejMoxn8Lu_X3FZaCXmi",
    "2OOjqQDMHDjRcXQCniWXnp",
    "27MTenki57DQsfMryX_1U0",
    "3KXtmVvejA78j_iAS$kydj",
}
PVC110_CONTROLLED_EXCEPTION_IDS = {
    "178mqyyzzFowLcbXcH6prO",
    "0bfVg4Ys1CevZs$qxhkXTo",
}
ORIGIN_HANDOFF_DECISIONS = {
    "COORD-HANDOFF-C003-A104": "A104",
    "COORD-HANDOFF-C003-A105": "A105",
    "COORD-HANDOFF-C003-WFIN": "WFIN",
    "COORD-HANDOFF-C003-PLUM": "PLUM",
    "COORD-HANDOFF-C003-ELEC": "ELEC",
    "COORD-HANDOFF-C003-RCP1": "RCP1",
    "COORD-HANDOFF-C003-INT1": "INT1",
    "COORD-HANDOFF-C003-DET1": "DET1",
}
REQUIRED_CONTROLLED_EXCEPTION_DECISIONS = {
    "COORD-DIRECTION-C003-001": "direction-noise-exception",
    "COORD-WALL-C005-CURVE": "wall-geometry",
    "COORD-WALL-C003-G3": "wall-relationship",
    "COORD-MATL-C003-GB01": "material-layer",
    "FLOOR-TILE-A105-GROUT-001": "flooring-joint",
    "FLOOR-TILE-A105-DRAIN-001": "flooring-joint",
    "FURNITURE-ANCHOR-C003-LOOSE": "furniture-origin-exception",
    "COORD-SOCKET-C003-CONTROLLED": "electric-appliance-origin-exception",
    "COORD-FLOW-C003-PVC110-BUNDLE": "flow-segment-bundle-exception",
}


def parse_standard_steps(text: str) -> list[dict[str, Any]]:
    steps: dict[int, dict[str, Any]] = {}
    pattern = re.compile(
        r"^\|\s*(\d{1,2})\s*\|\s*([^|]+?)\s*\|\s*([^|]+?)\s*\|\s*([^|]+?)\s*\|$"
    )
    for line in text.splitlines():
        match = pattern.match(line)
        if not match:
            continue
        number = int(match.group(1))
        if 1 <= number <= 15:
            steps[number] = {
                "step": number,
                "stage": match.group(2).strip(),
                "operation": match.group(3).strip(),
                "acceptance": match.group(4).strip(),
            }
    if sorted(steps) != list(range(1, 16)):
        raise RuntimeError("normalization standard must define exactly steps 1-15")
    return [steps[number] for number in range(1, 16)]


def read_decisions(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        rows = list(csv.DictReader(handle))
    required = {
        "decision_id",
        "scope",
        "object_guid",
        "current_value",
        "proposed_value",
        "basis",
        "confidence",
        "review_required",
        "status",
        "notes",
    }
    if not rows or set(rows[0]) != required:
        raise RuntimeError(f"invalid decision schema: {path}")
    for row in rows:
        if None in row:
            raise RuntimeError(f"decision row has extra columns: {row['decision_id']}")
        float(row["confidence"])
        if row["review_required"] not in {"yes", "no"}:
            raise RuntimeError(
                f"invalid review_required for {row['decision_id']}: "
                f"{row['review_required']}"
            )
    return rows


def is_normalization_decision(row: dict[str, str]) -> bool:
    return row["decision_id"].startswith(NORMALIZATION_DECISION_PREFIXES)


def parse_object_global_ids(value: str) -> set[str]:
    return {part.strip() for part in value.split(";") if part.strip()}


def parse_ifc_entity_ids(value: str) -> set[int]:
    ids: set[int] = set()
    for part in value.split(";"):
        token = part.strip().removeprefix("#")
        if token:
            ids.add(int(token))
    return ids


def controlled_exception_decision_gate(
    decisions: list[dict[str, str]],
) -> dict[str, Any]:
    matched: list[str] = []
    missing_or_invalid: list[str] = []
    for decision_id, expected_scope in REQUIRED_CONTROLLED_EXCEPTION_DECISIONS.items():
        rows = [row for row in decisions if row["decision_id"] == decision_id]
        if (
            len(rows) == 1
            and rows[0]["scope"] == expected_scope
            and rows[0]["status"] == "implemented"
        ):
            matched.append(decision_id)
        else:
            missing_or_invalid.append(decision_id)
    return {
        "pass": not missing_or_invalid,
        "recorded": sorted(matched),
        "missing_or_invalid": sorted(missing_or_invalid),
    }


def approved_direction_exception_ids(
    decisions: list[dict[str, str]],
) -> set[int]:
    matching = [
        row
        for row in decisions
        if row["decision_id"] == "COORD-DIRECTION-C003-001"
        and row["scope"] == "direction-noise-exception"
        and row["status"] == "implemented"
        and row["review_required"] == "no"
    ]
    if len(matching) != 1:
        return set()
    return parse_ifc_entity_ids(matching[0]["object_guid"])


def approved_furniture_origin_exceptions(
    decisions: list[dict[str, str]],
    furniture_records: list[dict[str, Any]],
) -> set[str]:
    candidate_ids = {
        record["global_id"]
        for record in furniture_records
        if record.get("normalization_disposition")
        == "controlled_exception_keep_complex_origin"
        and record.get("review_required") is False
        and record.get("automatic_write_allowed") is False
    }
    matching = [
        row
        for row in decisions
        if row["decision_id"] == "FURNITURE-ANCHOR-C003-LOOSE"
        and row["status"] == "implemented"
    ]
    if len(matching) != 1:
        return set()
    recorded_ids = parse_object_global_ids(matching[0]["object_guid"])
    return candidate_ids if recorded_ids == candidate_ids else set()


def approved_delegated_opening_origins(
    decisions: list[dict[str, str]],
    opening_records: list[dict[str, Any]],
) -> set[str]:
    matching = [
        row
        for row in decisions
        if row["decision_id"] == "COORD-OPENING-C003-HOSTREL"
        and row["scope"] == "opening-origin"
        and row["status"] == "implemented"
    ]
    if len(matching) != 1:
        return set()
    candidate_ids = {
        record["global_id"]
        for record in opening_records
        if record.get("normalization_disposition") == "delegated_to_host_placement"
        and record.get("review_required") is False
        and record.get("automatic_write_allowed") is False
    }
    recorded_ids = parse_object_global_ids(matching[0]["object_guid"])
    return candidate_ids if recorded_ids == candidate_ids else set()


def recorded_independent_opening_origins(
    decisions: list[dict[str, str]],
    opening_records: list[dict[str, Any]],
) -> set[str]:
    candidate_ids = {
        record["global_id"]
        for record in opening_records
        if record.get("normalization_disposition")
        == "independent_opening_anchor_review"
        and record.get("review_required") is True
        and record.get("automatic_write_allowed") is False
    }
    matching = [
        row
        for row in decisions
        if row["decision_id"] == "COORD-OPENING-C003-INDEPENDENT"
        and row["scope"] == "opening-origin"
        and row["status"] == "pending"
        and row["review_required"] == "no"
    ]
    if len(matching) != 1:
        return set()
    recorded_ids = parse_object_global_ids(matching[0]["object_guid"])
    return candidate_ids if recorded_ids == candidate_ids else set()


def implemented_independent_opening_batch(
    decisions: list[dict[str, str]],
) -> bool:
    """Require the exact implemented shared-Opening batch boundary."""
    expected_ids = {
        "0_eLd9WVfFLvk6I8goe8it",
        "0vxak9xlD8Z97XrPYRGunC",
        "2LCZDONiLBShgK9CpOImZJ",
    }
    matching = [
        row
        for row in decisions
        if row["decision_id"] == "COORD-OPENING-C003-INDEPENDENT"
        and row["scope"] == "opening-origin"
        and row["status"] == "implemented"
        and row["review_required"] == "no"
        and parse_object_global_ids(row["object_guid"]) == expected_ids
    ]
    return len(matching) == 1


def approved_socket_origin_exceptions(
    decisions: list[dict[str, str]],
    remaining_records: list[dict[str, Any]],
) -> set[str]:
    """Accept only the exact implemented socket contact-preservation decision."""
    candidate_ids = {
        record["global_id"]
        for record in remaining_records
        if not record.get("within_review_tolerance", True)
        and record.get("ifc_class") == "IfcElectricAppliance"
        and record.get("review_queue") == "service_installation_anchor_review"
    }
    matching = [
        row
        for row in decisions
        if row["decision_id"] == "COORD-SOCKET-C003-CONTROLLED"
        and row["scope"] == "electric-appliance-origin-exception"
        and row["status"] == "implemented"
        and row["review_required"] == "no"
    ]
    if len(matching) != 1:
        return set()
    recorded_ids = parse_object_global_ids(matching[0]["object_guid"])
    if recorded_ids != SOCKET_CONTROLLED_EXCEPTION_IDS:
        return set()
    return recorded_ids if recorded_ids <= candidate_ids else set()


def approved_pvc110_origin_exceptions(
    decisions: list[dict[str, str]],
    remaining_records: list[dict[str, Any]],
) -> set[str]:
    """Accept the two bundled PVC110 origins only inside the exact exception."""
    candidate_ids = {
        record["global_id"]
        for record in remaining_records
        if not record.get("within_review_tolerance", True)
        and record.get("ifc_class") == "IfcFlowSegment"
        and record.get("review_queue") == "service_installation_anchor_review"
    }
    matching = [
        row
        for row in decisions
        if row["decision_id"] == "COORD-FLOW-C003-PVC110-BUNDLE"
        and row["scope"] == "flow-segment-bundle-exception"
        and row["status"] == "implemented"
        and row["review_required"] == "no"
    ]
    if len(matching) != 1:
        return set()
    recorded_ids = parse_object_global_ids(matching[0]["object_guid"])
    if recorded_ids != PVC110_CONTROLLED_EXCEPTION_IDS:
        return set()
    return recorded_ids if recorded_ids <= candidate_ids else set()


def approved_origin_handoffs(
    decisions: list[dict[str, str]],
    eligible_origin_ids: set[str],
) -> dict[str, set[str]]:
    """Return exact downstream ownership without treating it as an exception.

    Every canonical responsibility row must exist exactly once. A delegated row
    must still cover only current over-tolerance origins; an implemented row
    must cover only objects that have left that queue after downstream QA.
    """
    approved: dict[str, set[str]] = {}
    claimed: set[str] = set()
    for decision_id, owner_task in ORIGIN_HANDOFF_DECISIONS.items():
        matching = [
            row
            for row in decisions
            if row["decision_id"] == decision_id
            and row["scope"] == "c003-origin-responsibility"
            and f"{owner_task} " in f"{row['proposed_value']} "
        ]
        if len(matching) != 1:
            return {}
        row = matching[0]
        recorded_ids = parse_object_global_ids(row["object_guid"])
        if not recorded_ids:
            return {}
        if row["status"] == "delegated" and row["review_required"] == "yes":
            if not recorded_ids <= eligible_origin_ids:
                return {}
        elif row["status"] == "implemented" and row["review_required"] == "no":
            if recorded_ids & eligible_origin_ids:
                return {}
            continue
        else:
            return {}
        if claimed & recorded_ids:
            return {}
        approved[owner_task] = recorded_ids
        claimed |= recorded_ids
    return approved


def load_fresh_report(path: Path, current_sha256: str) -> dict[str, Any]:
    report = json.loads(path.read_text(encoding="utf-8"))
    report_hash = report.get("source", {}).get("sha256")
    if report_hash != current_sha256:
        raise RuntimeError(
            f"stale report {path}: expected {current_sha256}, found {report_hash}"
        )
    return report


def git_baseline_state(source: Path) -> dict[str, Any]:
    repository = Path(
        subprocess.check_output(
            ["git", "-C", str(source.parent), "rev-parse", "--show-toplevel"],
            text=True,
        ).strip()
    )
    relative = source.resolve().relative_to(repository.resolve()).as_posix()
    baseline_exists = (
        subprocess.run(
            ["git", "-C", str(repository), "cat-file", "-e", f"HEAD:{relative}"],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            check=False,
        ).returncode
        == 0
    )
    status = subprocess.check_output(
        ["git", "-C", str(repository), "status", "--short", "--", relative],
        text=True,
    ).strip()
    return {
        "repository": str(repository),
        "relative_ifc_path": relative,
        "head_baseline_exists": baseline_exists,
        "ifc_worktree_status": status,
        "ifc_clean": status == "",
    }


def attach_statuses(
    standard_steps: list[dict[str, Any]],
    evidence: dict[str, Any],
) -> list[dict[str, Any]]:
    statuses = {
        1: ("pass", "15-step standard parsed from the canonical Markdown file"),
        2: (
            "pass" if evidence["git"]["head_baseline_exists"] else "blocked",
            "Git HEAD contains a recoverable IFC baseline",
        ),
        3: (
            "pass" if evidence["all_reports_fresh"] else "blocked",
            "current-hash origin, surface, coordinate and wall reports are present",
        ),
        4: (
            "in_progress" if evidence["remaining_origins"] else "pass",
            f"{evidence['remaining_origins']} product origins still require class-specific anchors",
        ),
        5: (
            "pass" if evidence["remaining_origins_classified"] else "blocked",
            "every remaining origin has a review queue and no generic write authorization",
        ),
        6: (
            "in_progress" if evidence["pending_decisions"] else "pass",
            f"{len(evidence['pending_decisions'])} normalization decisions remain unresolved",
        ),
        7: (
            "in_progress" if evidence["pending_human_decisions"] else "pass",
            f"{len(evidence['pending_human_decisions'])} human decisions remain unresolved",
        ),
        8: (
            "pass" if evidence["implemented_decisions"] else "pending",
            f"{len(evidence['implemented_decisions'])} normalization decisions record completed writes",
        ),
        9: (
            "in_progress" if evidence["remaining_origins"] else "pass",
            "remaining origins are anchored, controlled exceptions, or exact downstream responsibilities",
        ),
        10: (
            "pass" if evidence["current_write_evidence"] else "blocked",
            "current IFC hash is recorded by an implemented normalization decision",
        ),
        11: (
            "pass" if evidence["mechanical_qa_pass"] else "blocked",
            "wall relations and direction residual gates pass at 0.1 mm",
        ),
        12: (
            "in_progress" if evidence["pending_human_decisions"] else "pass",
            "no unresolved C003 design-intent review remains",
        ),
        13: (
            "pass" if evidence["git"]["ifc_clean"] else "pending",
            "IFC batch is committed" if evidence["git"]["ifc_clean"] else "IFC worktree change is not committed",
        ),
        14: (
            "in_progress"
            if evidence["remaining_origins"] or evidence["pending_decisions"]
            else "pass",
            "surface, slab, furniture, services and opening queues remain" if evidence["remaining_origins"] else "all C003 classes are closed or assigned downstream",
        ),
        15: (
            "pass" if evidence["ready"] else "pending",
            "all exit gates pass"
            if evidence["ready"]
            else (
                "technical gates pass; commit gate remains"
                if not evidence["remaining_origins"]
                and not evidence["pending_decisions"]
                and not evidence["git"]["ifc_clean"]
                else "red-origin queues or unresolved C003 decisions remain"
            ),
        ),
    }
    return [
        {**step, "status": statuses[step["step"]][0], "evidence": statuses[step["step"]][1]}
        for step in standard_steps
    ]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--standard", required=True, type=Path)
    parser.add_argument("--decisions", required=True, type=Path)
    parser.add_argument("--coordinate-report", required=True, type=Path)
    parser.add_argument("--remaining-origin-report", required=True, type=Path)
    parser.add_argument("--surface-report", required=True, type=Path)
    parser.add_argument("--furniture-report", required=True, type=Path)
    parser.add_argument("--opening-report", required=True, type=Path)
    parser.add_argument("--wall-report", required=True, type=Path)
    parser.add_argument("--report", required=True, type=Path)
    parser.add_argument(
        "--expected-direction-exception-id", action="append", type=int, default=[]
    )
    parser.add_argument("--direction-ratio-tolerance", type=float, default=1e-6)
    parser.add_argument("--fail-on-not-ready", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    source_path = args.input.resolve()
    current_hash = sha256(source_path)
    model = ifcopenshell.open(source_path)
    standard_steps = parse_standard_steps(
        args.standard.read_text(encoding="utf-8")
    )
    decisions = read_decisions(args.decisions)
    reports = {
        "coordinate": load_fresh_report(args.coordinate_report, current_hash),
        "remaining_origin": load_fresh_report(
            args.remaining_origin_report, current_hash
        ),
        "surface": load_fresh_report(args.surface_report, current_hash),
        "furniture": load_fresh_report(args.furniture_report, current_hash),
        "opening": load_fresh_report(args.opening_report, current_hash),
        "wall": load_fresh_report(args.wall_report, current_hash),
    }
    normalization_decisions = [
        row for row in decisions if is_normalization_decision(row)
    ]
    pending_decisions = [
        row
        for row in normalization_decisions
        if row["status"] not in RESOLVED_DECISION_STATUSES
    ]
    pending_human_decisions = [
        row for row in pending_decisions if row["review_required"] == "yes"
    ]
    implemented_decisions = [
        row for row in normalization_decisions if row["status"] == "implemented"
    ]
    remaining_records = reports["remaining_origin"]["records"]
    direction = direction_inventory(model, args.direction_ratio_tolerance)
    direction_residual_ids = sorted(
        record["direction_id"]
        for record in direction["records"]
        if record["category"] == "near_cardinal_eligible"
    )
    recorded_direction_exception_ids = approved_direction_exception_ids(decisions)
    cli_direction_exception_ids = set(args.expected_direction_exception_id)
    expected_direction_exception_ids = sorted(recorded_direction_exception_ids)
    wall_alignment = reports["wall"]["alignment"]
    wall_relations_pass = (
        wall_alignment["coplanar_edges"]["over_tolerance"] == 0
        and wall_alignment["junctions"]["over_tolerance"] == 0
    )
    direction_cli_consistent = (
        not cli_direction_exception_ids
        or cli_direction_exception_ids == recorded_direction_exception_ids
    )
    direction_gate_pass = (
        direction_residual_ids == expected_direction_exception_ids
        and direction_cli_consistent
    )
    current_write_evidence = any(
        current_hash in row["notes"] for row in implemented_decisions
    )
    remaining_origins_classified = all(
        record.get("review_queue")
        and record.get("automatic_write_allowed") is False
        for record in remaining_records
        if not record["within_review_tolerance"]
    )
    raw_remaining_origin_ids = {
        record["global_id"]
        for record in remaining_records
        if not record["within_review_tolerance"]
    }
    furniture_origin_exception_ids = approved_furniture_origin_exceptions(
        decisions, reports["furniture"]["records"]
    )
    delegated_opening_origin_ids = approved_delegated_opening_origins(
        decisions, reports["opening"]["records"]
    )
    independent_opening_origin_ids = recorded_independent_opening_origins(
        decisions, reports["opening"]["records"]
    )
    socket_origin_exception_ids = approved_socket_origin_exceptions(
        decisions, remaining_records
    )
    pvc110_origin_exception_ids = approved_pvc110_origin_exceptions(
        decisions, remaining_records
    )
    approved_origin_exception_ids = (
        furniture_origin_exception_ids
        | socket_origin_exception_ids
        | pvc110_origin_exception_ids
    )
    pre_handoff_exclusion_ids = (
        approved_origin_exception_ids | delegated_opening_origin_ids
    )
    handoff_eligible_ids = raw_remaining_origin_ids - pre_handoff_exclusion_ids
    origin_handoffs = approved_origin_handoffs(decisions, handoff_eligible_ids)
    origin_handoff_ids = set().union(*origin_handoffs.values()) if origin_handoffs else set()
    origin_handoff_gate_pass = origin_handoff_ids == handoff_eligible_ids
    origin_queue_exclusion_ids = pre_handoff_exclusion_ids | origin_handoff_ids
    remaining_origin_ids = raw_remaining_origin_ids - origin_queue_exclusion_ids
    remaining_origins = len(remaining_origin_ids)
    controlled_exception_decisions = controlled_exception_decision_gate(decisions)
    opening_summary = reports["opening"]["summary"]
    opening_boundary_gate_pass = (
        len(delegated_opening_origin_ids)
        == opening_summary["delegated_to_host_placement"]
        and opening_summary["delegated_to_host_placement"] > 0
        and opening_summary["independent_opening_anchor_review"] == 0
        and not independent_opening_origin_ids
        and implemented_independent_opening_batch(decisions)
    )
    git = git_baseline_state(source_path)
    ready = (
        remaining_origins == 0
        and not pending_decisions
        and wall_relations_pass
        and direction_gate_pass
        and controlled_exception_decisions["pass"]
        and opening_boundary_gate_pass
        and origin_handoff_gate_pass
        and current_write_evidence
        and git["ifc_clean"]
    )
    evidence = {
        "all_reports_fresh": True,
        "git": git,
        "remaining_origins": remaining_origins,
        "raw_remaining_origins": len(raw_remaining_origin_ids),
        "controlled_origin_exception_ids": sorted(approved_origin_exception_ids),
        "controlled_origin_exceptions": len(approved_origin_exception_ids),
        "origin_queue_exclusion_ids": sorted(origin_queue_exclusion_ids),
        "origin_queue_exclusions": len(origin_queue_exclusion_ids),
        "furniture_origin_exception_ids": sorted(furniture_origin_exception_ids),
        "delegated_opening_origin_ids": sorted(delegated_opening_origin_ids),
        "independent_opening_origin_ids": sorted(independent_opening_origin_ids),
        "socket_origin_exception_ids": sorted(socket_origin_exception_ids),
        "pvc110_origin_exception_ids": sorted(pvc110_origin_exception_ids),
        "origin_handoffs": {
            task: sorted(global_ids) for task, global_ids in sorted(origin_handoffs.items())
        },
        "origin_handoff_ids": sorted(origin_handoff_ids),
        "origin_handoff_count": len(origin_handoff_ids),
        "origin_handoff_gate_pass": origin_handoff_gate_pass,
        "remaining_origins_classified": remaining_origins_classified,
        "remaining_origin_queue_counts": reports["remaining_origin"]["summary"][
            "by_review_queue"
        ],
        "pending_decisions": [row["decision_id"] for row in pending_decisions],
        "pending_human_decisions": [
            row["decision_id"] for row in pending_human_decisions
        ],
        "implemented_decisions": [
            row["decision_id"] for row in implemented_decisions
        ],
        "direction_inventory_counts": direction["counts"],
        "direction_residual_ids": direction_residual_ids,
        "expected_direction_exception_ids": expected_direction_exception_ids,
        "direction_gate_pass": direction_gate_pass,
        "direction_cli_consistent": direction_cli_consistent,
        "controlled_exception_decisions": controlled_exception_decisions,
        "opening_boundary_gate_pass": opening_boundary_gate_pass,
        "wall_relations": {
            "coplanar_edges": {
                key: wall_alignment["coplanar_edges"][key]
                for key in ("total", "within_tolerance", "over_tolerance", "max_gap_mm")
            },
            "junctions": {
                key: wall_alignment["junctions"][key]
                for key in ("total", "within_tolerance", "over_tolerance", "max_gap_mm")
            },
        },
        "wall_relations_pass": wall_relations_pass,
        "current_write_evidence": current_write_evidence,
        "mechanical_qa_pass": (
            wall_relations_pass
            and direction_gate_pass
            and controlled_exception_decisions["pass"]
            and opening_boundary_gate_pass
            and origin_handoff_gate_pass
        ),
        "ready": ready,
    }
    steps = attach_statuses(standard_steps, evidence)
    report = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "mode": "read-only-normalization-exit-audit",
        "source": {
            "path": str(source_path),
            "sha256": current_hash,
            "schema": model.schema,
            "entity_count": len(list(model)),
            "root_count": len(model.by_type("IfcRoot")),
        },
        "canonical_inputs": {
            "standard": str(args.standard.resolve()),
            "decisions": str(args.decisions.resolve()),
        },
        "evidence": evidence,
        "steps": steps,
        "summary": {
            "pass": sum(step["status"] == "pass" for step in steps),
            "in_progress": sum(step["status"] == "in_progress" for step in steps),
            "pending": sum(step["status"] == "pending" for step in steps),
            "blocked": sum(step["status"] == "blocked" for step in steps),
            "ready": ready,
        },
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(
        json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(
        json.dumps(
            {
                "report": str(args.report),
                **report["summary"],
                "remaining_origins": remaining_origins,
                "pending_decisions": len(pending_decisions),
                "pending_human_decisions": len(pending_human_decisions),
                "direction_residual_ids": direction_residual_ids,
                "wall_relations_pass": wall_relations_pass,
                "ifc_clean": git["ifc_clean"],
            },
            ensure_ascii=False,
        )
    )
    if args.fail_on_not_ready and not ready:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
