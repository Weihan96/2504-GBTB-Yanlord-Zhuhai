#!/usr/bin/env python3
"""Compile developer HVAC-panel, video-intercom, and doorbell references for E-302/E-304."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

from e304_router_cad_evidence import (
    EXPECTED_PLAN_DXF_SHA256,
    entity_by_handle,
    first,
    paper_to_model,
    parse_dxf,
    point,
    source_to_ifc,
    text_value,
)


AC_TEXT_HANDLES = ("224298", "2242A9", "22435A", "224364", "25F3E1")
AC_TEXT = "空调开关(底边离地1300mm)"
AC_CANDIDATE_IDS = {"DEV-S002", "DEV-S004", "DEV-S009", "DEV-S010", "DEV-S011"}
VIEWPORT_HANDLE = "224238"
ACCESS_DEFINITIONS = (
    (
        "DEV-INTERCOM-001",
        "E-304",
        "video_intercom",
        "224383",
        "279257",
        "可视对讲机(底边离地1400mm)",
        1400.0,
        [3066.658536, -634.998805],
        "R01",
    ),
    (
        "DEV-DOORBELL-001",
        "E-304",
        "doorbell",
        "266571",
        "266342",
        "门铃(底边离地1300mm)",
        1300.0,
        [6053.880075, -664.974809],
        None,
    ),
)


def parse_args() -> argparse.Namespace:
    root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ifc", type=Path, default=root / "2504 GBTB Yanlord Zhuhai.ifc")
    parser.add_argument(
        "--developer-handover",
        type=Path,
        default=root / "build/mep-positioning/developer-handover-mep-candidate.json",
    )
    parser.add_argument("--plan-dxf", type=Path, default=root / "tmp/dwg/d1-handover-plan.dxf")
    parser.add_argument(
        "--output", type=Path, default=root / "build/elec/elec-developer-control-reference.json"
    )
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> int:
    args = parse_args()
    ifc_hash = sha256(args.ifc)
    if sha256(args.plan_dxf) != EXPECTED_PLAN_DXF_SHA256:
        raise RuntimeError("official handover plan DXF hash changed")
    handover = json.loads(args.developer_handover.read_text(encoding="utf-8"))
    if handover.get("source_ifc_sha256") != ifc_hash:
        raise RuntimeError("developer handover MEP report does not match the current formal IFC")

    entities = parse_dxf(args.plan_dxf)
    viewport = entity_by_handle(entities, VIEWPORT_HANDLE)
    if viewport["type"] != "VIEWPORT" or first(viewport, "69") != "4":
        raise RuntimeError("developer handover plan viewport identity changed")
    ac_text_entities = [entity_by_handle(entities, handle) for handle in AC_TEXT_HANDLES]
    if any(
        entity["type"] != "TEXT"
        or first(entity, "67", "0") != "1"
        or text_value(entity) != AC_TEXT
        for entity in ac_text_entities
    ):
        raise RuntimeError("developer AC-panel height annotations changed")

    hvac_rows = [
        row
        for row in handover.get("records", [])
        if row.get("candidate_role") == "ac_control_panel_reference"
    ]
    if {row.get("candidate_id") for row in hvac_rows} != AC_CANDIDATE_IDS:
        raise RuntimeError("developer handover must contain the five known AC control-panel references")
    hvac_references: list[dict[str, Any]] = []
    for row in sorted(hvac_rows, key=lambda item: item["candidate_id"]):
        if (
            row.get("source_status") != "developer_handover_existing_reference"
            or row.get("automatic_ifc_write_allowed") is not False
        ):
            raise RuntimeError(f"{row['candidate_id']}: developer reference safety status changed")
        hvac_references.append(
            {
                "candidate_id": row["candidate_id"],
                "sheet_id": "E-302",
                "device_role": "hvac_control_panel",
                "position_mm": [*row["ifc_position_mm"], 1300.0],
                "room_reference": row["candidate_space"]["candidate_reference"],
                "room_name": row["candidate_space"]["candidate_long_name"],
                "source_block": row["source_block"],
                "installation_height_mm": 1300.0,
                "height_evidence": {
                    "text": AC_TEXT,
                    "paper_space_text_handles": list(AC_TEXT_HANDLES),
                    "matching_basis": "five model-space AC blocks and five paper-space AC labels; count parity",
                },
                "source_status": "developer_handover_existing_reference_field_confirmation_pending",
                "review_required": True,
                "automatic_ifc_write_allowed": False,
            }
        )

    access_references: list[dict[str, Any]] = []
    for (
        candidate_id,
        sheet_id,
        role,
        text_handle,
        leader_handle,
        expected_text,
        height_mm,
        expected_ifc_xy,
        room_reference,
    ) in ACCESS_DEFINITIONS:
        text_entity = entity_by_handle(entities, text_handle)
        leader = entity_by_handle(entities, leader_handle)
        if (
            text_entity["type"] != "TEXT"
            or first(text_entity, "67", "0") != "1"
            or text_value(text_entity) != expected_text
        ):
            raise RuntimeError(f"{candidate_id}: developer annotation changed")
        if leader["type"] != "LEADER" or first(leader, "67", "0") != "1":
            raise RuntimeError(f"{candidate_id}: developer annotation leader changed")
        if first(text_entity, "330") != first(leader, "330") or first(text_entity, "8") != first(leader, "8"):
            raise RuntimeError(f"{candidate_id}: text and leader no longer share paper-space context")
        ifc_xy = source_to_ifc(paper_to_model(point(leader), viewport))
        if max(abs(actual - expected) for actual, expected in zip(ifc_xy, expected_ifc_xy)) > 0.001:
            raise RuntimeError(f"{candidate_id}: developer leader coordinate changed: {ifc_xy}")
        access_references.append(
            {
                "candidate_id": candidate_id,
                "sheet_id": sheet_id,
                "device_role": role,
                "position_mm": [round(ifc_xy[0], 6), round(ifc_xy[1], 6), height_mm],
                "room_reference": room_reference,
                "installation_height_mm": height_mm,
                "source_evidence": {
                    "text": expected_text,
                    "text_handle": text_handle,
                    "leader_handle": leader_handle,
                    "viewport_handle": VIEWPORT_HANDLE,
                },
                "source_status": "developer_handover_existing_reference_field_confirmation_pending",
                "review_required": True,
                "automatic_ifc_write_allowed": False,
            }
        )

    report = {
        "mode": "read_only_developer_control_reference",
        "source_ifc_sha256": ifc_hash,
        "source_dependencies": [
            {"path": str(args.developer_handover.resolve()), "sha256": sha256(args.developer_handover)},
            {"path": str(args.plan_dxf.resolve()), "sha256": sha256(args.plan_dxf)},
        ],
        "summary": {
            "hvac_control_panels": len(hvac_references),
            "video_intercoms": sum(row["device_role"] == "video_intercom" for row in access_references),
            "doorbells": sum(row["device_role"] == "doorbell" for row in access_references),
        },
        "hvac_control_panel_references": hvac_references,
        "access_control_references": access_references,
        "gates": {
            "five_hvac_control_panel_references_present": len(hvac_references) == 5,
            "five_hvac_height_annotations_present": len(ac_text_entities) == 5,
            "one_video_intercom_reference_present": sum(row["device_role"] == "video_intercom" for row in access_references) == 1,
            "one_doorbell_reference_present": sum(row["device_role"] == "doorbell" for row in access_references) == 1,
            "all_references_pending_field_confirmation": all(
                row["source_status"].endswith("field_confirmation_pending")
                for row in [*hvac_references, *access_references]
            ),
            "automatic_ifc_write_allowed": False,
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"summary": report["summary"], "gates": report["gates"]}, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
