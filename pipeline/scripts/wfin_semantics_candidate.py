#!/usr/bin/env python3
"""Create the read-only WFIN wall-finish intent IFC candidate."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import uuid
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import ifcopenshell
import ifcopenshell.api
import ifcopenshell.geom
import ifcopenshell.util.element
import numpy as np


GUID_NAMESPACE = uuid.UUID("af122a9c-3690-4c84-8ae5-d5fe14fa43de")
PSET_NAME = "Pset_WallFinishIntent"
WFIN_HANDOFF_IDS = {
    "0WdsBKMo545Ryicomy6Mqg",
    "0moBrr1cf5WgtBzqXWKtqn",
    "2ntxn4aYnB0xQOSuraf1r2",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--register", required=True, type=Path)
    parser.add_argument("--segments", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--issues", required=True, type=Path)
    parser.add_argument("--report", required=True, type=Path)
    parser.add_argument("--tolerance-mm", type=float, default=0.1)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def deterministic_guid(global_id: str, role: str) -> str:
    return ifcopenshell.guid.compress(uuid.uuid5(GUID_NAMESPACE, f"WFIN:{global_id}:{role}").hex)


def read_register(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))
    if len(rows) != 51:
        raise RuntimeError(f"expected 51 WFIN rows, got {len(rows)}")
    counts: dict[str, int] = {}
    for row in rows:
        counts[row["candidate_finish_code"]] = counts.get(row["candidate_finish_code"], 0) + 1
        expected_status = "deferred_material_review" if row["candidate_finish_code"] == "MULTI_FINISH_SPLIT_REQUIRED" else "confirmed_candidate"
        if row["status"] != expected_status or row["formal_ifc_write_allowed"] != "no":
            raise RuntimeError(f"unapproved WFIN row {row['covering_global_id']}")
    if counts != {"WHITE_WALL": 27, "TADELAKT": 23, "MULTI_FINISH_SPLIT_REQUIRED": 1}:
        raise RuntimeError(f"unexpected WFIN partition {counts}")
    return rows


def read_segments(path: Path) -> dict[str, list[dict[str, str]]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))
    counts: dict[str, int] = {}
    by_covering: dict[str, list[dict[str, str]]] = {}
    for row in rows:
        counts[row["candidate_finish_code"]] = counts.get(row["candidate_finish_code"], 0) + 1
        by_covering.setdefault(row["covering_global_id"], []).append(row)
        if row["formal_ifc_write_allowed"] != "no":
            raise RuntimeError(f"segment write boundary missing for {row['segment_id']}")
    if len(rows) != 52 or counts != {"WHITE_WALL": 28, "TADELAKT": 24} or len(by_covering) != 51:
        raise RuntimeError(f"unexpected WFIN segment partition: rows={len(rows)}, counts={counts}, objects={len(by_covering)}")
    return by_covering


def world_vertices(settings: ifcopenshell.geom.settings, product: ifcopenshell.entity_instance) -> np.ndarray:
    shape = ifcopenshell.geom.create_shape(settings, product)
    vertices = np.asarray(shape.geometry.verts, dtype=float).reshape((-1, 3)) * 1000.0
    return vertices[np.lexsort((vertices[:, 2], vertices[:, 1], vertices[:, 0]))]


def material_signature(product: ifcopenshell.entity_instance) -> Any:
    material = ifcopenshell.util.element.get_material(product)
    if material is None:
        return None
    if material.is_a("IfcMaterial"):
        return {"type": material.is_a(), "name": str(material.Name or "")}
    if material.is_a("IfcMaterialLayerSetUsage"):
        layer_set = material.ForLayerSet
        return {
            "type": material.is_a(),
            "layer_set_name": str(layer_set.LayerSetName or ""),
            "direction": str(material.LayerSetDirection),
            "sense": str(material.DirectionSense),
            "offset": float(material.OffsetFromReferenceLine),
            "layers": [
                {
                    "name": str(layer.Name or ""),
                    "thickness": float(layer.LayerThickness),
                    "material": str(layer.Material.Name or "") if layer.Material else "",
                }
                for layer in layer_set.MaterialLayers
            ],
        }
    return {"type": material.is_a(), "step": str(material)}


def write_issues(path: Path) -> None:
    fieldnames = [
        "issue_id", "scope", "object_guid", "current_evidence", "required_action_or_decision",
        "basis", "confidence", "review_required", "status", "stop_condition",
    ]
    rows = [
        {
            "issue_id": "WFIN-R01",
            "scope": "finish-boundary-and-missing-face-review",
            "object_guid": "06wFwLoDD6ie5iCTnc_yad; 3jha5L04zBjh0$pMl_tHLy; 04rs0EDjn2EvxytEQSxWRB; 3OVQygdDn17huGOgJJFTOY; 2nUpG$tzj0vR6tg5N9FJJ_; 2Rr8WvEy591Pj6dg4T6RQK",
            "current_evidence": "当前画面中的 0e0XOb$L18ZBVYJiQJrQ1p 与 3bMoS7bIT8wBmjqRvdPunu 已由用户确认为整件大白墙；剩余跨界对象及主卫干区返口/低墙/门垛/踢脚留待统一选材",
            "required_action_or_decision": "暂停 WFIN 几何和材料判断；最终统一选材时恢复 Blender 真深度审核，再决定剩余分界与缺失饰面",
            "basis": "用户明确要求当前先往后推进，最后统一选材后再判断；Space 方盒相交不再自动覆盖人审结论",
            "confidence": "1.00",
            "review_required": "yes",
            "status": "deferred_by_user",
            "stop_condition": "统一选材和人审恢复前，不得继续拆分、替换材料或新增 CLADDING",
        },
        {
            "issue_id": "WFIN-R02",
            "scope": "tadelakt-system",
            "object_guid": "",
            "current_evidence": "当前候选分段中 24 段为 Tadelakt 意图；最终产品选择由用户统一延后",
            "required_action_or_decision": "确认产品系统、颜色样板、完成面总厚、基层和湿区防水/收口节点",
            "basis": "参考链接只支持材质方向与湿区适用背景，不能替代本项目施工参数",
            "confidence": "1.00",
            "review_required": "yes",
            "status": "decision_required",
            "stop_condition": "参数未确认前只写候选意图 Pset，不替换材料层或发布节点",
        },
        {
            "issue_id": "WFIN-R03",
            "scope": "white-wall-system",
            "object_guid": "",
            "current_evidence": "当前候选分段中 28 段为大白墙；其中两件跨 Space 对象已由用户确认为整件大白墙，最终产品选择统一延后",
            "required_action_or_decision": "确认涂料体系、白色样板/光泽、基层处理和完成面总厚",
            "basis": "大白墙房间范围已确认，但产品与施工层次尚未确认",
            "confidence": "1.00",
            "review_required": "yes",
            "status": "decision_required",
            "stop_condition": "参数未确认前不发布材料表和收口节点",
        },
        {
            "issue_id": "WFIN-R04",
            "scope": "installation-anchor",
            "object_guid": "; ".join(sorted(WFIN_HANDOFF_IDS)),
            "current_evidence": "C003 精确移交 1 个 Painting Proxy 与 2 个卫生间飘窗 CLADDING；原点仍超过 0.1 mm",
            "required_action_or_decision": "在完成面与构件角色明确后生成纯原点重设候选，保持世界几何不变",
            "basis": "两个 CLADDING 当前原点位于自身顶点但不是整数；Painting 的安装基准尚未证明",
            "confidence": "1.00",
            "review_required": "yes",
            "status": "anchor_candidate_required",
            "stop_condition": "未证明锚点在实际几何且世界几何 0.0 mm 前不得写 placement",
        },
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    args = parse_args()
    rows = read_register(args.register)
    segments_by_covering = read_segments(args.segments)
    model = ifcopenshell.open(args.input)
    if model.schema != "IFC4":
        raise RuntimeError(f"expected IFC4, got {model.schema}")
    settings = ifcopenshell.geom.settings()
    settings.set(settings.USE_WORLD_COORDS, True)
    source_roots = {root.GlobalId for root in model.by_type("IfcRoot")}
    source_geometry = {}
    source_materials = {}
    source_long_names = {}
    for row in rows:
        covering = model.by_guid(row["covering_global_id"])
        if covering is None or not covering.is_a("IfcCovering"):
            raise RuntimeError(f"missing IfcCovering {row['covering_global_id']}")
        for space_id, space_name in zip(row["space_global_id"].split("; "), row["space_long_name"].split("; ")):
            space = model.by_guid(space_id)
            if space is None or not space.is_a("IfcSpace"):
                raise RuntimeError(f"missing IfcSpace {space_id}")
            if str(space.LongName or space.Name or "") != space_name:
                raise RuntimeError(f"Space LongName drift for {space.GlobalId}")
            reference = ifcopenshell.util.element.get_psets(space).get("Pset_SpaceCommon", {}).get("Reference")
            if not reference:
                raise RuntimeError(f"Space Reference missing for {space.GlobalId}")
        if PSET_NAME in ifcopenshell.util.element.get_psets(covering):
            raise RuntimeError(f"formal IFC already contains {PSET_NAME} on {covering.GlobalId}")
        source_geometry[covering.GlobalId] = world_vertices(settings, covering)
        source_materials[covering.GlobalId] = material_signature(covering)
        source_long_names[space.GlobalId] = str(space.LongName or "")

    applied = []
    for row in rows:
        covering = model.by_guid(row["covering_global_id"])
        spaces = [model.by_guid(space_id) for space_id in row["space_global_id"].split("; ")]
        references = [str(ifcopenshell.util.element.get_psets(space)["Pset_SpaceCommon"]["Reference"]) for space in spaces]
        covering_segments = segments_by_covering[covering.GlobalId]
        pset = ifcopenshell.api.run("pset.add_pset", model, product=covering, name=PSET_NAME)
        pset.GlobalId = deterministic_guid(covering.GlobalId, "PSET")
        relation = next(
            rel for rel in covering.IsDefinedBy
            if rel.is_a("IfcRelDefinesByProperties") and rel.RelatingPropertyDefinition == pset
        )
        relation.GlobalId = deterministic_guid(covering.GlobalId, "REL")
        contains_tadelakt = any(segment["candidate_finish_code"] == "TADELAKT" for segment in covering_segments)
        segment_summary = " | ".join(
            f'{segment["segment_id"]}:{segment["candidate_finish_code"]}:{segment["start_mm"]}-{segment["end_mm"]}mm'
            for segment in covering_segments
        )
        ifcopenshell.api.run(
            "pset.edit_pset",
            model,
            pset=pset,
            properties={
                "FinishCode": row["candidate_finish_code"],
                "FinishName": row["candidate_finish"],
                "AssignedSpaceReference": "; ".join(references),
                "AssignedSpaceLongName": row["space_long_name"],
                "AssignmentBasis": "USER_CONFIRMED_ROOM_SCOPE_AND_GRID_SPACE_BOUNDARY_SEGMENTS",
                "CandidateSegmentCount": len(covering_segments),
                "CandidateSegmentSummary": segment_summary,
                "ReviewStatus": "DEFERRED_MATERIAL_REVIEW" if row["candidate_finish_code"] == "MULTI_FINISH_SPLIT_REQUIRED" else "CONFIRMED_CURRENT_INTENT",
                "ExistingMaterialAssociationPreserved": True,
                "ProductSystemPending": True,
                "ColourSamplePending": True,
                "TotalThicknessPending": True,
                "SubstratePending": True,
                "WaterproofingDetailPending": contains_tadelakt,
                "FormalIfcWriteAllowed": False,
            },
        )
        applied.append(
            {
                "covering_global_id": covering.GlobalId,
                "finish_code": row["candidate_finish_code"],
                "space_reference": "; ".join(references),
                "space_long_name": row["space_long_name"],
                "pset_global_id": pset.GlobalId,
                "relation_global_id": relation.GlobalId,
            }
        )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    model.write(args.output)
    candidate = ifcopenshell.open(args.output)
    candidate_settings = ifcopenshell.geom.settings()
    candidate_settings.set(candidate_settings.USE_WORLD_COORDS, True)
    expected_new_roots = {
        value for record in applied for value in (record["pset_global_id"], record["relation_global_id"])
    }
    candidate_roots = {root.GlobalId for root in candidate.by_type("IfcRoot")}
    if candidate_roots != source_roots | expected_new_roots:
        raise RuntimeError("WFIN candidate Root boundary drift")

    maximum_change = 0.0
    pset_count = 0
    for row in rows:
        covering = candidate.by_guid(row["covering_global_id"])
        pset = ifcopenshell.util.element.get_psets(covering).get(PSET_NAME, {})
        if pset.get("FinishCode") != row["candidate_finish_code"]:
            raise RuntimeError(f"WFIN FinishCode mismatch for {covering.GlobalId}")
        if pset.get("AssignedSpaceLongName") != row["space_long_name"]:
            raise RuntimeError(f"WFIN Space mismatch for {covering.GlobalId}")
        if pset.get("FormalIfcWriteAllowed") is not False:
            raise RuntimeError(f"WFIN formal write boundary missing for {covering.GlobalId}")
        pset_count += 1
        before = source_geometry[covering.GlobalId]
        after = world_vertices(candidate_settings, covering)
        if before.shape != after.shape:
            raise RuntimeError(f"WFIN geometry topology drift for {covering.GlobalId}")
        maximum_change = max(maximum_change, float(np.max(np.abs(before - after), initial=0.0)))
        if material_signature(covering) != source_materials[covering.GlobalId]:
            raise RuntimeError(f"WFIN material association changed for {covering.GlobalId}")
    if maximum_change > args.tolerance_mm:
        raise RuntimeError(f"WFIN geometry changed by {maximum_change} mm")

    write_issues(args.issues)
    report = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "source": {"path": str(args.input), "sha256": sha256(args.input), "schema": model.schema},
        "candidate": {"path": str(args.output), "sha256": sha256(args.output)},
        "automatic_formal_ifc_write_allowed": False,
        "intent_pset": {
            "name": PSET_NAME,
            "standard_status": "project-specific design-intent Pset; not a buildingSMART standard Pset",
        },
        "cladding_count": len(rows),
        "intent_pset_count": pset_count,
        "single_finish_object_count": sum(row["candidate_finish_code"] != "MULTI_FINISH_SPLIT_REQUIRED" for row in rows),
        "mixed_finish_object_count": sum(row["candidate_finish_code"] == "MULTI_FINISH_SPLIT_REQUIRED" for row in rows),
        "finish_segment_count": sum(len(value) for value in segments_by_covering.values()),
        "tadelakt_segment_count": sum(segment["candidate_finish_code"] == "TADELAKT" for value in segments_by_covering.values() for segment in value),
        "white_wall_segment_count": sum(segment["candidate_finish_code"] == "WHITE_WALL" for value in segments_by_covering.values() for segment in value),
        "new_root_count": len(expected_new_roots),
        "maximum_world_vertex_change_mm": maximum_change,
        "material_associations_preserved": True,
        "open_issue_count": 4,
        "formal_write_blockers": ["WFIN-R01", "WFIN-R02", "WFIN-R03", "WFIN-R04"],
        "qa": {
            "segment_partition_complete": True,
            "mixed_finish_whole_object_assignment_blocked": True,
            "space_references_present": True,
            "root_boundary_exact": True,
            "geometry_within_tolerance": maximum_change <= args.tolerance_mm,
            "materials_unchanged": True,
            "formal_ifc_unchanged": True,
        },
        "applied": applied,
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(
        f"WFIN semantics candidate: {len(rows)} intent Psets, 52 finish segments, 1 deferred split-review object, "
        f"geometry max change {maximum_change:.6f} mm, materials preserved, formal IFC unchanged"
    )


if __name__ == "__main__":
    main()
