#!/usr/bin/env python3
"""Verify reloaded TRAP01 detail Drawing state and finalize audit evidence."""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import bpy
import ifcopenshell.util.element
from bonsai import tool


ROOT = Path(__file__).resolve().parents[2]
PRODUCT_DIR = ROOT / "output/review/highpoly-types/trap01"
FORMAL_IFC = ROOT / "2504 GBTB Yanlord Zhuhai.ifc"
FORMAL_SHA256 = "7a50b87e8f48a7c2bfbaa9f1dabdfca7155fa684b2325c66aa1ae15c8ab0a25c"
DERIVED_IFC = PRODUCT_DIR / "Geberit-151.116.11.1-TRAP01-derived-drawing.ifc"
OUTPUT_DIR = PRODUCT_DIR / "bonsai-drawings/cabinet-internal-detail"
PREPERSIST = OUTPUT_DIR / "TRAP01-cabinet-internal-detail-prepersist.json"
EVIDENCE = OUTPUT_DIR / "TRAP01-cabinet-internal-detail-create-drawing-evidence.json"
SESSION_BLEND = PRODUCT_DIR / "Geberit-151.116.11.1-TRAP01-cabinet-internal-detail.blend"
TARGET_GLOBAL_ID = "2Ak2ma0lvBEA49UpplzUqi"
COMPONENT_PSET = "Pset_Trap01DrawingComponent"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def persisted_path_count(annotation):
    return sum(
        len(item.Elements)
        for representation in annotation.Representation.Representations
        for item in representation.Items
        if item.is_a("IfcGeometricCurveSet")
    )


def main():
    if sha256(FORMAL_IFC) != FORMAL_SHA256:
        raise RuntimeError("formal IFC hash mismatch during TRAP01 detail reload verification")
    if Path(tool.Ifc.get_path()).resolve() != DERIVED_IFC.resolve():
        raise RuntimeError(f"wrong IFC loaded after save/reload: {tool.Ifc.get_path()}")
    data = json.loads(PREPERSIST.read_text(encoding="utf-8"))
    model = tool.Ifc.get()
    verified = []
    for record in data["outputs"]["views"]:
        drawing = model.by_guid(record["drawing_global_id"])
        annotation = model.by_guid(record["annotation_global_id"])
        if drawing is None or annotation is None:
            raise RuntimeError(f"persisted detail state missing for {record['view']}")
        if persisted_path_count(annotation) != record["path_count"]:
            raise RuntimeError(f"persisted detail path count drifted for {record['view']}")
        pset = ifcopenshell.util.element.get_pset(drawing, "EPset_Drawing")
        excluded = set(filter(None, pset.get("Exclude", "").split(",")))
        expected_excluded = {item["global_id"] for item in data["filters"]["excluded_occluders"]}
        if excluded != expected_excluded:
            raise RuntimeError(f"persisted Exclude filter drifted for {record['view']}: {excluded}")
        verified.append({
            "view": record["view"],
            "drawing_global_id": drawing.GlobalId,
            "annotation_global_id": annotation.GlobalId,
            "persisted_path_count": persisted_path_count(annotation),
            "persisted_include_count": len(list(filter(None, pset.get("Include", "").split(",")))),
            "persisted_excluded_global_ids": sorted(excluded),
        })
    components = [item for item in model.by_type("IfcAnnotation") if item.ObjectType == "DRAWING_COMPONENT"]
    roles = {
        ifcopenshell.util.element.get_pset(item, COMPONENT_PSET).get("ComponentRole"): item.GlobalId
        for item in components
    }
    if set(roles) != {"fixed_body", "horizontal_adjustable", "vertical_adjustable"}:
        raise RuntimeError("TRAP01 component semantics changed during Drawing-only batch")
    bpy.ops.wm.save_as_mainfile(filepath=str(SESSION_BLEND), check_existing=False)
    data.update({
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "status": "persisted_reloaded_verified",
        "persistence": {
            "provider_capability": "save_ifc_file",
            "method": "public bonsai-mcp save_ifc_file(overwrite=true, reload=true)",
            "reload_verified": True,
            "derived_ifc": str(DERIVED_IFC),
            "derived_ifc_bytes": DERIVED_IFC.stat().st_size,
            "derived_ifc_sha256": sha256(DERIVED_IFC),
            "session_blend": str(SESSION_BLEND),
            "session_blend_bytes": SESSION_BLEND.stat().st_size,
            "session_blend_sha256": sha256(SESSION_BLEND),
        },
        "postState": {
            "formal_ifc_sha256": sha256(FORMAL_IFC),
            "formal_ifc_bytes_unchanged": sha256(FORMAL_IFC) == FORMAL_SHA256,
            "detail_drawing_count": len(verified),
            "detail_annotation_count": len(verified),
            "component_count": len(components),
            "component_global_ids_by_role": roles,
            "reloaded_views": verified,
        },
        "tests": {
            "all_create_drawing_finished": all(item["create_result"] == ["FINISHED"] for item in data["outputs"]["views"]),
            "all_occluders_absent_from_svg": all(item["svg"]["all_occluders_absent"] for item in data["outputs"]["views"]),
            "all_target_body_projection_counts_zero": all(item["svg"]["target_ifc_body_group_count"] == 0 for item in data["outputs"]["views"]),
            "all_current_annotation_counts_one": all(item["svg"]["current_annotation_group_count"] == 1 for item in data["outputs"]["views"]),
            "all_persisted_path_counts_match": all(item["persisted_path_count"] == next(source["path_count"] for source in data["outputs"]["views"] if source["view"] == item["view"]) for item in verified),
            "blue_dashed_reference_displayed": False,
            "formal_ifc_bytes_unchanged": sha256(FORMAL_IFC) == FORMAL_SHA256,
        },
        "visual": {
            "preview_manifest": str(PRODUCT_DIR / "TRAP01-cabinet-internal-detail-manifest.json"),
            "rendered_preview_pending": True,
        },
        "verdict": "pass_pending_rendered_preview",
        "pass": False,
    })
    EVIDENCE.write_text(json.dumps(data, indent=2, ensure_ascii=False, default=str) + "\n", encoding="utf-8")
    print(json.dumps({"evidence": str(EVIDENCE), "derived_ifc_sha256": sha256(DERIVED_IFC), "session_blend": str(SESSION_BLEND), "verified": verified}, indent=2))


main()
