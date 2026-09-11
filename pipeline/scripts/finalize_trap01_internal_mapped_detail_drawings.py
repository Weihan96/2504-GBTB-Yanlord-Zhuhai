#!/usr/bin/env python3
"""Finalize the reloaded, project-axis-aware TRAP01 detail Drawing evidence."""

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
OUTPUT_DIR = PRODUCT_DIR / "bonsai-drawings/cabinet-internal-detail-v2"
PREPERSIST = OUTPUT_DIR / "TRAP01-cabinet-internal-detail-v2-prepersist.json"
EVIDENCE = OUTPUT_DIR / "TRAP01-cabinet-internal-detail-v2-create-drawing-evidence.json"
SESSION_BLEND = PRODUCT_DIR / "Geberit-151.116.11.1-TRAP01-cabinet-internal-detail-v2.blend"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def path_count(annotation):
    return sum(
        len(item.Elements)
        for representation in annotation.Representation.Representations
        for item in representation.Items
        if item.is_a("IfcGeometricCurveSet")
    )


def main():
    if sha256(FORMAL_IFC) != FORMAL_SHA256 or Path(tool.Ifc.get_path()).resolve() != DERIVED_IFC.resolve():
        raise RuntimeError("TRAP01 v2 reload save boundary failed")
    data = json.loads(PREPERSIST.read_text(encoding="utf-8"))
    model = tool.Ifc.get()
    verified = []
    for source in data["outputs"]["views"]:
        drawing = model.by_guid(source["drawing_global_id"])
        annotation = model.by_guid(source["annotation_global_id"])
        if drawing is None or annotation is None or path_count(annotation) != source["persisted_path_count_expected"]:
            raise RuntimeError(f"TRAP01 v2 persisted state missing for {source['view']}")
        pset = ifcopenshell.util.element.get_pset(drawing, "EPset_Drawing")
        expected = {item["global_id"] for item in data["filters"]["excluded_occluders"]}
        excluded = set(filter(None, pset.get("Exclude", "").split(",")))
        if excluded != expected:
            raise RuntimeError(f"TRAP01 v2 Exclude filter drifted for {source['view']}")
        verified.append({
            "view": source["view"], "source_product_view": source["source_product_view"],
            "drawing_global_id": drawing.GlobalId, "annotation_global_id": annotation.GlobalId,
            "persisted_path_count": path_count(annotation),
            "persisted_include_count": len(list(filter(None, pset.get("Include", "").split(",")))),
            "persisted_excluded_global_ids": sorted(excluded),
        })
    components = [item for item in model.by_type("IfcAnnotation") if item.ObjectType == "DRAWING_COMPONENT"]
    if len(components) != 3:
        raise RuntimeError("TRAP01 component count changed")
    bpy.ops.wm.save_as_mainfile(filepath=str(SESSION_BLEND), check_existing=False)
    data.update({
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "status": "persisted_reloaded_verified",
        "persistence": {
            "provider_capability": "save_ifc_file",
            "method": "public bonsai-mcp save_ifc_file(overwrite=true, reload=true)",
            "reload_verified": True,
            "derived_ifc": str(DERIVED_IFC), "derived_ifc_bytes": DERIVED_IFC.stat().st_size, "derived_ifc_sha256": sha256(DERIVED_IFC),
            "session_blend": str(SESSION_BLEND), "session_blend_bytes": SESSION_BLEND.stat().st_size, "session_blend_sha256": sha256(SESSION_BLEND),
        },
        "postState": {"formal_ifc_sha256": sha256(FORMAL_IFC), "formal_ifc_bytes_unchanged": True, "detail_drawing_count": len(verified), "detail_annotation_count": len(verified), "component_count": len(components), "reloaded_views": verified},
        "tests": {"all_create_drawing_finished": all(item["create_result"] == ["FINISHED"] for item in data["outputs"]["views"]), "all_occluders_absent": all(item["svg"]["all_occluders_absent"] for item in data["outputs"]["views"]), "all_target_bodies_suppressed": all(item["svg"]["target_ifc_body_group_count"] == 0 for item in data["outputs"]["views"]), "all_annotations_once": all(item["svg"]["current_annotation_group_count"] == 1 for item in data["outputs"]["views"]), "blue_dashed_reference_displayed": False, "formal_ifc_bytes_unchanged": True},
        "visual": {"preview_manifest": str(PRODUCT_DIR / "TRAP01-cabinet-internal-detail-v2-manifest.json"), "rendered_preview_pending": True},
        "verdict": "pass_pending_rendered_preview", "pass": False,
    })
    EVIDENCE.write_text(json.dumps(data, indent=2, ensure_ascii=False, default=str) + "\n", encoding="utf-8")
    print(json.dumps({"evidence": str(EVIDENCE), "derived_ifc_sha256": sha256(DERIVED_IFC), "session_blend": str(SESSION_BLEND), "verified": verified}, indent=2))


main()
