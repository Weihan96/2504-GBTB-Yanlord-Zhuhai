#!/usr/bin/env python3
"""Create corrected TRAP01 detail Drawings with project-axis-aware linework planes."""

from __future__ import annotations

import json
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import bpy
import ifcopenshell.api
import ifcopenshell.util.element
from bonsai import tool


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "pipeline/scripts"))
import create_trap01_internal_detail_drawings as base  # noqa: E402


PRODUCT_DIR = ROOT / "output/review/highpoly-types/trap01"
DERIVED_IFC = PRODUCT_DIR / "Geberit-151.116.11.1-TRAP01-derived-drawing.ifc"
FORMAL_IFC = ROOT / "2504 GBTB Yanlord Zhuhai.ifc"
OUTPUT_DIR = PRODUCT_DIR / "bonsai-drawings/cabinet-internal-detail-v2"
PREPERSIST = OUTPUT_DIR / "TRAP01-cabinet-internal-detail-v2-prepersist.json"
CANDIDATE = PRODUCT_DIR / "candidate-representations.json"
VIEW_DEFINITIONS = {
    "plan": {"drawing_name": "TRAP01-CABINET-INTERNAL-DETAIL-V2-PLAN", "target_view": "PLAN_VIEW", "location_hint": "PLAN"},
    "front": {"drawing_name": "TRAP01-CABINET-INTERNAL-DETAIL-V2-FRONT", "target_view": "ELEVATION_VIEW", "location_hint": "EAST"},
    "side": {"drawing_name": "TRAP01-CABINET-INTERNAL-DETAIL-V2-SIDE", "target_view": "ELEVATION_VIEW", "location_hint": "SOUTH"},
}
# The target has identity rotation.  Project EAST looks along local X and must
# therefore use the local YZ (side) projection; project SOUTH looks along local
# Y and must use the local XZ (front) projection.
SOURCE_VIEW_BY_PROJECT_VIEW = {"plan": "plan", "front": "side", "side": "front"}
base.VIEW_DEFINITIONS = VIEW_DEFINITIONS


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    if base.sha256(FORMAL_IFC) != base.FORMAL_SHA256:
        raise RuntimeError("formal IFC hash mismatch before TRAP01 mapped detail workflow")
    if Path(tool.Ifc.get_path()).resolve() != DERIVED_IFC.resolve():
        raise RuntimeError(f"wrong IFC loaded: {tool.Ifc.get_path()}")
    model = tool.Ifc.get()
    names = {definition["drawing_name"] for definition in VIEW_DEFINITIONS.values()}
    if any(item.ObjectType == "DRAWING" and item.Name in names for item in model.by_type("IfcAnnotation")):
        raise RuntimeError("TRAP01 mapped detail Drawings already exist")
    target = model.by_guid(base.TARGET_GLOBAL_ID)
    target_obj = tool.Ifc.get_object(target) if target else None
    if target is None or target_obj is None:
        raise RuntimeError("TRAP01 target missing")
    candidate = json.loads(CANDIDATE.read_text(encoding="utf-8"))
    target_bbox = base.world_bbox(target_obj)
    centre = [(target_bbox[0][axis] + target_bbox[1][axis]) / 2 for axis in range(3)]
    local_crop = ((centre[0] - 0.65, centre[1] - 0.65, 0.20), (centre[0] + 0.65, centre[1] + 0.65, 1.10))
    context_elements = []
    for element in model.by_type("IfcElement"):
        if element == target or element.GlobalId in base.OCCLUDERS:
            continue
        obj = tool.Ifc.get_object(element)
        if obj is None:
            continue
        try:
            if base.intersects(base.world_bbox(obj), local_crop):
                context_elements.append(element)
        except Exception:
            continue
    components = [item for item in model.by_type("IfcAnnotation") if item.ObjectType == "DRAWING_COMPONENT"]
    roles = {
        ifcopenshell.util.element.get_pset(item, base.COMPONENT_PSET).get("ComponentRole"): item
        for item in components
    }
    if set(roles) != {"fixed_body", "horizontal_adjustable", "vertical_adjustable"}:
        raise RuntimeError("TRAP01 component semantics drifted")
    override = base.view3d_override()
    records = []
    for project_view, definition in VIEW_DEFINITIONS.items():
        source_view = SOURCE_VIEW_BY_PROJECT_VIEW[project_view]
        output_svg = OUTPUT_DIR / f"{definition['drawing_name']}.svg"
        drawing, camera, width, height, clip_end = base.add_detail_drawing(
            model, definition, centre, context_elements, output_svg
        )
        with bpy.context.temp_override(**override):
            activate = bpy.ops.bim.activate_drawing(drawing=drawing.id(), should_view_from_camera=False)
        if activate != {"FINISHED"}:
            raise RuntimeError(f"TRAP01 mapped detail {project_view} activation failed")
        annotation = base.add_current_annotation(
            model, drawing, target, target_obj, project_view,
            candidate["views"][source_view]["proxy_paths_mm"],
            coordinate_view=source_view,
        )
        cprops = tool.Drawing.get_camera_props(camera)
        cprops.has_annotation = True
        pset = ifcopenshell.util.element.get_pset(drawing, "EPset_Drawing")
        ifcopenshell.api.pset.edit_pset(model, pset=model.by_id(pset["id"]), properties={"HasAnnotation": True})
        dprops = tool.Drawing.get_document_props()
        dprops.should_use_underlay_cache = False
        dprops.should_use_linework_cache = False
        dprops.should_use_annotation_cache = False
        with bpy.context.temp_override(**override):
            result = bpy.ops.bim.create_drawing(print_all=False, open_viewer=False, sync=False)
        if result != {"FINISHED"} or not output_svg.is_file() or not output_svg.stat().st_size:
            raise RuntimeError(f"TRAP01 mapped detail {project_view} Create Drawing failed: {result}")
        style = base.style_svg(output_svg, annotation.GlobalId)
        svg = base.inspect_svg(output_svg, base.TARGET_GLOBAL_ID, annotation.GlobalId)
        if not svg["no_target_or_annotation_duplicate"] or not svg["all_occluders_absent"]:
            raise RuntimeError(f"TRAP01 mapped detail {project_view} SVG gate failed: {svg}")
        records.append({
            "view": project_view,
            "source_product_view": source_view,
            "drawing_global_id": drawing.GlobalId,
            "annotation_global_id": annotation.GlobalId,
            "persisted_path_count_expected": len(candidate["views"][source_view]["proxy_paths_mm"]),
            "output_svg": str(output_svg),
            "output_svg_bytes": output_svg.stat().st_size,
            "output_svg_sha256": base.sha256(output_svg),
            "style": style,
            "svg": svg,
            "camera": {"type": camera.data.type, "width_m": width, "height_m": height, "clip_start_m": camera.data.clip_start, "clip_end_m": clip_end, "scale": "1:5", "matrix_world": [list(row) for row in camera.matrix_world]},
            "create_result": sorted(result),
        })
    data = {
        "schema_version": 1,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "status": "awaiting_save_ifc_file_and_reload_verification",
        "task": "TRAP01 cabinet-internal orthographic detail with project-axis-aware accepted proxy linework",
        "courseEvidence": {
            "evidence_mode": "embedded-course-index",
            "lesson": "085000 Introduction to Drawings",
            "timestamps": ["01:10 camera boundary/scale", "01:52 drawing depth", "01:59-02:13 Create Drawing/SVG", "02:36-03:16 Element Filters"],
            "course_fact": "Camera, depth and Include/Exclude filters are configured before Create Drawing; generated SVG is inspected.",
            "screen_observation": None,
            "provenance_note": "No private course screenshot was directly viewed.",
        },
        "plan": ["inspect", "map product-local views to project camera axes", "create filtered 1:5 Drawings", "Create Drawing", "save_ifc_file", "reload", "render and verify"],
        "preState": {"derived_ifc_sha256": base.sha256(DERIVED_IFC), "formal_ifc_sha256": base.sha256(FORMAL_IFC), "component_count": len(components)},
        "execution": {"provider": "public bonsai-mcp", "capability": "execute_blender_code", "generator": "bpy.ops.bim.create_drawing", "linework_mode": "OPENCASCADE", "blender": bpy.app.version_string},
        "save_boundary": {"derived_ifc": str(DERIVED_IFC), "formal_ifc": str(FORMAL_IFC), "geometry_mutation_allowed": False, "drawing_filter_annotation_only": True},
        "target": {"global_id": base.TARGET_GLOBAL_ID, "world_bbox_m": [list(target_bbox[0]), list(target_bbox[1])], "world_centre_m": centre, "placement_rotation": "identity"},
        "view_mapping": SOURCE_VIEW_BY_PROJECT_VIEW,
        "filters": {"include_count": len(context_elements), "include_global_ids": [item.GlobalId for item in context_elements], "include_ifc_class_counts": dict(sorted(Counter(item.is_a() for item in context_elements).items())), "excluded_occluders": [{"global_id": guid, "reason": reason} for guid, reason in base.OCCLUDERS.items()]},
        "adjustable_components": {role: {"global_id": item.GlobalId, "pset": ifcopenshell.util.element.get_pset(item, base.COMPONENT_PSET)} for role, item in roles.items()},
        "outputs": {"views": records},
        "formal_ifc_bytes_unchanged": base.sha256(FORMAL_IFC) == base.FORMAL_SHA256,
    }
    PREPERSIST.write_text(json.dumps(data, indent=2, ensure_ascii=False, default=str) + "\n", encoding="utf-8")
    print(json.dumps({"prepersist": str(PREPERSIST), "views": records, "ready_for_save_ifc_file": True}, indent=2))


if __name__ == "__main__":
    main()
