#!/usr/bin/env python3
"""Create corrected TRAP01 cabinet-internal detail Drawings from the current IFC Body.

The first detail attempt used view-specific review annotations whose local
planes collapsed in two project views.  This corrected set includes the actual
installed target Body exactly once and uses no review LINEWORK annotation.
"""

from __future__ import annotations

import hashlib
import json
import sys
import xml.etree.ElementTree as ET
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
FORMAL_IFC = ROOT / "2504 GBTB Yanlord Zhuhai.ifc"
FORMAL_SHA256 = "7a50b87e8f48a7c2bfbaa9f1dabdfca7155fa684b2325c66aa1ae15c8ab0a25c"
DERIVED_IFC = PRODUCT_DIR / "Geberit-151.116.11.1-TRAP01-derived-drawing.ifc"
OUTPUT_DIR = PRODUCT_DIR / "bonsai-drawings/cabinet-internal-body-detail"
PREPERSIST = OUTPUT_DIR / "TRAP01-cabinet-internal-body-detail-prepersist.json"
TARGET_GLOBAL_ID = base.TARGET_GLOBAL_ID
OCCLUDERS = base.OCCLUDERS
BLACK = base.BLACK
GREY = base.GREY
GEOMETRY_TAGS = base.GEOMETRY_TAGS
COMPONENT_PSET = base.COMPONENT_PSET
VIEW_DEFINITIONS = {
    "plan": {"drawing_name": "TRAP01-CABINET-INTERNAL-BODY-DETAIL-PLAN", "target_view": "PLAN_VIEW", "location_hint": "PLAN"},
    "front": {"drawing_name": "TRAP01-CABINET-INTERNAL-BODY-DETAIL-FRONT", "target_view": "ELEVATION_VIEW", "location_hint": "EAST"},
    "side": {"drawing_name": "TRAP01-CABINET-INTERNAL-BODY-DETAIL-SIDE", "target_view": "ELEVATION_VIEW", "location_hint": "SOUTH"},
}
base.VIEW_DEFINITIONS = VIEW_DEFINITIONS


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def local_name(value: str) -> str:
    return value.rsplit("}", 1)[-1].split(":")[-1]


def has_identity(element, guid):
    attrs = {local_name(key): value for key, value in element.attrib.items()}
    classes = element.attrib.get("class", "").split()
    return attrs.get("guid") == guid or f"GlobalId-{guid}" in classes


def style_and_inspect(svg_path: Path):
    raw_sha = sha256(svg_path)
    ET.register_namespace("", "http://www.w3.org/2000/svg")
    ET.register_namespace("ifc", "http://www.ifcopenshell.org/ns")
    tree = ET.parse(svg_path)
    root = tree.getroot()
    parents = {child: parent for parent in root.iter() for child in parent}
    target_groups = [element for element in root.iter() if has_identity(element, TARGET_GLOBAL_ID)]
    if not target_groups:
        raise RuntimeError("TRAP01 Body group missing from corrected detail SVG")
    target_geometry = {
        element
        for group in target_groups
        for element in group.iter()
        if local_name(element.tag) in GEOMETRY_TAGS
    }
    for group in target_groups:
        parent = parents.get(group)
        if parent is not None:
            parent.remove(group)
            parent.append(group)
    black = grey = 0
    for element in root.iter():
        if local_name(element.tag) not in GEOMETRY_TAGS:
            continue
        is_target = element in target_geometry
        colour = f"stroke:{BLACK};stroke-width:0.45;fill:none" if is_target else f"stroke:{GREY};stroke-width:0.20;fill:none;stroke-opacity:0.62"
        existing = element.attrib.get("style", "").rstrip(";")
        element.attrib["style"] = f"{existing};{colour}" if existing else colour
        black += int(is_target)
        grey += int(not is_target)
    root.attrib.update({
        "data-create-drawing-result": "FINISHED",
        "data-trap01-detail": "cabinet-internal-orthographic-current-body",
        "data-current-installed-body-global-id": TARGET_GLOBAL_ID,
        "data-blue-dashed-reference-displayed": "false",
        "data-review-annotation-displayed": "false",
        "data-occluder-global-ids": ",".join(OCCLUDERS),
    })
    tree.write(svg_path, encoding="utf-8", xml_declaration=True)
    root = ET.parse(svg_path).getroot()
    target_group_count = int(any(has_identity(element, TARGET_GLOBAL_ID) for element in root.iter()))
    annotation_group_count = sum(
        1 for element in root.iter()
        if "current-installed" in element.attrib.get("class", "")
    )
    excluded_counts = {
        guid: int(any(has_identity(element, guid) for element in root.iter()))
        for guid in OCCLUDERS
    }
    geometry_count = sum(local_name(element.tag) in GEOMETRY_TAGS for element in root.iter())
    result = {
        "root_tag": local_name(root.tag),
        "data_scale": root.attrib.get("data-scale"),
        "view_box": root.attrib.get("viewBox"),
        "geometry_element_count": geometry_count,
        "target_ifc_body_group_count": target_group_count,
        "current_annotation_group_count": annotation_group_count,
        "excluded_occluder_group_counts": excluded_counts,
        "body_once_annotation_zero": target_group_count == 1 and annotation_group_count == 0,
        "all_occluders_absent": all(value == 0 for value in excluded_counts.values()),
        "black_geometry_count": black,
        "grey_geometry_count": grey,
        "raw_sha256": raw_sha,
    }
    if not result["body_once_annotation_zero"] or not result["all_occluders_absent"] or black == 0 or grey == 0:
        raise RuntimeError(f"TRAP01 corrected detail SVG gate failed: {result}")
    return result


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    if sha256(FORMAL_IFC) != FORMAL_SHA256:
        raise RuntimeError("formal IFC hash mismatch before corrected TRAP01 detail Drawings")
    if Path(tool.Ifc.get_path()).resolve() != DERIVED_IFC.resolve():
        raise RuntimeError(f"wrong IFC loaded: {tool.Ifc.get_path()}")
    model = tool.Ifc.get()
    names = {definition["drawing_name"] for definition in VIEW_DEFINITIONS.values()}
    if any(item.ObjectType == "DRAWING" and item.Name in names for item in model.by_type("IfcAnnotation")):
        raise RuntimeError("corrected TRAP01 body detail Drawings already exist")
    target = model.by_guid(TARGET_GLOBAL_ID)
    target_obj = tool.Ifc.get_object(target) if target else None
    if target is None or target_obj is None:
        raise RuntimeError("TRAP01 target Body is missing")
    target_bbox = base.world_bbox(target_obj)
    centre = [(target_bbox[0][axis] + target_bbox[1][axis]) / 2 for axis in range(3)]
    local_crop = ((centre[0] - 0.65, centre[1] - 0.65, 0.20), (centre[0] + 0.65, centre[1] + 0.65, 1.10))
    retained_context = []
    for element in model.by_type("IfcElement"):
        if element == target or element.GlobalId in OCCLUDERS:
            continue
        obj = tool.Ifc.get_object(element)
        if obj is None:
            continue
        try:
            if base.intersects(base.world_bbox(obj), local_crop):
                retained_context.append(element)
        except Exception:
            continue
    include_elements = [target, *retained_context]
    components = [item for item in model.by_type("IfcAnnotation") if item.ObjectType == "DRAWING_COMPONENT"]
    roles = {
        ifcopenshell.util.element.get_pset(item, COMPONENT_PSET).get("ComponentRole"): item
        for item in components
    }
    if set(roles) != {"fixed_body", "horizontal_adjustable", "vertical_adjustable"}:
        raise RuntimeError("TRAP01 component semantics drifted")
    pre_state = {
        "derived_ifc_sha256": sha256(DERIVED_IFC),
        "formal_ifc_sha256": sha256(FORMAL_IFC),
        "drawing_count": len([item for item in model.by_type("IfcAnnotation") if item.ObjectType == "DRAWING"]),
        "component_count": len(components),
    }
    override = base.view3d_override()
    records = []
    for view, definition in VIEW_DEFINITIONS.items():
        output_svg = OUTPUT_DIR / f"{definition['drawing_name']}.svg"
        drawing, camera, width, height, clip_end = base.add_detail_drawing(
            model, definition, centre, include_elements, output_svg
        )
        with bpy.context.temp_override(**override):
            activate = bpy.ops.bim.activate_drawing(drawing=drawing.id(), should_view_from_camera=False)
        if activate != {"FINISHED"}:
            raise RuntimeError(f"TRAP01 corrected detail {view} activation failed")
        cprops = tool.Drawing.get_camera_props(camera)
        cprops.has_annotation = False
        cprops.linework_mode = "FREESTYLE"
        drawing_pset = ifcopenshell.util.element.get_pset(drawing, "EPset_Drawing")
        ifcopenshell.api.pset.edit_pset(
            model,
            pset=model.by_id(drawing_pset["id"]),
            properties={"HasAnnotation": False, "LineworkMode": "FREESTYLE"},
        )
        dprops = tool.Drawing.get_document_props()
        dprops.should_use_underlay_cache = False
        dprops.should_use_linework_cache = False
        dprops.should_use_annotation_cache = False
        with bpy.context.temp_override(**override):
            create_result = bpy.ops.bim.create_drawing(print_all=False, open_viewer=False, sync=False)
        if create_result != {"FINISHED"} or not output_svg.is_file() or output_svg.stat().st_size == 0:
            raise RuntimeError(f"TRAP01 corrected detail {view} Create Drawing failed: {create_result}")
        svg = style_and_inspect(output_svg)
        cache = output_svg.parent / "cache" / f"{output_svg.stem}-linework.svg"
        if not cache.is_file() or not cache.stat().st_size:
            raise RuntimeError(f"TRAP01 corrected detail {view} linework cache missing")
        records.append({
            "view": view,
            "drawing_global_id": drawing.GlobalId,
            "drawing_name": drawing.Name,
            "output_svg": str(output_svg),
            "output_svg_bytes": output_svg.stat().st_size,
            "output_svg_sha256": sha256(output_svg),
            "linework_cache": str(cache),
            "linework_cache_sha256": sha256(cache),
            "svg": svg,
            "camera": {
                "type": camera.data.type,
                "width_m": width,
                "height_m": height,
                "clip_start_m": camera.data.clip_start,
                "clip_end_m": clip_end,
                "scale": "1:5",
                "matrix_world": [list(row) for row in camera.matrix_world],
            },
            "create_result": sorted(create_result),
        })
    data = {
        "schema_version": 1,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "status": "awaiting_save_ifc_file_and_reload_verification",
        "task": "corrected TRAP01 cabinet-internal orthographic detail using the actual current installed IFC Body once per view",
        "courseEvidence": {
            "evidence_mode": "embedded-course-index",
            "lesson": "085000 Introduction to Drawings",
            "timestamps": ["01:10 camera boundary/scale", "01:52 drawing depth", "01:59-02:13 Create Drawing/SVG", "02:36-03:16 Include/Exclude filters"],
            "course_fact": "Configure orthographic camera, depth and Element Filters before Create Drawing; verify the generated SVG.",
            "screen_observation": None,
            "provenance_note": "No private course screenshot was directly viewed.",
        },
        "plan": ["inspect", "resolve occluders", "create 1:5 Body-only detail Drawings", "Create Drawing", "save_ifc_file", "reload", "inspect SVG and rendered previews"],
        "preState": pre_state,
        "execution": {"provider": "public bonsai-mcp", "capability": "execute_blender_code", "generator": "bpy.ops.bim.create_drawing", "linework_mode": "FREESTYLE", "blender": bpy.app.version_string},
        "save_boundary": {"derived_ifc": str(DERIVED_IFC), "formal_ifc": str(FORMAL_IFC), "geometry_mutation_allowed": False, "drawing_filter_only": True},
        "target": {"global_id": TARGET_GLOBAL_ID, "world_bbox_m": [list(target_bbox[0]), list(target_bbox[1])], "world_centre_m": centre},
        "filters": {
            "include_count": len(include_elements),
            "include_global_ids": [item.GlobalId for item in include_elements],
            "include_ifc_class_counts": dict(sorted(Counter(item.is_a() for item in include_elements).items())),
            "excluded_occluders": [{"global_id": guid, "reason": reason} for guid, reason in OCCLUDERS.items()],
        },
        "adjustable_components": {role: {"global_id": item.GlobalId, "pset": ifcopenshell.util.element.get_pset(item, COMPONENT_PSET)} for role, item in roles.items()},
        "outputs": {"views": records},
        "formal_ifc_bytes_unchanged": sha256(FORMAL_IFC) == FORMAL_SHA256,
    }
    PREPERSIST.write_text(json.dumps(data, indent=2, ensure_ascii=False, default=str) + "\n", encoding="utf-8")
    print(json.dumps({"prepersist": str(PREPERSIST), "views": records, "ready_for_save_ifc_file": True}, indent=2))


if __name__ == "__main__":
    main()
