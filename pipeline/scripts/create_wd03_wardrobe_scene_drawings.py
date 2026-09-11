#!/usr/bin/env python3
"""Create approved WD03 semantic linework in its actual project context."""

from __future__ import annotations

import contextlib
import hashlib
import json
import os
import sys
import traceback
import xml.etree.ElementTree as ET
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import addon_utils
import bpy
import ifcopenshell
import ifcopenshell.api
import ifcopenshell.util.element
from bonsai import tool
from bonsai.core import drawing as core_drawing


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "pipeline/scripts"))
import create_gessi316_54294_main_bathroom_drawing as shared  # noqa: E402


PRODUCT_DIR = ROOT / "output/review/highpoly-types/wd03"
FORMAL_IFC = ROOT / "2504 GBTB Yanlord Zhuhai.ifc"
FORMAL_SHA256 = "7a50b87e8f48a7c2bfbaa9f1dabdfca7155fa684b2325c66aa1ae15c8ab0a25c"
DERIVED_IFC = PRODUCT_DIR / "Poliform-Senzafine-WD03-derived-drawing.ifc"
WRITE_REPORT = PRODUCT_DIR / "Poliform-Senzafine-WD03-derived-drawing-report.json"
APPROVAL = ROOT / "pipeline/decisions/wd03-drawing-approval.json"
CANDIDATE = PRODUCT_DIR / "candidate-representations.json"
OUTPUT_DIR = PRODUCT_DIR / "bonsai-drawings/wardrobe"
EVIDENCE = OUTPUT_DIR / "WD03-WARDROBE-create-drawing-evidence.json"
SESSION_BLEND = PRODUCT_DIR / "Poliform-Senzafine-WD03-project-drawings.blend"
TARGET_GLOBAL_ID = "3cmikd9MTB$egM5KQaNgUf"
SOURCE_KIND = "geometry_derived_simplified_proxy"
SOURCE_LABEL_ZH = "基于原始高模几何生成的简化图纸表达"
BLACK = "#151a20"
GREY = "#a3abb3"
EXPECTED_PATH_COUNTS = {"plan": 5, "front": 9, "side": 3}
VIEW_DEFINITIONS = {
    "plan": {"drawing_name": "WD03-WARDROBE-PLAN", "target_view": "PLAN_VIEW", "location_hint": "PLAN"},
    "front": {"drawing_name": "WD03-WARDROBE-FRONT", "target_view": "ELEVATION_VIEW", "location_hint": "SOUTH"},
    "side": {"drawing_name": "WD03-WARDROBE-SIDE", "target_view": "ELEVATION_VIEW", "location_hint": "EAST"},
}
GEOMETRY_TAGS = {"path", "polyline", "polygon", "line", "circle", "ellipse", "rect"}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def local_name(value: str) -> str:
    return value.rsplit("}", 1)[-1].split(":")[-1]


def coordinates_mm(view, first, second):
    if view == "plan":
        return float(first), float(second), 0.0
    if view == "front":
        return float(first), 0.0, float(second)
    return 0.0, float(first), float(second)


def curve_representation(model, context, view, paths):
    polylines = []
    for path in paths:
        points = [
            model.create_entity("IfcCartesianPoint", Coordinates=coordinates_mm(view, *point))
            for point in path
        ]
        if len(points) >= 2:
            polylines.append(model.create_entity("IfcPolyline", Points=points))
    if len(polylines) != EXPECTED_PATH_COUNTS[view]:
        raise RuntimeError(f"WD03 {view} path count drifted")
    curve_set = model.create_entity("IfcGeometricCurveSet", Elements=polylines)
    colour = model.create_entity(
        "IfcColourRgb", Name="WD03 approved semantic black", Red=21 / 255, Green=26 / 255, Blue=32 / 255
    )
    style = model.create_entity(
        "IfcCurveStyle",
        Name="WD03 approved geometry-derived semantic linework",
        CurveFont=None,
        CurveWidth=model.create_entity("IfcPositiveLengthMeasure", 0.35),
        CurveColour=colour,
        ModelOrDraughting=True,
    )
    model.create_entity("IfcStyledItem", Item=curve_set, Styles=[style], Name="Black semantic LINEWORK")
    return model.create_entity(
        "IfcShapeRepresentation",
        ContextOfItems=context,
        RepresentationIdentifier="Annotation",
        RepresentationType="GeometricCurveSet",
        Items=[curve_set],
    )


def add_semantic_annotation(model, drawing, target, target_obj, view, paths):
    context = tool.Drawing.get_annotation_context(VIEW_DEFINITIONS[view]["target_view"])
    if context is None:
        context = tool.Drawing.create_annotation_context(VIEW_DEFINITIONS[view]["target_view"])
    obj = core_drawing.add_annotation(
        tool.Ifc, tool.Collector, tool.Drawing,
        drawing=drawing, object_type="LINEWORK", relating_type=None, enable_editing=False,
    )
    annotation = tool.Ifc.get_entity(obj)
    annotation.Name = f"WD03 approved semantic black linework / {view}"
    annotation.Description = f"{SOURCE_LABEL_ZH}; representative {TARGET_GLOBAL_ID}; no official CAD geometry."
    annotation.ObjectPlacement = target.ObjectPlacement
    obj.matrix_world = target_obj.matrix_world
    vertices, edges = [], []
    for path in paths:
        indices = []
        for first, second in path:
            indices.append(len(vertices))
            vertices.append(tuple(value / 1000.0 for value in coordinates_mm(view, first, second)))
        edges.extend(zip(indices, indices[1:]))
    obj.data.clear_geometry()
    obj.data.from_pydata(vertices, edges, [])
    obj.data.update()
    annotation.Representation = model.create_entity(
        "IfcProductDefinitionShape", Representations=[curve_representation(model, context, view, paths)]
    )
    pset = ifcopenshell.api.pset.add_pset(model, product=annotation, name="EPset_Annotation")
    ifcopenshell.api.pset.edit_pset(
        model,
        pset=pset,
        properties={
            "Classes": "review-target-wd03 geometry-derived semantic-linework approved",
            "TargetGlobalId": TARGET_GLOBAL_ID,
            "SourceKind": SOURCE_KIND,
            "SourceLabelZh": SOURCE_LABEL_ZH,
            "OfficialCadUsed": False,
            "LineRole": "approved-geometry-derived-semantic-linework",
            "IfcCurveStyleColour": BLACK,
        },
    )
    return annotation


def view3d_override():
    for window in bpy.context.window_manager.windows:
        for area in window.screen.areas:
            if area.type == "VIEW_3D":
                region = next((item for item in area.regions if item.type == "WINDOW"), None)
                if region:
                    return {"window": window, "screen": window.screen, "area": area, "region": region, "scene": bpy.context.scene}
    raise RuntimeError("Bonsai Drawing requires a real VIEW_3D area")


def element_matches(element, guid):
    attrs = {local_name(key): value for key, value in element.attrib.items()}
    classes = element.attrib.get("class", "").split()
    return attrs.get("guid") == guid or f"GlobalId-{guid}" in classes


def style_and_inspect(svg_path: Path, annotation_guid: str):
    raw_sha = sha256(svg_path)
    ET.register_namespace("", "http://www.w3.org/2000/svg")
    ET.register_namespace("ifc", "http://www.ifcopenshell.org/ns")
    tree = ET.parse(svg_path)
    root = tree.getroot()
    groups = [element for element in root.iter() if element_matches(element, annotation_guid)]
    if not groups:
        raise RuntimeError("WD03 semantic annotation group missing from SVG")
    semantic_geometry = {id(item) for group in groups for item in group.iter() if local_name(item.tag) in GEOMETRY_TAGS}
    black_count = grey_count = geometry_count = projection_count = 0
    for element in root.iter():
        classes = element.attrib.get("class", "").split()
        projection_count += int("projection" in classes)
        if local_name(element.tag) not in GEOMETRY_TAGS:
            continue
        geometry_count += 1
        is_black = id(element) in semantic_geometry or element_matches(element, annotation_guid)
        colour = f"stroke:{BLACK};stroke-width:0.35;fill:none" if is_black else f"stroke:{GREY};stroke-width:0.20;fill:none;stroke-opacity:0.70"
        current = element.attrib.get("style", "").rstrip(";")
        element.attrib["style"] = f"{current};{colour}" if current else colour
        black_count += int(is_black)
        grey_count += int(not is_black)
    target_present = any(element_matches(element, TARGET_GLOBAL_ID) for element in root.iter())
    root.attrib.update({
        "data-create-drawing-result": "FINISHED",
        "data-source-kind": SOURCE_KIND,
        "data-source-label-zh": SOURCE_LABEL_ZH,
        "data-official-cad-used": "false",
        "data-target-body-projection-suppressed": str(not target_present).lower(),
    })
    tree.write(svg_path, encoding="utf-8", xml_declaration=True)
    result = {
        "root_tag": local_name(root.tag),
        "width": root.attrib.get("width"),
        "height": root.attrib.get("height"),
        "view_box": root.attrib.get("viewBox"),
        "geometry_element_count": geometry_count,
        "projection_group_count": projection_count,
        "semantic_black_geometry_count": black_count,
        "grey_context_geometry_count": grey_count,
        "representative_body_present": target_present,
        "target_body_projection_count": int(target_present),
        "semantic_annotation_present": black_count > 0,
        "no_body_annotation_duplicate": not target_present and black_count > 0,
        "bonsai_generated_sha256_before_review_style": raw_sha,
        "post_style_only": True,
    }
    if result["root_tag"] != "svg" or geometry_count == 0 or black_count == 0 or grey_count == 0 or not result["no_body_annotation_duplicate"]:
        raise RuntimeError(f"WD03 SVG gate failed: {result}")
    return result


def persisted_path_count(annotation):
    return sum(
        len(item.Elements)
        for representation in annotation.Representation.Representations
        for item in representation.Items
        if item.is_a("IfcGeometricCurveSet")
    )


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    if sha256(FORMAL_IFC) != FORMAL_SHA256:
        raise RuntimeError("formal IFC hash mismatch")
    approval = json.loads(APPROVAL.read_text(encoding="utf-8"))
    candidate = json.loads(CANDIDATE.read_text(encoding="utf-8"))
    report = json.loads(WRITE_REPORT.read_text(encoding="utf-8"))
    if (
        approval.get("status") != "approved"
        or approval.get("derived_ifc_write_allowed") is not True
        or approval.get("formal_authoritative_ifc_write_allowed") is not False
        or candidate.get("source_kind") != SOURCE_KIND
        or candidate.get("representative_global_id") != TARGET_GLOBAL_ID
        or candidate.get("official_cad_used") is not False
        or report.get("representation_path_counts") != EXPECTED_PATH_COUNTS
        or report.get("pass") is not True
    ):
        raise RuntimeError("WD03 approval/write gate failed")
    for view in VIEW_DEFINITIONS:
        if len(candidate["views"][view]["proxy_paths_mm"]) != EXPECTED_PATH_COUNTS[view]:
            raise RuntimeError(f"WD03 {view} semantic candidate drifted")

    derived_before = sha256(DERIVED_IFC)
    os.chdir(ROOT)
    with contextlib.suppress(Exception):
        addon_utils.disable("bl_ext.user_default.project_control", default_set=False, handle_error=None)
    load_result = bpy.ops.bim.load_project(filepath=str(DERIVED_IFC), should_start_fresh_session=True, use_detailed_tooltip=True)
    if load_result != {"FINISHED"} or not tool.Ifc.get():
        raise RuntimeError(f"Bonsai failed to load WD03 derived IFC: {load_result}")
    model = tool.Ifc.get()
    target = model.by_guid(TARGET_GLOBAL_ID)
    target_obj = tool.Ifc.get_object(target) if target else None
    if target is None or target_obj is None:
        raise RuntimeError("WD03 target missing")
    target_bbox = shared.world_bbox(target_obj)
    centre = [(target_bbox[0][axis] + target_bbox[1][axis]) / 2 for axis in range(3)]
    scene_bbox = ((centre[0] - 2.2, centre[1] - 2.2, -0.2), (centre[0] + 2.2, centre[1] + 2.2, 3.0))
    records = shared.room_elements(scene_bbox)
    context_elements = [element for element, *_ in records if element.GlobalId != TARGET_GLOBAL_ID]
    context_counts = Counter(element.is_a() for element in context_elements)
    if len(context_elements) < 5 or not any(name.startswith("IfcWall") for name in context_counts):
        raise RuntimeError("WD03 project context crop is incomplete")

    override = view3d_override()
    view_records = []
    for view, definition in VIEW_DEFINITIONS.items():
        output_svg = OUTPUT_DIR / f"{definition['drawing_name']}.svg"
        drawing, camera, width, height, clip_end = shared.add_drawing(model, definition, scene_bbox, context_elements, output_svg)
        with bpy.context.temp_override(**override):
            activated = bpy.ops.bim.activate_drawing(drawing=drawing.id(), should_view_from_camera=False)
        if activated != {"FINISHED"}:
            raise RuntimeError(f"WD03 {view} activation failed")
        annotation = add_semantic_annotation(model, drawing, target, target_obj, view, candidate["views"][view]["proxy_paths_mm"])
        cprops = tool.Drawing.get_camera_props(camera)
        cprops.has_annotation = True
        drawing_pset = ifcopenshell.util.element.get_pset(drawing, "EPset_Drawing")
        ifcopenshell.api.pset.edit_pset(model, pset=model.by_id(drawing_pset["id"]), properties={"HasAnnotation": True})
        dprops = tool.Drawing.get_document_props()
        dprops.should_use_underlay_cache = False
        dprops.should_use_linework_cache = False
        dprops.should_use_annotation_cache = False
        with bpy.context.temp_override(**override):
            create_result = bpy.ops.bim.create_drawing(print_all=False, open_viewer=False, sync=False)
        if create_result != {"FINISHED"} or not output_svg.is_file() or output_svg.stat().st_size == 0:
            raise RuntimeError(f"WD03 {view} Create Drawing failed: {create_result}")
        svg = style_and_inspect(output_svg, annotation.GlobalId)
        cache = output_svg.parent / "cache" / f"{output_svg.stem}-linework.svg"
        if not cache.is_file() or cache.stat().st_size == 0:
            raise RuntimeError(f"WD03 {view} linework cache missing")
        view_records.append({
            "view": view,
            "drawing_global_id": drawing.GlobalId,
            "annotation_global_id": annotation.GlobalId,
            "path_count": EXPECTED_PATH_COUNTS[view],
            "output_svg": output_svg,
            "cache": cache,
            "svg": svg,
            "create_result": sorted(create_result),
            "camera": {"type": camera.data.type, "width_m": width, "height_m": height, "clip_end_m": clip_end, "matrix_world": [list(row) for row in camera.matrix_world]},
        })

    temporary = DERIVED_IFC.with_suffix(".ifc.drawings-next")
    model.write(temporary)
    os.replace(temporary, DERIVED_IFC)
    reload_result = bpy.ops.bim.load_project(filepath=str(DERIVED_IFC), should_start_fresh_session=True, use_detailed_tooltip=True)
    if reload_result != {"FINISHED"} or not tool.Ifc.get():
        raise RuntimeError("WD03 persisted derived IFC reload failed")
    reopened = tool.Ifc.get()
    for record in view_records:
        drawing = reopened.by_guid(record["drawing_global_id"])
        annotation = reopened.by_guid(record["annotation_global_id"])
        if drawing is None or annotation is None or persisted_path_count(annotation) != record["path_count"]:
            raise RuntimeError(f"WD03 persisted {record['view']} Drawing/LINEWORK missing")
        record["persisted_path_count"] = persisted_path_count(annotation)
    bpy.ops.wm.save_as_mainfile(filepath=str(SESSION_BLEND), check_existing=False)
    if sha256(FORMAL_IFC) != FORMAL_SHA256:
        raise RuntimeError("formal IFC changed during WD03 Drawing workflow")

    evidence = {
        "schema_version": 1,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "task": "approved WD03 product-level derived IFC plus Bonsai Create Drawing Plan/Front/Side in actual project context",
        "courseEvidence": {
            "mode": "embedded-course-index",
            "lesson": "085000 Introduction to Drawings",
            "timestamps": ["01:59 Create Drawing", "02:06 Drawing created", "02:13 SVG inspection", "02:41 Element Filters"],
            "course_fact": "Drawing camera, depth and element filters precede Create Drawing; generated SVG is then inspected.",
            "current_inference": "The workflow was adapted to Blender 4.5.3 LTS and Bonsai/IfcOpenShell 0.8.4 using semantic Drawing state and bpy.ops.bim.create_drawing.",
        },
        "preState": {"derived_ifc_sha256": derived_before, "formal_ifc_sha256": FORMAL_SHA256, "target_global_id": TARGET_GLOBAL_ID},
        "execution": {"generator": "bpy.ops.bim.create_drawing", "linework_mode": "OPENCASCADE", "source_kind": SOURCE_KIND, "source_label_zh": SOURCE_LABEL_ZH, "official_cad_used": False},
        "persistence": {"method": "temporary IFC write, atomic replace, Bonsai reload and Blend save", "reload_result": sorted(reload_result)},
        "postState": {"derived_ifc_sha256": sha256(DERIVED_IFC), "formal_ifc_sha256": sha256(FORMAL_IFC), "drawing_count": 3},
        "outputs": {"views": [
            {
                "view": record["view"], "drawing_global_id": record["drawing_global_id"], "annotation_global_id": record["annotation_global_id"],
                "path_count": record["path_count"], "persisted_path_count": record["persisted_path_count"], "camera": record["camera"],
                "create_drawing": {"operator": "bpy.ops.bim.create_drawing", "result": record["create_result"]},
                "svg": {"path": str(record["output_svg"]), "bytes": record["output_svg"].stat().st_size, "sha256": sha256(record["output_svg"]), **record["svg"]},
                "linework_cache": {"path": str(record["cache"]), "bytes": record["cache"].stat().st_size, "sha256": sha256(record["cache"])},
            } for record in view_records
        ], "session_blend": {"path": str(SESSION_BLEND), "bytes": SESSION_BLEND.stat().st_size, "sha256": sha256(SESSION_BLEND)}},
        "visual": {"representative_body_displayed": False, "semantic_annotation_displayed_once_per_view": True, "context_colour": GREY, "semantic_colour": BLACK},
        "context": {"include_count": len(context_elements), "include_ifc_class_counts": dict(sorted(context_counts.items())), "representative_body_suppressed_from_drawing_include": True},
        "tests": {"all_create_drawing_finished": True, "all_body_annotation_duplicates_absent": all(record["svg"]["no_body_annotation_duplicate"] for record in view_records), "all_path_counts_persisted": all(record["path_count"] == record["persisted_path_count"] for record in view_records), "formal_ifc_unchanged": True},
        "verdict": "pass",
        "pass": True,
    }
    EVIDENCE.write_text(json.dumps(evidence, indent=2, ensure_ascii=False, default=str) + "\n", encoding="utf-8")
    print(json.dumps({"evidence": str(EVIDENCE), "derived_ifc_sha256": sha256(DERIVED_IFC), "pass": True}, indent=2))


if __name__ == "__main__":
    try:
        main()
    except Exception:
        PRODUCT_DIR.mkdir(parents=True, exist_ok=True)
        (PRODUCT_DIR / "wd03-create-drawing-error.log").write_text(traceback.format_exc(), encoding="utf-8")
        raise
    finally:
        if not bpy.app.background:
            def quit_after_handlers():
                bpy.ops.wm.quit_blender()
                return None
            bpy.app.timers.register(quit_after_handlers, first_interval=1.0)
