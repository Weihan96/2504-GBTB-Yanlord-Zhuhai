#!/usr/bin/env python3
"""Create approved Gessi316 54146 Drawings in the actual main bathroom."""

from __future__ import annotations

import contextlib
import hashlib
import json
import os
import sys
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


PRODUCT_DIR = ROOT / "output/review/highpoly-types/gessi316-54146"
FORMAL_IFC = ROOT / "2504 GBTB Yanlord Zhuhai.ifc"
FORMAL_SHA256 = "7a50b87e8f48a7c2bfbaa9f1dabdfca7155fa684b2325c66aa1ae15c8ab0a25c"
APPROVAL = ROOT / "pipeline/decisions/gessi316-54146-drawing-approval.json"
CANDIDATE = PRODUCT_DIR / "candidate-representations.json"
TARGET_GLOBAL_ID = "04DLh1Jk9Dcu9ibcaE0id8"
ROOM_GLOBAL_ID = "3a4COIs5X7lgDirMBDT4Vs"
ROOM_NAME = "主卫湿区"
SOURCE_KIND = "native_dwg_review_simplification"
SOURCE_LABEL_ZH = "基于官方 Gessi 54146 G000 原生 DWG 轮廓的简化蓝线审核表达"
SOURCE_DWG_SHA256 = "c8ddb90f61565d5273a32574f777812bbb1d9ef33e2b98479f1e7c381e69950d"
BLUE = "#1677c8"
GREY = "#a3abb3"
EXPECTED_PATH_COUNTS = {"plan": 16, "front": 104, "side": 116}
VIEW_DEFINITIONS = {
    "plan": {"drawing_name": "GESSI316-54146-MAIN-BATH-PLAN", "target_view": "PLAN_VIEW", "location_hint": "PLAN"},
    "front": {"drawing_name": "GESSI316-54146-MAIN-BATH-FRONT", "target_view": "ELEVATION_VIEW", "location_hint": "EAST"},
    "side": {"drawing_name": "GESSI316-54146-MAIN-BATH-SIDE", "target_view": "ELEVATION_VIEW", "location_hint": "SOUTH"},
}
GEOMETRY_TAGS = {"path", "polyline", "polygon", "line", "circle", "ellipse", "rect"}


shared.PRODUCT_DIR = PRODUCT_DIR
shared.FORMAL_IFC = FORMAL_IFC
shared.FORMAL_SHA256 = FORMAL_SHA256
shared.TARGET_GLOBAL_ID = TARGET_GLOBAL_ID
shared.ROOM_GLOBAL_ID = ROOM_GLOBAL_ID
shared.ROOM_NAME = ROOM_NAME
shared.ARTICLE = "54146 G000"
shared.DRAWING_PRODUCT_CODE = "54146"
shared.SOURCE_DWG_SHA256 = SOURCE_DWG_SHA256
shared.SOURCE_KIND = SOURCE_KIND
shared.SOURCE_LABEL_ZH = SOURCE_LABEL_ZH
shared.BLUE = BLUE
shared.EXPECTED_PATH_COUNTS = EXPECTED_PATH_COUNTS
shared.VIEW_DEFINITIONS = VIEW_DEFINITIONS


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def local_name(value: str) -> str:
    return value.rsplit("}", 1)[-1].split(":")[-1]


def coordinates_mm(view: str, first: float, second: float):
    if view == "plan":
        return float(first), float(second), 0.0
    if view == "front":
        return float(first), 0.0, float(second)
    return 0.0, float(first), float(second)


def bounds_3d(points):
    return [
        [min(point[axis] for point in points) for axis in range(3)],
        [max(point[axis] for point in points) for axis in range(3)],
    ]


def annotation_representation(model, context, view, paths):
    polylines = []
    for path in paths:
        points = [
            model.create_entity("IfcCartesianPoint", Coordinates=coordinates_mm(view, first, second))
            for first, second in path
        ]
        if len(points) >= 2:
            polylines.append(model.create_entity("IfcPolyline", Points=points))
    if len(polylines) != EXPECTED_PATH_COUNTS[view]:
        raise RuntimeError(f"Gessi316 54146 {view} review path count drifted")
    curve_set = model.create_entity("IfcGeometricCurveSet", Elements=polylines)
    colour = model.create_entity(
        "IfcColourRgb", Name="Gessi 54146 review blue",
        Red=22.0 / 255.0, Green=119.0 / 255.0, Blue=200.0 / 255.0,
    )
    style = model.create_entity(
        "IfcCurveStyle", Name="Gessi 54146 official-outline review simplification",
        CurveFont=None, CurveWidth=model.create_entity("IfcPositiveLengthMeasure", 0.35),
        CurveColour=colour, ModelOrDraughting=True,
    )
    model.create_entity("IfcStyledItem", Item=curve_set, Styles=[style], Name="Blue review LINEWORK")
    representation = model.create_entity(
        "IfcShapeRepresentation", ContextOfItems=context,
        RepresentationIdentifier="Annotation", RepresentationType="GeometricCurveSet",
        Items=[curve_set],
    )
    return representation, len(polylines), sum(max(0, len(path) - 1) for path in paths)


def add_review_annotation(model, drawing, target, target_obj, view, paths):
    target_view = VIEW_DEFINITIONS[view]["target_view"]
    context = tool.Drawing.get_annotation_context(target_view) or tool.Drawing.create_annotation_context(target_view)
    obj = core_drawing.add_annotation(
        tool.Ifc, tool.Collector, tool.Drawing, drawing=drawing,
        object_type="LINEWORK", relating_type=None, enable_editing=False,
    )
    annotation = tool.Ifc.get_entity(obj)
    annotation.Name = f"Gessi 54146 official-outline review simplification / {view}"
    annotation.Description = (
        f"{SOURCE_LABEL_ZH}; target {TARGET_GLOBAL_ID}; original official DWG remains separate evidence; "
        "fine spray-nozzle detail intentionally removed."
    )
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
    representation, path_count, edge_count = annotation_representation(model, context, view, paths)
    annotation.Representation = model.create_entity("IfcProductDefinitionShape", Representations=[representation])
    pset = ifcopenshell.api.pset.add_pset(model, product=annotation, name="EPset_Annotation")
    ifcopenshell.api.pset.edit_pset(
        model, pset=pset,
        properties={
            "Classes": "review-target-gessi54146 native-dwg-review-simplification",
            "TargetGlobalId": TARGET_GLOBAL_ID,
            "ArticleNumber": "54146",
            "Configuration": "G000",
            "SourceKind": SOURCE_KIND,
            "SourceLabelZh": SOURCE_LABEL_ZH,
            "SourceDwgSha256": SOURCE_DWG_SHA256,
            "UnalteredOfficialDwg": False,
            "OriginalOfficialDwgEvidencePreserved": True,
            "FineSprayNozzleDetailPathCount": 0,
            "SourceScale": 1.0,
            "SideVerticalReflectionApplied": view == "side",
            "IfcCurveStyleColour": BLUE,
        },
    )
    return annotation, path_count, edge_count


def style_svg(svg_path: Path, annotation_guid: str, view: str):
    raw_sha = sha256(svg_path)
    ET.register_namespace("", "http://www.w3.org/2000/svg")
    ET.register_namespace("ifc", "http://www.ifcopenshell.org/ns")
    tree = ET.parse(svg_path)
    root = tree.getroot()
    groups = []
    for element in root.iter():
        attrs = {local_name(key): value for key, value in element.attrib.items()}
        classes = element.attrib.get("class", "").split()
        if attrs.get("guid") == annotation_guid or f"GlobalId-{annotation_guid}" in classes or "review-target-gessi54146" in classes:
            groups.append(element)
    if not groups:
        raise RuntimeError("Bonsai SVG has no Gessi316 54146 review Annotation")
    blue_ids = {
        id(element) for group in groups for element in group.iter()
        if local_name(element.tag) in GEOMETRY_TAGS
    }
    blue_count = grey_count = 0
    for element in root.iter():
        if local_name(element.tag) not in GEOMETRY_TAGS:
            continue
        existing = element.attrib.get("style", "").rstrip(";")
        if id(element) in blue_ids:
            colour = f"stroke:{BLUE};stroke-width:0.35;fill:none"
            blue_count += 1
        else:
            colour = f"stroke:{GREY};stroke-width:0.22;fill:none;stroke-opacity:0.72"
            grey_count += 1
        element.attrib["style"] = f"{existing};{colour}" if existing else colour
    for group in groups:
        classes = group.attrib.get("class", "").split()
        for name in ("review-target-gessi54146", "native-dwg-review-simplification"):
            if name not in classes:
                classes.append(name)
        group.attrib.update({
            "class": " ".join(classes), "data-source-kind": SOURCE_KIND,
            "data-source-dwg-sha256": SOURCE_DWG_SHA256,
            "data-unaltered-official-dwg": "false",
            "data-fine-spray-nozzle-detail-path-count": "0",
            "data-source-scale": "1.0", "data-reviewed-view": view,
        })
    root.attrib.update({
        "data-create-drawing-result": "FINISHED", "data-context-colour": GREY,
        "data-review-linework-colour": BLUE, "data-target-body-projection-suppressed": "true",
    })
    tree.write(svg_path, encoding="utf-8", xml_declaration=True)
    if blue_count == 0 or grey_count == 0:
        raise RuntimeError("Gessi316 54146 SVG colour-layer gate failed")
    return {
        "bonsai_generated_sha256_before_review_style": raw_sha,
        "blue_geometry_element_count": blue_count,
        "grey_context_geometry_element_count": grey_count,
        "post_style_only": True,
    }


def persisted_path_count(annotation):
    return sum(
        len(item.Elements)
        for representation in annotation.Representation.Representations
        for item in representation.Items
        if item.is_a("IfcGeometricCurveSet")
    )


def view3d_override():
    for window in bpy.context.window_manager.windows:
        for area in window.screen.areas:
            if area.type != "VIEW_3D":
                continue
            region = next((item for item in area.regions if item.type == "WINDOW"), None)
            if region:
                return {"window": window, "screen": window.screen, "area": area, "region": region, "scene": bpy.context.scene}
    raise RuntimeError("Bonsai Drawing requires a real VIEW_3D area")


def main():
    arguments = sys.argv[sys.argv.index("--") + 1 :]
    if len(arguments) != 3:
        raise SystemExit("expected: derived.ifc output-directory session.blend")
    derived_ifc = Path(arguments[0]).resolve()
    output_dir = Path(arguments[1]).resolve()
    session_blend = Path(arguments[2]).resolve()
    if sha256(FORMAL_IFC) != FORMAL_SHA256:
        raise RuntimeError("formal IFC hash mismatch before Gessi316 54146 Drawing workflow")
    approval = json.loads(APPROVAL.read_text(encoding="utf-8"))
    candidate = json.loads(CANDIDATE.read_text(encoding="utf-8"))
    if (
        approval.get("status") != "approved"
        or approval.get("approved_views") != ["plan", "front", "side"]
        or approval.get("derived_ifc_write_allowed") is not True
        or approval.get("formal_authoritative_ifc_write_allowed") is not False
        or candidate.get("source_kind") != SOURCE_KIND
        or candidate.get("source_dwg_sha256") != SOURCE_DWG_SHA256
        or candidate.get("unaltered_official_cad_used_as_review_representation") is not False
    ):
        raise RuntimeError("Gessi316 54146 product-level approval/source gate failed")
    paths_by_view = {}
    for view in VIEW_DEFINITIONS:
        item = candidate["views"][view]
        paths = item["review_simplified_official_outline_paths_mm"]
        audit = item["handle_line_texture_simplification"]
        if (
            len(paths) != EXPECTED_PATH_COUNTS[view]
            or item.get("unaltered_official_dwg") is not False
            or audit.get("envelope_delta_mm") != [0.0, 0.0]
            or audit.get("centre_delta_mm") != [0.0, 0.0]
            or audit.get("pass") is not True
        ):
            raise RuntimeError(f"Gessi316 54146 {view} review line gate failed")
        paths_by_view[view] = paths
    derived_before_sha = sha256(derived_ifc)
    os.chdir(ROOT)
    with contextlib.suppress(Exception):
        addon_utils.disable("bl_ext.user_default.project_control", default_set=False, handle_error=None)
    load_result = bpy.ops.bim.load_project(filepath=str(derived_ifc), should_start_fresh_session=True, use_detailed_tooltip=True)
    if load_result != {"FINISHED"} or not tool.Ifc.get():
        raise RuntimeError(f"Bonsai failed to load Gessi316 54146 derived IFC: {load_result}")
    model = tool.Ifc.get()
    existing_drawing_count = sum(
        1 for item in model.by_type("IfcAnnotation") if item.ObjectType == "DRAWING"
    )
    for definition in VIEW_DEFINITIONS.values():
        if any(item.Name == definition["drawing_name"] for item in model.by_type("IfcAnnotation")):
            raise RuntimeError(f"Gessi316 54146 Drawing already exists: {definition['drawing_name']}")
    target = model.by_guid(TARGET_GLOBAL_ID)
    room = model.by_guid(ROOM_GLOBAL_ID)
    target_obj = tool.Ifc.get_object(target) if target else None
    room_obj = tool.Ifc.get_object(room) if room else None
    if target is None or room is None or target_obj is None or room_obj is None or room.LongName != ROOM_NAME:
        raise RuntimeError("Gessi316 54146 target or main bathroom identity missing")
    room_bbox = shared.world_bbox(room_obj)
    target_bbox = shared.world_bbox(target_obj)
    records = shared.room_elements(room_bbox)
    original_context = [record[0] for record in records]
    if original_context.count(target) != 1:
        raise RuntimeError("expected Gessi316 54146 target exactly once in main bathroom")
    context_elements = [element for element in original_context if element != target]
    if len(context_elements) < 10:
        raise RuntimeError("main bathroom project context is unexpectedly empty")
    context_counts = Counter(element.is_a() for element in context_elements)
    output_dir.mkdir(parents=True, exist_ok=True)
    override = view3d_override()
    view_records = []
    for view, definition in VIEW_DEFINITIONS.items():
        output_svg = output_dir / f"{definition['drawing_name']}.svg"
        drawing, camera, width, height, clip_end = shared.add_drawing(model, definition, room_bbox, context_elements, output_svg)
        if view == "plan":
            camera.matrix_world.translation.z = room_bbox[1][2] - 0.05
            camera.data.clip_end = room_bbox[1][2] - room_bbox[0][2] + 0.40
            clip_end = camera.data.clip_end
        with bpy.context.temp_override(**override):
            activate_result = bpy.ops.bim.activate_drawing(drawing=drawing.id(), should_view_from_camera=False)
        if activate_result != {"FINISHED"}:
            raise RuntimeError(f"Gessi316 54146 {view} Drawing activation failed")
        annotation, path_count, edge_count = add_review_annotation(model, drawing, target, target_obj, view, paths_by_view[view])
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
            raise RuntimeError(f"Gessi316 54146 Bonsai Create Drawing failed for {view}: {create_result}")
        style = style_svg(output_svg, annotation.GlobalId, view)
        svg = shared.inspect_svg(output_svg, TARGET_GLOBAL_ID, annotation.GlobalId)
        if (
            svg["root_tag"] != "svg" or svg["geometry_element_count"] == 0
            or svg["projection_group_count"] == 0
            or svg["target_ifc_projection_group_count"] != 0
            or svg["review_annotation_geometry_count"] == 0
        ):
            raise RuntimeError(f"Gessi316 54146 {view} SVG structure/duplicate gate failed: {svg}")
        cache = output_svg.parent / "cache" / f"{output_svg.stem}-linework.svg"
        if not cache.is_file() or cache.stat().st_size == 0:
            raise RuntimeError(f"Gessi316 54146 {view} linework cache missing")
        view_records.append({
            "view": view, "drawing_global_id": drawing.GlobalId, "drawing_name": drawing.Name,
            "annotation_global_id": annotation.GlobalId, "path_count": path_count, "edge_count": edge_count,
            "expected_points": [coordinates_mm(view, first, second) for path in paths_by_view[view] for first, second in path],
            "output_svg": output_svg, "cache": cache, "style": style, "svg": svg,
            "camera": {"type": camera.data.type, "width_m": width, "height_m": height,
                       "clip_start_m": camera.data.clip_start, "clip_end_m": clip_end,
                       "matrix_world": [list(row) for row in camera.matrix_world],
                       "resolution": [bpy.context.scene.render.resolution_x, bpy.context.scene.render.resolution_y]},
            "create_result": sorted(create_result),
        })
    temporary = derived_ifc.with_suffix(".ifc.drawings-next")
    model.write(str(temporary))
    os.replace(temporary, derived_ifc)
    derived_after_sha = sha256(derived_ifc)
    reload_result = bpy.ops.bim.load_project(filepath=str(derived_ifc), should_start_fresh_session=True, use_detailed_tooltip=True)
    if reload_result != {"FINISHED"} or not tool.Ifc.get():
        raise RuntimeError("persisted Gessi316 54146 IFC reload failed")
    reopened = tool.Ifc.get()
    output_views = []
    for record in view_records:
        view = record["view"]
        drawing = reopened.by_guid(record["drawing_global_id"])
        annotation = reopened.by_guid(record["annotation_global_id"])
        if drawing is None or annotation is None or persisted_path_count(annotation) != record["path_count"]:
            raise RuntimeError(f"reloaded Gessi316 54146 {view} Drawing/Annotation missing")
        representation = next(item for item in annotation.Representation.Representations if item.RepresentationIdentifier == "Annotation")
        points = [tuple(float(value) for value in point.Coordinates) for item in representation.Items for path in item.Elements for point in path.Points]
        residual = max(
            abs(bounds_3d(points)[bound][axis] - bounds_3d(record["expected_points"])[bound][axis])
            for bound in range(2) for axis in range(3)
        )
        annotation_pset = ifcopenshell.util.element.get_pset(annotation, "EPset_Annotation")
        drawing_pset = ifcopenshell.util.element.get_pset(drawing, "EPset_Drawing")
        include = drawing_pset.get("Include", "").split(",")
        assignments = [rel for rel in reopened.by_type("IfcRelAssignsToGroup") if drawing in rel.RelatedObjects and annotation in rel.RelatedObjects]
        if (
            residual > 0.000001 or TARGET_GLOBAL_ID in include or len(include) != len(context_elements)
            or annotation_pset.get("UnalteredOfficialDwg") is not False
            or annotation_pset.get("FineSprayNozzleDetailPathCount") != 0
            or annotation_pset.get("SourceScale") != 1.0 or not assignments
        ):
            raise RuntimeError(f"persisted Gessi316 54146 {view} duplicate/metadata gate failed")
        output_views.append({
            "view": view, "drawing": {"global_id": drawing.GlobalId, "name": drawing.Name, "epset_drawing": drawing_pset},
            "linework_annotation": {"global_id": annotation.GlobalId, "name": annotation.Name, "epset_annotation": annotation_pset,
                                    "drawing_group_assignment_global_id": assignments[0].GlobalId},
            "review_path_count": record["path_count"], "persisted_review_path_count": persisted_path_count(annotation),
            "review_edge_count": record["edge_count"], "persisted_coordinate_residual_mm": residual,
            "camera": record["camera"],
            "create_drawing": {"operator": "bpy.ops.bim.create_drawing", "arguments": {"print_all": False, "open_viewer": False, "sync": False},
                               "result": record["create_result"], "linework_mode": "OPENCASCADE", "target_view": VIEW_DEFINITIONS[view]["target_view"]},
            "svg": {"path": str(record["output_svg"]), "bytes": record["output_svg"].stat().st_size,
                    "sha256": sha256(record["output_svg"]), **record["style"], **record["svg"]},
            "linework_cache": {"path": str(record["cache"]), "bytes": record["cache"].stat().st_size, "sha256": sha256(record["cache"])},
        })
    bpy.ops.wm.save_as_mainfile(filepath=str(session_blend), check_existing=False)
    if sha256(FORMAL_IFC) != FORMAL_SHA256:
        raise RuntimeError("formal IFC changed during Gessi316 54146 Drawing workflow")
    evidence = {
        "schema_version": 1, "generated_at": datetime.now(timezone.utc).isoformat(),
        "task": "approved Gessi316 54146 Bonsai Create Drawing Plan/Front/Side in actual main-bathroom context",
        "workflow": ["inspect", "plan", "execute", "persist", "reload", "verify"],
        "provider": {"preferred": "bonsai-mcp 1.1.0 at b0e67b1", "status": "unsupported",
                     "reason": "Blender bridge 127.0.0.1:9878 connection refused",
                     "fallback": "local Blender 4.5.3 LTS with installed Bonsai/IfcOpenShell 0.8.4"},
        "course_evidence": {"evidence_mode": "embedded-course-index", "lesson": "085000 Introduction to Drawings",
                            "timestamps": ["01:03 active Drawing camera", "01:59 Create Drawing", "02:13 SVG inspection"],
                            "course_fact": "Create Drawing refreshes SVG after the Drawing camera, linework, annotation and filters are configured.",
                            "current_inference": "The course workflow was adapted to Bonsai 0.8.4 semantic state and bpy.ops.bim.create_drawing."},
        "preState": {"derived_ifc": str(derived_ifc), "derived_ifc_sha256": derived_before_sha,
                     "formal_ifc_sha256": sha256(FORMAL_IFC), "drawing_count": existing_drawing_count},
        "execution": {"generator": "bpy.ops.bim.create_drawing", "linework_mode": "OPENCASCADE"},
        "persistence": {"method": "atomic IFC replace; bpy.ops.bim.load_project reload; save Blend",
                        "reload_result": sorted(reload_result), "post_reload_drawing_count": existing_drawing_count + len(output_views),
                        "post_reload_annotation_count": len(output_views)},
        "postState": {"derived_ifc": str(derived_ifc), "derived_ifc_sha256": derived_after_sha,
                      "formal_ifc_sha256": sha256(FORMAL_IFC), "drawing_count": existing_drawing_count + len(output_views)},
        "source_kind": SOURCE_KIND, "source_label_zh": SOURCE_LABEL_ZH,
        "source_dwg_sha256": SOURCE_DWG_SHA256, "unaltered_official_dwg": False,
        "original_official_dwg_evidence_preserved": True, "fine_spray_nozzle_detail_path_count": 0,
        "target": {"global_id": TARGET_GLOBAL_ID, "world_bbox_m": [list(target_bbox[0]), list(target_bbox[1])],
                   "target_include_count_before_suppression": original_context.count(target), "target_include_count_after_suppression": 0},
        "room": {"global_id": ROOM_GLOBAL_ID, "name": ROOM_NAME, "bbox_m": [list(room_bbox[0]), list(room_bbox[1])]},
        "context": {"project_context_retained": True, "include_count": len(context_elements),
                    "include_ifc_class_counts": dict(sorted(context_counts.items())), "target_body_suppressed": True},
        "views": output_views,
        "session_blend": {"path": str(session_blend), "bytes": session_blend.stat().st_size, "sha256": sha256(session_blend)},
        "versions": {"blender": bpy.app.version_string, "ifcopenshell": ifcopenshell.version,
                     "ifc_schema": reopened.schema, "bonsai_generator": "Bonsai 0.8.4 bim.create_drawing"},
        "formal_ifc_bytes_unchanged": True, "verdict": "pass", "pass": True,
    }
    evidence_path = output_dir / "GESSI316-54146-MAIN-BATH-create-drawing-evidence.json"
    evidence_path.write_text(json.dumps(evidence, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps({"derived_ifc": str(derived_ifc), "evidence": str(evidence_path), "pass": True}, ensure_ascii=False))


if __name__ == "__main__":
    main()
