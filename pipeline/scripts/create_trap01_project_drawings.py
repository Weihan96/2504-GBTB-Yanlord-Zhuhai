#!/usr/bin/env python3
"""Create approved TRAP01 current-configuration Drawings in actual project context.

Run inside the dedicated Blender 4.5.3 / Bonsai 0.8.4 session through the
public bonsai-mcp Provider. The formal IFC is read-only. The product-level
derived IFC persists three identifiable adjustable components plus the native
Bonsai Drawings and current-installed black LINEWORK annotations.
"""

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


PRODUCT_DIR = ROOT / "output/review/highpoly-types/trap01"
FORMAL_IFC = ROOT / "2504 GBTB Yanlord Zhuhai.ifc"
FORMAL_SHA256 = "7a50b87e8f48a7c2bfbaa9f1dabdfca7155fa684b2325c66aa1ae15c8ab0a25c"
DERIVED_IFC = PRODUCT_DIR / "Geberit-151.116.11.1-TRAP01-derived-drawing.ifc"
SESSION_BLEND = PRODUCT_DIR / "Geberit-151.116.11.1-TRAP01-project-drawings.blend"
OUTPUT_DIR = PRODUCT_DIR / "bonsai-drawings/project-context"
EVIDENCE = OUTPUT_DIR / "TRAP01-create-drawing-evidence.json"
APPROVAL = ROOT / "pipeline/decisions/trap01-drawing-approval.json"
CANDIDATE = PRODUCT_DIR / "candidate-representations.json"
WRITE_REPORT = PRODUCT_DIR / "Geberit-151.116.11.1-TRAP01-derived-drawing-report.json"
TARGET_GLOBAL_ID = "2Ak2ma0lvBEA49UpplzUqi"
SOURCE_KIND = "geometry_derived_simplified_proxy"
SOURCE_LABEL_ZH = "基于原始高模几何生成的简化图纸表达"
COMPONENT_PSET = "Pset_Trap01DrawingComponent"
BLACK = "#111820"
GREY = "#a3abb3"
BLUE = "#1677c8"
EXPECTED_PATH_COUNTS = {"plan": 16, "front": 1, "side": 8}
VIEW_DEFINITIONS = {
    "plan": {"drawing_name": "TRAP01-PROJECT-CONTEXT-PLAN", "target_view": "PLAN_VIEW", "location_hint": "PLAN"},
    "front": {"drawing_name": "TRAP01-PROJECT-CONTEXT-FRONT", "target_view": "ELEVATION_VIEW", "location_hint": "EAST"},
    "side": {"drawing_name": "TRAP01-PROJECT-CONTEXT-SIDE", "target_view": "ELEVATION_VIEW", "location_hint": "SOUTH"},
}
GEOMETRY_TAGS = {"path", "polyline", "polygon", "line", "circle", "ellipse", "rect"}


shared.ROOM_NAME = "TRAP01 representative installation context"
shared.ARTICLE = "151.116.11.1"
shared.SOURCE_LABEL_ZH = SOURCE_LABEL_ZH


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


def curve_representation(model, context, view, paths):
    polylines = []
    for path in paths:
        points = [model.create_entity("IfcCartesianPoint", Coordinates=coordinates_mm(view, *point)) for point in path]
        if len(points) >= 2:
            polylines.append(model.create_entity("IfcPolyline", Points=points))
    if len(polylines) != EXPECTED_PATH_COUNTS[view]:
        raise RuntimeError(f"TRAP01 {view} current path count drifted")
    curve_set = model.create_entity("IfcGeometricCurveSet", Elements=polylines)
    colour = model.create_entity(
        "IfcColourRgb", Name="TRAP01 current installed black", Red=17.0 / 255.0, Green=24.0 / 255.0, Blue=32.0 / 255.0
    )
    style = model.create_entity(
        "IfcCurveStyle", Name="TRAP01 current installed drawing linework", CurveFont=None,
        CurveWidth=model.create_entity("IfcPositiveLengthMeasure", 0.35), CurveColour=colour, ModelOrDraughting=True,
    )
    model.create_entity("IfcStyledItem", Item=curve_set, Styles=[style], Name="Current installed LINEWORK")
    return model.create_entity(
        "IfcShapeRepresentation", ContextOfItems=context, RepresentationIdentifier="Annotation",
        RepresentationType="GeometricCurveSet", Items=[curve_set],
    )


def add_current_annotation(model, drawing, target, target_obj, view, paths):
    target_view = VIEW_DEFINITIONS[view]["target_view"]
    context = tool.Drawing.get_annotation_context(target_view) or tool.Drawing.create_annotation_context(target_view)
    obj = core_drawing.add_annotation(
        tool.Ifc, tool.Collector, tool.Drawing, drawing=drawing,
        object_type="LINEWORK", relating_type=None, enable_editing=False,
    )
    annotation = tool.Ifc.get_entity(obj)
    annotation.Name = f"TRAP01 current installed configuration / {view}"
    annotation.Description = (
        f"{SOURCE_LABEL_ZH}; target {TARGET_GLOBAL_ID}; current installed black configuration only. "
        "Official blue dashed adjustable excess is intentionally excluded from scene Drawings."
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
    representation = curve_representation(model, context, view, paths)
    annotation.Representation = model.create_entity("IfcProductDefinitionShape", Representations=[representation])
    pset = ifcopenshell.api.pset.add_pset(model, product=annotation, name="EPset_Annotation")
    ifcopenshell.api.pset.edit_pset(
        model,
        pset=pset,
        properties={
            "Classes": "review-target-trap01 current-installed geometry-derived",
            "TargetGlobalId": TARGET_GLOBAL_ID,
            "SourceKind": SOURCE_KIND,
            "SourceLabelZh": SOURCE_LABEL_ZH,
            "LineRole": "current-installed",
            "IfcCurveStyleColour": BLACK,
            "ReferenceExtensionDisplayed": False,
            "BlueDashedReferenceExcluded": True,
        },
    )
    return annotation, len(paths), sum(max(0, len(path) - 1) for path in paths)


def style_svg(svg_path: Path, annotation_guid: str):
    raw_sha = sha256(svg_path)
    ET.register_namespace("", "http://www.w3.org/2000/svg")
    ET.register_namespace("ifc", "http://www.ifcopenshell.org/ns")
    tree = ET.parse(svg_path)
    root = tree.getroot()
    target_groups = []
    for element in root.iter():
        attrs = {local_name(key): value for key, value in element.attrib.items()}
        classes = element.attrib.get("class", "").split()
        if attrs.get("guid") == annotation_guid or f"GlobalId-{annotation_guid}" in classes:
            target_groups.append(element)
    if not target_groups:
        raise RuntimeError("TRAP01 current-installed annotation group missing from Bonsai SVG")
    black = grey = 0
    for element in root.iter():
        if local_name(element.tag) not in GEOMETRY_TAGS:
            continue
        attrs = {local_name(key): value for key, value in element.attrib.items()}
        classes = element.attrib.get("class", "").split()
        is_target = attrs.get("guid") == annotation_guid or f"GlobalId-{annotation_guid}" in classes
        colour = f"stroke:{BLACK};stroke-width:0.35;fill:none" if is_target else f"stroke:{GREY};stroke-width:0.22;fill:none;stroke-opacity:0.72"
        existing = element.attrib.get("style", "").rstrip(";")
        element.attrib["style"] = f"{existing};{colour}" if existing else colour
        black += int(is_target)
        grey += int(not is_target)
    root.attrib.update({
        "data-create-drawing-result": "FINISHED",
        "data-trap01-current-installed-colour": BLACK,
        "data-context-colour": GREY,
        "data-blue-dashed-reference-displayed": "false",
        "data-adjustable-component-mode": "fixed_body+horizontal_adjustable+vertical_adjustable",
    })
    tree.write(svg_path, encoding="utf-8", xml_declaration=True)
    content = svg_path.read_text(encoding="utf-8")
    if black == 0 or grey == 0 or BLUE.lower() in content.lower():
        raise RuntimeError("TRAP01 scene SVG colour/reference gate failed")
    return {
        "bonsai_generated_sha256_before_review_style": raw_sha,
        "current_black_geometry_count": black,
        "grey_context_geometry_count": grey,
        "blue_dashed_reference_displayed": False,
        "post_style_only": True,
    }


def inspect_svg(svg_path: Path, target_guid: str, annotation_guid: str, component_guids):
    root = ET.parse(svg_path).getroot()
    def has_ifc_identity(element, guid):
        attrs = {local_name(key): value for key, value in element.attrib.items()}
        classes = element.attrib.get("class", "").split()
        return attrs.get("guid") == guid or f"GlobalId-{guid}" in classes

    geometry_count = sum(local_name(element.tag) in GEOMETRY_TAGS for element in root.iter())
    projection_groups = sum("projection" in element.attrib.get("class", "").split() for element in root.iter())
    # Bonsai propagates the annotation's GlobalId class to every descendant
    # SVG path. Count distinct IFC identities, not descendant path occurrences.
    target_groups = int(any(has_ifc_identity(element, target_guid) for element in root.iter()))
    annotation_groups = int(any(has_ifc_identity(element, annotation_guid) for element in root.iter()))
    component_groups = sum(any(has_ifc_identity(element, item) for element in root.iter()) for item in component_guids)
    return {
        "root_tag": local_name(root.tag),
        "data_scale": root.attrib.get("data-scale"),
        "width": root.attrib.get("width"),
        "height": root.attrib.get("height"),
        "view_box": root.attrib.get("viewBox"),
        "geometry_element_count": geometry_count,
        "projection_group_count": projection_groups,
        "target_ifc_body_group_count": target_groups,
        "current_annotation_group_count": annotation_groups,
        "adjustable_component_reference_group_count": component_groups,
        "no_target_or_annotation_duplicate": target_groups == 0 and annotation_groups == 1 and component_groups == 0,
        "external_reference_count": sum("href" in local_name(key) and not value.startswith("#") for element in root.iter() for key, value in element.attrib.items()),
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
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    if sha256(FORMAL_IFC) != FORMAL_SHA256:
        raise RuntimeError("formal IFC hash mismatch before TRAP01 Drawing workflow")
    approval = json.loads(APPROVAL.read_text(encoding="utf-8"))
    candidate = json.loads(CANDIDATE.read_text(encoding="utf-8"))
    report = json.loads(WRITE_REPORT.read_text(encoding="utf-8"))
    if (
        approval.get("status") != "approved"
        or approval.get("approved_views") != ["plan", "front", "side"]
        or approval.get("derived_ifc_write_allowed") is not True
        or approval.get("formal_authoritative_ifc_write_allowed") is not False
        or set(approval.get("approved_component_scope", [])) != {"fixed_body", "horizontal_adjustable", "vertical_adjustable"}
        or report.get("component_model_pass") is not True
        or candidate.get("source_kind") != SOURCE_KIND
        or candidate.get("representative_global_id") != TARGET_GLOBAL_ID
    ):
        raise RuntimeError("TRAP01 product-level approval/component gate failed")
    for view in VIEW_DEFINITIONS:
        if len(candidate["views"][view]["proxy_paths_mm"]) != EXPECTED_PATH_COUNTS[view]:
            raise RuntimeError(f"TRAP01 {view} candidate path count drifted")
    derived_before_sha = sha256(DERIVED_IFC)
    with contextlib.suppress(Exception):
        addon_utils.disable("bl_ext.user_default.project_control", default_set=False, handle_error=None)
    load_result = bpy.ops.bim.load_project(filepath=str(DERIVED_IFC), should_start_fresh_session=True, use_detailed_tooltip=True)
    if load_result != {"FINISHED"} or not tool.Ifc.get():
        raise RuntimeError(f"Bonsai failed to load TRAP01 derived IFC: {load_result}")
    model = tool.Ifc.get()
    target = model.by_guid(TARGET_GLOBAL_ID)
    target_obj = tool.Ifc.get_object(target) if target else None
    if target is None or target_obj is None:
        raise RuntimeError("TRAP01 target is missing")
    target_bbox = shared.world_bbox(target_obj)
    centre = [(target_bbox[0][axis] + target_bbox[1][axis]) / 2 for axis in range(3)]
    scene_bbox = ((centre[0] - 1.5, centre[1] - 1.5, -0.2), (centre[0] + 1.5, centre[1] + 1.5, 3.0))
    records = shared.room_elements(scene_bbox)
    original = [record[0] for record in records]
    if original.count(target) != 1:
        raise RuntimeError("expected exactly one TRAP01 target in project crop")
    context_elements = [element for element in original if element != target]
    if len(context_elements) < 5:
        raise RuntimeError("TRAP01 actual project context crop is unexpectedly empty")
    context_counts = Counter(element.is_a() for element in context_elements)
    components = [item for item in model.by_type("IfcAnnotation") if item.ObjectType == "DRAWING_COMPONENT"]
    roles = {ifcopenshell.util.element.get_pset(item, COMPONENT_PSET).get("ComponentRole"): item for item in components}
    if set(roles) != {"fixed_body", "horizontal_adjustable", "vertical_adjustable"}:
        raise RuntimeError("TRAP01 persisted adjustable components missing before Drawing generation")
    component_guids = [item.GlobalId for item in components]

    view_records = []
    for view, definition in VIEW_DEFINITIONS.items():
        output_svg = OUTPUT_DIR / f"{definition['drawing_name']}.svg"
        drawing, camera, width, height, clip_end = shared.add_drawing(model, definition, scene_bbox, context_elements, output_svg)
        override = view3d_override()
        with bpy.context.temp_override(**override):
            activate = bpy.ops.bim.activate_drawing(drawing=drawing.id(), should_view_from_camera=False)
        if activate != {"FINISHED"}:
            raise RuntimeError(f"TRAP01 {view} Drawing activation failed")
        annotation, path_count, edge_count = add_current_annotation(model, drawing, target, target_obj, view, candidate["views"][view]["proxy_paths_mm"])
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
            raise RuntimeError(f"TRAP01 Bonsai Create Drawing failed for {view}: {create_result}")
        style = style_svg(output_svg, annotation.GlobalId)
        svg = inspect_svg(output_svg, TARGET_GLOBAL_ID, annotation.GlobalId, component_guids)
        if (
            svg["root_tag"] != "svg" or svg["geometry_element_count"] == 0
            or svg["projection_group_count"] == 0 or not svg["no_target_or_annotation_duplicate"]
        ):
            raise RuntimeError(f"TRAP01 {view} SVG structure/duplicate gate failed: {svg}")
        cache = output_svg.parent / "cache" / f"{output_svg.stem}-linework.svg"
        if not cache.is_file() or cache.stat().st_size == 0:
            raise RuntimeError(f"TRAP01 {view} Bonsai linework cache missing")
        view_records.append({
            "view": view,
            "drawing_global_id": drawing.GlobalId,
            "drawing_name": drawing.Name,
            "annotation_global_id": annotation.GlobalId,
            "path_count": path_count,
            "edge_count": edge_count,
            "output_svg": output_svg,
            "cache": cache,
            "style": style,
            "svg": svg,
            "camera": {
                "type": camera.data.type, "width_m": width, "height_m": height,
                "clip_start_m": camera.data.clip_start, "clip_end_m": clip_end,
                "matrix_world": [list(row) for row in camera.matrix_world],
                "resolution": [bpy.context.scene.render.resolution_x, bpy.context.scene.render.resolution_y],
            },
            "create_result": sorted(create_result),
        })

    temporary = DERIVED_IFC.with_suffix(".ifc.drawings-next")
    model.write(temporary)
    os.replace(temporary, DERIVED_IFC)
    derived_after_sha = sha256(DERIVED_IFC)
    reload_result = bpy.ops.bim.load_project(filepath=str(DERIVED_IFC), should_start_fresh_session=True, use_detailed_tooltip=True)
    if reload_result != {"FINISHED"} or not tool.Ifc.get():
        raise RuntimeError("persisted TRAP01 IFC reload failed")
    reopened = tool.Ifc.get()
    for record in view_records:
        drawing = reopened.by_guid(record["drawing_global_id"])
        annotation = reopened.by_guid(record["annotation_global_id"])
        if drawing is None or annotation is None or persisted_path_count(annotation) != record["path_count"]:
            raise RuntimeError(f"persisted TRAP01 {record['view']} Drawing/LINEWORK missing")
        record["persisted_path_count"] = persisted_path_count(annotation)
    reopened_components = [item for item in reopened.by_type("IfcAnnotation") if item.ObjectType == "DRAWING_COMPONENT"]
    if len(reopened_components) != 3:
        raise RuntimeError("TRAP01 component count drifted after Drawing reload")
    bpy.ops.wm.save_as_mainfile(filepath=str(SESSION_BLEND), check_existing=False)
    if sha256(FORMAL_IFC) != FORMAL_SHA256:
        raise RuntimeError("formal IFC changed during TRAP01 Drawing workflow")

    output_views = []
    for record in view_records:
        output_views.append({
            "view": record["view"],
            "drawing": {"global_id": record["drawing_global_id"], "name": record["drawing_name"]},
            "current_annotation_global_id": record["annotation_global_id"],
            "path_count": record["path_count"],
            "persisted_path_count": record["persisted_path_count"],
            "edge_count": record["edge_count"],
            "camera": record["camera"],
            "create_drawing": {
                "operator": "bpy.ops.bim.create_drawing",
                "arguments": {"print_all": False, "open_viewer": False, "sync": False},
                "result": record["create_result"], "linework_mode": "OPENCASCADE",
                "target_view": VIEW_DEFINITIONS[record["view"]]["target_view"],
            },
            "svg": {"path": str(record["output_svg"]), "bytes": record["output_svg"].stat().st_size, "sha256": sha256(record["output_svg"]), **record["style"], **record["svg"]},
            "linework_cache": {"path": str(record["cache"]), "bytes": record["cache"].stat().st_size, "sha256": sha256(record["cache"])},
        })
    evidence = {
        "schema_version": 1,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "task": "approved TRAP01 componentised derived IFC plus Bonsai Create Drawing Plan/Front/Side in actual project context",
        "courseEvidence": {
            "source_path": "/Users/jiaxinchen/.codex/plugins/cache/personal/bonsai-course-operator/0.1.1/evals/course-index.json",
            "lesson": "085000 Introduction to Drawings",
            "timestamps": ["01:59 Create Drawing", "02:06 Drawing created", "02:13 SVG inspection"],
            "course_fact": "Create Drawing generates or refreshes SVG after Drawing camera, scale, depth and filters are configured.",
            "current_inference": "Bonsai 0.8.4 semantic Drawing state and bpy.ops.bim.create_drawing were used in Blender 4.5.3 LTS.",
        },
        "plan": ["inspect approved derived IFC", "create three scoped Drawings", "add current-installed black LINEWORK", "run native Create Drawing", "persist", "reload", "inspect SVG and duplicate state"],
        "preState": {"derived_ifc": str(DERIVED_IFC), "derived_ifc_sha256": derived_before_sha, "formal_ifc_sha256": sha256(FORMAL_IFC), "component_count": len(components), "drawing_count": 0},
        "execution": {"provider": "public bonsai-mcp b0e67b1", "bridge": "127.0.0.1:9878", "capability": "execute_blender_code", "generator": "bpy.ops.bim.create_drawing", "linework_mode": "OPENCASCADE"},
        "persistence": {"method": "model.write to temporary then atomic replace; bpy.ops.bim.load_project reload; bpy.ops.wm.save_as_mainfile", "reload_result": sorted(reload_result)},
        "postState": {"derived_ifc": str(DERIVED_IFC), "derived_ifc_sha256": derived_after_sha, "formal_ifc_sha256": sha256(FORMAL_IFC), "component_count": len(reopened_components), "drawing_count": len(view_records), "current_linework_annotation_count": len(view_records)},
        "outputs": {"views": output_views, "session_blend": {"path": str(SESSION_BLEND), "bytes": SESSION_BLEND.stat().st_size, "sha256": sha256(SESSION_BLEND)}},
        "visual": {"preview_manifest": str(PRODUCT_DIR / "TRAP01-project-drawing-manifest.json"), "blue_dashed_reference_displayed": False},
        "target": {"global_id": TARGET_GLOBAL_ID, "world_bbox_m": [list(target_bbox[0]), list(target_bbox[1])], "world_centre_m": centre},
        "context": {"project_context_retained": True, "include_count": len(context_elements), "include_ifc_class_counts": dict(sorted(context_counts.items())), "actual_target_body_suppressed_from_drawing_include": True, "current_installed_linework_inserted_once": True},
        "adjustable_components": {role: {"global_id": item.GlobalId, "pset": ifcopenshell.util.element.get_pset(item, COMPONENT_PSET)} for role, item in roles.items()},
        "versions": {"blender": bpy.app.version_string, "bonsai": "0.8.4", "ifcopenshell": ifcopenshell.version},
        "tests": {"all_create_drawing_finished": all(item["create_drawing"]["result"] == ["FINISHED"] for item in output_views), "all_target_body_projection_counts_zero": all(item["svg"]["target_ifc_body_group_count"] == 0 for item in output_views), "all_current_annotation_counts_one": all(item["svg"]["current_annotation_group_count"] == 1 for item in output_views), "all_component_reference_counts_zero": all(item["svg"]["adjustable_component_reference_group_count"] == 0 for item in output_views), "blue_dashed_reference_displayed": False},
        "formal_ifc_bytes_unchanged": sha256(FORMAL_IFC) == FORMAL_SHA256,
        "verdict": "pass",
        "pass": True,
    }
    EVIDENCE.write_text(json.dumps(evidence, indent=2, ensure_ascii=False, default=str) + "\n", encoding="utf-8")
    print(json.dumps({"evidence": str(EVIDENCE), "derived_ifc_sha256": derived_after_sha, "views": [item["svg"] for item in output_views], "pass": True}, indent=2))


try:
    main()
except Exception:
    (PRODUCT_DIR / "trap01-create-drawing-error.log").write_text(traceback.format_exc(), encoding="utf-8")
    raise
