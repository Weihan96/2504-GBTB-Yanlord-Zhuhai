#!/usr/bin/env python3
"""Create persisted BED01 Bonsai Drawings in the actual master bedroom.

The input is a byte-identical copy of the formal IFC. The BED01 Body is
excluded from Drawing Include and replaced by the approved Baxter native 2D
DWG linework as one persistent IFC LINEWORK Annotation per view.
"""

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
from mathutils import Vector


ROOT = Path(__file__).resolve().parents[2]
SCRIPTS_DIR = ROOT / "pipeline/scripts"
sys.path.insert(0, str(SCRIPTS_DIR))
import create_gessi316_54294_main_bathroom_drawing as shared  # noqa: E402


PRODUCT_DIR = ROOT / "output/review/highpoly-types/bed01"
FORMAL_IFC = ROOT / "2504 GBTB Yanlord Zhuhai.ifc"
FORMAL_SHA256 = "7a50b87e8f48a7c2bfbaa9f1dabdfca7155fa684b2325c66aa1ae15c8ab0a25c"
APPROVAL = ROOT / "pipeline/decisions/bed01-drawing-approval.json"
REFERENCE = PRODUCT_DIR / "official-dwg-review-reference.json"
REVIEW_MANIFEST = PRODUCT_DIR / "official-dwg-review-manifest.json"
REVIEW_MANIFEST_SHA256 = "1063e4d07ddf902f8c0ba66a5ceeae7a3ca53cf18736c522b7600ce8e8b2138a"
TARGET_GLOBAL_ID = "3IQBEqO5vDI8Z9k1Ltge_N"
ROOM_GLOBAL_ID = "3gHz6U6BfFXgV6PnRzfOf$"
ROOM_NAME = "主卧"
SOURCE_DWG_SHA256 = "ab697db9c27448c62a9a77537d7cc4b6286577c335328113d703184be28d9c4a"
SOURCE_KIND = "native_dwg_review_reference"
SOURCE_LABEL_ZH = "Baxter Casablanca 官方独立 2D DWG 原生蓝线"
BLUE = "#1677c8"
GREY = "#a3abb3"
EXPECTED_PATH_COUNTS = {"plan": 25, "front": 22, "side": 36}
VIEW_DEFINITIONS = {
    "plan": {
        "drawing_name": "BED01-MASTER-BEDROOM-PLAN",
        "target_view": "PLAN_VIEW",
        "location_hint": "PLAN",
    },
    "front": {
        "drawing_name": "BED01-MASTER-BEDROOM-FRONT",
        "target_view": "ELEVATION_VIEW",
        "location_hint": "EAST",
    },
    "side": {
        "drawing_name": "BED01-MASTER-BEDROOM-SIDE",
        "target_view": "ELEVATION_VIEW",
        "location_hint": "SOUTH",
    },
}
GEOMETRY_TAGS = {"path", "polyline", "polygon", "line", "circle", "ellipse", "rect"}


shared.PRODUCT_DIR = PRODUCT_DIR
shared.FORMAL_IFC = FORMAL_IFC
shared.FORMAL_SHA256 = FORMAL_SHA256
shared.TARGET_GLOBAL_ID = TARGET_GLOBAL_ID
shared.ROOM_GLOBAL_ID = ROOM_GLOBAL_ID
shared.ROOM_NAME = ROOM_NAME
shared.ARTICLE = "Casablanca 180"
shared.DRAWING_PRODUCT_CODE = "BED01"
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


def bounds_2d(paths):
    points = [point for path in paths for point in path]
    return {
        "minimum": [min(point[axis] for point in points) for axis in range(2)],
        "maximum": [max(point[axis] for point in points) for axis in range(2)],
    }


def align_paths(view: str, paths: list, alignment: dict) -> list:
    source = [[list(point) for point in path] for path in paths]
    reflect = alignment["view_direction_reflection_x"]
    source_bounds = bounds_2d(source)
    if reflect:
        centre_x = (source_bounds["minimum"][0] + source_bounds["maximum"][0]) / 2
        source = [
            [[2 * centre_x - first, second] for first, second in path]
            for path in source
        ]
    translate_x, translate_y = alignment["translation_mm"]
    aligned = [
        [[first + translate_x, second + translate_y] for first, second in path]
        for path in source
    ]
    expected_mode = (
        "rigid_reflection_and_translation_only" if view == "side" else "translation_only"
    )
    if (
        alignment["mode"] != expected_mode
        or reflect is not (view == "side")
        or alignment["uniform_scale"] != 1.0
        or alignment["anisotropic_scale_used"] is not False
        or alignment["source_geometry_deformed"] is not False
    ):
        raise RuntimeError(f"{view} rigid alignment gate failed")
    return aligned


def coordinates_mm(view: str, first: float, second: float):
    if view == "plan":
        return float(first), float(second), 0.0
    if view == "front":
        return float(first), 0.0, float(second)
    return 0.0, float(first), float(second)


def annotation_representation(model, context, view, paths):
    polylines = []
    for path in paths:
        points = [
            model.create_entity(
                "IfcCartesianPoint", Coordinates=coordinates_mm(view, first, second)
            )
            for first, second in path
        ]
        if len(points) >= 2:
            polylines.append(model.create_entity("IfcPolyline", Points=points))
    if len(polylines) != EXPECTED_PATH_COUNTS[view]:
        raise RuntimeError(f"{view} official path count drifted")
    curve_set = model.create_entity("IfcGeometricCurveSet", Elements=polylines)
    colour = model.create_entity(
        "IfcColourRgb",
        Name="Baxter official DWG blue",
        Red=22.0 / 255.0,
        Green=119.0 / 255.0,
        Blue=200.0 / 255.0,
    )
    curve_style = model.create_entity(
        "IfcCurveStyle",
        Name="Baxter Casablanca official native 2D DWG linework",
        CurveFont=None,
        CurveWidth=model.create_entity("IfcPositiveLengthMeasure", 0.35),
        CurveColour=colour,
        ModelOrDraughting=True,
    )
    model.create_entity(
        "IfcStyledItem",
        Item=curve_set,
        Styles=[curve_style],
        Name="Blue official native 2D LINEWORK",
    )
    representation = model.create_entity(
        "IfcShapeRepresentation",
        ContextOfItems=context,
        RepresentationIdentifier="Annotation",
        RepresentationType="GeometricCurveSet",
        Items=[curve_set],
    )
    return representation, len(polylines), sum(max(0, len(path) - 1) for path in paths)


def add_official_annotation(model, drawing, target, target_obj, view, paths, alignment):
    target_view = VIEW_DEFINITIONS[view]["target_view"]
    context = tool.Drawing.get_annotation_context(target_view)
    if context is None:
        context = tool.Drawing.create_annotation_context(target_view)
    obj = core_drawing.add_annotation(
        tool.Ifc,
        tool.Collector,
        tool.Drawing,
        drawing=drawing,
        object_type="LINEWORK",
        relating_type=None,
        enable_editing=False,
    )
    annotation = tool.Ifc.get_entity(obj)
    annotation.Name = f"Baxter Casablanca official native 2D DWG / {view}"
    annotation.Description = (
        f"{SOURCE_LABEL_ZH}; approved product-level Drawing reference for {TARGET_GLOBAL_ID}; "
        "1:1 rigid placement only; prior BIM ACIS envelope excluded."
    )
    annotation.ObjectPlacement = target.ObjectPlacement
    obj.matrix_world = target_obj.matrix_world

    vertices = []
    edges = []
    for path in paths:
        indices = []
        for first, second in path:
            indices.append(len(vertices))
            vertices.append(tuple(value / 1000 for value in coordinates_mm(view, first, second)))
        edges.extend(zip(indices, indices[1:]))
    obj.data.clear_geometry()
    obj.data.from_pydata(vertices, edges, [])
    obj.data.update()

    representation, path_count, edge_count = annotation_representation(
        model, context, view, paths
    )
    annotation.Representation = model.create_entity(
        "IfcProductDefinitionShape", Representations=[representation]
    )
    pset = ifcopenshell.api.pset.add_pset(model, product=annotation, name="EPset_Annotation")
    ifcopenshell.api.pset.edit_pset(
        model,
        pset=pset,
        properties={
            "Classes": "review-target-bed01 official-native-dwg",
            "TargetGlobalId": TARGET_GLOBAL_ID,
            "SourceKind": SOURCE_KIND,
            "SourceLabelZh": SOURCE_LABEL_ZH,
            "SourceDwgSha256": SOURCE_DWG_SHA256,
            "OfficialNative2dDwg": True,
            "SourceScale": 1.0,
            "AlignmentMode": alignment["mode"],
            "TranslationMillimetres": ",".join(
                f"{value:.6f}" for value in alignment["translation_mm"]
            ),
            "ViewDirectionReflectionX": alignment["view_direction_reflection_x"],
            "SupersededAcisCandidateExcluded": True,
            "IfcCurveStyleColour": BLUE,
        },
    )
    return annotation, representation, path_count, edge_count


def style_svg(svg_path: Path, annotation_global_id: str, view: str):
    raw_sha = sha256(svg_path)
    ET.register_namespace("", "http://www.w3.org/2000/svg")
    ET.register_namespace("ifc", "http://www.ifcopenshell.org/ns")
    tree = ET.parse(svg_path)
    root = tree.getroot()
    annotation_groups = []
    for element in root.iter():
        attributes = {local_name(key): value for key, value in element.attrib.items()}
        classes = element.attrib.get("class", "").split()
        if (
            attributes.get("guid") == annotation_global_id
            or f"GlobalId-{annotation_global_id}" in classes
            or "review-target-bed01" in classes
        ):
            annotation_groups.append(element)
    if not annotation_groups:
        raise RuntimeError("Bonsai SVG has no BED01 official annotation group")
    blue_ids = {
        id(element)
        for group in annotation_groups
        for element in group.iter()
        if local_name(element.tag) in GEOMETRY_TAGS
    }
    blue_geometry = 0
    grey_geometry = 0
    for element in root.iter():
        if local_name(element.tag) not in GEOMETRY_TAGS:
            continue
        existing = element.attrib.get("style", "").rstrip(";")
        if id(element) in blue_ids:
            colour = f"stroke:{BLUE};stroke-width:0.35;fill:none"
            blue_geometry += 1
        else:
            colour = f"stroke:{GREY};stroke-width:0.22;fill:none;stroke-opacity:0.72"
            grey_geometry += 1
        element.attrib["style"] = f"{existing};{colour}" if existing else colour
    for group in annotation_groups:
        classes = group.attrib.get("class", "").split()
        for class_name in ("review-target-bed01", "official-native-dwg"):
            if class_name not in classes:
                classes.append(class_name)
        group.attrib["class"] = " ".join(classes)
        group.attrib["data-source-kind"] = SOURCE_KIND
        group.attrib["data-source-dwg-sha256"] = SOURCE_DWG_SHA256
        group.attrib["data-source-scale"] = "1.0"
        group.attrib["data-reviewed-view"] = view
    root.attrib["data-create-drawing-result"] = "FINISHED"
    root.attrib["data-context-colour"] = GREY
    root.attrib["data-official-linework-colour"] = BLUE
    tree.write(svg_path, encoding="utf-8", xml_declaration=True)
    if blue_geometry == 0 or grey_geometry == 0:
        raise RuntimeError("SVG colour-layer gate failed")
    return {
        "bonsai_generated_sha256_before_review_style": raw_sha,
        "annotation_group_count": len(annotation_groups),
        "blue_geometry_element_count": blue_geometry,
        "grey_context_geometry_element_count": grey_geometry,
        "post_style_only": True,
    }


def inspect_svg(svg_path: Path, target_global_id: str, annotation_global_id: str):
    result = shared.inspect_svg(svg_path, target_global_id, annotation_global_id)
    root = ET.parse(svg_path).getroot()
    styles = [element.attrib.get("style", "") for element in root.iter()]
    result["blue_style_count"] = sum(BLUE in style for style in styles)
    result["grey_style_count"] = sum(GREY in style for style in styles)
    return result


def bounds_3d(points):
    return [
        [min(point[axis] for point in points) for axis in range(3)],
        [max(point[axis] for point in points) for axis in range(3)],
    ]


def main():
    arguments = sys.argv[sys.argv.index("--") + 1 :]
    if len(arguments) != 3:
        raise SystemExit("expected: derived.ifc formal.ifc output-directory")
    derived_ifc = Path(arguments[0]).resolve()
    formal_ifc = Path(arguments[1]).resolve()
    output_dir = Path(arguments[2]).resolve()
    if sha256(formal_ifc) != FORMAL_SHA256 or sha256(FORMAL_IFC) != FORMAL_SHA256:
        raise RuntimeError("formal IFC hash mismatch")
    if sha256(derived_ifc) != FORMAL_SHA256:
        raise RuntimeError("fresh BED01 derived IFC must start byte-identical to formal IFC")
    derived_before_sha = sha256(derived_ifc)

    approval = json.loads(APPROVAL.read_text(encoding="utf-8"))
    reference = json.loads(REFERENCE.read_text(encoding="utf-8"))
    review_manifest = json.loads(REVIEW_MANIFEST.read_text(encoding="utf-8"))
    if (
        approval.get("status") != "approved"
        or approval.get("derived_ifc_write_allowed") is not True
        or approval.get("formal_authoritative_ifc_write_allowed") is not False
        or approval.get("approved_views") != ["plan", "front", "side"]
        or approval.get("candidate_manifest_sha256") != REVIEW_MANIFEST_SHA256
        or sha256(REVIEW_MANIFEST) != REVIEW_MANIFEST_SHA256
        or reference.get("source_kind") != SOURCE_KIND
        or reference.get("source_dwg_sha256") != SOURCE_DWG_SHA256
        or reference.get("source_geometry_scaled") is not False
        or reference.get("source_geometry_anisotropically_fitted") is not False
        or not reference.get("superseded_error_candidates")
    ):
        raise RuntimeError("BED01 approval/native-DWG source gate failed")

    alignments = {item["view"]: item["alignment"] for item in review_manifest["views"]}
    aligned_views = {}
    for view in VIEW_DEFINITIONS:
        raw_paths = reference["views"][view]["paths_mm"]
        if len(raw_paths) != EXPECTED_PATH_COUNTS[view]:
            raise RuntimeError(f"{view} source path count drifted")
        aligned_views[view] = align_paths(view, raw_paths, alignments[view])

    os.chdir(ROOT)
    with contextlib.suppress(Exception):
        addon_utils.disable(
            "bl_ext.user_default.project_control", default_set=False, handle_error=None
        )
    load_result = bpy.ops.bim.load_project(
        filepath=str(derived_ifc),
        should_start_fresh_session=True,
        use_detailed_tooltip=True,
    )
    if load_result != {"FINISHED"} or not tool.Ifc.get():
        raise RuntimeError(f"Bonsai failed to load BED01 derived IFC: {load_result}")
    model = tool.Ifc.get()
    target = model.by_guid(TARGET_GLOBAL_ID)
    room = model.by_guid(ROOM_GLOBAL_ID)
    target_obj = tool.Ifc.get_object(target)
    room_obj = tool.Ifc.get_object(room)
    if target is None or room is None or target_obj is None or room_obj is None:
        raise RuntimeError("BED01 target or master-bedroom space missing")
    if room.LongName != ROOM_NAME:
        raise RuntimeError("master-bedroom identity drifted")

    room_bbox = shared.world_bbox(room_obj)
    target_bbox = shared.world_bbox(target_obj)
    records = shared.room_elements(room_bbox)
    original_context = [record[0] for record in records]
    if original_context.count(target) != 1:
        raise RuntimeError("expected BED01 exactly once in master-bedroom context")
    context_elements = [element for element in original_context if element != target]
    context_counts = Counter(element.is_a() for element in context_elements)
    output_dir.mkdir(parents=True, exist_ok=True)

    view_records = []
    for view, definition in VIEW_DEFINITIONS.items():
        output_svg = output_dir / f"{definition['drawing_name']}.svg"
        drawing, camera, width, height, clip_end = shared.add_drawing(
            model, definition, room_bbox, context_elements, output_svg
        )
        activate_result = bpy.ops.bim.activate_drawing(
            drawing=drawing.id(), should_view_from_camera=False
        )
        if activate_result != {"FINISHED"}:
            raise RuntimeError(f"failed to activate {view} Drawing: {activate_result}")
        annotation, _, path_count, edge_count = add_official_annotation(
            model,
            drawing,
            target,
            target_obj,
            view,
            aligned_views[view],
            alignments[view],
        )
        cprops = tool.Drawing.get_camera_props(camera)
        cprops.has_annotation = True
        drawing_pset_data = ifcopenshell.util.element.get_pset(drawing, "EPset_Drawing")
        ifcopenshell.api.pset.edit_pset(
            model,
            pset=model.by_id(drawing_pset_data["id"]),
            properties={"HasAnnotation": True},
        )
        dprops = tool.Drawing.get_document_props()
        dprops.should_use_underlay_cache = False
        dprops.should_use_linework_cache = False
        dprops.should_use_annotation_cache = False
        create_result = bpy.ops.bim.create_drawing(
            print_all=False, open_viewer=False, sync=False
        )
        if (
            create_result != {"FINISHED"}
            or not output_svg.is_file()
            or output_svg.stat().st_size == 0
        ):
            raise RuntimeError(f"Bonsai Create Drawing failed for {view}: {create_result}")
        style = style_svg(output_svg, annotation.GlobalId, view)
        svg = inspect_svg(output_svg, TARGET_GLOBAL_ID, annotation.GlobalId)
        if (
            svg["root_tag"] != "svg"
            or svg["geometry_element_count"] == 0
            or svg["projection_group_count"] == 0
            or svg["target_ifc_projection_group_count"] != 0
            or svg["review_annotation_geometry_count"] == 0
            or svg["blue_style_count"] == 0
            or svg["grey_style_count"] == 0
        ):
            raise RuntimeError(f"BED01 {view} SVG structure/layer gate failed: {svg}")
        cache_path = output_svg.parent / "cache" / f"{output_svg.stem}-linework.svg"
        if not cache_path.is_file() or cache_path.stat().st_size == 0:
            raise RuntimeError(f"BED01 {view} linework cache missing")
        view_records.append(
            {
                "view": view,
                "drawing_global_id": drawing.GlobalId,
                "drawing_name": drawing.Name,
                "annotation_global_id": annotation.GlobalId,
                "path_count": path_count,
                "edge_count": edge_count,
                "expected_points": [
                    coordinates_mm(view, first, second)
                    for path in aligned_views[view]
                    for first, second in path
                ],
                "output_svg": output_svg,
                "cache_path": cache_path,
                "create_result": sorted(create_result),
                "style": style,
                "svg": svg,
                "camera": {
                    "type": camera.data.type,
                    "width_m": width,
                    "height_m": height,
                    "clip_start_m": camera.data.clip_start,
                    "clip_end_m": clip_end,
                    "resolution": [
                        bpy.context.scene.render.resolution_x,
                        bpy.context.scene.render.resolution_y,
                    ],
                },
            }
        )

    temporary = derived_ifc.with_suffix(".ifc.next")
    model.write(str(temporary))
    os.replace(temporary, derived_ifc)
    derived_after_sha = sha256(derived_ifc)
    reload_result = bpy.ops.bim.load_project(
        filepath=str(derived_ifc),
        should_start_fresh_session=True,
        use_detailed_tooltip=True,
    )
    if reload_result != {"FINISHED"} or not tool.Ifc.get():
        raise RuntimeError(f"Bonsai reload failed: {reload_result}")
    reopened = tool.Ifc.get()

    output_views = []
    for record in view_records:
        view = record["view"]
        drawing = reopened.by_guid(record["drawing_global_id"])
        annotation = reopened.by_guid(record["annotation_global_id"])
        if drawing is None or annotation is None:
            raise RuntimeError(f"reloaded IFC lost {view} Drawing or Annotation")
        representation = next(
            (
                item
                for item in annotation.Representation.Representations
                if item.RepresentationIdentifier == "Annotation"
                and item.RepresentationType == "GeometricCurveSet"
                and item.ContextOfItems.TargetView == VIEW_DEFINITIONS[view]["target_view"]
            ),
            None,
        )
        if representation is None:
            raise RuntimeError(f"reloaded IFC lost {view} Annotation representation")
        paths = [path for item in representation.Items for path in item.Elements]
        points = [
            tuple(float(value) for value in point.Coordinates)
            for path in paths
            for point in path.Points
        ]
        expected_bounds = bounds_3d(record["expected_points"])
        persisted_bounds = bounds_3d(points)
        residual = max(
            abs(persisted_bounds[bound][axis] - expected_bounds[bound][axis])
            for bound in range(2)
            for axis in range(3)
        )
        annotation_pset = ifcopenshell.util.element.get_pset(annotation, "EPset_Annotation")
        drawing_pset = ifcopenshell.util.element.get_pset(drawing, "EPset_Drawing")
        include = drawing_pset.get("Include", "").split(",")
        drawing_assignment = next(
            (
                relation
                for relation in reopened.by_type("IfcRelAssignsToGroup")
                if drawing in relation.RelatedObjects and annotation in relation.RelatedObjects
            ),
            None,
        )
        styled_items = [
            styled for item in representation.Items for styled in item.StyledByItem
        ]
        if (
            len(paths) != EXPECTED_PATH_COUNTS[view]
            or residual > 0.000001
            or annotation_pset.get("OfficialNative2dDwg") is not True
            or annotation_pset.get("SourceScale") != 1.0
            or not styled_items
            or drawing_assignment is None
            or TARGET_GLOBAL_ID in include
            or len(include) != len(context_elements)
        ):
            raise RuntimeError(f"persisted BED01 {view} gate failed")
        output_views.append(
            {
                "view": view,
                "drawing": {
                    "global_id": drawing.GlobalId,
                    "name": drawing.Name,
                    "object_type": drawing.ObjectType,
                    "epset_drawing": drawing_pset,
                },
                "linework_annotation": {
                    "global_id": annotation.GlobalId,
                    "name": annotation.Name,
                    "epset_annotation": annotation_pset,
                    "ifc_curve_style_persisted": True,
                    "ifc_curve_style_colour": BLUE,
                    "drawing_group_association_persisted": True,
                    "drawing_group_assignment_global_id": drawing_assignment.GlobalId,
                },
                "review_path_count": EXPECTED_PATH_COUNTS[view],
                "persisted_review_path_count": len(paths),
                "review_edge_count": record["edge_count"],
                "persisted_coordinate_bounds_mm": persisted_bounds,
                "expected_coordinate_bounds_mm": expected_bounds,
                "persisted_coordinate_residual_mm": residual,
                "alignment": alignments[view],
                "camera": record["camera"],
                "create_drawing": {
                    "operator": "bpy.ops.bim.create_drawing",
                    "arguments": {
                        "print_all": False,
                        "open_viewer": False,
                        "sync": False,
                    },
                    "result": record["create_result"],
                    "linework_mode": "OPENCASCADE",
                    "target_view": VIEW_DEFINITIONS[view]["target_view"],
                },
                "svg": {
                    "path": str(record["output_svg"]),
                    "bytes": record["output_svg"].stat().st_size,
                    "sha256": sha256(record["output_svg"]),
                    **record["style"],
                    **record["svg"],
                },
                "linework_cache": {
                    "path": str(record["cache_path"]),
                    "bytes": record["cache_path"].stat().st_size,
                    "sha256": sha256(record["cache_path"]),
                },
            }
        )

    if sha256(formal_ifc) != FORMAL_SHA256 or sha256(FORMAL_IFC) != FORMAL_SHA256:
        raise RuntimeError("formal IFC changed during BED01 Drawing creation")
    evidence = {
        "schema_version": 1,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "task": "Bonsai Create Drawing for approved Baxter Casablanca native 2D DWG linework in the actual master-bedroom context",
        "workflow": ["inspect", "plan", "execute", "persist", "reload", "verify"],
        "provider": {
            "preferred": "bonsai-mcp 1.1.0",
            "status": "unsupported",
            "reason": "Blender bridge 127.0.0.1:9878 connection refused",
            "fallback": "local Blender 4.5.3 LTS with installed Bonsai/IfcOpenShell 0.8.4",
        },
        "course_evidence": {
            "query": "085000 Create Drawing SVG persistent IFC annotation linework",
            "top_lesson": "085000 Introduction to Drawings",
            "timestamps": ["01:59 Create Drawing", "02:06 Drawing created", "02:13 SVG in browser"],
            "course_fact": "Create Drawing generates or refreshes SVG after Drawing camera, scale, depth and filters are configured.",
            "source_files_packaged": False,
            "source_files_note": "The installed skill contains indexed lesson records but not the referenced private source Markdown or screenshots.",
            "current_version_inference": "The recorded workflow was adapted to Blender 4.5.3 LTS and Bonsai/IfcOpenShell 0.8.4 using semantic Drawing state and bpy.ops.bim.create_drawing.",
        },
        "formal_ifc": str(formal_ifc),
        "formal_ifc_sha256": sha256(formal_ifc),
        "formal_ifc_bytes_unchanged": sha256(formal_ifc) == FORMAL_SHA256,
        "derived_ifc": str(derived_ifc),
        "derived_ifc_sha256_before": derived_before_sha,
        "derived_ifc_sha256_after": derived_after_sha,
        "save_boundary": "one BED01 product-level full-project IFC copy under the BED01 review package",
        "approval_record": str(APPROVAL),
        "approval_record_sha256": sha256(APPROVAL),
        "source_kind": SOURCE_KIND,
        "source_label_zh": SOURCE_LABEL_ZH,
        "source_dwg_sha256": SOURCE_DWG_SHA256,
        "source_scale": 1.0,
        "superseded_acis_candidate_excluded": True,
        "target": {
            "global_id": TARGET_GLOBAL_ID,
            "ifc_class": target.is_a(),
            "type_name": next(relation.RelatingType.Name for relation in target.IsTypedBy),
            "world_bbox_m": [list(target_bbox[0]), list(target_bbox[1])],
            "target_include_count_before_suppression": original_context.count(target),
            "target_include_count_after_suppression": 0,
        },
        "room": {
            "global_id": ROOM_GLOBAL_ID,
            "name": ROOM_NAME,
            "bbox_m": [list(room_bbox[0]), list(room_bbox[1])],
        },
        "context": {
            "project_context_retained": True,
            "include_count": len(context_elements),
            "include_global_ids": [element.GlobalId for element in context_elements],
            "include_ifc_class_counts": dict(sorted(context_counts.items())),
            "context_colour": GREY,
        },
        "views": output_views,
        "persistence": {
            "method": "model.write to temporary then atomic replace; bpy.ops.bim.load_project reload",
            "reload_result": sorted(reload_result),
            "post_reload_drawing_count": len(output_views),
            "post_reload_annotation_count": len(output_views),
        },
        "versions": {
            "blender": bpy.app.version_string,
            "ifcopenshell": ifcopenshell.version,
            "ifc_schema": reopened.schema,
            "bonsai_generator": "Bonsai 0.8.4 bim.create_drawing",
        },
        "pass": True,
    }
    evidence_path = output_dir / "BED01-MASTER-BEDROOM-create-drawing-evidence.json"
    evidence_path.write_text(
        json.dumps(evidence, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print(
        json.dumps(
            {
                "derived_ifc": str(derived_ifc),
                "evidence": str(evidence_path),
                "create_drawing_results": {
                    item["view"]: item["create_drawing"]["result"] for item in output_views
                },
                "pass": True,
            },
            ensure_ascii=False,
        )
    )


if __name__ == "__main__":
    main()
