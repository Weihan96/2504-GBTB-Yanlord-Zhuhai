#!/usr/bin/env python3
"""Create one native Bonsai Drawing for Gessi 54294 in the actual main bathroom.

Run inside Blender against a fresh per-view copy of the formal IFC.  The formal
IFC is never opened for write.  The target Body is suppressed from Drawing
Include and replaced by the approved de-textured official-outline candidate as
an IFC LINEWORK Annotation persisted in the drawing session.
"""

from __future__ import annotations

import contextlib
import hashlib
import json
import os
import re
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
import ifcopenshell.util.placement
from bonsai import tool
from bonsai.core import drawing as core_drawing
from mathutils import Matrix, Vector


ROOT = Path(__file__).resolve().parents[2]
PRODUCT_DIR = ROOT / "output/review/highpoly-types/gessi316-54294"
FORMAL_IFC = ROOT / "2504 GBTB Yanlord Zhuhai.ifc"
FORMAL_SHA256 = "7a50b87e8f48a7c2bfbaa9f1dabdfca7155fa684b2325c66aa1ae15c8ab0a25c"
CANDIDATE = PRODUCT_DIR / "candidate-representations.json"
TARGET_GLOBAL_ID = "2iKOL78$H0N9Yd9$ky3pW4"
ROOM_GLOBAL_ID = "3a4COIs5X7lgDirMBDT4Vs"
ROOM_NAME = "主卫湿区"
FINISHED_WALL_COVERING_GLOBAL_ID = "3WekjaeUn1qfL6oR1KYD_2"
SIDE_OCCLUDING_CONTEXT_GLOBAL_IDS = {
    "0nlHIdaVrV7xIKsYXraYCa",  # south boundary wall
    "04rs0EDjn2EvxytEQSxWRB",  # overlapping south boundary wall
    "26PwmC1AX6ZROC65DohUg5",  # south wall finished covering
}
ARTICLE = "45089_54294"
DRAWING_PRODUCT_CODE = "54294"
SOURCE_DWG_SHA256 = "dfe8bcf7e100c14187c387b75a35e3fe7f718c61595fadf66adfe92ca7648ff4"
SOURCE_KIND = "native_dwg_review_simplification"
SOURCE_LABEL_ZH = "基于官方 Gessi 54294 原生 DWG 轮廓的去纹审核简化表达"
BLUE = "#1677c8"
EXPECTED_PATH_COUNTS = {"plan": 69, "front": 1099, "side": 35}
VIEW_DEFINITIONS = {
    "plan": {
        "drawing_name": "GESSI316-54294-MAIN-BATH-PLAN",
        "target_view": "PLAN_VIEW",
        "location_hint": "PLAN",
    },
    "front": {
        "drawing_name": "GESSI316-54294-MAIN-BATH-FRONT",
        "target_view": "ELEVATION_VIEW",
        "location_hint": "EAST",
    },
    "side": {
        "drawing_name": "GESSI316-54294-MAIN-BATH-SIDE",
        "target_view": "ELEVATION_VIEW",
        "location_hint": "SOUTH",
    },
}
BOUNDARY_CLASSES = {
    "IfcWall", "IfcWallStandardCase", "IfcSlab", "IfcBeam", "IfcCovering",
    "IfcDoor", "IfcWindow",
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


def world_bbox(obj):
    corners = [obj.matrix_world @ Vector(corner) for corner in obj.bound_box]
    return (
        tuple(min(point[axis] for point in corners) for axis in range(3)),
        tuple(max(point[axis] for point in corners) for axis in range(3)),
    )


def room_elements(room_bbox):
    minimum, maximum = room_bbox
    records = []
    for obj in bpy.data.objects:
        element = tool.Ifc.get_entity(obj)
        if (
            element is None
            or not element.is_a("IfcElement")
            or element.is_a("IfcOpeningElement")
            or element.is_a("IfcVirtualElement")
            or obj.type != "MESH"
            or element.Representation is None
        ):
            continue
        obj_minimum, obj_maximum = world_bbox(obj)
        if not (
            obj_maximum[0] >= minimum[0]
            and obj_minimum[0] <= maximum[0]
            and obj_maximum[1] >= minimum[1]
            and obj_minimum[1] <= maximum[1]
            and obj_maximum[2] >= -0.20
            and obj_minimum[2] <= 3.00
        ):
            continue
        centre = tuple((obj_minimum[axis] + obj_maximum[axis]) / 2 for axis in range(3))
        centre_in_room = (
            minimum[0] - 0.05 <= centre[0] <= maximum[0] + 0.05
            and minimum[1] - 0.05 <= centre[1] <= maximum[1] + 0.05
        )
        if centre_in_room or element.is_a() in BOUNDARY_CLASSES:
            records.append((element, obj, obj_minimum, obj_maximum))
    records.sort(key=lambda item: item[0].GlobalId)
    return records


def candidate_coordinates_mm(view, first, second):
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
                "IfcCartesianPoint",
                Coordinates=candidate_coordinates_mm(view, first, second),
            )
            for first, second in path
        ]
        if len(points) >= 2:
            polylines.append(model.create_entity("IfcPolyline", Points=points))
    if len(polylines) != EXPECTED_PATH_COUNTS[view]:
        raise RuntimeError(f"{view} review path count drifted")
    curve_set = model.create_entity("IfcGeometricCurveSet", Elements=polylines)
    colour = model.create_entity(
        "IfcColourRgb",
        Name="Gessi review blue",
        Red=22.0 / 255.0,
        Green=119.0 / 255.0,
        Blue=200.0 / 255.0,
    )
    curve_style = model.create_entity(
        "IfcCurveStyle",
        Name="Gessi 54294 de-textured review linework",
        CurveFont=None,
        CurveWidth=model.create_entity("IfcPositiveLengthMeasure", 0.35),
        CurveColour=colour,
        ModelOrDraughting=True,
    )
    model.create_entity(
        "IfcStyledItem",
        Item=curve_set,
        Styles=[curve_style],
        Name="Blue review LINEWORK",
    )
    representation = model.create_entity(
        "IfcShapeRepresentation",
        ContextOfItems=context,
        RepresentationIdentifier="Annotation",
        RepresentationType="GeometricCurveSet",
        Items=[curve_set],
    )
    return representation, len(polylines), sum(max(0, len(path) - 1) for path in paths)


def add_review_annotation(model, drawing, target, target_obj, view, paths):
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
    annotation.Name = f"Gessi {DRAWING_PRODUCT_CODE} de-textured official-outline review / {view}"
    annotation.Description = (
        f"{SOURCE_LABEL_ZH}; target {TARGET_GLOBAL_ID}; original official DWG remains "
        "separate evidence and its knurl detail is intentionally excluded."
    )
    annotation.ObjectPlacement = target.ObjectPlacement
    obj.matrix_world = target_obj.matrix_world

    vertices = []
    edges = []
    for path in paths:
        path_indices = []
        for first, second in path:
            coordinates = tuple(value / 1000.0 for value in candidate_coordinates_mm(view, first, second))
            path_indices.append(len(vertices))
            vertices.append(coordinates)
        edges.extend(zip(path_indices, path_indices[1:]))
    obj.data.clear_geometry()
    obj.data.from_pydata(vertices, edges, [])
    obj.data.update()

    representation, path_count, edge_count = annotation_representation(
        model, context, view, paths
    )
    annotation.Representation = model.create_entity(
        "IfcProductDefinitionShape", Representations=[representation]
    )
    pset = ifcopenshell.api.pset.add_pset(
        model, product=annotation, name="EPset_Annotation"
    )
    ifcopenshell.api.pset.edit_pset(
        model,
        pset=pset,
        properties={
            "Classes": "review-target-gessi54294 native-dwg-review-simplification",
            "TargetGlobalId": TARGET_GLOBAL_ID,
            "ArticleNumber": ARTICLE,
            "DrawingProductCode": DRAWING_PRODUCT_CODE,
            "SourceKind": SOURCE_KIND,
            "SourceLabelZh": SOURCE_LABEL_ZH,
            "SourceDwgSha256": SOURCE_DWG_SHA256,
            "UnalteredOfficialDwg": False,
            "OriginalOfficialDwgEvidencePreserved": True,
            "HandleTextureDetailPathCount": 0,
            "IfcCurveStyleColour": BLUE,
        },
    )
    return annotation, representation, path_count, edge_count


def add_drawing(model, definition, room_bbox, include_elements, output_svg):
    minimum, maximum = room_bbox
    centre_x = (minimum[0] + maximum[0]) / 2
    centre_y = (minimum[1] + maximum[1]) / 2
    centre_z = 1.20
    target_view = definition["target_view"]
    if target_view == "PLAN_VIEW":
        storey = next(item for item in model.by_type("IfcBuildingStorey") if item.Name == "FFL")
        location_hint = storey.id()
        bpy.context.scene.cursor.location = (centre_x, centre_y, 1.20)
    elif definition["location_hint"] == "EAST":
        location_hint = "EAST"
        bpy.context.scene.cursor.location = (maximum[0] + 0.35, centre_y, centre_z)
    else:
        location_hint = "SOUTH"
        bpy.context.scene.cursor.location = (centre_x, minimum[1] - 0.35, centre_z)

    before = {item.id() for item in model.by_type("IfcAnnotation") if item.ObjectType == "DRAWING"}
    core_drawing.add_drawing(
        tool.Ifc, tool.Collector, tool.Drawing,
        target_view=target_view, location_hint=location_hint,
    )
    created = [
        item for item in model.by_type("IfcAnnotation")
        if item.ObjectType == "DRAWING" and item.id() not in before
    ]
    if len(created) != 1:
        raise RuntimeError(f"expected one drawing, got {len(created)}")
    drawing = created[0]
    core_drawing.update_drawing_name(
        tool.Ifc, tool.Drawing, drawing=drawing, name=definition["drawing_name"]
    )
    camera = tool.Ifc.get_object(drawing) or tool.Drawing.import_drawing(drawing)
    if target_view == "PLAN_VIEW":
        matrix = Matrix.Identity(4)
        matrix.translation = (centre_x, centre_y, 1.20)
        camera.matrix_world = matrix
        width = maximum[0] - minimum[0] + 0.60
        height = maximum[1] - minimum[1] + 0.60
        clip_end = 1.50
    elif definition["location_hint"] == "EAST":
        camera.matrix_world = tool.Drawing.generate_drawing_matrix("ELEVATION_VIEW", "EAST")
        width = maximum[1] - minimum[1] + 0.60
        height = 3.60
        clip_end = maximum[0] - minimum[0] + 0.80
    else:
        camera.matrix_world = tool.Drawing.generate_drawing_matrix("ELEVATION_VIEW", "SOUTH")
        width = maximum[0] - minimum[0] + 0.60
        height = 3.60
        clip_end = maximum[1] - minimum[1] + 0.80

    camera.data.type = "ORTHO"
    camera.data.clip_start = 0.002
    camera.data.clip_end = clip_end
    cprops = tool.Drawing.get_camera_props(camera)
    cprops.update_props = False
    cprops.camera_type = "ORTHO"
    cprops.target_view = target_view
    cprops.custom_scale_numerator = "1"
    cprops.custom_scale_denominator = "25"
    cprops.diagram_scale = "CUSTOM"
    cprops.has_underlay = False
    cprops.has_linework = True
    cprops.has_annotation = False
    cprops.linework_mode = "OPENCASCADE"
    cprops.fill_mode = "NONE"
    cprops.cut_mode = "BISECT"
    cprops.dpi = 300
    cprops.width = width
    cprops.height = height
    resolution_x, resolution_y = cprops.update_camera_resolution()
    cprops.update_props = True
    bpy.context.scene.render.resolution_x = resolution_x
    bpy.context.scene.render.resolution_y = resolution_y

    pset_data = ifcopenshell.util.element.get_pset(drawing, "EPset_Drawing")
    pset = model.by_id(pset_data["id"])
    ifcopenshell.api.pset.edit_pset(
        model,
        pset=pset,
        properties={
            "TargetView": target_view,
            "Scale": "1/25",
            "HumanScale": "1:25",
            "HasUnderlay": False,
            "HasLinework": True,
            "HasAnnotation": False,
            "GlobalReferencing": True,
            "DPI": 300,
            "LineworkMode": "OPENCASCADE",
            "FillMode": "NONE",
            "CutMode": "BISECT",
            "Include": ",".join(element.GlobalId for element in include_elements),
        },
    )
    drawing.Description = (
        f"Bonsai native Drawing / {ROOM_NAME} / Gessi {ARTICLE} / "
        f"{definition['drawing_name']} / {SOURCE_LABEL_ZH}"
    )
    reference = tool.Drawing.get_drawing_document(drawing)
    relative_output = os.path.relpath(
        output_svg, Path(tool.Ifc.get_path()).resolve().parent
    )
    ifcopenshell.api.document.edit_reference(
        model,
        reference=reference,
        attributes={"Location": Path(relative_output).as_posix()},
    )
    output_svg.parent.mkdir(parents=True, exist_ok=True)
    return drawing, camera, width, height, clip_end


def style_review_annotation(svg_path, annotation_global_id, view):
    raw_sha = sha256(svg_path)
    ET.register_namespace("", "http://www.w3.org/2000/svg")
    ET.register_namespace("ifc", "http://www.ifcopenshell.org/ns")
    tree = ET.parse(svg_path)
    root = tree.getroot()
    groups = []
    for element in root.iter():
        attributes = {local_name(key): value for key, value in element.attrib.items()}
        classes = element.attrib.get("class", "").split()
        if (
            attributes.get("guid") == annotation_global_id
            or f"GlobalId-{annotation_global_id}" in classes
            or "review-target-gessi54294" in classes
        ):
            groups.append(element)
    if not groups:
        raise RuntimeError("Bonsai SVG has no Gessi review annotation group")
    geometry_count = 0
    for group in groups:
        classes = group.attrib.get("class", "").split()
        for class_name in (
            "review-target-gessi54294",
            "native-dwg-review-simplification",
        ):
            if class_name not in classes:
                classes.append(class_name)
        group.attrib["class"] = " ".join(classes)
        group.attrib["data-source-kind"] = SOURCE_KIND
        group.attrib["data-source-label-zh"] = SOURCE_LABEL_ZH
        group.attrib["data-unaltered-official-dwg"] = "false"
        group.attrib["data-handle-texture-detail-path-count"] = "0"
        group.attrib["data-target-global-id"] = TARGET_GLOBAL_ID
        group.attrib["data-reviewed-view"] = view
        for element in group.iter():
            if local_name(element.tag) not in GEOMETRY_TAGS:
                continue
            geometry_count += 1
            existing = element.attrib.get("style", "").rstrip(";")
            highlight = f"stroke:{BLUE};stroke-width:0.35;fill:none"
            element.attrib["style"] = f"{existing};{highlight}" if existing else highlight
    root.attrib["data-create-drawing-result"] = "FINISHED"
    root.attrib["data-review-source-geometry-modified-after-create-drawing"] = "false"
    tree.write(svg_path, encoding="utf-8", xml_declaration=True)
    return {
        "bonsai_generated_sha256_before_review_style": raw_sha,
        "annotation_group_count": len(groups),
        "annotation_geometry_element_count": geometry_count,
        "post_style_only": True,
    }


def svg_numbers(value):
    return [float(item) for item in re.findall(r"[-+]?(?:\d*\.\d+|\d+)(?:[eE][-+]?\d+)?", value)]


def svg_group_x_coordinates(group):
    coordinates = []
    for element in group.iter():
        tag = local_name(element.tag)
        if tag == "line":
            coordinates.extend(float(element.attrib[key]) for key in ("x1", "x2"))
        elif tag in {"polyline", "polygon"}:
            coordinates.extend(svg_numbers(element.attrib.get("points", ""))[::2])
        elif tag == "path":
            coordinates.extend(svg_numbers(element.attrib.get("d", ""))[::2])
    return coordinates


def inspect_svg(
    svg_path,
    target_global_id,
    annotation_global_id,
    finished_wall_global_id=None,
    occluding_context_global_ids=None,
):
    root = ET.parse(svg_path).getroot()
    classes = Counter()
    geometry_count = 0
    target_groups = 0
    annotation_geometry = 0
    annotation_x_coordinates = []
    finished_wall_groups = []
    occluding_context_groups = []
    for element in root.iter():
        attributes = {local_name(key): value for key, value in element.attrib.items()}
        element_classes = element.attrib.get("class", "").split()
        classes.update(element_classes)
        if local_name(element.tag) in GEOMETRY_TAGS:
            geometry_count += 1
        if attributes.get("guid") == target_global_id:
            target_groups += 1
        if (
            f"GlobalId-{annotation_global_id}" in element_classes
            and local_name(element.tag) in GEOMETRY_TAGS
        ):
            annotation_geometry += 1
            if local_name(element.tag) == "line":
                annotation_x_coordinates.extend(
                    float(element.attrib[key]) for key in ("x1", "x2")
                )
        if attributes.get("guid") == finished_wall_global_id:
            finished_wall_groups.append(element)
        if attributes.get("guid") in (occluding_context_global_ids or set()):
            occluding_context_groups.append(element)
    finished_wall_x_coordinates = [
        value for group in finished_wall_groups for value in svg_group_x_coordinates(group)
    ]
    return {
        "root_tag": local_name(root.tag),
        "width": root.attrib.get("width"),
        "height": root.attrib.get("height"),
        "viewBox": root.attrib.get("viewBox"),
        "data_scale": root.attrib.get("data-scale"),
        "geometry_element_count": geometry_count,
        "projection_group_count": classes.get("projection", 0),
        "class_counts": dict(sorted(classes.items())),
        "target_ifc_projection_group_count": target_groups,
        "review_annotation_geometry_count": annotation_geometry,
        "finished_wall_projection_group_count": len(finished_wall_groups),
        "occluding_context_projection_group_count": len(occluding_context_groups),
        "finished_wall_x_coordinates_svg_mm": finished_wall_x_coordinates,
        "review_annotation_x_coordinates_svg_mm": annotation_x_coordinates,
    }


def bounds_3d(points):
    return [
        [min(point[axis] for point in points) for axis in range(3)],
        [max(point[axis] for point in points) for axis in range(3)],
    ]


def main():
    arguments = sys.argv[sys.argv.index("--") + 1 :]
    if len(arguments) != 4:
        raise SystemExit("expected: drawing-session.ifc formal.ifc VIEW output-directory")
    session_ifc = Path(arguments[0]).resolve()
    formal_ifc = Path(arguments[1]).resolve()
    view = arguments[2]
    output_dir = Path(arguments[3]).resolve()
    if view not in VIEW_DEFINITIONS:
        raise RuntimeError(f"unsupported view: {view}")
    if sha256(formal_ifc) != FORMAL_SHA256 or sha256(FORMAL_IFC) != FORMAL_SHA256:
        raise RuntimeError("formal IFC hash mismatch")
    if sha256(session_ifc) != FORMAL_SHA256:
        raise RuntimeError("fresh drawing session must start byte-identical to formal IFC")
    session_before_sha = sha256(session_ifc)

    candidate = json.loads(CANDIDATE.read_text(encoding="utf-8"))
    candidate_view = candidate["views"][view]
    paths = candidate_view["review_simplified_official_outline_paths_mm"]
    audit = candidate_view["handle_line_texture_simplification"]
    if (
        candidate.get("source_kind") != SOURCE_KIND
        or candidate_view.get("source_kind") != SOURCE_KIND
        or candidate_view.get("unaltered_official_dwg") is not False
        or candidate_view.get("source_dwg_sha256") != SOURCE_DWG_SHA256
        or len(paths) != EXPECTED_PATH_COUNTS[view]
        or audit.get("review_texture_detail_path_count") != 0
        or audit.get("handle_envelope_delta_mm") != [0.0, 0.0]
        or audit.get("three_hole_and_spout_positions_unchanged") is not True
        or audit.get("pass") is not True
    ):
        raise RuntimeError("Gessi de-textured review candidate gate failed")

    os.chdir(ROOT)
    with contextlib.suppress(Exception):
        addon_utils.disable("bl_ext.user_default.project_control", default_set=False, handle_error=None)
    load_result = bpy.ops.bim.load_project(
        filepath=str(session_ifc),
        should_start_fresh_session=True,
        use_detailed_tooltip=True,
    )
    if load_result != {"FINISHED"} or not tool.Ifc.get():
        raise RuntimeError(f"Bonsai failed to load session: {load_result}")
    model = tool.Ifc.get()
    target = model.by_guid(TARGET_GLOBAL_ID)
    room = model.by_guid(ROOM_GLOBAL_ID)
    target_obj = tool.Ifc.get_object(target)
    room_obj = tool.Ifc.get_object(room)
    if target is None or room is None or target_obj is None or room_obj is None:
        raise RuntimeError("Gessi target or main-bathroom space missing")
    if room.LongName != ROOM_NAME:
        raise RuntimeError("main-bathroom room identity drifted")

    room_bbox = world_bbox(room_obj)
    target_bbox = world_bbox(target_obj)
    records = room_elements(room_bbox)
    original_context = [record[0] for record in records]
    if original_context.count(target) != 1:
        raise RuntimeError("expected Gessi target exactly once in room context")
    context_elements = [element for element in original_context if element != target]
    finished_wall = model.by_guid(FINISHED_WALL_COVERING_GLOBAL_ID)
    finished_wall_obj = tool.Ifc.get_object(finished_wall) if finished_wall else None
    if finished_wall_obj is None or finished_wall not in context_elements:
        raise RuntimeError("main-bathroom finished wall covering is missing from context")
    if view == "side":
        available_context_ids = {element.GlobalId for element in context_elements}
        if not SIDE_OCCLUDING_CONTEXT_GLOBAL_IDS.issubset(available_context_ids):
            raise RuntimeError("expected all south elevation occluders before Side filtering")
        context_elements = [
            element
            for element in context_elements
            if element.GlobalId not in SIDE_OCCLUDING_CONTEXT_GLOBAL_IDS
        ]
    output_svg = output_dir / f"{VIEW_DEFINITIONS[view]['drawing_name']}.svg"
    evidence_path = output_dir / f"{VIEW_DEFINITIONS[view]['drawing_name']}-create-drawing-evidence.json"
    drawing, camera, width, height, clip_end = add_drawing(
        model, VIEW_DEFINITIONS[view], room_bbox, context_elements, output_svg
    )
    activate_result = bpy.ops.bim.activate_drawing(
        drawing=drawing.id(), should_view_from_camera=False
    )
    if activate_result != {"FINISHED"}:
        raise RuntimeError(f"failed to activate drawing: {activate_result}")
    annotation, representation, path_count, edge_count = add_review_annotation(
        model, drawing, target, target_obj, view, paths
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
    if create_result != {"FINISHED"} or not output_svg.is_file() or output_svg.stat().st_size == 0:
        raise RuntimeError(f"Bonsai Create Drawing failed: {create_result}")
    style = style_review_annotation(output_svg, annotation.GlobalId, view)
    svg = inspect_svg(
        output_svg,
        TARGET_GLOBAL_ID,
        annotation.GlobalId,
        FINISHED_WALL_COVERING_GLOBAL_ID if view == "side" else None,
        SIDE_OCCLUDING_CONTEXT_GLOBAL_IDS if view == "side" else None,
    )
    if (
        svg["root_tag"] != "svg"
        or svg["geometry_element_count"] == 0
        or svg["projection_group_count"] == 0
        or svg["target_ifc_projection_group_count"] != 0
        or svg["review_annotation_geometry_count"] == 0
    ):
        raise RuntimeError(f"Bonsai SVG structure/layer gate failed: {svg}")
    side_scene_wall_gate = None
    if view == "side":
        finished_wall_bbox = world_bbox(finished_wall_obj)
        finished_wall_face_x_m = finished_wall_bbox[1][0]
        target_body_mount_face_x_m = target_bbox[0][0]
        local_wall_anchor_y_mm = candidate[
            "review_only_geometry_simplification"
        ]["side_mechanical_gate"]["official_mapped_targets"]["wall_surface_y_mm"]
        annotation_wall_anchor_world = target_obj.matrix_world @ Vector(
            (0.0, local_wall_anchor_y_mm / 1000.0, 0.0)
        )
        annotation_to_finished_wall_residual_mm = (
            annotation_wall_anchor_world.x - finished_wall_face_x_m
        ) * 1000.0
        actual_body_to_finished_wall_residual_mm = (
            target_body_mount_face_x_m - finished_wall_face_x_m
        ) * 1000.0
        finished_wall_x_coordinates = svg["finished_wall_x_coordinates_svg_mm"]
        annotation_x_coordinates = svg["review_annotation_x_coordinates_svg_mm"]
        if not finished_wall_x_coordinates or not annotation_x_coordinates:
            raise RuntimeError("Side SVG lacks finished-wall or review Annotation coordinates")
        finished_wall_face_svg_x_mm = max(finished_wall_x_coordinates)
        annotation_wall_anchor_svg_x_mm = min(
            annotation_x_coordinates,
            key=lambda value: abs(value - finished_wall_face_svg_x_mm),
        )
        svg_residual_mm = (
            annotation_wall_anchor_svg_x_mm - finished_wall_face_svg_x_mm
        ) * 25.0
        side_scene_wall_gate = {
            "tolerance_mm": 0.2,
            "coordinate_space": "actual project world coordinates and generated Side SVG",
            "finished_wall_covering_global_id": FINISHED_WALL_COVERING_GLOBAL_ID,
            "finished_wall_face_x_m": finished_wall_face_x_m,
            "occluding_south_context_global_ids": sorted(SIDE_OCCLUDING_CONTEXT_GLOBAL_IDS),
            "occluding_south_context_excluded_from_side_drawing": True,
            "annotation_rigid_translation_applied_mm": [0.0, 0.0, 0.0],
            "annotation_wall_anchor_world_x_m": annotation_wall_anchor_world.x,
            "annotation_to_finished_wall_residual_mm": annotation_to_finished_wall_residual_mm,
            "actual_body_mount_face_x_m": target_body_mount_face_x_m,
            "actual_body_to_finished_wall_residual_mm": actual_body_to_finished_wall_residual_mm,
            "finished_wall_face_svg_x_mm": finished_wall_face_svg_x_mm,
            "annotation_wall_anchor_svg_x_mm": annotation_wall_anchor_svg_x_mm,
            "svg_residual_mm": svg_residual_mm,
            "finished_wall_projection_group_count": svg["finished_wall_projection_group_count"],
            "occluding_context_projection_group_count": svg["occluding_context_projection_group_count"],
            "root_cause": "the 20 mm wall finish was hidden by the south boundary wall layers in the exterior-side projection; the Annotation already matched the finished face",
            "pass": (
                abs(annotation_to_finished_wall_residual_mm) <= 0.2
                and abs(actual_body_to_finished_wall_residual_mm) <= 0.2
                and abs(svg_residual_mm) <= 0.2
                and svg["finished_wall_projection_group_count"] >= 1
                and svg["occluding_context_projection_group_count"] == 0
            ),
        }
        if side_scene_wall_gate["pass"] is not True:
            raise RuntimeError(f"Side scene finished-wall gate failed: {side_scene_wall_gate}")
    cache_path = output_svg.parent / "cache" / f"{output_svg.stem}-linework.svg"
    if not cache_path.is_file() or cache_path.stat().st_size == 0:
        raise RuntimeError("Bonsai linework cache missing")

    camera_record = {
        "type": camera.data.type,
        "width_m": width,
        "height_m": height,
        "clip_start_m": camera.data.clip_start,
        "clip_end_m": clip_end,
        "resolution": [
            bpy.context.scene.render.resolution_x,
            bpy.context.scene.render.resolution_y,
        ],
    }

    temporary = session_ifc.with_suffix(".ifc.next")
    model.write(str(temporary))
    os.replace(temporary, session_ifc)
    persisted_sha = sha256(session_ifc)
    reload_result = bpy.ops.bim.load_project(
        filepath=str(session_ifc),
        should_start_fresh_session=True,
        use_detailed_tooltip=True,
    )
    if reload_result != {"FINISHED"} or not tool.Ifc.get():
        raise RuntimeError(f"Bonsai reload failed: {reload_result}")
    reopened = tool.Ifc.get()
    reopened_drawing = reopened.by_guid(drawing.GlobalId)
    reopened_annotation = reopened.by_guid(annotation.GlobalId)
    if reopened_drawing is None or reopened_annotation is None:
        raise RuntimeError("reloaded session lost Drawing or LINEWORK Annotation")
    reopened_representation = next(
        (
            item for item in reopened_annotation.Representation.Representations
            if item.RepresentationIdentifier == "Annotation"
            and item.RepresentationType == "GeometricCurveSet"
            and item.ContextOfItems.TargetView == VIEW_DEFINITIONS[view]["target_view"]
        ),
        None,
    )
    if reopened_representation is None:
        raise RuntimeError("reloaded session lost review representation")
    persisted_paths = [
        path
        for item in reopened_representation.Items
        for path in item.Elements
    ]
    persisted_points = [
        tuple(float(value) for value in point.Coordinates)
        for path in persisted_paths
        for point in path.Points
    ]
    expected_points = [
        candidate_coordinates_mm(view, first, second)
        for path in paths
        for first, second in path
    ]
    persisted_bounds = bounds_3d(persisted_points)
    expected_bounds = bounds_3d(expected_points)
    coordinate_residual = max(
        abs(persisted_bounds[bound][axis] - expected_bounds[bound][axis])
        for bound in range(2)
        for axis in range(3)
    )
    if len(persisted_paths) != EXPECTED_PATH_COUNTS[view] or coordinate_residual > 0.000001:
        raise RuntimeError("persisted review coordinates drifted")
    annotation_pset = ifcopenshell.util.element.get_pset(
        reopened_annotation, "EPset_Annotation"
    )
    styled_items = [
        styled
        for item in reopened_representation.Items
        for styled in item.StyledByItem
    ]
    if (
        annotation_pset.get("HandleTextureDetailPathCount") != 0
        or annotation_pset.get("UnalteredOfficialDwg") is not False
        or not styled_items
    ):
        raise RuntimeError("persisted review metadata/style gate failed")

    drawing_pset = ifcopenshell.util.element.get_pset(reopened_drawing, "EPset_Drawing")
    include = drawing_pset.get("Include", "").split(",")
    if TARGET_GLOBAL_ID in include or len(include) != len(context_elements):
        raise RuntimeError("Drawing Include duplicate-target policy failed")
    if view == "side" and (
        FINISHED_WALL_COVERING_GLOBAL_ID not in include
        or SIDE_OCCLUDING_CONTEXT_GLOBAL_IDS.intersection(include)
    ):
        raise RuntimeError("Side Drawing finished-wall filter did not persist")
    reopened_target = reopened.by_guid(TARGET_GLOBAL_ID)
    reopened_annotation_matrix = ifcopenshell.util.placement.get_local_placement(
        reopened_annotation.ObjectPlacement
    )
    reopened_target_matrix = ifcopenshell.util.placement.get_local_placement(
        reopened_target.ObjectPlacement
    )
    placement_residual = max(
        abs(float(reopened_annotation_matrix[row][column]) - float(reopened_target_matrix[row][column]))
        for row in range(4)
        for column in range(4)
    )
    if placement_residual > 0.000001:
        raise RuntimeError("reloaded Annotation placement drifted from the target placement")
    if side_scene_wall_gate is not None:
        side_scene_wall_gate["post_reload_annotation_target_placement_residual"] = placement_residual
    class_counts = Counter(element.is_a() for element in context_elements)
    mechanical_gate = candidate["review_only_geometry_simplification"]["side_mechanical_gate"]
    evidence = {
        "schema_version": 1,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "task": "Bonsai Create Drawing for Gessi 45089_54294 de-textured official-outline review in actual main-bathroom context",
        "workflow": ["inspect", "plan", "execute", "persist", "reload", "verify"],
        "provider": {
            "preferred": "bonsai-mcp 1.1.0",
            "status": "supported_for_inspection",
            "inspection_project": "IFC4 / My Project",
            "inspection_note": "The active Provider held another product-level session, so isolated mutation used a separate local Blender process.",
            "execution": "local Blender 4.5.3 LTS with installed Bonsai/IfcOpenShell 0.8.4",
        },
        "course_evidence": {
            "query": "Bonsai create drawing SVG PDF annotation linework save reload verify",
            "top_lesson": "085000 Introduction to Drawings",
            "course_fact": "Create Drawing generates or refreshes SVG after a Drawing camera, scale, depth and filters are configured.",
            "source_files_packaged": False,
            "source_files_note": "The installed skill contains the indexed lesson records but not the referenced private source Markdown or screenshots.",
        },
        "view": view,
        "formal_ifc": str(formal_ifc),
        "formal_ifc_sha256": sha256(formal_ifc),
        "formal_ifc_bytes_unchanged": sha256(formal_ifc) == FORMAL_SHA256,
        "drawing_session_ifc": str(session_ifc),
        "drawing_session_sha256_before": session_before_sha,
        "drawing_session_sha256_after": persisted_sha,
        "drawing_session_is_safe_per_view_copy": True,
        "source_kind": SOURCE_KIND,
        "source_label_zh": SOURCE_LABEL_ZH,
        "source_dwg_sha256": SOURCE_DWG_SHA256,
        "unaltered_official_dwg": False,
        "original_official_dwg_evidence_preserved": True,
        "handle_texture_detail_path_count": 0,
        "review_path_count": EXPECTED_PATH_COUNTS[view],
        "persisted_review_path_count": len(persisted_paths),
        "review_edge_count": edge_count,
        "persisted_coordinate_bounds_mm": persisted_bounds,
        "expected_coordinate_bounds_mm": expected_bounds,
        "persisted_coordinate_residual_mm": coordinate_residual,
        "mechanical_gate": {
            "tolerance_mm": mechanical_gate["tolerance_mm"],
            "wall_anchor_residual_mm": mechanical_gate["review_proxy_residual_mm"]["wall_surface_y_mm"],
            "straight_spout_axis_residual_mm": mechanical_gate["review_proxy_residual_mm"]["straight_spout_axis_z_mm"],
            "lower_lip_residual_mm": mechanical_gate["review_proxy_residual_mm"]["lower_lip_z_mm"],
            "outlet_endpoint_residual_mm": mechanical_gate["review_proxy_residual_mm"]["outlet_endpoint_y_mm"],
            "visible_reach_residual_mm": mechanical_gate["review_proxy_residual_mm"]["visible_reach_mm"],
            "visible_envelope_residual_mm": mechanical_gate["review_proxy_residual_mm"]["visible_envelope_mm"],
            "pass": mechanical_gate["feature_gate_pass"],
        },
        "side_scene_finished_wall_gate": side_scene_wall_gate,
        "target": {
            "global_id": TARGET_GLOBAL_ID,
            "ifc_class": target.is_a(),
            "type_name": next(relation.RelatingType.Name for relation in target.IsTypedBy),
            "world_bbox_m": [list(target_bbox[0]), list(target_bbox[1])],
            "target_include_count_before_suppression": original_context.count(target),
            "target_include_count_after_suppression": include.count(TARGET_GLOBAL_ID),
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
            "include_ifc_class_counts": dict(sorted(class_counts.items())),
        },
        "drawing": {
            "global_id": drawing.GlobalId,
            "name": drawing.Name,
            "object_type": drawing.ObjectType,
            "epset_drawing": drawing_pset,
        },
        "linework_annotation": {
            "global_id": annotation.GlobalId,
            "ifc_class": annotation.is_a(),
            "name": annotation.Name,
            "predefined_type": ifcopenshell.util.element.get_predefined_type(annotation),
            "epset_annotation": annotation_pset,
            "ifc_curve_style_persisted": True,
            "ifc_curve_style_colour": BLUE,
        },
        "camera": camera_record,
        "create_drawing": {
            "operator": "bpy.ops.bim.create_drawing",
            "arguments": {"print_all": False, "open_viewer": False, "sync": False},
            "result": sorted(create_result),
            "linework_mode": "OPENCASCADE",
            "target_view": VIEW_DEFINITIONS[view]["target_view"],
        },
        "persistence": {
            "method": "model.write to temporary then atomic replace; bpy.ops.bim.load_project reload",
            "reload_result": sorted(reload_result),
            "post_reload_drawing_found": reopened_drawing is not None,
            "post_reload_annotation_found": reopened_annotation is not None,
        },
        "svg": {
            "path": str(output_svg),
            "bytes": output_svg.stat().st_size,
            "sha256": sha256(output_svg),
            **style,
            **svg,
        },
        "linework_cache": {
            "path": str(cache_path),
            "bytes": cache_path.stat().st_size,
            "sha256": sha256(cache_path),
        },
        "versions": {
            "blender": bpy.app.version_string,
            "ifcopenshell": ifcopenshell.version,
            "ifc_schema": reopened.schema,
            "bonsai_generator": "Bonsai 0.8.4 bim.create_drawing",
        },
        "pass": True,
    }
    evidence_path.write_text(
        json.dumps(evidence, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps({
        "view": view,
        "drawing": drawing.Name,
        "create_drawing_result": sorted(create_result),
        "svg": str(output_svg),
        "svg_sha256": sha256(output_svg),
        "session_ifc": str(session_ifc),
        "evidence": str(evidence_path),
        "pass": True,
    }, ensure_ascii=False))


if __name__ == "__main__":
    main()
