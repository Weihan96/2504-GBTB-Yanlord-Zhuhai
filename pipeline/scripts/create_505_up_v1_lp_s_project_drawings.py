#!/usr/bin/env python3
"""Create persisted 505 UP Bonsai Drawings in the actual entrance space.

The input is a product-level copy of the approved derived IFC.  The target
Body is excluded from each Drawing Include list and replaced by the approved
black geometry-derived proxy as a persistent IFC LINEWORK Annotation.
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
import ifcopenshell.util.placement
from bonsai import tool
from bonsai.core import drawing as core_drawing
from mathutils import Vector


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "pipeline/scripts"))
import create_gessi316_54294_main_bathroom_drawing as shared  # noqa: E402


PRODUCT_DIR = ROOT / "output/review/highpoly-types/505-up-v1-lp-s"
FORMAL_IFC = ROOT / "2504 GBTB Yanlord Zhuhai.ifc"
FORMAL_SHA256 = "7a50b87e8f48a7c2bfbaa9f1dabdfca7155fa684b2325c66aa1ae15c8ab0a25c"
APPROVAL = ROOT / "pipeline/decisions/505-up-v1-lp-s-drawing-approval.json"
CANDIDATE = PRODUCT_DIR / "candidate-representations.json"
MANIFEST = PRODUCT_DIR / "manifest.json"
TARGET_GLOBAL_ID = "19MpdkWqXC7uhUNhLQgrce"
ROOM_GLOBAL_ID = "36NNomu9TDP9xkhh6D2l8K"
ROOM_NAME = "玄关"
SOURCE_KIND = "geometry_derived_simplified_proxy"
SOURCE_LABEL_ZH = "基于原始高模几何生成的简化图纸表达"
BLACK = "#111111"
GREY = "#a3abb3"
EXPECTED_PATH_COUNTS = {"plan": 20, "front": 281, "side": 93}
VIEW_DEFINITIONS = {
    "plan": {
        "drawing_name": "MOLTENI-505-UP-ENTRANCE-PLAN",
        "target_view": "PLAN_VIEW",
        "location_hint": "PLAN",
    },
    "front": {
        "drawing_name": "MOLTENI-505-UP-ENTRANCE-FRONT",
        "target_view": "ELEVATION_VIEW",
        "location_hint": "NORTH",
    },
    "side": {
        "drawing_name": "MOLTENI-505-UP-ENTRANCE-SIDE",
        "target_view": "ELEVATION_VIEW",
        "location_hint": "EAST",
    },
}
GEOMETRY_TAGS = {"path", "polyline", "polygon", "line", "circle", "ellipse", "rect"}


shared.PRODUCT_DIR = PRODUCT_DIR
shared.FORMAL_IFC = FORMAL_IFC
shared.FORMAL_SHA256 = FORMAL_SHA256
shared.TARGET_GLOBAL_ID = TARGET_GLOBAL_ID
shared.ROOM_GLOBAL_ID = ROOM_GLOBAL_ID
shared.ROOM_NAME = ROOM_NAME
shared.ARTICLE = "505 UP V1.LP.S"
shared.DRAWING_PRODUCT_CODE = "505-UP"
shared.SOURCE_KIND = SOURCE_KIND
shared.SOURCE_LABEL_ZH = SOURCE_LABEL_ZH
shared.BLUE = BLACK
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


def union_bbox(first, second):
    return (
        tuple(min(first[0][axis], second[0][axis]) for axis in range(3)),
        tuple(max(first[1][axis], second[1][axis]) for axis in range(3)),
    )


def vector_record(value):
    return [round(float(component), 9) for component in value]


def matrix_record(value):
    return [
        [round(float(component), 9) for component in row]
        for row in value
    ]


def orient_front_camera_north(model, drawing, camera, target_obj, crop_bbox):
    """Place Front inside the room, looking south at the cabinet front face."""
    before = camera.matrix_world.copy()
    minimum, maximum = crop_bbox
    centre_x = (minimum[0] + maximum[0]) / 2.0
    bpy.context.scene.cursor.location = (centre_x, maximum[1] + 0.35, 1.20)
    camera.matrix_world = tool.Drawing.generate_drawing_matrix(
        "ELEVATION_VIEW", "NORTH"
    )
    bpy.context.view_layer.update()
    ifcopenshell.api.geometry.edit_object_placement(
        model,
        product=drawing,
        matrix=camera.matrix_world,
        is_si=True,
        should_transform_children=False,
    )
    target_axes = target_obj.matrix_world.to_3x3()
    before_axes = before.to_3x3()
    after_axes = camera.matrix_world.to_3x3()
    local_x_world = target_axes @ Vector((1.0, 0.0, 0.0))
    front_outward_world = target_axes @ Vector((0.0, 1.0, 0.0))
    screen_right_before = before_axes @ Vector((1.0, 0.0, 0.0))
    screen_right_after = after_axes @ Vector((1.0, 0.0, 0.0))
    view_direction_before = before_axes @ Vector((0.0, 0.0, -1.0))
    view_direction_after = after_axes @ Vector((0.0, 0.0, -1.0))
    audit = {
        "method": "product_world_axes_plus_room_facing_component_depth_plus_camera_basis",
        "camera_before_matrix_m": matrix_record(before),
        "camera_after_matrix_m": matrix_record(camera.matrix_world),
        "target_local_x_world": vector_record(local_x_world),
        "target_front_outward_world": vector_record(front_outward_world),
        "screen_right_before_world": vector_record(screen_right_before),
        "screen_right_after_world": vector_record(screen_right_after),
        "view_direction_before_world": vector_record(view_direction_before),
        "view_direction_after_world": vector_record(view_direction_after),
        "local_x_screen_dot_before": round(local_x_world.dot(screen_right_before), 9),
        "local_x_screen_dot_after": round(local_x_world.dot(screen_right_after), 9),
        "front_outward_view_dot_before": round(
            front_outward_world.dot(view_direction_before), 9
        ),
        "front_outward_view_dot_after": round(
            front_outward_world.dot(view_direction_after), 9
        ),
        "mechanical_anchors": {
            "back_wall_world_y_m": [-1.6, -1.4],
            "cabinet_world_y_m": [-1.4, -0.9635005],
            "entrance_room_world_y_m": [-1.5, 0.1],
            "front_components_local_y_mm": {
                "slats": [324.0, 349.0],
                "display_maximum": 420.499451,
            },
        },
        "original_pre_state_evidence": {
            "camera_matrix_m": [
                [1.0, 0.0, 0.0, 4.702590466],
                [0.0, 0.0, -1.0, -1.850000024],
                [0.0, 1.0, 0.0, 1.200000048],
                [0.0, 0.0, 0.0, 1.0],
            ],
            "project_front_svg_sha256": "f55f2ac9d385e2df6bfd0d1ff2df41451725d34db0afa66cfdc0e80a888f5165",
            "interpretation": "south camera looked north at the cabinet back and reversed local X",
        },
        "post_state": "north camera looks south at the room-facing cabinet front; local +X is screen right",
    }
    audit["input_was_mirrored_south_camera"] = (
        audit["local_x_screen_dot_before"] < -0.999
        and view_direction_before.y > 0.999
    )
    audit["input_was_already_corrected_north_camera"] = (
        audit["local_x_screen_dot_before"] > 0.999
        and view_direction_before.y < -0.999
    )
    if (
        not (
            audit["input_was_mirrored_south_camera"]
            or audit["input_was_already_corrected_north_camera"]
        )
        or audit["local_x_screen_dot_after"] <= 0.999
        or view_direction_after.y >= -0.999
    ):
        raise RuntimeError(f"505 Front handedness gate failed: {audit}")
    audit["pass"] = True
    return audit


def add_505_drawing(model, definition, crop_bbox, include_elements, output_svg):
    drawing, camera, width, height, clip_end = shared.add_drawing(
        model, definition, crop_bbox, include_elements, output_svg
    )
    if definition["location_hint"] == "NORTH":
        target_obj = tool.Ifc.get_object(model.by_guid(TARGET_GLOBAL_ID))
        orient_front_camera_north(model, drawing, camera, target_obj, crop_bbox)
    return drawing, camera, width, height, clip_end


def view3d_override():
    """Return a real GUI VIEW_3D context for Drawing operators called by MCP."""
    for window in bpy.context.window_manager.windows:
        for area in window.screen.areas:
            if area.type != "VIEW_3D":
                continue
            region = next((item for item in area.regions if item.type == "WINDOW"), None)
            if region is not None:
                return {
                    "window": window,
                    "screen": window.screen,
                    "area": area,
                    "region": region,
                    "scene": bpy.context.scene,
                }
    raise RuntimeError("Bonsai Drawing requires a real VIEW_3D area")


def annotation_representation(model, context, view, paths):
    polylines = []
    for path in paths:
        points = [
            model.create_entity(
                "IfcCartesianPoint",
                Coordinates=coordinates_mm(view, first, second),
            )
            for first, second in path
        ]
        if len(points) >= 2:
            polylines.append(model.create_entity("IfcPolyline", Points=points))
    if len(polylines) != EXPECTED_PATH_COUNTS[view]:
        raise RuntimeError(f"{view} approved proxy path count drifted")
    curve_set = model.create_entity("IfcGeometricCurveSet", Elements=polylines)
    colour = model.create_entity(
        "IfcColourRgb",
        Name="Approved geometry-derived black",
        Red=17.0 / 255.0,
        Green=17.0 / 255.0,
        Blue=17.0 / 255.0,
    )
    curve_style = model.create_entity(
        "IfcCurveStyle",
        Name="Molteni 505 UP approved black simplified linework",
        CurveFont=None,
        CurveWidth=model.create_entity("IfcPositiveLengthMeasure", 0.35),
        CurveColour=colour,
        ModelOrDraughting=True,
    )
    model.create_entity(
        "IfcStyledItem",
        Item=curve_set,
        Styles=[curve_style],
        Name="Black approved geometry-derived LINEWORK",
    )
    representation = model.create_entity(
        "IfcShapeRepresentation",
        ContextOfItems=context,
        RepresentationIdentifier="Annotation",
        RepresentationType="GeometricCurveSet",
        Items=[curve_set],
    )
    return representation, len(polylines), sum(max(0, len(path) - 1) for path in paths)


def add_proxy_annotation(model, drawing, target, target_obj, view, paths):
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
    annotation.Name = f"Molteni 505 UP approved black proxy / {view}"
    annotation.Description = (
        f"{SOURCE_LABEL_ZH}; target {TARGET_GLOBAL_ID}; official 505 UP family DWGs "
        "remain identity/component references and are not inserted as project linework."
    )
    annotation.ObjectPlacement = target.ObjectPlacement
    obj.matrix_world = target_obj.matrix_world

    vertices = []
    edges = []
    for path in paths:
        indices = []
        for first, second in path:
            indices.append(len(vertices))
            vertices.append(tuple(value / 1000.0 for value in coordinates_mm(view, first, second)))
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
    pset = ifcopenshell.api.pset.add_pset(
        model, product=annotation, name="EPset_Annotation"
    )
    ifcopenshell.api.pset.edit_pset(
        model,
        pset=pset,
        properties={
            "Classes": "review-target-505-up geometry-derived-simplified-proxy",
            "TargetGlobalId": TARGET_GLOBAL_ID,
            "SourceKind": SOURCE_KIND,
            "SourceLabelZh": SOURCE_LABEL_ZH,
            "OfficialCadUsed": False,
            "BlueOfficialAtomicComponentCompositionSelected": False,
            "IfcCurveStyleColour": BLACK,
        },
    )
    return annotation, path_count, edge_count


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
            or "review-target-505-up" in classes
        ):
            annotation_groups.append(element)
    if not annotation_groups:
        raise RuntimeError("Bonsai SVG has no 505 UP LINEWORK Annotation group")
    black_ids = {
        id(element)
        for group in annotation_groups
        for element in group.iter()
        if local_name(element.tag) in GEOMETRY_TAGS
    }
    black_count = 0
    grey_count = 0
    for element in root.iter():
        if local_name(element.tag) not in GEOMETRY_TAGS:
            continue
        existing = element.attrib.get("style", "").rstrip(";")
        if id(element) in black_ids:
            style = f"stroke:{BLACK};stroke-width:0.35;fill:none"
            black_count += 1
        else:
            style = f"stroke:{GREY};stroke-width:0.22;fill:none;stroke-opacity:0.72"
            grey_count += 1
        element.attrib["style"] = f"{existing};{style}" if existing else style
    for group in annotation_groups:
        classes = group.attrib.get("class", "").split()
        for class_name in ("review-target-505-up", SOURCE_KIND):
            if class_name not in classes:
                classes.append(class_name)
        group.attrib["class"] = " ".join(classes)
        group.attrib["data-source-kind"] = SOURCE_KIND
        group.attrib["data-source-label-zh"] = SOURCE_LABEL_ZH
        group.attrib["data-official-cad-used"] = "false"
        group.attrib["data-target-global-id"] = TARGET_GLOBAL_ID
        group.attrib["data-reviewed-view"] = view
    root.attrib["data-create-drawing-result"] = "FINISHED"
    root.attrib["data-context-colour"] = GREY
    root.attrib["data-approved-linework-colour"] = BLACK
    root.attrib["data-target-body-projection-excluded"] = "true"
    tree.write(svg_path, encoding="utf-8", xml_declaration=True)
    if black_count == 0 or grey_count == 0:
        raise RuntimeError("505 UP SVG colour-layer gate failed")
    return {
        "bonsai_generated_sha256_before_review_style": raw_sha,
        "annotation_group_count": len(annotation_groups),
        "black_geometry_element_count": black_count,
        "grey_context_geometry_element_count": grey_count,
        "post_style_only": True,
    }


def refresh_plan_only(session_ifc: Path, formal_ifc: Path, output_dir: Path):
    """Replace the persisted Plan LINEWORK and run Create Drawing for Plan only."""
    if sha256(formal_ifc) != FORMAL_SHA256 or sha256(FORMAL_IFC) != FORMAL_SHA256:
        raise RuntimeError("formal IFC hash mismatch")
    session_before_sha = sha256(session_ifc)
    approval = json.loads(APPROVAL.read_text(encoding="utf-8"))
    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    candidate = json.loads(CANDIDATE.read_text(encoding="utf-8"))
    if (
        approval.get("status") != "approved"
        or approval.get("derived_ifc_write_allowed") is not True
        or approval.get("formal_authoritative_ifc_write_allowed") is not False
        or approval.get("candidate_manifest_sha256") != sha256(MANIFEST)
        or candidate.get("source_kind") != SOURCE_KIND
        or candidate.get("official_cad_used") is not False
        or manifest.get("formal_ifc_sha256") != FORMAL_SHA256
    ):
        raise RuntimeError("505 UP approval/source gate failed")
    plan_paths = candidate["views"]["plan"]["proxy_paths_mm"]
    repair_audit = candidate["views"]["plan"].get("component_boundary_extraction")
    if len(plan_paths) != EXPECTED_PATH_COUNTS["plan"] or not repair_audit:
        raise RuntimeError("505 UP repaired Plan candidate gate failed")

    output_dir.mkdir(parents=True, exist_ok=True)
    front_svg = output_dir / f"{VIEW_DEFINITIONS['front']['drawing_name']}.svg"
    side_svg = output_dir / f"{VIEW_DEFINITIONS['side']['drawing_name']}.svg"
    unchanged_before = {"front": sha256(front_svg), "side": sha256(side_svg)}

    os.chdir(ROOT)
    with contextlib.suppress(Exception):
        addon_utils.disable(
            "bl_ext.user_default.project_control", default_set=False, handle_error=None
        )
    load_result = bpy.ops.bim.load_project(
        filepath=str(session_ifc),
        should_start_fresh_session=True,
        use_detailed_tooltip=True,
    )
    if load_result != {"FINISHED"} or not tool.Ifc.get():
        raise RuntimeError(f"Bonsai failed to load 505 UP drawing session: {load_result}")
    for handler in list(bpy.app.handlers.depsgraph_update_post):
        if getattr(handler, "__module__", "").startswith("project_control"):
            bpy.app.handlers.depsgraph_update_post.remove(handler)

    model = tool.Ifc.get()
    target = model.by_guid(TARGET_GLOBAL_ID)
    drawing = next(
        item
        for item in model.by_type("IfcAnnotation")
        if item.ObjectType == "DRAWING"
        and item.Name == VIEW_DEFINITIONS["plan"]["drawing_name"]
    )
    annotation = next(
        item
        for item in model.by_type("IfcAnnotation")
        if item.ObjectType == "LINEWORK"
        and item.Name == "Molteni 505 UP approved black proxy / plan"
    )
    camera = tool.Ifc.get_object(drawing) or tool.Drawing.import_drawing(drawing)
    drawing_group = next(
        inverse.RelatingGroup
        for inverse in model.get_inverse(drawing)
        if inverse.is_a("IfcRelAssignsToGroup")
    )
    tool.Drawing.import_annotations_in_group(drawing_group)
    annotation_obj = tool.Ifc.get_object(annotation)
    if target is None or camera is None or annotation_obj is None:
        raise RuntimeError("505 UP Plan Drawing objects did not reload")

    representation = next(
        item
        for item in annotation.Representation.Representations
        if item.RepresentationIdentifier == "Annotation"
        and item.RepresentationType == "GeometricCurveSet"
    )
    curve_set = next(
        item for item in representation.Items if item.is_a("IfcGeometricCurveSet")
    )
    polylines = []
    for path in plan_paths:
        points = [
            model.create_entity(
                "IfcCartesianPoint", Coordinates=coordinates_mm("plan", first, second)
            )
            for first, second in path
        ]
        polylines.append(model.create_entity("IfcPolyline", Points=points))
    curve_set.Elements = polylines

    vertices = []
    edges = []
    for path in plan_paths:
        indices = []
        for first, second in path:
            indices.append(len(vertices))
            vertices.append(
                tuple(
                    value / 1000.0
                    for value in coordinates_mm("plan", first, second)
                )
            )
        edges.extend(zip(indices, indices[1:]))
    annotation_obj.data.clear_geometry()
    annotation_obj.data.from_pydata(vertices, edges, [])
    annotation_obj.data.update()

    for pset in model.by_type("IfcPropertySet"):
        if pset.Name != "Pset_Molteni505UpDrawingSource":
            continue
        properties = {item.Name: item for item in pset.HasProperties}
        if "PlanPathCount" in properties:
            properties["PlanPathCount"].NominalValue = model.create_entity(
                "IfcText", str(EXPECTED_PATH_COUNTS["plan"])
            )

    override = view3d_override()
    with bpy.context.temp_override(**override):
        activate_result = bpy.ops.bim.activate_drawing(
            drawing=drawing.id(), should_view_from_camera=False
        )
    if activate_result != {"FINISHED"}:
        raise RuntimeError(f"failed to activate repaired Plan Drawing: {activate_result}")
    dprops = tool.Drawing.get_document_props()
    dprops.should_use_underlay_cache = False
    dprops.should_use_linework_cache = False
    dprops.should_use_annotation_cache = False
    with bpy.context.temp_override(**override):
        create_result = bpy.ops.bim.create_drawing(
            print_all=False, open_viewer=False, sync=False
        )
    output_svg = output_dir / f"{VIEW_DEFINITIONS['plan']['drawing_name']}.svg"
    if (
        create_result != {"FINISHED"}
        or not output_svg.is_file()
        or output_svg.stat().st_size == 0
    ):
        raise RuntimeError(f"Bonsai Plan Create Drawing failed: {create_result}")
    style = style_svg(output_svg, annotation.GlobalId, "plan")
    svg = shared.inspect_svg(output_svg, TARGET_GLOBAL_ID, annotation.GlobalId)
    if (
        svg["root_tag"] != "svg"
        or svg["geometry_element_count"] == 0
        or svg["projection_group_count"] == 0
        or svg["target_ifc_projection_group_count"] != 0
        or svg["review_annotation_geometry_count"] == 0
    ):
        raise RuntimeError(f"repaired 505 UP Plan SVG gate failed: {svg}")

    temporary = session_ifc.with_suffix(".ifc.next")
    model.write(str(temporary))
    os.replace(temporary, session_ifc)
    session_after_sha = sha256(session_ifc)
    reload_result = bpy.ops.bim.load_project(
        filepath=str(session_ifc),
        should_start_fresh_session=True,
        use_detailed_tooltip=True,
    )
    if reload_result != {"FINISHED"} or not tool.Ifc.get():
        raise RuntimeError(f"Bonsai Plan repair reload failed: {reload_result}")
    reopened = tool.Ifc.get()
    persisted_counts = {}
    drawing_count = 0
    annotation_count = 0
    for view, definition in VIEW_DEFINITIONS.items():
        reopened_drawing = next(
            item
            for item in reopened.by_type("IfcAnnotation")
            if item.ObjectType == "DRAWING" and item.Name == definition["drawing_name"]
        )
        reopened_annotation = next(
            item
            for item in reopened.by_type("IfcAnnotation")
            if item.ObjectType == "LINEWORK"
            and item.Name == f"Molteni 505 UP approved black proxy / {view}"
        )
        drawing_count += 1
        annotation_count += 1
        reopened_representation = next(
            item
            for item in reopened_annotation.Representation.Representations
            if item.RepresentationIdentifier == "Annotation"
            and item.RepresentationType == "GeometricCurveSet"
        )
        paths = [path for item in reopened_representation.Items for path in item.Elements]
        persisted_counts[view] = len(paths)
        include = ifcopenshell.util.element.get_pset(
            reopened_drawing, "EPset_Drawing"
        ).get("Include", "").split(",")
        if TARGET_GLOBAL_ID in include:
            raise RuntimeError(f"505 UP target Body re-entered {view} projection")
    if persisted_counts != EXPECTED_PATH_COUNTS:
        raise RuntimeError(f"505 UP post-reload path counts drifted: {persisted_counts}")

    reopened_annotation = reopened.by_guid(annotation.GlobalId)
    reopened_representation = next(
        item
        for item in reopened_annotation.Representation.Representations
        if item.RepresentationIdentifier == "Annotation"
    )
    persisted_paths = [
        path for item in reopened_representation.Items for path in item.Elements
    ]
    points = [
        tuple(float(value) for value in point.Coordinates)
        for path in persisted_paths
        for point in path.Points
    ]
    expected_points = [
        coordinates_mm("plan", first, second)
        for path in plan_paths
        for first, second in path
    ]
    persisted_bounds = bounds_3d(points)
    expected_bounds = bounds_3d(expected_points)
    residual = max(
        abs(persisted_bounds[bound][axis] - expected_bounds[bound][axis])
        for bound in range(2)
        for axis in range(3)
    )
    if residual > 0.000001:
        raise RuntimeError(f"505 UP Plan persisted coordinate residual: {residual}")

    unchanged_after = {"front": sha256(front_svg), "side": sha256(side_svg)}
    if unchanged_after != unchanged_before:
        raise RuntimeError("Front/Side SVG changed during Plan-only Create Drawing")
    blend_path = PRODUCT_DIR / "Molteni-505-UP-V1-LP-S-project-drawings.blend"
    bpy.ops.wm.save_as_mainfile(filepath=str(blend_path))
    if sha256(formal_ifc) != FORMAL_SHA256 or sha256(FORMAL_IFC) != FORMAL_SHA256:
        raise RuntimeError("formal IFC changed during 505 UP Plan repair")

    evidence_path = PRODUCT_DIR / "505-UP-ENTRANCE-create-drawing-evidence.json"
    evidence = json.loads(evidence_path.read_text(encoding="utf-8"))
    plan_record = next(item for item in evidence["views"] if item["view"] == "plan")
    cache_path = output_svg.parent / "cache" / f"{output_svg.stem}-linework.svg"
    drawing_pset = ifcopenshell.util.element.get_pset(
        reopened.by_guid(drawing.GlobalId), "EPset_Drawing"
    )
    annotation_pset = ifcopenshell.util.element.get_pset(
        reopened_annotation, "EPset_Annotation"
    )
    plan_record.update(
        {
            "review_path_count": EXPECTED_PATH_COUNTS["plan"],
            "persisted_review_path_count": len(persisted_paths),
            "review_edge_count": sum(max(0, len(path) - 1) for path in plan_paths),
            "persisted_coordinate_bounds_mm": persisted_bounds,
            "expected_coordinate_bounds_mm": expected_bounds,
            "persisted_coordinate_residual_mm": residual,
            "component_boundary_extraction": repair_audit,
            "drawing": {
                "global_id": drawing.GlobalId,
                "name": drawing.Name,
                "object_type": drawing.ObjectType,
                "epset_drawing": drawing_pset,
            },
            "linework_annotation": {
                **plan_record["linework_annotation"],
                "epset_annotation": annotation_pset,
            },
            "create_drawing": {
                "operator": "bpy.ops.bim.create_drawing",
                "arguments": {"print_all": False, "open_viewer": False, "sync": False},
                "result": sorted(create_result),
                "linework_mode": "OPENCASCADE",
                "target_view": VIEW_DEFINITIONS["plan"]["target_view"],
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
        }
    )
    evidence.update(
        {
            "generated_at": datetime.now(timezone.utc).isoformat(),
            "task": "Plan-only Bonsai Create Drawing refresh for the mechanically restored 505 UP DISPLAY component boundary",
            "drawing_session_sha256_before": session_before_sha,
            "drawing_session_sha256_after": session_after_sha,
            "approval_record_sha256": sha256(APPROVAL),
            "candidate_sha256": sha256(CANDIDATE),
            "bonsai_session": {"path": str(blend_path), "sha256": sha256(blend_path)},
            "plan_component_boundary_repair": repair_audit,
            "plan_only_refresh": {
                "create_result": sorted(create_result),
                "front_side_svg_sha256_before": unchanged_before,
                "front_side_svg_sha256_after": unchanged_after,
                "front_side_unchanged": unchanged_before == unchanged_after,
            },
            "persistence": {
                "method": "update existing Plan LINEWORK; model.write to temporary then atomic replace; bpy.ops.bim.load_project reload",
                "reload_result": sorted(reload_result),
                "post_reload_drawing_count": drawing_count,
                "post_reload_annotation_count": annotation_count,
                "post_reload_path_counts": persisted_counts,
            },
            "pass": True,
        }
    )
    evidence_path.write_text(
        json.dumps(evidence, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print(
        json.dumps(
            {
                "drawing_session_ifc": str(session_ifc),
                "evidence": str(evidence_path),
                "blend": str(blend_path),
                "create_drawing_result": sorted(create_result),
                "persisted_path_counts": persisted_counts,
                "front_side_unchanged": unchanged_before == unchanged_after,
                "target_body_projection_count": svg["target_ifc_projection_group_count"],
                "pass": True,
            },
            ensure_ascii=False,
        )
    )


def refresh_plan_front_revision(session_ifc: Path, formal_ifc: Path, output_dir: Path):
    """Persist the open Plan seam and correct Front to the north camera."""
    if sha256(formal_ifc) != FORMAL_SHA256 or sha256(FORMAL_IFC) != FORMAL_SHA256:
        raise RuntimeError("formal IFC hash mismatch")
    session_before_sha = sha256(session_ifc)
    approval = json.loads(APPROVAL.read_text(encoding="utf-8"))
    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    candidate = json.loads(CANDIDATE.read_text(encoding="utf-8"))
    if (
        approval.get("status") != "revision_pending_review"
        or approval.get("pending_reapproval_views") != ["plan", "front"]
        or approval.get("derived_ifc_write_allowed") is not True
        or approval.get("formal_authoritative_ifc_write_allowed") is not False
        or approval.get("candidate_manifest_sha256") != sha256(MANIFEST)
        or candidate.get("review_status")
        not in ("visual_review_pending", "revision_pending_review")
        or candidate.get("source_kind") != SOURCE_KIND
        or candidate.get("official_cad_used") is not False
        or manifest.get("formal_ifc_sha256") != FORMAL_SHA256
    ):
        raise RuntimeError("505 UP revision/source gate failed")
    approved_paths = {
        view: candidate["views"][view]["proxy_paths_mm"] for view in VIEW_DEFINITIONS
    }
    if {view: len(paths) for view, paths in approved_paths.items()} != EXPECTED_PATH_COUNTS:
        raise RuntimeError("505 UP revision path counts drifted")
    component_audit = candidate["views"]["plan"].get("component_semantics")
    if (
        not component_audit
        or component_audit.get("interface", {}).get("closed") is not False
        or component_audit.get("interface", {}).get("full_display_footprint_added")
        is not False
    ):
        raise RuntimeError("505 UP Plan open-interface audit missing")

    output_dir.mkdir(parents=True, exist_ok=True)
    output_svgs = {
        view: output_dir / f"{definition['drawing_name']}.svg"
        for view, definition in VIEW_DEFINITIONS.items()
    }
    before_hashes = {view: sha256(path) for view, path in output_svgs.items()}

    os.chdir(ROOT)
    with contextlib.suppress(Exception):
        addon_utils.disable(
            "bl_ext.user_default.project_control", default_set=False, handle_error=None
        )
    load_result = bpy.ops.bim.load_project(
        filepath=str(session_ifc),
        should_start_fresh_session=True,
        use_detailed_tooltip=True,
    )
    if load_result != {"FINISHED"} or not tool.Ifc.get():
        raise RuntimeError(f"Bonsai failed to load 505 UP drawing session: {load_result}")
    for handler in list(bpy.app.handlers.depsgraph_update_post):
        if getattr(handler, "__module__", "").startswith("project_control"):
            bpy.app.handlers.depsgraph_update_post.remove(handler)

    model = tool.Ifc.get()
    target = model.by_guid(TARGET_GLOBAL_ID)
    room = model.by_guid(ROOM_GLOBAL_ID)
    target_obj = tool.Ifc.get_object(target)
    room_obj = tool.Ifc.get_object(room)
    if target is None or room is None or target_obj is None or room_obj is None:
        raise RuntimeError("505 UP target or entrance room did not reload")
    crop_bbox = union_bbox(shared.world_bbox(room_obj), shared.world_bbox(target_obj))

    states = {}
    for view in ("plan", "front"):
        definition = VIEW_DEFINITIONS[view]
        drawing = next(
            item
            for item in model.by_type("IfcAnnotation")
            if item.ObjectType == "DRAWING" and item.Name == definition["drawing_name"]
        )
        annotation = next(
            item
            for item in model.by_type("IfcAnnotation")
            if item.ObjectType == "LINEWORK"
            and item.Name == f"Molteni 505 UP approved black proxy / {view}"
        )
        camera = tool.Ifc.get_object(drawing) or tool.Drawing.import_drawing(drawing)
        drawing_group = next(
            inverse.RelatingGroup
            for inverse in model.get_inverse(drawing)
            if inverse.is_a("IfcRelAssignsToGroup")
        )
        tool.Drawing.import_annotations_in_group(drawing_group)
        annotation_obj = tool.Ifc.get_object(annotation)
        if camera is None or annotation_obj is None:
            raise RuntimeError(f"505 UP {view} Drawing objects did not import")
        states[view] = {
            "drawing": drawing,
            "annotation": annotation,
            "camera": camera,
            "annotation_obj": annotation_obj,
        }

    plan_annotation = states["plan"]["annotation"]
    plan_representation = next(
        item
        for item in plan_annotation.Representation.Representations
        if item.RepresentationIdentifier == "Annotation"
        and item.RepresentationType == "GeometricCurveSet"
    )
    plan_curve_set = next(
        item for item in plan_representation.Items if item.is_a("IfcGeometricCurveSet")
    )
    plan_polylines = []
    for path in approved_paths["plan"]:
        points = [
            model.create_entity(
                "IfcCartesianPoint", Coordinates=coordinates_mm("plan", first, second)
            )
            for first, second in path
        ]
        plan_polylines.append(model.create_entity("IfcPolyline", Points=points))
    plan_curve_set.Elements = plan_polylines

    vertices = []
    edges = []
    for path in approved_paths["plan"]:
        indices = []
        for first, second in path:
            indices.append(len(vertices))
            vertices.append(
                tuple(
                    value / 1000.0
                    for value in coordinates_mm("plan", first, second)
                )
            )
        edges.extend(zip(indices, indices[1:]))
    plan_obj = states["plan"]["annotation_obj"]
    plan_obj.data.clear_geometry()
    plan_obj.data.from_pydata(vertices, edges, [])
    plan_obj.data.update()

    front_orientation = orient_front_camera_north(
        model,
        states["front"]["drawing"],
        states["front"]["camera"],
        target_obj,
        crop_bbox,
    )

    override = view3d_override()
    generated = {}
    for view in ("plan", "front"):
        drawing = states[view]["drawing"]
        annotation = states[view]["annotation"]
        with bpy.context.temp_override(**override):
            activate_result = bpy.ops.bim.activate_drawing(
                drawing=drawing.id(), should_view_from_camera=False
            )
        if activate_result != {"FINISHED"}:
            raise RuntimeError(f"failed to activate revised {view} Drawing")
        dprops = tool.Drawing.get_document_props()
        dprops.should_use_underlay_cache = False
        dprops.should_use_linework_cache = False
        dprops.should_use_annotation_cache = False
        with bpy.context.temp_override(**override):
            create_result = bpy.ops.bim.create_drawing(
                print_all=False, open_viewer=False, sync=False
            )
        output_svg = output_svgs[view]
        if (
            create_result != {"FINISHED"}
            or not output_svg.is_file()
            or output_svg.stat().st_size == 0
        ):
            raise RuntimeError(f"Bonsai Create Drawing failed for revised {view}")
        style = style_svg(output_svg, annotation.GlobalId, view)
        svg = shared.inspect_svg(output_svg, TARGET_GLOBAL_ID, annotation.GlobalId)
        if (
            svg["root_tag"] != "svg"
            or svg["geometry_element_count"] == 0
            or svg["target_ifc_projection_group_count"] != 0
            or svg["review_annotation_geometry_count"] == 0
        ):
            raise RuntimeError(f"revised 505 UP {view} SVG gate failed: {svg}")
        cache_path = output_svg.parent / "cache" / f"{output_svg.stem}-linework.svg"
        generated[view] = {
            "create_result": sorted(create_result),
            "style": style,
            "svg": svg,
            "cache_path": cache_path,
        }

    temporary = session_ifc.with_suffix(".ifc.next")
    model.write(str(temporary))
    os.replace(temporary, session_ifc)
    session_after_sha = sha256(session_ifc)
    reload_result = bpy.ops.bim.load_project(
        filepath=str(session_ifc),
        should_start_fresh_session=True,
        use_detailed_tooltip=True,
    )
    if reload_result != {"FINISHED"} or not tool.Ifc.get():
        raise RuntimeError(f"Bonsai 505 revision reload failed: {reload_result}")
    reopened = tool.Ifc.get()
    persisted_counts = {}
    for view, definition in VIEW_DEFINITIONS.items():
        drawing = next(
            item
            for item in reopened.by_type("IfcAnnotation")
            if item.ObjectType == "DRAWING" and item.Name == definition["drawing_name"]
        )
        annotation = next(
            item
            for item in reopened.by_type("IfcAnnotation")
            if item.ObjectType == "LINEWORK"
            and item.Name == f"Molteni 505 UP approved black proxy / {view}"
        )
        persisted_counts[view] = sum(
            len(item.Elements or [])
            for representation in annotation.Representation.Representations
            for item in representation.Items
            if item.is_a("IfcGeometricCurveSet")
        )
        include = ifcopenshell.util.element.get_pset(drawing, "EPset_Drawing").get(
            "Include", ""
        ).split(",")
        if TARGET_GLOBAL_ID in include:
            raise RuntimeError(f"505 UP target Body re-entered {view} projection")
    if persisted_counts != EXPECTED_PATH_COUNTS:
        raise RuntimeError(f"505 UP post-reload path counts drifted: {persisted_counts}")

    reopened_front = next(
        item
        for item in reopened.by_type("IfcAnnotation")
        if item.ObjectType == "DRAWING"
        and item.Name == VIEW_DEFINITIONS["front"]["drawing_name"]
    )
    front_matrix = ifcopenshell.util.placement.get_local_placement(
        reopened_front.ObjectPlacement
    )
    front_rotation = [[round(float(front_matrix[row][column]), 9) for column in range(3)] for row in range(3)]
    expected_rotation = [[-1.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, 1.0, 0.0]]
    rotation_residual = max(
        abs(front_rotation[row][column] - expected_rotation[row][column])
        for row in range(3)
        for column in range(3)
    )
    if rotation_residual > 0.000001:
        raise RuntimeError(f"505 Front camera did not persist NORTH: {front_rotation}")
    front_orientation["persisted_camera_matrix_project_units"] = matrix_record(front_matrix)
    front_orientation["persisted_rotation_residual"] = rotation_residual
    front_orientation["persisted_north_rotation_pass"] = True

    after_hashes = {view: sha256(path) for view, path in output_svgs.items()}
    if after_hashes["side"] != before_hashes["side"]:
        raise RuntimeError("505 Side SVG changed during Plan/Front revision")
    if (
        front_orientation["input_was_mirrored_south_camera"]
        and after_hashes["front"] == before_hashes["front"]
    ):
        raise RuntimeError("505 Front SVG did not change after handedness correction")
    blend_path = PRODUCT_DIR / "Molteni-505-UP-V1-LP-S-project-drawings.blend"
    bpy.ops.wm.save_as_mainfile(filepath=str(blend_path))
    if sha256(formal_ifc) != FORMAL_SHA256 or sha256(FORMAL_IFC) != FORMAL_SHA256:
        raise RuntimeError("formal IFC changed during 505 revision")

    evidence_path = PRODUCT_DIR / "505-UP-ENTRANCE-create-drawing-evidence.json"
    evidence = json.loads(evidence_path.read_text(encoding="utf-8"))
    for view in ("plan", "front"):
        record = next(item for item in evidence["views"] if item["view"] == view)
        annotation = states[view]["annotation"]
        output_svg = output_svgs[view]
        cache_path = generated[view]["cache_path"]
        record.update(
            {
                "review_path_count": EXPECTED_PATH_COUNTS[view],
                "persisted_review_path_count": EXPECTED_PATH_COUNTS[view],
                "create_drawing": {
                    "operator": "bpy.ops.bim.create_drawing",
                    "arguments": {
                        "print_all": False,
                        "open_viewer": False,
                        "sync": False,
                    },
                    "result": generated[view]["create_result"],
                    "linework_mode": "OPENCASCADE",
                    "target_view": VIEW_DEFINITIONS[view]["target_view"],
                },
                "svg": {
                    "path": str(output_svg),
                    "bytes": output_svg.stat().st_size,
                    "sha256": sha256(output_svg),
                    **generated[view]["style"],
                    **generated[view]["svg"],
                },
                "linework_cache": {
                    "path": str(cache_path),
                    "bytes": cache_path.stat().st_size,
                    "sha256": sha256(cache_path),
                },
            }
        )
        if view == "plan":
            record.pop("component_boundary_extraction", None)
            record["component_semantics"] = component_audit

    evidence.update(
        {
            "generated_at": datetime.now(timezone.utc).isoformat(),
            "task": "Bonsai Plan open-interface and Front north-camera revision for 505 UP",
            "review_status": "revision_pending_review",
            "drawing_session_sha256_before": session_before_sha,
            "drawing_session_sha256_after": session_after_sha,
            "approval_record_sha256": sha256(APPROVAL),
            "candidate_sha256": sha256(CANDIDATE),
            "bonsai_session": {"path": str(blend_path), "sha256": sha256(blend_path)},
            "plan_component_semantics": component_audit,
            "front_orientation_revision": front_orientation,
            "plan_front_revision": {
                "create_results": {
                    view: generated[view]["create_result"] for view in generated
                },
                "svg_sha256_before": before_hashes,
                "svg_sha256_after": after_hashes,
                "side_unchanged": before_hashes["side"] == after_hashes["side"],
                "target_body_projection_counts": {
                    view: generated[view]["svg"]["target_ifc_projection_group_count"]
                    for view in generated
                },
            },
            "persistence": {
                "method": "update Plan LINEWORK and Front Drawing placement; atomic IFC replace; Bonsai reload",
                "reload_result": sorted(reload_result),
                "post_reload_drawing_count": 3,
                "post_reload_annotation_count": 3,
                "post_reload_path_counts": persisted_counts,
            },
            "pass": True,
        }
    )
    evidence.pop("plan_component_boundary_repair", None)
    evidence.pop("plan_only_refresh", None)
    evidence_path.write_text(
        json.dumps(evidence, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print(
        json.dumps(
            {
                "drawing_session_ifc": str(session_ifc),
                "evidence": str(evidence_path),
                "create_results": {
                    view: generated[view]["create_result"] for view in generated
                },
                "persisted_path_counts": persisted_counts,
                "front_orientation": front_orientation,
                "side_unchanged": before_hashes["side"] == after_hashes["side"],
                "formal_ifc_sha256": sha256(formal_ifc),
                "pass": True,
            },
            ensure_ascii=False,
        )
    )


def main():
    arguments = sys.argv[sys.argv.index("--") + 1 :]
    if len(arguments) == 4 and arguments[3] == "--plan-front-revision":
        return refresh_plan_front_revision(
            Path(arguments[0]).resolve(),
            Path(arguments[1]).resolve(),
            Path(arguments[2]).resolve(),
        )
    if len(arguments) == 4 and arguments[3] == "--plan-only":
        return refresh_plan_only(
            Path(arguments[0]).resolve(),
            Path(arguments[1]).resolve(),
            Path(arguments[2]).resolve(),
        )
    if len(arguments) != 3:
        raise SystemExit(
            "expected: drawing-session.ifc formal.ifc output-directory [--plan-only]"
        )
    session_ifc = Path(arguments[0]).resolve()
    formal_ifc = Path(arguments[1]).resolve()
    output_dir = Path(arguments[2]).resolve()
    if sha256(formal_ifc) != FORMAL_SHA256 or sha256(FORMAL_IFC) != FORMAL_SHA256:
        raise RuntimeError("formal IFC hash mismatch")
    session_before_sha = sha256(session_ifc)

    approval = json.loads(APPROVAL.read_text(encoding="utf-8"))
    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    candidate = json.loads(CANDIDATE.read_text(encoding="utf-8"))
    if (
        approval.get("status") != "approved"
        or approval.get("derived_ifc_write_allowed") is not True
        or approval.get("formal_authoritative_ifc_write_allowed") is not False
        or approval.get("approved_views") != ["plan", "front", "side"]
        or approval.get("candidate_manifest_sha256") != sha256(MANIFEST)
        or approval.get("approved_candidate", {}).get("source_kind") != SOURCE_KIND
        or approval.get("approved_candidate", {}).get("official_cad_used") is not False
        or manifest.get("formal_ifc_sha256") != FORMAL_SHA256
        or candidate.get("source_kind") != SOURCE_KIND
        or candidate.get("official_cad_used") is not False
    ):
        raise RuntimeError("505 UP approval/source gate failed")
    approved_paths = {}
    for view in VIEW_DEFINITIONS:
        item = candidate["views"][view]
        paths = item["proxy_paths_mm"]
        if (
            item.get("source_kind") != SOURCE_KIND
            or item.get("official_cad_paths_mm") != []
            or len(paths) != EXPECTED_PATH_COUNTS[view]
        ):
            raise RuntimeError(f"505 UP {view} approved candidate gate failed")
        approved_paths[view] = paths

    os.chdir(ROOT)
    with contextlib.suppress(Exception):
        addon_utils.disable(
            "bl_ext.user_default.project_control", default_set=False, handle_error=None
        )
    load_result = bpy.ops.bim.load_project(
        filepath=str(session_ifc),
        should_start_fresh_session=True,
        use_detailed_tooltip=True,
    )
    if load_result != {"FINISHED"} or not tool.Ifc.get():
        raise RuntimeError(f"Bonsai failed to load 505 UP drawing session: {load_result}")
    for handler in list(bpy.app.handlers.depsgraph_update_post):
        if getattr(handler, "__module__", "").startswith("project_control"):
            bpy.app.handlers.depsgraph_update_post.remove(handler)

    model = tool.Ifc.get()
    target = model.by_guid(TARGET_GLOBAL_ID)
    room = model.by_guid(ROOM_GLOBAL_ID)
    target_obj = tool.Ifc.get_object(target)
    room_obj = tool.Ifc.get_object(room)
    if target is None or room is None or target_obj is None or room_obj is None:
        raise RuntimeError("505 UP target or entrance space missing")
    if room.LongName != ROOM_NAME:
        raise RuntimeError("entrance-space identity drifted")

    room_bbox = shared.world_bbox(room_obj)
    target_bbox = shared.world_bbox(target_obj)
    crop_bbox = union_bbox(room_bbox, target_bbox)
    records = shared.room_elements(crop_bbox)
    original_context = [record[0] for record in records]
    if original_context.count(target) != 1:
        raise RuntimeError("expected 505 UP exactly once in entrance context")
    context_elements = [element for element in original_context if element != target]
    context_counts = Counter(element.is_a() for element in context_elements)
    output_dir.mkdir(parents=True, exist_ok=True)

    view_records = []
    override = view3d_override()
    for view, definition in VIEW_DEFINITIONS.items():
        output_svg = output_dir / f"{definition['drawing_name']}.svg"
        drawing, camera, width, height, clip_end = add_505_drawing(
            model, definition, crop_bbox, context_elements, output_svg
        )
        with bpy.context.temp_override(**override):
            activate_result = bpy.ops.bim.activate_drawing(
                drawing=drawing.id(), should_view_from_camera=False
            )
        if activate_result != {"FINISHED"}:
            raise RuntimeError(f"failed to activate {view} Drawing: {activate_result}")
        annotation, path_count, edge_count = add_proxy_annotation(
            model, drawing, target, target_obj, view, approved_paths[view]
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
        with bpy.context.temp_override(**override):
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
        svg = shared.inspect_svg(output_svg, TARGET_GLOBAL_ID, annotation.GlobalId)
        if (
            svg["root_tag"] != "svg"
            or svg["geometry_element_count"] == 0
            or svg["projection_group_count"] == 0
            or svg["target_ifc_projection_group_count"] != 0
            or svg["review_annotation_geometry_count"] == 0
        ):
            raise RuntimeError(f"505 UP {view} SVG structure/layer gate failed: {svg}")
        cache_path = output_svg.parent / "cache" / f"{output_svg.stem}-linework.svg"
        if not cache_path.is_file() or cache_path.stat().st_size == 0:
            raise RuntimeError(f"505 UP {view} linework cache missing")
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
                    for path in approved_paths[view]
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

    temporary = session_ifc.with_suffix(".ifc.next")
    model.write(str(temporary))
    os.replace(temporary, session_ifc)
    session_after_sha = sha256(session_ifc)
    reload_result = bpy.ops.bim.load_project(
        filepath=str(session_ifc),
        should_start_fresh_session=True,
        use_detailed_tooltip=True,
    )
    if reload_result != {"FINISHED"} or not tool.Ifc.get():
        raise RuntimeError(f"Bonsai reload failed: {reload_result}")
    reopened = tool.Ifc.get()

    outputs = []
    for record in view_records:
        view = record["view"]
        drawing = reopened.by_guid(record["drawing_global_id"])
        annotation = reopened.by_guid(record["annotation_global_id"])
        if drawing is None or annotation is None:
            raise RuntimeError(f"reloaded IFC lost {view} Drawing or LINEWORK Annotation")
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
        persisted_bounds = bounds_3d(points)
        expected_bounds = bounds_3d(record["expected_points"])
        residual = max(
            abs(persisted_bounds[bound][axis] - expected_bounds[bound][axis])
            for bound in range(2)
            for axis in range(3)
        )
        annotation_pset = ifcopenshell.util.element.get_pset(annotation, "EPset_Annotation")
        drawing_pset = ifcopenshell.util.element.get_pset(drawing, "EPset_Drawing")
        include = drawing_pset.get("Include", "").split(",")
        assignment = next(
            (
                relation
                for relation in reopened.by_type("IfcRelAssignsToGroup")
                if drawing in relation.RelatedObjects and annotation in relation.RelatedObjects
            ),
            None,
        )
        styled_items = [styled for item in representation.Items for styled in item.StyledByItem]
        if (
            len(paths) != EXPECTED_PATH_COUNTS[view]
            or residual > 0.000001
            or annotation_pset.get("SourceKind") != SOURCE_KIND
            or annotation_pset.get("OfficialCadUsed") is not False
            or not styled_items
            or assignment is None
            or TARGET_GLOBAL_ID in include
            or len(include) != len(context_elements)
        ):
            raise RuntimeError(f"persisted 505 UP {view} gate failed")
        outputs.append(
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
                    "ifc_curve_style_colour": BLACK,
                    "drawing_group_assignment_global_id": assignment.GlobalId,
                },
                "review_path_count": EXPECTED_PATH_COUNTS[view],
                "persisted_review_path_count": len(paths),
                "review_edge_count": record["edge_count"],
                "persisted_coordinate_bounds_mm": persisted_bounds,
                "expected_coordinate_bounds_mm": expected_bounds,
                "persisted_coordinate_residual_mm": residual,
                "camera": record["camera"],
                "create_drawing": {
                    "operator": "bpy.ops.bim.create_drawing",
                    "arguments": {"print_all": False, "open_viewer": False, "sync": False},
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

    blend_path = PRODUCT_DIR / "Molteni-505-UP-V1-LP-S-project-drawings.blend"
    bpy.ops.wm.save_as_mainfile(filepath=str(blend_path))
    if sha256(formal_ifc) != FORMAL_SHA256 or sha256(FORMAL_IFC) != FORMAL_SHA256:
        raise RuntimeError("formal IFC changed during 505 UP Drawing creation")
    evidence = {
        "schema_version": 1,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "task": "Bonsai Create Drawing for approved black Molteni 505 UP proxy in the actual entrance-space context",
        "workflow": ["inspect", "plan", "execute", "persist", "reload", "verify"],
        "provider": {
            "name": "bonsai-mcp",
            "status": "supported",
            "bridge": "127.0.0.1:9878",
            "execution": "execute_blender_code",
        },
        "course_evidence": {
            "lesson": "085000 Introduction to Drawings",
            "source": "/Users/jiaxinchen/.codex/worktrees/fb76/2504 GBTB Yanlord Zhuhai/research/bonsai-course/lessons/085000/README.md",
            "screenshots": {
                "active_drawing_layers": "/Users/jiaxinchen/.codex/worktrees/fb76/2504 GBTB Yanlord Zhuhai/research/bonsai-course/lessons/085000/screenshots/085000-01m03s-active-drawing-camera.png",
                "create_drawing": "/Users/jiaxinchen/.codex/worktrees/fb76/2504 GBTB Yanlord Zhuhai/research/bonsai-course/lessons/085000/screenshots/085000-01m59s-create-drawing-button.png",
                "svg_result": "/Users/jiaxinchen/.codex/worktrees/fb76/2504 GBTB Yanlord Zhuhai/research/bonsai-course/lessons/085000/screenshots/085000-02m13s-svg-in-browser.png",
            },
            "course_fact": "Create Drawing generates or refreshes SVG after Drawing camera, scale, depth, layers and filters are configured.",
            "current_version_inference": "The recorded 4.2.1/0.8.2 workflow is adapted to Blender 4.5.3 LTS and Bonsai/IfcOpenShell 0.8.4 using semantic Drawing state and bpy.ops.bim.create_drawing.",
        },
        "formal_ifc": str(formal_ifc),
        "formal_ifc_sha256": sha256(formal_ifc),
        "formal_ifc_bytes_unchanged": sha256(formal_ifc) == FORMAL_SHA256,
        "drawing_session_ifc": str(session_ifc),
        "drawing_session_sha256_before": session_before_sha,
        "drawing_session_sha256_after": session_after_sha,
        "save_boundary": "one 505 UP product-level full-project IFC copy under its review package",
        "bonsai_session": {"path": str(blend_path), "sha256": sha256(blend_path)},
        "approval_record": str(APPROVAL),
        "approval_record_sha256": sha256(APPROVAL),
        "candidate": str(CANDIDATE),
        "candidate_sha256": sha256(CANDIDATE),
        "source_kind": SOURCE_KIND,
        "source_label_zh": SOURCE_LABEL_ZH,
        "official_cad_used": False,
        "blue_official_atomic_component_composition_selected": False,
        "target": {
            "global_id": TARGET_GLOBAL_ID,
            "ifc_class": target.is_a(),
            "name": target.Name,
            "world_bbox_m": [list(target_bbox[0]), list(target_bbox[1])],
            "target_include_count_before_suppression": original_context.count(target),
            "target_include_count_after_suppression": 0,
        },
        "room": {
            "global_id": ROOM_GLOBAL_ID,
            "name": ROOM_NAME,
            "bbox_m": [list(room_bbox[0]), list(room_bbox[1])],
            "crop_bbox_m": [list(crop_bbox[0]), list(crop_bbox[1])],
        },
        "context": {
            "project_context_retained": True,
            "include_count": len(context_elements),
            "include_global_ids": [element.GlobalId for element in context_elements],
            "include_ifc_class_counts": dict(sorted(context_counts.items())),
            "context_colour": GREY,
            "target_colour": BLACK,
        },
        "views": outputs,
        "persistence": {
            "method": "model.write to temporary then atomic replace; bpy.ops.bim.load_project reload",
            "reload_result": sorted(reload_result),
            "post_reload_drawing_count": len(outputs),
            "post_reload_annotation_count": len(outputs),
        },
        "versions": {
            "blender": bpy.app.version_string,
            "ifcopenshell": ifcopenshell.version,
            "ifc_schema": reopened.schema,
            "bonsai_generator": "Bonsai 0.8.4 bpy.ops.bim.create_drawing",
        },
        "pass": True,
    }
    evidence_path = PRODUCT_DIR / "505-UP-ENTRANCE-create-drawing-evidence.json"
    evidence_path.write_text(
        json.dumps(evidence, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print(
        json.dumps(
            {
                "drawing_session_ifc": str(session_ifc),
                "evidence": str(evidence_path),
                "blend": str(blend_path),
                "create_drawing_results": {
                    item["view"]: item["create_drawing"]["result"] for item in outputs
                },
                "pass": True,
            },
            ensure_ascii=False,
        )
    )


if __name__ == "__main__":
    main()
