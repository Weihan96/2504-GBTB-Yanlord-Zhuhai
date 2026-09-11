#!/usr/bin/env python3
"""Create one native Bonsai Drawing for the approved sxb010 dining-bay blind.

Run inside Blender/Bonsai with a per-view copy of the approved product-derived
IFC.  The script keeps the dining-bay and dining-room project context, removes
the target Body from the ordinary projection to prevent a double expression,
and persists the approved merged-line proxy as a blue IFC LINEWORK Annotation.
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
PRODUCT_DIR = ROOT / "output/review/highpoly-types/sxb010"
FORMAL_IFC = ROOT / "2504 GBTB Yanlord Zhuhai.ifc"
FORMAL_SHA256 = "7a50b87e8f48a7c2bfbaa9f1dabdfca7155fa684b2325c66aa1ae15c8ab0a25c"
APPROVED_DERIVED = PRODUCT_DIR / "sxb010-derived-drawing.ifc"
CANDIDATE = PRODUCT_DIR / "candidate-representations.json"
SIMPLIFICATION_AUDIT = PRODUCT_DIR / "line-simplification-audit.json"
TARGET_GLOBAL_ID = "1O9JRXCI56VRUbpuLJy86Z"
DINING_BAY_GLOBAL_ID = "3AlgzRX3jCvPS420cah9Qn"
DINING_GLOBAL_ID = "2JnxeB$or6D8nffRg7GQbf"
ROOM_NAMES = {"餐厅飘窗", "餐厅"}
SOURCE_LABEL_ZH = "基于原始高模几何生成的简化图纸表达"
BLUE = "#1677c8"
CONTEXT_GREY = "#aab4be"
EXPECTED_PATH_COUNTS = {"plan": 5, "front": 51, "side": 55}
VIEW_DEFINITIONS = {
    "plan": {
        "drawing_name": "SXB010-DINING-BAY-PLAN",
        "target_view": "PLAN_VIEW",
        "location_hint": "PLAN",
    },
    "front": {
        "drawing_name": "SXB010-DINING-BAY-FRONT",
        "target_view": "ELEVATION_VIEW",
        "location_hint": "NORTH",
    },
    "side": {
        "drawing_name": "SXB010-DINING-BAY-SIDE",
        "target_view": "ELEVATION_VIEW",
        "location_hint": "EAST",
    },
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


def union_bbox(*bounds):
    return (
        tuple(min(bound[0][axis] for bound in bounds) for axis in range(3)),
        tuple(max(bound[1][axis] for bound in bounds) for axis in range(3)),
    )


def overlaps(first, second, padding=0.0):
    return all(
        first[1][axis] >= second[0][axis] - padding
        and first[0][axis] <= second[1][axis] + padding
        for axis in range(3)
    )


def context_elements(context_bbox):
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
            or not getattr(element, "GlobalId", None)
        ):
            continue
        bounds = world_bbox(obj)
        if overlaps(bounds, context_bbox, padding=0.35):
            records.append((element, obj, bounds))
    records.sort(key=lambda item: item[0].GlobalId)
    return records


def path_bounds(paths):
    first = [point[0] for path in paths for point in path]
    second = [point[1] for path in paths for point in path]
    return (min(first), min(second)), (max(first), max(second))


def source_coordinate_centres(candidate):
    front = candidate["views"]["front"]["proxy_paths_mm"]
    plan = candidate["views"]["plan"]["proxy_paths_mm"]
    front_bounds = path_bounds(front)
    plan_bounds = path_bounds(plan)
    return {
        "x": (front_bounds[0][0] + front_bounds[1][0]) / 2,
        "y": (front_bounds[0][1] + front_bounds[1][1]) / 2,
        "z": (plan_bounds[0][1] + plan_bounds[1][1]) / 2,
    }


def coordinates_mm(view, first, second, centres):
    if view == "plan":
        return (float(first), centres["y"], float(second))
    if view == "front":
        return (float(first), float(second), centres["z"])
    return (centres["x"], float(second), float(first))


def annotation_representation(model, context, view, paths, centres):
    polylines = []
    for path in paths:
        points = [
            model.create_entity(
                "IfcCartesianPoint",
                Coordinates=coordinates_mm(view, first, second, centres),
            )
            for first, second in path
        ]
        if len(points) >= 2:
            polylines.append(model.create_entity("IfcPolyline", Points=points))
    if len(polylines) != EXPECTED_PATH_COUNTS[view]:
        raise RuntimeError(f"{view} merged-line path count drifted")
    curve_set = model.create_entity("IfcGeometricCurveSet", Elements=polylines)
    return model.create_entity(
        "IfcShapeRepresentation",
        ContextOfItems=context,
        RepresentationIdentifier="Annotation",
        RepresentationType="GeometricCurveSet",
        Items=[curve_set],
    )


def add_linework_annotation(model, drawing, target, target_obj, view, paths, centres):
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
    annotation.Name = f"sxb010 approved merged-line blind / {view}"
    annotation.Description = (
        f"{SOURCE_LABEL_ZH}; 45 slat identities retained; close neighbouring lines "
        "merged without changing the outer envelope."
    )
    annotation.ObjectPlacement = target.ObjectPlacement
    obj.matrix_world = target_obj.matrix_world

    vertices = []
    edges = []
    for path in paths:
        indices = []
        for first, second in path:
            point = tuple(value / 1000 for value in coordinates_mm(view, first, second, centres))
            indices.append(len(vertices))
            vertices.append(point)
        edges.extend(zip(indices, indices[1:]))
    obj.data.clear_geometry()
    obj.data.from_pydata(vertices, edges, [])
    obj.data.update()

    representation = annotation_representation(model, context, view, paths, centres)
    annotation.Representation = model.create_entity(
        "IfcProductDefinitionShape", Representations=[representation]
    )
    pset = ifcopenshell.api.pset.add_pset(model, product=annotation, name="EPset_Annotation")
    ifcopenshell.api.pset.edit_pset(
        model,
        pset=pset,
        properties={
            "Classes": "geometry-derived-simplified-proxy review-target-sxb010",
            "TargetGlobalId": TARGET_GLOBAL_ID,
            "SourceKind": "geometry_derived_simplified_proxy",
            "SourceLabelZh": SOURCE_LABEL_ZH,
            "LineSimplification": "near-neighbour centreline merge; 45 slats retained",
        },
    )
    return annotation, representation, len(edges)


def camera_pose(camera, location, target):
    camera.location = Vector(location)
    camera.rotation_euler = (Vector(target) - camera.location).to_track_quat("-Z", "Y").to_euler()


def view3d_override():
    """Return a real GUI VIEW_3D override for operators called by the MCP bridge."""
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


def add_drawing(model, definition, crop_bbox, include_elements, output_svg):
    name = definition["drawing_name"]
    target_view = definition["target_view"]
    if any(
        item.ObjectType == "DRAWING" and item.Name == name
        for item in model.by_type("IfcAnnotation")
    ):
        raise RuntimeError(f"drawing already exists in fresh session: {name}")

    minimum, maximum = crop_bbox
    centre = tuple((minimum[axis] + maximum[axis]) / 2 for axis in range(3))
    if target_view == "PLAN_VIEW":
        storey = next(item for item in model.by_type("IfcBuildingStorey") if item.Name == "FFL")
        location_hint = storey.id()
        bpy.context.scene.cursor.location = (centre[0], centre[1], 3.25)
    else:
        location_hint = definition["location_hint"]
        bpy.context.scene.cursor.location = centre

    before = {
        item.id()
        for item in model.by_type("IfcAnnotation")
        if item.ObjectType == "DRAWING"
    }
    core_drawing.add_drawing(
        tool.Ifc,
        tool.Collector,
        tool.Drawing,
        target_view=target_view,
        location_hint=location_hint,
    )
    created = [
        item
        for item in model.by_type("IfcAnnotation")
        if item.ObjectType == "DRAWING" and item.id() not in before
    ]
    if len(created) != 1:
        raise RuntimeError(f"expected one Drawing, got {len(created)}")
    drawing = created[0]
    core_drawing.update_drawing_name(tool.Ifc, tool.Drawing, drawing=drawing, name=name)
    camera = tool.Ifc.get_object(drawing) or tool.Drawing.import_drawing(drawing)

    if target_view == "PLAN_VIEW":
        camera_pose(camera, (centre[0], centre[1], 4.20), (centre[0], centre[1], 0.0))
        width = maximum[0] - minimum[0] + 0.70
        height = maximum[1] - minimum[1] + 0.70
        clip_end = 4.50
    elif definition["location_hint"] == "NORTH":
        camera_pose(camera, (centre[0], maximum[1] + 1.20, 1.45), (centre[0], minimum[1], 1.45))
        width = maximum[0] - minimum[0] + 0.70
        height = 3.40
        clip_end = maximum[1] - minimum[1] + 2.00
    else:
        camera_pose(camera, (maximum[0] + 1.20, centre[1], 1.45), (minimum[0], centre[1], 1.45))
        width = maximum[1] - minimum[1] + 0.70
        height = 3.40
        clip_end = maximum[0] - minimum[0] + 2.00

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
    ifcopenshell.api.run(
        "attribute.edit_attributes",
        model,
        product=drawing,
        attributes={
            "Description": (
                f"Bonsai native Drawing / 餐厅飘窗+餐厅 / sxb010 / {name} / "
                f"{SOURCE_LABEL_ZH}"
            )
        },
    )
    reference = tool.Drawing.get_drawing_document(drawing)
    location = os.path.relpath(output_svg, Path(tool.Ifc.get_path()).resolve().parent)
    ifcopenshell.api.document.edit_reference(
        model,
        reference=reference,
        attributes={"Location": Path(location).as_posix()},
    )
    output_svg.parent.mkdir(parents=True, exist_ok=True)
    return drawing, camera, width, height, clip_end


def style_review_svg(svg_path, annotation_global_id, view):
    raw_sha = sha256(svg_path)
    ET.register_namespace("", "http://www.w3.org/2000/svg")
    ET.register_namespace("ifc", "http://www.ifcopenshell.org/ns")
    tree = ET.parse(svg_path)
    root = tree.getroot()
    target_groups = []
    geometry_count = 0
    for element in root.iter():
        if local_name(element.tag) in GEOMETRY_TAGS:
            geometry_count += 1
            existing = element.attrib.get("style", "").rstrip(";")
            style = f"stroke:{CONTEXT_GREY};stroke-width:0.22;fill:none"
            element.attrib["style"] = f"{existing};{style}" if existing else style
        attributes = {local_name(key): value for key, value in element.attrib.items()}
        classes = element.attrib.get("class", "").split()
        if (
            attributes.get("guid") == annotation_global_id
            or f"GlobalId-{annotation_global_id}" in classes
            or "review-target-sxb010" in classes
        ):
            target_groups.append(element)
    if not target_groups:
        raise RuntimeError("Bonsai SVG has no sxb010 LINEWORK Annotation group")
    target_geometry_count = 0
    for group in target_groups:
        classes = group.attrib.get("class", "").split()
        for class_name in ("review-target-sxb010", "geometry-derived-simplified-proxy"):
            if class_name not in classes:
                classes.append(class_name)
        group.attrib["class"] = " ".join(classes)
        group.attrib["data-source-kind"] = "geometry_derived_simplified_proxy"
        group.attrib["data-source-label-zh"] = SOURCE_LABEL_ZH
        group.attrib["data-target-global-id"] = TARGET_GLOBAL_ID
        group.attrib["data-reviewed-view"] = view
        for element in group.iter():
            if local_name(element.tag) not in GEOMETRY_TAGS:
                continue
            target_geometry_count += 1
            existing = element.attrib.get("style", "").rstrip(";")
            style = f"stroke:{BLUE};stroke-width:0.35;fill:none"
            element.attrib["style"] = f"{existing};{style}" if existing else style
    root.attrib["data-create-drawing-result"] = "FINISHED"
    root.attrib["data-context-colour"] = "light-grey"
    root.attrib["data-review-highlight-geometry-modified"] = "false"
    tree.write(svg_path, encoding="utf-8", xml_declaration=True)
    return {
        "bonsai_generated_sha256_before_review_style": raw_sha,
        "geometry_element_count": geometry_count,
        "target_group_count": len(target_groups),
        "target_geometry_element_count": target_geometry_count,
        "review_style_only": True,
        "geometry_modified": False,
    }


def inspect_svg(svg_path):
    root = ET.parse(svg_path).getroot()
    classes = Counter()
    geometry_count = 0
    for element in root.iter():
        classes.update(element.attrib.get("class", "").split())
        if local_name(element.tag) in GEOMETRY_TAGS:
            geometry_count += 1
    return {
        "root_tag": local_name(root.tag),
        "width": root.attrib.get("width"),
        "height": root.attrib.get("height"),
        "viewBox": root.attrib.get("viewBox"),
        "data_scale": root.attrib.get("data-scale"),
        "geometry_element_count": geometry_count,
        "projection_group_count": classes.get("projection", 0),
        "class_counts": dict(sorted(classes.items())),
    }


def persisted_annotation_state(model, annotation_global_id, view):
    annotation = model.by_guid(annotation_global_id)
    target_view = VIEW_DEFINITIONS[view]["target_view"]
    representation = next(
        (
            item
            for item in annotation.Representation.Representations
            if item.RepresentationIdentifier == "Annotation"
            and item.ContextOfItems.TargetView == target_view
            and item.RepresentationType == "GeometricCurveSet"
        ),
        None,
    )
    if representation is None:
        raise RuntimeError("persisted session lost sxb010 Annotation representation")
    paths = [path for item in representation.Items for path in item.Elements]
    return annotation, representation, paths


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
    approved_derived_sha = sha256(APPROVED_DERIVED)
    if sha256(session_ifc) != approved_derived_sha:
        raise RuntimeError("fresh drawing session must start byte-identical to approved derived IFC")
    session_before_sha = sha256(session_ifc)

    os.chdir(ROOT)
    with contextlib.suppress(Exception):
        addon_utils.disable("bl_ext.user_default.project_control", default_set=False, handle_error=None)
    load_result = bpy.ops.bim.load_project(
        filepath=str(session_ifc),
        should_start_fresh_session=True,
        use_detailed_tooltip=True,
    )
    if load_result != {"FINISHED"} or not tool.Ifc.get():
        raise RuntimeError(f"Bonsai failed to load drawing session: {load_result}")
    for handler in list(bpy.app.handlers.depsgraph_update_post):
        if getattr(handler, "__module__", "").startswith("project_control"):
            bpy.app.handlers.depsgraph_update_post.remove(handler)

    model = tool.Ifc.get()
    target = model.by_guid(TARGET_GLOBAL_ID)
    dining_bay = model.by_guid(DINING_BAY_GLOBAL_ID)
    dining = model.by_guid(DINING_GLOBAL_ID)
    target_obj = tool.Ifc.get_object(target)
    dining_bay_obj = tool.Ifc.get_object(dining_bay)
    dining_obj = tool.Ifc.get_object(dining)
    if any(item is None for item in (target, dining_bay, dining, target_obj, dining_bay_obj, dining_obj)):
        raise RuntimeError("sxb010 target or dining spaces are missing")
    if {dining_bay.LongName, dining.LongName} != ROOM_NAMES:
        raise RuntimeError("dining-space identities drifted")

    candidate = json.loads(CANDIDATE.read_text(encoding="utf-8"))
    audit = json.loads(SIMPLIFICATION_AUDIT.read_text(encoding="utf-8"))
    paths = candidate["views"][view]["proxy_paths_mm"]
    if len(paths) != EXPECTED_PATH_COUNTS[view]:
        raise RuntimeError("approved candidate merged-line path count mismatch")
    view_audit = audit["views"][view]
    if (
        view_audit["slat_count_before"] != 45
        or view_audit["slat_centerline_count_after"] != 45
        or view_audit["outer_envelope_delta_mm"] != [0.0, 0.0, 0.0, 0.0]
        or not view_audit["pass"]
    ):
        raise RuntimeError("line-simplification audit gate failed")
    if view == "front":
        slats = paths[1:46]
        if len(slats) != 45 or any(
            len(path) != 2
            or abs(path[0][1] - path[1][1]) > 0.000001
            or abs(abs(path[1][0] - path[0][0]) - 1804.4116) > 0.01
            for path in slats
        ):
            raise RuntimeError("45 retained front slat centrelines drifted")

    crop_bbox = union_bbox(world_bbox(dining_bay_obj), world_bbox(dining_obj))
    records = context_elements(crop_bbox)
    originals = [record[0] for record in records]
    target_count_before = originals.count(target)
    include_elements = [element for element in originals if element != target]
    if target_count_before != 1:
        raise RuntimeError(f"expected target once before suppression, got {target_count_before}")
    centres = source_coordinate_centres(candidate)
    definition = VIEW_DEFINITIONS[view]
    output_svg = output_dir / f"{definition['drawing_name']}.svg"
    evidence_path = output_dir / f"{definition['drawing_name']}-create-drawing-evidence.json"
    drawing, camera, width, height, clip_end = add_drawing(
        model, definition, crop_bbox, include_elements, output_svg
    )
    override = view3d_override()
    with bpy.context.temp_override(**override):
        activate_result = bpy.ops.bim.activate_drawing(
            drawing=drawing.id(), should_view_from_camera=False
        )
    if activate_result != {"FINISHED"}:
        raise RuntimeError(f"failed to activate Drawing: {activate_result}")
    annotation, _, annotation_edge_count = add_linework_annotation(
        model, drawing, target, target_obj, view, paths, centres
    )
    cprops = tool.Drawing.get_camera_props(camera)
    cprops.has_annotation = True
    pset_data = ifcopenshell.util.element.get_pset(drawing, "EPset_Drawing")
    ifcopenshell.api.pset.edit_pset(
        model,
        pset=model.by_id(pset_data["id"]),
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
    if create_result != {"FINISHED"} or not output_svg.is_file() or output_svg.stat().st_size == 0:
        raise RuntimeError(f"Bonsai Create Drawing failed: {create_result}")

    review_style = style_review_svg(output_svg, annotation.GlobalId, view)
    svg = inspect_svg(output_svg)
    if svg["root_tag"] != "svg" or svg["geometry_element_count"] == 0:
        raise RuntimeError("Bonsai SVG structure gate failed")
    cache_path = output_svg.parent / "cache" / f"{output_svg.stem}-linework.svg"
    if not cache_path.is_file() or cache_path.stat().st_size == 0:
        raise RuntimeError("Bonsai linework cache is missing")

    temporary = session_ifc.with_suffix(".ifc.next")
    model.write(str(temporary))
    reopened = ifcopenshell.open(temporary)
    reopened_drawing = reopened.by_guid(drawing.GlobalId)
    reopened_annotation, _, reopened_paths = persisted_annotation_state(
        reopened, annotation.GlobalId, view
    )
    if reopened_drawing is None or reopened_annotation is None:
        raise RuntimeError("persisted session lost Drawing or LINEWORK Annotation")
    if len(reopened_paths) != EXPECTED_PATH_COUNTS[view]:
        raise RuntimeError("persisted merged-line path count drifted")
    os.replace(temporary, session_ifc)
    reloaded = ifcopenshell.open(session_ifc)
    reloaded_drawing = reloaded.by_guid(drawing.GlobalId)
    _, _, reloaded_paths = persisted_annotation_state(reloaded, annotation.GlobalId, view)
    if reloaded_drawing is None or len(reloaded_paths) != EXPECTED_PATH_COUNTS[view]:
        raise RuntimeError("post-persist reload verification failed")

    drawing_pset = ifcopenshell.util.element.get_pset(drawing, "EPset_Drawing")
    drawing_include = drawing_pset.get("Include", "").split(",")
    if drawing_include.count(TARGET_GLOBAL_ID) != 0:
        raise RuntimeError("Drawing Include still contains the target Body")
    class_counts = Counter(element.is_a() for element in include_elements)
    evidence = {
        "schema_version": 1,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "task": "Bonsai Create Drawing for approved sxb010 in dining-bay and dining-room context",
        "view": view,
        "formal_ifc": str(formal_ifc),
        "formal_ifc_sha256": sha256(formal_ifc),
        "formal_ifc_bytes_unchanged": sha256(formal_ifc) == FORMAL_SHA256,
        "approved_derived_ifc": str(APPROVED_DERIVED),
        "approved_derived_ifc_sha256": approved_derived_sha,
        "drawing_session_ifc": str(session_ifc),
        "drawing_session_sha256_before": session_before_sha,
        "drawing_session_sha256_after": sha256(session_ifc),
        "drawing_session_is_safe_per_view_copy": True,
        "source_kind": "geometry_derived_simplified_proxy",
        "source_label_zh": SOURCE_LABEL_ZH,
        "candidate_sha256": sha256(CANDIDATE),
        "line_simplification_audit_sha256": sha256(SIMPLIFICATION_AUDIT),
        "line_simplification": {
            "slat_count_before": view_audit["slat_count_before"],
            "slat_centerline_count_after": view_audit["slat_centerline_count_after"],
            "outer_envelope_delta_mm": view_audit["outer_envelope_delta_mm"],
            "path_count": EXPECTED_PATH_COUNTS[view],
            "persisted_path_count": len(reloaded_paths),
        },
        "target": {
            "global_id": TARGET_GLOBAL_ID,
            "ifc_class": target.is_a(),
            "name": target.Name,
            "object_type": target.ObjectType,
            "bbox_m": [list(world_bbox(target_obj)[0]), list(world_bbox(target_obj)[1])],
            "include_count_before_suppression": target_count_before,
            "include_count_after_suppression": drawing_include.count(TARGET_GLOBAL_ID),
            "duplicate_target_is_real_second_instance": False,
        },
        "spaces": [
            {"global_id": dining_bay.GlobalId, "long_name": dining_bay.LongName},
            {"global_id": dining.GlobalId, "long_name": dining.LongName},
        ],
        "context": {
            "project_context_retained": True,
            "crop_bbox_m": [list(crop_bbox[0]), list(crop_bbox[1])],
            "include_count": len(include_elements),
            "include_ifc_class_counts": dict(sorted(class_counts.items())),
            "context_colour": "light_grey",
            "target_colour": "blue",
        },
        "drawing": {
            "id": drawing.id(),
            "global_id": drawing.GlobalId,
            "name": drawing.Name,
            "object_type": drawing.ObjectType,
            "epset_drawing": drawing_pset,
            "document_reference_location": tool.Drawing.get_drawing_document(drawing).Location,
        },
        "linework_annotation": {
            "global_id": annotation.GlobalId,
            "ifc_class": annotation.is_a(),
            "predefined_type": ifcopenshell.util.element.get_predefined_type(annotation),
            "name": annotation.Name,
            "associated_to_drawing_group": True,
            "target_global_id": TARGET_GLOBAL_ID,
            "geometry_source": "approved proxy_paths_mm",
            "edge_count": annotation_edge_count,
        },
        "camera": {
            "type": camera.data.type,
            "matrix_world": [[float(value) for value in row] for row in camera.matrix_world],
            "width_m": width,
            "height_m": height,
            "clip_start_m": camera.data.clip_start,
            "clip_end_m": clip_end,
            "resolution": [
                bpy.context.scene.render.resolution_x,
                bpy.context.scene.render.resolution_y,
            ],
        },
        "create_drawing": {
            "operator": "bpy.ops.bim.create_drawing",
            "arguments": {"print_all": False, "open_viewer": False, "sync": False},
            "result": sorted(create_result),
            "linework_mode": "OPENCASCADE",
            "target_view": definition["target_view"],
        },
        "svg": {
            "path": str(output_svg),
            "bytes": output_svg.stat().st_size,
            "sha256": sha256(output_svg),
            **review_style,
            **svg,
        },
        "linework_cache": {
            "path": str(cache_path),
            "bytes": cache_path.stat().st_size,
            "sha256": sha256(cache_path),
        },
        "persistence": {
            "method": "IfcOpenShell write to per-view session IFC",
            "reloaded_from_persisted_session": True,
            "drawing_found_after_reload": True,
            "annotation_path_count_after_reload": len(reloaded_paths),
        },
        "versions": {
            "blender": bpy.app.version_string,
            "ifcopenshell": ifcopenshell.version,
            "ifc_schema": model.schema,
            "bonsai_generator": "Bonsai 0.8.4 bpy.ops.bim.create_drawing",
        },
        "pass": True,
    }
    evidence_path.write_text(
        json.dumps(evidence, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print(
        json.dumps(
            {
                "view": view,
                "drawing": drawing.Name,
                "create_drawing_result": sorted(create_result),
                "svg": str(output_svg),
                "svg_sha256": sha256(output_svg),
                "session_ifc": str(session_ifc),
                "session_ifc_sha256": sha256(session_ifc),
                "evidence": str(evidence_path),
                "pass": True,
            },
            ensure_ascii=False,
        )
    )


if __name__ == "__main__":
    main()
