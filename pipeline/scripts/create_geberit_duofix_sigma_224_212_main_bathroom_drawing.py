#!/usr/bin/env python3
"""Create one native Bonsai Drawing for the main-bathroom Duofix instance.

Run inside Blender with a per-view copy of the approved derived IFC.  The copy
keeps the full project context while mapping exactly one approved native-DWG
view to Bonsai's Body + TargetView selection contract.  The formal IFC and the
approved derived IFC remain read-only.
"""

from __future__ import annotations

import hashlib
import json
import os
import sys
import xml.etree.ElementTree as ET
import contextlib
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import bpy
import addon_utils
import ifcopenshell
import ifcopenshell.api
import ifcopenshell.util.element
import ifcopenshell.util.placement
from bonsai import tool
from bonsai.core import drawing as core_drawing
from mathutils import Matrix, Vector


ROOT = Path(__file__).resolve().parents[2]
PRODUCT_DIR = ROOT / "output/review/highpoly-types/geberit-duofix-sigma-224-212"
FORMAL_IFC = ROOT / "2504 GBTB Yanlord Zhuhai.ifc"
FORMAL_SHA256 = "7a50b87e8f48a7c2bfbaa9f1dabdfca7155fa684b2325c66aa1ae15c8ab0a25c"
APPROVED_DERIVED = PRODUCT_DIR / "Geberit-Duofix-Sigma-224-212-derived-drawing.ifc"
APPROVED_DERIVED_SHA256 = "c1a0d6af0ef0897fb6a53c3b70dd2a8ae35d8a63deca06bbc9379e3decbe8897"
CANDIDATE = PRODUCT_DIR / "candidate-representations.json"
TARGET_GLOBAL_ID = "3pvAlH5C14v8uVEJ1LmK8M"
SOURCE_REPRESENTATIVE_GLOBAL_ID = "3hgNkx97vCTOC2eewCpMNk"
ROOM_GLOBAL_ID = "3a4COIs5X7lgDirMBDT4Vs"
ROOM_NAME = "主卫湿区"
ARTICLE = "224.212.00.2"
BLUE = "#1677c8"
EXPECTED_PATH_COUNTS = {"plan": 253, "front": 790, "side": 225}
VIEW_DEFINITIONS = {
    "plan": {
        "drawing_name": "GEBERIT-224212-MAIN-BATH-PLAN",
        "target_view": "PLAN_VIEW",
        "location_hint": "PLAN",
    },
    "front": {
        "drawing_name": "GEBERIT-224212-MAIN-BATH-FRONT",
        "target_view": "ELEVATION_VIEW",
        "location_hint": "EAST",
    },
    "side": {
        "drawing_name": "GEBERIT-224212-MAIN-BATH-LEFT-SIDE",
        "target_view": "ELEVATION_VIEW",
        "location_hint": "SOUTH",
    },
}
BOUNDARY_CLASSES = {
    "IfcWall",
    "IfcWallStandardCase",
    "IfcSlab",
    "IfcBeam",
    "IfcCovering",
    "IfcDoor",
    "IfcWindow",
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
        boundary = element.is_a() in BOUNDARY_CLASSES
        if centre_in_room or boundary:
            records.append((element, obj, obj_minimum, obj_maximum))
    records.sort(key=lambda item: item[0].GlobalId)
    return records


def body_target_context(model, target_view):
    matches = [
        context
        for context in model.by_type("IfcGeometricRepresentationSubContext")
        if context.ContextType == "Model"
        and context.ContextIdentifier == "Body"
        and context.TargetView == target_view
    ]
    if len(matches) != 1:
        raise RuntimeError(f"expected one Model/Body/{target_view} context, got {len(matches)}")
    return matches[0]


def mapped_body_translation(element):
    for representation in element.Representation.Representations:
        if (
            representation.RepresentationIdentifier == "Body"
            and representation.RepresentationType == "MappedRepresentation"
        ):
            mapped_item = next(
                (item for item in representation.Items if item.is_a("IfcMappedItem")),
                None,
            )
            if mapped_item is not None:
                matrix = ifcopenshell.util.placement.get_mappeditem_transformation(mapped_item)
                return tuple(float(matrix[axis][3]) for axis in range(3))
    raise RuntimeError(f"{element.GlobalId} has no mapped Body translation")


def instance_mapping_offset(model, target):
    source = model.by_guid(SOURCE_REPRESENTATIVE_GLOBAL_ID)
    source_translation = mapped_body_translation(source)
    target_translation = mapped_body_translation(target)
    delta = tuple(target_translation[axis] - source_translation[axis] for axis in range(3))
    expected = (61.16010284, -65.40448331, 326.0)
    if any(abs(delta[axis] - expected[axis]) > 0.001 for axis in range(3)):
        raise RuntimeError(f"main-bathroom mapped Body offset drifted: {delta}")
    return source, source_translation, target_translation, delta


def adjusted_coordinates_mm(view, first, second, mapping_offset_mm):
    offset_x, offset_y, offset_z = mapping_offset_mm
    if view == "plan":
        return (float(first) + offset_x, float(second) + offset_y, 0.0)
    if view == "front":
        return (float(first) + offset_x, 0.0, float(second) + offset_z)
    return (0.0, float(first) + offset_y, float(second) + offset_z)


def official_curve_representation(model, context, identifier, view, paths, mapping_offset_mm):
    polylines = []
    for path in paths:
        points = []
        for first, second in path:
            coordinates = adjusted_coordinates_mm(view, first, second, mapping_offset_mm)
            points.append(model.create_entity("IfcCartesianPoint", Coordinates=coordinates))
        if len(points) >= 2:
            polylines.append(model.create_entity("IfcPolyline", Points=points))
    if len(polylines) != EXPECTED_PATH_COUNTS[view]:
        raise RuntimeError(f"{view} native-DWG path count drifted")
    curve_set = model.create_entity("IfcGeometricCurveSet", Elements=polylines)
    return model.create_entity(
        "IfcShapeRepresentation",
        ContextOfItems=context,
        RepresentationIdentifier=identifier,
        RepresentationType="GeometricCurveSet",
        Items=[curve_set],
    )


def add_official_linework_annotation(
    model, drawing, target, target_obj, view, paths, mapping_offset_mm
):
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
    annotation.Name = f"Geberit Duofix {ARTICLE} official native DWG / {view}"
    annotation.Description = (
        f"Approved official native-DWG linework overlay for {TARGET_GLOBAL_ID}; "
        "generated as a Bonsai Drawing LINEWORK annotation because the concealed "
        "cistern is removed by architectural hidden-line processing."
    )
    annotation.ObjectPlacement = target.ObjectPlacement
    obj.matrix_world = target_obj.matrix_world

    vertices = []
    edges = []
    for path in paths:
        path_indices = []
        for first, second in path:
            coordinates_mm = adjusted_coordinates_mm(
                view, first, second, mapping_offset_mm
            )
            coordinates = tuple(value / 1000 for value in coordinates_mm)
            path_indices.append(len(vertices))
            vertices.append(coordinates)
        edges.extend(zip(path_indices, path_indices[1:]))
    mesh = obj.data
    mesh.clear_geometry()
    mesh.from_pydata(vertices, edges, [])
    mesh.update()

    representation = official_curve_representation(
        model, context, "Annotation", view, paths, mapping_offset_mm
    )
    annotation.Representation = model.create_entity(
        "IfcProductDefinitionShape", Representations=[representation]
    )
    pset = ifcopenshell.api.pset.add_pset(model, product=annotation, name="EPset_Annotation")
    ifcopenshell.api.pset.edit_pset(
        model,
        pset=pset,
        properties={
            "Classes": "official-native-dwg review-target-duofix",
            "TargetGlobalId": TARGET_GLOBAL_ID,
            "ArticleNumber": ARTICLE,
            "SourceKind": "native_dwg",
            "InstanceMappingOffsetMillimetres": ",".join(
                f"{value:.6f}" for value in mapping_offset_mm
            ),
        },
    )
    return annotation, representation, len(edges)


def path_bounds(paths):
    first = [point[0] for path in paths for point in path]
    second = [point[1] for path in paths for point in path]
    return (min(first), min(second)), (max(first), max(second))


def toilet_outlet_reference(front_paths):
    candidates = []
    for path in front_paths:
        if len(path) < 8:
            continue
        first = [point[0] for point in path]
        second = [point[1] for point in path]
        width = max(first) - min(first)
        height = max(second) - min(second)
        centre = ((min(first) + max(first)) / 2, (min(second) + max(second)) / 2)
        closed = Vector(path[0]).xy.copy()
        if (
            (closed - Vector(path[-1]).xy).length < 0.5
            and 80 <= width <= 120
            and abs(width - height) < 0.01
            and centre[1] < -600
        ):
            candidates.append((width, centre))
    if not candidates:
        raise RuntimeError("official front DWG toilet outlet reference not found")
    return max(candidates)[1]


def flush_panel_reference(front_paths):
    rectangle_lines = []
    for path in front_paths:
        if len(path) != 2:
            continue
        first = [point[0] for point in path]
        second = [point[1] for point in path]
        width = max(first) - min(first)
        height = max(second) - min(second)
        if (abs(width - 344.0) < 0.01 and height < 0.01) or (
            abs(height - 290.0) < 0.01 and width < 0.01
        ):
            rectangle_lines.extend(path)
    if len(rectangle_lines) != 8:
        raise RuntimeError("official front DWG flush-panel reference rectangle drifted")
    first = [point[0] for point in rectangle_lines]
    second = [point[1] for point in rectangle_lines]
    return {
        "bounds_mm": [[min(first), min(second)], [max(first), max(second)]],
        "centre_mm": [(min(first) + max(first)) / 2, (min(second) + max(second)) / 2],
    }


def mapped_world_point(target_obj, view, source_point, mapping_offset_mm):
    coordinates_mm = adjusted_coordinates_mm(
        view, source_point[0], source_point[1], mapping_offset_mm
    )
    local_metres = Vector(tuple(value / 1000 for value in coordinates_mm))
    return tuple(float(value) for value in (target_obj.matrix_world @ local_metres))


def elevation_alignment(model, target_obj, context_records, all_views, mapping_offset_mm):
    front_paths = all_views["front"]["official_native_dwg_paths_mm"]
    (_, source_minimum_z), (_, source_maximum_z) = path_bounds(front_paths)
    target_minimum, target_maximum = world_bbox(target_obj)
    target_xy = target_obj.matrix_world.translation.xy

    def containing_floor_record(ifc_class, maximum_height):
        matches = []
        for element, _, minimum, maximum in context_records:
            if (
                element.is_a() == ifc_class
                and minimum[0] <= target_xy.x <= maximum[0]
                and minimum[1] <= target_xy.y <= maximum[1]
                and maximum[2] - minimum[2] <= maximum_height
            ):
                matches.append((element, minimum, maximum))
        return min(matches, key=lambda record: abs(record[2][2])) if matches else None

    finish = containing_floor_record("IfcCovering", 0.10)
    slab = containing_floor_record("IfcSlab", 1.0)
    ffl = next(storey for storey in model.by_type("IfcBuildingStorey") if storey.Name == "FFL")
    ffl_z = float(
        ifcopenshell.util.placement.get_local_placement(ffl.ObjectPlacement)[2][3]
        / 1000
    )
    support_foot = mapped_world_point(
        target_obj, "front", (0.0, source_minimum_z), mapping_offset_mm
    )[2]
    frame_top = mapped_world_point(
        target_obj, "front", (0.0, source_maximum_z), mapping_offset_mm
    )[2]
    outlet_source = toilet_outlet_reference(front_paths)
    outlet_world = mapped_world_point(
        target_obj, "front", outlet_source, mapping_offset_mm
    )
    flush_reference = flush_panel_reference(front_paths)
    flush_world = mapped_world_point(
        target_obj, "front", flush_reference["centre_mm"], mapping_offset_mm
    )
    return {
        "ffl_elevation_m": ffl_z,
        "nearby_finish": {
            "global_id": finish[0].GlobalId if finish else None,
            "ifc_class": finish[0].is_a() if finish else None,
            "top_elevation_m": finish[2][2] if finish else None,
        },
        "nearby_structural_slab": {
            "global_id": slab[0].GlobalId if slab else None,
            "top_elevation_m": slab[2][2] if slab else None,
            "bottom_elevation_m": slab[1][2] if slab else None,
        },
        "project_target_body_bbox_m": [list(target_minimum), list(target_maximum)],
        "support_foot_elevation_m": support_foot,
        "support_foot_residual_to_target_body_minimum_mm": (
            support_foot - target_minimum[2]
        )
        * 1000,
        "support_foot_offset_from_ffl_mm": (support_foot - ffl_z) * 1000,
        "support_foot_offset_from_nearby_finish_mm": (
            (support_foot - finish[2][2]) * 1000 if finish else None
        ),
        "frame_top_elevation_m": frame_top,
        "frame_top_residual_to_target_body_maximum_mm": (
            frame_top - target_maximum[2]
        )
        * 1000,
        "toilet_outlet_centre": {
            "source_local_mm": list(outlet_source),
            "project_world_m": list(outlet_world),
            "elevation_from_ffl_mm": (outlet_world[2] - ffl_z) * 1000,
        },
        "flush_panel_reference": {
            "source_bounds_mm": flush_reference["bounds_mm"],
            "source_centre_mm": flush_reference["centre_mm"],
            "project_world_centre_m": list(flush_world),
            "elevation_from_ffl_mm": (flush_world[2] - ffl_z) * 1000,
        },
    }


def add_drawing(model, definition, room_bbox, include_elements, output_svg):
    name = definition["drawing_name"]
    target_view = definition["target_view"]
    matches = [
        drawing
        for drawing in model.by_type("IfcAnnotation")
        if drawing.ObjectType == "DRAWING" and drawing.Name == name
    ]
    if matches:
        raise RuntimeError(f"drawing already exists in fresh view session: {name}")

    minimum, maximum = room_bbox
    centre_x = (minimum[0] + maximum[0]) / 2
    centre_y = (minimum[1] + maximum[1]) / 2
    centre_z = 1.20
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
        raise RuntimeError(f"expected one drawing, got {len(created)}")
    drawing = created[0]
    core_drawing.update_drawing_name(tool.Ifc, tool.Drawing, drawing=drawing, name=name)
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
    ifcopenshell.api.run(
        "attribute.edit_attributes",
        model,
        product=drawing,
        attributes={
            "Description": (
                f"Bonsai native Drawing / {ROOM_NAME} / Geberit Duofix {ARTICLE} / "
                f"{definition['drawing_name']} / official native DWG selected by project owner"
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


def add_review_highlight(svg_path, target_global_id, annotation_global_id, view):
    raw_sha = sha256(svg_path)
    ET.register_namespace("", "http://www.w3.org/2000/svg")
    ET.register_namespace("ifc", "http://www.ifcopenshell.org/ns")
    tree = ET.parse(svg_path)
    root = tree.getroot()
    target_groups = []
    for element in root.iter():
        attributes = {local_name(key): value for key, value in element.attrib.items()}
        classes = element.attrib.get("class", "").split()
        if (
            attributes.get("guid") == annotation_global_id
            or f"GlobalId-{annotation_global_id}" in classes
            or "review-target-duofix" in classes
        ):
            target_groups.append(element)
    if not target_groups:
        raise RuntimeError("Bonsai SVG has no main-bathroom Duofix group")
    target_geometry_count = 0
    for group in target_groups:
        classes = group.attrib.get("class", "").split()
        for class_name in ("review-target-duofix", "official-native-dwg"):
            if class_name not in classes:
                classes.append(class_name)
        group.attrib["class"] = " ".join(classes)
        group.attrib["data-source-kind"] = "native_dwg"
        group.attrib["data-article-number"] = ARTICLE
        group.attrib["data-target-global-id"] = target_global_id
        group.attrib["data-reviewed-view"] = view
        for element in group.iter():
            if local_name(element.tag) not in GEOMETRY_TAGS:
                continue
            target_geometry_count += 1
            existing = element.attrib.get("style", "").rstrip(";")
            highlight = f"stroke:{BLUE};stroke-width:0.35;fill:none"
            element.attrib["style"] = f"{existing};{highlight}" if existing else highlight
    root.attrib["data-create-drawing-result"] = "FINISHED"
    root.attrib["data-review-highlight-geometry-modified"] = "false"
    tree.write(svg_path, encoding="utf-8", xml_declaration=True)
    return {
        "bonsai_generated_sha256_before_review_style": raw_sha,
        "target_group_count": len(target_groups),
        "target_geometry_element_count": target_geometry_count,
        "review_style_only": True,
        "geometry_modified": False,
    }


def inspect_target_layers(svg_path, target_global_id, annotation_global_id):
    root = ET.parse(svg_path).getroot()
    target_ifc_groups = []
    annotation_geometry = []
    for element in root.iter():
        attributes = {local_name(key): value for key, value in element.attrib.items()}
        classes = element.attrib.get("class", "").split()
        if attributes.get("guid") == target_global_id:
            target_ifc_groups.append(element)
        if (
            f"GlobalId-{annotation_global_id}" in classes
            and local_name(element.tag) in GEOMETRY_TAGS
        ):
            annotation_geometry.append(element)
    target_ifc_geometry_count = sum(
        1
        for group in target_ifc_groups
        for element in group.iter()
        if local_name(element.tag) in GEOMETRY_TAGS
    )
    return {
        "target_ifc_projection_group_count": len(target_ifc_groups),
        "target_ifc_projection_geometry_count": target_ifc_geometry_count,
        "official_annotation_geometry_count": len(annotation_geometry),
        "black_target_body_policy": "excluded_from_drawing_include",
        "blue_target_layer_policy": "approved_official_native_dwg_ifc_linework_annotation_only",
    }


def inspect_svg(svg_path):
    root = ET.parse(svg_path).getroot()
    classes = Counter()
    geometry_count = 0
    references = []
    for element in root.iter():
        for class_name in element.attrib.get("class", "").split():
            classes[class_name] += 1
        if local_name(element.tag) in GEOMETRY_TAGS:
            geometry_count += 1
        for key, value in element.attrib.items():
            if local_name(key) == "href":
                references.append(value)
    return {
        "root_tag": local_name(root.tag),
        "width": root.attrib.get("width"),
        "height": root.attrib.get("height"),
        "viewBox": root.attrib.get("viewBox"),
        "data_scale": root.attrib.get("data-scale"),
        "geometry_element_count": geometry_count,
        "projection_group_count": classes.get("projection", 0),
        "class_counts": dict(sorted(classes.items())),
        "reference_count": len(references),
        "references": references,
    }


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
    if sha256(APPROVED_DERIVED) != APPROVED_DERIVED_SHA256:
        raise RuntimeError("approved derived IFC hash mismatch")
    if sha256(session_ifc) != APPROVED_DERIVED_SHA256:
        raise RuntimeError("fresh drawing session must start byte-identical to approved derived IFC")
    session_before_sha = sha256(session_ifc)

    os.chdir(ROOT)
    with contextlib.suppress(Exception):
        addon_utils.disable("bl_ext.user_default.project_control", default_set=False, handle_error=None)
    result = bpy.ops.bim.load_project(
        filepath=str(session_ifc),
        should_start_fresh_session=True,
        use_detailed_tooltip=True,
    )
    if result != {"FINISHED"} or not tool.Ifc.get():
        raise RuntimeError(f"Bonsai failed to load drawing session: {result}")
    model = tool.Ifc.get()
    target = model.by_guid(TARGET_GLOBAL_ID)
    room = model.by_guid(ROOM_GLOBAL_ID)
    target_obj = tool.Ifc.get_object(target)
    room_obj = tool.Ifc.get_object(room)
    if target is None or room is None or target_obj is None or room_obj is None:
        raise RuntimeError("main-bathroom target or space missing")
    if room.LongName != ROOM_NAME:
        raise RuntimeError("main-bathroom space identity drifted")

    candidate = json.loads(CANDIDATE.read_text(encoding="utf-8"))
    candidate_view = candidate["views"][view]
    paths = candidate_view["official_native_dwg_paths_mm"]
    if len(paths) != EXPECTED_PATH_COUNTS[view]:
        raise RuntimeError("approved candidate native-DWG path count mismatch")
    source_representative, source_mapping, target_mapping, mapping_offset = (
        instance_mapping_offset(model, target)
    )
    room_bbox = world_bbox(room_obj)
    context_records = room_elements(room_bbox)
    original_context_elements = [record[0] for record in context_records]
    original_target_include_count = original_context_elements.count(target)
    context_elements = [
        element for element in original_context_elements if element != target
    ]
    if original_target_include_count != 1:
        raise RuntimeError(
            f"expected main-bathroom target once before duplicate suppression, got {original_target_include_count}"
        )
    type_object = next(relation.RelatingType for relation in target.IsTypedBy)
    same_type_instances = sorted(
        {
            related.GlobalId
            for relation in type_object.Types
            for related in relation.RelatedObjects
        }
    )
    same_type_context_instances = sorted(
        element.GlobalId
        for element in original_context_elements
        if any(relation.RelatingType == type_object for relation in element.IsTypedBy)
    )
    if same_type_context_instances != [TARGET_GLOBAL_ID]:
        raise RuntimeError(
            f"unexpected second Duofix instance in main-bath context: {same_type_context_instances}"
        )
    alignment = elevation_alignment(
        model, target_obj, context_records, candidate["views"], mapping_offset
    )
    definition = VIEW_DEFINITIONS[view]
    output_svg = output_dir / f"{definition['drawing_name']}.svg"
    evidence_path = output_dir / f"{definition['drawing_name']}-create-drawing-evidence.json"
    drawing, camera, width, height, clip_end = add_drawing(
        model,
        definition,
        room_bbox,
        context_elements,
        output_svg,
    )
    result = bpy.ops.bim.activate_drawing(drawing=drawing.id(), should_view_from_camera=False)
    if result != {"FINISHED"}:
        raise RuntimeError(f"failed to activate drawing: {result}")
    annotation, official_representation, official_edge_count = add_official_linework_annotation(
        model,
        drawing,
        target,
        target_obj,
        view,
        paths,
        mapping_offset,
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
    create_result = bpy.ops.bim.create_drawing(
        print_all=False,
        open_viewer=False,
        sync=False,
    )
    if create_result != {"FINISHED"} or not output_svg.is_file() or output_svg.stat().st_size == 0:
        raise RuntimeError(f"Bonsai Create Drawing failed: {create_result}")

    highlight = add_review_highlight(
        output_svg, TARGET_GLOBAL_ID, annotation.GlobalId, view
    )
    svg = inspect_svg(output_svg)
    target_layers = inspect_target_layers(
        output_svg, TARGET_GLOBAL_ID, annotation.GlobalId
    )
    if svg["root_tag"] != "svg" or svg["geometry_element_count"] == 0 or svg["projection_group_count"] == 0:
        raise RuntimeError("Bonsai SVG structure gate failed")
    if (
        target_layers["target_ifc_projection_group_count"] != 0
        or target_layers["official_annotation_geometry_count"] == 0
    ):
        raise RuntimeError(f"duplicate target layer policy failed: {target_layers}")
    cache_path = output_svg.parent / "cache" / f"{output_svg.stem}-linework.svg"
    if not cache_path.is_file() or cache_path.stat().st_size == 0:
        raise RuntimeError("Bonsai linework cache is missing")

    temporary = session_ifc.with_suffix(".ifc.next")
    model.write(str(temporary))
    reopened = ifcopenshell.open(temporary)
    reopened_drawing = reopened.by_guid(drawing.GlobalId)
    reopened_annotation = reopened.by_guid(annotation.GlobalId)
    reopened_representation = next(
        (
            representation
            for representation in reopened_annotation.Representation.Representations
            if representation.RepresentationIdentifier == "Annotation"
            and representation.ContextOfItems.TargetView == definition["target_view"]
            and representation.RepresentationType == "GeometricCurveSet"
        ),
        None,
    )
    if reopened_drawing is None or reopened_annotation is None or reopened_representation is None:
        raise RuntimeError("persisted drawing session lost the Drawing or official representation")
    persisted_path_count = sum(len(item.Elements) for item in reopened_representation.Items)
    if persisted_path_count != EXPECTED_PATH_COUNTS[view]:
        raise RuntimeError("persisted official representation path count drifted")
    persisted_points = [
        tuple(float(value) for value in point.Coordinates)
        for item in reopened_representation.Items
        for path in item.Elements
        for point in path.Points
    ]
    expected_points = [
        adjusted_coordinates_mm(view, first, second, mapping_offset)
        for path in paths
        for first, second in path
    ]
    persisted_bounds = [
        [min(point[axis] for point in persisted_points) for axis in range(3)],
        [max(point[axis] for point in persisted_points) for axis in range(3)],
    ]
    expected_bounds = [
        [min(point[axis] for point in expected_points) for axis in range(3)],
        [max(point[axis] for point in expected_points) for axis in range(3)],
    ]
    persisted_coordinate_residual_mm = max(
        abs(persisted_bounds[bound][axis] - expected_bounds[bound][axis])
        for bound in range(2)
        for axis in range(3)
    )
    if persisted_coordinate_residual_mm > 0.000001:
        raise RuntimeError("persisted official representation coordinate mapping drifted")
    os.replace(temporary, session_ifc)

    reference = tool.Drawing.get_drawing_document(drawing)
    pset = ifcopenshell.util.element.get_pset(drawing, "EPset_Drawing")
    class_counts = Counter(element.is_a() for element in context_elements)
    drawing_include = pset.get("Include", "").split(",")
    if drawing_include.count(TARGET_GLOBAL_ID) != 0 or len(drawing_include) != len(context_elements):
        raise RuntimeError("Drawing Include still contains the suppressed target Body")
    evidence = {
        "schema_version": 2,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "task": "Bonsai Create Drawing from approved Geberit 224.212.00.2 native-DWG source in the actual main-bathroom context",
        "view": view,
        "formal_ifc": str(formal_ifc),
        "formal_ifc_sha256": sha256(formal_ifc),
        "formal_ifc_bytes_unchanged": sha256(formal_ifc) == FORMAL_SHA256,
        "approved_derived_ifc": str(APPROVED_DERIVED),
        "approved_derived_ifc_sha256": sha256(APPROVED_DERIVED),
        "drawing_session_ifc": str(session_ifc),
        "drawing_session_sha256_before": session_before_sha,
        "drawing_session_sha256_after": sha256(session_ifc),
        "drawing_session_is_safe_per_view_copy": True,
        "source_kind": "native_dwg",
        "article_number": ARTICLE,
        "native_dwg_code": candidate_view["native_dwg_code"],
        "native_dwg_sha256": candidate_view["native_dwg_sha256"],
        "official_native_dwg_path_count": EXPECTED_PATH_COUNTS[view],
        "persisted_native_dwg_path_count": persisted_path_count,
        "official_native_dwg_edge_count": official_edge_count,
        "persisted_adjusted_coordinate_bounds_mm": persisted_bounds,
        "expected_adjusted_coordinate_bounds_mm": expected_bounds,
        "persisted_coordinate_residual_mm": persisted_coordinate_residual_mm,
        "target": {
            "global_id": TARGET_GLOBAL_ID,
            "ifc_class": target.is_a(),
            "type_name": next(relation.RelatingType.Name for relation in target.IsTypedBy),
            "same_type_project_instance_count": len(same_type_instances),
            "same_type_project_global_ids": same_type_instances,
            "same_type_main_bath_context_count_before_suppression": len(
                same_type_context_instances
            ),
            "same_type_main_bath_context_global_ids_before_suppression": same_type_context_instances,
            "target_include_count_before_suppression": original_target_include_count,
            "target_include_count_after_suppression": drawing_include.count(
                TARGET_GLOBAL_ID
            ),
            "duplicate_target_is_real_second_instance": False,
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
            "non_target_project_context_element_count": len(context_elements),
            "drawing_element_policy": (
                "65 non-target project elements plus the current Drawing's official "
                "LINEWORK annotation; target Body excluded to prevent double expression"
            ),
        },
        "instance_mapping_alignment": {
            "source_representative_global_id": source_representative.GlobalId,
            "source_mapped_body_translation_mm": list(source_mapping),
            "target_mapped_body_translation_mm": list(target_mapping),
            "applied_delta_mm": list(mapping_offset),
            "reason": (
                "Official DWG paths were normalized in the source representative's "
                "mapped Body coordinates and must inherit the main-bath instance's "
                "MappingTarget translation delta."
            ),
        },
        "elevation_alignment": alignment,
        "drawing": {
            "id": drawing.id(),
            "global_id": drawing.GlobalId,
            "name": drawing.Name,
            "object_type": drawing.ObjectType,
            "epset_drawing": pset,
            "document_reference_location": reference.Location,
        },
        "official_linework_annotation": {
            "global_id": annotation.GlobalId,
            "ifc_class": annotation.is_a(),
            "predefined_type": ifcopenshell.util.element.get_predefined_type(annotation),
            "name": annotation.Name,
            "associated_to_drawing_group": True,
            "target_global_id": TARGET_GLOBAL_ID,
            "geometry_source": "official_native_dwg_paths_mm",
            "instance_mapping_offset_mm": list(mapping_offset),
            "hidden_line_boundary": (
                "The concealed cistern body is removed by architectural hidden-line "
                "processing; the approved source is emitted through Bonsai's native "
                "Drawing LINEWORK annotation layer."
            ),
        },
        "camera": {
            "type": camera.data.type,
            "matrix_world": [[float(value) for value in row] for row in camera.matrix_world],
            "width_m": width,
            "height_m": height,
            "clip_start_m": camera.data.clip_start,
            "clip_end_m": clip_end,
            "resolution": [bpy.context.scene.render.resolution_x, bpy.context.scene.render.resolution_y],
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
            **highlight,
            **svg,
            **target_layers,
        },
        "linework_cache": {
            "path": str(cache_path),
            "bytes": cache_path.stat().st_size,
            "sha256": sha256(cache_path),
        },
        "versions": {
            "blender": bpy.app.version_string,
            "ifcopenshell": ifcopenshell.version,
            "ifc_schema": model.schema,
            "bonsai_generator": "Bonsai 0.8.4 bim.create_drawing",
        },
        "pass": True,
    }
    evidence_path.write_text(json.dumps(evidence, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps({
        "view": view,
        "drawing": drawing.Name,
        "create_drawing_result": sorted(create_result),
        "svg": str(output_svg),
        "svg_sha256": sha256(output_svg),
        "evidence": str(evidence_path),
        "pass": True,
    }, ensure_ascii=False))


if __name__ == "__main__":
    main()
