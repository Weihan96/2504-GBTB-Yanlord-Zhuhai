#!/usr/bin/env python3
"""Create three orthographic cabinet-internal detail Drawings for TRAP01.

This script is executed inside Blender/Bonsai through the public bonsai-mcp
Provider.  It only adds Drawing/filter/LINEWORK state to the already approved
product-level derived IFC.  The formal project IFC and TRAP01 component
geometry remain read-only.
"""

from __future__ import annotations

import hashlib
import json
import os
import sys
import traceback
import xml.etree.ElementTree as ET
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import bpy
import ifcopenshell
import ifcopenshell.api
import ifcopenshell.util.element
from bonsai import tool
from bonsai.core import drawing as core_drawing
from mathutils import Matrix, Vector


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "pipeline/scripts"))
import create_gessi316_54294_main_bathroom_drawing as shared  # noqa: E402


PRODUCT_DIR = ROOT / "output/review/highpoly-types/trap01"
FORMAL_IFC = ROOT / "2504 GBTB Yanlord Zhuhai.ifc"
FORMAL_SHA256 = "7a50b87e8f48a7c2bfbaa9f1dabdfca7155fa684b2325c66aa1ae15c8ab0a25c"
DERIVED_IFC = PRODUCT_DIR / "Geberit-151.116.11.1-TRAP01-derived-drawing.ifc"
OUTPUT_DIR = PRODUCT_DIR / "bonsai-drawings/cabinet-internal-detail"
PREPERSIST = OUTPUT_DIR / "TRAP01-cabinet-internal-detail-prepersist.json"
CANDIDATE = PRODUCT_DIR / "candidate-representations.json"
TARGET_GLOBAL_ID = "2Ak2ma0lvBEA49UpplzUqi"
COMPONENT_PSET = "Pset_Trap01DrawingComponent"
SOURCE_KIND = "geometry_derived_simplified_proxy"
SOURCE_LABEL_ZH = "基于原始高模几何生成的简化图纸表达"
BLACK = "#111820"
GREY = "#a3abb3"
EXPECTED_PATH_COUNTS = {"plan": 16, "front": 1, "side": 8}
GEOMETRY_TAGS = {"path", "polyline", "polygon", "line", "circle", "ellipse", "rect"}

# These project objects conceal the trap in at least one orthographic view.
# Their GlobalIds were resolved from the derived IFC immediately before this
# Drawing batch.  They are excluded only from the detail Drawings.
OCCLUDERS = {
    "1FgLPMw$5B4wBH2ySMkXE1": "washbasin/counter assembly above the trap",
    "2YnGmmoYbF28Vy25$aNKT8": "vanity/cabinet volume intersecting the trap envelope",
    "3OVQygdDn17huGOgJJFTOY": "low cabinet partition enclosing the trap",
}

VIEW_DEFINITIONS = {
    "plan": {
        "drawing_name": "TRAP01-CABINET-INTERNAL-DETAIL-PLAN",
        "target_view": "PLAN_VIEW",
        "location_hint": "PLAN",
    },
    "front": {
        "drawing_name": "TRAP01-CABINET-INTERNAL-DETAIL-FRONT",
        "target_view": "ELEVATION_VIEW",
        "location_hint": "EAST",
    },
    "side": {
        "drawing_name": "TRAP01-CABINET-INTERNAL-DETAIL-SIDE",
        "target_view": "ELEVATION_VIEW",
        "location_hint": "SOUTH",
    },
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def local_name(value: str) -> str:
    return value.rsplit("}", 1)[-1].split(":")[-1]


def world_bbox(obj):
    points = [obj.matrix_world @ Vector(corner) for corner in obj.bound_box]
    return (
        tuple(min(point[axis] for point in points) for axis in range(3)),
        tuple(max(point[axis] for point in points) for axis in range(3)),
    )


def intersects(first, second):
    return all(first[1][axis] >= second[0][axis] and first[0][axis] <= second[1][axis] for axis in range(3))


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
        raise RuntimeError(f"TRAP01 {view} path count drifted")
    curve_set = model.create_entity("IfcGeometricCurveSet", Elements=polylines)
    colour = model.create_entity(
        "IfcColourRgb", Name="TRAP01 current installed black",
        Red=17.0 / 255.0, Green=24.0 / 255.0, Blue=32.0 / 255.0,
    )
    style = model.create_entity(
        "IfcCurveStyle", Name="TRAP01 current installed internal detail",
        CurveFont=None, CurveWidth=model.create_entity("IfcPositiveLengthMeasure", 0.35),
        CurveColour=colour, ModelOrDraughting=True,
    )
    model.create_entity("IfcStyledItem", Item=curve_set, Styles=[style], Name="Current installed internal detail")
    return model.create_entity(
        "IfcShapeRepresentation", ContextOfItems=context,
        RepresentationIdentifier="Annotation", RepresentationType="GeometricCurveSet",
        Items=[curve_set],
    )


def add_current_annotation(model, drawing, target, target_obj, view, paths, coordinate_view=None):
    coordinate_view = coordinate_view or view
    target_view = VIEW_DEFINITIONS[view]["target_view"]
    context = tool.Drawing.get_annotation_context(target_view) or tool.Drawing.create_annotation_context(target_view)
    obj = core_drawing.add_annotation(
        tool.Ifc, tool.Collector, tool.Drawing, drawing=drawing,
        object_type="LINEWORK", relating_type=None, enable_editing=False,
    )
    annotation = tool.Ifc.get_entity(obj)
    annotation.Name = f"TRAP01 current installed cabinet-internal detail / {view}"
    annotation.Description = (
        f"{SOURCE_LABEL_ZH}; current installed configuration only; cabinet occluders excluded by Drawing filter."
    )
    annotation.ObjectPlacement = target.ObjectPlacement
    obj.matrix_world = target_obj.matrix_world
    vertices, edges = [], []
    for path in paths:
        indices = []
        for first, second in path:
            indices.append(len(vertices))
            vertices.append(tuple(value / 1000.0 for value in coordinates_mm(coordinate_view, first, second)))
        edges.extend(zip(indices, indices[1:]))
    obj.data.clear_geometry()
    obj.data.from_pydata(vertices, edges, [])
    obj.data.update()
    annotation.Representation = model.create_entity(
        "IfcProductDefinitionShape",
        Representations=[curve_representation(model, context, coordinate_view, paths)],
    )
    pset = ifcopenshell.api.pset.add_pset(model, product=annotation, name="EPset_Annotation")
    ifcopenshell.api.pset.edit_pset(
        model,
        pset=pset,
        properties={
            "Classes": "review-target-trap01 current-installed cabinet-internal-detail",
            "TargetGlobalId": TARGET_GLOBAL_ID,
            "SourceKind": SOURCE_KIND,
            "SourceLabelZh": SOURCE_LABEL_ZH,
            "LineRole": "current-installed-internal-detail",
            "IfcCurveStyleColour": BLACK,
            "ReferenceExtensionDisplayed": False,
            "BlueDashedReferenceExcluded": True,
        },
    )
    return annotation


def view3d_override():
    for window in bpy.context.window_manager.windows:
        for area in window.screen.areas:
            if area.type != "VIEW_3D":
                continue
            region = next((item for item in area.regions if item.type == "WINDOW"), None)
            if region:
                return {"window": window, "screen": window.screen, "area": area, "region": region, "scene": bpy.context.scene}
    raise RuntimeError("Bonsai Drawing requires a real VIEW_3D area")


def add_detail_drawing(model, definition, centre, include_elements, output_svg):
    target_view = definition["target_view"]
    if target_view == "PLAN_VIEW":
        storey = next(item for item in model.by_type("IfcBuildingStorey") if item.Name == "FFL")
        location_hint = storey.id()
        bpy.context.scene.cursor.location = (centre[0], centre[1], 0.98)
    elif definition["location_hint"] == "EAST":
        location_hint = "EAST"
        bpy.context.scene.cursor.location = (centre[0] + 0.55, centre[1], centre[2])
    else:
        location_hint = "SOUTH"
        bpy.context.scene.cursor.location = (centre[0], centre[1] - 0.55, centre[2])

    before = {item.id() for item in model.by_type("IfcAnnotation") if item.ObjectType == "DRAWING"}
    core_drawing.add_drawing(tool.Ifc, tool.Collector, tool.Drawing, target_view=target_view, location_hint=location_hint)
    created = [
        item for item in model.by_type("IfcAnnotation")
        if item.ObjectType == "DRAWING" and item.id() not in before
    ]
    if len(created) != 1:
        raise RuntimeError(f"expected one detail Drawing, got {len(created)}")
    drawing = created[0]
    core_drawing.update_drawing_name(tool.Ifc, tool.Drawing, drawing=drawing, name=definition["drawing_name"])
    camera = tool.Ifc.get_object(drawing) or tool.Drawing.import_drawing(drawing)
    if target_view == "PLAN_VIEW":
        matrix = Matrix.Identity(4)
        matrix.translation = (centre[0], centre[1], 0.98)
        camera.matrix_world = matrix
        width, height, clip_end = 0.90, 0.80, 0.72
    elif definition["location_hint"] == "EAST":
        camera.matrix_world = tool.Drawing.generate_drawing_matrix("ELEVATION_VIEW", "EAST")
        width, height, clip_end = 0.80, 0.80, 0.95
    else:
        camera.matrix_world = tool.Drawing.generate_drawing_matrix("ELEVATION_VIEW", "SOUTH")
        width, height, clip_end = 0.90, 0.80, 0.95

    camera.data.type = "ORTHO"
    camera.data.clip_start = 0.002
    camera.data.clip_end = clip_end
    cprops = tool.Drawing.get_camera_props(camera)
    cprops.update_props = False
    cprops.camera_type = "ORTHO"
    cprops.target_view = target_view
    cprops.custom_scale_numerator = "1"
    cprops.custom_scale_denominator = "5"
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

    drawing_pset = ifcopenshell.util.element.get_pset(drawing, "EPset_Drawing")
    ifcopenshell.api.pset.edit_pset(
        model,
        pset=model.by_id(drawing_pset["id"]),
        properties={
            "TargetView": target_view,
            "Scale": "1/5",
            "HumanScale": "1:5",
            "HasUnderlay": False,
            "HasLinework": True,
            "HasAnnotation": False,
            "GlobalReferencing": True,
            "DPI": 300,
            "LineworkMode": "OPENCASCADE",
            "FillMode": "NONE",
            "CutMode": "BISECT",
            "Include": ",".join(element.GlobalId for element in include_elements),
            "Exclude": ",".join(OCCLUDERS),
        },
    )
    drawing.Description = (
        "Bonsai native orthographic cabinet-internal detail / TRAP01 current installed configuration / "
        "occluding basin, vanity volume and low partition excluded"
    )
    reference = tool.Drawing.get_drawing_document(drawing)
    relative_output = os.path.relpath(output_svg, Path(tool.Ifc.get_path()).resolve().parent)
    ifcopenshell.api.document.edit_reference(model, reference=reference, attributes={"Location": Path(relative_output).as_posix()})
    output_svg.parent.mkdir(parents=True, exist_ok=True)
    return drawing, camera, width, height, clip_end


def style_svg(svg_path: Path, annotation_guid: str):
    raw_sha = sha256(svg_path)
    ET.register_namespace("", "http://www.w3.org/2000/svg")
    ET.register_namespace("ifc", "http://www.ifcopenshell.org/ns")
    tree = ET.parse(svg_path)
    root = tree.getroot()
    parents = {child: parent for parent in root.iter() for child in parent}
    target_groups = []
    for element in root.iter():
        attrs = {local_name(key): value for key, value in element.attrib.items()}
        classes = element.attrib.get("class", "").split()
        if attrs.get("guid") == annotation_guid or f"GlobalId-{annotation_guid}" in classes:
            target_groups.append(element)
    if not target_groups:
        raise RuntimeError("TRAP01 detail annotation group missing from Bonsai SVG")
    # Keep current-installed black linework visibly above retained grey context.
    for group in target_groups:
        parent = parents.get(group)
        if parent is not None:
            parent.remove(group)
            parent.append(group)
    black = grey = 0
    for element in root.iter():
        if local_name(element.tag) not in GEOMETRY_TAGS:
            continue
        attrs = {local_name(key): value for key, value in element.attrib.items()}
        classes = element.attrib.get("class", "").split()
        is_target = attrs.get("guid") == annotation_guid or f"GlobalId-{annotation_guid}" in classes
        colour = f"stroke:{BLACK};stroke-width:0.45;fill:none" if is_target else f"stroke:{GREY};stroke-width:0.20;fill:none;stroke-opacity:0.62"
        existing = element.attrib.get("style", "").rstrip(";")
        element.attrib["style"] = f"{existing};{colour}" if existing else colour
        black += int(is_target)
        grey += int(not is_target)
    root.attrib.update({
        "data-create-drawing-result": "FINISHED",
        "data-trap01-detail": "cabinet-internal-orthographic",
        "data-trap01-current-installed-colour": BLACK,
        "data-context-colour": GREY,
        "data-blue-dashed-reference-displayed": "false",
        "data-occluder-global-ids": ",".join(OCCLUDERS),
    })
    tree.write(svg_path, encoding="utf-8", xml_declaration=True)
    if black == 0 or grey == 0:
        raise RuntimeError("TRAP01 detail SVG colour gate failed")
    return {"raw_sha256": raw_sha, "black_geometry_count": black, "grey_geometry_count": grey}


def inspect_svg(svg_path: Path, target_guid: str, annotation_guid: str):
    root = ET.parse(svg_path).getroot()

    def has_identity(element, guid):
        attrs = {local_name(key): value for key, value in element.attrib.items()}
        classes = element.attrib.get("class", "").split()
        return attrs.get("guid") == guid or f"GlobalId-{guid}" in classes

    target_groups = int(any(has_identity(element, target_guid) for element in root.iter()))
    annotation_groups = int(any(has_identity(element, annotation_guid) for element in root.iter()))
    excluded_groups = {guid: int(any(has_identity(element, guid) for element in root.iter())) for guid in OCCLUDERS}
    geometry_count = sum(local_name(element.tag) in GEOMETRY_TAGS for element in root.iter())
    return {
        "root_tag": local_name(root.tag),
        "data_scale": root.attrib.get("data-scale"),
        "view_box": root.attrib.get("viewBox"),
        "geometry_element_count": geometry_count,
        "target_ifc_body_group_count": target_groups,
        "current_annotation_group_count": annotation_groups,
        "excluded_occluder_group_counts": excluded_groups,
        "no_target_or_annotation_duplicate": target_groups == 0 and annotation_groups == 1,
        "all_occluders_absent": all(count == 0 for count in excluded_groups.values()),
        "external_reference_count": sum(
            "href" in local_name(key) and not value.startswith("#")
            for element in root.iter() for key, value in element.attrib.items()
        ),
    }


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    if sha256(FORMAL_IFC) != FORMAL_SHA256:
        raise RuntimeError("formal IFC hash mismatch before TRAP01 detail Drawing workflow")
    if Path(tool.Ifc.get_path()).resolve() != DERIVED_IFC.resolve():
        raise RuntimeError(f"wrong IFC loaded: {tool.Ifc.get_path()}")
    model = tool.Ifc.get()
    existing_names = {
        item.Name for item in model.by_type("IfcAnnotation")
        if item.ObjectType == "DRAWING" and item.Name in {definition["drawing_name"] for definition in VIEW_DEFINITIONS.values()}
    }
    if existing_names:
        raise RuntimeError(f"detail Drawings already exist: {sorted(existing_names)}")
    target = model.by_guid(TARGET_GLOBAL_ID)
    target_obj = tool.Ifc.get_object(target) if target else None
    if target is None or target_obj is None:
        raise RuntimeError("TRAP01 target is missing")
    candidate = json.loads(CANDIDATE.read_text(encoding="utf-8"))
    if candidate.get("representative_global_id") != TARGET_GLOBAL_ID or candidate.get("source_kind") != SOURCE_KIND:
        raise RuntimeError("TRAP01 candidate identity/source gate failed")
    for view in VIEW_DEFINITIONS:
        if len(candidate["views"][view]["proxy_paths_mm"]) != EXPECTED_PATH_COUNTS[view]:
            raise RuntimeError(f"TRAP01 {view} candidate path count drifted")
    target_bbox = world_bbox(target_obj)
    centre = [(target_bbox[0][axis] + target_bbox[1][axis]) / 2 for axis in range(3)]
    local_crop = (
        (centre[0] - 0.65, centre[1] - 0.65, 0.20),
        (centre[0] + 0.65, centre[1] + 0.65, 1.10),
    )
    context_elements = []
    for element in model.by_type("IfcElement"):
        if element == target or element.GlobalId in OCCLUDERS:
            continue
        obj = tool.Ifc.get_object(element)
        if obj is None:
            continue
        try:
            if intersects(world_bbox(obj), local_crop):
                context_elements.append(element)
        except Exception:
            continue
    if len(context_elements) < 5:
        raise RuntimeError("TRAP01 retained detail context is unexpectedly sparse")
    missing_occluders = [guid for guid in OCCLUDERS if model.by_guid(guid) is None]
    if missing_occluders:
        raise RuntimeError(f"resolved occluders missing: {missing_occluders}")
    components = [item for item in model.by_type("IfcAnnotation") if item.ObjectType == "DRAWING_COMPONENT"]
    roles = {
        ifcopenshell.util.element.get_pset(item, COMPONENT_PSET).get("ComponentRole"): item
        for item in components
    }
    if set(roles) != {"fixed_body", "horizontal_adjustable", "vertical_adjustable"}:
        raise RuntimeError("TRAP01 adjustable component semantics drifted")
    pre_state = {
        "derived_ifc": str(DERIVED_IFC),
        "derived_ifc_sha256": sha256(DERIVED_IFC),
        "formal_ifc_sha256": sha256(FORMAL_IFC),
        "drawing_count": len([item for item in model.by_type("IfcAnnotation") if item.ObjectType == "DRAWING"]),
        "component_count": len(components),
    }
    override = view3d_override()
    view_records = []
    for view, definition in VIEW_DEFINITIONS.items():
        output_svg = OUTPUT_DIR / f"{definition['drawing_name']}.svg"
        drawing, camera, width, height, clip_end = add_detail_drawing(
            model, definition, centre, context_elements, output_svg
        )
        with bpy.context.temp_override(**override):
            activate = bpy.ops.bim.activate_drawing(drawing=drawing.id(), should_view_from_camera=False)
        if activate != {"FINISHED"}:
            raise RuntimeError(f"TRAP01 detail {view} activation failed")
        annotation = add_current_annotation(
            model, drawing, target, target_obj, view, candidate["views"][view]["proxy_paths_mm"]
        )
        cprops = tool.Drawing.get_camera_props(camera)
        cprops.has_annotation = True
        drawing_pset = ifcopenshell.util.element.get_pset(drawing, "EPset_Drawing")
        ifcopenshell.api.pset.edit_pset(
            model, pset=model.by_id(drawing_pset["id"]), properties={"HasAnnotation": True}
        )
        dprops = tool.Drawing.get_document_props()
        dprops.should_use_underlay_cache = False
        dprops.should_use_linework_cache = False
        dprops.should_use_annotation_cache = False
        with bpy.context.temp_override(**override):
            create_result = bpy.ops.bim.create_drawing(print_all=False, open_viewer=False, sync=False)
        if create_result != {"FINISHED"} or not output_svg.is_file() or output_svg.stat().st_size == 0:
            raise RuntimeError(f"TRAP01 detail {view} Create Drawing failed: {create_result}")
        style = style_svg(output_svg, annotation.GlobalId)
        svg = inspect_svg(output_svg, TARGET_GLOBAL_ID, annotation.GlobalId)
        if (
            svg["root_tag"] != "svg" or svg["geometry_element_count"] == 0
            or not svg["no_target_or_annotation_duplicate"] or not svg["all_occluders_absent"]
            or svg["external_reference_count"] != 0
        ):
            raise RuntimeError(f"TRAP01 detail {view} SVG gate failed: {svg}")
        cache = output_svg.parent / "cache" / f"{output_svg.stem}-linework.svg"
        if not cache.is_file() or cache.stat().st_size == 0:
            raise RuntimeError(f"TRAP01 detail {view} linework cache missing")
        view_records.append({
            "view": view,
            "drawing_global_id": drawing.GlobalId,
            "drawing_name": drawing.Name,
            "annotation_global_id": annotation.GlobalId,
            "path_count": EXPECTED_PATH_COUNTS[view],
            "output_svg": str(output_svg),
            "output_svg_bytes": output_svg.stat().st_size,
            "output_svg_sha256": sha256(output_svg),
            "linework_cache": str(cache),
            "linework_cache_sha256": sha256(cache),
            "style": style,
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
    evidence = {
        "schema_version": 1,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "status": "awaiting_save_ifc_file_and_reload_verification",
        "task": "TRAP01 cabinet-internal orthographic Plan/Front/Side detail Drawings with occluders filtered",
        "courseEvidence": {
            "evidence_mode": "embedded-course-index",
            "lesson": "085000 Introduction to Drawings",
            "timestamps": [
                "01:10 camera boundary/scale",
                "01:52 drawing depth",
                "01:59-02:13 Create Drawing and SVG inspection",
                "02:36-03:16 Element Filters and Include whitelist result",
            ],
            "course_fact": "Configure the orthographic Drawing camera, depth and element filters before Create Drawing; an Include filter constrains the SVG result.",
            "screen_observation": None,
            "provenance_note": "Private course screenshots were not directly viewed; only the packaged embedded-course-index and hashes were used.",
            "current_inference": "Bonsai 0.8.4 uses semantic EPset_Drawing Include/Exclude state and bpy.ops.bim.create_drawing in Blender 4.5.3 LTS.",
        },
        "plan": [
            "inspect approved derived IFC and component semantics",
            "resolve cabinet occluders by GlobalId",
            "create tight 1:5 orthographic detail cameras",
            "apply Include whitelist and explicit Exclude filter",
            "add current-installed LINEWORK once per view",
            "run native Create Drawing",
            "persist with Provider save_ifc_file",
            "reload and verify IFC/SVG/rendered output",
        ],
        "preState": pre_state,
        "execution": {
            "provider": "public bonsai-mcp",
            "provider_contract_status": "supported",
            "capability": "execute_blender_code",
            "generator": "bpy.ops.bim.create_drawing",
            "linework_mode": "OPENCASCADE",
            "blender": bpy.app.version_string,
            "bonsai": "0.8.4",
            "ifcopenshell": ifcopenshell.version,
        },
        "save_boundary": {
            "derived_ifc": str(DERIVED_IFC),
            "formal_ifc": str(FORMAL_IFC),
            "geometry_mutation_allowed": False,
            "drawing_filter_annotation_only": True,
        },
        "target": {
            "global_id": TARGET_GLOBAL_ID,
            "world_bbox_m": [list(target_bbox[0]), list(target_bbox[1])],
            "world_centre_m": centre,
        },
        "filters": {
            "include_count": len(context_elements),
            "include_global_ids": [element.GlobalId for element in context_elements],
            "include_ifc_class_counts": dict(sorted(Counter(element.is_a() for element in context_elements).items())),
            "excluded_occluders": [{"global_id": guid, "reason": reason} for guid, reason in OCCLUDERS.items()],
        },
        "adjustable_components": {
            role: {"global_id": item.GlobalId, "pset": ifcopenshell.util.element.get_pset(item, COMPONENT_PSET)}
            for role, item in roles.items()
        },
        "outputs": {"views": view_records},
        "formal_ifc_bytes_unchanged": sha256(FORMAL_IFC) == FORMAL_SHA256,
    }
    PREPERSIST.write_text(json.dumps(evidence, indent=2, ensure_ascii=False, default=str) + "\n", encoding="utf-8")
    print(json.dumps({"prepersist": str(PREPERSIST), "views": view_records, "ready_for_save_ifc_file": True}, indent=2))


if __name__ == "__main__":
    try:
        main()
    except Exception:
        (PRODUCT_DIR / "trap01-internal-detail-error.log").write_text(traceback.format_exc(), encoding="utf-8")
        raise
