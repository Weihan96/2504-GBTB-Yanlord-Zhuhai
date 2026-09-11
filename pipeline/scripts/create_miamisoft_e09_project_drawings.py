#!/usr/bin/env python3
"""Create approved Miami Soft E09 Drawings in the actual project context."""

from __future__ import annotations

import contextlib
import hashlib
import json
import os
import shutil
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


PRODUCT_DIR = ROOT / "output/review/highpoly-types/miamisoft-e09"
FORMAL_IFC = ROOT / "2504 GBTB Yanlord Zhuhai.ifc"
FORMAL_SHA256 = "7a50b87e8f48a7c2bfbaa9f1dabdfca7155fa684b2325c66aa1ae15c8ab0a25c"
DERIVED_IFC = PRODUCT_DIR / "Baxter-Miami-Soft-E09-derived-drawing.ifc"
SESSION_BLEND = PRODUCT_DIR / "Baxter-Miami-Soft-E09-project-drawings.blend"
OUTPUT_DIR = PRODUCT_DIR / "bonsai-drawings/living-area"
EVIDENCE = OUTPUT_DIR / "MIAMISOFT-E09-create-drawing-evidence.json"
APPROVAL = ROOT / "pipeline/decisions/miamisoft-e09-drawing-approval.json"
CANDIDATE = PRODUCT_DIR / "candidate-representations.json"
TARGET_GLOBAL_ID = "3osoWufdD1mhDDAM6lcix4"
SOURCE_DWG_SHA256 = "efc3315feed0a2dac8938aebbd7336856dabba6721395e4520783c39ba0cc51a"
SOURCE_KIND = "native_dwg"
SOURCE_LABEL_ZH = "Baxter 官方 Miami Soft E09 原生二维 DWG"
BLUE = "#1677c8"
GREY = "#a3abb3"
HIDDEN_DASH = "2.4,1.5"
EXPECTED_PATH_COUNTS = {"plan": 64, "front": 130, "side": 85}
VIEW_DEFINITIONS = {
    "plan": {
        "drawing_name": "MIAMISOFT-E09-LIVING-AREA-PLAN",
        "target_view": "PLAN_VIEW",
        "location_hint": "PLAN",
    },
    "front": {
        "drawing_name": "MIAMISOFT-E09-LIVING-AREA-FRONT",
        "target_view": "ELEVATION_VIEW",
        "location_hint": "EAST",
    },
    "side": {
        "drawing_name": "MIAMISOFT-E09-LIVING-AREA-SIDE",
        "target_view": "ELEVATION_VIEW",
        "location_hint": "SOUTH",
    },
}
GEOMETRY_TAGS = {"path", "polyline", "polygon", "line", "circle", "ellipse", "rect"}


shared.PRODUCT_DIR = PRODUCT_DIR
shared.FORMAL_IFC = FORMAL_IFC
shared.FORMAL_SHA256 = FORMAL_SHA256
shared.TARGET_GLOBAL_ID = TARGET_GLOBAL_ID
shared.ROOM_GLOBAL_ID = "target-centred-project-context"
shared.ROOM_NAME = "起居区（Miami E09 代表实例周边）"
shared.ARTICLE = "Miami Soft E09 dx/r"
shared.DRAWING_PRODUCT_CODE = "MIAMISOFT-E09"
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


def path_bbox(path):
    return (
        [min(point[0] for point in path), min(point[1] for point in path)],
        [max(point[0] for point in path), max(point[1] for point in path)],
    )


def split_visible_hidden(view: str, paths: list):
    if view != "side":
        return paths, [], None
    visible, hidden = [], []
    for path in paths:
        minimum, maximum = path_bbox(path)
        centre = [
            (minimum[0] + maximum[0]) / 2,
            (minimum[1] + maximum[1]) / 2,
        ]
        size = [maximum[0] - minimum[0], maximum[1] - minimum[1]]
        is_hidden_roller = (
            abs(centre[0] - 347.3536) <= 0.01
            and abs(centre[1] - 665.0000) <= 0.01
            and abs(size[0] - 269.9901) <= 0.02
            and abs(size[1] - 269.9950) <= 0.02
        )
        (hidden if is_hidden_roller else visible).append(path)
    if len(hidden) != 1 or len(visible) != EXPECTED_PATH_COUNTS[view] - 1:
        raise RuntimeError("official hidden roller path identification drifted")
    hidden_path_index = paths.index(hidden[0])
    first_hidden_edge = sum(max(0, len(path) - 1) for path in paths[:hidden_path_index])
    hidden_edge_count = max(0, len(hidden[0]) - 1)
    return paths, hidden, [first_hidden_edge, first_hidden_edge + hidden_edge_count]


def curve_representation(model, context, view, paths, hidden: bool):
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
        Name=(
            "Baxter Miami Soft official hidden native-DWG linework"
            if hidden
            else "Baxter Miami Soft official visible native-DWG linework"
        ),
        CurveFont=None,
        CurveWidth=model.create_entity("IfcPositiveLengthMeasure", 0.35),
        CurveColour=colour,
        ModelOrDraughting=True,
    )
    model.create_entity(
        "IfcStyledItem",
        Item=curve_set,
        Styles=[curve_style],
        Name="Blue official hidden LINEWORK" if hidden else "Blue official visible LINEWORK",
    )
    representation = model.create_entity(
        "IfcShapeRepresentation",
        ContextOfItems=context,
        RepresentationIdentifier="Annotation",
        RepresentationType="GeometricCurveSet",
        Items=[curve_set],
    )
    return representation, len(polylines), sum(max(0, len(path) - 1) for path in paths)


def add_annotation(model, drawing, target, target_obj, view, paths, hidden: bool):
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
    line_role = "hidden" if hidden else "visible"
    annotation.Name = f"Baxter Miami Soft E09 official native DWG / {view} / {line_role}"
    annotation.Description = (
        f"{SOURCE_LABEL_ZH}; approved product-level linework for {TARGET_GLOBAL_ID}; "
        "rigid reviewed placement only; high-poly component offsets are accepted and not force-fitted; "
        f"line role={line_role}."
    )
    annotation.ObjectPlacement = target.ObjectPlacement
    obj.matrix_world = target_obj.matrix_world
    vertices, edges = [], []
    for path in paths:
        indices = []
        for first, second in path:
            indices.append(len(vertices))
            vertices.append(
                tuple(value / 1000.0 for value in coordinates_mm(view, first, second))
            )
        edges.extend(zip(indices, indices[1:]))
    obj.data.clear_geometry()
    obj.data.from_pydata(vertices, edges, [])
    obj.data.update()
    representation, path_count, edge_count = curve_representation(
        model, context, view, paths, hidden
    )
    annotation.Representation = model.create_entity(
        "IfcProductDefinitionShape", Representations=[representation]
    )
    pset = ifcopenshell.api.pset.add_pset(model, product=annotation, name="EPset_Annotation")
    ifcopenshell.api.pset.edit_pset(
        model,
        pset=pset,
        properties={
            "Classes": f"review-target-miamisoft-e09 official-native-dwg {line_role}-line",
            "TargetGlobalId": TARGET_GLOBAL_ID,
            "SourceKind": SOURCE_KIND,
            "SourceLabelZh": SOURCE_LABEL_ZH,
            "SourceDwgSha256": SOURCE_DWG_SHA256,
            "OfficialNative2dDwg": True,
            "SourceScale": 1.0,
            "RigidReviewedPlacementOnly": True,
            "ProjectHighPolyComponentOffsetsAccepted": True,
            "ComponentForceFitApplied": False,
            "LineRole": line_role,
            "IfcCurveStyleColour": BLUE,
            "HiddenLinePattern": "DASHED" if hidden else "CONTINUOUS",
            "OfficialHiddenLineSourceHandle": "1115B" if view == "side" else None,
            "OfficialHiddenLinePathIndex": 70 if view == "side" else None,
            "OfficialHiddenLinePathCount": 1 if view == "side" else 0,
            "OfficialHiddenLinePattern": "DASHED" if view == "side" else "NONE",
        },
    )
    return annotation, path_count, edge_count


def style_svg(svg_path: Path, visible_guid: str, hidden_edge_range, view: str):
    raw_sha = sha256(svg_path)
    ET.register_namespace("", "http://www.w3.org/2000/svg")
    ET.register_namespace("ifc", "http://www.ifcopenshell.org/ns")
    tree = ET.parse(svg_path)
    root = tree.getroot()
    visible_groups = []
    for element in root.iter():
        attributes = {local_name(key): value for key, value in element.attrib.items()}
        classes = element.attrib.get("class", "").split()
        if attributes.get("guid") == visible_guid or f"GlobalId-{visible_guid}" in classes:
            visible_groups.append(element)
    if not visible_groups:
        raise RuntimeError("Bonsai SVG official annotation group gate failed")
    official_geometry = [
        element
        for element in root.iter()
        if local_name(element.tag) in GEOMETRY_TAGS
        and (
            {local_name(key): value for key, value in element.attrib.items()}.get("guid")
            == visible_guid
            or f"GlobalId-{visible_guid}" in element.attrib.get("class", "").split()
        )
    ]
    if not official_geometry:
        raise RuntimeError("Bonsai SVG official annotation geometry is missing")
    hidden_ids = set()
    if hidden_edge_range is not None:
        first, last = hidden_edge_range
        if last > len(official_geometry) or last <= first:
            raise RuntimeError("Bonsai SVG hidden edge range drifted")
        hidden_ids = {id(element) for element in official_geometry[first:last]}
    visible_ids = {id(element) for element in official_geometry}
    blue_visible = blue_hidden = grey = 0
    for element in root.iter():
        if local_name(element.tag) not in GEOMETRY_TAGS:
            continue
        existing = element.attrib.get("style", "").rstrip(";")
        if id(element) in hidden_ids:
            colour = (
                f"stroke:{BLUE};stroke-width:0.35;fill:none;"
                f"stroke-dasharray:{HIDDEN_DASH}"
            )
            blue_hidden += 1
        elif id(element) in visible_ids:
            colour = f"stroke:{BLUE};stroke-width:0.35;fill:none"
            blue_visible += 1
        else:
            colour = f"stroke:{GREY};stroke-width:0.22;fill:none;stroke-opacity:0.72"
            grey += 1
        element.attrib["style"] = f"{existing};{colour}" if existing else colour
    for group in visible_groups:
        classes = group.attrib.get("class", "").split()
        for class_name in (
            "review-target-miamisoft-e09",
            "official-native-dwg",
        ):
            if class_name not in classes:
                classes.append(class_name)
        group.attrib["class"] = " ".join(classes)
        group.attrib["data-source-kind"] = SOURCE_KIND
        group.attrib["data-source-dwg-sha256"] = SOURCE_DWG_SHA256
    for element in official_geometry:
        element.attrib["data-line-role"] = "hidden" if id(element) in hidden_ids else "visible"
        if id(element) in hidden_ids:
            element.attrib["data-hidden-line-pattern"] = "DASHED"
    root.attrib["data-create-drawing-result"] = "FINISHED"
    root.attrib["data-official-linework-colour"] = BLUE
    root.attrib["data-context-colour"] = GREY
    root.attrib["data-component-force-fit-applied"] = "false"
    tree.write(svg_path, encoding="utf-8", xml_declaration=True)
    if blue_visible == 0 or grey == 0 or (view == "side" and blue_hidden == 0):
        raise RuntimeError("SVG colour/hidden-line gate failed")
    return {
        "bonsai_generated_sha256_before_review_style": raw_sha,
        "visible_annotation_group_count": len(visible_groups),
        "hidden_annotation_group_count": 0,
        "blue_visible_geometry_count": blue_visible,
        "blue_hidden_geometry_count": blue_hidden,
        "grey_context_geometry_count": grey,
        "hidden_dasharray": HIDDEN_DASH if hidden_ids else None,
        "post_style_only": True,
    }


def inspect_svg(svg_path: Path, target_guid: str, annotation_guids: list[str]):
    root = ET.parse(svg_path).getroot()
    geometry_count = sum(
        local_name(element.tag) in GEOMETRY_TAGS for element in root.iter()
    )
    projection_groups = 0
    target_groups = 0
    annotation_groups = Counter()
    for element in root.iter():
        attrs = {local_name(key): value for key, value in element.attrib.items()}
        classes = element.attrib.get("class", "").split()
        guid = attrs.get("guid")
        if "projection" in classes:
            projection_groups += 1
        if guid == target_guid or f"GlobalId-{target_guid}" in classes:
            target_groups += 1
        for annotation_guid in annotation_guids:
            if guid == annotation_guid or f"GlobalId-{annotation_guid}" in classes:
                annotation_groups[annotation_guid] += 1
    return {
        "root_tag": local_name(root.tag),
        "geometry_element_count": geometry_count,
        "projection_group_count": projection_groups,
        "target_ifc_projection_group_count": target_groups,
        "annotation_group_counts": dict(annotation_groups),
        "annotation_presence_all": all(annotation_groups[guid] >= 1 for guid in annotation_guids),
        "duplicate_target_or_annotation": target_groups != 0
        or any(annotation_groups[guid] == 0 for guid in annotation_guids),
    }


def persisted_paths(annotation):
    representations = annotation.Representation.Representations
    curve_sets = [
        item
        for representation in representations
        for item in representation.Items
        if item.is_a("IfcGeometricCurveSet")
    ]
    return [curve for curve_set in curve_sets for curve in curve_set.Elements]


def view3d_override():
    """Return a real GUI VIEW_3D override for operators called by MCP."""
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


def main():
    (PRODUCT_DIR / "miamisoft-e09-create-drawing-error.log").unlink(missing_ok=True)
    if sha256(FORMAL_IFC) != FORMAL_SHA256:
        raise RuntimeError("formal IFC hash mismatch before Miami E09 workflow")
    approval = json.loads(APPROVAL.read_text(encoding="utf-8"))
    candidate = json.loads(CANDIDATE.read_text(encoding="utf-8"))
    if (
        approval.get("status") != "approved"
        or approval.get("approved_views") != ["plan", "front", "side"]
        or approval.get("derived_ifc_write_allowed") is not True
        or approval.get("formal_authoritative_ifc_write_allowed") is not False
        or candidate.get("derived_ifc_write_allowed") is not True
        or candidate.get("formal_ifc_write_allowed") is not False
        or candidate.get("source_kind") != SOURCE_KIND
        or candidate.get("representative_global_id") != TARGET_GLOBAL_ID
    ):
        raise RuntimeError("Miami E09 product-level approval/write gate failed")
    for view in VIEW_DEFINITIONS:
        if len(candidate["views"][view]["candidate_paths_mm"]) != EXPECTED_PATH_COUNTS[view]:
            raise RuntimeError(f"{view} candidate path count drifted")
        alignment = candidate["views"][view]["alignment"]
        if (
            alignment.get("uniform_scale") != 1.0
            or alignment.get("anisotropic_scale_used") is not False
            or alignment.get("source_geometry_deformed") is not False
        ):
            raise RuntimeError(f"{view} rigid candidate gate failed")

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    derived_previous_sha = sha256(DERIVED_IFC) if DERIVED_IFC.is_file() else None
    shutil.copy2(FORMAL_IFC, DERIVED_IFC)
    if sha256(DERIVED_IFC) != FORMAL_SHA256:
        raise RuntimeError("derived IFC did not start byte-identical to formal IFC")
    derived_before_sha = sha256(DERIVED_IFC)

    os.chdir(ROOT)
    with contextlib.suppress(Exception):
        addon_utils.disable("bl_ext.user_default.project_control", default_set=False, handle_error=None)
    load_result = bpy.ops.bim.load_project(
        filepath=str(DERIVED_IFC),
        should_start_fresh_session=True,
        use_detailed_tooltip=True,
    )
    if load_result != {"FINISHED"} or not tool.Ifc.get():
        raise RuntimeError(f"Bonsai failed to load derived IFC: {load_result}")
    model = tool.Ifc.get()
    target = model.by_guid(TARGET_GLOBAL_ID)
    target_obj = tool.Ifc.get_object(target) if target else None
    if target is None or target_obj is None:
        raise RuntimeError("Miami E09 representative is missing")
    target_bbox = shared.world_bbox(target_obj)
    target_matrix = [list(row) for row in target_obj.matrix_world]
    centre = [
        (target_bbox[0][axis] + target_bbox[1][axis]) / 2 for axis in range(3)
    ]
    scene_bbox = (
        (centre[0] - 3.0, centre[1] - 3.0, -0.2),
        (centre[0] + 3.0, centre[1] + 3.0, 3.0),
    )
    records = shared.room_elements(scene_bbox)
    original_context = [record[0] for record in records]
    if original_context.count(target) != 1:
        raise RuntimeError("expected exactly one Miami E09 target in scene crop")
    context_elements = [element for element in original_context if element != target]
    if len(context_elements) < 5:
        raise RuntimeError("actual project context crop is unexpectedly empty")
    context_counts = Counter(element.is_a() for element in context_elements)

    view_records = []
    for view, definition in VIEW_DEFINITIONS.items():
        paths = candidate["views"][view]["candidate_paths_mm"]
        official_paths, hidden_paths, hidden_edge_range = split_visible_hidden(view, paths)
        output_svg = OUTPUT_DIR / f"{definition['drawing_name']}.svg"
        drawing, camera, width, height, clip_end = shared.add_drawing(
            model, definition, scene_bbox, context_elements, output_svg
        )
        override = view3d_override()
        with bpy.context.temp_override(**override):
            activate_result = bpy.ops.bim.activate_drawing(
                drawing=drawing.id(), should_view_from_camera=False
            )
        if activate_result != {"FINISHED"}:
            raise RuntimeError(f"failed to activate {view} Drawing")
        visible, visible_count, visible_edges = add_annotation(
            model, drawing, target, target_obj, view, official_paths, False
        )
        hidden = None
        hidden_count = len(hidden_paths)
        hidden_edges = sum(max(0, len(path) - 1) for path in hidden_paths)
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
        style = style_svg(
            output_svg,
            visible.GlobalId,
            hidden_edge_range,
            view,
        )
        annotation_guids = [visible.GlobalId]
        svg = inspect_svg(output_svg, TARGET_GLOBAL_ID, annotation_guids)
        if (
            svg["root_tag"] != "svg"
            or svg["geometry_element_count"] == 0
            or svg["projection_group_count"] == 0
            or svg["duplicate_target_or_annotation"]
        ):
            raise RuntimeError(f"{view} SVG structure/duplicate gate failed: {svg}")
        cache_path = output_svg.parent / "cache" / f"{output_svg.stem}-linework.svg"
        if not cache_path.is_file() or cache_path.stat().st_size == 0:
            raise RuntimeError(f"{view} Bonsai linework cache missing")
        view_records.append(
            {
                "view": view,
                "drawing_global_id": drawing.GlobalId,
                "drawing_name": drawing.Name,
                "visible_annotation_global_id": visible.GlobalId,
                "hidden_annotation_global_id": hidden.GlobalId if hidden else None,
                "visible_path_count": visible_count,
                "hidden_path_count": hidden_count,
                "edge_count": visible_edges,
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
                    "matrix_world": [list(row) for row in camera.matrix_world],
                    "resolution": [
                        bpy.context.scene.render.resolution_x,
                        bpy.context.scene.render.resolution_y,
                    ],
                },
            }
        )

    temporary = DERIVED_IFC.with_suffix(".ifc.next")
    model.write(str(temporary))
    os.replace(temporary, DERIVED_IFC)
    derived_after_sha = sha256(DERIVED_IFC)
    reload_result = bpy.ops.bim.load_project(
        filepath=str(DERIVED_IFC),
        should_start_fresh_session=True,
        use_detailed_tooltip=True,
    )
    if reload_result != {"FINISHED"} or not tool.Ifc.get():
        raise RuntimeError("persisted Miami E09 IFC reload failed")
    reopened = tool.Ifc.get()
    for record in view_records:
        drawing = reopened.by_guid(record["drawing_global_id"])
        visible = reopened.by_guid(record["visible_annotation_global_id"])
        hidden = None
        if drawing is None or visible is None:
            raise RuntimeError(f"persisted {record['view']} Drawing/Annotation missing")
        persisted_visible = len(persisted_paths(visible))
        annotation_pset = ifcopenshell.util.element.get_pset(visible, "EPset_Annotation")
        persisted_hidden = int(annotation_pset.get("OfficialHiddenLinePathCount", 0))
        if (
            persisted_visible != record["visible_path_count"]
            or persisted_hidden != record["hidden_path_count"]
        ):
            raise RuntimeError(f"persisted {record['view']} path count drifted")
        record["persisted_visible_path_count"] = persisted_visible
        record["persisted_hidden_path_count"] = persisted_hidden

    bpy.ops.wm.save_as_mainfile(filepath=str(SESSION_BLEND), check_existing=False)
    if sha256(FORMAL_IFC) != FORMAL_SHA256:
        raise RuntimeError("formal IFC changed during Miami E09 workflow")

    output_views = []
    for record in view_records:
        output_views.append(
            {
                "view": record["view"],
                "drawing": {
                    "global_id": record["drawing_global_id"],
                    "name": record["drawing_name"],
                },
                "visible_annotation_global_id": record["visible_annotation_global_id"],
                "hidden_annotation_global_id": record["hidden_annotation_global_id"],
                "visible_path_count": record["visible_path_count"],
                "hidden_path_count": record["hidden_path_count"],
                "persisted_visible_path_count": record["persisted_visible_path_count"],
                "persisted_hidden_path_count": record["persisted_hidden_path_count"],
                "edge_count": record["edge_count"],
                "alignment": candidate["views"][record["view"]]["alignment"],
                "component_force_fit_applied": False,
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
                    "target_view": VIEW_DEFINITIONS[record["view"]]["target_view"],
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
    evidence = {
        "schema_version": 1,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "task": "Bonsai Create Drawing for approved Baxter Miami Soft E09 native-DWG linework in the actual project living-area context",
        "workflow": ["inspect", "plan", "execute", "persist", "reload", "verify"],
        "provider": {
            "name": "public bonsai-mcp",
            "status": "supported",
            "bridge": "127.0.0.1:9878",
            "execution": "mcp__bonsai_mcp__execute_blender_code",
        },
        "course_evidence": {
            "mode": "embedded-course-index",
            "lesson": "085000 Introduction to Drawings",
            "timestamps": ["01:59 Create Drawing", "02:06 Drawing created", "02:13 SVG inspection"],
            "course_fact": "Create Drawing generates or refreshes SVG after the Drawing camera, scale, depth and filters are configured.",
            "current_version_inference": "Bonsai 0.8.4 semantic Drawing state and bpy.ops.bim.create_drawing were used in Blender 4.5.3 LTS.",
        },
        "formal_ifc": str(FORMAL_IFC),
        "formal_ifc_sha256": sha256(FORMAL_IFC),
        "formal_ifc_bytes_unchanged": sha256(FORMAL_IFC) == FORMAL_SHA256,
        "derived_ifc": str(DERIVED_IFC),
        "derived_ifc_sha256_before": derived_before_sha,
        "derived_ifc_previous_sha256": derived_previous_sha,
        "derived_ifc_sha256_after": derived_after_sha,
        "session_blend": {
            "path": str(SESSION_BLEND),
            "bytes": SESSION_BLEND.stat().st_size,
            "sha256": sha256(SESSION_BLEND),
        },
        "save_boundary": "one Miami Soft E09 product-level full-project derived IFC copy",
        "approval_record": str(APPROVAL),
        "approval_record_sha256": sha256(APPROVAL),
        "candidate_record": str(CANDIDATE),
        "candidate_record_sha256": sha256(CANDIDATE),
        "source_kind": SOURCE_KIND,
        "source_label_zh": SOURCE_LABEL_ZH,
        "source_dwg_sha256": SOURCE_DWG_SHA256,
        "target": {
            "global_id": TARGET_GLOBAL_ID,
            "ifc_class": target.is_a(),
            "world_bbox_m": [list(target_bbox[0]), list(target_bbox[1])],
            "world_centre_m": centre,
            "placement_matrix": target_matrix,
            "component_force_fit_applied": False,
            "project_highpoly_component_offsets_accepted": True,
        },
        "scene_crop_bbox_m": [list(scene_bbox[0]), list(scene_bbox[1])],
        "context": {
            "project_context_retained": True,
            "include_count": len(context_elements),
            "include_ifc_class_counts": dict(sorted(context_counts.items())),
            "actual_target_body_suppressed_from_drawing_include": True,
            "context_colour": GREY,
        },
        "views": output_views,
        "persistence": {
            "method": "model.write to temporary then atomic replace; bpy.ops.bim.load_project reload; bpy.ops.wm.save_as_mainfile",
            "reload_result": sorted(reload_result),
            "post_reload_drawing_count": len(output_views),
            "post_reload_visible_annotation_count": len(output_views),
            "post_reload_hidden_annotation_count": 0,
            "post_reload_hidden_line_semantic_path_count": 1,
        },
        "versions": {
            "blender": bpy.app.version_string,
            "ifcopenshell": ifcopenshell.version,
            "ifc_schema": reopened.schema,
            "bonsai_generator": "Bonsai 0.8.4 bpy.ops.bim.create_drawing",
        },
        "tests": {
            "all_create_drawing_finished": all(
                view["create_drawing"]["result"] == ["FINISHED"] for view in output_views
            ),
            "all_target_body_projection_counts_zero": all(
                view["svg"]["target_ifc_projection_group_count"] == 0 for view in output_views
            ),
            "all_annotation_groups_unique": all(
                not view["svg"]["duplicate_target_or_annotation"] for view in output_views
            ),
            "side_hidden_path_count": next(
                view["hidden_path_count"] for view in output_views if view["view"] == "side"
            ),
            "side_hidden_dasharray": next(
                view["svg"]["hidden_dasharray"] for view in output_views if view["view"] == "side"
            ),
            "formal_ifc_unchanged": True,
        },
        "pass": True,
    }
    EVIDENCE.write_text(json.dumps(evidence, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(
        json.dumps(
            {
                "derived_ifc": str(DERIVED_IFC),
                "session_blend": str(SESSION_BLEND),
                "evidence": str(EVIDENCE),
                "views": {view["view"]: view["create_drawing"]["result"] for view in output_views},
                "pass": True,
            },
            ensure_ascii=False,
        )
    )


if __name__ == "__main__":
    try:
        main()
    except Exception:
        PRODUCT_DIR.mkdir(parents=True, exist_ok=True)
        (PRODUCT_DIR / "miamisoft-e09-create-drawing-error.log").write_text(
            traceback.format_exc(), encoding="utf-8"
        )
        raise
