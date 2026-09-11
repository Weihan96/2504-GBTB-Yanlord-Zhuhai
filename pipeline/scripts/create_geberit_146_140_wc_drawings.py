#!/usr/bin/env python3
"""Create approved Geberit 146.140 Drawings in its actual WC context.

Run in the Blender GUI after the approval-gated writer has created the
product-level derived IFC.  The representative Body and the second product
instance are excluded from the Drawing Include set; the exact approved native
DWG linework is inserted once as IFC LINEWORK Annotation for each view.
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


PRODUCT_DIR = ROOT / "output/review/highpoly-types/geberit-146-140"
FORMAL_IFC = ROOT / "2504 GBTB Yanlord Zhuhai.ifc"
FORMAL_SHA256 = "7a50b87e8f48a7c2bfbaa9f1dabdfca7155fa684b2325c66aa1ae15c8ab0a25c"
DERIVED_IFC = PRODUCT_DIR / "Geberit-146-140-derived-drawing.ifc"
SESSION_BLEND = PRODUCT_DIR / "Geberit-146-140-WC-project-drawings.blend"
OUTPUT_DIR = PRODUCT_DIR / "bonsai-drawings/wc"
EVIDENCE = OUTPUT_DIR / "GEBERIT-146-140-create-drawing-evidence.json"
APPROVAL = ROOT / "pipeline/decisions/geberit-146-140-drawing-approval.json"
CANDIDATE = PRODUCT_DIR / "candidate-representations.json"
WRITE_REPORT = PRODUCT_DIR / "Geberit-146-140-derived-drawing-report.json"
TARGET_GLOBAL_ID = "1rhZG98PPCSxaLeMFLTYb9"
SIBLING_GLOBAL_ID = "0UtU7yPb10ku4gsbGoM_sp"
ARTICLE = "146.140.11.1"
SOURCE_KIND = "native_dwg"
SOURCE_LABEL_ZH = "Geberit 官方归档 146.140.11.1 原生二维 DWG"
BLUE = "#1677c8"
GREY = "#a3abb3"
EXPECTED_PATH_COUNTS = {"plan": 50, "front": 82, "side": 65}
REPRESENTATION_IDS = {
    "plan": "Geberit146140Plan",
    "front": "Geberit146140Front",
    "side": "Geberit146140Side",
}
VIEW_DEFINITIONS = {
    "plan": {
        "drawing_name": "GEBERIT-146-140-WC-PLAN",
        "target_view": "PLAN_VIEW",
        "location_hint": "PLAN",
    },
    "front": {
        "drawing_name": "GEBERIT-146-140-WC-FRONT",
        "target_view": "ELEVATION_VIEW",
        "location_hint": "EAST",
    },
    "side": {
        "drawing_name": "GEBERIT-146-140-WC-SIDE",
        "target_view": "ELEVATION_VIEW",
        "location_hint": "SOUTH",
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
            model.create_entity(
                "IfcCartesianPoint",
                Coordinates=coordinates_mm(view, first, second),
            )
            for first, second in path
        ]
        if len(points) >= 2:
            polylines.append(model.create_entity("IfcPolyline", Points=points))
    if len(polylines) != EXPECTED_PATH_COUNTS[view]:
        raise RuntimeError(f"Geberit {view} approved path count drifted")
    curve_set = model.create_entity("IfcGeometricCurveSet", Elements=polylines)
    colour = model.create_entity(
        "IfcColourRgb",
        Name="Geberit official native DWG blue",
        Red=22.0 / 255.0,
        Green=119.0 / 255.0,
        Blue=200.0 / 255.0,
    )
    style = model.create_entity(
        "IfcCurveStyle",
        Name="Geberit 146.140 approved native DWG linework",
        CurveFont=None,
        CurveWidth=model.create_entity("IfcPositiveLengthMeasure", 0.35),
        CurveColour=colour,
        ModelOrDraughting=True,
    )
    model.create_entity(
        "IfcStyledItem", Item=curve_set, Styles=[style], Name="Blue approved LINEWORK"
    )
    return model.create_entity(
        "IfcShapeRepresentation",
        ContextOfItems=context,
        RepresentationIdentifier="Annotation",
        RepresentationType="GeometricCurveSet",
        Items=[curve_set],
    )


def add_official_annotation(model, drawing, target, target_obj, view, paths):
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
    annotation.Name = f"Geberit {ARTICLE} approved official linework / {view}"
    annotation.Description = (
        f"{SOURCE_LABEL_ZH}; representative {TARGET_GLOBAL_ID}; exact approved "
        "Plan/Front/Side native-DWG path set, with no scaling or deformation."
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
    annotation.Representation = model.create_entity(
        "IfcProductDefinitionShape",
        Representations=[curve_representation(model, context, view, paths)],
    )
    pset = ifcopenshell.api.pset.add_pset(
        model, product=annotation, name="EPset_Annotation"
    )
    ifcopenshell.api.pset.edit_pset(
        model,
        pset=pset,
        properties={
            "Classes": "review-target-geberit146140 native-dwg approved-linework",
            "TargetGlobalId": TARGET_GLOBAL_ID,
            "ExcludedSiblingGlobalId": SIBLING_GLOBAL_ID,
            "ArticleNumber": ARTICLE,
            "SourceKind": SOURCE_KIND,
            "SourceLabelZh": SOURCE_LABEL_ZH,
            "LineRole": "approved-official-native-dwg",
            "UniformScale": 1.0,
            "SourceGeometryDeformed": False,
            "IfcCurveStyleColour": BLUE,
        },
    )
    return annotation, len(paths), sum(max(0, len(path) - 1) for path in paths)


def view3d_override():
    for window in bpy.context.window_manager.windows:
        for area in window.screen.areas:
            if area.type != "VIEW_3D":
                continue
            region = next((item for item in area.regions if item.type == "WINDOW"), None)
            if region:
                return {
                    "window": window,
                    "screen": window.screen,
                    "area": area,
                    "region": region,
                    "scene": bpy.context.scene,
                }
    raise RuntimeError("Bonsai Drawing requires a real VIEW_3D area")


def element_matches(element, guid):
    attrs = {local_name(key): value for key, value in element.attrib.items()}
    classes = element.attrib.get("class", "").split()
    return attrs.get("guid") == guid or f"GlobalId-{guid}" in classes


def style_svg(svg_path: Path, annotation_guid: str):
    raw_sha = sha256(svg_path)
    ET.register_namespace("", "http://www.w3.org/2000/svg")
    ET.register_namespace("ifc", "http://www.ifcopenshell.org/ns")
    tree = ET.parse(svg_path)
    root = tree.getroot()
    annotation_groups = [element for element in root.iter() if element_matches(element, annotation_guid)]
    if not annotation_groups:
        raise RuntimeError("Geberit approved annotation group missing from Bonsai SVG")
    blue_geometry = set()
    for group in annotation_groups:
        for element in group.iter():
            if local_name(element.tag) in GEOMETRY_TAGS:
                blue_geometry.add(id(element))
    grey_count = blue_count = 0
    for element in root.iter():
        if local_name(element.tag) not in GEOMETRY_TAGS:
            continue
        is_blue = id(element) in blue_geometry or element_matches(element, annotation_guid)
        colour = (
            f"stroke:{BLUE};stroke-width:0.35;fill:none"
            if is_blue
            else f"stroke:{GREY};stroke-width:0.22;fill:none;stroke-opacity:0.72"
        )
        existing = element.attrib.get("style", "").rstrip(";")
        element.attrib["style"] = f"{existing};{colour}" if existing else colour
        blue_count += int(is_blue)
        grey_count += int(not is_blue)
    root.attrib.update(
        {
            "data-create-drawing-result": "FINISHED",
            "data-official-native-dwg-colour": BLUE,
            "data-context-colour": GREY,
            "data-target-body-projection-suppressed": "true",
            "data-sibling-product-projection-suppressed": "true",
            "data-uniform-scale": "1.0",
            "data-source-geometry-deformed": "false",
        }
    )
    tree.write(svg_path, encoding="utf-8", xml_declaration=True)
    if blue_count == 0 or grey_count == 0:
        raise RuntimeError("Geberit SVG colour gate failed")
    return {
        "bonsai_generated_sha256_before_review_style": raw_sha,
        "blue_official_geometry_count": blue_count,
        "grey_context_geometry_count": grey_count,
        "post_style_only": True,
    }


def inspect_svg(svg_path: Path, annotation_guid: str):
    root = ET.parse(svg_path).getroot()
    target_present = any(element_matches(element, TARGET_GLOBAL_ID) for element in root.iter())
    sibling_present = any(element_matches(element, SIBLING_GLOBAL_ID) for element in root.iter())
    annotation_present = any(element_matches(element, annotation_guid) for element in root.iter())
    geometry_count = sum(
        local_name(element.tag) in GEOMETRY_TAGS for element in root.iter()
    )
    projection_groups = sum(
        "projection" in element.attrib.get("class", "").split() for element in root.iter()
    )
    return {
        "root_tag": local_name(root.tag),
        "data_scale": root.attrib.get("data-scale"),
        "width": root.attrib.get("width"),
        "height": root.attrib.get("height"),
        "view_box": root.attrib.get("viewBox"),
        "geometry_element_count": geometry_count,
        "projection_group_count": projection_groups,
        "representative_body_present": target_present,
        "sibling_product_present": sibling_present,
        "official_annotation_present": annotation_present,
        "no_body_annotation_or_instance_duplicate": (
            not target_present and not sibling_present and annotation_present
        ),
        "external_reference_count": sum(
            "href" in local_name(key) and not value.startswith("#")
            for element in root.iter()
            for key, value in element.attrib.items()
        ),
    }


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
        raise RuntimeError("formal IFC hash mismatch before Geberit Drawing workflow")
    approval = json.loads(APPROVAL.read_text(encoding="utf-8"))
    candidate = json.loads(CANDIDATE.read_text(encoding="utf-8"))
    report = json.loads(WRITE_REPORT.read_text(encoding="utf-8"))
    if (
        approval.get("status") != "approved"
        or approval.get("approved_views") != ["plan", "front", "side"]
        or approval.get("derived_ifc_write_allowed") is not True
        or approval.get("formal_authoritative_ifc_write_allowed") is not False
        or candidate.get("source_kind") != SOURCE_KIND
        or candidate.get("representative_global_id") != TARGET_GLOBAL_ID
        or report.get("pass") is not True
        or report.get("representation_path_counts") != EXPECTED_PATH_COUNTS
    ):
        raise RuntimeError("Geberit product-level approval/write gate failed")
    for view in VIEW_DEFINITIONS:
        paths = candidate["views"][view]["official_native_dwg_paths_mm"]
        if len(paths) != EXPECTED_PATH_COUNTS[view]:
            raise RuntimeError(f"Geberit {view} approved candidate drifted")

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
        raise RuntimeError(f"Bonsai failed to load Geberit derived IFC: {load_result}")
    model = tool.Ifc.get()
    target = model.by_guid(TARGET_GLOBAL_ID)
    sibling = model.by_guid(SIBLING_GLOBAL_ID)
    target_obj = tool.Ifc.get_object(target) if target else None
    if target is None or target_obj is None or sibling is None:
        raise RuntimeError("Geberit representative or sibling identity is missing")
    target_bbox = shared.world_bbox(target_obj)
    centre = [(target_bbox[0][axis] + target_bbox[1][axis]) / 2 for axis in range(3)]
    scene_bbox = (
        (centre[0] - 1.6, centre[1] - 1.6, -0.2),
        (centre[0] + 1.6, centre[1] + 1.6, 3.0),
    )
    records = shared.room_elements(scene_bbox)
    original = [record[0] for record in records]
    if original.count(target) != 1:
        raise RuntimeError("expected exactly one representative Geberit in WC crop")
    same_type_in_crop = [
        element.GlobalId
        for element in original
        if element.IsTypedBy and element.IsTypedBy[0].RelatingType.Name == "Geberit 146.140"
    ]
    context_elements = [
        element for element in original if element.GlobalId not in {TARGET_GLOBAL_ID, SIBLING_GLOBAL_ID}
    ]
    if len(context_elements) < 5:
        raise RuntimeError("Geberit actual WC context crop is unexpectedly empty")
    context_counts = Counter(element.is_a() for element in context_elements)
    if not any(name.startswith("IfcWall") for name in context_counts):
        raise RuntimeError("Geberit WC context has no wall")
    if not any(name in context_counts for name in ("IfcSlab", "IfcCovering")):
        raise RuntimeError("Geberit WC context has no floor/covering")

    view_records = []
    for view, definition in VIEW_DEFINITIONS.items():
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
            raise RuntimeError(f"Geberit {view} Drawing activation failed")
        paths = candidate["views"][view]["official_native_dwg_paths_mm"]
        annotation, path_count, edge_count = add_official_annotation(
            model, drawing, target, target_obj, view, paths
        )
        cprops = tool.Drawing.get_camera_props(camera)
        cprops.has_annotation = True
        drawing_pset = ifcopenshell.util.element.get_pset(drawing, "EPset_Drawing")
        ifcopenshell.api.pset.edit_pset(
            model,
            pset=model.by_id(drawing_pset["id"]),
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
            raise RuntimeError(f"Geberit Bonsai Create Drawing failed for {view}: {create_result}")
        style = style_svg(output_svg, annotation.GlobalId)
        svg = inspect_svg(output_svg, annotation.GlobalId)
        if (
            svg["root_tag"] != "svg"
            or svg["geometry_element_count"] == 0
            or svg["projection_group_count"] == 0
            or not svg["no_body_annotation_or_instance_duplicate"]
        ):
            raise RuntimeError(f"Geberit {view} SVG/duplicate gate failed: {svg}")
        cache = output_svg.parent / "cache" / f"{output_svg.stem}-linework.svg"
        if not cache.is_file() or cache.stat().st_size == 0:
            raise RuntimeError(f"Geberit {view} Bonsai linework cache missing")
        view_records.append(
            {
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
                "create_result": sorted(create_result),
            }
        )

    temporary = DERIVED_IFC.with_suffix(".ifc.drawings-next")
    model.write(temporary)
    os.replace(temporary, DERIVED_IFC)
    derived_after_sha = sha256(DERIVED_IFC)
    reload_result = bpy.ops.bim.load_project(
        filepath=str(DERIVED_IFC),
        should_start_fresh_session=True,
        use_detailed_tooltip=True,
    )
    if reload_result != {"FINISHED"} or not tool.Ifc.get():
        raise RuntimeError("persisted Geberit IFC reload failed")
    reopened = tool.Ifc.get()
    reopened_target = reopened.by_guid(TARGET_GLOBAL_ID)
    persisted_representations = {
        item.RepresentationIdentifier: sum(
            len(curve_set.Elements)
            for curve_set in item.Items
            if curve_set.is_a("IfcGeometricCurveSet")
        )
        for item in reopened_target.Representation.Representations
        if item.RepresentationIdentifier in REPRESENTATION_IDS.values()
    }
    expected_persisted = {
        REPRESENTATION_IDS[view]: count for view, count in EXPECTED_PATH_COUNTS.items()
    }
    if persisted_representations != expected_persisted:
        raise RuntimeError("Geberit approved product representations drifted after reload")
    for record in view_records:
        drawing = reopened.by_guid(record["drawing_global_id"])
        annotation = reopened.by_guid(record["annotation_global_id"])
        if drawing is None or annotation is None or persisted_path_count(annotation) != record["path_count"]:
            raise RuntimeError(f"persisted Geberit {record['view']} Drawing/LINEWORK missing")
        record["persisted_path_count"] = persisted_path_count(annotation)
    bpy.ops.wm.save_as_mainfile(filepath=str(SESSION_BLEND), check_existing=False)
    if sha256(FORMAL_IFC) != FORMAL_SHA256:
        raise RuntimeError("formal IFC changed during Geberit Drawing workflow")

    output_views = []
    for record in view_records:
        output_views.append(
            {
                "view": record["view"],
                "drawing": {
                    "global_id": record["drawing_global_id"],
                    "name": record["drawing_name"],
                },
                "official_annotation_global_id": record["annotation_global_id"],
                "path_count": record["path_count"],
                "persisted_path_count": record["persisted_path_count"],
                "edge_count": record["edge_count"],
                "camera": record["camera"],
                "create_drawing": {
                    "operator": "bpy.ops.bim.create_drawing",
                    "arguments": {"print_all": False, "open_viewer": False, "sync": False},
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
                    "path": str(record["cache"]),
                    "bytes": record["cache"].stat().st_size,
                    "sha256": sha256(record["cache"]),
                },
            }
        )
    evidence = {
        "schema_version": 1,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "task": "approved Geberit AquaClean Sela 146.140 product-level derived IFC plus Bonsai Create Drawing Plan/Front/Side in actual WC context",
        "courseEvidence": {
            "mode": "embedded-course-index",
            "source_path": "/Users/jiaxinchen/.codex/plugins/cache/personal/bonsai-course-operator/0.1.1+codex.20260827111347/skills/operate-bonsai-from-course/evals/course-index.json",
            "lesson": "085000 Introduction to Drawings",
            "timestamps": ["01:59 Create Drawing", "02:06 Drawing created", "02:13 SVG inspection"],
            "course_fact": "Create Drawing generates or refreshes SVG after Drawing camera, scale, depth and filters are configured.",
            "current_inference": "Bonsai 0.8.4 semantic Drawing state and bpy.ops.bim.create_drawing were used in Blender 4.5.3 LTS.",
        },
        "plan": [
            "inspect approved derived IFC",
            "scope representative WC instance and project context",
            "create three Drawings",
            "insert exact official LINEWORK once",
            "run native Create Drawing",
            "persist and reload",
            "inspect SVG structure and duplicate state",
        ],
        "preState": {
            "derived_ifc": str(DERIVED_IFC),
            "derived_ifc_sha256": derived_before_sha,
            "formal_ifc_sha256": sha256(FORMAL_IFC),
            "drawing_count": 0,
            "representative_global_id": TARGET_GLOBAL_ID,
            "sibling_global_id": SIBLING_GLOBAL_ID,
        },
        "execution": {
            "provider": "public bonsai-mcp b0e67b1 contract; direct Blender GUI adapter because bridge was initially unreachable",
            "bridge_preflight": "127.0.0.1:9878 connection refused",
            "capability": "execute_blender_code adapter",
            "generator": "bpy.ops.bim.create_drawing",
            "linework_mode": "OPENCASCADE",
        },
        "persistence": {
            "method": "model.write to temporary then atomic replace; bpy.ops.bim.load_project reload; bpy.ops.wm.save_as_mainfile",
            "reload_result": sorted(reload_result),
        },
        "postState": {
            "derived_ifc": str(DERIVED_IFC),
            "derived_ifc_sha256": derived_after_sha,
            "formal_ifc_sha256": sha256(FORMAL_IFC),
            "drawing_count": len(view_records),
            "official_linework_annotation_count": len(view_records),
            "approved_product_representation_path_counts": persisted_representations,
        },
        "outputs": {
            "views": output_views,
            "session_blend": {
                "path": str(SESSION_BLEND),
                "bytes": SESSION_BLEND.stat().st_size,
                "sha256": sha256(SESSION_BLEND),
            },
        },
        "visual": {
            "preview_manifest": str(PRODUCT_DIR / "Geberit-146-140-WC-drawing-manifest.json"),
            "representative_body_displayed": False,
            "sibling_instance_displayed": False,
            "official_annotation_displayed_once_per_view": True,
        },
        "target": {
            "global_id": TARGET_GLOBAL_ID,
            "world_bbox_m": [list(target_bbox[0]), list(target_bbox[1])],
            "world_centre_m": centre,
            "same_type_instances_in_crop_before_filter": same_type_in_crop,
        },
        "context": {
            "project_context_retained": True,
            "include_count": len(context_elements),
            "include_ifc_class_counts": dict(sorted(context_counts.items())),
            "representative_body_suppressed_from_drawing_include": True,
            "sibling_product_suppressed_from_drawing_include": True,
            "official_linework_inserted_once": True,
        },
        "versions": {
            "blender": bpy.app.version_string,
            "bonsai": "0.8.4",
            "ifcopenshell": ifcopenshell.version,
            "ifc_schema": reopened.schema,
        },
        "tests": {
            "all_create_drawing_finished": all(item["create_drawing"]["result"] == ["FINISHED"] for item in output_views),
            "all_body_and_sibling_duplicates_absent": all(item["svg"]["no_body_annotation_or_instance_duplicate"] for item in output_views),
            "all_path_counts_persisted": all(item["path_count"] == item["persisted_path_count"] for item in output_views),
            "formal_ifc_unchanged": True,
        },
        "formal_ifc_bytes_unchanged": sha256(FORMAL_IFC) == FORMAL_SHA256,
        "verdict": "pass",
        "pass": True,
    }
    EVIDENCE.write_text(
        json.dumps(evidence, indent=2, ensure_ascii=False, default=str) + "\n",
        encoding="utf-8",
    )
    print(
        json.dumps(
            {
                "evidence": str(EVIDENCE),
                "derived_ifc_sha256": derived_after_sha,
                "views": {item["view"]: item["svg"] for item in output_views},
                "pass": True,
            },
            indent=2,
        )
    )


try:
    main()
except Exception:
    PRODUCT_DIR.mkdir(parents=True, exist_ok=True)
    (PRODUCT_DIR / "geberit-146-140-create-drawing-error.log").write_text(
        traceback.format_exc(), encoding="utf-8"
    )
    raise
finally:
    if not bpy.app.background:
        # Let Bonsai's load-post handlers settle before asking the GUI process
        # to exit; immediate quit during a completed reload can race Blender's
        # delayed window-manager exit handler.
        def quit_after_handlers():
            bpy.ops.wm.quit_blender()
            return None

        bpy.app.timers.register(quit_after_handlers, first_interval=1.0)
