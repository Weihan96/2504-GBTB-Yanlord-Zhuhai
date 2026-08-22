#!/usr/bin/env python3
"""Build native Bonsai interior-elevation Drawings from the DWG view register.

Run this file inside the already-open Blender/Bonsai session.  It creates real
IfcAnnotation/ObjectType=DRAWING cameras and asks Bonsai's ``bim.create_drawing``
operator to produce vector SVG linework.  It never renders a Workbench image.

The project contains deeply nested furniture BReps which can overflow the
OpenCASCADE SVG serialiser stack, while Blender 4.5 Freestyle is unstable for
some rooms.  The generator supports both native Bonsai linework modes.  In
OpenCASCADE mode it mechanically excludes only representations above a stated
node-count safety threshold and records every exclusion in the source report.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import os
import xml.etree.ElementTree as ET
from collections import Counter, deque
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

import bpy
import ifcopenshell
import ifcopenshell.api
import ifcopenshell.util.element
from bonsai import tool
from bonsai.core import drawing as core_drawing
from mathutils import Vector


REGISTER = Path("pipeline/decisions/int1-elevation-view-register.csv")
OUTPUT_DIR = Path(
    os.environ.get("INT1_BONSAI_OUTPUT_DIR", "drawings/elevations/native")
)
REPORT_DIR = Path(
    os.environ.get("INT1_BONSAI_REPORT_DIR", "build/int1/native-bonsai")
)
PROJECT_ROOT = Path(__file__).resolve().parents[2]
CAMERA_Z_M = 1.25
VIEW_BOTTOM_M = -0.10
VIEW_TOP_M = 2.60
HORIZONTAL_MARGIN_M = 0.30
WALL_BACK_MARGIN_M = 0.30
CLIP_START_M = 0.002
INTEGER_TOLERANCE_MM = 0.01
MAJOR_INTEGER_RESIDUAL_MM = 0.10
LINEWORK_MODE = os.environ.get("INT1_BONSAI_LINEWORK_MODE", "FREESTYLE").upper()
MAX_OCC_REPRESENTATION_NODES = 18_000
REPRESENTATION_NODE_COUNTS: dict[str, int] = {}

DIRECTION = {
    "+Y": {"hint": "SOUTH", "suffix": "PY", "axis": 1, "sign": 1},
    "+X": {"hint": "WEST", "suffix": "PX", "axis": 0, "sign": 1},
    "-Y": {"hint": "NORTH", "suffix": "NY", "axis": 1, "sign": -1},
    "-X": {"hint": "EAST", "suffix": "NX", "axis": 0, "sign": -1},
}

BOUNDARY_CLASSES = {
    "IfcWall",
    "IfcSlab",
    "IfcBeam",
    "IfcCovering",
    "IfcDoor",
    "IfcWindow",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def source_ifc_lineage() -> dict[str, Any]:
    """Describe the authoritative IFC and the IFC actually used for linework.

    A drawing-only candidate may add target-view Curve3D representations while
    preserving the detailed MODEL_VIEW bodies from the formal IFC.  Recording
    both files prevents the candidate from being mistaken for the SSOT.
    """

    drawing_source = Path(tool.Ifc.get_path()).resolve()
    formal_source = Path(
        os.environ.get("INT1_BONSAI_FORMAL_IFC", str(drawing_source))
    ).resolve()
    if not formal_source.is_file():
        raise FileNotFoundError(formal_source)

    def display(path: Path) -> str:
        try:
            return path.relative_to(PROJECT_ROOT).as_posix()
        except ValueError:
            return str(path)

    return {
        "formal_ifc": display(formal_source),
        "formal_ifc_sha256": sha256(formal_source),
        "drawing_source_ifc": display(drawing_source),
        "drawing_source_ifc_sha256": sha256(drawing_source),
        "drawing_source_is_derived": drawing_source != formal_source,
    }


def world_bbox(obj: bpy.types.Object) -> tuple[tuple[float, float, float], tuple[float, float, float]]:
    points = [obj.matrix_world @ Vector(corner) for corner in obj.bound_box]
    return (
        tuple(min(point[i] for point in points) for i in range(3)),
        tuple(max(point[i] for point in points) for i in range(3)),
    )


def bbox_corners(
    minimum: tuple[float, float, float], maximum: tuple[float, float, float]
) -> list[Vector]:
    return [
        Vector((x, y, z))
        for x in (minimum[0], maximum[0])
        for y in (minimum[1], maximum[1])
        for z in (minimum[2], maximum[2])
    ]


def read_rows(
    project_root: Path,
    sheet_ids: set[str] | None,
    view_ids: set[str] | None = None,
) -> list[dict[str, str]]:
    with (project_root / REGISTER).open(newline="", encoding="utf-8-sig") as handle:
        rows = list(csv.DictReader(handle))
    if sheet_ids:
        rows = [row for row in rows if row["sheet_id"] in sheet_ids]
    if view_ids:
        rows = [row for row in rows if row["view_id"] in view_ids]
    rows.sort(key=lambda row: int(row["view_id"]))
    if not rows:
        raise RuntimeError("no elevation views selected")
    return rows


def drawing_name(row: dict[str, str]) -> str:
    if row.get("drawing_name"):
        return row["drawing_name"]
    suffix = DIRECTION[row["direction"]]["suffix"]
    room = row["space_reference"].split(" ", 1)[0]
    return f'{row["sheet_id"]}-{row["view_id"]}-{room}-{suffix}'


def get_space_bbox(row: dict[str, str]) -> tuple[tuple[float, float, float], tuple[float, float, float]]:
    global_ids = row.get("space_global_ids", row["space_global_id"]).split(";")
    boxes = []
    for global_id in global_ids:
        space = tool.Ifc.get().by_guid(global_id)
        obj = tool.Ifc.get_object(space)
        if not obj:
            raise RuntimeError(f"space is not loaded: {global_id}")
        boxes.append(world_bbox(obj))
    minimum = [min(box[0][axis] for box in boxes) for axis in range(3)]
    maximum = [max(box[1][axis] for box in boxes) for axis in range(3)]
    # Public unfolded elevations use an explicit plan strip.  This makes the
    # view scope auditable and avoids clipping to one room's Space boundary.
    overrides = (
        ("scope_min_x_mm", minimum, 0),
        ("scope_min_y_mm", minimum, 1),
        ("scope_max_x_mm", maximum, 0),
        ("scope_max_y_mm", maximum, 1),
    )
    for key, target, axis in overrides:
        if row.get(key):
            target[axis] = float(row[key]) / 1000
    return tuple(minimum), tuple(maximum)


def is_demolish_wall(element: ifcopenshell.entity_instance) -> bool:
    if not element.is_a("IfcWall"):
        return False
    return (
        ifcopenshell.util.element.get_pset(element, "Pset_WallCommon", "Status")
        == "DEMOLISH"
    )


def room_elements(
    room_bbox: tuple[tuple[float, float, float], tuple[float, float, float]]
) -> list[tuple[ifcopenshell.entity_instance, bpy.types.Object, tuple, tuple]]:
    minimum, maximum = room_bbox
    results = []
    for obj in bpy.data.objects:
        element = tool.Ifc.get_entity(obj)
        if (
            not element
            or not element.is_a("IfcElement")
            or element.is_a("IfcOpeningElement")
            or element.is_a("IfcVirtualElement")
            or is_demolish_wall(element)
            or obj.type != "MESH"
            or element.Representation is None
        ):
            continue
        if element.is_a("IfcBuildingElementProxy") and (element.Name or "") == "CAD":
            continue
        obj_minimum, obj_maximum = world_bbox(obj)
        overlaps = (
            obj_maximum[0] >= minimum[0]
            and obj_minimum[0] <= maximum[0]
            and obj_maximum[1] >= minimum[1]
            and obj_minimum[1] <= maximum[1]
            and obj_maximum[2] >= VIEW_BOTTOM_M
            and obj_minimum[2] <= max(VIEW_TOP_M, maximum[2])
        )
        if not overlaps:
            continue
        centre = tuple((obj_minimum[i] + obj_maximum[i]) / 2 for i in range(3))
        centre_in_room = (
            minimum[0] - 0.05 <= centre[0] <= maximum[0] + 0.05
            and minimum[1] - 0.05 <= centre[1] <= maximum[1] + 0.05
        )
        boundary = element.is_a() in BOUNDARY_CLASSES or (
            element.is_a("IfcBuildingElementProxy") and "Ceiling" in (element.Name or "")
        )
        if centre_in_room or boundary:
            results.append((element, obj, obj_minimum, obj_maximum))
    results.sort(key=lambda item: item[0].GlobalId)
    return results


def representation_node_count(element: ifcopenshell.entity_instance) -> int:
    cached = REPRESENTATION_NODE_COUNTS.get(element.GlobalId)
    if cached is not None:
        return cached
    representation = element.Representation
    seen: set[int] = set()
    queue = deque([representation])
    while queue:
        current = queue.popleft()
        if not hasattr(current, "id") or current.id() in seen:
            continue
        seen.add(current.id())
        for value in current:
            if hasattr(value, "id"):
                queue.append(value)
            elif isinstance(value, (tuple, list)):
                queue.extend(item for item in value if hasattr(item, "id"))
    count = len(seen)
    REPRESENTATION_NODE_COUNTS[element.GlobalId] = count
    return count


def linework_elements(
    elements: list[tuple[ifcopenshell.entity_instance, bpy.types.Object, tuple, tuple]],
) -> tuple[list[tuple[ifcopenshell.entity_instance, bpy.types.Object, tuple, tuple]], list[dict[str, Any]]]:
    if LINEWORK_MODE != "OPENCASCADE":
        return elements, []
    included = []
    excluded = []
    for item in elements:
        element = item[0]
        if ifcopenshell.util.representation.get_representation(
            element, "Model", "Body", "ELEVATION_VIEW"
        ):
            included.append(item)
            continue
        node_count = representation_node_count(element)
        if node_count > MAX_OCC_REPRESENTATION_NODES:
            excluded.append(
                {
                    "global_id": element.GlobalId,
                    "ifc_class": element.is_a(),
                    "name": element.Name,
                    "representation_node_count": node_count,
                    "reason": "OpenCASCADE serializer safety threshold",
                }
            )
        else:
            included.append(item)
    return included, excluded


def camera_dimensions(
    row: dict[str, str], room_bbox: tuple[tuple[float, float, float], tuple[float, float, float]]
) -> tuple[float, float, float]:
    minimum, maximum = room_bbox
    direction = DIRECTION[row["direction"]]
    anchor = (float(row["ifc_x_mm"]) / 1000, float(row["ifc_y_mm"]) / 1000)
    axis = direction["axis"]
    target = maximum[axis] if direction["sign"] > 0 else minimum[axis]
    wall_back_margin = float(row.get("wall_back_margin_m", WALL_BACK_MARGIN_M))
    horizontal_margin = float(row.get("horizontal_margin_m", HORIZONTAL_MARGIN_M))
    view_bottom = float(row.get("view_bottom_m", VIEW_BOTTOM_M))
    view_top = float(row.get("view_top_m", VIEW_TOP_M))
    clip_end = abs(target - anchor[axis]) + wall_back_margin
    projected_axis = 0 if axis == 1 else 1
    width = maximum[projected_axis] - minimum[projected_axis] + horizontal_margin
    height = view_top - view_bottom
    return width, height, clip_end


def camera_z(row: dict[str, str]) -> float:
    view_bottom = float(row.get("view_bottom_m", VIEW_BOTTOM_M))
    view_top = float(row.get("view_top_m", VIEW_TOP_M))
    return (view_bottom + view_top) / 2


def find_or_create_drawing(row: dict[str, str]) -> ifcopenshell.entity_instance:
    model = tool.Ifc.get()
    name = drawing_name(row)
    matches = [
        drawing
        for drawing in model.by_type("IfcAnnotation")
        if getattr(drawing, "ObjectType", None) == "DRAWING" and drawing.Name == name
    ]
    if len(matches) > 1:
        raise RuntimeError(f"duplicate native drawing: {name}")
    if matches:
        return matches[0]

    bpy.context.scene.cursor.location = (
        float(row["ifc_x_mm"]) / 1000,
        float(row["ifc_y_mm"]) / 1000,
        camera_z(row),
    )
    before = {
        drawing.id()
        for drawing in model.by_type("IfcAnnotation")
        if getattr(drawing, "ObjectType", None) == "DRAWING"
    }
    core_drawing.add_drawing(
        tool.Ifc,
        tool.Collector,
        tool.Drawing,
        target_view="ELEVATION_VIEW",
        location_hint=DIRECTION[row["direction"]]["hint"],
    )
    created = [
        drawing
        for drawing in model.by_type("IfcAnnotation")
        if getattr(drawing, "ObjectType", None) == "DRAWING" and drawing.id() not in before
    ]
    if len(created) != 1:
        raise RuntimeError(f"expected one native drawing, got {len(created)}")
    drawing = created[0]
    core_drawing.update_drawing_name(
        tool.Ifc, tool.Drawing, drawing=drawing, name=name
    )
    return drawing


def configure_drawing(
    project_root: Path,
    row: dict[str, str],
    drawing: ifcopenshell.entity_instance,
    room_bbox: tuple[tuple[float, float, float], tuple[float, float, float]],
    elements: list[tuple[ifcopenshell.entity_instance, bpy.types.Object, tuple, tuple]],
) -> tuple[bpy.types.Object, Path]:
    model = tool.Ifc.get()
    name = drawing_name(row)
    camera = tool.Ifc.get_object(drawing)
    if not camera:
        camera = tool.Drawing.import_drawing(drawing)

    bpy.context.scene.cursor.location = (
        float(row["ifc_x_mm"]) / 1000,
        float(row["ifc_y_mm"]) / 1000,
        camera_z(row),
    )
    camera.matrix_world = tool.Drawing.generate_drawing_matrix(
        "ELEVATION_VIEW", DIRECTION[row["direction"]]["hint"]
    )
    width, height, clip_end = camera_dimensions(row, room_bbox)
    camera.data.type = "ORTHO"
    camera.data.clip_start = CLIP_START_M
    camera.data.clip_end = clip_end

    cprops = tool.Drawing.get_camera_props(camera)
    cprops.update_props = False
    cprops.camera_type = "ORTHO"
    cprops.target_view = "ELEVATION_VIEW"
    human_scale = row.get("human_scale", "1:50")
    scale = row.get("scale", "1/50")
    enum_scale = f"{human_scale}|{scale}"
    if enum_scale == "1:30|1/30":
        cprops.custom_scale_numerator = "1"
        cprops.custom_scale_denominator = "30"
        cprops.diagram_scale = "CUSTOM"
    else:
        cprops.diagram_scale = enum_scale
    cprops.has_underlay = False
    cprops.has_linework = True
    cprops.has_annotation = False
    cprops.linework_mode = LINEWORK_MODE
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
    include = ",".join(item[0].GlobalId for item in elements)
    ifcopenshell.api.pset.edit_pset(
        model,
        pset=pset,
        properties={
            "TargetView": "ELEVATION_VIEW",
            "Scale": scale,
            "HumanScale": human_scale,
            "HasUnderlay": False,
            "HasLinework": True,
            "HasAnnotation": False,
            "GlobalReferencing": True,
            "DPI": 300,
            "LineworkMode": LINEWORK_MODE,
            "FillMode": "NONE",
            "CutMode": "BISECT",
            "Include": include,
        },
    )
    ifcopenshell.api.run(
        "attribute.edit_attributes",
        model,
        product=drawing,
        attributes={
            "Description": (
                f'{row["sheet_id"]} {row["official_title"]} / '
                f'{row["view_id"]} / {row["space_reference"]} / '
                f'{row["direction"]} / DWG handle {row["source_handle"]}'
            )
        },
    )

    relative_svg = OUTPUT_DIR / f"{name}.svg"
    output_svg = project_root / relative_svg
    current_ifc_parent = Path(tool.Ifc.get_path()).resolve().parent
    location_from_current_ifc = os.path.relpath(output_svg, current_ifc_parent)
    reference = tool.Drawing.get_drawing_document(drawing)
    ifcopenshell.api.document.edit_reference(
        model,
        reference=reference,
        attributes={"Location": Path(location_from_current_ifc).as_posix()},
    )
    output_svg.parent.mkdir(parents=True, exist_ok=True)
    return camera, output_svg


def integer_residual_mm(
    minimum: tuple[float, float, float], maximum: tuple[float, float, float]
) -> float:
    values = [coordinate * 1000 for coordinate in (*minimum, *maximum)]
    return max(abs(value - round(value)) for value in values)


def project_bbox_to_svg(
    camera: bpy.types.Object,
    minimum: tuple[float, float, float],
    maximum: tuple[float, float, float],
    width_m: float,
    height_m: float,
) -> tuple[float, float, float, float] | None:
    local = [camera.matrix_world.inverted() @ point for point in bbox_corners(minimum, maximum)]
    if all(-point.z < camera.data.clip_start or -point.z > camera.data.clip_end for point in local):
        return None
    xs = [(point.x + width_m / 2) * 20 for point in local]
    ys = [(height_m / 2 - point.y) * 20 for point in local]
    x0, x1 = max(0.0, min(xs)), min(width_m * 20, max(xs))
    y0, y1 = max(0.0, min(ys)), min(height_m * 20, max(ys))
    if x1 <= x0 or y1 <= y0:
        return None
    return x0, y0, x1 - x0, y1 - y0


def add_integer_highlights(
    svg_path: Path,
    camera: bpy.types.Object,
    elements: list[tuple[ifcopenshell.entity_instance, bpy.types.Object, tuple, tuple]],
    width_m: float,
    height_m: float,
) -> list[dict[str, Any]]:
    records = []
    parts = [
        '<g id="noninteger-highlights" fill="none">',
        '<style>.noninteger-major{stroke:#B88A5A;stroke-width:.30;}'
        '.noninteger-minor{stroke:#D8C7A1;stroke-width:.22;stroke-dasharray:1.2,.7;}'
        '.noninteger-label{font-family:Arial,sans-serif;font-size:2px;fill:#8A6746;stroke:none;}</style>',
    ]
    for element, _obj, minimum, maximum in elements:
        residual = integer_residual_mm(minimum, maximum)
        if residual <= INTEGER_TOLERANCE_MM:
            continue
        projected = project_bbox_to_svg(camera, minimum, maximum, width_m, height_m)
        if not projected:
            continue
        x, y, width, height = projected
        severity = "major" if residual > MAJOR_INTEGER_RESIDUAL_MM else "minor"
        parts.append(
            f'<rect class="noninteger-{severity}" x="{x:.3f}" y="{y:.3f}" '
            f'width="{width:.3f}" height="{height:.3f}" data-guid="{element.GlobalId}" '
            f'data-residual-mm="{residual:.6f}"/>'
        )
        records.append(
            {
                "global_id": element.GlobalId,
                "ifc_class": element.is_a(),
                "name": element.Name,
                "residual_mm": residual,
                "severity": severity,
                "projected_bbox_mm": [x, y, width, height],
            }
        )
    parts.append("</g>")
    source = svg_path.read_text(encoding="utf-8")
    source = source.replace("</svg>", "".join(parts) + "</svg>")
    svg_path.write_text(source, encoding="utf-8")
    return records


def valid_linework_cache(
    project_root: Path,
    svg_path: Path,
    drawing_name_value: str,
) -> bool:
    """Accept a Freestyle cache only when it belongs to the loaded IFC.

    A syntactically valid SVG is not sufficient: after an IFC product rename
    or geometry edit, Bonsai's cache can still contain the previous linework.
    The source report provides the IFC content hash, while GlobalId/name pairs
    in the cached SVG catch legacy reports created before this gate existed.
    """

    cache_path = svg_path.parent / "cache" / f"{svg_path.stem}-linework.svg"
    if not cache_path.is_file() or cache_path.stat().st_size < 1000:
        return False
    report_path = project_root / REPORT_DIR / f"{drawing_name_value}-source.json"
    if not report_path.is_file():
        return False
    try:
        report = json.loads(report_path.read_text(encoding="utf-8"))
    except (json.JSONDecodeError, OSError):
        return False
    if report.get("formal_ifc_sha256_before_save") != sha256(Path(tool.Ifc.get_path())):
        return False
    try:
        root = ET.parse(cache_path).getroot()
    except ET.ParseError:
        return False
    if not any(element.tag.endswith("path") for element in root.iter()):
        return False
    for svg_element in root.iter():
        attributes = {
            key.rsplit("}", 1)[-1].split(":")[-1]: value
            for key, value in svg_element.attrib.items()
        }
        global_id = attributes.get("guid")
        cached_name = attributes.get("name")
        if not global_id or cached_name is None:
            continue
        entity = tool.Ifc.get().by_guid(global_id)
        if entity and cached_name != (entity.Name or ""):
            return False
    return True


def build_view(project_root: Path, row: dict[str, str]) -> dict[str, Any]:
    room_bbox = get_space_bbox(row)
    room_element_candidates = room_elements(room_bbox)
    elements, complexity_exclusions = linework_elements(room_element_candidates)
    drawing = find_or_create_drawing(row)
    camera, svg_path = configure_drawing(project_root, row, drawing, room_bbox, elements)

    result = bpy.ops.bim.activate_drawing(
        drawing=drawing.id(), should_view_from_camera=False
    )
    if result != {"FINISHED"}:
        raise RuntimeError(f"failed to activate {drawing.Name}: {result}")
    dprops = tool.Drawing.get_document_props()
    dprops.should_use_underlay_cache = False
    dprops.should_use_linework_cache = LINEWORK_MODE == "FREESTYLE" and valid_linework_cache(
        project_root, svg_path, drawing.Name
    )
    dprops.should_use_annotation_cache = False
    result = bpy.ops.bim.create_drawing(
        print_all=False, open_viewer=False, sync=False
    )
    if result != {"FINISHED"} or not svg_path.is_file() or not svg_path.stat().st_size:
        raise RuntimeError(f"failed to create {drawing.Name}: {result}")

    width, height, clip_end = camera_dimensions(row, room_bbox)
    highlights = add_integer_highlights(svg_path, camera, elements, width, height)
    lineage = source_ifc_lineage()
    report = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "generator": "Bonsai 0.8.4 bim.create_drawing",
        "linework_mode": LINEWORK_MODE,
        "target_view": "ELEVATION_VIEW",
        "sheet_id": row["sheet_id"],
        "view_id": row["view_id"],
        "direction": row["direction"],
        "source_handle": row["source_handle"],
        "source_locator": row["source_locator"],
        "space_reference": row["space_reference"],
        "space_global_id": row["space_global_id"],
        "space_global_ids": row.get("space_global_ids", row["space_global_id"]).split(";"),
        "drawing": {
            "id": drawing.id(),
            "global_id": drawing.GlobalId,
            "name": drawing.Name,
            "object_type": drawing.ObjectType,
        },
        "camera": {
            "location_m": [float(value) for value in camera.matrix_world.translation],
            "matrix_world": [
                [float(value) for value in matrix_row]
                for matrix_row in camera.matrix_world
            ],
            "width_m": width,
            "height_m": height,
            "clip_start_m": camera.data.clip_start,
            "clip_end_m": clip_end,
            "ortho_scale_m": camera.data.ortho_scale,
            "resolution": [bpy.context.scene.render.resolution_x, bpy.context.scene.render.resolution_y],
        },
        "scale": row.get("scale", "1/50"),
        "human_scale": row.get("human_scale", "1:50"),
        "room_bbox_m": [list(room_bbox[0]), list(room_bbox[1])],
        "include_count": len(elements),
        "include_global_ids": [item[0].GlobalId for item in elements],
        "include_class_counts": dict(Counter(item[0].is_a() for item in elements)),
        "lightweight_elevation_global_ids": [
            item[0].GlobalId
            for item in elements
            if ifcopenshell.util.representation.get_representation(
                item[0], "Model", "Body", "ELEVATION_VIEW"
            )
        ],
        "demolish_wall_count": 0,
        "complexity_exclusion_count": len(complexity_exclusions),
        "complexity_exclusions": complexity_exclusions,
        "opencascade_representation_node_limit": (
            MAX_OCC_REPRESENTATION_NODES if LINEWORK_MODE == "OPENCASCADE" else None
        ),
        "integer_highlight_count": len(highlights),
        "major_integer_highlight_count": sum(item["severity"] == "major" for item in highlights),
        "maximum_integer_residual_mm": max((item["residual_mm"] for item in highlights), default=0.0),
        "integer_highlights": highlights,
        "svg": str(svg_path.relative_to(project_root)),
        "svg_sha256": sha256(svg_path),
        **lineage,
        # Legacy cache key: this intentionally follows the IFC loaded into
        # Bonsai, which may be the drawing-only derivative rather than SSOT.
        "formal_ifc_sha256_before_save": sha256(Path(tool.Ifc.get_path())),
        "pass": True,
    }
    report_path = project_root / REPORT_DIR / f"{drawing.Name}-source.json"
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return report


def run(
    sheet_ids: Iterable[str] | None = None,
    view_ids: Iterable[str] | None = None,
) -> list[dict[str, Any]]:
    if not tool.Ifc.get():
        raise RuntimeError("no IFC project is loaded")
    if LINEWORK_MODE == "FREESTYLE" and not hasattr(bpy.context.scene, "svg_export"):
        bpy.ops.preferences.addon_enable(
            module="bl_ext.blender_org.freestyle_svg_exporter"
        )
    project_root = PROJECT_ROOT
    selected = set(sheet_ids) if sheet_ids else None
    selected_views = {str(view_id).zfill(2) for view_id in view_ids} if view_ids else None
    rows = read_rows(project_root, selected, selected_views)
    reports = [build_view(project_root, row) for row in rows]
    summary_path = project_root / REPORT_DIR / "native-elevation-summary.json"
    summary_path.write_text(
        json.dumps(
            {
                "generated_at": datetime.now(timezone.utc).isoformat(),
                "view_count": len(reports),
                "sheet_ids": sorted({report["sheet_id"] for report in reports}),
                "views": reports,
                "pass": True,
            },
            ensure_ascii=False,
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    print(
        json.dumps(
            {
                "view_count": len(reports),
                "sheet_ids": sorted({report["sheet_id"] for report in reports}),
                "drawing_ids": [report["drawing"]["id"] for report in reports],
                "pass": True,
            },
            ensure_ascii=False,
        )
    )
    return reports


if __name__ == "__main__":
    run()
