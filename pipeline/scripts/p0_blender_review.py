"""Prepare one true-depth Blender scene for the accumulated A-104/A-105 review."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import bpy
import bonsai.tool as tool
import ifcopenshell.util.element
from mathutils import Vector


LABEL_COLLECTION = "P0_REVIEW_LABELS"
LEGACY_LABEL_COLLECTIONS = ("A104_REVIEW_LABELS",)
COLORS = {
    "existing_wall": (0.48, 0.52, 0.57, 1.0),
    "new_wall": (0.58, 0.26, 0.72, 1.0),
    "door": (0.00, 0.72, 0.78, 1.0),
    "window": (0.10, 0.36, 0.95, 1.0),
    "a104_review": (1.00, 0.03, 0.02, 1.0),
    "wet_floor": (0.00, 0.38, 0.88, 1.0),
    "material_pending": (1.00, 0.34, 0.00, 1.0),
    "a105_slab": (0.72, 0.15, 0.82, 1.0),
    "furniture": (1.00, 0.72, 0.00, 1.0),
    "grid": (0.12, 0.78, 0.20, 1.0),
    "floor_context": (0.66, 0.68, 0.72, 1.0),
    "sanitary": (0.88, 0.90, 0.93, 1.0),
    "drainage": (0.12, 0.65, 0.78, 1.0),
    "slope_arrow": (1.00, 1.00, 1.00, 1.0),
}


def entity_for(obj: bpy.types.Object):
    try:
        return tool.Ifc.get_entity(obj)
    except Exception:
        return None


def remove_label_collection(name: str) -> None:
    collection = bpy.data.collections.get(name)
    if collection is None:
        return
    for obj in list(collection.objects):
        bpy.data.objects.remove(obj, do_unlink=True)
    bpy.data.collections.remove(collection)


def remove_labels() -> None:
    remove_label_collection(LABEL_COLLECTION)
    for name in LEGACY_LABEL_COLLECTIONS:
        remove_label_collection(name)


def world_bounds(obj: bpy.types.Object) -> tuple[Vector, Vector]:
    corners = [obj.matrix_world @ Vector(corner) for corner in obj.bound_box]
    minimum = Vector(tuple(min(point[axis] for point in corners) for axis in range(3)))
    maximum = Vector(tuple(max(point[axis] for point in corners) for axis in range(3)))
    return minimum, maximum


def project_paths() -> tuple[Path, Path, Path, Path]:
    root = Path(str(tool.Ifc.get_path())).resolve().parent
    return (
        root / "pipeline/decisions/a104-door-window-review.csv",
        root / "pipeline/decisions/a105-floor-review.csv",
        root / "build/a105/a105-report.json",
        root / "pipeline/decisions/p0-review.csv",
    )


def load_review_data() -> tuple[dict[str, dict[str, str]], dict[str, dict[str, str]], dict, set[str]]:
    a104_path, a105_path, report_path, decisions_path = project_paths()
    with a104_path.open(encoding="utf-8-sig", newline="") as handle:
        a104 = {row["global_id"]: row for row in csv.DictReader(handle)}
    with a105_path.open(encoding="utf-8-sig", newline="") as handle:
        a105 = {row["global_id"]: row for row in csv.DictReader(handle)}
    report = json.loads(report_path.read_text(encoding="utf-8"))
    confirmed_a104_ids: set[str] = set()
    with decisions_path.open(encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle):
            if row["decision_id"].startswith("A104-") and row["status"] == "confirmed":
                confirmed_a104_ids.update(value.strip() for value in row["object_guid"].split(";") if value.strip())
    if len(a104) != 19 or len(a105) != 21:
        raise RuntimeError(f"review register count drift: A104={len(a104)}, A105={len(a105)}")
    return a104, a105, report, confirmed_a104_ids


def switch_to_model_body(global_ids: set[str]) -> int:
    switched = 0
    model = tool.Ifc.get()
    for global_id in sorted(global_ids):
        entity = model.by_guid(global_id)
        obj = tool.Ifc.get_object(entity) if entity else None
        if obj is None or entity.Representation is None:
            raise RuntimeError(f"cannot find Blender object/representation for {global_id}")
        body = next(
            (
                representation
                for representation in entity.Representation.Representations
                if representation.RepresentationIdentifier == "Body"
                and representation.ContextOfItems.ContextType == "Model"
            ),
            None,
        )
        if body is None:
            raise RuntimeError(f"no Model/Body representation for {global_id}")
        active = tool.Geometry.get_active_representation(obj)
        if active is None or active.id() != body.id():
            bpy.context.view_layer.objects.active = obj
            obj.select_set(True)
            result = bpy.ops.bim.switch_representation(
                obj=obj.name,
                ifc_definition_id=body.id(),
                disable_opening_subtractions=False,
            )
            obj.select_set(False)
            if "FINISHED" not in result:
                raise RuntimeError(f"failed to show Model/Body for {global_id}: {result}")
            switched += 1
    return switched


def add_label(
    collection: bpy.types.Collection,
    obj: bpy.types.Object,
    body: str,
    color: tuple[float, float, float, float],
    size: float = 0.13,
) -> None:
    minimum, maximum = world_bounds(obj)
    curve = bpy.data.curves.new(f"P0_LABEL_{body}", "FONT")
    curve.body = body
    curve.align_x = "CENTER"
    curve.align_y = "CENTER"
    curve.size = size
    curve.extrude = 0.003
    label = bpy.data.objects.new(f"P0_LABEL_{body}", curve)
    collection.objects.link(label)
    label.location = (
        (minimum.x + maximum.x) / 2.0,
        (minimum.y + maximum.y) / 2.0,
        maximum.z + 0.045,
    )
    label.color = color
    label.show_in_front = False


def add_bbox_outline(collection: bpy.types.Collection, obj: bpy.types.Object) -> None:
    minimum, maximum = world_bounds(obj)
    curve = bpy.data.curves.new(f"P0_OUTLINE_{obj.name}", "CURVE")
    curve.dimensions = "3D"
    curve.bevel_depth = 0.012
    curve.bevel_resolution = 1
    spline = curve.splines.new("POLY")
    spline.points.add(3)
    z = maximum.z + 0.006
    points = [
        (minimum.x, minimum.y, z, 1.0),
        (maximum.x, minimum.y, z, 1.0),
        (maximum.x, maximum.y, z, 1.0),
        (minimum.x, maximum.y, z, 1.0),
    ]
    for point, value in zip(spline.points, points):
        point.co = value
    spline.use_cyclic_u = True
    outline = bpy.data.objects.new(f"P0_OUTLINE_{obj.name}", curve)
    collection.objects.link(outline)
    outline.color = COLORS["material_pending"]
    outline.show_in_front = False


def add_slope_arrow(collection: bpy.types.Collection, obj: bpy.types.Object, row: dict[str, str]) -> None:
    plane = json.loads(row["top_plane"])
    direction = Vector((-float(plane["a_dz_dx"]), -float(plane["b_dz_dy"]), 0.0))
    if direction.length <= 1e-12:
        return
    direction.normalize()
    minimum, maximum = world_bounds(obj)
    centre = (minimum + maximum) / 2.0
    centre.z = (
        float(plane["a_dz_dx"]) * centre.x * 1000.0
        + float(plane["b_dz_dy"]) * centre.y * 1000.0
        + float(plane["c_mm"])
    ) / 1000.0 + 0.012
    shortest_side = min(maximum.x - minimum.x, maximum.y - minimum.y)
    length = min(0.28, max(0.12, shortest_side * 0.42))
    head_length = min(0.065, length * 0.28)
    perpendicular = Vector((-direction.y, direction.x, 0.0))
    start = centre - direction * length / 2.0
    end = centre + direction * length / 2.0
    head_base = end - direction * head_length
    points = (
        (start, end),
        (head_base + perpendicular * head_length * 0.55, end),
        (head_base - perpendicular * head_length * 0.55, end),
    )
    curve = bpy.data.curves.new(f"P0_SLOPE_{row['candidate_id']}", "CURVE")
    curve.dimensions = "3D"
    curve.bevel_depth = 0.008
    curve.bevel_resolution = 1
    for first, second in points:
        spline = curve.splines.new("POLY")
        spline.points.add(1)
        spline.points[0].co = (*first, 1.0)
        spline.points[1].co = (*second, 1.0)
    arrow = bpy.data.objects.new(f"P0_SLOPE_{row['candidate_id']}", curve)
    collection.objects.link(arrow)
    arrow.color = COLORS["slope_arrow"]
    arrow.show_in_front = False


def configure_objects(
    a104: dict[str, dict[str, str]],
    a105: dict[str, dict[str, str]],
    slab_global_id: str,
    confirmed_a104_ids: set[str],
) -> tuple[
    dict[str, int],
    list[tuple[bpy.types.Object, str, tuple[float, float, float, float], float]],
    list[bpy.types.Object],
    list[tuple[bpy.types.Object, dict[str, str]]],
]:
    counts: dict[str, int] = {}
    labels: list[tuple[bpy.types.Object, str, tuple[float, float, float, float], float]] = []
    outlines: list[bpy.types.Object] = []
    arrows: list[tuple[bpy.types.Object, dict[str, str]]] = []
    for obj in bpy.context.scene.objects:
        obj.show_in_front = False
        entity = entity_for(obj)
        if entity is None:
            continue
        visible = False
        key = ""
        if entity.is_a("IfcWall"):
            common = ifcopenshell.util.element.get_psets(entity).get("Pset_WallCommon", {})
            if str(common.get("Status") or "").upper() != "DEMOLISH":
                visible = True
                key = "new_wall" if str(common.get("Status") or "").upper() == "NEW" else "existing_wall"
        elif entity.is_a("IfcDoor") or entity.is_a("IfcWindow"):
            visible = True
            row = a104.get(entity.GlobalId)
            if row is None:
                raise RuntimeError(f"door/window missing from A-104 register: {entity.GlobalId}")
            if row["review_required"] == "yes" and entity.GlobalId not in confirmed_a104_ids:
                key = "a104_review"
                labels.append((obj, f"{row['candidate_id']} {row['review_group']}", COLORS[key], 0.15))
            else:
                key = "door" if entity.is_a("IfcDoor") else "window"
                if entity.is_a("IfcDoor"):
                    labels.append((obj, row["candidate_id"], COLORS[key], 0.11))
        elif entity.is_a("IfcFurniture") or entity.is_a("IfcFurnishingElement"):
            visible = True
            key = "furniture"
        elif entity.is_a("IfcGridAxis"):
            visible = True
            key = "grid"
        elif entity.GlobalId in a105:
            visible = True
            row = a105[entity.GlobalId]
            key = "wet_floor" if row["kind"] == "sloped_wet_tile" else "material_pending"
            if row["kind"] == "sloped_wet_tile":
                plane = json.loads(row["top_plane"])
                text = f"{row['candidate_id']} {float(plane['slope_percent']):.3f}% {plane['downhill_direction']}"
                arrows.append((obj, row))
            else:
                text = f"{row['candidate_id']} {row.get('confirmed_material') or 'MATERIAL?'} / 50 mm"
                outlines.append(obj)
            labels.append((obj, text, COLORS[key], 0.105))
        elif entity.GlobalId == slab_global_id:
            visible = True
            key = "a105_slab"
            labels.append((obj, "A105-R04 AIRCRETE SLAB?", COLORS[key], 0.14))
        elif entity.is_a("IfcSlab") or entity.is_a("IfcCovering"):
            _minimum, maximum = world_bounds(obj)
            visible = maximum.z <= 0.15
            key = "floor_context"
        elif entity.is_a("IfcSanitaryTerminal"):
            visible = True
            key = "sanitary"
        elif entity.is_a("IfcFlowSegment"):
            visible = True
            key = "drainage"
        obj.hide_set(not visible)
        if visible:
            obj.color = COLORS[key]
            obj.display_type = "SOLID"
            # Hide the zero-thickness source plane after measuring it; a
            # true-depth orange solid curve below shows its extents without
            # obscuring the real sloped floor solids beneath Z=0.
            if key == "material_pending":
                obj.hide_set(True)
            counts[key] = counts.get(key, 0) + 1
    return counts, labels, outlines, arrows


def configure_viewport() -> None:
    screen = bpy.context.screen
    if screen is None:
        return
    for area in screen.areas:
        if area.type != "VIEW_3D":
            continue
        space = area.spaces.active
        space.shading.type = "SOLID"
        space.shading.light = "STUDIO"
        space.shading.color_type = "OBJECT"
        space.shading.background_type = "VIEWPORT"
        space.shading.background_color = (0.055, 0.055, 0.055)
        space.shading.show_xray = False
        space.overlay.show_overlays = True
        space.overlay.show_outline_selected = True
        space.region_3d.view_perspective = "PERSP"


def main() -> None:
    model = tool.Ifc.get()
    if model is None:
        raise RuntimeError("no IFC is loaded")
    a104, a105, report, confirmed_a104_ids = load_review_data()
    slab_global_id = report["delegated_slab"]["global_id"]
    remove_labels()
    if bpy.context.object is not None and bpy.context.object.mode != "OBJECT":
        bpy.ops.object.mode_set(mode="OBJECT")
    bpy.ops.object.select_all(action="DESELECT")
    switched = switch_to_model_body(set(a104) | set(a105) | {slab_global_id})
    counts, label_rows, outline_rows, arrow_rows = configure_objects(a104, a105, slab_global_id, confirmed_a104_ids)
    collection = bpy.data.collections.new(LABEL_COLLECTION)
    bpy.context.scene.collection.children.link(collection)
    for obj in outline_rows:
        add_bbox_outline(collection, obj)
    for obj, row in arrow_rows:
        add_slope_arrow(collection, obj, row)
    for args in label_rows:
        add_label(collection, *args)
    bpy.ops.object.select_all(action="DESELECT")
    configure_viewport()
    print(
        {
            "ifc_path": str(tool.Ifc.get_path()),
            "schema": model.schema,
            "a104_review_objects": sum(
                row["review_required"] == "yes" and global_id not in confirmed_a104_ids
                for global_id, row in a104.items()
            ),
            "a105_floor_objects": len(a105),
            "a105_slab_global_id": slab_global_id,
            "labels": len(label_rows),
            "material_pending_outlines": len(outline_rows),
            "slope_arrows": len(arrow_rows),
            "switched_to_model_body": switched,
            "context": counts,
            "selected": len(bpy.context.selected_objects),
            "show_in_front": sum(obj.show_in_front for obj in bpy.context.scene.objects),
            "solid_view": True,
            "xray": False,
        }
    )


main()
