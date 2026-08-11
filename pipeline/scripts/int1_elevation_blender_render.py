"""Render the 36 INT1 CAD-indexed elevations from the open Bonsai IFC.

Run this file inside the already-open Blender/Bonsai session.  It is deliberately
render-only: it never saves the .blend or writes the IFC.  Every scene and viewport
property changed for the batch is restored before the manifest is written.
"""

from __future__ import annotations

import csv
import hashlib
import json
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable


VIEW_REGISTER_REL = Path("pipeline/decisions/int1-elevation-view-register.csv")
A102_WALL_REVIEW_REL = Path("pipeline/decisions/a102-wall-review.csv")
A102_DEMOLITION_REVIEW_REL = Path("pipeline/decisions/a102-demolition-review.csv")
OUTPUT_DIR_REL = Path("build/int1/elevations/raw")
MANIFEST_REL = Path("build/int1/elevation-render-manifest.json")

EXPECTED_VIEW_COUNT = 36
RENDER_WIDTH_PX = 1800
RENDER_HEIGHT_PX = 1200
SPACE_MARGIN_M = 0.35
VERTICAL_MARGIN_M = 0.25
DEPTH_MARGIN_M = 0.50
MIN_ORTHO_HEIGHT_M = 2.50
DEMOLITION_STATUSES = {"DEMOLISH", "DEMOLISHED"}
EXCLUDED_IFC_CLASSES = {"IfcSpace", "IfcGrid", "IfcAnnotation", "IfcOpeningElement"}
DIRECTION_VECTORS = {
    "+X": (1.0, 0.0),
    "-X": (-1.0, 0.0),
    "+Y": (0.0, 1.0),
    "-Y": (0.0, -1.0),
}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def first_value(row: dict[str, str], *names: str) -> str:
    for name in names:
        value = str(row.get(name, "")).strip()
        if value:
            return value
    return ""


def normalize_direction(raw: str) -> str:
    token = raw.strip().upper().replace(" ", "")
    aliases = {
        "E": "+X",
        "EAST": "+X",
        "东": "+X",
        "W": "-X",
        "WEST": "-X",
        "西": "-X",
        "N": "+Y",
        "NORTH": "+Y",
        "北": "+Y",
        "S": "-Y",
        "SOUTH": "-Y",
        "南": "-Y",
    }
    token = aliases.get(token, token)
    if token not in DIRECTION_VECTORS:
        raise ValueError(f"unsupported elevation direction: {raw!r}")
    return token


def _view_sort_key(view_id: str) -> tuple[int, str]:
    digits = "".join(character for character in view_id if character.isdigit())
    return (int(digits) if digits else 10**9, view_id)


def load_view_rows(path: Path) -> list[dict[str, Any]]:
    if not path.exists():
        raise FileNotFoundError(f"missing elevation view register: {path}")
    with path.open(newline="", encoding="utf-8-sig") as handle:
        source_rows = list(csv.DictReader(handle))

    rows: list[dict[str, Any]] = []
    seen: set[str] = set()
    for source in source_rows:
        view_id = first_value(source, "view_id", "view_number", "cad_view_id", "index")
        if not view_id:
            raise ValueError("each elevation register row requires view_id")
        normalized_view_id = view_id.zfill(2) if view_id.isdigit() else view_id
        if normalized_view_id in seen:
            raise ValueError(f"duplicate elevation view_id: {normalized_view_id}")
        seen.add(normalized_view_id)
        anchor_x = first_value(source, "anchor_x_mm", "ifc_anchor_x_mm", "ifc_x_mm")
        anchor_y = first_value(source, "anchor_y_mm", "ifc_anchor_y_mm", "ifc_y_mm")
        if not anchor_x or not anchor_y:
            raise ValueError(f"view {view_id} requires IFC anchor_x_mm and anchor_y_mm")
        direction = normalize_direction(first_value(source, "direction", "view_direction", "ifc_direction"))
        sheet_id = first_value(source, "sheet_id", "drawing_id", "el_sheet", "sheet")
        anchor_id = first_value(source, "anchor_id", "anchor", "anchor_cluster")
        if not sheet_id:
            raise ValueError(f"view {view_id} requires sheet_id")
        if not anchor_id:
            raise ValueError(f"view {view_id} requires anchor_id")
        rows.append(
            {
                "view_id": normalized_view_id,
                "sheet_id": sheet_id,
                "anchor_id": anchor_id,
                "anchor_x_mm": float(anchor_x),
                "anchor_y_mm": float(anchor_y),
                "direction": direction,
                "space_global_id": first_value(source, "space_global_id", "space_guid", "ifc_space_global_id"),
                "space_name": first_value(source, "space_name", "room_name"),
                "space_reference": first_value(source, "space_reference", "space_ref", "room_id"),
                "review_status": first_value(source, "review_status", "status"),
                "source": dict(source),
            }
        )
    rows.sort(key=lambda row: _view_sort_key(row["view_id"]))
    if len(rows) != EXPECTED_VIEW_COUNT:
        raise ValueError(f"expected {EXPECTED_VIEW_COUNT} unique elevation views, found {len(rows)}")
    return rows


def resolve_demolition_register(root: Path) -> Path:
    preferred = root / A102_WALL_REVIEW_REL
    if preferred.exists():
        return preferred
    fallback = root / A102_DEMOLITION_REVIEW_REL
    if fallback.exists():
        return fallback
    raise FileNotFoundError(
        f"missing demolition register: expected {preferred} or current fallback {fallback}"
    )


def load_demolition_ids(path: Path) -> set[str]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        rows = list(csv.DictReader(handle))
    ids = {
        first_value(row, "candidate_global_id", "global_id", "ifc_global_id")
        for row in rows
        if first_value(row, "ifc_status_candidate", "status", "ifc_status").upper()
        in DEMOLITION_STATUSES
    }
    ids.discard("")
    if not ids:
        raise ValueError(f"no DEMOLISH wall IDs found in {path}")
    return ids


def project_root_from_ifc_path(ifc_path: str) -> Path:
    path = Path(ifc_path).expanduser().resolve()
    if not path.is_file() or path.suffix.lower() != ".ifc":
        raise RuntimeError(f"Bonsai does not report a valid open IFC path: {path}")
    return path.parent


@dataclass(frozen=True)
class Bounds:
    minimum: tuple[float, float, float]
    maximum: tuple[float, float, float]
    corners: tuple[tuple[float, float, float], ...]


def object_world_bounds(obj: Any) -> Bounds | None:
    if not getattr(obj, "bound_box", None):
        return None
    corners = tuple(tuple(obj.matrix_world @ _vector(corner)) for corner in obj.bound_box)
    if not corners:
        return None
    return Bounds(
        minimum=tuple(min(point[axis] for point in corners) for axis in range(3)),
        maximum=tuple(max(point[axis] for point in corners) for axis in range(3)),
        corners=corners,
    )


def _vector(values: Iterable[float]) -> Any:
    from mathutils import Vector

    return Vector(values)


def bbox_json_mm(bounds: Bounds) -> dict[str, list[float]]:
    return {
        "min": [round(value * 1000.0, 6) for value in bounds.minimum],
        "max": [round(value * 1000.0, 6) for value in bounds.maximum],
    }


def bbox_min_mm(bounds: Bounds) -> list[float]:
    return [round(value * 1000.0, 6) for value in bounds.minimum]


def bbox_max_mm(bounds: Bounds) -> list[float]:
    return [round(value * 1000.0, 6) for value in bounds.maximum]


def entity_container_name(entity: Any) -> str:
    import ifcopenshell.util.element

    container = ifcopenshell.util.element.get_container(entity)
    return str(getattr(container, "Name", "") or "").strip()


def is_typed_ceiling(entity: Any) -> bool:
    """Return true only for IFC coverings explicitly typed as ceilings."""
    import ifcopenshell.util.element

    if not entity.is_a("IfcCovering"):
        return False
    covering_type = ifcopenshell.util.element.get_type(entity)
    occurrence_predefined = str(getattr(entity, "PredefinedType", "") or "").upper()
    type_predefined = str(getattr(covering_type, "PredefinedType", "") or "").upper()
    return "CEILING" in {occurrence_predefined, type_predefined}


def is_renderable_entity(entity: Any, demolition_ids: set[str]) -> bool:
    if entity is None or not entity.is_a("IfcProduct"):
        return False
    if any(entity.is_a(name) for name in EXCLUDED_IFC_CLASSES):
        return False
    if str(getattr(entity, "GlobalId", "")) in demolition_ids:
        return False
    if any(
        entity.is_a(name)
        for name in ("IfcWall", "IfcDoor", "IfcWindow", "IfcFurniture", "IfcSanitaryTerminal")
    ):
        return True
    if entity.is_a("IfcDistributionElement") or entity.is_a("IfcElementAssembly"):
        return True
    if is_typed_ceiling(entity):
        return True
    if entity.is_a("IfcBuildingElementProxy"):
        name = str(getattr(entity, "Name", "") or "").lower()
        return entity_container_name(entity).upper() == "DCL" and "ceiling" in name
    return False


def bounds_intersect_xy(first: Bounds, second: Bounds, margin: float) -> bool:
    return not (
        first.maximum[0] < second.minimum[0] - margin
        or first.minimum[0] > second.maximum[0] + margin
        or first.maximum[1] < second.minimum[1] - margin
        or first.minimum[1] > second.maximum[1] + margin
    )


def point_in_bounds_xy(x_m: float, y_m: float, bounds: Bounds, tolerance_m: float = 0.01) -> bool:
    return (
        bounds.minimum[0] - tolerance_m <= x_m <= bounds.maximum[0] + tolerance_m
        and bounds.minimum[1] - tolerance_m <= y_m <= bounds.maximum[1] + tolerance_m
    )


def resolve_space(model: Any, tool: Any, row: dict[str, Any]) -> tuple[Any, Any, Bounds]:
    requested_guid = row["space_global_id"]
    if requested_guid:
        entity = model.by_guid(requested_guid)
        if entity is None or not entity.is_a("IfcSpace"):
            raise RuntimeError(f"view {row['view_id']} references missing IfcSpace {requested_guid}")
        obj = tool.Ifc.get_object(entity)
        bounds = object_world_bounds(obj) if obj is not None else None
        if bounds is None:
            raise RuntimeError(f"view {row['view_id']} has no loaded geometry for IfcSpace {requested_guid}")
        return entity, obj, bounds

    x_m = row["anchor_x_mm"] / 1000.0
    y_m = row["anchor_y_mm"] / 1000.0
    candidates: list[tuple[float, Any, Any, Bounds]] = []
    for entity in model.by_type("IfcSpace"):
        obj = tool.Ifc.get_object(entity)
        bounds = object_world_bounds(obj) if obj is not None else None
        if bounds is None or not point_in_bounds_xy(x_m, y_m, bounds):
            continue
        area = (bounds.maximum[0] - bounds.minimum[0]) * (bounds.maximum[1] - bounds.minimum[1])
        candidates.append((area, entity, obj, bounds))
    if not candidates:
        raise RuntimeError(
            f"view {row['view_id']} anchor ({row['anchor_x_mm']}, {row['anchor_y_mm']}) mm "
            "does not match loaded IfcSpace geometry"
        )
    _, entity, obj, bounds = min(candidates, key=lambda item: item[0])
    return entity, obj, bounds


def projection_extents(
    bounds: Bounds,
    anchor_m: tuple[float, float],
    forward_xy: tuple[float, float],
) -> dict[str, float]:
    lateral_xy = (-forward_xy[1], forward_xy[0])
    lateral = []
    depth = []
    for x, y, _z in bounds.corners:
        offset_x = x - anchor_m[0]
        offset_y = y - anchor_m[1]
        lateral.append(offset_x * lateral_xy[0] + offset_y * lateral_xy[1])
        depth.append(offset_x * forward_xy[0] + offset_y * forward_xy[1])
    return {
        "lateral_min": min(lateral),
        "lateral_max": max(lateral),
        "depth_min": min(depth),
        "depth_max": max(depth),
        "z_min": bounds.minimum[2],
        "z_max": bounds.maximum[2],
    }


def projected_bbox_json(
    bounds: Bounds,
    anchor_m: tuple[float, float],
    forward_xy: tuple[float, float],
) -> dict[str, list[float]]:
    extents = projection_extents(bounds, anchor_m, forward_xy)
    return {
        "u_min_mm": round(extents["lateral_min"] * 1000.0, 6),
        "u_max_mm": round(extents["lateral_max"] * 1000.0, 6),
        "z_min_mm": round(extents["z_min"] * 1000.0, 6),
        "z_max_mm": round(extents["z_max"] * 1000.0, 6),
        "depth_min_mm": round(extents["depth_min"] * 1000.0, 6),
        "depth_max_mm": round(extents["depth_max"] * 1000.0, 6),
    }


def capture_viewport_states(bpy: Any) -> list[tuple[Any, dict[str, Any], dict[str, Any]]]:
    shading_names = (
        "type",
        "light",
        "color_type",
        "single_color",
        "show_xray",
        "show_shadows",
        "show_cavity",
        "cavity_type",
        "show_specular_highlight",
        "show_outline",
        "background_type",
        "background_color",
    )
    overlay_names = ("show_wireframes", "show_relationship_lines")
    states = []
    for screen in bpy.data.screens:
        for area in screen.areas:
            if area.type != "VIEW_3D":
                continue
            space = area.spaces.active
            states.append(
                (
                    space,
                    {name: getattr(space.shading, name) for name in shading_names if hasattr(space.shading, name)},
                    {name: getattr(space.overlay, name) for name in overlay_names if hasattr(space.overlay, name)},
                )
            )
    return states


def capture_scene_state(bpy: Any) -> dict[str, Any]:
    scene = bpy.context.scene
    image = scene.render.image_settings
    shading = scene.display.shading
    shading_names = (
        "light",
        "color_type",
        "single_color",
        "show_xray",
        "show_shadows",
        "show_cavity",
        "cavity_type",
        "show_specular_highlight",
        "show_outline",
        "outline_color",
        "background_type",
        "background_color",
    )
    return {
        "camera": scene.camera,
        "render": {
            "engine": scene.render.engine,
            "filepath": scene.render.filepath,
            "resolution_x": scene.render.resolution_x,
            "resolution_y": scene.render.resolution_y,
            "resolution_percentage": scene.render.resolution_percentage,
            "film_transparent": scene.render.film_transparent,
        },
        "image": {
            "file_format": image.file_format,
            "color_mode": image.color_mode,
            "color_depth": image.color_depth,
            "compression": image.compression,
        },
        "world_color": tuple(scene.world.color) if scene.world is not None else None,
        "display_shading": {
            name: getattr(shading, name) for name in shading_names if hasattr(shading, name)
        },
        "viewport_states": capture_viewport_states(bpy),
        "objects": {
            obj.name_full: {
                "object": obj,
                "hide_render": obj.hide_render,
                "hide_viewport": obj.hide_viewport,
                "hide_set": obj.hide_get(),
                "show_in_front": obj.show_in_front,
            }
            for obj in bpy.context.scene.objects
        },
        "selected": list(bpy.context.selected_objects),
        "active": bpy.context.view_layer.objects.active,
    }


def assign_if_present(target: Any, name: str, value: Any) -> None:
    if hasattr(target, name):
        setattr(target, name, value)


def restore_scene_state(bpy: Any, state: dict[str, Any], temporary_camera: Any | None) -> dict[str, bool]:
    scene = bpy.context.scene
    for name, value in state["render"].items():
        setattr(scene.render, name, value)
    for name, value in state["image"].items():
        setattr(scene.render.image_settings, name, value)
    for name, value in state["display_shading"].items():
        assign_if_present(scene.display.shading, name, value)
    if scene.world is not None and state["world_color"] is not None:
        scene.world.color = state["world_color"]
    for space, shading_state, overlay_state in state["viewport_states"]:
        for name, value in shading_state.items():
            assign_if_present(space.shading, name, value)
        for name, value in overlay_state.items():
            assign_if_present(space.overlay, name, value)
    for item in state["objects"].values():
        obj = item["object"]
        if obj.name not in bpy.data.objects:
            continue
        obj.hide_render = item["hide_render"]
        obj.hide_viewport = item["hide_viewport"]
        obj.hide_set(item["hide_set"])
        obj.show_in_front = item["show_in_front"]
    scene.camera = state["camera"]
    for obj in list(bpy.context.selected_objects):
        obj.select_set(False)
    for obj in state["selected"]:
        if obj.name in bpy.data.objects:
            obj.select_set(True)
    if state["active"] is not None and state["active"].name in bpy.data.objects:
        bpy.context.view_layer.objects.active = state["active"]
    if temporary_camera is not None and temporary_camera.name in bpy.data.objects:
        camera_data = temporary_camera.data
        bpy.data.objects.remove(temporary_camera, do_unlink=True)
        if camera_data is not None and camera_data.users == 0:
            bpy.data.cameras.remove(camera_data)
    return {
        "render_state": scene.render.engine == state["render"]["engine"]
        and scene.render.filepath == state["render"]["filepath"],
        "camera_state": scene.camera == state["camera"],
        "object_visibility_state": all(
            item["object"].name not in bpy.data.objects
            or (
                item["object"].hide_render == item["hide_render"]
                and item["object"].hide_viewport == item["hide_viewport"]
                and item["object"].hide_get() == item["hide_set"]
                and item["object"].show_in_front == item["show_in_front"]
            )
            for item in state["objects"].values()
        ),
        "viewport_shading_state": all(
            all(getattr(space.shading, name) == value for name, value in shading_state.items())
            and all(getattr(space.overlay, name) == value for name, value in overlay_state.items())
            for space, shading_state, overlay_state in state["viewport_states"]
        ),
    }


def configure_cad_like_render(bpy: Any) -> str:
    scene = bpy.context.scene
    selected_engine = ""
    errors = []
    for engine in ("BLENDER_WORKBENCH_NEXT", "BLENDER_WORKBENCH"):
        try:
            scene.render.engine = engine
            selected_engine = engine
            break
        except (TypeError, ValueError) as error:
            errors.append(f"{engine}: {error}")
    if not selected_engine:
        raise RuntimeError(f"no Blender 4.x Workbench render engine available: {'; '.join(errors)}")

    scene.render.resolution_x = RENDER_WIDTH_PX
    scene.render.resolution_y = RENDER_HEIGHT_PX
    scene.render.resolution_percentage = 100
    scene.render.film_transparent = False
    scene.render.image_settings.file_format = "PNG"
    scene.render.image_settings.color_mode = "RGBA"
    scene.render.image_settings.color_depth = "8"
    scene.render.image_settings.compression = 15
    if scene.world is not None:
        scene.world.color = (1.0, 1.0, 1.0)

    shading = scene.display.shading
    assign_if_present(shading, "light", "STUDIO")
    assign_if_present(shading, "color_type", "SINGLE")
    assign_if_present(shading, "single_color", (0.78, 0.80, 0.82))
    assign_if_present(shading, "show_xray", False)
    assign_if_present(shading, "show_shadows", True)
    assign_if_present(shading, "show_cavity", True)
    assign_if_present(shading, "cavity_type", "WORLD")
    assign_if_present(shading, "show_specular_highlight", False)
    assign_if_present(shading, "show_outline", True)
    assign_if_present(shading, "outline_color", (0.01, 0.01, 0.01))
    assign_if_present(shading, "background_type", "WORLD")

    for screen in bpy.data.screens:
        for area in screen.areas:
            if area.type != "VIEW_3D":
                continue
            space = area.spaces.active
            space.shading.type = "SOLID"
            space.shading.color_type = "SINGLE"
            assign_if_present(space.shading, "single_color", (0.78, 0.80, 0.82))
            space.shading.show_xray = False
            assign_if_present(space.shading, "show_outline", True)
            assign_if_present(space.shading, "outline_color", (0.01, 0.01, 0.01))
            space.overlay.show_wireframes = False
            space.overlay.show_relationship_lines = False
    for obj in bpy.context.scene.objects:
        obj.show_in_front = False
    return selected_engine


def create_temporary_camera(bpy: Any) -> Any:
    camera_data = bpy.data.cameras.new("INT1_ELEVATION_RENDER_CAMERA_DATA")
    camera_data.type = "ORTHO"
    camera = bpy.data.objects.new("INT1_ELEVATION_RENDER_CAMERA", camera_data)
    bpy.context.scene.collection.objects.link(camera)
    bpy.context.scene.camera = camera
    return camera


def collect_render_objects(
    bpy: Any,
    model: Any,
    tool: Any,
    space_bounds: Bounds,
    demolition_ids: set[str],
) -> list[dict[str, Any]]:
    records = []
    seen_objects: set[Any] = set()
    for entity in model.by_type("IfcProduct"):
        if not is_renderable_entity(entity, demolition_ids):
            continue
        obj = tool.Ifc.get_object(entity)
        if obj is None or obj in seen_objects:
            continue
        bounds = object_world_bounds(obj)
        if bounds is None or not bounds_intersect_xy(bounds, space_bounds, SPACE_MARGIN_M):
            continue
        seen_objects.add(obj)
        records.append({"entity": entity, "object": obj, "bounds": bounds})

    return records


def create_sanitized_render_copies(
    bpy: Any, render_scene: Any, records: list[dict[str, Any]]
) -> tuple[Any, list[Any]]:
    """Create plain meshes without Bonsai IDProperties for Workbench rendering.

    Blender 4.5 can segfault in IDP_CopyProperty_ex when Workbench evaluates
    Bonsai meshes carrying BIMMeshProperties.  Rebuilding only vertex/face data
    keeps the rendered geometry and world transform while isolating the renderer
    from IFC/Bonsai custom properties.  Originals are never edited.
    """
    collection = bpy.data.collections.new("INT1_SANITIZED_RENDER")
    render_scene.collection.children.link(collection)
    copies = []
    for index, record in enumerate(records):
        source = record["object"]
        if source.type != "MESH" or source.data is None:
            continue
        vertices = [tuple(vertex.co) for vertex in source.data.vertices]
        faces = [tuple(polygon.vertices) for polygon in source.data.polygons]
        mesh = bpy.data.meshes.new(f"INT1_SANITIZED_MESH_{index:03d}")
        mesh.from_pydata(vertices, [], faces)
        mesh.update()
        copy = bpy.data.objects.new(f"INT1_SANITIZED_OBJECT_{index:03d}", mesh)
        copy.matrix_world = source.matrix_world.copy()
        copy.show_in_front = False
        collection.objects.link(copy)
        copies.append(copy)
    return collection, copies


def create_isolated_render_scene(bpy: Any, source_scene: Any, camera: Any) -> Any:
    """Build a scene whose dependency graph contains no Bonsai source meshes."""
    render_scene = bpy.data.scenes.new("INT1_ISOLATED_RENDER_SCENE")
    render_scene.render.engine = source_scene.render.engine
    render_scene.render.resolution_x = RENDER_WIDTH_PX
    render_scene.render.resolution_y = RENDER_HEIGHT_PX
    render_scene.render.resolution_percentage = 100
    render_scene.render.film_transparent = False
    render_scene.render.image_settings.file_format = "PNG"
    render_scene.render.image_settings.color_mode = "RGBA"
    render_scene.render.image_settings.color_depth = "8"
    render_scene.render.image_settings.compression = 15
    shading = render_scene.display.shading
    assign_if_present(shading, "light", "STUDIO")
    assign_if_present(shading, "color_type", "SINGLE")
    assign_if_present(shading, "single_color", (0.78, 0.80, 0.82))
    assign_if_present(shading, "show_xray", False)
    assign_if_present(shading, "show_shadows", True)
    assign_if_present(shading, "show_cavity", True)
    assign_if_present(shading, "cavity_type", "WORLD")
    assign_if_present(shading, "show_specular_highlight", False)
    assign_if_present(shading, "show_outline", True)
    assign_if_present(shading, "outline_color", (0.01, 0.01, 0.01))
    assign_if_present(shading, "background_type", "WORLD")
    render_scene.collection.objects.link(camera)
    render_scene.camera = camera
    return render_scene


def remove_sanitized_render_copies(bpy: Any, collection: Any | None, copies: list[Any]) -> None:
    for copy in copies:
        mesh = copy.data
        bpy.data.objects.remove(copy, do_unlink=True)
        if mesh is not None and mesh.users == 0:
            bpy.data.meshes.remove(mesh)
    if collection is not None and collection.name in bpy.data.collections:
        bpy.data.collections.remove(collection)


def camera_for_view(
    camera: Any,
    row: dict[str, Any],
    space_bounds: Bounds,
    object_records: list[dict[str, Any]],
) -> tuple[dict[str, Any], dict[str, float]]:
    forward_xy = DIRECTION_VECTORS[row["direction"]]
    lateral_xy = (-forward_xy[1], forward_xy[0])
    anchor_m = (row["anchor_x_mm"] / 1000.0, row["anchor_y_mm"] / 1000.0)
    space_extents = projection_extents(space_bounds, anchor_m, forward_xy)

    lateral_min = space_extents["lateral_min"] - SPACE_MARGIN_M
    lateral_max = space_extents["lateral_max"] + SPACE_MARGIN_M
    depth_min = space_extents["depth_min"] - DEPTH_MARGIN_M
    depth_max = space_extents["depth_max"] + DEPTH_MARGIN_M
    z_min = min(0.0, space_extents["z_min"]) - 0.05
    z_max = max(space_extents["z_max"], 2.4)
    for record in object_records:
        extents = projection_extents(record["bounds"], anchor_m, forward_xy)
        z_min = min(z_min, extents["z_min"] - 0.05)
        z_max = max(z_max, extents["z_max"] + VERTICAL_MARGIN_M)

    lateral_center = (lateral_min + lateral_max) / 2.0
    z_center = (z_min + z_max) / 2.0
    target_x = anchor_m[0] + lateral_xy[0] * lateral_center + forward_xy[0] * space_extents["depth_min"]
    target_y = anchor_m[1] + lateral_xy[1] * lateral_center + forward_xy[1] * space_extents["depth_min"]
    camera_distance = 1.0
    camera.location = (
        target_x - forward_xy[0] * camera_distance,
        target_y - forward_xy[1] * camera_distance,
        z_center,
    )
    camera.rotation_euler = _vector((forward_xy[0], forward_xy[1], 0.0)).to_track_quat("-Z", "Y").to_euler()

    width = lateral_max - lateral_min
    height = max(MIN_ORTHO_HEIGHT_M, z_max - z_min)
    aspect = RENDER_WIDTH_PX / RENDER_HEIGHT_PX
    camera.data.ortho_scale = max(height, width / aspect)
    visible_height = camera.data.ortho_scale
    visible_width = visible_height * aspect
    frame_u_min = lateral_center - visible_width / 2.0
    frame_u_max = lateral_center + visible_width / 2.0
    frame_z_min = z_center - visible_height / 2.0
    frame_z_max = z_center + visible_height / 2.0
    camera.data.clip_start = 0.01
    camera.data.clip_end = camera_distance + (depth_max - space_extents["depth_min"]) + 1.0
    return {
        "location_mm": [round(value * 1000.0, 6) for value in camera.location],
        "rotation_euler_rad": [round(value, 12) for value in camera.rotation_euler],
        "ortho_scale_mm": round(camera.data.ortho_scale * 1000.0, 6),
        "clip_start_mm": round(camera.data.clip_start * 1000.0, 6),
        "clip_end_mm": round(camera.data.clip_end * 1000.0, 6),
        "resolution_px": [RENDER_WIDTH_PX, RENDER_HEIGHT_PX],
    }, {
        "u_min_mm": round(frame_u_min * 1000.0, 6),
        "u_max_mm": round(frame_u_max * 1000.0, 6),
        "z_min_mm": round(frame_z_min * 1000.0, 6),
        "z_max_mm": round(frame_z_max * 1000.0, 6),
    }


def render_view(
    bpy: Any,
    root: Path,
    row: dict[str, Any],
    space: Any,
    space_bounds: Bounds,
    records: list[dict[str, Any]],
    camera: Any,
) -> dict[str, Any]:
    forward_xy = DIRECTION_VECTORS[row["direction"]]
    anchor_m = (row["anchor_x_mm"] / 1000.0, row["anchor_y_mm"] / 1000.0)
    camera_state, frame = camera_for_view(camera, row, space_bounds, records)
    output_path = root / OUTPUT_DIR_REL / f"INT1-ELEV-{row['view_id']}.png"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    render_scene = None
    collection = None
    copies: list[Any] = []
    try:
        render_scene = create_isolated_render_scene(bpy, bpy.context.scene, camera)
        collection, copies = create_sanitized_render_copies(bpy, render_scene, records)
        render_scene.render.filepath = str(output_path)
        bpy.ops.render.render(write_still=True, scene=render_scene.name)
    finally:
        remove_sanitized_render_copies(bpy, collection, copies)
        if render_scene is not None and render_scene.name in bpy.data.scenes:
            bpy.data.scenes.remove(render_scene)
    if not output_path.is_file() or output_path.stat().st_size == 0:
        raise RuntimeError(f"Blender did not create elevation PNG: {output_path}")

    objects = []
    for record in records:
        entity = record["entity"]
        objects.append(
            {
                "global_id": str(getattr(entity, "GlobalId", "")),
                "ifc_class": entity.is_a(),
                "name": str(getattr(entity, "Name", "") or ""),
                "container": entity_container_name(entity),
                "bbox_min_mm": bbox_min_mm(record["bounds"]),
                "bbox_max_mm": bbox_max_mm(record["bounds"]),
                "projected": projected_bbox_json(record["bounds"], anchor_m, forward_xy),
                "source": "world_bbox",
            }
        )
    objects.sort(key=lambda item: (item["ifc_class"], item["global_id"]))
    return {
        "view_id": row["view_id"],
        "sheet_id": row["sheet_id"],
        "direction": row["direction"],
        "anchor_id": row["anchor_id"],
        "space_reference": row["space_reference"] or row["space_name"] or str(getattr(space, "Name", "") or ""),
        "space_global_id": str(space.GlobalId),
        "png": {
            "path": str(output_path.relative_to(root)),
            "sha256": sha256_file(output_path),
            "width_px": RENDER_WIDTH_PX,
            "height_px": RENDER_HEIGHT_PX,
        },
        "camera": camera_state,
        "frame": frame,
        "objects": objects,
        "demolish_visible_count": 0,
    }


def validate_manifest_schema(manifest: dict[str, Any]) -> None:
    top_required = {"source_ifc_sha256", "view_register_sha256", "generated_at", "views"}
    view_required = {
        "view_id",
        "sheet_id",
        "direction",
        "anchor_id",
        "space_reference",
        "space_global_id",
        "png",
        "camera",
        "frame",
        "objects",
    }
    png_required = {"path", "sha256", "width_px", "height_px"}
    camera_required = {
        "location_mm",
        "rotation_euler_rad",
        "ortho_scale_mm",
        "clip_start_mm",
        "clip_end_mm",
        "resolution_px",
    }
    frame_required = {"u_min_mm", "u_max_mm", "z_min_mm", "z_max_mm"}
    object_required = {
        "global_id",
        "ifc_class",
        "name",
        "container",
        "bbox_min_mm",
        "bbox_max_mm",
        "projected",
        "source",
    }
    projected_required = {
        "u_min_mm",
        "u_max_mm",
        "z_min_mm",
        "z_max_mm",
        "depth_min_mm",
        "depth_max_mm",
    }
    missing = top_required - manifest.keys()
    if missing:
        raise ValueError(f"manifest missing top-level fields: {sorted(missing)}")
    if len(manifest["views"]) != EXPECTED_VIEW_COUNT:
        raise ValueError(f"manifest must contain {EXPECTED_VIEW_COUNT} views")
    if manifest.get("demolish_visible_count") != 0:
        raise ValueError("manifest demolish_visible_count must be 0")
    for view in manifest["views"]:
        if missing := view_required - view.keys():
            raise ValueError(f"view {view.get('view_id')} missing fields: {sorted(missing)}")
        if view.get("demolish_visible_count") != 0:
            raise ValueError(f"view {view['view_id']} includes visible DEMOLISH walls")
        for label, value, required in (
            ("png", view["png"], png_required),
            ("camera", view["camera"], camera_required),
            ("frame", view["frame"], frame_required),
        ):
            if missing := required - value.keys():
                raise ValueError(f"view {view['view_id']} {label} missing fields: {sorted(missing)}")
        for item in view["objects"]:
            if missing := object_required - item.keys():
                raise ValueError(f"view {view['view_id']} object missing fields: {sorted(missing)}")
            if item["source"] != "world_bbox":
                raise ValueError(f"view {view['view_id']} object source must be world_bbox")
            if missing := projected_required - item["projected"].keys():
                raise ValueError(f"view {view['view_id']} projected object missing fields: {sorted(missing)}")


def run() -> dict[str, Any]:
    import bpy
    import bonsai.tool as tool

    model = tool.Ifc.get()
    if model is None:
        raise RuntimeError("No IFC is loaded in Bonsai")
    ifc_path = Path(str(tool.Ifc.get_path())).resolve()
    root = project_root_from_ifc_path(str(ifc_path))
    view_register_path = root / VIEW_REGISTER_REL
    rows = load_view_rows(view_register_path)
    demolition_register = resolve_demolition_register(root)
    demolition_ids = load_demolition_ids(demolition_register)
    ifc_hash = sha256_file(ifc_path)

    state = capture_scene_state(bpy)
    temporary_camera = None
    views: list[dict[str, Any]] = []
    engine = ""
    restoration_checks: dict[str, bool] = {}
    try:
        engine = configure_cad_like_render(bpy)
        temporary_camera = create_temporary_camera(bpy)
        for row in rows:
            space, _space_obj, space_bounds = resolve_space(model, tool, row)
            records = collect_render_objects(bpy, model, tool, space_bounds, demolition_ids)
            views.append(render_view(bpy, root, row, space, space_bounds, records, temporary_camera))
    finally:
        restoration_checks = restore_scene_state(bpy, state, temporary_camera)

    if len(views) != EXPECTED_VIEW_COUNT:
        raise RuntimeError(f"rendered {len(views)} of {EXPECTED_VIEW_COUNT} elevation views")
    if not all(restoration_checks.values()):
        raise RuntimeError(f"scene state restoration failed: {restoration_checks}")

    manifest = {
        "schema_version": 1,
        "generator": "pipeline/scripts/int1_elevation_blender_render.py",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "ifc_path": str(ifc_path),
        "source_ifc_sha256": ifc_hash,
        "view_register_sha256": sha256_file(view_register_path),
        "view_register": str(VIEW_REGISTER_REL),
        "demolition_register": str(demolition_register.relative_to(root)),
        "render_engine": engine,
        "render_style": "solid true-depth grey-white entities with black Workbench outlines",
        "view_count": len(views),
        "output_directory": str(OUTPUT_DIR_REL),
        "excluded_ifc_classes": sorted(EXCLUDED_IFC_CLASSES),
        "excluded_demolition_global_ids": sorted(demolition_ids),
        "demolish_visible_count": 0,
        "state_restored": True,
        "restoration_checks": restoration_checks,
        "formal_ifc_write": False,
        "blend_save": False,
        "views": views,
    }
    validate_manifest_schema(manifest)
    manifest_path = root / MANIFEST_REL
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return manifest


if __name__ == "__main__":
    RESULT = run()
    print(json.dumps({"view_count": RESULT["view_count"], "manifest": str(MANIFEST_REL)}, ensure_ascii=False))
