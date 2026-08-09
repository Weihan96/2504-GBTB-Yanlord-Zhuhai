"""Prepare a three-colour Blender review for the C003→INT1 handoff."""

from __future__ import annotations

import csv
from pathlib import Path

import bpy
import bonsai.tool as tool
import ifcopenshell.util.element


ROOT_COLLECTION = "INT1_HANDOFF_REVIEW"
GROUPS = {
    "named": ("01_NAMED_GREEN", (0.10, 0.80, 0.20, 1.0)),
    "geometry": ("02_GEOMETRY_ORANGE", (1.00, 0.35, 0.05, 1.0)),
    "legacy": ("03_LEGACY_RED", (1.00, 0.08, 0.08, 1.0)),
}
FURNITURE_YELLOW = (1.0, 0.65, 0.0, 1.0)


def project_root() -> Path:
    return Path(str(tool.Ifc.get_path())).resolve().parent


def read_rows() -> list[dict[str, str]]:
    path = project_root() / "pipeline/decisions/int1-handoff-review.csv"
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def ensure_collection(name: str, parent: bpy.types.Collection) -> bpy.types.Collection:
    collection = bpy.data.collections.get(name)
    if collection is None:
        collection = bpy.data.collections.new(name)
    if collection not in parent.children.values():
        parent.children.link(collection)
    return collection


def clear_review_links(collection: bpy.types.Collection) -> None:
    for obj in list(collection.objects):
        collection.objects.unlink(obj)


def group_key(role: str) -> str:
    if role.startswith("named_"):
        return "named"
    if role == "legacy_cad_reference_needs_disposition":
        return "legacy"
    return "geometry"


def configure_viewport() -> None:
    for screen in bpy.data.screens:
        for area in screen.areas:
            if area.type != "VIEW_3D":
                continue
            area.spaces.active.shading.type = "SOLID"
            area.spaces.active.shading.color_type = "OBJECT"
            area.spaces.active.shading.show_xray = False
            area.spaces.active.overlay.show_wireframes = False
            area.spaces.active.overlay.show_relationship_lines = False


def build_review() -> dict:
    model = tool.Ifc.get()
    if model is None:
        raise RuntimeError("No IFC is loaded in Bonsai")
    rows = read_rows()
    if len(rows) != 16:
        raise RuntimeError(f"expected 16 INT1 handoff records, found {len(rows)}")

    root = ensure_collection(ROOT_COLLECTION, bpy.context.scene.collection)
    review_groups = {}
    for key, (name, _colour) in GROUPS.items():
        collection = ensure_collection(name, root)
        clear_review_links(collection)
        review_groups[key] = collection

    for furniture in model.by_type("IfcFurniture"):
        obj = tool.Ifc.get_object(furniture)
        if obj is not None:
            obj.color = FURNITURE_YELLOW
            obj.show_in_front = False

    counts = {key: 0 for key in GROUPS}
    review_objects = []
    for row in rows:
        entity = model.by_guid(row["global_id"])
        obj = tool.Ifc.get_object(entity) if entity is not None else None
        if obj is None:
            raise RuntimeError(f"Blender object not found for {row['global_id']}")
        key = group_key(row["candidate_role"])
        collection = review_groups[key]
        if obj not in collection.objects.values():
            collection.objects.link(obj)
        obj.color = GROUPS[key][1]
        obj.show_in_front = False
        obj.hide_set(False)
        obj.hide_viewport = False
        obj["int1_candidate_role"] = row["candidate_role"]
        obj["int1_review_status"] = row["review_status"]
        counts[key] += 1
        review_objects.append(obj)

    for wall in model.by_type("IfcWall"):
        status = ifcopenshell.util.element.get_psets(wall).get("Pset_WallCommon", {}).get("Status")
        if status == "DEMOLISH":
            obj = tool.Ifc.get_object(wall)
            if obj is not None:
                obj.hide_set(True)

    for name in ("A102_DEMOLITION_REFERENCE", "RCP1_HVAC_CONSTRAINT_AUTHORING"):
        collection = bpy.data.collections.get(name)
        if collection is not None:
            collection.hide_viewport = True
            collection.hide_render = True

    for obj in list(bpy.context.selected_objects):
        obj.select_set(False)
    configure_viewport()
    bpy.context.scene["int1_review_legend"] = "green=named; orange=geometry-only; red=legacy CAD"
    bpy.context.scene["int1_review_source"] = str(project_root() / "pipeline/decisions/int1-handoff-review.csv")
    return {
        "review_objects": len(review_objects),
        "named_green": counts["named"],
        "geometry_orange": counts["geometry"],
        "legacy_red": counts["legacy"],
        "show_in_front": False,
        "solid_view": True,
        "xray": False,
        "wireframes": False,
        "formal_ifc_write": False,
        "blend_save": False,
    }


RESULT = build_review()
print(RESULT)
