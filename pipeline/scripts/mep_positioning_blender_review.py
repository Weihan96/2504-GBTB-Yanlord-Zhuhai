"""Prepare a depth-aware Blender review of current water and electrical positioning."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import bpy
import bonsai.tool as tool


ROOT_COLLECTION = "MEP_POSITIONING_REVIEW"
GROUPS = {
    "plum": ("01_PLUM_BLUE", (0.05, 0.30, 1.00, 1.0)),
    "elec": ("02_ELEC_PURPLE", (0.75, 0.08, 1.00, 1.0)),
    "review": ("03_REVIEW_RED", (1.00, 0.03, 0.03, 1.0)),
}
FURNITURE_YELLOW = (1.0, 0.65, 0.0, 1.0)
OFF_PLAN_WASTE_TERMINAL_ID = "3IVqCnhGr51hY4LrOq_5G_"


def project_root() -> Path:
    return Path(str(tool.Ifc.get_path())).resolve().parent


def read_json(relative: str) -> dict:
    return json.loads((project_root() / relative).read_text(encoding="utf-8"))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


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


def link_product(model, global_id: str, collection: bpy.types.Collection, colour, role: str):
    entity = model.by_guid(global_id)
    obj = tool.Ifc.get_object(entity) if entity is not None else None
    if obj is None:
        raise RuntimeError(f"Blender object not found for {global_id}")
    if obj not in collection.objects.values():
        collection.objects.link(obj)
    obj.color = colour
    obj.show_in_front = False
    obj.hide_set(False)
    obj.hide_viewport = False
    obj["mep_review_role"] = role
    return obj


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
    current_hash = sha256(Path(str(tool.Ifc.get_path())))
    plum = read_json("build/plum/p201-demand-endpoints.json")
    p202 = read_json("build/plum/p202-existing-location-register.json")
    plum_connections = read_json("build/plum/plum-connection-candidate.json")
    elec = read_json("build/elec/elec-existing-candidate.json")
    elec_positions = read_json("build/elec/elec-positioning-candidate.json")
    source = read_json("build/mep-positioning/source-audit.json")
    hashes = {
        plum["source_ifc_sha256"], p202["source_ifc_sha256"],
        plum_connections["source_ifc_sha256"], elec["source"]["sha256"],
        elec_positions["source_ifc_sha256"], source["source_ifc_sha256"],
    }
    if hashes != {current_hash}:
        raise RuntimeError("MEP review inputs are stale relative to loaded IFC")

    root = ensure_collection(ROOT_COLLECTION, bpy.context.scene.collection)
    collections = {}
    for key, (name, _colour) in GROUPS.items():
        collection = ensure_collection(name, root)
        clear_review_links(collection)
        collection.hide_viewport = False
        collection.hide_render = False
        collections[key] = collection

    service_ids = {
        row["global_id"] for row in plum["demand_endpoints"] if row["service_demand_candidate"]
    }
    non_service = {
        row["global_id"]: row["candidate_role"]
        for row in plum["demand_endpoints"] if not row["service_demand_candidate"]
    }
    p202_ids = {row["global_id"] for row in p202["objects"]}
    plum_ids = (p202_ids - non_service.keys()) - {OFF_PLAN_WASTE_TERMINAL_ID}
    for global_id in sorted(plum_ids):
        role = "service_demand_candidate" if global_id in service_ids else "existing_drainage_or_assembly"
        link_product(model, global_id, collections["plum"], GROUPS["plum"][1], role)

    e301 = elec["sheets"]["E-301"]
    e303 = elec["sheets"]["E-303"]
    elec_roles = {
        row["global_id"]: row["candidate_role"] for row in elec_positions["socket_candidates"]
    }
    established_elec = e301["lights"] + e303["sockets"] + e303["typed_equipment"]
    for row in established_elec:
        role = elec_roles.get(row["global_id"], "existing_light_or_typed_equipment")
        link_product(model, row["global_id"], collections["elec"], GROUPS["elec"][1], role)

    review_roles = dict(non_service)
    review_roles[OFF_PLAN_WASTE_TERMINAL_ID] = "off_plan_waste_terminal"
    review_roles.update({
        row["global_id"]: row["candidate_role"]
        for row in elec_positions["proxy_identity_candidates"]
    })
    for global_id, role in sorted(review_roles.items()):
        link_product(model, global_id, collections["review"], GROUPS["review"][1], role)

    for furniture in model.by_type("IfcFurniture"):
        obj = tool.Ifc.get_object(furniture)
        if obj is not None:
            obj.color = FURNITURE_YELLOW
            obj.show_in_front = False

    for name in (
        "A102_DEMOLITION_REFERENCE",
        "RCP1_HVAC_CONSTRAINT_AUTHORING",
        "INT1_HANDOFF_REVIEW",
    ):
        collection = bpy.data.collections.get(name)
        if collection is not None:
            collection.hide_viewport = True
            collection.hide_render = True

    for obj in list(bpy.context.selected_objects):
        obj.select_set(False)
    configure_viewport()
    bpy.context.scene["mep_review_legend"] = "blue=PLUM; purple=established ELEC; red=semantic/off-plan/unresolved"
    bpy.context.scene["mep_review_source_pdf"] = source["source_drawings"]["latest"]["path"]
    return {
        "plum_blue": len(plum_ids),
        "elec_purple": len(established_elec),
        "review_red": len(review_roles),
        "service_demand_candidates": len(service_ids),
        "show_in_front": False,
        "solid_view": True,
        "xray": False,
        "wireframes": False,
        "formal_ifc_write": False,
        "blend_save": False,
    }


RESULT = build_review()
print(RESULT)
