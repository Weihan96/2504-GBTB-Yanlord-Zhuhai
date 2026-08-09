"""Install a rebuildable Blender HVAC constraint preview without writing IFC."""

from __future__ import annotations

import json
import hashlib
from pathlib import Path

import bpy
import bonsai.tool as tool
from bpy.props import StringProperty
from mathutils import Matrix, Vector


COLLECTION_NAME = "RCP1_HVAC_CONSTRAINT_AUTHORING"
NODE_GROUP_NAME = "RCP1_HVAC_CONSTRAINT_GEOMETRY"
MATERIAL_NAME = "RCP1 HVAC Constraint Preview"
ROUTE_PREFIX = "RCP1_HVAC_ROUTE_"
ANCHOR_PREFIX = "RCP1_HVAC_ANCHOR_"
DRIVER_CLASSES_KEY = "rcp1_hvac_constraint_preview_classes"
DRIVER_HANDLER_KEY = "rcp1_hvac_constraint_preview_handler"
DRIVER_STATE_KEY = "rcp1_hvac_constraint_preview_state"
PREVIEW_RADIUS_M = 0.015
FILLET_RADIUS_M = 0.12
LEGACY_REVIEW_COLLECTIONS = (
    "RCP1_HVAC_REMODEL_REVIEW",
    "RCP1_PAIR_REVIEW_TEMP",
    "RCP1_ROUTE_READINESS_REVIEW",
    "RCP1_ROUTE_CONSTRAINT_REVIEW",
)


def project_root() -> Path:
    return Path(str(tool.Ifc.get_path())).resolve().parent


def load_report() -> dict:
    return json.loads((project_root() / "build/rcp1/route-readiness-candidate.json").read_text(encoding="utf-8"))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def ensure_collection() -> bpy.types.Collection:
    collection = bpy.data.collections.get(COLLECTION_NAME)
    if collection is None:
        collection = bpy.data.collections.new(COLLECTION_NAME)
        bpy.context.scene.collection.children.link(collection)
    return collection


def remove_generated_objects(collection: bpy.types.Collection) -> None:
    for obj in list(collection.objects):
        bpy.data.objects.remove(obj, do_unlink=True)


def ensure_material() -> bpy.types.Material:
    material = bpy.data.materials.get(MATERIAL_NAME)
    if material is None:
        material = bpy.data.materials.new(MATERIAL_NAME)
    material.diffuse_color = (0.0, 0.75, 1.0, 1.0)
    return material


def ensure_node_group() -> bpy.types.GeometryNodeTree:
    group = bpy.data.node_groups.get(NODE_GROUP_NAME)
    if group is not None:
        bpy.data.node_groups.remove(group, do_unlink=True)
    group = bpy.data.node_groups.new(NODE_GROUP_NAME, "GeometryNodeTree")
    group.interface.new_socket(name="Geometry", in_out="INPUT", socket_type="NodeSocketGeometry")
    group.interface.new_socket(name="Geometry", in_out="OUTPUT", socket_type="NodeSocketGeometry")
    nodes = group.nodes
    links = group.links
    group_input = nodes.new("NodeGroupInput")
    group_output = nodes.new("NodeGroupOutput")
    fillet = nodes.new("GeometryNodeFilletCurve")
    fillet.mode = "POLY"
    fillet.inputs["Radius"].default_value = FILLET_RADIUS_M
    fillet.inputs["Count"].default_value = 6
    fillet.inputs["Limit Radius"].default_value = True
    profile = nodes.new("GeometryNodeCurvePrimitiveCircle")
    profile.mode = "RADIUS"
    profile.inputs["Resolution"].default_value = 12
    profile.inputs["Radius"].default_value = PREVIEW_RADIUS_M
    curve_to_mesh = nodes.new("GeometryNodeCurveToMesh")
    curve_to_mesh.inputs["Fill Caps"].default_value = True
    set_material = nodes.new("GeometryNodeSetMaterial")
    set_material.inputs["Material"].default_value = ensure_material()
    links.new(group_input.outputs["Geometry"], fillet.inputs["Curve"])
    links.new(fillet.outputs["Curve"], curve_to_mesh.inputs["Curve"])
    links.new(profile.outputs["Curve"], curve_to_mesh.inputs["Profile Curve"])
    links.new(curve_to_mesh.outputs["Mesh"], set_material.inputs["Geometry"])
    links.new(set_material.outputs["Geometry"], group_output.inputs["Geometry"])
    return group


def ifc_object(global_id: str | None) -> bpy.types.Object | None:
    if not global_id:
        return None
    model = tool.Ifc.get()
    entity = model.by_guid(global_id) if model is not None else None
    return tool.Ifc.get_object(entity) if entity is not None else None


def create_anchor(
    collection: bpy.types.Collection,
    anchor_id: str,
    anchor_kind: str,
    global_id: str | None,
    location_mm: list[float],
) -> bpy.types.Object:
    name = f"{ANCHOR_PREFIX}{anchor_id}"
    anchor = bpy.data.objects.get(name)
    if anchor is None:
        anchor = bpy.data.objects.new(name, None)
        collection.objects.link(anchor)
    anchor.empty_display_type = "SPHERE"
    anchor.empty_display_size = 0.065
    if anchor_kind == "blender_bend":
        anchor.color = (1.0, 0.28, 0.02, 1.0)
    elif anchor_id == "H01":
        anchor.color = (1.0, 0.16, 0.72, 1.0)
    elif anchor_id == "H07":
        anchor.color = (0.08, 0.42, 1.0, 1.0)
    else:
        anchor.color = (0.0, 0.75, 1.0, 1.0)
    anchor.show_in_front = False
    anchor["rcp1_anchor_id"] = anchor_id
    anchor["rcp1_anchor_kind"] = anchor_kind
    anchor["rcp1_global_id"] = global_id or ""
    world = Matrix.Translation(Vector(tuple(value / 1000.0 for value in location_mm)))
    parent = ifc_object(global_id)
    anchor.parent = parent
    anchor.matrix_world = world
    return anchor


def create_curve(
    collection: bpy.types.Collection,
    route_id: str,
    anchors: list[bpy.types.Object],
    node_group: bpy.types.GeometryNodeTree,
) -> bpy.types.Object:
    curve = bpy.data.curves.new(f"{ROUTE_PREFIX}{route_id}_DATA", "CURVE")
    curve.dimensions = "3D"
    curve.resolution_u = 2
    obj = bpy.data.objects.new(f"{ROUTE_PREFIX}{route_id}", curve)
    collection.objects.link(obj)
    obj["rcp1_route_id"] = route_id
    obj["rcp1_anchor_names"] = json.dumps([anchor.name for anchor in anchors])
    obj["rcp1_status"] = "constraint_skeleton_not_fabrication_geometry"
    obj.show_in_front = False
    modifier = obj.modifiers.new("RCP1 HVAC Constraint Geometry", "NODES")
    modifier.node_group = node_group
    update_curve(obj, anchors)
    return obj


def update_curve(route: bpy.types.Object, anchors: list[bpy.types.Object]) -> None:
    curve = route.data
    curve.splines.clear()
    if len(anchors) < 2:
        return
    spline = curve.splines.new("POLY")
    spline.points.add(len(anchors) - 1)
    inverse = route.matrix_world.inverted_safe()
    for point, anchor in zip(spline.points, anchors):
        local = inverse @ anchor.matrix_world.translation
        point.co = (*local, 1.0)


def anchor_signature() -> tuple:
    rows = []
    for obj in sorted(bpy.data.objects, key=lambda item: item.name):
        if not obj.name.startswith(ANCHOR_PREFIX):
            continue
        rows.append((obj.name, *(round(value, 9) for value in obj.matrix_world.translation)))
    return tuple(rows)


def rebuild_preview() -> int:
    updated = 0
    for route in bpy.data.objects:
        if not route.name.startswith(ROUTE_PREFIX) or route.type != "CURVE":
            continue
        names = json.loads(route.get("rcp1_anchor_names", "[]"))
        anchors = [bpy.data.objects.get(name) for name in names]
        if any(anchor is None for anchor in anchors):
            continue
        update_curve(route, anchors)
        updated += 1
    return updated


def automatic_preview_handler(_scene: bpy.types.Scene, _depsgraph: bpy.types.Depsgraph) -> None:
    state = bpy.app.driver_namespace.setdefault(DRIVER_STATE_KEY, {})
    signature = anchor_signature()
    if signature == state.get("signature"):
        return
    state["signature"] = signature
    rebuild_preview()


def install_handler() -> None:
    previous = bpy.app.driver_namespace.get(DRIVER_HANDLER_KEY)
    if previous in bpy.app.handlers.depsgraph_update_post:
        bpy.app.handlers.depsgraph_update_post.remove(previous)
    bpy.app.handlers.depsgraph_update_post.append(automatic_preview_handler)
    bpy.app.driver_namespace[DRIVER_HANDLER_KEY] = automatic_preview_handler
    bpy.app.driver_namespace[DRIVER_STATE_KEY] = {"signature": anchor_signature()}


class RCP1HVAC_OT_rebuild_preview(bpy.types.Operator):
    bl_idname = "rcp1_hvac.rebuild_preview"
    bl_label = "重建路线预览"
    bl_description = "从当前空调、洞口和弯点位置重建 Blender 预览；不写 IFC"

    def execute(self, context):
        count = rebuild_preview()
        context.scene["rcp1_hvac_last_action"] = f"已重建 {count} 条路线；IFC 未修改"
        self.report({"INFO"}, context.scene["rcp1_hvac_last_action"])
        return {"FINISHED"}


class RCP1HVAC_OT_add_bend(bpy.types.Operator):
    bl_idname = "rcp1_hvac.add_bend"
    bl_label = "在光标处添加弯点"
    bl_description = "在所选路线终点前添加临时弯点并自动更新预览；不写 CSV 或 IFC"

    def execute(self, context):
        route_id = context.scene.rcp1_hvac_active_route.strip()
        route = bpy.data.objects.get(f"{ROUTE_PREFIX}{route_id}")
        if route is None:
            self.report({"ERROR"}, f"找不到路线 {route_id}")
            return {"CANCELLED"}
        names = json.loads(route.get("rcp1_anchor_names", "[]"))
        index = 1 + sum(1 for name in names if "_BEND_" in name)
        anchor_id = f"{route_id}_BEND_{index:02d}"
        anchor = create_anchor(
            ensure_collection(),
            anchor_id,
            "blender_bend",
            None,
            [value * 1000.0 for value in context.scene.cursor.location],
        )
        insert_at = max(1, len(names) - 1)
        names.insert(insert_at, anchor.name)
        route["rcp1_anchor_names"] = json.dumps(names)
        rebuild_preview()
        context.view_layer.objects.active = anchor
        anchor.select_set(True)
        context.scene["rcp1_hvac_last_action"] = f"已添加临时弯点 {anchor_id}；尚未写入决策表或 IFC"
        self.report({"INFO"}, context.scene["rcp1_hvac_last_action"])
        return {"FINISHED"}


class RCP1HVAC_OT_delete_bend(bpy.types.Operator):
    bl_idname = "rcp1_hvac.delete_bend"
    bl_label = "删除所选弯点"
    bl_description = "删除所选临时弯点并重建预览；不写 IFC"

    def execute(self, context):
        bend = context.active_object
        if bend is None or bend.get("rcp1_anchor_kind") != "blender_bend":
            self.report({"ERROR"}, "请先选择一个橙色临时弯点")
            return {"CANCELLED"}
        for route in bpy.data.objects:
            if not route.name.startswith(ROUTE_PREFIX) or route.type != "CURVE":
                continue
            names = json.loads(route.get("rcp1_anchor_names", "[]"))
            if bend.name in names:
                names.remove(bend.name)
                route["rcp1_anchor_names"] = json.dumps(names)
        name = bend.name
        bpy.data.objects.remove(bend, do_unlink=True)
        rebuild_preview()
        context.scene["rcp1_hvac_last_action"] = f"已删除临时弯点 {name}；IFC 未修改"
        self.report({"INFO"}, context.scene["rcp1_hvac_last_action"])
        return {"FINISHED"}


class RCP1HVAC_OT_export_candidate(bpy.types.Operator):
    bl_idname = "rcp1_hvac.export_candidate"
    bl_label = "导出路线候选"
    bl_description = "把当前 Blender 锚点顺序和位置导出到 build 供机械检查；不写决策表或 IFC"

    def execute(self, context):
        ifc_path = Path(str(tool.Ifc.get_path())).resolve()
        routes = []
        for route in sorted(bpy.data.objects, key=lambda item: item.name):
            if not route.name.startswith(ROUTE_PREFIX) or route.type != "CURVE":
                continue
            anchors = []
            for name in json.loads(route.get("rcp1_anchor_names", "[]")):
                anchor = bpy.data.objects.get(name)
                if anchor is None:
                    continue
                anchors.append({
                    "anchor_id": anchor.get("rcp1_anchor_id", anchor.name),
                    "anchor_kind": anchor.get("rcp1_anchor_kind", "unknown"),
                    "global_id": anchor.get("rcp1_global_id", "") or None,
                    "world_mm": [round(value * 1000.0, 6) for value in anchor.matrix_world.translation],
                })
            routes.append({
                "route_id": route.get("rcp1_route_id", route.name),
                "status": "blender_preview_candidate_human_review_required",
                "anchors": anchors,
            })
        output = project_root() / "build/rcp1/hvac-route-preview-candidate.json"
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(json.dumps({
            "mode": "blender_hvac_route_preview_candidate",
            "source_ifc": str(ifc_path),
            "source_ifc_sha256": sha256(ifc_path),
            "source_of_truth": "IFC plus approved route and waypoint registers",
            "routes": routes,
            "preview_parameters": {
                "geometry_nodes": NODE_GROUP_NAME,
                "preview_radius_m": PREVIEW_RADIUS_M,
                "fillet_radius_m": FILLET_RADIUS_M,
                "fabrication_geometry": False,
            },
            "formal_ifc_write_allowed": False,
        }, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
        context.scene["rcp1_hvac_last_action"] = f"已导出 {len(routes)} 条路线候选；决策表和 IFC 未修改"
        self.report({"INFO"}, context.scene["rcp1_hvac_last_action"])
        return {"FINISHED"}


class RCP1HVAC_PT_constraint_authoring(bpy.types.Panel):
    bl_label = "RCP1 HVAC 约束路线"
    bl_idname = "RCP1HVAC_PT_constraint_authoring"
    bl_space_type = "VIEW_3D"
    bl_region_type = "UI"
    bl_category = "RCP1 HVAC"

    def draw(self, context):
        layout = self.layout
        layout.label(text="青色：已确认锚点；橙色：临时弯点")
        layout.label(text="预览会自动更新；不会自动写 IFC")
        layout.operator("rcp1_hvac.rebuild_preview")
        layout.prop(context.scene, "rcp1_hvac_active_route", text="路线 ID")
        row = layout.row(align=True)
        row.operator("rcp1_hvac.add_bend")
        row.operator("rcp1_hvac.delete_bend")
        layout.operator("rcp1_hvac.export_candidate")
        layout.label(text=context.scene.get("rcp1_hvac_last_action", "等待操作"))


CLASSES = (
    RCP1HVAC_OT_rebuild_preview,
    RCP1HVAC_OT_add_bend,
    RCP1HVAC_OT_delete_bend,
    RCP1HVAC_OT_export_candidate,
    RCP1HVAC_PT_constraint_authoring,
)


def register_ui() -> None:
    previous_classes = bpy.app.driver_namespace.get(DRIVER_CLASSES_KEY, ())
    for cls in reversed(previous_classes):
        try:
            bpy.utils.unregister_class(cls)
        except (RuntimeError, ValueError):
            pass
    for cls in CLASSES:
        bpy.utils.register_class(cls)
    if not hasattr(bpy.types.Scene, "rcp1_hvac_active_route"):
        bpy.types.Scene.rcp1_hvac_active_route = StringProperty(default="RCP1-SERVICE-A02")
    bpy.app.driver_namespace[DRIVER_CLASSES_KEY] = CLASSES


def build_authoring_scene() -> dict:
    report = load_report()
    collection = ensure_collection()
    remove_generated_objects(collection)
    node_group = ensure_node_group()
    anchors_by_id: dict[str, bpy.types.Object] = {}
    routes_created = []
    hidden_review_collections = []
    for name in LEGACY_REVIEW_COLLECTIONS:
        legacy_collection = bpy.data.collections.get(name)
        if legacy_collection is None:
            continue
        legacy_collection.hide_viewport = True
        legacy_collection.hide_render = True
        hidden_review_collections.append(name)
    for route in report["confirmed_route_graph"]:
        anchors = []
        for waypoint in route["waypoints"]:
            anchor_id = waypoint["anchor_id"]
            anchor = anchors_by_id.get(anchor_id)
            if anchor is None:
                anchor = create_anchor(
                    collection,
                    anchor_id,
                    waypoint["anchor_kind"],
                    waypoint.get("anchor_global_id") or None,
                    waypoint["centre_mm"],
                )
                anchors_by_id[anchor_id] = anchor
            anchors.append(anchor)
        if len(anchors) < 2:
            continue
        create_curve(collection, route["route_id"], anchors, node_group)
        routes_created.append(route["route_id"])
    register_ui()
    install_handler()
    for obj in bpy.context.scene.objects:
        obj.show_in_front = False
    if bpy.context.screen is not None:
        for area in bpy.context.screen.areas:
            if area.type != "VIEW_3D":
                continue
            area.spaces.active.shading.type = "SOLID"
            area.spaces.active.shading.color_type = "MATERIAL"
            area.spaces.active.shading.show_xray = False
            area.spaces.active.overlay.show_wireframes = False
    return {
        "collection": COLLECTION_NAME,
        "routes": routes_created,
        "anchor_count": len(anchors_by_id),
        "endpoint_anchors": [anchor_id for anchor_id in ("H01", "H07") if anchor_id in anchors_by_id],
        "hidden_legacy_review_collections": hidden_review_collections,
        "geometry_nodes": NODE_GROUP_NAME,
        "automatic_preview_update": True,
        "ifc_write": False,
        "blend_save": False,
        "show_in_front": False,
        "xray": False,
        "wireframes": False,
    }


RESULT = build_authoring_scene()
print(RESULT)
