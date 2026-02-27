# pyright: reportInvalidTypeForm=false
# ============================================================
# Outline2D (Batch IFC) + Stable UI Logs (CollectionProperty)
# + Optional BonsaiBIM: Add Representation (From Object)
#
# v0.5.0  (single-file final)
# - Avoids MULTILINE StringProperty registration completely
# - Logs shown in N-panel via UIList, supports Copy/Clear
# - Supports multi-select IFC objects (process each)
# - Optional preprocess: split each IFC obj by MATERIAL then LOOSE
# - Optional preprocess: remove unused materials
# - Optional: place filled outline at bbox center
# - Optional: join per IFC object (parts -> one outline)
# - Optional: set joined outline origin = original IFC object
# - Optional: BonsaiBIM add_representation(method='OBJECT'), auto-disabled if Bonsai missing
# ============================================================

bl_info = {
    "name": "Outline2D (Batch IFC) + Stable Logs + Bonsai Representation",
    "author": "ChatGPT",
    "version": (0, 5, 0),
    "blender": (4, 0, 0),
    "location": "View3D > Sidebar (N) > Outline2D",
    "category": "Object",
}

import bpy
import bmesh
from mathutils import Vector, Matrix
from datetime import datetime


# ============================================================
# Availability
# ============================================================
def bonsai_is_available() -> bool:
    try:
        return hasattr(bpy.ops, "bim") and hasattr(bpy.ops.bim, "add_representation")
    except Exception:
        return False


def _ensure_shapely():
    try:
        from shapely.geometry import Polygon
        from shapely.ops import unary_union
        return Polygon, unary_union, None
    except Exception as e:
        return None, None, e


# ============================================================
# Stable Logs (NO MULTILINE StringProperty)
# ============================================================
def _ts() -> str:
    return datetime.now().strftime("%H:%M:%S")


class OUTLINE2D_LogLine(bpy.types.PropertyGroup):
    line: bpy.props.StringProperty(name="Line", default="")


class OUTLINE2D_LogBuffer(bpy.types.PropertyGroup):
    lines: bpy.props.CollectionProperty(type=OUTLINE2D_LogLine)
    active_index: bpy.props.IntProperty(name="Active Index", default=0)


def _logbuf(context: bpy.types.Context) -> OUTLINE2D_LogBuffer:
    return context.scene.outline2d_logbuf


def log_ui(context: bpy.types.Context, msg: str):
    """Append one log line to buffer + print to console."""
    props = context.scene.outline2d_props
    buf = _logbuf(context)

    line = f"[{_ts()}] {msg}"
    print(line)

    it = buf.lines.add()
    it.line = line

    max_lines = max(20, int(getattr(props, "ui_log_max_lines", 200)))
    extra = len(buf.lines) - max_lines
    if extra > 0:
        for _ in range(extra):
            buf.lines.remove(0)

    buf.active_index = max(0, len(buf.lines) - 1)


def clear_ui_logs(context: bpy.types.Context):
    buf = _logbuf(context)
    buf.lines.clear()
    buf.active_index = 0


def get_logs_as_text(context: bpy.types.Context) -> str:
    buf = _logbuf(context)
    return "\n".join([it.line for it in buf.lines])


class OUTLINE2D_OT_clear_logs(bpy.types.Operator):
    bl_idname = "outline2d.clear_logs"
    bl_label = "Clear Logs"
    bl_options = {"INTERNAL"}

    def execute(self, context):
        clear_ui_logs(context)
        return {"FINISHED"}


class OUTLINE2D_OT_copy_logs(bpy.types.Operator):
    bl_idname = "outline2d.copy_logs"
    bl_label = "Copy Logs"
    bl_options = {"INTERNAL"}

    def execute(self, context):
        context.window_manager.clipboard = get_logs_as_text(context)
        self.report({"INFO"}, "Logs copied to clipboard.")
        return {"FINISHED"}


class OUTLINE2D_UL_logs(bpy.types.UIList):
    def draw_item(self, context, layout, data, item, icon, active_data, active_propname, index):
        layout.label(text=getattr(item, "line", ""))


# ============================================================
# Materials
# ============================================================
def _copy_material_slots(src_obj: bpy.types.Object, dst_mesh: bpy.types.Mesh):
    dst_mesh.materials.clear()
    for slot in src_obj.material_slots:
        if slot.material:
            dst_mesh.materials.append(slot.material)


def _remove_unused_material_slots(obj: bpy.types.Object):
    """Remove unused material slots; remap material indices."""
    if not obj or obj.type != "MESH":
        return

    mesh = obj.data
    mats = mesh.materials
    if not mats or len(mats) == 0:
        return

    used = set(int(p.material_index) for p in mesh.polygons)
    used = {i for i in used if 0 <= i < len(mats)}

    if len(used) == 0:
        mats.clear()
        obj.active_material_index = 0
        return

    keep = [i for i in range(len(mats)) if i in used]
    if len(keep) == len(mats):
        return

    mapping = {old_i: new_i for new_i, old_i in enumerate(keep)}
    new_mats = [mats[i] for i in keep]

    mats.clear()
    for m in new_mats:
        mats.append(m)

    for p in mesh.polygons:
        p.material_index = mapping.get(int(p.material_index), 0)

    obj.active_material_index = mapping.get(int(obj.active_material_index), 0)


# ============================================================
# Geometry helpers
# ============================================================
def _set_selection(context: bpy.types.Context, active_obj, selected_objs):
    vl = context.view_layer
    for o in vl.objects:
        o.select_set(False)
    for o in selected_objs:
        if o and o.name in bpy.data.objects:
            o.select_set(True)
    vl.objects.active = active_obj


def _bbox_center_world(obj: bpy.types.Object) -> Vector:
    pts = [obj.matrix_world @ Vector(corner) for corner in obj.bound_box]
    c = Vector((0.0, 0.0, 0.0))
    for p in pts:
        c += p
    return c / 8.0


def _face_loop_world_xy(face: bmesh.types.BMFace, world_mat: Matrix):
    pts = []
    for v in face.verts:
        p = world_mat @ v.co
        pts.append((float(p.x), float(p.y)))

    if len(pts) < 3:
        return None

    cleaned = []
    for p in pts:
        if not cleaned or p != cleaned[-1]:
            cleaned.append(p)

    if len(cleaned) < 3:
        return None

    if cleaned[0] != cleaned[-1]:
        cleaned.append(cleaned[0])

    if len(cleaned) < 4:
        return None

    return cleaned


def _set_origin_matrix_keep_world_geometry(obj: bpy.types.Object, target_matrix: Matrix):
    """Set object origin/transform to target_matrix without moving geometry in world."""
    if not obj or obj.type != "MESH":
        return

    M0 = obj.matrix_world.copy()
    M1 = target_matrix.copy()

    try:
        T = M1.inverted() @ M0
    except Exception:
        obj.matrix_world = M1
        return

    obj.data.transform(T)
    obj.data.update()
    obj.matrix_world = M1


def _join_objects(context: bpy.types.Context, objs, joined_name: str):
    if not objs:
        return None
    if len(objs) == 1:
        objs[0].name = joined_name
        return objs[0]

    try:
        if context.mode != "OBJECT":
            bpy.ops.object.mode_set(mode="OBJECT")
    except Exception:
        pass

    vl = context.view_layer
    for o in vl.objects:
        o.select_set(False)

    for o in objs:
        if o and o.name in bpy.data.objects:
            o.select_set(True)

    vl.objects.active = objs[-1]

    try:
        bpy.ops.object.join()
    except Exception as e:
        log_ui(context, f"[Outline2D] Join failed: {e}")
        return None

    joined = vl.objects.active
    if joined:
        joined.name = joined_name
    return joined


# ============================================================
# Create filled outline mesh
# ============================================================
def _create_filled_outline_mesh(
    context: bpy.types.Context,
    src_obj: bpy.types.Object,
    exterior_coords_xy,
    name: str,
    triangulate: bool,
    place_at_bbox_center: bool,
):
    if not exterior_coords_xy or len(exterior_coords_xy) < 4:
        return None

    pts = list(exterior_coords_xy)
    if pts[0] == pts[-1]:
        pts = pts[:-1]
    if len(pts) < 3:
        return None

    center = _bbox_center_world(src_obj) if place_at_bbox_center else Vector((0.0, 0.0, 0.0))

    mesh = bpy.data.meshes.new(name + "_MESH")
    _copy_material_slots(src_obj, mesh)

    obj = bpy.data.objects.new(name, mesh)
    context.collection.objects.link(obj)

    obj.matrix_world.identity()
    obj.location = center if place_at_bbox_center else (0.0, 0.0, 0.0)

    bm = bmesh.new()

    if place_at_bbox_center:
        bmv = [bm.verts.new((x - center.x, y - center.y, 0.0)) for (x, y) in pts]
    else:
        bmv = [bm.verts.new((x, y, 0.0)) for (x, y) in pts]

    bm.verts.ensure_lookup_table()

    try:
        face = bm.faces.new(bmv)
    except ValueError:
        bm.free()
        return None

    face.normal_update()
    if face.normal.z < 0:
        bmesh.ops.reverse_faces(bm, faces=[face])

    if len(mesh.materials) > 0:
        mi = min(max(int(src_obj.active_material_index), 0), len(mesh.materials) - 1)
        face.material_index = mi

    if triangulate:
        bmesh.ops.triangulate(bm, faces=[face])

    bm.to_mesh(mesh)
    bm.free()
    mesh.update()

    if place_at_bbox_center:
        obj.location.z = center.z

    return obj


# ============================================================
# Preprocess: split by material then loose
# ============================================================
def split_by_material_then_loose(context: bpy.types.Context, src: bpy.types.Object, duplicate_source: bool):
    if not src or src.type != "MESH":
        return [], None, "Object must be a mesh"

    try:
        if context.mode != "OBJECT":
            bpy.ops.object.mode_set(mode="OBJECT")
    except Exception:
        pass

    temp_root = src
    if duplicate_source:
        _set_selection(context, src, [src])
        try:
            bpy.ops.object.duplicate(linked=False)
        except Exception as e:
            return [], None, f"Duplicate failed: {e}"
        temp_root = context.view_layer.objects.active

    _set_selection(context, temp_root, [temp_root])

    # Separate by MATERIAL
    try:
        bpy.ops.object.mode_set(mode="EDIT")
        bpy.ops.mesh.separate(type="MATERIAL")
        bpy.ops.object.mode_set(mode="OBJECT")
    except Exception as e:
        try:
            bpy.ops.object.mode_set(mode="OBJECT")
        except Exception:
            pass
        return [], temp_root, f"Separate by MATERIAL failed: {e}"

    material_objs = [o for o in context.selected_objects if o and o.type == "MESH"]
    all_parts = []

    # For each material object, separate by LOOSE
    for mo in material_objs:
        _set_selection(context, mo, [mo])
        try:
            bpy.ops.object.mode_set(mode="EDIT")
            bpy.ops.mesh.separate(type="LOOSE")
            bpy.ops.object.mode_set(mode="OBJECT")
        except Exception as e:
            try:
                bpy.ops.object.mode_set(mode="OBJECT")
            except Exception:
                pass
            return [], temp_root, f"Separate by LOOSE failed on {mo.name}: {e}"

        loose_objs = [o for o in context.selected_objects if o and o.type == "MESH"]
        all_parts.extend(loose_objs)

    # Unique by name
    uniq, seen = [], set()
    for o in all_parts:
        if o.name not in seen:
            seen.add(o.name)
            uniq.append(o)

    return uniq, temp_root, None


# ============================================================
# Core: process one mesh object => filled outline
# ============================================================
def process_one_object(
    context: bpy.types.Context,
    src: bpy.types.Object,
    *,
    Polygon,
    unary_union,
    pick_largest: bool,
    triangulate: bool,
    place_at_bbox_center: bool,
    name_prefix: str,
):
    if not src or src.type != "MESH":
        return None, "Object is not a mesh"

    depsgraph = context.evaluated_depsgraph_get()
    src_eval = src.evaluated_get(depsgraph)
    mesh_eval = src_eval.to_mesh()
    if not mesh_eval:
        return None, "Failed to get evaluated mesh"

    bm = bmesh.new()
    bm.from_mesh(mesh_eval)
    bm.faces.ensure_lookup_table()

    world_mat = src.matrix_world.copy()
    polys = []

    for f in bm.faces:
        coords = _face_loop_world_xy(f, world_mat)
        if not coords:
            continue
        try:
            poly = Polygon(coords)
            if not poly.is_valid:
                poly = poly.buffer(0)
            if not poly.is_empty:
                polys.append(poly)
        except Exception:
            continue

    bm.free()
    src_eval.to_mesh_clear()

    if not polys:
        return None, "No valid faces to union"

    try:
        u = unary_union(polys)
    except Exception as e:
        return None, f"unary_union failed: {e}"

    if u.is_empty:
        return None, "Union is empty"

    geom = u
    if geom.geom_type == "MultiPolygon":
        geoms = list(geom.geoms)
        geom = max(geoms, key=lambda g: g.area) if pick_largest else geoms[0]

    if geom.geom_type != "Polygon":
        return None, f"Unexpected geometry type: {geom.geom_type}"

    exterior = list(geom.exterior.coords)
    if len(exterior) < 4:
        return None, "Exterior too small"

    out_name = f"{name_prefix}_{src.name}_Outline2D_Fill"
    out_obj = _create_filled_outline_mesh(
        context=context,
        src_obj=src,
        exterior_coords_xy=exterior,
        name=out_name,
        triangulate=triangulate,
        place_at_bbox_center=place_at_bbox_center,
    )
    if not out_obj:
        return None, "Failed to create filled outline mesh"

    return out_obj, None


# ============================================================
# Bonsai: add representation (From Object)
# ============================================================
def bonsai_add_representation_pipeline(
    context: bpy.types.Context,
    *,
    ifc_obj: bpy.types.Object,
    outline_obj: bpy.types.Object,
    context_id: str | int | None = None,
) -> bool:
    def log(m: str):
        log_ui(context, f"[Bonsai] {m}")

    if not bonsai_is_available():
        log("BonsaiBIM not available (operator missing)")
        return False

    if not ifc_obj or ifc_obj.type != "MESH":
        log("ifc_obj invalid or not MESH")
        return False
    if not outline_obj or outline_obj.type != "MESH":
        log("outline_obj invalid or not MESH")
        return False

    try:
        if bpy.context.mode != "OBJECT":
            bpy.ops.object.mode_set(mode="OBJECT")
    except Exception as e:
        log(f"mode_set OBJECT failed (non-fatal): {e}")

    scene = bpy.context.scene
    if not hasattr(scene, "BIMGeometryProperties"):
        log("scene has no BIMGeometryProperties")
        return False
    gprops = scene.BIMGeometryProperties

    try:
        gprops.representation_from_object = None
    except Exception:
        pass

    try:
        gprops.representation_from_object = outline_obj
        log(f"representation_from_object = {outline_obj.name}")
    except Exception as e:
        log(f"set representation_from_object failed: {e}")
        return False

    if context_id is not None:
        if not hasattr(ifc_obj, "BIMGeometryProperties"):
            log("ifc_obj has no BIMGeometryProperties")
            return False
        try:
            ifc_obj.BIMGeometryProperties.contexts = str(context_id)
            log(f"contexts = {ifc_obj.BIMGeometryProperties.contexts}")
        except Exception as e:
            log(f"set contexts failed: {e}")
            return False
    else:
        log("context_id not provided; using existing contexts")

    try:
        for o in bpy.context.selected_objects:
            o.select_set(False)
        ifc_obj.select_set(True)
        bpy.context.view_layer.objects.active = ifc_obj
        log(f"active_object = {ifc_obj.name}")
    except Exception as e:
        log(f"selection/active setup failed: {e}")
        return False

    try:
        bpy.ops.bim.add_representation(representation_conversion_method="OBJECT")
        log("add_representation executed")
        return True
    except Exception as e:
        log(f"operator failed: {e}")
        return False


# ============================================================
# Per-IFC pipeline
# ============================================================
def run_pipeline_for_ifc_object(
    context: bpy.types.Context,
    ifc_obj: bpy.types.Object,
    props,
    Polygon,
    unary_union,
) -> tuple[list[bpy.types.Object], list[tuple[str, str]]]:
    created: list[bpy.types.Object] = []
    failed: list[tuple[str, str]] = []

    initial_matrix = ifc_obj.matrix_world.copy()
    name_prefix = ifc_obj.name

    log_ui(context, f"[Outline2D] === Start: {ifc_obj.name} ===")

    # Determine candidates (split or not)
    if props.preprocess_split_active:
        parts, _temp_root, err = split_by_material_then_loose(
            context,
            ifc_obj,
            duplicate_source=props.split_duplicate_source,
        )
        if err:
            failed.append((ifc_obj.name, err))
            log_ui(context, f"[Outline2D] FAIL {ifc_obj.name}: {err}")
            return [], failed

        candidates = parts
        log_ui(context, f"[Outline2D] Split parts: {len(candidates)}")

        if not candidates:
            msg = "Split produced no parts."
            failed.append((ifc_obj.name, msg))
            log_ui(context, f"[Outline2D] FAIL {ifc_obj.name}: {msg}")
            return [], failed
    else:
        candidates = [ifc_obj]
        log_ui(context, "[Outline2D] No split (process object directly)")

    # Process each candidate -> outline mesh
    outlines: list[bpy.types.Object] = []
    for src in candidates:
        if props.preprocess_remove_unused_materials:
            try:
                _remove_unused_material_slots(src)
            except Exception as e:
                msg = f"Remove unused materials failed: {e}"
                failed.append((src.name, msg))
                log_ui(context, f"[Outline2D] FAIL {src.name}: {msg}")
                continue

        out_obj, msg = process_one_object(
            context,
            src,
            Polygon=Polygon,
            unary_union=unary_union,
            pick_largest=props.pick_largest,
            triangulate=props.triangulate,
            place_at_bbox_center=props.place_at_bbox_center,
            name_prefix=name_prefix,
        )
        if out_obj:
            outlines.append(out_obj)
        else:
            failed.append((src.name, msg or "Unknown error"))
            log_ui(context, f"[Outline2D] FAIL {src.name}: {msg or 'Unknown error'}")

    if not outlines:
        log_ui(context, f"[Outline2D] === End: {ifc_obj.name} (no outlines) ===")
        return [], failed

    # Join per IFC object
    if props.post_join_created and len(outlines) >= 2:
        joined = _join_objects(context, outlines, f"{ifc_obj.name}_Outline2D_Joined")
        if joined:
            if props.post_join_origin_from_initial_active:
                _set_origin_matrix_keep_world_geometry(joined, initial_matrix)
            created = [joined]
            log_ui(context, f"[Outline2D] Joined -> {joined.name}")
        else:
            created = outlines
            log_ui(context, f"[Outline2D] Join failed; kept {len(outlines)} outlines")
    else:
        created = outlines
        log_ui(context, f"[Outline2D] Created outlines: {len(created)} (no join)")

    # Delete split parts
    if props.preprocess_split_active and props.post_delete_split_parts:
        try:
            if context.mode != "OBJECT":
                bpy.ops.object.mode_set(mode="OBJECT")
        except Exception:
            pass
        try:
            vl = context.view_layer
            for o in vl.objects:
                o.select_set(False)
            for o in candidates:
                if o and o.name in bpy.data.objects:
                    o.select_set(True)
            vl.objects.active = candidates[-1]
            bpy.ops.object.delete()
            log_ui(context, "[Outline2D] Deleted split parts")
        except Exception as e:
            log_ui(context, f"[Outline2D] Delete split parts failed: {e}")

    # Bonsai add representation
    if props.bonsai_add_representation and created:
        outline_obj = created[-1]
        ctx_id = props.bonsai_context_id.strip() or None
        ok = bonsai_add_representation_pipeline(
            context,
            ifc_obj=ifc_obj,
            outline_obj=outline_obj,
            context_id=ctx_id,
        )
        if not ok:
            failed.append((ifc_obj.name, "Bonsai add_representation failed"))
            log_ui(context, f"[Outline2D] FAIL {ifc_obj.name}: Bonsai add_representation failed")
        else:
            log_ui(context, f"[Outline2D] Bonsai rep added from {outline_obj.name}")

    log_ui(context, f"[Outline2D] === End: {ifc_obj.name} (created {len(created)}) ===")
    return created, failed


# ============================================================
# UI Props
# ============================================================
class OUTLINE2D_Props(bpy.types.PropertyGroup):
    # Defaults (match your screenshot intent)
    preprocess_split_active: bpy.props.BoolProperty(
        name="Preprocess: Split (Material → Loose)",
        default=True,
    )
    split_duplicate_source: bpy.props.BoolProperty(
        name="Split: Duplicate Source",
        default=True,
    )
    post_delete_split_parts: bpy.props.BoolProperty(
        name="Post: Delete Split Parts",
        default=True,
    )

    preprocess_remove_unused_materials: bpy.props.BoolProperty(
        name="Preprocess: Remove Unused Materials",
        default=True,
    )

    place_at_bbox_center: bpy.props.BoolProperty(
        name="Place at BBox Center",
        default=True,
    )

    pick_largest: bpy.props.BoolProperty(
        name="Pick Largest Island",
        default=True,
    )

    triangulate: bpy.props.BoolProperty(
        name="Triangulate",
        default=False,
    )

    post_join_created: bpy.props.BoolProperty(
        name="Post: Join Result (per IFC)",
        default=True,
    )

    post_join_origin_from_initial_active: bpy.props.BoolProperty(
        name="Post: Join Origin = IFC Object",
        default=True,
    )

    bonsai_add_representation: bpy.props.BoolProperty(
        name="Bonsai: Add Representation (From Object)",
        default=True,
    )
    bonsai_context_id: bpy.props.StringProperty(
        name="Bonsai Context ID",
        default="",
    )

    # Logging controls
    ui_log_max_lines: bpy.props.IntProperty(
        name="Max Log Lines",
        default=250,
        min=20,
        max=5000,
    )
    ui_log_clear_on_run: bpy.props.BoolProperty(
        name="Clear Logs on Run",
        default=True,
    )


# ============================================================
# Main Operator
# ============================================================
class OUTLINE2D_OT_run(bpy.types.Operator):
    bl_idname = "outline2d.run_pipeline"
    bl_label = "Run Outline2D Pipeline"
    bl_options = {"REGISTER", "UNDO"}

    def execute(self, context):
        props = context.scene.outline2d_props

        if props.ui_log_clear_on_run:
            clear_ui_logs(context)

        log_ui(context, "[Outline2D] ===== RUN =====")

        Polygon, unary_union, err = _ensure_shapely()
        if err is not None:
            log_ui(context, f"[Outline2D] ERROR: Shapely import failed: {err}")
            self.report({"ERROR"}, f"Shapely import failed: {err}")
            return {"CANCELLED"}

        # Multi-select targets
        selected_mesh = [o for o in context.selected_objects if o and o.type == "MESH"]
        if not selected_mesh:
            if context.active_object and context.active_object.type == "MESH":
                selected_mesh = [context.active_object]
            else:
                log_ui(context, "[Outline2D] ERROR: Select at least one mesh/IFC object.")
                self.report({"ERROR"}, "Select at least one mesh/IFC object.")
                return {"CANCELLED"}

        active = context.active_object if context.active_object in selected_mesh else None
        targets = [active] + [o for o in selected_mesh if o != active] if active else selected_mesh
        log_ui(context, f"[Outline2D] Targets: {len(targets)}")

        # Bonsai auto-disable if missing
        if not bonsai_is_available():
            if props.bonsai_add_representation:
                log_ui(context, "[Outline2D] Bonsai missing -> disabling Bonsai options")
            props.bonsai_add_representation = False

        all_created: list[bpy.types.Object] = []
        all_failed: list[tuple[str, str]] = []

        # Process each IFC object
        for ifc_obj in targets:
            try:
                # Make current IFC active (helps ops)
                for o in context.selected_objects:
                    o.select_set(False)
                ifc_obj.select_set(True)
                context.view_layer.objects.active = ifc_obj
            except Exception:
                pass

            created, failed = run_pipeline_for_ifc_object(
                context,
                ifc_obj,
                props,
                Polygon,
                unary_union,
            )
            all_created.extend(created)
            all_failed.extend(failed)

        # Select created results
        if all_created:
            for o in context.view_layer.objects:
                o.select_set(False)
            for o in all_created:
                if o and o.name in bpy.data.objects:
                    o.select_set(True)
            context.view_layer.objects.active = all_created[-1]

        log_ui(context, f"[Outline2D] SUMMARY: Created {len(all_created)} | Failed {len(all_failed)}")
        if all_failed:
            for name, msg in all_failed:
                log_ui(context, f"[Outline2D] FAIL {name}: {msg}")

        if all_failed and not all_created:
            self.report({"ERROR"}, f"All failed ({len(all_failed)}). Check Logs.")
            return {"CANCELLED"}

        if all_failed:
            self.report({"WARNING"}, f"Created {len(all_created)}. Failed {len(all_failed)} (see Logs).")
        else:
            self.report({"INFO"}, f"Batch complete. Created {len(all_created)} object(s).")

        return {"FINISHED"}


# ============================================================
# Panel
# ============================================================
class OUTLINE2D_PT_panel(bpy.types.Panel):
    bl_label = "Outline2D"
    bl_idname = "OUTLINE2D_PT_panel"
    bl_space_type = "VIEW_3D"
    bl_region_type = "UI"
    bl_category = "Outline2D"

    def draw(self, context):
        layout = self.layout
        props = context.scene.outline2d_props
        bonsai_ok = bonsai_is_available()

        layout.prop(props, "preprocess_split_active")
        col = layout.column(align=True)
        col.enabled = props.preprocess_split_active
        col.prop(props, "split_duplicate_source")
        col.prop(props, "post_delete_split_parts")

        layout.separator()
        layout.prop(props, "preprocess_remove_unused_materials")
        layout.prop(props, "place_at_bbox_center")
        layout.prop(props, "pick_largest")
        layout.prop(props, "triangulate")

        layout.separator()
        layout.prop(props, "post_join_created")
        row = layout.row()
        row.enabled = props.post_join_created
        row.prop(props, "post_join_origin_from_initial_active")

        layout.separator()

        brow = layout.row()
        brow.enabled = bonsai_ok
        brow.prop(props, "bonsai_add_representation")

        bcol = layout.column(align=True)
        bcol.enabled = bonsai_ok and props.bonsai_add_representation
        bcol.prop(props, "bonsai_context_id")

        if not bonsai_ok:
            layout.label(text="BonsaiBIM not installed: Bonsai options disabled", icon="INFO")

        layout.separator()
        layout.operator("outline2d.run_pipeline", icon="FACESEL")

        layout.separator()
        layout.prop(props, "ui_log_clear_on_run")
        layout.prop(props, "ui_log_max_lines")

        box = layout.box()
        box.label(text="Logs")

        r = box.row(align=True)
        r.operator("outline2d.clear_logs", icon="TRASH", text="")
        r.operator("outline2d.copy_logs", icon="COPYDOWN", text="Copy")

        box.template_list(
            "OUTLINE2D_UL_logs",
            "",
            context.scene.outline2d_logbuf,
            "lines",
            context.scene.outline2d_logbuf,
            "active_index",
            rows=10,
        )


# ============================================================
# Register / Unregister (safe)
# ============================================================
classes = (
    OUTLINE2D_LogLine,
    OUTLINE2D_LogBuffer,
    OUTLINE2D_UL_logs,
    OUTLINE2D_Props,
    OUTLINE2D_OT_clear_logs,
    OUTLINE2D_OT_copy_logs,
    OUTLINE2D_OT_run,
    OUTLINE2D_PT_panel,
)


def register():
    # Do not “force unregister” here to avoid Blender RNA edge cases.
    # This file is intended as an addon: install+enable once.
    for c in classes:
        bpy.utils.register_class(c)

    bpy.types.Scene.outline2d_props = bpy.props.PointerProperty(type=OUTLINE2D_Props)
    bpy.types.Scene.outline2d_logbuf = bpy.props.PointerProperty(type=OUTLINE2D_LogBuffer)

    # Auto-disable Bonsai checkbox if missing
    try:
        if not bonsai_is_available():
            bpy.context.scene.outline2d_props.bonsai_add_representation = False
    except Exception:
        pass


def unregister():
    # Remove scene props first
    try:
        if hasattr(bpy.types.Scene, "outline2d_props"):
            del bpy.types.Scene.outline2d_props
    except Exception:
        pass

    try:
        if hasattr(bpy.types.Scene, "outline2d_logbuf"):
            del bpy.types.Scene.outline2d_logbuf
    except Exception:
        pass

    for c in reversed(classes):
        bpy.utils.unregister_class(c)


if __name__ == "__main__":
    register()