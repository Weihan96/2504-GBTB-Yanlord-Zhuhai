# pyright: reportInvalidTypeForm=false

bl_info = {
    "name": "Silhouette Core: Flatten + Outer Boundary Rebuild (World Axis)",
    "author": "ChatGPT",
    "version": (0, 2, 1),
    "blender": (4, 0, 0),
    "location": "View3D > Sidebar (N) > Silhouette",
    "category": "Object",
}

import bpy
import bmesh
import math
from mathutils import Vector
from typing import Literal, List, Tuple, Dict


Vec2 = Tuple[float, float]


# ----------------------------
# Geometry helpers
# ----------------------------

def bbox_center_world(obj: bpy.types.Object) -> Vector:
    corners = [obj.matrix_world @ Vector(c) for c in obj.bound_box]
    minv = Vector((min(c.x for c in corners), min(c.y for c in corners), min(c.z for c in corners)))
    maxv = Vector((max(c.x for c in corners), max(c.y for c in corners), max(c.z for c in corners)))
    return (minv + maxv) * 0.5


def copy_materials(src_obj: bpy.types.Object, dst_obj: bpy.types.Object):
    dst_obj.data.materials.clear()
    for slot in src_obj.material_slots:
        if slot.material:
            dst_obj.data.materials.append(slot.material)


def make_projection_body(
    context: bpy.types.Context,
    src: bpy.types.Object,
    axis: Literal["X", "Y", "Z"],
    apply_modifiers: bool = True,
    apply_scale: bool = True,
) -> bpy.types.Object:
    """
    Duplicate src -> (optionally) bake modifiers on the duplicate -> apply scale -> return proj object.
    Materials synced from src.
    """
    if axis not in {"X", "Y", "Z"}:
        raise ValueError("axis must be 'X','Y','Z'")

    # Create a lightweight duplicate object (we may replace its mesh with evaluated mesh)
    proj = src.copy()
    proj.data = src.data.copy()
    proj.animation_data_clear()
    proj.name = f"{src.name}_PROJ_{axis}"

    col = src.users_collection[0] if src.users_collection else context.scene.collection
    col.objects.link(proj)

    # Sync materials
    copy_materials(src, proj)

    # Bake modifiers ON proj (not on src)
    if apply_modifiers and len(proj.modifiers) > 0:
        depsgraph = context.evaluated_depsgraph_get()
        eval_obj = proj.evaluated_get(depsgraph)
        # new_from_object bakes modifiers into a new mesh datablock
        baked_mesh = bpy.data.meshes.new_from_object(
            eval_obj,
            preserve_all_data_layers=True,
            depsgraph=depsgraph,
        )
        old_mesh = proj.data
        proj.data = baked_mesh
        try:
            bpy.data.meshes.remove(old_mesh)
        except Exception:
            pass

        # Keep materials after mesh replacement
        copy_materials(src, proj)

        # Clear modifiers to avoid double effects later
        proj.modifiers.clear()

    # Apply scale if requested (on proj only)
    if apply_scale:
        if context.mode != "OBJECT":
            bpy.ops.object.mode_set(mode="OBJECT")
        bpy.ops.object.select_all(action="DESELECT")
        proj.select_set(True)
        context.view_layer.objects.active = proj
        bpy.ops.object.transform_apply(location=False, rotation=False, scale=True)

    return proj


def flatten_to_bbox_center_plane_world(obj: bpy.types.Object, axis: Literal["X", "Y", "Z"]):
    """
    Flatten obj vertices onto world axis plane at bbox center:
      axis="Z": z = center.z
      axis="Y": y = center.y
      axis="X": x = center.x
    """
    center = bbox_center_world(obj)

    target_val = center.x if axis == "X" else center.y if axis == "Y" else center.z

    mw = obj.matrix_world
    imw = mw.inverted_safe()
    me = obj.data

    for v in me.vertices:
        pw = mw @ v.co
        if axis == "X":
            pw.x = target_val
        elif axis == "Y":
            pw.y = target_val
        else:
            pw.z = target_val
        v.co = imw @ pw

    me.update()


# ----------------------------
# Cleanup after flatten (kept, but NOT used in crash-safe outline)
# ----------------------------

def cleanup_after_flatten(
    obj: bpy.types.Object,
    dist_ratio: float = 1e-6,
    area_ratio: float = 1e-12,
):
    """
    Run in EDIT mode on obj.
    - merge-by-distance (remove doubles)
    - remove degenerate geometry (zero-area faces, etc.)
    Ratios are relative to bbox diagonal/area so it works in cm/m scales.

    NOTE: This can be a crash trigger on some degenerate flattened meshes.
    It's kept for reference but NOT called by the new outline selector.
    """
    if bpy.context.mode != "EDIT_MESH":
        bpy.ops.object.mode_set(mode="EDIT")

    bm = bmesh.from_edit_mesh(obj.data)
    bm.verts.ensure_lookup_table()
    bm.edges.ensure_lookup_table()
    bm.faces.ensure_lookup_table()

    minv = Vector((1e30, 1e30, 1e30))
    maxv = Vector((-1e30, -1e30, -1e30))
    for v in bm.verts:
        co = v.co
        minv.x = min(minv.x, co.x); minv.y = min(minv.y, co.y); minv.z = min(minv.z, co.z)
        maxv.x = max(maxv.x, co.x); maxv.y = max(maxv.y, co.y); maxv.z = max(maxv.z, co.z)
    diag = (maxv - minv).length
    if diag <= 0.0:
        diag = 1.0

    dist = max(diag * dist_ratio, 1e-12)
    area_eps = max((diag * diag) * area_ratio, 1e-18)

    try:
        bmesh.ops.remove_doubles(bm, verts=bm.verts, dist=dist)
    except Exception:
        pass

    bm.verts.ensure_lookup_table()
    bm.edges.ensure_lookup_table()
    bm.faces.ensure_lookup_table()

    tiny_faces = [f for f in bm.faces if f.calc_area() <= area_eps]
    if tiny_faces:
        bmesh.ops.delete(bm, geom=tiny_faces, context="FACES")

    try:
        bmesh.ops.dissolve_degenerate(bm, edges=bm.edges, dist=dist)
    except Exception:
        pass

    bmesh.update_edit_mesh(obj.data, loop_triangles=False, destructive=True)


# ----------------------------
# 2D hull utilities (pure python)
# ----------------------------

def _project_world_to_2d(pw: Vector, axis: Literal["X", "Y", "Z"]) -> Vec2:
    # axis=Z -> (x,y), axis=Y -> (x,z), axis=X -> (y,z)
    if axis == "Z":
        return (pw.x, pw.y)
    if axis == "Y":
        return (pw.x, pw.z)
    return (pw.y, pw.z)


def _cross(o: Vec2, a: Vec2, b: Vec2) -> float:
    return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])


def convex_hull_2d(points: List[Vec2]) -> List[Vec2]:
    """
    Monotonic chain convex hull.
    Return hull vertices in CCW order without repeating the first point.
    """
    pts = sorted(set(points))
    if len(pts) <= 2:
        return pts

    lower: List[Vec2] = []
    for p in pts:
        while len(lower) >= 2 and _cross(lower[-2], lower[-1], p) <= 0:
            lower.pop()
        lower.append(p)

    upper: List[Vec2] = []
    for p in reversed(pts):
        while len(upper) >= 2 and _cross(upper[-2], upper[-1], p) <= 0:
            upper.pop()
        upper.append(p)

    hull = lower[:-1] + upper[:-1]
    return hull


def _ensure_vgroup(obj: bpy.types.Object, name: str) -> bpy.types.VertexGroup:
    vg = obj.vertex_groups.get(name)
    if vg is None:
        vg = obj.vertex_groups.new(name=name)
    return vg


def _dist_point_to_seg(p: Vec2, a: Vec2, b: Vec2) -> float:
    # distance from point p to segment ab (2D)
    ax, ay = a
    bx, by = b
    px, py = p
    abx = bx - ax
    aby = by - ay
    apx = px - ax
    apy = py - ay
    denom = abx * abx + aby * aby
    if denom <= 1e-30:
        dx = px - ax
        dy = py - ay
        return math.sqrt(dx * dx + dy * dy)
    t = (apx * abx + apy * aby) / denom
    t = 0.0 if t < 0.0 else 1.0 if t > 1.0 else t
    cx = ax + t * abx
    cy = ay + t * aby
    dx = px - cx
    dy = py - cy
    return math.sqrt(dx * dx + dy * dy)


# ----------------------------
# Outer boundary selection (CRASH-SAFE)
#   - no cleanup
#   - no delete
#   - no new verts/edges
#   - only select EXISTING edges near hull boundary
# ----------------------------

def select_outer_boundary_edges_to_vgroup(
    obj: bpy.types.Object,
    axis: Literal["X", "Y", "Z"] = "Z",
    group_name: str = "__SIL_OUTER__",
    keep_selected_in_edit: bool = True,
) -> bool:
    """
    Crash-safe outline selector:
      1) project existing verts to 2D
      2) compute convex hull (2D)
      3) mark existing verts close to hull boundary
      4) select existing edges whose endpoints are marked and whose midpoint lies near hull boundary
      5) write marked verts into vertex group (debug)

    IMPORTANT:
      - Does NOT call cleanup_after_flatten()
      - Does NOT delete or create geometry
      - Uses destructive=False update
    """
    axis = axis.upper()
    if axis not in {"X", "Y", "Z"}:
        raise ValueError("axis must be 'X','Y','Z'")

    if bpy.context.mode != "EDIT_MESH":
        bpy.ops.object.mode_set(mode="EDIT")

    bm = bmesh.from_edit_mesh(obj.data)
    bm.verts.ensure_lookup_table()
    bm.edges.ensure_lookup_table()

    mw = obj.matrix_world

    # 1) project verts
    v2p: Dict[int, Vec2] = {}
    pts_raw: List[Vec2] = []

    min_u = 1e30
    min_v = 1e30
    max_u = -1e30
    max_v = -1e30

    for v in bm.verts:
        pw = mw @ v.co
        u, w = _project_world_to_2d(pw, axis)
        if not (math.isfinite(u) and math.isfinite(w)):
            continue
        v2p[v.index] = (u, w)
        pts_raw.append((u, w))
        min_u = min(min_u, u); min_v = min(min_v, w)
        max_u = max(max_u, u); max_v = max(max_v, w)

    if len(pts_raw) < 3:
        return False

    span = max(max_u - min_u, max_v - min_v)
    if span <= 0.0:
        return False

    # 2) quantize unique points (stability)
    q = max(span * 1e-5, 1e-6)
    uniq: Dict[Tuple[int, int], Vec2] = {}
    for (u, v) in pts_raw:
        ku = int(round(u / q))
        kv = int(round(v / q))
        uniq[(ku, kv)] = (u, v)

    uniq_pts = list(uniq.values())
    if len(uniq_pts) < 3:
        return False

    # 3) hull
    hull = convex_hull_2d(uniq_pts)
    if len(hull) < 3:
        return False

    hull_segs = []
    for i in range(len(hull)):
        a = hull[i]
        b = hull[(i + 1) % len(hull)]
        hull_segs.append((a, b))

    eps = q * 3.0

    # 4) mark existing verts near hull boundary
    hull_vert_idx = set()
    for vidx, p in v2p.items():
        dmin = 1e30
        for (a, b) in hull_segs:
            d = _dist_point_to_seg(p, a, b)
            if d < dmin:
                dmin = d
        if dmin <= eps:
            hull_vert_idx.add(vidx)

    if len(hull_vert_idx) < 3:
        return False

    # 5) select edges near hull boundary
    chosen_edges = []
    for e in bm.edges:
        i0 = e.verts[0].index
        i1 = e.verts[1].index
        if i0 not in hull_vert_idx or i1 not in hull_vert_idx:
            continue

        p0 = v2p.get(i0)
        p1 = v2p.get(i1)
        if p0 is None or p1 is None:
            continue

        mid = ((p0[0] + p1[0]) * 0.5, (p0[1] + p1[1]) * 0.5)

        dmin = 1e30
        for (a, b) in hull_segs:
            d = _dist_point_to_seg(mid, a, b)
            if d < dmin:
                dmin = d

        if dmin <= eps:
            chosen_edges.append(e)

    if not chosen_edges:
        return False

    for e in bm.edges:
        e.select = False
    for e in chosen_edges:
        e.select = True

    # 6) write marked verts into vgroup (deform layer)
    vg = _ensure_vgroup(obj, group_name)
    deform = bm.verts.layers.deform.verify()
    gid = vg.index

    # clear old weights in this group
    for v in bm.verts:
        if gid in v[deform]:
            del v[deform][gid]

    for v in bm.verts:
        if v.index in hull_vert_idx:
            v[deform][gid] = 1.0

    bmesh.update_edit_mesh(obj.data, loop_triangles=False, destructive=False)

    if keep_selected_in_edit:
        bpy.ops.mesh.select_mode(type="EDGE")
    else:
        bpy.ops.object.mode_set(mode="OBJECT")

    return True


# ----------------------------
# Rebuild clean face from selected boundary
# ----------------------------

def rebuild_clean_face_from_selected_boundary(obj: bpy.types.Object):
    """
    Assumes in EDIT mode and outer boundary EDGES are selected.
    Steps:
      1) duplicate selected edges
      2) fill (edge_face_add)
      3) invert selection -> delete old flattened projection faces
      4) delete loose
    """
    if bpy.context.mode != "EDIT_MESH":
        bpy.ops.object.mode_set(mode="EDIT")

    bpy.ops.mesh.select_mode(type="EDGE")

    # 1) duplicate boundary edges
    bpy.ops.mesh.duplicate()

    # 2) fill to create face
    bpy.ops.mesh.edge_face_add()  # 'F'

    # 3) delete old projection faces
    bpy.ops.mesh.select_mode(type="FACE")
    bpy.ops.mesh.select_all(action="INVERT")
    bpy.ops.mesh.delete(type="FACE")

    # 4) cleanup loose edges/verts
    bpy.ops.mesh.select_all(action="SELECT")
    bpy.ops.mesh.delete_loose(use_verts=True, use_edges=True, use_faces=False)
    bpy.ops.mesh.select_all(action="DESELECT")


# ----------------------------
# Operator / UI
# ----------------------------

class OBJECT_OT_silhouette_inner_core(bpy.types.Operator):
    bl_idname = "object.silhouette_inner_core"
    bl_label = "Silhouette Core (One Axis)"
    bl_options = {"REGISTER", "UNDO"}

    axis: bpy.props.EnumProperty(
        name="Axis",
        items=[("X", "X", "World X"), ("Y", "Y", "World Y"), ("Z", "Z", "World Z")],
        default="Z",
    )
    apply_modifiers: bpy.props.BoolProperty(
        name="Apply Modifiers on Projection",
        default=True,
    )
    apply_scale: bpy.props.BoolProperty(
        name="Apply Scale on Projection",
        default=True,
    )

    def execute(self, context):
        src = context.view_layer.objects.active
        if not src or src.type != "MESH":
            self.report({"ERROR"}, "Please set an active mesh object.")
            return {"CANCELLED"}

        if context.mode != "OBJECT":
            bpy.ops.object.mode_set(mode="OBJECT")

        # 1) projection body (dup + modifiers baked + scale applied)
        proj = make_projection_body(
            context,
            src,
            self.axis,
            apply_modifiers=self.apply_modifiers,
            apply_scale=self.apply_scale,
        )

        # 2) flatten onto bbox-center world-axis plane
        flatten_to_bbox_center_plane_world(proj, self.axis)

        # 3) select outer boundary + rebuild clean face
        bpy.ops.object.select_all(action="DESELECT")
        proj.select_set(True)
        context.view_layer.objects.active = proj
        bpy.ops.object.mode_set(mode="EDIT")

        ok = select_outer_boundary_edges_to_vgroup(proj, axis=self.axis, group_name="__SIL_OUTER__")
        if not ok:
            bpy.ops.object.mode_set(mode="OBJECT")
            self.report({"WARNING"}, "Hull selection failed (not enough points / degenerate / no matching edges).")
            return {"CANCELLED"}

        # rebuild_clean_face_from_selected_boundary(proj)

        bpy.ops.object.mode_set(mode="OBJECT")
        self.report({"INFO"}, f"Done: {src.name} -> {proj.name} ({self.axis})")
        return {"FINISHED"}


class VIEW3D_PT_silhouette_core_panel(bpy.types.Panel):
    bl_label = "Silhouette (Core)"
    bl_idname = "VIEW3D_PT_silhouette_core_panel"
    bl_space_type = "VIEW_3D"
    bl_region_type = "UI"
    bl_category = "Silhouette"

    def draw(self, context):
        layout = self.layout
        layout.operator("object.silhouette_inner_core", icon="MOD_SOLIDIFY")


classes = (
    OBJECT_OT_silhouette_inner_core,
    VIEW3D_PT_silhouette_core_panel,
)


def register():
    for c in classes:
        bpy.utils.register_class(c)


def unregister():
    for c in reversed(classes):
        bpy.utils.unregister_class(c)


if __name__ == "__main__":
    register()
