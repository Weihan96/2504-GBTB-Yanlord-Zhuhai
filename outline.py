# pyright: reportInvalidTypeForm=false

bl_info = {
    "name": "Outline Mesh (IfcOpenShell-style, with Faces & Materials)",
    "author": "ChatGPT",
    "version": (0, 1, 7),
    "blender": (4, 0, 0),
    "location": "View3D > Sidebar (N) > Silhouette",
    "category": "Object",
}

import bpy
import bmesh
from math import radians
from mathutils import Vector
from typing import Literal


Axis = Literal["+X", "-X", "+Y", "-Y", "+Z", "-Z"]


# ----------------------------
# Materials
# ----------------------------

def _copy_material_slots(src_obj: bpy.types.Object, dst_mesh: bpy.types.Mesh):
    """
    Copy material slots from src_obj to dst_mesh.
    IMPORTANT: Must exist BEFORE bmesh.to_mesh if you want material_index preserved.
    """
    dst_mesh.materials.clear()
    for slot in src_obj.material_slots:
        dst_mesh.materials.append(slot.material if slot.material else None)


# ----------------------------
# Core: visible faces -> copy -> flatten
# ----------------------------

def generate_outline_mesh_with_faces(obj: bpy.types.Object, axis: Axis = "+Z") -> bpy.types.Mesh:
    def get_visible_faces(obj: bpy.types.Object, bm: bmesh.types.BMesh, axis: Axis) -> list[bmesh.types.BMFace]:
        """
        Visible face heuristic (IfcOpenShell-style):
        - reject faces whose normal points along the ray direction (back-facing to the "view")
        - raycast from a point on a far plane toward the object along the ray direction
        - if first hit is this face (or not this object), mark visible
        """
        distance = max(obj.dimensions.xyz)
        if distance <= 0:
            distance = 1.0

        eps = 0.002
        depsgraph = bpy.context.evaluated_depsgraph_get()

        # direction_local = ray direction in OBJ LOCAL space
        # far_val = coordinate on the far plane (also in OBJ LOCAL space, using obj.bound_box)
        if axis == "+Z":
            far_val = max(co[2] for co in obj.bound_box) + eps
            direction_local = Vector((0, 0, -1))
        elif axis == "-Z":
            far_val = min(co[2] for co in obj.bound_box) - eps
            direction_local = Vector((0, 0, 1))
        elif axis == "+Y":
            far_val = max(co[1] for co in obj.bound_box) + eps
            direction_local = Vector((0, -1, 0))
        elif axis == "-Y":
            far_val = min(co[1] for co in obj.bound_box) - eps
            direction_local = Vector((0, 1, 0))
        elif axis == "+X":
            far_val = max(co[0] for co in obj.bound_box) + eps
            direction_local = Vector((-1, 0, 0))
        else:  # "-X"
            far_val = min(co[0] for co in obj.bound_box) - eps
            direction_local = Vector((1, 0, 0))

        global_direction = obj.matrix_world.to_quaternion() @ direction_local

        visible_faces: list[bmesh.types.BMFace] = []
        for face in bm.faces:
            # If face normal points along ray direction, it's "facing away" from the view side -> skip
            if direction_local.dot(face.normal) > 0:
                continue

            c = face.calc_center_median()

            if axis in {"+Z", "-Z"}:
                start_local = Vector((c.x, c.y, far_val))
            elif axis in {"+Y", "-Y"}:
                start_local = Vector((c.x, far_val, c.z))
            else:  # {"+X","-X"}
                start_local = Vector((far_val, c.y, c.z))

            start_world = obj.matrix_world @ start_local

            hit, loc, norm, idx, o, mw = bpy.context.scene.ray_cast(
                depsgraph, start_world, global_direction, distance=distance
            )

            if o != obj or idx == face.index:
                visible_faces.append(face)

        return visible_faces

    def get_contour_edges(visible_faces: list[bmesh.types.BMFace]) -> list[bmesh.types.BMEdge]:
        contour_edges: list[bmesh.types.BMEdge] = []
        vis_set = set(visible_faces)
        for face in visible_faces:
            for edge in face.edges:
                total_linked_faces = len(edge.link_faces)
                if total_linked_faces == 1:
                    contour_edges.append(edge)
                elif total_linked_faces == 2:
                    other_face = edge.link_faces[0] if edge.link_faces[1] == face else edge.link_faces[1]
                    if other_face not in vis_set:
                        contour_edges.append(edge)
        return contour_edges

    def get_crease_edges(visible_faces: list[bmesh.types.BMFace], threshold: float) -> list[bmesh.types.BMEdge]:
        crease_edges: list[bmesh.types.BMEdge] = []
        for face in visible_faces:
            for edge in face.edges:
                if len(edge.link_faces) == 2:
                    angle = edge.link_faces[0].normal.angle(edge.link_faces[1].normal)
                    if abs(angle) > threshold:
                        crease_edges.append(edge)
        return crease_edges

    # --- Source bmesh ---
    bm = bmesh.new()
    bm.from_mesh(obj.data)
    bm.faces.ensure_lookup_table()

    visible_faces = get_visible_faces(obj, bm, axis=axis)
    if not visible_faces:
        bm.free()
        raise RuntimeError("No visible faces found.")

    # keep the outline computation (same as before; currently not used downstream)
    outline_edges = set(get_contour_edges(visible_faces))
    outline_edges.update(get_crease_edges(visible_faces, radians(60)))
    _ = outline_edges

    # --- Copy VISIBLE FACES to bm_new, preserve material_index ---
    bm_new = bmesh.new()
    vert_map: dict[bmesh.types.BMVert, bmesh.types.BMVert] = {}

    max_mi = max(len(obj.material_slots) - 1, 0)

    for f in visible_faces:
        new_face_verts = []
        for v in f.verts:
            nv = vert_map.get(v)
            if nv is None:
                nv = bm_new.verts.new(v.co.copy())
                vert_map[v] = nv
            new_face_verts.append(nv)

        try:
            nf = bm_new.faces.new(new_face_verts)
        except ValueError:
            nf = None

        if nf is not None:
            mi = f.material_index
            if mi < 0:
                mi = 0
            elif mi > max_mi:
                mi = max_mi
            nf.material_index = mi

    bm_new.verts.ensure_lookup_table()

    # --- Flatten along axis in new bmesh ---
    for vert in bm_new.verts:
        if axis in {"+Z", "-Z"}:
            vert.co.z = 0.0
        elif axis in {"+Y", "-Y"}:
            vert.co.y = 0.0
        else:  # {"+X","-X"}
            vert.co.x = 0.0

    # CRITICAL FIX: copy materials BEFORE to_mesh
    new_mesh = bpy.data.meshes.new("outline_tmp")
    _copy_material_slots(obj, new_mesh)
    bm_new.to_mesh(new_mesh)

    bm_new.free()
    bm.free()

    new_mesh.update()
    return new_mesh


# ----------------------------
# Operator / UI
# ----------------------------

class OBJECT_OT_generate_outline_mesh(bpy.types.Operator):
    bl_idname = "object.generate_outline_mesh_ifc_faces"
    bl_label = "Generate Outline Mesh (Keep Faces+Mats)"
    bl_options = {"REGISTER", "UNDO"}

    def execute(self, context: bpy.types.Context):
        src = context.view_layer.objects.active
        if not src or src.type != "MESH":
            self.report({"ERROR"}, "Please select an active mesh object.")
            return {"CANCELLED"}

        axis: Axis = context.scene.sil_axis_dropdown  # type: ignore[assignment]

        tmp = src # TODO add duplicate and apply rotation+scale if needed in the future, currently we just use the original mesh as is for raycasting and vertex positions, so that we can preserve edit mode selection and avoid depsgraph issues with evaluated meshes. This means the operator will work best on objects with no rotation and scale of 1, but it should still produce correct results even if those transforms are present (just not perfectly flattened along the axis).
        try:
            # ALWAYS apply rotation+scale on a temp copy

            new_mesh = generate_outline_mesh_with_faces(tmp, axis=axis)

            suffix = axis.replace("+", "p").replace("-", "m")
            new_obj = bpy.data.objects.new(f"{src.name}_OUTLINE_FACES_{suffix}", new_mesh)

            col = src.users_collection[0] if src.users_collection else context.scene.collection
            col.objects.link(new_obj)

            # Use the temp object's matrix_world (after apply rot+scale)
            new_obj.matrix_world = tmp.matrix_world.copy()

        except Exception as e:
            self.report({"ERROR"}, f"Generate failed: {e}")
            return {"CANCELLED"}

        bpy.ops.object.select_all(action="DESELECT")
        new_obj.select_set(True)
        context.view_layer.objects.active = new_obj

        self.report({"INFO"}, f"Outline (faces+mats) created: {new_obj.name}")
        return {"FINISHED"}


class VIEW3D_PT_outline_mesh_panel(bpy.types.Panel):
    bl_label = "Silhouette"
    bl_idname = "VIEW3D_PT_outline_mesh_panel"
    bl_space_type = "VIEW_3D"
    bl_region_type = "UI"
    bl_category = "Silhouette"

    def draw(self, context):
        layout = self.layout
        layout.label(text="IfcOpenShell-style Outline")
        layout.prop(context.scene, "sil_axis_dropdown", text="Axis")
        layout.operator("object.generate_outline_mesh_ifc_faces", icon="MESH_DATA")


classes = (
    OBJECT_OT_generate_outline_mesh,
    VIEW3D_PT_outline_mesh_panel,
)


def register():
    for c in classes:
        bpy.utils.register_class(c)

    # N-panel dropdown (6 options)
    bpy.types.Scene.sil_axis_dropdown = bpy.props.EnumProperty(
        name="Axis",
        items=[
            ("+X", "+X", "Project from +X side (ray -X), flatten X=0"),
            ("-X", "-X", "Project from -X side (ray +X), flatten X=0"),
            ("+Y", "+Y", "Project from +Y side (ray -Y), flatten Y=0"),
            ("-Y", "-Y", "Project from -Y side (ray +Y), flatten Y=0"),
            ("+Z", "+Z", "Project from +Z side (ray -Z), flatten Z=0"),
            ("-Z", "-Z", "Project from -Z side (ray +Z), flatten Z=0"),
        ],
        default="+Z",
    )


def unregister():
    if hasattr(bpy.types.Scene, "sil_axis_dropdown"):
        delattr(bpy.types.Scene, "sil_axis_dropdown")

    for c in reversed(classes):
        bpy.utils.unregister_class(c)


if __name__ == "__main__":
    register()
