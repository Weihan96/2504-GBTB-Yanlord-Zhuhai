# pyright: reportInvalidTypeForm=false
bl_info = {
    "name": "CAD 2D (HLR + Material Boundaries + Fills)",
    "author": "ChatGPT",
    "version": (0, 1, 5),
    "blender": (4, 0, 0),
    "location": "View3D > Sidebar (N) > CAD2D",
    "description": "Clean CAD-like 2D outline + material seams + per-material fill mesh from Mesh using face-level HLR. Minimal UI; multi-line tooltips.",
    "category": "Object",
}

import bpy
import bmesh
from mathutils import Vector
from mathutils.geometry import tessellate_polygon
from math import isfinite, radians, hypot
from bpy.props import (
    BoolProperty,
    EnumProperty,
    FloatProperty,
    IntProperty,
    StringProperty,
    PointerProperty,
    CollectionProperty,
)
from bpy.types import Operator, Panel, PropertyGroup, UIList

# -----------------------------
# Shapely (required)
# -----------------------------
try:
    from shapely.geometry import Polygon, LineString, MultiLineString
    from shapely.ops import unary_union
    _SHAPELY_OK = True
    _SHAPELY_ERR = ""
except Exception as e:
    _SHAPELY_OK = False
    _SHAPELY_ERR = repr(e)


# ============================================================
# Logs
# ============================================================

class CAD2D_LogItem(PropertyGroup):
    level: StringProperty(name="Level", default="INFO")
    message: StringProperty(name="Message", default="")


class CAD2D_Settings(PropertyGroup):
    # --- Basic ---
    view: EnumProperty(
        name="View",
        items=[
            ("TOP", "Top", "Top (XY)"),
            ("FRONT", "Front", "Front (XZ)"),
            ("RIGHT", "Right", "Right (YZ)"),
        ],
        default="TOP",
        description=(
            "Projection View\n"
            "TOP: XY (looking -Z)\n"
            "FRONT: XZ (looking -Y)\n"
            "RIGHT: YZ (looking -X)"
        ),
    )

    # Debris cleanup (unit-aware UI via subtype DISTANCE)
    min_seg_len: FloatProperty(
        name="Min Segment Length",
        default=0.001,
        min=0.0,
        soft_max=0.1,
        precision=6,
        subtype="DISTANCE",
        description=(
            "Min Segment Length\n"
            "Drop very short polyline segments (removes debris).\n"
            "Uses scene length units."
        ),
    )

    # Material smoothing as ppm (unit-agnostic, ratio of bbox diagonal)
    mat_smooth_ppm: IntProperty(
        name="Smooth (ppm)",
        default=50,
        min=0,
        soft_max=500,
        description=(
            "Material Smooth (ppm)\n"
            "Smoothing distance = bbox_diagonal * ppm / 1,000,000\n"
            "Higher: merges tiny surface noise\n"
            "Lower: preserves small details\n"
            "Unit-agnostic (relative to model size)"
        ),
    )

    # Fill output
    generate_fills: BoolProperty(
        name="Generate Fills (Mesh)",
        default=True,
        description=(
            "Generate Fills (Mesh)\n"
            "Create 2D fill meshes for each visible material region.\n"
            "Faces keep material_index and copy material slots from source."
        ),
    )
    fill_z_offset: FloatProperty(
        name="Fill Z Offset",
        default=-0.001,
        soft_min=-0.1,
        soft_max=0.1,
        precision=6,
        subtype="DISTANCE",
        description=(
            "Fill Z Offset\n"
            "Offset fills along +Z/-Z so linework renders above fills.\n"
            "Uses scene length units."
        ),
    )
    fill_uv_scale: FloatProperty(
        name="Fill UV Scale",
        default=1.0,
        min=1e-8,
        soft_max=10.0,
        precision=4,
        description=(
            "Fill UV Scale\n"
            "UV = (x * scale, y * scale) in projected 2D coordinates.\n"
            "Adjust to control texture scale on fills."
        ),
    )

    # Output
    output_collection: StringProperty(
        name="Output Collection",
        default="CAD_2D",
        description="Collection to store generated outputs (curves + optional fills).",
    )
    unlink_from_scene_root: BoolProperty(
        name="Unlink From Scene Root",
        default=True,
        description=(
            "Unlink From Scene Root\n"
            "Keep outputs only under the output collection (clean Outliner)."
        ),
    )

    # --- Advanced (hidden by default) ---
    show_advanced: BoolProperty(
        name="Show Advanced",
        default=False,
        description="Show advanced parameters (usually not needed).",
    )

    triangulate: BoolProperty(
        name="Triangulate",
        default=True,
        description="Triangulate faces for robust polygon projection (recommended).",
    )
    per_island: BoolProperty(
        name="Per Island",
        default=True,
        description="Process connected face islands separately.",
    )

    simplify_tol: FloatProperty(
        name="Simplify Tol",
        default=1e-5,
        min=0.0,
        soft_max=1e-2,
        precision=8,
        description="Simplify tolerance for linework. Small recommended.",
    )
    snap_decimals: IntProperty(
        name="Snap Decimals",
        default=6,
        min=0,
        max=10,
        description="Round projected coordinates to this number of decimals.",
    )
    min_visible_area: FloatProperty(
        name="Min Visible Area",
        default=0.0,
        min=0.0,
        soft_max=1e-2,
        precision=10,
        description="Drop tiny visible polygon fragments (0 disables).",
    )

    detail_mode: EnumProperty(
        name="Detail Mode",
        items=[
            ("NONE", "None", "No detail lines"),
            ("MATERIAL", "Material Boundaries", "Use material visible regions to generate seams"),
            ("HARD", "Hard Edges", "Hard edges filtered by visible faces (advanced)"),
            ("BOTH", "Both", "Combine Material boundaries + Hard edges (advanced)"),
        ],
        default="MATERIAL",
        description=(
            "Detail Mode\n"
            "Material Boundaries: seams from visible material regions (recommended)\n"
            "Hard Edges: dihedral edges filtered by visibility (advanced)\n"
            "Both: combine both"
        ),
    )
    detail_angle_deg: FloatProperty(
        name="Hard Edge Angle (deg)",
        default=30.0,
        min=0.0,
        max=180.0,
        precision=1,
        description="Hard edge threshold (dihedral angle >= this).",
    )

    clip_details_to_visible: BoolProperty(
        name="Clip Details To Visible Mask",
        default=True,
        description=(
            "Clip Details To Visible Mask\n"
            "Clip detail lines by visible mask in 2D (reduces strays)."
        ),
    )
    clip_buffer: FloatProperty(
        name="Clip Buffer",
        default=0.0,
        min=0.0,
        soft_max=1e-2,
        precision=8,
        description="Optional buffer applied to visible mask before clipping details.",
    )

    # Logs
    logs: CollectionProperty(type=CAD2D_LogItem)
    log_index: IntProperty(name="Log Index", default=0)


def _st(context) -> CAD2D_Settings:
    return context.scene.cad2d_settings


def _log(context, level: str, msg: str):
    st = _st(context)
    it = st.logs.add()
    it.level = level
    it.message = msg
    st.log_index = max(0, len(st.logs) - 1)


# ============================================================
# View basis & projection
# ============================================================

def _get_view_basis(view: str):
    v = view.upper()
    if v == "TOP":
        view_dir = Vector((0, 0, -1))
        right = Vector((1, 0, 0))
        up = Vector((0, 1, 0))
    elif v == "FRONT":
        view_dir = Vector((0, -1, 0))
        right = Vector((1, 0, 0))
        up = Vector((0, 0, 1))
    elif v == "RIGHT":
        view_dir = Vector((-1, 0, 0))
        right = Vector((0, -1, 0))
        up = Vector((0, 0, 1))
    else:
        raise ValueError("View must be TOP/FRONT/RIGHT")
    return right, up, view_dir


def _project_to_2d(v: Vector, right: Vector, up: Vector, snap_decimals: int):
    x = v.dot(right)
    y = v.dot(up)
    if not (isfinite(x) and isfinite(y)):
        return None
    return (round(x, snap_decimals), round(y, snap_decimals))


def _ensure_collection(name: str):
    col = bpy.data.collections.get(name)
    if not col:
        col = bpy.data.collections.new(name)
        bpy.context.scene.collection.children.link(col)
    return col


# ============================================================
# Curve building
# ============================================================

def _add_lines_to_curve_data(curve_data: bpy.types.Curve, lines, min_seg_len: float):
    def add_polyline(coords):
        if len(coords) < 2:
            return
        filtered = [coords[0]]
        for p in coords[1:]:
            prev = filtered[-1]
            if ((p[0] - prev[0]) ** 2 + (p[1] - prev[1]) ** 2) ** 0.5 >= min_seg_len:
                filtered.append(p)
        if len(filtered) < 2:
            return
        spline = curve_data.splines.new("POLY")
        spline.points.add(len(filtered) - 1)
        for i, (x, y) in enumerate(filtered):
            spline.points[i].co = (x, y, 0.0, 1.0)

    if lines is None:
        return

    if isinstance(lines, LineString):
        add_polyline(list(lines.coords))
    elif isinstance(lines, MultiLineString):
        for g in lines.geoms:
            add_polyline(list(g.coords))
    else:
        if hasattr(lines, "geoms"):
            for g in lines.geoms:
                if isinstance(g, LineString):
                    add_polyline(list(g.coords))
                elif isinstance(g, MultiLineString):
                    for gg in g.geoms:
                        add_polyline(list(gg.coords))


def _make_curve_object(name: str, lines, min_seg_len: float):
    curve_data = bpy.data.curves.new(name, type="CURVE")
    curve_data.dimensions = "2D"
    _add_lines_to_curve_data(curve_data, lines, min_seg_len)
    return bpy.data.objects.new(name, curve_data)


# ============================================================
# Face islands
# ============================================================

def _split_face_islands(bm: bmesh.types.BMesh):
    bm.faces.ensure_lookup_table()
    visited = set()
    islands = []
    for f in bm.faces:
        if f.index in visited:
            continue
        stack = [f]
        visited.add(f.index)
        comp = []
        while stack:
            cur = stack.pop()
            comp.append(cur)
            for e in cur.edges:
                for nf in e.link_faces:
                    if nf.index not in visited:
                        visited.add(nf.index)
                        stack.append(nf)
        islands.append(comp)
    return islands


# ============================================================
# HLR island: visible mask + visible parts by material
# ============================================================

def _hlr_island(
    faces,
    right: Vector,
    up: Vector,
    view_dir: Vector,
    snap_decimals: int,
    min_visible_area: float,
):
    """
    Returns:
      visible_mask (Polygon/MultiPolygon),
      visible_faces (set[int]),
      visible_parts_by_mat (dict[int, list[geom]]),
      stats
    """
    items = []
    face_total = 0
    face_kept = 0

    for f in faces:
        face_total += 1
        if f.normal.dot(view_dir) >= 0:
            continue

        pts2d = [_project_to_2d(v.co, right, up, snap_decimals) for v in f.verts]
        if any(p is None for p in pts2d):
            continue

        poly = Polygon(pts2d)
        if poly.is_empty or (not poly.is_valid) or poly.area <= 1e-12:
            continue

        depth = f.calc_center_median().dot(view_dir)
        mat = int(getattr(f, "material_index", 0))
        items.append((depth, f.index, mat, poly))
        face_kept += 1

    if not items:
        return None, set(), {}, (face_kept, face_total, 0, 0)

    items.sort(key=lambda x: x[0])  # far -> near

    visible_region = None
    visible_parts = []
    visible_faces = set()
    visible_parts_by_mat = {}

    polys_total = len(items)
    visible_parts_count = 0

    for _, fidx, mat, poly in items:
        if visible_region is None:
            visible_part = poly
            visible_region = poly
        else:
            visible_part = poly.difference(visible_region)
            if visible_part.is_empty:
                continue
            visible_region = visible_region.union(poly)

        if visible_part.is_empty:
            continue

        if min_visible_area > 0:
            try:
                if visible_part.area < min_visible_area:
                    continue
            except Exception:
                pass

        visible_faces.add(fidx)
        visible_parts.append(visible_part)
        visible_parts_by_mat.setdefault(mat, []).append(visible_part)
        visible_parts_count += 1

    if not visible_parts:
        return None, visible_faces, visible_parts_by_mat, (face_kept, face_total, polys_total, 0)

    visible_mask = unary_union(visible_parts)
    if visible_mask.is_empty:
        return None, visible_faces, visible_parts_by_mat, (face_kept, face_total, polys_total, visible_parts_count)

    return visible_mask, visible_faces, visible_parts_by_mat, (face_kept, face_total, polys_total, visible_parts_count)


# ============================================================
# Detail: hard edges filtered by visible faces (advanced)
# ============================================================

def _detail_lines_hard_edges(
    bm: bmesh.types.BMesh,
    visible_faces: set,
    right: Vector,
    up: Vector,
    view_dir: Vector,
    snap_decimals: int,
    angle_threshold_rad: float,
):
    segs = []
    bm.edges.ensure_lookup_table()

    for e in bm.edges:
        if len(e.link_faces) != 2:
            continue
        f1, f2 = e.link_faces[0], e.link_faces[1]

        ang = f1.normal.angle(f2.normal)
        if ang < angle_threshold_rad:
            continue

        if (f1.index not in visible_faces) and (f2.index not in visible_faces):
            continue

        if (f1.normal.dot(view_dir) >= 0) and (f2.normal.dot(view_dir) >= 0):
            continue

        p0 = _project_to_2d(e.verts[0].co, right, up, snap_decimals)
        p1 = _project_to_2d(e.verts[1].co, right, up, snap_decimals)
        if p0 is None or p1 is None or p0 == p1:
            continue

        segs.append(LineString([p0, p1]))

    if not segs:
        return None
    return unary_union(segs)


# ============================================================
# Material regions (union visible parts per mat) + optional smooth (ppm)
# ============================================================

def _build_regions_by_material_ppm(
    visible_parts_by_mat_all: dict,
    smooth_ppm: int,
):
    regions = {}
    for mat, parts in visible_parts_by_mat_all.items():
        if not parts:
            continue
        try:
            reg = unary_union(parts)
        except Exception:
            continue
        if reg.is_empty:
            continue
        regions[mat] = reg

    if not regions:
        return {}

    ppm = max(0, int(smooth_ppm))
    if ppm <= 0:
        return regions

    try:
        all_union = unary_union(list(regions.values()))
        if all_union.is_empty:
            return regions
        minx, miny, maxx, maxy = all_union.bounds
        diag = hypot(maxx - minx, maxy - miny)
        eps = diag * (ppm / 1_000_000.0)
    except Exception:
        eps = 0.0

    if eps <= 0:
        return regions

    for k, reg in list(regions.items()):
        try:
            regions[k] = reg.buffer(eps).buffer(-eps)
        except Exception:
            pass

    return regions


# ============================================================
# Detail: material boundaries (from regions)
# ============================================================

def _detail_lines_material_boundaries(
    regions_by_mat: dict,
    outline_lines,
    simplify_tol: float,
):
    if not regions_by_mat:
        return None

    boundaries = []
    for reg in regions_by_mat.values():
        b = reg.boundary
        if not b.is_empty:
            boundaries.append(b)

    if not boundaries:
        return None

    seam = unary_union(boundaries)

    try:
        outline_buf = outline_lines.buffer(max(simplify_tol, 1e-10))
        seam = seam.difference(outline_buf)
    except Exception:
        pass

    if seam.is_empty:
        return None

    if simplify_tol > 0:
        try:
            seam = seam.simplify(simplify_tol)
        except Exception:
            pass

    return seam


# ============================================================
# Fill Mesh: build 2D triangles from shapely polygons (supports holes)
# ============================================================

def _iter_polygons(geom):
    if geom is None or geom.is_empty:
        return
    gtype = getattr(geom, "geom_type", "")
    if gtype == "Polygon":
        yield geom
    elif hasattr(geom, "geoms"):
        for gg in geom.geoms:
            yield from _iter_polygons(gg)


def _ring_coords_no_close(ring):
    coords = list(ring.coords)
    if len(coords) >= 2 and coords[0] == coords[-1]:
        coords = coords[:-1]
    return coords


def _build_fill_mesh_object(
    name: str,
    regions_by_mat: dict,
    src_obj: bpy.types.Object,
    z: float,
    uv_scale: float,
):
    me = bpy.data.meshes.new(name)
    obj = bpy.data.objects.new(name, me)

    # Copy material slots from source object
    me.materials.clear()
    if src_obj and src_obj.type == "MESH":
        for slot in src_obj.material_slots:
            if slot.material:
                me.materials.append(slot.material)

    bm = bmesh.new()
    uv_layer = bm.loops.layers.uv.new("UVMap")

    def add_polygon(poly: Polygon, mat_index: int):
        loops = []
        ext = _ring_coords_no_close(poly.exterior)
        if len(ext) < 3:
            return 0
        loops.append([Vector((x, y)) for (x, y) in ext])

        for interior in poly.interiors:
            hole = _ring_coords_no_close(interior)
            if len(hole) >= 3:
                loops.append([Vector((x, y)) for (x, y) in hole])

        tris = tessellate_polygon(loops)
        if not tris:
            return 0

        flat2d = []
        for loop in loops:
            flat2d.extend(loop)

        bm_verts = []
        for v2 in flat2d:
            bm_verts.append(bm.verts.new((v2.x, v2.y, z)))

        created = 0
        max_mat = max(0, len(me.materials) - 1)
        mat_index = max(0, min(int(mat_index), max_mat)) if len(me.materials) > 0 else 0

        for (i1, i2, i3) in tris:
            try:
                f = bm.faces.new((bm_verts[i1], bm_verts[i2], bm_verts[i3]))
            except ValueError:
                continue
            f.material_index = mat_index
            for loop in f.loops:
                co = loop.vert.co
                loop[uv_layer].uv = (co.x * uv_scale, co.y * uv_scale)
            created += 1
        return created

    tri_count = 0
    for mat, region in regions_by_mat.items():
        for poly in _iter_polygons(region):
            tri_count += add_polygon(poly, mat)

    bm.normal_update()
    bm.to_mesh(me)
    bm.free()

    return obj, tri_count


# ============================================================
# Generate
# ============================================================

def generate_cad2d(context, src_obj: bpy.types.Object):
    if not _SHAPELY_OK:
        raise RuntimeError(f"Shapely import failed: {_SHAPELY_ERR}")
    if not src_obj or src_obj.type != "MESH":
        raise RuntimeError("Active object must be a Mesh.")

    st = _st(context)
    right, up, view_dir = _get_view_basis(st.view)
    angle_threshold_rad = radians(max(0.0, st.detail_angle_deg))

    # Duplicate mesh (do not touch original)
    tmp = src_obj.copy()
    tmp.data = src_obj.data.copy()
    context.scene.collection.objects.link(tmp)

    bm = bmesh.new()
    bm.from_mesh(tmp.data)
    bm.verts.ensure_lookup_table()

    # Bake world transform into vertices (equivalent to apply transforms on the duplicate, operator-free)
    bm.transform(tmp.matrix_world)

    if st.triangulate:
        bm.faces.ensure_lookup_table()
        bmesh.ops.triangulate(bm, faces=bm.faces[:], quad_method="BEAUTY", ngon_method="BEAUTY")
        bm.faces.ensure_lookup_table()

    # Islands
    if st.per_island:
        islands = _split_face_islands(bm)
    else:
        bm.faces.ensure_lookup_table()
        islands = [list(bm.faces)]

    outline_boundaries = []
    union_masks = []
    union_visible_faces = set()
    visible_parts_by_mat_all = {}

    face_total_sum = 0
    face_kept_sum = 0
    polys_sum = 0
    visible_parts_sum = 0
    solved_islands = 0

    for faces in islands:
        visible_mask, visible_faces, visible_parts_by_mat, stats = _hlr_island(
            faces=faces,
            right=right,
            up=up,
            view_dir=view_dir,
            snap_decimals=st.snap_decimals,
            min_visible_area=st.min_visible_area,
        )

        union_visible_faces |= visible_faces

        for mat, parts in visible_parts_by_mat.items():
            if parts:
                visible_parts_by_mat_all.setdefault(mat, []).extend(parts)

        face_kept, face_total, polys_total, visible_parts = stats
        face_total_sum += face_total
        face_kept_sum += face_kept
        polys_sum += polys_total
        visible_parts_sum += visible_parts

        if visible_mask is None:
            continue

        solved_islands += 1
        union_masks.append(visible_mask)

        bnd = visible_mask.boundary
        if not bnd.is_empty:
            outline_boundaries.append(bnd)

    if not outline_boundaries:
        bm.free()
        bpy.data.objects.remove(tmp, do_unlink=True)
        raise RuntimeError("No outline extracted. Try different view or check normals/scale.")

    outline = unary_union(outline_boundaries)
    if st.simplify_tol > 0:
        outline = outline.simplify(st.simplify_tol)

    visible_union_mask = unary_union(union_masks) if union_masks else None

    # Regions by material (ppm smoothing)
    regions_by_mat = _build_regions_by_material_ppm(
        visible_parts_by_mat_all=visible_parts_by_mat_all,
        smooth_ppm=st.mat_smooth_ppm,
    )

    # Detail
    details = None
    if st.detail_mode in {"MATERIAL", "BOTH"}:
        details = _detail_lines_material_boundaries(
            regions_by_mat=regions_by_mat,
            outline_lines=outline,
            simplify_tol=st.simplify_tol,
        )

    if st.detail_mode in {"HARD", "BOTH"}:
        hard = _detail_lines_hard_edges(
            bm=bm,
            visible_faces=union_visible_faces,
            right=right,
            up=up,
            view_dir=view_dir,
            snap_decimals=st.snap_decimals,
            angle_threshold_rad=angle_threshold_rad,
        )
        if hard is not None and (not hard.is_empty):
            details = hard if details is None else unary_union([details, hard])

    # Clip details to visible mask
    if details is not None and (not details.is_empty) and visible_union_mask is not None and st.clip_details_to_visible:
        mask = visible_union_mask
        if st.clip_buffer > 0:
            try:
                mask = mask.buffer(st.clip_buffer)
            except Exception:
                pass
        try:
            details = details.intersection(mask)
        except Exception:
            pass

        if st.simplify_tol > 0:
            try:
                details = details.simplify(st.simplify_tol)
            except Exception:
                pass

    # Cleanup temp
    bm.free()
    bpy.data.objects.remove(tmp, do_unlink=True)

    # Outputs
    base = f"CAD_2D_{src_obj.name}_{st.view}"
    outline_name = f"{base}_OUTLINE"
    detail_name = f"{base}_DETAIL"
    fill_name = f"{base}_FILL"

    outline_obj = _make_curve_object(outline_name, outline, st.min_seg_len)

    detail_obj = None
    if details is not None and (not getattr(details, "is_empty", False)):
        detail_obj = _make_curve_object(detail_name, details, st.min_seg_len)

    fill_obj = None
    fill_tris = 0
    if st.generate_fills and regions_by_mat:
        fill_obj, fill_tris = _build_fill_mesh_object(
            name=fill_name,
            regions_by_mat=regions_by_mat,
            src_obj=src_obj,
            z=st.fill_z_offset,
            uv_scale=st.fill_uv_scale,
        )

    col = _ensure_collection(st.output_collection)
    col.objects.link(outline_obj)
    if detail_obj:
        col.objects.link(detail_obj)
    if fill_obj:
        col.objects.link(fill_obj)

    if st.unlink_from_scene_root:
        sc = bpy.context.scene.collection
        if outline_obj.name in sc.objects:
            sc.objects.unlink(outline_obj)
        if detail_obj and detail_obj.name in sc.objects:
            sc.objects.unlink(detail_obj)
        if fill_obj and fill_obj.name in sc.objects:
            sc.objects.unlink(fill_obj)

    return {
        "outline_obj": outline_obj,
        "detail_obj": detail_obj,
        "fill_obj": fill_obj,
        "fill_tris": fill_tris,
        "face_kept": face_kept_sum,
        "face_total": face_total_sum,
        "polys_total": polys_sum,
        "visible_parts": visible_parts_sum,
        "islands_total": len(islands),
        "islands_solved": solved_islands,
        "materials_seen": len(regions_by_mat),
    }


# ============================================================
# Operators
# ============================================================

class CAD2D_OT_Generate(Operator):
    bl_idname = "cad2d.generate"
    bl_label = "Generate CAD 2D"
    bl_options = {"REGISTER", "UNDO"}

    def execute(self, context):
        st = _st(context)
        obj = context.view_layer.objects.active

        _log(context, "INFO", "v0.1.5 Start")
        _log(context, "INFO", f"View={st.view} Detail={st.detail_mode} Smooth={st.mat_smooth_ppm}ppm MinSeg={st.min_seg_len}")
        _log(context, "INFO", f"Fill={st.generate_fills} Z={st.fill_z_offset} UVScale={st.fill_uv_scale}")

        if not obj:
            _log(context, "ERROR", "No active object.")
            self.report({"ERROR"}, "No active object.")
            return {"CANCELLED"}

        try:
            stats = generate_cad2d(context, obj)
            _log(context, "INFO", f"OK | Outline: {stats['outline_obj'].name}")
            _log(context, "INFO", f"OK | Detail: {stats['detail_obj'].name if stats['detail_obj'] else '(none)'}")
            _log(context, "INFO", f"OK | Fill: {stats['fill_obj'].name if stats['fill_obj'] else '(none)'} tris={stats['fill_tris']}")
            _log(context, "INFO", f"Islands: {stats['islands_solved']}/{stats['islands_total']} Mats: {stats['materials_seen']}")
            self.report({"INFO"}, "CAD2D generated.")
            return {"FINISHED"}
        except Exception as e:
            msg = f"{type(e).__name__}: {e}"
            _log(context, "ERROR", msg)
            self.report({"ERROR"}, msg)
            return {"CANCELLED"}


class CAD2D_OT_ClearLogs(Operator):
    bl_idname = "cad2d.clear_logs"
    bl_label = "Clear Logs"
    bl_options = {"REGISTER"}

    def execute(self, context):
        st = _st(context)
        st.logs.clear()
        st.log_index = 0
        return {"FINISHED"}


class CAD2D_OT_DeleteOutputs(Operator):
    bl_idname = "cad2d.delete_outputs"
    bl_label = "Delete Outputs"
    bl_options = {"REGISTER", "UNDO"}

    def execute(self, context):
        st = _st(context)
        col = bpy.data.collections.get(st.output_collection)
        if not col:
            _log(context, "INFO", f"No collection '{st.output_collection}'.")
            return {"FINISHED"}

        def is_out(o: bpy.types.Object):
            if not o:
                return False
            if not o.name.startswith("CAD_2D_"):
                return False
            return o.name.endswith("_OUTLINE") or o.name.endswith("_DETAIL") or o.name.endswith("_FILL")

        to_delete = [o for o in list(col.objects) if is_out(o)]
        for o in to_delete:
            bpy.data.objects.remove(o, do_unlink=True)

        _log(context, "INFO", f"Deleted {len(to_delete)} outputs.")
        return {"FINISHED"}


# ============================================================
# UIList + Panel
# ============================================================

class CAD2D_UL_Logs(UIList):
    def draw_item(self, context, layout, data, item, icon, active_data, active_propname, index):
        lvl = getattr(item, "level", "INFO")
        msg = getattr(item, "message", "")
        icon_map = {"INFO": "INFO", "WARN": "ERROR", "ERROR": "CANCEL"}
        row = layout.row(align=True)
        row.label(text=lvl, icon=icon_map.get(lvl, "DOT"))
        row.label(text=msg)


class CAD2D_PT_Panel(Panel):
    bl_label = "CAD 2D"
    bl_idname = "CAD2D_PT_panel"
    bl_space_type = "VIEW_3D"
    bl_region_type = "UI"
    bl_category = "CAD2D"

    def draw(self, context):
        layout = self.layout
        st = _st(context)

        if not _SHAPELY_OK:
            box = layout.box()
            box.label(text="Shapely not available", icon="ERROR")
            box.label(text=_SHAPELY_ERR[:160])

        # Basic
        layout.label(text="Basic", icon="PREFERENCES")
        layout.prop(st, "view")

        # Cleanup
        layout.separator()
        layout.label(text="Cleanup", icon="MODIFIER")
        layout.prop(st, "min_seg_len")

        # Detail
        layout.separator()
        layout.label(text="Detail", icon="MATERIAL")
        layout.prop(st, "detail_mode")
        layout.prop(st, "mat_smooth_ppm")

        # Fill
        layout.separator()
        layout.label(text="Fill", icon="MESH_GRID")
        layout.prop(st, "generate_fills")
        col = layout.column()
        col.enabled = st.generate_fills
        col.prop(st, "fill_z_offset")
        col.prop(st, "fill_uv_scale")

        # Output
        layout.separator()
        layout.label(text="Output", icon="OUTLINER_COLLECTION")
        layout.prop(st, "output_collection")
        layout.prop(st, "unlink_from_scene_root")

        layout.separator()
        row = layout.row(align=True)
        row.operator("cad2d.generate", icon="OUTLINER_OB_CURVE")
        row.operator("cad2d.delete_outputs", icon="TRASH")

        # Advanced
        layout.separator()
        layout.prop(st, "show_advanced", toggle=True)
        if st.show_advanced:
            adv = layout.box()
            adv.label(text="Advanced", icon="TOOL_SETTINGS")
            adv.prop(st, "triangulate")
            adv.prop(st, "per_island")
            adv.prop(st, "simplify_tol")
            adv.prop(st, "snap_decimals")
            adv.prop(st, "min_visible_area")
            adv.prop(st, "clip_details_to_visible")
            sub = adv.column()
            sub.enabled = st.clip_details_to_visible
            sub.prop(st, "clip_buffer")

            adv.separator()
            adv.label(text="Hard Edges", icon="EDGESEL")
            adv.prop(st, "detail_angle_deg")

        # Logs
        layout.separator()
        layout.label(text="Logs", icon="TEXT")
        layout.operator("cad2d.clear_logs", icon="X")
        layout.template_list("CAD2D_UL_Logs", "", st, "logs", st, "log_index", rows=6)


# ============================================================
# Register
# ============================================================

classes = (
    CAD2D_LogItem,
    CAD2D_Settings,
    CAD2D_OT_Generate,
    CAD2D_OT_ClearLogs,
    CAD2D_OT_DeleteOutputs,
    CAD2D_UL_Logs,
    CAD2D_PT_Panel,
)

def register():
    for c in classes:
        bpy.utils.register_class(c)
    bpy.types.Scene.cad2d_settings = PointerProperty(type=CAD2D_Settings)

def unregister():
    if hasattr(bpy.types.Scene, "cad2d_settings"):
        del bpy.types.Scene.cad2d_settings
    for c in reversed(classes):
        bpy.utils.unregister_class(c)

if __name__ == "__main__":
    try:
        unregister()
    except Exception:
        pass
    register()