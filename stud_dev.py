import bpy
import bmesh
import math
from mathutils import Vector, Matrix
from mathutils.kdtree import KDTree
import ifcopenshell
from bonsai.bim.ifc import IfcStore
import ifcopenshell.util.element as elem_util
# =================================================================
#  常量：龙骨规格 & 几何参数
# =================================================================

STUD_EDGE_WIDTH = 0.034       # 边骨 34mm
STUD_MAIN_WIDTH = 0.028       # 主骨宽度 28mm（短边方向）
STUD_SEC_WIDTH = 0.049        # 副骨 49mm（长度不影响计算，只在扣 1mm 时用）
STUD_MAIN_EXCESS = 0.100      # 主骨余量 100mm
STUD_ALIGN_OFFSET = 0.0395    # 原点偏移 39.5mm
STUD_SEC_CUT = 0.001          # 副骨两端各扣 0.5mm，总共 1mm


# =================================================================
#  Log 工具
# =================================================================

def stud_log_set(context, s):
    context.scene.stud_dev_props.log = s


def stud_log_append(context, s):
    log = context.scene.stud_dev_props.log
    context.scene.stud_dev_props.log = f"{log}\n{s}" if log else s


# =================================================================
#  获取 IFC
# =================================================================

def get_ifc_model():
    try:
        return IfcStore.get_file()
    except:
        return None


def update_stud_type_enum(self, context):
    """从 IFC 模型动态加载 IfcMemberType"""
    model = get_ifc_model()
    if not model:
        return [('NONE', 'No IFC Loaded', '')]
    return [
        (t.GlobalId, f"{t.Name or '(Unnamed)'} ({t.GlobalId})", "")
        for t in model.by_type("IfcMemberType")
    ]


def find_member_type(model, guid):
    if not guid or guid == "NONE":
        return None
    for t in model.by_type("IfcMemberType"):
        if t.GlobalId == guid:
            return t
    return None


# =================================================================
#  Mesh 主轴检测 + 向量对齐
# =================================================================

def detect_mesh_axis(obj):
    bb = [Vector(c) for c in obj.bound_box]
    xs = [v.x for v in bb]
    ys = [v.y for v in bb]
    zs = [v.z for v in bb]

    len_x = max(xs) - min(xs)
    len_y = max(ys) - min(ys)
    len_z = max(zs) - min(zs)

    if len_x >= len_y and len_x >= len_z:
        return Vector((1, 0, 0))
    elif len_y >= len_x and len_y >= len_z:
        return Vector((0, 1, 0))
    else:
        return Vector((0, 0, 1))


def rotation_from_vector_to_vector(a: Vector, b: Vector):
    a = a.normalized()
    b = b.normalized()

    axis = a.cross(b)
    if axis.length < 1e-6:
        return Matrix.Identity(3) if a.dot(b) > 0 else Matrix.Rotation(
            3.14159265, 3, Vector((1, 0, 0))
        )
    return Matrix.Rotation(a.angle(b), 3, axis)


# =================================================================
#  判断：是不是 Profile（MaterialProfileSet）
# =================================================================

def is_profile_based_type(type_obj):
    """
    多数真实场景中，IfcMaterialProfileSet 是最可靠的 Profile 类型判定方式
    """
    matset = elem_util.get_material(type_obj)
    if (
        matset
        and matset.is_a("IfcMaterialProfileSet")
        and hasattr(matset, "MaterialProfiles")
        and len(matset.MaterialProfiles) > 0
        and matset.MaterialProfiles[0].Profile is not None
    ):
        return True
    return False


# =================================================================
#  姿态矩阵构造（Profile 专用：local Z = extrusion）
# =================================================================

def calc_profile_transform(start: Vector, end: Vector, roll_rad: float):
    direction = end - start
    length = direction.length
    if length < 1e-6:
        raise ValueError("两点太近")

    z_axis = direction.normalized()

    # 默认世界 Z 为 up
    world_up = Vector((0, 0, 1))
    if abs(z_axis.dot(world_up)) > 0.999:
        world_up = Vector((0, 1, 0))

    x_axis = world_up.cross(z_axis).normalized()
    y_axis = z_axis.cross(x_axis).normalized()

    # 绕挤出轴做 Roll
    R_roll = Matrix.Rotation(roll_rad, 4, z_axis)

    x_axis = (R_roll @ x_axis).normalized()
    y_axis = (R_roll @ y_axis).normalized()

    rot = Matrix((x_axis, y_axis, z_axis)).transposed()

    mat = rot.to_4x4()
    mat.translation = start
    return mat, length


# =================================================================
#  姿态矩阵构造（Mesh 专用：用 mesh 主轴对齐）
# =================================================================

def calc_mesh_transform(mesh_obj, start: Vector, end: Vector, roll_rad: float):
    target_dir = (end - start).normalized()
    mesh_axis = detect_mesh_axis(mesh_obj)

    R_align = rotation_from_vector_to_vector(mesh_axis, target_dir)
    R_roll = Matrix.Rotation(roll_rad, 4, target_dir)

    mat = R_roll @ R_align.to_4x4()
    mat.translation = start
    return mat


# =================================================================
#  ★ 新增：add_ifc_array（含内部私有 find_array_owner）
# =================================================================

def add_ifc_array(obj, axis_world: Vector, spacing: float, count: int, context):
    """
    在指定对象 obj 上创建 IFC Array 阵列。
    使用世界坐标 axis_world 方向，以 spacing 为间距，生成 count 个实例。
    不创建新的 IfcProduct，只对 obj 本身添加 IFC Array。
    """

    if count <= 1:
        return

    # ==============================================
    # 内部函数：找到最新 IFC Array 控制对象（私有）
    # ==============================================
    def _find_array_owner(_obj):
        ao = bpy.context.active_object
        if ao and hasattr(ao, "BIMArrayProperties"):
            return ao

        if hasattr(_obj, "BIMArrayProperties"):
            return _obj

        for child in _obj.children:
            if hasattr(child, "BIMArrayProperties"):
                return child

        if hasattr(_obj, "BIMObjectProperties"):
            iid = _obj.BIMObjectProperties.ifc_definition_id
            for other in bpy.data.objects:
                if (
                    hasattr(other, "BIMObjectProperties")
                    and other.BIMObjectProperties.ifc_definition_id == iid
                    and hasattr(other, "BIMArrayProperties")
                ):
                    return other

        return None

    # ==============================================
    # 1. 创建 IFC Array
    # ==============================================
    try:
        bpy.ops.bim.add_array()
    except Exception as e:
        stud_log_append(context, f"❌ add_array 失败: {e}")
        return

    # ==============================================
    # 2. 获取最新 Array 控制对象
    # ==============================================
    arr_owner = _find_array_owner(obj)
    if not arr_owner:
        stud_log_append(context, "❌ 找不到 IFC Array 控制对象")
        return

    # ==============================================
    # 3. 启用编辑 IFC Array
    # ==============================================
    try:
        bpy.ops.bim.enable_editing_array(item=-1)
    except Exception:
        pass

    arr = arr_owner.BIMArrayProperties
    axis = axis_world.normalized()

    arr.x = axis.x * spacing
    arr.y = axis.y * spacing
    arr.z = axis.z * spacing

    arr.count = count
    arr.use_local_space = False
    arr.sync_children = True

    try:
        bpy.ops.bim.edit_array(item=-1)
    except Exception:
        pass

    stud_log_append(
        context,
        f"✔ IFC Array: spacing={spacing:.4f}, count={count}, axis={axis}"
    )


# =================================================================
#  通用检查：参考面 mesh 是否合法
# =================================================================

def validate_reference_mesh(obj):
    """
    多面 Mesh 的合法性判断：
    - 允许垂直（法向 dot ≈ 0）
    - 允许平行（法向 dot ≈ ±1）
    - 其他角度一律不允许
    """

    if obj.type != "MESH":
        return False, "参考对象必须是 Mesh"

    mesh = obj.data

    if len(mesh.polygons) == 0:
        return False, "参考 Mesh 没有 polygon"

    # 收集世界空间法向量
    normals = []
    for poly in mesh.polygons:
        n = obj.matrix_world.to_3x3() @ poly.normal
        n.normalize()
        normals.append(n)

    # 单面永远合法
    if len(normals) == 1:
        return True, ""

    # 多面：允许以下情况：
    #   dot ≈ 0     → 垂直
    #   dot ≈ 1/-1  → 平行
    #   其他情况    → 非法
    for i in range(len(normals)):
        for j in range(i + 1, len(normals)):
            dot = normals[i].dot(normals[j])
            if abs(dot) < 1e-4:
                # 垂直 → 合法
                continue
            if abs(abs(dot) - 1.0) < 1e-4:
                # 平行 or 反向平行 → 合法
                continue

            # 其余角度非法
            return (
                False,
                f"参考 Mesh 中 polygon {i} 与 {j} 的法向夹角非法，dot={dot:.4f}（必须垂直或平行）"
            )

    return True, ""


# =================================================================
#  创建实例（最终调用逻辑完全一致）
# =================================================================

def create_stud_instance(context, model, type_obj, start, end, roll_rad):
    try:
        bpy.ops.bim.add_occurrence(
            relating_type_id=type_obj.id(),
            from_invoke=False,
            representation_template="EXTRUSION",
        )
    except Exception as e:
        stud_log_append(context, f"❌ add_occurrence 失败: {e}")
        return None

    obj = bpy.context.active_object
    if not obj:
        stud_log_append(context, "⚠ add_occurrence 后未找到 active_object")
        return None

    if is_profile_based_type(type_obj):
        mat, length = calc_profile_transform(start, end, roll_rad)
        try:
            bpy.ops.bim.change_profile_depth(depth=length)
            stud_log_append(context, f"✔ Profile 挤出深度 = {length:.3f}")
        except Exception as e:
            stud_log_append(context, f"⚠ 挤出深度更新失败: {e}")
    else:
        # =================================================================
        # Mesh 类型：姿态 + 自动轴向 IFC Array 拼接
        # =================================================================
        stud_log_append(context, "✔ Mesh 类型")

        # 先做姿态变换
        mat = calc_mesh_transform(obj, start, end, roll_rad)
        obj.matrix_world = mat

        # ---------------------------------------------------------------
        # 🔧 新增：自动轴向拼接（使用 add_ifc_array）
        # ---------------------------------------------------------------
        target_vec = end - start
        target_len = target_vec.length
        if target_len > 1e-6:

            mesh_axis_local = detect_mesh_axis(obj)
            bb = [Vector(c) for c in obj.bound_box]
            proj = [v.dot(mesh_axis_local) for v in bb]
            unit_len = max(proj) - min(proj)

            if unit_len > 1e-6:
                count = int(math.ceil(target_len / unit_len))
                if count > 1:
                    axis_world = target_vec.normalized()
                    spacing = unit_len
                    add_ifc_array(obj, axis_world, spacing, count, context)
                    stud_log_append(
                        context,
                        f"✔ 自动轴向拼接：unit={unit_len:.4f}, target={target_len:.4f}, count={count}"
                    )

    obj.matrix_world = mat
    return obj


# =================================================================
#  获取参考面顶点（按顺时针排序）
# =================================================================

def get_ordered_face_vertices(ref_obj):
    """返回按顺时针排序的参考面顶点（世界坐标）"""

    mesh = ref_obj.data
    verts = [ref_obj.matrix_world @ v.co for v in mesh.vertices]

    if not verts or len(verts) < 3:
        return []

    # 计算中心点
    center = Vector((0, 0, 0))
    for v in verts:
        center += v
    center /= len(verts)

    # 排序：按 atan2
    ordered = sorted(
        verts,
        key=lambda p: math.atan2((p - center).y, (p - center).x)
    )

    return ordered


# =================================================================
#  获取参考面长边和短边方向（局部坐标）
# =================================================================

def get_long_short_axis_in_local(ref_obj):
    """
    根据参考面的局部 bound_box 判断：
        local_long_axis   —— 长边方向的局部单位向量
        local_short_axis  —— 短边方向的局部单位向量

    此函数仅负责“轴归类”，不负责 extrusion 起止点。
    """

    bb = ref_obj.bound_box
    xs = [co[0] for co in bb]
    ys = [co[1] for co in bb]
    zs = [co[2] for co in bb]

    len_x = max(xs) - min(xs)
    len_y = max(ys) - min(ys)
    len_z = max(zs) - min(zs)

    lengths = {"x": len_x, "y": len_y, "z": len_z}
    axes = ["x", "y", "z"]

    # 厚度轴 = 最短的轴
    thickness_axis = min(axes, key=lambda a: lengths[a])

    # 平面轴 = 其余两个
    plane_axes = [a for a in axes if a != thickness_axis]
    a1, a2 = plane_axes

    # 长边 / 短边
    if lengths[a1] >= lengths[a2]:
        long_axis = a1
        short_axis = a2
    else:
        long_axis = a2
        short_axis = a1

    # 映射为向量
    axis_to_vec = {
        "x": Vector((1, 0, 0)),
        "y": Vector((0, 1, 0)),
        "z": Vector((0, 0, 1)),
    }

    local_long_axis = axis_to_vec[long_axis]
    local_short_axis = axis_to_vec[short_axis]

    return local_long_axis, local_short_axis


def create_offset_object_from_ref(ref_obj, offset_dist=0.02, epsilon=1e-5):
    """
    基于 ref_obj 创建 offset 后的 mesh：
        1. Solidify
        2. 删除与 ref_obj 重合顶点（旧壳）
        3. Flip normals（统一外向）
        4. 清理孤立面与边
    """
    # ---------------------------------------
    # 1. Duplicate object
    # ---------------------------------------
    offset_obj = ref_obj.copy()
    offset_obj.data = ref_obj.data.copy()
    ref_obj.users_collection[0].objects.link(offset_obj)

    # ---------------------------------------
    # 2. Solidify modifier
    # ---------------------------------------
    # 修正方向：始终向内偏移
    actual_thickness = -abs(offset_dist)

    mod = offset_obj.modifiers.new("OffsetTemp", "SOLIDIFY")
    mod.thickness = actual_thickness
    mod.offset = 1.0
    mod.use_even_offset = True

    bpy.context.view_layer.objects.active = offset_obj
    bpy.ops.object.modifier_apply(modifier=mod.name)

    # ---------------------------------------
    # 3. Build KDTree from ref_obj verts
    # ---------------------------------------
    ref_mesh = ref_obj.data
    size = len(ref_mesh.vertices)
    kd = KDTree(size)

    for i, v in enumerate(ref_mesh.vertices):
        world_co = ref_obj.matrix_world @ v.co
        kd.insert(world_co, i)

    kd.balance()

    # ---------------------------------------
    # 4. BMesh: 删除与 ref_obj 重合的顶点（旧壳）
    # ---------------------------------------
    bm = bmesh.new()
    bm.from_mesh(offset_obj.data)

    verts_to_delete = []

    for v in bm.verts:
        world_v = offset_obj.matrix_world @ v.co

        co, index, dist = kd.find(world_v)
        if dist < epsilon:
            verts_to_delete.append(v)

    bmesh.ops.delete(bm, geom=verts_to_delete, context='VERTS')

    # ---------------------------------------
    # 5. Flip normals（统一翻面）
    # ---------------------------------------
    for f in bm.faces:
        f.normal_flip()

    # ---------------------------------------
    # 6. 清理孤立面与无效 edge
    # ---------------------------------------
    invalid_faces = [f for f in bm.faces if not f.is_valid]
    if invalid_faces:
        bmesh.ops.delete(bm, geom=invalid_faces, context='FACES')

    invalid_edges = [e for e in bm.edges if not e.is_valid]
    if invalid_edges:
        bmesh.ops.delete(bm, geom=invalid_edges, context='EDGES')

    # ---------------------------------------
    # 7. 输出结果
    # ---------------------------------------
    bm.to_mesh(offset_obj.data)
    bm.free()

    return offset_obj


# =================================================================
#  参考面分析：厚度轴 / 平面轴 / 长边 / 短边
# =================================================================
#
#   2D 平面简图（忽略厚度轴）：
#
#       long_axis
#     <------------>
#   +------------------+
#   |                  |
#   |                  |  short_axis
#   |                  |
#   +------------------+
#
# 厚度轴 thickness_axis = bound_box 最小尺寸方向（通常是墙厚/板厚）
# 平面轴 plane_axes = 其余两个轴
# long_axis / short_axis 由平面两轴的长度大小决定
#
# 约定：
#   - 主龙骨 extrusion：沿 long_axis
#   - 副龙骨 extrusion：沿 short_axis

def analyse_reference_panel(ref_obj):
    """
    解析参考面（Mesh）的局部包围盒 bound_box，自动识别坐标轴意义，用于龙骨排布。

    返回一个 tuple（严格的顺序）：

        (
            local_main_start,        # 主龙骨 extrusion 起点（沿 long_axis）
            local_main_end,          # 主龙骨 extrusion 终点

            local_short_axis_vec,    # 短边方向（主龙骨排布方向）
            short_length,            # 短边长度

            local_sec_start,         # 副龙骨 extrusion 起点（沿 short_axis）
            local_sec_end,           # 副龙骨 extrusion 终点

            local_long_axis_vec,     # 长边方向（副龙骨排布方向）
            long_length,             # 长边长度
        )
    """

    bb = ref_obj.bound_box
    xs = [co[0] for co in bb]
    ys = [co[1] for co in bb]
    zs = [co[2] for co in bb]

    min_x, max_x = min(xs), max(xs)
    min_y, max_y = min(ys), max(ys)
    min_z, max_z = min(zs), max(zs)

    len_x = max_x - min_x
    len_y = max_y - min_y
    len_z = max_z - min_z

    lengths = {"x": len_x, "y": len_y, "z": len_z}
    mins = {"x": min_x, "y": min_y, "z": min_z}
    maxs = {"x": max_x, "y": max_y, "z": max_z}
    centers = {
        "x": 0.5 * (min_x + max_x),
        "y": 0.5 * (min_y + max_y),
        "z": 0.5 * (min_z + max_z),
    }

    local_long_axis_vec, local_short_axis_vec = get_long_short_axis_in_local(ref_obj)

    # 反查 axis 名称（保持你原本代码兼容）
    vec_to_axis = {
        (1,0,0): "x",
        (0,1,0): "y",
        (0,0,1): "z",
    }
    long_axis = vec_to_axis[tuple(local_long_axis_vec)]
    short_axis = vec_to_axis[tuple(local_short_axis_vec)]
    # =============================================================

    # 主龙骨 extrusion start/end（沿长边）
    coords_start = dict(centers)
    coords_end = dict(centers)

    # 厚度轴 = 非长非短的轴
    thickness_axis = [a for a in ["x","y","z"] if a not in (long_axis, short_axis)][0]

    coords_start[thickness_axis] = centers[thickness_axis]
    coords_end[thickness_axis] = centers[thickness_axis]

    coords_start[short_axis] = mins[short_axis]
    coords_end[short_axis] = mins[short_axis]

    coords_start[long_axis] = mins[long_axis]
    coords_end[long_axis] = maxs[long_axis]

    local_main_start = Vector((coords_start["x"], coords_start["y"], coords_start["z"]))
    local_main_end   = Vector((coords_end["x"],   coords_end["y"],   coords_end["z"]))

    # 副龙骨 extrusion start/end（沿短边）
    sec_coords_start = dict(centers)
    sec_coords_end   = dict(centers)

    sec_coords_start[thickness_axis] = centers[thickness_axis]
    sec_coords_end[thickness_axis]   = centers[thickness_axis]

    sec_coords_start[long_axis] = mins[long_axis]
    sec_coords_end[long_axis]   = mins[long_axis]

    sec_coords_start[short_axis] = mins[short_axis]
    sec_coords_end[short_axis]   = maxs[short_axis]

    local_sec_start = Vector((sec_coords_start["x"], sec_coords_start["y"], sec_coords_start["z"]))
    local_sec_end   = Vector((sec_coords_end["x"],   sec_coords_end["y"],   sec_coords_end["z"]))


    return (
        local_main_start,
        local_main_end,
        local_short_axis_vec,
        lengths[short_axis],
        local_sec_start,
        local_sec_end,
        local_long_axis_vec,
        lengths[long_axis],
    )


# =================================================================
#  IFC 阵列：根据传入 count 执行 IFC Array 排布（几何逻辑已全部外移）
# =================================================================

def array_studs_on_reference(
    context,
    model,
    type_obj,
    ref_obj,
    local_start,
    local_end,
    local_offset,
    local_axis_vec,
    count,
    spacing,
    roll_rad=0.0,
):
    """
    使用 IFC Array 在参考面上沿指定方向排布龙骨。

    ⭐ 本函数仅负责阵列，不参与数量计算或偏移计算。

    参数说明（按人类习惯排序）：
    -------------------------------------------------------------
    context           Blender 上下文
    model             IfcOpenShell 模型
    type_obj          IfcMemberType 对象
    ref_obj           参考面（Mesh），用于 matrix_world 变换

    local_start       基准龙骨的局部起点
    local_end         基准龙骨的局部终点
    local_offset      基于参考面的局部偏移（local）

    local_axis_vec    阵列方向（local），需为单位向量

    count             阵列数量（由 Operator 预先计算）
    spacing           阵列间距（Operator 提供）

    roll_rad          基准龙骨的旋转角度（默认为 0）
    -------------------------------------------------------------
    """

    if count <= 0:
        stud_log_append(context, "❌ count 必须 > 0")
        return

    mw = ref_obj.matrix_world

    # 1️⃣ 创建第一根龙骨（基准对象）
    world_start = mw @ (local_start + local_offset)
    world_end   = mw @ (local_end   + local_offset)

    base_obj = create_stud_instance(
        context, model, type_obj,
        world_start, world_end, roll_rad
    )
    if not base_obj:
        stud_log_append(context, "❌ 创建基准龙骨失败")
        return

    # =================================================================
    # 🚀 新逻辑：使用 add_ifc_array 完成 IFC 阵列
    # =================================================================
    axis_world = (mw.to_3x3() @ local_axis_vec).normalized()

    try:
        add_ifc_array(
            base_obj,
            axis_world,
            spacing,
            count,
            context
        )
    except Exception as e:
        stud_log_append(context, f"❌ IFC 阵列失败: {e}")
        return

    stud_log_append(context, f"🎉 IFC Array 完成，共 {count} 根")
    return base_obj
# =================================================================
#  描边：沿参考面四周生成龙骨（基于 mesh 顶点）
# =================================================================

def outline_studs_on_reference(
    context,
    model,
    type_obj,
    ref_obj,
    local_offset,
    roll_rad,
):
    """沿参考面顶点顺序生成描边龙骨"""

    verts = get_ordered_face_vertices(ref_obj)
    if len(verts) < 3:
        stud_log_append(context, "⚠ 参考面顶点不足 3 个，无法描边")
        return

    mw = ref_obj.matrix_world
    inv_mw = mw.inverted()

    # 按顺序连接：v1→v2, v2→v3, ..., vn→v1
    count = len(verts)
    studs = []
    for i in range(count):
        world_start = verts[i]
        world_end   = verts[(i + 1) % count]

        local_start = inv_mw @ world_start + local_offset
        local_end   = inv_mw @ world_end   + local_offset

        studs.append(create_stud_instance(
            context,
            model,
            type_obj,
            mw @ local_start,
            mw @ local_end,
            roll_rad,
        ))

    stud_log_append(context, f"✔ 描边龙骨已生成，共 {count} 条")
    return studs

# =================================================================
#  布局计算：主骨 / 副骨数量 & 偏移 & 挤出长度
# =================================================================

def compute_stud_layout(
    short_length: float,
    long_length: float,
    spacing: float,
    sec_spacing: float,
    local_main_start: Vector,
    local_main_end: Vector,
    local_sec_start: Vector,
    local_sec_end: Vector,
):
    """
    根据参考面尺寸和龙骨规则，计算：
      - 主骨 / 副骨数量
      - 主骨挤出长度
      - 短边 / 长边方向偏移
      - 副骨扣减 1mm 后的新起止点
      - 主骨 extrusion 新终点
    """

    # 计算副骨数量（沿长边）
    base_len_sec = STUD_SEC_WIDTH + 2 * STUD_MAIN_EXCESS + 2 * STUD_EDGE_WIDTH
    if long_length <= base_len_sec:
        sec_count = 1
    else:
        sec_count = int((long_length - base_len_sec) // sec_spacing) + 1

    # 主骨长度（沿长边 extrusion）
    main_extrude_len = (sec_count - 1) * sec_spacing + STUD_MAIN_EXCESS * 2 + STUD_MAIN_EXCESS

    # 主龙骨数量（沿短边）
    base_len_main = STUD_MAIN_WIDTH + 2 * STUD_EDGE_WIDTH
    if short_length <= base_len_main:
        main_count = 1
    else:
        main_count = int((short_length - base_len_main) // spacing) + 1

    # 主骨短边方向居中偏移
    main_short_offset = STUD_EDGE_WIDTH + (
        short_length
        - 2 * STUD_EDGE_WIDTH
        - ((main_count - 1) * spacing + STUD_MAIN_WIDTH)
    ) / 2

    # 主骨长边方向偏移
    main_long_offset = (
        STUD_EDGE_WIDTH
        + (long_length - 2 * STUD_EDGE_WIDTH - main_extrude_len) / 2
        + STUD_ALIGN_OFFSET
    )

    # 副骨长边方向偏移
    sec_long_offset = main_long_offset + STUD_MAIN_EXCESS

    # 副骨：扣除 1mm
    sec_vec = local_sec_end - local_sec_start
    sec_dir = sec_vec.normalized()
    sec_len = sec_vec.length
    new_sec_len = sec_len - STUD_SEC_CUT
    shrink = (sec_len - new_sec_len) / 2

    adjusted_sec_start = local_sec_start + sec_dir * shrink
    adjusted_sec_end   = local_sec_end   - sec_dir * shrink

    # 主骨 extrusion 重设为 main_extrude_len
    main_dir = (local_main_end - local_main_start).normalized()
    adjusted_main_start = local_main_start
    adjusted_main_end   = local_main_start + main_dir * main_extrude_len

    return {
        "sec_count": sec_count,
        "main_count": main_count,
        "main_extrude_len": main_extrude_len,
        "main_short_offset": main_short_offset,
        "main_long_offset": main_long_offset,
        "sec_long_offset": sec_long_offset,
        "adjusted_sec_start": adjusted_sec_start,
        "adjusted_sec_end": adjusted_sec_end,
        "adjusted_main_start": adjusted_main_start,
        "adjusted_main_end": adjusted_main_end,
        "sec_len_original": sec_len,
        "sec_len_new": new_sec_len,
    }


# =================================================================
#  在 canonical 面（法向 +Z，XY 为面内轴）上生成所有龙骨
# =================================================================

def generate_studs_on_canonical_panel(context, model, props, ref_panel):
    """
    在 canonical 面（ref_panel）上生成龙骨阵列。

    要求 ref_panel 满足：
        - 已 canonical 化（法向 = +Z）
        - local X/Y 为面内两个正交方向
        - scale 已应用
        - local 空间中直接可用于几何分析

    生成内容：
        - 边龙骨 outline
        - 主龙骨阵列（沿 long axis）
        - 副龙骨阵列（沿 short axis）

    返回：
        - 返回所有生成的 stud 对象（local / world 均可）
        - 这些对象将在 generate_studs_on_mesh 中被外部应用矩阵变换 T
    """
    generated_studs = []

    def _append(studs):
        if isinstance(studs, list):
            generated_studs.extend(studs)
        elif studs:
            generated_studs.append(studs)

    # 1. 主龙骨类型
    type_obj = find_member_type(model, props.selected_type)
    if not type_obj:
        stud_log_set(context, "❌ 未选择主龙骨类型")
        return generated_studs

    # 2. 分析 canonical 面的几何
    try:
        (
            local_main_start,
            local_main_end,
            local_short_axis_vec,
            short_length,
            local_sec_start,
            local_sec_end,
            local_long_axis_vec,
            long_length,
        ) = analyse_reference_panel(ref_panel)
    except Exception as e:
        stud_log_set(context, f"❌ 参考面分析失败: {e}")
        return generated_studs

    # 3. 布局逻辑
    layout = compute_stud_layout(
        short_length=short_length,
        long_length=long_length,
        spacing=props.spacing,
        sec_spacing=props.secondary_spacing,
        local_main_start=local_main_start,
        local_main_end=local_main_end,
        local_sec_start=local_sec_start,
        local_sec_end=local_sec_end,
    )

    sec_count           = layout["sec_count"]
    main_count          = layout["main_count"]
    main_short_offset   = layout["main_short_offset"]
    main_long_offset    = layout["main_long_offset"]
    sec_long_offset     = layout["sec_long_offset"]
    adjusted_sec_start  = layout["adjusted_sec_start"]
    adjusted_sec_end    = layout["adjusted_sec_end"]
    adjusted_main_start = layout["adjusted_main_start"]
    adjusted_main_end   = layout["adjusted_main_end"]
    sec_len_original    = layout["sec_len_original"]
    sec_len_new         = layout["sec_len_new"]

    # 日志
    stud_log_append(context, f"ℹ 副龙骨数量 = {sec_count}")
    stud_log_append(context, f"ℹ 主龙骨数量 = {main_count}")
    stud_log_append(context, f"ℹ 副龙骨调整: {sec_len_original:.4f} → {sec_len_new:.4f}")

    # 4. 偏移
    local_main_offset = (
        local_short_axis_vec * main_short_offset +
        local_long_axis_vec * main_long_offset +
        Vector((props.offset_x, props.offset_y, props.offset_z))
    )

    local_sec_offset = (
        local_long_axis_vec * sec_long_offset +
        Vector((
            props.secondary_offset_x,
            props.secondary_offset_y,
            props.secondary_offset_z,
        ))
    )

    # 5. 边龙骨
    edge_type_obj = find_member_type(model, props.edge_type)
    if edge_type_obj:
        edge_offset = Vector((
            props.edge_offset_x,
            props.edge_offset_y,
            props.edge_offset_z,
        ))
        studs = outline_studs_on_reference(
            context, model, edge_type_obj, ref_panel,
            edge_offset, props.edge_roll_rad
        )
        _append(studs)

    # 6. 主龙骨阵列
    base_main_stud = array_studs_on_reference(
        context, model, type_obj, ref_panel,
        adjusted_main_start, adjusted_main_end,
        local_main_offset,
        local_short_axis_vec,
        main_count, props.spacing, props.roll_rad
    )
    _append(base_main_stud)

    # 7. 副龙骨阵列
    secondary_type_obj = find_member_type(model, props.secondary_type)
    if secondary_type_obj:
        base_sec_stud = array_studs_on_reference(
            context, model, secondary_type_obj, ref_panel,
            adjusted_sec_start, adjusted_sec_end,
            local_sec_offset,
            local_long_axis_vec,
            sec_count, props.secondary_spacing, props.secondary_roll_rad
        )
        _append(base_sec_stud)

    # ======================================================
    # 8. 生成完成后，绕参考面中心 X 轴旋转 180°（翻到另一侧）
    # ======================================================
    if generated_studs and ref_panel and ref_panel.type == "MESH":
        # 计算参考面中心（世界空间）
        verts = ref_panel.data.vertices
        if len(verts) > 0:
            center_world = sum(
                (ref_panel.matrix_world @ v.co for v in verts),
                Vector()
            ) / len(verts)

            # 绕 X 轴旋转 180°
            R = Matrix.Rotation(math.pi, 4, 'X')
            T_to_center     = Matrix.Translation(center_world)
            T_from_center   = Matrix.Translation(-center_world)
            M_flip = T_to_center @ R @ T_from_center

            for obj in generated_studs:
                if obj:
                    obj.matrix_world = M_flip @ obj.matrix_world

    return generated_studs


def create_canonical_panel_from_polygon(src_obj, poly):
    """
    生成 canonical 面 + canonical→original 的变换矩阵 T。

    返回：
      panel_canonical   —— 世界坐标下法向 = (0, 0, 1)，几何在 canonical 空间
      T                 —— canonical → original polygon 的世界变换
    """
    # ==========================================================
    # 1. 提取 polygon 顶点（世界空间）
    # ==========================================================
    verts_world = [
        src_obj.matrix_world @ src_obj.data.vertices[i].co
        for i in poly.vertices
    ]

    # polygon 中心
    C = sum(verts_world, Vector()) / len(verts_world)

    # polygon 法向（世界空间）
    n = poly.normal.copy()
    n = (src_obj.matrix_world.to_3x3() @ n).normalized()   # 作为 +Z 方向

    # ==========================================================
    # 2. 构造 polygon 的局部坐标系 (t, b, n)
    # ==========================================================
    # 取最长边方向作为 t（面内某一方向）
    edges = []
    for i in range(len(verts_world)):
        v0 = verts_world[i]
        v1 = verts_world[(i + 1) % len(verts_world)]
        edges.append(v1 - v0)

    t = max(edges, key=lambda e: e.length).normalized()  # +X
    b = n.cross(t).normalized()                          # +Y，与 t、n 右手系

    # ==========================================================
    # 3. 构造 canonical→original 的世界变换矩阵 T
    # ==========================================================
    # canonical:
    #   +X → t
    #   +Y → b
    #   +Z → n
    #
    T = Matrix((
        (t.x,  b.x,  n.x,  C.x),
        (t.y,  b.y,  n.y,  C.y),
        (t.z,  b.z,  n.z,  C.z),
        (0.0,  0.0,  0.0,  1.0),
    ))

    # ==========================================================
    # 4. 在 canonical 空间构造 panel_canonical 的几何
    # ==========================================================
    mesh = bpy.data.meshes.new(f"{src_obj.name}_canonical_face_{poly.index}")
    panel_canonical = bpy.data.objects.new(mesh.name, mesh)
    src_obj.users_collection[0].objects.link(panel_canonical)

    T_inv = T.inverted()
    verts_canonical = [T_inv @ v for v in verts_world]

    # 顶点顺序可以保持原顺序，此时在 canonical 中法向大致为 (0, 0, 1)
    face_indices = tuple(range(len(verts_canonical)))

    mesh.from_pydata(
        [v.to_tuple() for v in verts_canonical],
        [],
        [face_indices],
    )
    mesh.update()

    # canonical 面保持 world_matrix = Identity：
    # 此时面在世界坐标下几何已经 canonical 化，
    # 法向约为 (0, 0, 1)
    panel_canonical.matrix_world = Matrix.Identity(4)

    return panel_canonical, T


# =================================================================
#  将单面 Object 设置成 IfcVirtualElement
# =================================================================

def assign_virtual_element(obj):
    """将 obj 标记为 IfcVirtualElement（需在 OBJECT 模式下调用）"""
    if not obj:
        return

    # 设为 active & 选中
    bpy.ops.object.select_all(action='DESELECT')
    obj.select_set(True)
    bpy.context.view_layer.objects.active = obj

    # 先进入 EDIT，全选面，再回到 OBJECT
    bpy.ops.object.mode_set(mode='EDIT')
    bpy.ops.mesh.select_all(action='SELECT')
    bpy.ops.object.mode_set(mode='OBJECT')

    # 指定 IFC 类
    bpy.ops.bim.assign_class(
        ifc_class="IfcVirtualElement",
        predefined_type="",
        userdefined_type="",
        props_to_pset=False,
    )

def parent_objects(parent_obj, children):
    """将 children 全部 parent 到 parent_obj （Blender + IFC）"""
    for child in children:
        if not child:
            continue
        child.parent = parent_obj

        # IFC parent-child 关系
        try:
            elem_util.assign_parent(child, parent_obj)
        except:
            pass

# =================================================================
#  多面 Mesh：循环生成各面龙骨 + 建立 IfcVirtualElement
# =================================================================

def generate_studs_on_mesh(context, model, props, ref_obj):
    """
    正确顺序：
      1. create canonical panel
      2. generate studs in canonical
      3. apply T to panel
      4. apply T to studs
      5. parent studs to panel
      6. assign IfcVirtualElement
    """

    all_studs = []

    if ref_obj.mode != "OBJECT":
        bpy.ops.object.mode_set(mode="OBJECT")

    mesh = ref_obj.data

    for poly in mesh.polygons:
        panel_canonical, T = create_canonical_panel_from_polygon(ref_obj, poly)
        studs_local = generate_studs_on_canonical_panel(
            context, model, props, panel_canonical
        )
        parent_objects(panel_canonical, studs_local)
        assign_virtual_element(panel_canonical)

        # move panel to original position
        panel_canonical.matrix_world = T
        all_studs.extend(studs_local)

        stud_log_append(context, f"✔ polygon {poly.index} 完成，生成 {len(studs_local)} 根龙骨")

    return all_studs


# =================================================================
#  UI 属性
# =================================================================

class StudDevProps(bpy.types.PropertyGroup):
    # ------------------------------
    # 主龙骨选定
    # ------------------------------

    # 选择 IfcMemberType
    selected_type: bpy.props.EnumProperty(
        name="主龙骨",
        items=update_stud_type_enum,
    )

    # Offset：在参考面局部坐标下的微调（让位石膏板等）
    offset_x: bpy.props.FloatProperty(
        name="Offset X",
        default=0.0,
        description="参考面局部 X 方向偏移",
        unit="LENGTH",
    )
    offset_y: bpy.props.FloatProperty(
        name="Offset Y",
        default=0.0,
        description="参考面局部 Y 方向偏移",
        unit="LENGTH",
    )
    offset_z: bpy.props.FloatProperty(
        name="Offset Z",
        default=0.019,           # 19mm 双层8.5mm石膏板厚度
        description="参考面局部 Z 方向偏移（常用于让位厚度）",
        unit="LENGTH",
    )

    # 绕挤出轴的 Roll，用来控制“覆面方向”
    roll_rad: bpy.props.FloatProperty(
        name="Roll",
        default=0.0,
        description="绕主轴旋转角度（控制覆面朝向）",
        unit="ROTATION",
    )

    # ------------------------------
    # 副龙骨选定
    # ------------------------------

    secondary_type: bpy.props.EnumProperty(
        name="副龙骨",
        items=update_stud_type_enum,
    )

    secondary_offset_x: bpy.props.FloatProperty(
        name="Offset X",
        default=0.0,
        description="副龙骨在参考面局部 X 方向偏移",
        unit="LENGTH",
    )
    secondary_offset_y: bpy.props.FloatProperty(
        name="Offset Y",
        default=0.0,
        description="副龙骨在参考面局部 Y 方向偏移",
        unit="LENGTH",
    )
    secondary_offset_z: bpy.props.FloatProperty(
        name="Offset Z",
        default=0.0005,
        description="副龙骨在参考面局部 Z 方向偏移",
        unit="LENGTH",
    )

    secondary_roll_rad: bpy.props.FloatProperty(
        name="Roll",
        default=0.0,
        description="副龙骨绕主轴旋转角度（覆面朝向）",
        unit="ROTATION",
    )

    # ------------------------------
    # 边龙骨选定
    # ------------------------------

    edge_type: bpy.props.EnumProperty(
        name="边龙骨",
        items=update_stud_type_enum,
    )

    edge_offset_x: bpy.props.FloatProperty(
        name="Offset X",
        default=0.0,
        description="边龙骨在参考面局部 X 方向偏移",
        unit="LENGTH",
    )
    edge_offset_y: bpy.props.FloatProperty(
        name="Offset Y",
        default=0.0,
        description="边龙骨在参考面局部 Y 方向偏移",
        unit="LENGTH",
    )
    edge_offset_z: bpy.props.FloatProperty(
        name="Offset Z",
        default=0.0,
        description="边龙骨在参考面局部 Z 方向偏移",
        unit="LENGTH",
    )

    edge_roll_rad: bpy.props.FloatProperty(
        name="Roll",
        default=0.0,
        description="边龙骨绕主轴旋转角度（覆面朝向）",
        unit="ROTATION",
    )

    # ------------------------------
    # 转角龙骨选定
    # ------------------------------

    corner_type: bpy.props.EnumProperty(
        name="转角龙骨",
        items=update_stud_type_enum,
    )

    corner_offset_x: bpy.props.FloatProperty(
        name="Offset X",
        default=0.0,
        description="转角龙骨在参考面局部 X 方向偏移",
        unit="LENGTH",
    )
    corner_offset_y: bpy.props.FloatProperty(
        name="Offset Y",
        default=0.0,
        description="转角龙骨在参考面局部 Y 方向偏移",
        unit="LENGTH",
    )
    corner_offset_z: bpy.props.FloatProperty(
        name="Offset Z",
        default=0.019,
        description="转角龙骨在参考面局部 Z 方向偏移",
        unit="LENGTH",
    )

    corner_roll_rad: bpy.props.FloatProperty(
        name="Roll",
        default=0.0,
        description="转角龙骨绕主轴旋转角度（覆面朝向）",
        unit="ROTATION",
    )

    # ------------------------------
    # 参考面排布
    # ------------------------------

    # 参考面对象
    ref_obj: bpy.props.PointerProperty(
        name="参考面",
        type=bpy.types.Object,
        description="用于排布龙骨的参考面（Mesh）",
    )

    # 阵列间距 = duplication 平移距离
    spacing: bpy.props.FloatProperty(
        name="主龙骨间距",
        default=0.6,  # 例：600mm 龙骨间距
        min=0.001,
        description="沿短边方向的排布间距（轴线间距）",
        unit="LENGTH",
    )

    # 副龙骨排布间距
    secondary_spacing: bpy.props.FloatProperty(
        name="副龙骨间距",
        default=0.3,
        min=0.001,
        description="副龙骨沿与主龙骨垂直方向的排布间距",
        unit="LENGTH",
    )

    log: bpy.props.StringProperty(default="")


# =================================================================
#  Operator：参考面 Scale ≠ 1，是否 Apply？
# =================================================================

class IFC_OT_ConfirmApplyScale(bpy.types.Operator):
    bl_idname = "ifc.confirm_apply_scale"
    bl_label = "参考面 Scale ≠ 1，是否 Apply？"

    ref_obj_name: bpy.props.StringProperty()
    original_operator: bpy.props.StringProperty(default="ifc.array_stud_from_multiref")

    def execute(self, context):
        obj = bpy.data.objects.get(self.ref_obj_name)
        if obj:
            bpy.ops.object.select_all(action='DESELECT')
            obj.select_set(True)
            context.view_layer.objects.active = obj
            bpy.ops.object.transform_apply(location=False, rotation=False, scale=True)
            stud_log_append(context, f"✔ 已应用参考面 Scale：{obj.name}")

        # 自动继续执行排布（无需再点一次按钮）
        bpy.ops.ifc.array_stud_from_multiref(bypass_scale_check=True)
        return {"FINISHED"}

    def invoke(self, context, event):
        return context.window_manager.invoke_props_dialog(self, width=300)

    def draw(self, context):
        layout = self.layout
        layout.label(text="参考面缩放不是 1，是否自动应用？")
        layout.label(text=f"对象：{self.ref_obj_name}")


# =================================================================
#  Operator：测试 offset_obj 生成（不排布龙骨）
# =================================================================

class IFC_OT_PolygonOffset(bpy.types.Operator):
    bl_idname = "ifc.polygon_offset"
    bl_label  = "测试 Polygon Offset"

    def execute(self, context):
        props = context.scene.stud_dev_props
        stud_log_set(context, "")

        # ----------------------------------------
        # 0. 获取参考对象
        # ----------------------------------------
        ref_obj = props.ref_obj or context.active_object
        if not ref_obj or ref_obj.type != "MESH":
            stud_log_set(context, "❌ 请指定参考面（Mesh）")
            return {"FINISHED"}

        # ----------------------------------------
        # 1. 创建 offset 副本（缩小一圈）
        # ----------------------------------------
        offset_value = props.offset_z
        try:
            offset_obj = create_offset_object_from_ref(ref_obj, offset_value)
        except Exception as e:
            stud_log_set(context, f"❌ 创建 offset_obj 失败：{e}")
            return {"FINISHED"}

        stud_log_append(context, 
            f"🎉 已成功基于 {ref_obj.name} 生成 offset 对象：{offset_obj.name}\n"
            f"   使用偏移量 offset = {offset_value:.4f} m"
        )

        return {"FINISHED"}



# =================================================================
#  Operator：为多面参考面生成龙骨（支持 scale 检查 + offset 预处理）
# =================================================================

class IFC_OT_ArrayStud_FromMultiRef(bpy.types.Operator):
    bl_idname = "ifc.array_stud_from_multiref"
    bl_label  = "为多面参考面生成龙骨"

    bypass_scale_check: bpy.props.BoolProperty(default=False)

    def execute(self, context):
        props = context.scene.stud_dev_props
        stud_log_set(context, "")

        model = get_ifc_model()
        if not model:
            stud_log_set(context, "❌ 无 IFC 模型")
            return {"CANCELLED"}

        # ------------------------------
        # 0. 获取参考对象
        # ------------------------------
        ref_obj = props.ref_obj or context.active_object
        if not ref_obj or ref_obj.type != "MESH":
            stud_log_set(context, "❌ 请选择一个 Mesh 作为参考面")
            return {"CANCELLED"}

        # ------------------------------
        # 0b. 验证参考 mesh（支持多面）
        # ------------------------------
        ok, reason = validate_reference_mesh(ref_obj)
        if not ok:
            stud_log_set(context, f"❌ 无法作为多面参考面：{reason}")
            return {"CANCELLED"}

        # ------------------------------
        # 1. Scale 检查（必须最前）
        # ------------------------------
        if not self.bypass_scale_check:
            sx, sy, sz = ref_obj.scale
            if (abs(sx - 1.0) > 1e-6) or (abs(sy - 1.0) > 1e-6) or (abs(sz - 1.0) > 1e-6):
                return bpy.ops.ifc.confirm_apply_scale(
                    'INVOKE_DEFAULT',
                    ref_obj_name=ref_obj.name
                )

        # ------------------------------
        # 2. 先创建 offset 对象（关键）
        # ------------------------------
        offset_dist = props.offset_z  # 默认 19 mm，双层石膏板
        offset_obj = create_offset_object_from_ref(ref_obj, offset_dist)

        if offset_obj is None:
            stud_log_set(context, "❌ 创建 offset_obj 失败（请检查参考面是否封闭、几何是否异常）")
            return {"CANCELLED"}

        stud_log_append(context, f"✔ 创建 offset_obj：{offset_obj.name}")

        # ------------------------------
        # 3. 执行多面排布（主龙骨、副龙骨、边龙骨）
        # ------------------------------
        studs = generate_studs_on_mesh(
            context, model, props,
            offset_obj,
        )

        if not studs:
            stud_log_append(context, "⚠ 未生成任何龙骨")
        else:
            stud_log_append(context, f"🎉 多面龙骨生成完成，共 {len(studs)} 根")

        return {"FINISHED"}


# =================================================================
#  UI
# =================================================================

class IFC_PT_StudDevPanel(bpy.types.Panel):
    bl_idname = "IFC_PT_StudDevPanel"
    bl_label = "Stud Dev Tools"
    bl_space_type = "VIEW_3D"
    bl_region_type = "UI"
    bl_category = "IFC"

    def draw(self, context):
        layout = self.layout
        props = context.scene.stud_dev_props

        col = layout.column(align=True)
        col.label(text="主龙骨：")
        col.prop(props, "selected_type", text="")

        col.separator()
        col.label(text="主龙骨偏移：")
        col.prop(props, "offset_x")
        col.prop(props, "offset_y")
        col.prop(props, "offset_z")
        col.prop(props, "roll_rad")

        col.separator()
        col.label(text="副龙骨：")
        col.prop(props, "secondary_type", text="")
        col.prop(props, "secondary_offset_x")
        col.prop(props, "secondary_offset_y")
        col.prop(props, "secondary_offset_z")
        col.prop(props, "secondary_roll_rad")

        col.separator()
        col.label(text="边龙骨：")
        col.prop(props, "edge_type", text="")
        col.prop(props, "edge_offset_x")
        col.prop(props, "edge_offset_y")
        col.prop(props, "edge_offset_z")
        col.prop(props, "edge_roll_rad")

        col.separator()
        col.label(text="转角龙骨：")
        col.prop(props, "corner_type", text="")
        col.prop(props, "corner_offset_x")
        col.prop(props, "corner_offset_y")
        col.prop(props, "corner_offset_z")
        col.prop(props, "corner_roll_rad")

        col.separator()
        col.label(text="参考面：")
        col.prop(props, "ref_obj", text="")
        col.prop(props, "spacing")
        col.prop(props, "secondary_spacing")
        col.separator()
        col.operator("ifc.polygon_offset", text="测试 Polygon Offset")
        col.operator("ifc.array_stud_from_multiref", text="生成多面龙骨")

        col.separator()
        col.label(text="日志：")
        col.prop(props, "log")


# =================================================================
#  注册逻辑
# =================================================================

classes = (
    StudDevProps,
    IFC_PT_StudDevPanel,
    IFC_OT_ConfirmApplyScale,
    IFC_OT_PolygonOffset,
    IFC_OT_ArrayStud_FromMultiRef
)

def register():
    for cls in classes:
        bpy.utils.register_class(cls)
    bpy.types.Scene.stud_dev_props = bpy.props.PointerProperty(type=StudDevProps)

def unregister():
    for cls in reversed(classes):
        bpy.utils.unregister_class(cls)
    del bpy.types.Scene.stud_dev_props


if __name__ == "__main__":
    register()
