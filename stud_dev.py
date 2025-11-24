import bpy
import math
from mathutils import Vector, Matrix
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
                import math as _math
                count = int(_math.ceil(target_len / unit_len))
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

    axes = ["x", "y", "z"]

    # 厚度轴
    thickness_axis = min(axes, key=lambda a: lengths[a])

    # 平面轴 = 剩下两个
    plane_axes = [a for a in axes if a != thickness_axis]

    # 长边/短边判断
    a1, a2 = plane_axes
    if lengths[a1] >= lengths[a2]:
        long_axis = a1
        short_axis = a2
    else:
        long_axis = a2
        short_axis = a1

    # ---------------------------
    # 主龙骨 extrusion start/end（沿长边）
    # ---------------------------
    coords_start = dict(centers)
    coords_end = dict(centers)

    # 厚度轴取中心（避免偏移）
    coords_start[thickness_axis] = centers[thickness_axis]
    coords_end[thickness_axis] = centers[thickness_axis]

    # 短边取 MIN（一侧靠齐）
    coords_start[short_axis] = mins[short_axis]
    coords_end[short_axis] = mins[short_axis]

    # 长边 extrusion：min → max
    coords_start[long_axis] = mins[long_axis]
    coords_end[long_axis] = maxs[long_axis]

    local_main_start = Vector((coords_start["x"], coords_start["y"], coords_start["z"]))
    local_main_end   = Vector((coords_end["x"],   coords_end["y"],   coords_end["z"]))

    # ---------------------------
    # 副龙骨 extrusion start/end（沿短边）
    # ---------------------------
    sec_coords_start = dict(centers)
    sec_coords_end   = dict(centers)

    # 厚度轴中心
    sec_coords_start[thickness_axis] = centers[thickness_axis]
    sec_coords_end[thickness_axis]   = centers[thickness_axis]

    # 长边取 MIN（与主龙骨一致）
    sec_coords_start[long_axis] = mins[long_axis]
    sec_coords_end[long_axis]   = mins[long_axis]

    # **短边 extrusion：min → max**
    sec_coords_start[short_axis] = mins[short_axis]
    sec_coords_end[short_axis]   = maxs[short_axis]

    local_sec_start = Vector((sec_coords_start["x"], sec_coords_start["y"], sec_coords_start["z"]))
    local_sec_end   = Vector((sec_coords_end["x"],   sec_coords_end["y"],   sec_coords_end["z"]))

    # ---------------------------
    # 构造方向向量
    # ---------------------------

    if short_axis == "x":
        local_short_axis_vec = Vector((1, 0, 0))
    elif short_axis == "y":
        local_short_axis_vec = Vector((0, 1, 0))
    else:
        local_short_axis_vec = Vector((0, 0, 1))

    if long_axis == "x":
        local_long_axis_vec = Vector((1, 0, 0))
    elif long_axis == "y":
        local_long_axis_vec = Vector((0, 1, 0))
    else:
        local_long_axis_vec = Vector((0, 0, 1))

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
    for i in range(count):
        world_start = verts[i]
        world_end   = verts[(i + 1) % count]

        local_start = inv_mw @ world_start + local_offset
        local_end   = inv_mw @ world_end   + local_offset

        create_stud_instance(
            context,
            model,
            type_obj,
            mw @ local_start,
            mw @ local_end,
            roll_rad,
        )

    stud_log_append(context, f"✔ 描边龙骨已生成，共 {count} 条")


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
        default=0.0,
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
        default=0.0,
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
    original_operator: bpy.props.StringProperty(default="ifc.array_stud_from_ref")

    def execute(self, context):
        obj = bpy.data.objects.get(self.ref_obj_name)
        if obj:
            bpy.ops.object.select_all(action='DESELECT')
            obj.select_set(True)
            context.view_layer.objects.active = obj
            bpy.ops.object.transform_apply(location=False, rotation=False, scale=True)
            stud_log_append(context, f"✔ 已应用参考面 Scale：{obj.name}")

        # 自动继续执行排布（无需再点一次按钮）
        bpy.ops.ifc.array_stud_from_ref(bypass_scale_check=True)
        return {"FINISHED"}

    def invoke(self, context, event):
        return context.window_manager.invoke_props_dialog(self, width=300)

    def draw(self, context):
        layout = self.layout
        layout.label(text="参考面缩放不是 1，是否自动应用？")
        layout.label(text=f"对象：{self.ref_obj_name}")


# =================================================================
#  Operator：参考面 → 排布龙骨（阵列）
# =================================================================

class IFC_OT_ArrayStud_FromRef(bpy.types.Operator):
    bl_idname = "ifc.array_stud_from_ref"
    bl_label = "参考面 → 排布龙骨"

    bypass_scale_check: bpy.props.BoolProperty(default=False)

    def execute(self, context):
        props = context.scene.stud_dev_props
        stud_log_set(context, "")

        # ----------------------------------------
        # 0. 取得参考面对象
        # ----------------------------------------
        ref_obj = props.ref_obj or context.active_object
        if not ref_obj or ref_obj.type != "MESH":
            stud_log_set(context, "❌ 请指定参考面（Mesh）")
            return {"FINISHED"}

        # ----------------------------------------
        # 1. Scale 检查（必须最前）
        # ----------------------------------------
        if not self.bypass_scale_check:
            sx, sy, sz = ref_obj.scale
            if (abs(sx - 1.0) > 1e-6) or (abs(sy - 1.0) > 1e-6) or (abs(sz - 1.0) > 1e-6):
                return bpy.ops.ifc.confirm_apply_scale(
                    'INVOKE_DEFAULT',
                    ref_obj_name=ref_obj.name
                )

        # ----------------------------------------
        # 2. IFC 模型检查
        # ----------------------------------------
        model = get_ifc_model()
        if not model:
            stud_log_set(context, "❌ 无 IFC 模型")
            return {"FINISHED"}

        # ----------------------------------------
        # 3. 主龙骨类型检查
        # ----------------------------------------
        type_obj = find_member_type(model, props.selected_type)
        if not type_obj:
            stud_log_set(context, "❌ 未选择类型")
            return {"FINISHED"}

        # ----------------------------------------
        # 4. 解析参考面
        # ----------------------------------------
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
            ) = analyse_reference_panel(ref_obj)
        except Exception as e:
            stud_log_set(context, f"❌ 参考面分析失败: {e}")
            return {"FINISHED"}

        # ===================================================================
        # 5~14. 布局计算：调用独立函数（compute_stud_layout）
        # ===================================================================
        layout_info = compute_stud_layout(
            short_length=short_length,
            long_length=long_length,
            spacing=props.spacing,
            sec_spacing=props.secondary_spacing,
            local_main_start=local_main_start,
            local_main_end=local_main_end,
            local_sec_start=local_sec_start,
            local_sec_end=local_sec_end,
        )

        sec_count = layout_info["sec_count"]
        main_count = layout_info["main_count"]
        main_extrude_len = layout_info["main_extrude_len"]
        main_short_offset = layout_info["main_short_offset"]
        main_long_offset = layout_info["main_long_offset"]
        sec_long_offset = layout_info["sec_long_offset"]
        adjusted_sec_start = layout_info["adjusted_sec_start"]
        adjusted_sec_end = layout_info["adjusted_sec_end"]
        adjusted_main_start = layout_info["adjusted_main_start"]
        adjusted_main_end = layout_info["adjusted_main_end"]
        sec_len_original = layout_info["sec_len_original"]
        sec_len_new = layout_info["sec_len_new"]

        stud_log_append(context, f"ℹ 副龙骨数量 = {sec_count}")
        stud_log_append(context, f"ℹ 主龙骨挤出长度 = {main_extrude_len:.4f}")
        stud_log_append(context, f"ℹ 主龙骨数量 = {main_count}")
        stud_log_append(context, f"ℹ 主龙骨短边偏移 = {main_short_offset:.4f}")
        stud_log_append(context, f"ℹ 主龙骨长边偏移 = {main_long_offset:.4f}")
        stud_log_append(context, f"ℹ 副龙骨长边偏移 = {sec_long_offset:.4f}")
        stud_log_append(
            context,
            f"✔ 副龙骨调整: 原长度={sec_len_original:.4f} → 新长度={sec_len_new:.4f}"
        )

        # ===================================================================
        # 12. 准备偏移（local 坐标）
        # ===================================================================
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

        # ===================================================================
        # 15. 边龙骨描边
        # ===================================================================
        edge_type_obj = find_member_type(model, props.edge_type)
        if edge_type_obj:
            local_edge_offset = Vector((
                props.edge_offset_x,
                props.edge_offset_y,
                props.edge_offset_z,
            ))
            outline_studs_on_reference(
                context,
                model,
                edge_type_obj,
                ref_obj,
                local_edge_offset,
                props.edge_roll_rad,
            )

        # ===================================================================
        # 16. 主骨排布（正确数量、正确位置）
        # ===================================================================
        array_studs_on_reference(
            context,
            model,
            type_obj,
            ref_obj,
            adjusted_main_start,
            adjusted_main_end,
            local_main_offset,
            local_short_axis_vec,
            main_count,          # 由布局计算函数直接给出数量
            props.spacing,       # 由属性保证必须 > 0
            props.roll_rad,
        )

        # ===================================================================
        # 17. 副骨排布
        # ===================================================================
        secondary_type_obj = find_member_type(model, props.secondary_type)

        if secondary_type_obj:
            array_studs_on_reference(
                context,
                model,
                secondary_type_obj,
                ref_obj,
                adjusted_sec_start,
                adjusted_sec_end,
                local_sec_offset,
                local_long_axis_vec,
                sec_count,              # 由布局计算函数直接给出数量
                props.secondary_spacing,
                props.secondary_roll_rad,
            )

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
        col.operator("ifc.array_stud_from_ref", text="生成龙骨")

        col.separator()
        col.label(text="日志：")
        col.prop(props, "log")


# =================================================================
#  注册逻辑
# =================================================================

classes = (
    StudDevProps,
    IFC_PT_StudDevPanel,
    IFC_OT_ArrayStud_FromRef,
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
