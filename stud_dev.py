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
#  Corner 工具：识别阴阳角 + 侧面 + 20mm 偏移
# =================================================================

def find_corner_edges(offset_obj,
                      include_side_side_corners: bool = False,
                      side_z_eps: float = 0.01):
    """
    基于 offset_obj 分析所有 polygon，识别：
      - 每个面的世界坐标顶点
      - 每条边是否为转角边（is_corner）
      - 阴角 / 阳角 corner_type = "inner" / "outer"
      - 邻接面的世界法向 neighbor_normal

    参数：
      include_side_side_corners:
        False（默认）时：两个“侧面”之间的交线不视为 corner；
        True 时：保持原逻辑，侧面-侧面交线也算 corner。

      side_z_eps:
        判断侧面用的阈值，和 classify_side_faces 保持一致：
          |normal.z| < side_z_eps → 侧面

    返回 panels_raw 结构：
        {
          poly_index: {
            "verts_world": [Vector, ...]   # 按 polygon 顶点顺序
            "edges": [
              {
                "v1_index": int,
                "v2_index": int,
                "v1_world": Vector,
                "v2_world": Vector,
                "is_corner": bool,
                "corner_type": "inner"/"outer"/None,
                "neighbor_normal": Vector or None,
              },
              ...
            ],
            "normal_world": Vector,
          },
          ...
        }
    """
    mesh = offset_obj.data
    mw = offset_obj.matrix_world

    # 多边形世界法向
    poly_normals = {}
    for poly in mesh.polygons:
        n_world = (mw.to_3x3() @ poly.normal).normalized()
        poly_normals[poly.index] = n_world

    # 与 classify_side_faces 同一规则：|normal.z| < 0.01 → 侧面
    is_side = {
        idx: (abs(n.z) < side_z_eps)
        for idx, n in poly_normals.items()
    }

    # 边与多边形的关联：key = (min(v1,v2), max(v1,v2))
    edge_map = {}
    for poly in mesh.polygons:
        verts = poly.vertices
        n = len(verts)
        for i in range(n):
            v1 = verts[i]
            v2 = verts[(i + 1) % n]
            key = tuple(sorted((v1, v2)))
            edge_map.setdefault(key, []).append((poly.index, i))

    panels_raw = {}

    for poly in mesh.polygons:
        verts = poly.vertices
        n = len(verts)

        data = {
            "verts_world": [mw @ mesh.vertices[i].co for i in verts],
            "edges": [],
            "normal_world": poly_normals[poly.index],
        }
        panels_raw[poly.index] = data

        for i in range(n):
            v1_idx = verts[i]
            v2_idx = verts[(i + 1) % n]
            v1_world = mw @ mesh.vertices[v1_idx].co
            v2_world = mw @ mesh.vertices[v2_idx].co

            key = tuple(sorted((v1_idx, v2_idx)))
            attached = edge_map.get(key, [])

            is_corner = False
            corner_type = None
            neighbor_normal = None

            if len(attached) == 2:
                (pA, _), (pB, _) = attached
                other_poly_idx = pB if pA == poly.index else pA

                n_self = poly_normals[poly.index]
                n_other = poly_normals[other_poly_idx]

                # 如果不希望把 “两个侧面之间” 的交线当成 corner，则直接跳过
                if (
                    not include_side_side_corners
                    and is_side.get(poly.index, False)
                    and is_side.get(other_poly_idx, False)
                ):
                    # 保持 is_corner=False, neighbor_normal=None
                    pass
                else:
                    dot = n_self.dot(n_other)

                    # 垂直：认为是转角边
                    if abs(dot) < 1e-3:
                        is_corner = True
                        neighbor_normal = n_other

                        e_dir = (v2_world - v1_world)
                        if e_dir.length > 1e-6:
                            e_dir.normalize()
                            cross_n = n_self.cross(n_other)
                            sign = cross_n.dot(e_dir)
                            if sign > 0:
                                corner_type = "outer"
                            else:
                                corner_type = "inner"

            data["edges"].append({
                "v1_index": v1_idx,
                "v2_index": v2_idx,
                "v1_world": v1_world,
                "v2_world": v2_world,
                "vert_key": key,
                "is_corner": is_corner,
                "corner_type": corner_type,
                "neighbor_normal": neighbor_normal,
            })

    return panels_raw


# ===============================================================
#  BUTT INNER CORNER：基于 offset_obj（正确偏移后）
# ===============================================================

def find_butt_inner_corners_final(offset_obj,
                                  plane_eps=1e-3,
                                  aabb_eps=1e-3,
                                  dot_orth_eps=1e-3):
    """
    检测 “非侧面 B 的边贴到 侧面 A 的偏移平面（offset_obj）内部”的 butt inner corner。

    返回：
        butt_map[(polyB_idx, edge_key)] = {
            "polyA": polyA_idx,         # 侧面
            "neighbor_normal": normal_of_A,
            "v1_world": p1_proj,        # 投影后的点
            "v2_world": p2_proj,
            "dist": average_plane_dist,
        }
    """

    mesh = offset_obj.data
    mw   = offset_obj.matrix_world

    # --- 构建缓存 ---
    poly_count = len(mesh.polygons)

    poly_normals = [ (mw.to_3x3() @ poly.normal).normalized()
                      for poly in mesh.polygons ]

    poly_centers = [
        sum((mw @ mesh.vertices[i].co for i in poly.vertices), Vector()) / len(poly.vertices)
        for poly in mesh.polygons
    ]

    # AABB（用于判断投影是否在面内）
    poly_aabb_min = []
    poly_aabb_max = []
    for poly in mesh.polygons:
        ws = [mw @ mesh.vertices[i].co for i in poly.vertices]
        poly_aabb_min.append(Vector((min(p.x for p in ws),
                                     min(p.y for p in ws),
                                     min(p.z for p in ws))))
        poly_aabb_max.append(Vector((max(p.x for p in ws),
                                     max(p.y for p in ws),
                                     max(p.z for p in ws))))

    # 侧面判断：法向不接近 Z 方向 → 侧面
    side_mask = [abs(n.dot(Vector((0,0,1)))) < 0.9 for n in poly_normals]

    side_polys    = [i for i in range(poly_count) if side_mask[i]]
    nonside_polys = [i for i in range(poly_count) if not side_mask[i]]

    butt_map = {}

    # ==========================================================
    # 遍历：非侧面 B 的边是否贴在 侧面 A 上
    # ==========================================================
    for pb in nonside_polys:
        polyB = mesh.polygons[pb]
        vertsB = list(polyB.vertices)
        nB = poly_normals[pb]

        vcount = len(vertsB)

        for i in range(vcount):
            v1 = vertsB[i]
            v2 = vertsB[(i + 1) % vcount]

            p1 = mw @ mesh.vertices[v1].co
            p2 = mw @ mesh.vertices[v2].co

            edge_key = tuple(sorted((v1, v2)))

            # 遍历所有侧面 A
            for pa in side_polys:
                polyA = mesh.polygons[pa]
                nA = poly_normals[pa]
                cA = poly_centers[pa]

                # 1. 法向必须有明显夹角（即 A/B 应该接近垂直）
                if abs(nA.dot(nB)) > dot_orth_eps:
                    continue

                # 2. 判断 p1/p2 是否靠近 A 的平面
                d1 = abs((p1 - cA).dot(nA))
                d2 = abs((p2 - cA).dot(nA))

                if d1 > plane_eps or d2 > plane_eps:
                    continue

                # 3. 投影到 A 平面
                p1_proj = p1 - nA * ( (p1 - cA).dot(nA) )
                p2_proj = p2 - nA * ( (p2 - cA).dot(nA) )

                # 4. 投影点是否落在 A 的 AABB 内
                bb_min = poly_aabb_min[pa]
                bb_max = poly_aabb_max[pa]

                def in_aabb(q):
                    return (
                        (bb_min.x - aabb_eps) <= q.x <= (bb_max.x + aabb_eps) and
                        (bb_min.y - aabb_eps) <= q.y <= (bb_max.y + aabb_eps) and
                        (bb_min.z - aabb_eps) <= q.z <= (bb_max.z + aabb_eps)
                    )

                if not (in_aabb(p1_proj) and in_aabb(p2_proj)):
                    continue

                # --- 找到 BUTT INNER CORNER ---
                butt_map[(pb, i)] = {
                    "polyA": pa,
                    "neighbor_normal": nA,
                    "v1_world": p1_proj,
                    "v2_world": p2_proj,
                    "dist": (d1 + d2) * 0.5,
                }
                break   # 当前边无需继续找其它 A

    return butt_map


# ===============================================================
#  将新版 butt-corner 信息合并进 panels_raw
#  （基于 poly_idx + edge_idx）
# ===============================================================

def merge_butt_inner_corners_into_panels(panels_raw, butt_map):
    """
    适配新版 butt_map（key = (poly_idx, edge_idx)）。

    但不能覆盖已有的拓扑 corner（特别是 outer），
    只对 non-corner 或 inner 类型的边进行补充。
    """

    for (poly_idx, edge_idx), info in butt_map.items():

        if poly_idx not in panels_raw:
            continue

        edges = panels_raw[poly_idx].get("edges")
        if not edges:
            continue

        if edge_idx < 0 or edge_idx >= len(edges):
            continue

        edge = edges[edge_idx]

        # ======================================================
        # 🚫 不能覆盖已有 corner（特别是拓扑 outer）
        # ======================================================
        if edge.get("is_corner"):
            # 已经是 outer → 绝对不能覆盖
            if edge.get("corner_type") == "outer":
                continue

            # 已经是 inner，但来自拓扑检测 → 不覆盖
            # 保留你需要的行为，可以根据需要决定
            if edge.get("is_butt") is not True:
                continue

        # ======================================================
        # ✔ 覆盖 / 添加 butt inner corner 信息
        # ======================================================
        edge["is_corner"] = True
        edge["corner_type"] = "inner"
        edge["neighbor_normal"] = info.get("neighbor_normal")
        edge["is_butt"] = True


def classify_side_faces(panels_raw):
    """
    标记侧面 / 非侧面：
      - |normal.z| < 阈值 → 视为侧面（墙面等）
      - 其余视为非侧面（顶棚 / 地面）
    """
    for poly_idx, data in panels_raw.items():
        n_world = data["normal_world"]
        data["is_side"] = abs(n_world.z) < 0.01


def apply_corner_offset(panels_raw, offset_dist=0.02):
    """
    仅对【侧面 + 阳角】的 corner edge 做 20mm 偏移（沿邻面法向反向）。

    修复点：
      - 以前是“按 edge 单独偏移”，会导致同一顶点在不同 edge 上有不同 final 位置；
      - 现在改为“按顶点聚合偏移”，保证同一 polygon 内共享顶点的所有边使用同一 final 坐标。
    """
    def vkey(v):
        # 使用坐标做 key，防止浮点误差
        return (round(v.x, 6), round(v.y, 6), round(v.z, 6))

    for poly_idx, data in panels_raw.items():
        is_side = data.get("is_side", False)
        edges = data.get("edges", [])

        if not edges:
            continue

        # 1️⃣ 建立顶点表：world → final（初始 = world）
        vert_map = {}
        for e in edges:
            for vw in (e["v1_world"], e["v2_world"]):
                k = vkey(vw)
                if k not in vert_map:
                    # 存 world & final，两者起初相同
                    vert_map[k] = {
                        "world": vw.copy(),
                        "final": vw.copy(),
                    }

        # 2️⃣ 对满足条件的 corner edge，将其两个端点的 final 一起偏移
        for e in edges:
            if not (is_side and e.get("is_corner") and e.get("corner_type") == "outer"):
                continue
            neigh_n = e.get("neighbor_normal")
            if neigh_n is None:
                continue

            offset_vec = neigh_n * (-offset_dist)

            for vw in (e["v1_world"], e["v2_world"]):
                k = vkey(vw)
                vert_map[k]["final"] = vert_map[k]["final"] + offset_vec

        # 3️⃣ 回写到每条 edge 的 v*_final（从同一 vert_map 读，保证连续性）
        for e in edges:
            k1 = vkey(e["v1_world"])
            k2 = vkey(e["v2_world"])
            e["v1_final"] = vert_map[k1]["final"].copy()
            e["v2_final"] = vert_map[k2]["final"].copy()


# ===============================================================
#  Corner + 侧面 + 阴角 + 20mm 偏移 一站式封装
# ===============================================================

def build_corner_panels_data(context, ref_obj, offset_obj, offset_dist=0.019):
    """
    一站式封装：

      1. 基于 offset_obj 识别拓扑 corner（原有逻辑）
      2. 基于 ref_obj 额外识别 butt 阴角（非拓扑、贴合），并合并到 panels_raw
      3. 标记侧面 is_side
      4. 仅对【侧面 + 阳角】corner 边做 20mm 偏移（沿邻面法向反向）
      5. 返回 panels_raw

    最终 panels_raw 结构：
      {
        poly_idx: {
          "normal_world": Vector,
          "is_side": bool,
          "edges": [
            {
              "v1_world": Vector,
              "v2_world": Vector,
              "v1_final": Vector,
              "v2_final": Vector,
              "vert_key": (vi, vj),
              "is_corner": bool,
              "corner_type": "inner"/"outer"/None,
              "neighbor_normal": Vector|None,
              ...
            },
            ...
          ]
        },
        ...
      }
    """
    # 1. 原有 corner 检测（基于 offset_obj）
    panels_raw = find_corner_edges(offset_obj)

    # 2. 新增：在 ref_obj 上检测“贴在顶面的侧边阴角”
    butt_map = find_butt_inner_corners_final(offset_obj, offset_dist)
    merge_butt_inner_corners_into_panels(panels_raw, butt_map)

    # 3. 标记侧面（你原来的函数）
    classify_side_faces(panels_raw)

    # 4. 对 侧面+阳角 corner 边做 20mm 偏移
    apply_corner_offset(panels_raw, offset_dist=offset_dist)

    return panels_raw


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


def build_final_offset_mesh_from_panels_raw(panels_raw, base_name="OffsetFinal"):
    """
    基于 panels_raw 构建最终 offset mesh（方案 A：所有 polygon 顶点独立）
    仅使用 e["v1_final"] 作为 polygon 的顶点顺序来源。

    返回：
        Object（新创建的 mesh 对象）
    """
    # 创建 mesh & object
    mesh = bpy.data.meshes.new(f"{base_name}_mesh")
    obj = bpy.data.objects.new(f"{base_name}_obj", mesh)
    bpy.context.collection.objects.link(obj)

    bm = bmesh.new()

    # 逐 polygon 构建几何
    for poly_idx, pdata in panels_raw.items():
        edges = pdata["edges"]

        # 按 edges 顺序得到 polygon 顶点列表
        verts_world = [e["v1_final"] for e in edges]

        bm_face_verts = []
        for v in verts_world:
            bm_face_verts.append(bm.verts.new(v))

        # 创建面
        try:
            bm.faces.new(bm_face_verts)
        except ValueError:
            # 如果 face 已经存在（理论上不会发生，因为全部独立顶点）
            pass

    # 输出 mesh
    bm.to_mesh(mesh)
    mesh.update()

    return obj


# =================================================================
#  获取参考面顶点（按顺时针排序）
# =================================================================

def get_ordered_face_edges(ref_obj):
    """
    返回单个面的边列表：[(world_v1, world_v2, edge_index), ...]
    约定：
      - 只处理单 polygon（canonical panel 正是这种情况）
      - edge_index 按 polygon 顶点顺序从 0..n-1
    """
    mesh = ref_obj.data
    if not mesh.polygons:
        return []

    poly = mesh.polygons[0]
    verts = poly.vertices
    n = len(verts)

    edges = []
    for i in range(n):
        v1_idx = verts[i]
        v2_idx = verts[(i + 1) % n]
        v1_world = ref_obj.matrix_world @ mesh.vertices[v1_idx].co
        v2_world = ref_obj.matrix_world @ mesh.vertices[v2_idx].co
        edges.append((v1_world, v2_world, i))

    return edges


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
    edge_type_obj,
    ref_obj,
    edge_offset,
    edge_roll_rad,
    is_side_face=False,
    corner_edge_indices=None,
    corner_edge_types=None,
    corner_inner_type_obj=None,
    corner_outer_type_obj=None,
    corner_offset=None,
    corner_roll_rad=None,
):
    """沿参考面边缘生成描边龙骨，支持转角（阴/阳角）判断"""

    edges = get_ordered_face_edges(ref_obj)
    if len(edges) < 1:
        stud_log_append(context, "⚠ 参考面边数不足，无法描边")
        return

    mw = ref_obj.matrix_world
    inv_mw = mw.inverted()

    studs = []

    for world_start, world_end, edge_idx in edges:
        # 角信息
        is_corner = corner_edge_indices is not None and edge_idx in corner_edge_indices
        ctype = None
        if is_corner and corner_edge_types:
            ctype = corner_edge_types.get(edge_idx)

        # 决定本条边用什么类型 & 偏移 & roll
        stud_type = edge_type_obj
        offset_vec = edge_offset
        roll = edge_roll_rad

        if is_corner:
            # 侧面上的 corner：不生成转角龙骨，直接跳过
            if is_side_face:
                continue

            # 非侧面：按阴角/阳角用不同类型
            if ctype == "inner" and corner_inner_type_obj:
                stud_type = corner_inner_type_obj
            elif ctype == "outer" and corner_outer_type_obj:
                stud_type = corner_outer_type_obj
            # 使用 corner offset / roll（如果有）
            if corner_offset is not None:
                offset_vec = corner_offset
            if corner_roll_rad is not None:
                roll = corner_roll_rad

        # 若未设置类型则跳过
        if not stud_type:
            continue

        # 变换到 local，再加偏移
        local_start = inv_mw @ world_start + offset_vec
        local_end   = inv_mw @ world_end   + offset_vec

        world_start_off = mw @ local_start
        world_end_off   = mw @ local_end

        stud_obj = create_stud_instance(
            context,
            model,
            stud_type,
            world_start_off,
            world_end_off,
            roll,
        )
        studs.append(stud_obj)

    stud_log_append(context, f"✔ 描边龙骨已生成，共 {len(studs)} 条")
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

def generate_studs_on_canonical_panel(
    context,
    model,
    props,
    ref_panel,
    is_side_face=False,
    corner_edge_indices=None,
    corner_edge_types=None,
):
    """
    在 canonical 面（ref_panel）上生成龙骨阵列。

    新增：
      - is_side_face         ：该面是否为侧面
      - corner_edge_indices  ：哪些边是 corner
      - corner_edge_types    ：corner 边是 "inner"/"outer"
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

    # 5. 边龙骨（含转角龙骨）
    edge_type_obj         = find_member_type(model, props.edge_type)
    corner_inner_type_obj = find_member_type(model, props.corner_inner_type)
    corner_outer_type_obj = find_member_type(model, props.corner_outer_type)

    if edge_type_obj:
        edge_offset = Vector((
            props.edge_offset_x,
            props.edge_offset_y,
            props.edge_offset_z,
        ))
        corner_offset = Vector((
            props.corner_offset_x,
            props.corner_offset_y,
            props.corner_offset_z,
        ))

        studs = outline_studs_on_reference(
            context, model,
            edge_type_obj,
            ref_panel,
            edge_offset,
            props.edge_roll_rad,
            is_side_face=is_side_face,
            corner_edge_indices=corner_edge_indices,
            corner_edge_types=corner_edge_types,
            corner_inner_type_obj=corner_inner_type_obj,
            corner_outer_type_obj=corner_outer_type_obj,
            corner_offset=corner_offset,
            corner_roll_rad=props.corner_roll_rad,
        )
        _append(studs)

    # 6. 主龙骨阵列
    if is_side_face:
        pass
    else:
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

    return generated_studs


def create_canonical_panel_from_polygon(
    src_obj,
    poly,
    verts_world_override=None,
    corner_edge_flags=None,
    corner_edge_types=None,
):
    """
    生成 canonical 面 + canonical→original 的变换矩阵 T。

    可选：
      - verts_world_override: 使用外部提供的世界坐标顶点（按 polygon 顶点顺序）
      - corner_edge_flags/types: 每条边是否为 corner 以及阴/阳角类型
    """

    # ==========================================================
    # 1. 提取 polygon 顶点（世界空间）
    # ==========================================================
    if verts_world_override is not None:
        verts_world = list(verts_world_override)
    else:
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
    edges_vec = []
    for i in range(len(verts_world)):
        v0 = verts_world[i]
        v1 = verts_world[(i + 1) % len(verts_world)]
        edges_vec.append(v1 - v0)

    t = max(edges_vec, key=lambda e: e.length).normalized()  # +X
    b = n.cross(t).normalized()                              # +Y，与 t、n 右手系

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

    face_indices = tuple(range(len(verts_canonical)))

    mesh.from_pydata(
        [v.to_tuple() for v in verts_canonical],
        [],
        [face_indices],
    )
    mesh.update()

    # canonical 面保持 world_matrix = Identity：
    panel_canonical.matrix_world = Matrix.Identity(4)

    # ==========================================================
    # 5. Corner 边信息映射到 canonical 边 index
    # ==========================================================
    corner_edge_indices = set()
    corner_edge_type_map = {}

    if corner_edge_flags is not None:
        for i, is_corner in enumerate(corner_edge_flags):
            if is_corner:
                corner_edge_indices.add(i)
                if corner_edge_types and i < len(corner_edge_types):
                    ctype = corner_edge_types[i]
                    if ctype:
                        corner_edge_type_map[i] = ctype

    return panel_canonical, T, corner_edge_indices, corner_edge_type_map


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

def generate_studs_on_mesh(context, model, props, ref_obj, panels_raw):
    """
    ref_obj: 已经 offset 过的 offset_obj
    panels_raw: 外部预先构建好的 corner + side + v*_final 数据
    """

    all_studs = []
    if ref_obj.mode != "OBJECT":
        bpy.ops.object.mode_set(mode="OBJECT")

    mesh = ref_obj.data

    for poly in mesh.polygons:
        panel_data = panels_raw.get(poly.index)
        if not panel_data:
            continue

        edges_data = panel_data["edges"]
        is_side    = panel_data.get("is_side", False)

        # canonical 顶点序列 = v1_final 列表
        verts_world_override = [e["v1_final"] for e in edges_data]

        corner_flags = [e["is_corner"] for e in edges_data]
        corner_types = [e["corner_type"] for e in edges_data]

        panel_canonical, T, corner_edge_indices, corner_edge_types = create_canonical_panel_from_polygon(
            ref_obj,
            poly,
            verts_world_override=verts_world_override,
            corner_edge_flags=corner_flags,
            corner_edge_types=corner_types,
        )

        studs_local = generate_studs_on_canonical_panel(
            context,
            model,
            props,
            panel_canonical,
            is_side_face=is_side,
            corner_edge_indices=corner_edge_indices,
            corner_edge_types=corner_edge_types,
        )

        parent_objects(panel_canonical, studs_local)
        assign_virtual_element(panel_canonical)

        panel_canonical.matrix_world = T
        all_studs.extend(studs_local)

        stud_log_append(context, f"✔ polygon {poly.index} 完成，生成 {len(studs_local)} 根龙骨")

    return all_studs


# ============================================================
#  Utility: 安全删除 Blender 对象（不影响当前选择）
# ============================================================
def delete_objects_safely(objs):
    if not objs:
        return

    # 支持单个对象传入
    if not isinstance(objs, (list, tuple, set)):
        objs = [objs]

    # 过滤掉 None 或已被移除的对象
    objs = [o for o in objs if o and o.name in bpy.data.objects]

    if not objs:
        return

    # 记录当前选中对象集合
    prev_selection = bpy.context.selected_objects.copy()
    prev_active = bpy.context.view_layer.objects.active

    # 取消所有选中
    bpy.ops.object.select_all(action='DESELECT')

    # 选择要删除的对象
    for o in objs:
        o.select_set(True)

    # 设置 active obj（为 delete 操作所需）
    bpy.context.view_layer.objects.active = objs[0]

    # 删除
    try:
        bpy.ops.object.delete()
    except Exception as e:
        print(f"[WARN] 删除对象失败: {e}")

    # 恢复原选中状态
    bpy.ops.object.select_all(action='DESELECT')
    for o in prev_selection:
        if o and o.name in bpy.data.objects:
            o.select_set(True)
    bpy.context.view_layer.objects.active = prev_active


def draw_panels_raw_debug(context, panels_raw, name_prefix="PDBG"):
    """
    可视化 panels_raw（只画 final 边）：
      outer 阳角 = 绿
      inner 阴角 = 蓝
      non-corner = 黄
    object.name 极简格式：
      PDBG_p{poly}_e{edge}_{OUT/IN/N}_nn(x,y,z)_F
    """

    # 删除旧的 debug 对象
    for obj in list(bpy.data.objects):
        if obj.name.startswith(name_prefix):
            bpy.data.objects.remove(obj, do_unlink=True)

    # 创建集合
    col_name = f"{name_prefix}_COL"
    if col_name in bpy.data.collections:
        debug_col = bpy.data.collections[col_name]
    else:
        debug_col = bpy.data.collections.new(col_name)
        context.scene.collection.children.link(debug_col)

    # Object Color
    C_GREEN  = (0.0, 1.0, 0.0, 1.0)   # outer
    C_BLUE   = (0.0, 0.4, 1.0, 1.0)   # inner
    C_YELLOW = (1.0, 1.0, 0.0, 1.0)   # non-corner

    def mk(name, p1, p2, col):
        mesh = bpy.data.meshes.new(name+"_M")
        obj  = bpy.data.objects.new(name, mesh)
        mesh.from_pydata([p1, p2], [(0,1)], [])
        debug_col.objects.link(obj)
        obj.display_type = 'WIRE'
        obj.show_in_front = True
        obj.show_wire = True
        obj.color = col
        return obj

    # ------ 遍历每个 polygon ------
    for poly_idx, pdata in panels_raw.items():
        edges = pdata.get("edges", [])
        if not edges:
            continue

        for ei, e in enumerate(edges):
            v1f = e.get("v1_final")
            v2f = e.get("v2_final")
            is_corner = e.get("is_corner")
            ctype = e.get("corner_type")
            neigh = e.get("neighbor_normal")

            # 简短 corner 表示
            if not is_corner:
                ctag = "N"    # non-corner
                col = C_YELLOW
            else:
                if ctype == "outer":
                    ctag = "OUT"
                    col = C_GREEN
                elif ctype == "inner":
                    ctag = "IN"
                    col = C_BLUE
                else:
                    ctag = "N"
                    col = C_YELLOW

            # 简短 neighbor normal
            if neigh is None:
                nn = "None"
            else:
                nn = f"{neigh.x:.2f},{neigh.y:.2f},{neigh.z:.2f}"

            # 极简 object 名称
            name = f"{name_prefix}_p{poly_idx}_e{ei}_{ctag}_nn({nn})_F"

            mk(name, v1f, v2f, col)

    stud_log_append(context, "🎨 Debug：已绘制 final 边（无 world 边）")

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
    # 转角龙骨选定（阴 / 阳）
    # ------------------------------

    corner_inner_type: bpy.props.EnumProperty(
        name="内角龙骨",
        items=update_stud_type_enum,
    )

    corner_outer_type: bpy.props.EnumProperty(
        name="外角龙骨",
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
        # 3. 生成 panels_raw（可调试）
        panels_raw = build_corner_panels_data(context, ref_obj, offset_obj, offset_dist=offset_value)

        draw_panels_raw_debug(context, panels_raw)
        
        final_obj = build_final_offset_mesh_from_panels_raw(panels_raw)

        stud_log_append(context, "✔ 已可视化 final_obj，检查侧面20mm偏移是否正确")


        # ============================
        # DEBUG：打印 panels_raw 结构
        # ============================
        stud_log_append(context, f"\n===== DEBUG: panels_raw，共 {len(panels_raw)} 个 polygon =====")

        ref_me = ref_obj.data
        off_me = offset_obj.data
        M_ref = ref_obj.matrix_world
        M_off = offset_obj.matrix_world

        for poly in ref_me.polygons:
            pi = poly.index
            pdata = panels_raw.get(pi)

            ref_normal = (M_ref.to_3x3() @ poly.normal).normalized()
            if pi < len(off_me.polygons):
                off_normal = (M_off.to_3x3() @ off_me.polygons[pi].normal).normalized()
            else:
                off_normal = None

            stud_log_append(
                context,
                f"\n--- 面 #{pi} ---\n"
                f"ref_normal   = {ref_normal}\n"
                f"offset_normal= {off_normal}"
            )

            if not pdata:
                stud_log_append(context, "  (⚠ panels_raw 中无该面的数据)")
                continue

            is_side = pdata.get("is_side")
            stud_log_append(context, f"  is_side = {is_side}")

            edges = pdata.get("edges", [])
            stud_log_append(context, f"  edges 数量 = {len(edges)}")

            for ei, e in enumerate(edges):
                v1_world = e.get("v1_world")
                v2_world = e.get("v2_world")
                v1_final = e.get("v1_final")
                v2_final = e.get("v2_final")

                is_corner = e.get("is_corner")
                corner_type = e.get("corner_type")
                neighbor_normal = e.get("neighbor_normal")

                # dot 用来辅助看内外角（如果存在邻面法向）
                if neighbor_normal is not None:
                    try:
                        dp = ref_normal.dot(neighbor_normal)
                    except:
                        dp = "ERR"
                else:
                    dp = "N/A"

                stud_log_append(
                    context,
                    f"    Edge #{ei}:\n"
                    f"      is_corner     = {is_corner}\n"
                    f"      corner_type   = {corner_type}\n"
                    f"      neighbor_norm = {neighbor_normal}\n"
                    f"      dot(ref, neigh)= {dp}\n"
                    f"      v1_world      = {v1_world}\n"
                    f"      v2_world      = {v2_world}\n"
                    f"      v1_final      = {v1_final}\n"
                    f"      v2_final      = {v2_final}"
                )
        
        return {"FINISHED"}



# =================================================================
#  Operator：为多面参考面生成龙骨
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

        ref_obj = props.ref_obj or context.active_object
        if not ref_obj or ref_obj.type != "MESH":
            stud_log_set(context, "❌ 请选择一个 Mesh 作为参考面")
            return {"CANCELLED"}

        ok, reason = validate_reference_mesh(ref_obj)
        if not ok:
            stud_log_set(context, f"❌ 无法作为多面参考面：{reason}")
            return {"CANCELLED"}

        if not self.bypass_scale_check:
            sx, sy, sz = ref_obj.scale
            if (abs(sx - 1.0) > 1e-6) or (abs(sy - 1.0) > 1e-6) or (abs(sz - 1.0) > 1e-6):
                return bpy.ops.ifc.confirm_apply_scale(
                    'INVOKE_DEFAULT',
                    ref_obj_name=ref_obj.name
                )

        # 2. 创建 offset_obj
        offset_value = props.offset_z
        offset_obj = create_offset_object_from_ref(ref_obj, offset_value)
        if offset_obj is None:
            stud_log_set(context, "❌ 创建 offset_obj 失败")
            return {"CANCELLED"}

        stud_log_append(context, f"✔ 创建 offset_obj：{offset_obj.name}")

        # 3. 生成 panels_raw（可调试）
        panels_raw = build_corner_panels_data(context, ref_obj, offset_obj, offset_dist=offset_value)
        stud_log_append(context, "✔ 完成 corner/side 预处理")

        # 4. 生成龙骨
        final_obj = build_final_offset_mesh_from_panels_raw(panels_raw)
        
        studs = generate_studs_on_mesh(
            context, model, props,
            final_obj,
            panels_raw,
        )

        delete_objects_safely([final_obj, offset_obj])

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
        col.prop(props, "corner_inner_type", text="内角")
        col.prop(props, "corner_outer_type", text="外角")
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
