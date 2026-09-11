"""Approved IFC Asset Browser library, with a separate read-only review array."""
bl_info = {"name": "IFC Assets", "author": "Project review tools", "version": (0, 8, 3),
           "blender": (4, 5, 0), "location": "3D View > Sidebar > Assets", "category": "3D View"}

import json
import hashlib
import math
import textwrap
from pathlib import Path
import bpy
import bpy.utils.previews
from bpy.props import StringProperty, EnumProperty, IntProperty, FloatProperty, BoolProperty
from mathutils import Matrix, Vector

ROOT = Path(__file__).resolve().parents[3]
DEFAULT_CATALOG = ROOT / "output/review/approved-product-library/catalog.json"
if DEFAULT_CATALOG.with_name("all-review-catalog.json").is_file():
    DEFAULT_CATALOG = DEFAULT_CATALOG.with_name("all-review-catalog.json")
runtime_catalog = ROOT / "output/review/approved-product-library/runtime/catalog.json"
portable_catalog = Path(__file__).resolve().parents[2] / "catalog.json"
if runtime_catalog.is_file():
    DEFAULT_CATALOG = runtime_catalog
if portable_catalog.is_file():
    DEFAULT_CATALOG = portable_catalog
if not DEFAULT_CATALOG.is_file():
    DEFAULT_CATALOG = Path(bpy.data.filepath).parent / "catalog.json"
_catalog = {"products": []}
_preview = None
_enum_items = []
_catalog_path = DEFAULT_CATALOG


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def resolve_path(value):
    path = Path(value)
    return path if path.is_absolute() else (_catalog_path.parent / path).resolve()


def entry(key):
    return next((p for p in _catalog["products"] if p["id"] == key), None)


def enum_items(self, context):
    global _enum_items
    query = self.search.strip().casefold()
    _enum_items = [(p["id"], p["name"], p.get("approval_label", "已验收"),
                    _preview[p["id"] + ":iso"].icon_id if _preview and p["id"] + ":iso" in _preview else 0, i)
                   for i, p in enumerate(_catalog["products"]) if not query or query in (p["id"] + " " + p["name"]).casefold()]
    return _enum_items


def asset_library_enum(self, context):
    return [("LOCAL", "IFC 单品卡片", "仅内存缩略图，模型直接来自 IFC", 1)]


def asset_library_value(self):
    return 1


def load_catalog(path):
    global _catalog, _catalog_path, _preview
    candidate_path = Path(path).resolve()
    data = json.loads(candidate_path.read_text())
    def candidate_asset(value):
        asset = Path(value)
        return asset if asset.is_absolute() else (candidate_path.parent / asset).resolve()
    assert data.get("schema_version") in (2, 3, 4), "不支持的图库目录版本"
    assert len({p["id"] for p in data["products"]}) == len(data["products"])
    for product in data["products"]:
        if data["schema_version"] in (3, 4):
            from .catalog_rules import validate_entry
            validate_entry(product, candidate_asset, digest)
            continue
        assert product["single_product_approval_status"] == "approved"
        scene_status = product["scene_approval_status"]
        assert scene_status in ("approved", "pending")
        assert product["approval_label"] == {"approved": "单品及场景已验收", "pending": "单品已通过、场景待验收"}[scene_status]
        assert product["display_front_basis"] == {"approved": "approved_scene_camera", "pending": "candidate_scene_camera_pending_review"}[scene_status]
        if scene_status == "pending":
            authority = product["library_authorization"]
            authority_path = candidate_asset(authority["path"])
            assert digest(authority_path) == authority["sha256"]
            inclusion = json.loads(authority_path.read_text())
            assert inclusion["library_inclusion_authorized"] is True
            assert product["id"] in inclusion["products"]
            assert inclusion["scene_approval_status"] == "pending"
        assert product["package_validation"] == "pass"
        assert candidate_asset(product["ifc_path"]).is_file()
        for view in ("plan", "front", "side", "iso"):
            assert candidate_asset(product["previews"][view]).is_file()
    candidate_preview = bpy.utils.previews.new()
    try:
        for p in data["products"]:
            for view, path in p["previews"].items():
                candidate_preview.load(p["id"] + ":" + view, str(candidate_asset(path)), "IMAGE")
    except Exception:
        bpy.utils.previews.remove(candidate_preview)
        raise
    if _preview:
        bpy.utils.previews.remove(_preview)
    _preview = candidate_preview
    _catalog = data
    _catalog_path = candidate_path


def review_scene(context):
    if context.scene.get("review_library_display") and context.scene.get("display_schema_version") == 2:
        return context.scene
    scene = bpy.data.scenes.new("单品审核库 · XY 阵列")
    scene["review_library_display"] = True
    scene["display_schema_version"] = 2
    scene["ifc_write_allowed"] = False
    scene.unit_settings.system = "METRIC"
    scene.unit_settings.length_unit = "METERS"
    context.window.scene = scene
    return scene


def body_collection(scene, product):
    """All Body representations are read, in SI units. No IFC is changed."""
    import ifcopenshell
    import ifcopenshell.geom
    import ifcopenshell.util.placement
    import ifcopenshell.util.unit
    existing = next((c for c in scene.collection.children if c.get("review_product_id") == product["id"]), None)
    if existing:
        if existing.get("ifc_sha256") != product["ifc_sha256"]:
            raise RuntimeError("已载入的展示副本版本不同；请保留旧场景并新建验收场景")
        existing["scene_approval_status"] = product["scene_approval_status"]
        return existing
    path = resolve_path(product["ifc_path"])
    if digest(path) != product["ifc_sha256"]:
        raise RuntimeError("单品 IFC 哈希已变化，拒绝载入未验证版本")
    model = ifcopenshell.open(str(path))
    target = model.by_guid(product["global_id"])
    assert len([e for e in model.by_type("IfcElement") if e.Representation and not e.is_a("IfcOpeningElement")]) == 1
    assert not model.by_type("IfcAnnotation"), "不是纯单品 IFC"
    unit = ifcopenshell.util.unit.calculate_unit_scale(model)
    original_matrix = ifcopenshell.util.placement.get_local_placement(target.ObjectPlacement)
    original_matrix[:3, 3] *= unit
    collection = bpy.data.collections.new(product["name"])
    collection["review_product_id"] = product["id"]
    collection["ifc_path"] = str(path)
    collection["ifc_sha256"] = product["ifc_sha256"]
    collection["global_id"] = product["global_id"]
    collection["scene_approval_status"] = product["scene_approval_status"]
    collection["product_local_to_project_matrix_m"] = json.dumps(original_matrix.tolist())
    front = product["display_front_normal_project"]
    assert abs(front[2]) < .01
    yaw = -math.pi / 2 - math.atan2(front[1], front[0])
    # Preserve project vertical, then yaw only to face the recipe Front toward
    # -Y. Pending scene cameras remain explicitly unapproved in the catalogue.
    project_rotation = Matrix(original_matrix.tolist()).to_3x3().to_4x4()
    display_rotation = Matrix.Rotation(yaw, 4, "Z") @ project_rotation
    collection["product_local_to_display_orientation"] = json.dumps([list(r) for r in display_rotation])
    collection["display_only"] = True
    scene.collection.children.link(collection)
    settings = ifcopenshell.geom.settings()
    # Curves in Body are retained as edges; triangle bodies retain faces.
    settings.set("dimensionality", ifcopenshell.ifcopenshell_wrapper.CURVES_SURFACES_AND_SOLIDS)
    points = []
    try:
        for i, rep in enumerate(target.Representation.Representations):
            if rep.RepresentationIdentifier != "Body" or rep.ContextOfItems.TargetView != "MODEL_VIEW":
                continue
            shape = ifcopenshell.geom.create_shape(settings, rep)
            geometry = getattr(shape, "geometry", shape)
            vertices = list(zip(*[iter(geometry.verts)] * 3))
            faces = list(zip(*[iter(geometry.faces)] * 3))
            edges = [] if faces else list(zip(*[iter(geometry.edges)] * 2))
            mesh = bpy.data.meshes.new(f"{product['id']} Body {i}")
            mesh.from_pydata(vertices, edges, faces)
            mesh.update()
            for polygon in mesh.polygons:
                polygon.use_smooth = True
            if hasattr(mesh, "set_sharp_from_angle"):
                mesh.set_sharp_from_angle(angle=.7)
            obj = bpy.data.objects.new(mesh.name, mesh)
            obj["review_only"] = True
            obj["source_ifc_global_id"] = product["global_id"]
            obj["source_body_representation_index"] = i
            collection.objects.link(obj)
            material = bpy.data.materials.get("单品库 · 中性灰") or bpy.data.materials.new("单品库 · 中性灰")
            material.diffuse_color = (.58, .64, .69, 1)
            obj.data.materials.append(material)
            points.extend(display_rotation @ Vector(v) for v in vertices)
        assert points, "产品 Body 无可视几何"
        bounds = [[min(p[i] for p in points) for i in range(3)], [max(p[i] for p in points) for i in range(3)]]
        collection["display_oriented_bounds_m"] = json.dumps(bounds)
        collection["body_representation_count"] = sum(r.RepresentationIdentifier == "Body" and r.ContextOfItems.TargetView == "MODEL_VIEW" for r in target.Representation.Representations)
        assert digest(path) == product["ifc_sha256"]
        return collection
    except Exception:
        # Remove only newly created display objects, never source or user data.
        for obj in list(collection.objects):
            bpy.data.objects.remove(obj, do_unlink=True)
        bpy.data.collections.remove(collection)
        raise


def place(collection, x, y):
    minimum, maximum = json.loads(collection["display_oriented_bounds_m"])
    translation = Vector((x - (minimum[0] + maximum[0]) / 2, y - (minimum[1] + maximum[1]) / 2, -minimum[2]))
    rotation = Matrix(json.loads(collection["product_local_to_display_orientation"]))
    for obj in collection.objects:
        if obj.get("review_only"):
            obj.matrix_world = Matrix.Translation(translation) @ rotation
    original = Matrix(json.loads(collection["product_local_to_project_matrix_m"]))
    project_to_display = Matrix.Translation(translation) @ rotation @ original.inverted()
    assert (project_to_display.to_3x3() @ Vector((0, 0, 1)) - Vector((0, 0, 1))).length < 1e-5
    collection["project_to_display_matrix_m"] = json.dumps([list(row) for row in project_to_display])
    collection["display_slot_xy_m"] = (x, y)
    label = next((o for o in collection.objects if o.get("review_label")), None)
    if label is None:
        curve = bpy.data.curves.new(collection.name + " label", "FONT")
        label = bpy.data.objects.new(curve.name, curve)
        label["review_label"] = True
        collection.objects.link(label)
    key = collection["review_product_id"]
    label.data.body = {"miamisoft-e09": "MIAMI E09", "gessi316-54294": "GESSI 54294",
                       "geberit-duofix-sigma-224-212": "DUOFIX", "gessi316-54145": "GESSI 54145",
                       "geberit-154-446-ks-1": "SHOWER CHANNEL", "falper-sorgente": "FALPER WFB",
                       "street-h": "STREET-H SUPPORT"}.get(key, key.upper())
    if collection["scene_approval_status"] == "pending":
        label.data.body += "\nScene pending"
    label.data.size = .16
    label.data.align_x = "CENTER"
    label.location = (x, y - (maximum[1] - minimum[1]) / 2 - .30, .01)


def arrange_all(context, columns=3, gap=1.):
    scene = review_scene(context)
    collections = [body_collection(scene, p) for p in _catalog["products"]]
    if not collections:
        return []
    bounds = [json.loads(c["display_oriented_bounds_m"]) for c in collections]
    pitch_x = max(b[1][0] - b[0][0] for b in bounds) + gap
    pitch_y = max(b[1][1] - b[0][1] for b in bounds) + gap
    for i, collection in enumerate(collections):
        place(collection, (i % columns) * pitch_x, (i // columns) * pitch_y)
    scene["review_catalog_path"] = str(_catalog_path)
    scene["array_columns"] = columns
    scene["array_pitch_m"] = (pitch_x, pitch_y)
    return collections


class REVIEWLIB_Properties(bpy.types.PropertyGroup):
    asset_library: EnumProperty(items=asset_library_enum, get=asset_library_value)
    assets: bpy.props.CollectionProperty(type=bpy.types.AssetHandle)
    active_asset: IntProperty(default=0)
    show_catalog_settings: BoolProperty(name="目录设置", default=False)
    catalog_path: StringProperty(name="目录文件", subtype="FILE_PATH", default=str(DEFAULT_CATALOG))
    search: StringProperty(name="搜索")
    product: EnumProperty(name="单品", items=enum_items)
    columns: IntProperty(name="每行数量", default=3, min=1, max=12)
    gap: FloatProperty(name="单品净间距", default=1., min=.1, subtype="DISTANCE", unit="LENGTH")


class REVIEWLIB_OT_Refresh(bpy.types.Operator):
    bl_idname = "review_library.refresh"
    bl_label = "刷新单品目录"
    def execute(self, context):
        try:
            load_catalog(bpy.path.abspath(context.scene.review_library.catalog_path))
            from . import asset_cards
            asset_cards.rebuild(__import__(__name__, fromlist=["*"]))
            if _catalog["products"]:
                context.scene.review_library.product = _catalog["products"][0]["id"]
            return {"FINISHED"}
        except Exception as e:
            self.report({"ERROR"}, str(e)); return {"CANCELLED"}


class REVIEWLIB_OT_Load(bpy.types.Operator):
    bl_idname = "review_library.load"
    bl_label = "载入此单品展示副本"
    bl_options = {"INTERNAL"}  # legacy read-only array worker, not an insertion UI
    bl_description = "只创建 Blender 展示网格；不修改 IFC，不写回项目坐标"
    def execute(self, context):
        product = entry(context.scene.review_library.product)
        if not product:
            return {"CANCELLED"}
        try:
            scene = review_scene(context)
            collection = body_collection(scene, product)
            if "display_slot_xy_m" not in collection:
                # A new item gets a free slot; loading an existing grid member
                # must never silently move it back over the first product.
                occupied = [c for c in scene.collection.children if "display_slot_xy_m" in c]
                bounds = json.loads(collection["display_oriented_bounds_m"])
                half_width = (bounds[1][0] - bounds[0][0]) / 2
                right_edge = max((float(c["display_slot_xy_m"][0]) +
                    (json.loads(c["display_oriented_bounds_m"])[1][0] -
                     json.loads(c["display_oriented_bounds_m"])[0][0]) / 2 for c in occupied), default=0.)
                x = right_edge + scene.review_library.gap + half_width if occupied else 0.
                place(collection, x, 0.)
            scene.review_library.product = product["id"]
            return {"FINISHED"}
        except Exception as e:
            self.report({"ERROR"}, str(e)); return {"CANCELLED"}


class REVIEWLIB_OT_Arrange(bpy.types.Operator):
    bl_idname = "review_library.arrange"
    bl_label = "全部载入并排列到 XY 平面"
    bl_options = {"INTERNAL"}
    def execute(self, context):
        props = context.scene.review_library
        try:
            arrange_all(context, props.columns, props.gap)
            return {"FINISHED"}
        except Exception as e:
            self.report({"ERROR"}, str(e)); return {"CANCELLED"}


class REVIEWLIB_OT_Image(bpy.types.Operator):
    bl_idname = "review_library.image"
    bl_label = "查看大图"
    view: StringProperty()
    def execute(self, context):
        product = entry(context.scene.review_library.product)
        if not product:
            return {"CANCELLED"}
        path = resolve_path(product["previews"][self.view])
        area = next((a for a in context.screen.areas if a.type == "IMAGE_EDITOR"), None)
        if area is None:
            viewport = next((a for a in context.screen.areas if a.type == "VIEW_3D"), None)
            if viewport is None:
                return {"CANCELLED"}
            previous = {a.as_pointer() for a in context.screen.areas}
            with context.temp_override(area=viewport):
                bpy.ops.screen.area_split(direction="VERTICAL", factor=.65)
            area = next(a for a in context.screen.areas if a.as_pointer() not in previous)
            area.type = "IMAGE_EDITOR"
        area.spaces.active.image = bpy.data.images.load(str(path), check_existing=True)
        with context.temp_override(area=area, region=next(r for r in area.regions if r.type == "WINDOW")):
            bpy.ops.image.view_all(fit_view=True)
        return {"FINISHED"}


CLASSES = (REVIEWLIB_Properties, REVIEWLIB_OT_Refresh, REVIEWLIB_OT_Load,
           REVIEWLIB_OT_Arrange, REVIEWLIB_OT_Image)


def register():
    for cls in CLASSES:
        bpy.utils.register_class(cls)
    bpy.types.Scene.review_library = bpy.props.PointerProperty(type=REVIEWLIB_Properties)
    if DEFAULT_CATALOG.is_file():
        load_catalog(DEFAULT_CATALOG)
    from . import native_assets
    native_assets.register()
    from . import drawing_integration
    drawing_integration.register()
    from . import asset_ui
    asset_ui.register()


def unregister():
    global _preview
    from . import asset_ui
    asset_ui.unregister()
    from . import drawing_integration
    drawing_integration.unregister()
    from . import native_assets
    native_assets.unregister()
    del bpy.types.Scene.review_library
    for cls in reversed(CLASSES):
        bpy.utils.unregister_class(cls)
    if _preview:
        bpy.utils.previews.remove(_preview)
        _preview = None
