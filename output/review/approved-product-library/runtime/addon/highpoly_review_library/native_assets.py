"""Native Asset Browser catalogue; release-to-place is one Bonsai transaction."""
from pathlib import Path
import json
import sys
import logging
import itertools
import gpu
from gpu_extras.batch import batch_for_shader
import ifcopenshell.geom
import bpy
import bonsai.tool as tool
from bpy.props import StringProperty, FloatVectorProperty
from mathutils import Matrix, Vector
from bpy_extras import view3d_utils
from bonsai.bim.module.project.operator import AppendLibraryElement
from bonsai.bim import import_ifc
from .ifc_asset import prepare_source, insert_prepared
from .loading_feedback import LoadingFeedback, notify

LIBRARY_NAME = "高模 IFC 单品审核库"


class PlacementPreview:
    """Viewport-only wire box. Never creates a Blender object or an IFC entity."""
    def __init__(self, bounds):
        self.bounds = bounds
        self.point = None
        self.area = None
        self.outline_visible = True
        self.shader = gpu.shader.from_builtin("UNIFORM_COLOR")
        self.handle = bpy.types.SpaceView3D.draw_handler_add(self.draw, (), "WINDOW", "POST_VIEW")

    def vertices(self):
        low, high = self.bounds
        corners = [Vector(tuple((high if bit else low)[i] for i, bit in enumerate(bits))) + self.point
                   for bits in itertools.product((0, 1), repeat=3)]
        return [corners[i] for a in range(8) for b in range(a+1, 8)
                if (a ^ b) in (1, 2, 4) for i in (a, b)]

    def update(self, target, point):
        previous = self.area
        self.area = target[0] if target else None
        self.point = Vector(point) if target else None
        for area in (previous, self.area):
            if area:
                area.tag_redraw()

    def draw(self):
        if not self.outline_visible or self.point is None or bpy.context.area != self.area:
            return
        depth = gpu.state.depth_test_get()
        width = gpu.state.line_width_get()
        try:
            gpu.state.depth_test_set("NONE")
            gpu.state.line_width_set(2.0)
            self.shader.bind()
            self.shader.uniform_float("color", (0.1, 0.65, 1.0, 1.0))
            batch_for_shader(self.shader, "LINES", {"pos": self.vertices()}).draw(self.shader)
        finally:
            gpu.state.depth_test_set(depth)
            gpu.state.line_width_set(width)

    def close(self):
        if self.handle is not None:
            bpy.types.SpaceView3D.draw_handler_remove(self.handle, "WINDOW")
            self.handle = None
        self.update(None, None)


class ProductImporter:
    # Reuse the installed Bonsai importer methods, without constructing an RNA
    # Operator directly (Blender only constructs those via bpy.ops).
    import_materials = AppendLibraryElement.import_materials
    import_styles = AppendLibraryElement.import_styles
    import_material_styles = AppendLibraryElement.import_material_styles

    def import_product_from_ifc(self, element, context, product_id=None):
        self.file = tool.Ifc.get()
        settings = import_ifc.IfcImportSettings.factory(context, tool.Ifc.get_path(), logging.getLogger("ImportIFC"))
        importer = import_ifc.IfcImporter(settings)
        importer.file = self.file
        # Each approved product legitimately has dedicated drawing contexts.
        # Native full-file import's >100-context merge can invalidate settings;
        # refresh the Loader directly, without rewriting any existing contexts.
        tool.Loader.load_settings()
        importer.material_creator.load_existing_materials()
        self.import_materials(element, importer)
        self.import_styles(element, importer)
        notify(product_id, .65, "材质就绪，正在生成几何")
        # A parent context iterator can emit Body + all three drawing shapes
        # for the same product. Import exactly one primary 3D representation;
        # all approved representations remain in the IFC product unchanged.
        bodies = [r for r in element.Representation.Representations if r.RepresentationIdentifier == "Body"]
        body = next((r for r in bodies if getattr(r.ContextOfItems, "TargetView", None) == "MODEL_VIEW"), bodies[0])
        geom_settings = ifcopenshell.geom.settings()
        geom_settings.set("dimensionality", ifcopenshell.ifcopenshell_wrapper.CURVES_SURFACES_AND_SOLIDS)
        shape = ifcopenshell.geom.create_shape(geom_settings, element, body)
        notify(product_id, .85, "几何就绪，正在建立场景对象")
        importer.create_product(element, shape)
        importer.place_objects_in_collections()
        notify(product_id, .95, "正在完成 Bonsai 关联")


def library():
    return sys.modules[__package__]


def asset_product(context):
    from .asset_cards import product_for_id
    asset = getattr(context, "asset", None)
    if asset:
        return product_for_id(getattr(asset, "local_id", None), library())
    # template_asset_view supplies active_file, not Asset Browser's asset.
    active_file = getattr(context, "active_file", None)
    if active_file:
        path = active_file.relative_path.replace("\\", "/")
        if path.startswith("Collection/"):
            return product_for_id(bpy.data.collections.get(path[len("Collection/"):]), library())
    return None


def drop_matrix(product, point):
    record = json.loads((library()._catalog_path.parent / "native-assets" / "placements.json").read_text())
    if record[product["id"]]["source_sha256"] != product["ifc_sha256"]:
        raise ValueError("资产摆放数据与单品 IFC 版本不一致，请重新构建资产库")
    matrix = Matrix(record[product["id"]]["local_to_asset_m"])
    matrix.translation += Vector(point)
    return matrix


class REVIEWLIB_OT_InsertIFC(bpy.types.Operator, tool.Ifc.Operator):
    bl_idname = "review_library.insert_ifc"
    bl_label = "放置 IFC 单品"
    bl_options = {"REGISTER", "UNDO"}
    product_id: StringProperty()
    # Absolute IFC project placement, not the browser thumbnail/display slot.
    project_matrix: FloatVectorProperty(size=16)

    @classmethod
    def poll(cls, context):
        return bool(tool.Ifc.get()) and not context.scene.get("review_library_display", False)

    def _execute(self, context):
        product = library().entry(self.product_id)
        if not product:
            raise ValueError("未知的审核单品")
        target = tool.Ifc.get()
        container = tool.Root.get_default_container()
        if container is None:
            raise ValueError("请先在 Bonsai 选择默认楼层")
        notify(self.product_id, .05, "正在读取单品 IFC")
        source, source_product = prepare_source(library().resolve_path(product["ifc_path"]),
            product["ifc_sha256"], product["global_id"], target,
            require_drawings=product.get("insertion_content", "approved_3d_2d") == "approved_3d_2d")
        notify(self.product_id, .25, "单品校验通过，正在接入 IFC")
        rows = [list(self.project_matrix[i:i+4]) for i in range(0, 16, 4)]
        element = insert_prepared(target, source, source_product, rows, container, {
            "ProductId": product["id"], "SourceGlobalId": product["global_id"],
            "SourceSHA256": product["ifc_sha256"], "SceneApproval": product["scene_approval_status"],
            "SingleProductApproval": product["single_product_approval_status"],
            "ReviewCategory": product.get("review_category", "approved"),
            "InsertionContent": product.get("insertion_content", "approved_3d_2d")})
        notify(self.product_id, .5, "IFC 接入完成，正在载入材质")
        importer = ProductImporter()
        importer.import_product_from_ifc(element, context, self.product_id)
        obj = tool.Ifc.get_object(element)
        if obj is None:
            raise RuntimeError("IFC 已创建，但 Bonsai 几何导入失败；请撤销本次放置")
        for selected in context.selected_objects:
            selected.select_set(False)
        obj.select_set(True)
        context.view_layer.objects.active = obj
        context.scene.review_library.product = product["id"]
        self.report({"INFO"}, "已加入当前 IFC 内存；Ctrl+S 使用 Bonsai 保存")


class REVIEWLIB_OT_DragIFC(bpy.types.Operator):
    bl_idname = "review_library.drag_ifc"
    bl_label = "拖放 IFC 单品"
    bl_options = {"INTERNAL", "UNDO"}

    @classmethod
    def poll(cls, context):
        return context.area is not None and bool(asset_product(context))

    def invoke(self, context, event):
        self.product_id = asset_product(context)["id"]
        self._target = None
        self._point = None
        if not tool.Ifc.get() or context.scene.get("review_library_display"):
            self.report({"WARNING"}, "请切换到已载入 IFC 的项目 Scene，再拖入单品")
            return {"CANCELLED"}
        if not tool.Root.get_default_container():
            self.report({"WARNING"}, "请先在 Bonsai 选择默认楼层")
            return {"CANCELLED"}
        try:
            product = library().entry(self.product_id)
            drop_matrix(product, (0, 0, 0))  # Fail before entering modal if stale.
            records = json.loads((library()._catalog_path.parent / "native-assets/placements.json").read_text())
            self._preview = PlacementPreview(records[self.product_id]["asset_bounds_m"])
        except Exception as error:
            self.report({"ERROR"}, str(error))
            return {"CANCELLED"}
        context.window.cursor_modal_set("CROSSHAIR")
        context.workspace.status_text_set("拖到三维视图放置 · 左键释放确认 · Esc / 右键取消 · 使用默认楼层标高")
        context.window_manager.modal_handler_add(self)
        return {"RUNNING_MODAL"}

    def finish(self, context):
        preview = getattr(self, "_preview", None)
        if preview:
            preview.close()
        context.window.cursor_modal_restore()
        context.workspace.status_text_set(None)

    def modal(self, context, event):
        try:
            return self._modal(context, event)
        except Exception as error:
            self.finish(context)
            self.report({"ERROR"}, str(error))
            return {"CANCELLED"}

    def cancel(self, context):
        self.finish(context)

    def _modal(self, context, event):
        if event.type in {"ESC", "RIGHTMOUSE", "WINDOW_DEACTIVATE"}:
            self.finish(context)
            return {"CANCELLED"}
        self._target = None
        for area in context.screen.areas:
            if area.type != "VIEW_3D":
                continue
            # Sidebar/toolbar regions overlap WINDOW in Blender's coordinates.
            # Releasing back on the asset grid must cancel, not place behind it.
            if any(r.type != "WINDOW" and r.width > 1 and r.height > 1 and
                   r.x <= event.mouse_x < r.x+r.width and r.y <= event.mouse_y < r.y+r.height
                   for r in area.regions):
                continue
            region = next((r for r in area.regions if r.type == "WINDOW" and
                r.x <= event.mouse_x < r.x+r.width and r.y <= event.mouse_y < r.y+r.height), None)
            if region is None:
                continue
            rv3d = area.spaces.active.region_3d
            xy = (event.mouse_x-region.x, event.mouse_y-region.y)
            ray = view3d_utils.region_2d_to_vector_3d(region, rv3d, xy)
            origin = view3d_utils.region_2d_to_origin_3d(region, rv3d, xy)
            plane_z = tool.Root.get_default_container_elevation()
            if abs(ray.z) < 1e-7:
                continue
            distance = (plane_z-origin.z)/ray.z
            if rv3d.is_perspective and distance < 0:
                continue
            self._point = origin + ray*distance
            self._target = (area, region)
            context.workspace.status_text_set(f"{self.product_id} · 放置点 {self._point.x:.3f}, {self._point.y:.3f}, {self._point.z:.3f} m · 释放确认 / Esc 取消")
            break
        if getattr(self, "_preview", None):
            self._preview.update(self._target, self._point)
        if event.type == "LEFTMOUSE" and event.value == "RELEASE":
            if self._target is None:
                self.finish(context)
                return {"CANCELLED"}
            area, region = self._target
            matrix = drop_matrix(library().entry(self.product_id), self._point)
            # Surveyor is Bonsai's own Blender-offset -> project-coordinate path.
            probe = bpy.data.objects.new("IFC placement coordinate probe", None)
            try:
                probe.matrix_world = matrix
                project_matrix = tool.Surveyor.get_absolute_matrix(probe)
            finally:
                bpy.data.objects.remove(probe)
            try:
                with context.temp_override(area=area, region=region):
                    with LoadingFeedback(self._preview, self.product_id, context) as loading:
                        loading.advance(0., "正在加载单品")
                        result = bpy.ops.review_library.insert_ifc(product_id=self.product_id,
                            project_matrix=[float(v) for row in project_matrix for v in row])
                        if result != {"FINISHED"}:
                            raise RuntimeError("单品载入未完成")
                        # Only show full height after the native transaction returns.
                        loading.advance(1., "单品加载完成")
            finally:
                self.finish(context)
            return {"FINISHED"}
        return {"RUNNING_MODAL"}


class REVIEWLIB_OT_SelectAsset(bpy.types.Operator):
    bl_idname = "review_library.select_asset"
    bl_label = "查看单品图"
    bl_description = "点击松手查看名称、验收状态和三视图；按住拖动放置 IFC"

    def invoke(self, context, event):
        result = self.execute(context)
        # The native preview list activates its row on LEFTMOUSE PRESS,
        # before Blender can recognise CLICK_DRAG. Opening a modal popup on
        # that press steals the rest of the gesture from the drag operator.
        # A completed click invokes us again on RELEASE (or CLICK). Native
        # custom dragging suppresses activation, so its release stays a drop.
        if event.type == 'LEFTMOUSE' and event.value in {'RELEASE', 'CLICK'}:
            from .asset_ui import open_details
            open_details(context)
        return result

    def execute(self, context):
        product = asset_product(context)
        if product:
            context.scene.review_library.product = product["id"]
        return {"FINISHED"}


CLASSES = (REVIEWLIB_OT_InsertIFC, REVIEWLIB_OT_DragIFC, REVIEWLIB_OT_SelectAsset)


def register():
    for cls in CLASSES:
        bpy.utils.register_class(cls)
    # Blender's supported custom drag hook belongs to template_asset_view.
    # A File Browser keymap cannot intercept native thumbnail drag buttons.
    # Do not override Blender's global collection/object drop operators.
    from . import asset_cards
    asset_cards.register(library())


def unregister():
    from . import asset_cards
    asset_cards.unregister()
    for cls in reversed(CLASSES):
        bpy.utils.unregister_class(cls)
