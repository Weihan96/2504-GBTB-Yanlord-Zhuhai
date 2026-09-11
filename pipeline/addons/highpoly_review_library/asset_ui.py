"""Compact asset browser; review details belong to an on-demand popup."""
from collections import Counter
import sys
import textwrap
import bpy


def library():
    return sys.modules[__package__]


def open_details(context):
    if library().entry(context.scene.review_library.product):
        return bpy.ops.review_library.details('INVOKE_DEFAULT')
    return {'CANCELLED'}


def draw_details(layout, product, lib):
    for line in textwrap.wrap(product['name'], width=42):
        layout.label(text=line)
    icon = 'CHECKMARK' if product['scene_approval_status'] == 'approved' else 'TIME'
    layout.label(text=product['approval_label'], icon=icon)
    if product.get('insertion_content') == 'body_only':
        layout.label(text='仅含 3D · 二维图纸尚未接入', icon='INFO')
    elif product.get('insertion_content') == 'review_body_views':
        layout.label(text='3D / Plan / Front / Side · Body', icon='OUTLINER_DATA_MESH')
    for line in textwrap.wrap(product.get('source_label', ''), width=28):
        layout.label(text=line)
    layout.separator()
    row = layout.row()
    for view in ('plan', 'front', 'side'):
        col = row.column(align=True)
        key = product['id'] + ':' + view
        if lib._preview and key in lib._preview:
            col.template_icon(icon_value=lib._preview[key].icon_id, scale=4.0)
        op = col.operator('review_library.image', text=view.title())
        op.view = view


class REVIEWLIB_OT_Details(bpy.types.Operator):
    bl_idname = 'review_library.details'
    bl_label = 'Asset Details'
    bl_description = '查看当前资产的名称、验收状态、来源和三视图；不会插入或保存 IFC'
    bl_options = {'INTERNAL'}

    @classmethod
    def poll(cls, context):
        return bool(library().entry(context.scene.review_library.product))

    def invoke(self, context, event):
        return context.window_manager.invoke_popup(self, width=360)

    def draw(self, context):
        lib = library()
        product = lib.entry(context.scene.review_library.product)
        if product:
            draw_details(self.layout, product, lib)

    def execute(self, context):
        return {'FINISHED'}


class REVIEWLIB_PT_Library(bpy.types.Panel):
    bl_label = 'Assets'
    bl_idname = 'REVIEWLIB_PT_library'
    bl_space_type = 'VIEW_3D'
    bl_region_type = 'UI'
    bl_category = 'Assets'

    def draw(self, context):
        lib = library()
        layout = self.layout
        props = context.scene.review_library
        header = layout.row(align=True)
        header.label(text=f"{len(lib._catalog['products'])} assets")
        header.operator('review_library.details', text='', icon='INFO')
        header.prop(props, 'show_catalog_settings', text='', icon='FILE_FOLDER')
        header.operator('review_library.refresh', text='', icon='FILE_REFRESH')
        counts = Counter(p.get('category_label', p['approval_label']) for p in lib._catalog['products'])
        for label, count in counts.items():
            layout.label(text=f'{label}：{count}')
        if props.show_catalog_settings:
            layout.prop(props, 'catalog_path', text='')
        layout.label(text='点击查看 · 拖动放置', icon='ASSET_MANAGER')
        if not lib._catalog['products']:
            layout.label(text='暂无可用资产', icon='INFO')
            return
        layout.template_asset_view('approved_ifc_products', props, 'asset_library',
            props, 'assets', props, 'active_asset', filter_id_types={'filter_group'},
            display_options={'NO_LIBRARY'}, activate_operator='review_library.select_asset',
            drag_operator='review_library.drag_ifc')
        if context.scene.get('review_library_display'):
            layout.label(text='当前为审核阵列；请切回 IFC 项目', icon='INFO')
        else:
            layout.label(text='仅改内存 · Ctrl+S 使用 Bonsai 保存', icon='FILE_TICK')


CLASSES = (REVIEWLIB_OT_Details, REVIEWLIB_PT_Library)


def register():
    for cls in CLASSES:
        bpy.utils.register_class(cls)


def unregister():
    for cls in reversed(CLASSES):
        bpy.utils.unregister_class(cls)
