"""Run from Blender's Text Editor to show the repo-local product library.

No addon installation, preferences save, IFC write or Scene replacement.
The current Bonsai project stays active; dragging creates an in-memory element.
"""
import bpy
import importlib.util
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
path = ROOT / "pipeline/addons/highpoly_review_library/__init__.py"
name = "highpoly_review_library"
if name in sys.modules:
    library = sys.modules[name]
    if library.bl_info["version"] < (0, 8, 0):
        raise RuntimeError("旧版展示库已加载：请在保存自己的工作后重启 Blender，再运行此脚本。不要热替换原生资产面板。")
else:
    spec = importlib.util.spec_from_file_location(name, path)
    library = importlib.util.module_from_spec(spec)
    sys.modules[name] = library
    spec.loader.exec_module(library)
    library.register()
props = bpy.context.scene.review_library
if library._catalog["products"] and not props.is_property_set("active_asset"):
    props.active_asset = 0
    props.product = min(library._catalog["products"], key=lambda p: p["name"].casefold())["id"]
# A freshly opened window has not built the new sidebar's panel categories
# yet. Materialize them before asking RNA to activate our category.
for area in bpy.context.screen.areas:
    if area.type == "VIEW_3D":
        area.spaces.active.show_region_ui = True
        area.tag_redraw()
bpy.ops.wm.redraw_timer(type="DRAW_WIN_SWAP", iterations=2)
for area in bpy.context.screen.areas:
    if area.type == "VIEW_3D":
        area.spaces.active.show_region_ui = True
        area.spaces.active.overlay.show_extras = False
        for region in area.regions:
            if region.type == "UI":
                # This RNA setter needs the owning UI region in context;
                # Text Editor / bridge context otherwise reports read-only.
                with bpy.context.temp_override(area=area, region=region):
                    region.active_panel_category = "Assets"
        area.tag_redraw()
print(f"已加载 {len(library._catalog['products'])} 个可拖放 IFC 资产；未插入模型、未保存 IFC。")
