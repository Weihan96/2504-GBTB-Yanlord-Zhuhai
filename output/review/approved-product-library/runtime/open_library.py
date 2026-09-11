"""Run in Blender after opening a target IFC with Bonsai. Does not save."""
from pathlib import Path
import importlib.util
import sys
import bpy

root=Path(__file__).resolve().parent
name="highpoly_review_library"
if name in sys.modules:
    library=sys.modules[name]
    if library.bl_info["version"] < (0,8,0) or library._catalog_path != root/"catalog.json":
        raise RuntimeError("当前窗口载入的是其他版本或路径；请保存自己的工作后，在新窗口加载此图库。")
else:
    spec=importlib.util.spec_from_file_location(name,root/"addon/highpoly_review_library/__init__.py")
    library=importlib.util.module_from_spec(spec)
    sys.modules[name]=library
    spec.loader.exec_module(library)
    library.register()
for area in bpy.context.screen.areas:
    if area.type=="VIEW_3D":
        area.spaces.active.show_region_ui=True
        area.tag_redraw()
print("已加载 IFC Assets；3D 视图按 N → Assets。点击查看，拖动放置。未插入模型，未保存 IFC。")
