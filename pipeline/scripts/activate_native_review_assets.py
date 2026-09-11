"""Session-local UI upgrade. Preserve existing IFC memory and review Scene.

Redirect this task's package inspection session to a copied test file, so a
user's Ctrl+S during acceptance cannot overwrite an approved library package.
Never use this script to redirect a normal project session.
"""
from pathlib import Path
import importlib.util
import sys
import json
import tempfile
import shutil
import bpy
import bonsai.tool as tool

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "output/review/approved-product-library"
source_path = Path(tool.Ifc.get_path()).resolve()
assert source_path == (ROOT / "output/review/highpoly-types/bed01/product-package/BED01-product.ifc").resolve(), source_path
assert len(tool.Ifc.get().by_type("IfcElement")) == 1
original_scene = next(s for s in bpy.data.scenes if not s.get("review_library_display") and any(tool.Ifc.get_entity(o) for o in s.objects))
bpy.context.window.scene = original_scene
fixture = Path(tempfile.mkdtemp(prefix="native-library-user-review-")) / "library-review.ifc"
shutil.copy2(source_path, fixture)
tool.Ifc.set_path(str(fixture))
old = sys.modules.get("highpoly_review_library")
if old:
    old.unregister()
for name in list(sys.modules):
    if name == "highpoly_review_library" or name.startswith("highpoly_review_library."):
        del sys.modules[name]
spec = importlib.util.spec_from_file_location("highpoly_review_library", ROOT / "pipeline/addons/highpoly_review_library/__init__.py")
addon = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = addon
spec.loader.exec_module(addon)
addon.register()
storey = tool.Ifc.get().by_type("IfcBuildingStorey")[0]
tool.Spatial.set_default_container(storey)
viewport = next(a for a in bpy.context.screen.areas if a.type == "VIEW_3D")
with bpy.context.temp_override(area=viewport, region=next(r for r in viewport.regions if r.type == "WINDOW")):
    bpy.ops.view3d.view_all(center=False)
viewport.spaces.active.show_region_ui = True
bpy.ops.wm.redraw_timer(type="DRAW_WIN_SWAP", iterations=2)
with bpy.context.temp_override(area=viewport, region=next(r for r in viewport.regions if r.type == "UI")):
    bpy.context.region.active_panel_category = "单品库"
(OUT / "native-live-session.json").write_text(json.dumps({"fixture": str(fixture),
    "source_package": str(source_path), "original_scene": original_scene.name,
    "source_sha256": addon.digest(source_path), "fixture_initial_sha256": addon.digest(fixture)}, indent=2)+"\n")
print("NATIVE_LIBRARY_UI_READY", fixture)
