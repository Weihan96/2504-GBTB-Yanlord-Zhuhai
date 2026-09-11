"""Fresh task-owned IFC acceptance bootstrap; no IFC or Blend save."""
from pathlib import Path
import runpy
import bpy
import bonsai.tool as tool
import ifcopenshell.util.element

scripts = Path(__file__).resolve().parent
runpy.run_path(str(scripts / "open_approved_library_review.py"), run_name="__main__")

def show_library():
    runpy.run_path(str(scripts / "activate_approved_library.py"), run_name="__main__")
    targets = [e for e in tool.Ifc.get().by_type("IfcElement")
               if ifcopenshell.util.element.get_pset(e, "ReviewLibrarySource")]
    if targets:
        obj = tool.Ifc.get_object(targets[-1])
        for other in bpy.context.selected_objects:
            other.select_set(False)
        obj.select_set(True)
        bpy.context.view_layer.objects.active = obj
        for area in bpy.context.screen.areas:
            if area.type == "VIEW_3D":
                with bpy.context.temp_override(area=area, region=next(r for r in area.regions if r.type == "WINDOW")):
                    bpy.ops.view3d.view_selected(use_all_regions=False)
    return None

bpy.app.timers.register(show_library, first_interval=.5)
