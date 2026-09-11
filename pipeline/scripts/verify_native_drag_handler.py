"""Exercise the production drag handler with synthetic events in the owned UI.

This is NOT a physical-pointer/keymap test. Every insertion is immediately
undone; only a fixture named library-review.ifc is accepted.
"""
from pathlib import Path
from types import SimpleNamespace
import json
import bpy
import bonsai.tool as tool
import highpoly_review_library as library
from highpoly_review_library.native_assets import REVIEWLIB_OT_DragIFC, PlacementPreview
from mathutils import Vector

assert Path(tool.Ifc.get_path()).name == "library-review.ifc"
fixture = Path(tool.Ifc.get_path())
before_hash = library.digest(fixture)
before_ids = {p.GlobalId for p in tool.Ifc.get().by_type("IfcElement")}

class DragHarness:
    product_id = "street-h"
    finish = REVIEWLIB_OT_DragIFC.finish
    _modal = REVIEWLIB_OT_DragIFC._modal
    def __init__(self):
        self._point = None
        self._target = None
        records = json.loads((library._catalog_path.parent / "native-assets/placements.json").read_text())
        self._preview = PlacementPreview(records[self.product_id]["asset_bounds_m"])
    def report(self, kind, message):
        raise AssertionError(message)

viewport = next(a for a in bpy.context.screen.areas if a.type == "VIEW_3D")
region = next(r for r in viewport.regions if r.type == "WINDOW")
with bpy.context.temp_override(area=viewport, region=region):
    preview_test = DragHarness()
    preview_test._preview.update((viewport, region), (2., 3., 0.))
    vertices = preview_test._preview.vertices()
    assert len(vertices) == 24
    assert abs(min(v.z for v in vertices)) < 1e-6
    preview_test._preview.update((viewport, region), (4., 5., 1.))
    assert all((b-a-Vector((2, 2, 1))).length < 1e-5 for a, b in zip(vertices, preview_test._preview.vertices()))
    preview_test._preview.draw()  # Exercise actual GPU shader and batch.
    # Cancellation/outside-drop must never begin an IFC transaction.
    assert REVIEWLIB_OT_DragIFC.modal(preview_test, bpy.context,
        SimpleNamespace(type="ESC", value="PRESS")) == {"CANCELLED"}
    assert preview_test._preview.handle is None
    assert REVIEWLIB_OT_DragIFC.modal(DragHarness(), bpy.context,
        SimpleNamespace(type="LEFTMOUSE", value="RELEASE", mouse_x=-10, mouse_y=-10)) == {"CANCELLED"}
    assert before_ids == {p.GlobalId for p in tool.Ifc.get().by_type("IfcElement")}
    bpy.ops.ed.undo_push(message="Before synthetic IFC drag test")
    release_test = DragHarness()
    result = REVIEWLIB_OT_DragIFC.modal(release_test, bpy.context,
        SimpleNamespace(type="LEFTMOUSE", value="RELEASE", mouse_x=region.x+region.width*.25,
            mouse_y=region.y+region.height*.75))
    assert result == {"FINISHED"}
    assert release_test._preview.handle is None
    created = [p for p in tool.Ifc.get().by_type("IfcElement") if p.GlobalId not in before_ids]
    assert len(created) == 1
    element = created[0]
    assert len([o for o in bpy.context.scene.objects if tool.Ifc.get_entity(o) == element]) == 1
    bpy.ops.ed.undo_push(message="Synthetic IFC drag test")
    assert bpy.ops.ed.undo() == {"FINISHED"}
assert before_ids == {p.GlobalId for p in tool.Ifc.get().by_type("IfcElement")}
assert library.digest(fixture) == before_hash
out = Path(__file__).resolve().parents[2] / "output/review/approved-product-library/ifc-card-drag-validation.json"
out.write_text(json.dumps({"status": "pass", "cancel_no_mutation": True, "outside_drop_no_mutation": True,
    "bounding_box_gpu_draw": True, "bounding_box_translation_matches_pointer": True,
    "preview_handler_removed_on_cancel_and_release": True,
    "synthetic_release_creates_one_real_ifc_element": True, "single_blender_object": True,
    "native_undo_restores_previous_elements": True, "disk_unchanged": True,
    "physical_pointer_test": "not_verified; synthetic handler events only"}, indent=2)+"\n")
print("NATIVE_DRAG_HANDLER_PASS")
