"""Verify registered UI routing and prepare a timed style-only loading preview."""
import hashlib
import json
import os
from pathlib import Path
import bpy
import bonsai.tool as tool
import highpoly_review_library as lib
from highpoly_review_library import native_assets as native, asset_ui, loading_feedback as feedback
from mathutils import Vector

assert os.getpid() == 33367 and feedback._active is None
ROOT = Path(lib.__file__).resolve().parents[3]
OUT = ROOT / 'output/review/approved-product-library'
assert Path(tool.Ifc.get_path()).resolve() == (OUT / 'ifc-card-ui/library-review.ifc').resolve()


def state():
    return {'memory': hashlib.sha256(tool.Ifc.get().to_string().encode()).hexdigest(),
            'disk': hashlib.sha256(Path(tool.Ifc.get_path()).read_bytes()).hexdigest(),
            'objects': [(o.as_pointer(), [list(r) for r in o.matrix_world]) for o in bpy.data.objects],
            'selection': [o.as_pointer() for o in bpy.context.selected_objects],
            'product': bpy.context.scene.review_library.product,
            'elements': len(tool.Ifc.get().by_type('IfcElement'))}


before = state()
calls = []
original_product = native.asset_product
original_details = asset_ui.open_details
area = next(a for a in bpy.context.screen.areas if a.type == 'VIEW_3D')
region = next(r for r in area.regions if r.type == 'WINDOW')
try:
    # Exercise the actual registered activation operator, with a controlled
    # native-card context; no mouse events, drag or insert are simulated here.
    native.asset_product = lambda context: lib.entry(before['product'])
    asset_ui.open_details = lambda context: calls.append(context.scene.review_library.product)
    with bpy.context.temp_override(area=area, region=region):
        assert bpy.ops.review_library.select_asset() == {'FINISHED'}
    assert calls == [before['product']]
finally:
    native.asset_product = original_product
    asset_ui.open_details = original_details
assert state() == before
report = {'status': 'routing_pass_pending_style_capture', 'version': '0.8.0',
    'task': 'Assets click popup and borderless loading feedback',
    'preState': before, 'postState': state(), 'registered_click_routes_to_popup': True,
    'physical_mouse_test': 'Not automated; registered operator tested with controlled card context',
    'migration_note': 'First category activation needed a sidebar redraw; subsequently activated Assets. No IFC API mutation, save or reload was used.',
    'courseEvidence': {'mode':'embedded-course-index', 'lesson':'052000', 'fact':'UI/Blender state and persisted IFC state are separate'},
    'visual': {'popup':str(OUT / 'assets-popup-details.png')},
    'persistence': 'UI code persisted separately; current IFC remains untouched'}
(OUT / 'assets-popup-ui-validation.json').write_text(json.dumps(report, ensure_ascii=False, indent=2)+'\n')
# Real draw code with a clearly labelled synthetic fraction: purely a style
# check, not a claimed import-progress measurement. Restore automatically.
with bpy.context.temp_override(area=area, region=region):
    d = area.spaces.active.region_3d.view_distance * .22
    preview = native.PlacementPreview([[-d/2,-d/2,0],[d/2,d/2,d*2]])
    preview.update((area,region), area.spaces.active.region_3d.view_location - Vector((0,0,d*.65)))
    loading = feedback.LoadingFeedback(preview, 'style-only-test', bpy.context)
    loading.__enter__()
    loading.advance(.65, '样式预览 65%（不导入模型）')
    assert preview.outline_visible is False


def cleanup():
    with bpy.context.temp_override(area=area, region=region):
        loading.__exit__(None, None, None)
        preview.close()
    assert state() == before
    report.update(status='pass', postState=state(), borderless_style_only=True,
        loading_handler_removed=loading.handle is None, placement_handler_removed=preview.handle is None)
    (OUT / 'assets-popup-ui-validation.json').write_text(json.dumps(report, ensure_ascii=False, indent=2)+'\n')
    asset_ui._style_cleanup = None


asset_ui._style_cleanup = cleanup
bpy.app.timers.register(cleanup, first_interval=45)
print('Registered activation routes correctly; style-only preview active for up to 45 seconds')
