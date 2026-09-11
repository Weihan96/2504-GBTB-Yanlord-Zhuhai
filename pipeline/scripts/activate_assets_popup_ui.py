"""Narrow UI-only update of the owned review session; no IFC/RNA property reload."""
import ast
import hashlib
import importlib
import json
import os
from pathlib import Path
import bpy
import bonsai.tool as tool
import highpoly_review_library as lib
from highpoly_review_library import native_assets as native, loading_feedback as feedback

ROOT = Path(lib.__file__).resolve().parents[3]
OUT = ROOT / 'output/review/approved-product-library'
assert os.getpid() == 33367
assert Path(tool.Ifc.get_path()).resolve() == (OUT / 'ifc-card-ui/library-review.ifc').resolve()
assert feedback._active is None, 'Wait for insertion to finish'
assert not any(op.bl_idname == 'REVIEW_LIBRARY_OT_drag_ifc' for op in bpy.context.window.modal_operators)


def state():
    f = tool.Ifc.get()
    return {'pid': os.getpid(), 'file': tool.Ifc.get_path(), 'identity': id(f),
        'memory_hash': hashlib.sha256(f.to_string().encode()).hexdigest(),
        'disk_hash': hashlib.sha256(Path(tool.Ifc.get_path()).read_bytes()).hexdigest(),
        'objects': [(o.as_pointer(), o.name, [list(r) for r in o.matrix_world]) for o in bpy.data.objects],
        'selected': [o.as_pointer() for o in bpy.context.selected_objects],
        'active': bpy.context.view_layer.objects.active.as_pointer() if bpy.context.view_layer.objects.active else None,
        'product': bpy.context.scene.review_library.product,
        'property_pointer': bpy.context.scene.review_library.as_pointer(),
        'catalog_pointer': id(lib._catalog), 'preview_pointer': id(lib._preview),
        'cards': [(c.as_pointer(), c.name) for c in bpy.data.collections if c.get('review_library_card')]}


before = state()
(OUT / 'assets-popup-ui-prestate.json').write_text(json.dumps(before, ensure_ascii=False, indent=2) + '\n')
old_panel = getattr(lib, 'REVIEWLIB_PT_Library', None)
assert old_panel is not None, 'UI already updated; do not replay a live migration'
spaces = [(a, a.spaces.active.show_region_ui) for a in bpy.context.screen.areas if a.type == 'VIEW_3D']
try:
    # Detach the asset-view widget before replacing its Panel, retaining all
    # Collection cards, properties, native operators and undo history.
    for area, shown in spaces:
        area.spaces.active.show_region_ui = False
    bpy.ops.wm.redraw_timer(type='DRAW_WIN_SWAP', iterations=2)
    bpy.utils.unregister_class(old_panel)
    lib.CLASSES = tuple(c for c in lib.CLASSES if c is not old_panel)
    del lib.REVIEWLIB_PT_Library
    from highpoly_review_library import asset_ui
    asset_ui.register()
    importlib.reload(feedback)
    native.LoadingFeedback = feedback.LoadingFeedback
    tree = ast.parse(Path(native.__file__).read_text())
    preview = next(n for n in tree.body if isinstance(n, ast.ClassDef) and n.name == 'PlacementPreview')
    exec(compile(ast.Module(body=[preview], type_ignores=[]), native.__file__, 'exec'), native.__dict__)
    select = next(n for n in tree.body if isinstance(n, ast.ClassDef) and n.name == 'REVIEWLIB_OT_SelectAsset')
    execute = next(n for n in select.body if isinstance(n, ast.FunctionDef) and n.name == 'execute')
    exec(compile(ast.Module(body=[execute], type_ignores=[]), native.__file__, 'exec'), native.__dict__)
    native.REVIEWLIB_OT_SelectAsset.execute = native.__dict__.pop('execute')
    tree = ast.parse(Path(lib.__file__).read_text())
    nodes = [n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name in ('register', 'unregister')]
    exec(compile(ast.Module(body=nodes, type_ignores=[]), lib.__file__, 'exec'), lib.__dict__)
    lib.bl_info.update(name='IFC Assets', version=(0,8,0), location='3D View > Sidebar > Assets')
finally:
    for area, shown in spaces:
        area.spaces.active.show_region_ui = shown
        area.tag_redraw()
    # Panel categories are rebuilt only after the restored sidebar is drawn.
    bpy.ops.wm.redraw_timer(type='DRAW_WIN_SWAP', iterations=2)
    for area, shown in spaces:
        for region in area.regions:
            if region.type == 'UI' and hasattr(region, 'active_panel_category'):
                region.active_panel_category = 'Assets'
        area.tag_redraw()
bpy.ops.wm.redraw_timer(type='DRAW_WIN_SWAP', iterations=2)
after = state()
assert before == after, 'User state changed during interface update'
report = {'status': 'ui_updated_pending_visual', 'preState': before, 'postState': after,
    'version': '0.8.0', 'preserved_ifc_and_undo': True, 'old_panel_removed': True,
    'new_tab': 'Assets', 'images': [], 'save_boundary': 'UI only; no IFC save, reload or insertion'}
(OUT / 'assets-popup-ui-validation.json').write_text(json.dumps(report, ensure_ascii=False, indent=2) + '\n')
bpy.ops.screen.screenshot(filepath=str(OUT / 'assets-popup-panel.png'))
print(json.dumps({'status': report['status'], 'model_unchanged': before == after}))
