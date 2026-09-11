"""Refresh only preview data in the owned acceptance window; never reload IFC/RNA."""
import hashlib
import importlib
import json
import os
from pathlib import Path
import bpy
import bonsai.tool as tool
import highpoly_review_library as library
from highpoly_review_library import catalog_rules

ROOT = Path(library.__file__).resolve().parents[3]
OUT = ROOT / 'output/review/approved-product-library'
RUNTIME = OUT / 'runtime'
EXPECTED_IFC = OUT / 'ifc-card-ui/library-review.ifc'
assert os.getpid() == 33367, 'Wrong task window; inspect ownership before running'
assert Path(tool.Ifc.get_path()).resolve() == EXPECTED_IFC.resolve()
assert not getattr(library.native_assets, '_active_drag', None), 'Do not refresh during a drag'
digest = lambda path: hashlib.sha256(Path(path).read_bytes()).hexdigest()


def state():
    model = tool.Ifc.get()
    return {'pid': os.getpid(), 'ifc': tool.Ifc.get_path(), 'disk_sha256': digest(EXPECTED_IFC),
            'memory_ifc_sha256': hashlib.sha256(model.to_string().encode()).hexdigest(),
            'ifc_identity': id(model), 'elements': len(model.by_type('IfcElement')),
            'objects': [(o.as_pointer(), o.name, [list(r) for r in o.matrix_world]) for o in bpy.data.objects],
            'selected': [o.as_pointer() for o in bpy.context.selected_objects],
            'active': bpy.context.view_layer.objects.active.as_pointer() if bpy.context.view_layer.objects.active else None,
            'scene': bpy.context.scene.as_pointer(), 'product': bpy.context.scene.review_library.product}


before = state()
importlib.reload(catalog_rules)
library.load_catalog(RUNTIME / 'catalog.json')
library.bl_info['version'] = (0, 7, 1)
# Existing Image Editor datablocks cache file contents independently of thumbnail caches.
paths = {library.resolve_path(p['previews'][v]) for p in library._catalog['products'] for v in ('plan', 'front', 'side')}
for image in bpy.data.images:
    if image.filepath and Path(bpy.path.abspath(image.filepath)).resolve() in paths:
        image.reload()
screens = [a for a in bpy.context.screen.areas if a.type == 'IMAGE_EDITOR']
assert screens, 'Use the existing preview editor; do not change the user layout'
product = library.entry(before['product'])
buttons = []
for view in ('plan', 'front', 'side', 'iso', 'plan'):
    assert bpy.ops.review_library.image(view=view) == {'FINISHED'}
    actual = screens[0].spaces.active.image
    assert Path(bpy.path.abspath(actual.filepath)).resolve() == library.resolve_path(product['previews'][view])
    assert tuple(actual.size) == ((1600, 1200) if view != 'iso' else tuple(actual.size))
    buttons.append({'product': product['id'], 'view': view, 'image': actual.filepath,
                    'size': list(actual.size), 'sha256': digest(actual.filepath)})
for area in bpy.context.screen.areas:
    area.tag_redraw()
bpy.ops.wm.redraw_timer(type='DRAW_WIN_SWAP', iterations=2)
after = state()
assert before == after, 'Preview refresh changed model, transform, selection or persistence state'
report = {'status': 'pass', 'task': 'Repair drawing preview entry points without changing IFC state',
          'preState': before, 'postState': after, 'buttons': buttons,
          'persistence': 'Only preview data refreshed; no IFC save/reload and no RNA replacement',
          'preview_product_count': len(library._catalog['products']),
          'preview_cache_count': len(library._preview), 'screenshot': str(OUT / 'drawing-preview-ui.png')}
(OUT / 'drawing-preview-ui-validation.json').write_text(json.dumps(report, ensure_ascii=False, indent=2) + '\n')
bpy.ops.screen.screenshot(filepath=report['screenshot'])
print(json.dumps({'status': 'pass', 'products': len(library._catalog['products']),
                  'tested_buttons': len(buttons), 'model_unchanged': before == after}, ensure_ascii=False))
