"""Data/function-only update: preserve native card IDs, IFC state and registered RNA."""
import ast
import hashlib
import json
import os
from pathlib import Path
import bpy
import bonsai.tool as tool
import highpoly_review_library as lib
from highpoly_review_library import asset_ui, asset_cards, loading_feedback

assert os.getpid() == 33367 and loading_feedback._active is None
ROOT = Path(lib.__file__).resolve().parents[3]
OUT = ROOT / 'output/review/approved-product-library'
assert Path(tool.Ifc.get_path()).resolve() == (OUT / 'ifc-card-ui/library-review.ifc').resolve()


def state():
    return {'pid':os.getpid(),'memory':hashlib.sha256(tool.Ifc.get().to_string().encode()).hexdigest(),
        'disk':hashlib.sha256(Path(tool.Ifc.get_path()).read_bytes()).hexdigest(),
        'objects':[(o.as_pointer(),o.name,[list(r) for r in o.matrix_world]) for o in bpy.data.objects],
        'selected':[o.as_pointer() for o in bpy.context.selected_objects],
        'product':bpy.context.scene.review_library.product,
        'cards':[(c.as_pointer(),c.name) for c in bpy.data.collections if c.get(asset_cards.KEY)==asset_cards.OWNER]}


before=state()
report={'status':'updating','preState':before,'version':'0.8.1'}
report_path=OUT/'thumbnail-live-validation.json'
report_path.write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
tree=ast.parse(Path(asset_ui.__file__).read_text())
fn=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='draw_details')
exec(compile(ast.Module(body=[fn],type_ignores=[]),asset_ui.__file__,'exec'),asset_ui.__dict__)
lib.load_catalog(OUT/'runtime/catalog.json')
lib.bl_info['version']=(0,8,1)
count=0
for card in bpy.data.collections:
    product=asset_cards.product_for_id(card,lib)
    if not product: continue
    with bpy.context.temp_override(id=card):
        assert bpy.ops.ed.lib_id_load_custom_preview(filepath=str(lib.resolve_path(product['previews']['iso'])))=={'FINISHED'}
    assert card.preview and all(card.preview.image_size)
    count+=1
assert count==39
for area in bpy.context.screen.areas: area.tag_redraw()
bpy.ops.wm.redraw_timer(type='DRAW_WIN_SWAP',iterations=2)
after=state()
assert before==after
report.update(status='pass',postState=after,thumbnail_cards_refreshed=count,
    persistence='Only thumbnail buffers, catalog and popup draw function updated; no IFC write/reload or RNA replacement',
    courseEvidence={'mode':'embedded-course-index','lesson':'052000','fact':'UI/Blender data is separate from IFC persistence'})
report_path.write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps({'status':'pass','cards':count,'user_state_unchanged':before==after}))
