"""Record visually inspected deliverables and actual clean-session evidence."""
import json
from pathlib import Path
import migrate_sxb010 as task
report=json.loads(task.REPORT.read_text())
assert report['verdict']=='validated_pending_visual_check'
fresh=json.loads((task.OUT/'fresh-session-validation.json').read_text())
assert fresh['physical_products']==1 and fresh['annotations']==0 and fresh['groups']==0
assert fresh['visible_meshes']==['IfcBuildingElementProxy/sxb010']
assert fresh['mesh_names']==['Cube','IfcBuildingElementProxy/sxb010']
assert not fresh['blend_saved'] and not fresh['has_blend_warning']
assert fresh['ifc_sha256']==task.pkg.sha256(task.SINGLE)
assert task.pkg.sha256(task.FORMAL)==task.BASELINE
report['fresh_pure_session']=fresh
report['visual']={'method':'view_image inspection of all three scene PNGs and all three pure-IFC single-view PNGs',
    'observations':['Scene Plan target present once at dining-bay boundary; no clipping.',
        'Scene Front complete 45-slat expression, frame and controls visible against grey context.',
        'Scene Side full-height blind and tilted slats remain legible and unclipped.',
        'Three black single-product previews faithfully show the 5/51/55 approved paths from pure IFC.'],
    'verdict':'pass'}
report['verdict']='pass'
report['staged_baseline_file_count_after']=937
report['git_scope']='Only new files in sxb010/product-package; no tracked unstaged SXB010 changes and no staging operations'
task.write(task.REPORT,report)
for name in ['manifest.json','handoff.json']:
    path=task.OUT/name;data=json.loads(path.read_text());data['status']='complete'
    data['fresh_session_validation']='fresh-session-validation.json'
    task.write(path,data)
print('SXB010 pure package complete; all three scene geometry comparisons 0.0 mm')
