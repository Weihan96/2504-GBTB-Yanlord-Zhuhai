"""Package the drag gesture fix without touching models, approvals or index."""
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import zipfile
import ast

root = Path(__file__).resolve().parents[2]
out = root / 'output/review/approved-product-library'
runtime = out / 'runtime'
addon = root / 'pipeline/addons/highpoly_review_library'
read = lambda p: json.loads(p.read_text())
digest = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
plan = read(out/'portable-cleanup-plan.json')
old = read(out/'portable-runtime-manifest.json')
report = read(out/'drag-gesture-validation.json')
assert report['status'] == 'pass'
info = next(n for n in ast.parse((addon/'__init__.py').read_text()).body
            if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='bl_info' for t in n.targets))
version = ast.literal_eval(info.value)['version']
if version >= (0,8,3):
    bounds = read(out/'loading-bounds-validation.json')
    assert bounds['status']=='pass' and bounds['before']==bounds['after']
    assert bounds['blue_outline_visible'] and bounds['handlers_removed']
version = '.'.join(map(str,version))
index = subprocess.check_output(['git','ls-files','--stage','-z'],cwd=root)
assert hashlib.sha256(index).hexdigest() == plan['protected_index_sha256']
assert digest(root/'2504 GBTB Yanlord Zhuhai.ifc') == plan['formal_ifc_sha256']
for relative, sha in plan['protected_files'].items():
    assert digest(root/relative) == sha, relative
for item in old['files']:
    if item['path'].startswith(('ifc/','approvals/','previews/','native-assets/')) or item['path']=='catalog.json':
        assert digest(runtime/item['path']) == item['sha256'], item['path']
for source in addon.glob('*.py'):
    shutil.copy2(source,runtime/'addon/highpoly_review_library'/source.name)
with zipfile.ZipFile(out/'highpoly_review_library.zip','w',zipfile.ZIP_DEFLATED) as archive:
    for source in sorted(addon.glob('*.py')):
        assert digest(source)==digest(runtime/'addon/highpoly_review_library'/source.name)
        archive.write(source,'highpoly_review_library/'+source.name)
files = [{'path':str(p.relative_to(runtime)),'bytes':p.stat().st_size,'sha256':digest(p)}
         for p in sorted(runtime.rglob('*')) if p.is_file() and '__pycache__' not in p.parts]
manifest = dict(status='pass',version=version,files=files,count=len(files),
                bytes=sum(p['bytes'] for p in files),review_categories=old['review_categories'])
(out/'portable-runtime-manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2)+'\n')
assert subprocess.check_output(['git','ls-files','--stage','-z'],cwd=root)==index
print(json.dumps({'status':'pass','version':version,'runtime_bytes':manifest['bytes'],
                  'protected_sources':len(plan['protected_files']),'index_unchanged':True}))
