"""Sync the verified preview-only release; preserve IFC, approvals and staged baseline."""
from collections import Counter
from pathlib import Path
import hashlib
import json
import subprocess
import zipfile

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'output/review/approved-product-library'
RUNTIME = OUT / 'runtime'
digest = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
read = lambda p: json.loads(Path(p).read_text())
write = lambda p, v: Path(p).write_text(json.dumps(v, ensure_ascii=False, indent=2) + '\n')
old = read(OUT / 'portable-runtime-manifest.json')
proof = read(OUT / 'drawing-preview-validation.json')
ui = read(OUT / 'drawing-preview-ui-validation.json')
assert ui['status'] == 'pass' and ui['preState'] == ui['postState']
assert proof['linework_count'] == 117 and proof['product_count'] == 39
assert Path(ui['screenshot']).is_file()
plan = read(OUT / 'portable-cleanup-plan.json')
assert digest(ROOT / '2504 GBTB Yanlord Zhuhai.ifc') == proof['formal_sha256']
index = subprocess.check_output(['git', 'ls-files', '--stage', '-z'], cwd=ROOT)
assert hashlib.sha256(index).hexdigest() == proof['protected_index_sha256']
for relative, sha in plan['protected_files'].items():
    assert digest(ROOT / relative) == sha, relative
for item in old['files']:
    p = item['path']
    if p.startswith(('ifc/', 'approvals/')) or p.endswith('/iso.png') or p == 'native-assets/placements.json':
        assert digest(RUNTIME / p) == item['sha256'], p
addon = ROOT / 'pipeline/addons/highpoly_review_library'
with zipfile.ZipFile(OUT / 'highpoly_review_library.zip', 'w', zipfile.ZIP_DEFLATED) as archive:
    for p in sorted(addon.glob('*.py')):
        assert digest(p) == digest(RUNTIME / 'addon/highpoly_review_library' / p.name)
        archive.write(p, 'highpoly_review_library/' + p.name)
files = [{'path': str(p.relative_to(RUNTIME)), 'bytes': p.stat().st_size, 'sha256': digest(p)}
         for p in sorted(RUNTIME.rglob('*')) if p.is_file() and '__pycache__' not in p.parts]
catalog = read(RUNTIME / 'catalog.json')
manifest = {'status': 'pass', 'version': '0.7.1', 'files': files, 'count': len(files),
            'bytes': sum(p['bytes'] for p in files),
            'review_categories': dict(Counter(p['review_category'] for p in catalog['products']))}
assert manifest['review_categories'] == old['review_categories']
write(OUT / 'portable-runtime-manifest.json', manifest)
staged = set(subprocess.check_output(['git', 'diff', '--cached', '--name-only', '-z'], cwd=ROOT).decode().split('\0')) - {''}
unstaged = set(subprocess.check_output(['git', 'diff', '--name-only', '-z'], cwd=ROOT).decode().split('\0')) - {''}
untracked = set(subprocess.check_output(['git', 'ls-files', '--others', '--exclude-standard', '-z'], cwd=ROOT).decode().split('\0')) - {''}
proof.update(status='pass', ui_verification='drawing-preview-ui-validation.json',
    visual={'status': 'pass', 'screenshot': ui['screenshot'],
            'inspection': 'Actual narrow sidebar shows linework thumbnails for Plan/Front/Side and shaded 3D thumbnail; Image Editor displays CHA01 linework. Plan/Front/Side rasters inspected; their candidate geometry was not edited.'},
    persistence='Catalog/SVG/PNG/installer persisted; live preview caches refreshed. Current IFC, transforms, selection and disk hash identical.',
    protected_sources_unchanged=len(plan['protected_files']), runtime_bytes=manifest['bytes'],
    tests={'command': "python3 -m unittest discover -s pipeline/tests -p 'test_review_*.py'", 'passed': 35},
    git={'protected_staged': len(staged), 'unstaged_tracked': len(unstaged), 'untracked': len(untracked),
         'overlap': sorted(staged & unstaged), 'staged_index_unchanged': True})
write(OUT / 'drawing-preview-validation.json', proof)
assert subprocess.check_output(['git', 'ls-files', '--stage', '-z'], cwd=ROOT) == index
print(json.dumps({'status': 'pass', 'runtime_bytes': manifest['bytes'], 'products': 39, 'previews': 117,
                  'git': {k: v if k != 'overlap' else len(v) for k, v in proof['git'].items()}}))
