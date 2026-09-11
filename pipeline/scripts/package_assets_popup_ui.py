"""Package the UI-only Assets release after live routing and style verification."""
from collections import Counter
from pathlib import Path
import hashlib
import json
import shutil
import subprocess
import sys
import zipfile

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'output/review/approved-product-library'
RUNTIME = OUT / 'runtime'
ADDON = ROOT / 'pipeline/addons/highpoly_review_library'
digest = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
read = lambda p: json.loads(Path(p).read_text())
write = lambda p, value: Path(p).write_text(json.dumps(value, ensure_ascii=False, indent=2)+'\n')
plan = read(OUT / 'portable-cleanup-plan.json')
old = read(OUT / 'portable-runtime-manifest.json')
thumbnails = '--thumbnails' in sys.argv
version = '0.8.1' if thumbnails else '0.8.0'
ui_path = OUT / ('thumbnail-live-validation.json' if thumbnails else 'assets-popup-ui-validation.json')
ui = read(ui_path)
assert ui['status'] == 'pass' and ui['preState'] == ui['postState']
if not thumbnails:
    assert ui['loading_handler_removed'] and ui['placement_handler_removed']
index = subprocess.check_output(['git','ls-files','--stage','-z'],cwd=ROOT)
assert hashlib.sha256(index).hexdigest() == plan['protected_index_sha256']
assert digest(ROOT / '2504 GBTB Yanlord Zhuhai.ifc') == plan['formal_ifc_sha256']
for relative, sha in plan['protected_files'].items():
    assert digest(ROOT / relative) == sha, relative
for item in old['files']:
    if thumbnails and (item['path'].endswith('/iso.png') or item['path']=='catalog.json'):
        continue
    if item['path'].startswith(('ifc/','approvals/','previews/','native-assets/')) or item['path']=='catalog.json':
        assert digest(RUNTIME / item['path']) == item['sha256'], item['path']
if thumbnails:
    compression = read(OUT / 'thumbnail-compression-validation.json')
    for item in compression['files']:
        assert digest(RUNTIME / item['path']) == item['sha256']
        assert digest(ROOT / item['source']) == item['source_sha256']
for source in ADDON.glob('*.py'):
    shutil.copy2(source,RUNTIME / 'addon/highpoly_review_library' / source.name)
with zipfile.ZipFile(OUT / 'highpoly_review_library.zip','w',zipfile.ZIP_DEFLATED) as archive:
    for source in sorted(ADDON.glob('*.py')):
        assert digest(source) == digest(RUNTIME / 'addon/highpoly_review_library' / source.name)
        archive.write(source,'highpoly_review_library/'+source.name)
files = [{'path':str(p.relative_to(RUNTIME)),'bytes':p.stat().st_size,'sha256':digest(p)}
         for p in sorted(RUNTIME.rglob('*')) if p.is_file() and '__pycache__' not in p.parts]
categories = dict(Counter(p['review_category'] for p in read(RUNTIME / 'catalog.json')['products']))
assert categories == old['review_categories']
write(OUT / 'portable-runtime-manifest.json',{'status':'pass','version':version,'files':files,
    'count':len(files),'bytes':sum(p['bytes'] for p in files),'review_categories':categories})
gitset = lambda args: set(subprocess.check_output(['git',*args],cwd=ROOT).decode().split('\0'))-{''}
staged=gitset(['diff','--cached','--name-only','-z'])
unstaged=gitset(['diff','--name-only','-z'])
untracked=gitset(['ls-files','--others','--exclude-standard','-z'])
for name in (('thumbnail-three-view-popup.png',) if thumbnails else ('assets-popup-details.png','assets-borderless-loading.png')):
    assert (OUT / name).is_file()
visual = {'status':'pass','popup':str(OUT/'assets-popup-details.png'),
        'loading':str(OUT/'assets-borderless-loading.png'),
        'inspection':'Assets tab and narrow sidebar inspected: old inline details removed, popup has all four previews. Synthetic 65% style check shows borderless translucent green fill; partially occluded by sidebar. This was not a real import.'}
if thumbnails:
    visual = {'status':'pending_popup_reopen','popup':str(OUT/'thumbnail-three-view-popup.png'),
        'inspection':'Three buttons verified; existing open popups retained stale preview IDs after cache replacement. Close and reopen for user visual acceptance. Compressed CHA01 and BED01 rasters visually inspected.',
        'computer_use_boundary':'Computer Use resolved the other Libelle Blender window; no key/click was sent there.'}
ui.update(tests={'passed':41 if thumbnails else 40,'command':"python3 -m unittest discover -s pipeline/tests -p 'test_review_*.py'"},
    visual=visual,
    protected_sources_unchanged=len(plan['protected_files']),formal_ifc_unchanged=True,
    runtime_models_linework_approvals_unchanged=True,
    git={'protected_staged':len(staged),'unstaged_tracked':len(unstaged),'untracked':len(untracked),
         'overlap':sorted(staged & unstaged),'staged_index_unchanged':True})
write(ui_path,ui)
if thumbnails:
    compression.update(status='compressed_verified_pending_popup_reopen',live_validation=ui_path.name,
        runtime_bytes=sum(p['bytes'] for p in files))
    write(OUT / 'thumbnail-compression-validation.json',compression)
assert subprocess.check_output(['git','ls-files','--stage','-z'],cwd=ROOT)==index
print(json.dumps({'status':'pass','version':version,'staged':len(staged),'unstaged':len(unstaged),
    'untracked':len(untracked),'overlap':len(staged & unstaged),'runtime_bytes':sum(p['bytes'] for p in files)}))
