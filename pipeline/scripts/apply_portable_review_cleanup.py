"""Move precisely verified legacy review files to a recoverable outside-git archive."""
from pathlib import Path
from datetime import datetime
import json
import hashlib
import shutil
import subprocess

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/"output/review/approved-product-library"
digest=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
plan=json.loads((OUT/"portable-cleanup-plan.json").read_text())
validation=json.loads((OUT/"portable-isolation-validation.json").read_text())
visual=json.loads((OUT/"portable-raster-validation.json").read_text())
assert validation["status"]=="structured_pass" and len(validation["results"])==39 and validation["svg_count"]==117
assert visual["status"]=="pass" and visual["images"]==117
assert digest(ROOT/"2504 GBTB Yanlord Zhuhai.ifc")==plan["formal_ifc_sha256"]=="7a50b87e8f48a7c2bfbaa9f1dabdfca7155fa684b2325c66aa1ae15c8ab0a25c"
index=subprocess.check_output(["git","ls-files","--stage","-z"])
assert hashlib.sha256(index).hexdigest()==plan["protected_index_sha256"]
capsule=Path(validation["capsule"])/"runtime"
for p in capsule.rglob("*"):
    if p.is_file() and "__pycache__" not in p.parts:
        target=OUT/"runtime"/p.relative_to(capsule)
        assert target.is_file() and digest(target)==digest(p), str(target)
for rel,sha in plan["protected_files"].items():
    assert digest(ROOT/rel)==sha,rel
# Live users of old review files block cleanup, not merely another task's Blender.
processes=subprocess.check_output(["ps","-axo","pid,command"]).decode()
for row in processes.splitlines():
    if "/Applications/Blender.app/Contents/MacOS/Blender" in row and str(ROOT) in row:
        raise RuntimeError("A Blender process still references this worktree; inspect before cleanup: "+row)
for item in plan["files"]:
    source=ROOT/item["path"]
    assert source.is_file() and not source.is_symlink() and digest(source)==item["sha256"],item["path"]
archive=Path("/Users/jiaxinchen/.codex/archives")/("highpoly-library-cleanup-"+datetime.now().strftime("%Y%m%d-%H%M%S"))
assert not archive.exists()
archive.mkdir(parents=True)
report={"status":"moving","archive":str(archive),"recoverable":True,"files":[],"protected_index_sha256":plan["protected_index_sha256"]}
report_path=OUT/"portable-cleanup-result.json"
def save():
    data=json.dumps(report,ensure_ascii=False,indent=2)+"\n"
    report_path.write_text(data)
    (archive/"recovery-manifest.json").write_text(data)
save()
try:
    for item in plan["files"]:
        source=ROOT/item["path"]
        target=archive/item["path"]
        target.parent.mkdir(parents=True,exist_ok=True)
        assert not target.exists()
        shutil.move(str(source),str(target))
        assert digest(target)==item["sha256"] and not source.exists()
        report["files"].append(item)
        save()
    assert subprocess.check_output(["git","ls-files","--stage","-z"])==index
    assert digest(ROOT/"2504 GBTB Yanlord Zhuhai.ifc")==plan["formal_ifc_sha256"]
    for rel,sha in plan["protected_files"].items():
        assert digest(ROOT/rel)==sha
    report.update(status="moved_verified",count=len(report["files"]),bytes=sum(f["bytes"] for f in report["files"]),formal_ifc_unchanged=True,source_and_approval_records_unchanged=True)
finally:
    save()
print(json.dumps({k:v for k,v in report.items() if k!="files"}),flush=True)
