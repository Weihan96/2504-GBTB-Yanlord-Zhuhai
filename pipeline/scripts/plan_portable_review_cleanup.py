"""Resolve exact obsolete review files; this script never removes anything."""
from pathlib import Path
import json
import hashlib
import subprocess

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/"output/review/approved-product-library"
catalog=json.loads((OUT/"runtime/catalog.json").read_text())
ids={p["id"] for p in catalog["products"]}
digest=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
targets={}
def add(value,reason):
    p=Path(value)
    p=p if p.is_absolute() else ROOT/p
    if not p.is_file():return
    assert not p.is_symlink()
    rel=p.relative_to(ROOT)
    assert "libelle" not in str(rel).lower()
    assert rel.parts[:2]==("output","review") or rel.parts[:2]==("pipeline","scripts")
    assert "runtime" not in rel.parts and "official-source" not in rel.parts
    targets[str(rel)]={"path":str(rel),"bytes":p.stat().st_size,"sha256":digest(p),"reason":reason}
audit=json.loads((OUT/"cleanup-audit.json").read_text())
classes={"blend_cache_candidate_requires_dependency_and_session_check","isolated_legacy_candidate_requires_dependency_and_render_check",
         "full_scene_legacy_candidate_blocked","keep_full_scene_pending_product_or_unresolved_scope","keep_active_library_ifc",
         "keep_active_material_or_library_asset"}
for f in audit["files"]:
    p=Path(f["path"])
    within=p.is_relative_to(Path("output/review/approved-product-library")) or (len(p.parts)>3 and p.parts[:3]==("output","review","highpoly-types") and p.parts[3] in ids)
    if within and f["classification"] in classes and p.suffix.lower() in (".ifc",".blend",".blend1"):
        add(p,"Obsolete review binary; current catalog uses the isolated runtime replacement")
# Duplicated runtime inputs created after the older cleanup audit.
for folder in ("body-packages","candidate-packages","native-assets"):
    for p in (OUT/folder).glob("*"):
        if p.suffix.lower() in (".ifc",".blend",".blend1") or p.name=="placements.json":
            add(p,"Copied into the verified runtime; no longer the active library path")
# Retire only obsolete version-bound packaging and hardcoded-session tests.
for name in ("package_library_v05.py","package_library_v06.py","finalize_library_v04.py","verify_loading_feedback_visual.py","verify_live_body_library.py",
             "render_review_library.py","verify_saved_review_library.py","verify_body_view_drawings.py","normalize_library_body_views.py"):
    add(ROOT/"pipeline/scripts"/name,"Historical version/session-specific helper replaced by portable validation")
slugs={slug.replace("-","_") for slug in ids}
for prefix,suffix in (("extract_","_bonsai_review_ifc.py"),("render_","_bonsai_review.py")):
    for p in (ROOT/"pipeline/scripts").glob(prefix+"*"+suffix):
        if p.name[len(prefix):-len(suffix)] in slugs:
            add(p,"Legacy per-product isolated IFC or Blender camera generator, superseded by current Body runtime")
# Never remove CAD, source records, approval records, or SVG/PNG acceptance evidence.
protected={}
for slug in sorted(ids):
    for p in (ROOT/"output/review/highpoly-types"/slug).rglob("*"):
        if p.is_file() and (p.suffix.lower() in (".dwg",".dxf",".pdf",".zip") or "official-source" in p.parts):
            protected[str(p.relative_to(ROOT))]=digest(p)
for p in catalog["products"]:
    record=OUT/"runtime"/p["approval_record"]["path"]
    protected[str(record.relative_to(ROOT))]=digest(record)
original=json.loads((OUT/"all-review-catalog.json").read_text())
for p in original["products"]:
    record=Path(p["approval_record"]["path"])
    if not record.is_absolute():record=(OUT/record).resolve()
    protected[str(record.relative_to(ROOT))]=digest(record)
assert not set(targets)&set(protected)
report={"status":"resolved_pending_all_117_drawing_and_visual_gates","files":list(targets.values()),
        "count":len(targets),"bytes":sum(f["bytes"] for f in targets.values()),"protected_files":protected,
        "formal_ifc_sha256":digest(ROOT/"2504 GBTB Yanlord Zhuhai.ifc"),
        "protected_index_sha256":hashlib.sha256(subprocess.check_output(["git","ls-files","--stage","-z"])).hexdigest(),
        "operation":"recoverable move outside git, only after verification; no Git history or stash operations"}
(OUT/"portable-cleanup-plan.json").write_text(json.dumps(report,ensure_ascii=False,indent=2)+"\n")
print(json.dumps({"count":report["count"],"bytes":report["bytes"],"protected_files":len(protected)}))
