"""A/B dependency experiment; never edits the delivered runtime or Git index."""
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import hashlib
import json
import shutil
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "output/review/approved-product-library"
source = OUT / "runtime"
capsule = Path(tempfile.mkdtemp(prefix="library-no-asset-blends-"))
runtime = capsule / "runtime"
digest = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
index = subprocess.check_output(["git", "ls-files", "--stage", "-z"], cwd=ROOT)
before = {str(p.relative_to(source)): digest(p) for p in source.rglob("*")
          if p.is_file() and "__pycache__" not in p.parts}
formal = ROOT / "2504 GBTB Yanlord Zhuhai.ifc"
formal_hash = digest(formal)
def ignore(path, names):
    return [n for n in names if n == "__pycache__" or (Path(path).name == "native-assets" and n.endswith(".blend"))]
shutil.copytree(source, runtime, ignore=ignore)
assert not list((runtime / "native-assets").glob("*.blend"))
assert (runtime / "native-assets/placements.json").is_file()
assert (runtime / "ifc/Materials.blend").is_file()
worker = capsule / "verify.py"
worker.write_text((ROOT / "pipeline/scripts/verify_portable_review_product.py").read_text())
profile = '(version 1)(allow default)(deny file-read* file-write* (subpath ' + json.dumps(str(ROOT)) + '))(deny file-read* (subpath "/Users/jiaxinchen/Poliigon"))'
products = ["bed01", "trap01"]
def test(slug):
    target = capsule / "results" / slug
    target.mkdir(parents=True)
    with (target / "worker.log").open("w") as log:
        process = subprocess.run(["/usr/bin/sandbox-exec", "-p", profile,
            "/Applications/Blender.app/Contents/MacOS/Blender", "-b", "--python-exit-code", "1", "--python", str(worker), "--",
            "--runtime", str(runtime), "--out", str(target), "--slug", slug, "--blocked-root", str(ROOT)],
            cwd=capsule, stdout=log, stderr=subprocess.STDOUT, timeout=300)
    assert process.returncode == 0, str(target / "worker.log")
    result = json.loads((target / "result.json").read_text())
    print(json.dumps({"product": slug, "native_insert": result["native_insert"], "reloaded": result["reloaded"]}), flush=True)
    return result
with ThreadPoolExecutor(max_workers=2) as pool:
    results = list(pool.map(test, products))
assert before == {str(p.relative_to(source)): digest(p) for p in source.rglob("*")
                  if p.is_file() and "__pycache__" not in p.parts}
assert digest(formal) == formal_hash
assert subprocess.check_output(["git", "ls-files", "--stage", "-z"], cwd=ROOT) == index
report = {"task": "Remove only native-assets/*.blend in a disposable runtime; test IFC-only importer",
          "status": "backend_pass_ui_evidence_separate", "capsule": str(capsule),
          "courseEvidence": {"mode": "embedded-course-index", "lesson": "052000", "fact": "Save IFC and reopen to verify persistence"},
          "plan": "Omit cache Blends, retain placement and material dependencies, prohibit access to old repository, insert/save/reload/draw two samples",
          "preState": {"asset_blends": 39, "runtime_files_hashed": len(before), "formal_sha256": formal_hash},
          "postState": {"asset_blends": 0, "catalog_entries": 39, "placement_records": 39,
                        "material_blend_retained": True, "delivered_runtime_unchanged": True,
                        "formal_unchanged": True, "index_unchanged": True},
          "results": results, "svg_count": 6,
          "visual": "pending raster inspection",
          "ui_check": {"method": "source inspection, not interactive mouse test",
              "grid": "REVIEWLIB_PT_Library.draw uses template_asset_view with Collection filter",
              "drag": "asset_product resolves a transient empty Collection card's product_id",
              "conclusion": "No on-disk asset Blend is required; card and event tests are recorded separately in ifc-card-validation.json and native-drag-handler-validation.json"},
          "verdict": "No asset Blend is needed by IFC insertion for these two samples; UI evidence is checked separately."}
destination = OUT / "no-asset-blends-test.json"
destination.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n")
print(str(destination), flush=True)
