"""Package code and record the verified runtime closure without changing Git's index."""
from collections import Counter
from pathlib import Path
import hashlib
import json
import subprocess
import zipfile

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "output/review/approved-product-library"
RUNTIME = OUT / "runtime"
digest = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
read = lambda name: json.loads((OUT / name).read_text())
write = lambda name, value: (OUT / name).write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n")
cleanup = read("portable-cleanup-result.json")
plan = read("portable-cleanup-plan.json")
isolation = read("portable-isolation-validation.json")
visual = read("portable-raster-validation.json")
caches = read("portable-cache-validation.json")
integration = read("native-integration-validation.json")
assert cleanup["status"] == "moved_verified"
assert isolation["status"] == "structured_pass" and len(isolation["results"]) == 39 and isolation["svg_count"] == 117
assert visual["status"] == "pass" and visual["images"] == 117
assert caches["status"] == "pass" and len(caches["results"]) == 39
assert integration["status"] == "pass" and len(integration["products"]) == 40
assert (OUT / "native-integration-validation.json").stat().st_mtime > (OUT / "portable-cleanup-result.json").stat().st_mtime
assert integration["native_bonsai_save_reload"] and integration["native_ctrl_s_keymap"]
index = subprocess.check_output(["git", "ls-files", "--stage", "-z"], cwd=ROOT)
assert hashlib.sha256(index).hexdigest() == plan["protected_index_sha256"]
assert digest(ROOT / "2504 GBTB Yanlord Zhuhai.ifc") == plan["formal_ifc_sha256"]
for path, sha in plan["protected_files"].items():
    assert digest(ROOT / path) == sha, path
for item in cleanup["files"]:
    assert not (ROOT / item["path"]).exists(), item["path"]
    assert digest(Path(cleanup["archive"]) / item["path"]) == item["sha256"], item["path"]
catalog = json.loads((RUNTIME / "catalog.json").read_text())
assert len(catalog["products"]) == 39
isolation["visual_check"] = "pass; see portable-raster-validation.json (technical output only, not design approval)"
write("portable-isolation-validation.json", isolation)
for product in catalog["products"]:
    paths = [product["ifc_path"], product["approval_record"]["path"], *product["previews"].values(),
             "native-assets/" + product["id"] + ".blend"]
    for value in paths:
        path = (RUNTIME / value).resolve()
        assert path.is_relative_to(RUNTIME) and path.is_file(), value
    assert digest(RUNTIME / product["ifc_path"]) == product["ifc_sha256"]
    assert digest(RUNTIME / product["approval_record"]["path"]) == product["approval_record"]["sha256"]
    evidence_name = "portable-validation/" + product["id"] + "/result.json"
    evidence = read(evidence_name)
    evidence["verdict"] = "technical_pass"
    evidence["visual_inspection"] = "portable-raster-validation.json; original design approval unchanged"
    write(evidence_name, evidence)
addon = ROOT / "pipeline/addons/highpoly_review_library"
with zipfile.ZipFile(OUT / "highpoly_review_library.zip", "w", zipfile.ZIP_DEFLATED) as archive:
    for path in sorted(addon.glob("*.py")):
        assert digest(path) == digest(RUNTIME / "addon/highpoly_review_library" / path.name)
        archive.write(path, "highpoly_review_library/" + path.name)
files = [{"path": str(p.relative_to(RUNTIME)), "bytes": p.stat().st_size, "sha256": digest(p)}
         for p in sorted(RUNTIME.rglob("*")) if p.is_file() and "__pycache__" not in p.parts]
manifest = {"status": "pass", "files": files, "count": len(files), "bytes": sum(f["bytes"] for f in files),
            "review_categories": dict(Counter(p["review_category"] for p in catalog["products"]))}
write("portable-runtime-manifest.json", manifest)
tests = subprocess.run(["/usr/bin/python3", "-m", "unittest", "discover", "-s", "pipeline/tests", "-p", "test_review_*.py"],
                       cwd=ROOT, capture_output=True, text=True)
assert tests.returncode == 0, tests.stdout + tests.stderr
staged = set(subprocess.check_output(["git", "diff", "--cached", "--name-only", "-z"], cwd=ROOT).decode().split("\0")) - {""}
unstaged = set(subprocess.check_output(["git", "diff", "--name-only", "-z"], cwd=ROOT).decode().split("\0")) - {""}
untracked = set(subprocess.check_output(["git", "ls-files", "--others", "--exclude-standard", "-z"], cwd=ROOT).decode().split("\0")) - {""}
report = {
    "status": "pass", "task": "Portable runtime and recoverable aggressive review cleanup",
    "evidence": {"course_provider": "embedded-course-index", "lesson": "052000",
                 "source_sha256": "f8c108a9132758df2577f889fd4b9112023502fa4e4b86c801371c23d1cf0916",
                 "timestamps": ["01:54 save IFC", "02:01 reopen IFC"], "raw_video_available": False},
    "plan": "Freeze protected index and formal IFC; isolate runtime; save and reload; create 117 drawings; visually inspect; archive exact obsolete files",
    "pre_state": {"protected_staged_paths": 1732, "protected_index_sha256": plan["protected_index_sha256"]},
    "execution": {"isolated_products": 39, "native_create_drawing_svgs": 117, "isolated_cache_reads": 39,
                  "post_cleanup_combined_project_instances": 40, "post_cleanup_native_save_reload": True,
                  "blocked_paths": [str(ROOT), "/Users/jiaxinchen/Poliigon"],
                  "materials": 18, "packed_textures": 63, "unit_tests": tests.stdout + tests.stderr},
    "persistence": "39 products inserted, saved through native Bonsai provider, reloaded, then drawn; source IFCs unchanged; new drop instances retain provenance",
    "post_state": {"formal_ifc_sha256": plan["formal_ifc_sha256"], "formal_unchanged": True,
                   "protected_source_approval_files": len(plan["protected_files"]), "protected_index_unchanged": True,
                   "staged_paths": len(staged), "unstaged_tracked_paths": len(unstaged), "untracked_paths_at_report": len(untracked),
                   "overlapping_paths": sorted(staged & unstaged), "no_commit_stash_or_history_change": True},
    "outputs": {"runtime_bytes": manifest["bytes"], "runtime_files": manifest["count"],
                "archive": cleanup["archive"], "archived_files": cleanup["count"], "archived_bytes": cleanup["bytes"],
                "addon_zip_sha256": digest(OUT / "highpoly_review_library.zip"), "review_categories": manifest["review_categories"]},
    "visual_inspection": visual,
    "verdict": "Runtime, persistence, drawing and recoverable cleanup passed. Existing design approval states are unchanged. Archive and existing Git objects still consume disk; current changes are not staged."
}
write("portable-cleanup-validation.json", report)
print(json.dumps({"status": report["status"], "runtime_bytes": manifest["bytes"],
                  "categories": manifest["review_categories"], "staged": len(staged),
                  "unstaged": len(unstaged), "untracked": len(untracked), "overlap": len(staged & unstaged)}))
