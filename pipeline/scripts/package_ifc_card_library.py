"""Publish the verified geometry-free card library and archive exact asset caches."""
from pathlib import Path
from collections import Counter
import hashlib
import json
import shutil
import subprocess
import tempfile
import zipfile

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "output/review/approved-product-library"
RUNTIME = OUT / "runtime"
digest = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
read = lambda name: json.loads((OUT / name).read_text())
write = lambda name, value: (OUT / name).write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n")
old = read("portable-runtime-manifest.json")
plan = read("portable-cleanup-plan.json")
catalog = read("runtime/catalog.json")
cards = read("ifc-card-validation.json")
drag = read("ifc-card-drag-validation.json")
integration = read("native-integration-validation.json")
drawings = read("no-asset-blends-test.json")
assert cards["status"] == drag["status"] == integration["status"] == "pass"
assert len(integration["products"]) == 40 and integration["native_bonsai_save_reload"]
assert drawings["svg_count"] == 6 and drawings["visual"]["status"] == "pass"
assert not (OUT / "ifc-card-cache-archive.json").exists(), "Already archived: inspect receipt before rerunning"
index = subprocess.check_output(["git", "ls-files", "--stage", "-z"], cwd=ROOT)
assert hashlib.sha256(index).hexdigest() == plan["protected_index_sha256"]
formal = ROOT / "2504 GBTB Yanlord Zhuhai.ifc"
assert digest(formal) == plan["formal_ifc_sha256"]
for relative, sha in plan["protected_files"].items():
    assert digest(ROOT / relative) == sha
old_files = {item["path"]: item for item in old["files"]}
targets = [RUNTIME / "native-assets" / (product["id"] + ".blend") for product in catalog["products"]]
assert len(targets) == 39 and set(targets) == set((RUNTIME / "native-assets").glob("*.blend"))
for path in targets:
    assert path.is_file() and not path.is_symlink()
    assert digest(path) == old_files[str(path.relative_to(RUNTIME))]["sha256"]
for item in old["files"]:
    if item["path"].startswith(("ifc/", "previews/", "approvals/")) or item["path"] == "native-assets/placements.json":
        assert digest(RUNTIME / item["path"]) == item["sha256"], item["path"]
archive = Path(tempfile.mkdtemp(prefix="ifc-asset-caches-", dir="/Users/jiaxinchen/.codex/archives"))
receipt = {"status": "moving", "archive": str(archive), "recoverable": True, "files": []}
def checkpoint():
    write("ifc-card-cache-archive.json", receipt)
    (archive / "recovery-manifest.json").write_text(json.dumps(receipt, ensure_ascii=False, indent=2) + "\n")
checkpoint()
for source in targets:
    item = old_files[str(source.relative_to(RUNTIME))]
    target = archive / source.name
    shutil.move(str(source), str(target))
    assert digest(target) == item["sha256"] and not source.exists()
    receipt["files"].append({**item, "restore_to": str(source), "archived_to": str(target)})
    checkpoint()
receipt.update(status="moved_verified", count=39, bytes=sum(i["bytes"] for i in receipt["files"]))
checkpoint()
addon = ROOT / "pipeline/addons/highpoly_review_library"
with zipfile.ZipFile(OUT / "highpoly_review_library.zip", "w", zipfile.ZIP_DEFLATED) as zipped:
    for path in sorted(addon.glob("*.py")):
        assert digest(path) == digest(RUNTIME / "addon/highpoly_review_library" / path.name)
        zipped.write(path, "highpoly_review_library/" + path.name)
files = [{"path": str(p.relative_to(RUNTIME)), "bytes": p.stat().st_size, "sha256": digest(p)}
         for p in sorted(RUNTIME.rglob("*")) if p.is_file() and "__pycache__" not in p.parts]
manifest = {"status": "pass", "version": "0.7.0", "files": files, "count": len(files),
            "bytes": sum(f["bytes"] for f in files),
            "review_categories": dict(Counter(p["review_category"] for p in catalog["products"]))}
write("portable-runtime-manifest.json", manifest)
assert digest(formal) == plan["formal_ifc_sha256"]
assert subprocess.check_output(["git", "ls-files", "--stage", "-z"], cwd=ROOT) == index
report = {"status": "implemented_verified_pending_user_mouse_acceptance", "version": "0.7.0",
          "task": "Remove all persisted product asset Blends; native cards directly dispatch IFC insertion",
          "courseEvidence": {"mode": "embedded-course-index", "lesson": "052000", "fact": "Save and reopen IFC to prove persistence"},
          "plan": "Keep native drag grid with image/ID-only cards; retain original IFC insertion and feedback; validate; archive old caches",
          "preState": {"runtime_bytes": old["bytes"], "asset_blends": 39, "formal_sha256": plan["formal_ifc_sha256"]},
          "execution": {"cards": cards, "drag_handler": drag, "combined_instances_saved_reloaded": 40, "regression_svgs": 6},
          "persistence": "No card Blend, preferences or formal IFC saved. Native integration tests save disposable IFCs only. Ctrl+S remains Bonsai-owned.",
          "postState": {"runtime_bytes": manifest["bytes"], "asset_blends": 0, "placement_bytes": (RUNTIME / "native-assets/placements.json").stat().st_size,
                        "formal_unchanged": True, "protected_index_unchanged": True, "source_ifcs_materials_previews_approvals_unchanged": True},
          "outputs": {"archive": str(archive), "archived_bytes": receipt["bytes"], "net_bytes_reduced": old["bytes"] - manifest["bytes"]},
          "visual": {"grid": str(OUT / "ifc-card-ui/grid.png"), "real_loading_stage": str(OUT / "ifc-card-ui/green-real-stage.png"),
                     "svg_contact": drawings["visual"]["contact"], "inspection": "Native grid and six SVG rasters inspected; screenshot shows actual 65% import-stage translucent green volume. Camera was reframed for legibility."},
          "boundaries": ["Physical mouse launch was not automated: Computer Use selected another task's Blender, which was not operated.",
                         "Native LOCAL asset grid may also list the user's other Collection assets; only owned cards dispatch IFC insertion.",
                         "39 product cache Blends removed; material dependency Materials.blend deliberately retained.",
                         "No new design approvals; no staging, commit, stash or history edits."],
          "verdict": "Dependency removed and runtime behavior verified at card, handler, persistence and drawing layers; physical mouse acceptance remains with user."}
write("ifc-card-release-validation.json", report)
print(json.dumps({"status": report["status"], "runtime_bytes": manifest["bytes"], "archived_bytes": receipt["bytes"], "archive": str(archive)}))
