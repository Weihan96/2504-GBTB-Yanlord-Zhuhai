"""Record completed human-agent visual inspection of the migration evidence."""
import json
from pathlib import Path
import migrate_trap01 as m
import review_product_package as pkg

report = json.loads(m.REPORT.read_text())
assert report["verdict"] == "validated_pending_visual_inspection"
viewport = json.loads((m.OUT / "bonsai-product-viewport.json").read_text())
assert viewport["viewport"]["objects_in_view_total"] == 1
assert viewport["viewport"]["objects_in_view"][0]["global_id"] == m.GUID
assert viewport["active_ifc"] == str(m.SINGLE)
assert viewport["ifc_sha256"] == pkg.sha256(m.SINGLE)
assert not viewport["blend_saved"]
visual = {
    "inspected_by": "TRAP01 migration subagent",
    "method": "Opened all three regenerated PNGs and the cleanly reloaded pure-IFC Bonsai viewport image with view_image",
    "plan": "Complete blue installed assembly visible: left straight outlet, central joint and concentric distal inlet circles; gray project context retained; no cabinet occluder covers product.",
    "front": "Complete blue end-on assembly visible from upper collar to lower trap bowl and support tabs; gray context retained.",
    "side": "Complete blue curved upper inlet tube, coupling, outlet and lower trap body visible; not substituted by a generic silhouette.",
    "pure_body": "Only the selected TRAP01 3D product is visible after public bridge clear-and-reload; complete upper elbow and lower body fit the view.",
    "scene_previews": [r["preview"] for r in report["approved_scene_comparison"]],
    "pure_viewport": viewport["image"],
    "pass": True,
}
report.update({"visual": visual, "verdict": "pass", "stage": "pure_product_and_scene_validation_complete",
    "clean_pure_ui_reload": {"method": "bonsai_bridge._reload_ifc_project", "visible_physical_products": 1, "guid": m.GUID, "persisted_evidence": "bonsai-product-viewport.json"},
    "execution_notes": [
        "Initial raw bridge command execute_blender_code was rejected without mutation; installed protocol uses execute_code.",
        "Long scene request exceeded provider response deadline, but completed on Blender main thread. Completion was established from saved reload report, independent IFC/SVG checks and a fresh live pure-project query; it was not inferred from the timeout.",
        "Non-fresh native load retained orphan Blender scene objects when returning from temporary project. No such objects were saved into pure IFC. Final accepted display was cleared and reloaded through the public provider; visible physical product count is one.",
        "Native Bonsai serializer consolidated microscopic SVG segments. Maximum physical drawing deviation is 0.001666 mm; persisted approved curve geometry is unchanged.",
    ], "formal_sha256_final": pkg.sha256(m.FORMAL), "cleanup_performed": False})
assert report["formal_sha256_final"] == m.FORMAL_SHA
m.write(m.REPORT, report)
manifest_path = m.OUT / "package-manifest.json"
manifest = json.loads(manifest_path.read_text())
manifest["status"] = "complete_verified"
manifest["visual"] = visual
manifest["validation_sha256"] = pkg.sha256(m.REPORT)
m.write(manifest_path, manifest)
m.write(m.OUT / "handoff.json", {
    "profile_key": "trap01", "status": "complete_verified", "approved_source_version": "official-detail-v4",
    "persistent_ifc": pkg.record(m.SINGLE), "scene_recipe": pkg.record(m.RECIPE),
    "package_manifest": pkg.record(manifest_path), "validation": pkg.record(m.REPORT),
    "formal_ifc_unchanged": True, "formal_sha256": m.FORMAL_SHA,
    "only_product_directory_modified": str(m.OUT), "staged_baseline_untouched": True,
    "files_staged_by_this_agent": 0, "old_files_deleted": 0, "formal_ifc_written": False,
    "bonsai_session": {"task_id": "migrate-trap01-pure-product", "pid": 61461, "port": 9884, "active_ifc": str(m.SINGLE), "clean_reload_verified": True},
    "cleanup_requires_confirmation": "cleanup-proposal.json", "scene_svg_names": [f"TRAP01-SCENE-{v.upper()}.svg" for v in m.VIEWS],
})
print(json.dumps({"status": "complete_verified", "pure_bytes": m.SINGLE.stat().st_size, "formal_unchanged": True}))
