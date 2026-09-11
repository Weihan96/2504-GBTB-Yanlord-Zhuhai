"""Record performed visual inspection; never modifies IFC or approval records."""
import json
import migrate_bed01 as task

report = json.loads(task.REPORT.read_text())
assert report["verdict"] == "validated_pending_visual_check"
assert report["attachment"]["linework_geometry_source"] == "pure_product_ifc_representations"
assert task.pkg.sha256(task.FORMAL) == task.FORMAL_HASH
report["visual"] = {
    "method": "inspected three Inkscape-rendered scene PNGs with view_image",
    "observations": ["Plan blue bed linework and grey master-bedroom context visible without target clipping.",
                     "Front blue headboard and base visible in the same approved elevation.",
                     "Side blue headboard, mattress, base and feet visible; approved family linework unchanged."],
    "verdict": "pass"}
report["fresh_pure_session"] = {
    "launcher": "bonsai-launcher", "pid": 62301, "port": 9885,
    "log": "/var/folders/rz/d8p6s4y50ws53rd150nm2s3r0000gn/T/codex-task-blender-501/1ef36077eac34a30b1d4-1788804811626.log",
    "readiness": "TASK_BLENDER_IFC_READY", "blend_saved": False, "has_blend_warning": False,
    "ifc": task.pkg.record(task.SINGLE), "visible_meshes": ["IfcFurniture/Furniture"],
    "material_dependencies_resolve": True,
    "display_only_changes": "Default startup Cube and furniture type preview hidden; no IFC changes."}
report["transport_note"] = "Second scene request exceeded bridge 120-second response deadline after starting. The completed disk report, pure-IFC linework-source marker, saved/reloaded scene, independent schema/SVG checks and subsequent live bridge query prove completion. It was not blindly retried."
report["verdict"] = "pass"
task.write(task.REPORT, report)
for name in ("manifest.json", "handoff.json"):
    path = task.OUT / name
    data = json.loads(path.read_text())
    data["status"] = "complete"
    if name == "manifest.json":
        data["single_product_views"] = {v: str(task.PRODUCT / f"official-dwg-{v}.svg") for v in ("plan", "front", "side")}
        data["existing_body_iso_preview"] = task.pkg.record(task.PRODUCT / "bonsai-camera-iso.png")
        data["note_preview"] = "Existing actual Body saved-camera preview retained; Body fingerprint unchanged. Library-specific thumbnails may be regenerated without changing IFC."
    task.write(path, data)
print("BED01 pure package: pass")
