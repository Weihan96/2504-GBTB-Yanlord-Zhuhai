"""Finalize independently observed validation; no IFC or approval writes."""
import json
import migrate_miami_e09 as task

report = json.loads(task.REPORT.read_text())
assert report["verdict"] == "validated_pending_visual_check"
assert report["attachment"]["linework_geometry_source"] == "pure_product_ifc_representations"
assert task.pkg.sha256(task.FORMAL) == task.FORMAL_HASH
report["visual"] = {
    "method": "view_image of three new Inkscape-rendered scene PNGs, isolated Side thumbnail and existing actual Body camera render",
    "observations": [
        "Plan blue E09 remains in the approved living-room sofa arrangement; target is not clipped.",
        "Front blue cushions, arm and feet occupy the original approved elevation in grey context.",
        "Side blue arm, cushions and base retain accepted positions; no component force-fit applied.",
        "The source hidden circle uses 259 independently styled short SVG lines. Dash attributes are preserved, but per-segment dash-phase reset makes the circle appear nearly continuous at thumbnail scale; inherited, not corrected by this storage migration."
    ], "verdict": "pass_storage_equivalence"}
report["fresh_pure_session"] = {
    "launcher": "bonsai-launcher", "pid": 64908, "port": 9886,
    "owner": "5c976bee992602cbe91d", "task_id": "migrate_miami_e09_pure_product",
    "log": "/var/folders/rz/d8p6s4y50ws53rd150nm2s3r0000gn/T/codex-task-blender-501/5c976bee992602cbe91d-1788805492419.log",
    "readiness": "TASK_BLENDER_IFC_READY", "blend_saved": False, "has_blend_warning": False,
    "physical_elements": 1, "annotations": 0, "visible_meshes": ["IfcFurniture/MiamiSoft E09"],
    "session_cleanup": "Only task-owned PID 63910 was replaced after persisted scene validation; known temporary display state discarded. No unrelated window touched."
}
report["approval_completion_evidence"] = {
    "path": "pipeline/decisions/highpoly-subtask-completion-2026-09-05.json",
    "product": "miamisoft-e09", "user_evidence": "通过，可以关闭负责它的子代理了"}
report["verdict"] = "pass"
task.write(task.REPORT, report)
for name in ("manifest.json", "handoff.json"):
    path = task.OUT / name
    data = json.loads(path.read_text())
    data["status"] = "complete"
    if name == "manifest.json":
        data["single_product_views"] = {v: task.pkg.record(task.OUT / f"MIAMI-E09-SINGLE-{v.upper()}.svg") for v in ("plan", "front", "side")}
        data["library_previews"] = {v: task.pkg.record(task.OUT / f"MIAMI-E09-SINGLE-{v.upper()}.png") for v in ("plan", "front", "side")}
        data["library_previews"]["iso"] = task.pkg.record(task.PRODUCT / "bonsai-camera-iso.png")
        data["existing_body_iso_preview"] = data["library_previews"]["iso"]
        data["preview_provenance"] = {
            "2d": "Target linework isolated from newly generated real Bonsai scene SVGs; no geometry change; reframed only.",
            "3d": "Existing actual IFC Body saved-camera render; target Body fingerprint unchanged.",
            "3d_manifest": task.pkg.record(task.PRODUCT / "bonsai-review-manifest.json")}
        data["known_limitations"].append(report["visual"]["observations"][-1])
    task.write(path, data)
print("Miami E09 pure package: pass")
