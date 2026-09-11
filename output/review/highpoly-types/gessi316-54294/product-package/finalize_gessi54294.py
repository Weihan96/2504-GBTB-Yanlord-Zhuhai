"""Record independently observed storage-migration acceptance evidence."""
import json
import ifcopenshell
import numpy as np
import migrate_gessi54294 as task


def finish(session):
    report = json.loads(task.REPORT.read_text())
    assert report["verdict"] in ("validated_pending_visual_check", "pass")
    assert task.pkg.sha256(task.FORMAL) == task.FORMAL_HASH
    formal = ifcopenshell.open(str(task.FORMAL))
    scene = ifcopenshell.open(report["temporary_project"]["path"])
    formal_positions = {e.GlobalId: task.pkg.placement(e).tolist() for e in formal.by_type("IfcElement")}
    scene_positions = {e.GlobalId: task.pkg.placement(e).tolist() for e in scene.by_type("IfcElement")}
    assert formal_positions == scene_positions
    report["independent_validation"]["all_formal_product_placements_unchanged"] = True
    report["independent_validation"]["placement_verified_product_count"] = len(formal_positions)
    report["visual"] = {
        "method": "view_image of new scene PNGs and isolated product PNGs",
        "observations": [
            "Plan retains the approved three-hole group in the actual main-bathroom context.",
            "Front retains the de-textured blue outlines with no extra positional shift.",
            "Side retains the approved finished-wall face, with the occluding substrate excluded only by its Drawing settings.",
            "Single-product previews isolate the persisted approved linework; no redraw or geometry repair."
        ], "verdict": "pass_storage_equivalence"}
    report["fresh_pure_session"] = session
    report["verdict"] = "pass"
    task.write(task.REPORT, report)
    for name in ("manifest.json", "handoff.json"):
        path = task.OUT / name
        data = json.loads(path.read_text())
        data["status"] = "complete"
        if name == "manifest.json":
            data["single_product_views"] = {v: task.pkg.record(task.OUT / f"GESSI54294-SINGLE-{v.upper()}.svg") for v in ("plan", "front", "side")}
            data["library_previews"] = {v: task.pkg.record(task.OUT / f"GESSI54294-SINGLE-{v.upper()}.png") for v in ("plan", "front", "side")}
            iso = task.PRODUCT / "bonsai-camera-iso.png"
            data["library_previews"]["iso"] = task.pkg.record(iso)
            data["preview_provenance"] = {"2d": "Approved product lines isolated from regenerated real Bonsai scene SVGs; reframing plus thumbnail-only stroke width/round-cap adaptation. No geometry change; approved scene SVG styling unchanged.",
                                           "3d": "Existing actual Body camera render; all original Body representations preserved unchanged."}
            data["handle_texture_detail_path_count"] = 0
            data["side_annotation_translation_mm"] = [0, 0, 0]
            data["approved_finished_wall_residual_mm"] = 0.02276718446054815
        task.write(path, data)
    cleanup = json.loads((task.OUT / "cleanup-proposal.json").read_text())
    cleanup["temporary_directories"] = sorted(set(cleanup["temporary_directories"] + [
        "/var/folders/rz/d8p6s4y50ws53rd150nm2s3r0000gn/T/gessi54294-pure-package-scene-31y6j538"]))
    task.write(task.OUT / "cleanup-proposal.json", cleanup)
    print("Gessi54294 pure package: pass")
