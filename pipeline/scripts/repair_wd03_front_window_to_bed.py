#!/usr/bin/env python3
"""Regenerate only WD03 Front from the guest-bedroom window toward the bed."""
import contextlib
import json
import sys
import traceback
from datetime import datetime, timezone
from pathlib import Path

import addon_utils
import bpy
import ifcopenshell.api.geometry
import ifcopenshell.util.placement
from mathutils import Matrix, Vector

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "pipeline/scripts"))
import create_wd03_wardrobe_scene_drawings as wd
from bonsai import tool

OUT = wd.PRODUCT_DIR / "bonsai-drawings/wardrobe"
AUDIT = OUT / "WD03-FRONT-window-to-bed-evidence.json"
WINDOW = "2nisU4bAr9MfEP9b7LwE2g"
BED = "1i_pqgLv9A7uuV7MjaArBW"
DRAWING = "1P1egc73v5kOpl_BXKYZ62"
ANNOTATION = "2uAkszkbj8YBoncHNaFrYK"


def record(path):
    return {"path": str(path), "sha256": wd.sha256(path), "bytes": path.stat().st_size}


def main():
    approval = json.loads(wd.APPROVAL.read_text())
    assert approval["derived_ifc_write_allowed"] and not approval["formal_authoritative_ifc_write_allowed"]
    assert wd.sha256(wd.FORMAL_IFC) == wd.FORMAL_SHA256
    frozen = [wd.CANDIDATE, wd.PRODUCT_DIR / "plan.svg", wd.PRODUCT_DIR / "front.svg", wd.PRODUCT_DIR / "side.svg"]
    frozen += [OUT / f"WD03-WARDROBE-{v}.svg" for v in ("PLAN", "SIDE")]
    before = {str(p): wd.sha256(p) for p in frozen}
    derived_before = wd.sha256(wd.DERIVED_IFC)
    with contextlib.suppress(Exception):
        addon_utils.disable("bl_ext.user_default.project_control", default_set=False, handle_error=None)
    assert bpy.ops.bim.load_project(filepath=str(wd.DERIVED_IFC), should_start_fresh_session=True) == {"FINISHED"}
    model = tool.Ifc.get()
    drawing = model.by_guid(DRAWING)
    camera = tool.Ifc.get_object(drawing) or tool.Drawing.import_drawing(drawing)
    annotation = model.by_guid(ANNOTATION)
    assert wd.persisted_path_count(annotation) == 9
    references = {}
    for role, guid in (("guest_bed", BED), ("guest_window", WINDOW), ("wd03", wd.TARGET_GLOBAL_ID)):
        element = model.by_guid(guid)
        obj = tool.Ifc.get_object(element)
        assert obj is not None, role
        references[role] = {"global_id": guid, "type_name": getattr(__import__("ifcopenshell").util.element.get_type(element), "Name", None), "world_bbox_m": wd.shared.world_bbox(obj), "world_placement_mm": ifcopenshell.util.placement.get_local_placement(element.ObjectPlacement).tolist()}
    old_matrix = [list(row) for row in camera.matrix_world]
    # The actual guest window is north of BED02, which is north of WD03.
    window_y = references["guest_window"]["world_placement_mm"][1][3] / 1000
    bed_y = references["guest_bed"]["world_placement_mm"][1][3] / 1000
    wardrobe_y = references["wd03"]["world_placement_mm"][1][3] / 1000
    assert window_y > bed_y > wardrobe_y
    matrix = Matrix(((-1, 0, 0, camera.location.x), (0, 0, 1, window_y - 0.30), (0, 1, 0, camera.location.z), (0, 0, 0, 1)))
    camera.matrix_world = matrix
    camera.data.clip_end = max(5.2, matrix.translation.y - wardrobe_y + 0.8)
    bpy.context.view_layer.update()
    drawing.Description = "WD03 guest-bedroom Front: camera on the window side, looking toward BED02 and the wardrobe; approved semantic black geometry."
    # Persist the camera placement through the normal IFC geometry API.
    import numpy as np
    ifcopenshell.api.geometry.edit_object_placement(model, product=drawing, matrix=np.array(matrix), is_si=True, should_transform_children=False)
    override = wd.view3d_override()
    with bpy.context.temp_override(**override):
        assert bpy.ops.bim.activate_drawing(drawing=drawing.id(), should_view_from_camera=False) == {"FINISHED"}
    dprops = tool.Drawing.get_document_props()
    dprops.should_use_underlay_cache = False
    dprops.should_use_linework_cache = False
    dprops.should_use_annotation_cache = False
    svg = OUT / "WD03-WARDROBE-FRONT.svg"
    with bpy.context.temp_override(**override):
        result = bpy.ops.bim.create_drawing(print_all=False, open_viewer=False, sync=False)
    assert result == {"FINISHED"} and svg.is_file()
    validation = wd.style_and_inspect(svg, ANNOTATION)
    assert validation["no_body_annotation_duplicate"]
    import bonsai_bridge
    persistence = bonsai_bridge._h_save_ifc_file({"output_path": str(wd.DERIVED_IFC), "overwrite": True, "reload": True})
    reopened = tool.Ifc.get()
    count = wd.persisted_path_count(reopened.by_guid(ANNOTATION))
    assert count == 9
    persisted_matrix = ifcopenshell.util.placement.get_local_placement(reopened.by_guid(DRAWING).ObjectPlacement)
    assert np.allclose(persisted_matrix[:3, :3], np.array(matrix)[:3, :3], atol=1e-6)
    assert np.allclose(persisted_matrix[:3, 3] / 1000, matrix.translation, atol=1e-6)
    bpy.ops.wm.save_as_mainfile(filepath=str(wd.SESSION_BLEND), check_existing=False)
    assert before == {str(p): wd.sha256(p) for p in frozen}
    assert wd.sha256(wd.FORMAL_IFC) == wd.FORMAL_SHA256
    audit = {
        "task": "WD03 Front camera corrected from guest-bedroom window toward BED02 and wardrobe",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "versions": {"blender": bpy.app.version_string, "ifcopenshell": ifcopenshell.version, "bridge": bonsai_bridge.bl_info.get("version")},
        "courseEvidence": {"mode": "embedded-course-index", "lesson": "085000 Introduction to Drawings", "timestamps": ["01:03 active drawing camera", "01:24 move camera boundary", "01:59 Create Drawing"]},
        "plan": ["inspect real IFC world coordinates", "reverse Front camera from window side", "native Create Drawing", "provider save_ifc_file and reload", "verify frozen geometry and Plan/Side"],
        "preState": {"derived_ifc_sha256": derived_before, "front_camera_matrix_m": old_matrix, "frozen_files": before},
        "world_direction_evidence": references,
        "execution": {"generator": "bpy.ops.bim.create_drawing", "result": sorted(result), "mode": "OPENCASCADE", "adapter": "isolated Blender process using installed public bridge save_ifc_file handler; unrelated live MCP scene preserved"},
        "persistence": persistence,
        "postState": {"front_camera_matrix_m": [list(row) for row in matrix], "view_direction_world": [0, -1, 0], "persisted_drawing_matrix_mm": persisted_matrix.tolist(), "persisted_annotation_path_count": count, "formal_ifc_sha256": wd.FORMAL_SHA256},
        "outputs": {"front_svg": {**record(svg), **validation}, "derived_ifc": record(wd.DERIVED_IFC), "session_blend": record(wd.SESSION_BLEND)},
        "tests": {"window_bed_wardrobe_world_order_verified": True, "front_points_from_window_to_bed": True, "candidate_and_plan_side_bytes_unchanged": True, "annotation_path_count_9": True, "target_body_no_duplicate": True},
        "visual": {"status": "awaiting_render_inspection"},
        "review_status": "front_scene_svg_pending_user_review",
        "geometry_approval_quote": "但是这个几何的看上去是OK的。",
        "verdict": "awaiting_visual_verification",
    }
    AUDIT.write_text(json.dumps(audit, ensure_ascii=False, indent=2, default=str) + "\n")
    print(json.dumps({"evidence": str(AUDIT), "svg": str(svg), "generated": True}))


if __name__ == "__main__":
    try:
        main()
    except Exception:
        (OUT / "WD03-FRONT-window-to-bed-error.log").write_text(traceback.format_exc())
        raise
    finally:
        bpy.app.timers.register(lambda: bpy.ops.wm.quit_blender() and None, first_interval=2.0)
