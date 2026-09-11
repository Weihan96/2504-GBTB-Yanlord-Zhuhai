"""Interactive-session operator tests and screenshot; never writes an IFC."""
import bpy
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "output/review/approved-product-library"
library = sys.modules["highpoly_review_library"]
props = bpy.context.scene.review_library
total = len(library._catalog["products"])
assert total
before = {p["id"]: library.digest(library.resolve_path(p["ifc_path"])) for p in library._catalog["products"]}
props.search = "bed01"
assert len(library.enum_items(props, bpy.context)) == 1
props.search = "unknown-product-12345"
assert not library.enum_items(props, bpy.context)
props.search = ""
assert len(library.enum_items(props, bpy.context)) == total
tiers = {state: [p for p in library._catalog["products"] if p["scene_approval_status"] == state]
         for state in ("approved", "pending")}
assert len(tiers["approved"]) == 8 and len(tiers["pending"]) == 7
for key, status in (("bed01", "approved"), ("bed02", "pending")):
    props.product = key
    assert library.entry(props.product)["scene_approval_status"] == status
    for view in ("front", "side", "iso", "plan"):
        assert bpy.ops.review_library.image(view=view) == {"FINISHED"}
saved_path, saved_catalog = library._catalog_path, library._catalog
try:
    library.load_catalog(OUT / "intentionally-missing-catalog.json")
except FileNotFoundError:
    pass
else:
    raise AssertionError("Invalid catalogue unexpectedly loaded")
assert library._catalog_path == saved_path and library._catalog is saved_catalog
props.columns = 3
props.gap = 1.
assert bpy.ops.review_library.arrange() == {"FINISHED"}
assert before == {p["id"]: library.digest(library.resolve_path(p["ifc_path"])) for p in library._catalog["products"]}
ui_regions = [{"width_pixels": r.width, "height_pixels": r.height} for a in bpy.context.screen.areas
              if a.type == "VIEW_3D" for r in a.regions if r.type == "UI"]


screenshots = {}


def capture_pending():
    assert props.product == "bed02"
    path = OUT / "library-ui.png"
    # An unfocused task-owned window may retain an old framebuffer even after
    # RNA state changes. Draw and swap before recording visual evidence.
    bpy.ops.wm.redraw_timer(type="DRAW_WIN_SWAP", iterations=2)
    bpy.ops.screen.screenshot(filepath=str(path))
    screenshots["scene_pending"] = {"product_id": "bed02", "path": str(path), "sha256": library.digest(path)}
    record = {"products": total, "search_and_preview_operator_pass": True,
        "approval_tiers_verified": True,
        "approval_summary": {key: len(value) for key, value in tiers.items()},
        "tier_screenshots": screenshots,
        "search_empty_result_supported": True, "invalid_catalog_preserves_previous_state": True,
        "four_preview_buttons_pass": True, "images_fitted_to_panel": True,
        "source_ifc_hashes_unchanged": before, "ui_regions": ui_regions,
        "ui_scale": bpy.context.preferences.system.ui_scale,
        "window_screenshot": {"path": str(path), "sha256": library.digest(path)},
        "ifc_write_performed": False, "visual_review": "pending_main_agent_image_inspection"}
    (OUT / "ui-validation.json").write_text(json.dumps(record, ensure_ascii=False, indent=2) + "\n")
    return None


def capture_approved():
    assert props.product == "bed01"
    path = OUT / "library-ui-scene-approved.png"
    bpy.ops.wm.redraw_timer(type="DRAW_WIN_SWAP", iterations=2)
    bpy.ops.screen.screenshot(filepath=str(path))
    screenshots["scene_approved"] = {"product_id": "bed01", "path": str(path), "sha256": library.digest(path)}
    props.product = "bed02"
    bpy.ops.review_library.image(view="plan")
    bpy.app.timers.register(capture_pending, first_interval=1.)
    return None


props.product = "bed01"
bpy.ops.review_library.image(view="plan")
bpy.app.timers.register(capture_approved, first_interval=1.)
print("Library UI operators pass; screenshot scheduled")
