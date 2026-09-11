"""Regression contracts for the native, geometry-free card drag surface."""
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[2]
ADDON = ROOT / "pipeline/addons/highpoly_review_library"


class CardContracts(unittest.TestCase):
    def test_ifc_drag_has_no_blend_path_dispatch(self):
        source = (ADDON / "native_assets.py").read_text()
        self.assertNotIn(".blend/Collection/", source)
        self.assertIn("product_for_id", source)
        self.assertIn("prepare_source", source)
        self.assertNotIn("libs.new", source)

    def test_native_drag_grid_is_kept(self):
        source = (ADDON / "asset_ui.py").read_text()
        self.assertIn("template_asset_view", source)
        self.assertIn("drag_operator='review_library.drag_ifc'", source)

    def test_cards_never_generate_geometry_or_save(self):
        source = (ADDON / "asset_cards.py").read_text()
        for forbidden in ("meshes.new", "objects.new", "save_as_mainfile", "save_project", "libraries.load"):
            self.assertNotIn(forbidden, source)
        self.assertIn("undo_post", source)
        self.assertIn("redo_post", source)

    def test_placement_builder_no_longer_recreates_blend_cache(self):
        source = (ROOT / "pipeline/scripts/build_native_review_assets.py").read_text()
        self.assertNotIn("save_as_mainfile", source)
        self.assertIn("placements.json", source)
