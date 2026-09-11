"""Background checks for geometry-free transient cards; no project writes."""
from pathlib import Path
from types import SimpleNamespace
import importlib.util
import json
import sys
import bpy

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location("highpoly_review_library", ROOT / "pipeline/addons/highpoly_review_library/__init__.py")
library = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = library
spec.loader.exec_module(library)
before_meshes = len(bpy.data.meshes)
before_objects = len(bpy.data.objects)
before_collections = len(bpy.data.collections)
preferences = [(p.name, p.path) for p in bpy.context.preferences.filepaths.asset_libraries]
library.register()
from highpoly_review_library import asset_cards
from highpoly_review_library.native_assets import asset_product, drop_matrix
cards = [c for c in bpy.data.collections if c.get(asset_cards.KEY) == asset_cards.OWNER]
assert len(cards) == 39
for card in cards:
    assert not card.objects and not card.children
    assert card.preview and all(card.preview.image_size)
    entry = library.entry(card["product_id"])
    assert asset_product(SimpleNamespace(asset=SimpleNamespace(local_id=card))) == entry
    assert asset_product(SimpleNamespace(active_file=SimpleNamespace(relative_path="Collection/" + card.name))) == entry
    drop_matrix(entry, (1, 2, 3))
other = bpy.data.collections.new("Other user's collection")
other.asset_mark()
assert asset_product(SimpleNamespace(asset=SimpleNamespace(local_id=other))) is None
assert len(bpy.data.objects) == before_objects and len(bpy.data.meshes) == before_meshes
asset_cards.rebuild(library)
assert len([c for c in bpy.data.collections if c.get(asset_cards.KEY) == asset_cards.OWNER]) == 39
assert len(bpy.data.collections) == before_collections + 40
library.unregister()
assert not [c for c in bpy.data.collections if c.get(asset_cards.KEY) == asset_cards.OWNER]
assert other.asset_data
assert len(bpy.data.collections) == before_collections + 1
assert preferences == [(p.name, p.path) for p in bpy.context.preferences.filepaths.asset_libraries]
report = {"status": "pass", "cards": 39, "objects_added": 0, "meshes_added": 0,
          "previews": 39, "direct_product_id_resolution": True, "placements": 39,
          "other_asset_not_insertable": True, "preferences_unchanged": True,
          "rebuild_and_unregister": True}
(ROOT / "output/review/approved-product-library/ifc-card-validation.json").write_text(json.dumps(report, indent=2) + "\n")
print(json.dumps(report))
