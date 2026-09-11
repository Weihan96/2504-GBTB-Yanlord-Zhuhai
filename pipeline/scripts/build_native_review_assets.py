"""Refresh IFC-derived placement metadata only; never persist asset Blends or IFC."""
from pathlib import Path
import importlib.util
import sys
import json
import bpy
from mathutils import Matrix, Vector

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location("highpoly_review_library", ROOT / "pipeline/addons/highpoly_review_library/__init__.py")
addon = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = addon
spec.loader.exec_module(addon)
addon.register()
out = addon._catalog_path.parent / "native-assets"
out.mkdir(exist_ok=True)
placements = {}
for product in addon._catalog["products"]:
    scene = bpy.data.scenes.new("Asset builder " + product["id"])
    collection = addon.body_collection(scene, product)
    addon.place(collection, 0, 0)
    # Remove display labels only from this worker-owned asset, never user scenes.
    for obj in list(collection.objects):
        if obj.get("review_label"):
            bpy.data.objects.remove(obj, do_unlink=True)
    collection.name = "[" + product.get("category_label", product["approval_label"]) + "] " + product["name"]
    matrix = next(iter(collection.objects)).matrix_world.copy()
    points = [obj.matrix_world @ v.co for obj in collection.objects if obj.type == "MESH" for v in obj.data.vertices]
    placements[product["id"]] = {"local_to_asset_m": [list(r) for r in matrix],
        "asset_bounds_m": [[min(p[i] for p in points) for i in range(3)], [max(p[i] for p in points) for i in range(3)]],
        "anchor": "oriented Body bounding-box bottom centre", "source_sha256": product["ifc_sha256"]}
    for obj in list(collection.objects):
        bpy.data.objects.remove(obj, do_unlink=True)
    bpy.data.collections.remove(collection)
    bpy.data.scenes.remove(scene)
    print("IFC_PLACEMENT", product["id"], flush=True)
(out / "placements.json").write_text(json.dumps(placements, indent=2) + "\n")
print("IFC_PLACEMENTS_READY", len(placements), flush=True)
