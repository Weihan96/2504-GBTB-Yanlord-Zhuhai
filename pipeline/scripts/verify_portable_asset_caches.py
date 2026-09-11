"""Run in background Blender against a relocated portable runtime."""
from pathlib import Path
import argparse
import json
import sys
import bpy

parser = argparse.ArgumentParser()
parser.add_argument("--runtime", required=True)
parser.add_argument("--output", required=True)
args = parser.parse_args(sys.argv[sys.argv.index("--") + 1:])
runtime = Path(args.runtime)
catalog = json.loads((runtime / "catalog.json").read_text())
placements = json.loads((runtime / "native-assets/placements.json").read_text())
results = []
for product in catalog["products"]:
    slug = product["id"]
    path = runtime / "native-assets" / (slug + ".blend")
    with bpy.data.libraries.load(str(path), link=False, assets_only=True) as (source, target):
        assert len(source.collections) == 1, (slug, source.collections)
        target.collections = source.collections
    collection = target.collections[0]
    assert collection and collection.asset_data and len(collection.all_objects) > 0, slug
    assert collection.preview and all(collection.preview.image_size), (slug, "missing asset thumbnail")
    assert placements[slug]["source_sha256"] == product["ifc_sha256"]
    results.append({"id": slug, "objects": len(collection.all_objects),
                    "thumbnail_size": list(collection.preview.image_size), "status": "pass"})
    bpy.data.collections.remove(collection)
    bpy.data.orphans_purge(do_recursive=True)
assert len(results) == 39
Path(args.output).write_text(json.dumps({"status": "pass", "results": results}, indent=2) + "\n")
print("PORTABLE_ASSET_CACHES_PASS", flush=True)
