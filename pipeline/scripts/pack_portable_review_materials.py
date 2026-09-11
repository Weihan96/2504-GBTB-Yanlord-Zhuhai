"""Pack external material images into the runtime copy, never original assets."""
from pathlib import Path
import json
import hashlib
import bpy
import ifcopenshell

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT/"output/review/approved-product-library"
target = OUT/"runtime/ifc/Materials.blend"
before = hashlib.sha256(target.read_bytes()).hexdigest()
source_library = ROOT/"Materials.blend"
assert hashlib.sha256(source_library.read_bytes()).hexdigest() == "18a1e92531bc56c3938a16c64189bc61da83354529c22a9c4d63f0f05cfd69a2"
required = set()
for product in json.loads((OUT/"runtime/catalog.json").read_text())["products"]:
    model = ifcopenshell.open(str(OUT/"runtime"/product["ifc_path"]))
    required.update(s.Identification.removeprefix("materials/") for s in model.by_type("IfcExternallyDefinedSurfaceStyle"))
# Read the hash-identical original library so //Textures paths resolve against
# their actual original directory. Only its referenced materials are retained.
with bpy.data.libraries.load(str(source_library),link=False) as (src,dst):
    assert required.issubset(src.materials), required-set(src.materials)
    dst.materials = sorted(required)
materials = set(dst.materials)
images = set()
seen = set()
def visit(tree):
    if tree is None or tree in seen:
        return
    seen.add(tree)
    for node in tree.nodes:
        image = getattr(node,"image",None)
        if image and image.source in ("FILE","TILED","SEQUENCE","MOVIE"):
            images.add(image)
        visit(getattr(node,"node_tree",None))
for material in materials:
    visit(material.node_tree)
records,missing = [],[]
for image in images:
    path = Path(bpy.path.abspath(image.filepath)).resolve()
    if image.packed_file or list(image.packed_files):
        records.append({"image":image.name,"already_packed":True})
    elif path.is_file():
        image.pack()
        assert image.packed_file or list(image.packed_files)
        records.append({"image":image.name,"source":str(path),"bytes":path.stat().st_size,"sha256":hashlib.sha256(path.read_bytes()).hexdigest()})
    else:
        missing.append({"image":image.name,"source":str(path)})
report = {"before_sha256":before,"images":records,"missing":missing,"materials":len(materials)}
if not missing:
    bpy.data.libraries.write(str(target),materials,fake_user=True,compress=True)
    report["after_sha256"] = hashlib.sha256(target.read_bytes()).hexdigest()
    report["status"] = "packed"
else:
    report["status"] = "missing_external_textures"
(OUT/"portable-materials.json").write_text(json.dumps(report,indent=2)+"\n")
print(json.dumps({"status":report["status"],"materials":len(materials),"images":len(records),"missing":len(missing)}),flush=True)
