#!/usr/bin/env python3
"""Extract the single electric flue check-valve proxy for Bonsai rendering."""

import hashlib
import subprocess
from pathlib import Path

import ifcopenshell

try:
    import ifcpatch
except ModuleNotFoundError:
    script = Path(__file__).resolve()
    blender = Path("/Applications/Blender.app/Contents/MacOS/Blender")
    if not blender.is_file():
        raise RuntimeError("ifcpatch is unavailable and Blender/Bonsai fallback was not found")
    expression = (
        "import addon_utils,runpy; "
        "addon_utils.enable('bonsai_bridge', default_set=False, persistent=False); "
        f"runpy.run_path({str(script)!r}, run_name='__main__')"
    )
    raise SystemExit(subprocess.run([str(blender), "--background", "--python-expr", expression], check=False).returncode)


ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "2504 GBTB Yanlord Zhuhai.ifc"
OUTPUT = ROOT / "output/review/highpoly-types/electric-flue-check-valve/Electric-flue-check-valve-bonsai-isolated.ifc"
GLOBAL_ID = "1faflkXXH6M9cnYPE9Liir"
FORMAL_SHA256 = "7a50b87e8f48a7c2bfbaa9f1dabdfca7155fa684b2325c66aa1ae15c8ab0a25c"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


if sha256(SOURCE) != FORMAL_SHA256:
    raise RuntimeError("formal IFC hash mismatch")
source = ifcopenshell.open(SOURCE)
isolated = ifcpatch.execute({"input": str(SOURCE), "file": source, "recipe": "ExtractElements", "arguments": [GLOBAL_ID]})
product = isolated.by_guid(GLOBAL_ID)
if product is None:
    raise RuntimeError("electric flue check-valve proxy was not extracted")
model_body = [
    representation
    for representation in product.Representation.Representations
    if representation.RepresentationIdentifier == "Body"
    and representation.ContextOfItems.TargetView == "MODEL_VIEW"
]
if len(model_body) != 1:
    raise RuntimeError("isolated proxy lost its actual MODEL_VIEW Body representation")
represented_products = [
    element
    for element in isolated.by_type("IfcElement")
    if element.Representation and not element.is_a("IfcOpeningElement")
]
represented_openings = [
    element
    for element in isolated.by_type("IfcOpeningElement")
    if element.Representation
]
if represented_products != [product] or len(represented_openings) != 2:
    raise RuntimeError("isolated review IFC must retain one product and its two Boolean openings")
OUTPUT.parent.mkdir(parents=True, exist_ok=True)
isolated.write(OUTPUT)
if sha256(SOURCE) != FORMAL_SHA256:
    OUTPUT.unlink(missing_ok=True)
    raise RuntimeError("formal IFC bytes changed during isolation")
print(OUTPUT)
