#!/usr/bin/env python3
"""Extract one Gessi316 54038 IFC Body for Bonsai camera rendering."""

import hashlib
import subprocess
import sys
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
OUTPUT = ROOT / "output/review/highpoly-types/gessi316-54038/Gessi316-54038-bonsai-isolated.ifc"
GLOBAL_ID = "245NU$zZL0d9tYTVwwBdk$"
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
    raise RuntimeError("Gessi representative was not extracted")
identifiers = {representation.RepresentationIdentifier for representation in product.Representation.Representations}
if "Body" not in identifiers:
    raise RuntimeError("isolated Gessi IFC lost its actual Body representation")
represented = [element for element in isolated.by_type("IfcElement") if element.Representation]
if len(represented) != 1:
    raise RuntimeError("Bonsai review IFC must contain exactly one represented element")
OUTPUT.parent.mkdir(parents=True, exist_ok=True)
isolated.write(OUTPUT)
if sha256(SOURCE) != FORMAL_SHA256:
    OUTPUT.unlink(missing_ok=True)
    raise RuntimeError("formal IFC bytes changed during isolation")
print(OUTPUT)
