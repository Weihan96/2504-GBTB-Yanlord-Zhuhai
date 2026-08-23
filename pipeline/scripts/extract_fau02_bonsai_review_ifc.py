#!/usr/bin/env python3
"""Extract one FAU02 IFC Body for the actual Bonsai camera review."""

import hashlib
import sys
from pathlib import Path

import ifcopenshell

try:
    import ifcpatch
except ModuleNotFoundError:
    bonsai_packages = (
        Path.home()
        / "Library/Application Support/Blender/4.5/extensions/.local/lib/python3.11/site-packages"
    )
    sys.path.append(str(bonsai_packages))
    import ifcpatch


ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "2504 GBTB Yanlord Zhuhai.ifc"
OUTPUT = ROOT / "output/review/highpoly-types/fau02/FAU02-bonsai-isolated.ifc"
GLOBAL_ID = "36ZX3QPyD7SvlXsDKMP8rY"
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
isolated = ifcpatch.execute(
    {"input": str(SOURCE), "file": source, "recipe": "ExtractElements", "arguments": [GLOBAL_ID]}
)
product = isolated.by_guid(GLOBAL_ID)
if product is None:
    raise RuntimeError("FAU02 representative was not extracted")
model_body = [
    representation
    for representation in product.Representation.Representations
    if representation.RepresentationIdentifier == "Body"
    and representation.ContextOfItems.TargetView == "MODEL_VIEW"
]
if len(model_body) != 1:
    raise RuntimeError("isolated FAU02 IFC lost its actual MODEL_VIEW Body representation")
represented = [element for element in isolated.by_type("IfcElement") if element.Representation]
if len(represented) != 1:
    raise RuntimeError("Bonsai review IFC must contain exactly one represented element")
OUTPUT.parent.mkdir(parents=True, exist_ok=True)
isolated.write(OUTPUT)
if sha256(SOURCE) != FORMAL_SHA256:
    OUTPUT.unlink(missing_ok=True)
    raise RuntimeError("formal IFC bytes changed during isolation")
print(OUTPUT)
