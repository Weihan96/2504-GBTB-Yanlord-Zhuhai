#!/usr/bin/env python3
"""Extract the single Molteni 505 UP project instance for Bonsai review."""

import hashlib
from pathlib import Path

import ifcopenshell
import ifcpatch


ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "2504 GBTB Yanlord Zhuhai.ifc"
OUTPUT = ROOT / "output/review/highpoly-types/505-up-v1-lp-s/Molteni-505-UP-V1-LP-S-bonsai-isolated.ifc"
GLOBAL_ID = "19MpdkWqXC7uhUNhLQgrce"
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
    raise RuntimeError("505 UP representative was not extracted")
body_representations = [
    representation
    for representation in product.Representation.Representations
    if representation.RepresentationIdentifier == "Body"
]
if not any(representation.ContextOfItems.TargetView == "MODEL_VIEW" for representation in body_representations):
    raise RuntimeError("isolated 505 UP IFC lost its actual MODEL_VIEW Body representation")
represented = [element for element in isolated.by_type("IfcElement") if element.Representation]
if len(represented) != 1:
    raise RuntimeError("Bonsai review IFC must contain exactly one represented element")
OUTPUT.parent.mkdir(parents=True, exist_ok=True)
isolated.write(OUTPUT)
if sha256(SOURCE) != FORMAL_SHA256:
    OUTPUT.unlink(missing_ok=True)
    raise RuntimeError("formal IFC bytes changed during isolation")
print(OUTPUT)
