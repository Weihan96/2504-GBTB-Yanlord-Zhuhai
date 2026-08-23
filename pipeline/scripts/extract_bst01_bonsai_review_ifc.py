#!/usr/bin/env python3
"""Extract one Baxter Ninfea BST01 IFC Body for Bonsai camera rendering."""

import hashlib
from pathlib import Path
import ifcopenshell
import ifcpatch

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "2504 GBTB Yanlord Zhuhai.ifc"
OUTPUT = ROOT / "output/review/highpoly-types/bst01/Baxter-Ninfea-BST01-bonsai-isolated.ifc"
GLOBAL_ID = "3eic1dzkn5heTIn4PhF37v"
FORMAL_SHA256 = "7a50b87e8f48a7c2bfbaa9f1dabdfca7155fa684b2325c66aa1ae15c8ab0a25c"

def sha256(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''): h.update(block)
    return h.hexdigest()

if sha256(SOURCE) != FORMAL_SHA256: raise RuntimeError("formal IFC hash mismatch")
source=ifcopenshell.open(SOURCE)
isolated=ifcpatch.execute({"input":str(SOURCE),"file":source,"recipe":"ExtractElements","arguments":[GLOBAL_ID]})
product=isolated.by_guid(GLOBAL_ID)
if product is None: raise RuntimeError("BST01 representative was not extracted")
bodies=[r for r in product.Representation.Representations if r.RepresentationIdentifier=="Body"]
if not any(r.ContextOfItems.TargetView=="MODEL_VIEW" for r in bodies): raise RuntimeError("isolated BST01 IFC lost its actual MODEL_VIEW Body representation")
if len([e for e in isolated.by_type("IfcElement") if e.Representation]) != 1: raise RuntimeError("Bonsai review IFC must contain exactly one represented element")
OUTPUT.parent.mkdir(parents=True,exist_ok=True)
isolated.write(OUTPUT)
if sha256(SOURCE) != FORMAL_SHA256:
    OUTPUT.unlink(missing_ok=True)
    raise RuntimeError("formal IFC bytes changed during isolation")
print(OUTPUT)
