#!/usr/bin/env python3
"""Extract the approved Falper representative into a Bonsai review IFC."""

from pathlib import Path

import ifcopenshell
import ifcpatch


ROOT = Path(__file__).resolve().parents[2]
SOURCE = (
    ROOT
    / "output/review/highpoly-types/falper-sorgente/"
    "Falper-Sorgente-WFB-derived-drawing-corrected.ifc"
)
OUTPUT = (
    ROOT
    / "output/review/highpoly-types/falper-sorgente/"
    "Falper-Sorgente-WFB-bonsai-isolated.ifc"
)
GLOBAL_ID = "350tdaubr8QP3Cu2YMQZIN"


source = ifcopenshell.open(SOURCE)
isolated = ifcpatch.execute(
    {
        "input": str(SOURCE),
        "file": source,
        "recipe": "ExtractElements",
        "arguments": [GLOBAL_ID],
    }
)
product = isolated.by_guid(GLOBAL_ID)
if product is None:
    raise RuntimeError("Falper representative was not extracted")
identifiers = {
    representation.RepresentationIdentifier
    for representation in product.Representation.Representations
}
if not {"Body", "FalperWFBPlan", "FalperWFBFront", "FalperWFBSide"}.issubset(identifiers):
    raise RuntimeError("isolated Bonsai IFC lost approved representations")
if len([element for element in isolated.by_type("IfcElement") if element.Representation]) != 1:
    raise RuntimeError("Bonsai review IFC must contain exactly one represented element")
OUTPUT.parent.mkdir(parents=True, exist_ok=True)
isolated.write(OUTPUT)
print(OUTPUT)
