#!/usr/bin/env python3
"""Mechanically verify the native Drawing batch without tessellation churn."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import xml.etree.ElementTree as ET
from collections import Counter, deque
from pathlib import Path

import ifcopenshell
import ifcopenshell.util.element
import ifcopenshell.util.placement


PROJECT_ROOT = Path(__file__).resolve().parents[2]
REGISTER = PROJECT_ROOT / "pipeline/decisions/int1-elevation-view-register.csv"
NATIVE_DIR = PROJECT_ROOT / "drawings/elevations/native"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def forward_entities(roots) -> dict[int, str]:
    queue = deque(root for root in roots if root is not None)
    entities: dict[int, str] = {}
    while queue:
        entity = queue.popleft()
        if not hasattr(entity, "id") or entity.id() in entities:
            continue
        entities[entity.id()] = str(entity)
        for value in entity:
            if hasattr(value, "id"):
                queue.append(value)
            elif isinstance(value, (tuple, list)):
                queue.extend(item for item in value if hasattr(item, "id"))
    return entities


def product_geometry_graph(product) -> dict[int, str]:
    return forward_entities((product.ObjectPlacement, product.Representation))


def matrix_max_delta(first, second) -> float:
    return max(
        abs(float(first[row][column]) - float(second[row][column]))
        for row in range(4)
        for column in range(4)
    )


def drawing_document(model: ifcopenshell.file, drawing):
    relations = [
        inverse
        for inverse in model.get_inverse(drawing)
        if inverse.is_a("IfcRelAssociatesDocument")
    ]
    if len(relations) != 1:
        raise RuntimeError(f"{drawing.Name}: expected one document relation")
    return relations[0].RelatingDocument


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--candidate", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()
    source = ifcopenshell.open(arguments.source)
    candidate = ifcopenshell.open(arguments.candidate)
    candidate_by_guid = {
        product.GlobalId: product for product in candidate.by_type("IfcProduct")
    }

    product_records = []
    for source_product in source.by_type("IfcProduct"):
        candidate_product = candidate_by_guid.get(source_product.GlobalId)
        if candidate_product is None:
            product_records.append(
                {"global_id": source_product.GlobalId, "status": "removed"}
            )
            continue
        source_graph = product_geometry_graph(source_product)
        candidate_graph = product_geometry_graph(candidate_product)
        graph_equal = source_graph == candidate_graph
        placement_delta = matrix_max_delta(
            ifcopenshell.util.placement.get_local_placement(
                source_product.ObjectPlacement
            ),
            ifcopenshell.util.placement.get_local_placement(
                candidate_product.ObjectPlacement
            ),
        )
        product_records.append(
            {
                "global_id": source_product.GlobalId,
                "ifc_class": source_product.is_a(),
                "name": source_product.Name,
                "represented": source_product.Representation is not None,
                "geometry_graph_entity_count": len(source_graph),
                "geometry_graph_exact": graph_equal,
                "placement_matrix_max_delta": placement_delta,
                "world_geometry_delta_mm": 0.0 if graph_equal and placement_delta == 0.0 else math.inf,
                "status": "unchanged" if graph_equal and placement_delta == 0.0 else "changed",
            }
        )

    drawings = [
        drawing
        for drawing in candidate.by_type("IfcAnnotation")
        if getattr(drawing, "ObjectType", None) == "DRAWING"
        and (drawing.Name or "").startswith("EL-")
    ]
    drawing_records = []
    for drawing in sorted(drawings, key=lambda item: item.Name):
        pset = ifcopenshell.util.element.get_pset(drawing, "EPset_Drawing")
        document = drawing_document(candidate, drawing)
        svg_path = PROJECT_ROOT / document.Location
        root = ET.parse(svg_path).getroot()
        image_count = sum(element.tag.endswith("image") for element in root.iter())
        drawing_records.append(
            {
                "name": drawing.Name,
                "global_id": drawing.GlobalId,
                "target_view": pset.get("TargetView"),
                "scale": pset.get("Scale"),
                "linework_mode": pset.get("LineworkMode"),
                "has_underlay": pset.get("HasUnderlay"),
                "has_annotation": pset.get("HasAnnotation"),
                "document_location": document.Location,
                "svg_sha256": sha256(svg_path),
                "svg_image_element_count": image_count,
                "svg_has_noninteger_highlights": any(
                    element.attrib.get("id") == "noninteger-highlights"
                    for element in root.iter()
                ),
            }
        )

    changed = [record for record in product_records if record["status"] != "unchanged"]
    represented = [record for record in product_records if record.get("represented")]
    drawing_failures = [
        record
        for record in drawing_records
        if record["target_view"] != "ELEVATION_VIEW"
        or record["scale"] != "1/50"
        or record["has_underlay"] is not False
        or record["svg_image_element_count"] != 0
        or not record["svg_has_noninteger_highlights"]
    ]
    result = {
        "source_ifc_sha256": sha256(arguments.source),
        "candidate_ifc_sha256": sha256(arguments.candidate),
        "source_product_count": len(product_records),
        "source_represented_product_count": len(represented),
        "original_product_geometry_graph_exact_count": sum(
            record.get("geometry_graph_exact", False) for record in product_records
        ),
        "original_product_changed_count": len(changed),
        "original_product_maximum_world_geometry_delta_mm": (
            0.0 if not changed else math.inf
        ),
        "original_product_maximum_placement_matrix_delta": max(
            (record.get("placement_matrix_max_delta", 0.0) for record in product_records),
            default=0.0,
        ),
        "native_elevation_drawing_count": len(drawing_records),
        "linework_mode_counts": dict(
            Counter(record["linework_mode"] for record in drawing_records)
        ),
        "drawing_failure_count": len(drawing_failures),
        "changed_original_products": changed,
        "drawing_failures": drawing_failures,
        "drawings": drawing_records,
        "proof": (
            "Exact equality of every original product ObjectPlacement and "
            "Representation forward STEP subgraph proves 0.0 mm world-geometry change."
        ),
        "pass": not changed and len(drawing_records) == 36 and not drawing_failures,
    }
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(
        json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(
        json.dumps(
            {key: value for key, value in result.items() if key not in {"drawings", "changed_original_products", "drawing_failures"}},
            ensure_ascii=False,
        )
    )
    if not result["pass"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
