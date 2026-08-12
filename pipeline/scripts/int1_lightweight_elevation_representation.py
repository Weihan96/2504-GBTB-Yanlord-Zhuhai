#!/usr/bin/env python3
"""Add curve-only ELEVATION_VIEW bodies for deep-BRep drawing elements.

The detailed MODEL_VIEW Body remains untouched.  Bonsai prioritises a
Model/Body/ELEVATION_VIEW representation when producing elevation drawings,
so these small wire profiles prevent OpenCASCADE from serialising every edge
of the original curved or tessellated body.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
from typing import Iterable

import ifcopenshell
import ifcopenshell.api.context
import ifcopenshell.geom
import ifcopenshell.util.representation


EXPECTED_SOURCE_SHA256 = (
    "e17a3b6a068508140df555960f5bb72d2b7cc00f3c40adb0ccd10754ecd388be"
)
TARGETS = {
    "0V9CnYT0n3GvkLmXhxvv$u": ("DRA01", "pipe_flange"),
    "0YgqsahOH4FvZ_nLwl_Tlf": ("DRA01", "pipe_flange"),
    "2wViHGrwX2hOYCb3MDNb$9": ("DRA01", "pipe_flange"),
    "2S2c498tb7$gzdukjhCGVQ": ("Geberit 154.446.KS.1", "linear_drain"),
    "2iKOL78$H0N9Yd9$ky3pW4": ("Gessi316 54294", "sanitary_fitting"),
    "3hgNkx97vCTOC2eewCpMNk": ("Geberit Duofix Sigma", "duofix_frame"),
    "3pvAlH5C14v8uVEJ1LmK8M": ("Geberit Duofix Sigma", "duofix_frame"),
    "2xmcLzu1rDTeMzRuNxPDyE": ("HIMA01", "cabinet"),
    "3IQBEqO5vDI8Z9k1Ltge_N": ("BED01", "bed"),
    "1MzM8Ms2vFo8KEm503j9w2": ("SIS04", "cabinet"),
    "1O9JRXCI56VRUbpuLJy86Z": (None, "ceiling_panel"),
    "1i_pqgLv9A7uuV7MjaArBW": ("BED02", "bed"),
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def elevation_context(model: ifcopenshell.file):
    existing = ifcopenshell.util.representation.get_context(
        model, "Model", "Body", "ELEVATION_VIEW"
    )
    if existing:
        return existing
    parent = ifcopenshell.util.representation.get_context(model, "Model")
    if not parent:
        raise RuntimeError("IFC has no Model representation context")
    return ifcopenshell.api.context.add_context(
        model,
        context_type="Model",
        context_identifier="Body",
        target_view="ELEVATION_VIEW",
        parent=parent,
    )


def local_bbox_mm(product) -> tuple[float, float, float, float, float, float]:
    settings = ifcopenshell.geom.settings()
    settings.set(settings.USE_WORLD_COORDS, False)
    shape = ifcopenshell.geom.create_shape(settings, product)
    vertices = shape.geometry.verts
    minimum = [min(vertices[index::3]) * 1000 for index in range(3)]
    maximum = [max(vertices[index::3]) * 1000 for index in range(3)]
    return (*minimum, *maximum)


def line(a, b):
    return (tuple(float(value) for value in a), tuple(float(value) for value in b))


def box_edges(bounds) -> list[tuple[tuple[float, ...], tuple[float, ...]]]:
    x0, y0, z0, x1, y1, z1 = bounds
    points = {
        (ix, iy, iz): (x, y, z)
        for ix, x in enumerate((x0, x1))
        for iy, y in enumerate((y0, y1))
        for iz, z in enumerate((z0, z1))
    }
    edges = []
    for key, point in points.items():
        for axis in range(3):
            neighbour = list(key)
            neighbour[axis] += 1
            neighbour = tuple(neighbour)
            if neighbour in points:
                edges.append(line(point, points[neighbour]))
    return edges


def rectangle_xz(x0, x1, y, z0, z1):
    return [
        line((x0, y, z0), (x1, y, z0)),
        line((x1, y, z0), (x1, y, z1)),
        line((x1, y, z1), (x0, y, z1)),
        line((x0, y, z1), (x0, y, z0)),
    ]


def rectangle_yz(x, y0, y1, z0, z1):
    return [
        line((x, y0, z0), (x, y1, z0)),
        line((x, y1, z0), (x, y1, z1)),
        line((x, y1, z1), (x, y0, z1)),
        line((x, y0, z1), (x, y0, z0)),
    ]


def profile_edges(bounds, profile: str):
    x0, y0, z0, x1, y1, z1 = bounds
    dx, dy, dz = x1 - x0, y1 - y0, z1 - z0
    cx, cy = (x0 + x1) / 2, (y0 + y1) / 2
    edges = box_edges(bounds)

    if profile == "pipe_flange":
        # Two square flange rings and an axial centreline.  No circles or
        # surfaces are used: Bonsai therefore exports only these controlled
        # polylines instead of thousands of tessellation edges.
        for z in (z0 + dz * 0.18, z1 - dz * 0.18):
            edges += [
                line((x0, y0, z), (x1, y0, z)),
                line((x1, y0, z), (x1, y1, z)),
                line((x1, y1, z), (x0, y1, z)),
                line((x0, y1, z), (x0, y0, z)),
            ]
        edges.append(line((cx, cy, z0), (cx, cy, z1)))
    elif profile == "linear_drain":
        y = cy
        edges += rectangle_xz(x0 + dx * 0.04, x1 - dx * 0.04, y, z0 + dz * 0.18, z1 - dz * 0.18)
        edges.append(line((cx, y, z0), (cx, y, z1)))
    elif profile == "sanitary_fitting":
        y = cy
        edges += rectangle_xz(x0 + dx * 0.15, x1 - dx * 0.15, y, z0 + dz * 0.18, z1 - dz * 0.18)
        edges += [line((cx, y, z0), (cx, y, z1)), line((x0, y, z0 + dz * 0.5), (x1, y, z0 + dz * 0.5))]
    elif profile == "duofix_frame":
        y = cy
        edges += rectangle_xz(x0 + dx * 0.12, x1 - dx * 0.12, y, z0 + dz * 0.08, z1 - dz * 0.08)
        edges += [
            line((x0, y, z0 + dz * 0.32), (x1, y, z0 + dz * 0.32)),
            line((x0, y, z0 + dz * 0.68), (x1, y, z0 + dz * 0.68)),
            line((cx, y, z0), (cx, y, z1)),
        ]
    elif profile == "cabinet":
        y = cy
        edges += [
            line((x0, y, z0 + dz / 3), (x1, y, z0 + dz / 3)),
            line((x0, y, z0 + 2 * dz / 3), (x1, y, z0 + 2 * dz / 3)),
            line((cx, y, z0), (cx, y, z1)),
        ]
    elif profile == "bed":
        y = cy
        x = cx
        edges += [
            line((x0, y, z0 + dz * 0.48), (x1, y, z0 + dz * 0.48)),
            line((x0, y, z0 + dz * 0.68), (x1, y, z0 + dz * 0.68)),
            line((x, y0, z0 + dz * 0.48), (x, y1, z0 + dz * 0.48)),
        ]
    elif profile == "ceiling_panel":
        # Thin overhead proxy: retain its perimeter and two panel joints in
        # each principal elevation.  A planar wire is sufficient here; a
        # curved or surfaced substitute would only create unwanted SVG edges.
        edges += [
            line((x0, cy, z0 + dz / 3), (x1, cy, z0 + dz / 3)),
            line((x0, cy, z0 + 2 * dz / 3), (x1, cy, z0 + 2 * dz / 3)),
            line((cx, y0, z0), (cx, y1, z0)),
            line((cx, y0, z1), (cx, y1, z1)),
        ]
    else:
        raise RuntimeError(f"unsupported lightweight profile: {profile}")

    # A centre-plane outline in each principal elevation keeps the object
    # legible even when its bounding edges overlap after projection.
    edges += rectangle_xz(x0, x1, cy, z0, z1)
    edges += rectangle_yz(cx, y0, y1, z0, z1)
    return edges


def representation_for_context(product, context):
    if not product.Representation:
        return None
    matches = [
        representation
        for representation in product.Representation.Representations
        if representation.ContextOfItems == context
    ]
    if len(matches) > 1:
        raise RuntimeError(f"{product.GlobalId}: duplicate ELEVATION_VIEW bodies")
    return matches[0] if matches else None


def add_representation(model, product, context, profile: str):
    existing = representation_for_context(product, context)
    bounds = local_bbox_mm(product)
    if existing:
        if existing.RepresentationType != "Curve3D" or any(
            not item.is_a("IfcPolyline") for item in existing.Items
        ):
            raise RuntimeError(
                f"{product.GlobalId}: existing ELEVATION_VIEW is not controlled Curve3D"
            )
        return existing, bounds, len(existing.Items), "existing"
    edges = profile_edges(bounds, profile)
    items = [
        model.createIfcPolyline(
            [model.createIfcCartesianPoint(start), model.createIfcCartesianPoint(end)]
        )
        for start, end in edges
    ]
    representation = model.createIfcShapeRepresentation(
        context, "Body", "Curve3D", items
    )
    definition = product.Representation
    definition.Representations = tuple(definition.Representations) + (representation,)
    return representation, bounds, len(edges), "created"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--allow-source-hash", default=EXPECTED_SOURCE_SHA256)
    args = parser.parse_args()
    source = args.input.resolve()
    output = args.output.resolve()
    if source == output:
        raise RuntimeError("candidate output must not overwrite source IFC")
    source_hash = sha256(source)
    if source_hash != args.allow_source_hash:
        raise RuntimeError(
            f"source IFC hash mismatch: expected {args.allow_source_hash}, got {source_hash}"
        )

    model = ifcopenshell.open(source)
    context = elevation_context(model)
    records = []
    for global_id, (expected_type, profile) in TARGETS.items():
        product = model.by_guid(global_id)
        actual_type = next(
            (relation.RelatingType.Name for relation in product.IsTypedBy), None
        )
        if actual_type != expected_type:
            raise RuntimeError(
                f"{global_id}: expected type {expected_type}, got {actual_type}"
            )
        representation, bounds, edge_count, action = add_representation(
            model, product, context, profile
        )
        records.append(
            {
                "global_id": global_id,
                "ifc_class": product.is_a(),
                "type_name": actual_type,
                "profile": profile,
                "representation_id": representation.id(),
                "representation_type": representation.RepresentationType,
                "context": [context.ContextType, context.ContextIdentifier, context.TargetView],
                "local_bbox_mm": [round(value, 6) for value in bounds],
                "edge_count": edge_count,
                "surface_count": 0,
                "action": action,
                "source": "model-derived controlled outline; official manufacturer CAD may supersede after evidence review",
            }
        )

    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_name(output.name + ".next")
    model.write(temporary)
    verified = ifcopenshell.open(temporary)
    verified_context = ifcopenshell.util.representation.get_context(
        verified, "Model", "Body", "ELEVATION_VIEW"
    )
    if not verified_context:
        raise RuntimeError("candidate reload lost ELEVATION_VIEW context")
    for global_id in TARGETS:
        representation = representation_for_context(
            verified.by_guid(global_id), verified_context
        )
        if not representation or representation.RepresentationType != "Curve3D":
            raise RuntimeError(f"{global_id}: candidate reload lost lightweight body")
        if any(not item.is_a("IfcPolyline") for item in representation.Items):
            raise RuntimeError(f"{global_id}: non-polyline item entered lightweight body")
    os.replace(temporary, output)

    report = {
        "source_ifc": str(source),
        "source_ifc_sha256": source_hash,
        "candidate_ifc": str(output),
        "candidate_ifc_sha256": sha256(output),
        "target_count": len(records),
        "total_edge_count": sum(record["edge_count"] for record in records),
        "surface_count": 0,
        "records": records,
        "pass": len(records) == 12 and all(record["edge_count"] < 40 for record in records),
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps({key: value for key, value in report.items() if key != "records"}, ensure_ascii=False))
    if not report["pass"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
