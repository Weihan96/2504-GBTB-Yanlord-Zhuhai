#!/usr/bin/env python3
"""Repair the BED02 Viktor native-DWG review reference and SVG package.

The official Viktor family DWG is valid.  The superseded review extraction was
incomplete because it selected only ARREDO.  This BED02-only repair selects the
160x200 group from ARREDO plus geometric entities on layer 0, while excluding
the _QUOTE dimension layer and layer-0 TEXT/MTEXT labels.  It never writes IFC.

Run with:
  uv run --with 'ezdxf>=1.4,<2' --with 'ifcopenshell==0.8.4' \
    python pipeline/scripts/repair_bed02_viktor_official_dwg_review.py
"""

from __future__ import annotations

import hashlib
import json
import math
import re
import subprocess
import sys
import tempfile
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import ezdxf
from ezdxf import disassemble


ROOT = Path(__file__).resolve().parents[2]
FORMAL_IFC = ROOT / "2504 GBTB Yanlord Zhuhai.ifc"
FORMAL_SHA256 = "7a50b87e8f48a7c2bfbaa9f1dabdfca7155fa684b2325c66aa1ae15c8ab0a25c"
FOLDER = ROOT / "output/review/highpoly-types/bed02"
SOURCE_DWG = FOLDER / "official-source/official-download/Viktor_Letto.dwg"
SOURCE_DWG_RELATIVE = "output/review/highpoly-types/bed02/official-source/official-download/Viktor_Letto.dwg"
SOURCE_DWG_SHA256 = "3f0676a004f3e744779093188d75d58238153332d9d1810f28821cafc79a9bbd"
OFFICIAL_DOWNLOAD_URL = "https://dam.baxter.it/asset/6c3331bf-d64d-4bf6-83f9-40c1ffc41344/Baxter_Viktor_Bed_2D_3D.zip"
REFERENCE_PATH = FOLDER / "official-dwg-review-reference.json"
PACKAGE_PATH = FOLDER / "official-dwg-review-manifest.json"
SOURCE_RECORD_PATH = FOLDER / "official-source/source-access-record.json"
PRIMARY_MANIFEST_PATH = FOLDER / "manifest.json"
PRIMARY_INDEX_PATH = FOLDER / "index.html"
OFFICIAL_INDEX_PATH = FOLDER / "official-dwg-review-index.html"
SHARED_GENERATOR = ROOT / "pipeline/scripts/generate_official_dwg_family_review_packages.py"

SELECTED_GEOMETRY_LAYERS = ["ARREDO", "0"]
ALLOWED_ENTITY_TYPES = {
    "ARC",
    "CIRCLE",
    "ELLIPSE",
    "LINE",
    "LWPOLYLINE",
    "POLYLINE",
    "SPLINE",
}
VIEW_BOXES = {
    "plan": [98300.0, -84500.0, 100100.0, -81800.0],
    "front": [98300.0, -81600.0, 100100.0, -80000.0],
    "side": [100700.0, -81600.0, 103200.0, -80000.0],
}
EXPECTED = {
    "plan": {"path_count": 51, "path_count_by_layer": {"ARREDO": 28, "0": 23}, "size": [1723.392281, 2340.708974]},
    "front": {"path_count": 99, "path_count_by_layer": {"ARREDO": 91, "0": 8}, "size": [1733.666161, 1060.0]},
    "side": {"path_count": 84, "path_count_by_layer": {"ARREDO": 77, "0": 7}, "size": [2351.640392, 1065.443243]},
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, value: dict) -> None:
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def relative(path: Path) -> str:
    return path.resolve().relative_to(ROOT).as_posix()


def path_bounds(paths: list[list[list[float]]]) -> dict:
    points = [point for path in paths for point in path]
    minimum = [min(point[axis] for point in points) for axis in range(2)]
    maximum = [max(point[axis] for point in points) for axis in range(2)]
    return {
        "minimum": [round(value, 6) for value in minimum],
        "maximum": [round(value, 6) for value in maximum],
        "size": [round(maximum[axis] - minimum[axis], 6) for axis in range(2)],
    }


def normalise(paths: list[list[list[float]]]) -> list[list[list[float]]]:
    minimum_x = min(point[0] for path in paths for point in path)
    minimum_y = min(point[1] for path in paths for point in path)
    return [
        [[round(x - minimum_x, 6), round(y - minimum_y, 6)] for x, y in path]
        for path in paths
    ]


def inventory(doc: ezdxf.document.Drawing) -> tuple[list[dict], dict]:
    records: list[dict] = []
    excluded_entity_counts: Counter[str] = Counter()
    selected_entity_counts: Counter[str] = Counter()
    for entity in doc.modelspace():
        layer = entity.dxf.layer
        entity_type = entity.dxftype()
        if layer not in SELECTED_GEOMETRY_LAYERS:
            continue
        if entity_type not in ALLOWED_ENTITY_TYPES:
            excluded_entity_counts[f"{layer}:{entity_type}"] += 1
            continue
        selected_entity_counts[f"{layer}:{entity_type}"] += 1
        for primitive in disassemble.to_primitives(disassemble.recursive_decompose([entity])):
            try:
                vertices = [
                    vertex
                    for vertex in primitive.vertices()
                    if math.isfinite(vertex.x) and math.isfinite(vertex.y)
                ]
            except Exception:
                continue
            if len(vertices) < 2:
                continue
            path = [[float(vertex.x), float(vertex.y)] for vertex in vertices]
            bounds = path_bounds([path])
            minimum_x, minimum_y = bounds["minimum"]
            maximum_x, maximum_y = bounds["maximum"]
            records.append(
                {
                    "path": path,
                    "source_layer": layer,
                    "source_entity_type": entity_type,
                    "source_handle": entity.dxf.handle,
                    "bounds": [minimum_x, minimum_y, maximum_x, maximum_y],
                    "centre": [(minimum_x + maximum_x) / 2.0, (minimum_y + maximum_y) / 2.0],
                }
            )
    return records, {
        "selected_entity_counts": dict(sorted(selected_entity_counts.items())),
        "excluded_entity_counts": dict(sorted(excluded_entity_counts.items())),
    }


def select_view(records: list[dict], box: list[float]) -> tuple[list[list[list[float]]], dict[str, int]]:
    minimum_x, minimum_y, maximum_x, maximum_y = box
    selected = [
        record
        for record in records
        if minimum_x <= record["centre"][0] <= maximum_x
        and minimum_y <= record["centre"][1] <= maximum_y
    ]
    paths = normalise([record["path"] for record in selected])
    count_by_layer = Counter(record["source_layer"] for record in selected)
    return paths, {layer: count_by_layer[layer] for layer in SELECTED_GEOMETRY_LAYERS}


def close_enough(actual: list[float], expected: list[float], tolerance: float = 0.001) -> bool:
    return all(abs(actual[index] - expected[index]) <= tolerance for index in range(len(expected)))


def build_superseded_record(existing_reference: dict | None, existing_package: dict | None) -> list[dict]:
    if existing_reference and existing_reference.get("superseded_extractions"):
        return existing_reference["superseded_extractions"]
    if not existing_reference:
        return []
    old_views = {}
    package_views = {
        record["view"]: record for record in (existing_package or {}).get("views", [])
    }
    for view in ("plan", "front", "side"):
        source_view = existing_reference["views"][view]
        package_view = package_views.get(view, {})
        old_views[view] = {
            "path_count": source_view["path_count"],
            "bounds_mm": source_view["bounds_mm"],
            "svg_sha256": package_view.get("svg_sha256"),
        }
    return [
        {
            "status": "superseded_incomplete_extraction_excluded_from_approval",
            "source_dwg_valid": True,
            "source_dwg_sha256": SOURCE_DWG_SHA256,
            "generator": existing_reference.get("generator"),
            "reference_sha256_before_repair": sha256(REFERENCE_PATH),
            "selected_geometry_layers": ["ARREDO"],
            "excluded_geometry_layer_in_error": "0",
            "reason": "The DWG is correct, but the ARREDO-only extraction omitted required native bed geometry on layer 0. The superseded SVGs are extraction errors, not manufacturer-DWG errors.",
            "views": old_views,
        }
    ]


def update_source_record(reference: dict, package: dict) -> None:
    record = load_json(SOURCE_RECORD_PATH)
    record["official_native_dwg_review_reference"] = {
        "status": "corrected_native_dwg_reference_pending_visual_approval",
        "source_dwg": SOURCE_DWG_RELATIVE,
        "source_dwg_sha256": SOURCE_DWG_SHA256,
        "configuration": reference["configuration"],
        "selected_geometry_layers": SELECTED_GEOMETRY_LAYERS,
        "allowed_geometry_entity_types": sorted(ALLOWED_ENTITY_TYPES),
        "excluded_layers": reference["excluded_layers"],
        "excluded_entity_types": reference["excluded_entity_types"],
        "views": {
            view: {
                "path_count": reference["views"][view]["path_count"],
                "path_count_by_source_layer": reference["views"][view]["path_count_by_source_layer"],
                "bounds_mm": reference["views"][view]["bounds_mm"],
            }
            for view in ("plan", "front", "side")
        },
        "review_package_manifest": relative(PACKAGE_PATH),
        "review_package_manifest_sha256": sha256(PACKAGE_PATH),
        "primary_candidate_replaced": False,
        "formal_ifc_write_performed": False,
        "superseded_extractions": reference["superseded_extractions"],
    }
    record["scope"] = (
        "Manufacturer-authenticated six-variant family CAD. The corrected 160x200 native-DWG "
        "Plan/Front/Side is rendered as a blue 1:1 review reference; it is not the primary "
        "project drawing geometry and is not an exact project Body configuration."
    )
    write_json(SOURCE_RECORD_PATH, record)


def update_primary_manifest(reference: dict, package: dict) -> None:
    manifest = load_json(PRIMARY_MANIFEST_PATH)
    manifest["blue_line_present"] = True
    manifest["blue_line_role"] = "official_native_dwg_review_reference_only_not_primary_candidate"
    manifest["official_cad_used"] = False
    manifest["official_identity_evidence_only"] = True
    manifest["drawing_source"]["official_product_cad_status"] = (
        "manufacturer_authenticated_family_cad_acquired_and_rendered_as_review_reference_not_primary"
    )
    manifest["official_reference"]["official_cad_used"] = False
    manifest["official_reference"]["official_review_reference_rendered"] = True
    manifest["official_reference"]["scope"] = (
        "Manufacturer family identity, nominal dimensions, and native 160x200 review linework; "
        "not an exact project configuration, not a project shop drawing, and not the primary candidate."
    )
    manifest["official_dwg_review_supplement"] = {
        "status": "pending_visual_approval",
        "package": relative(PACKAGE_PATH),
        "package_sha256": sha256(PACKAGE_PATH),
        "reference": relative(REFERENCE_PATH),
        "reference_sha256": sha256(REFERENCE_PATH),
        "selected_geometry_layers": SELECTED_GEOMETRY_LAYERS,
        "excluded_layers": reference["excluded_layers"],
        "views": {
            record["view"]: {
                "svg": record["svg"],
                "svg_sha256": record["svg_sha256"],
                "path_count": record["official_reference_path_count"],
                "bounds_mm": reference["views"][record["view"]]["bounds_mm"],
                "uniform_scale": record["alignment"]["uniform_scale"],
                "view_direction_reflection_x": record["alignment"]["view_direction_reflection_x"],
            }
            for record in package["views"]
        },
        "superseded_extractions": reference["superseded_extractions"],
        "derived_ifc_write_performed": False,
        "formal_authoritative_ifc_write_performed": False,
    }
    source_record_hash = sha256(SOURCE_RECORD_PATH)
    manifest["official_source_access_record_sha256"] = source_record_hash
    manifest["drawing_source"]["official_source_access_record_sha256"] = source_record_hash
    manifest["formal_ifc_sha256"] = FORMAL_SHA256
    manifest["formal_ifc_bytes_unchanged"] = sha256(FORMAL_IFC) == FORMAL_SHA256
    write_json(PRIMARY_MANIFEST_PATH, manifest)


def update_indexes() -> None:
    official = OFFICIAL_INDEX_PATH.read_text(encoding="utf-8")
    marker = "BED02 extraction correction"
    correction = (
        f'<p><strong>{marker}:</strong> the official DWG is valid. The superseded ARREDO-only SVG '
        "extraction omitted required layer-0 bed geometry. This package selects geometric entities "
        "from ARREDO + 0, excludes _QUOTE and TEXT/MTEXT, and keeps native scale 1:1.</p>"
    )
    if marker not in official:
        official = official.replace("<nav>", f"{correction}<nav>", 1)
    OFFICIAL_INDEX_PATH.write_text(official, encoding="utf-8")

    primary = PRIMARY_INDEX_PATH.read_text(encoding="utf-8")
    paragraph = (
        "<p>Black line = <strong>基于原始高模几何生成的简化图纸表达</strong>. "
        "The supplemental blue line is the corrected 1:1 native Viktor 160×200 DWG review reference "
        "from ARREDO + layer 0 geometry; _QUOTE dimensions and TEXT/MTEXT are excluded. "
        "It remains family/reference evidence only and has not replaced the black primary candidate.</p>"
    )
    primary = re.sub(
        r"(<h1>Baxter Viktor / project BED02</h1>)<p>.*?</p><nav>",
        rf"\1{paragraph}<nav>",
        primary,
        count=1,
        flags=re.DOTALL,
    )
    PRIMARY_INDEX_PATH.write_text(primary, encoding="utf-8")


def main() -> None:
    if sha256(FORMAL_IFC) != FORMAL_SHA256:
        raise RuntimeError("formal IFC hash mismatch before BED02 repair")
    if sha256(SOURCE_DWG) != SOURCE_DWG_SHA256:
        raise RuntimeError("official Viktor DWG hash mismatch")

    existing_reference = load_json(REFERENCE_PATH) if REFERENCE_PATH.is_file() else None
    existing_package = load_json(PACKAGE_PATH) if PACKAGE_PATH.is_file() else None
    superseded = build_superseded_record(existing_reference, existing_package)

    with tempfile.TemporaryDirectory(prefix="bed02-viktor-official-dwg-") as temporary:
        dxf = Path(temporary) / "Viktor_Letto.dxf"
        conversion = subprocess.run(
            ["dwg2dxf", "-y", "-o", str(dxf), str(SOURCE_DWG)],
            check=False,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        if conversion.returncode != 0 or not dxf.is_file():
            raise RuntimeError(f"DWG conversion failed: {conversion.stderr[-2000:]}")
        document = ezdxf.readfile(dxf)
        records, inventory_audit = inventory(document)
        available_layers = sorted(layer.dxf.name for layer in document.layers)
        views = {}
        for view, box in VIEW_BOXES.items():
            paths, count_by_layer = select_view(records, box)
            bounds = path_bounds(paths)
            expected = EXPECTED[view]
            if len(paths) != expected["path_count"]:
                raise RuntimeError(f"{view} path count {len(paths)} != {expected['path_count']}")
            if count_by_layer != expected["path_count_by_layer"]:
                raise RuntimeError(f"{view} layer counts {count_by_layer} != {expected['path_count_by_layer']}")
            if not close_enough(bounds["size"], expected["size"]):
                raise RuntimeError(f"{view} bounds {bounds['size']} != {expected['size']}")
            views[view] = {
                "geometry_kind": "native_dwg_2d_linework",
                "paths_mm": paths,
                "path_count": len(paths),
                "path_count_by_source_layer": count_by_layer,
                "bounds_mm": bounds,
                "selection_box_source_coordinates": box,
                "selection_rule": "allowed_geometry_entity_path_bbox_centre_inside_160x200_view_box",
                "selected_geometry_layers": SELECTED_GEOMETRY_LAYERS,
                "excluded_dimension_layer": "_QUOTE",
                "excluded_entity_types": ["TEXT", "MTEXT"],
                "blue_stroke_style": "solid",
                "dedicated_native_2d_view_present": True,
                "reflect_x_for_project_view_direction": view == "side",
                "source_scaled": False,
            }

    screenshots = []
    for screenshot in (
        Path.home() / "Desktop/Viktor-official-DWG-full-model-space-2026-08-25.png",
        Path.home() / "Desktop/Viktor-official-DWG-160x200-views-2026-08-25.png",
    ):
        if screenshot.is_file():
            screenshots.append({"absolute_path": str(screenshot), "sha256": sha256(screenshot)})

    reference = {
        "schema_version": 2,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "generator": "pipeline/scripts/repair_bed02_viktor_official_dwg_review.py",
        "manufacturer": "Baxter",
        "family": "Viktor",
        "configuration": "160x200 / official dimension labels 1720 x 2340 x 1060 mm",
        "reference_role": "official_family_identity_and_native_160x200_linework_reference_not_exact_project_configuration",
        "source_kind": "native_dwg_review_reference",
        "source_dwg": SOURCE_DWG_RELATIVE,
        "source_dwg_absolute_path": str(SOURCE_DWG.resolve()),
        "source_dwg_sha256": SOURCE_DWG_SHA256,
        "official_download_url": OFFICIAL_DOWNLOAD_URL,
        "third_party_cad_used": False,
        "source_geometry_scaled": False,
        "source_geometry_anisotropically_fitted": False,
        "primary_candidate_replaced": False,
        "source_layers_available": available_layers,
        "selected_geometry_layers": SELECTED_GEOMETRY_LAYERS,
        "allowed_geometry_entity_types": sorted(ALLOWED_ENTITY_TYPES),
        "excluded_layers": {
            "_QUOTE": "Official red dimensions remain source evidence and are excluded from blue approval linework.",
            "Defpoints": "Non-printing definition layer; excluded.",
        },
        "excluded_entity_types": {
            "TEXT": "Family-size labels are source evidence, not product geometry.",
            "MTEXT": "The layer-0 160x200 title would otherwise become a false Front path; excluded.",
        },
        "inventory_audit": inventory_audit,
        "official_dimension_cross_check_mm": {
            "nominal_overall": [1720.0, 2340.0, 1060.0],
            "native_soft_linework_bounds": {
                view: views[view]["bounds_mm"]["size"] for view in ("plan", "front", "side")
            },
            "note": "Small width/height overhangs are native upholstered curves; no source path is scaled or deformed.",
        },
        "autocad_visual_verification": {
            "application": "AutoCAD 2024",
            "trusted_dwg_reported": True,
            "screenshots": screenshots,
        },
        "superseded_extractions": superseded,
        "note": (
            "The authenticated DWG contains six family variants. Only the labelled 160x200 group is selected. "
            "The corrected extraction includes native geometric entities from ARREDO and layer 0, excludes "
            "_QUOTE and text labels, and remains 1:1. The project Body is a width-reduced, non-uniformly changed "
            "model, so the blue linework remains family/reference evidence only."
        ),
        "views": views,
        "formal_ifc_sha256": FORMAL_SHA256,
        "formal_ifc_bytes_unchanged": True,
        "pass": True,
    }
    write_json(REFERENCE_PATH, reference)

    subprocess.run([sys.executable, str(SHARED_GENERATOR), "--slug", "bed02"], check=True, cwd=ROOT)
    update_indexes()

    package = load_json(PACKAGE_PATH)
    package["repair_generator"] = "pipeline/scripts/repair_bed02_viktor_official_dwg_review.py"
    package["selected_geometry_layers"] = SELECTED_GEOMETRY_LAYERS
    package["allowed_geometry_entity_types"] = sorted(ALLOWED_ENTITY_TYPES)
    package["excluded_layers"] = reference["excluded_layers"]
    package["excluded_entity_types"] = reference["excluded_entity_types"]
    package["official_dimension_cross_check_mm"] = reference["official_dimension_cross_check_mm"]
    package["superseded_extractions"] = superseded
    package["index_sha256"] = sha256(OFFICIAL_INDEX_PATH)
    package["formal_ifc_bytes_unchanged"] = sha256(FORMAL_IFC) == FORMAL_SHA256
    package["pass"] = True
    write_json(PACKAGE_PATH, package)

    update_source_record(reference, package)
    update_primary_manifest(reference, package)

    if sha256(FORMAL_IFC) != FORMAL_SHA256:
        raise RuntimeError("formal IFC changed during BED02 repair")
    print(relative(REFERENCE_PATH))
    print(relative(PACKAGE_PATH))
    for view in ("plan", "front", "side"):
        print(relative(FOLDER / f"official-dwg-{view}.svg"))


if __name__ == "__main__":
    main()
