#!/usr/bin/env python3
"""Mechanically select review-only references from four acquired family DWGs.

This extractor never writes IFC.  It converts each immutable DWG to a temporary
DXF, selects the documented configuration/view regions, and stores only review
reference paths.  Native 2D paths are never scaled.  Where a source has no
dedicated 2D view, a dashed audit envelope is emitted and labelled as such.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import subprocess
import tempfile
import zipfile
from datetime import datetime, timezone
from pathlib import Path

import ezdxf
from ezdxf import disassemble


ROOT = Path(__file__).resolve().parents[2]
FORMAL_IFC = ROOT / "2504 GBTB Yanlord Zhuhai.ifc"
FORMAL_SHA256 = "7a50b87e8f48a7c2bfbaa9f1dabdfca7155fa684b2325c66aa1ae15c8ab0a25c"

PRODUCTS = {
    "bed01": {
        "manufacturer": "Baxter",
        "family": "Casablanca",
        "configuration": "180x200 / official native 2D group labelled 2100 x 2520 x 900 mm overall",
        "dwg": "output/review/highpoly-types/bed01/official-source/2d-3d-download/CASABLANCA/2D/Casablanca_Letto.dwg",
        "sha256": "ab697db9c27448c62a9a77537d7cc4b6286577c335328113d703184be28d9c4a",
        "url": "https://dam.baxter.it/asset/76c053fa-2ac5-4a60-85e3-d9df4ceb9cf0/Baxter_Casablanca_Bed_2D_3D.zip",
        "reference_role": "official_native_2d_180x200_configuration_reference_not_exact_project_body",
        "layers": ["_ARREDO", "_MATERASSO", "_PIEDINI", "_CUCITURE"],
        "view_boxes": {
            "plan": [-300.0, 5000.0, 1900.0, 7800.0],
            "front": [-300.0, 8300.0, 1900.0, 9400.0],
            "side": [2500.0, 8300.0, 5300.0, 9400.0],
        },
        "view_reflect_x": {"side": True},
        "archive": "output/review/highpoly-types/bed01/official-source/2d-3d-download/Baxter_Casablanca_Bed_2D_3D.zip",
        "archive_sha256": "acf67fddcfaa5e5bb50eb90f78a5eac262bfa318a767a938fc29f73931bb5970",
        "archive_member": "CASABLANCA/2D/Casablanca_Letto.dwg",
        "byte_identical_archives": [
            "/Users/jiaxinchen/Downloads/Baxter_Casablanca_Bed_2D_3D.zip",
            "/Users/jiaxinchen/Downloads/Baxter_Casablanca_Bed_2D_3D (1).zip",
            "/Users/jiaxinchen/Downloads/Baxter_Casablanca_Bed_2D_3D (2).zip",
        ],
        "available_layers": ["0", "Defpoints", "_QUOTE", "_MATERASSO", "_PIEDINI", "_CUCITURE", "_ARREDO"],
        "excluded_layers": {
            "_QUOTE": "dimension and annotation evidence retained in the AutoCAD screenshots, excluded from approval linework",
            "0": "six family-size labels retained in the AutoCAD screenshots; this file has no geometric entities on layer 0",
            "Defpoints": "non-printing definition layer with no selected product geometry",
        },
        "autocad_screenshots": [
            "output/review/highpoly-types/bed01/official-source/autocad-casablanca-2d-full-modelspace.png",
            "output/review/highpoly-types/bed01/official-source/autocad-casablanca-180x200-three-view.png",
        ],
        "note": "The authenticated Baxter 2D/3D archive contains the independent native 2D file Casablanca_Letto.dwg. Only the labelled 180x200 Plan/Front/Side group is selected from the four product-geometry layers. _QUOTE dimensions and layer-0 family labels remain visible in the AutoCAD evidence but are excluded from approval linework. No source path is scaled or fitted.",
    },
    "bed02": {
        "manufacturer": "Baxter",
        "family": "Viktor",
        "configuration": "160x200 / official 1720 x 2340 x 1060 mm family variant",
        "dwg": "output/review/highpoly-types/bed02/official-source/official-download/Viktor_Letto.dwg",
        "sha256": "3f0676a004f3e744779093188d75d58238153332d9d1810f28821cafc79a9bbd",
        "url": "https://dam.baxter.it/asset/6c3331bf-d64d-4bf6-83f9-40c1ffc41344/Baxter_Viktor_Bed_2D_3D.zip",
        "reference_role": "official_family_identity_reference_not_exact_project_configuration",
        "layers": ["ARREDO"],
        "view_boxes": {
            "plan": [98300.0, -84500.0, 100100.0, -81800.0],
            "front": [98300.0, -81600.0, 100100.0, -80000.0],
            "side": [100700.0, -81600.0, 103200.0, -80000.0],
        },
        "view_reflect_x": {"side": True},
        "note": "The authenticated DWG contains six family variants. Only the labelled 160x200 group is selected. The project Body is a width-reduced, non-uniformly changed model, so these paths are family/identity reference only and are not fitted to it.",
    },
    "sis04": {
        "manufacturer": "Molteni&C",
        "family": "Sistema 7 Wall Unit",
        "configuration": "1962 x 370 x 722 mm / project 4 Doors",
        "dwg": "output/review/highpoly-types/sis04/official-source/official-download/0000_2D_Sistema-7_Wall-Units_Kitchens.dwg",
        "sha256": "77f240f65f1ab68a168b363b18f8294122ec913259f93745a01d19d534678263",
        "url": "https://res.cloudinary.com/molteni/raw/upload/v1752749168/0000_2D_Sistema-7_Wall-Units_Kitchens.dwg",
        "reference_role": "exact_official_wall_unit_configuration_reference",
        "layers": ["DADA_FURNITURES", "_DADA FURNITURES", "_DADA TEXTURES"],
        "contained_view_boxes": {
            "plan": [-7735.945, 2835.027, -5753.943, 3205.029],
        },
        "envelopes": {
            "front": [1962.0, 722.0],
            "side": [370.0, 722.0],
        },
        "note": "The native Plan paths are selected from the exact 1962 x 370 mm group. The DWG does not provide a dedicated closed Front/Side path set for the selected project pose; those SVG views therefore use dashed dimension envelopes from the same DWG and say so explicitly. Catalogue depth 331 mm is internal; 370 mm is the DWG overall depth.",
    },
    "hima01": {
        "manufacturer": "Poliform",
        "family": "HIMA",
        "configuration": "PVA11 / 3 elements / 2330 x 115 x 1000 mm",
        "dwg": "output/review/highpoly-types/hima01/official-source/official-download/Poliform-HIMA-screen.dwg",
        "sha256": "f0e1b980aa102f3e06fe602f9e7a9b40db45d76e45200c08e4b8b4a0eec1d523",
        "url": "https://s3.poliform.it/2025/09/Poliform-HIMA-screen.dwg",
        "reference_role": "official_pva11_family_and_height_reference_not_exact_project_folded_pose",
        "layers": ["drawing", "ropes"],
        "view_boxes": {
            "plan": [-12700.0, -1450.0, -10000.0, -1200.0],
            "front": [-12700.0, -800.0, -10000.0, 300.0],
        },
        "envelopes": {"side": [115.0, 1000.0]},
        "note": "PVA11 is the official three-element 1000 mm variant. The project Body is folded to a different Plan envelope; no scaling or pose deformation is applied. Side is a dashed DWG dimension envelope because no dedicated PVA11 side path set is published in this file.",
    },
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def rectangle(width: float, height: float) -> list[list[list[float]]]:
    return [[[0.0, 0.0], [width, 0.0], [width, height], [0.0, height], [0.0, 0.0]]]


def bounds(path: list[list[float]]) -> tuple[float, float, float, float]:
    xs = [point[0] for point in path]
    ys = [point[1] for point in path]
    return min(xs), min(ys), max(xs), max(ys)


def normalise(paths: list[list[list[float]]]) -> list[list[list[float]]]:
    minimum_x = min(point[0] for path in paths for point in path)
    minimum_y = min(point[1] for path in paths for point in path)
    return [
        [[round(x - minimum_x, 6), round(y - minimum_y, 6)] for x, y in path]
        for path in paths
    ]


def path_inventory(doc: ezdxf.document.Drawing, layers: list[str]) -> list[dict]:
    entities = [entity for entity in doc.modelspace() if entity.dxf.layer in layers]
    records = []
    for primitive in disassemble.to_primitives(disassemble.recursive_decompose(entities)):
        try:
            vertices = [
                vertex for vertex in primitive.vertices()
                if math.isfinite(vertex.x) and math.isfinite(vertex.y)
            ]
        except Exception:
            continue
        if len(vertices) < 2:
            continue
        path = [[float(vertex.x), float(vertex.y)] for vertex in vertices]
        minimum_x, minimum_y, maximum_x, maximum_y = bounds(path)
        records.append({
            "path": path,
            "bounds": [minimum_x, minimum_y, maximum_x, maximum_y],
            "centre": [(minimum_x + maximum_x) / 2.0, (minimum_y + maximum_y) / 2.0],
        })
    return records


def select_by_centre(records: list[dict], box: list[float]) -> list[list[list[float]]]:
    minimum_x, minimum_y, maximum_x, maximum_y = box
    return [
        record["path"] for record in records
        if minimum_x <= record["centre"][0] <= maximum_x
        and minimum_y <= record["centre"][1] <= maximum_y
    ]


def select_contained(records: list[dict], box: list[float]) -> list[list[list[float]]]:
    minimum_x, minimum_y, maximum_x, maximum_y = box
    tolerance = 0.01
    return [
        record["path"] for record in records
        if record["bounds"][0] >= minimum_x - tolerance
        and record["bounds"][1] >= minimum_y - tolerance
        and record["bounds"][2] <= maximum_x + tolerance
        and record["bounds"][3] <= maximum_y + tolerance
    ]


def path_bounds(paths: list[list[list[float]]]) -> dict:
    points = [point for path in paths for point in path]
    minimum = [min(point[axis] for point in points) for axis in range(2)]
    maximum = [max(point[axis] for point in points) for axis in range(2)]
    return {
        "minimum": [round(value, 6) for value in minimum],
        "maximum": [round(value, 6) for value in maximum],
        "size": [round(maximum[axis] - minimum[axis], 6) for axis in range(2)],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--slug", choices=PRODUCTS)
    args = parser.parse_args()
    if sha256(FORMAL_IFC) != FORMAL_SHA256:
        raise RuntimeError("formal IFC hash mismatch before DWG reference extraction")
    selected_products = PRODUCTS.items() if args.slug is None else [(args.slug, PRODUCTS[args.slug])]
    for slug, config in selected_products:
        source = ROOT / config["dwg"]
        if sha256(source) != config["sha256"]:
            raise RuntimeError(f"{slug} official DWG hash mismatch")
        archive_inventory_path = None
        archive_inventory_sha256 = None
        archive_member_sha256 = None
        byte_identical_archives = []
        if config.get("archive"):
            archive = ROOT / config["archive"]
            if sha256(archive) != config["archive_sha256"]:
                raise RuntimeError(f"{slug} official archive hash mismatch")
            with zipfile.ZipFile(archive) as package:
                members = [
                    {
                        "name": info.filename,
                        "is_directory": info.is_dir(),
                        "uncompressed_bytes": info.file_size,
                        "compressed_bytes": info.compress_size,
                        "crc32": f"{info.CRC:08x}",
                    }
                    for info in package.infolist()
                ]
                archive_member_sha256 = hashlib.sha256(package.read(config["archive_member"])).hexdigest()
            if archive_member_sha256 != config["sha256"]:
                raise RuntimeError(f"{slug} archive member and extracted DWG differ")
            inventory = {
                "schema_version": 1,
                "archive": config["archive"],
                "archive_sha256": config["archive_sha256"],
                "official_download_url": config["url"],
                "member_count": len(members),
                "directory_count": sum(member["is_directory"] for member in members),
                "dwg_members": [member["name"] for member in members if member["name"].lower().endswith(".dwg")],
                "members": members,
                "pass": True,
            }
            archive_inventory_path = archive.parent / "archive-inventory.json"
            archive_inventory_path.write_text(json.dumps(inventory, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
            archive_inventory_sha256 = sha256(archive_inventory_path)
            for copy_path_text in config.get("byte_identical_archives", []):
                copy_path = Path(copy_path_text)
                copy_hash = sha256(copy_path)
                if copy_hash != config["archive_sha256"]:
                    raise RuntimeError(f"{slug} local archive copy differs: {copy_path}")
                byte_identical_archives.append({
                    "absolute_path": str(copy_path),
                    "sha256": copy_hash,
                    "byte_identical_to_preserved_archive": True,
                })
        views = {}
        with tempfile.TemporaryDirectory(prefix=f"{slug}-official-dwg-") as temporary:
            dxf = Path(temporary) / f"{slug}.dxf"
            conversion = subprocess.run(
                ["dwg2dxf", "-y", "-o", str(dxf), str(source)],
                check=False,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
            )
            if conversion.returncode != 0 or not dxf.is_file():
                raise RuntimeError(
                    f"{slug} DWG conversion failed: {conversion.stderr[-2000:]}"
                )
            doc = ezdxf.readfile(dxf)
            records = path_inventory(doc, config["layers"])
            for view, box in config.get("view_boxes", {}).items():
                paths = normalise(select_by_centre(records, box))
                if not paths:
                    raise RuntimeError(f"{slug} {view} selection is empty")
                views[view] = {
                    "geometry_kind": "native_dwg_2d_linework",
                    "paths_mm": paths,
                    "path_count": len(paths),
                    "bounds_mm": path_bounds(paths),
                    "selection_box_source_coordinates": box,
                    "selection_rule": "path_bbox_centre_inside_documented_configuration_view_box",
                    "blue_stroke_style": "solid",
                    "dedicated_native_2d_view_present": True,
                    "reflect_x_for_project_view_direction": config.get("view_reflect_x", {}).get(view, False),
                    "source_scaled": False,
                }
            for view, box in config.get("contained_view_boxes", {}).items():
                paths = normalise(select_contained(records, box))
                if not paths:
                    raise RuntimeError(f"{slug} {view} contained selection is empty")
                views[view] = {
                    "geometry_kind": "native_dwg_2d_linework",
                    "paths_mm": paths,
                    "path_count": len(paths),
                    "bounds_mm": path_bounds(paths),
                    "selection_box_source_coordinates": box,
                    "selection_rule": "path_bbox_fully_contained_in_exact_configuration_envelope",
                    "blue_stroke_style": "solid",
                    "dedicated_native_2d_view_present": True,
                    "reflect_x_for_project_view_direction": False,
                    "source_scaled": False,
                }
            for view, size in config.get("envelopes", {}).items():
                paths = rectangle(*size)
                views[view] = {
                    "geometry_kind": "native_dwg_dimension_envelope_not_dedicated_2d_view",
                    "paths_mm": paths,
                    "path_count": 1,
                    "bounds_mm": path_bounds(paths),
                    "blue_stroke_style": "dashed",
                    "dedicated_native_2d_view_present": False,
                    "source_scaled": False,
                }
        if set(views) != {"plan", "front", "side"}:
            raise RuntimeError(f"{slug} does not have a complete review view record")
        output = ROOT / "output/review/highpoly-types" / slug / "official-dwg-review-reference.json"
        record = {
            "schema_version": 1,
            "generated_at": datetime.now(timezone.utc).isoformat(),
            "generator": "pipeline/scripts/extract_official_dwg_family_review_references.py",
            "manufacturer": config["manufacturer"],
            "family": config["family"],
            "configuration": config["configuration"],
            "reference_role": config["reference_role"],
            "source_kind": "native_dwg_review_reference",
            "source_dwg": config["dwg"],
            "source_dwg_absolute_path": str(source.resolve()),
            "source_dwg_sha256": config["sha256"],
            "official_download_url": config["url"],
            "third_party_cad_used": False,
            "source_geometry_scaled": False,
            "source_geometry_anisotropically_fitted": False,
            "primary_candidate_replaced": False,
            "note": config["note"],
            "source_archive": config.get("archive"),
            "source_archive_sha256": config.get("archive_sha256"),
            "source_archive_member": config.get("archive_member"),
            "source_archive_member_sha256": archive_member_sha256,
            "source_archive_inventory": (str(archive_inventory_path.relative_to(ROOT)) if archive_inventory_path else None),
            "source_archive_inventory_sha256": archive_inventory_sha256,
            "local_byte_identical_archive_copies": byte_identical_archives,
            "source_layers_available": config.get("available_layers"),
            "selected_geometry_layers": config.get("layers"),
            "excluded_layers": config.get("excluded_layers"),
            "autocad_source_screenshots": config.get("autocad_screenshots", []),
            "superseded_error_candidates": ([{
                "status": "superseded_error_candidate_excluded_from_approval",
                "source_dwg": "output/review/highpoly-types/bed01/official-source/bim-download/CASABLANCA_BED BASE 180 X 200 H90.dwg",
                "source_dwg_sha256": "1d3fb7931d2252f301498acdea30a145303f0a2aadbeb8e1638fdf7b800642cb",
                "geometry_kind": "native_dwg_3d_solid_projected_envelope",
                "reason": "ACIS 3DSOLID projection envelopes are not the independently published manufacturer 2D Plan/Front/Side linework and must never enter approval.",
            }] if slug == "bed01" else []),
            "views": views,
            "formal_ifc_sha256": FORMAL_SHA256,
            "formal_ifc_bytes_unchanged": sha256(FORMAL_IFC) == FORMAL_SHA256,
            "pass": True,
        }
        output.write_text(json.dumps(record, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
        print(output.relative_to(ROOT))
    if sha256(FORMAL_IFC) != FORMAL_SHA256:
        raise RuntimeError("formal IFC changed during DWG reference extraction")


if __name__ == "__main__":
    main()
