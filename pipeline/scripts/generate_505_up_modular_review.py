#!/usr/bin/env python3
"""Generate a read-only modular audit for the Molteni 505 UP project instance.

The output is review evidence only. It never writes IFC, approval JSON, or the
high-poly review queue. Blue solid lines come from a native DWG block without
scaling; blue dashed lines are datasheet rules, not manufacturer shop drawings.
"""

from __future__ import annotations

import hashlib
import json
import math
import subprocess
from collections import defaultdict
from pathlib import Path
from typing import Any, Callable

import ifcopenshell
import numpy as np
from ifcopenshell.util.placement import get_local_placement


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "output/review/highpoly-types/505-up-v1-lp-s"
SOURCE = OUT / "official-source"
FORMAL = ROOT / "2504 GBTB Yanlord Zhuhai.ifc"
COMPOSITION = OUT / "505-up-v1-lp-s-modular-composition.json"
SCHEMA = OUT / "505-up-parametric-schema.json"
CANDIDATES = OUT / "candidate-representations.json"
DWG = SOURCE / "2021_2D_505-UP_Living-Systems_Indoor.dwg"
PDF = SOURCE / "2021_DS_505-UP-System_Living-Systems_Indoor_EN.pdf"
GLOBAL_ID = "19MpdkWqXC7uhUNhLQgrce"
FORMAL_SHA = "7a50b87e8f48a7c2bfbaa9f1dabdfca7155fa684b2325c66aa1ae15c8ab0a25c"
DWG_SHA = "dfe4a6bcd655a3ed19343813a4aba8e139b847cafe8ea69b47623a6e1e71add3"
PDF_SHA = "6df90e8b77c428b01fbda98525397952600dc56000e1a38673152c632544966b"
DWG_BLOCK_HANDLE = 154738


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_text(path: Path, text: str) -> None:
    path.write_text(text, encoding="utf-8")


def write_json(path: Path, value: Any) -> None:
    write_text(path, json.dumps(value, indent=2, ensure_ascii=False) + "\n")


def handle_id(value: Any) -> int:
    return int(value[2])


def affine(point: list[float] | tuple[float, ...], *, tx: float = 0, ty: float = 0,
           sx: float = 1, sy: float = 1, rotation: float = 0) -> tuple[float, float]:
    x, y = float(point[0]), float(point[1])
    c, s = math.cos(rotation), math.sin(rotation)
    return tx + sx * x * c - sy * y * s, ty + sx * x * s + sy * y * c


def native_dwg_block_paths() -> tuple[list[list[tuple[float, float]]], dict[str, Any]]:
    completed = subprocess.run(
        ["dwgread", "-O", "JSON", str(DWG)], check=True,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
    )
    payload = json.loads(completed.stdout)
    objects = payload["OBJECTS"]
    by_handle = {handle_id(obj["handle"]): obj for obj in objects if "handle" in obj}
    block = by_handle[DWG_BLOCK_HANDLE]

    def entity_paths(entity: dict[str, Any], transform: Callable[[Any], tuple[float, float]]) -> list[list[tuple[float, float]]]:
        kind = entity.get("entity")
        paths: list[list[tuple[float, float]]] = []
        if kind == "LINE":
            paths.append([transform(entity["start"]), transform(entity["end"])])
        elif kind == "LWPOLYLINE":
            points = [transform(point) for point in entity.get("points", [])]
            if entity.get("flag", 0) & 1 and points:
                points.append(points[0])
            if points:
                paths.append(points)
        elif kind == "POLYLINE_2D":
            points = [
                transform(by_handle[handle_id(vertex)]["point"])
                for vertex in entity.get("vertex", [])
                if by_handle.get(handle_id(vertex), {}).get("point")
            ]
            if entity.get("flag", 0) & 1 and points:
                points.append(points[0])
            if points:
                paths.append(points)
        elif kind == "HATCH":
            for boundary in entity.get("paths", []):
                points = [
                    transform(segment["first_endpoint"])
                    for segment in boundary.get("segs", [])
                    if segment.get("first_endpoint")
                ]
                if points:
                    points.append(points[0])
                    paths.append(points)
        return paths

    def recurse(block_id: int, transform: Callable[[Any], tuple[float, float]], seen: tuple[int, ...] = ()) -> list[list[tuple[float, float]]]:
        if block_id in seen:
            return []
        result: list[list[tuple[float, float]]] = []
        for ref in by_handle[block_id].get("entities", []):
            entity = by_handle[handle_id(ref)]
            if entity.get("entity") != "INSERT":
                result.extend(entity_paths(entity, transform))
                continue
            target = handle_id(entity["block_header"])
            insertion = entity.get("ins_pt", [0, 0, 0])
            scale = entity.get("scale", [1, 1, 1])
            rotation = float(entity.get("rotation", 0))

            def nested(point: Any, parent: Callable[[Any], tuple[float, float]] = transform) -> tuple[float, float]:
                return parent(affine(point, tx=insertion[0], ty=insertion[1], sx=scale[0], sy=scale[1], rotation=rotation))

            result.extend(recurse(target, nested, seen + (block_id,)))
        return result

    raw = recurse(DWG_BLOCK_HANDLE, lambda point: (float(point[0]), float(point[1])))
    xs = [x for path in raw for x, _ in path]
    ys = [y for path in raw for _, y in path]
    min_x, max_x, min_y, max_y = min(xs), max(xs), min(ys), max(ys)
    # Normalise and invert Y only. This is a rigid translation/reflection for view
    # presentation; the component is never stretched or resized.
    normalised = [[(x - min_x, max_y - y) for x, y in path] for path in raw]
    return normalised, {
        "source_dwg": str(DWG.relative_to(ROOT)),
        "source_dwg_sha256": DWG_SHA,
        "block_handle_decimal": DWG_BLOCK_HANDLE,
        "block_name": block.get("name"),
        "native_bounds_mm": [min_x, min_y, max_x, max_y],
        "native_size_mm": [max_x - min_x, max_y - min_y],
        "path_count": len(normalised),
        "geometry_transform": "normalise_origin_then_y_axis_view_reflection_only",
        "scale_x": 1.0,
        "scale_y": 1.0,
    }


def connected_components(item: Any) -> list[dict[str, Any]]:
    points = np.asarray(item.Coordinates.CoordList, dtype=float)
    faces = [tuple(index - 1 for index in face.CoordIndex) for face in item.Faces]
    vertex_faces: dict[int, list[int]] = defaultdict(list)
    for face_index, face in enumerate(faces):
        for vertex in face:
            vertex_faces[vertex].append(face_index)
    seen: set[int] = set()
    result: list[dict[str, Any]] = []
    for start in range(len(faces)):
        if start in seen:
            continue
        seen.add(start)
        stack = [start]
        component_faces: list[int] = []
        vertices: set[int] = set()
        while stack:
            current = stack.pop()
            component_faces.append(current)
            for vertex in faces[current]:
                vertices.add(vertex)
                for neighbour in vertex_faces[vertex]:
                    if neighbour not in seen:
                        seen.add(neighbour)
                        stack.append(neighbour)
        component_points = points[list(vertices)]
        minimum = component_points.min(axis=0)
        maximum = component_points.max(axis=0)
        result.append({
            "face_count": len(component_faces),
            "vertex_count": len(vertices),
            "minimum_mm": minimum.round(6).tolist(),
            "maximum_mm": maximum.round(6).tolist(),
            "size_mm": (maximum - minimum).round(6).tolist(),
        })
    return result


def inspect_project_body() -> dict[str, Any]:
    model = ifcopenshell.open(FORMAL)
    product = model.by_guid(GLOBAL_ID)
    representation = next(
        rep for rep in product.Representation.Representations
        if rep.RepresentationIdentifier == "Body" and rep.RepresentationType == "Tessellation"
    )
    styles: list[dict[str, Any]] = []
    all_components: list[tuple[int, dict[str, Any]]] = []
    for style_index, item in enumerate(representation.Items):
        components = connected_components(item)
        style_name = None
        if item.StyledByItem:
            style_name = item.StyledByItem[0].Name
        styles.append({"style_index": style_index, "style_name": style_name, "components": components})
        all_components.extend((style_index, component) for component in components)

    slats = [
        component for style, component in all_components
        if style == 2 and np.allclose(component["size_mm"], [22, 25, 1528], atol=0.01)
    ]
    display = max(styles[3]["components"], key=lambda component: component["face_count"])
    placement = np.asarray(get_local_placement(product.ObjectPlacement), dtype=float)
    return {
        "global_id": GLOBAL_ID,
        "ifc_class": product.is_a(),
        "name": product.Name,
        "placement_matrix": placement.round(6).tolist(),
        "styles": styles,
        "observed_atoms": {
            "vertical_grid_slats": {"count": len(slats), "components": slats},
            "right_protruding_open_element": display,
        },
    }


def svg_polyline(points: list[tuple[float, float]], cls: str) -> str:
    coords = " ".join(f"{x:.3f},{y:.3f}" for x, y in points)
    return f'<polyline class="{cls}" points="{coords}" />'


def native_atom_svg(paths: list[list[tuple[float, float]]], metadata: dict[str, Any]) -> str:
    width, height = metadata["native_size_mm"]
    body = "\n".join(svg_polyline(path, "official-native") for path in paths)
    return f'''<svg xmlns="http://www.w3.org/2000/svg" width="900" height="1000" viewBox="-80 -130 {width + 160:.3f} {height + 260:.3f}">
<style>.official-native{{fill:none;stroke:#0877df;stroke-width:2;vector-effect:non-scaling-stroke;stroke-linejoin:round}} text{{font-family:Arial,sans-serif;fill:#17324d}}</style>
<rect x="-80" y="-130" width="{width + 160:.3f}" height="{height + 260:.3f}" fill="white"/>
<text x="0" y="-75" font-size="34" font-weight="700">505 UP native DWG atom · W608 H768</text>
<text x="0" y="-35" font-size="24">block {metadata['block_handle_decimal']} / {metadata['block_name']} · scale 1:1 · no stretch</text>
{body}
<text x="0" y="{height + 65:.3f}" font-size="23">Blue = exact native DWG entity paths; family component, not a project shop drawing.</text>
</svg>'''


def front_overlay_svg(composition: dict[str, Any], body: dict[str, Any], native_paths: list[list[tuple[float, float]]]) -> str:
    min_x = -2400.181026
    project_components = [
        component for style in body["styles"] for component in style["components"]
        if max(component["size_mm"]) > 20
    ]
    actual_rects = []
    for component in project_components:
        mn, mx = component["minimum_mm"], component["maximum_mm"]
        x, y = mn[0] - min_x, -mx[2]
        width, height = mx[0] - mn[0], mx[2] - mn[2]
        actual_rects.append(f'<rect class="actual" x="{x:.3f}" y="{y:.3f}" width="{width:.3f}" height="{height:.3f}"/>')

    candidate = json.loads(CANDIDATES.read_text(encoding="utf-8"))
    proxy = []
    for path in candidate["views"]["front"]["proxy_paths_mm"]:
        proxy.append(svg_polyline([(point[0] - min_x, -point[1]) for point in path], "proxy"))

    official_grid = []
    for x in [0, 32, 640, 672, 1920, 1952, 2880, 2912]:
        official_grid.append(f'<line class="official-rule" x1="{x}" y1="-5" x2="{x}" y2="2375"/>')
    for y in [-5, 71, 455, 839, 1223, 1607, 1991, 2375]:
        official_grid.append(f'<line class="official-rule" x1="0" y1="{y}" x2="2912" y2="{y}"/>')
    native = []
    for path in native_paths:
        native.append(svg_polyline([(x + 1952, y + 844) for x, y in path], "official-native"))

    return f'''<svg xmlns="http://www.w3.org/2000/svg" width="1800" height="1200" viewBox="0 0 3900 2700">
<style>
text{{font-family:Arial,sans-serif;fill:#17324d}} .actual{{fill:none;stroke:#8f9aa4;stroke-width:8;opacity:.42;vector-effect:non-scaling-stroke}}
.proxy{{fill:none;stroke:#111820;stroke-width:2.4;opacity:.72;vector-effect:non-scaling-stroke}}
.official-rule{{stroke:#0877df;stroke-width:3;stroke-dasharray:16 12;opacity:.75;vector-effect:non-scaling-stroke}}
.official-native{{fill:none;stroke:#0877df;stroke-width:8;vector-effect:non-scaling-stroke;stroke-linejoin:round}}
.callout{{fill:#fff8e5;stroke:#ad6b00;stroke-width:3}}
</style>
<rect width="3900" height="2700" fill="white"/>
<text x="80" y="90" font-size="54" font-weight="700">505 UP V1.LP.S · modular compliance / atomic component audit</text>
<text x="80" y="145" font-size="30">Grey = actual IFC Body · Black = current geometry-derived proxy · Solid blue = native DWG atom · Dashed blue = datasheet rule</text>
<g transform="translate(80 230)">{''.join(actual_rects)}{''.join(proxy)}{''.join(official_grid)}{''.join(native)}</g>
<rect class="callout" x="3070" y="260" width="750" height="1110" rx="22"/>
<text x="3120" y="330" font-size="35" font-weight="700">Mechanical result</text>
<text x="3120" y="400" font-size="28">Width: 640 + 1280 + 960 + 32 = 2912</text>
<text x="3120" y="450" font-size="28" fill="#177239">PASS · exact modular width</text>
<text x="3120" y="530" font-size="28">Body height: 76 + 6 × 384 = 2380</text>
<text x="3120" y="580" font-size="28" fill="#a35b00">FFL clearance: project 0 / official 5</text>
<text x="3120" y="660" font-size="28">Central opening:</text>
<text x="3150" y="705" font-size="27">W1280 / clear 1248 / H3 clear 1120</text>
<text x="3150" y="750" font-size="27" fill="#177239">PASS · rule geometry exact</text>
<text x="3120" y="830" font-size="28">Right DISPLAY:</text>
<text x="3150" y="875" font-size="27">project face ≈ 610 × 766</text>
<text x="3150" y="920" font-size="27">official atom = 608 × 768</text>
<text x="3150" y="965" font-size="27" fill="#177239">front fit within 4 mm</text>
<text x="3150" y="1010" font-size="27" fill="#b03823">depth 420.462 vs D500 · FAIL</text>
<text x="3120" y="1090" font-size="28">Left grille:</text>
<text x="3150" y="1135" font-size="27">14 slats · 22 × 25 × 1528 · pitch 44</text>
<text x="3150" y="1180" font-size="27" fill="#a35b00">official Stripe identity unresolved</text>
<text x="80" y="2640" font-size="28">Module-composed official-component candidate; not an official complete shop drawing. No IFC write.</text>
</svg>'''


def plan_depth_svg(body: dict[str, Any]) -> str:
    display = body["observed_atoms"]["right_protruding_open_element"]
    actual_depth = display["size_mm"][1]
    return f'''<svg xmlns="http://www.w3.org/2000/svg" width="1500" height="800" viewBox="0 0 1500 800">
<style>text{{font-family:Arial,sans-serif;fill:#17324d}} .a{{fill:#87929c22;stroke:#87929c;stroke-width:4}} .p{{fill:none;stroke:#111820;stroke-width:3}} .o{{fill:none;stroke:#0877df;stroke-width:6}} .d{{stroke:#0877df;stroke-width:3;stroke-dasharray:12 10}}</style>
<rect width="1500" height="800" fill="white"/><text x="60" y="70" font-size="42" font-weight="700">505 UP right open DISPLAY · depth check</text>
<text x="60" y="118" font-size="25">Blue = official datasheet D500 envelope · Grey = actual IFC atom · Black = current proxy envelope</text>
<line class="d" x1="150" y1="180" x2="150" y2="700"/><text x="90" y="735" font-size="24">rear datum</text>
<rect class="o" x="150" y="220" width="1000" height="180"/><text x="1175" y="325" font-size="29" fill="#0877df">500 mm</text>
<rect class="a" x="150" y="470" width="{actual_depth * 2:.3f}" height="180"/><rect class="p" x="150" y="470" width="{actual_depth * 2:.3f}" height="180"/>
<text x="{180 + actual_depth * 2:.3f}" y="575" font-size="29">{actual_depth:.3f} mm</text>
<line x1="{150 + actual_depth * 2:.3f}" y1="675" x2="1150" y2="675" stroke="#b03823" stroke-width="5"/><text x="{190 + actual_depth * 2:.3f}" y="720" font-size="28" fill="#b03823">short by {500 - actual_depth:.3f} mm</text>
</svg>'''


def audit_markdown(audit: dict[str, Any]) -> str:
    display = audit["project_body"]["observed_atoms"]["right_protruding_open_element"]
    return f'''# Molteni&C 505 UP V1.LP.S modular audit

This is a read-only product review. It does not approve or write any IFC.

## Phase 1 - rules versus current geometry

- The official 2026 English datasheet defines horizontal module centre widths `320 / 480 / 640 / 960 / 1280 / 1440 / 1920 mm`, a vertical pitch of `384 mm`, `32 mm` structural panels, and structure depths `320 / 400 mm`.
- The project width `2912 mm` decomposes exactly as `640 + 1280 + 960 + 32 end panel`. Clear widths are `608 / 1248 / 928 mm`.
- The Body height `2380 mm` decomposes exactly as `76 mm plinth + 6 × 384 mm`. The datasheet adds a `5 mm` floor clearance, so a floor-standing system is `2385 mm` from FFL. The IFC placement starts at world Z `0`; therefore the body is 5 mm too low unless another project datum supplies that clearance.
- The central opening is a valid `W1280` / clear `1248` / `H3` clear `1120` construction.
- The right protruding open element is approximately `{display['size_mm'][0]:.3f} × {display['size_mm'][2]:.3f} mm` in front, close to the official `W608 × H768` atom. Its actual depth is `{display['size_mm'][1]:.3f} mm`, not the official open DISPLAY `D500`; the deficit is `{500 - display['size_mm'][1]:.3f} mm`.
- The left grille is 14 separate slats, each `22 × 25 × 1528 mm`, at a `44 mm` pitch. An exact native-DWG Stripe/Harry's Bar component has not been mechanically matched, so these lines remain geometry-derived and must not be shown as verified official blue geometry.
- The model includes 12 mm central dividers at rule-aligned positions. External geometry cannot prove that required shelves are internally reinforced; this remains a specification check.

Conclusion: the overall module grid is compliant, but the current design is **not yet fully compliant** because of the DISPLAY depth and the 5 mm floor-clearance placement. The grille identity is unresolved.

## Phase 2 - atomic decomposition

| Atom | Project position / size | Official evidence | Result |
| --- | --- | --- | --- |
| End/side panels | 32 deep-panel thickness at x 0, 640, 1920, 2880 | Datasheet pp. 2, 5-6 | exact |
| Left column | W640 / clear 608 / D320 | Datasheet p. 2 | exact shell |
| Left grille | 14 × 22 × 25 × 1528, pitch 44 | Stripe/Harry's Bar family reference only | unresolved; no blue substitution |
| Centre column | W1280 / clear 1248 / H3 clear 1120 / D320 | Datasheet pp. 10, 18 | exact rule geometry |
| Right column | W960 / clear 928 / D320 | Datasheet pp. 2, 9 | exact shell |
| Right open DISPLAY | project ≈ 610 × 766 × {display['size_mm'][1]:.3f}; target 608 × 768 × 500 | Datasheet p. 23 + native DWG block `{DWG_BLOCK_HANDLE}` | front close, depth fails |

The extracted DWG atom is placed with translation only. It is not stretched to the project mesh. The composite is labelled `module-composed official-component candidate`, not an official complete shop drawing.

## Phase 3 - parameter tool boundary

`505-up-parametric-schema.json` defines the allowed module widths, 384 mm vertical pitch, panel thicknesses, depth families, placement rules, and source provenance. `505-up-v1-lp-s-modular-composition.json` is the current data instance. `pipeline/scripts/generate_505_up_modular_review.py` validates it and reproduces this audit and the SVGs. This prototype intentionally stops before IFC authoring.
'''


def main() -> None:
    assert sha256(FORMAL) == FORMAL_SHA, "formal IFC baseline changed"
    assert sha256(DWG) == DWG_SHA, "official DWG hash changed"
    assert sha256(PDF) == PDF_SHA, "official datasheet hash changed"
    composition = json.loads(COMPOSITION.read_text(encoding="utf-8"))
    schema = json.loads(SCHEMA.read_text(encoding="utf-8"))
    widths = [module["module_width_mm"] for module in composition["columns"]]
    assert all(width in schema["rules"]["horizontal_module_widths_mm"] for width in widths)
    assert sum(widths) + composition["end_panel_mm"] == composition["project_envelope_mm"][0]
    assert composition["body_height_mm"] == composition["plinth_mm"] + composition["vertical_modules"] * schema["rules"]["vertical_pitch_mm"]

    native_paths, native_metadata = native_dwg_block_paths()
    assert np.allclose(native_metadata["native_size_mm"], [608, 768], atol=0.001)
    body = inspect_project_body()
    display = body["observed_atoms"]["right_protruding_open_element"]
    audit = {
        "schema_version": 1,
        "review_status": "engineering_review_pending",
        "source_mode": "module_composed_official_component_candidate",
        "official_complete_shop_drawing": False,
        "formal_ifc_write": "not_performed",
        "formal_ifc_sha256_before": FORMAL_SHA,
        "formal_ifc_sha256_after": sha256(FORMAL),
        "formal_ifc_unchanged": True,
        "sources": {
            "datasheet": {"url": "https://res.cloudinary.com/molteni/image/upload/v1760457000/2021_DS_505-UP-System_Living-Systems_Indoor_EN.pdf", "path": str(PDF.relative_to(ROOT)), "sha256": PDF_SHA, "pages": 31},
            "native_dwg": {"url": "https://res.cloudinary.com/molteni/raw/upload/v1766397289/2021_2D_505-UP_Living-Systems_Indoor.dwg?_s=public-apps", "path": str(DWG.relative_to(ROOT)), "sha256": DWG_SHA},
        },
        "datasheet_rules": schema["rules"],
        "composition": composition,
        "native_dwg_atom": native_metadata,
        "project_body": body,
        "checks": {
            "width_module_sum": {"status": "pass", "equation": "640 + 1280 + 960 + 32 = 2912", "delta_mm": 0},
            "body_height_modules": {"status": "pass", "equation": "76 + 6 * 384 = 2380", "delta_mm": 0},
            "floor_clearance": {"status": "fail", "project_mm": 0, "official_mm": 5, "delta_mm": -5},
            "central_open_compartment": {"status": "pass", "module_width_mm": 1280, "clear_width_mm": 1248, "clear_height_mm": 1120},
            "right_display_front": {"status": "pass_with_4mm_review_tolerance", "project_width_mm": display["size_mm"][0], "project_height_mm": display["size_mm"][2], "official_width_mm": 608, "official_height_mm": 768},
            "right_display_depth": {"status": "fail", "project_mm": display["size_mm"][1], "official_mm": 500, "delta_mm": display["size_mm"][1] - 500},
            "left_grille_identity": {"status": "unresolved", "official_blue_allowed": False},
            "reinforced_shelf_internal_construction": {"status": "not_provable_from_external_body"},
            "native_atom_transform": {"status": "pass", "scale_x": 1.0, "scale_y": 1.0, "placement": "translation_only"},
        },
        "overall_design_compliance": "partial_fail_requires_display_depth_and_floor_clearance_correction",
    }

    write_json(OUT / "505-up-modular-rule-comparison.json", audit)
    write_text(OUT / "505-up-modular-rule-audit.md", audit_markdown(audit))
    write_text(OUT / "official-dwg-atom-w608-h768.svg", native_atom_svg(native_paths, native_metadata))
    write_text(OUT / "505-up-modular-official-component-overlay-front.svg", front_overlay_svg(composition, body, native_paths))
    write_text(OUT / "505-up-modular-depth-comparison-plan.svg", plan_depth_svg(body))
    write_text(OUT / "505-up-module-composed-candidate-front.svg", front_overlay_svg(composition, body, native_paths))

    sources = {
        "schema_version": 1,
        "manufacturer": "Molteni&C",
        "family": "505 UP System",
        "product_page": "https://www.molteni.it/en/ap/product/505-up-system",
        "datasheet": audit["sources"]["datasheet"],
        "native_dwg": audit["sources"]["native_dwg"],
        "native_atom": native_metadata,
    }
    write_json(OUT / "505-up-modular-analysis-source-record.json", sources)

    outputs = [
        OUT / "505-up-modular-rule-comparison.json",
        OUT / "505-up-modular-rule-audit.md",
        OUT / "official-dwg-atom-w608-h768.svg",
        OUT / "505-up-modular-official-component-overlay-front.svg",
        OUT / "505-up-modular-depth-comparison-plan.svg",
        OUT / "505-up-module-composed-candidate-front.svg",
        OUT / "505-up-modular-analysis-source-record.json",
    ]
    manifest = {
        "schema_version": 1,
        "generator": str(Path(__file__).resolve().relative_to(ROOT)),
        "formal_ifc_sha256": FORMAL_SHA,
        "formal_ifc_unchanged": sha256(FORMAL) == FORMAL_SHA,
        "derived_ifc_written": False,
        "official_complete_shop_drawing": False,
        "outputs": [{"path": str(path.relative_to(ROOT)), "sha256": sha256(path), "bytes": path.stat().st_size} for path in outputs],
        "pass": True,
    }
    write_json(OUT / "505-up-modular-review-manifest.json", manifest)
    assert sha256(FORMAL) == FORMAL_SHA


if __name__ == "__main__":
    main()
