#!/usr/bin/env python3
"""Generate the Molteni 505 UP project-composition review.

The Plan proxy retains the whole-product silhouette and mechanically restores
only the interface seam between the main cabinet envelope and the independently
modelled protruding DISPLAY component. It never closes the entire DISPLAY
footprint over the main cabinet.
"""

import argparse
import json
import math
import re
from datetime import datetime, timezone
from pathlib import Path

import ifcopenshell
from shapely.geometry import LineString, Polygon
from shapely.ops import unary_union

from falper_sorgente_linework import ROOT, load_json, relative, sha256, write_json
from int1_highpoly_type_review import (
    display_edge_sample,
    mesh_for_one_product,
    product_type_name,
    projected_raw_edges,
    projected_silhouette,
    svg_path,
)
import gessi316_54294_review as shared


shared.PROFILE_KEY = "505-up-v1-lp-s"
shared.ARTICLE_NUMBER = "505 UP System / project V1.LP.S"
shared.GENERATOR = "pipeline/scripts/505_up_v1_lp_s_review.py"
shared.OFFICIAL_CAD_STATUS = "official_native_family_dwg_archived_exact_project_configuration_not_matched"
shared.OFFICIAL_CAD_ACQUIRED = True
shared.OFFICIAL_CAD_EXACT_PROJECT_CONFIGURATION_MATCH = False
shared.OUTPUT_DIR = ROOT / "output/review/highpoly-types/505-up-v1-lp-s"
shared.REGISTER = shared.OUTPUT_DIR / "profile.json"
shared.ACCESS_RECORD = shared.OUTPUT_DIR / "official-source/source-access-record.json"
shared.EXPECTED_DESCRIPTION = None

FORMAL_IFC = ROOT / "2504 GBTB Yanlord Zhuhai.ifc"
FORMAL_SHA256 = "7a50b87e8f48a7c2bfbaa9f1dabdfca7155fa684b2325c66aa1ae15c8ab0a25c"
APPROVAL = ROOT / "pipeline/decisions/505-up-v1-lp-s-drawing-approval.json"
SOURCE_KIND = "geometry_derived_simplified_proxy"
DISPLAY_WIDTH_MM = 610.0
DISPLAY_DEPTH_MM = 420.462
DISPLAY_HEIGHT_MM = 766.0
DISPLAY_SIZE_TOLERANCE_MM = 0.75
BASE_PATH_COUNTS = {"plan": 19, "front": 281, "side": 93}
REPAIRED_PATH_COUNTS = {"plan": 20, "front": 281, "side": 93}


def mesh_components(vertices, faces):
    """Return connected mesh components without relying on object names."""
    parent = list(range(len(vertices)))

    def find(index):
        while parent[index] != index:
            parent[index] = parent[parent[index]]
            index = parent[index]
        return index

    for face in faces:
        anchor = face[0]
        for index in face[1:]:
            root_a, root_b = find(anchor), find(index)
            if root_a != root_b:
                parent[root_b] = root_a
    grouped = {}
    for index in range(len(vertices)):
        grouped.setdefault(find(index), []).append(index)
    records = []
    for root, indices in grouped.items():
        points = [vertices[index] for index in indices]
        minimum = [min(point[axis] for point in points) for axis in range(3)]
        maximum = [max(point[axis] for point in points) for axis in range(3)]
        records.append(
            {
                "root": root,
                "indices": indices,
                "minimum": minimum,
                "maximum": maximum,
                "size": [maximum[axis] - minimum[axis] for axis in range(3)],
            }
        )
    return records, find


def polygon_area(path):
    return abs(
        sum(
            start[0] * end[1] - end[0] * start[1]
            for start, end in zip(path, path[1:])
        )
    ) / 2.0


def component_plan_geometry(component, vertices, faces, find, simplify_mm):
    """Return the mechanically projected Plan polygon for one mesh component."""
    remap = {old: new for new, old in enumerate(component["indices"])}
    component_vertices = [vertices[index] for index in component["indices"]]
    component_faces = [
        tuple(remap[index] for index in face)
        for face in faces
        if find(face[0]) == component["root"]
    ]
    paths = projected_silhouette(
        component_vertices, component_faces, (0, 1), simplify_mm
    )
    polygons = []
    for path in paths:
        if len(path) < 4:
            continue
        polygon = Polygon(path)
        if polygon.is_valid and polygon.area >= 1.0:
            polygons.append(polygon)
    if not polygons:
        raise RuntimeError(f"component {component['root']} has no Plan polygon")
    return unary_union(polygons)


def plan_component_interface(vertices, faces, simplify_mm):
    """Extract the true DISPLAY/main-envelope interface and component evidence."""
    components, find = mesh_components(vertices, faces)
    display_matches = [
        component
        for component in components
        if abs(component["size"][0] - DISPLAY_WIDTH_MM) <= DISPLAY_SIZE_TOLERANCE_MM
        and abs(component["size"][1] - DISPLAY_DEPTH_MM) <= DISPLAY_SIZE_TOLERANCE_MM
        and abs(component["size"][2] - DISPLAY_HEIGHT_MM) <= DISPLAY_SIZE_TOLERANCE_MM
    ]
    if len(display_matches) != 1:
        raise RuntimeError(
            "expected one mechanically identified 505 DISPLAY component, found "
            f"{len(display_matches)}"
        )
    display = display_matches[0]
    main = max(
        (component for component in components if component is not display),
        key=lambda component: component["size"][0] * component["size"][1],
    )
    covers = [
        component
        for component in components
        if abs(component["maximum"][2]) <= DISPLAY_SIZE_TOLERANCE_MM
        and 75.0 <= component["size"][2] <= 77.0
        and abs(component["size"][1] - main["size"][1])
        <= DISPLAY_SIZE_TOLERANCE_MM
        and component["size"][0] >= 500.0
    ]
    slats = [
        component
        for component in components
        if abs(component["size"][0] - 22.0) <= DISPLAY_SIZE_TOLERANCE_MM
        and abs(component["size"][1] - 25.0) <= DISPLAY_SIZE_TOLERANCE_MM
        and abs(component["size"][2] - 1528.0) <= DISPLAY_SIZE_TOLERANCE_MM
    ]
    if len(covers) != 2 or len(slats) != 14:
        raise RuntimeError(
            f"505 Plan component semantics drifted: covers={len(covers)} slats={len(slats)}"
        )
    slat_centres = sorted(
        (component["minimum"][0] + component["maximum"][0]) / 2.0
        for component in slats
    )
    slat_spacing = [
        slat_centres[index + 1] - slat_centres[index]
        for index in range(len(slat_centres) - 1)
    ]
    if any(abs(value - 44.0) > DISPLAY_SIZE_TOLERANCE_MM for value in slat_spacing):
        raise RuntimeError(f"505 slat spacing drifted: {slat_spacing}")

    display_geometry = component_plan_geometry(
        display, vertices, faces, find, simplify_mm
    )
    interface_y = main["maximum"][1]
    probe = LineString(
        [
            (display["minimum"][0] - 1000.0, interface_y),
            (display["maximum"][0] + 1000.0, interface_y),
        ]
    )
    interface = display_geometry.intersection(probe)
    if interface.geom_type != "LineString" or interface.is_empty:
        raise RuntimeError(
            f"DISPLAY/main envelope must yield one interface line, got {interface.wkt}"
        )
    seam = [(float(x), float(y)) for x, y in interface.coords]
    if len(seam) != 2 or abs(interface.length - DISPLAY_WIDTH_MM) > 1.0:
        raise RuntimeError(f"DISPLAY interface length drifted: {interface.length}")

    def component_record(component):
        return {
            "root_vertex_index": component["root"],
            "minimum_mm": [round(value, 6) for value in component["minimum"]],
            "maximum_mm": [round(value, 6) for value in component["maximum"]],
            "size_mm": [round(value, 6) for value in component["size"]],
        }

    audit = {
        "method": "connected_components_then_display_main_envelope_intersection",
        "main_envelope_component": component_record(main),
        "display_component": component_record(display),
        "cover_plate_components": [
            component_record(component)
            for component in sorted(covers, key=lambda item: item["minimum"][0])
        ],
        "slat_components": {
            "count": len(slats),
            "centres_x_mm": [round(value, 6) for value in slat_centres],
            "spacing_mm": [round(value, 6) for value in slat_spacing],
            "component_size_mm": [22.0, 25.0, 1528.0],
        },
        "interface": {
            "main_envelope_front_y_mm": round(interface_y, 6),
            "path_mm": [[round(x, 6), round(y, 6)] for x, y in seam],
            "length_mm": round(interface.length, 6),
            "closed": False,
            "full_display_footprint_added": False,
        },
        "plan_screen_semantics_after_world_placement": {
            "rear_top_cover": "local y 0..320 mm; world y -1384..-1064 mm",
            "right_slats": "world x 5095..5697 mm; world y -1060..-1035 mm",
            "front_left_display_protrusion": "world x 3153..3763 mm; world y -1383.963..-963.501 mm",
        },
        "pass": True,
    }
    return seam, audit


def structural_proxy_builder(vertices, faces, axes, simplify_mm, view):
    """Keep the Body silhouette plus true mesh creases, without triangle diagonals."""
    silhouette = projected_silhouette(vertices, faces, axes, simplify_mm)
    if view == "plan":
        seam, _ = plan_component_interface(vertices, faces, simplify_mm)
        return silhouette + [seam]
    edge_normals = {}
    for face in faces:
        a, b, c = (vertices[index] for index in face)
        ab = tuple(b[index] - a[index] for index in range(3))
        ac = tuple(c[index] - a[index] for index in range(3))
        normal = (
            ab[1] * ac[2] - ab[2] * ac[1],
            ab[2] * ac[0] - ab[0] * ac[2],
            ab[0] * ac[1] - ab[1] * ac[0],
        )
        length = math.sqrt(sum(value * value for value in normal))
        if length <= 1e-9:
            continue
        normal = tuple(value / length for value in normal)
        for start, end in ((face[0], face[1]), (face[1], face[2]), (face[2], face[0])):
            edge_normals.setdefault(tuple(sorted((start, end))), []).append(normal)

    crease_limit = math.cos(math.radians(3.0))
    axis_intervals = {"h": {}, "v": {}}
    for (start, end), normals in edge_normals.items():
        sharp = len(normals) == 1 or any(
            abs(sum(a * b for a, b in zip(normals[i], normals[j]))) < crease_limit
            for i in range(len(normals))
            for j in range(i + 1, len(normals))
        )
        if not sharp:
            continue
        p1 = (vertices[start][axes[0]], vertices[start][axes[1]])
        p2 = (vertices[end][axes[0]], vertices[end][axes[1]])
        dx, dy = p2[0] - p1[0], p2[1] - p1[1]
        if abs(dy) <= 0.8 and abs(dx) >= 24.0:
            coordinate = round(((p1[1] + p2[1]) / 2.0) * 2.0) / 2.0
            axis_intervals["h"].setdefault(coordinate, []).append(sorted((p1[0], p2[0])))
        elif abs(dx) <= 0.8 and abs(dy) >= 24.0:
            coordinate = round(((p1[0] + p2[0]) / 2.0) * 2.0) / 2.0
            axis_intervals["v"].setdefault(coordinate, []).append(sorted((p1[1], p2[1])))

    structural = []
    for direction, groups in axis_intervals.items():
        for coordinate, intervals in groups.items():
            intervals.sort()
            merged = []
            for start, end in intervals:
                if merged and start <= merged[-1][1] + 1.0:
                    merged[-1][1] = max(merged[-1][1], end)
                else:
                    merged.append([start, end])
            for start, end in merged:
                if end - start < 30.0:
                    continue
                structural.append(
                    [(start, coordinate), (end, coordinate)]
                    if direction == "h"
                    else [(coordinate, start), (coordinate, end)]
                )
    return silhouette + structural


shared.PROXY_BUILDER = structural_proxy_builder


def rounded(paths):
    return [
        [[round(float(x), 6), round(float(y), 6)] for x, y in path]
        for path in paths
    ]


def update_plan_svg(path: Path, raw_edges, proxy):
    """Replace only the black proxy path while preserving the existing review sheet."""
    points = [point for edge in raw_edges for point in edge]
    points.extend(point for item in proxy for point in item)
    min_x, max_x = min(point[0] for point in points), max(point[0] for point in points)
    min_y, max_y = min(point[1] for point in points), max(point[1] for point in points)
    padding = max(max_x - min_x, max_y - min_y) * 0.07
    min_x, max_x = min_x - padding, max_x + padding
    min_y, max_y = min_y - padding, max_y + padding
    scale = min(980.0 / (max_x - min_x), 750.0 / (max_y - min_y))

    def transform(point):
        return (
            60.0 + (point[0] - min_x) * scale,
            170.0 + 750.0 - (point[1] - min_y) * scale,
        )

    data = " ".join(
        svg_path([item], transform, close=item[0] == item[-1])
        for item in proxy
    )
    text = path.read_text(encoding="utf-8")
    text = re.sub(r' data-component-boundary-extraction="[^"]*"', "", text)
    text = re.sub(r' data-component-interface-extraction="[^"]*"', "", text)
    text = re.sub(r' data-full-display-footprint-added="[^"]*"', "", text)
    pattern = re.compile(
        r'(<path class="simplified-proxy-silhouette geometry-derived"[^>]*?) d="[^"]*"'
    )
    replacement = (
        r'\1 data-component-interface-extraction="mechanical-connected-components"'
        r' data-display-component-size-mm="610,420.462,766"'
        r' data-full-display-footprint-added="false"'
        f' d="{data}"'
    )
    updated, count = pattern.subn(replacement, text, count=1)
    if count != 1:
        raise RuntimeError("505 Plan SVG proxy path was not uniquely located")
    updated = re.sub(
        r"derived proxy paths: \d+",
        f"derived proxy paths: {len(proxy)}",
        updated,
        count=1,
    )
    path.write_text(updated, encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=FORMAL_IFC)
    parser.add_argument("--register", type=Path, default=shared.REGISTER)
    parser.add_argument("--output", type=Path, default=shared.OUTPUT_DIR)
    args = parser.parse_args()
    source = args.input.resolve()
    output = args.output.resolve()
    if sha256(source) != FORMAL_SHA256:
        raise RuntimeError("formal IFC hash mismatch")
    profile = load_json(args.register.resolve())["profiles"][shared.PROFILE_KEY]
    candidate_path = output / "candidate-representations.json"
    manifest_path = output / "manifest.json"
    candidate = load_json(candidate_path)
    manifest = load_json(manifest_path)
    product, vertices, faces = mesh_for_one_product(
        ifcopenshell.open(source), profile["representative_global_id"]
    )
    actual_identity = product_type_name(product) or product.Name
    if actual_identity != profile["ifc_type_name"]:
        raise RuntimeError("505 representative identity drifted")
    old_hashes = {
        name: sha256(output / f"{name}.svg") for name in ("front", "side")
    }
    old_paths = {
        name: candidate["views"][name]["proxy_paths_mm"]
        for name in ("front", "side")
    }
    old_counts = {
        name: len(candidate["views"][name]["proxy_paths_mm"])
        for name in ("plan", "front", "side")
    }
    if old_counts not in (BASE_PATH_COUNTS, REPAIRED_PATH_COUNTS):
        raise RuntimeError(f"505 candidate pre-state path counts drifted: {old_counts}")

    simplify_mm = float(profile["silhouette_simplify_mm"])
    seam, audit = plan_component_interface(vertices, faces, simplify_mm)
    proxy = projected_silhouette(vertices, faces, (0, 1), simplify_mm) + [seam]
    if len(proxy) != REPAIRED_PATH_COUNTS["plan"]:
        raise RuntimeError(f"505 repaired Plan path count drifted: {len(proxy)}")
    raw_edges = display_edge_sample(
        projected_raw_edges(vertices, faces, (0, 1)), maximum=1000
    )
    plan_svg = output / "plan.svg"
    update_plan_svg(plan_svg, raw_edges, proxy)

    candidate["views"]["plan"]["proxy_paths_mm"] = rounded(proxy)
    candidate["views"]["plan"].pop("component_boundary_extraction", None)
    candidate["views"]["plan"]["component_semantics"] = audit
    candidate["review_status"] = "visual_review_pending"
    candidate["views"]["plan"]["semantic_path_count"] = 1
    candidate["views"]["plan"]["silhouette_path_count"] = len(proxy) - 1
    write_json(candidate_path, candidate)

    plan_record = next(item for item in manifest["views"] if item["view"] == "plan")
    plan_record["svg_sha256"] = sha256(plan_svg)
    plan_record["silhouette_path_count"] = len(proxy) - 1
    plan_record["semantic_path_count"] = 1
    plan_record["derived_proxy_path_count"] = len(proxy)
    plan_record.pop("component_boundary_extraction", None)
    plan_record["component_semantics"] = audit
    manifest["generated_at"] = datetime.now(timezone.utc).isoformat()
    manifest["candidate_representations_sha256"] = sha256(candidate_path)
    manifest.pop("plan_component_boundary_extraction", None)
    manifest["plan_component_semantics"] = audit
    manifest["review_status"] = "visual_review_pending"
    manifest["approved_for_drawing_ifc"] = False
    manifest["formal_ifc_sha256_after_generation"] = sha256(source)
    manifest["formal_ifc_bytes_unchanged"] = sha256(source) == FORMAL_SHA256
    manifest["pass"] = True
    write_json(manifest_path, manifest)

    approval = load_json(APPROVAL)
    if (
        approval.get("status") not in ("approved", "revision_pending_review")
        or approval.get("derived_ifc_write_allowed") is not True
        or approval.get("formal_authoritative_ifc_write_allowed") is not False
    ):
        raise RuntimeError("505 approval gate drifted")
    approval["candidate_manifest_sha256"] = sha256(manifest_path)
    approval["status"] = "revision_pending_review"
    approval["pending_reapproval_views"] = ["plan", "front"]
    approval["latest_review"] = {
        "outcome": "rejected",
        "front_user_evidence": "这儿的我怎么感觉左右颠倒啊？",
        "plan_user_evidence": "这里你加的线在哪儿呢？加错了吧。我能不能给你画一个边？不应该是上面一整块盖板，然后右边是几个格栅，然后前面左边是一个简单的凸出来吧。你怎么突然反倒变成一个整体了？",
    }
    approval["plan_component_boundary_repair"] = {
        "status": "previous_closed_footprint_rejected_replaced_by_interface_seam",
        "user_authorization": "返工 505 Plan；撤销错误整体闭合轮廓并按真实分件恢复",
        "path_count_before": old_counts["plan"],
        "path_count_after": len(proxy),
        "mechanical_extraction": audit,
    }
    write_json(APPROVAL, approval)

    if any(sha256(output / f"{name}.svg") != old_hashes[name] for name in old_hashes):
        raise RuntimeError("505 Front/Side SVG changed during Plan-only repair")
    if any(candidate["views"][name]["proxy_paths_mm"] != old_paths[name] for name in old_paths):
        raise RuntimeError("505 Front/Side proxy paths changed during Plan-only repair")
    print(
        json.dumps(
            {
                "plan_paths_before": old_counts["plan"],
                "plan_paths_after": len(proxy),
                "display_component": audit,
                "front_svg_sha256": old_hashes["front"],
                "side_svg_sha256": old_hashes["side"],
                "formal_ifc_sha256": sha256(source),
                "manifest": relative(manifest_path),
                "approval": relative(APPROVAL),
                "pass": True,
            },
            ensure_ascii=False,
            indent=2,
        )
    )


def write_index(manifest):
    cards = "".join(
        f'<article><h2>{item["view"].title()}</h2><a href="{item["view"]}.svg"><img src="{item["view"]}.svg"></a></article>'
        for item in manifest["views"]
    )
    bonsai_cards = "".join(
        f'<article><h2>Bonsai {view.title()}</h2><a href="bonsai-camera-{filename}.png"><img src="bonsai-camera-{filename}.png"></a></article>'
        for view, filename in (("plan", "plan"), ("front", "front-elevation"), ("side", "side-elevation"), ("iso", "iso"))
    )
    (shared.OUTPUT_DIR / "index.html").write_text(
        f'''<!doctype html><html><meta charset="utf-8"><title>Molteni 505 UP review</title><style>body{{font:16px Arial;margin:30px;background:#f5f3ef;color:#1f2d3d}}main{{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:18px}}article{{background:#fff;padding:14px;border-radius:10px}}img{{width:100%}}a{{margin-right:18px}}</style><h1>Molteni&C 505 UP System / project 505 UP V1.LP.S</h1><p>Official native family DWGs are archived below. No published catalogue composition matches the project-specific slatted-left / display-right arrangement, so the candidate drawing line remains {shared.SOURCE_LABEL_EN}; official CAD is not rearranged or falsely shown as blue project geometry.</p><nav><a href="manifest.json">Manifest</a><a href="candidate-representations.json">Candidate</a><a href="profile.json">Profile</a><a href="official-source/source-access-record.json">Source record</a><a href="project-context-furniture-plan-review.svg">Project furniture plan</a><a href="bonsai-review-manifest.json">Bonsai evidence</a><a href="official-source/molteni-505-up-technical-library-full-preview.svg">Official 2021 CAD SVG</a><a href="official-source/molteni-505-up-inspiring-solution-full-preview.svg">Official inspiring solutions CAD SVG</a><a href="https://www.molteni.it/en/ap/product/505-up-system">Official page</a></nav><h2>Project drawing context</h2><main><article><h2>Furniture Plan</h2><a href="project-context-furniture-plan-review.svg"><img src="project-context-furniture-plan-review-preview.png"></a></article></main><h2>Three-view project candidate</h2><main>{cards}</main><h2>Actual Bonsai camera renders from isolated IFC Body</h2><main>{bonsai_cards}</main><code>{manifest["formal_ifc_sha256"]}</code></html>''',
        encoding="utf-8",
    )


shared.write_index = write_index


if __name__ == "__main__":
    main()
