#!/usr/bin/env python3
"""Generate the review-only 505 UP Plan from first-visible semantic faces.

This batch intentionally stops at the single-product SVG. It does not touch
the persisted product-level drawing IFC or the Bonsai project Drawing.
"""

import json
import math
import sys
from collections import Counter
from hashlib import sha256 as hashlib_sha256
from importlib import import_module
from pathlib import Path

import ifcopenshell
from shapely.geometry import LineString, Point, Polygon
from shapely.ops import linemerge, unary_union

SCRIPT_DIR = Path(__file__).resolve().parent
ROOT = SCRIPT_DIR.parent.parent
sys.path.insert(0, str(SCRIPT_DIR))

from int1_highpoly_type_review import mesh_for_one_product  # noqa: E402

diagnose = import_module("diagnose_505_up_plan_visibility")
review = import_module("505_up_v1_lp_s_review")

FORMAL_IFC = ROOT / "2504 GBTB Yanlord Zhuhai.ifc"
FORMAL_SHA256 = "7a50b87e8f48a7c2bfbaa9f1dabdfca7155fa684b2325c66aa1ae15c8ab0a25c"
PRODUCT = ROOT / "output/review/highpoly-types/505-up-v1-lp-s"
GLOBAL_ID = "19MpdkWqXC7uhUNhLQgrce"
SOURCE_KIND = "geometry_derived_simplified_proxy"
SOURCE_LABEL_ZH = "基于原始高模几何生成的简化图纸表达"
USER_AUTHORIZATION = "没问题，那么基于此给我新的单品svg"
FRONT_SVG_SHA256 = "0f942b3799c9854119765445258a8eb422274c5cb198582e2fe369c5b13c975b"
SIDE_SVG_SHA256 = "8f4224ef5fb3594bd369cc850ff49bd728c266bba8ac966528a1c361ed3eb156"
DERIVED_IFC_SHA256 = "a49d1b9d2859fe596761ecc70c845cf888022f32c988a71946e8bfa4b0338417"
SESSION_IFC_SHA256 = "b8e9c5825225bd5f1df3a5ca5b5c9d0945e29acc9f1435c1e4dc98c3525b7687"
PLAN_SCENE_SHA256 = "f67a0a3867e105e6ffc2a4507ba07f40d233c72142d46c5c06e2c9555ccf29d5"

SEMANTIC_KEYS = [
    "A_main_top",
    "B_right_slats",
    "C_front_left_display",
    "D_mid_shelf_lip",
    "E_lower_front_rails",
]
SEMANTIC_LABELS = {
    "A_main_top": "A 柜体与背板最高顶面 / 2380 mm",
    "B_right_slats": "B 右侧 14 片格栅顶面 / 1976 mm",
    "C_front_left_display": "C 前方左侧 DISPLAY 顶面 / 1610.705 mm",
    "D_mid_shelf_lip": "D 中部前沿可见窄带 / 470 mm",
    "E_lower_front_rails": "E 下部前沿横向构件 / 440 mm",
}
SEMANTIC_COLOURS = {
    "A_main_top": "#4c78a8",
    "B_right_slats": "#17a398",
    "C_front_left_display": "#f28e2b",
    "D_mid_shelf_lip": "#edc948",
    "E_lower_front_rails": "#b279a2",
}


def sha256(path: Path) -> str:
    digest = hashlib_sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def round_point(point):
    return [round(float(point[0]), 6), round(float(point[1]), 6)]


def collect_lines(geometry):
    lines = []

    def visit(part):
        if part.geom_type in {"LineString", "LinearRing"}:
            lines.append(LineString(part.coords))
        elif hasattr(part, "geoms"):
            for child in part.geoms:
                visit(child)

    visit(geometry)
    return lines


def component_groups(components):
    main = max(components, key=lambda item: math.prod(item["size"]))
    back = next(
        item for item in components
        if abs(item["size"][0] - 2912.0) <= 0.5 and abs(item["size"][1] - 16.0) <= 0.5
    )
    slats = {
        item["root"] for item in components
        if abs(item["size"][0] - 22.0) <= 0.5
        and abs(item["size"][1] - 25.0) <= 0.5
        and abs(item["size"][2] - 1528.0) <= 0.5
    }
    display = next(
        item for item in components
        if abs(item["size"][0] - 610.0) <= 0.75
        and abs(item["size"][1] - 420.462) <= 0.75
        and abs(item["size"][2] - 766.0) <= 0.75
    )
    lip = next(
        item for item in components
        if abs(item["size"][0] - 1248.0) <= 0.5
        and abs(item["size"][1] - 330.0) <= 0.5
        and abs(item["size"][2] - 42.0) <= 0.5
    )
    rails = {
        item["root"] for item in components
        if abs(item["size"][1] - 32.0) <= 0.5
        and abs(item["size"][2] - 376.0) <= 0.5
    }
    if len(slats) != 14 or len(rails) != 2:
        raise RuntimeError(f"component classification drifted: slats={len(slats)} rails={len(rails)}")
    return main, back, slats, display, lip, rails


def visible_semantic_regions(vertices, faces, components, find):
    main, back, slats, display, lip, rails = component_groups(components)
    parts = {key: [] for key in SEMANTIC_KEYS}
    surfaces = diagnose.horizontal_top_surfaces(vertices, faces, components, find)
    for item in surfaces:
        key = diagnose.semantic_key(
            item, main["root"], back["root"], slats, display["root"], lip["root"], rails
        )
        if key in parts and item["visible_geometry"].area >= 10.0:
            parts[key].append(item["visible_geometry"])

    raw = {key: unary_union(items).buffer(0) for key, items in parts.items()}
    # Half a millimetre closes tessellation-only cracks while retaining the
    # model's 4 mm and larger component gaps. Repartitioning in descending
    # world height makes every Plan point belong to its first visible face.
    regions = {}
    covered = Polygon()
    for key in SEMANTIC_KEYS:
        cleaned = (
            raw[key]
            .buffer(0.5, join_style=2)
            .buffer(-0.5, join_style=2)
            .simplify(0.5, preserve_topology=True)
        )
        visible = cleaned.difference(covered).buffer(0)
        regions[key] = visible
        covered = unary_union([covered, visible]).buffer(0)
    return regions


def build_semantic_linework(regions):
    probes = {key: geometry.buffer(1e-6) for key, geometry in regions.items()}

    def classify(point):
        matches = [key for key in SEMANTIC_KEYS if probes[key].contains(point)]
        if len(matches) > 1:
            raise RuntimeError(f"semantic partition overlaps at {point.wkt}: {matches}")
        return matches[0] if matches else "background"

    network = unary_union([geometry.boundary for geometry in regions.values()])
    records = []
    retained = []
    for line in collect_lines(network):
        coordinates = list(line.coords)
        for start, end in zip(coordinates, coordinates[1:]):
            dx, dy = end[0] - start[0], end[1] - start[1]
            length = math.hypot(dx, dy)
            if length <= 1e-6:
                continue
            midpoint = ((start[0] + end[0]) / 2.0, (start[1] + end[1]) / 2.0)
            normal = (-dy / length, dx / length)
            side_a = classify(Point(midpoint[0] + normal[0] * 2.0, midpoint[1] + normal[1] * 2.0))
            side_b = classify(Point(midpoint[0] - normal[0] * 2.0, midpoint[1] - normal[1] * 2.0))
            valid = side_a != side_b
            record = {
                "start_mm": round_point(start),
                "end_mm": round_point(end),
                "length_mm": round(length, 6),
                "side_a": side_a,
                "side_b": side_b,
                "internal": side_a != "background" and side_b != "background",
                "semantic_separation_valid": valid,
            }
            records.append(record)
            if valid:
                retained.append(LineString([start, end]))

    if any(not record["semantic_separation_valid"] for record in records):
        raise RuntimeError("same-semantic or unexplained boundary survived normalization")
    merged = collect_lines(linemerge(unary_union(retained)))
    merged.sort(key=lambda line: (round(line.bounds[0], 6), round(line.bounds[1], 6), round(line.length, 6)))
    paths = [[round_point(point) for point in line.coords] for line in merged]
    pairs = Counter(
        tuple(sorted((record["side_a"], record["side_b"])))
        for record in records
        if record["internal"]
    )
    required = {
        ("A_main_top", "C_front_left_display"),
        ("A_main_top", "D_mid_shelf_lip"),
        ("C_front_left_display", "E_lower_front_rails"),
        ("D_mid_shelf_lip", "E_lower_front_rails"),
    }
    missing = required.difference(pairs)
    if missing:
        raise RuntimeError(f"required semantic interfaces missing: {sorted(missing)}")
    return paths, records, pairs


def geometry_path(geometry, transform):
    commands = []
    for polygon in diagnose.polygons(geometry):
        for ring in [polygon.exterior, *polygon.interiors]:
            points = [transform(float(x), float(y)) for x, y in ring.coords]
            commands.append("M " + " L ".join(f"{x:.3f},{y:.3f}" for x, y in points) + " Z")
    return " ".join(commands)


def linework_path(paths, transform):
    commands = []
    for path in paths:
        points = [transform(float(point[0]), float(point[1])) for point in path]
        command = "M " + " L ".join(f"{x:.3f},{y:.3f}" for x, y in points)
        if path[0] == path[-1]:
            command += " Z"
        commands.append(command)
    return " ".join(commands)


def transforms(regions, plot):
    geometry = unary_union(list(regions.values()))
    min_x, min_y, max_x, max_y = geometry.bounds
    screen_min_x, screen_max_x = -max_x, -min_x
    screen_min_y, screen_max_y = -max_y, -min_y
    span_x = screen_max_x - screen_min_x
    span_y = screen_max_y - screen_min_y
    scale = min(plot[2] / span_x, plot[3] / span_y)
    offset_x = plot[0] + (plot[2] - span_x * scale) / 2.0
    offset_y = plot[1] + (plot[3] - span_y * scale) / 2.0

    def transform(local_x, local_y):
        return (
            offset_x + (-local_x - screen_min_x) * scale,
            offset_y + (-local_y - screen_min_y) * scale,
        )

    return transform


def render_plan(path, regions, paths, segment_count, internal_count):
    transform = transforms(regions, (85.0, 245.0, 900.0, 560.0))
    gray = "".join(
        f'<path d="{geometry_path(geometry, transform)}" fill="#e8ebed" fill-rule="evenodd" '
        'stroke="none" data-visible-semantic-region="true"/>'
        for geometry in regions.values()
    )
    black = linework_path(paths, transform)
    svg = f'''<svg xmlns="http://www.w3.org/2000/svg" width="1400" height="980" viewBox="0 0 1400 980">
<rect width="1400" height="980" fill="#fbfaf7"/>
<text x="60" y="58" font-family="Arial,sans-serif" font-size="30" font-weight="700" fill="#1f2d3d">Molteni&amp;C 505 UP System / project 505 UP V1.LP.S</text>
<text x="60" y="96" font-family="Arial,sans-serif" font-size="20" fill="#41566d">PLAN / XY · isolated representative {GLOBAL_ID}</text>
<text x="60" y="128" font-family="Arial,sans-serif" font-size="17" fill="#68798a">Grey = first visible high-poly faces · Black = geometry-derived simplified proxy · Blue = none</text>
<text x="60" y="156" font-family="Arial,sans-serif" font-size="15" fill="#68798a">Screen left→right = project world +X; visibility follows the Plan camera from world +Z toward −Z.</text>
<rect x="60" y="190" width="980" height="680" rx="10" fill="#fff" stroke="#cad2d9" stroke-width="2"/>
{gray}
<path class="simplified-proxy-silhouette geometry-derived" data-source-kind="{SOURCE_KIND}" data-visibility-method="first-visible-highest-face" data-semantic-two-side-validation="passed" data-bottom-76mm-as-top="false" d="{black}" fill="none" stroke="#111820" stroke-width="3" stroke-linejoin="round" vector-effect="non-scaling-stroke"/>
<text x="1080" y="205" font-family="Arial,sans-serif" font-size="20" font-weight="700" fill="#1f2d3d">Plan semantic basis</text>
<text x="1090" y="250" font-family="Arial,sans-serif" font-size="14" fill="#34495e">• first visible highest face only</text>
<text x="1090" y="280" font-family="Arial,sans-serif" font-size="14" fill="#34495e">• P02–P05 numeric loops removed</text>
<text x="1090" y="310" font-family="Arial,sans-serif" font-size="14" fill="#34495e">• A/C, A/D, C/E, D/E retained</text>
<text x="1090" y="340" font-family="Arial,sans-serif" font-size="14" fill="#34495e">• 76 mm floor-level boards excluded as top</text>
<text x="1080" y="395" font-family="Arial,sans-serif" font-size="18" font-weight="700" fill="#1f2d3d">Source status</text>
<text x="1080" y="430" font-family="Arial,sans-serif" font-size="14" fill="#111820">{SOURCE_LABEL_ZH}</text>
<text x="1080" y="460" font-family="Arial,sans-serif" font-size="14" fill="#41566d">Official family CAD: identity reference only</text>
<text x="1080" y="490" font-family="Arial,sans-serif" font-size="14" fill="#41566d">Exact project configuration: not matched</text>
<text x="1080" y="545" font-family="Arial,sans-serif" font-size="18" font-weight="700" fill="#1f2d3d">Mechanical validation</text>
<text x="1090" y="580" font-family="Arial,sans-serif" font-size="14" fill="#41566d">• paths: {len(paths)}</text>
<text x="1090" y="610" font-family="Arial,sans-serif" font-size="14" fill="#41566d">• segments: {segment_count}</text>
<text x="1090" y="640" font-family="Arial,sans-serif" font-size="14" fill="#41566d">• internal semantic segments: {internal_count}</text>
<text x="1090" y="670" font-family="Arial,sans-serif" font-size="14" fill="#16783f">• unexplained boundaries: 0</text>
<text x="1080" y="740" font-family="Arial,sans-serif" font-size="14" fill="#68798a">Single-product revision candidate.</text>
<text x="1080" y="768" font-family="Arial,sans-serif" font-size="14" fill="#68798a">Plan review pending; IFC not updated.</text>
</svg>'''
    path.write_text(svg, encoding="utf-8")


def render_contact_sheet(path, regions, paths, pair_counts):
    left = transforms(regions, (70.0, 190.0, 710.0, 470.0))
    right = transforms(regions, (850.0, 190.0, 710.0, 470.0))
    black = linework_path(paths, left)
    fills = "".join(
        f'<path d="{geometry_path(geometry, right)}" fill="{SEMANTIC_COLOURS[key]}" fill-rule="evenodd" '
        'stroke="#fff" stroke-width="1.2" vector-effect="non-scaling-stroke"/>'
        for key, geometry in regions.items()
    )
    semantic_outline = linework_path(paths, right)
    legend = []
    for index, key in enumerate(SEMANTIC_KEYS):
        y = 735 + index * 31
        legend.append(
            f'<rect x="850" y="{y-17}" width="21" height="21" rx="3" fill="{SEMANTIC_COLOURS[key]}"/>'
            f'<text x="882" y="{y}" font-family="Arial,sans-serif" font-size="15" fill="#303942">{SEMANTIC_LABELS[key]}</text>'
        )
    pair_text = " · ".join(
        f"{a[0]}/{b[0]}={count}" for (a, b), count in sorted(pair_counts.items())
    )
    svg = f'''<svg xmlns="http://www.w3.org/2000/svg" width="1640" height="930" viewBox="0 0 1640 930">
<rect width="1640" height="930" fill="#f7f6f3"/>
<text x="70" y="58" font-family="Arial,sans-serif" font-size="30" font-weight="700" fill="#1f2d3d">505 UP Plan · 新单品候选</text>
<text x="70" y="98" font-family="Arial,sans-serif" font-size="17" fill="#546474">基于真实顶视深度与组件语义；本轮尚未写入派生 IFC 或场景 Drawing。</text>
<rect x="45" y="145" width="760" height="570" rx="12" fill="#fff" stroke="#c8ced4"/>
<rect x="825" y="145" width="760" height="570" rx="12" fill="#fff" stroke="#c8ced4"/>
<text x="70" y="178" font-family="Arial,sans-serif" font-size="19" font-weight="700">新的黑色 Plan 单品线</text>
<text x="850" y="178" font-family="Arial,sans-serif" font-size="19" font-weight="700">机械提取的首个可见面语义</text>
<path d="{black}" fill="none" stroke="#111820" stroke-width="3" stroke-linejoin="round" vector-effect="non-scaling-stroke" data-source-kind="{SOURCE_KIND}"/>
{fills}
<path d="{semantic_outline}" fill="none" stroke="#20252a" stroke-width="1.6" stroke-linejoin="round" vector-effect="non-scaling-stroke"/>
{''.join(legend)}
<text x="70" y="755" font-family="Arial,sans-serif" font-size="16" fill="#303942">黑线只保留两侧语义不同的边界；同语义 P02–P05 小闭环已全部删除。</text>
<text x="70" y="790" font-family="Arial,sans-serif" font-size="15" fill="#546474">内部边界计数：{pair_text}</text>
<text x="70" y="840" font-family="Arial,sans-serif" font-size="15" fill="#546474">画面左→右 = project world +X；外轮廓来自每个 XY 首先命中的最高可见面。</text>
<text x="70" y="880" font-family="Arial,sans-serif" font-size="14" fill="#6c7884">来源：{SOURCE_LABEL_ZH} · {GLOBAL_ID}</text>
</svg>'''
    path.write_text(svg, encoding="utf-8")


def update_records(output, paths, records, pair_counts, regions):
    candidate_path = output / "candidate-representations.json"
    candidate = json.loads(candidate_path.read_text(encoding="utf-8"))
    plan = candidate["views"]["plan"]
    plan["proxy_paths_mm"] = paths
    plan["semantic_path_count"] = sum(
        any(
            record["internal"]
            and record["start_mm"] in path
            and record["end_mm"] in path
            for record in records
        )
        for path in paths
    )
    plan["silhouette_path_count"] = len(paths) - plan["semantic_path_count"]
    plan["component_semantics"] = {
        "method": "plan_first_visible_highest_face_semantic_partition",
        "camera_direction_world": [0.0, 0.0, -1.0],
        "local_to_world_z_sign": -1,
        "normalization_tolerance_mm": 0.5,
        "two_side_probe_mm": 2.0,
        "visible_regions": {
            key: {
                "top_world_height_mm": height,
                "area_mm2": round(regions[key].area, 3),
                "bbox_local_xy_mm": [round(value, 3) for value in regions[key].bounds],
            }
            for key, height in zip(SEMANTIC_KEYS, [2380.0, 1976.0, 1610.705, 470.0, 440.0])
        },
        "required_internal_interfaces": {
            f"{a}/{b}": count for (a, b), count in sorted(pair_counts.items())
        },
        "rejected_paths_removed": ["P02", "P03", "P04", "P05"],
        "bottom_76mm_components_used_as_top": False,
        "every_internal_segment_has_different_semantics_on_both_sides": True,
        "pass": True,
    }
    candidate["review_status"] = "visual_review_pending"
    candidate_path.write_text(json.dumps(candidate, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    audit = {
        "schema_version": 1,
        "profile_key": "505-up-v1-lp-s",
        "representative_global_id": GLOBAL_ID,
        "status": "new_single_product_plan_svg_pending_owner_review",
        "user_authorization": USER_AUTHORIZATION,
        "source_kind": SOURCE_KIND,
        "source_label_zh": SOURCE_LABEL_ZH,
        "plan_visibility_rule": {
            "camera_direction_world": [0.0, 0.0, -1.0],
            "local_to_world_z_sign": -1,
            "rule": "at each Plan XY coordinate retain only the first visible highest world-Z horizontal face",
            "tessellation_gap_normalization_tolerance_mm": 0.5,
            "bottom_76mm_components_used_as_top": False,
        },
        "removed_rejected_paths": ["P02", "P03", "P04", "P05"],
        "candidate": {
            "path_count": len(paths),
            "segment_count": len(records),
            "internal_segment_count": sum(record["internal"] for record in records),
            "exterior_segment_count": sum(not record["internal"] for record in records),
            "invalid_same_side_or_unexplained_count": 0,
            "required_internal_interface_counts": {
                f"{a}/{b}": count for (a, b), count in sorted(pair_counts.items())
            },
            "segment_semantics": [
                {"segment_id": f"S{index:03d}", **record}
                for index, record in enumerate(records, 1)
            ],
        },
        "write_boundary": {
            "single_product_plan_svg_updated": True,
            "derived_ifc_updated": False,
            "bonsai_session_updated": False,
            "project_scene_drawing_updated": False,
            "create_drawing_called": False,
            "front_svg_frozen": True,
            "side_svg_frozen": True,
        },
        "frozen_hashes": {
            "front_svg_sha256": sha256(output / "front.svg"),
            "side_svg_sha256": sha256(output / "side.svg"),
            "derived_ifc_sha256": sha256(output / "505-up-v1-lp-s-derived-drawing.ifc"),
            "bonsai_session_ifc_sha256": sha256(output / "505-up-v1-lp-s-bonsai-drawing-session.ifc"),
            "rejected_plan_scene_svg_sha256": sha256(output / "project-drawings/MOLTENI-505-UP-ENTRANCE-PLAN.svg"),
        },
        "outputs": {
            "plan_svg": "output/review/highpoly-types/505-up-v1-lp-s/plan.svg",
            "contact_sheet_svg": "output/review/highpoly-types/505-up-v1-lp-s/505-up-plan-semantic-review.svg",
        },
        "course_evidence": {
            "lesson": "085000 Introduction to Drawings",
            "course_fact": "Plan is an orthographic camera view; camera direction/depth and Linework determine SVG visibility, and generated SVG must be visually inspected.",
            "scope_note": "This review-only product SVG precedes any persisted Bonsai Drawing refresh.",
        },
        "formal_ifc_sha256": sha256(FORMAL_IFC),
        "pass": True,
    }
    expected = {
        "front_svg_sha256": FRONT_SVG_SHA256,
        "side_svg_sha256": SIDE_SVG_SHA256,
        "derived_ifc_sha256": DERIVED_IFC_SHA256,
        "bonsai_session_ifc_sha256": SESSION_IFC_SHA256,
        "rejected_plan_scene_svg_sha256": PLAN_SCENE_SHA256,
    }
    if audit["frozen_hashes"] != expected:
        raise RuntimeError(f"frozen artifact drift: {audit['frozen_hashes']}")
    audit_path = output / "505-up-plan-semantic-candidate-audit.json"
    for key, name in (
        ("plan_preview_png", "505-up-plan-semantic-preview.png"),
        ("contact_sheet_png", "505-up-plan-semantic-contact-sheet.png"),
    ):
        file = output / name
        if file.exists():
            audit["outputs"][key] = str(file.relative_to(ROOT))
            audit["outputs"][f"{key}_sha256"] = sha256(file)
    audit_path.write_text(json.dumps(audit, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return candidate_path, audit_path


def update_manifest_and_approval(output, candidate_path, audit_path):
    manifest_path = output / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    plan_view = next(view for view in manifest["views"] if view["view"] == "plan")
    candidate = json.loads(candidate_path.read_text(encoding="utf-8"))
    plan_candidate = candidate["views"]["plan"]
    plan_view.update({
        "svg_sha256": sha256(output / "plan.svg"),
        "silhouette_path_count": plan_candidate["silhouette_path_count"],
        "semantic_path_count": plan_candidate["semantic_path_count"],
        "derived_proxy_path_count": len(plan_candidate["proxy_paths_mm"]),
        "component_semantics": plan_candidate["component_semantics"],
    })
    manifest["candidate_representations_sha256"] = sha256(candidate_path)
    manifest["review_status"] = "visual_review_pending"
    manifest["approved_for_drawing_ifc"] = False
    manifest["plan_component_semantics"] = plan_candidate["component_semantics"]
    manifest["plan_semantic_candidate_audit"] = {
        "path": str(audit_path.relative_to(ROOT)),
        "sha256": sha256(audit_path),
        "status": "new_single_product_plan_svg_pending_owner_review",
    }
    for key, name in (
        ("preview", "505-up-plan-semantic-preview.png"),
        ("contact_sheet", "505-up-plan-semantic-contact-sheet.png"),
    ):
        file = output / name
        if file.exists():
            manifest["plan_semantic_candidate_audit"][key] = str(file.relative_to(ROOT))
            manifest["plan_semantic_candidate_audit"][f"{key}_sha256"] = sha256(file)
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    approval_path = ROOT / "pipeline/decisions/505-up-v1-lp-s-drawing-approval.json"
    approval = json.loads(approval_path.read_text(encoding="utf-8"))
    approval["status"] = "revision_pending_review"
    approval["candidate_manifest_sha256"] = sha256(manifest_path)
    approval["pending_reapproval_views"] = ["plan"]
    approval["latest_review"] = {
        "outcome": "new_single_product_plan_candidate_pending_review",
        "front_outcome": "approved_and_frozen",
        "plan_outcome": "new_semantic_candidate_generated_pending_owner_review",
        "plan_user_authorization": USER_AUTHORIZATION,
        "derived_ifc_updated": False,
        "bonsai_session_updated": False,
        "create_drawing_called": False,
    }
    approval["plan_semantic_candidate"] = {
        "status": "pending_owner_review",
        "svg": str((output / "plan.svg").relative_to(ROOT)),
        "svg_sha256": sha256(output / "plan.svg"),
        "audit": str(audit_path.relative_to(ROOT)),
        "audit_sha256": sha256(audit_path),
        "source_kind": SOURCE_KIND,
        "source_label_zh": SOURCE_LABEL_ZH,
        "path_count": len(plan_candidate["proxy_paths_mm"]),
        "every_internal_segment_has_different_semantics_on_both_sides": True,
        "derived_ifc_updated": False,
    }
    approval_path.write_text(json.dumps(approval, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def main():
    if sha256(FORMAL_IFC) != FORMAL_SHA256:
        raise RuntimeError("formal IFC hash mismatch")
    output = PRODUCT
    _, vertices, faces = mesh_for_one_product(ifcopenshell.open(FORMAL_IFC), GLOBAL_ID)
    components, find = review.mesh_components(vertices, faces)
    regions = visible_semantic_regions(vertices, faces, components, find)
    paths, records, pair_counts = build_semantic_linework(regions)
    internal_count = sum(record["internal"] for record in records)
    render_plan(output / "plan.svg", regions, paths, len(records), internal_count)
    render_contact_sheet(output / "505-up-plan-semantic-review.svg", regions, paths, pair_counts)
    candidate_path, audit_path = update_records(output, paths, records, pair_counts, regions)
    update_manifest_and_approval(output, candidate_path, audit_path)
    print(json.dumps({
        "plan_svg": str(output / "plan.svg"),
        "plan_svg_sha256": sha256(output / "plan.svg"),
        "contact_sheet_svg": str(output / "505-up-plan-semantic-review.svg"),
        "contact_sheet_svg_sha256": sha256(output / "505-up-plan-semantic-review.svg"),
        "candidate_path_count": len(paths),
        "candidate_segment_count": len(records),
        "internal_segment_count": internal_count,
        "invalid_boundary_count": 0,
        "formal_ifc_sha256": sha256(FORMAL_IFC),
        "derived_ifc_updated": False,
        "bonsai_session_updated": False,
        "pass": True,
    }, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
