#!/usr/bin/env python3
"""Place exact Miami Soft E09 native-DWG linework in project review context."""

from __future__ import annotations

import json
import re
import xml.etree.ElementTree as ET
from pathlib import Path

import ifcopenshell
import numpy as np
from ifcopenshell.util.placement import get_local_placement

from falper_sorgente_linework import ROOT, load_json, relative, sha256, write_json
from render_hima01_project_context import rebase_external_resources
from render_sis04_project_context import (
    add_white_review_background,
    suppress_plan_elevation_markers,
    write_uncached_png_preview,
)


PRODUCT_DIR = ROOT / "output/review/highpoly-types/miamisoft-e09"
CANDIDATE_PATH = PRODUCT_DIR / "candidate-representations.json"
MANIFEST_PATH = PRODUCT_DIR / "manifest.json"
REFERENCE_PATH = PRODUCT_DIR / "official-dwg-review-reference.json"
APPROVAL_PATH = ROOT / "pipeline/decisions/miamisoft-e09-drawing-approval.json"
INDEX_PATH = PRODUCT_DIR / "index.html"
FORMAL_IFC = ROOT / "2504 GBTB Yanlord Zhuhai.ifc"
FORMAL_IFC_SHA256 = "7a50b87e8f48a7c2bfbaa9f1dabdfca7155fa684b2325c66aa1ae15c8ab0a25c"
SOURCE_DWG_SHA256 = "efc3315feed0a2dac8938aebbd7336856dabba6721395e4520783c39ba0cc51a"
REPRESENTATIVE_GLOBAL_ID = "3osoWufdD1mhDDAM6lcix4"
PLAN_DRAWING_GLOBAL_ID = "3mtEdbDYn1d9iNK0Iw113w"
SIDE_DRAWING_GLOBAL_ID = "3oQZvKp$DBYQZCNExowkE1"
BLUE = "#1677c8"
PAPER_SCALE_SVG_UNITS_PER_MM = 0.02
NUMBER = re.compile(r"[-+]?(?:\d*\.\d+|\d+)(?:[eE][-+]?\d+)?")

VIEWS = (
    {
        "view": "plan",
        "candidate_view": "plan",
        "source": ROOT / "drawings/Furniture Plan.svg",
        "output": PRODUCT_DIR / "project-context-furniture-plan.svg",
        "review": PRODUCT_DIR / "project-context-furniture-plan-review.svg",
        "preview": PRODUCT_DIR / "project-context-furniture-plan-review-preview.png",
        "drawing_global_id": PLAN_DRAWING_GLOBAL_ID,
        "local_axes": (0, 1),
        "padding": (12.0, 12.0),
        "expected_path_count": 64,
        "recorded_view_reflection": False,
    },
    {
        "view": "side",
        "candidate_view": "side",
        "source": ROOT / "drawings/elevations/native/EL-01-02-R20-PY.svg",
        "output": PRODUCT_DIR / "project-context-r20-side-elevation.svg",
        "review": PRODUCT_DIR / "project-context-r20-side-elevation-review.svg",
        "preview": PRODUCT_DIR / "project-context-r20-side-elevation-review-preview.png",
        "drawing_global_id": SIDE_DRAWING_GLOBAL_ID,
        "local_axes": (1, 2),
        "padding": (8.0, 6.0),
        "expected_path_count": 85,
        "recorded_view_reflection": True,
    },
)


def canonical_sha256(value: object) -> str:
    import hashlib

    payload = json.dumps(value, ensure_ascii=False, separators=(",", ":"), sort_keys=True)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def guid(element: ET.Element) -> str | None:
    return next((value for key, value in element.attrib.items() if key.endswith("guid")), None)


def path_points(element: ET.Element) -> list[tuple[float, float]]:
    values = [float(value) for value in NUMBER.findall(element.attrib.get("d", ""))]
    return list(zip(values[0::2], values[1::2]))


def projection_bbox(source: Path) -> tuple[tuple[float, float], tuple[float, float]]:
    groups: list[list[tuple[float, float]]] = []
    for element in ET.parse(source).getroot().iter():
        if not element.tag.endswith("g") or guid(element) != REPRESENTATIVE_GLOBAL_ID:
            continue
        if "projection" not in element.attrib.get("class", "").split():
            continue
        points = [
            point
            for child in element.iter()
            if child.tag.endswith("path")
            for point in path_points(child)
        ]
        if points:
            groups.append(points)
    if len(groups) != 1:
        raise RuntimeError(f"expected one complete target projection in {source}, found {len(groups)}")
    points = groups[0]
    return (
        (min(point[0] for point in points), min(point[1] for point in points)),
        (max(point[0] for point in points), max(point[1] for point in points)),
    )


def svg_viewbox(source: Path) -> tuple[float, float, float, float]:
    match = re.search(r'viewBox="([^"]+)"', source.read_text(encoding="utf-8"))
    if match is None:
        raise RuntimeError(f"missing viewBox in {source}")
    values = tuple(float(value) for value in match.group(1).split())
    if len(values) != 4:
        raise RuntimeError(f"invalid viewBox in {source}")
    return values  # type: ignore[return-value]


def point_to_svg(
    point: list[float],
    axes: tuple[int, int],
    object_matrix: np.ndarray,
    camera_inverse: np.ndarray,
    viewbox: tuple[float, float, float, float],
) -> tuple[float, float]:
    local = np.array([0.0, 0.0, 0.0, 1.0])
    local[axes[0]] = point[0]
    local[axes[1]] = point[1]
    camera = camera_inverse @ object_matrix @ local
    centre_x = viewbox[0] + viewbox[2] / 2.0
    centre_y = viewbox[1] + viewbox[3] / 2.0
    return (
        centre_x + float(camera[0]) * PAPER_SCALE_SVG_UNITS_PER_MM,
        centre_y - float(camera[1]) * PAPER_SCALE_SVG_UNITS_PER_MM,
    )


def official_svg_paths(
    paths: list[list[list[float]]],
    axes: tuple[int, int],
    object_matrix: np.ndarray,
    camera_inverse: np.ndarray,
    viewbox: tuple[float, float, float, float],
) -> tuple[str, tuple[tuple[float, float], tuple[float, float]]]:
    elements: list[str] = []
    all_points: list[tuple[float, float]] = []
    for index, path in enumerate(paths):
        transformed = [point_to_svg(point, axes, object_matrix, camera_inverse, viewbox) for point in path]
        if len(transformed) < 2:
            continue
        all_points.extend(transformed)
        command = f"M {transformed[0][0]:.6f},{transformed[0][1]:.6f} " + " ".join(
            f"L {point[0]:.6f},{point[1]:.6f}" for point in transformed[1:]
        )
        elements.append(
            f'<path class="official-reference native-dwg" data-official-source-path-index="{index}" d="{command}"/>'
        )
    if len(elements) != len(paths):
        raise RuntimeError("official source contains a degenerate path")
    bounds = (
        (min(point[0] for point in all_points), min(point[1] for point in all_points)),
        (max(point[0] for point in all_points), max(point[1] for point in all_points)),
    )
    return "\n  ".join(elements), bounds


def mark_actual_ifc_body(content: str) -> str:
    pattern = rf'(<g\b(?=[^>]*\bifc:guid="{re.escape(REPRESENTATIVE_GLOBAL_ID)}")[^>]*)(>)'
    content, count = re.subn(pattern, r'\1 data-context-role="actual-ifc-body"\2', content)
    if count < 1:
        raise RuntimeError("target IFC Body group not found in project drawing")
    return content


def write_context_svg(source: Path, target: Path, view: str, path_elements: str) -> None:
    content = rebase_external_resources(source.read_text(encoding="utf-8"), source, target)
    content = mark_actual_ifc_body(content)
    style = f"""<style id="miamisoft-e09-project-context-review-style">
svg &gt; *:not(defs):not(style):not(#miamisoft-e09-{view}-official-overlay) {{ opacity:0.34; filter:grayscale(1); }}
g[data-context-role="actual-ifc-body"] {{ opacity:0.10 !important; }}
#miamisoft-e09-{view}-official-overlay {{ opacity:1 !important; }}
#miamisoft-e09-{view}-official-overlay path {{ fill:none !important; stroke:{BLUE} !important; stroke-width:0.075; stroke-linecap:round; stroke-linejoin:round; }}
</style>"""
    overlay = f"""<g id="miamisoft-e09-{view}-official-overlay" class="official-reference native-dwg project-context-review-overlay" data-source-kind="native_dwg" data-source-dwg-sha256="{SOURCE_DWG_SHA256}" data-ifc-guid="{REPRESENTATIVE_GLOBAL_ID}" data-review-only="true" data-bonsai-scene-svg="false" data-ifc-write-performed="false" data-uniform-geometry-scale="1.0" fill="none">
  {path_elements}
</g>"""
    if "</svg>" not in content:
        raise RuntimeError(f"invalid project SVG: {source}")
    content = content.rsplit("</svg>", 1)[0] + style + "\n" + overlay + "\n</svg>\n"
    target.write_text(content, encoding="utf-8")


def write_review_crop(
    source: Path,
    target: Path,
    bounds: tuple[tuple[float, float], tuple[float, float]],
    padding: tuple[float, float],
) -> None:
    minimum, maximum = bounds
    x = minimum[0] - padding[0]
    y = minimum[1] - padding[1]
    width = maximum[0] - minimum[0] + 2.0 * padding[0]
    height = maximum[1] - minimum[1] + 2.0 * padding[1]
    content = source.read_text(encoding="utf-8")
    content, count = re.subn(
        r'viewBox="[^"]+"',
        f'viewBox="{x:.6f} {y:.6f} {width:.6f} {height:.6f}"',
        content,
        count=1,
    )
    if count != 1:
        raise RuntimeError(f"could not crop {source}")
    target.write_text(content, encoding="utf-8")


def rounded_matrix(matrix: np.ndarray) -> list[list[float]]:
    return [[round(float(value), 9) for value in row] for row in matrix]


def within_viewbox(
    bounds: tuple[tuple[float, float], tuple[float, float]],
    viewbox: tuple[float, float, float, float],
) -> bool:
    minimum, maximum = bounds
    return (
        minimum[0] >= viewbox[0]
        and minimum[1] >= viewbox[1]
        and maximum[0] <= viewbox[0] + viewbox[2]
        and maximum[1] <= viewbox[1] + viewbox[3]
    )


def update_index() -> None:
    content = INDEX_PATH.read_text(encoding="utf-8")
    content = content.replace(
        "<h2>Project drawing context</h2>",
        "<h2>Project-context review overlays · official native DWG blue line</h2>",
        1,
    )
    note = (
        '<p class="project-context-review-note"><strong>Blue</strong> is the exact official E09 native DWG placed at the representative instance world transform; '
        'project drawing context is light grey and the current actual IFC Body is lighter grey. These two files are review-only project-context overlays, '
        'not Bonsai scene Drawing SVGs, and no IFC was written.</p>'
    )
    heading = "<h2>Project-context review overlays · official native DWG blue line</h2>"
    if note not in content:
        content = content.replace(heading, heading + note, 1)
    INDEX_PATH.write_text(content, encoding="utf-8")


def main() -> None:
    if sha256(FORMAL_IFC) != FORMAL_IFC_SHA256:
        raise RuntimeError("formal IFC baseline hash changed")
    candidate = load_json(CANDIDATE_PATH)
    reference = load_json(REFERENCE_PATH)
    manifest = load_json(MANIFEST_PATH)
    approval = load_json(APPROVAL_PATH)
    if candidate.get("source_kind") != "native_dwg" or candidate.get("official_cad_used") is not True:
        raise RuntimeError("Miami candidate is not the verified official native DWG")
    if reference.get("source_dwg_sha256") != SOURCE_DWG_SHA256:
        raise RuntimeError("official DWG source hash changed")
    if approval.get("status") != "pending" or approval.get("derived_ifc_write_allowed") is not False:
        raise RuntimeError("review write gate is not closed")

    model = ifcopenshell.open(str(FORMAL_IFC))
    product = model.by_guid(REPRESENTATIVE_GLOBAL_ID)
    object_matrix = np.asarray(get_local_placement(product.ObjectPlacement), dtype=float)
    records = []
    for spec in VIEWS:
        view = spec["view"]
        source = spec["source"]
        viewbox = svg_viewbox(source)
        drawing = model.by_guid(spec["drawing_global_id"])
        camera_matrix = np.asarray(get_local_placement(drawing.ObjectPlacement), dtype=float)
        camera_inverse = np.linalg.inv(camera_matrix)
        candidate_view = candidate["views"][spec["candidate_view"]]
        paths = candidate_view["official_cad_paths_mm"]
        alignment = candidate_view["alignment"]
        if len(paths) != spec["expected_path_count"]:
            raise RuntimeError(f"{view} official path count changed")
        if alignment.get("uniform_scale") != 1.0 or alignment.get("source_geometry_deformed") is not False:
            raise RuntimeError(f"{view} official geometry is not 1:1")
        if alignment.get("view_direction_reflection_x") is not spec["recorded_view_reflection"]:
            raise RuntimeError(f"{view} recorded view reflection changed")
        elements, official_bounds = official_svg_paths(
            paths,
            spec["local_axes"],
            object_matrix,
            camera_inverse,
            viewbox,
        )
        if not within_viewbox(official_bounds, viewbox):
            raise RuntimeError(f"{view} official blue line is clipped by the project viewBox")
        target_bbox = projection_bbox(source)
        crop_bounds = (
            (min(target_bbox[0][0], official_bounds[0][0]), min(target_bbox[0][1], official_bounds[0][1])),
            (max(target_bbox[1][0], official_bounds[1][0]), max(target_bbox[1][1], official_bounds[1][1])),
        )
        write_context_svg(source, spec["output"], view, elements)
        write_review_crop(spec["output"], spec["review"], crop_bounds, spec["padding"])
        if view == "plan":
            suppress_plan_elevation_markers(spec["review"])
        add_white_review_background(spec["review"])
        write_uncached_png_preview(spec["review"], spec["preview"])
        records.append(
            {
                "view": view,
                "candidate_view": spec["candidate_view"],
                "source": relative(source),
                "source_sha256": sha256(source),
                "output": relative(spec["output"]),
                "output_sha256": sha256(spec["output"]),
                "review_crop": relative(spec["review"]),
                "review_crop_sha256": sha256(spec["review"]),
                "review_preview": relative(spec["preview"]),
                "review_preview_sha256": sha256(spec["preview"]),
                "overlay": {
                    "ifc_guid": REPRESENTATIVE_GLOBAL_ID,
                    "source_kind": "native_dwg",
                    "official_cad_used": True,
                    "third_party_cad_used": False,
                    "source_dwg_sha256": SOURCE_DWG_SHA256,
                    "source_path_count": len(paths),
                    "source_path_semantic_sha256": canonical_sha256(paths),
                    "blue_line_present": True,
                    "stroke": BLUE,
                    "uniform_geometry_scale": 1.0,
                    "recorded_view_direction_reflection_x": spec["recorded_view_reflection"],
                    "reflection_already_applied_in_candidate": spec["recorded_view_reflection"],
                    "additional_reflection_applied_in_project_context": False,
                    "alignment_by_bbox_fit": False,
                    "rigid_world_placement_only": True,
                    "paper_scale_svg_units_per_mm": PAPER_SCALE_SVG_UNITS_PER_MM,
                    "local_axes": list(spec["local_axes"]),
                    "object_placement_matrix_mm": rounded_matrix(object_matrix),
                    "drawing_camera_placement_matrix_mm": rounded_matrix(camera_matrix),
                    "drawing_camera_inverse_matrix": rounded_matrix(camera_inverse),
                    "project_viewbox": [round(value, 9) for value in viewbox],
                    "project_svg_bounds": [[round(value, 6) for value in pair] for pair in official_bounds],
                    "within_project_viewbox": True,
                    "no_clipping": True,
                },
                "project_context_style": {
                    "surrounding_project_elements": "light_grey",
                    "actual_ifc_body": "lighter_grey",
                    "official_native_dwg": "blue",
                },
            }
        )

    context_manifest = {
        "schema_version": 2,
        "generator": "pipeline/scripts/render_miamisoft_e09_project_context.py",
        "artifact_kind": "review_only_project_context_overlay",
        "not_bonsai_scene_drawing_svg": True,
        "ifc_write_performed": False,
        "derived_ifc_write_allowed": False,
        "formal_ifc_write_allowed": False,
        "formal_ifc_sha256": FORMAL_IFC_SHA256,
        "formal_ifc_bytes_unchanged": sha256(FORMAL_IFC) == FORMAL_IFC_SHA256,
        "source_kind": "native_dwg",
        "source_label_zh": "Baxter 官方 Miami Soft E09 原生二维 DWG 项目上下文审核蓝线",
        "official_cad_used": True,
        "third_party_cad_used": False,
        "blue_line_present": True,
        "project_context_retained": True,
        "walls_and_surrounding_project_elements_retained": True,
        "overlay_top_layer_with_white_mask": False,
        "geometry_scale": 1.0,
        "alignment_by_bbox_fit": False,
        "context_view_scope": {
            "included": ["plan", "side"],
            "excluded": ["front"],
            "reason": "Only Furniture Plan and the complete R20 +Y side projection are available as mechanically complete project contexts; no front context is invented.",
        },
        "provider_pre_state": {
            "provider": "bonsai_mcp",
            "online": True,
            "ifc_schema": "IFC4",
            "project_name": "My Project",
            "query_only": True,
            "provider_mutation_used": False,
        },
        "course_evidence": {
            "lesson": "085000 Underlay / Linework / Annotation",
            "source": "research/bonsai-course/lessons/085000/README.md (indexed course evidence)",
            "fact": "The drawing workflow separates Underlay, Linework and Annotation and requires visual inspection after SVG generation.",
            "current_version_adaptation": "Review-only analogue: project context is light-grey underlay, current Body is lighter-grey comparison, and official native DWG is blue linework; generated SVG and PNG are visually inspected.",
        },
        "review_annotation_suppression": {
            "plan": "official-elevation-anchor groups only in review crop",
            "side": "none",
            "full_project_svgs_keep_original_annotations": True,
        },
        "views": records,
        "pass": all(record["overlay"]["no_clipping"] for record in records),
    }
    context_manifest_path = PRODUCT_DIR / "project-context-manifest.json"
    write_json(context_manifest_path, context_manifest)

    reference["project_context_official_dwg_overlay"] = {
        "manifest": relative(context_manifest_path),
        "manifest_sha256": sha256(context_manifest_path),
        "views": ["plan", "side"],
        "source_kind": "native_dwg",
        "blue_line_present": True,
        "uniform_geometry_scale": 1.0,
        "alignment_by_bbox_fit": False,
        "representative_global_id": REPRESENTATIVE_GLOBAL_ID,
        "artifact_kind": "review_only_project_context_overlay",
        "not_bonsai_scene_drawing_svg": True,
        "ifc_write_performed": False,
        "pass": True,
    }
    write_json(REFERENCE_PATH, reference)
    candidate["official_dwg_review_reference_sha256"] = sha256(REFERENCE_PATH)
    candidate["derived_ifc_write_allowed"] = False
    candidate["formal_ifc_write_allowed"] = False
    candidate["review_status"] = "visual_review_pending"
    write_json(CANDIDATE_PATH, candidate)

    manifest["official_reference"]["mechanical_selection_record_sha256"] = sha256(REFERENCE_PATH)
    manifest["candidate_representations_sha256"] = sha256(CANDIDATE_PATH)
    manifest["project_context"] = {
        "manifest": relative(context_manifest_path),
        "manifest_sha256": sha256(context_manifest_path),
        "artifact_kind": "review_only_project_context_overlay",
        "not_bonsai_scene_drawing_svg": True,
        "ifc_write_performed": False,
        "source_kind": "native_dwg",
        "official_cad_used": True,
        "blue_line_present": True,
        "uniform_geometry_scale": 1.0,
        "alignment_by_bbox_fit": False,
        "walls_and_surrounding_project_elements_retained": True,
        "pass": True,
    }
    manifest["review_status"] = "visual_review_pending"
    manifest["approved_for_drawing_ifc"] = False
    manifest["derived_ifc_write_allowed"] = False
    manifest["formal_ifc_write"] = "not performed"
    manifest["formal_ifc_sha256_after_generation"] = sha256(FORMAL_IFC)
    manifest["formal_ifc_bytes_unchanged"] = sha256(FORMAL_IFC) == FORMAL_IFC_SHA256
    write_json(MANIFEST_PATH, manifest)

    approval["candidate_manifest_sha256"] = sha256(MANIFEST_PATH)
    approval["status"] = "pending"
    approval["approved_views"] = []
    approval["derived_ifc_write_allowed"] = False
    approval["formal_authoritative_ifc_write_allowed"] = False
    write_json(APPROVAL_PATH, approval)
    update_index()

    if sha256(FORMAL_IFC) != FORMAL_IFC_SHA256:
        raise RuntimeError("formal IFC changed during review overlay generation")
    print(
        json.dumps(
            {
                "manifest": relative(context_manifest_path),
                "manifest_sha256": sha256(context_manifest_path),
                "views": [
                    {
                        "view": record["view"],
                        "path_count": record["overlay"]["source_path_count"],
                        "output": record["output"],
                        "review": record["review_crop"],
                        "preview": record["review_preview"],
                    }
                    for record in records
                ],
                "formal_ifc_sha256": sha256(FORMAL_IFC),
                "ifc_write_performed": False,
                "pass": context_manifest["pass"],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
