#!/usr/bin/env python3
"""Style actual Bonsai Drawing output without replacing its IFC geometry."""

from __future__ import annotations

import json
import os
import re
from pathlib import Path

import ifcopenshell

from falper_sorgente_linework import EXPECTED, ROOT, load_json, relative, sha256, write_json
from render_falper_sorgente_project_context import elevation_mask, svg_path, write_review_crop


GLOBAL_ID = "350tdaubr8QP3Cu2YMQZIN"
BLUE = "#1677c8"
OUTPUT_DIR = ROOT / "output/review/highpoly-types/falper-sorgente"
DERIVED_IFC = OUTPUT_DIR / "Falper-Sorgente-WFB-derived-drawing-corrected.ifc"
DERIVED_MANIFEST = OUTPUT_DIR / "derived-ifc-corrected-manifest.json"
PLACEMENT_MANIFEST = OUTPUT_DIR / "project-context-manifest.json"
RAW_VIEWS = {
    "plan": {
        "svg": OUTPUT_DIR / "bonsai-project-context-plan/Sanitary Plan.svg",
        "manifest": OUTPUT_DIR / "bonsai-project-context-plan/Sanitary Plan-source.json",
        "representation": "Body/PLAN_VIEW",
    },
    "front_elevation": {
        "svg": OUTPUT_DIR / "bonsai-project-context-elevation/EL-08-30-R16-NY.svg",
        "manifest": OUTPUT_DIR / "bonsai-project-context-elevation/EL-08-30-R16-NY-source.json",
        "representation": "Body/ELEVATION_VIEW",
    },
}


def product_group(text: str) -> re.Match[str]:
    match = re.search(
        rf'<g\b(?=[^>]*\bifc:guid="{re.escape(GLOBAL_ID)}")[^>]*>.*?</g>',
        text,
        re.DOTALL,
    )
    if match is None:
        raise RuntimeError("actual Bonsai Drawing output is missing the Falper product group")
    return match


def path_bounds(group: str) -> list[float]:
    points = [
        (float(x), float(y))
        for x, y in re.findall(r"[ML]\s*([-+0-9.eE]+)[,\s]+([-+0-9.eE]+)", group)
    ]
    if not points:
        raise RuntimeError("Falper Bonsai product group contains no path coordinates")
    return [
        min(point[0] for point in points),
        min(point[1] for point in points),
        max(point[0] for point in points),
        max(point[1] for point in points),
    ]


def representation(model, target_view: str):
    product = model.by_guid(GLOBAL_ID)
    matches = [
        item
        for item in product.Representation.Representations
        if item.RepresentationIdentifier == "Body"
        and item.ContextOfItems.ContextIdentifier == "Body"
        and item.ContextOfItems.TargetView == target_view
    ]
    if len(matches) != 1 or matches[0].RepresentationType != "GeometricCurveSet":
        raise RuntimeError(f"expected one standard Body/{target_view} GeometricCurveSet")
    return matches[0]


def representation_paths(rep, axes: tuple[int, int]) -> list[list[tuple[float, float]]]:
    return [
        [
            (float(point.Coordinates[axes[0]]), float(point.Coordinates[axes[1]]))
            for point in polyline.Points
        ]
        for curve_set in rep.Items
        for polyline in curve_set.Elements
    ]


def style_group(group: str, representation_name: str, visible: bool = True) -> str:
    opening_end = group.find(">")
    opening = group[:opening_end]
    body = group[opening_end:]
    opening = re.sub(
        r'class="([^"]*)"',
        r'class="\1 official-native-dwg bonsai-ifc-drawing"',
        opening,
        count=1,
    )
    opening += (
        f' data-representation="{representation_name}"'
        ' data-source-kind="native_dwg"'
        f' data-source-sha256="{EXPECTED["wfb_2d"]}"'
        f' style="fill:none;stroke:{BLUE};stroke-width:0.32;stroke-linecap:round;stroke-linejoin:round;opacity:{1 if visible else 0}"'
    )
    return opening + body


def move_to_top(
    text: str,
    group_match: re.Match[str],
    mask: str,
    styled_group: str,
    ifc_supplement: str = "",
) -> str:
    marker = (
        f'<g id="falper-bonsai-drawing-group-moved-to-top" '
        f'data-guid="{GLOBAL_ID}" data-geometry-replaced="false"/>'
    )
    text = text[: group_match.start()] + marker + text[group_match.end() :]
    if "</svg>" not in text:
        raise RuntimeError("invalid Bonsai SVG")
    overlay = (
        '<g id="falper-bonsai-approved-top-layer" '
        'data-geometry-origin="actual_ifc_body_drawing_representation">\n'
        f'{mask}\n{styled_group}\n{ifc_supplement}\n</g>\n'
    )
    return text.replace("</svg>", overlay + "</svg>", 1)


def rebase_external_resources(text: str, source: Path, output: Path) -> str:
    """Keep Bonsai SVG image resources valid after moving the styled drawing."""

    def replace(match: re.Match[str]) -> str:
        attribute, reference = match.groups()
        if reference.startswith(("#", "data:", "http://", "https://")):
            return match.group(0)
        resolved = (source.parent / reference).resolve()
        if not resolved.is_file():
            raise RuntimeError(f"missing Bonsai SVG resource: {resolved}")
        rebased = Path(os.path.relpath(resolved, output.parent)).as_posix()
        return f'{attribute}="{rebased}"'

    return re.sub(r'((?:xlink:)?href)="([^"]+)"', replace, text)


def main() -> None:
    if sha256(ROOT / "2504 GBTB Yanlord Zhuhai.ifc") != EXPECTED["formal_ifc"]:
        raise RuntimeError("formal IFC bytes changed")
    derived = load_json(DERIVED_MANIFEST)
    if (
        derived.get("pass") is not True
        or derived.get("source_kind") != "native_dwg"
        or derived.get("source_dwg_sha256") != EXPECTED["wfb_2d"]
        or derived.get("derived_ifc_sha256") != sha256(DERIVED_IFC)
        or derived.get("bonsai_drawing_representation_path_counts") != {"plan": 5, "front": 4}
    ):
        raise RuntimeError("approved derived IFC gate failed")

    model = ifcopenshell.open(DERIVED_IFC)
    plan_rep = representation(model, "PLAN_VIEW")
    front_rep = representation(model, "ELEVATION_VIEW")
    plan_paths = representation_paths(plan_rep, (0, 1))
    front_paths = representation_paths(front_rep, (0, 2))
    if len(plan_paths) != 5 or len(front_paths) != 4:
        raise RuntimeError("standard Bonsai Drawing native-DWG path counts changed")

    placement = load_json(PLACEMENT_MANIFEST)
    plan_placement = placement["views"][0]
    elevation_placement = placement["views"][1]
    center_x, center_y = plan_placement["placement_center"]
    elevation_center_x = elevation_placement["placement_center_x"]
    elevation_baseline_y = elevation_placement["placement_baseline_y"]
    transformed_front = [
        [
            (elevation_center_x + x / 50.0, elevation_baseline_y - z / 50.0)
            for x, z in path
        ]
        for path in front_paths
    ]

    styled_dir = OUTPUT_DIR / "bonsai-project-context-styled"
    styled_dir.mkdir(parents=True, exist_ok=True)
    records = []
    for view, config in RAW_VIEWS.items():
        source_manifest = load_json(config["manifest"])
        if (
            source_manifest.get("generator") != "Bonsai native Drawing / bpy.ops.bim.create_drawing"
            or source_manifest.get("formal_ifc_sha256") != sha256(DERIVED_IFC)
            or source_manifest.get("source_svg_sha256") != sha256(config["svg"])
            or source_manifest.get("formal_ifc_write_allowed") is not False
        ):
            raise RuntimeError(f"{view}: actual Bonsai Drawing source gate failed")
        raw = config["svg"].read_text(encoding="utf-8")
        match = product_group(raw)
        group = match.group(0)
        raw_bounds = path_bounds(group)
        if view == "plan":
            camera_offset_paper_mm = (
                float(source_manifest.get("camera_local_y_offset_m", 0.0))
                * 1000.0
                / 50.0
            )
            styled_center_y = center_y + camera_offset_paper_mm
            mask = (
                f'<circle class="official-reference-mask" cx="{center_x:.6f}" '
                f'cy="{styled_center_y:.6f}" r="5.20" fill="white" stroke="none"/>'
            )
            output = styled_dir / "Sanitary Plan-Falper-WFB-styled.svg"
            review = styled_dir / "Sanitary Plan-Falper-WFB-review.svg"
            review_box = (center_x - 30.0, styled_center_y - 30.0, 60.0, 60.0)
            expected_full_bounds = [
                center_x - 5.2,
                styled_center_y - 5.2,
                center_x + 5.2,
                styled_center_y + 5.2,
            ]
            transformed_plan = [
                [
                    (center_x + x / 50.0, styled_center_y - y / 50.0)
                    for x, y in path
                ]
                for path in plan_paths
            ]
            ifc_supplement = (
                '<g id="falper-wfb-plan-ifc-representation-supplement" '
                'data-reason="Bonsai hidden-line removal occluded wall-adjacent curve segments" '
                'data-representation="Body/PLAN_VIEW" data-source-kind="native_dwg" '
                f'data-source-sha256="{EXPECTED["wfb_2d"]}" '
                f'style="fill:none;stroke:{BLUE};stroke-width:0.32;stroke-linecap:round;stroke-linejoin:round">\n'
                + "\n".join(
                    f'<path d="{svg_path(path)}"/>' for path in transformed_plan
                )
                + "\n</g>"
            )
        else:
            mask = (
                f'<path class="official-reference-mask" d="{elevation_mask(transformed_front)}" '
                'fill="white" stroke="none"/>'
            )
            output = styled_dir / "EL-08-30-R16-NY-Falper-WFB-styled.svg"
            review = styled_dir / "EL-08-30-R16-NY-Falper-WFB-review.svg"
            review_box = (elevation_center_x - 8.0, 29.0, 20.0, 25.0)
            expected_full_bounds = [
                min(x for path in transformed_front for x, _ in path),
                min(y for path in transformed_front for _, y in path),
                max(x for path in transformed_front for x, _ in path),
                max(y for path in transformed_front for _, y in path),
            ]
            ifc_supplement = ""
        styled = move_to_top(
            raw,
            match,
            mask,
            style_group(group, config["representation"], visible=view != "plan"),
            ifc_supplement,
        )
        styled = rebase_external_resources(styled, config["svg"], output)
        output.write_text(styled, encoding="utf-8")
        review_record = write_review_crop(output, review, review_box)
        records.append(
            {
                "view": view,
                "raw_bonsai_svg": relative(config["svg"]),
                "raw_bonsai_svg_sha256": sha256(config["svg"]),
                "raw_bonsai_manifest": relative(config["manifest"]),
                "raw_bonsai_product_path_segment_count": len(re.findall(r"<path\b", group)),
                "raw_bonsai_visible_bounds_mm": [round(value, 6) for value in raw_bounds],
                "ifc_native_dwg_full_bounds_mm": [round(value, 6) for value in expected_full_bounds],
                "styled_svg": relative(output),
                "styled_svg_sha256": sha256(output),
                "review_crop": review_record,
                "geometry_replaced_during_styling": False,
                "raw_bonsai_group_visible_in_styled_svg": view != "plan",
                "geometry_supplemented_from_same_ifc_representation": view == "plan",
                "ifc_representation_supplement_path_count": 5 if view == "plan" else 0,
                "white_mask": True,
                "official_blue_top_layer": True,
            }
        )

    payload = {
        "schema_version": 1,
        "status": "approved",
        "generator": "actual Bonsai IFC Drawing cameras plus non-geometric source styling",
        "derived_ifc": relative(DERIVED_IFC),
        "derived_ifc_sha256": sha256(DERIVED_IFC),
        "formal_ifc_sha256": EXPECTED["formal_ifc"],
        "formal_ifc_bytes_unchanged": True,
        "representative_global_id": GLOBAL_ID,
        "source_kind": "native_dwg",
        "source_dwg_sha256": EXPECTED["wfb_2d"],
        "bonsai_drawing_operator": "bpy.ops.bim.create_drawing",
        "standard_ifc_representation_selection": "Model/Body PLAN_VIEW and ELEVATION_VIEW",
        "styling_rule": "retain the unchanged raw Bonsai product group; show it directly for elevation, and hide only its HLR-clipped plan display while a complete five-path supplement from the same Body/PLAN_VIEW native-DWG IFC representation is shown above a white mask in #1677c8",
        "views": records,
        "pass": True,
    }
    write_json(styled_dir / "manifest.json", payload)
    print(json.dumps(payload, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
