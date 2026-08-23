#!/usr/bin/env python3
"""Archive FAU02's existing IFC projection in complete project context."""

from __future__ import annotations

import json
import shutil
from pathlib import Path

from falper_sorgente_linework import ROOT, load_json, relative, sha256, write_json
import render_hima01_project_context as shared


PRODUCT_DIR = ROOT / "output/review/highpoly-types/fau02"
CANDIDATE = PRODUCT_DIR / "candidate-representations.json"
GLOBAL_ID = "36ZX3QPyD7SvlXsDKMP8rY"
SOURCE_KIND = "geometry_derived_simplified_proxy"
SOURCE_LABEL_ZH = "基于原始高模几何生成的简化图纸表达"
GENERATOR = "pipeline/scripts/render_fau02_project_context.py"
CONTEXT_VIEWS = (
    (
        "plan",
        ROOT / "drawings/FFL PLAN.svg",
        "project-context-furniture-plan.svg",
        (12.0, 12.0),
        "actual IFC plan projection retained with basin/floor-plan occlusion",
    ),
    (
        "front",
        ROOT / "drawings/elevations/native/EL-P02-B-PUBLIC-PX.svg",
        "project-context-front-elevation.svg",
        (8.0, 6.0),
        "actual IFC elevation projection retained with wall and furniture occlusion",
    ),
    (
        "side",
        ROOT / "drawings/elevations/native/EL-P01-02-LIVING-SOUTH-NY.svg",
        "project-context-side-elevation.svg",
        (8.0, 6.0),
        "actual IFC elevation projection retained with wall and furniture occlusion",
    ),
)


def main() -> None:
    candidate = load_json(CANDIDATE)
    if (
        candidate.get("source_kind") != SOURCE_KIND
        or candidate.get("official_cad_used") is not False
        or candidate.get("third_party_cad_used") is not False
        or candidate.get("review_status") != "visual_review_pending"
    ):
        raise RuntimeError("FAU02 pending geometry-derived candidate gate failed")
    shared.GLOBAL_ID = GLOBAL_ID
    PRODUCT_DIR.mkdir(parents=True, exist_ok=True)
    records = []
    for view, source, filename, padding, visibility_note in CONTEXT_VIEWS:
        source_content = source.read_text(encoding="utf-8")
        if "IfcWall" not in source_content or "IfcSanitaryTerminal" not in source_content:
            raise RuntimeError(f"FAU02 context source is missing wall or sanitary context: {source}")
        if view != "plan" and "IfcFurniture" not in source_content:
            raise RuntimeError(f"FAU02 elevation context is missing furniture: {source}")
        points = shared.group_points(source)
        target_bbox = shared.bbox(points)
        target = PRODUCT_DIR / filename
        shutil.copyfile(source, target)
        crop = PRODUCT_DIR / filename.replace(".svg", "-review.svg")
        shared.write_review_crop(target, crop, target_bbox, *padding)
        preview = PRODUCT_DIR / filename.replace(".svg", "-review-preview.png")
        shared.write_png_preview(crop, preview)
        records.append({
            "view": view,
            "source": relative(source),
            "output": relative(target),
            "output_sha256": sha256(target),
            "review_crop": relative(crop),
            "review_crop_sha256": sha256(crop),
            "review_preview": relative(preview),
            "review_preview_sha256": sha256(preview),
            "existing_project_projection": {
                "ifc_guid": GLOBAL_ID,
                "geometry_source": "actual_project_IFC_representation",
                "visible_projection_bbox_svg_units": [
                    [round(value, 6) for value in target_bbox[0]],
                    [round(value, 6) for value in target_bbox[1]],
                ],
                "visibility_note": visibility_note,
                "occlusion_retained": True,
            },
        })
    manifest = {
        "schema_version": 1,
        "generator": GENERATOR,
        "source_kind": SOURCE_KIND,
        "source_label_zh": SOURCE_LABEL_ZH,
        "candidate_source_kind": SOURCE_KIND,
        "candidate_source_label_zh": SOURCE_LABEL_ZH,
        "official_cad_acquired": True,
        "official_cad_exact_project_configuration_match": False,
        "official_cad_used": False,
        "third_party_cad_used": False,
        "blue_line_present": False,
        "project_context_retained": True,
        "walls_and_surrounding_project_elements_retained": True,
        "furniture_retained_in_elevations": True,
        "existing_actual_ifc_projection_retained": True,
        "overlay_top_layer_with_white_mask": False,
        "top_layer_substitution_performed": False,
        "occlusion_retained": True,
        "reason_no_top_layer_substitution": "GH2 failed the exact-project match gate; duplicating a full silhouette on top would erase valid basin, wall and furniture occlusion before human approval",
        "separate_candidate_views": {
            view: relative(PRODUCT_DIR / f"{view}.svg") for view in ("plan", "front", "side")
        },
        "views": records,
        "pass": True,
    }
    target = PRODUCT_DIR / "project-context-manifest.json"
    write_json(target, manifest)
    print(json.dumps({"manifest": relative(target), "outputs": [item["output"] for item in records], "pass": True}, indent=2))


if __name__ == "__main__":
    main()
