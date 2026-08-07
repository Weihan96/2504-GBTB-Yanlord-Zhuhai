#!/usr/bin/env python3
"""Extract the A-102 demolition-wall review boundary without writing IFC.

The source PDF is a plotted demolition markup over the handed-over apartment
plan.  Its coloured vector rectangles are converted to project world
coordinates through explicit grid anchors.  The output remains a review
register: every inferred wall requires human confirmation before any
IfcWall with Status=DEMOLISH may be created.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import uuid
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Sequence


EXPECTED_PDF_SHA256 = "84a18d211fe8e6c3e05f56b35027a96799d7c54cd4bc75cd602c72f65ea327f0"
EXPECTED_DWG_SHA256 = "ba355f6a90732ad07f843d59e8bab5e1da9daffe7aed74889d1d265b2ce22d7e"
PROTECTED_WALL_GUID = "0hKdvAZkn1TejLgJhK_vDp"
PROTECTED_OPENING_GUID = "1YxMx6s0r3ZPPohkRKXWbl"
GUID_NAMESPACE = uuid.UUID("7d5fc9ce-f4c1-4adf-8a34-f7ef0dd9a102")

# Rect indices and page-space bounds are frozen evidence extracted from the
# hash-guarded PDF.  If pdfplumber is available, they are re-extracted and
# compared rather than trusted silently.
CANONICAL_RECTS: tuple[tuple[int, str, float, float, float, float], ...] = (
    (191, "ALREADY_REMOVED", 445.965, 454.990, 482.996, 459.849),
    (192, "ALREADY_REMOVED", 446.004, 520.660, 478.035, 525.519),
    (193, "ALREADY_REMOVED", 532.960, 519.888, 573.961, 525.552),
    (194, "ALREADY_REMOVED", 623.332, 519.946, 640.348, 526.055),
    (195, "ALREADY_REMOVED", 689.480, 520.155, 790.480, 526.265),
    (196, "ALREADY_REMOVED", 1069.208, 443.235, 1124.192, 449.345),
    (197, "ALREADY_REMOVED", 1008.657, 443.235, 1019.743, 449.345),
    (198, "ALREADY_REMOVED", 625.967, 526.564, 637.052, 623.556),
    (199, "ALREADY_REMOVED", 790.677, 520.910, 801.764, 661.910),
    (200, "ALREADY_REMOVED", 565.108, 460.106, 572.092, 519.114),
    (201, "ALREADY_REMOVED", 1008.708, 285.124, 1015.692, 443.116),
    (202, "PLANNED_DEMOLITION", 756.660, 454.075, 789.660, 460.185),
    (203, "PLANNED_DEMOLITION", 751.806, 385.115, 789.814, 391.225),
)

# Plot-space coordinate versus the current project grid coordinate in mm.
X_ANCHORS: tuple[tuple[float, float], ...] = (
    (-6600.0, 439.831),
    (-4800.0, 538.411),
    (-1800.0, 702.691),
    (-100.0, 795.811),
    (3900.0, 1014.871),
    (6600.0, 1162.711),
)
Y_ANCHORS: tuple[tuple[float, float], ...] = (
    (4800.0, 198.710),
    (3900.0, 247.970),
    (3300.0, 280.850),
    (2700.0, 313.730),
    (100.0, 456.110),
    (-1200.0, 527.270),
    (-4500.0, 707.990),
    (-5300.0, 751.790),
)

# Nominal review boxes are deliberately explicit.  They use the construction
# drawing's 50 mm control grid and are not IFC write authority.
NOMINAL_BBOXES: tuple[tuple[float, float, float, float], ...] = (
    (-6500.0, 0.0, -5800.0, 100.0),
    (-6500.0, -1200.0, -5900.0, -1100.0),
    (-4900.0, -1150.0, -4150.0, -1050.0),
    (-3250.0, -1150.0, -2950.0, -1050.0),
    (-2050.0, -1150.0, -200.0, -1050.0),
    (4900.0, 250.0, 5900.0, 350.0),
    (3800.0, 250.0, 4000.0, 350.0),
    (-3200.0, -3000.0, -3000.0, -1200.0),
    (-200.0, -3700.0, 0.0, -1100.0),
    (-4300.0, -1050.0, -4200.0, 0.0),
    (3800.0, 300.0, 3900.0, 3200.0),
    (-800.0, 0.0, -200.0, 100.0),
    (-900.0, 1300.0, -200.0, 1400.0),
)

CONFIDENCE: tuple[float, ...] = (
    0.75,
    0.75,
    0.85,
    0.70,
    0.75,
    0.65,
    0.65,
    0.90,
    0.90,
    0.75,
    0.80,
    0.80,
    0.60,
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def linear_fit(anchors: Sequence[tuple[float, float]]) -> tuple[float, float, list[float]]:
    world = [item[0] for item in anchors]
    paper = [item[1] for item in anchors]
    world_mean = sum(world) / len(world)
    paper_mean = sum(paper) / len(paper)
    denominator = sum((value - world_mean) ** 2 for value in world)
    slope = sum(
        (world_value - world_mean) * (paper_value - paper_mean)
        for world_value, paper_value in anchors
    ) / denominator
    intercept = paper_mean - slope * world_mean
    residuals = [
        (paper_value - intercept) / slope - world_value
        for world_value, paper_value in anchors
    ]
    return slope, intercept, residuals


def color_close(value: Any, expected: tuple[float, float, float], tolerance: float = 1e-5) -> bool:
    if not isinstance(value, (tuple, list)) or len(value) != 3:
        return False
    return all(abs(float(actual) - target) <= tolerance for actual, target in zip(value, expected))


def extract_pdf_rects(pdf_path: Path) -> tuple[list[tuple[int, str, float, float, float, float]], str]:
    try:
        import pdfplumber  # type: ignore
    except ImportError:
        return list(CANONICAL_RECTS), "hash-guarded-canonical"

    yellow = (1.0, 0.92549, 0.227451)
    red = (0.956863, 0.262745, 0.211765)
    records: list[tuple[int, str, float, float, float, float]] = []
    with pdfplumber.open(pdf_path) as pdf:
        if len(pdf.pages) != 1:
            raise RuntimeError(f"expected one PDF page, found {len(pdf.pages)}")
        page = pdf.pages[0]
        for index, rect in enumerate(page.rects):
            if rect["x0"] <= 400.0 or rect["top"] >= 800.0:
                continue
            stroke = rect.get("stroking_color")
            if color_close(stroke, yellow):
                status = "ALREADY_REMOVED"
            elif color_close(stroke, red):
                status = "PLANNED_DEMOLITION"
            else:
                continue
            records.append(
                (
                    index,
                    status,
                    float(rect["x0"]),
                    float(rect["top"]),
                    float(rect["x1"]),
                    float(rect["bottom"]),
                )
            )
    return records, "pdfplumber-vector-reextract"


def compare_rectangles(
    actual: Sequence[tuple[int, str, float, float, float, float]],
    expected: Sequence[tuple[int, str, float, float, float, float]],
    tolerance: float = 0.01,
) -> float:
    if len(actual) != len(expected):
        raise RuntimeError(f"expected {len(expected)} demolition rectangles, found {len(actual)}")
    maximum = 0.0
    for found, frozen in zip(actual, expected):
        if found[:2] != frozen[:2]:
            raise RuntimeError(f"demolition rectangle identity changed: {found[:2]} != {frozen[:2]}")
        delta = max(abs(a - b) for a, b in zip(found[2:], frozen[2:]))
        maximum = max(maximum, delta)
        if delta > tolerance:
            raise RuntimeError(f"PDF rectangle {found[0]} changed by {delta:.6f} pt")
    return maximum


def world_bbox_from_paper(
    rect: tuple[int, str, float, float, float, float],
    x_fit: tuple[float, float],
    y_fit: tuple[float, float],
) -> tuple[float, float, float, float]:
    _, _, x0, top, x1, bottom = rect
    sx, bx = x_fit
    sy, by = y_fit
    wx0 = (x0 - bx) / sx
    wx1 = (x1 - bx) / sx
    wy_top = (top - by) / sy
    wy_bottom = (bottom - by) / sy
    return wx0, min(wy_top, wy_bottom), wx1, max(wy_top, wy_bottom)


def basis_for(index: int, status: str) -> str:
    basis = (
        "hash-verified PDF vector rectangle; 14 Grid anchors; "
        "50 mm nominal review grid; DWG handover base hash recorded"
    )
    if index == 12:
        basis += "; adjacent to current final wall at X=-200 mm"
    if index == 13:
        basis += (
            f"; overlaps protected adjusted wall {PROTECTED_WALL_GUID}; "
            f"must not move current Opening {PROTECTED_OPENING_GUID}"
        )
    if status == "PLANNED_DEMOLITION":
        basis += "; source PDF marks this segment as planned rather than already removed"
    return basis


def candidate_global_id(candidate_id: str) -> str:
    import ifcopenshell

    return ifcopenshell.guid.compress(uuid.uuid5(GUID_NAMESPACE, candidate_id).hex)


def build_records(
    rects: Sequence[tuple[int, str, float, float, float, float]],
    x_fit: tuple[float, float],
    y_fit: tuple[float, float],
    pdf_hash: str,
    dwg_hash: str,
) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    for index, (rect, nominal, confidence) in enumerate(
        zip(rects, NOMINAL_BBOXES, CONFIDENCE), start=1
    ):
        rect_index, source_status, x0, top, x1, bottom = rect
        raw = world_bbox_from_paper(rect, x_fit, y_fit)
        nx0, ny0, nx1, ny1 = nominal
        records.append(
            {
                "candidate_id": f"D{index:02d}",
                "candidate_global_id": candidate_global_id(f"D{index:02d}"),
                "source_status": source_status,
                "ifc_status_candidate": "DEMOLISH",
                "pdf_page": 1,
                "pdf_rect_index": rect_index,
                "paper_x0_pt": x0,
                "paper_top_pt": top,
                "paper_x1_pt": x1,
                "paper_bottom_pt": bottom,
                "raw_x_min_mm": raw[0],
                "raw_y_min_mm": raw[1],
                "raw_x_max_mm": raw[2],
                "raw_y_max_mm": raw[3],
                "candidate_x_min_mm": nx0,
                "candidate_y_min_mm": ny0,
                "candidate_x_max_mm": nx1,
                "candidate_y_max_mm": ny1,
                "candidate_z_min_mm": 0.0,
                "candidate_z_max_mm": 3000.0,
                "nominal_length_mm": max(nx1 - nx0, ny1 - ny0),
                "nominal_thickness_mm": min(nx1 - nx0, ny1 - ny0),
                "maximum_edge_inference_mm": max(
                    abs(raw[0] - nx0),
                    abs(raw[1] - ny0),
                    abs(raw[2] - nx1),
                    abs(raw[3] - ny1),
                ),
                "basis": basis_for(index, source_status),
                "confidence": confidence,
                "review_required": "yes",
                "review_status": "pending",
                "formal_ifc_write_allowed": "no",
                "pdf_sha256": pdf_hash,
                "dwg_sha256": dwg_hash,
            }
        )
    return records


def candidate_wall_matrix(record: dict[str, Any]):
    import numpy

    x0 = float(record["candidate_x_min_mm"]) / 1000.0
    y0 = float(record["candidate_y_min_mm"]) / 1000.0
    x1 = float(record["candidate_x_max_mm"]) / 1000.0
    y1 = float(record["candidate_y_max_mm"]) / 1000.0
    dx = x1 - x0
    dy = y1 - y0
    matrix = numpy.eye(4)
    if dx >= dy:
        matrix[0, 3] = x0
        matrix[1, 3] = y0
        return matrix, dx, dy
    matrix[0, 0] = 0.0
    matrix[0, 1] = -1.0
    matrix[1, 0] = 1.0
    matrix[1, 1] = 0.0
    matrix[0, 3] = x1
    matrix[1, 3] = y0
    return matrix, dy, dx


def create_candidate_ifc(
    source_path: Path,
    output_path: Path,
    records: Sequence[dict[str, Any]],
    tolerance_mm: float = 0.1,
) -> dict[str, Any]:
    import ifcopenshell
    import ifcopenshell.api
    import ifcopenshell.util.element
    import ifcopenshell.util.representation

    from a103_wall_plan_candidate import geometry_settings, world_bbox_mm
    from geometry_alignment_audit import geometry_difference_audit

    if source_path.resolve() == output_path.resolve():
        raise RuntimeError("A-102 candidate must not overwrite the formal IFC")
    source = ifcopenshell.open(source_path)
    source_root_ids = {root.GlobalId for root in source.by_type("IfcRoot")}
    candidate = ifcopenshell.open(source_path)
    body_context = ifcopenshell.util.representation.get_context(
        candidate, "Model", "Body", "MODEL_VIEW"
    )
    if body_context is None:
        raise RuntimeError("IFC Body/MODEL_VIEW context is missing")
    ffl = next((storey for storey in candidate.by_type("IfcBuildingStorey") if storey.Name == "FFL"), None)
    if ffl is None:
        raise RuntimeError("FFL storey is missing")

    created = []
    for record in records:
        candidate_id = str(record["candidate_id"])
        length_mm = float(record["nominal_length_mm"])
        thickness_mm = float(record["nominal_thickness_mm"])
        height_mm = float(record["candidate_z_max_mm"]) - float(record["candidate_z_min_mm"])
        source_label = "已经拆除" if record["source_status"] == "ALREADY_REMOVED" else "计划拆除"
        wall = ifcopenshell.api.run(
            "root.create_entity",
            candidate,
            ifc_class="IfcWall",
            predefined_type="NOTDEFINED",
            name=f"A102 {candidate_id} {source_label}候选 {length_mm:.0f}×{thickness_mm:.0f}×{height_mm:.0f} mm",
        )
        wall.GlobalId = record["candidate_global_id"]
        wall.Tag = candidate_id
        representation = ifcopenshell.api.run(
            "geometry.add_wall_representation",
            candidate,
            context=body_context,
            length=length_mm / 1000.0,
            height=height_mm / 1000.0,
            thickness=thickness_mm / 1000.0,
        )
        ifcopenshell.api.run(
            "geometry.assign_representation", candidate, product=wall, representation=representation
        )
        matrix, matrix_length, matrix_thickness = candidate_wall_matrix(record)
        if abs(matrix_length * 1000.0 - length_mm) > 1e-6 or abs(matrix_thickness * 1000.0 - thickness_mm) > 1e-6:
            raise RuntimeError(f"A-102 wall matrix dimensions disagree for {candidate_id}")
        ifcopenshell.api.run(
            "geometry.edit_object_placement",
            candidate,
            product=wall,
            matrix=matrix,
            is_si=True,
            should_transform_children=False,
        )
        ifcopenshell.api.run("spatial.assign_container", candidate, products=[wall], relating_structure=ffl)
        common = ifcopenshell.api.run("pset.add_pset", candidate, product=wall, name="Pset_WallCommon")
        ifcopenshell.api.run(
            "pset.edit_pset", candidate, pset=common, properties={"Status": "DEMOLISH"}
        )
        quantities = ifcopenshell.api.run(
            "pset.add_qto", candidate, product=wall, name="Qto_WallBaseQuantities"
        )
        ifcopenshell.api.run(
            "pset.edit_qto",
            candidate,
            qto=quantities,
            properties={
                "Length": length_mm,
                "Width": thickness_mm,
                "Height": height_mm,
            },
        )
        review = ifcopenshell.api.run(
            "pset.add_pset", candidate, product=wall, name="Pset_A102DemolitionReview"
        )
        ifcopenshell.api.run(
            "pset.edit_pset",
            candidate,
            pset=review,
            properties={
                "CandidateId": candidate_id,
                "SourceStatus": str(record["source_status"]),
                "Confidence": float(record["confidence"]),
                "ReviewStatus": "PENDING",
                "FormalIfcWriteAllowed": False,
                "SourcePdfSha256": str(record["pdf_sha256"]),
                "SourceDwgSha256": str(record["dwg_sha256"]),
                "InferenceBasis": str(record["basis"]),
            },
        )
        created.append(wall.GlobalId)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    candidate.write(output_path)
    reopened = ifcopenshell.open(output_path)
    settings = geometry_settings()
    wall_checks: list[dict[str, Any]] = []
    for record in records:
        wall = reopened.by_guid(record["candidate_global_id"])
        bbox = world_bbox_mm(settings, wall)
        expected_min = [
            float(record["candidate_x_min_mm"]),
            float(record["candidate_y_min_mm"]),
            float(record["candidate_z_min_mm"]),
        ]
        expected_max = [
            float(record["candidate_x_max_mm"]),
            float(record["candidate_y_max_mm"]),
            float(record["candidate_z_max_mm"]),
        ]
        bbox_delta = max(
            abs(actual - expected)
            for actual, expected in zip(
                bbox["min_mm"] + bbox["max_mm"], expected_min + expected_max
            )
        )
        psets = ifcopenshell.util.element.get_psets(wall)
        quantities = psets.get("Qto_WallBaseQuantities", {})
        wall_checks.append(
            {
                "candidate_id": record["candidate_id"],
                "global_id": wall.GlobalId,
                "ifc_class": wall.is_a(),
                "name": wall.Name,
                "tag": wall.Tag,
                "bbox_mm": bbox,
                "maximum_bbox_delta_mm": bbox_delta,
                "status": psets.get("Pset_WallCommon", {}).get("Status"),
                "review_status": psets.get("Pset_A102DemolitionReview", {}).get("ReviewStatus"),
                "quantities_mm": {
                    "Length": quantities.get("Length"),
                    "Width": quantities.get("Width"),
                    "Height": quantities.get("Height"),
                },
                "within_tolerance": bbox_delta <= tolerance_mm,
            }
        )

    compared_classes = ("IfcWall", "IfcOpeningElement", "IfcDoor", "IfcWindow")
    original_ids = [
        product.GlobalId
        for ifc_class in compared_classes
        for product in source.by_type(ifc_class)
        if getattr(product, "GlobalId", None)
    ]
    original_geometry = geometry_difference_audit(
        reopened,
        source,
        str(source_path.resolve()),
        tolerance_mm=tolerance_mm,
        classes=(),
        global_ids=original_ids,
    )
    reopened_root_ids = {root.GlobalId for root in reopened.by_type("IfcRoot")}
    status_counts = {
        status: sum(
            ifcopenshell.util.element.get_psets(wall).get("Pset_WallCommon", {}).get("Status") == status
            for wall in reopened.by_type("IfcWall")
        )
        for status in ("EXISTING", "NEW", "DEMOLISH")
    }
    passed = (
        len(created) == 13
        and len(set(created)) == 13
        and len(reopened.by_type("IfcWall")) == 101
        and status_counts == {"EXISTING": 84, "NEW": 4, "DEMOLISH": 13}
        and all(check["ifc_class"] == "IfcWall" for check in wall_checks)
        and all(check["status"] == "DEMOLISH" for check in wall_checks)
        and all(check["review_status"] == "PENDING" for check in wall_checks)
        and all(
            check["quantities_mm"]
            == {
                "Length": float(record["nominal_length_mm"]),
                "Width": float(record["nominal_thickness_mm"]),
                "Height": float(record["candidate_z_max_mm"]) - float(record["candidate_z_min_mm"]),
            }
            for check, record in zip(wall_checks, records)
        )
        and all(check["within_tolerance"] for check in wall_checks)
        and original_geometry["total"] == len(original_ids)
        and original_geometry["over_tolerance"] == 0
        and source_root_ids <= reopened_root_ids
    )
    return {
        "path": str(output_path.resolve()),
        "sha256": sha256(output_path),
        "schema": reopened.schema,
        "wall_count": len(reopened.by_type("IfcWall")),
        "status_counts": status_counts,
        "created_global_ids": created,
        "wall_checks": wall_checks,
        "original_product_geometry": original_geometry,
        "source_root_ids_preserved": source_root_ids <= reopened_root_ids,
        "pass": passed,
    }


def write_register(path: Path, records: Sequence[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(records)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pdf", required=True, type=Path)
    parser.add_argument("--dwg", required=True, type=Path)
    parser.add_argument("--ifc", required=True, type=Path)
    parser.add_argument("--output-ifc", required=True, type=Path)
    parser.add_argument("--output-register", required=True, type=Path)
    parser.add_argument("--report", required=True, type=Path)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    pdf_hash = sha256(args.pdf)
    dwg_hash = sha256(args.dwg)
    if pdf_hash != EXPECTED_PDF_SHA256:
        raise RuntimeError(f"unexpected demolition PDF hash: {pdf_hash}")
    if dwg_hash != EXPECTED_DWG_SHA256:
        raise RuntimeError(f"unexpected handover DWG hash: {dwg_hash}")
    ifc_hash_before = sha256(args.ifc)
    ifc_bytes = args.ifc.read_bytes()
    for global_id in (PROTECTED_WALL_GUID, PROTECTED_OPENING_GUID):
        if global_id.encode() not in ifc_bytes:
            raise RuntimeError(f"protected IFC GlobalId is missing: {global_id}")

    rects, vector_mode = extract_pdf_rects(args.pdf)
    rectangle_delta = compare_rectangles(rects, CANONICAL_RECTS)
    sx, bx, x_residuals = linear_fit(X_ANCHORS)
    sy, by, y_residuals = linear_fit(Y_ANCHORS)
    records = build_records(rects, (sx, bx), (sy, by), pdf_hash, dwg_hash)
    write_register(args.output_register, records)
    candidate_ifc = create_candidate_ifc(args.ifc, args.output_ifc, records)

    ifc_hash_after = sha256(args.ifc)
    phase_counts = {
        status: sum(record["source_status"] == status for record in records)
        for status in ("ALREADY_REMOVED", "PLANNED_DEMOLITION")
    }
    candidate_grid_errors = [
        abs(float(record[key]) / 50.0 - round(float(record[key]) / 50.0))
        for record in records
        for key in (
            "candidate_x_min_mm",
            "candidate_y_min_mm",
            "candidate_x_max_mm",
            "candidate_y_max_mm",
        )
    ]
    gates = {
        "candidate_count": len(records),
        "source_status_counts": phase_counts,
        "maximum_grid_anchor_residual_mm": max(
            max(abs(value) for value in x_residuals),
            max(abs(value) for value in y_residuals),
        ),
        "maximum_pdf_rectangle_reextract_delta_pt": rectangle_delta,
        "candidate_edges_on_50mm_grid": max(candidate_grid_errors) <= 1e-9,
        "formal_ifc_sha256_before": ifc_hash_before,
        "formal_ifc_sha256_after": ifc_hash_after,
        "formal_ifc_unchanged": ifc_hash_before == ifc_hash_after,
        "protected_wall_global_id": PROTECTED_WALL_GUID,
        "protected_opening_global_id": PROTECTED_OPENING_GUID,
        "automatic_ifc_write_allowed": False,
        "candidate_ifc_pass": candidate_ifc["pass"],
        "candidate_ifc_wall_count": candidate_ifc["wall_count"],
        "candidate_ifc_status_counts": candidate_ifc["status_counts"],
        "candidate_ifc_original_products_over_tolerance": candidate_ifc["original_product_geometry"]["over_tolerance"],
    }
    passed = (
        len(records) == 13
        and phase_counts == {"ALREADY_REMOVED": 11, "PLANNED_DEMOLITION": 2}
        and gates["maximum_grid_anchor_residual_mm"] <= 0.5
        and gates["candidate_edges_on_50mm_grid"]
        and gates["formal_ifc_unchanged"]
        and gates["candidate_ifc_pass"]
        and all(record["review_required"] == "yes" for record in records)
    )
    report = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "mode": "read-only-a102-demolition-review-candidate",
        "sources": {
            "pdf": {"path": str(args.pdf.resolve()), "sha256": pdf_hash, "vector_mode": vector_mode},
            "dwg": {"path": str(args.dwg.resolve()), "sha256": dwg_hash, "use": "handover geometry reference"},
            "ifc": {"path": str(args.ifc.resolve()), "sha256": ifc_hash_before, "use": "protected final-built model"},
            "candidate_ifc": candidate_ifc,
        },
        "mapping": {
            "paper_x_equals_slope_times_world_x_plus_intercept": {"slope": sx, "intercept": bx},
            "paper_top_equals_slope_times_world_y_plus_intercept": {"slope": sy, "intercept": by},
            "x_anchor_residuals_mm": x_residuals,
            "y_anchor_residuals_mm": y_residuals,
        },
        "records": records,
        "gates": gates,
        "pass": passed,
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps({"report": str(args.report), "register": str(args.output_register), "gates": gates, "pass": passed}, ensure_ascii=False))
    if not passed:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
