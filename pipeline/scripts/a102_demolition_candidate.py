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
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Sequence


EXPECTED_PDF_SHA256 = "84a18d211fe8e6c3e05f56b35027a96799d7c54cd4bc75cd602c72f65ea327f0"
EXPECTED_DWG_SHA256 = "ba355f6a90732ad07f843d59e8bab5e1da9daffe7aed74889d1d265b2ce22d7e"
PROTECTED_WALL_GUID = "0hKdvAZkn1TejLgJhK_vDp"
PROTECTED_OPENING_GUID = "1YxMx6s0r3ZPPohkRKXWbl"

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
    }
    passed = (
        len(records) == 13
        and phase_counts == {"ALREADY_REMOVED": 11, "PLANNED_DEMOLITION": 2}
        and gates["maximum_grid_anchor_residual_mm"] <= 0.5
        and gates["candidate_edges_on_50mm_grid"]
        and gates["formal_ifc_unchanged"]
        and all(record["review_required"] == "yes" for record in records)
    )
    report = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "mode": "read-only-a102-demolition-review-candidate",
        "sources": {
            "pdf": {"path": str(args.pdf.resolve()), "sha256": pdf_hash, "vector_mode": vector_mode},
            "dwg": {"path": str(args.dwg.resolve()), "sha256": dwg_hash, "use": "handover geometry reference"},
            "ifc": {"path": str(args.ifc.resolve()), "sha256": ifc_hash_before, "use": "protected final-built model"},
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
