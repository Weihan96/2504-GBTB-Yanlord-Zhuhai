#!/usr/bin/env python3
"""Shared native-DWG and vector-PDF geometry helpers for Falper Sorgente."""

from __future__ import annotations

import hashlib
import json
import math
import shutil
import subprocess
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SCOPE = "family_reference_not_project_shop_drawing"
PRODUCT_URL = "https://falper.it/en/freestanding-washbasins/"
AUTOCAD_2D_URL = (
    "https://falper.it/wp-content/uploads/2025/11/"
    "Lavabi-Freestanding-Autocad-2D.zip"
)
TECHNICAL_DWG_URL = (
    "https://falper.it/wp-content/uploads/2025/11/"
    "Lavabi-Freestanding-Scheda-tecnica-DWG.zip"
)
AUTOCAD_3D_URL = (
    "https://falper.it/wp-content/uploads/2025/11/"
    "Lavabi-Freestanding-Autocad-3D.zip"
)
EXPECTED = {
    "formal_ifc": "7a50b87e8f48a7c2bfbaa9f1dabdfca7155fa684b2325c66aa1ae15c8ab0a25c",
    "pdf": "b8bc94982bdd4305f2ffe5240daf429146470677157022860a0ef591514373c2",
    "wfa_2d": "48ffe94e5a692cd03e8a09e67c658ebb42362817df03007bd90ff7a912567441",
    "wfb_2d": "4b3d8900f8ae249d7e9634339f3529d7d38a940dc2d0c909a7e90e31d3946831",
    "wfa_3d": "6a9d4b45607ae7915b2b5925dfa818693984df4b3474bbf6b17ab3d53b72e956",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def relative(path: Path) -> str:
    try:
        return path.resolve().relative_to(ROOT).as_posix()
    except ValueError:
        return path.resolve().as_posix()


def load_json(path: Path) -> dict:
    with path.open(encoding="utf-8") as handle:
        return json.load(handle)


def write_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )


def _de_boor_point(entity: dict, parameter: float) -> tuple[float, float]:
    degree = int(entity["degree"])
    knots = [float(value) for value in entity["knots"]]
    controls = [
        (float(point["x"]), float(point["y"])) for point in entity["ctrl_pts"]
    ]
    count = len(controls) - 1
    if parameter >= knots[count + 1]:
        return controls[-1]
    span = degree
    for index in range(degree, count + 1):
        if knots[index] <= parameter < knots[index + 1]:
            span = index
            break
    values = [list(controls[span - degree + index]) for index in range(degree + 1)]
    for level in range(1, degree + 1):
        for index in range(degree, level - 1, -1):
            knot_index = span - degree + index
            denominator = knots[knot_index + degree - level + 1] - knots[knot_index]
            alpha = 0.0 if abs(denominator) < 1e-15 else (
                (parameter - knots[knot_index]) / denominator
            )
            values[index][0] = (
                (1.0 - alpha) * values[index - 1][0] + alpha * values[index][0]
            )
            values[index][1] = (
                (1.0 - alpha) * values[index - 1][1] + alpha * values[index][1]
            )
    return values[degree][0], values[degree][1]


def sample_spline(entity: dict, steps: int = 768) -> list[list[float]]:
    degree = int(entity["degree"])
    controls = entity["ctrl_pts"]
    knots = [float(value) for value in entity["knots"]]
    start = knots[degree]
    end = knots[len(controls)]
    return [
        [round(value, 9) for value in _de_boor_point(entity, start + (end - start) * index / steps)]
        for index in range(steps + 1)
    ]


def sample_circle(center: list[float], radius: float, steps: int = 256) -> list[list[float]]:
    cx, cy = float(center[0]), float(center[1])
    return [
        [
            round(cx + radius * math.cos(2.0 * math.pi * index / steps), 9),
            round(cy + radius * math.sin(2.0 * math.pi * index / steps), 9),
        ]
        for index in range(steps + 1)
    ]


def _dwg_json(path: Path) -> tuple[dict, str]:
    executable = shutil.which("dwgread")
    if not executable:
        raise RuntimeError(
            "LibreDWG dwgread is required to parse the native DWG; install libredwg"
        )
    with tempfile.TemporaryDirectory(prefix="falper-native-dwg-") as temporary:
        output = Path(temporary) / "source.json"
        completed = subprocess.run(
            [executable, "-O", "JSON", "-o", str(output), str(path)],
            capture_output=True,
            text=True,
        )
        if completed.returncode != 0 or not output.is_file():
            raise RuntimeError(
                f"dwgread failed for {path}: {completed.stderr.strip()}"
            )
        # LibreDWG preserves legacy MText codepage bytes in its JSON strings.
        document = json.loads(output.read_bytes().decode("latin-1"))
    version = subprocess.run(
        [executable, "--version"], capture_output=True, text=True, check=True
    ).stdout.splitlines()[0]
    return document, version


def _same_center(left: list[float], right: list[float], tolerance: float = 0.01) -> bool:
    return math.dist((float(left[0]), float(left[1])), (float(right[0]), float(right[1]))) <= tolerance


def extract_native_dwg_variant(path: Path, code: str) -> dict:
    expected_key = f"{code.lower()}_2d"
    actual_hash = sha256(path)
    if actual_hash != EXPECTED[expected_key]:
        raise RuntimeError(
            f"{code} native DWG hash mismatch: expected {EXPECTED[expected_key]}, got {actual_hash}"
        )
    document, parser_version = _dwg_json(path)
    entities = [item for item in document["OBJECTS"] if item.get("entmode") == 2]

    circles = [item for item in entities if item.get("entity") == "CIRCLE"]
    plan_group = None
    for outer in circles:
        if abs(float(outer.get("radius", 0.0)) - 260.0) > 0.01:
            continue
        group = [item for item in circles if _same_center(item["center"], outer["center"])]
        radii = [float(item["radius"]) for item in group]
        if all(any(abs(value - expected) <= 0.01 for value in radii) for expected in (260, 250, 50, 45)):
            plan_group = sorted(group, key=lambda item: float(item["radius"]), reverse=True)
            break
    if plan_group is None or len(plan_group) != 5:
        raise RuntimeError(f"{code}: expected one five-circle Sorgente plan group")
    plan_center = plan_group[0]["center"]
    plan_paths = []
    for circle in plan_group:
        raw = sample_circle(circle["center"], float(circle["radius"]))
        plan_paths.append(
            [
                [
                    round(point[0] - float(plan_center[0]), 6),
                    round(point[1] - float(plan_center[1]), 6),
                ]
                for point in raw
            ]
        )

    splines = []
    for entity in entities:
        if entity.get("entity") != "SPLINE" or int(entity.get("degree", 0)) != 3:
            continue
        controls = entity.get("ctrl_pts", [])
        if not controls:
            continue
        minimum_y = min(float(point["y"]) for point in controls)
        maximum_y = max(float(point["y"]) for point in controls)
        minimum_x = min(float(point["x"]) for point in controls)
        maximum_x = max(float(point["x"]) for point in controls)
        if (
            850.0 <= maximum_y - minimum_y <= 950.0
            and float(plan_center[0]) - 350.0 <= minimum_x
            and maximum_x <= float(plan_center[0]) + 350.0
        ):
            splines.append(entity)
    if len(splines) != 2:
        raise RuntimeError(f"{code}: expected two native elevation SPLINE entities, found {len(splines)}")
    sampled_splines = [(entity, sample_spline(entity)) for entity in splines]
    sampled_splines.sort(key=lambda item: min(point[0] for point in item[1]))
    curve_left = min(point[0] for _, path_points in sampled_splines for point in path_points)
    curve_right = max(point[0] for _, path_points in sampled_splines for point in path_points)
    curve_bottom = min(point[1] for _, path_points in sampled_splines for point in path_points)
    curve_top = max(point[1] for _, path_points in sampled_splines for point in path_points)

    top_candidates = []
    base_candidates = []
    for entity in entities:
        if entity.get("entity") != "LINE":
            continue
        start, end = entity.get("start"), entity.get("end")
        if not start or not end or abs(float(start[1]) - float(end[1])) > 0.01:
            continue
        width = abs(float(start[0]) - float(end[0]))
        y = float(start[1])
        if 450.0 <= width <= 550.0 and abs(y - curve_top) <= 5.0:
            top_candidates.append(entity)
        if 180.0 <= width <= 260.0 and abs(y - curve_bottom) <= 10.0:
            base_candidates.append(entity)
    if not top_candidates or not base_candidates:
        raise RuntimeError(
            f"{code}: native elevation top/base line candidates are incomplete"
        )
    # Dimension witness and rim lines sit within a few millimetres of the
    # outline.  The exterior top is the highest wide line; floor contact is
    # the lowest short line.
    top_line = max(top_candidates, key=lambda item: float(item["start"][1]))
    base_line = min(base_candidates, key=lambda item: float(item["start"][1]))
    left = curve_left
    right = curve_right
    top = float(top_line["start"][1])
    bottom = float(base_line["start"][1])
    center_x = (left + right) / 2.0
    scale_x = 520.0 / (right - left)
    scale_z = 900.0 / (top - bottom)

    def normalize_elevation(points: list[list[float]]) -> list[list[float]]:
        return [
            [
                round((float(point[0]) - center_x) * scale_x, 6),
                round((float(point[1]) - bottom) * scale_z, 6),
            ]
            for point in points
        ]

    elevation_paths = [normalize_elevation(points) for _, points in sampled_splines]
    for line in (top_line, base_line):
        elevation_paths.append(normalize_elevation([line["start"], line["end"]]))

    return {
        "model_code": code,
        "source_dwg": relative(path),
        "source_dwg_sha256": actual_hash,
        "source_kind": "native_dwg",
        "dwg_file_version": document["FILEHEADER"]["version"],
        "parser": parser_version,
        "units": "mm",
        "scope": SCOPE,
        "views": {
            "plan": {
                "paths_mm": plan_paths,
                "native_entity_types": ["CIRCLE"] * 5,
                "native_entity_indices": [int(item["index"]) for item in plan_group],
                "normalisation": {
                    "origin": "outer_circle_center",
                    "scale": 1.0,
                    "outer_diameter_mm": 520.0,
                },
            },
            "elevation": {
                "paths_mm": elevation_paths,
                "native_entity_types": ["SPLINE", "SPLINE", "LINE", "LINE"],
                "native_entity_indices": [
                    int(sampled_splines[0][0]["index"]),
                    int(sampled_splines[1][0]["index"]),
                    int(top_line["index"]),
                    int(base_line["index"]),
                ],
                "normalisation": {
                    "published_envelope_mm": [520.0, 900.0],
                    "native_exterior_bbox_mm": [left, bottom, right, top],
                    "horizontal_scale": round(scale_x, 12),
                    "vertical_scale": round(scale_z, 12),
                },
            },
        },
    }


def paths_geometry(paths: list[list[list[float]]]):
    from shapely.geometry import LineString, MultiLineString

    lines = [LineString(path) for path in paths if len(path) >= 2]
    return MultiLineString(lines)


def pathset_diff_mm(
    expected_paths: list[list[list[float]]],
    observed_paths: list[list[list[float]]],
) -> dict:
    expected = paths_geometry(expected_paths)
    observed = paths_geometry(observed_paths)
    return {
        "hausdorff_mm": round(expected.hausdorff_distance(observed), 6),
        "expected_length_mm": round(expected.length, 6),
        "observed_length_mm": round(observed.length, 6),
        "length_delta_mm": round(abs(expected.length - observed.length), 6),
    }
