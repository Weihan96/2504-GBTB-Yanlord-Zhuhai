#!/usr/bin/env python3
"""Build the P-202 existing drainage and sanitary location candidate.

The existing Sanitary Plan remains the underlay. Observable drainage geometry
is projected from the frozen IFC, while small placement markers register all
49 P-202 objects. Each placement marker is not a connector or rough-in point,
and it never provides routing authority.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import subprocess
import sys
import xml.etree.ElementTree as ET
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

import ifcopenshell
import ifcopenshell.util.placement
from shapely.geometry import Polygon
from shapely.ops import unary_union

from geometry_alignment_audit import geometry_settings, world_mesh_mm
from plum_existing_candidate import (
    ASSEMBLY_IDS,
    PVC110_IDS,
    drainage_products,
)


SVG_NS = "http://www.w3.org/2000/svg"
IFC_NS = "http://www.ifcopenshell.org/ns"
SVG = f"{{{SVG_NS}}}"
IFC = f"{{{IFC_NS}}}"
SCALE = 50.0
ORIGIN_X = 200.0
ORIGIN_Y = 200.0
PAGE_WIDTH_MM = 500.0
PAGE_HEIGHT_MM = 400.0
GEOMETRY_SETTINGS = geometry_settings()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument("--ifc", type=Path, default=Path("2504 GBTB Yanlord Zhuhai.ifc"))
    parser.add_argument("--underlay", type=Path, default=Path("drawings/Sanitary Plan.svg"))
    parser.add_argument("--registry", type=Path, default=Path("build/plum/p202-existing-location-register.json"))
    parser.add_argument("--plum-report", type=Path, default=Path("build/plum/plum-report.json"))
    parser.add_argument("--output-svg", type=Path, default=Path("drawings/Sanitary Plan-P202-candidate.svg"))
    parser.add_argument("--output-pdf", type=Path, default=Path("output/pdf/P-202-drainage-sanitary-candidate.pdf"))
    parser.add_argument("--report", type=Path, default=Path("build/plum/p202-candidate-report.json"))
    parser.add_argument("--render-report", type=Path, default=Path("build/plum/p202-pdf-render-report.json"))
    parser.add_argument("--proof-png", type=Path, default=Path("build/plum/p202-proof.png"))
    parser.add_argument(
        "--expected-ifc-sha256",
        help="Optional caller-frozen source hash; omit to use the hash-checked PLUM registry",
    )
    parser.add_argument("--render-script", type=Path, default=Path("pipeline/scripts/render_svg_pdf.py"))
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def svg_point(point_mm: Iterable[float]) -> tuple[float, float]:
    x, y = list(point_mm)[:2]
    return ORIGIN_X + float(x) / SCALE, ORIGIN_Y - float(y) / SCALE


def world_origin_mm(product: Any) -> list[float]:
    matrix = ifcopenshell.util.placement.get_local_placement(product.ObjectPlacement)
    return [float(value) for value in matrix[:3, 3]]


def path_from_geometry(geometry: Any) -> str:
    polygons: list[str] = []
    if geometry.is_empty:
        return ""
    shapes = list(geometry.geoms) if hasattr(geometry, "geoms") else [geometry]
    for shape in shapes:
        if not isinstance(shape, Polygon):
            continue
        rings = [shape.exterior, *shape.interiors]
        for ring in rings:
            points = [svg_point(point) for point in ring.coords]
            if len(points) < 3:
                continue
            polygons.append(
                "M" + " L".join(f"{x:.4f},{y:.4f}" for x, y in points) + " Z"
            )
    return " ".join(polygons)


def projected_product_path(product: Any) -> str:
    vertices, faces = world_mesh_mm(GEOMETRY_SETTINGS, product)
    triangles = []
    for face in faces:
        polygon = Polygon([(vertices[index][0], vertices[index][1]) for index in face])
        if polygon.is_valid and polygon.area > 1e-6:
            triangles.append(polygon)
    if not triangles:
        return ""
    projected = unary_union(triangles).buffer(0).simplify(0.05, preserve_topology=True)
    return path_from_geometry(projected)


def add_text(parent: Any, x: float, y: float, text: str, css_class: str, **attributes: str) -> Any:
    node = ET.SubElement(parent, SVG + "text", {"x": str(x), "y": str(y), "class": css_class, **attributes})
    node.text = text
    return node


def add_marker(group: Any, product: Any, marker_class: str, label: str) -> None:
    x, y = svg_point(world_origin_mm(product))
    marker = ET.SubElement(
        group,
        SVG + "g",
        {
            "class": f"p202-location-marker {marker_class}",
            "data-guid": product.GlobalId,
            "data-role": "existing-object-placement-not-connector",
        },
    )
    ET.SubElement(marker, SVG + "circle", {"cx": f"{x:.4f}", "cy": f"{y:.4f}", "r": "1.35"})
    ET.SubElement(marker, SVG + "line", {"x1": f"{x-1.8:.4f}", "y1": f"{y:.4f}", "x2": f"{x+1.8:.4f}", "y2": f"{y:.4f}"})
    ET.SubElement(marker, SVG + "line", {"x1": f"{x:.4f}", "y1": f"{y-1.8:.4f}", "x2": f"{x:.4f}", "y2": f"{y+1.8:.4f}"})
    marker.set("aria-label", f"{label} {product.GlobalId}")


def source_sanitary_guids(root: Any) -> set[str]:
    return {
        group.get(IFC + "guid")
        for group in root.iter(SVG + "g")
        if "IfcSanitaryTerminal" in group.get("class", "") and group.get(IFC + "guid")
    }


def parse_pgm(path: Path) -> tuple[int, int, bytes]:
    data = path.read_bytes()
    tokens: list[bytes] = []
    index = 0
    while len(tokens) < 4:
        while index < len(data) and data[index:index + 1].isspace():
            index += 1
        if data[index:index + 1] == b"#":
            index = data.index(b"\n", index) + 1
            continue
        end = index
        while end < len(data) and not data[end:end + 1].isspace():
            end += 1
        tokens.append(data[index:end])
        index = end
    while index < len(data) and data[index:index + 1].isspace():
        index += 1
    if tokens[0] != b"P5" or tokens[3] != b"255":
        raise RuntimeError(f"unsupported PGM header in {path}")
    width, height = int(tokens[1]), int(tokens[2])
    pixels = data[index:index + width * height]
    if len(pixels) != width * height:
        raise RuntimeError(f"truncated PGM pixels in {path}")
    return width, height, pixels


def raster_metrics(pdf_path: Path, pgm_path: Path) -> dict[str, Any]:
    subprocess.run(
        ["pdftoppm", "-f", "1", "-singlefile", "-r", "20", "-gray", str(pdf_path), str(pgm_path.with_suffix(""))],
        check=True,
        capture_output=True,
    )
    width, height, pixels = parse_pgm(pgm_path)
    ink = [index for index, value in enumerate(pixels) if value < 245]
    if not ink:
        return {"width_px": width, "height_px": height, "ink_ratio": 0.0, "ink_bbox_page_ratio": 0.0}
    xs = [index % width for index in ink]
    ys = [index // width for index in ink]
    return {
        "width_px": width,
        "height_px": height,
        "ink_ratio": len(ink) / len(pixels),
        "ink_bbox_page_ratio": ((max(xs) - min(xs) + 1) * (max(ys) - min(ys) + 1)) / (width * height),
    }


def main() -> None:
    args = parse_args()
    root_path = args.root.resolve()
    resolve = lambda path: path if path.is_absolute() else root_path / path
    ifc_path = resolve(args.ifc)
    underlay_path = resolve(args.underlay)
    registry_path = resolve(args.registry)
    plum_report_path = resolve(args.plum_report)
    output_svg = resolve(args.output_svg)
    output_pdf = resolve(args.output_pdf)
    report_path = resolve(args.report)
    render_report_path = resolve(args.render_report)
    proof_png = resolve(args.proof_png)
    render_script = resolve(args.render_script)

    ifc_sha = sha256(ifc_path)
    if args.expected_ifc_sha256 and ifc_sha != args.expected_ifc_sha256:
        raise RuntimeError(f"formal IFC hash changed: expected {args.expected_ifc_sha256}, got {ifc_sha}")
    registry = json.loads(registry_path.read_text(encoding="utf-8"))
    plum_report = json.loads(plum_report_path.read_text(encoding="utf-8"))
    if registry["source_ifc_sha256"] != ifc_sha or plum_report["source"]["ifc_sha256"] != ifc_sha:
        raise RuntimeError("PLUM source registry/report is stale for the formal IFC")
    if len(registry["objects"]) != 49 or not plum_report["qa"]["candidate_registry_pass"]:
        raise RuntimeError("P-202 existing-object registry is incomplete or failed")

    model = ifcopenshell.open(ifc_path)
    sanitary = sorted(model.by_type("IfcSanitaryTerminal"), key=lambda item: item.GlobalId)
    waste = sorted(model.by_type("IfcWasteTerminal"), key=lambda item: item.GlobalId)
    assemblies = [model.by_guid(global_id) for global_id in ASSEMBLY_IDS]
    drainage = drainage_products(model)
    if (len(sanitary), len(waste), len(assemblies), len(drainage)) != (27, 3, 3, 16):
        raise RuntimeError("P-202 observable object counts changed")

    ET.register_namespace("", SVG_NS)
    ET.register_namespace("ifc", IFC_NS)
    tree = ET.parse(underlay_path)
    svg_root = tree.getroot()
    base_sanitary = source_sanitary_guids(svg_root)
    svg_root.set("width", f"{PAGE_WIDTH_MM}mm")
    svg_root.set("height", f"{PAGE_HEIGHT_MM}mm")
    svg_root.set("viewBox", f"0 0 {PAGE_WIDTH_MM:.0f} {PAGE_HEIGHT_MM:.0f}")
    svg_root.set("data-sheet", "P-202")
    svg_root.set("data-source-ifc-sha256", ifc_sha)
    svg_root.set("data-status", "candidate-existing-location-connectivity-missing")

    definitions = svg_root.find(SVG + "defs")
    if definitions is None:
        definitions = ET.SubElement(svg_root, SVG + "defs")
    style = ET.SubElement(definitions, SVG + "style", {"type": "text/css"})
    style.text = """
@page{size:500mm 400mm;margin:0}html,body{margin:0;width:500mm;height:400mm;overflow:hidden}
.p202-overlay{pointer-events:none}.p202-drainage-geometry{fill:#228be6;fill-opacity:.18;stroke:#1864ab;stroke-width:.45;fill-rule:evenodd}
.p202-waste-geometry{fill:#e64980;fill-opacity:.28;stroke:#a61e4d;stroke-width:.55;fill-rule:evenodd}
.p202-location-marker circle{fill:#fff;stroke-width:.45}.p202-location-marker line{stroke-width:.35}
.p202-sanitary-marker circle,.p202-sanitary-marker line{stroke:#2b8a3e}.p202-waste-marker circle,.p202-waste-marker line{stroke:#a61e4d}
.p202-drainage-marker circle,.p202-drainage-marker line{stroke:#1864ab}.p202-assembly-marker circle,.p202-assembly-marker line{stroke:#6741d9}
.p202-panel-bg{fill:#fff;stroke:#102f43;stroke-width:.55}.p202-panel-rule{stroke:#ced4da;stroke-width:.35}
.p202-title{font:700 7px Arial,"Noto Sans CJK SC",sans-serif;fill:#102f43}.p202-sub{font:3.1px Arial,"Noto Sans CJK SC",sans-serif;fill:#526777}
.p202-head{font:700 3.4px Arial,"Noto Sans CJK SC",sans-serif;fill:#102f43}.p202-note{font:2.8px Arial,"Noto Sans CJK SC",sans-serif;fill:#243b53}
.p202-small{font:2.45px Arial,"Noto Sans CJK SC",sans-serif;fill:#526777}.p202-badge{font:700 2.7px Arial,"Noto Sans CJK SC",sans-serif}
"""

    overlay = ET.SubElement(svg_root, SVG + "g", {"class": "p202-overlay", "data-layer": "observable-ifc-overlays"})
    drain_geometry_count = 0
    for product in drainage:
        path_data = projected_product_path(product)
        if path_data:
            ET.SubElement(overlay, SVG + "path", {
                "d": path_data, "class": "p202-drainage-geometry", "data-guid": product.GlobalId,
            })
            drain_geometry_count += 1
    waste_geometry_count = 0
    for product in waste:
        path_data = projected_product_path(product)
        if path_data:
            ET.SubElement(overlay, SVG + "path", {
                "d": path_data, "class": "p202-waste-geometry", "data-guid": product.GlobalId,
            })
            waste_geometry_count += 1

    marker_group = ET.SubElement(overlay, SVG + "g", {"data-layer": "existing-object-placement-markers"})
    for product in sanitary:
        add_marker(marker_group, product, "p202-sanitary-marker", "洁具既有 ObjectPlacement")
    for product in waste:
        add_marker(marker_group, product, "p202-waste-marker", "地漏既有 ObjectPlacement")
    for product in drainage:
        add_marker(marker_group, product, "p202-drainage-marker", "排水对象既有 ObjectPlacement")
    for product in assemblies:
        add_marker(marker_group, product, "p202-assembly-marker", "组合件既有 ObjectPlacement")

    panel = ET.SubElement(svg_root, SVG + "g", {"data-layer": "p202-title-panel"})
    ET.SubElement(panel, SVG + "rect", {"x": "400", "y": "0", "width": "100", "height": "400", "class": "p202-panel-bg"})
    add_text(panel, 407, 18, "P-202", "p202-title")
    add_text(panel, 407, 27, "既有排水及洁具定位候选", "p202-head")
    add_text(panel, 407, 35, "非施工发布版", "p202-sub")
    ET.SubElement(panel, SVG + "rect", {"x": "407", "y": "43", "width": "39", "height": "7", "rx": "2", "fill": "#d3f9d8"})
    add_text(panel, 426.5, 48, "既有定位", "p202-badge", **{"text-anchor": "middle", "fill": "#2b8a3e"})
    ET.SubElement(panel, SVG + "rect", {"x": "450", "y": "43", "width": "43", "height": "7", "rx": "2", "fill": "#ffe3e3"})
    add_text(panel, 471.5, 48, "连接数据缺失", "p202-badge", **{"text-anchor": "middle", "fill": "#c92a2a"})
    add_text(panel, 407, 63, "当前 IFC 登记", "p202-head")
    for index, text in enumerate((
        "27 件洁具 / 卫生终端",
        "3 个线性地漏",
        "16 个可识别排水对象",
        "3 个洁具组合件",
        "PVC110：2 产品 / 6 分支",
    )):
        add_text(panel, 410, 72 + index * 7, text, "p202-note")
    ET.SubElement(panel, SVG + "line", {"x1": "407", "y1": "110", "x2": "493", "y2": "110", "class": "p202-panel-rule"})
    add_text(panel, 407, 122, "图例", "p202-head")
    legend = (
        ("#2b8a3e", "绿圈：洁具 ObjectPlacement"),
        ("#e64980", "粉色：线性地漏几何 / 原点"),
        ("#228be6", "蓝色：既有排水几何 / 原点"),
        ("#6741d9", "紫圈：洁具组合件原点"),
    )
    for index, (color, text) in enumerate(legend):
        y = 133 + index * 9
        ET.SubElement(panel, SVG + "circle", {"cx": "411", "cy": str(y-1), "r": "2", "fill": "#fff", "stroke": color, "stroke-width": ".55"})
        add_text(panel, 417, y, text, "p202-note")
    add_text(panel, 407, 174, "重要边界", "p202-head")
    notes = (
        "• 所有十字圈仅为既有 ObjectPlacement。",
        "• 不是给排水接口、粗装点或放线点。",
        "• 当前没有 IfcDistributionPort。",
        "• 当前没有供排水 System / 连通关系。",
        "• 蓝色投影保留既有几何，不代表新路径。",
        "• PVC110 维持 2 个产品，不拆分正式 IFC。",
    )
    for index, text in enumerate(notes):
        add_text(panel, 407, 184 + index * 8, text, "p202-small")
    ET.SubElement(panel, SVG + "line", {"x1": "407", "y1": "237", "x2": "493", "y2": "237", "class": "p202-panel-rule"})
    add_text(panel, 407, 249, "施工前停止条件", "p202-head")
    stops = (
        "1. 厂家粗装尺寸 / 接口中心未确认。",
        "2. 排水连接、坡度、通气及检修未确认。",
        "3. 厨房排水转接口需保留可检修条件。",
        "4. 无 Port / System 时不得冻结连接方案。",
    )
    for index, text in enumerate(stops):
        add_text(panel, 407, 260 + index * 9, text, "p202-small")
    add_text(panel, 407, 310, "底图可见性", "p202-head")
    add_text(panel, 407, 320, f"源 Sanitary Plan 可见洁具：{len(base_sanitary)} / 27", "p202-small")
    add_text(panel, 407, 328, "本候选用绿圈登记全部 27 件。", "p202-small")
    add_text(panel, 407, 336, "地漏 3IVq… 当前几何位于户型上方，待人审。", "p202-small")
    add_text(panel, 407, 348, f"IFC SHA {ifc_sha[:16]}…", "p202-small")
    add_text(panel, 407, 356, "状态：候选 / CONNECTION DATA MISSING", "p202-small")
    add_text(panel, 407, 383, "P-202 | 1:50 | 1 / 1", "p202-head")

    output_svg.parent.mkdir(parents=True, exist_ok=True)
    tree.write(output_svg, encoding="utf-8", xml_declaration=True)
    subprocess.run([
        sys.executable, str(render_script),
        "--input-svg", str(output_svg), "--output-pdf", str(output_pdf),
        "--proof-png", str(proof_png), "--report", str(render_report_path),
        "--width-mm", str(PAGE_WIDTH_MM), "--height-mm", str(PAGE_HEIGHT_MM),
        "--proof-dpi", "100",
    ], cwd=root_path, check=True)
    render_report = json.loads(render_report_path.read_text(encoding="utf-8"))
    metrics = raster_metrics(output_pdf, report_path.parent / "p202-raster-audit.pgm")
    marker_counts = Counter()
    for group in svg_root.iter(SVG + "g"):
        classes = set(group.get("class", "").split())
        for marker_name in ("sanitary", "waste", "drainage", "assembly"):
            if f"p202-{marker_name}-marker" in classes:
                marker_counts[marker_name] += 1
    qa = {
        "source_ifc_same_as_registry": True,
        "source_underlay_sanitary_count": len(base_sanitary),
        "sanitary_marker_count": marker_counts["sanitary"],
        "waste_marker_count": marker_counts["waste"],
        "drainage_marker_count": marker_counts["drainage"],
        "assembly_marker_count": marker_counts["assembly"],
        "drainage_projected_geometry_count": drain_geometry_count,
        "waste_projected_geometry_count": waste_geometry_count,
        "pvc110_product_count": plum_report["qa"]["pvc110_product_count"],
        "pvc110_branch_count": plum_report["qa"]["pvc110_branch_count"],
        "pvc110_world_geometry_unchanged": plum_report["qa"]["pvc110_world_geometry_unchanged"],
        "connectivity_status": plum_report["qa"]["distribution_data"]["status"],
        "connectivity_qa_passed": False,
        "off_plan_waste_terminal_ids": [
            product.GlobalId for product in waste if svg_point(world_origin_mm(product))[1] < 60.0
        ],
        "pdf_page_pass": render_report["pass"],
        "raster_metrics": metrics,
    }
    qa["candidate_mechanical_pass"] = all((
        qa["source_underlay_sanitary_count"] == 25,
        qa["sanitary_marker_count"] == 27,
        qa["waste_marker_count"] == 3,
        qa["drainage_marker_count"] == 16,
        qa["assembly_marker_count"] == 3,
        qa["drainage_projected_geometry_count"] == 16,
        qa["waste_projected_geometry_count"] == 3,
        qa["pvc110_product_count"] == 2,
        qa["pvc110_branch_count"] == 6,
        qa["pvc110_world_geometry_unchanged"],
        qa["connectivity_status"] == "data_missing",
        qa["off_plan_waste_terminal_ids"] == ["3IVqCnhGr51hY4LrOq_5G_"],
        qa["pdf_page_pass"], metrics["ink_ratio"] >= 0.002,
        metrics["ink_bbox_page_ratio"] >= 0.50,
    ))
    qa["construction_release_pass"] = False
    report = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "source": {
            "ifc": str(ifc_path), "ifc_sha256": ifc_sha,
            "underlay": str(underlay_path), "underlay_sha256": sha256(underlay_path),
            "registry": str(registry_path),
        },
        "outputs": {
            "svg": str(output_svg), "svg_sha256": sha256(output_svg),
            "pdf": str(output_pdf), "pdf_sha256": sha256(output_pdf),
            "proof_png": str(proof_png),
        },
        "qa": qa,
    }
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(qa, ensure_ascii=False))
    if not qa["candidate_mechanical_pass"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
