"""Build a read-only library from explicitly accepted, verified pure packages.

Thumbnail vectors are product LINEWORK cropped from the newly regenerated
Bonsai scene SVG. They are not a new drawing algorithm or old proxy previews.
"""
import copy
import json
import os
from pathlib import Path
import subprocess
import xml.etree.ElementTree as ET
import ifcopenshell
import review_product_package as pkg
from pure_product_package import graph_from_json

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "output/review/approved-product-library"
PRODUCTS = ROOT / "output/review/highpoly-types"
FORMAL = ROOT / "2504 GBTB Yanlord Zhuhai.ifc"
FORMAL_HASH = "7a50b87e8f48a7c2bfbaa9f1dabdfca7155fa684b2325c66aa1ae15c8ab0a25c"
# Whole-product scene acceptance confirmed in the user's review history.
# Other single-view/partial approvals are not silently promoted to this list.
ACCEPTED = {
    "bed01": "Baxter Casablanca 180 / BED01",
    "wd03": "Poliform Senzafine / WD03",
    "trap01": "Geberit 151.116.11.1 / TRAP01",
    "miamisoft-e09": "Baxter Miami Soft / E09",
    "wd02": "Poliform Senzafine glass wardrobe / WD02",
    "sxb010": "Hunter Douglas blinds / SXB010",
    "gessi316-54294": "Gessi316 Meccanica / 54294",
    "geberit-duofix-sigma-224-212": "Geberit Duofix Sigma / 224.212",
}
SINGLE_ACCEPTED = {
    "hima01": "Poliform Hima / HIMA01",
    "geberit-154-446-ks-1": "Geberit shower channel / 154.446.KS.1",
    "bed02": "Baxter Viktor 160 × 200 / BED02",
    "gessi316-54145": "Gessi316 Meccanica / 54145 G000",
    "street-h": "antoniolupi STREET-H / 支撑子构件",
    "tab02": "RODA Bernardo 367 / TAB02",
    "falper-sorgente": "Falper Sorgente / WFB",
}
APPROVAL_LABELS = {"approved": "单品及场景已验收", "pending": "单品已通过、场景待验收"}
SOURCE_LABELS = {
    "bed01": "官方家族 DWG · 非项目精确配置",
    "wd03": "基于原始高模几何生成的简化图纸表达",
    "wd02": "基于原始高模几何生成的简化图纸表达",
    "trap01": "官方 DWG · 按项目安装长度配置",
    "miamisoft-e09": "官方 E09 DWG · 保留项目高模差异",
    "sxb010": "高模几何简化 · 蓝色仅为审核高亮",
    "gessi316-54294": "官方 DWG 去纹审核简化表达",
    "geberit-duofix-sigma-224-212": "官方 224.212.00.2 原生 DWG",
    "hima01": "基于原始高模几何生成的简化图纸表达",
    "geberit-154-446-ks-1": "官方 154.446.KS.1 原生 DWG",
    "bed02": "官方 160×200 家族 DWG · 非项目精确配置",
    "gessi316-54145": "基于官方 54145 G000 DWG 的简化表达",
    "street-h": "支撑子构件 · 基于原始高模几何生成的简化图纸表达",
    "tab02": "基于原始高模几何生成的简化图纸表达",
    "falper-sorgente": "官方 WFB 家族 CAD · 非项目施工图",
}
SVG = "http://www.w3.org/2000/svg"
ET.register_namespace("", SVG)


def authorized_scope(scope, authorization):
    """Inclusion authority is distinct from scene approval and technical pass."""
    assert set(scope["whole_product_accepted"]) == set(ACCEPTED)
    assert set(scope["single_product_accepted_scene_not_final"]) == set(SINGLE_ACCEPTED)
    assert authorization["pure_product_ifc_write_authorized"] is True
    assert authorization["library_inclusion_authorized"] is True
    assert authorization["formal_ifc_write_authorized"] is False
    assert authorization["scene_approval_status"] == "pending"
    assert authorization["single_product_approval_status"] == "approved"
    assert authorization["user_response"] == "是"
    assert set(authorization["products"]) == set(SINGLE_ACCEPTED)
    assert not set(ACCEPTED) & set(SINGLE_ACCEPTED)
    return {**ACCEPTED, **SINGLE_ACCEPTED}


def approval_fields(key):
    assert key in ACCEPTED or key in SINGLE_ACCEPTED
    scene = "approved" if key in ACCEPTED else "pending"
    return {"single_product_approval_status": "approved",
            "scene_approval_status": scene, "approval_label": APPROVAL_LABELS[scene]}


def validate_approval_fields(product):
    for field, expected in approval_fields(product["id"]).items():
        assert product.get(field) == expected, (product["id"], field, expected)


def relative(path):
    return os.path.relpath(Path(path).resolve(), OUT)


def cropped_linework(source, guid, destination):
    root = ET.parse(source).getroot()
    lines = [e for e in root.iter() if f"GlobalId-{guid}" in e.get("class", "").split()]
    assert lines and all(e.tag == f"{{{SVG}}}line" for e in lines), "Unexpected non-line SVG representation"
    parents = {e: p for p in root.iter() for e in p}
    for line in lines:
        node = line
        while node is not root:
            assert not node.get("transform"), "Thumbnail extraction needs explicit transformed-SVG support"
            node = parents[node]
    xs = [float(e.get(k)) for e in lines for k in ("x1", "x2")]
    ys = [float(e.get(k)) for e in lines for k in ("y1", "y2")]
    pad = max(max(xs) - min(xs), max(ys) - min(ys)) * .065 + .5
    x, y = min(xs) - pad, min(ys) - pad
    width, height = max(xs) - x + pad, max(ys) - y + pad
    cropped = ET.Element(f"{{{SVG}}}svg", {"version": "1.1", "viewBox": f"{x} {y} {width} {height}",
        "width": f"{width}mm", "height": f"{height}mm", "data-source": "Bonsai scene product LINEWORK crop",
        "data-source-svg-sha256": pkg.sha256(source), "data-annotation-guid": guid})
    # Keep the source style rules and every actual approved line unchanged.
    for item in root:
        if item.tag == f"{{{SVG}}}defs":
            cropped.append(copy.deepcopy(item))
    for line in lines:
        cropped.append(copy.deepcopy(line))
    ET.ElementTree(cropped).write(destination, encoding="utf-8", xml_declaration=True)
    png = destination.with_suffix(".png")
    subprocess.run(["/Applications/Inkscape.app/Contents/MacOS/inkscape", str(destination),
        "--export-area-page", "--export-background=white", "--export-background-opacity=1",
        "--export-width=640", f"--export-filename={png}"], check=True, capture_output=True)
    return {"source_svg": pkg.record(source), "annotation_global_id": guid, "line_count": len(lines),
            "single_svg": pkg.record(destination), "preview": pkg.record(png), "line_geometry_changed": False}


def build():
    assert pkg.sha256(FORMAL) == FORMAL_HASH
    OUT.mkdir(parents=True, exist_ok=True)
    authorization_path = OUT / "migration-authorization.json"
    authorization = json.loads(authorization_path.read_text())
    scope = json.loads((OUT / "acceptance-scope.json").read_text())
    authorized = authorized_scope(scope, authorization)
    products, pending = [], []
    for key, name in authorized.items():
        directory = PRODUCTS / key / "product-package"
        handoff_path = directory / "handoff.json"
        if not handoff_path.is_file():
            pending.append({"id": key, "reason": "migration_not_complete"})
            continue
        handoff = json.loads(handoff_path.read_text())
        complete = (handoff.get("verdict") == "pass" or handoff.get("validation_verdict") == "pass"
                    or handoff.get("status") in ("complete", "complete_verified"))
        if not complete:
            pending.append({"id": key, "reason": "migration_validation_pending"})
            continue
        if key in SINGLE_ACCEPTED:
            assert handoff.get("scene_approval_status") == "pending", "Technical migration cannot approve scenes"
        paths = list(directory.glob("*-product.ifc"))
        assert len(paths) == 1, (key, paths)
        path = paths[0]
        model = ifcopenshell.open(str(path))
        assert len(model.by_type("IfcElement")) == 1 and not model.by_type("IfcAnnotation")
        assert not model.by_type("IfcGroup")
        assert not [s for s in model.by_type("IfcSpatialElement") if s.Representation]
        target = model.by_type("IfcElement")[0]
        reps = target.Representation.Representations
        assert len([r for r in reps if r.RepresentationIdentifier != "Body"]) == 3
        recipe_path = directory / "scene-recipe.json"
        recipe = json.loads(recipe_path.read_text())
        assert recipe["target_global_id"] == target.GlobalId
        assert recipe["linework_geometry_included"] is False and recipe["body_geometry_included"] is False
        scene_metadata = graph_from_json(recipe["graph"])
        front_camera = scene_metadata.by_guid(recipe["views"]["front"]["drawing_guid"])
        front_normal = pkg.placement(front_camera)[:3, 2]
        assert abs(float(front_normal[2])) < .01, "Front camera is not horizontal"
        approval = ROOT / f"pipeline/decisions/{key}-drawing-approval.json"
        assert approval.is_file()
        target_dir = OUT / "previews" / key
        target_dir.mkdir(parents=True, exist_ok=True)
        previews, vectors, scene_vectors, evidence = {}, {}, {}, {}
        for view in ("plan", "front", "side"):
            scenes = list(directory.glob(f"*-SCENE-{view.upper()}.svg"))
            assert len(scenes) == 1, (key, view, scenes)
            vector = target_dir / f"{view}.svg"
            evidence[view] = cropped_linework(scenes[0], recipe["views"][view]["annotation_guid"], vector)
            vectors[view] = relative(vector)
            scene_vectors[view] = relative(scenes[0])
            previews[view] = relative(vector.with_suffix(".png"))
        iso = PRODUCTS / key / "bonsai-camera-iso.png"
        assert iso.is_file(), "Missing actual Body camera preview"
        previews["iso"] = relative(iso)
        products.append({"id": key, "name": name, "global_id": target.GlobalId,
            "source_label": SOURCE_LABELS[key],
            "display_front_normal_project": front_normal.tolist(),
            "display_front_basis": "approved_scene_camera" if key in ACCEPTED else "candidate_scene_camera_pending_review",
            **approval_fields(key), "package_validation": "pass",
            "ifc_path": relative(path), "ifc_sha256": pkg.sha256(path), "ifc_bytes": path.stat().st_size,
            "body_fingerprint": pkg.body_fingerprint(target),
            "single_svg": vectors, "scene_svg": scene_vectors, "previews": previews,
            "thumbnail_evidence": evidence, "iso_evidence": {**pkg.record(iso),
                "kind": "existing_actual_Body_saved_camera_render", "basis": "migration preserves Body fingerprint"},
            "approval_record": pkg.record(approval), "validation_record": pkg.record(handoff_path),
            "recipe": pkg.record(recipe_path), "ifc_write_allowed": False,
            "library_authorization": pkg.record(authorization_path) if key in SINGLE_ACCEPTED else None})
    catalog = {"schema_version": 2, "products": products, "pending": pending,
        "count_policy": "whole-product acceptance plus explicitly authorized single-approved products; never a fixed 34",
        "approval_summary": {"scene_approved": sum(p["scene_approval_status"] == "approved" for p in products),
                             "scene_pending": sum(p["scene_approval_status"] == "pending" for p in products)},
        "authorized_product_count": len(authorized),
        "display_policy": "Blender copies only; source IFC retains project coordinates",
        "formal_ifc_sha256": FORMAL_HASH, "formal_write_allowed": False,
        "partial_approvals_not_promoted": True}
    (OUT / "catalog.json").write_text(json.dumps(catalog, ensure_ascii=False, indent=2) + "\n")
    assert pkg.sha256(FORMAL) == FORMAL_HASH
    print(json.dumps({"ready": [p["id"] for p in products], "pending": pending}, ensure_ascii=False))


if __name__ == "__main__":
    build()
