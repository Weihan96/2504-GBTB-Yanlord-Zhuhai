"""Include every in-scope queue product without promoting any approval.

Existing validated packages stay byte-identical. New candidates extract only
the original Body, not unapproved drawing expressions, into separate files.
"""
from pathlib import Path
import importlib.util
import json
import os
import shutil
import ifcopenshell
from review_product_package import ScopedCopy, selected_scope, sha256, body_fingerprint, placement

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "output/review/approved-product-library"
spec = importlib.util.spec_from_file_location("catalog_rules", ROOT / "pipeline/addons/highpoly_review_library/catalog_rules.py")
rules = importlib.util.module_from_spec(spec)
spec.loader.exec_module(rules)


def main():
    formal = ROOT / "2504 GBTB Yanlord Zhuhai.ifc"
    formal_hash = sha256(formal)
    assert formal_hash == "7a50b87e8f48a7c2bfbaa9f1dabdfca7155fa684b2325c66aa1ae15c8ab0a25c"
    source = ifcopenshell.open(str(formal))
    old = json.loads((OUT / "catalog.json").read_text())
    existing = {p["id"]: p for p in old["products"]}
    queue = json.loads((ROOT / "output/review/highpoly-types/highpoly-review-queue.json").read_text())
    products, excluded = [], []
    for item in queue["products"]:
        slug = item["slug"]
        if item.get("do_not_touch") or item.get("excluded"):
            excluded.append({"id": slug, "reason": item["status"]})
            continue
        folder = ROOT / item["folder"]
        approval_path = ROOT / item["approval"]["record"]
        approval = json.loads(approval_path.read_text())
        if slug in existing:
            p = existing[slug]
            p["insertion_content"] = "approved_3d_2d"
        else:
            manifest = json.loads((folder / "manifest.json").read_text())
            guid = manifest.get("representative_global_id", item["representative_global_id"])
            target = source.by_guid(guid)
            # Physical assembly ancestors are not part of this single product.
            # Forward ObjectPlacement dependencies still preserve their frame.
            scope = {e for e in selected_scope(source, target, []) if not e.is_a("IfcElement") or e == target}
            voids = list(getattr(target, "HasOpenings", ()))
            scope.update(r.RelatedOpeningElement for r in voids)
            package = ifcopenshell.file(schema=source.schema)
            clone = ScopedCopy(source, package, scope)
            for root in sorted(scope, key=lambda e: e.id()):
                clone.copy(root)
            for relation in voids:
                clone.copy(relation)
            clone.inverses()
            element = package.by_guid(guid)
            element.Representation.Representations = tuple(r for r in element.Representation.Representations if r.RepresentationIdentifier == "Body")
            destination = OUT / "candidate-packages" / (slug + ".ifc")
            destination.parent.mkdir(exist_ok=True)
            for style in package.by_type("IfcExternallyDefinedSurfaceStyle"):
                location = style.Location
                if location and "://" not in location and not Path(location).is_absolute():
                    dependency = (formal.parent / location).resolve()
                    copied = (destination.parent / location).resolve()
                    assert copied.is_relative_to(destination.parent.resolve()), "材质路径不能越出候选包目录"
                    assert dependency.is_file(), dependency
                    copied.parent.mkdir(parents=True, exist_ok=True)
                    if not copied.exists():
                        shutil.copy2(dependency, copied)
                    assert sha256(copied) == sha256(dependency)
            package.write(str(destination))
            reopened = ifcopenshell.open(str(destination))
            assert len([e for e in reopened.by_type("IfcElement") if not e.is_a("IfcOpeningElement")]) == 1
            assert len(reopened.by_type("IfcOpeningElement")) == len(voids)
            assert not reopened.by_type("IfcAnnotation")
            assert body_fingerprint(target) == body_fingerprint(reopened.by_guid(guid))
            assert (placement(target) == placement(reopened.by_guid(guid))).all()
            renders = json.loads((folder / "bonsai-review-manifest.json").read_text())["renders"]
            previews = {r["view"].replace("-elevation", ""): str(ROOT / r["path"]) for r in renders}
            assert all(Path(v).is_file() for v in previews.values()), slug
            p = {"id": slug, "name": manifest.get("display_name", item["type_name"]),
                 "global_id": guid, "ifc_path": os.path.relpath(destination, OUT), "ifc_sha256": sha256(destination),
                 "previews": previews, "single_svg": {view: str(folder / (view + ".svg")) for view in ("plan", "front", "side")},
                 "display_front_normal_project": [0, -1, 0], "display_front_basis": "project_axes_unreviewed",
                 "single_product_approval_status": approval["status"], "scene_approval_status": "pending",
                 "approval_label": "单品已通过、场景待验收" if approval["status"] == "approved" else ("部分视图通过、其余待验收" if approval["status"] == "partially_approved" else "尚未验收"),
                 "source_label": manifest.get("source_label_zh", item["source_kind"]),
                 "package_validation": "pass", "insertion_content": "body_only",
                 "candidate_review_index": str(folder / "index.html"),
                 "preview_note": "原始 Bonsai Body 相机图；待审二维 SVG 单独列出，未嵌入 IFC",
                 "approval_record": {"path": str(approval_path), "sha256": sha256(approval_path)}}
        p["review_category"] = rules.classify(p["single_product_approval_status"], p["scene_approval_status"])
        p["category_label"] = rules.CATEGORIES[p["review_category"]]
        p["formal_ifc_write_authorized"] = False
        rules.validate_entry(p, lambda value: (OUT / value).resolve(), sha256)
        products.append(p)
        print("CATALOG_PRODUCT", slug, p["category_label"], flush=True)
    assert sha256(formal) == formal_hash
    result = {"schema_version": 3, "products": products, "excluded": excluded,
              "formal_sha256_unchanged": formal_hash,
              "inclusion_authority": "用户要求将尚未通过验证的模型全部加入图库并标注验收类别；纳入不代表批准"}
    (OUT / "all-review-catalog.json").write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n")


if __name__ == "__main__":
    main()
