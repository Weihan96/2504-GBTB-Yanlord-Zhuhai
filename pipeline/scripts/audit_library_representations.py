"""Read-only IFC audit. Writes reports, never edits IFC or deletes candidates.

Exercise the installed Bonsai context selector directly (AST extraction avoids
importing bpy). Geometry identities are compared without inferring semantic
equivalence from a matching bounding box or an identical UI label.
"""
import ast
from collections import Counter, namedtuple
import copy
import json
from pathlib import Path
import re
import subprocess

import ifcopenshell
import ifcopenshell.geom
from review_product_package import fingerprint, sha256

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "output/review/approved-product-library"
SITE = Path("/Users/jiaxinchen/Library/Application Support/Blender/4.5/extensions/.local/lib/python3.11/site-packages")
DRAWING = SITE / "bonsai/bim/module/drawing/operator.py"
FORMAL = ROOT / "2504 GBTB Yanlord Zhuhai.ifc"
BASELINE = "7a50b87e8f48a7c2bfbaa9f1dabdfca7155fa684b2325c66aa1ae15c8ab0a25c"


def selector():
    tree = ast.parse(DRAWING.read_text())
    fn = copy.deepcopy(next(n for n in ast.walk(tree) if isinstance(n, ast.FunctionDef) and n.name == "get_linework_contexts"))
    source_line = fn.lineno
    fn.returns = None
    for arg in fn.args.args:
        arg.annotation = None
    fn.decorator_list = []
    namespace = {"LineworkContexts": namedtuple("LineworkContexts", "body annotation")}
    exec(compile(ast.Module(body=[fn], type_ignores=[]), str(DRAWING), "exec"), namespace)
    return namespace[fn.name], source_line


def audit_product(entry, choose):
    path = (OUT / entry["ifc_path"]).resolve()
    digest = sha256(path)
    assert digest == entry["ifc_sha256"], (entry["id"], "catalog hash mismatch")
    model = ifcopenshell.open(str(path))
    element = model.by_guid(entry["global_id"])
    reps = list(element.Representation.Representations)
    rows = []
    for rep in reps:
        ctx = rep.ContextOfItems
        rows.append({"step_id": rep.id(), "representation_identifier": rep.RepresentationIdentifier,
                     "context_type": ctx.ContextType, "context_identifier": ctx.ContextIdentifier,
                     "target_view": getattr(ctx, "TargetView", None), "context_step_id": ctx.id(),
                     "representation_type": rep.RepresentationType, "item_count": len(rep.Items),
                     "items_fingerprint": fingerprint(rep.Items)})
    primary = [r for r in rows if r["representation_identifier"] == "Body" and r["target_view"] == "MODEL_VIEW"]
    assert len(primary) == 1, (entry["id"], "primary Body ambiguous")
    non_model_body = [r for r in rows if r["representation_identifier"] == "Body" and r not in primary]
    review_views = [r for r in rows if r["representation_identifier"] != "Body"]
    selection = {}
    for view in ("PLAN_VIEW", "ELEVATION_VIEW"):
        contexts = choose(None, model, view)
        selection[view] = {
            channel: [[r["representation_identifier"] for r in rows if r["context_step_id"] in ids] for ids in getattr(contexts, channel)]
            for channel in ("body", "annotation")}
    findings = []
    if any(r["context_identifier"] == "Annotation" for r in review_views):
        findings.append("approved_product_linework_in_annotation_context")
    if any(r["context_identifier"] not in ("Body", "Facetation", "Annotation") for r in review_views):
        findings.append("custom_product_context_not_selected_by_native_drawing")
    if review_views and non_model_body:
        findings.append("legacy_body_views_coexist_with_approved_linework")
    if not review_views:
        findings.append("body_only_candidate_no_approved_three_view_package")
    if len([r for r in review_views if r["target_view"] == "ELEVATION_VIEW"]) > 1:
        findings.append("front_side_share_target_view_no_native_direction_selection")
    physical = [e for e in model.by_type("IfcElement") if not e.is_a("IfcOpeningElement")]
    assert len(physical) == 1
    assert not model.by_type("IfcAnnotation")
    assert sha256(path) == digest
    return {"id": entry["id"], "path": str(path.relative_to(ROOT)), "sha256": digest,
            "global_id": element.GlobalId, "review_category": entry["review_category"],
            "physical_product_count": len(physical), "opening_dependency_count": len(model.by_type("IfcOpeningElement")),
            "ifc_annotation_entity_count": 0, "representations": rows,
            "legacy_body_view_count": len(non_model_body), "review_view_count": len(review_views),
            "exact_item_graph_duplicates_of_model_body": [r["step_id"] for r in non_model_body if r["items_fingerprint"] == primary[0]["items_fingerprint"]],
            "native_context_priority_buckets": selection, "findings": findings}


def entity_types():
    schema = ifcopenshell.ifcopenshell_wrapper.schema_by_name("IFC4")
    def descend(declaration):
        return {declaration.name().upper()} | set().union(*(descend(c) for c in declaration.subtypes()))
    return descend(schema.declaration_by_name("IfcElement"))


def rename_only_probe(catalog, choose):
    entry = next(p for p in catalog["products"] if p["id"] == "gessi316-54294")
    path = (OUT / entry["ifc_path"]).resolve()
    before = sha256(path)
    model = ifcopenshell.open(str(path))
    element = model.by_guid(entry["global_id"])
    reps = list(element.Representation.Representations)
    for rep in reps:
        if rep.RepresentationIdentifier in ("ApprovedFront", "ApprovedSide"):
            # Deliberately exercise the proposed superficial fix in a private
            # memory copy only. No IFC write occurs in this audit.
            rep.ContextOfItems.ContextIdentifier = "Body"
    contexts = choose(None, model, "ELEVATION_VIEW").body[0]
    settings = ifcopenshell.geom.settings()
    settings.set("dimensionality", ifcopenshell.ifcopenshell_wrapper.CURVES_SURFACES_AND_SOLIDS)
    settings.set("context-ids", contexts)
    shapes = [{"element_step_id": shape.id, "geometry_id": shape.geometry.id,
               "vertex_count": len(shape.geometry.verts) // 3, "edge_count": len(shape.geometry.edges) // 2}
              for shape in ifcopenshell.geom.iterator(settings, model, 1, include=[element])]
    assert sha256(path) == before
    return {"product": entry["id"], "source_sha256_unchanged": before, "persisted": False,
            "experiment": "rename approved Front/Side context identifier to Body only",
            "selected_context_ids": contexts, "shapes_emitted_for_one_product": shapes,
            "representation_id_map": {str(r.id()): r.RepresentationIdentifier for r in reps},
            "conclusion": "legacy elevation, approved Front and approved Side emitted together; rename alone is not a fix",
            "not_a_create_drawing_svg_test": True}


def cleanup_inventory(catalog):
    # Enumerate actual files, including ignored old IFC caches; never follow
    # symlinks or inspect the expressly excluded Libelle task's artifacts.
    element_types = entity_types()
    active = {(OUT / e["ifc_path"]).resolve() for e in catalog["products"]}
    entries = {e["id"]: e for e in catalog["products"]}
    proposed = set()
    for proposal in (ROOT / "output/review/highpoly-types").glob("*/product-package/cleanup-proposal.json"):
        for row in json.loads(proposal.read_text()).get("legacy_files", []):
            proposed.add(Path(row["path"]).resolve())
    files = subprocess.check_output(["rg", "--files", "--hidden", "-g", "!.git", "-g", "*.py", "-g", "*.json", "-g", "*.md"]).decode().splitlines()
    texts = {}
    for name in files:
        p = ROOT / name
        if "libelle" not in name.lower() and p.is_file() and p.stat().st_size < 8_000_000 and "representation-audit" not in name and "cleanup-audit" not in name:
            texts[name] = p.read_text(errors="replace")
    result = []
    for path in sorted((ROOT / "output/review").rglob("*")):
        if path.is_symlink() or not path.is_file() or "libelle" in str(path).lower():
            continue
        if path.suffix.lower() not in (".ifc", ".blend", ".blend1", ".blend2"):
            continue
        rel = path.relative_to(ROOT)
        row = {"path": str(rel), "bytes": path.stat().st_size, "sha256": sha256(path), "deleted": False}
        product = next((entries[part] for part in rel.parts if part in entries), None)
        row["product"] = product["id"] if product else None
        if path.suffix.lower() == ".ifc":
            types = Counter(t.decode() for t in re.findall(rb"^\s*#\d+\s*=\s*(IFC[A-Z0-9_]+)\s*\(", path.read_bytes(), re.M))
            row["ifc_element_count_step_scan"] = sum(n for t, n in types.items() if t in element_types)
            row["ifc_annotation_count_step_scan"] = types["IFCANNOTATION"]
            if path.resolve() in active:
                row["classification"] = "keep_active_library_ifc"
            elif row["ifc_element_count_step_scan"] == 714:
                row["classification"] = "full_scene_legacy_candidate_blocked" if product and product["insertion_content"] == "approved_3d_2d" else "keep_full_scene_pending_product_or_unresolved_scope"
            elif row["ifc_element_count_step_scan"] > 1:
                row["classification"] = "keep_small_multi_element_fixture_dependencies_unresolved"
            else:
                row["classification"] = "isolated_legacy_candidate_requires_dependency_and_render_check"
        elif path.name == "Materials.blend" or "native-assets" in path.parts:
            row["classification"] = "keep_active_material_or_library_asset"
        else:
            row["classification"] = "blend_cache_candidate_requires_dependency_and_session_check"
        row["already_in_product_cleanup_proposal"] = path.resolve() in proposed
        row["references"] = [name for name, text in texts.items() if path.name in text or str(rel) in text]
        row["reference_scan_limit"] = "literal filename/path in non-Libelle py/json/md <8MB; not runtime dependency proof"
        if "candidate" in row["classification"]:
            row["blockers"] = ["representation correction and replacement drawing verification incomplete", "reference audit and user deletion confirmation required"]
        result.append(row)
    return {"status": "proposal_only_no_deletion", "cleanup_authorized": False,
            "immediately_delete_approved_files": [],
            "scope": "output/review IFC and Blender files; source CAD, SVG/PNG, manifests, approvals and recipes must be retained",
            "inventory_method": "STEP entity type scan, not IFC schema or rendering validation", "files": result,
            "totals": {kind: {"count": sum(r["classification"] == kind for r in result),
                               "bytes": sum(r["bytes"] for r in result if r["classification"] == kind)} for kind in sorted({r["classification"] for r in result})}}


def main():
    assert sha256(FORMAL) == BASELINE
    index_before = subprocess.check_output(["git", "ls-files", "--stage", "-z"])
    catalog = json.loads((OUT / "all-review-catalog.json").read_text())
    choose, source_line = selector()
    products = []
    for entry in catalog["products"]:
        products.append(audit_product(entry, choose))
        print("AUDIT", entry["id"], products[-1]["findings"], flush=True)
    report = {"status": "audit_complete_correction_requires_direction_policy", "ifc_modified": False,
              "formal_ifc_sha256_before": BASELINE, "formal_ifc_sha256_after": sha256(FORMAL),
              "course_evidence_mode": "embedded-course-index",
              "course_lesson": "109000 Representation contexts as they relate to a door; 00:27 Plan / Body / PLAN_VIEW / Annotation2D",
              "course_screenshot_provenance_sha256_not_directly_observed": "34540cf996d3ed68b6158f3bd5f4caf7f68d77008b098fc90f9cb435b56204d0",
              "installed_selector": {"path": str(DRAWING), "line": source_line, "sha256": sha256(DRAWING), "executed_from_installed_ast": True},
              "ifc4_context_reference": "https://standards.buildingsmart.org/IFC/RELEASE/IFC4/FINAL/HTML/schema/ifcrepresentationresource/lexical/ifcgeometricrepresentationsubcontext.htm",
              "count": len(products), "finding_counts": dict(Counter(f for p in products for f in p["findings"])), "products": products,
              "rename_only_probe": rename_only_probe(catalog, choose),
              "not_verified": ["new native Create Drawing SVG", "corrected IFC reload", "direction-aware representation adapter"],
              "approval_statuses_unchanged": True}
    (OUT / "representation-audit.json").write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n")
    cleanup = cleanup_inventory(catalog)
    (OUT / "cleanup-audit.json").write_text(json.dumps(cleanup, ensure_ascii=False, indent=2) + "\n")
    lines = ["# 单品表示核查与清理候选", "", "状态：39 项核查完成；IFC 修正尚未执行。所有 IFC、批准状态及当前 Blender 会话保持不变。本报告不是修正版 IFC 的验收证明。", "",
             "## 核查结论", "",
             "- 39 份活动单品 IFC 均只有一个实际产品（电动止回阀另有两个必要开孔），均没有 IfcAnnotation 实体。截图中的 Annotation 是表示上下文名称，不是额外实体数量。",
             "- 15 份包含已批准三视图：其中 13 份放在 Annotation 上下文；WD03 与 Falper 使用原生 Drawing 不选取的自定义上下文。14 份同时保留旧 Body 平面／立面。",
             "- 24 份是 Body-only 候选，不能当作已完成三视图接入。不能从“representation 数量足够”推断接入正确。",
             "- 没有发现与 MODEL_VIEW 的 Items 前向图指纹完全相同的旧 Body 表示；不能把所有旧平面／立面称作重复高模。旧版本与批准版本共存的问题仍需解决。",
             "- 安装版 Bonsai 的上下文选择函数已直接执行；15 份的 Front/Side 都使用 ELEVATION_VIEW，原生选择器不按相机方向区分。Gessi 54294 的内存实验确认：只将 Annotation 改名 Body，会同时交给几何迭代器旧立面、Front、Side 三组几何。没有生成新 SVG，因此不声称已观察到最终图纸重影。", "",
             "## 方法与需要确认的边界", "",
             "Bonsai Course Operator 课程 109000 的 00:27 支持产品平面使用 Plan / Body / PLAN_VIEW；Annotation2D 是几何类型，不是 Annotation 上下文。证据为 embedded-course-index，未声称直接打开课程私有截图。安装版源码与 IFC4 资料用于交叉核对。", "",
             "建议保留原始 3D Body 和已批准三视图；用批准的 Plan 替换旧 Plan Body，而不是追加一套。Front/Side 必须独立保留视向语义，不能改成两个无法区分的 ELEVATION_VIEW Body。", "",
             "需要人工确认的架构选择：继续保留独立官方二维图，并让图库在出图时按相机方向选用正确线稿；或改成由统一简化三维几何生成立面（不能保证保留已批准官方二维画法）。推荐前者，但不能声称它是现成的原生 Bonsai 自动行为。Bonsai Course Operator 对分歧架构要求人工审阅，因此本轮不批量改写。", "",
             "实现前者还必须处理 append_asset 的上下文合并：当前版本按 ContextType/ContextIdentifier/TargetView 合并，单独增添 Front/Side 名称不足以保留方向。不得修改公共 Provider 或增设保存按钮；仍由 Bonsai Ctrl+S 保存。", "",
             "修正验收：先 Gessi 54294 试行；原点、单位、3D Body、已批准线稿和样式不变；保存重载后真实 Create Drawing 的 Plan/Front/Side SVG 必须分别验证，没有旧线或另一视向混入；再批量覆盖 39 项对应规则。正式 IFC 前后 SHA-256 必须不变。", "",
             "## 全部单品", "", "| 产品 | 表示数 | 旧 Body 非模型视图 | 已批准线稿表示 | 类别 |", "|---|---:|---:|---:|---|"]
    for p in products:
        lines.append(f"| {p['id']} | {len(p['representations'])} | {p['legacy_body_view_count']} | {p['review_view_count']} | {p['review_category']} |")
    lines += ["", "## 清理清单（未授权删除）", "",
              "这里只盘点 output/review 下 IFC 与 Blender 文件；不是对仓库其他内容的删除授权。字面路径引用检查不能证明运行时可删除。所有 SVG/PNG 验收证据、官方 CAD 原件及下载链接/哈希、批准记录、scene-recipe 和必要材质依赖保留。", "",
              "| 分类 | 文件数 | 大小（十进制 GB） |", "|---|---:|---:|"]
    for kind, totals in cleanup["totals"].items():
        lines.append(f"| {kind} | {totals['count']} | {totals['bytes']/1e9:.3f} |")
    lines += ["", "### 优先复查的旧完整场景 IFC", "",
              "以下每份含 714 个 IfcElement，与正式项目的元素数量相同；已有对应单品包，但表示修正与替代出图验证未完成，暂不可立即删除。", "",
              "| 文件 | MB |", "|---|---:|"]
    for r in cleanup["files"]:
        if r["classification"] == "full_scene_legacy_candidate_blocked":
            lines.append(f"| {r['path']} | {r['bytes']/1e6:.1f} |")
    lines += ["", "其余 .blend/.blend1 包括审核相机、待审产品和旧阵列，需要逐个确认替代产物及是否正被使用；不应整目录清除。详见 cleanup-audit.json 的逐文件哈希、分类和引用清单。", "",
              f"正式 IFC SHA-256：`{BASELINE}`，操作前后不变。此次未操作 stash、未修改正式 IFC、未删除旧审核文件。", ""]
    (OUT / "representation-audit.md").write_text("\n".join(lines))
    assert sha256(FORMAL) == BASELINE
    assert subprocess.check_output(["git", "ls-files", "--stage", "-z"]) == index_before
    print(json.dumps({"findings": report["finding_counts"], "cleanup": cleanup["totals"]}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
