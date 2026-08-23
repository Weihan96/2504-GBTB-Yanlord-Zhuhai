#!/usr/bin/env python3
"""Audit the approval-gated high-poly drawing review queue.

This audit is read-only with respect to IFC files. It verifies that every
in-scope product has a self-contained human-review package, real Bonsai camera
evidence, and project drawings which retain surrounding project context.
"""

from __future__ import annotations

import argparse
import hashlib
import html
import json
import subprocess
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
INVENTORY = ROOT / "pipeline/decisions/highpoly-product-inventory.json"
FORMAL_IFC = ROOT / "2504 GBTB Yanlord Zhuhai.ifc"
EXPECTED_FORMAL_HASH = (
    "7a50b87e8f48a7c2bfbaa9f1dabdfca7155fa684b2325c66aa1ae15c8ab0a25c"
)
DEFAULT_OUTPUT = ROOT / "output/review/highpoly-types/highpoly-review-queue.json"
STASH_RECOVERY_AUDIT = ROOT / "pipeline/decisions/highpoly-stash-recovery-audit.json"
FALLBACK_SOURCE_KIND = "geometry_derived_simplified_proxy"
FALLBACK_SOURCE_LABEL_ZH = "基于原始高模几何生成的简化图纸表达"
OFFICIAL_SOURCE_KINDS = {
    "native_dwg",
    "native_dxf",
    "official_native_dwg",
    "official_native_dxf",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def run_git(*arguments: str) -> str:
    result = subprocess.run(
        ["git", *arguments],
        cwd=ROOT,
        check=True,
        capture_output=True,
        text=True,
    )
    return result.stdout.strip()


def relative(path: Path) -> str:
    return path.relative_to(ROOT).as_posix()


def drawing_source(manifest: dict[str, Any]) -> tuple[str, bool]:
    source_kind = manifest.get("source_kind")
    drawing = manifest.get("drawing_source") or {}
    official = manifest.get("official_reference") or {}
    if not source_kind:
        source_kind = drawing.get("source_kind")
    if not source_kind:
        source_kind = official.get("source_kind")
    source_kind = source_kind or "unresolved"

    official_cad_used = manifest.get("official_cad_used")
    if not isinstance(official_cad_used, bool):
        official_cad_used = drawing.get("official_cad_used")
    if not isinstance(official_cad_used, bool):
        official_cad_used = source_kind in {
            "native_dwg",
            "native_dxf",
            "official_native_dwg",
            "official_native_dxf",
        }
    return source_kind, official_cad_used


def context_output_paths(context: dict[str, Any]) -> list[Path]:
    paths: list[Path] = []
    for view in context.get("views", []):
        value = view.get("styled_svg") or view.get("output")
        if value:
            paths.append(ROOT / value)
    return paths


def approval_status_is_pending(status: str) -> bool:
    return status in {"pending", "visual_review_pending"}


def source_label_zh(manifest: dict[str, Any]) -> str | None:
    drawing = manifest.get("drawing_source") or {}
    return manifest.get("source_label_zh") or drawing.get("source_label_zh")


def candidate_views(candidate: dict[str, Any]) -> list[dict[str, Any]]:
    views = candidate.get("views") or {}
    return list(views.values()) if isinstance(views, dict) else list(views)


def source_provenance_gate(
    *,
    completed: bool,
    folder: Path,
    manifest: dict[str, Any],
    candidate: dict[str, Any],
    context: dict[str, Any],
    source_kind: str,
    official_cad_used: bool,
) -> dict[str, Any]:
    access_path = folder / "official-source/source-access-record.json"
    access = load_json(access_path) if access_path.is_file() else {}
    views = candidate_views(candidate)

    if completed:
        falper_source = ROOT / "pipeline/decisions/falper-sorgente-official-source-manifest.json"
        checks = {
            "manifest_source_is_official_native_cad": source_kind in OFFICIAL_SOURCE_KINDS,
            "official_cad_used": official_cad_used is True,
            "context_source_is_official_native_cad": context.get("source_kind") in OFFICIAL_SOURCE_KINDS,
            "official_source_manifest_exists": falper_source.is_file(),
        }
        mode = "approved_official_native_cad"
    elif official_cad_used:
        official_paths_present = bool(views) and all(
            bool(view.get("official_native_dwg_paths_mm") or view.get("official_native_dxf_paths_mm"))
            for view in views
        )
        context_blue = bool(
            context.get("blue_line_present") is True
            or context.get("blue_line_top_layer_with_white_mask") is True
        )
        checks = {
            "manifest_source_is_official_native_cad": source_kind in OFFICIAL_SOURCE_KINDS,
            "candidate_source_is_official_native_cad": candidate.get("source_kind") in OFFICIAL_SOURCE_KINDS,
            "candidate_marks_official_cad_used": candidate.get("official_cad_used") is True,
            "candidate_has_native_paths_for_all_views": official_paths_present,
            "context_source_is_official_native_cad": context.get("source_kind") in OFFICIAL_SOURCE_KINDS,
            "context_has_blue_official_top_layer": context_blue,
            "source_access_record_exists": access_path.is_file(),
            "source_access_record_passes": access.get("pass") is True,
        }
        mode = "official_native_cad"
    else:
        official_project_paths_absent = bool(views) and all(
            not view.get("official_cad_paths_mm") for view in views
        )
        separate_family_references_not_substituted = all(
            view.get("official_family_paths_used_as_project_representation") is not True
            for view in views
        )
        context_blue_absent = (
            context.get("blue_line_present") is not True
            and context.get("blue_product_cad_line_present") is not True
            and context.get("blue_line_top_layer_with_white_mask") is not True
        )
        checks = {
            "manifest_uses_fallback_kind": source_kind == FALLBACK_SOURCE_KIND,
            "manifest_uses_exact_fallback_label": source_label_zh(manifest) == FALLBACK_SOURCE_LABEL_ZH,
            "manifest_marks_official_cad_unused": manifest.get("official_cad_used") is False,
            "candidate_uses_fallback_kind": candidate.get("source_kind") == FALLBACK_SOURCE_KIND,
            "candidate_uses_exact_fallback_label": candidate.get("source_label_zh") == FALLBACK_SOURCE_LABEL_ZH,
            "candidate_marks_official_cad_unused": candidate.get("official_cad_used") is False,
            "candidate_has_no_project_official_paths": official_project_paths_absent,
            "separate_family_references_not_substituted": separate_family_references_not_substituted,
            "context_uses_fallback_kind": context.get("source_kind") == FALLBACK_SOURCE_KIND,
            "context_uses_exact_fallback_label": context.get("source_label_zh") == FALLBACK_SOURCE_LABEL_ZH,
            "context_marks_official_cad_unused": context.get("official_cad_used") is False,
            "context_has_no_blue_project_cad_line": context_blue_absent,
            "source_access_record_exists": access_path.is_file(),
            "source_access_record_passes": access.get("pass") is True,
        }
        mode = "geometry_derived_fallback"
    return {
        "mode": mode,
        "source_kind": source_kind,
        "source_label_zh": source_label_zh(manifest),
        "official_cad_used": official_cad_used,
        "checks": checks,
        "pass": all(checks.values()),
    }


def audit_product(product: dict[str, Any], formal_hash: str) -> dict[str, Any]:
    status = product["status"]
    folder = ROOT / product["folder"]
    slug = folder.name
    if status == "excluded_other_task_do_not_touch":
        return {
            "type_name": product["type_name"],
            "slug": slug,
            "status": status,
            "excluded": True,
            "do_not_touch": True,
            "folder_exists": folder.exists(),
            "mechanical_pass": True,
        }

    pending = status == "review_ready_pending_approval"
    completed = status == "completed"
    required_names = [
        "plan.svg",
        "front.svg",
        "side.svg",
        "manifest.json",
        "candidate-representations.json",
        "bonsai-review-manifest.json",
        "project-context-manifest.json",
        "index.html",
    ]
    if pending:
        required_names.extend(["profile.json", "review-contact-sheet.png"])
    missing = [name for name in required_names if not (folder / name).is_file()]

    manifest_path = folder / "manifest.json"
    approval_path = ROOT / f"pipeline/decisions/{slug}-drawing-approval.json"
    bonsai_path = folder / "bonsai-review-manifest.json"
    context_path = folder / "project-context-manifest.json"
    manifest = load_json(manifest_path) if manifest_path.is_file() else {}
    candidate_path = folder / "candidate-representations.json"
    candidate = load_json(candidate_path) if candidate_path.is_file() else {}
    approval = load_json(approval_path) if approval_path.is_file() else {}
    bonsai = load_json(bonsai_path) if bonsai_path.is_file() else {}
    context = load_json(context_path) if context_path.is_file() else {}
    source_kind, official_cad_used = drawing_source(manifest)
    source_gate = source_provenance_gate(
        completed=completed,
        folder=folder,
        manifest=manifest,
        candidate=candidate,
        context=context,
        source_kind=source_kind,
        official_cad_used=official_cad_used,
    )

    manifest_hash = sha256(manifest_path) if manifest_path.is_file() else None
    approval_hash_matches = bool(
        manifest_hash
        and approval.get("candidate_manifest_sha256") == manifest_hash
    )
    approval_status = approval.get("status", "missing")
    approval_gate_pass = (
        approval_hash_matches
        and (
            (pending and approval_status_is_pending(approval_status))
            or (completed and approval_status == "approved")
        )
    )

    saved_camera_count = bonsai.get("bonsai_session", {}).get("saved_camera_count")
    actual_bonsai = (
        bonsai.get("mode") == "actual_bonsai_ifc_body_camera_render"
        or (
            bonsai.get("camera_renders_are_render_results") is True
            and bonsai.get("camera_renders", {}).get("render_source")
            == "actual_ifc_body_representation"
        )
    )
    bonsai_pass = bool(
        bonsai.get("pass") is True
        and actual_bonsai
        and saved_camera_count == 4
        and bonsai.get("geometry_product_count") == 1
        and bonsai.get("whole_model_render") is False
    )

    context_outputs = context_output_paths(context)
    context_outputs_exist = bool(context_outputs) and all(
        path.is_file() for path in context_outputs
    )
    if pending:
        context_pass = bool(
            context.get("pass") is True
            and context.get("project_context_retained") is True
            and context.get("walls_and_surrounding_project_elements_retained")
            is True
            and context_outputs_exist
        )
    else:
        # Falper predates the shared context flags. Its two approved project
        # drawings are proven directly to retain the IfcWall groups.
        context_pass = bool(
            completed
            and context.get("status") == "approved"
            and len(context_outputs) >= 2
            and context_outputs_exist
            and all(
                "class=\"IfcWall" in path.read_text(encoding="utf-8")
                for path in context_outputs
            )
        )

    derived_ifcs = sorted(
        relative(path) for path in folder.glob("*derived*.ifc") if path.is_file()
    )
    derived_gate_pass = (pending and not derived_ifcs) or (
        completed and bool(derived_ifcs)
    )
    writer_path = ROOT / f"pipeline/scripts/{slug.replace('-', '_')}_drawing_ifc.py"
    test_path = ROOT / f"pipeline/tests/{slug}-highpoly-review.test.ts"
    commits = [
        commit
        for commit in run_git(
            "log", "--all", "--format=%H", "--", product["folder"]
        ).splitlines()
        if commit
    ]
    commit_gate_pass = (pending and not commits) or (completed and bool(commits))

    formal_hashes = {
        manifest.get("formal_ifc_sha256"),
        bonsai.get("formal_ifc_sha256"),
    }
    formal_hashes.discard(None)
    formal_gate_pass = bool(
        formal_hash == EXPECTED_FORMAL_HASH
        and formal_hashes
        and formal_hashes == {EXPECTED_FORMAL_HASH}
    )
    package_pass = not missing
    mechanical_pass = all(
        [
            package_pass,
            formal_gate_pass,
            approval_gate_pass,
            bonsai_pass,
            context_pass,
            source_gate["pass"],
            derived_gate_pass,
            commit_gate_pass,
            writer_path.is_file(),
            test_path.is_file(),
        ]
    )

    return {
        "type_name": product["type_name"],
        "slug": slug,
        "status": status,
        "folder": product["folder"],
        "representative_global_id": product.get("representative_global_id"),
        "source_kind": source_kind,
        "official_cad_used": official_cad_used,
        "drawing_source_gate": source_gate,
        "package": {
            "required_files": required_names,
            "missing_files": missing,
            "pass": package_pass,
            "index": relative(folder / "index.html"),
            "contact_sheet": relative(folder / "review-contact-sheet.png")
            if (folder / "review-contact-sheet.png").is_file()
            else None,
        },
        "formal_ifc": {
            "sha256": formal_hash,
            "manifest_hashes": sorted(formal_hashes),
            "bytes_unchanged": formal_gate_pass,
        },
        "bonsai": {
            "mode": bonsai.get("mode") or "legacy_actual_ifc_body_camera_render",
            "actual_ifc_body_camera_render": actual_bonsai,
            "saved_camera_count": saved_camera_count,
            "one_product_only": bonsai.get("geometry_product_count") == 1,
            "whole_model_render": bonsai.get("whole_model_render"),
            "pass": bonsai_pass,
        },
        "project_context": {
            "view_count": len(context.get("views", [])),
            "project_context_retained": context.get(
                "project_context_retained", completed
            ),
            "walls_and_surrounding_project_elements_retained": context.get(
                "walls_and_surrounding_project_elements_retained", completed
            ),
            "outputs": [relative(path) for path in context_outputs],
            "pass": context_pass,
            "acceptance_rule": (
                "project Plan/Elevation retains walls and surrounding furniture/equipment; "
                "the product drawing representation is overlaid in-place"
            ),
        },
        "approval": {
            "record": relative(approval_path),
            "status": approval_status,
            "candidate_manifest_sha256": approval.get(
                "candidate_manifest_sha256"
            ),
            "manifest_sha256": manifest_hash,
            "hash_matches": approval_hash_matches,
            "write_allowed": bool(
                approval.get("derived_ifc_write_allowed")
                or approval.get("formal_ifc_write_allowed")
            ),
            "pass": approval_gate_pass,
        },
        "derived_ifcs": derived_ifcs,
        "derived_ifc_gate_pass": derived_gate_pass,
        "writer": relative(writer_path),
        "writer_exists": writer_path.is_file(),
        "test": relative(test_path),
        "test_exists": test_path.is_file(),
        "commits": commits,
        "commit_gate_pass": commit_gate_pass,
        "mechanical_pass": mechanical_pass,
    }


def write_html(report: dict[str, Any], path: Path) -> None:
    rows: list[str] = []
    for product in report["products"]:
        if product.get("excluded"):
            rows.append(
                "<tr class='excluded'>"
                f"<td>{html.escape(product['type_name'])}</td>"
                f"<td>{html.escape(product['status'])}</td>"
                "<td colspan='7'>Other task: do not touch</td></tr>"
            )
            continue
        package = product["package"]
        context = product["project_context"]
        approval = product["approval"]
        links = [
            f"<a href='{html.escape(product['slug'])}/index.html'>folder</a>",
            f"<a href='{html.escape(product['slug'])}/project-context-manifest.json'>context</a>",
        ]
        if package["contact_sheet"]:
            links.append(
                f"<a href='{html.escape(product['slug'])}/review-contact-sheet.png'>contact sheet</a>"
            )
        rows.append(
            f"<tr class='{'pass' if product['mechanical_pass'] else 'fail'}'>"
            f"<td>{html.escape(product['type_name'])}</td>"
            f"<td>{html.escape(product['status'])}</td>"
            f"<td>{html.escape(product['source_kind'])}</td>"
            f"<td>{'yes' if product['official_cad_used'] else 'no'}</td>"
            f"<td>{context['view_count']} / {'pass' if context['pass'] else 'fail'}</td>"
            f"<td>{product['bonsai']['saved_camera_count']} / {'pass' if product['bonsai']['pass'] else 'fail'}</td>"
            f"<td>{html.escape(approval['status'])} / {'hash ok' if approval['hash_matches'] else 'hash fail'}</td>"
            f"<td>{len(product['derived_ifcs'])} / {len(product['commits'])}</td>"
            f"<td>{' · '.join(links)}</td></tr>"
        )

    summary = report["summary"]
    path.write_text(
        "<!doctype html><html lang='zh-CN'><meta charset='utf-8'>"
        "<title>High-poly drawing review queue</title>"
        "<style>body{font:15px system-ui;margin:32px;color:#17202a}"
        "table{border-collapse:collapse;width:100%}th,td{padding:9px;border:1px solid #ccd1d1;text-align:left}"
        "th{background:#edf2f7}.pass{background:#f4fbf5}.fail{background:#fff0f0}.excluded{background:#f5f5f5}"
        "code{font-size:12px}a{color:#1677c8}</style>"
        "<h1>High-poly drawing review queue</h1>"
        f"<p>已完成 {summary['completed_count']}；待人眼批准 {summary['pending_human_approval_count']}；"
        f"其他任务隔离 {summary['excluded_other_task_count']}。机械门："
        f"<strong>{'PASS' if summary['package_mechanical_pass'] else 'FAIL'}</strong>。</p>"
        "<p>项目上下文验收：Plan/Elevation 必须保留墙体及周边家具/设备，产品 representation 原位叠加；"
        "单体 SVG 或孤立相机渲染不能替代项目图纸。</p>"
        f"<p>正式 IFC SHA-256：<code>{html.escape(report['formal_ifc_sha256'])}</code></p>"
        "<table><thead><tr><th>Product</th><th>Status</th><th>Drawing source</th>"
        "<th>Official CAD</th><th>Project context</th><th>Bonsai cameras</th>"
        "<th>Approval</th><th>Derived IFC / commits</th><th>Review links</th></tr></thead>"
        f"<tbody>{''.join(rows)}</tbody></table></html>",
        encoding="utf-8",
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    output = args.output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)

    inventory = load_json(INVENTORY)
    formal_hash = sha256(FORMAL_IFC)
    products = [audit_product(product, formal_hash) for product in inventory["products"]]
    active = [product for product in products if not product.get("excluded")]
    pending = [
        product
        for product in active
        if product["status"] == "review_ready_pending_approval"
    ]
    completed = [product for product in active if product["status"] == "completed"]
    excluded = [product for product in products if product.get("excluded")]
    stashes = run_git("stash", "list").splitlines()
    stash_recovery = load_json(STASH_RECOVERY_AUDIT)
    stash_commits = {
        line.split("\t", 1)[0]
        for line in run_git("stash", "list", "--format=%H%x09%gd%x09%gs").splitlines()
    }
    stash_recovery_pass = bool(
        stash_recovery.get("pass") is True
        and stash_recovery.get("source_stash", {}).get("commit") in stash_commits
        and stash_recovery.get("stash_apply_performed") is False
        and stash_recovery.get("stash_pop_performed") is False
        and stash_recovery.get("stash_drop_performed") is False
        and stash_recovery.get("stash_mutation_performed") is False
        and stash_recovery.get("protected_other_task_stash", {}).get("content_inspected") is False
        and stash_recovery.get("protected_other_task_stash", {}).get("mutated") is False
        and all(
            record.get("current_exists") is True
            for record in stash_recovery.get("untracked_snapshot", [])
        )
    )
    report = {
        "schema_version": 1,
        "formal_ifc": relative(FORMAL_IFC),
        "formal_ifc_sha256": formal_hash,
        "formal_ifc_bytes_unchanged": formal_hash == EXPECTED_FORMAL_HASH,
        "stashes_observed_read_only": stashes,
        "stash_mutation_performed": False,
        "stash_recovery": {
            "audit": relative(STASH_RECOVERY_AUDIT),
            "source_stash_commit": stash_recovery.get("source_stash", {}).get("commit"),
            "source_untracked_tree": stash_recovery.get("source_untracked_tree"),
            "recovered_file_count": len(stash_recovery.get("untracked_snapshot", [])),
            "candidate_product_slugs": stash_recovery.get("candidate_product_slugs", []),
            "protected_other_task_stash": stash_recovery.get("protected_other_task_stash"),
            "pass": stash_recovery_pass,
        },
        "products": products,
        "summary": {
            "product_count": len(products),
            "completed_count": len(completed),
            "fully_packaged_pending_count": sum(
                product["mechanical_pass"] for product in pending
            ),
            "pending_human_approval_count": len(pending),
            "excluded_other_task_count": len(excluded),
            "derived_ifc_product_count": sum(
                bool(product["derived_ifcs"]) for product in active
            ),
            "independently_committed_product_count": sum(
                bool(product["commits"]) for product in active
            ),
            "package_mechanical_pass": bool(active)
            and all(product["mechanical_pass"] for product in active)
            and stash_recovery_pass,
            "stash_recovery_pass": stash_recovery_pass,
            "goal_complete": False,
            "remaining_gate": (
                "Each pending product requires explicit human approval before its "
                "derived IFC is written and that product folder is committed."
            ),
        },
    }
    output.write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    write_html(report, output.with_suffix(".html"))
    print(json.dumps(report["summary"], ensure_ascii=False))
    if not report["summary"]["package_mechanical_pass"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
