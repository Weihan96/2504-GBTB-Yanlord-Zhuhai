import { expect, test } from "bun:test";
import { createHash } from "node:crypto";
import { existsSync, readFileSync } from "node:fs";
import { join, resolve } from "node:path";

const root = resolve(import.meta.dir, "../..");
const product = join(root, "output/review/highpoly-types/bed01");
const formal = join(root, "2504 GBTB Yanlord Zhuhai.ifc");
const formalHash = "7a50b87e8f48a7c2bfbaa9f1dabdfca7155fa684b2325c66aa1ae15c8ab0a25c";
const sourceHash = "ab697db9c27448c62a9a77537d7cc4b6286577c335328113d703184be28d9c4a";
const pathCounts = { plan: 25, front: 22, side: 36 } as const;

function sha256(path: string) {
  return createHash("sha256").update(readFileSync(path)).digest("hex");
}

test("BED01 approval opens only the product-level official-DWG derived write gate", () => {
  const approvalPath = join(root, "pipeline/decisions/bed01-drawing-approval.json");
  const approval = JSON.parse(readFileSync(approvalPath, "utf8"));
  expect(approval).toMatchObject({
    profile_key: "bed01",
    status: "approved",
    reviewer: "project_owner",
    approved_views: ["plan", "front", "side"],
    derived_ifc_write_allowed: true,
    formal_authoritative_ifc_write_allowed: false,
  });
  expect(approval.approved_candidate).toMatchObject({
    source_kind: "native_dwg_review_reference",
    source_dwg_sha256: sourceHash,
    alignment: "Plan and Front translation only; Side centre reflection plus translation; scale 1.0",
    superseded_acis_candidate_excluded: true,
  });
  expect(sha256(formal)).toBe(formalHash);
});

test("BED01 persisted IFC retains three real Bonsai Drawings and native-DWG Annotations", () => {
  const manifestPath = join(product, "master-bedroom-bonsai-drawing-manifest.json");
  const manifest = JSON.parse(readFileSync(manifestPath, "utf8"));
  const evidence = JSON.parse(
    readFileSync(join(root, manifest.create_drawing_evidence.path), "utf8"),
  );
  expect(manifest).toMatchObject({
    pass: true,
    formal_ifc_sha256: formalHash,
    formal_ifc_bytes_unchanged: true,
    source_kind: "native_dwg_review_reference",
    source_dwg_sha256: sourceHash,
    source_scale: 1,
    superseded_acis_candidate_excluded: true,
    room_global_id: "3gHz6U6BfFXgV6PnRzfOf$",
    room_name: "主卧",
    target_global_id: "3IQBEqO5vDI8Z9k1Ltge_N",
    actual_target_body_suppressed_from_drawing_include: true,
  });
  expect(manifest.tests).toMatchObject({
    all_views_create_drawing_finished: true,
    all_views_opencascade: true,
    all_persisted_path_counts_match: true,
    all_previews_blue_and_grey: true,
    all_previews_not_clipped: true,
    actual_bed_body_projection_group_count: 0,
    formal_ifc_hash_preserved: true,
  });
  expect(sha256(join(root, manifest.derived_ifc.path))).toBe(manifest.derived_ifc.sha256);
  expect(evidence.persistence).toMatchObject({
    reload_result: ["FINISHED"],
    post_reload_drawing_count: 3,
    post_reload_annotation_count: 3,
  });
  expect(evidence.formal_ifc_bytes_unchanged).toBeTrue();
  expect(evidence.views).toHaveLength(3);

  for (const view of manifest.views as any[]) {
    const key = view.view as keyof typeof pathCounts;
    expect(view.create_drawing).toMatchObject({
      operator: "bpy.ops.bim.create_drawing",
      result: ["FINISHED"],
      linework_mode: "OPENCASCADE",
    });
    expect(view.review_path_count).toBe(pathCounts[key]);
    expect(view.persisted_review_path_count).toBe(pathCounts[key]);
    expect(view.persisted_coordinate_residual_mm).toBeLessThanOrEqual(0.000001);
    const persistedView = evidence.views.find((item: any) => item.view === view.view);
    expect(persistedView.linework_annotation.drawing_group_association_persisted).toBeTrue();
    expect(persistedView.linework_annotation.drawing_group_assignment_global_id).toHaveLength(22);
    expect(view.alignment.uniform_scale).toBe(1);
    expect(view.alignment.anisotropic_scale_used).toBeFalse();
    expect(view.alignment.source_geometry_deformed).toBeFalse();
    expect(view.pdf_preview_gate).toMatchObject({
      blue_line_present: true,
      grey_context_present: true,
      no_page_edge_clipping: true,
      pass: true,
    });
    for (const artifact of [view.svg, view.page_pdf, view.pdf_rendered_preview_png]) {
      expect(existsSync(join(root, artifact.path))).toBeTrue();
      expect(sha256(join(root, artifact.path))).toBe(artifact.sha256);
    }
    const svg = readFileSync(join(root, view.svg.path), "utf8");
    expect(svg).toContain('data-create-drawing-result="FINISHED"');
    expect(svg).toContain('data-context-colour="#a3abb3"');
    expect(svg).toContain('data-official-linework-colour="#1677c8"');
    expect(svg).toContain("official-native-dwg");
  }
  expect(sha256(formal)).toBe(formalHash);
});
