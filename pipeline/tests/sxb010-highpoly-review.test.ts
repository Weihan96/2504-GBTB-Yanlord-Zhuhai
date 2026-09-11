import { expect, test } from "bun:test";
import { createHash } from "node:crypto";
import { existsSync, mkdtempSync, readFileSync, rmSync, writeFileSync } from "node:fs";
import { tmpdir } from "node:os";
import { join, resolve } from "node:path";

const root = resolve(import.meta.dir, "../..");
const product = join(root, "output/review/highpoly-types/sxb010");
const formal = join(root, "2504 GBTB Yanlord Zhuhai.ifc");
const formalHash = "7a50b87e8f48a7c2bfbaa9f1dabdfca7155fa684b2325c66aa1ae15c8ab0a25c";
const sourceKind = "geometry_derived_simplified_proxy";
const sourceLabelZh = "基于原始高模几何生成的简化图纸表达";

function sha256(path: string) {
  return createHash("sha256").update(readFileSync(path)).digest("hex");
}

test("owner direction and official Hunter Douglas sources do not masquerade as project CAD", () => {
  const access = JSON.parse(readFileSync(join(product, "official-source/source-access-record.json"), "utf8"));
  expect(access.identity_status).toBe("owner_confirmed_family_direction_exact_sku_finish_and_operation_pending");
  expect(access.owner_direction_evidence.source).toBe(
    "stash@{1}:drawings/evidence/OWNER-INT1-全屋窗饰方案-20260818.md",
  );
  expect(access.owner_direction_evidence.git_blob_sha1).toBe("1915cc02e5d992571a9a451a52ed4d8ad8216769");
  expect(access.owner_direction_evidence.sha256).toBe("319afadd3ec1308eadee50d7422f977621b05d540411620de8c26a25df07f49a");
  expect(access.owner_direction_evidence.stash_operation).toBe("read_only_git_show_only_no_pop_no_apply");
  expect(access.official_product_cad.acquired).toBeFalse();
  expect(access.official_product_cad.local_cad_files).toEqual([]);
  expect(access.official_product_cad.official_eu_filtered_download_formats).toEqual(["PDF", "PDF"]);
  const brochure = access.official_identity_sources.find(
    (item: any) => item.kind === "manufacturer_technical_brochure",
  );
  expect(brochure.printed_pages_visually_checked).toEqual([2, 3]);
  expect(brochure.pdf_file_pages_visually_checked).toEqual([4, 5]);
  expect(access.drawing_geometry_source).toMatchObject({
    source_kind: sourceKind,
    source_label_zh: sourceLabelZh,
    official_cad_used: false,
    third_party_cad_used: false,
  });
  for (const evidence of access.official_identity_sources) {
    expect(sha256(join(root, evidence.local_path))).toBe(evidence.sha256);
    if (evidence.filtered_response_local_path) {
      expect(sha256(join(root, evidence.filtered_response_local_path))).toBe(evidence.filtered_response_sha256);
      const filtered = readFileSync(join(root, evidence.filtered_response_local_path), "utf8");
      expect((filtered.match(/\.pdf/g) || []).length).toBe(2);
      expect(filtered.toLowerCase()).not.toContain(".dwg");
      expect(filtered.toLowerCase()).not.toContain(".dxf");
    }
  }
  expect(access.dimension_cross_check.project_width_within_official_limit).toBeTrue();
  expect(access.dimension_cross_check.project_height_within_official_limit).toBeTrue();
  expect(access.dimension_cross_check.status).toBe(
    "family_and_size_envelope_consistent_exact_product_and_shop_geometry_unconfirmed",
  );
});

test("single untyped representative uses semantic Plan Front Side axes and no blue CAD line", () => {
  const manifest = JSON.parse(readFileSync(join(product, "manifest.json"), "utf8"));
  const candidate = JSON.parse(readFileSync(join(product, "candidate-representations.json"), "utf8"));
  expect(manifest.representative_global_id).toBe("1O9JRXCI56VRUbpuLJy86Z");
  expect(manifest.registered_instance_global_ids).toEqual(["1O9JRXCI56VRUbpuLJy86Z"]);
  expect(manifest.ifc_type_name).toBe("sxb010");
  expect(manifest.ifc_type_description).toBeNull();
  expect(manifest.geometry_product_count).toBe(1);
  expect(manifest.whole_model_render).toBeFalse();
  expect(manifest.formal_ifc_bytes_unchanged).toBeTrue();
  expect(manifest.bounds_mm.size).toEqual([1811.999451, 1596.075623, 63.393311]);
  expect(candidate.source_kind).toBe(sourceKind);
  expect(candidate.source_label_zh).toBe(sourceLabelZh);
  expect(candidate.official_cad_used).toBeFalse();
  expect(candidate.third_party_cad_used).toBeFalse();
  expect(Object.fromEntries(manifest.views.map((view: any) => [view.view, view.projection_axes]))).toEqual({
    plan: [0, 2],
    front: [0, 1],
    side: [2, 1],
  });
  expect(Object.fromEntries(manifest.views.map((view: any) => [view.view, view.silhouette_path_count]))).toEqual({
    plan: 5,
    front: 51,
    side: 55,
  });
  for (const view of manifest.views) {
    expect(view.drawing_line_source_kind).toBe(sourceKind);
    expect(view.drawing_line_source_label_zh).toBe(sourceLabelZh);
    expect(view.official_cad_path_count).toBe(0);
    expect(view.blue_line_present).toBeFalse();
    const svg = readFileSync(join(root, view.svg), "utf8");
    expect(svg).toContain('class="simplified-proxy-silhouette geometry-derived"');
    expect(svg).toContain('data-source-kind="geometry_derived_simplified_proxy"');
    expect(svg).not.toContain('class="official-reference native-dwg"');
    expect(svg).not.toContain("#1677c8");
  }
});

test("22 mm near-line merge keeps all slats, controls and the exact projected envelope", () => {
  const audit = JSON.parse(readFileSync(join(product, "line-simplification-audit.json"), "utf8"));
  const manifest = JSON.parse(readFileSync(join(product, "manifest.json"), "utf8"));
  const candidate = JSON.parse(readFileSync(join(product, "candidate-representations.json"), "utf8"));
  expect(audit.source_kind).toBe(sourceKind);
  expect(audit.source_label_zh).toBe(sourceLabelZh);
  expect(audit.merge_threshold_mm).toBe(22);
  expect(audit.review_status).toBe("visual_review_pending");
  expect(audit.derived_ifc_write_allowed).toBeFalse();
  expect(audit.formal_ifc_bytes_unchanged).toBeTrue();
  expect(audit.pass).toBeTrue();
  expect(manifest.line_simplification_audit_sha256).toBe(sha256(join(product, "line-simplification-audit.json")));
  expect(candidate.near_line_merge_threshold_mm).toBe(22);

  const expected = {
    plan: { beforePaths: 8, afterPaths: 5, beforeSegments: 63, afterSegments: 8 },
    front: { beforePaths: 241, afterPaths: 51, beforeSegments: 1444, afterSegments: 60 },
    side: { beforePaths: 84, afterPaths: 55, beforeSegments: 696, afterSegments: 58 },
  };
  for (const [view, record] of Object.entries(audit.views) as any) {
    expect(record.before_path_count).toBe(expected[view as keyof typeof expected].beforePaths);
    expect(record.after_path_count).toBe(expected[view as keyof typeof expected].afterPaths);
    expect(record.before_segment_count).toBe(expected[view as keyof typeof expected].beforeSegments);
    expect(record.after_segment_count).toBe(expected[view as keyof typeof expected].afterSegments);
    expect(record.merge_threshold_at_final_svg_scale_px).toBeGreaterThan(8.9);
    expect(record.merge_threshold_at_final_svg_scale_px).toBeLessThan(10.5);
    expect(Math.max(...record.outer_envelope_delta_mm)).toBe(0);
    expect(record.before_bounds_mm.minimum).toEqual(record.after_bounds_mm.minimum);
    expect(record.before_bounds_mm.maximum).toEqual(record.after_bounds_mm.maximum);
    expect(record.slat_count_before).toBe(45);
    expect(record.slat_centerline_count_after).toBe(45);
    expect(record.slat_center_pitch_after_mm.minimum).toBeGreaterThan(record.merge_threshold_mm);
    expect(record.headrail_count_preserved).toBe(1);
    expect(record.bottom_rail_count_preserved).toBe(1);
    expect(record.guide_and_cord_component_count_before).toBe(7);
    expect(record.guide_and_cord_axis_count_after).toBe(3);
    expect(record.control_component_count_preserved).toBe(9);
    expect(record.outer_envelope_preserved).toBeTrue();
    expect(record.slat_rhythm_preserved).toBeTrue();
  }
});

test("project plan and R07 elevations retain context with mechanically located semantic views", () => {
  const context = JSON.parse(readFileSync(join(product, "project-context-manifest.json"), "utf8"));
  expect(context.source_kind).toBe(sourceKind);
  expect(context.official_cad_used).toBeFalse();
  expect(context.third_party_cad_used).toBeFalse();
  expect(context.project_context_retained).toBeTrue();
  expect(context.walls_and_surrounding_project_elements_retained).toBeTrue();
  expect(context.overlay_top_layer_with_white_mask).toBeTrue();
  expect(context.blue_line_present).toBeFalse();
  expect(context.legacy_plan_source_group_count).toBe(2);
  expect(context.selected_legacy_plan_source_group).toBe("densest_Body_projection_not_Curve3D");
  expect(context.views.map((view: any) => view.view)).toEqual(["plan", "front", "side"]);
  expect(context.semantic_view_mapping).toEqual({
    plan: { candidate_axes: [0, 2], source: "drawings/FFL PLAN.svg" },
    front: {
      candidate_axes: [0, 1],
      source: "drawings/elevations/native/EL-04-10-R07-PY.svg",
    },
    side: {
      candidate_axes: [2, 1],
      source: "drawings/elevations/native/EL-04-11-R07-PX.svg",
    },
  });
  expect(context.project_world_bbox_m).toEqual([
    [-2.578700412869453, -4.9666923713684, 0.833773986816406],
    [-0.766700962185859, -4.90329906082153, 2.4298496093749997],
  ]);
  expect(context.context_target_derivation.no_nonuniform_fit_or_blank_space_guessing).toBeTrue();
  expect(context.rejected_context_evidence.used_as_alignment_target).toBeFalse();
  expect(context.review_annotation_suppression.geometry_removed).toBeFalse();
  const expectedPathCounts = { plan: 5, front: 51, side: 55 };
  for (const view of context.views) {
    expect(view.overlay.path_count).toBe(expectedPathCounts[view.view as keyof typeof expectedPathCounts]);
    expect(Math.max(...view.overlay.fit.bbox_absolute_delta_svg_units)).toBeLessThanOrEqual(
      view.overlay.fit.bbox_tolerance_svg_units,
    );
    expect(view.overlay.fit.uniform_scale_preserved).toBeTrue();
    expect(sha256(join(root, view.review_preview))).toBe(view.review_preview_sha256);
    const svg = readFileSync(join(root, view.output), "utf8");
    expect(svg).toContain("IfcFurniture");
    expect(svg).toMatch(/IfcWall|IfcSlab|IfcCovering/);
    expect(svg).toContain('class="geometry-derived-simplified-proxy project-context-overlay"');
    expect(svg).toContain('class="geometry-derived-proxy-mask"');
    expect(svg).toContain('class="geometry-derived-proxy"');
    expect(svg.indexOf('class="geometry-derived-proxy-mask"')).toBeLessThan(
      svg.indexOf('class="geometry-derived-proxy"'),
    );
    expect(svg).not.toContain('class="official-reference native-dwg"');
  }
  const planSvg = readFileSync(join(root, context.views[0].output), "utf8");
  expect(planSvg).toContain("IfcWall");
  expect(planSvg).not.toContain('class="official-elevation-anchor"');
  expect(context.pass).toBeTrue();
});

test("Bonsai evidence maps local-axis cameras to semantic drawing views without relabelling geometry", () => {
  const evidence = JSON.parse(readFileSync(join(product, "bonsai-review-manifest.json"), "utf8"));
  expect(evidence.mode).toBe("actual_bonsai_ifc_body_camera_render");
  expect(evidence.geometry_product_count).toBe(1);
  expect(evidence.whole_model_render).toBeFalse();
  expect(evidence.formal_ifc_bytes_unchanged).toBeTrue();
  expect(evidence.bonsai_session.saved_active_representation).toBe("Body");
  expect(evidence.bonsai_session.ifc_context_identifier).toBe("Body");
  expect(evidence.bonsai_session.saved_camera_count).toBe(4);
  expect(evidence.semantic_view_mapping).toEqual({
    plan: {
      render_view_record: "front",
      source_camera_path: "output/review/highpoly-types/sxb010/bonsai-camera-front-elevation.png",
      path: "output/review/highpoly-types/sxb010/bonsai-camera-semantic-plan.png",
      project_axes: [0, 2],
    },
    front: {
      render_view_record: "plan",
      source_camera_path: "output/review/highpoly-types/sxb010/bonsai-camera-plan.png",
      path: "output/review/highpoly-types/sxb010/bonsai-camera-semantic-front-elevation.png",
      project_axes: [0, 1],
    },
    side: {
      render_view_record: "side",
      source_camera_path: "output/review/highpoly-types/sxb010/bonsai-camera-side-elevation.png",
      path: "output/review/highpoly-types/sxb010/bonsai-camera-semantic-side-elevation.png",
      project_axes: [2, 1],
    },
  });
  for (const mapping of Object.values(evidence.semantic_view_mapping) as any[]) {
    expect(sha256(join(root, mapping.path))).toBe(sha256(join(root, mapping.source_camera_path)));
  }
  for (const render of evidence.renders) {
    expect(render.camera_type).toBe("ORTHO");
    expect(render.render_operation).toBe("bpy.ops.render.render(write_still=True)");
    expect(render.camera_ortho_scale_m + 1e-6).toBeGreaterThanOrEqual(
      render.projected_width_m * render.framing_margin_factor,
    );
    expect(render.camera_ortho_scale_m + 1e-6).toBeGreaterThanOrEqual(
      render.projected_height_m * render.image_aspect * render.framing_margin_factor,
    );
    expect(sha256(join(root, render.path))).toBe(render.sha256);
  }
  expect(sha256(join(root, evidence.bonsai_session.path))).toBe(evidence.bonsai_session.sha256);
  expect(sha256(join(root, evidence.isolated_ifc))).toBe(evidence.isolated_ifc_sha256);
});

test("pending sxb010 review leaves the formal IFC byte-identical and creates no derived IFC", () => {
  const candidate = JSON.parse(readFileSync(join(product, "candidate-representations.json"), "utf8"));
  expect(candidate.formal_ifc_write_allowed).toBeFalse();
  expect(candidate.review_status).toBe("visual_review_pending");
  expect(existsSync(join(product, "Hunter-Douglas-sxb010-derived-drawing.ifc"))).toBeFalse();
  expect(sha256(formal)).toBe(formalHash);
});

test("approved sxb010 record still requires explicit apply and forbids formal IFC writes", () => {
  const temporary = mkdtempSync(join(tmpdir(), "sxb010-gate-"));
  const output = join(temporary, "forbidden.ifc");
  const script = join(root, "pipeline/scripts/sxb010_drawing_ifc.py");
  const approval = JSON.parse(readFileSync(join(root, "pipeline/decisions/sxb010-drawing-approval.json"), "utf8"));
  const manifest = join(product, "manifest.json");
  expect(approval.status).toBe("approved");
  expect(approval.approved_views).toEqual(["plan", "front", "side"]);
  expect(approval.derived_ifc_write_allowed).toBeTrue();
  expect(approval.formal_authoritative_ifc_write_allowed).toBeFalse();
  expect(approval.candidate_manifest_sha256).toBe(sha256(manifest));

  const withoutApply = Bun.spawnSync(["python3", script, "--input", formal, "--output", output], {
    cwd: root,
    stderr: "pipe",
  });
  expect(withoutApply.exitCode).not.toBe(0);
  expect(withoutApply.stderr.toString()).toContain("IFC write requires the explicit --apply flag");
  expect(existsSync(output)).toBeFalse();

  expect(sha256(formal)).toBe(formalHash);
  rmSync(temporary, { recursive: true, force: true });
});

test("synthetic approval writes only a temporary derived IFC with traceable fallback sources", () => {
  const temporary = mkdtempSync(join(tmpdir(), "sxb010-approved-"));
  const output = join(temporary, "derived.ifc");
  const approvalPath = join(temporary, "approval.json");
  const script = join(root, "pipeline/scripts/sxb010_drawing_ifc.py");
  const manifest = join(product, "manifest.json");
  const pending = JSON.parse(
    readFileSync(join(root, "pipeline/decisions/sxb010-drawing-approval.json"), "utf8"),
  );
  writeFileSync(
    approvalPath,
    JSON.stringify(
      {
        ...pending,
        candidate_manifest_sha256: sha256(manifest),
        status: "approved",
        reviewer: "Automated synthetic gate test",
        review_date: "2026-08-23",
        approved_views: ["plan", "front", "side"],
        derived_ifc_write_allowed: true,
        formal_authoritative_ifc_write_allowed: false,
        approval_evidence: "Synthetic test only; output is deleted after verification.",
      },
      null,
      2,
    ) + "\n",
  );
  const result = Bun.spawnSync(
    [
      "python3",
      script,
      "--input",
      formal,
      "--approval",
      approvalPath,
      "--output",
      output,
      "--apply",
    ],
    { cwd: root, stdout: "pipe", stderr: "pipe" },
  );
  expect(result.exitCode, result.stderr.toString()).toBe(0);
  const report = JSON.parse(result.stdout.toString());
  expect(report.representations).toEqual({
    plan: "Sxb010Plan",
    front: "Sxb010Front",
    side: "Sxb010Side",
  });
  expect(report.representation_path_counts).toEqual({ plan: 5, front: 51, side: 55 });
  expect(report.source_kind).toBe(sourceKind);
  expect(report.source_label_zh).toBe(sourceLabelZh);
  expect(report.official_cad_geometry_included).toBeFalse();
  expect(report.source_document_associations).toEqual([
    "HUNTER-DOUGLAS-SXB010-OFFICIAL-DOWNLOAD-FILTER",
    "HUNTER-DOUGLAS-SXB010-OFFICIAL-DOWNLOADS",
    "HUNTER-DOUGLAS-SXB010-OFFICIAL-PRODUCT-FLYER",
    "HUNTER-DOUGLAS-SXB010-OFFICIAL-PRODUCT-PAGE",
    "HUNTER-DOUGLAS-SXB010-OFFICIAL-TECHNICAL-BROCHURE",
    "HUNTER-DOUGLAS-SXB010-OWNER-DIRECTION-EVIDENCE",
    "HUNTER-DOUGLAS-SXB010-SOURCE-ACCESS-RECORD",
  ]);
  expect(existsSync(output)).toBeTrue();
  expect(sha256(formal)).toBe(formalHash);
  rmSync(temporary, { recursive: true, force: true });
}, 20000);
