import { expect, test } from "bun:test";
import { createHash } from "node:crypto";
import { existsSync, mkdtempSync, readFileSync, rmSync, writeFileSync } from "node:fs";
import { tmpdir } from "node:os";
import { dirname, join, resolve } from "node:path";

const root = resolve(import.meta.dir, "../..");
const product = join(root, "output/review/highpoly-types/gessi316-54294");
const sourceDir = join(product, "official-source");
const formal = join(root, "2504 GBTB Yanlord Zhuhai.ifc");
const formalHash = "7a50b87e8f48a7c2bfbaa9f1dabdfca7155fa684b2325c66aa1ae15c8ab0a25c";
const sourceKind = "native_dwg";
const sourceLabelZh = "Gessi 官方精确型号 54294 原生 DWG 图纸表达";
const reviewSourceKind = "native_dwg_review_simplification";
const reviewSourceLabelZh = "基于官方 Gessi 54294 原生 DWG 轮廓的去纹审核简化表达";
const dwgHash = "dfe8bcf7e100c14187c387b75a35e3fe7f718c61595fadf66adfe92ca7648ff4";
const zipHash = "fad0d98e94a83bb488b5f0703472862d628759f21700c1a87c3af343f4c1bdb9";
const pdfHash = "82058bb27752f27e6fc92029783d67c15fd443befcd393f7d572fd0147a581e3";
const pathCounts = { plan: 1481, front: 1099, side: 745 };
const reviewPathCounts = { plan: 69, front: 1099, side: 35 };
const scope = "official Gessi exact 54294 family reference and 45089_54294 article combination; not a project shop drawing";

function sha256(path: string) {
  return createHash("sha256").update(readFileSync(path)).digest("hex");
}

test("exact public Gessi 54294 DWG and PDF sources are archived and cross-verified", () => {
  const access = JSON.parse(readFileSync(join(sourceDir, "source-access-record.json"), "utf8"));
  const revalidation = JSON.parse(readFileSync(join(sourceDir, "official-source-revalidation.json"), "utf8"));
  const pdfVerification = JSON.parse(readFileSync(join(sourceDir, "official-pdf-verification.json"), "utf8"));
  expect(access).toMatchObject({
    manufacturer: "Gessi",
    project_ifc_type_name: "Gessi316 54294",
    resolved_article_number: "45089_54294",
    external_visible_product_code: "54294",
    companion_built_in_product_code: "45089",
    scope,
    pass: true,
  });
  expect(access.official_product_cad).toMatchObject({
    authentication_required: false,
    acquired: true,
    exact_article_match: true,
    exact_project_configuration_match: false,
  });
  expect(access.dimension_cross_check).toMatchObject({
    official_dwg_wall_surface_x_mm: -18.818959,
    official_dwg_visible_wall_to_outlet_reach_mm: 209.812362,
    project_ifc_visible_wall_to_outlet_reach_mm: 202.846329,
    official_vs_project_visible_reach_delta_mm: 6.966033,
    technical_pdf_adjustable_depth_range_mm: [190, 210],
    project_body_inside_adjustable_depth_range: true,
    previous_side_overlay_transform_error_mm: 18.797845,
  });
  expect(access.drawing_geometry_source).toMatchObject({
    source_kind: sourceKind,
    source_label_zh: sourceLabelZh,
    official_cad_used: true,
    third_party_cad_used: false,
    adjacent_product_cad_used: false,
    source_dwg_sha256: dwgHash,
  });
  expect(access.identity_and_geometry_policy).toMatchObject({
    "54294_native_dwg_used_for_three_views": true,
    "45089_companion_dwg_used_as_54294_geometry": false,
    adjacent_gessi_product_cad_used: false,
    third_party_cad_used: false,
  });
  expect(sha256(join(sourceDir, "GPF5429400000G000_3.dwg"))).toBe(dwgHash);
  expect(sha256(join(sourceDir, "GPF5429400000G000_arc.zip"))).toBe(zipHash);
  expect(sha256(join(sourceDir, "GPF5429400000G000_1.pdf"))).toBe(pdfHash);
  expect(revalidation).toMatchObject({
    all_downloaded_bytes_match_local_archive: true,
    all_zip_members_match_extracted_files: true,
    pass: true,
  });
  expect(revalidation.zip_members["54294"].member_sha256).toBe(dwgHash);
  expect(revalidation.drawing_geometry_policy.companion_45089_dwg_used_as_54294_geometry).toBeFalse();
  expect(pdfVerification.pass).toBeTrue();
  expect(pdfVerification.mechanical_cross_check.pass).toBeTrue();
  expect(pdfVerification.visual_review.status).toBe("codex_visual_qa_pass");
});

test("native DWG linework extraction is repeatable and dimensionally exact", () => {
  const temporary = mkdtempSync(join(tmpdir(), "gessi316-54294-linework-"));
  const output = join(temporary, "linework.json");
  const run = Bun.spawnSync([
    "python3", join(root, "pipeline/scripts/gessi316_54294_linework.py"), "--output", output,
  ], { cwd: root, stdout: "pipe", stderr: "pipe" });
  if (run.exitCode !== 0) throw new Error(run.stderr.toString());
  const linework = JSON.parse(readFileSync(output, "utf8"));
  expect(linework.source_kind).toBe(sourceKind);
  expect(linework.official_sources.native_dwg.sha256).toBe(dwgHash);
  expect(linework.official_sources.native_dwg_zip.sha256).toBe(zipHash);
  expect(linework.official_sources.technical_vector_pdf.sha256).toBe(pdfHash);
  expect(linework.identity_gates).toMatchObject({
    "54294_native_dwg_used": true,
    "45089_companion_dwg_used_as_54294_geometry": false,
    adjacent_gessi_product_cad_used: false,
    third_party_cad_used: false,
  });
  expect(linework.nominal_dimension_cross_check).toMatchObject({
    native_dwg_plan_envelope_mm: [361.7, 209.812412],
    absolute_delta_mm: [0.3, 0.187588],
    tolerance_mm: 0.5,
    pass: true,
  });
  for (const [view, count] of Object.entries(pathCounts)) {
    expect(linework.views[view].source_dwg_sha256).toBe(dwgHash);
    expect(linework.views[view].paths_mm).toHaveLength(count);
  }
  rmSync(temporary, { recursive: true, force: true });
}, 60_000);

test("review generator retains actual Body evidence and rebuilds the proxy to exact Side anchors", () => {
  const temporary = mkdtempSync(join(tmpdir(), "gessi316-54294-review-"));
  const run = Bun.spawnSync([
    "python3", join(root, "pipeline/scripts/gessi316_54294_review.py"), "--output", temporary,
  ], { cwd: root, stdout: "pipe", stderr: "pipe" });
  if (run.exitCode !== 0) throw new Error(run.stderr.toString());
  const generated = JSON.parse(readFileSync(join(temporary, "manifest.json"), "utf8"));
  expect(generated).toMatchObject({
    formal_ifc_bytes_unchanged: true,
    original_mesh_vertex_count: 21299,
    original_mesh_face_count: 41570,
    mesh_vertex_count: 1822,
    mesh_face_count: 3616,
  });
  expect(generated.review_only_geometry_simplification).toMatchObject({
    mode: "review_only_smooth_handles_and_official_max_reach_spout",
    removed_component_count: 3,
    removed_vertex_count: 20387,
    removed_face_count: 39762,
    handle_replacement_component_count: 2,
    handle_replacement_vertex_count: 260,
    handle_replacement_face_count: 512,
    spout_replacement_component_count: 1,
    spout_replacement_vertex_count: 650,
    spout_replacement_face_count: 1296,
    optimised_vertex_count: 1822,
    optimised_face_count: 3616,
    main_handle_cylinder_count_preserved: 2,
    control_levers_preserved: true,
    spout_component_count_preserved: 1,
    three_hole_trim_disc_count_preserved: 3,
    plan_front_identity_and_three_hole_layout_preserved: true,
    pass: true,
  });
  const mechanical = generated.review_only_geometry_simplification.side_mechanical_gate;
  expect(mechanical).toMatchObject({
    tolerance_mm: 0.2,
    actual_body_residual_mm: {
      wall_surface_y_mm: 0,
      straight_spout_axis_z_mm: 0,
      lower_lip_z_mm: 3.711117,
      outlet_endpoint_y_mm: 6.966033,
      visible_reach_mm: -6.966033,
      visible_envelope_mm: [-6.966033, -4.211117],
    },
    review_proxy_residual_mm: {
      wall_surface_y_mm: 0,
      straight_spout_axis_z_mm: 0,
      lower_lip_z_mm: 0,
      outlet_endpoint_y_mm: 0,
      visible_reach_mm: 0,
      visible_envelope_mm: [0, 0],
    },
    previous_axis_only_overlay: {
      wall_surface_residual_mm: -18.797845,
      outlet_endpoint_residual_mm: -25.763878,
      lower_lip_residual_mm: -3.711117,
      straight_spout_axis_residual_mm: 0,
    },
    correct_wall_and_axis_translation_mm: [28.613938, -21.051758],
    official_line_scaled: false,
    official_line_mirrored: false,
    wrong_dwg_view_selected: false,
    latest_exact_54294_dwg_used: true,
    project_body_within_official_adjustable_depth_range_190_210_mm: true,
    feature_gate_pass: true,
  });
  const side = generated.views.find((view: any) => view.view === "side");
  expect(side.translation_alignment).toMatchObject({
    scale_applied: false,
    mirror_applied: false,
    prior_wall_surface_residual_mm: -18.797845,
    prior_outlet_endpoint_residual_mm: -25.763878,
    corrected_anchor_residual_mm: [0, 0],
    pass: true,
  });
  expect(side.ifc_dwg_compatibility).toMatchObject({
    review_proxy_visible_side_envelope_mm: [209.812362, 71.551758],
    official_native_dwg_visible_side_envelope_mm: [209.812362, 71.551758],
    absolute_delta_mm: [0, 0],
    tolerance_mm: 0.2,
    pass: true,
  });
  expect(existsSync(join(temporary, "plan.svg"))).toBeTrue();
  expect(existsSync(join(temporary, "front.svg"))).toBeTrue();
  expect(existsSync(join(temporary, "side.svg"))).toBeTrue();
  rmSync(temporary, { recursive: true, force: true });
}, 60_000);

test("single-instance three-view SVGs place exact native DWG blue linework above IFC comparisons", () => {
  const manifest = JSON.parse(readFileSync(join(product, "manifest.json"), "utf8"));
  const candidate = JSON.parse(readFileSync(join(product, "candidate-representations.json"), "utf8"));
  expect(manifest.representative_global_id).toBe("2iKOL78$H0N9Yd9$ky3pW4");
  expect(manifest.registered_instance_global_ids).toEqual(["2iKOL78$H0N9Yd9$ky3pW4"]);
  expect(manifest.geometry_product_count).toBe(1);
  expect(manifest.whole_model_render).toBeFalse();
  expect(manifest.formal_ifc_bytes_unchanged).toBeTrue();
  expect(manifest).toMatchObject({
    original_mesh_vertex_count: 21299,
    original_mesh_face_count: 41570,
    mesh_vertex_count: 1822,
    mesh_face_count: 3616,
  });
  expect(manifest.review_only_geometry_simplification).toMatchObject({
    mode: "review_only_smooth_handles_and_official_max_reach_spout",
    removed_component_count: 3,
    removed_vertex_count: 20387,
    removed_face_count: 39762,
    handle_replacement_component_count: 2,
    handle_replacement_vertex_count: 260,
    handle_replacement_face_count: 512,
    spout_replacement_component_count: 1,
    spout_replacement_vertex_count: 650,
    spout_replacement_face_count: 1296,
    optimised_vertex_count: 1822,
    optimised_face_count: 3616,
    vertex_reduction_percent: 91.446,
    face_reduction_percent: 91.301,
    smooth_segments_per_sleeve: 64,
    main_handle_cylinder_count_preserved: 2,
    control_levers_preserved: true,
    spout_component_count_preserved: 1,
    three_hole_trim_disc_count_preserved: 3,
    formal_ifc_modified: false,
    pass: true,
  });
  expect(manifest.review_only_geometry_simplification.original_bounds_mm.minimum[0]).toBe(manifest.review_only_geometry_simplification.optimised_bounds_mm.minimum[0]);
  expect(manifest.review_only_geometry_simplification.original_bounds_mm.maximum[0]).toBe(manifest.review_only_geometry_simplification.optimised_bounds_mm.maximum[0]);
  expect(candidate.review_only_geometry_simplification).toEqual(manifest.review_only_geometry_simplification);
  expect(manifest.drawing_source).toMatchObject({ source_kind: reviewSourceKind, original_source_kind: sourceKind, source_dwg_sha256: dwgHash, official_cad_used: true, third_party_cad_used: false });
  expect(manifest.companion_45089_dwg_used_as_54294_geometry).toBeFalse();
  expect(candidate).toMatchObject({ source_kind: reviewSourceKind, original_source_kind: sourceKind, source_dwg_sha256: dwgHash, official_cad_used: true, unaltered_official_cad_used_as_review_representation: false, original_official_cad_evidence_preserved: true, third_party_cad_used: false, review_status: "visual_review_pending" });
  for (const view of manifest.views) {
    expect(view.drawing_line_source_kind).toBe(reviewSourceKind);
    expect(view.original_official_cad_path_count).toBe(pathCounts[view.view as keyof typeof pathCounts]);
    expect(view.review_simplified_path_count).toBe(reviewPathCounts[view.view as keyof typeof reviewPathCounts]);
    expect(view.blue_line_present).toBeTrue();
    expect(view.white_mask_present).toBeTrue();
    expect(view.ifc_dwg_compatibility.pass).toBeTrue();
    expect(candidate.views[view.view].original_official_native_dwg_paths_mm).toHaveLength(pathCounts[view.view as keyof typeof pathCounts]);
    expect(candidate.views[view.view].review_simplified_official_outline_paths_mm).toHaveLength(reviewPathCounts[view.view as keyof typeof reviewPathCounts]);
    expect(candidate.views[view.view].handle_line_texture_simplification.review_texture_detail_path_count).toBe(0);
    expect(candidate.views[view.view].handle_line_texture_simplification.handle_envelope_delta_mm).toEqual([0, 0]);
    const svg = readFileSync(join(root, view.svg), "utf8");
    expect(svg).toContain('class="original-highpoly"');
    expect(svg).toContain('class="simplified-proxy-silhouette geometry-derived"');
    expect(svg).toContain('class="official-reference-mask"');
    expect(svg).toContain('class="review-simplified-reference official-outline-derived"');
    expect(svg).toContain('data-unaltered-official-dwg="false"');
    expect(svg).toContain("#1677c8");
    expect(svg).toContain("Blue = de-textured review simplification based on official 54294 outline");
    expect(svg.indexOf('class="official-reference-mask"')).toBeLessThan(svg.indexOf('class="review-simplified-reference official-outline-derived"'));
    const evidence = readFileSync(join(root, view.original_official_evidence_svg), "utf8");
    expect(evidence).toContain('class="official-reference native-dwg original-unaltered-evidence"');
    expect(evidence).toContain('data-unaltered-official-dwg="true"');
  }
  expect(candidate.views.plan.handle_line_texture_simplification.original_texture_detail_path_count).toBe(1410);
  expect(candidate.views.front.handle_line_texture_simplification.original_texture_detail_path_count).toBe(0);
  expect(candidate.views.side.handle_line_texture_simplification.original_texture_detail_path_count).toBe(709);
  const side = manifest.views.find((view: any) => view.view === "side");
  expect(side.translation_alignment).toMatchObject({
    mode: "wall_plane_and_straight_spout_axis_translation_only",
    scale_applied: false,
    mirror_applied: false,
    official_straight_spout_wall_z_mm: [29.051758, 49.051758],
    official_straight_spout_axis_z_mm: 39.051758,
    official_wall_surface_x_mm: -18.818959,
    prior_axis_only_translation_mm: [9.816093, -21.051758],
    corrected_wall_and_axis_translation_mm: [28.613938, -21.051758],
    prior_wall_surface_residual_mm: -18.797845,
    prior_outlet_endpoint_residual_mm: -25.763878,
    corrected_anchor_residual_mm: [0, 0],
    pass: true,
  });
});

test("project drawings retain walls and furniture while de-textured official-outline review stays topmost", () => {
  const context = JSON.parse(readFileSync(join(product, "project-context-manifest.json"), "utf8"));
  expect(context).toMatchObject({ source_kind: reviewSourceKind, official_cad_used: true, unaltered_official_cad_used_as_review_representation: false, original_official_cad_evidence_preserved: true, review_texture_detail_path_count: 0, third_party_cad_used: false, project_context_retained: true, walls_and_surrounding_project_elements_retained: true, overlay_top_layer_with_white_mask: true, blue_line_present: true, pass: true });
  expect(context.views.map((view: any) => [view.view, view.candidate_view])).toEqual([["plan", "plan"], ["front", "side"], ["side", "front"]]);
  for (const view of context.views) {
    expect(view.overlay.source_dwg_sha256).toBe(dwgHash);
    expect(view.overlay.path_count).toBe(reviewPathCounts[view.candidate_view as keyof typeof reviewPathCounts]);
    expect(view.overlay.fit.uniform_scale_preserved).toBeTrue();
    expect(Math.max(...view.overlay.fit.bbox_absolute_delta_svg_units)).toBeLessThanOrEqual(view.overlay.fit.bbox_tolerance_svg_units);
    expect(sha256(join(root, view.review_preview))).toBe(view.review_preview_sha256);
    const svgPath = join(root, view.output);
    const svg = readFileSync(svgPath, "utf8");
    expect(svg).toContain("IfcWall");
    expect(svg).toContain("IfcSanitaryTerminal");
    expect(svg).toContain('class="official-reference-envelope-mask"');
    expect(svg).toContain('class="official-reference-mask"');
    expect(svg).toContain('class="review-simplified-reference official-outline-derived project-context-overlay"');
    const overlay = svg.slice(svg.indexOf('class="review-simplified-reference official-outline-derived project-context-overlay"'));
    expect(overlay.indexOf('class="official-reference-envelope-mask"')).toBeLessThan(overlay.indexOf('class="official-reference-mask"'));
    expect(overlay.indexOf('class="official-reference-mask"')).toBeLessThan(overlay.indexOf('class="review-simplified-reference official-outline-derived"'));
    for (const match of svg.matchAll(/(?:href|xlink:href)="([^"#]+)"/g)) {
      if (/^(?:https?:|data:)/.test(match[1])) continue;
      expect(existsSync(resolve(dirname(svgPath), match[1]))).toBeTrue();
    }
  }
});

test("Bonsai evidence is actual isolated IFC Body rendering from four saved cameras", () => {
  const evidence = JSON.parse(readFileSync(join(product, "bonsai-review-manifest.json"), "utf8"));
  expect(evidence).toMatchObject({ mode: "actual_bonsai_ifc_body_camera_render", representative_global_id: "2iKOL78$H0N9Yd9$ky3pW4", geometry_product_count: 1, whole_model_render: false, formal_ifc_bytes_unchanged: true });
  expect(evidence.bonsai_session.saved_active_representation).toBe("Body");
  expect(evidence.bonsai_session.saved_camera_count).toBe(4);
  expect(evidence.review_only_geometry_simplification).toMatchObject({
    mode: "review_only_smooth_handles_and_official_max_reach_spout",
    removed_component_count: 3,
    original_vertex_count: 21299,
    original_face_count: 41570,
    optimised_vertex_count: 1822,
    optimised_face_count: 3616,
    vertex_reduction_percent: 91.446,
    face_reduction_percent: 91.301,
    main_cylinder_outline_preserved: true,
    control_levers_preserved: true,
    spout_preserved: true,
    spout_rebuilt_to_official_54294_maximum_reach: true,
    three_hole_layout_preserved: true,
    formal_ifc_modified: false,
    pass: true,
  });
  expect(evidence.renders.map((render: any) => render.view)).toEqual(["plan", "front", "side", "iso"]);
  for (const render of evidence.renders) {
    expect(render.camera_type).toBe("ORTHO");
    expect(render.render_operation).toBe("bpy.ops.render.render(write_still=True)");
    expect(sha256(join(root, render.path))).toBe(render.sha256);
  }
  expect(sha256(join(root, evidence.isolated_ifc))).toBe(evidence.isolated_ifc_sha256);
  expect(evidence.isolated_ifc).toContain("Gessi316-45089-54294-bonsai-review-simplified.ifc");
  expect(evidence.actual_project_body_comparison).toMatchObject({
    mode: "unchanged_actual_project_ifc_body_side_comparison",
    pass: true,
  });
  expect(sha256(join(root, evidence.actual_project_body_comparison.side_render.path))).toBe(evidence.actual_project_body_comparison.side_render.sha256);
  const simplification = JSON.parse(readFileSync(join(product, "review-geometry-simplification.json"), "utf8"));
  expect(simplification).toMatchObject({
    mode: "review_only_smooth_handles_and_official_max_reach_spout",
    formal_ifc_bytes_unchanged: true,
    removed_vertex_count: 20387,
    removed_face_count: 39762,
    replacement_vertex_count: 910,
    replacement_face_count: 1808,
    spout_replacement_vertex_count: 650,
    spout_replacement_face_count: 1296,
    optimised_vertex_count: 1822,
    optimised_face_count: 3616,
    pass: true,
  });
  expect(sha256(join(root, simplification.review_ifc))).toBe(simplification.review_ifc_sha256);
});

test("scene SVG approval cannot write an IFC and formal IFC remains byte-identical", () => {
  const temporary = mkdtempSync(join(tmpdir(), "gessi316-54294-gate-"));
  const output = join(temporary, "forbidden.ifc");
  const script = join(root, "pipeline/scripts/gessi316_54294_drawing_ifc.py");
  const approval = JSON.parse(readFileSync(join(root, "pipeline/decisions/gessi316-54294-drawing-approval.json"), "utf8"));
  const manifest = join(product, "manifest.json");
  expect(approval).toMatchObject({
    status: "approved",
    reviewer: "project owner",
    review_date: "2026-08-27",
    approved_views: ["plan", "front", "side"],
    derived_ifc_write_allowed: false,
    formal_authoritative_ifc_write_allowed: false,
    review_stage: "main_bathroom_scene_svg_approved_ifc_write_closed",
    approval_scope_kind: "scene_svg_only",
    scope,
  });
  expect(approval.approved_scene_svgs).toMatchObject({
    manifest_sha256: "19ad23e5b8ac5a43cf4d2291f9e6c36e7e5c02ae807f1761aebbbce00ed31de8",
    source_kind: reviewSourceKind,
    source_label_zh: reviewSourceLabelZh,
    unaltered_official_cad: false,
    handle_texture_detail_path_count: 0,
  });
  expect(approval.approved_scene_svgs.views.plan.sha256).toBe("3724d63ac60ec6f2d6d4b2fda23f5a135e6113ad5a6455c4043dca27f6d9b0cc");
  expect(approval.approved_scene_svgs.views.front.sha256).toBe("bff75eef3805946758d803be6c54b6acd2d4330a48c6c1357a621f4969a568ba");
  expect(approval.approved_scene_svgs.views.side).toMatchObject({
    sha256: "b705eedd1a2a850a5a99ca41902459b55995881d8bec65925446526557dca286",
    annotation_translation_mm: [0, 0, 0],
  });
  expect(Math.abs(approval.approved_scene_svgs.views.side.finished_wall_residual_mm)).toBeLessThanOrEqual(0.2);
  expect(approval.candidate_manifest_sha256).toBe(sha256(manifest));
  const withoutApply = Bun.spawnSync(["python3", script, "--input", formal, "--output", output], { cwd: root, stderr: "pipe" });
  expect(withoutApply.exitCode).not.toBe(0);
  expect(withoutApply.stderr.toString()).toContain("IFC write requires the explicit --apply flag");
  const pending = Bun.spawnSync(["python3", script, "--input", formal, "--output", output, "--apply"], { cwd: root, stderr: "pipe" });
  expect(pending.exitCode).not.toBe(0);
  expect(pending.stderr.toString()).toContain("approval gate rejected IFC write");
  expect(existsSync(output)).toBeFalse();
  expect(sha256(formal)).toBe(formalHash);
  rmSync(temporary, { recursive: true, force: true });
});

test("temporary scoped approval writes open native-DWG representations and source associations", () => {
  const manifest = join(product, "manifest.json");
  const temporary = mkdtempSync(join(tmpdir(), "gessi316-54294-approved-"));
  const approvalPath = join(temporary, "approval.json");
  const output = join(temporary, "derived.ifc");
  const report = join(temporary, "report.json");
  writeFileSync(approvalPath, JSON.stringify({
    schema_version: 1, profile_key: "gessi316-54294", ifc_type_name: "Gessi316 54294",
    candidate_manifest_sha256: sha256(manifest), status: "approved", reviewer: "automated gate fixture",
    review_date: "2026-08-23", approved_views: ["plan", "front", "side"], derived_ifc_write_allowed: true,
    formal_authoritative_ifc_write_allowed: false, scope, approval_evidence: "temporary automated writer verification only",
  }, null, 2));
  const run = Bun.spawnSync(["python3", join(root, "pipeline/scripts/gessi316_54294_drawing_ifc.py"), "--input", formal, "--manifest", manifest, "--approval", approvalPath, "--output", output, "--report", report, "--apply"], { cwd: root, stdout: "pipe", stderr: "pipe" });
  if (run.exitCode !== 0) throw new Error(run.stderr.toString());
  const result = JSON.parse(readFileSync(report, "utf8"));
  expect(result.pass).toBeTrue();
  expect(result.formal_ifc_bytes_unchanged).toBeTrue();
  expect(result.representation_path_counts).toEqual(reviewPathCounts);
  expect(result.source_kind).toBe(reviewSourceKind);
  expect(result.source_label_zh).toBe(reviewSourceLabelZh);
  expect(result.official_cad_geometry_included).toBeFalse();
  expect(result.representation_geometry_source).toBe("review_simplified_official_outline_paths_mm based on GPF5429400000G000_3.dwg");
  expect(result.source_document_associations).toContain("GESSI316-45089-54294-OFFICIAL-NATIVE-DWG");
  expect(result.source_document_associations).toContain("GESSI316-45089-54294-OFFICIAL-TECHNICAL-PDF");
  expect(result.source_document_associations).toContain("GESSI316-45089-54294-SOURCE-ACCESS-RECORD");
  expect(sha256(formal)).toBe(formalHash);
  rmSync(temporary, { recursive: true, force: true });
}, 120_000);

test("actual main-bathroom Bonsai Drawings persist de-textured blue IFC LINEWORK and PDF previews", () => {
  const manifestPath = join(product, "main-bathroom-bonsai-drawing-manifest.json");
  const manifest = JSON.parse(readFileSync(manifestPath, "utf8"));
  expect(manifest).toMatchObject({
    pass: true,
    formal_ifc_sha256: formalHash,
    formal_ifc_bytes_unchanged: true,
    source_kind: reviewSourceKind,
    unaltered_official_dwg_used_as_drawing_linework: false,
    original_official_dwg_evidence_preserved: true,
    handle_texture_detail_path_count: 0,
    target_global_id: "2iKOL78$H0N9Yd9$ky3pW4",
    room_global_id: "3a4COIs5X7lgDirMBDT4Vs",
    room_name: "主卫湿区",
    project_context_retained: true,
  });
  expect(manifest.workflow).toEqual(["inspect", "plan", "execute", "persist", "reload", "verify"]);
  expect(manifest.provider).toMatchObject({
    status: "supported_for_inspection",
    inspection_project: "IFC4 / My Project",
    execution: "isolated local Blender 4.5.3 LTS / Bonsai-IfcOpenShell 0.8.4",
  });
  expect(manifest.combined_pdf.page_count).toBe(3);
  expect(manifest.combined_pdf.page_order).toEqual(["plan", "front", "side"]);
  expect(sha256(join(root, manifest.combined_pdf.path))).toBe(manifest.combined_pdf.sha256);
  expect(manifest.views.map((item: any) => item.view)).toEqual(["plan", "front", "side"]);
  const expectedCounts = { plan: 69, front: 1099, side: 35 };
  for (const item of manifest.views) {
    expect(item.create_drawing).toMatchObject({ operator: "bpy.ops.bim.create_drawing", result: ["FINISHED"], linework_mode: "OPENCASCADE" });
    expect(item.persistence).toMatchObject({ reload_result: ["FINISHED"], post_reload_drawing_found: true, post_reload_annotation_found: true });
    expect(item.review_path_count).toBe(expectedCounts[item.view as keyof typeof expectedCounts]);
    expect(item.persisted_review_path_count).toBe(expectedCounts[item.view as keyof typeof expectedCounts]);
    expect(item.persisted_coordinate_residual_mm).toBe(0);
    expect(item.mechanical_gate).toMatchObject({
      wall_anchor_residual_mm: 0,
      straight_spout_axis_residual_mm: 0,
      lower_lip_residual_mm: 0,
      outlet_endpoint_residual_mm: 0,
      visible_reach_residual_mm: 0,
      visible_envelope_residual_mm: [0, 0],
      pass: true,
    });
    for (const artifact of [item.drawing_session_ifc, item.svg, item.page_pdf, item.create_drawing_evidence, item.pdf_rendered_preview_png]) {
      expect(sha256(join(root, artifact.path))).toBe(artifact.sha256);
      expect(artifact.bytes).toBeGreaterThan(0);
    }
    expect(item.pdf_preview_colour_gate).toMatchObject({ blue_line_present: true, handle_texture_detail_path_count: 0, pass: true });
    const svg = readFileSync(join(root, item.svg.path), "utf8");
    expect(svg).toContain("review-target-gessi54294");
    expect(svg).toContain("native-dwg-review-simplification");
    expect(svg).toContain('data-unaltered-official-dwg="false"');
    expect(svg).toContain('data-handle-texture-detail-path-count="0"');
    expect(svg).toContain("stroke:#1677c8");
    expect(svg).not.toContain("original-unaltered-evidence");
  }
  const side = manifest.views.find((item: any) => item.view === "side");
  expect(side.side_scene_finished_wall_gate).toMatchObject({
    tolerance_mm: 0.2,
    coordinate_space: "actual project world coordinates and generated Side SVG",
    finished_wall_covering_global_id: "3WekjaeUn1qfL6oR1KYD_2",
    occluding_south_context_global_ids: [
      "04rs0EDjn2EvxytEQSxWRB",
      "0nlHIdaVrV7xIKsYXraYCa",
      "26PwmC1AX6ZROC65DohUg5",
    ],
    occluding_south_context_excluded_from_side_drawing: true,
    annotation_rigid_translation_applied_mm: [0, 0, 0],
    finished_wall_projection_group_count: 1,
    occluding_context_projection_group_count: 0,
    post_reload_annotation_target_placement_residual: 0,
    pass: true,
  });
  expect(Math.abs(side.side_scene_finished_wall_gate.annotation_to_finished_wall_residual_mm)).toBeLessThanOrEqual(0.2);
  expect(Math.abs(side.side_scene_finished_wall_gate.actual_body_to_finished_wall_residual_mm)).toBeLessThanOrEqual(0.2);
  expect(Math.abs(side.side_scene_finished_wall_gate.svg_residual_mm)).toBeLessThanOrEqual(0.2);
  expect(side.svg_validation.finished_wall_projection_group_count).toBe(1);
  expect(side.svg_validation.occluding_context_projection_group_count).toBe(0);
  expect(side.svg_validation.target_ifc_projection_group_count).toBe(0);
  const sideSvg = readFileSync(join(root, side.svg.path), "utf8");
  expect(sideSvg).toContain('ifc:guid="3WekjaeUn1qfL6oR1KYD_2"');
  expect(sideSvg).not.toContain('ifc:guid="04rs0EDjn2EvxytEQSxWRB"');
  expect(sideSvg).not.toContain('ifc:guid="0nlHIdaVrV7xIKsYXraYCa"');
  expect(sideSvg).not.toContain('ifc:guid="26PwmC1AX6ZROC65DohUg5"');
  expect(manifest.tests).toMatchObject({
    all_views_create_drawing_finished: true,
    all_sessions_reloaded: true,
    all_pngs_blue_line_present: true,
    all_pngs_handle_texture_detail_path_count: 0,
    all_persisted_coordinate_residual_mm: 0,
    side_finished_wall_scene_gate_pass: true,
    formal_ifc_hash_preserved: true,
  });
  expect(Math.abs(manifest.tests.side_finished_wall_scene_residual_mm)).toBeLessThanOrEqual(0.2);
  expect(sha256(formal)).toBe(formalHash);
}, 120_000);
