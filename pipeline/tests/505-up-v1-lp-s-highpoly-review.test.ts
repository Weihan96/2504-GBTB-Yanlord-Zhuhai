import { expect, test } from "bun:test";
import { createHash } from "node:crypto";
import { existsSync, mkdtempSync, readFileSync, rmSync, writeFileSync } from "node:fs";
import { tmpdir } from "node:os";
import { join, resolve } from "node:path";

const root = resolve(import.meta.dir, "../..");
const product = join(root, "output/review/highpoly-types/505-up-v1-lp-s");
const formal = join(root, "2504 GBTB Yanlord Zhuhai.ifc");
const formalHash = "7a50b87e8f48a7c2bfbaa9f1dabdfca7155fa684b2325c66aa1ae15c8ab0a25c";
const sourceKind = "geometry_derived_simplified_proxy";
const sourceLabelZh = "基于原始高模几何生成的简化图纸表达";
const scope = "exact manufacturer family and native family CAD; project-specific V1.LP.S composition is not mechanically matched to an official catalogue composition and is not a project shop drawing";

function sha256(path: string) {
  return createHash("sha256").update(readFileSync(path)).digest("hex");
}

test("official native 505 UP family DWGs are archived but do not masquerade as the project composition", () => {
  const access = JSON.parse(readFileSync(join(product, "official-source/source-access-record.json"), "utf8"));
  expect(access).toMatchObject({ manufacturer: "Molteni&C", family: "505 UP System", designer: "Nicola Gallizia", scope, pass: true });
  expect(access.official_product_cad).toMatchObject({ acquired: true, exact_project_configuration_match: false, official_cad_used_as_candidate_geometry: false });
  expect(access.configuration_cross_check).toMatchObject({ exact_match_count: 0, pass: true });
  expect(access.project_identity_evidence).toMatchObject({ family_identity_status: "confirmed", project_suffix_status: "unverified_project_or_asset_composition_suffix" });
  const dwgs = access.official_identity_sources.filter((source: any) => source.kind === "manufacturer_native_family_dwg");
  expect(dwgs.map((source: any) => source.sha256)).toEqual([
    "dfe4a6bcd655a3ed19343813a4aba8e139b847cafe8ea69b47623a6e1e71add3",
    "1e7988d5c7f00522736144fbecde1814032e5673af8cd13a47979e4347231646",
  ]);
  for (const source of dwgs) expect(sha256(join(root, source.local_path))).toBe(source.sha256);
  expect(access.native_dwg_conversion_audit.official_reference_svgs.map((item: any) => sha256(join(root, item.path)))).toEqual(
    access.native_dwg_conversion_audit.official_reference_svgs.map((item: any) => item.sha256),
  );
  expect(access.drawing_geometry_source).toMatchObject({ source_kind: sourceKind, source_label_zh: sourceLabelZh, official_cad_used: false, third_party_cad_used: false });
});

test("one 505 UP Body produces three geometry-derived candidates with zero blue CAD paths", () => {
  const manifest = JSON.parse(readFileSync(join(product, "manifest.json"), "utf8"));
  const candidate = JSON.parse(readFileSync(join(product, "candidate-representations.json"), "utf8"));
  expect(manifest).toMatchObject({ representative_global_id: "19MpdkWqXC7uhUNhLQgrce", geometry_product_count: 1, whole_model_render: false, formal_ifc_bytes_unchanged: true, review_status: "visual_review_pending", approved_for_drawing_ifc: false, pass: true });
  expect(manifest.project_context.manifest_sha256).toBe(sha256(join(root, manifest.project_context.manifest)));
  expect(manifest.bonsai_review.manifest_sha256).toBe(sha256(join(root, manifest.bonsai_review.manifest)));
  expect(candidate).toMatchObject({ article_number: "505 UP System / project V1.LP.S", source_kind: sourceKind, source_label_zh: sourceLabelZh, official_cad_used: false, third_party_cad_used: false, formal_ifc_write_allowed: false });
  expect(Object.fromEntries(Object.entries(candidate.views).map(([view, item]: any) => [view, item.proxy_paths_mm.length]))).toEqual({ plan: 24, front: 281, side: 93 });
  expect(candidate.views.plan.component_semantics).toMatchObject({
    method: "plan_first_visible_highest_face_semantic_partition",
    camera_direction_world: [0, 0, -1],
    local_to_world_z_sign: -1,
    normalization_tolerance_mm: 0.5,
    two_side_probe_mm: 2,
    required_internal_interfaces: {
      "A_main_top/C_front_left_display": 2,
      "A_main_top/D_mid_shelf_lip": 1,
      "C_front_left_display/E_lower_front_rails": 2,
      "D_mid_shelf_lip/E_lower_front_rails": 6,
    },
    rejected_paths_removed: ["P02", "P03", "P04", "P05"],
    bottom_76mm_components_used_as_top: false,
    every_internal_segment_has_different_semantics_on_both_sides: true,
    pass: true,
  });
  expect(sha256(join(product, "front.svg"))).toBe("0f942b3799c9854119765445258a8eb422274c5cb198582e2fe369c5b13c975b");
  expect(sha256(join(product, "side.svg"))).toBe("8f4224ef5fb3594bd369cc850ff49bd728c266bba8ac966528a1c361ed3eb156");
  for (const view of ["plan", "front", "side"]) {
    expect(candidate.views[view].official_cad_paths_mm).toEqual([]);
    const svg = readFileSync(join(product, `${view}.svg`), "utf8");
    expect(svg).toContain(`data-source-kind="${sourceKind}"`);
    expect(svg).not.toContain('class="official-reference native-dwg"');
    expect(svg).not.toContain("#1677c8");
  }
});

test("project plan retains real walls and furniture while no project elevation is invented", async () => {
  const context = JSON.parse(readFileSync(join(product, "project-context-manifest.json"), "utf8"));
  expect(context.views.map((view: any) => view.view)).toEqual(["plan"]);
  expect(context.context_view_scope).toMatchObject({ included: ["plan"], excluded: ["front", "side"] });
  expect(context).toMatchObject({ project_context_retained: true, walls_and_surrounding_project_elements_retained: true, overlay_top_layer_with_white_mask: true, blue_line_present: false, pass: true });
  const view = context.views[0];
  expect(view.source).toBe("drawings/Furniture Plan.svg");
  expect(view.overlay.fit).toMatchObject({ rotate_quarter_turns: 0, uniform_scale_preserved: true, scale_svg_units_per_mm: 0.02, pass: true });
  expect(Math.max(...view.overlay.fit.bbox_absolute_delta_svg_units)).toBeLessThanOrEqual(view.overlay.fit.bbox_tolerance_svg_units);
  expect(sha256(join(root, view.review_preview))).toBe(view.review_preview_sha256);
  const svg = readFileSync(join(root, view.output), "utf8");
  expect(svg).toContain('class="geometry-derived-simplified-proxy project-context-overlay"');
  expect(svg.indexOf('class="geometry-derived-proxy-mask"')).toBeLessThan(svg.indexOf('class="geometry-derived-proxy"'));
  const nativeElevations = new Bun.Glob("drawings/elevations/native/*.svg");
  for await (const path of nativeElevations.scan(root)) expect(readFileSync(join(root, path), "utf8")).not.toContain("19MpdkWqXC7uhUNhLQgrce");
});

test("actual Bonsai session contains four camera renders of the isolated MODEL_VIEW Body", () => {
  const evidence = JSON.parse(readFileSync(join(product, "bonsai-review-manifest.json"), "utf8"));
  expect(evidence).toMatchObject({ mode: "actual_bonsai_ifc_body_camera_render", representative_global_id: "19MpdkWqXC7uhUNhLQgrce", geometry_product_count: 1, whole_model_render: false, formal_ifc_bytes_unchanged: true, pass: true });
  expect(evidence.bonsai_session).toMatchObject({ saved_active_representation: "Body", ifc_context_identifier: "Body", ifc_target_view: "MODEL_VIEW", saved_camera_count: 4, front_camera_local_y_sign: 1 });
  expect(evidence.renders.map((render: any) => render.view)).toEqual(["plan", "front", "side", "iso"]);
  for (const render of evidence.renders) {
    expect(render.render_operation).toBe("bpy.ops.render.render(write_still=True)");
    expect(sha256(join(root, render.path))).toBe(render.sha256);
  }
  expect(sha256(join(root, evidence.bonsai_session.path))).toBe(evidence.bonsai_session.sha256);
  expect(existsSync(join(product, "Molteni-505-UP-V1-LP-S-bonsai-review.blend1"))).toBeFalse();
});

test("approved Front is frozen while rejected Plan keeps only the existing derived IFC gate", () => {
  const approval = JSON.parse(readFileSync(join(root, "pipeline/decisions/505-up-v1-lp-s-drawing-approval.json"), "utf8"));
  expect(approval).toMatchObject({
    status: "revision_pending_review",
    reviewer: "project_owner",
    approved_views: ["plan", "front", "side"],
    approval_evidence: "先黑色的505，下一个吧",
    derived_ifc_write_authorization_evidence: "505 必须写入 IFC",
    derived_ifc_write_allowed: true,
    formal_authoritative_ifc_write_allowed: false,
    pending_reapproval_views: ["plan"],
    latest_review: {
      outcome: "new_single_product_plan_candidate_pending_review",
      front_outcome: "approved_and_frozen",
      plan_outcome: "new_semantic_candidate_generated_pending_owner_review",
      plan_user_authorization: "没问题，那么基于此给我新的单品svg",
      derived_ifc_updated: false,
      bonsai_session_updated: false,
      create_drawing_called: false,
    },
    front_revision_freeze: {
      status: "approved_and_frozen",
      scene_svg_sha256: "85331ec1c29c42e053f06dce7fcf58a68f4cba8d594cb0044d28c0af1a6f10bf",
      rendered_png_sha256: "8ff7aa9aab8ae33789f9b504ebbea63cbbf9910f75f08a8c581be08fd37361d9",
      linework_path_count: 281,
    },
    approved_candidate: {
      source_kind: sourceKind,
      source_label_zh: sourceLabelZh,
      official_cad_used: false,
      blue_official_atomic_component_composition_selected: false,
    },
  });
  expect(approval.candidate_manifest_sha256).toBe(sha256(join(product, "manifest.json")));
  expect(approval.plan_semantic_candidate).toMatchObject({
    status: "pending_owner_review",
    source_kind: sourceKind,
    path_count: 24,
    every_internal_segment_has_different_semantics_on_both_sides: true,
    derived_ifc_updated: false,
  });
  expect(sha256(formal)).toBe(formalHash);
});

test("new Plan single-product SVG contains only mechanically valid first-visible boundaries", () => {
  const audit = JSON.parse(readFileSync(join(product, "505-up-plan-semantic-candidate-audit.json"), "utf8"));
  expect(audit).toMatchObject({
    status: "new_single_product_plan_svg_pending_owner_review",
    user_authorization: "没问题，那么基于此给我新的单品svg",
    source_kind: sourceKind,
    plan_visibility_rule: {
      camera_direction_world: [0, 0, -1],
      local_to_world_z_sign: -1,
      tessellation_gap_normalization_tolerance_mm: 0.5,
      bottom_76mm_components_used_as_top: false,
    },
    removed_rejected_paths: ["P02", "P03", "P04", "P05"],
    candidate: {
      path_count: 24,
      segment_count: 87,
      internal_segment_count: 11,
      exterior_segment_count: 76,
      invalid_same_side_or_unexplained_count: 0,
      required_internal_interface_counts: {
        "A_main_top/C_front_left_display": 2,
        "A_main_top/D_mid_shelf_lip": 1,
        "C_front_left_display/E_lower_front_rails": 2,
        "D_mid_shelf_lip/E_lower_front_rails": 6,
      },
    },
    write_boundary: {
      single_product_plan_svg_updated: true,
      derived_ifc_updated: false,
      bonsai_session_updated: false,
      project_scene_drawing_updated: false,
      create_drawing_called: false,
      front_svg_frozen: true,
      side_svg_frozen: true,
    },
    pass: true,
  });
  expect(audit.candidate.segment_semantics.every((segment: any) => segment.side_a !== segment.side_b)).toBeTrue();
  expect(audit.candidate.segment_semantics.filter((segment: any) => segment.internal).every(
    (segment: any) => segment.side_a !== "background" && segment.side_b !== "background",
  )).toBeTrue();
  for (const [pathKey, hashKey] of [
    ["plan_preview_png", "plan_preview_png_sha256"],
    ["contact_sheet_png", "contact_sheet_png_sha256"],
  ]) expect(sha256(join(root, audit.outputs[pathKey]))).toBe(audit.outputs[hashKey]);
  const svg = readFileSync(join(product, "plan.svg"), "utf8");
  expect(svg).toContain('data-visibility-method="first-visible-highest-face"');
  expect(svg).toContain('data-semantic-two-side-validation="passed"');
  expect(svg).toContain('data-bottom-76mm-as-top="false"');
  expect(svg).not.toContain("#1677c8");
  expect(sha256(join(product, "front.svg"))).toBe("0f942b3799c9854119765445258a8eb422274c5cb198582e2fe369c5b13c975b");
  expect(sha256(join(product, "side.svg"))).toBe("8f4224ef5fb3594bd369cc850ff49bd728c266bba8ac966528a1c361ed3eb156");
  expect(sha256(join(product, "505-up-v1-lp-s-derived-drawing.ifc"))).toBe("a49d1b9d2859fe596761ecc70c845cf888022f32c988a71946e8bfa4b0338417");
  expect(sha256(join(product, "505-up-v1-lp-s-bonsai-drawing-session.ifc"))).toBe("b8e9c5825225bd5f1df3a5ca5b5c9d0945e29acc9f1435c1e4dc98c3525b7687");
  expect(sha256(formal)).toBe(formalHash);
});

test("Plan diagnostic exposes depth-visible semantic regions without changing final Drawing artifacts", () => {
  const audit = JSON.parse(readFileSync(join(product, "505-up-plan-visibility-semantics.json"), "utf8"));
  expect(audit).toMatchObject({
    status: "diagnostic_only_plan_revision_required",
    component_count: 38,
    plan_camera_visibility_rule: {
      camera_direction_world: [0, 0, -1],
      local_to_world_z_sign: -1,
      finding: "local z near 0 is world floor level, not the product top",
    },
    decision: {
      current_plan_candidate_pass: false,
      final_plan_updated: false,
      derived_ifc_updated: false,
      bonsai_session_updated: false,
      create_drawing_called: false,
    },
    frozen_views: {
      plan_rejected_scene_svg_sha256: "f67a0a3867e105e6ffc2a4507ba07f40d233c72142d46c5c06e2c9555ccf29d5",
      front_scene_svg_sha256: "85331ec1c29c42e053f06dce7fcf58a68f4cba8d594cb0044d28c0af1a6f10bf",
      side_scene_svg_sha256: "e5c16902267f3ba56157230aca3ff7b234b446f7e47b5a7b7cecc248c3643c84",
      derived_ifc_sha256: "a49d1b9d2859fe596761ecc70c845cf888022f32c988a71946e8bfa4b0338417",
      bonsai_session_ifc_sha256: "b8e9c5825225bd5f1df3a5ca5b5c9d0945e29acc9f1435c1e4dc98c3525b7687",
    },
    pass: true,
  });
  const regions = Object.fromEntries(
    audit.visible_surface_records.map((record: any) => [record.semantic_region, record.top_world_height_mm]),
  );
  expect(regions).toMatchObject({
    A_main_top: 2380,
    B_right_slats: 1976,
    C_front_left_display: 1610.705,
    D_mid_shelf_lip: 470,
    E_lower_front_rails: 440,
  });
  expect(audit.candidate_path_semantics.slice(1, 5).every((item: any) => item.classification === "unexplained_numeric_residual_loop")).toBeTrue();
  expect(audit.candidate_path_semantics.at(-1)).toMatchObject({ path_id: "P20", classification: "main_top_to_display_depth_boundary" });
  expect(audit.candidate_segment_summary).toEqual({
    segment_count: 97,
    valid_two_side_semantic_count: 81,
    invalid_same_side_or_unexplained_count: 16,
    invalid_segment_ids: [
      "P02-S01", "P02-S02", "P02-S03", "P02-S04",
      "P03-S01", "P03-S02", "P03-S03", "P03-S04",
      "P04-S01", "P04-S02", "P04-S03", "P04-S04",
      "P05-S01", "P05-S02", "P05-S03", "P05-S04",
    ],
  });
  expect(audit.candidate_segment_semantics.filter((item: any) => !item.semantic_separation_valid).every(
    (item: any) => item.side_a === "A_main_top" && item.side_b === "A_main_top",
  )).toBeTrue();
  for (const [pathKey, hashKey] of [
    ["component_visibility_svg", "component_visibility_svg_sha256"],
    ["depth_visibility_svg", "depth_visibility_svg_sha256"],
    ["component_visibility_png", "component_visibility_png_sha256"],
    ["depth_visibility_png", "depth_visibility_png_sha256"],
  ]) expect(sha256(join(root, audit.outputs[pathKey]))).toBe(audit.outputs[hashKey]);
  expect(sha256(formal)).toBe(formalHash);
});

test("persisted product-level IFC and revised Bonsai drawings reload with review path counts", () => {
  const derived = join(product, "505-up-v1-lp-s-derived-drawing.ifc");
  const session = join(product, "505-up-v1-lp-s-bonsai-drawing-session.ifc");
  const report = JSON.parse(readFileSync(join(product, "derived-drawing-write-report.json"), "utf8"));
  const evidence = JSON.parse(readFileSync(join(product, "505-UP-ENTRANCE-create-drawing-evidence.json"), "utf8"));
  const drawingManifest = JSON.parse(readFileSync(join(product, "505-up-project-drawing-manifest.json"), "utf8"));
  expect(report).toMatchObject({
    pass: true,
    formal_ifc_bytes_unchanged: true,
    representation_path_counts: { plan: 20, front: 281, side: 93 },
    official_cad_geometry_included: false,
    source_kind: sourceKind,
  });
  expect(sha256(derived)).toBe(report.derived_ifc_sha256);
  expect(evidence).toMatchObject({
    provider: { name: "bonsai-mcp", status: "supported", execution: "execute_blender_code" },
    formal_ifc_bytes_unchanged: true,
    source_kind: sourceKind,
    official_cad_used: false,
    blue_official_atomic_component_composition_selected: false,
    target: { global_id: "19MpdkWqXC7uhUNhLQgrce", target_include_count_after_suppression: 0 },
    room: { name: "玄关" },
    pass: true,
  });
  expect(sha256(session)).toBe(evidence.drawing_session_sha256_after);
  expect(drawingManifest.approval).toMatchObject({ status: "revision_pending_review", pending_reapproval_views: ["plan"] });
  expect(sha256(join(root, drawingManifest.approval.record))).not.toBe(drawingManifest.approval.record_sha256);
  expect(drawingManifest.approval.status).toBe("revision_pending_review");
  expect(sha256(join(root, drawingManifest.derived_ifc.path))).toBe(drawingManifest.derived_ifc.sha256);
  expect(sha256(join(root, drawingManifest.bonsai_session.ifc))).toBe(drawingManifest.bonsai_session.ifc_sha256);
  expect(sha256(join(root, drawingManifest.revision_audit.path))).toBe(drawingManifest.revision_audit.sha256);
  for (const view of drawingManifest.views) {
    expect(sha256(join(root, view.svg))).toBe(view.svg_sha256);
    expect(sha256(join(root, view.rendered_pdf_page))).toBe(view.rendered_pdf_page_sha256);
  }
  const reload = Bun.spawnSync([
    "python3",
    "-c",
    [
      "import ifcopenshell,json,sys",
      "f=ifcopenshell.open(sys.argv[1])",
      "drawings=[x for x in f.by_type('IfcAnnotation') if getattr(x,'ObjectType',None)=='DRAWING' and (getattr(x,'Name',None) or '').startswith('MOLTENI-505-UP-ENTRANCE-')]",
      "linework=[x for x in f.by_type('IfcAnnotation') if getattr(x,'ObjectType',None)=='LINEWORK' and (getattr(x,'Name',None) or '').startswith('Molteni 505 UP approved black proxy')]",
      "counts={x.Name.rsplit('/',1)[-1].strip():sum(len(item.Elements or []) for rep in x.Representation.Representations or [] for item in rep.Items or [] if item.is_a('IfcGeometricCurveSet')) for x in linework}",
      "print(json.dumps({'schema':f.schema,'drawings':len(drawings),'linework':len(linework),'counts':counts}))",
    ].join(";"),
    session,
  ], { cwd: root, stdout: "pipe", stderr: "pipe" });
  if (reload.exitCode !== 0) throw new Error(reload.stderr.toString());
  expect(JSON.parse(reload.stdout.toString())).toEqual({
    schema: "IFC4",
    drawings: 3,
    linework: 3,
    counts: { plan: 20, front: 281, side: 93 },
  });
  expect(evidence.views.map((view: any) => [view.view, view.persisted_review_path_count])).toEqual([
    ["plan", 20],
    ["front", 281],
    ["side", 93],
  ]);
  for (const view of evidence.views) {
    expect(view.create_drawing).toMatchObject({ operator: "bpy.ops.bim.create_drawing", result: ["FINISHED"] });
    expect(view.persisted_coordinate_residual_mm).toBe(0);
    expect(view.svg.post_style_only).toBe(true);
    expect(view.svg.black_geometry_element_count).toBeGreaterThan(0);
    expect(view.svg.grey_context_geometry_element_count).toBeGreaterThan(0);
    expect(sha256(view.svg.path)).toBe(view.svg.sha256);
  }
  expect(evidence).toMatchObject({ review_status: "revision_pending_review" });
  expect(evidence.plan_component_semantics).toMatchObject({
    method: "connected_components_then_display_main_envelope_intersection",
    slat_components: { count: 14 },
    interface: {
      path_mm: [[-446.18847, 320], [163.643711, 320]],
      closed: false,
      full_display_footprint_added: false,
    },
    pass: true,
  });
  expect(evidence.front_orientation_revision).toMatchObject({
    method: "product_world_axes_plus_room_facing_component_depth_plus_camera_basis",
    local_x_screen_dot_after: 1,
    front_outward_view_dot_after: -1,
    persisted_north_rotation_pass: true,
    pass: true,
  });
  expect(evidence.plan_front_revision).toMatchObject({
    create_results: { plan: ["FINISHED"], front: ["FINISHED"] },
    side_unchanged: true,
    target_body_projection_counts: { plan: 0, front: 0 },
    svg_sha256_after: {
      plan: "f67a0a3867e105e6ffc2a4507ba07f40d233c72142d46c5c06e2c9555ccf29d5",
      front: "85331ec1c29c42e053f06dce7fcf58a68f4cba8d594cb0044d28c0af1a6f10bf",
      side: "e5c16902267f3ba56157230aca3ff7b234b446f7e47b5a7b7cecc248c3643c84",
    },
  });
  expect(existsSync(evidence.bonsai_session.path)).toBeTrue();
  expect(sha256(evidence.bonsai_session.path)).toBe(evidence.bonsai_session.sha256);
  expect(sha256(formal)).toBe(formalHash);
});

test("scoped temporary approval writes verified representations and native-DWG source associations", () => {
  const manifest = join(product, "manifest.json");
  const temporary = mkdtempSync(join(tmpdir(), "molteni-505-approved-"));
  const approvalPath = join(temporary, "approval.json");
  const output = join(temporary, "derived.ifc");
  const report = join(temporary, "report.json");
  writeFileSync(approvalPath, JSON.stringify({ schema_version: 1, profile_key: "505-up-v1-lp-s", ifc_type_name: "505 UP V1.LP.S", candidate_manifest_sha256: sha256(manifest), status: "approved", reviewer: "automated gate fixture", review_date: "2026-08-23", approved_views: ["plan", "front", "side"], derived_ifc_write_allowed: true, formal_authoritative_ifc_write_allowed: false, scope, approval_evidence: "temporary automated writer verification only" }, null, 2));
  const run = Bun.spawnSync(["python3", join(root, "pipeline/scripts/505_up_v1_lp_s_drawing_ifc.py"), "--input", formal, "--manifest", manifest, "--approval", approvalPath, "--output", output, "--report", report, "--apply"], { cwd: root, stdout: "pipe", stderr: "pipe" });
  if (run.exitCode !== 0) throw new Error(run.stderr.toString());
  const result = JSON.parse(readFileSync(report, "utf8"));
  expect(result).toMatchObject({ pass: true, formal_ifc_bytes_unchanged: true, representation_path_counts: { plan: 24, front: 281, side: 93 }, official_cad_geometry_included: false, source_kind: sourceKind, source_label_zh: sourceLabelZh });
  expect(result.representations).toEqual({ plan: "Molteni505UpPlan", front: "Molteni505UpFront", side: "Molteni505UpSide" });
  expect(result.source_document_associations).toContain("MOLTENI-505-UP-OFFICIAL-TECHNICAL-DWG");
  expect(result.source_document_associations).toContain("MOLTENI-505-UP-OFFICIAL-INSPIRING-DWG");
  expect(sha256(formal)).toBe(formalHash);
  rmSync(temporary, { recursive: true, force: true });
}, 30_000);
