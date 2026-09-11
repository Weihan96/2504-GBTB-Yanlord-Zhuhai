import { expect, test } from "bun:test";
import { createHash } from "node:crypto";
import { existsSync, readFileSync } from "node:fs";
import { join, resolve } from "node:path";

const root = resolve(import.meta.dir, "../..");
const product = join(root, "output/review/highpoly-types/miamisoft-e09");
const formal = join(root, "2504 GBTB Yanlord Zhuhai.ifc");
const formalHash = "7a50b87e8f48a7c2bfbaa9f1dabdfca7155fa684b2325c66aa1ae15c8ab0a25c";
const sourceKind = "native_dwg";
const sourceLabelZh = "Baxter 官方 Miami Soft E09 原生二维 DWG 候选";
const scope = "Exact Baxter Miami Soft E09 dx/r native Plan/Front/Side selected from the official family DWG; 1:1, translation only, with Side view-direction reflection; pending review and not yet written to IFC.";
const zipHash = "f422e836ed7d666616340d696af93002fac6c9deff5c9e4a714ea9b8e1cbafd8";
const dwgHash = "efc3315feed0a2dac8938aebbd7336856dabba6721395e4520783c39ba0cc51a";

function sha256(path: string) {
  return createHash("sha256").update(readFileSync(path)).digest("hex");
}

function sha256Text(value: string) {
  return createHash("sha256").update(value).digest("hex");
}

test("exact Baxter Miami Soft E09 native DWG is archived with source hashes", () => {
  const access = JSON.parse(readFileSync(join(product, "official-source/source-access-record.json"), "utf8"));
  expect(access).toMatchObject({
    manufacturer: "Baxter",
    family: "Miami Soft",
    designer: "Paola Navone",
    project_ifc_type_name: "MiamiSoft E09",
    project_ifc_type_description: "R terminal module 130 x 108 h70/80 cm",
    resolved_variant: "E09 - dx/r - right terminal module - 130 x 108 h70/80 cm",
    scope,
    pass: true,
  });
  expect(access.project_identity_evidence).toMatchObject({
    exact_model_code_status: "confirmed",
    handed_variant_status: "confirmed_E09_dx_r_right_terminal",
  });
  expect(access.official_product_cad).toMatchObject({
    authentication_required_on_user_download_page: true,
    official_asset_url_returned_by_product_drawings_api: true,
    acquired: true,
    exact_project_configuration_match: true,
    family_archive_sha256: zipHash,
    third_party_cad_used: false,
    official_vector_pdf_archived: true,
    official_measurement_svg_archived: true,
    official_vector_references_used_as_cad_geometry: false,
  });
  expect(access.official_product_cad.local_cad_files).toEqual([expect.objectContaining({ sha256: dwgHash, format: "AutoCAD 2018/2019/2020 DWG" })]);
  expect(access.dimension_cross_check).toMatchObject({
    official_nominal_width_depth_height_mm: [1300, 1080, 800],
    project_ifc_body_local_xyz_mm: [1343.172913, 1144.75415, 797.479919],
    maximum_absolute_delta_mm: 64.75415,
    tolerance_mm: 70,
    pass: true,
  });
  const pdf = access.official_identity_sources.find((source: any) => source.kind === "manufacturer_vector_technical_sheet");
  expect(pdf.page).toBe(6);
  expect(pdf.vector_object_audit).toMatchObject({ drawing_count: 8000, drawing_item_count: 44568, word_count: 438, character_count: 2035 });
  expect(pdf.vector_object_audit.required_text_found).toEqual(["E09 - dx/r", "R terminal module.", "130 x 108 h70/80 cm"]);
  const svg = access.official_identity_sources.find((source: any) => source.kind === "manufacturer_measurement_svg");
  expect(svg.vector_object_audit).toEqual({ path_count: 466, use_count: 76, group_count: 66 });
  for (const source of access.official_identity_sources) {
    expect(sha256(join(root, source.local_path))).toBe(source.sha256);
    if (source.preview_path) expect(sha256(join(root, source.preview_path))).toBe(source.preview_sha256);
  }
  const inventory = JSON.parse(readFileSync(join(product, "official-source/official-download/archive-inventory.json"), "utf8"));
  expect(inventory).toMatchObject({ archive_sha256: zipHash, selected_native_2d_sha256: dwgHash, member_count: 14 });
  expect(sha256(join(root, inventory.archive))).toBe(zipHash);
  expect(sha256(join(root, inventory.selected_native_2d_local_path))).toBe(dwgHash);
});

test("one Miami Soft E09 Body is compared with exact native-DWG blue candidates", () => {
  const manifest = JSON.parse(readFileSync(join(product, "manifest.json"), "utf8"));
  const candidate = JSON.parse(readFileSync(join(product, "candidate-representations.json"), "utf8"));
  const inventory = JSON.parse(readFileSync(join(root, "pipeline/decisions/highpoly-product-inventory.json"), "utf8"));
  const inventoryEntry = inventory.products.find((entry: any) => entry.type_name === "MiamiSoft E09");
  expect(inventoryEntry.status).toBe("review_ready_pending_approval");
  expect(inventory.review_ready_pending_approval_product_count).toBe(38);
  expect(manifest).toMatchObject({
    representative_global_id: "3osoWufdD1mhDDAM6lcix4",
    geometry_product_count: 1,
    whole_model_render: false,
    source_kind: sourceKind,
    official_cad_acquired: true,
    official_cad_exact_project_configuration_match: true,
    official_cad_used: true,
    blue_line_present: true,
    pass: true,
  });
  expect(manifest.bounds_mm.size).toEqual([1343.172913, 1144.75415, 797.479919]);
  expect(candidate).toMatchObject({
    article_number: "Baxter Miami Soft E09 dx/r",
    source_kind: sourceKind,
    source_label_zh: sourceLabelZh,
    official_cad_used: true,
    third_party_cad_used: false,
    formal_ifc_write_allowed: false,
    derived_ifc_write_allowed: true,
    review_status: "approved_for_product_level_derived_ifc",
  });
  expect(Object.fromEntries(Object.entries(candidate.views).map(([view, value]: any) => [view, value.candidate_paths_mm.length]))).toEqual({ plan: 64, front: 130, side: 85 });
  expect(candidate.views.plan.comparison_geometry_derived_proxy_paths_mm).toHaveLength(25);
  expect(candidate.views.plan.comparison_geometry_derived_proxy_provenance).toMatchObject({
    source_kind: "geometry_derived_simplified_proxy",
    source_revision: "139d9ec178889e4f6bc0b81095b9bc79f33b3b79",
    source_candidate_sha256: "bfd6da993c754be473ee0296bd356ba7597c2918fac1e4372445268a3b90f305",
    path_count: 25,
    recomputed: false,
    official_blue_paths_reused_as_black: false,
  });
  expect(manifest.views.find((view: any) => view.view === "plan")).toMatchObject({
    comparison_proxy_path_count: 25,
    official_cad_path_count: 64,
    layer_path_counts: {
      actual_ifc_body_displayed_edges: 989,
      geometry_derived_simplified_proxy: 25,
      official_native_dwg: 64,
    },
  });
  for (const view of ["plan", "front", "side"]) {
    expect(candidate.views[view].official_cad_paths_mm).toEqual(candidate.views[view].candidate_paths_mm);
    expect(candidate.views[view].alignment.uniform_scale).toBe(1);
    expect(candidate.views[view].alignment.anisotropic_scale_used).toBeFalse();
    expect(candidate.views[view].alignment.source_geometry_deformed).toBeFalse();
    const svgText = readFileSync(join(product, `${view}.svg`), "utf8");
    expect(svgText).toContain(`data-source-kind="${sourceKind}"`);
    expect(svgText).toContain('class="official-candidate native-dwg"');
    expect(svgText).toContain('data-source-scaled="false"');
    expect(svgText).toContain("#1677c8");
  }
  expect(candidate.views.plan.alignment.view_direction_reflection_x).toBeFalse();
  expect(candidate.views.front.alignment.view_direction_reflection_x).toBeFalse();
  expect(candidate.views.side.alignment.view_direction_reflection_x).toBeTrue();
  const planSvg = readFileSync(join(product, "plan.svg"), "utf8");
  const layerPath = (className: string) => planSvg.match(new RegExp(`<path class="${className}"[^>]*?\\sd="([^"]*)"`))?.[1] ?? "";
  expect(sha256Text(layerPath("actual-ifc-body"))).toBe("268e5010a56c43438f246fd9c6bf174a9920616f4d5b78588887b22c922b49a6");
  expect(sha256Text(layerPath("geometry-derived-simplified-proxy comparison-only"))).toBe("2f60dd288c55b9cb842940f0f5e31a7bdf1e752b2fbdc8889fcbd9ca5474c000");
  expect(sha256Text(layerPath("official-candidate native-dwg"))).toBe("be78291dda48329b490abb01d4823c58f1fc9aa178cd1028a867fe64b5eeff14");
});

test("Miami Soft project context places exact official blue line at the IFC world transform", () => {
  const context = JSON.parse(readFileSync(join(product, "project-context-manifest.json"), "utf8"));
  const officialPathSemanticHashes: Record<string, string> = {
    plan: "d820eaf8c4f637d597be8eb9ad5afaeec715c1226eef959e297cf22ccc65ea0a",
    side: "3a026cb5f0fc1c247a3c7aadc2cf1b72c2654c0a50dce8650492dafba08efcb5",
  };
  expect(context.views.map((view: any) => view.view)).toEqual(["plan", "side"]);
  expect(context.context_view_scope).toMatchObject({ included: ["plan", "side"], excluded: ["front"] });
  expect(context).toMatchObject({
    schema_version: 2,
    artifact_kind: "review_only_project_context_overlay",
    not_bonsai_scene_drawing_svg: true,
    ifc_write_performed: false,
    derived_ifc_write_allowed: false,
    formal_ifc_write_allowed: false,
    formal_ifc_sha256: formalHash,
    formal_ifc_bytes_unchanged: true,
    source_kind: "native_dwg",
    official_cad_used: true,
    project_context_retained: true,
    walls_and_surrounding_project_elements_retained: true,
    overlay_top_layer_with_white_mask: false,
    blue_line_present: true,
    geometry_scale: 1,
    alignment_by_bbox_fit: false,
    pass: true,
  });
  expect(context.provider_pre_state).toMatchObject({ online: true, ifc_schema: "IFC4", project_name: "My Project", query_only: true, provider_mutation_used: false });
  const plan = context.views.find((view: any) => view.view === "plan");
  const side = context.views.find((view: any) => view.view === "side");
  expect(plan.overlay).toMatchObject({
    source_kind: "native_dwg",
    source_path_count: 64,
    blue_line_present: true,
    uniform_geometry_scale: 1,
    recorded_view_direction_reflection_x: false,
    alignment_by_bbox_fit: false,
    rigid_world_placement_only: true,
    within_project_viewbox: true,
    no_clipping: true,
  });
  expect(side.overlay).toMatchObject({
    source_kind: "native_dwg",
    source_path_count: 85,
    blue_line_present: true,
    uniform_geometry_scale: 1,
    recorded_view_direction_reflection_x: true,
    reflection_already_applied_in_candidate: true,
    additional_reflection_applied_in_project_context: false,
    alignment_by_bbox_fit: false,
    rigid_world_placement_only: true,
    within_project_viewbox: true,
    no_clipping: true,
  });
  expect(plan.overlay.object_placement_matrix_mm).toEqual([
    [0, 1, 0, 2287.878274918],
    [-1, 0, 0, 803.948163986],
    [0, 0, 1, 10],
    [0, 0, 0, 1],
  ]);
  expect(plan.overlay.drawing_camera_placement_matrix_mm[2][3]).toBe(1700);
  expect(side.overlay.drawing_camera_placement_matrix_mm).toEqual([
    [1, 0, 0, 1820.961952209],
    [0, 0, -1, 626.980006695],
    [0, 1, 0, 1250],
    [0, 0, 0, 1],
  ]);
  for (const view of context.views) {
    expect(sha256(join(root, view.output))).toBe(view.output_sha256);
    expect(sha256(join(root, view.review_crop))).toBe(view.review_crop_sha256);
    expect(sha256(join(root, view.review_preview))).toBe(view.review_preview_sha256);
    expect(view.overlay.source_path_semantic_sha256).toBe(officialPathSemanticHashes[view.view]);
    const full = readFileSync(join(root, view.output), "utf8");
    const crop = readFileSync(join(root, view.review_crop), "utf8");
    expect(full).toContain('class="official-reference native-dwg project-context-review-overlay"');
    expect(full).toContain('data-review-only="true"');
    expect(full).toContain('data-bonsai-scene-svg="false"');
    expect(full).toContain('data-ifc-write-performed="false"');
    expect(full).toContain('data-uniform-geometry-scale="1.0"');
    expect(full).toContain('data-context-role="actual-ifc-body"');
    expect(full).toContain(`#miamisoft-e09-${view.view}-official-overlay path { fill:none !important; stroke:#1677c8 !important;`);
    expect((full.match(/class="official-reference native-dwg"/g) ?? [])).toHaveLength(view.overlay.source_path_count);
    expect(full).not.toContain("geometry-derived-proxy-mask");
    expect(crop).toContain('class="review-paper-background"');
  }
});

test("Miami Soft Bonsai evidence uses four cameras on the actual MODEL_VIEW IFC Body", () => {
  const evidence = JSON.parse(readFileSync(join(product, "bonsai-review-manifest.json"), "utf8"));
  expect(evidence).toMatchObject({
    mode: "actual_bonsai_ifc_body_camera_render",
    representative_global_id: "3osoWufdD1mhDDAM6lcix4",
    geometry_product_count: 1,
    whole_model_render: false,
    formal_ifc_bytes_unchanged: true,
    pass: true,
  });
  expect(evidence.bonsai_session).toMatchObject({ saved_active_representation: "Body", ifc_context_identifier: "Body", ifc_target_view: "MODEL_VIEW", saved_camera_count: 4, front_camera_local_y_sign: -1 });
  expect(evidence.renders.map((render: any) => render.view)).toEqual(["plan", "front", "side", "iso"]);
  for (const render of evidence.renders) {
    expect(render.camera_type).toBe("ORTHO");
    expect(render.render_operation).toBe("bpy.ops.render.render(write_still=True)");
    expect(sha256(join(root, render.path))).toBe(render.sha256);
  }
  expect(sha256(join(root, evidence.bonsai_session.path))).toBe(evidence.bonsai_session.sha256);
  expect(sha256(join(root, evidence.isolated_ifc))).toBe(evidence.isolated_ifc_sha256);
});

test("approved Miami Soft review opens only the product-level derived IFC gate", () => {
  const candidate = JSON.parse(readFileSync(join(product, "candidate-representations.json"), "utf8"));
  const approval = JSON.parse(readFileSync(join(root, "pipeline/decisions/miamisoft-e09-drawing-approval.json"), "utf8"));
  expect(candidate.formal_ifc_write_allowed).toBeFalse();
  expect(candidate.derived_ifc_write_allowed).toBeTrue();
  expect(candidate.review_status).toBe("approved_for_product_level_derived_ifc");
  expect(approval).toMatchObject({
    status: "approved",
    approved_views: ["plan", "front", "side"],
    derived_ifc_write_allowed: true,
    formal_authoritative_ifc_write_allowed: false,
    approval_evidence: { message: "批准通过，让 子代理 完成。" },
    write_execution: {
      component_force_fit_applied: false,
      side_hidden_line_source_handle: "1115B",
      side_hidden_line_pattern: "DASHED",
      formal_authoritative_ifc_written: false,
      scene_drawings_pending_visual_acceptance: true,
    },
  });
  const derived = join(product, "Baxter-Miami-Soft-E09-derived-drawing.ifc");
  expect(existsSync(derived)).toBeTrue();
  expect(sha256(derived)).toBe(approval.write_execution.derived_ifc_sha256);
  expect(existsSync(join(product, "Baxter-Miami-Soft-E09-bonsai-review.blend1"))).toBeFalse();
  expect(sha256(formal)).toBe(formalHash);
});

test("Bonsai Create Drawing persists three scene SVGs without target ghosting", () => {
  const evidence = JSON.parse(readFileSync(join(product, "bonsai-drawings/living-area/MIAMISOFT-E09-create-drawing-evidence.json"), "utf8"));
  const productManifest = JSON.parse(readFileSync(join(product, "manifest.json"), "utf8"));
  const drawingManifestPath = join(product, "living-area-bonsai-drawing-manifest.json");
  expect(evidence).toMatchObject({
    pass: true,
    formal_ifc_sha256: formalHash,
    formal_ifc_bytes_unchanged: true,
    save_boundary: "one Miami Soft E09 product-level full-project derived IFC copy",
    tests: {
      all_create_drawing_finished: true,
      all_target_body_projection_counts_zero: true,
      all_annotation_groups_unique: true,
      side_hidden_path_count: 1,
      side_hidden_dasharray: "2.4,1.5",
      formal_ifc_unchanged: true,
    },
  });
  expect(evidence.views.map((view: any) => view.view)).toEqual(["plan", "front", "side"]);
  for (const view of evidence.views) {
    expect(view.create_drawing).toMatchObject({ operator: "bpy.ops.bim.create_drawing", result: ["FINISHED"], linework_mode: "OPENCASCADE" });
    expect(view.component_force_fit_applied).toBeFalse();
    expect(view.svg.target_ifc_projection_group_count).toBe(0);
    expect(view.svg.duplicate_target_or_annotation).toBeFalse();
    expect(sha256(view.svg.path)).toBe(view.svg.sha256);
  }
  const side = evidence.views.find((view: any) => view.view === "side");
  const sideSvg = readFileSync(side.svg.path, "utf8");
  expect(side.persisted_hidden_path_count).toBe(1);
  expect((sideSvg.match(/data-line-role="hidden"/g) ?? [])).toHaveLength(259);
  expect((sideSvg.match(/stroke-dasharray:2.4,1.5/g) ?? [])).toHaveLength(259);
  expect(productManifest.derived_drawing_workflow).toMatchObject({
    status: "generated_for_scene_review",
    actual_target_body_suppressed_from_drawing_include: true,
    component_force_fit_applied: false,
    side_hidden_line_source_handle: "1115B",
    side_hidden_line_pattern: "DASHED",
    formal_authoritative_ifc_written: false,
    pass: true,
  });
  expect(sha256(drawingManifestPath)).toBe(productManifest.derived_drawing_workflow.manifest_sha256);
  expect(sha256(join(product, "Baxter-Miami-Soft-E09-living-area-drawings.pdf"))).toBe(productManifest.derived_drawing_workflow.combined_pdf_sha256);
  expect(sha256(formal)).toBe(formalHash);
});
