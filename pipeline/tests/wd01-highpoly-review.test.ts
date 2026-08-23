import { expect, test } from "bun:test";
import { createHash } from "node:crypto";
import { existsSync, mkdtempSync, readFileSync, rmSync, writeFileSync } from "node:fs";
import { tmpdir } from "node:os";
import { join, resolve } from "node:path";

const root = resolve(import.meta.dir, "../..");
const product = join(root, "output/review/highpoly-types/wd01");
const formal = join(root, "2504 GBTB Yanlord Zhuhai.ifc");
const formalHash = "7a50b87e8f48a7c2bfbaa9f1dabdfca7155fa684b2325c66aa1ae15c8ab0a25c";
const sourceKind = "geometry_derived_simplified_proxy";
const sourceLabelZh = "基于原始高模几何生成的简化图纸表达";
const scope = "official manufacturer Pivot and Senzafine system identity and remote technical-document evidence only; the project 1202.879 x 673.524 x 2390.023 mm arrangement is a custom configuration, not exact official CAD geometry or a project shop drawing";

function sha256(path: string) {
  return createHash("sha256").update(readFileSync(path)).digest("hex");
}

test("official Poliform sources establish the Pivot and Senzafine system but not exact WD01 CAD", () => {
  const access = JSON.parse(readFileSync(join(product, "official-source/source-access-record.json"), "utf8"));
  const evidencePath = join(product, "official-source/official-product-page-evidence.json");
  const evidence = JSON.parse(readFileSync(evidencePath, "utf8"));
  expect(access).toMatchObject({
    manufacturer: "Poliform",
    family: "Pivot + Senzafine Wardrobe",
    project_type_code: "WD01",
    ifc_type_description: "Poliform Pivot Senzafine",
    scope,
    pass: true,
  });
  expect(access.official_product_identity).toMatchObject({
    pivot_product_name: "Pivot",
    senzafine_product_name: "Senzafine Wardrobe",
    integrated_system_relationship_supported: true,
    customised_compositions_supported: true,
    project_ifc_type_description_exactly_names_manufacturer_and_both_systems: true,
  });
  expect(access.official_product_cad).toMatchObject({
    registration_form_and_captcha_required: true,
    public_exact_native_dwg_url_located: false,
    acquired: false,
    local_cad_files: [],
    exact_project_configuration_match: false,
  });
  expect(access.official_remote_documents).toMatchObject({
    identity_and_system_relationship_observed: true,
    binary_archived: false,
    binary_archive_blocker: "HTTP 403 Cloudflare challenge on direct retrieval",
    used_as_representation_geometry: false,
  });
  expect(access.dimension_cross_check).toMatchObject({
    project_ifc_body_bounds_mm: [1202.878845, 673.52434, 2390.022879],
    exact_dimension_match_claimed: false,
    geometry_scaled_to_catalogue_example: false,
  });
  expect(access.official_identity_sources).toHaveLength(1);
  expect(access.official_identity_sources[0].sha256).toBe(sha256(evidencePath));
  expect(evidence.official_page_observations.pivot_technical_pdf).toMatchObject({
    exact_project_1202_879_x_673_524_x_2390_023_mm_configuration_identified: false,
  });
  expect(evidence.drawing_geometry_conclusion).toMatchObject({
    official_native_dwg_acquired: false,
    official_native_dwg_used: false,
    third_party_cad_used: false,
    source_kind: sourceKind,
  });
});

test("one WD01 Body produces three geometry-derived views with zero blue CAD paths", () => {
  const manifest = JSON.parse(readFileSync(join(product, "manifest.json"), "utf8"));
  const candidate = JSON.parse(readFileSync(join(product, "candidate-representations.json"), "utf8"));
  const profile = JSON.parse(readFileSync(join(product, "profile.json"), "utf8"));
  const register = JSON.parse(readFileSync(join(root, "pipeline/decisions/highpoly-drawing-profile-register.json"), "utf8"));
  const inventory = JSON.parse(readFileSync(join(root, "pipeline/decisions/highpoly-product-inventory.json"), "utf8"));
  expect(manifest.representative_global_id).toBe("2mPTt7$nvBSAn28I9NWt6T");
  expect(manifest.registered_instance_global_ids).toEqual(["2mPTt7$nvBSAn28I9NWt6T"]);
  expect(manifest.geometry_product_count).toBe(1);
  expect(manifest.whole_model_render).toBeFalse();
  expect(manifest.formal_ifc_bytes_unchanged).toBeTrue();
  expect(manifest.ifc_type_description).toBe("Poliform Pivot Senzafine");
  expect(manifest.bounds_mm.size).toEqual([1202.878845, 673.52434, 2390.022879]);
  expect(profile.profiles.wd01.representative_global_id).toBe("2mPTt7$nvBSAn28I9NWt6T");
  expect(register.profiles.wd01.representative_global_id).toBe("2mPTt7$nvBSAn28I9NWt6T");
  expect(inventory.review_ready_pending_approval_product_count).toBe(38);
  expect(inventory.products.find((item: any) => item.type_name === "WD01")).toMatchObject({
    folder: "output/review/highpoly-types/wd01",
    status: "review_ready_pending_approval",
    source_strategy: "poliform_Pivot_Senzafine_official_system_identity_and_remote_technical_document_verified_gated_resources_exact_project_native_cad_not_acquired_geometry_derived_proxy",
  });
  expect(candidate.source_kind).toBe(sourceKind);
  expect(candidate.source_label_zh).toBe(sourceLabelZh);
  expect(candidate.official_cad_used).toBeFalse();
  expect(candidate.third_party_cad_used).toBeFalse();
  expect(manifest.views.map((view: any) => view.silhouette_path_count)).toEqual([7, 3, 7]);
  for (const view of manifest.views) {
    expect(view.official_cad_path_count).toBe(0);
    expect(view.blue_line_present).toBeFalse();
    const svg = readFileSync(join(root, view.svg), "utf8");
    expect(svg).toContain('class="simplified-proxy-silhouette geometry-derived"');
    expect(svg).toContain('data-source-kind="geometry_derived_simplified_proxy"');
    expect(svg).not.toContain('class="official-reference native-dwg"');
    expect(svg).not.toContain("#1677c8");
  }
});

test("WD01 project context keeps real walls and furniture in plan and both elevations", () => {
  const context = JSON.parse(readFileSync(join(product, "project-context-manifest.json"), "utf8"));
  expect(context).toMatchObject({
    source_kind: sourceKind,
    source_label_zh: sourceLabelZh,
    official_cad_used: false,
    third_party_cad_used: false,
    project_context_retained: true,
    walls_and_surrounding_project_elements_retained: true,
    overlay_top_layer_with_white_mask: true,
    blue_line_present: false,
    pass: true,
  });
  expect(context.context_view_scope).toMatchObject({ included: ["plan", "front", "side"], excluded: [] });
  expect(context.context_view_scope.reason).toContain("alternate R09 -Y elevation");
  expect(context.views.map((view: any) => view.view)).toEqual(["plan", "front", "side"]);
  expect(context.semantic_view_mapping).toEqual({
    plan: { candidate_axes: [0, 1], source: "drawings/Furniture Plan.svg", rotate_quarter_turns: 1, project_direction: "+Z" },
    front: { candidate_axes: [0, 2], source: "drawings/elevations/native/EL-05-15-R09-PX.svg", rotate_quarter_turns: 0, project_direction: "+X" },
    side: { candidate_axes: [1, 2], source: "drawings/elevations/native/EL-05-14-R09-PY.svg", rotate_quarter_turns: 0, project_direction: "+Y" },
  });
  expect(context.project_projection_note).toMatchObject({
    side_original_visible_projection_bbox_mm: [659.48965, 2389.9528],
    complete_ifc_body_side_bbox_mm: [673.52434, 2390.022879],
    geometry_stretched: false,
  });
  for (const view of context.views) {
    expect(view.overlay.fit.uniform_scale_preserved).toBeTrue();
    expect(view.overlay.fit.scale_svg_units_per_mm).toBe(0.02);
    expect(view.overlay.fit.transformation_mode).toBe("axis_swap_rigid_reflection_and_translation_only");
    expect(Math.max(...view.overlay.fit.bbox_absolute_delta_svg_units)).toBeLessThanOrEqual(view.overlay.fit.bbox_tolerance_svg_units);
    expect(sha256(join(root, view.review_preview))).toBe(view.review_preview_sha256);
    const svg = readFileSync(join(root, view.output), "utf8");
    expect(svg).toContain("IfcFurniture");
    expect(svg).toContain('class="geometry-derived-simplified-proxy project-context-overlay"');
    expect(svg).toContain('class="geometry-derived-proxy-mask"');
    expect(svg.indexOf('class="geometry-derived-proxy-mask"')).toBeLessThan(svg.indexOf('class="geometry-derived-proxy"'));
  }
});

test("WD01 Bonsai evidence uses four cameras on the actual IFC Body", () => {
  const evidence = JSON.parse(readFileSync(join(product, "bonsai-review-manifest.json"), "utf8"));
  expect(evidence.mode).toBe("actual_bonsai_ifc_body_camera_render");
  expect(evidence.geometry_product_count).toBe(1);
  expect(evidence.whole_model_render).toBeFalse();
  expect(evidence.formal_ifc_bytes_unchanged).toBeTrue();
  expect(evidence.bonsai_session.saved_active_representation).toBe("Body");
  expect(evidence.bonsai_session.ifc_context_identifier).toBe("Body");
  expect(evidence.bonsai_session.saved_camera_count).toBe(4);
  expect(evidence.bonsai_session.local_axis_extents_m).toMatchObject({ x: 1.202879003, y: 0.673524947, z: 2.390022993 });
  expect(evidence.renders.map((render: any) => render.view)).toEqual(["plan", "front", "side", "iso"]);
  for (const render of evidence.renders) {
    expect(render.camera_type).toBe("ORTHO");
    expect(render.render_operation).toBe("bpy.ops.render.render(write_still=True)");
    expect(sha256(join(root, render.path))).toBe(render.sha256);
  }
  expect(sha256(join(root, evidence.bonsai_session.path))).toBe(evidence.bonsai_session.sha256);
  expect(sha256(join(root, evidence.isolated_ifc))).toBe(evidence.isolated_ifc_sha256);
});

test("pending WD01 approval blocks every IFC write and preserves the formal IFC", () => {
  const candidate = JSON.parse(readFileSync(join(product, "candidate-representations.json"), "utf8"));
  const approval = JSON.parse(readFileSync(join(root, "pipeline/decisions/wd01-drawing-approval.json"), "utf8"));
  expect(candidate.formal_ifc_write_allowed).toBeFalse();
  expect(candidate.review_status).toBe("visual_review_pending");
  expect(approval.status).toBe("pending");
  expect(approval.derived_ifc_write_allowed).toBeFalse();
  expect(approval.formal_authoritative_ifc_write_allowed).toBeFalse();
  expect(approval.candidate_manifest_sha256).toBe(sha256(join(product, "manifest.json")));
  expect(existsSync(join(product, "Poliform-Pivot-Senzafine-WD01-derived-drawing.ifc"))).toBeFalse();
  expect(sha256(formal)).toBe(formalHash);
});

test("WD01 writer rejects missing apply and the pending human approval record", () => {
  const temporary = mkdtempSync(join(tmpdir(), "wd01-gate-"));
  const output = join(temporary, "forbidden.ifc");
  const script = join(root, "pipeline/scripts/wd01_drawing_ifc.py");
  const withoutApply = Bun.spawnSync(["python3", script, "--input", formal, "--output", output], { cwd: root, stderr: "pipe" });
  expect(withoutApply.exitCode).not.toBe(0);
  expect(withoutApply.stderr.toString()).toContain("IFC write requires the explicit --apply flag");
  expect(existsSync(output)).toBeFalse();
  const pendingApproval = Bun.spawnSync(["python3", script, "--input", formal, "--output", output, "--apply"], { cwd: root, stderr: "pipe" });
  expect(pendingApproval.exitCode).not.toBe(0);
  expect(pendingApproval.stderr.toString()).toContain("approval gate rejected IFC write");
  expect(existsSync(output)).toBeFalse();
  expect(sha256(formal)).toBe(formalHash);
  rmSync(temporary, { recursive: true, force: true });
});

test("a scoped temporary approval writes verified WD01 representations and source documents", () => {
  const manifest = join(product, "manifest.json");
  const temporary = mkdtempSync(join(tmpdir(), "wd01-approved-"));
  const approvalPath = join(temporary, "approval.json");
  const output = join(temporary, "derived.ifc");
  const report = join(temporary, "report.json");
  writeFileSync(approvalPath, JSON.stringify({
    schema_version: 1,
    profile_key: "wd01",
    ifc_type_name: "WD01",
    candidate_manifest_sha256: sha256(manifest),
    status: "approved",
    reviewer: "automated gate fixture",
    review_date: "2026-08-23",
    approved_views: ["plan", "front", "side"],
    derived_ifc_write_allowed: true,
    formal_authoritative_ifc_write_allowed: false,
    scope,
    approval_evidence: "temporary automated writer verification only",
  }, null, 2));
  const run = Bun.spawnSync([
    "python3", join(root, "pipeline/scripts/wd01_drawing_ifc.py"),
    "--input", formal, "--manifest", manifest, "--approval", approvalPath,
    "--output", output, "--report", report, "--apply",
  ], { cwd: root, stdout: "pipe", stderr: "pipe" });
  expect(run.exitCode).toBe(0);
  if (run.exitCode !== 0) throw new Error(run.stderr.toString());
  const result = JSON.parse(readFileSync(report, "utf8"));
  expect(result.pass).toBeTrue();
  expect(result.formal_ifc_bytes_unchanged).toBeTrue();
  expect(result.representations).toEqual({ plan: "Wd01Plan", front: "Wd01Front", side: "Wd01Side" });
  expect(result.representation_path_counts).toEqual({ plan: 7, front: 3, side: 7 });
  expect(result.source_kind).toBe(sourceKind);
  expect(result.source_label_zh).toBe(sourceLabelZh);
  expect(result.official_cad_geometry_included).toBeFalse();
  expect(result.source_document_associations).toEqual([
    "POLIFORM-PIVOT-SENZAFINE-WD01-OFFICIAL-PAGE-EVIDENCE",
    "POLIFORM-PIVOT-SENZAFINE-WD01-OFFICIAL-PIVOT-PAGE",
    "POLIFORM-PIVOT-SENZAFINE-WD01-OFFICIAL-SENZAFINE-PAGE",
    "POLIFORM-PIVOT-SENZAFINE-WD01-OFFICIAL-TECHNICAL-PDF",
    "POLIFORM-PIVOT-SENZAFINE-WD01-SOURCE-ACCESS-RECORD",
  ]);
  expect(sha256(formal)).toBe(formalHash);
  rmSync(temporary, { recursive: true, force: true });
}, 30_000);
