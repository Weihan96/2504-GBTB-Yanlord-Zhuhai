import { expect, test } from "bun:test";
import { createHash } from "node:crypto";
import { existsSync, mkdtempSync, readFileSync, rmSync, writeFileSync } from "node:fs";
import { tmpdir } from "node:os";
import { join, resolve } from "node:path";

const root = resolve(import.meta.dir, "../..");
const product = join(root, "output/review/highpoly-types/trap01");
const formal = join(root, "2504 GBTB Yanlord Zhuhai.ifc");
const formalHash = "7a50b87e8f48a7c2bfbaa9f1dabdfca7155fa684b2325c66aa1ae15c8ab0a25c";
const scope = "exact Geberit 151.116.11.1 adjustable family reference; project instance is a shortened installation configuration; not a project shop drawing";
const nativeHashes: Record<string, string> = {
  A: "01fa50da7cf39c6a6ab2a89339d1e4523270cb0943e8dce148f3d471f1a4fff2",
  G: "af55992e60aa1d84c9834da995f6ef8269593f43f46a6311cdf3e2e59530dd3f",
  L: "761eec37c809fc4319e9beaba6eaf472e19eefcbd16adf2815a4644d0a36fdc8",
  P: "c99cee2a0c9ae94b3a5535db4973982795eba64cedc5da135d5ca8f12264a8e1",
};

function sha256(path: string) {
  return createHash("sha256").update(readFileSync(path)).digest("hex");
}

test("exact Geberit 151.116.11.1 identity and official evidence are archived", () => {
  const access = JSON.parse(readFileSync(join(product, "official-source/source-access-record.json"), "utf8"));
  expect(access).toMatchObject({ manufacturer: "Geberit", resolved_article: "151.116.11.1", project_ifc_type_name: "TRAP01", project_ifc_type_description: "Space Saving Dip Tube Trap", scope, pass: true });
  expect(access.article_resolution).toMatchObject({ project_ifc_body_local_xyz_mm: [252.502579, 76.540974, 191.93998], selected_151_116_fixed_body_width_mm: 76.541707, selected_fixed_width_absolute_delta_mm: 0.000733, excluded_alternative: "151.117.11.1 / d40", selected_article_pass: true });
  expect(access.configuration_cross_check).toMatchObject({ project_instance_is_shortened_configuration: true, official_default_family_paths_used_as_project_representation: false, geometry_scaled_or_stretched_to_match: false, pass: true });
  const sourceHashes: Record<string, string> = {
    manufacturer_product_page: "566b4286b131cba21435fbb94e22ac46845886a6cc354750bd9340a2e0cd7864",
    manufacturer_product_data_sheet: "e331a53269493b86908611c1ab2f5d05c401aff2530c464546342a1dece0d06b",
    manufacturer_installation_instructions: "edbf40f45788262b056f4d1d366731c0242cd523b70a97215f00091a0ebe5e9a",
    manufacturer_maintenance_manual: "1bbafe09d372ba2c52058fb1a795c98086d85f6fb30625b36670e23dd7da6edf",
  };
  for (const source of access.official_sources) {
    expect(source.sha256).toBe(sourceHashes[source.kind]);
    expect(sha256(join(root, source.local_path))).toBe(source.sha256);
    if (source.preview_path) expect(sha256(join(root, source.preview_path))).toBe(source.preview_sha256);
  }
  expect(access.official_sources.find((source: any) => source.kind === "manufacturer_product_data_sheet").vector_object_audit).toMatchObject({ drawing_count: 55, drawing_item_count: 391, word_count: 142, character_count: 785 });
  expect(access.official_sources.find((source: any) => source.kind === "manufacturer_installation_instructions").page_2_text_audit.required_text_found).toEqual(["113-365", "0-252", "32", "40", "85-334"]);
});

test("direct native DWG parsing locks article width without force-fitting the adjustable default", () => {
  const linework = JSON.parse(readFileSync(join(product, "official-native-dwg-linework.json"), "utf8"));
  const access = JSON.parse(readFileSync(join(product, "official-source/source-access-record.json"), "utf8"));
  expect(linework).toMatchObject({ article_number: "151.116.11.1", source_kind: "native_dwg", pass: true });
  expect(Object.fromEntries(Object.entries(linework.views).map(([view, value]: any) => [view, value.paths_mm.length]))).toEqual({ plan: 92, front: 77, side: 117 });
  expect(access.native_cad_selection).toMatchObject({ article_number: "151.116.11.1", source_kind: "native_dwg", sha256: nativeHashes, path_counts: { plan: 92, front: 77, side: 117 }, official_cad_acquired: true, third_party_cad_used: false, pass: true });
  for (const [code, hash] of Object.entries(nativeHashes)) expect(sha256(join(product, `official-source/151.116.11.1_${code}.dwg`))).toBe(hash);
  expect(linework.configuration_cross_check).toMatchObject({ project_instance_is_shortened_configuration: true, official_default_family_paths_used_as_project_representation: false, fixed_diameter_match_pass: true, pass: true });
});

test("candidate aligns black configuration with solid and dashed official family ranges", () => {
  const manifest = JSON.parse(readFileSync(join(product, "manifest.json"), "utf8"));
  const candidate = JSON.parse(readFileSync(join(product, "candidate-representations.json"), "utf8"));
  const inventory = JSON.parse(readFileSync(join(root, "pipeline/decisions/highpoly-product-inventory.json"), "utf8"));
  expect(inventory.review_ready_pending_approval_product_count).toBe(38);
  expect(inventory.products.find((entry: any) => entry.type_name === "TRAP01")).toMatchObject({ status: "review_ready_pending_approval", representative_global_id: "2Ak2ma0lvBEA49UpplzUqi" });
  expect(manifest).toMatchObject({ representative_global_id: "2Ak2ma0lvBEA49UpplzUqi", geometry_product_count: 1, whole_model_render: false, review_status: "visual_review_pending", approved_for_drawing_ifc: false, pass: true });
  expect(manifest).toMatchObject({ source_kind: "geometry_derived_simplified_proxy", source_label_zh: "基于原始高模几何生成的简化图纸表达", official_cad_acquired: true, official_cad_used: false, official_cad_used_as_representation: false, third_party_cad_used: false, blue_line_present: true });
  expect(manifest.drawing_source).toMatchObject({ source_kind: "geometry_derived_simplified_proxy", source_label_zh: "基于原始高模几何生成的简化图纸表达", official_cad_acquired: true, official_cad_used_as_representation: false, third_party_cad_used: false });
  expect(manifest.configuration_note).toContain("shortened project installation configuration");
  expect(candidate).toMatchObject({ source_kind: "geometry_derived_simplified_proxy", source_label_zh: "基于原始高模几何生成的简化图纸表达", official_cad_acquired: true, official_cad_used: false, official_cad_used_as_representation: false, third_party_cad_used: false, blue_line_present: true, formal_ifc_write_allowed: false, review_status: "visual_review_pending" });
  for (const [view, counts] of Object.entries({ plan: [16, 92], front: [1, 77], side: [8, 117] })) {
    expect(candidate.views[view].proxy_paths_mm.length).toBe(counts[0]);
    expect(candidate.views[view].official_cad_paths_mm).toEqual([]);
    expect(candidate.views[view].official_family_native_dwg_paths_mm.length).toBe(counts[1]);
    expect(candidate.views[view].adjustable_overlay.retained_solid_fragment_count).toBeGreaterThan(0);
    expect(candidate.views[view].adjustable_overlay.clipped_dashed_fragment_count).toBeGreaterThan(0);
    expect(candidate.views[view].adjustable_overlay).toMatchObject({ alignment_method: "translation_only_from_fixed_body_and_connection_axis_datums", scale: 1, reflection: false, horizontal_and_vertical_classified_independently: true, geometry_scaled_or_stretched: false, pass: true });
    expect(candidate.views[view].official_family_paths_used_as_project_representation).toBeFalse();
    const svg = readFileSync(join(product, `${view}.svg`), "utf8");
    expect(svg).toContain('class="configured-simplified-proxy geometry-derived"');
    expect(svg).toContain('class="official-family-native-dwg retained-configured-range"');
    expect(svg).toContain('class="official-family-native-dwg clipped-adjustable-extension"');
    expect(svg).toContain('stroke-dasharray="12 8"');
    expect(svg).toContain('data-installation-status="not-project-actual"');
    expect(svg).toContain('data-used-as-project-representation="false"');
  }
  const audit = JSON.parse(readFileSync(join(product, "adjustable-overlay-audit.json"), "utf8"));
  expect(audit).toMatchObject({ official_adjustment_ranges_mm: { horizontal_extension: [0, 252], h: [85, 334] }, geometry_scaled_or_stretched: false, formal_ifc_write: "not performed", pass: true });
  expect(audit.views.plan.clipped_extension_mm).toMatchObject({ horizontal: 159.52933, vertical: 0 });
  expect(audit.views.front.clipped_extension_mm.horizontal).toBeGreaterThan(180);
  expect(audit.views.front.clipped_extension_mm.vertical).toBeGreaterThan(200);
  expect(audit.views.side.clipped_extension_mm).toMatchObject({ horizontal: 0, vertical: 208.00002 });
});

test("project sanitary plan retains context and mechanically restores the occluded configured shape", () => {
  const context = JSON.parse(readFileSync(join(product, "project-context-manifest.json"), "utf8"));
  expect(context).toMatchObject({ source_kind: "geometry_derived_simplified_proxy", source_label_zh: "基于原始高模几何生成的简化图纸表达", official_cad_acquired: true, official_cad_used: false, official_cad_used_as_representation: false, blue_line_present: false, project_context_retained: true, walls_and_surrounding_project_elements_retained: true, overlay_top_layer_with_white_mask: true, blue_product_cad_line_present: false, pass: true });
  expect(context.context_view_scope).toMatchObject({ included: ["plan"], excluded: ["front", "side"] });
  expect(context.views).toHaveLength(1);
  expect(context.views[0].overlay).toMatchObject({ path_count: 16, fixed_project_scale_preserved: true, geometry_scaled_or_stretched: false });
  expect(context.views[0].overlay.mechanical_translation).toMatchObject({ project_visible_segment_count: 613, matching_edge_vote_count: 853, translation_svg_units: [71.64, 210.999997], scale_svg_units_per_mm: 0.02, original_projection_is_occluded: true, pass: true });
  expect(sha256(join(root, context.views[0].review_preview))).toBe(context.views[0].review_preview_sha256);
  const full = readFileSync(join(root, context.views[0].output), "utf8");
  expect(full).toContain('class="configured-proxy-mask"');
  expect(full).toContain('class="configured-proxy"');
  expect(full.indexOf('class="configured-proxy-mask"')).toBeLessThan(full.indexOf('class="configured-proxy"'));
});

test("actual Bonsai evidence saves four cameras and renders the isolated MODEL_VIEW Body", () => {
  const evidence = JSON.parse(readFileSync(join(product, "bonsai-review-manifest.json"), "utf8"));
  expect(evidence).toMatchObject({ mode: "actual_bonsai_ifc_body_camera_render", representative_global_id: "2Ak2ma0lvBEA49UpplzUqi", geometry_product_count: 1, whole_model_render: false, formal_ifc_bytes_unchanged: true, pass: true });
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

test("approved TRAP01 opens only its componentised derived IFC write gate", () => {
  const manifest = join(product, "manifest.json");
  const approval = JSON.parse(readFileSync(join(root, "pipeline/decisions/trap01-drawing-approval.json"), "utf8"));
  expect(approval).toMatchObject({ status: "approved", reviewer: "project_owner", approved_views: ["plan", "front", "side"], approved_component_scope: ["fixed_body", "horizontal_adjustable", "vertical_adjustable"], candidate_manifest_sha256: sha256(manifest), derived_ifc_write_allowed: true, formal_authoritative_ifc_write_allowed: false });
  expect(approval.approved_line_semantics).toMatchObject({ black: "current installed project configuration", blue_solid: "official native DWG linework retained within the current configured range", blue_dashed: "official adjustable-family excess beyond the current configured endpoint; reference only, not installed project geometry" });
  expect(existsSync(join(product, "Geberit-151.116.11.1-TRAP01-derived-drawing.ifc"))).toBeTrue();
  expect(existsSync(join(product, "Geberit-151.116.11.1-TRAP01-bonsai-review.blend1"))).toBeFalse();
  expect(sha256(formal)).toBe(formalHash);
  const temporary = mkdtempSync(join(tmpdir(), "trap01-gate-"));
  const output = join(temporary, "forbidden.ifc");
  const script = join(root, "pipeline/scripts/trap01_drawing_ifc.py");
  const withoutApply = Bun.spawnSync(["python3", script, "--input", formal, "--output", output], { cwd: root, stderr: "pipe" });
  expect(withoutApply.exitCode).not.toBe(0);
  expect(withoutApply.stderr.toString()).toContain("IFC write requires the explicit --apply flag");
  const approved = Bun.spawnSync(["python3", script, "--input", formal, "--output", output, "--apply"], { cwd: root, stderr: "pipe" });
  expect(approved.exitCode).toBe(0);
  expect(existsSync(output)).toBeTrue();
  expect(sha256(formal)).toBe(formalHash);
  rmSync(temporary, { recursive: true, force: true });
}, 30_000);

test("actual Bonsai scene Drawings show one current configuration and no adjustable excess", () => {
  const manifest = JSON.parse(readFileSync(join(product, "TRAP01-project-drawing-manifest.json"), "utf8"));
  const evidence = JSON.parse(readFileSync(join(product, "bonsai-drawings/project-context/TRAP01-create-drawing-evidence.json"), "utf8"));
  expect(manifest).toMatchObject({ mode: "actual_bonsai_create_drawing_project_context", formal_ifc_sha256: formalHash, formal_ifc_bytes_unchanged: true, component_model_pass: true, current_installation_only: true, blue_dashed_reference_displayed: false, no_duplicate_target_or_annotation: true, pass: true });
  expect(manifest.derived_ifc.sha256).toBe(sha256(join(root, manifest.derived_ifc.path)));
  expect(evidence).toMatchObject({ verdict: "pass", formal_ifc_bytes_unchanged: true, tests: { all_create_drawing_finished: true, all_target_body_projection_counts_zero: true, all_current_annotation_counts_one: true, all_component_reference_counts_zero: true, blue_dashed_reference_displayed: false }, pass: true });
  expect(evidence.postState).toMatchObject({ component_count: 3, drawing_count: 3, current_linework_annotation_count: 3 });
  expect(manifest.views.map((entry: any) => entry.view)).toEqual(["plan", "front", "side"]);
  for (const view of manifest.views) {
    expect(view.create_drawing).toMatchObject({ operator: "bpy.ops.bim.create_drawing", result: ["FINISHED"], linework_mode: "OPENCASCADE" });
    expect(view.svg_validation).toMatchObject({ target_ifc_body_group_count: 0, current_annotation_group_count: 1, adjustable_component_reference_group_count: 0, no_target_or_annotation_duplicate: true, blue_dashed_reference_displayed: false });
    expect(view.visual).toMatchObject({ black_present: true, grey_context_present: true, blue_dashed_reference_absent: true, pass: true });
    expect(sha256(join(root, view.svg.path))).toBe(view.svg.sha256);
    expect(sha256(join(root, view.page_pdf.path))).toBe(view.page_pdf.sha256);
    expect(sha256(join(root, view.preview.path))).toBe(view.preview.sha256);
  }
  expect(sha256(formal)).toBe(formalHash);
});

test("temporary scoped approval writes configured proxies plus official source associations", () => {
  const manifest = join(product, "manifest.json");
  const temporary = mkdtempSync(join(tmpdir(), "trap01-approved-"));
  const approvalPath = join(temporary, "approval.json");
  const output = join(temporary, "derived.ifc");
  const report = join(temporary, "report.json");
  writeFileSync(approvalPath, JSON.stringify({ schema_version: 1, profile_key: "trap01", ifc_type_name: "TRAP01", candidate_manifest_sha256: sha256(manifest), status: "approved", reviewer: "automated gate fixture", review_date: "2026-08-23", approved_views: ["plan", "front", "side"], derived_ifc_write_allowed: true, formal_authoritative_ifc_write_allowed: false, scope, approval_evidence: "temporary automated writer verification only" }, null, 2));
  const run = Bun.spawnSync(["python3", join(root, "pipeline/scripts/trap01_drawing_ifc.py"), "--input", formal, "--manifest", manifest, "--approval", approvalPath, "--output", output, "--report", report, "--apply"], { cwd: root, stdout: "pipe", stderr: "pipe" });
  expect(run.exitCode).toBe(0);
  if (run.exitCode !== 0) throw new Error(run.stderr.toString());
  const result = JSON.parse(readFileSync(report, "utf8"));
  expect(result).toMatchObject({ pass: true, formal_ifc_bytes_unchanged: true, representations: { plan: "Trap01Plan", front: "Trap01Front", side: "Trap01Side" }, representation_path_counts: { plan: 16, front: 1, side: 8 }, representation_geometry_source: "proxy_paths_mm derived from the configured isolated representative IFC Body", official_cad_geometry_included: false, source_property_set: "Pset_Trap01DrawingSource", source_kind: "geometry_derived_simplified_proxy", source_label_zh: "基于原始高模几何生成的简化图纸表达", component_model_pass: true });
  expect(Object.keys(result.adjustable_components).sort()).toEqual(["fixed_body", "horizontal_adjustable", "vertical_adjustable"]);
  expect(result.adjustment_parameters_mm).toMatchObject({ horizontal: { installed: 113.402666, minimum: 0, maximum: 252, official_default_overflow: 159.52933 }, vertical_h: { installed: 191.939979, minimum: 85, maximum: 334, official_default_overflow: 208.00002 } });
  expect(result).toMatchObject({ blue_solid_semantics: "official native DWG linework retained within the current configured range", blue_dashed_semantics: "official adjustable-family excess beyond the current configured endpoint; reference only, not installed project geometry", reference_extension_displayed_in_scene_by_default: false });
  expect(result.source_document_associations).toContain("GEBERIT-151-116-11-1-TRAP01-OFFICIAL-G-DWG");
  expect(result.source_document_associations).toContain("GEBERIT-151-116-11-1-TRAP01-OFFICIAL-PRODUCT-PAGE");
  expect(result.source_document_associations).toContain("GEBERIT-151-116-11-1-TRAP01-NATIVE-DWG-LINEWORK");
  expect(result.source_document_associations).toHaveLength(11);
  expect(sha256(formal)).toBe(formalHash);
  rmSync(temporary, { recursive: true, force: true });
}, 30_000);
