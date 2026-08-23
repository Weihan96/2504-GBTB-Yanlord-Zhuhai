import { expect, test } from "bun:test";
import { createHash } from "node:crypto";
import { existsSync, mkdtempSync, readFileSync, rmSync, writeFileSync } from "node:fs";
import { tmpdir } from "node:os";
import { join, resolve } from "node:path";

const root = resolve(import.meta.dir, "../..");
const product = join(root, "output/review/highpoly-types/marilyn-01");
const formal = join(root, "2504 GBTB Yanlord Zhuhai.ifc");
const formalHash = "7a50b87e8f48a7c2bfbaa9f1dabdfca7155fa684b2325c66aa1ae15c8ab0a25c";
const dwgHash = "cb9825eecab7c78ee28e7668e933c2fe413395a63bf75f3afab8cc4959864724";
const scope = "exact Baxter Marilyn bergere 86 x 100 x 94 cm family CAD reference; not a project shop drawing";

function sha256(path: string) {
  return createHash("sha256").update(readFileSync(path)).digest("hex");
}

test("exact Baxter Marilyn bergere identity and every official source hash are archived", () => {
  const access = JSON.parse(readFileSync(join(product, "official-source/source-access-record.json"), "utf8"));
  expect(access).toMatchObject({ manufacturer: "Baxter", family: "Marilyn", designer: "Draga & Aurel", project_ifc_type_name: "Marilyn 01", project_ifc_type_description: "Bergère armchair with swivel base W86D100H78", resolved_variant: "Marilyn bergere armchair with swivel base - 86 x 100 x 94 cm", scope, pass: true });
  expect(access.project_identity_evidence).toMatchObject({ variant_status: "confirmed_marilyn_bergere_86x100x94", description_status: "height_conflicts_with_official_bergere_and_matches_standard_armchair_height", body_geometry_status: "matches_official_bergere_within_35mm" });
  expect(access.dimension_cross_check).toMatchObject({ official_nominal_width_depth_height_mm: [860, 1000, 940], project_ifc_body_local_xyz_mm: [877.174957, 1013.911865, 970.476962], maximum_absolute_delta_mm: 30.476962, tolerance_mm: 35, pass: true });
  expect(access.drawing_geometry_source.source_label_zh).toBe(
    "基于 Baxter 精确型号原生 DWG 的官方图纸表达",
  );
  expect(access.official_source_revalidation).toMatchObject({
    marilyn_02_geometry_used: false,
    pass: true,
  });
  expect(existsSync(join(root, access.official_source_revalidation.path))).toBeTrue();
  const expected: Record<string, string> = {
    manufacturer_product_page: "e9f630d77aff716b01d21736117df23a491984183c2e4ea0f00c87a61df85140",
    manufacturer_current_vector_technical_sheet: "8210ae5ce84b467fd6cb46a93c665fbfaed1606e806210a751ef24ea026ff581",
    manufacturer_measurement_svg_matte: "4c2b4e94f1814bab8ffd6a154481c332577669d4afed43781c6a34fa509cc7c4",
    manufacturer_measurement_svg_glossy: "d5e0b04cce3465e60948c5ec6b33748666221744eac83642e71ee5e0f71e204c",
    manufacturer_native_2d_3d_zip: "db054a193c9c5b8dc45f0df39199b2d96dc3712a5a7e48821fe94dce1566fedb",
    manufacturer_native_2d_dwg: dwgHash,
    manufacturer_exact_model_3ds: "171dbd72b2298f201d834fb54a5d39bca7ca5b86a552b97a5cbc64e23dd5e385",
  };
  for (const source of access.official_sources) {
    expect(source.sha256).toBe(expected[source.kind]);
    expect(sha256(join(root, source.local_path))).toBe(source.sha256);
    if (source.preview_path) expect(sha256(join(root, source.preview_path))).toBe(source.preview_sha256);
  }
  expect(access.official_sources.find((source: any) => source.kind === "manufacturer_current_vector_technical_sheet").vector_object_audit).toMatchObject({ drawing_count: 730, drawing_item_count: 6956, word_count: 81, character_count: 390 });
  expect(access.official_sources.filter((source: any) => source.kind.includes("measurement_svg")).map((source: any) => source.vector_object_audit)).toEqual([{ path_count: 247, use_count: 54, group_count: 50 }, { path_count: 247, use_count: 54, group_count: 50 }]);
});

test("live Baxter sources preserve the exact bergere package and exclude Marilyn 02 geometry", () => {
  const evidence = JSON.parse(
    readFileSync(join(product, "official-source/official-source-revalidation.json"), "utf8"),
  );
  expect(evidence).toMatchObject({
    manufacturer: "Baxter",
    family: "Marilyn",
    project_ifc_type_name: "Marilyn 01",
    resolved_variant: "bergere armchair with swivel base, 86 x 100 x 94 cm",
    marilyn_02_geometry_used: false,
    pass: true,
  });
  expect(evidence.product_page).toMatchObject({
    http_status: 200,
    all_required_identity_tokens_found: true,
    pass: true,
  });
  expect(Object.values(evidence.product_page.required_identity_tokens).every(Boolean)).toBeTrue();
  expect(evidence.native_zip).toMatchObject({
    http_status: 200,
    expected_sha256: "db054a193c9c5b8dc45f0df39199b2d96dc3712a5a7e48821fe94dce1566fedb",
    downloaded_bytes_match_local_archive: true,
    dwg_member_sha256: dwgHash,
    dwg_member_matches_local_extraction: true,
    bergere_3ds_member_sha256: "171dbd72b2298f201d834fb54a5d39bca7ca5b86a552b97a5cbc64e23dd5e385",
    bergere_3ds_member_matches_local_extraction: true,
    pass: true,
  });
  expect(evidence.technical_pdf).toMatchObject({
    http_status: 200,
    identity_page: 13,
    expected_identity_page_render_sha256: "977e90a62a151b8a68c4abbdc41ba4beecbefd989bd41a220770debb9fb3cfc5",
    current_identity_page_render_sha256: "977e90a62a151b8a68c4abbdc41ba4beecbefd989bd41a220770debb9fb3cfc5",
    archived_identity_page_render_sha256: "977e90a62a151b8a68c4abbdc41ba4beecbefd989bd41a220770debb9fb3cfc5",
    current_identity_page_render_matches_archive: true,
    pass: true,
  });
  expect(Object.values(evidence.measurement_svgs).every((item: any) => item.pass)).toBeTrue();
});

test("direct native-DWG parsing supplies exactly 24/68/54 approved blue paths", () => {
  const linework = JSON.parse(readFileSync(join(product, "official-native-dwg-linework.json"), "utf8"));
  expect(linework).toMatchObject({ source_kind: "native_dwg", source_dwg_sha256: dwgHash, source_zip_sha256: "db054a193c9c5b8dc45f0df39199b2d96dc3712a5a7e48821fe94dce1566fedb", exact_model_3ds_sha256: "171dbd72b2298f201d834fb54a5d39bca7ca5b86a552b97a5cbc64e23dd5e385", pass: true });
  expect(Object.fromEntries(Object.entries(linework.views).map(([view, value]: any) => [view, value.paths_mm.length]))).toEqual({ plan: 24, front: 68, side: 54 });
  expect(linework.native_dimension_labels_in_exact_cluster).toEqual(["1000", "450", "860", "940"]);
  const counterexample = readFileSync(join(product, "official-source/counterexamples/README.md"), "utf8");
  expect(counterexample).toContain("negative");
  expect(counterexample).toContain("must never be used");
});

test("candidate SVGs put white masks below native-DWG blue linework and inventory includes Marilyn", () => {
  const manifest = JSON.parse(readFileSync(join(product, "manifest.json"), "utf8"));
  const candidate = JSON.parse(readFileSync(join(product, "candidate-representations.json"), "utf8"));
  const inventory = JSON.parse(readFileSync(join(root, "pipeline/decisions/highpoly-product-inventory.json"), "utf8"));
  expect(inventory.review_ready_pending_approval_product_count).toBe(38);
  expect(inventory.products.find((entry: any) => entry.type_name === "Marilyn 01")).toMatchObject({ status: "review_ready_pending_approval", representative_global_id: "3l2Ji4k2H9oOTDGWYUq7uV" });
  expect(manifest).toMatchObject({ representative_global_id: "3l2Ji4k2H9oOTDGWYUq7uV", geometry_product_count: 1, whole_model_render: false, source_kind: "native_dwg", official_cad_used: true, third_party_cad_used: false, review_status: "visual_review_pending", approved_for_drawing_ifc: false, pass: true });
  expect(manifest.bounds_mm.size).toEqual([877.174957, 1013.911865, 970.476962]);
  expect(candidate).toMatchObject({ source_kind: "native_dwg", source_dwg_sha256: dwgHash, official_cad_used: true, third_party_cad_used: false, formal_ifc_write_allowed: false, review_status: "visual_review_pending" });
  for (const [view, count] of Object.entries({ plan: 24, front: 68, side: 54 })) {
    expect(candidate.views[view].official_native_dwg_paths_mm.length).toBe(count);
    const svg = readFileSync(join(product, `${view}.svg`), "utf8");
    expect(svg).toContain('class="official-reference-mask"');
    expect(svg).toContain('class="official-reference native-dwg"');
    expect(svg.indexOf('class="official-reference-mask"')).toBeLessThan(svg.indexOf('class="official-reference native-dwg"'));
    expect(svg).toContain(`data-source-sha256="${dwgHash}"`);
  }
});

test("project plan and two elevations retain walls and furniture at fixed native scale", () => {
  const context = JSON.parse(readFileSync(join(product, "project-context-manifest.json"), "utf8"));
  expect(context.views.map((view: any) => view.view)).toEqual(["plan", "front", "side"]);
  expect(context).toMatchObject({ source_kind: "native_dwg", official_cad_used: true, third_party_cad_used: false, project_context_retained: true, walls_and_surrounding_project_elements_retained: true, overlay_top_layer_with_white_mask: true, blue_line_present: true, project_scale_svg_units_per_mm: 0.02, pass: true });
  expect(context.context_view_scope).toMatchObject({ included: ["plan", "front", "side"], excluded: [] });
  expect(context.project_projection_note.geometry_stretched).toBeFalse();
  expect(context.review_annotation_suppression.walls_furniture_and_ifc_geometry_removed).toBeFalse();
  for (const view of context.views) {
    expect(view.overlay.fit.uniform_scale_preserved).toBeTrue();
    expect(view.overlay.fit.transformation_mode).toContain("native_1_50_scale_only");
    expect(view.overlay.path_count).toBe(({ plan: 24, front: 68, side: 54 } as any)[view.view]);
    expect(sha256(join(root, view.review_preview))).toBe(view.review_preview_sha256);
    const full = readFileSync(join(root, view.output), "utf8");
    expect(full).toContain('class="official-native-dwg-mask"');
    expect(full).toContain('class="official-native-dwg-blue"');
    expect(full.indexOf('class="official-native-dwg-mask"')).toBeLessThan(full.indexOf('class="official-native-dwg-blue"'));
  }
});

test("actual Bonsai evidence saves four cameras and renders the isolated MODEL_VIEW Body", () => {
  const evidence = JSON.parse(readFileSync(join(product, "bonsai-review-manifest.json"), "utf8"));
  expect(evidence).toMatchObject({ mode: "actual_bonsai_ifc_body_camera_render", representative_global_id: "3l2Ji4k2H9oOTDGWYUq7uV", geometry_product_count: 1, whole_model_render: false, formal_ifc_bytes_unchanged: true, pass: true });
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

test("pending Marilyn approval forbids every derived/formal IFC write", () => {
  const manifest = join(product, "manifest.json");
  const approval = JSON.parse(readFileSync(join(root, "pipeline/decisions/marilyn-01-drawing-approval.json"), "utf8"));
  expect(approval).toMatchObject({ status: "pending", candidate_manifest_sha256: sha256(manifest), derived_ifc_write_allowed: false, formal_authoritative_ifc_write_allowed: false });
  expect(existsSync(join(product, "Baxter-Marilyn-01-derived-drawing.ifc"))).toBeFalse();
  expect(existsSync(join(product, "Baxter-Marilyn-01-bonsai-review.blend1"))).toBeFalse();
  expect(sha256(formal)).toBe(formalHash);
  const temporary = mkdtempSync(join(tmpdir(), "marilyn-01-gate-"));
  const output = join(temporary, "forbidden.ifc");
  const script = join(root, "pipeline/scripts/marilyn_01_drawing_ifc.py");
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

test("temporary scoped approval writes verified official-DWG representations and source associations", () => {
  const manifest = join(product, "manifest.json");
  const temporary = mkdtempSync(join(tmpdir(), "marilyn-01-approved-"));
  const approvalPath = join(temporary, "approval.json");
  const output = join(temporary, "derived.ifc");
  const report = join(temporary, "report.json");
  writeFileSync(approvalPath, JSON.stringify({ schema_version: 1, profile_key: "marilyn-01", ifc_type_name: "Marilyn 01", candidate_manifest_sha256: sha256(manifest), status: "approved", reviewer: "automated gate fixture", review_date: "2026-08-23", approved_views: ["plan", "front", "side"], derived_ifc_write_allowed: true, formal_authoritative_ifc_write_allowed: false, scope, approval_evidence: "temporary automated writer verification only" }, null, 2));
  const run = Bun.spawnSync(["python3", join(root, "pipeline/scripts/marilyn_01_drawing_ifc.py"), "--input", formal, "--manifest", manifest, "--approval", approvalPath, "--output", output, "--report", report, "--apply"], { cwd: root, stdout: "pipe", stderr: "pipe" });
  expect(run.exitCode).toBe(0);
  if (run.exitCode !== 0) throw new Error(run.stderr.toString());
  const result = JSON.parse(readFileSync(report, "utf8"));
  expect(result).toMatchObject({ pass: true, formal_ifc_bytes_unchanged: true, representations: { plan: "Marilyn01Plan", front: "Marilyn01Front", side: "Marilyn01Side" }, representation_path_counts: { plan: 24, front: 68, side: 54 }, representation_geometry_source: "official_native_dwg_paths_mm", official_cad_geometry_included: true, proxy_geometry_included: false, source_property_set: "Pset_Marilyn01DrawingSource", source_kind: "native_dwg", source_label_zh: "基于 Baxter 精确型号原生 DWG 的官方图纸表达" });
  expect(result.source_document_associations).toContain("BAXTER-MARILYN-01-NATIVE-DWG");
  expect(result.source_document_associations).toContain("BAXTER-MARILYN-01-OFFICIAL-PRODUCT-PAGE");
  expect(result.source_document_associations).toContain("BAXTER-MARILYN-01-NATIVE-LINEWORK-REGISTER");
  expect(result.source_document_associations.length).toBe(12);
  expect(sha256(formal)).toBe(formalHash);
  rmSync(temporary, { recursive: true, force: true });
}, 30_000);
