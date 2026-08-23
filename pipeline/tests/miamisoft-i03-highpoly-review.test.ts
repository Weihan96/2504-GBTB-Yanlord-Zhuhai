import { expect, test } from "bun:test";
import { createHash } from "node:crypto";
import { existsSync, mkdtempSync, readFileSync, rmSync, writeFileSync } from "node:fs";
import { tmpdir } from "node:os";
import { join, resolve } from "node:path";

const root = resolve(import.meta.dir, "../..");
const product = join(root, "output/review/highpoly-types/miamisoft-i03");
const formal = join(root, "2504 GBTB Yanlord Zhuhai.ifc");
const formalHash = "7a50b87e8f48a7c2bfbaa9f1dabdfca7155fa684b2325c66aa1ae15c8ab0a25c";
const sourceKind = "geometry_derived_simplified_proxy";
const sourceLabelZh = "基于原始高模几何生成的简化图纸表达";
const scope = "exact manufacturer family, I03 roll model and nominal dimensions; project IFC Description incorrectly states the I01 108 cm length while type code and Body geometry match I03; authenticated native 2D/3D/BIM not acquired; not a project shop drawing and not official CAD geometry";

function sha256(path: string) {
  return createHash("sha256").update(readFileSync(path)).digest("hex");
}

test("exact Baxter Miami Soft I03 roll identity is archived while the stale IFC Description is disclosed", () => {
  const access = JSON.parse(readFileSync(join(product, "official-source/source-access-record.json"), "utf8"));
  expect(access).toMatchObject({
    manufacturer: "Baxter",
    family: "Miami Soft",
    designer: "Paola Navone",
    project_ifc_type_name: "MiamiSoft I03",
    project_ifc_type_description: "Roll 108 x Ø27 cm",
    resolved_variant: "I03 - roll cushion - 170 x Ø27 cm",
    scope,
    pass: true,
  });
  expect(access.project_identity_evidence).toMatchObject({
    exact_model_code_status: "confirmed",
    module_status: "confirmed_I03_roll_170cm",
    description_status: "conflicts_with_official_I03_and_matches_adjacent_I01_length",
    body_geometry_status: "matches_official_I03_within_15mm",
  });
  expect(access.official_product_cad).toMatchObject({
    authentication_required: true,
    acquired: false,
    exact_project_configuration_match: false,
    local_cad_files: [],
    third_party_cad_used: false,
    official_vector_pdf_archived: true,
    official_measurement_svg_archived: true,
    official_vector_references_used_as_cad_geometry: false,
  });
  expect(access.dimension_cross_check).toMatchObject({
    official_nominal_length_diameter_mm: [1700, 270],
    project_ifc_body_local_xyz_mm: [1711.749878, 274.261993, 272.80899],
    body_minus_official_length_y_diameter_z_diameter_mm: [11.749878, 4.261993, 2.80899],
    maximum_absolute_delta_mm: 11.749878,
    tolerance_mm: 15,
    pass: true,
  });
  const pdf = access.official_identity_sources.find((source: any) => source.kind === "manufacturer_vector_technical_sheet");
  expect(pdf.page).toBe(8);
  expect(pdf.vector_object_audit).toMatchObject({ drawing_count: 56, drawing_item_count: 1856, word_count: 72, character_count: 255 });
  expect(pdf.vector_object_audit.required_text_found).toEqual(["I03 - Rullo. Roll.", "170 x Ø27 cm"]);
  expect(pdf.vector_object_audit.adjacent_models_excluded).toEqual(["I01 - 108 x Ø27 cm", "I02 - 140 x Ø27 cm"]);
  const svg = access.official_identity_sources.find((source: any) => source.kind === "manufacturer_measurement_svg");
  expect(svg.vector_object_audit).toEqual({ path_count: 22, use_count: 20, group_count: 28 });
  for (const source of access.official_identity_sources) {
    expect(sha256(join(root, source.local_path))).toBe(source.sha256);
    if (source.preview_path) expect(sha256(join(root, source.preview_path))).toBe(source.preview_sha256);
  }
});

test("one Miami Soft I03 Body produces three geometry-derived candidates with zero blue CAD paths", () => {
  const manifest = JSON.parse(readFileSync(join(product, "manifest.json"), "utf8"));
  const candidate = JSON.parse(readFileSync(join(product, "candidate-representations.json"), "utf8"));
  const inventory = JSON.parse(readFileSync(join(root, "pipeline/decisions/highpoly-product-inventory.json"), "utf8"));
  const inventoryEntry = inventory.products.find((entry: any) => entry.type_name === "MiamiSoft I03");
  expect(inventoryEntry.status).toBe("review_ready_pending_approval");
  expect(inventory.review_ready_pending_approval_product_count).toBe(38);
  expect(manifest).toMatchObject({ representative_global_id: "1hdNtfQRPCGfHLfqzdp6w0", geometry_product_count: 1, whole_model_render: false, official_cad_acquired: false, official_cad_used: false, blue_line_present: false, pass: true });
  expect(manifest.bounds_mm.size).toEqual([1711.749878, 274.261993, 272.80899]);
  expect(candidate).toMatchObject({ article_number: "Baxter Miami Soft I03 roll cushion", source_kind: sourceKind, source_label_zh: sourceLabelZh, official_cad_used: false, third_party_cad_used: false, formal_ifc_write_allowed: false, review_status: "visual_review_pending" });
  expect(Object.fromEntries(Object.entries(candidate.views).map(([view, value]: any) => [view, value.proxy_paths_mm.length]))).toEqual({ plan: 1, front: 1, side: 1 });
  for (const view of ["plan", "front", "side"]) {
    expect(candidate.views[view].official_cad_paths_mm).toEqual([]);
    const svgText = readFileSync(join(product, `${view}.svg`), "utf8");
    expect(svgText).toContain(`data-source-kind="${sourceKind}"`);
    expect(svgText).not.toContain('class="official-reference native-dwg"');
    expect(svgText).not.toContain("#1677c8");
  }
});

test("Miami Soft I03 project context retains walls and furniture with an exact complete side fit", () => {
  const context = JSON.parse(readFileSync(join(product, "project-context-manifest.json"), "utf8"));
  expect(context.views.map((view: any) => view.view)).toEqual(["plan", "side"]);
  expect(context.context_view_scope).toMatchObject({ included: ["plan", "side"], excluded: ["front"] });
  expect(context).toMatchObject({ source_kind: sourceKind, official_cad_used: false, project_context_retained: true, walls_and_surrounding_project_elements_retained: true, overlay_top_layer_with_white_mask: true, blue_line_present: false, mechanical_side_fit_tolerance_svg_units: 0.08, mechanical_side_fit_pass: true, pass: true });
  expect(context.project_projection_note.geometry_stretched).toBeFalse();
  const plan = context.views.find((view: any) => view.view === "plan");
  const side = context.views.find((view: any) => view.view === "side");
  expect(plan.overlay.fit.rotate_quarter_turns).toBe(1);
  expect(plan.overlay.fit.uniform_scale_preserved).toBeTrue();
  expect(Math.max(...plan.overlay.fit.bbox_absolute_delta_svg_units)).toBeLessThanOrEqual(0.25);
  expect(Math.max(...side.overlay.fit.bbox_absolute_delta_svg_units)).toBeLessThanOrEqual(0.08);
  for (const view of context.views) {
    expect(sha256(join(root, view.review_preview))).toBe(view.review_preview_sha256);
    const full = readFileSync(join(root, view.output), "utf8");
    expect(full).toContain('class="geometry-derived-proxy-mask" fill="#ffffff" fill-rule="evenodd"');
    expect(full.indexOf('class="geometry-derived-proxy-mask"')).toBeLessThan(full.indexOf('class="geometry-derived-proxy"'));
  }
});

test("Miami Soft I03 Bonsai evidence uses four cameras on the actual MODEL_VIEW IFC Body", () => {
  const evidence = JSON.parse(readFileSync(join(product, "bonsai-review-manifest.json"), "utf8"));
  expect(evidence).toMatchObject({ mode: "actual_bonsai_ifc_body_camera_render", representative_global_id: "1hdNtfQRPCGfHLfqzdp6w0", geometry_product_count: 1, whole_model_render: false, formal_ifc_bytes_unchanged: true, pass: true });
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

test("pending Miami Soft I03 review leaves formal IFC unchanged and creates no derived IFC", () => {
  const candidate = JSON.parse(readFileSync(join(product, "candidate-representations.json"), "utf8"));
  expect(candidate.formal_ifc_write_allowed).toBeFalse();
  expect(candidate.review_status).toBe("visual_review_pending");
  expect(existsSync(join(product, "Baxter-Miami-Soft-I03-derived-drawing.ifc"))).toBeFalse();
  expect(existsSync(join(product, "Baxter-Miami-Soft-I03-bonsai-review.blend1"))).toBeFalse();
  expect(sha256(formal)).toBe(formalHash);
});

test("Miami Soft I03 writer rejects missing apply and pending human approval", () => {
  const temporary = mkdtempSync(join(tmpdir(), "miamisoft-i03-gate-"));
  const output = join(temporary, "forbidden.ifc");
  const script = join(root, "pipeline/scripts/miamisoft_i03_drawing_ifc.py");
  const approval = JSON.parse(readFileSync(join(root, "pipeline/decisions/miamisoft-i03-drawing-approval.json"), "utf8"));
  const manifest = join(product, "manifest.json");
  expect(approval).toMatchObject({ status: "pending", derived_ifc_write_allowed: false, formal_authoritative_ifc_write_allowed: false });
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

test("scoped temporary approval writes verified Miami Soft I03 representations and source documents", () => {
  const manifest = join(product, "manifest.json");
  const temporary = mkdtempSync(join(tmpdir(), "miamisoft-i03-approved-"));
  const approvalPath = join(temporary, "approval.json");
  const output = join(temporary, "derived.ifc");
  const report = join(temporary, "report.json");
  writeFileSync(approvalPath, JSON.stringify({ schema_version: 1, profile_key: "miamisoft-i03", ifc_type_name: "MiamiSoft I03", candidate_manifest_sha256: sha256(manifest), status: "approved", reviewer: "automated gate fixture", review_date: "2026-08-23", approved_views: ["plan", "front", "side"], derived_ifc_write_allowed: true, formal_authoritative_ifc_write_allowed: false, scope, approval_evidence: "temporary automated writer verification only" }, null, 2));
  const run = Bun.spawnSync(["python3", join(root, "pipeline/scripts/miamisoft_i03_drawing_ifc.py"), "--input", formal, "--manifest", manifest, "--approval", approvalPath, "--output", output, "--report", report, "--apply"], { cwd: root, stdout: "pipe", stderr: "pipe" });
  expect(run.exitCode).toBe(0);
  if (run.exitCode !== 0) throw new Error(run.stderr.toString());
  const result = JSON.parse(readFileSync(report, "utf8"));
  expect(result).toMatchObject({ pass: true, formal_ifc_bytes_unchanged: true, representations: { plan: "MiamiSoftI03Plan", front: "MiamiSoftI03Front", side: "MiamiSoftI03Side" }, representation_path_counts: { plan: 1, front: 1, side: 1 }, source_kind: sourceKind, source_label_zh: sourceLabelZh, official_cad_geometry_included: false, source_property_set: "Pset_MiamiSoftI03DrawingSource" });
  expect(result.source_document_associations).toEqual([
    "BAXTER-MIAMI-SOFT-I03-OFFICIAL-MEASUREMENT-SVG",
    "BAXTER-MIAMI-SOFT-I03-OFFICIAL-MEASUREMENT-SVG-ARCHIVE",
    "BAXTER-MIAMI-SOFT-I03-OFFICIAL-PRODUCT-PAGE",
    "BAXTER-MIAMI-SOFT-I03-OFFICIAL-PRODUCT-PAGE-ARCHIVE",
    "BAXTER-MIAMI-SOFT-I03-OFFICIAL-TECHNICAL-SHEET",
    "BAXTER-MIAMI-SOFT-I03-OFFICIAL-TECHNICAL-SHEET-ARCHIVE",
    "BAXTER-MIAMI-SOFT-I03-SOURCE-ACCESS-RECORD",
  ]);
  expect(sha256(formal)).toBe(formalHash);
  rmSync(temporary, { recursive: true, force: true });
}, 30_000);
