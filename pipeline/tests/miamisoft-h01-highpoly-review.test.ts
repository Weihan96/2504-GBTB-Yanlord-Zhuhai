import { expect, test } from "bun:test";
import { createHash } from "node:crypto";
import { existsSync, mkdtempSync, readFileSync, rmSync, writeFileSync } from "node:fs";
import { tmpdir } from "node:os";
import { join, resolve } from "node:path";

const root = resolve(import.meta.dir, "../..");
const product = join(root, "output/review/highpoly-types/miamisoft-h01");
const formal = join(root, "2504 GBTB Yanlord Zhuhai.ifc");
const formalHash = "7a50b87e8f48a7c2bfbaa9f1dabdfca7155fa684b2325c66aa1ae15c8ab0a25c";
const sourceKind = "geometry_derived_simplified_proxy";
const sourceLabelZh = "基于原始高模几何生成的简化图纸表达";
const scope = "exact manufacturer family, H01 cushion model and nominal flexible face dimensions; authenticated native 2D/3D/BIM not acquired; not a project shop drawing and not official CAD geometry";

function sha256(path: string) {
  return createHash("sha256").update(readFileSync(path)).digest("hex");
}

test("exact Baxter Miami Soft H01 flexible cushion identity is archived without claiming native CAD", () => {
  const access = JSON.parse(readFileSync(join(product, "official-source/source-access-record.json"), "utf8"));
  expect(access).toMatchObject({
    manufacturer: "Baxter",
    family: "Miami Soft",
    designer: "Paola Navone",
    project_ifc_type_name: "MiamiSoft H01",
    project_ifc_type_description: "Cushion 80 x 80 cm",
    resolved_variant: "H01 - cushion - 80 x 80 cm",
    scope,
    pass: true,
  });
  expect(access.project_identity_evidence).toMatchObject({ exact_model_code_status: "confirmed", module_status: "confirmed_H01_cushion" });
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
    official_nominal_flexible_face_width_height_mm: [800, 800],
    project_ifc_body_local_xyz_mm: [809.757111, 564.975739, 403.963976],
    width_delta_mm: 9.757111,
    width_tolerance_mm: 15,
    width_check_pass: true,
    full_rigid_box_check_applicable: false,
    pass: true,
  });
  const pdf = access.official_identity_sources.find((source: any) => source.kind === "manufacturer_vector_technical_sheet");
  expect(pdf.page).toBe(8);
  expect(pdf.vector_object_audit).toMatchObject({ method: "pdfplumber page 8 vector object audit", drawing_count: 56, drawing_item_count: 1856, word_count: 72, character_count: 255 });
  expect(pdf.vector_object_audit.required_text_found).toEqual(["H01 - Cuscino. Cushion.", "80 x 80 cm"]);
  const svg = access.official_identity_sources.find((source: any) => source.kind === "manufacturer_measurement_svg");
  expect(svg.vector_object_audit).toEqual({ path_count: 50, use_count: 19, group_count: 28 });
  for (const source of access.official_identity_sources) {
    expect(sha256(join(root, source.local_path))).toBe(source.sha256);
    if (source.preview_path) expect(sha256(join(root, source.preview_path))).toBe(source.preview_sha256);
  }
});

test("one Miami Soft H01 Body produces three geometry-derived candidates with zero blue CAD paths", () => {
  const manifest = JSON.parse(readFileSync(join(product, "manifest.json"), "utf8"));
  const candidate = JSON.parse(readFileSync(join(product, "candidate-representations.json"), "utf8"));
  const inventory = JSON.parse(readFileSync(join(root, "pipeline/decisions/highpoly-product-inventory.json"), "utf8"));
  const inventoryEntry = inventory.products.find((entry: any) => entry.type_name === "MiamiSoft H01");
  expect(inventoryEntry.status).toBe("review_ready_pending_approval");
  expect(inventory.review_ready_pending_approval_product_count).toBe(38);
  expect(manifest).toMatchObject({ representative_global_id: "0B8rO47U14sgJ_ydZxqvs6", registered_instance_global_ids: ["0B8rO47U14sgJ_ydZxqvs6", "0b9rDAHDD2qAIQqZ8m9INg"], geometry_product_count: 1, whole_model_render: false, official_cad_acquired: false, official_cad_used: false, blue_line_present: false, pass: true });
  expect(manifest.bounds_mm.size).toEqual([809.757111, 564.975739, 403.963976]);
  expect(candidate).toMatchObject({ article_number: "Baxter Miami Soft H01 cushion", source_kind: sourceKind, source_label_zh: sourceLabelZh, official_cad_used: false, third_party_cad_used: false, formal_ifc_write_allowed: false, review_status: "visual_review_pending" });
  expect(Object.fromEntries(Object.entries(candidate.views).map(([view, value]: any) => [view, value.proxy_paths_mm.length]))).toEqual({ plan: 3, front: 3, side: 5 });
  for (const view of ["plan", "front", "side"]) {
    expect(candidate.views[view].official_cad_paths_mm).toEqual([]);
    const svgText = readFileSync(join(product, `${view}.svg`), "utf8");
    expect(svgText).toContain(`data-source-kind="${sourceKind}"`);
    expect(svgText).not.toContain('class="official-reference native-dwg"');
    expect(svgText).not.toContain("#1677c8");
  }
});

test("Miami Soft H01 project context replaces the stale representative PLAN_VIEW without removing other furniture", () => {
  const context = JSON.parse(readFileSync(join(product, "project-context-manifest.json"), "utf8"));
  expect(context.views.map((view: any) => view.view)).toEqual(["plan", "side"]);
  expect(context.context_view_scope).toMatchObject({ included: ["plan", "side"], excluded: ["front"] });
  expect(context).toMatchObject({ source_kind: sourceKind, official_cad_used: false, project_context_retained: true, walls_and_surrounding_project_elements_retained: true, overlay_top_layer_with_white_mask: true, blue_line_present: false, mechanical_side_fit_tolerance_svg_units: 0.15, mechanical_side_fit_pass: true, pass: true });
  expect(context.project_projection_note.geometry_stretched).toBeFalse();
  expect(context.plan_source_representation_comparison).toMatchObject({ source_representation: "existing Bonsai PLAN_VIEW projection group", replacement_representation: "actual MODEL_VIEW Body projection", representative_source_group_hidden_in_review_copy: true, other_project_elements_retained: true, geometry_stretched_to_match_source: false });
  expect(context.review_product_projection_replacement).toMatchObject({ representative_global_id: "0B8rO47U14sgJ_ydZxqvs6", source_projection_groups_hidden_per_view: 1, replacement_proxy_top_layer: true, second_type_instance_retained: "0b9rDAHDD2qAIQqZ8m9INg" });
  const plan = context.views.find((view: any) => view.view === "plan");
  const side = context.views.find((view: any) => view.view === "side");
  expect(plan.overlay.fit.rotate_quarter_turns).toBe(0);
  expect(plan.overlay.fit.uniform_scale_preserved).toBeTrue();
  expect(Math.max(...plan.overlay.fit.bbox_absolute_delta_svg_units)).toBeLessThanOrEqual(3.1);
  expect(Math.max(...side.overlay.fit.bbox_absolute_delta_svg_units)).toBeLessThanOrEqual(0.15);
  for (const view of context.views) {
    expect(sha256(join(root, view.review_preview))).toBe(view.review_preview_sha256);
    const full = readFileSync(join(root, view.output), "utf8");
    expect(full).toContain('class="geometry-derived-proxy-mask" fill="#ffffff" fill-rule="evenodd"');
    expect(full).toContain('data-review-replaced-by="geometry-derived-body-proxy"');
    expect(full).toContain('ifc:guid="0b9rDAHDD2qAIQqZ8m9INg"');
    expect(full.indexOf('class="geometry-derived-proxy-mask"')).toBeLessThan(full.indexOf('class="geometry-derived-proxy"'));
  }
});

test("Miami Soft H01 Bonsai evidence uses four cameras on the actual MODEL_VIEW IFC Body", () => {
  const evidence = JSON.parse(readFileSync(join(product, "bonsai-review-manifest.json"), "utf8"));
  expect(evidence).toMatchObject({ mode: "actual_bonsai_ifc_body_camera_render", representative_global_id: "0B8rO47U14sgJ_ydZxqvs6", geometry_product_count: 1, whole_model_render: false, formal_ifc_bytes_unchanged: true, pass: true });
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

test("pending Miami Soft H01 review leaves formal IFC unchanged and creates no derived IFC", () => {
  const candidate = JSON.parse(readFileSync(join(product, "candidate-representations.json"), "utf8"));
  expect(candidate.formal_ifc_write_allowed).toBeFalse();
  expect(candidate.review_status).toBe("visual_review_pending");
  expect(existsSync(join(product, "Baxter-Miami-Soft-H01-derived-drawing.ifc"))).toBeFalse();
  expect(existsSync(join(product, "Baxter-Miami-Soft-H01-bonsai-review.blend1"))).toBeFalse();
  expect(sha256(formal)).toBe(formalHash);
});

test("Miami Soft H01 writer rejects missing apply and pending human approval", () => {
  const temporary = mkdtempSync(join(tmpdir(), "miamisoft-h01-gate-"));
  const output = join(temporary, "forbidden.ifc");
  const script = join(root, "pipeline/scripts/miamisoft_h01_drawing_ifc.py");
  const approval = JSON.parse(readFileSync(join(root, "pipeline/decisions/miamisoft-h01-drawing-approval.json"), "utf8"));
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

test("scoped temporary approval writes verified Miami Soft H01 representations and source documents", () => {
  const manifest = join(product, "manifest.json");
  const temporary = mkdtempSync(join(tmpdir(), "miamisoft-h01-approved-"));
  const approvalPath = join(temporary, "approval.json");
  const output = join(temporary, "derived.ifc");
  const report = join(temporary, "report.json");
  writeFileSync(approvalPath, JSON.stringify({ schema_version: 1, profile_key: "miamisoft-h01", ifc_type_name: "MiamiSoft H01", candidate_manifest_sha256: sha256(manifest), status: "approved", reviewer: "automated gate fixture", review_date: "2026-08-23", approved_views: ["plan", "front", "side"], derived_ifc_write_allowed: true, formal_authoritative_ifc_write_allowed: false, scope, approval_evidence: "temporary automated writer verification only" }, null, 2));
  const run = Bun.spawnSync(["python3", join(root, "pipeline/scripts/miamisoft_h01_drawing_ifc.py"), "--input", formal, "--manifest", manifest, "--approval", approvalPath, "--output", output, "--report", report, "--apply"], { cwd: root, stdout: "pipe", stderr: "pipe" });
  expect(run.exitCode).toBe(0);
  if (run.exitCode !== 0) throw new Error(run.stderr.toString());
  const result = JSON.parse(readFileSync(report, "utf8"));
  expect(result).toMatchObject({ pass: true, formal_ifc_bytes_unchanged: true, representations: { plan: "MiamiSoftH01Plan", front: "MiamiSoftH01Front", side: "MiamiSoftH01Side" }, representation_path_counts: { plan: 3, front: 3, side: 5 }, source_kind: sourceKind, source_label_zh: sourceLabelZh, official_cad_geometry_included: false, source_property_set: "Pset_MiamiSoftH01DrawingSource" });
  expect(result.source_document_associations).toEqual([
    "BAXTER-MIAMI-SOFT-H01-OFFICIAL-MEASUREMENT-SVG",
    "BAXTER-MIAMI-SOFT-H01-OFFICIAL-MEASUREMENT-SVG-ARCHIVE",
    "BAXTER-MIAMI-SOFT-H01-OFFICIAL-PRODUCT-PAGE",
    "BAXTER-MIAMI-SOFT-H01-OFFICIAL-PRODUCT-PAGE-ARCHIVE",
    "BAXTER-MIAMI-SOFT-H01-OFFICIAL-TECHNICAL-SHEET",
    "BAXTER-MIAMI-SOFT-H01-OFFICIAL-TECHNICAL-SHEET-ARCHIVE",
    "BAXTER-MIAMI-SOFT-H01-SOURCE-ACCESS-RECORD",
  ]);
  expect(sha256(formal)).toBe(formalHash);
  rmSync(temporary, { recursive: true, force: true });
}, 30_000);
