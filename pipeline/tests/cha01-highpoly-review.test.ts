import { expect, test } from "bun:test";
import { createHash } from "node:crypto";
import { existsSync, mkdtempSync, readFileSync, rmSync, writeFileSync } from "node:fs";
import { tmpdir } from "node:os";
import { join, resolve } from "node:path";

const root = resolve(import.meta.dir, "../..");
const product = join(root, "output/review/highpoly-types/cha01");
const formal = join(root, "2504 GBTB Yanlord Zhuhai.ifc");
const formalHash = "7a50b87e8f48a7c2bfbaa9f1dabdfca7155fa684b2325c66aa1ae15c8ab0a25c";
const sourceKind = "geometry_derived_simplified_proxy";
const sourceLabelZh = "基于原始高模几何生成的简化图纸表达";
const scope = "exact Baxter Colette armchair 57 x 60 x 73 cm manufacturer identity; authenticated 2D/3D/BIM not acquired; technical PDFs and measurement SVG are identity evidence only; not a project shop drawing or official CAD geometry";

function sha256(path: string) {
  return createHash("sha256").update(readFileSync(path)).digest("hex");
}

test("exact Baxter Colette identity sources are archived without claiming CAD geometry", () => {
  const access = JSON.parse(readFileSync(join(product, "official-source/source-access-record.json"), "utf8"));
  expect(access).toMatchObject({
    manufacturer: "Baxter",
    family: "Colette",
    model: "Armchair / Poltroncina 57 x 60 x 73 cm",
    designer: "Roberto Lazzeroni",
    project_type_code: "CHA01",
    project_ifc_type_description: "Baxter Colette",
    scope,
    pass: true,
  });
  expect(access.official_product_cad).toMatchObject({
    authentication_required: true,
    acquired: false,
    exact_project_configuration_match: false,
    local_cad_files: [],
    third_party_cad_used: false,
  });
  expect(access.drawing_geometry_source).toEqual({
    source_kind: sourceKind,
    source_label_zh: sourceLabelZh,
    source_label_en: "simplified drawing representation derived from the original high-poly geometry",
    official_cad_used: false,
    third_party_cad_used: false,
  });
  expect(access.dimension_cross_check).toMatchObject({
    official_model_width_depth_height_mm: [570, 600, 730],
    project_ifc_body_local_xyz_mm: [568.466003, 593.482941, 734.464966],
    body_minus_official_mm: [-1.533997, -6.517059, 4.464966],
    maximum_absolute_delta_mm: 6.517059,
    status: "exact_model_geometry_consistent_no_scaling_or_fitting",
  });
  const expectedHashes: Record<string, string> = {
    "baxter-colette-product-page.html": "7aa860785daff5f259c8eeadea793d46d71fba90e964d573d6d481f35de57b9b",
    "baxter-login-page.html": "a40f10c226066c93c534522d8306ec225f6231ad16fbf33732c98f8af2472976",
    "Baxter_Colette_current-technical-sheet.pdf": "6f22cfa8feeb4b1cece737a1123f22dce6674cd39659b630ea5f88b6ab1c8239",
    "Baxter_Colette_current-technical-sheet-page-37-preview.png": "7b9c94987d7de04f6d2f19e1d9406b733d50baafddf273d3397cc4376f45d823",
    "Baxter_Colette_sedia_Indoor.pdf": "6e179cedc5689addbdcaffbd015fed45cbd811b5dfd6d03d9afc84734c8c4a5c",
    "Baxter_Colette_technical-sheet-page-6-preview.png": "b552095147e8f6362032c7cd643d017c164ec2d4bc9a864c72aff98251f452a3",
    "COLEPOCO57.svg": "637678ebd969fdb7600b11bf185fb7ebeb18d079c69eec42d591047bd22d0cb6",
  };
  for (const [name, hash] of Object.entries(expectedHashes)) {
    expect(sha256(join(product, "official-source", name))).toBe(hash);
  }
  expect(access.official_identity_sources[1]).toMatchObject({
    pdf_page_count: 40,
    exact_variant_page: 37,
    vector_object_audit: { drawing_count: 615, line_count: 202, curve_count: 289, word_count: 179, character_count: 957 },
  });
  expect(access.official_identity_sources[2]).toMatchObject({
    pdf_page_count: 6,
    exact_variant_page: 6,
    vector_object_audit: { line_count: 1270, curve_count: 950, word_count: 369, character_count: 1550 },
  });
  expect(access.official_identity_sources[3].vector_object_audit).toEqual({ path_count: 94, use_count: 54, group_count: 52 });
  for (const source of access.official_identity_sources) {
    expect(sha256(join(root, source.local_path))).toBe(source.sha256);
  }
});

test("CHA01 candidate isolates one Body and exposes three non-CAD views", () => {
  const manifest = JSON.parse(readFileSync(join(product, "manifest.json"), "utf8"));
  const candidate = JSON.parse(readFileSync(join(product, "candidate-representations.json"), "utf8"));
  const inventory = JSON.parse(readFileSync(join(root, "pipeline/decisions/highpoly-product-inventory.json"), "utf8"));
  expect(inventory.review_ready_pending_approval_product_count).toBe(38);
  expect(inventory.products.find((entry: any) => entry.type_name === "CHA01")).toMatchObject({
    status: "review_ready_pending_approval",
    representative_global_id: "1luHljRzDAhPNXxTDNu7qB",
    instance_count: 2,
  });
  expect(manifest).toMatchObject({
    ifc_type_name: "CHA01",
    ifc_type_description: "Baxter Colette",
    representative_global_id: "1luHljRzDAhPNXxTDNu7qB",
    registered_instance_global_ids: ["1luHljRzDAhPNXxTDNu7qB", "2tCmGkmFz2_Pg0jRjf3Nup"],
    geometry_product_count: 1,
    whole_model_render: false,
    source_kind: sourceKind,
    official_cad_used: false,
    third_party_cad_used: false,
    blue_line_present: false,
    formal_ifc_bytes_unchanged: true,
    review_status: "visual_review_pending",
    approved_for_drawing_ifc: false,
    pass: true,
  });
  expect(manifest.bounds_mm.size).toEqual([568.466003, 593.482941, 734.464966]);
  expect(candidate).toMatchObject({ source_kind: sourceKind, source_label_zh: sourceLabelZh, official_cad_used: false, third_party_cad_used: false, formal_ifc_write_allowed: false, review_status: "visual_review_pending" });
  expect(Object.fromEntries(Object.entries(candidate.views).map(([view, value]: any) => [view, value.proxy_paths_mm.length]))).toEqual({ plan: 1, front: 6, side: 4 });
  for (const view of manifest.views) {
    expect(view.official_cad_path_count).toBe(0);
    expect(view.blue_line_present).toBeFalse();
    const svg = readFileSync(join(root, view.svg), "utf8");
    expect(svg).toContain('class="simplified-proxy-silhouette geometry-derived"');
    expect(svg).not.toContain('class="official-reference native-dwg"');
    expect(svg).not.toContain("#1677c8");
  }
});

test("project plan and R07 elevations keep context and replace both stale projections at 1:50", () => {
  const context = JSON.parse(readFileSync(join(product, "project-context-manifest.json"), "utf8"));
  expect(context).toMatchObject({ source_kind: sourceKind, source_label_zh: sourceLabelZh, official_cad_used: false, third_party_cad_used: false, project_context_retained: true, walls_and_surrounding_project_elements_retained: true, overlay_top_layer_with_white_mask: true, blue_line_present: false, project_scale_svg_units_per_mm: 0.02, pass: true });
  expect(context.views.map((view: any) => [view.view, view.candidate_view])).toEqual([["plan", "plan"], ["front", "front"], ["side", "side"]]);
  expect(context.review_annotation_suppression.walls_furniture_and_ifc_geometry_removed).toBeFalse();
  const expectedCounts: Record<string, number> = { plan: 1, front: 6, side: 4 };
  for (const view of context.views) {
    expect(view.overlays.map((overlay: any) => overlay.ifc_guid)).toEqual(["1luHljRzDAhPNXxTDNu7qB", "2tCmGkmFz2_Pg0jRjf3Nup"]);
    expect(sha256(join(root, view.output))).toBe(view.output_sha256);
    expect(sha256(join(root, view.review_crop))).toBe(view.review_crop_sha256);
    expect(sha256(join(root, view.review_preview))).toBe(view.review_preview_sha256);
    for (const overlay of view.overlays) {
      expect(overlay.path_count).toBe(expectedCounts[view.view]);
      expect(overlay.original_project_projection_suppressed_in_review_copy).toBeTrue();
      expect(overlay.fit.uniform_scale_preserved).toBeTrue();
      expect(overlay.fit.geometry_stretched).toBeFalse();
      expect(overlay.fit.scale_svg_units_per_mm).toBe(0.02);
      expect(overlay.fit.transformation_mode).toContain("native_1_50_scale_only");
    }
    const full = readFileSync(join(root, view.output), "utf8");
    expect(full).toContain("IfcFurniture");
    expect(full.match(/class="geometry-derived-simplified-proxy project-context-overlay"/g)?.length).toBe(2);
    expect(full.match(/class="geometry-derived-proxy-mask"/g)?.length).toBe(2);
    expect(full.match(/class="geometry-derived-proxy"/g)?.length).toBe(2);
    expect(full).not.toContain('class="official-reference native-dwg"');
  }
  const plan = context.views.find((view: any) => view.view === "plan");
  const front = context.views.find((view: any) => view.view === "front");
  const side = context.views.find((view: any) => view.view === "side");
  expect(plan.overlays[0].fit.bbox_absolute_delta_svg_units).toEqual([2.623389, 0.03068]);
  expect(front.overlays[0].fit.bbox_absolute_delta_svg_units).toEqual([0, 0.5493]);
  expect(side.overlays[0].fit.bbox_absolute_delta_svg_units).toEqual([0, 0]);
  expect(readFileSync(join(root, plan.review_crop), "utf8")).not.toContain('class="official-elevation-marker"');
  expect(readFileSync(join(root, front.review_crop), "utf8")).not.toContain('id="noninteger-highlights"');
  expect(readFileSync(join(root, side.review_crop), "utf8")).not.toContain('id="noninteger-highlights"');
});

test("actual Bonsai evidence saves four cameras and renders the isolated MODEL_VIEW Body", () => {
  const evidence = JSON.parse(readFileSync(join(product, "bonsai-review-manifest.json"), "utf8"));
  expect(evidence).toMatchObject({ mode: "actual_bonsai_ifc_body_camera_render", representative_global_id: "1luHljRzDAhPNXxTDNu7qB", geometry_product_count: 1, whole_model_render: false, formal_ifc_bytes_unchanged: true, pass: true });
  expect(evidence.bonsai_session).toMatchObject({ saved_active_representation: "Body", ifc_context_identifier: "Body", ifc_target_view: "MODEL_VIEW", saved_camera_count: 4, front_camera_local_y_sign: -1 });
  expect(evidence.renders.map((render: any) => render.view)).toEqual(["plan", "front", "side", "iso"]);
  for (const render of evidence.renders) {
    expect(render.camera_type).toBe("ORTHO");
    expect(render.render_operation).toBe("bpy.ops.render.render(write_still=True)");
    expect(render.camera_ortho_scale_m + 2e-7).toBeGreaterThanOrEqual(render.projected_width_m * render.framing_margin_factor);
    expect(render.camera_ortho_scale_m + 2e-7).toBeGreaterThanOrEqual(render.projected_height_m * render.image_aspect * render.framing_margin_factor);
    expect(sha256(join(root, render.path))).toBe(render.sha256);
  }
  expect(sha256(join(root, evidence.bonsai_session.path))).toBe(evidence.bonsai_session.sha256);
  expect(sha256(join(root, evidence.isolated_ifc))).toBe(evidence.isolated_ifc_sha256);
});

test("pending CHA01 approval forbids every derived and formal IFC write", () => {
  const manifest = join(product, "manifest.json");
  const approval = JSON.parse(readFileSync(join(root, "pipeline/decisions/cha01-drawing-approval.json"), "utf8"));
  expect(approval).toMatchObject({ status: "pending", candidate_manifest_sha256: sha256(manifest), derived_ifc_write_allowed: false, formal_authoritative_ifc_write_allowed: false });
  expect(existsSync(join(product, "Baxter-Colette-CHA01-derived-drawing.ifc"))).toBeFalse();
  expect(existsSync(join(product, "Baxter-Colette-CHA01-bonsai-review.blend1"))).toBeFalse();
  expect(sha256(formal)).toBe(formalHash);
  const temporary = mkdtempSync(join(tmpdir(), "cha01-gate-"));
  const output = join(temporary, "forbidden.ifc");
  const script = join(root, "pipeline/scripts/cha01_drawing_ifc.py");
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

test("temporary scoped approval writes verified geometry-derived representations and source associations", () => {
  const manifest = join(product, "manifest.json");
  const temporary = mkdtempSync(join(tmpdir(), "cha01-approved-"));
  const approvalPath = join(temporary, "approval.json");
  const output = join(temporary, "derived.ifc");
  const report = join(temporary, "report.json");
  writeFileSync(approvalPath, JSON.stringify({ schema_version: 1, profile_key: "cha01", ifc_type_name: "CHA01", candidate_manifest_sha256: sha256(manifest), status: "approved", reviewer: "automated gate fixture", review_date: "2026-08-23", approved_views: ["plan", "front", "side"], derived_ifc_write_allowed: true, formal_authoritative_ifc_write_allowed: false, scope, approval_evidence: "temporary automated writer verification only" }, null, 2));
  const run = Bun.spawnSync(["python3", join(root, "pipeline/scripts/cha01_drawing_ifc.py"), "--input", formal, "--manifest", manifest, "--approval", approvalPath, "--output", output, "--report", report, "--apply"], { cwd: root, stdout: "pipe", stderr: "pipe" });
  if (run.exitCode !== 0) throw new Error(run.stderr.toString());
  expect(run.exitCode).toBe(0);
  const result = JSON.parse(readFileSync(report, "utf8"));
  expect(result).toMatchObject({ pass: true, formal_ifc_bytes_unchanged: true, representations: { plan: "Cha01Plan", front: "Cha01Front", side: "Cha01Side" }, representation_path_counts: { plan: 1, front: 6, side: 4 }, representation_geometry_source: "proxy_paths_mm derived from the isolated representative IFC Body", official_cad_geometry_included: false, source_property_set: "Pset_Cha01DrawingSource", source_kind: sourceKind, source_label_zh: sourceLabelZh });
  expect(result.source_document_associations).toContain("BAXTER-COLETTE-CHA01-OFFICIAL-PRODUCT-PAGE");
  expect(result.source_document_associations).toContain("BAXTER-COLETTE-CHA01-MEASUREMENT-SVG-ARCHIVE");
  expect(result.source_document_associations).toContain("BAXTER-COLETTE-CHA01-SOURCE-ACCESS-RECORD");
  expect(result.source_document_associations.length).toBe(11);
  expect(sha256(formal)).toBe(formalHash);
  rmSync(temporary, { recursive: true, force: true });
}, 30_000);
