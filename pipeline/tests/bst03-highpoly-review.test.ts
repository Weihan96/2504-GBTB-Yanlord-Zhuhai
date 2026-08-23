import { expect, test } from "bun:test";
import { createHash } from "node:crypto";
import { existsSync, mkdtempSync, readFileSync, rmSync, writeFileSync } from "node:fs";
import { tmpdir } from "node:os";
import { join, resolve } from "node:path";

const root = resolve(import.meta.dir, "../..");
const product = join(root, "output/review/highpoly-types/bst03");
const formal = join(root, "2504 GBTB Yanlord Zhuhai.ifc");
const formalHash = "7a50b87e8f48a7c2bfbaa9f1dabdfca7155fa684b2325c66aa1ae15c8ab0a25c";
const sourceKind = "geometry_derived_simplified_proxy";
const sourceLabelZh = "基于原始高模几何生成的简化图纸表达";

function sha256(path: string) {
  return createHash("sha256").update(readFileSync(path)).digest("hex");
}

test("exact Baxter Stone left-drawer identity is archived without claiming acquired CAD", () => {
  const access = JSON.parse(readFileSync(join(product, "official-source/source-access-record.json"), "utf8"));
  expect(access.manufacturer).toBe("Baxter");
  expect(access.family).toBe("Stone");
  expect(access.designer).toBe("Federico Peri");
  expect(access.project_ifc_type_name).toBe("BST03");
  expect(access.project_ifc_type_description).toBe("Stone Bedside Table with Drawer D45");
  expect(access.resolved_variant).toBe("freestanding bedside table with L drawer, 45 x 45 x 46 cm");
  expect(access.project_identity_evidence).toMatchObject({
    register_status: "confirmed",
    register_confidence: 1,
    confirmed_variant: "freestanding_bedside_table_with_L_drawer",
  });
  expect(access.official_product_cad.authentication_required).toBeTrue();
  expect(access.official_product_cad.acquired).toBeFalse();
  expect(access.official_product_cad.local_cad_files).toEqual([]);
  expect(access.official_product_cad.official_vector_pdf_archived).toBeTrue();
  expect(access.official_product_cad.official_vector_pdf_used_as_cad_geometry).toBeFalse();
  expect(access.dimension_cross_check).toMatchObject({
    official_variant_overall_mm: [450, 450, 460],
    pass_by_axis: [true, true, true],
    pass: true,
  });
  const pdf = access.official_identity_sources.find((source: any) => source.kind === "manufacturer_vector_technical_sheet");
  expect(pdf.page).toBe(10);
  expect(pdf.vector_object_audit).toMatchObject({ line_count: 810, curve_count: 75, rect_count: 63, char_count: 1860 });
  expect(pdf.vector_object_audit.required_text_found).toContain("Freestanding night table with r/l drawer.");
  for (const source of access.official_identity_sources) {
    expect(sha256(join(root, source.local_path))).toBe(source.sha256);
    if (source.preview_path) expect(sha256(join(root, source.preview_path))).toBe(source.preview_sha256);
  }
  expect(access.drawing_geometry_source).toMatchObject({
    source_kind: sourceKind,
    source_label_zh: sourceLabelZh,
    official_cad_used: false,
    third_party_cad_used: false,
  });
});

test("single BST03 representative has three clean geometry-derived views and zero blue CAD paths", () => {
  const manifest = JSON.parse(readFileSync(join(product, "manifest.json"), "utf8"));
  const candidate = JSON.parse(readFileSync(join(product, "candidate-representations.json"), "utf8"));
  const inventory = JSON.parse(readFileSync(join(root, "pipeline/decisions/highpoly-product-inventory.json"), "utf8"));
  const inventoryEntry = inventory.products.find((entry: any) => entry.type_name === "BST03");
  expect(inventoryEntry.status).toBe("review_ready_pending_approval");
  expect(inventory.review_ready_pending_approval_product_count).toBe(38);
  expect(manifest.representative_global_id).toBe("1HZoxe$df4cBb4UXXH5J2S");
  expect(manifest.geometry_product_count).toBe(1);
  expect(manifest.whole_model_render).toBeFalse();
  expect(manifest.bounds_mm.size).toEqual([462.460968, 462.47998, 466.218994]);
  expect(candidate.article_number).toBe("Baxter Stone L drawer 45 / BST03");
  expect(candidate.review_status).toBe("visual_review_pending");
  expect(candidate.formal_ifc_write_allowed).toBeFalse();
  expect(Object.fromEntries(Object.entries(candidate.views).map(([view, value]: any) => [view, value.proxy_paths_mm.length]))).toEqual({
    plan: 1,
    front: 4,
    side: 2,
  });
  for (const view of ["plan", "front", "side"]) {
    expect(candidate.views[view].source_kind).toBe(sourceKind);
    expect(candidate.views[view].official_cad_paths_mm).toEqual([]);
    const svg = readFileSync(join(product, `${view}.svg`), "utf8");
    expect(svg).toContain(`data-source-kind="${sourceKind}"`);
    expect(svg).not.toContain('class="official-reference native-dwg"');
    expect(svg).not.toContain("#1677c8");
  }
});

test("BST03 project plan retains walls and furniture without inventing project elevations", () => {
  const context = JSON.parse(readFileSync(join(product, "project-context-manifest.json"), "utf8"));
  expect(context.views.map((view: any) => view.view)).toEqual(["plan"]);
  expect(context.source_kind).toBe(sourceKind);
  expect(context.official_cad_used).toBeFalse();
  expect(context.project_context_retained).toBeTrue();
  expect(context.walls_and_surrounding_project_elements_retained).toBeTrue();
  expect(context.overlay_top_layer_with_white_mask).toBeTrue();
  expect(context.blue_line_present).toBeFalse();
  expect(context.project_projection_note.geometry_stretched).toBeFalse();
  expect(context.review_annotation_suppression.walls_furniture_and_ifc_geometry_removed).toBeFalse();
  const view = context.views[0];
  expect(view.overlay.fit.rotate_quarter_turns).toBe(1);
  expect(view.overlay.fit.uniform_scale_preserved).toBeTrue();
  expect(view.overlay.fit.scale_svg_units_per_mm).toBe(0.02);
  expect(Math.max(...view.overlay.fit.bbox_absolute_delta_svg_units)).toBeLessThanOrEqual(
    view.overlay.fit.bbox_tolerance_svg_units,
  );
  expect(sha256(join(root, view.review_preview))).toBe(view.review_preview_sha256);
  const full = readFileSync(join(root, view.output), "utf8");
  const crop = readFileSync(join(root, view.review_crop), "utf8");
  expect(full).toContain('class="geometry-derived-simplified-proxy project-context-overlay"');
  expect(full).toContain('class="geometry-derived-proxy-mask"');
  expect(full).toContain('class="geometry-derived-proxy"');
  expect(full.indexOf('class="geometry-derived-proxy-mask"')).toBeLessThan(full.indexOf('class="geometry-derived-proxy"'));
  expect(full).not.toContain('class="official-reference native-dwg"');
  expect(crop).toContain('class="review-paper-background"');
  expect(crop).not.toContain('class="official-elevation-anchor"');
});

test("BST03 Bonsai evidence uses the actual MODEL_VIEW IFC Body", () => {
  const evidence = JSON.parse(readFileSync(join(product, "bonsai-review-manifest.json"), "utf8"));
  expect(evidence.mode).toBe("actual_bonsai_ifc_body_camera_render");
  expect(evidence.representative_global_id).toBe("1HZoxe$df4cBb4UXXH5J2S");
  expect(evidence.geometry_product_count).toBe(1);
  expect(evidence.whole_model_render).toBeFalse();
  expect(evidence.formal_ifc_bytes_unchanged).toBeTrue();
  expect(evidence.bonsai_session.saved_active_representation).toBe("Body");
  expect(evidence.bonsai_session.ifc_context_identifier).toBe("Body");
  expect(evidence.bonsai_session.ifc_target_view).toBe("MODEL_VIEW");
  expect(evidence.bonsai_session.saved_camera_count).toBe(4);
  expect(evidence.renders.map((render: any) => render.view)).toEqual(["plan", "front", "side", "iso"]);
  for (const render of evidence.renders) {
    expect(render.camera_type).toBe("ORTHO");
    expect(render.render_operation).toBe("bpy.ops.render.render(write_still=True)");
    expect(render.camera_ortho_scale_m + 1e-8).toBeGreaterThanOrEqual(render.projected_width_m * render.framing_margin_factor);
    expect(render.camera_ortho_scale_m + 1e-8).toBeGreaterThanOrEqual(
      render.projected_height_m * render.image_aspect * render.framing_margin_factor,
    );
    expect(sha256(join(root, render.path))).toBe(render.sha256);
  }
  expect(sha256(join(root, evidence.bonsai_session.path))).toBe(evidence.bonsai_session.sha256);
  expect(sha256(join(root, evidence.isolated_ifc))).toBe(evidence.isolated_ifc_sha256);
});

test("pending BST03 review leaves formal IFC unchanged and creates no derived IFC", () => {
  const candidate = JSON.parse(readFileSync(join(product, "candidate-representations.json"), "utf8"));
  expect(candidate.formal_ifc_write_allowed).toBeFalse();
  expect(candidate.review_status).toBe("visual_review_pending");
  expect(existsSync(join(product, "Baxter-Stone-BST03-derived-drawing.ifc"))).toBeFalse();
  expect(existsSync(join(product, "Baxter-Stone-BST03-bonsai-review.blend1"))).toBeFalse();
  expect(sha256(formal)).toBe(formalHash);
});

test("BST03 writer rejects missing apply and pending human approval", () => {
  const temporary = mkdtempSync(join(tmpdir(), "bst03-gate-"));
  const output = join(temporary, "forbidden.ifc");
  const script = join(root, "pipeline/scripts/bst03_drawing_ifc.py");
  const approval = JSON.parse(readFileSync(join(root, "pipeline/decisions/bst03-drawing-approval.json"), "utf8"));
  const manifest = join(product, "manifest.json");
  expect(approval.status).toBe("pending");
  expect(approval.derived_ifc_write_allowed).toBeFalse();
  expect(approval.formal_authoritative_ifc_write_allowed).toBeFalse();
  expect(approval.candidate_manifest_sha256).toBe(sha256(manifest));

  const withoutApply = Bun.spawnSync(["python3", script, "--input", formal, "--output", output], {
    cwd: root,
    stderr: "pipe",
  });
  expect(withoutApply.exitCode).not.toBe(0);
  expect(withoutApply.stderr.toString()).toContain("IFC write requires the explicit --apply flag");
  expect(existsSync(output)).toBeFalse();

  const pending = Bun.spawnSync(["python3", script, "--input", formal, "--output", output, "--apply"], {
    cwd: root,
    stderr: "pipe",
  });
  expect(pending.exitCode).not.toBe(0);
  expect(pending.stderr.toString()).toContain("approval gate rejected IFC write");
  expect(existsSync(output)).toBeFalse();
  expect(sha256(formal)).toBe(formalHash);
  rmSync(temporary, { recursive: true, force: true });
});

test("scoped temporary approval writes verified BST03 representations and source documents", () => {
  const manifest = join(product, "manifest.json");
  const temporary = mkdtempSync(join(tmpdir(), "bst03-approved-"));
  const approvalPath = join(temporary, "approval.json");
  const output = join(temporary, "derived.ifc");
  const report = join(temporary, "report.json");
  writeFileSync(approvalPath, JSON.stringify({
    schema_version: 1,
    profile_key: "bst03",
    ifc_type_name: "BST03",
    candidate_manifest_sha256: sha256(manifest),
    status: "approved",
    reviewer: "automated gate fixture",
    review_date: "2026-08-23",
    approved_views: ["plan", "front", "side"],
    derived_ifc_write_allowed: true,
    formal_authoritative_ifc_write_allowed: false,
    scope: "exact manufacturer family, left-drawer variant and nominal dimensions; not a project shop drawing and not official CAD geometry",
    approval_evidence: "temporary automated writer verification only",
  }, null, 2));
  const run = Bun.spawnSync([
    "python3",
    join(root, "pipeline/scripts/bst03_drawing_ifc.py"),
    "--input", formal,
    "--manifest", manifest,
    "--approval", approvalPath,
    "--output", output,
    "--report", report,
    "--apply",
  ], { cwd: root, stdout: "pipe", stderr: "pipe" });
  expect(run.exitCode).toBe(0);
  if (run.exitCode !== 0) throw new Error(run.stderr.toString());
  const result = JSON.parse(readFileSync(report, "utf8"));
  expect(result.pass).toBeTrue();
  expect(result.formal_ifc_bytes_unchanged).toBeTrue();
  expect(result.representations).toEqual({ plan: "Bst03Plan", front: "Bst03Front", side: "Bst03Side" });
  expect(result.representation_path_counts).toEqual({ plan: 1, front: 4, side: 2 });
  expect(result.source_kind).toBe(sourceKind);
  expect(result.source_label_zh).toBe(sourceLabelZh);
  expect(result.official_cad_geometry_included).toBeFalse();
  expect(result.source_document_associations).toEqual([
    "BAXTER-STONE-BST03-OFFICIAL-PRODUCT-PAGE",
    "BAXTER-STONE-BST03-OFFICIAL-PRODUCT-PAGE-ARCHIVE",
    "BAXTER-STONE-BST03-OFFICIAL-TECHNICAL-SHEET",
    "BAXTER-STONE-BST03-OFFICIAL-TECHNICAL-SHEET-ARCHIVE",
    "BAXTER-STONE-BST03-SOURCE-ACCESS-RECORD",
  ]);
  expect(sha256(formal)).toBe(formalHash);
  rmSync(temporary, { recursive: true, force: true });
}, 30_000);
