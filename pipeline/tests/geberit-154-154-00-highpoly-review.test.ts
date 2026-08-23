import { expect, test } from "bun:test";
import { createHash } from "node:crypto";
import { existsSync, mkdtempSync, readFileSync, rmSync, writeFileSync } from "node:fs";
import { tmpdir } from "node:os";
import { join, resolve } from "node:path";

const root = resolve(import.meta.dir, "../..");
const product = join(root, "output/review/highpoly-types/geberit-154-154-00");
const formal = join(root, "2504 GBTB Yanlord Zhuhai.ifc");
const formalHash = "7a50b87e8f48a7c2bfbaa9f1dabdfca7155fa684b2325c66aa1ae15c8ab0a25c";
const sourceKind = "geometry_derived_simplified_proxy";
const sourceLabelZh = "基于原始高模几何生成的简化图纸表达";

function sha256(path: string) {
  return createHash("sha256").update(readFileSync(path)).digest("hex");
}

test("exact Geberit identity and vector evidence are archived without claiming CAD", () => {
  const access = JSON.parse(readFileSync(join(product, "official-source/source-access-record.json"), "utf8"));
  const revalidation = JSON.parse(
    readFileSync(join(product, "official-source/official-source-revalidation.json"), "utf8"),
  );
  expect(access.manufacturer).toBe("Geberit");
  expect(access.family).toBe("CleanLine shower channel installation set");
  expect(access.project_ifc_type_name).toBe("Geberit 154.154.00");
  expect(access.project_ifc_type_description).toBe(
    "Installation set for CleanLine shower channel, for screed height at inlet 90–220 mm",
  );
  expect(access.resolved_article_number).toBe("154.154.00.1");
  expect(access.official_product_cad.cad_drawings_field).toBe("undefined");
  expect(access.official_product_cad.native_dwg_http_statuses).toEqual({ A: 404, G: 404, L: 404, P: 404 });
  expect(access.official_product_cad.acquired).toBeFalse();
  expect(access.official_product_cad.local_cad_files).toEqual([]);
  expect(access.official_product_cad.official_vector_eps_archived).toBeTrue();
  expect(access.official_product_cad.official_vector_eps_used_as_cad_geometry).toBeFalse();
  expect(revalidation.pass).toBeTrue();
  expect(revalidation.article_number).toBe("154.154.00.1");
  expect(revalidation.product_page.http_status).toBe(200);
  expect(revalidation.product_page.downloaded_exact_article_with_cadDrawings_undefined).toBeTrue();
  expect(revalidation.product_page.archived_exact_article_with_cadDrawings_undefined).toBeTrue();
  expect([200, 401, 403]).toContain(revalidation.article_api.http_status);
  expect(revalidation.article_api.note).toContain("not used to infer CAD availability");
  expect(Object.fromEntries(
    Object.entries(revalidation.native_dwg_results).map(([code, item]: any) => [code, item.http_status]),
  )).toEqual({ A: 404, G: 404, L: 404, P: 404 });
  expect(Object.values(revalidation.native_dwg_results).every((item: any) => item.native_dwg_acquired === false)).toBeTrue();
  expect(Object.values(revalidation.official_eps_results).every((item: any) => item.pass)).toBeTrue();
  expect(revalidation.official_eps_used_as_cad_geometry).toBeFalse();
  expect(revalidation.third_party_cad_used).toBeFalse();
  expect(revalidation.drawing_geometry_source).toEqual({
    source_kind: sourceKind,
    source_label_zh: sourceLabelZh,
  });
  expect(access.dimension_cross_check.pass).toBeTrue();
  expect(access.dimension_cross_check.compared_dimensions.every((item: any) => item.pass)).toBeTrue();
  expect(access.drawing_geometry_source).toMatchObject({
    source_kind: sourceKind,
    source_label_zh: sourceLabelZh,
    official_cad_used: false,
    third_party_cad_used: false,
  });
  for (const source of access.official_identity_sources) {
    expect(sha256(join(root, source.local_path))).toBe(source.sha256);
    if (source.preview_path) expect(sha256(join(root, source.preview_path))).toBe(source.preview_sha256);
  }
});

test("single representative has three black geometry-derived views and zero blue CAD paths", () => {
  const manifest = JSON.parse(readFileSync(join(product, "manifest.json"), "utf8"));
  const candidate = JSON.parse(readFileSync(join(product, "candidate-representations.json"), "utf8"));
  const profile = JSON.parse(readFileSync(join(product, "profile.json"), "utf8"));
  const inventory = JSON.parse(readFileSync(join(root, "pipeline/decisions/highpoly-product-inventory.json"), "utf8"));
  const inventoryEntry = inventory.products.find((entry: any) => entry.type_name === "Geberit 154.154.00");
  expect(inventoryEntry.status).toBe("review_ready_pending_approval");
  expect(inventory.review_ready_pending_approval_product_count).toBe(38);
  expect(profile.profiles["geberit-154-154-00"].representative_global_id).toBe("14EazrLgP8whYZWY_yCuKy");
  expect(manifest.geometry_product_count).toBe(1);
  expect(manifest.whole_model_render).toBeFalse();
  expect(manifest.representative_global_id).toBe("14EazrLgP8whYZWY_yCuKy");
  expect(manifest.drawing_source).toMatchObject({ source_kind: sourceKind, official_cad_used: false });
  expect(candidate.article_number).toBe("154.154.00.1");
  expect(candidate.review_status).toBe("visual_review_pending");
  expect(candidate.formal_ifc_write_allowed).toBeFalse();
  expect(Object.fromEntries(Object.entries(candidate.views).map(([view, value]: any) => [view, value.proxy_paths_mm.length]))).toEqual({
    plan: 1,
    front: 3,
    side: 1,
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

test("project drawings retain context and distinguish project pipe blue from product CAD", () => {
  const context = JSON.parse(readFileSync(join(product, "project-context-manifest.json"), "utf8"));
  expect(context.source_kind).toBe(sourceKind);
  expect(context.official_cad_used).toBeFalse();
  expect(context.project_context_retained).toBeTrue();
  expect(context.walls_and_surrounding_project_elements_retained).toBeTrue();
  expect(context.overlay_top_layer_with_white_mask).toBeTrue();
  expect(context.blue_line_present).toBeFalse();
  expect(context.review_annotation_suppression).toMatchObject({
    plan: ["p202-location-marker", "official-elevation-anchor"],
    walls_furniture_and_ifc_geometry_removed: false,
  });
  expect(context.project_projection_note.geometry_stretched).toBeFalse();
  for (const view of context.views) {
    expect(view.overlay.fit.uniform_scale_preserved).toBeTrue();
    expect(view.overlay.fit.scale_svg_units_per_mm).toBe(0.02);
    expect(view.overlay.fit.transformation_mode).toBe("axis_swap_rigid_reflection_and_translation_only");
    expect(Math.max(...view.overlay.fit.bbox_absolute_delta_svg_units)).toBeLessThanOrEqual(
      view.overlay.fit.bbox_tolerance_svg_units,
    );
    expect(sha256(join(root, view.review_preview))).toBe(view.review_preview_sha256);
    const svg = readFileSync(join(root, view.output), "utf8");
    expect(svg).toContain('class="geometry-derived-simplified-proxy project-context-overlay"');
    expect(svg).toContain('class="geometry-derived-proxy-mask"');
    expect(svg).toContain('class="geometry-derived-proxy"');
    expect(svg.indexOf('class="geometry-derived-proxy-mask"')).toBeLessThan(svg.indexOf('class="geometry-derived-proxy"'));
    expect(svg).not.toContain('class="official-reference native-dwg"');
  }
});

test("Bonsai evidence contains actual IFC Body camera renders", () => {
  const evidence = JSON.parse(readFileSync(join(product, "bonsai-review-manifest.json"), "utf8"));
  expect(evidence.mode).toBe("actual_bonsai_ifc_body_camera_render");
  expect(evidence.representative_global_id).toBe("14EazrLgP8whYZWY_yCuKy");
  expect(evidence.geometry_product_count).toBe(1);
  expect(evidence.whole_model_render).toBeFalse();
  expect(evidence.formal_ifc_bytes_unchanged).toBeTrue();
  expect(evidence.bonsai_session.saved_active_representation).toBe("Body");
  expect(evidence.bonsai_session.ifc_context_identifier).toBe("Body");
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

test("pending review leaves the formal IFC byte-identical and creates no derived IFC", () => {
  const candidate = JSON.parse(readFileSync(join(product, "candidate-representations.json"), "utf8"));
  expect(candidate.formal_ifc_write_allowed).toBeFalse();
  expect(candidate.review_status).toBe("visual_review_pending");
  expect(existsSync(join(product, "Geberit-154-154-00-1-derived-drawing.ifc"))).toBeFalse();
  expect(sha256(formal)).toBe(formalHash);
});

test("writer rejects missing apply and the pending human approval record", () => {
  const temporary = mkdtempSync(join(tmpdir(), "geberit-154-154-00-gate-"));
  const output = join(temporary, "forbidden.ifc");
  const script = join(root, "pipeline/scripts/geberit_154_154_00_drawing_ifc.py");
  const approval = JSON.parse(
    readFileSync(join(root, "pipeline/decisions/geberit-154-154-00-drawing-approval.json"), "utf8"),
  );
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

  const pendingApproval = Bun.spawnSync(
    ["python3", script, "--input", formal, "--output", output, "--apply"],
    { cwd: root, stderr: "pipe" },
  );
  expect(pendingApproval.exitCode).not.toBe(0);
  expect(pendingApproval.stderr.toString()).toContain("approval gate rejected IFC write");
  expect(existsSync(output)).toBeFalse();
  expect(sha256(formal)).toBe(formalHash);
  rmSync(temporary, { recursive: true, force: true });
});

test("scoped temporary approval writes verified representations and source associations", () => {
  const manifest = join(product, "manifest.json");
  const temporary = mkdtempSync(join(tmpdir(), "geberit-154-154-00-approved-"));
  const approvalPath = join(temporary, "approval.json");
  const output = join(temporary, "derived.ifc");
  const report = join(temporary, "report.json");
  writeFileSync(approvalPath, JSON.stringify({
    schema_version: 1,
    profile_key: "geberit-154-154-00",
    ifc_type_name: "Geberit 154.154.00",
    candidate_manifest_sha256: sha256(manifest),
    status: "approved",
    reviewer: "automated gate fixture",
    review_date: "2026-08-23",
    approved_views: ["plan", "front", "side"],
    derived_ifc_write_allowed: true,
    formal_authoritative_ifc_write_allowed: false,
    scope: "exact manufacturer article identity and standard dimensions; not a project shop drawing and not official CAD geometry",
    approval_evidence: "temporary automated writer verification only",
  }, null, 2));
  const run = Bun.spawnSync([
    "python3",
    join(root, "pipeline/scripts/geberit_154_154_00_drawing_ifc.py"),
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
  expect(result.representations).toEqual({
    plan: "Geberit154154Plan",
    front: "Geberit154154Front",
    side: "Geberit154154Side",
  });
  expect(result.representation_path_counts).toEqual({ plan: 1, front: 3, side: 1 });
  expect(result.source_kind).toBe(sourceKind);
  expect(result.source_label_zh).toBe(sourceLabelZh);
  expect(result.official_cad_geometry_included).toBeFalse();
  expect(result.source_document_associations).toEqual([
    "GEBERIT-154-154-00-1-OFFICIAL-PRODUCT-PAGE",
    "GEBERIT-154-154-00-1-PRODUCT-PAGE-ARCHIVE",
    "GEBERIT-154-154-00-1-SOURCE-ACCESS-RECORD",
    "GEBERIT-154-154-00-1-VECTOR-EPS-FRONT",
    "GEBERIT-154-154-00-1-VECTOR-EPS-PERSPECTIVE",
    "GEBERIT-154-154-00-1-VECTOR-EPS-TOP",
  ]);
  expect(sha256(formal)).toBe(formalHash);
  rmSync(temporary, { recursive: true, force: true });
}, 30_000);
