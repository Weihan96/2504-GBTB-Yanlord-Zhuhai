import { expect, test } from "bun:test";
import { createHash } from "node:crypto";
import { existsSync, mkdtempSync, readFileSync, rmSync, writeFileSync } from "node:fs";
import { tmpdir } from "node:os";
import { join, resolve } from "node:path";

const root = resolve(import.meta.dir, "../..");
const product = join(root, "output/review/highpoly-types/bst02");
const formal = join(root, "2504 GBTB Yanlord Zhuhai.ifc");
const formalHash = "7a50b87e8f48a7c2bfbaa9f1dabdfca7155fa684b2325c66aa1ae15c8ab0a25c";
const sourceKind = "geometry_derived_simplified_proxy";
const sourceLabelZh = "基于原始高模几何生成的简化图纸表达";
const scope = "exact Baxter Beside 55 x 58 x 35 cm manufacturer identity and public measurement evidence; native 2D/3D/BIM requires login and was not acquired; public PDF/SVG are not representation CAD geometry or a project shop drawing; the 100 mm project Body height discrepancy is retained";

function sha256(path: string) {
  return createHash("sha256").update(readFileSync(path)).digest("hex");
}

test("exact Baxter Beside identity is archived while native CAD remains unavailable", () => {
  const access = JSON.parse(readFileSync(join(product, "official-source/source-access-record.json"), "utf8"));
  expect(access).toMatchObject({
    manufacturer: "Baxter",
    family: "Beside",
    designer: "Studiopepe",
    project_ifc_type_name: "BST02",
    project_ifc_type_description: "Beside 55x58xh35",
    resolved_variant: "Beside bedside table, 55 x 58 x 35 cm",
    scope,
    pass: true,
  });
  expect(access.official_product_cad).toMatchObject({
    authentication_required: true,
    acquired: false,
    exact_project_configuration_match: false,
    local_cad_files: [],
    official_vector_pdf_archived: true,
    official_measurement_svg_archived: true,
    official_vector_evidence_used_as_cad_geometry: false,
  });
  expect(access.dimension_cross_check).toMatchObject({
    official_variant_overall_mm: [550, 580, 350],
    project_ifc_body_local_xyz_mm: [560, 577.889252, 250],
    absolute_delta_mm: [10, 2.110748, 100],
    pass_by_axis: [true, true, false],
    pass: false,
    review_required: true,
  });
  const measurement = access.official_identity_sources.find((source: any) => source.kind === "manufacturer_public_measurement_svg");
  expect(measurement.vector_object_audit).toEqual({ path_count: 34, use_count: 45, group_count: 41 });
  for (const source of access.official_identity_sources) {
    expect(sha256(join(root, source.local_path))).toBe(source.sha256);
    if (source.preview_path) expect(sha256(join(root, source.preview_path))).toBe(source.preview_sha256);
  }
  expect(access.drawing_geometry_source).toMatchObject({
    source_kind: sourceKind,
    source_label_zh: sourceLabelZh,
    official_cad_used: false,
    official_public_vector_evidence_used_as_cad_geometry: false,
    third_party_cad_used: false,
  });
});

test("single BST02 representative has three geometry-derived views and zero blue CAD paths", () => {
  const manifest = JSON.parse(readFileSync(join(product, "manifest.json"), "utf8"));
  const candidate = JSON.parse(readFileSync(join(product, "candidate-representations.json"), "utf8"));
  const profile = JSON.parse(readFileSync(join(product, "profile.json"), "utf8"));
  const inventory = JSON.parse(readFileSync(join(root, "pipeline/decisions/highpoly-product-inventory.json"), "utf8"));
  const entry = inventory.products.find((item: any) => item.type_name === "BST02");
  expect(entry.status).toBe("review_ready_pending_approval");
  expect(inventory.review_ready_pending_approval_product_count).toBe(38);
  expect(manifest.representative_global_id).toBe("3hA0vKpcn44u4Tsx4tqiUz");
  expect(manifest.geometry_product_count).toBe(1);
  expect(manifest.whole_model_render).toBeFalse();
  expect(manifest.bounds_mm.size).toEqual([560, 577.889252, 250]);
  expect(manifest.dimension_cross_check.pass).toBeFalse();
  expect(manifest.dimension_cross_check.review_required).toBeTrue();
  expect(profile.profiles.bst02).toMatchObject({
    ifc_type_name: "BST02",
    representative_global_id: "3hA0vKpcn44u4Tsx4tqiUz",
    drawing_source: { source_kind: sourceKind, official_cad_used: false, third_party_cad_used: false },
  });
  expect(candidate.article_number).toBe("Baxter Beside 55x58xh35 / BST02");
  expect(candidate.review_status).toBe("visual_review_pending");
  expect(candidate.formal_ifc_write_allowed).toBeFalse();
  expect(Object.fromEntries(Object.entries(candidate.views).map(([view, value]: any) => [view, value.proxy_paths_mm.length]))).toEqual({ plan: 6, front: 1, side: 1 });
  for (const view of ["plan", "front", "side"]) {
    expect(candidate.views[view].source_kind).toBe(sourceKind);
    expect(candidate.views[view].official_cad_paths_mm).toEqual([]);
    const svg = readFileSync(join(product, `${view}.svg`), "utf8");
    expect(svg).toContain(`data-source-kind="${sourceKind}"`);
    expect(svg).not.toContain('class="official-reference native-dwg"');
    expect(svg).not.toContain("#1677c8");
  }
});

test("BST02 project context uses the real plan and exact R09 side projection", () => {
  const context = JSON.parse(readFileSync(join(product, "project-context-manifest.json"), "utf8"));
  expect(context.views.map((view: any) => view.view)).toEqual(["plan", "side"]);
  expect(context.context_view_scope.excluded).toEqual(["front"]);
  expect(context.project_context_retained).toBeTrue();
  expect(context.walls_and_surrounding_project_elements_retained).toBeTrue();
  expect(context.overlay_top_layer_with_white_mask).toBeTrue();
  expect(context.blue_line_present).toBeFalse();
  expect(context.project_projection_note.geometry_stretched).toBeFalse();
  expect(context.project_projection_note.official_height_discrepancy_mm).toBe(100);
  const [plan, side] = context.views;
  expect(plan.overlay.fit.rotate_quarter_turns).toBe(1);
  expect(Math.max(...plan.overlay.fit.bbox_absolute_delta_svg_units)).toBeLessThanOrEqual(plan.overlay.fit.bbox_tolerance_svg_units);
  expect(side.overlay.fit.rotate_quarter_turns).toBe(0);
  expect(side.overlay.fit.bbox_absolute_delta_svg_units).toEqual([0, 0]);
  for (const view of context.views) {
    expect(view.overlay.fit.uniform_scale_preserved).toBeTrue();
    expect(view.overlay.fit.scale_svg_units_per_mm).toBe(0.02);
    expect(sha256(join(root, view.review_preview))).toBe(view.review_preview_sha256);
    const full = readFileSync(join(root, view.output), "utf8");
    expect(full).toContain('class="geometry-derived-proxy-mask"');
    expect(full).toContain('class="geometry-derived-proxy"');
    expect(full.indexOf('class="geometry-derived-proxy-mask"')).toBeLessThan(full.indexOf('class="geometry-derived-proxy"'));
  }
});

test("BST02 Bonsai evidence renders the actual isolated MODEL_VIEW Body with four cameras", () => {
  const evidence = JSON.parse(readFileSync(join(product, "bonsai-review-manifest.json"), "utf8"));
  expect(evidence.mode).toBe("actual_bonsai_ifc_body_camera_render");
  expect(evidence.representative_global_id).toBe("3hA0vKpcn44u4Tsx4tqiUz");
  expect(evidence.geometry_product_count).toBe(1);
  expect(evidence.whole_model_render).toBeFalse();
  expect(evidence.formal_ifc_bytes_unchanged).toBeTrue();
  expect(evidence.bonsai_session).toMatchObject({ saved_active_representation: "Body", ifc_context_identifier: "Body", ifc_target_view: "MODEL_VIEW", saved_camera_count: 4 });
  expect(evidence.renders.map((render: any) => render.view)).toEqual(["plan", "front", "side", "iso"]);
  for (const render of evidence.renders) {
    expect(render.camera_type).toBe("ORTHO");
    expect(render.render_operation).toBe("bpy.ops.render.render(write_still=True)");
    expect(sha256(join(root, render.path))).toBe(render.sha256);
  }
  expect(sha256(join(root, evidence.bonsai_session.path))).toBe(evidence.bonsai_session.sha256);
  expect(sha256(join(root, evidence.isolated_ifc))).toBe(evidence.isolated_ifc_sha256);
});

test("pending BST02 review blocks every IFC write and preserves the formal IFC", () => {
  const temporary = mkdtempSync(join(tmpdir(), "bst02-gate-"));
  const output = join(temporary, "forbidden.ifc");
  const script = join(root, "pipeline/scripts/bst02_drawing_ifc.py");
  const approval = JSON.parse(readFileSync(join(root, "pipeline/decisions/bst02-drawing-approval.json"), "utf8"));
  const manifest = join(product, "manifest.json");
  expect(approval.status).toBe("pending");
  expect(approval.derived_ifc_write_allowed).toBeFalse();
  expect(approval.candidate_manifest_sha256).toBe(sha256(manifest));
  expect(existsSync(join(product, "Baxter-Beside-BST02-derived-drawing.ifc"))).toBeFalse();

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

test("scoped temporary approval writes verified BST02 representations and source documents", () => {
  const manifest = join(product, "manifest.json");
  const temporary = mkdtempSync(join(tmpdir(), "bst02-approved-"));
  const approvalPath = join(temporary, "approval.json");
  const output = join(temporary, "derived.ifc");
  const report = join(temporary, "report.json");
  writeFileSync(approvalPath, JSON.stringify({
    schema_version: 1,
    profile_key: "bst02",
    ifc_type_name: "BST02",
    candidate_manifest_sha256: sha256(manifest),
    status: "approved",
    reviewer: "automated gate fixture",
    review_date: "2026-08-23",
    approved_views: ["plan", "front", "side"],
    derived_ifc_write_allowed: true,
    formal_authoritative_ifc_write_allowed: false,
    scope,
    approval_evidence: "temporary automated writer verification only; 100 mm height discrepancy explicitly retained",
  }, null, 2));
  const run = Bun.spawnSync([
    "python3", join(root, "pipeline/scripts/bst02_drawing_ifc.py"),
    "--input", formal, "--manifest", manifest, "--approval", approvalPath,
    "--output", output, "--report", report, "--apply",
  ], { cwd: root, stdout: "pipe", stderr: "pipe" });
  expect(run.exitCode).toBe(0);
  if (run.exitCode !== 0) throw new Error(run.stderr.toString());
  const result = JSON.parse(readFileSync(report, "utf8"));
  expect(result.pass).toBeTrue();
  expect(result.formal_ifc_bytes_unchanged).toBeTrue();
  expect(result.representations).toEqual({ plan: "Bst02Plan", front: "Bst02Front", side: "Bst02Side" });
  expect(result.representation_path_counts).toEqual({ plan: 6, front: 1, side: 1 });
  expect(result.source_kind).toBe(sourceKind);
  expect(result.official_cad_geometry_included).toBeFalse();
  expect(result.source_document_associations).toEqual([
    "BAXTER-BESIDE-BST02-OFFICIAL-MEASUREMENT-SVG",
    "BAXTER-BESIDE-BST02-OFFICIAL-MEASUREMENT-SVG-ARCHIVE",
    "BAXTER-BESIDE-BST02-OFFICIAL-PRODUCT-PAGE",
    "BAXTER-BESIDE-BST02-OFFICIAL-PRODUCT-PAGE-ARCHIVE",
    "BAXTER-BESIDE-BST02-OFFICIAL-TECHNICAL-SHEET",
    "BAXTER-BESIDE-BST02-OFFICIAL-TECHNICAL-SHEET-ARCHIVE",
    "BAXTER-BESIDE-BST02-SOURCE-ACCESS-RECORD",
  ]);
  expect(sha256(formal)).toBe(formalHash);
  rmSync(temporary, { recursive: true, force: true });
}, 30_000);
