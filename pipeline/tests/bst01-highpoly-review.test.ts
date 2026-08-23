import { expect, test } from "bun:test";
import { createHash } from "node:crypto";
import { existsSync, mkdtempSync, readFileSync, rmSync, writeFileSync } from "node:fs";
import { tmpdir } from "node:os";
import { join, resolve } from "node:path";

const root = resolve(import.meta.dir, "../..");
const product = join(root, "output/review/highpoly-types/bst01");
const formal = join(root, "2504 GBTB Yanlord Zhuhai.ifc");
const formalHash = "7a50b87e8f48a7c2bfbaa9f1dabdfca7155fa684b2325c66aa1ae15c8ab0a25c";
const sourceKind = "geometry_derived_simplified_proxy";
const sourceLabelZh = "基于原始高模几何生成的简化图纸表达";
const scope = "exact Baxter Ninfea bedside-table family and nominal 42 x 42 x 45 cm dimensions; official right- and left-opening public vectors both archived because the project type does not record opening side; native 2D/3D/BIM requires login and was not acquired; public PDF/SVG are not representation CAD geometry or a project shop drawing";

function sha256(path: string) {
  return createHash("sha256").update(readFileSync(path)).digest("hex");
}

test("exact Baxter Ninfea family and both opening variants are archived without selecting one", () => {
  const access = JSON.parse(readFileSync(join(product, "official-source/source-access-record.json"), "utf8"));
  expect(access).toMatchObject({
    manufacturer: "Baxter",
    family: "Ninfea",
    designer: "Pietro Russo",
    project_ifc_type_name: "BST01",
    project_ifc_type_description: "Nifea Comodino diam42xh45",
    resolved_variant: "Ninfea bedside-table family, 42 x 42 x 45 cm; opening side unresolved",
    scope,
    pass: true,
  });
  expect(access.official_product_cad).toMatchObject({
    authentication_required: true,
    acquired: false,
    exact_project_configuration_match: false,
    project_opening_side_resolved: false,
    local_cad_files: [],
    official_vector_pdf_archived: true,
    official_left_and_right_measurement_svgs_archived: true,
    official_vector_evidence_used_as_cad_geometry: false,
  });
  expect(access.dimension_cross_check).toMatchObject({
    official_variant_overall_mm: [420, 420, 450],
    project_ifc_body_local_xyz_mm: [399.030914, 399.031906, 450],
    absolute_delta_mm: [20.969086, 20.968094, 0],
    pass_by_axis: [true, true, true],
    pass: true,
  });
  const variants = access.official_identity_sources.filter((source: any) => source.kind.includes("measurement_svg"));
  expect(variants.map((source: any) => source.kind)).toEqual([
    "manufacturer_public_measurement_svg_right_opening",
    "manufacturer_public_measurement_svg_left_opening",
  ]);
  for (const variant of variants) expect(variant.vector_object_audit).toEqual({ path_count: 507, use_count: 66, group_count: 57 });
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

test("single BST01 representative has three geometry-derived views and zero blue CAD paths", () => {
  const manifest = JSON.parse(readFileSync(join(product, "manifest.json"), "utf8"));
  const candidate = JSON.parse(readFileSync(join(product, "candidate-representations.json"), "utf8"));
  const profile = JSON.parse(readFileSync(join(product, "profile.json"), "utf8"));
  const inventory = JSON.parse(readFileSync(join(root, "pipeline/decisions/highpoly-product-inventory.json"), "utf8"));
  const entry = inventory.products.find((item: any) => item.type_name === "BST01");
  expect(entry.status).toBe("review_ready_pending_approval");
  expect(inventory.review_ready_pending_approval_product_count).toBe(38);
  expect(manifest.representative_global_id).toBe("3eic1dzkn5heTIn4PhF37v");
  expect(manifest.geometry_product_count).toBe(1);
  expect(manifest.whole_model_render).toBeFalse();
  expect(manifest.bounds_mm.size).toEqual([399.030914, 399.031906, 450]);
  expect(manifest.dimension_cross_check.pass).toBeTrue();
  expect(profile.profiles.bst01).toMatchObject({
    ifc_type_name: "BST01",
    representative_global_id: "3eic1dzkn5heTIn4PhF37v",
    drawing_source: { source_kind: sourceKind, official_cad_used: false, third_party_cad_used: false },
  });
  expect(candidate.article_number).toBe("Baxter Ninfea bedside table opening side unresolved / BST01");
  expect(candidate.review_status).toBe("visual_review_pending");
  expect(candidate.formal_ifc_write_allowed).toBeFalse();
  expect(Object.fromEntries(Object.entries(candidate.views).map(([view, value]: any) => [view, value.proxy_paths_mm.length]))).toEqual({ plan: 21, front: 1, side: 5 });
  for (const view of ["plan", "front", "side"]) {
    expect(candidate.views[view].source_kind).toBe(sourceKind);
    expect(candidate.views[view].official_cad_paths_mm).toEqual([]);
    const svg = readFileSync(join(product, `${view}.svg`), "utf8");
    expect(svg).toContain(`data-source-kind="${sourceKind}"`);
    expect(svg).not.toContain('class="official-reference native-dwg"');
    expect(svg).not.toContain("#1677c8");
  }
});

test("BST01 project context mechanically maps plan, front and side without stretching", () => {
  const context = JSON.parse(readFileSync(join(product, "project-context-manifest.json"), "utf8"));
  expect(context.views.map((view: any) => view.view)).toEqual(["plan", "front", "side"]);
  expect(context.context_view_scope.excluded).toEqual([]);
  expect(context.project_context_retained).toBeTrue();
  expect(context.walls_and_surrounding_project_elements_retained).toBeTrue();
  expect(context.overlay_top_layer_with_white_mask).toBeTrue();
  expect(context.blue_line_present).toBeFalse();
  expect(context.project_projection_note.geometry_stretched).toBeFalse();
  expect(context.project_projection_note.opening_side_claimed).toBeFalse();
  const [plan, front, side] = context.views;
  expect(plan.overlay.fit.rotate_quarter_turns).toBe(1);
  expect(Math.max(...plan.overlay.fit.bbox_absolute_delta_svg_units)).toBeLessThanOrEqual(plan.overlay.fit.bbox_tolerance_svg_units);
  expect(front.overlay.fit.bbox_absolute_delta_svg_units).toEqual([0, 0]);
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

test("BST01 Bonsai evidence renders the actual isolated MODEL_VIEW Body with four cameras", () => {
  const evidence = JSON.parse(readFileSync(join(product, "bonsai-review-manifest.json"), "utf8"));
  expect(evidence.mode).toBe("actual_bonsai_ifc_body_camera_render");
  expect(evidence.representative_global_id).toBe("3eic1dzkn5heTIn4PhF37v");
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

test("pending BST01 review blocks every IFC write and preserves the formal IFC", () => {
  const temporary = mkdtempSync(join(tmpdir(), "bst01-gate-"));
  const output = join(temporary, "forbidden.ifc");
  const script = join(root, "pipeline/scripts/bst01_drawing_ifc.py");
  const approval = JSON.parse(readFileSync(join(root, "pipeline/decisions/bst01-drawing-approval.json"), "utf8"));
  const manifest = join(product, "manifest.json");
  expect(approval.status).toBe("pending");
  expect(approval.derived_ifc_write_allowed).toBeFalse();
  expect(approval.candidate_manifest_sha256).toBe(sha256(manifest));
  expect(existsSync(join(product, "Baxter-Ninfea-BST01-derived-drawing.ifc"))).toBeFalse();

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

test("scoped temporary approval writes verified BST01 representations and unresolved-side provenance", () => {
  const manifest = join(product, "manifest.json");
  const temporary = mkdtempSync(join(tmpdir(), "bst01-approved-"));
  const approvalPath = join(temporary, "approval.json");
  const output = join(temporary, "derived.ifc");
  const report = join(temporary, "report.json");
  writeFileSync(approvalPath, JSON.stringify({
    schema_version: 1,
    profile_key: "bst01",
    ifc_type_name: "BST01",
    candidate_manifest_sha256: sha256(manifest),
    status: "approved",
    reviewer: "automated gate fixture",
    review_date: "2026-08-23",
    approved_views: ["plan", "front", "side"],
    derived_ifc_write_allowed: true,
    formal_authoritative_ifc_write_allowed: false,
    scope,
    approval_evidence: "temporary automated writer verification only; project opening side remains unresolved",
  }, null, 2));
  const run = Bun.spawnSync([
    "python3", join(root, "pipeline/scripts/bst01_drawing_ifc.py"),
    "--input", formal, "--manifest", manifest, "--approval", approvalPath,
    "--output", output, "--report", report, "--apply",
  ], { cwd: root, stdout: "pipe", stderr: "pipe" });
  expect(run.exitCode).toBe(0);
  if (run.exitCode !== 0) throw new Error(run.stderr.toString());
  const result = JSON.parse(readFileSync(report, "utf8"));
  expect(result.pass).toBeTrue();
  expect(result.formal_ifc_bytes_unchanged).toBeTrue();
  expect(result.representations).toEqual({ plan: "Bst01Plan", front: "Bst01Front", side: "Bst01Side" });
  expect(result.representation_path_counts).toEqual({ plan: 21, front: 1, side: 5 });
  expect(result.source_kind).toBe(sourceKind);
  expect(result.official_cad_geometry_included).toBeFalse();
  expect(result.source_document_associations).toEqual([
    "BAXTER-NINFEA-BST01-OFFICIAL-LEFT-OPENING-SVG",
    "BAXTER-NINFEA-BST01-OFFICIAL-LEFT-OPENING-SVG-ARCHIVE",
    "BAXTER-NINFEA-BST01-OFFICIAL-PRODUCT-PAGE",
    "BAXTER-NINFEA-BST01-OFFICIAL-PRODUCT-PAGE-ARCHIVE",
    "BAXTER-NINFEA-BST01-OFFICIAL-RIGHT-OPENING-SVG",
    "BAXTER-NINFEA-BST01-OFFICIAL-RIGHT-OPENING-SVG-ARCHIVE",
    "BAXTER-NINFEA-BST01-OFFICIAL-TECHNICAL-SHEET",
    "BAXTER-NINFEA-BST01-OFFICIAL-TECHNICAL-SHEET-ARCHIVE",
    "BAXTER-NINFEA-BST01-SOURCE-ACCESS-RECORD",
  ]);
  const inspect = Bun.spawnSync(["python3", "-c", "import ifcopenshell,ifcopenshell.util.element,json,sys; m=ifcopenshell.open(sys.argv[1]); p=m.by_guid('3eic1dzkn5heTIn4PhF37v'); print(json.dumps(ifcopenshell.util.element.get_pset(p,'Pset_Bst01DrawingSource')))" , output], { cwd: root, stdout: "pipe", stderr: "pipe" });
  expect(inspect.exitCode).toBe(0);
  const pset = JSON.parse(inspect.stdout.toString());
  expect(pset).toMatchObject({
    SourceKind: sourceKind,
    OfficialCadUsed: "false",
    OfficialOpeningSideResolved: "false",
    OfficialRightOpeningVectorEvidence: "archived_not_representation_geometry",
    OfficialLeftOpeningVectorEvidence: "archived_not_representation_geometry",
    RepresentationGeometrySource: "proxy_paths_mm derived from the isolated representative IFC MODEL_VIEW Body",
  });
  expect(sha256(formal)).toBe(formalHash);
  rmSync(temporary, { recursive: true, force: true });
}, 30_000);
