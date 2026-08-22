import { expect, test } from "bun:test";
import {
  existsSync,
  mkdtempSync,
  readFileSync,
  writeFileSync,
} from "node:fs";
import { createHash } from "node:crypto";
import { tmpdir } from "node:os";
import { join, resolve } from "node:path";

const root = resolve(import.meta.dir, "../..");
const formalIfc = join(root, "2504 GBTB Yanlord Zhuhai.ifc");
const formalHash = "7a50b87e8f48a7c2bfbaa9f1dabdfca7155fa684b2325c66aa1ae15c8ab0a25c";
const wfaHash = "48ffe94e5a692cd03e8a09e67c658ebb42362817df03007bd90ff7a912567441";
const wfbHash = "4b3d8900f8ae249d7e9634339f3529d7d38a940dc2d0c909a7e90e31d3946831";
const pdfHash = "b8bc94982bdd4305f2ffe5240daf429146470677157022860a0ef591514373c2";

function sha256(path: string): string {
  return createHash("sha256").update(readFileSync(path)).digest("hex");
}

test("archived Falper sources preserve native WFA/WFB identity", () => {
  expect(sha256(join(root, "drawings/evidence/FALPER-official-Sorgente-WFA-2D.dwg"))).toBe(wfaHash);
  expect(sha256(join(root, "drawings/evidence/FALPER-official-Sorgente-WFB-2D.dwg"))).toBe(wfbHash);
  expect(sha256(join(root, "drawings/evidence/FALPER-official-Sorgente-WFA-WFB.pdf"))).toBe(pdfHash);
  expect(sha256(formalIfc)).toBe(formalHash);
});

test("native DWG extractor decodes five plan and four elevation paths per model", () => {
  const temporary = mkdtempSync(join(tmpdir(), "falper-dwg-"));
  const output = join(temporary, "linework.json");
  const run = Bun.spawnSync([
    "python3",
    "pipeline/scripts/extract_falper_sorgente_dwg_linework.py",
    "--output",
    output,
  ], { cwd: root });
  expect(run.exitCode, run.stderr.toString()).toBe(0);
  const payload = JSON.parse(readFileSync(output, "utf8"));
  expect(payload.source_kind).toBe("native_dwg");
  expect(payload.scope).toBe("family_reference_not_project_shop_drawing");
  expect(payload.identity_guards.wfa_3d_used_as_wfb).toBe(false);
  expect(payload.variants.WFA.source_dwg_sha256).toBe(wfaHash);
  expect(payload.variants.WFB.source_dwg_sha256).toBe(wfbHash);
  for (const model of ["WFA", "WFB"]) {
    expect(payload.variants[model].views.plan.paths_mm).toHaveLength(5);
    expect(payload.variants[model].views.elevation.paths_mm).toHaveLength(4);
    expect(payload.variants[model].views.elevation.native_entity_types).toEqual([
      "SPLINE", "SPLINE", "LINE", "LINE",
    ]);
  }
});

test("official PDF mechanically cross-validates each WFA/WFB DWG view", () => {
  const report = JSON.parse(
    readFileSync(join(root, "pipeline/decisions/falper-sorgente-dwg-pdf-verification.json"), "utf8"),
  );
  expect(report.pass).toBe(true);
  expect(report.mode).toBe("mechanical_native_dwg_vs_official_vector_pdf");
  expect(report.identity_gates.blue_line_source_kind).toBe("native_dwg");
  expect(report.identity_gates.wfa_3d_used_as_wfb).toBe(false);
  expect(report.records).toHaveLength(4);
  expect(report.records.every((record: any) => record.pass)).toBe(true);
  expect(report.records.every((record: any) => record.hausdorff_mm <= record.tolerance_mm)).toBe(true);
  expect(report.wfa_wfb_native_family_geometry_diff.plan.hausdorff_mm).toBe(0);
  expect(report.wfa_wfb_native_family_geometry_diff.elevation.hausdorff_mm).toBe(0);
});

test("high-poly review generates only one instance and layers native WFB blue on top", () => {
  const temporary = mkdtempSync(join(tmpdir(), "falper-highpoly-"));
  const before = sha256(formalIfc);
  const run = Bun.spawnSync([
    "python3",
    "pipeline/scripts/int1_highpoly_type_review.py",
    "--input",
    formalIfc,
    "--profile-key",
    "falper-sorgente",
    "--output-root",
    temporary,
  ], { cwd: root });
  expect(run.exitCode, run.stderr.toString()).toBe(0);
  const directory = join(temporary, "falper-sorgente");
  const manifest = JSON.parse(readFileSync(join(directory, "manifest.json"), "utf8"));
  expect(manifest.geometry_product_count).toBe(1);
  expect(manifest.whole_model_render).toBe(false);
  expect(manifest.formal_ifc_write).toBe(false);
  expect(manifest.formal_ifc_bytes_unchanged).toBe(true);
  expect(manifest.representative_global_id).toBe("350tdaubr8QP3Cu2YMQZIN");
  expect(manifest.official_reference.source_kind).toBe("native_dwg");
  expect(manifest.official_reference.source_dwg_sha256).toBe(wfbHash);
  expect(manifest.official_reference.verification_pass).toBe(true);
  expect(manifest.review_status).toBe("visual_review_pending");
  expect(manifest.approved_for_drawing_ifc).toBe(false);
  expect(manifest.views.map((view: any) => view.view)).toEqual(["plan", "front", "side"]);
  const candidate = JSON.parse(readFileSync(join(directory, "candidate-representations.json"), "utf8"));
  expect(candidate.schema_version).toBe(2);
  expect(candidate.views.plan.proxy_paths_mm).toHaveLength(110);
  expect(candidate.views.front.proxy_paths_mm).toHaveLength(117);
  expect(candidate.views.side.proxy_paths_mm).toHaveLength(120);
  expect(candidate.views.plan.official_native_dwg_paths_mm).toHaveLength(5);
  expect(candidate.views.front.official_native_dwg_paths_mm).toHaveLength(4);
  expect(candidate.views.side.official_native_dwg_paths_mm).toHaveLength(4);
  expect(candidate.views.front.paths_mm).toBeUndefined();
  for (const view of manifest.views) {
    expect(view.blue_line_source_kind).toBe("native_dwg");
    const svg = readFileSync(join(directory, `${view.view}.svg`), "utf8");
    expect(svg).toContain('class="original-highpoly"');
    expect(svg).toContain('class="simplified-proxy-silhouette"');
    expect(svg).toContain('class="official-reference-mask"');
    expect(svg).toContain('class="official-reference native-dwg"');
    expect(svg).toContain('data-source-kind="native_dwg"');
    expect(svg.indexOf('class="official-reference-mask"')).toBeLessThan(
      svg.indexOf('class="official-reference native-dwg"'),
    );
  }
  expect(existsSync(join(directory, "index.html"))).toBe(true);
  expect(sha256(formalIfc)).toBe(before);
}, 30_000);

test("pending approval cannot write any IFC bytes", () => {
  const temporary = mkdtempSync(join(tmpdir(), "falper-gate-"));
  const approvalPath = join(temporary, "approval.json");
  const output = join(temporary, "derived.ifc");
  writeFileSync(approvalPath, JSON.stringify({
    schema_version: 1,
    profile_key: "falper-sorgente",
    model_code: "WFB",
    candidate_manifest_sha256: null,
    status: "visual_review_pending",
    reviewer: null,
    review_date: null,
    approved_views: [],
    formal_ifc_write_allowed: false,
    scope: "family_reference_not_project_shop_drawing",
  }));
  const before = sha256(formalIfc);
  const run = Bun.spawnSync([
    "python3",
    "pipeline/scripts/falper_sorgente_drawing_ifc.py",
    "--input",
    formalIfc,
    "--approval",
    approvalPath,
    "--output",
    output,
    "--apply",
  ], { cwd: root });
  expect(run.exitCode).not.toBe(0);
  expect(run.stderr.toString()).toContain("approval gate rejected IFC write");
  expect(existsSync(output)).toBe(false);
  expect(sha256(formalIfc)).toBe(before);
});

test("an exact explicit approval writes only a temporary derived IFC with source associations", () => {
  const temporary = mkdtempSync(join(tmpdir(), "falper-approved-"));
  const manifest = join(root, "output/review/highpoly-types/falper-sorgente/manifest.json");
  const approvalPath = join(temporary, "approval.json");
  const output = join(temporary, "derived.ifc");
  const report = join(temporary, "derived-ifc-manifest.json");
  writeFileSync(approvalPath, JSON.stringify({
    schema_version: 1,
    profile_key: "falper-sorgente",
    model_code: "WFB",
    candidate_manifest_sha256: sha256(manifest),
    status: "approved",
    reviewer: "mechanical-test-reviewer",
    review_date: "2026-08-22",
    approved_views: ["plan", "front", "side"],
    formal_ifc_write_allowed: true,
    scope: "family_reference_not_project_shop_drawing",
  }));
  const before = sha256(formalIfc);
  const run = Bun.spawnSync([
    "python3",
    "pipeline/scripts/falper_sorgente_drawing_ifc.py",
    "--input",
    formalIfc,
    "--manifest",
    manifest,
    "--approval",
    approvalPath,
    "--output",
    output,
    "--report",
    report,
    "--apply",
  ], { cwd: root });
  expect(run.exitCode, run.stderr.toString()).toBe(0);
  const result = JSON.parse(run.stdout.toString());
  expect(result.pass).toBe(true);
  expect(result.source_document_association).toBe(
    "IfcDocumentReference/IfcRelAssociatesDocument",
  );
  expect(result.source_property_set).toBe("Pset_FalperSorgenteDrawingSource");
  expect(result.formal_ifc_bytes_unchanged).toBe(true);
  expect(result.source_kind).toBe("native_dwg");
  expect(result.source_dwg_sha256).toBe(wfbHash);
  expect(result.representation_path_counts).toEqual({ plan: 5, front: 4, side: 4 });
  expect(result.bonsai_drawing_representation_path_counts).toEqual({ plan: 5, front: 4 });
  expect(result.bonsai_drawing_representation_selection).toBe(
    "Model/Body PLAN_VIEW and ELEVATION_VIEW",
  );
  expect(result.representation_geometry_source).toBe("official_native_dwg_paths_mm");
  expect(result.proxy_geometry_included).toBe(false);
  expect(result.representation_source_mapping).toBe(
    "FalperWFBPlan=plan;FalperWFBFront=elevation;FalperWFBSide=elevation",
  );
  expect(result.representations).toEqual(expect.arrayContaining([
    "FalperWFBPlan", "FalperWFBFront", "FalperWFBSide",
  ]));
  expect(existsSync(output)).toBe(true);
  expect(JSON.parse(readFileSync(report, "utf8"))).toEqual(result);
  expect(sha256(formalIfc)).toBe(before);
}, 30_000);

test("Bonsai evidence distinguishes camera renders from representation selection", () => {
  const manifest = JSON.parse(readFileSync(
    join(root, "output/review/highpoly-types/falper-sorgente/bonsai-review-manifest.json"),
    "utf8",
  ));
  expect(manifest.review_status).toBe("approved");
  expect(manifest.geometry_product_count).toBe(1);
  expect(manifest.whole_model_render).toBe(false);
  expect(manifest.bonsai_session.saved_active_representation).toBe("Body");
  expect(manifest.bonsai_session.saved_camera_count).toBe(4);
  expect(manifest.viewport_body_evidence_is_render_result).toBe(false);
  expect(manifest.representation_selection_evidence_is_render_result).toBe(false);
  expect(manifest.camera_renders_are_render_results).toBe(true);
  expect(manifest.camera_renders_accepted_as_project_context_drawing_evidence).toBe(false);
  expect(manifest.project_context_bonsai_drawing_manifest).toBe(
    "output/review/highpoly-types/falper-sorgente/bonsai-project-context-styled/manifest.json",
  );
  expect(manifest.camera_renders.render_source).toBe("actual_ifc_body_representation");
  expect(manifest.camera_renders.render_engine).toBe("BLENDER_WORKBENCH");
  expect(manifest.camera_renders.resolution_px).toEqual([1200, 1200]);
  for (const view of ["plan", "front_elevation", "side_elevation", "isometric"]) {
    const render = manifest.camera_renders[view];
    const path = join(root, render.path);
    expect(existsSync(path)).toBe(true);
    expect(sha256(path)).toBe(render.sha256);
    expect(render.camera).toStartWith("FALPER_CAM_");
  }
  expect(sha256(join(root, manifest.isolated_ifc.path))).toBe(manifest.isolated_ifc.sha256);
  expect(sha256(join(root, manifest.bonsai_session.path))).toBe(
    manifest.bonsai_session.sha256,
  );
  expect(sha256(formalIfc)).toBe(formalHash);
});

test("actual Bonsai Drawing cameras select official WFB geometry in project context", () => {
  const directory = join(root, "output/review/highpoly-types/falper-sorgente");
  const result = JSON.parse(readFileSync(
    join(directory, "bonsai-project-context-styled/manifest.json"),
    "utf8",
  ));
  expect(result.pass).toBe(true);
  expect(result.bonsai_drawing_operator).toBe("bpy.ops.bim.create_drawing");
  expect(result.standard_ifc_representation_selection).toBe(
    "Model/Body PLAN_VIEW and ELEVATION_VIEW",
  );
  expect(result.source_kind).toBe("native_dwg");
  expect(result.source_dwg_sha256).toBe(wfbHash);
  expect(result.formal_ifc_bytes_unchanged).toBe(true);
  expect(result.views.map((view: any) => view.view)).toEqual(["plan", "front_elevation"]);
  expect(result.views[0].raw_bonsai_product_path_segment_count).toBeGreaterThan(1000);
  expect(result.views[1].raw_bonsai_product_path_segment_count).toBeGreaterThan(500);
  const planSource = JSON.parse(readFileSync(
    join(directory, "bonsai-project-context-plan/Sanitary Plan-source.json"),
    "utf8",
  ));
  expect(planSource.camera_local_y_offset_m).toBe(-0.02);
  expect(Math.abs(
    result.views[0].raw_bonsai_visible_bounds_mm[3]
    - result.views[0].ifc_native_dwg_full_bounds_mm[3]
  )).toBeLessThanOrEqual(0.2);
  expect(result.views[0].geometry_supplemented_from_same_ifc_representation).toBe(true);
  expect(result.views[0].ifc_representation_supplement_path_count).toBe(5);
  expect(result.views[0].raw_bonsai_group_visible_in_styled_svg).toBe(false);
  expect(result.views[1].raw_bonsai_group_visible_in_styled_svg).toBe(true);
  expect(result.views.every((view: any) => view.geometry_replaced_during_styling === false)).toBe(true);
  expect(result.views.every((view: any) => view.white_mask && view.official_blue_top_layer)).toBe(true);
  for (const view of result.views) {
    expect(sha256(join(root, view.raw_bonsai_svg))).toBe(view.raw_bonsai_svg_sha256);
    expect(sha256(join(root, view.styled_svg))).toBe(view.styled_svg_sha256);
    const svg = readFileSync(join(root, view.styled_svg), "utf8");
    expect(svg).toContain('id="falper-bonsai-approved-top-layer"');
    expect(svg).toContain('data-geometry-origin="actual_ifc_body_drawing_representation"');
    expect(svg).toContain('data-source-kind="native_dwg"');
    expect(svg).toContain('data-geometry-replaced="false"');
    expect(svg).toContain("stroke:#1677c8");
    if (view.view === "plan") {
      expect(svg).toContain('id="falper-wfb-plan-ifc-representation-supplement"');
      expect(svg).toContain('data-reason="Bonsai hidden-line removal occluded wall-adjacent curve segments"');
    }
    expect(svg.indexOf('class="official-reference-mask"')).toBeLessThan(
      svg.indexOf('data-representation="Body/'),
    );
    expect(svg).toContain('class="IfcWall');
  }
  const counterexample = JSON.parse(readFileSync(
    join(directory, "INVALIDATED-isolated-orange-render-evidence.json"),
    "utf8",
  ));
  expect(counterexample.status).toBe("counterexample");
  expect(counterexample.replacement_manifest).toBe(
    "output/review/highpoly-types/falper-sorgente/bonsai-project-context-styled/manifest.json",
  );
  expect(sha256(formalIfc)).toBe(formalHash);
});

test("project plan and elevation retain context with native-DWG Falper overlays", () => {
  const temporary = mkdtempSync(join(tmpdir(), "falper-project-context-"));
  const sanitaryPlan = join(root, "drawings/Sanitary Plan.svg");
  const elevation = join(root, "drawings/elevations/native/EL-08-30-R16-NY.svg");
  const planBefore = sha256(sanitaryPlan);
  const elevationBefore = sha256(elevation);
  const run = Bun.spawnSync([
    "python3",
    "pipeline/scripts/render_falper_sorgente_project_context.py",
    "--output-dir",
    temporary,
  ], { cwd: root });
  expect(run.exitCode, run.stderr.toString()).toBe(0);

  const report = JSON.parse(readFileSync(join(temporary, "project-context-manifest.json"), "utf8"));
  expect(report.formal_ifc_write).toBe(false);
  expect(report.formal_ifc_sha256).toBe(formalHash);
  expect(report.source_kind).toBe("native_dwg");
  expect(report.source_dwg_sha256).toBe(wfbHash);
  expect(report.views.map((view: any) => view.native_dwg_path_count)).toEqual([5, 4]);
  expect(report.views.every((view: any) => view.white_mask && view.blue_top_layer)).toBe(true);

  const plan = readFileSync(join(temporary, "project-context-sanitary-plan.svg"), "utf8");
  const front = readFileSync(
    join(temporary, "project-context-elevation-EL-08-30-R16-NY.svg"),
    "utf8",
  );
  expect(plan).toContain('id="falper-wfb-plan-approved-overlay"');
  expect(plan).toContain('data-representation="FalperWFBPlan"');
  expect(plan).toContain('data-source-kind="native_dwg"');
  expect(plan).toContain('class="IfcWall');
  expect(front).toContain('id="falper-wfb-front-approved-overlay"');
  expect(front).toContain('data-representation="FalperWFBFront"');
  expect(front).toContain('data-source-kind="native_dwg"');
  expect(front).toContain('class="IfcWall');
  expect(front).toMatch(/M [-+0-9.]+,[-+0-9.]+ L [-+0-9.]+,[-+0-9.]+/);
  expect(plan.indexOf('id="falper-original-body-projection-removed"')).toBeLessThan(
    plan.indexOf('id="falper-wfb-plan-approved-overlay"'),
  );
  expect(front.indexOf('id="falper-original-body-projection-removed"')).toBeLessThan(
    front.indexOf('id="falper-wfb-front-approved-overlay"'),
  );
  expect(sha256(sanitaryPlan)).toBe(planBefore);
  expect(sha256(elevation)).toBe(elevationBefore);
  expect(sha256(formalIfc)).toBe(formalHash);
}, 30_000);
