import { expect, test } from "bun:test";
import { createHash } from "node:crypto";
import { existsSync, mkdtempSync, readFileSync, rmSync, writeFileSync } from "node:fs";
import { tmpdir } from "node:os";
import { join, resolve } from "node:path";

const root = resolve(import.meta.dir, "../..");
const product = join(root, "output/review/highpoly-types/geberit-154-154-00-1-f");
const formal = join(root, "2504 GBTB Yanlord Zhuhai.ifc");
const formalHash = "7a50b87e8f48a7c2bfbaa9f1dabdfca7155fa684b2325c66aa1ae15c8ab0a25c";
const sourceKind = "geometry_derived_simplified_proxy";
const sourceLabelZh = "基于原始高模几何生成的简化图纸表达";
const scope = "project-model flange subcomponent of the exact Geberit 154.154.00.1 parent installation set; .F is a project authoring decomposition label, not a separately published manufacturer article; parent EPS is not component CAD geometry or a project shop drawing";

function sha256(path: string) {
  return createHash("sha256").update(readFileSync(path)).digest("hex");
}

test("official evidence preserves the parent/component boundary and never claims .F as an article", () => {
  const access = JSON.parse(readFileSync(join(product, "official-source/source-access-record.json"), "utf8"));
  expect(access).toMatchObject({
    manufacturer: "Geberit",
    family: "CleanLine shower channel installation set",
    project_ifc_type_name: "Geberit 154.154.00.1.F",
    parent_manufacturer_article_number: "154.154.00.1",
    project_component_code: "154.154.00.1.F",
    project_component_role: "Flange",
    scope,
    pass: true,
  });
  expect(access.component_identity_boundary).toContain("not presented as a separately published Geberit article number");
  expect(access.official_component_cad).toMatchObject({
    published_as_separate_manufacturer_article: false,
    acquired: false,
    exact_project_component_match: false,
    local_cad_files: [],
    official_parent_vector_eps_archived: true,
    official_parent_vector_eps_used_as_component_cad_geometry: false,
    parent_native_dwg_http_statuses: { A: 404, G: 404, L: 404, P: 404 },
  });
  expect(access.project_component_geometry).toMatchObject({
    body_local_xyz_mm: [88, 210, 104.999998],
    official_standalone_component_dimensions_available: false,
    mechanical_dimension_cross_check_applicable: false,
  });
  expect(access.project_decomposition_evidence).toMatchObject({
    parent_installation_set_type: "Geberit 154.154.00",
    sibling_cover_type: "Geberit 154.154.00.1.C",
    component_type: "Geberit 154.154.00.1.F",
  });
  for (const source of access.official_identity_sources) {
    expect(sha256(join(root, source.local_path))).toBe(source.sha256);
    if (source.preview_path) expect(sha256(join(root, source.preview_path))).toBe(source.preview_sha256);
  }
});

test("one .F Body produces 2/7/6 black paths and zero parent-EPS blue paths", () => {
  const manifest = JSON.parse(readFileSync(join(product, "manifest.json"), "utf8"));
  const candidate = JSON.parse(readFileSync(join(product, "candidate-representations.json"), "utf8"));
  const profile = JSON.parse(readFileSync(join(product, "profile.json"), "utf8"));
  const inventory = JSON.parse(readFileSync(join(root, "pipeline/decisions/highpoly-product-inventory.json"), "utf8"));
  expect(inventory.review_ready_pending_approval_product_count).toBe(38);
  expect(inventory.products.find((entry: any) => entry.type_name === "Geberit 154.154.00.1.F")).toMatchObject({
    representative_global_id: "2jjNIn9gHBYwNSWwlM5T_i",
    status: "review_ready_pending_approval",
    folder: "output/review/highpoly-types/geberit-154-154-00-1-f",
  });
  expect(profile.profiles["geberit-154-154-00-1-f"].official_reference).toMatchObject({
    parent_article_number: "154.154.00.1",
    project_component_code: "154.154.00.1.F",
    project_component_role: "Flange",
  });
  expect(manifest).toMatchObject({
    representative_global_id: "2jjNIn9gHBYwNSWwlM5T_i",
    registered_instance_global_ids: ["2jjNIn9gHBYwNSWwlM5T_i"],
    geometry_product_count: 1,
    whole_model_render: false,
    formal_ifc_bytes_unchanged: true,
    bounds_mm: { size: [88, 210, 104.999998] },
  });
  expect(candidate).toMatchObject({
    article_number: "154.154.00.1 / project component .F",
    source_kind: sourceKind,
    source_label_zh: sourceLabelZh,
    official_cad_used: false,
    third_party_cad_used: false,
    review_status: "visual_review_pending",
  });
  expect(Object.fromEntries(Object.entries(candidate.views).map(([view, value]: any) => [view, value.proxy_paths_mm.length]))).toEqual({ plan: 2, front: 7, side: 6 });
  for (const view of manifest.views) {
    expect(view.official_cad_path_count).toBe(0);
    expect(view.blue_line_present).toBeFalse();
    const svg = readFileSync(join(root, view.svg), "utf8");
    expect(svg).toContain(`data-source-kind="${sourceKind}"`);
    expect(svg).not.toContain('class="official-reference native-dwg"');
    expect(svg).not.toContain("#1677c8");
  }
});

test("real project plan and R12 elevations retain context while disclosing clipped native elevations", () => {
  const context = JSON.parse(readFileSync(join(product, "project-context-manifest.json"), "utf8"));
  expect(context).toMatchObject({
    source_kind: sourceKind,
    official_cad_used: false,
    third_party_cad_used: false,
    project_context_retained: true,
    walls_and_surrounding_project_elements_retained: true,
    overlay_top_layer_with_white_mask: true,
    blue_line_present: false,
    pass: true,
  });
  expect(context.views.map((view: any) => [view.view, view.source])).toEqual([
    ["plan", "drawings/Sanitary Plan.svg"],
    ["front", "drawings/elevations/native/EL-06-20-R12-NY.svg"],
    ["side", "drawings/elevations/native/EL-06-21-R12-PX.svg"],
  ]);
  expect(context.review_annotation_suppression).toMatchObject({
    plan: ["official-elevation-anchor"],
    walls_furniture_and_ifc_geometry_removed: false,
  });
  expect(context.project_projection_note).toMatchObject({
    native_plan_projection_bbox_mm: [88, 210],
    native_front_projection_bbox_mm: [63.100047, 57.55335],
    native_side_projection_bbox_mm: [210, 18.3831],
    body_projection_bbox_mm: { plan: [88, 210], front: [88, 104.999998], side: [210, 104.999998] },
    geometry_stretched: false,
    native_project_plan_and_elevations_found: true,
  });
  for (const view of context.views) {
    expect(view.overlay.fit.scale_svg_units_per_mm).toBe(0.02);
    expect(view.overlay.fit.uniform_scale_preserved).toBeTrue();
    expect(view.overlay.fit.transformation_mode).toBe("axis_swap_rigid_reflection_and_translation_only");
    expect(Math.max(...view.overlay.fit.bbox_absolute_delta_svg_units)).toBeLessThanOrEqual(view.overlay.fit.bbox_tolerance_svg_units);
    expect(sha256(join(root, view.review_preview))).toBe(view.review_preview_sha256);
    const svg = readFileSync(join(root, view.output), "utf8");
    expect(svg).toContain('class="geometry-derived-proxy-mask"');
    expect(svg).toContain('class="geometry-derived-proxy"');
    expect(svg.indexOf('class="geometry-derived-proxy-mask"')).toBeLessThan(svg.indexOf('class="geometry-derived-proxy"'));
    expect(svg).not.toContain('class="official-reference native-dwg"');
  }
  const reviewPlan = readFileSync(join(product, "project-context-sanitary-plan-review.svg"), "utf8");
  expect(reviewPlan).not.toContain('class="official-elevation-anchor"');
  expect(readFileSync(join(product, "project-context-sanitary-plan.svg"), "utf8")).toContain('class="official-elevation-anchor"');
});

test("actual Bonsai session renders only the isolated .F MODEL_VIEW Body with four cameras", () => {
  const evidence = JSON.parse(readFileSync(join(product, "bonsai-review-manifest.json"), "utf8"));
  expect(evidence).toMatchObject({
    mode: "actual_bonsai_ifc_body_camera_render",
    representative_global_id: "2jjNIn9gHBYwNSWwlM5T_i",
    geometry_product_count: 1,
    whole_model_render: false,
    formal_ifc_bytes_unchanged: true,
  });
  expect(evidence.bonsai_session).toMatchObject({ saved_active_representation: "Body", ifc_context_identifier: "Body", saved_camera_count: 4 });
  expect(evidence.renders.map((render: any) => render.view)).toEqual(["plan", "front", "side", "iso"]);
  for (const render of evidence.renders) {
    expect(render.camera_type).toBe("ORTHO");
    expect(render.render_operation).toBe("bpy.ops.render.render(write_still=True)");
    expect(render.camera_ortho_scale_m + 1e-8).toBeGreaterThanOrEqual(render.projected_width_m * render.framing_margin_factor);
    expect(render.camera_ortho_scale_m + 1e-8).toBeGreaterThanOrEqual(render.projected_height_m * render.image_aspect * render.framing_margin_factor);
    expect(sha256(join(root, render.path))).toBe(render.sha256);
  }
  expect(sha256(join(root, evidence.bonsai_session.path))).toBe(evidence.bonsai_session.sha256);
  expect(sha256(join(root, evidence.isolated_ifc))).toBe(evidence.isolated_ifc_sha256);
});

test("pending .F review blocks derived IFC and keeps the formal IFC byte-identical", () => {
  const approval = JSON.parse(readFileSync(join(root, "pipeline/decisions/geberit-154-154-00-1-f-drawing-approval.json"), "utf8"));
  expect(approval).toMatchObject({ status: "pending", approved_views: [], derived_ifc_write_allowed: false, formal_authoritative_ifc_write_allowed: false, scope });
  expect(approval.candidate_manifest_sha256).toBe(sha256(join(product, "manifest.json")));
  expect(existsSync(join(product, "Geberit-154-154-00-1-F-derived-drawing.ifc"))).toBeFalse();
  expect(sha256(formal)).toBe(formalHash);
});

test(".F writer rejects missing apply and pending human approval", () => {
  const temporary = mkdtempSync(join(tmpdir(), "geberit-flange-gate-"));
  const output = join(temporary, "forbidden.ifc");
  const script = join(root, "pipeline/scripts/geberit_154_154_00_1_f_drawing_ifc.py");
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

test("temporary approval writes Body-derived .F representations and parent-evidence associations", () => {
  const manifest = join(product, "manifest.json");
  const temporary = mkdtempSync(join(tmpdir(), "geberit-flange-approved-"));
  const approval = join(temporary, "approval.json");
  const output = join(temporary, "derived.ifc");
  const report = join(temporary, "report.json");
  writeFileSync(approval, JSON.stringify({
    schema_version: 1, profile_key: "geberit-154-154-00-1-f", ifc_type_name: "Geberit 154.154.00.1.F",
    candidate_manifest_sha256: sha256(manifest), status: "approved", reviewer: "automated gate fixture",
    review_date: "2026-08-23", approved_views: ["plan", "front", "side"],
    derived_ifc_write_allowed: true, formal_authoritative_ifc_write_allowed: false,
    scope, approval_evidence: "temporary automated writer verification only",
  }, null, 2));
  const run = Bun.spawnSync([
    "python3", join(root, "pipeline/scripts/geberit_154_154_00_1_f_drawing_ifc.py"),
    "--input", formal, "--manifest", manifest, "--approval", approval,
    "--output", output, "--report", report, "--apply",
  ], { cwd: root, stdout: "pipe", stderr: "pipe" });
  if (run.exitCode !== 0) throw new Error(run.stderr.toString());
  const result = JSON.parse(readFileSync(report, "utf8"));
  expect(result).toMatchObject({
    pass: true,
    formal_ifc_bytes_unchanged: true,
    representations: { plan: "Geberit154154FlangePlan", front: "Geberit154154FlangeFront", side: "Geberit154154FlangeSide" },
    representation_path_counts: { plan: 2, front: 7, side: 6 },
    source_kind: sourceKind,
    source_label_zh: sourceLabelZh,
    official_cad_geometry_included: false,
    source_property_set: "Pset_Geberit154154FlangeDrawingSource",
  });
  expect(result.source_document_associations).toHaveLength(6);
  expect(result.source_document_associations).toContain("GEBERIT-154-154-00-1-F-OFFICIAL-PARENT-PRODUCT-PAGE");
  const inspect = Bun.spawnSync(["python3", "-c", [
    "import ifcopenshell,ifcopenshell.util.element,json,sys",
    "m=ifcopenshell.open(sys.argv[1])",
    "p=m.by_guid('2jjNIn9gHBYwNSWwlM5T_i')",
    "print(json.dumps(ifcopenshell.util.element.get_pset(p,'Pset_Geberit154154FlangeDrawingSource')))"
  ].join(";"), output], { cwd: root, stdout: "pipe", stderr: "pipe" });
  if (inspect.exitCode !== 0) throw new Error(inspect.stderr.toString());
  expect(JSON.parse(inspect.stdout.toString())).toMatchObject({
    ParentArticleNumber: "154.154.00.1",
    ProjectComponentCode: "154.154.00.1.F",
    ProjectComponentRole: "Flange",
    ComponentCodeIsManufacturerArticle: "false",
    OfficialComponentCadUsed: "false",
    ParentVectorEPSUsedAsComponentGeometry: "false",
    OfficialStandaloneComponentDimensionsAvailable: "false",
    NativeProjectElevationClippingRecorded: "true",
    EvidenceScope: scope,
  });
  expect(sha256(formal)).toBe(formalHash);
  rmSync(temporary, { recursive: true, force: true });
}, 30_000);
