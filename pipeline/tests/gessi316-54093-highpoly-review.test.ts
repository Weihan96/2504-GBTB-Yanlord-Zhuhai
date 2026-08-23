import { expect, test } from "bun:test";
import { createHash } from "node:crypto";
import { existsSync, mkdtempSync, readFileSync, rmSync, writeFileSync } from "node:fs";
import { tmpdir } from "node:os";
import { join, resolve } from "node:path";

const root = resolve(import.meta.dir, "../..");
const product = join(root, "output/review/highpoly-types/gessi316-54093");
const sourceDir = join(product, "official-source");
const formal = join(root, "2504 GBTB Yanlord Zhuhai.ifc");
const formalHash = "7a50b87e8f48a7c2bfbaa9f1dabdfca7155fa684b2325c66aa1ae15c8ab0a25c";
const sourceKind = "native_dwg";
const sourceLabelZh = "Gessi 官方精确型号 54093 G000 原生 DWG 图纸表达";
const scope = "official Gessi exact 54093 G000 family reference; not a project shop drawing";
const dwgHash = "a9308fc8498d34c8cf2f68fd28aaf90b11e62bfc59424e0a5a3ac7f0d28d5047";
const zipHash = "615d4747639b6acfa1daca28a894124fef26e74239de89ef68455b92db31d0f9";
const pdfHash = "2eb67af1be75093a12ab8bc3a739e5bdd6e4b1c4c544948795c87cb9598515c6";
const pathCounts = { plan: 7, front: 914, side: 21 };

function sha256(path: string) {
  return createHash("sha256").update(readFileSync(path)).digest("hex");
}

test("exact public Gessi 54093 G000 DWG and PDF are archived and cross-verified", () => {
  const access = JSON.parse(readFileSync(join(product, "official-source/source-access-record.json"), "utf8"));
  const revalidation = JSON.parse(readFileSync(join(sourceDir, "official-source-revalidation.json"), "utf8"));
  const pdfVerification = JSON.parse(readFileSync(join(sourceDir, "official-pdf-verification.json"), "utf8"));
  expect(access).toMatchObject({
    manufacturer: "Gessi",
    family: "Gessi316 Meccanica",
    project_ifc_type_name: "Gessi316 54093",
    resolved_article_number: "54093",
    resolved_configuration: "G000",
    scope,
    pass: true,
  });
  expect(access.official_product_cad).toMatchObject({
    authentication_required: false,
    acquired: true,
    exact_project_configuration_match: true,
  });
  expect(access.dimension_cross_check).toMatchObject({
    official_api_width_depth_height_mm: [50, 190, 273],
    native_dwg_plan_envelope_mm: [49.973479, 190.312335],
    native_dwg_front_envelope_mm: [50, 273.210406],
    native_dwg_side_envelope_mm: [190.335804, 273.202791],
    project_ifc_body_local_xyz_mm: [48.974943, 190.281433, 273.203123],
    maximum_official_nominal_delta_mm: 0.335804,
    tolerance_mm: 0.5,
    pass: true,
    geometry_stretched: false,
  });
  expect(access.drawing_geometry_source).toMatchObject({ source_kind: sourceKind, source_label_zh: sourceLabelZh, source_dwg_sha256: dwgHash, official_cad_used: true, third_party_cad_used: false });
  expect(access.identity_and_geometry_policy).toMatchObject({ "54093_g000_native_dwg_used": true, g001_or_a004_variant_cad_used: false, adjacent_gessi_product_cad_used: false, third_party_cad_used: false });
  expect(sha256(join(sourceDir, "GPF5409300000G000_3.dwg"))).toBe(dwgHash);
  expect(sha256(join(sourceDir, "GPF5409300000G000_arc.zip"))).toBe(zipHash);
  expect(sha256(join(sourceDir, "GPF5409300000G000_1.pdf"))).toBe(pdfHash);
  expect(revalidation.public_access.exact_54093_g000_native_dwg_publicly_downloadable).toBeTrue();
  expect(revalidation.drawing_geometry_policy.g001_or_a004_variant_cad_used).toBeFalse();
  expect(pdfVerification.mechanical_cross_check.pass).toBeTrue();
  expect(pdfVerification.visual_review.status).toBe("codex_visual_qa_pass");
  for (const source of access.official_identity_sources) {
    expect(sha256(join(root, source.local_path))).toBe(source.sha256);
    if (source.preview_path) expect(sha256(join(root, source.preview_path))).toBe(source.preview_sha256);
  }
});

test("54093 G000 native DWG linework extraction is repeatable", () => {
  const temporary = mkdtempSync(join(tmpdir(), "gessi316-54093-linework-"));
  const output = join(temporary, "linework.json");
  const run = Bun.spawnSync(["python3", join(root, "pipeline/scripts/gessi316_54093_linework.py"), "--output", output], { cwd: root, stdout: "pipe", stderr: "pipe" });
  if (run.exitCode !== 0) throw new Error(run.stderr.toString());
  const linework = JSON.parse(readFileSync(output, "utf8"));
  expect(linework.source_kind).toBe(sourceKind);
  expect(linework.official_sources.native_dwg.sha256).toBe(dwgHash);
  expect(linework.nominal_dimension_cross_check.pass).toBeTrue();
  expect(linework.identity_gates.g001_or_a004_variant_cad_used).toBeFalse();
  for (const [view, count] of Object.entries(pathCounts)) expect(linework.views[view].paths_mm).toHaveLength(count);
  rmSync(temporary, { recursive: true, force: true });
}, 60_000);

test("single-instance 54093 review has exact G000 blue DWG paths in all three views", () => {
  const manifest = JSON.parse(readFileSync(join(product, "manifest.json"), "utf8"));
  const candidate = JSON.parse(readFileSync(join(product, "candidate-representations.json"), "utf8"));
  const profile = JSON.parse(readFileSync(join(product, "profile.json"), "utf8"));
  const inventory = JSON.parse(readFileSync(join(root, "pipeline/decisions/highpoly-product-inventory.json"), "utf8"));
  const entry = inventory.products.find((item: any) => item.type_name === "Gessi316 54093");
  expect(inventory.review_ready_pending_approval_product_count).toBe(38);
  expect(entry).toMatchObject({
    representative_global_id: "1jM_suNMPAw8_nvZejQSnp",
    folder: "output/review/highpoly-types/gessi316-54093",
    status: "review_ready_pending_approval",
  });
  expect(manifest.representative_global_id).toBe("1jM_suNMPAw8_nvZejQSnp");
  expect(manifest.registered_instance_global_ids).toEqual(["1jM_suNMPAw8_nvZejQSnp"]);
  expect(manifest.geometry_product_count).toBe(1);
  expect(manifest.whole_model_render).toBeFalse();
  expect(manifest.formal_ifc_bytes_unchanged).toBeTrue();
  expect(manifest.article_number).toBe("54093");
  expect(profile.profiles["gessi316-54093"].official_reference.article_number).toBe("54093");
  expect(candidate).toMatchObject({
    profile_key: "gessi316-54093",
    article_number: "54093",
    source_kind: sourceKind,
    source_label_zh: sourceLabelZh,
    source_dwg_sha256: dwgHash,
    official_cad_used: true,
    third_party_cad_used: false,
  });
  expect(manifest.views.map((view: any) => view.silhouette_path_count)).toEqual([1, 2, 1]);
  for (const view of manifest.views) {
    expect(view.drawing_line_source_kind).toBe(sourceKind);
    expect(view.official_cad_path_count).toBe(pathCounts[view.view as keyof typeof pathCounts]);
    expect(view.blue_line_present).toBeTrue();
    expect(view.white_mask_present).toBeTrue();
    expect(view.ifc_dwg_compatibility.pass).toBeTrue();
    expect(candidate.views[view.view].official_native_dwg_paths_mm).toHaveLength(pathCounts[view.view as keyof typeof pathCounts]);
    const svg = readFileSync(join(root, view.svg), "utf8");
    expect(svg).toContain('class="simplified-proxy-silhouette geometry-derived"');
    expect(svg).toContain('class="official-reference-mask"');
    expect(svg).toContain('class="official-reference native-dwg"');
    expect(svg).toContain('data-source-kind="native_dwg"');
    expect(svg).toContain("#1677c8");
  }
});

test("54093 project context uses only the real clean sanitary plan and invents no elevation", () => {
  const context = JSON.parse(readFileSync(join(product, "project-context-manifest.json"), "utf8"));
  expect(context).toMatchObject({
    source_kind: sourceKind,
    source_label_zh: sourceLabelZh,
    official_cad_used: true,
    third_party_cad_used: false,
    project_context_retained: true,
    walls_and_surrounding_project_elements_retained: true,
    overlay_top_layer_with_white_mask: true,
    blue_line_present: true,
    pass: true,
  });
  expect(context.context_view_scope).toMatchObject({
    included: ["plan"],
    excluded: ["front", "side"],
  });
  expect(context.project_projection_note.native_project_elevation_found).toBeFalse();
  expect(context.views).toHaveLength(1);
  const view = context.views[0];
  expect(view.view).toBe("plan");
  expect(view.candidate_view).toBe("plan");
  expect(view.source).toBe("drawings/Sanitary Plan.svg");
  expect(view.overlay.fit).toMatchObject({
    scale_svg_units_per_mm: 0.02,
    uniform_scale_preserved: true,
    transformation_mode: "axis_swap_rigid_reflection_and_translation_only",
    pass: true,
  });
  expect(sha256(join(root, view.review_preview))).toBe(view.review_preview_sha256);
  const svg = readFileSync(join(root, view.output), "utf8");
  expect(svg).toContain("IfcWall");
  expect(svg).toContain("IfcSanitaryTerminal");
  expect(view.overlay.source_dwg_sha256).toBe(dwgHash);
  expect(view.overlay.path_count).toBe(pathCounts.plan);
  expect(Math.max(...view.overlay.fit.bbox_absolute_delta_svg_units)).toBeLessThanOrEqual(view.overlay.fit.bbox_tolerance_svg_units);
  expect(svg).toContain('class="official-reference-envelope-mask"');
  expect(svg).toContain('class="official-reference-mask"');
  expect(svg).toContain('class="official-reference native-dwg project-context-overlay"');
  expect(svg).not.toContain("p202-location-marker");
});

test("54093 Bonsai evidence contains four renders from the actual isolated IFC Body", () => {
  const manifest = JSON.parse(readFileSync(join(product, "manifest.json"), "utf8"));
  const evidence = JSON.parse(readFileSync(join(product, "bonsai-review-manifest.json"), "utf8"));
  expect(manifest.bonsai_review.manifest_sha256).toBe(sha256(join(product, "bonsai-review-manifest.json")));
  expect(evidence).toMatchObject({
    mode: "actual_bonsai_ifc_body_camera_render",
    representative_global_id: "1jM_suNMPAw8_nvZejQSnp",
    geometry_product_count: 1,
    whole_model_render: false,
    formal_ifc_bytes_unchanged: true,
  });
  expect(evidence.bonsai_session.saved_active_representation).toBe("Body");
  expect(evidence.bonsai_session.ifc_context_identifier).toBe("Body");
  expect(evidence.bonsai_session.saved_camera_count).toBe(4);
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

test("pending 54093 approval leaves formal IFC byte-identical and creates no derived IFC", () => {
  const candidate = JSON.parse(readFileSync(join(product, "candidate-representations.json"), "utf8"));
  const approval = JSON.parse(readFileSync(join(root, "pipeline/decisions/gessi316-54093-drawing-approval.json"), "utf8"));
  expect(candidate.formal_ifc_write_allowed).toBeFalse();
  expect(candidate.review_status).toBe("visual_review_pending");
  expect(approval).toMatchObject({
    status: "pending",
    approved_views: [],
    derived_ifc_write_allowed: false,
    formal_authoritative_ifc_write_allowed: false,
    scope,
  });
  expect(approval.candidate_manifest_sha256).toBe(sha256(join(product, "manifest.json")));
  expect(existsSync(join(product, "Gessi316-54093-derived-drawing.ifc"))).toBeFalse();
  expect(sha256(formal)).toBe(formalHash);
});

test("54093 writer rejects missing apply and the pending approval record", () => {
  const temporary = mkdtempSync(join(tmpdir(), "gessi316-54093-gate-"));
  const output = join(temporary, "forbidden.ifc");
  const script = join(root, "pipeline/scripts/gessi316_54093_drawing_ifc.py");
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

test("temporary scoped approval writes open 54093 G000 native-DWG representations and source links", () => {
  const manifest = join(product, "manifest.json");
  const temporary = mkdtempSync(join(tmpdir(), "gessi316-54093-approved-"));
  const approval = join(temporary, "approval.json");
  const output = join(temporary, "derived.ifc");
  const report = join(temporary, "report.json");
  writeFileSync(approval, JSON.stringify({
    schema_version: 1,
    profile_key: "gessi316-54093",
    ifc_type_name: "Gessi316 54093",
    candidate_manifest_sha256: sha256(manifest),
    status: "approved",
    reviewer: "automated gate fixture",
    review_date: "2026-08-23",
    approved_views: ["plan", "front", "side"],
    derived_ifc_write_allowed: true,
    formal_authoritative_ifc_write_allowed: false,
    scope,
    approval_evidence: "temporary automated writer verification only",
  }, null, 2));
  const run = Bun.spawnSync([
    "python3", join(root, "pipeline/scripts/gessi316_54093_drawing_ifc.py"),
    "--input", formal, "--manifest", manifest, "--approval", approval,
    "--output", output, "--report", report, "--apply",
  ], { cwd: root, stdout: "pipe", stderr: "pipe" });
  if (run.exitCode !== 0) throw new Error(run.stderr.toString());
  const result = JSON.parse(readFileSync(report, "utf8"));
  expect(result).toMatchObject({
    pass: true,
    formal_ifc_bytes_unchanged: true,
    representations: { plan: "Gessi54093Plan", front: "Gessi54093Front", side: "Gessi54093Side" },
    representation_path_counts: pathCounts,
    source_kind: sourceKind,
    source_label_zh: sourceLabelZh,
    official_cad_geometry_included: true,
    source_property_set: "Pset_Gessi31654093DrawingSource",
  });
  expect(result.representation_geometry_source).toBe("official_native_dwg_paths_mm extracted from GPF5409300000G000_3.dwg");
  expect(result.source_document_associations).toHaveLength(6);
  expect(result.source_document_associations).toContain("GESSI316-54093-OFFICIAL-NATIVE-DWG-ZIP");
  expect(result.source_document_associations).toContain("GESSI316-54093-OFFICIAL-TECHNICAL-PDF");
  const inspect = Bun.spawnSync(["python3", "-c", [
    "import ifcopenshell, ifcopenshell.util.element, json, sys",
    "m=ifcopenshell.open(sys.argv[1])",
    "p=m.by_guid('1jM_suNMPAw8_nvZejQSnp')",
    "print(json.dumps(ifcopenshell.util.element.get_pset(p, 'Pset_Gessi31654093DrawingSource')))"
  ].join(";"), output], { cwd: root, stdout: "pipe", stderr: "pipe" });
  if (inspect.exitCode !== 0) throw new Error(inspect.stderr.toString());
  const pset = JSON.parse(inspect.stdout.toString());
  expect(pset).toMatchObject({
    ArticleNumber: "54093",
    TechnicalDrawingNumber: "GPF5409300000G000",
    SourceDwgSha256: dwgHash,
    OfficialTechnicalDrawingLocalArchive: "true",
    OfficialProductCadStatus: "public_official_api_exact_54093_g000_native_dwg_acquired",
    OfficialCadUsed: "true",
    OfficialVectorEvidenceUsedAsCadGeometry: "false",
    DimensionCrossCheckPass: "true",
    NativeProjectElevationFound: "false",
    EvidenceScope: scope,
  });
  expect(sha256(formal)).toBe(formalHash);
  rmSync(temporary, { recursive: true, force: true });
}, 120_000);
