import { expect, test } from "bun:test";
import { createHash } from "node:crypto";
import { existsSync, mkdtempSync, readFileSync, rmSync, writeFileSync } from "node:fs";
import { tmpdir } from "node:os";
import { join, resolve } from "node:path";

const root = resolve(import.meta.dir, "../..");
const product = join(root, "output/review/highpoly-types/gessi316-54145");
const sourceDir = join(product, "official-source");
const formal = join(root, "2504 GBTB Yanlord Zhuhai.ifc");
const formalHash = "7a50b87e8f48a7c2bfbaa9f1dabdfca7155fa684b2325c66aa1ae15c8ab0a25c";
const sourceKind = "native_dwg";
const sourceLabelZh = "Gessi 官方精确型号 54145 G000 原生 DWG 图纸表达";
const scope = "official Gessi exact 54145 G000 family reference; not a project shop drawing";
const dwgHash = "9978b68468a61875acb0736aab45a62a98efadfd6be61d347e7ecaa94e08fd09";
const zipHash = "9678c2de6a276c0d1f6e3b29e6c39763c1be0763bccf1b8974e68408ede699fd";
const pdfHash = "148a16717193cbc066d0314c6d5369492ed582da0f026dbbe6a58d06f8279025";
const pathCounts = { plan: 22, front: 615, side: 653 };

function sha256(path: string) {
  return createHash("sha256").update(readFileSync(path)).digest("hex");
}

test("exact public Gessi 54145 G000 DWG and PDF are archived and cross-verified", () => {
  const access = JSON.parse(readFileSync(join(sourceDir, "source-access-record.json"), "utf8"));
  const revalidation = JSON.parse(readFileSync(join(sourceDir, "official-source-revalidation.json"), "utf8"));
  const pdfVerification = JSON.parse(readFileSync(join(sourceDir, "official-pdf-verification.json"), "utf8"));
  expect(access).toMatchObject({
    manufacturer: "Gessi", family: "Gessi316 Meccanica", project_ifc_type_name: "Gessi316 54145",
    resolved_article_number: "54145", resolved_configuration: "G000", scope, pass: true,
  });
  expect(access.official_product_cad).toMatchObject({ authentication_required: false, acquired: true, exact_project_configuration_match: true });
  expect(access.dimension_cross_check).toMatchObject({
    official_api_width_depth_height_mm: [300, 600, 115],
    native_dwg_plan_envelope_mm: [600, 299.993081],
    native_dwg_front_envelope_mm: [600, 118.95],
    native_dwg_side_envelope_mm: [300, 118.9206],
    project_ifc_body_local_xyz_mm: [599.818665, 299.818832, 119.07444],
    maximum_ifc_projection_delta_mm: 0.181335, tolerance_mm: 0.5, pass: true, geometry_stretched: false,
  });
  expect(access.drawing_geometry_source).toMatchObject({ source_kind: sourceKind, source_label_zh: sourceLabelZh, source_dwg_sha256: dwgHash, official_cad_used: true, third_party_cad_used: false });
  expect(access.identity_and_geometry_policy).toMatchObject({ "54145_g000_native_dwg_used": true, g001_variant_cad_used: false, adjacent_gessi_product_cad_used: false, third_party_cad_used: false });
  expect(sha256(join(sourceDir, "GPF5414500000G000_3.dwg"))).toBe(dwgHash);
  expect(sha256(join(sourceDir, "GPF5414500000G000_arc.zip"))).toBe(zipHash);
  expect(sha256(join(sourceDir, "GPF5414500000G000_1.pdf"))).toBe(pdfHash);
  expect(revalidation.public_access.exact_54145_g000_native_dwg_publicly_downloadable).toBeTrue();
  expect(revalidation.drawing_geometry_policy.g001_variant_cad_used).toBeFalse();
  expect(pdfVerification.mechanical_cross_check.pass).toBeTrue();
  expect(pdfVerification.visual_review.status).toBe("codex_visual_qa_pass");
  for (const source of access.official_identity_sources) expect(sha256(join(root, source.local_path))).toBe(source.sha256);
});

test("54145 G000 native DWG extraction partitions all paths repeatably", () => {
  const temporary = mkdtempSync(join(tmpdir(), "gessi316-54145-linework-"));
  const output = join(temporary, "linework.json");
  const run = Bun.spawnSync(["python3", join(root, "pipeline/scripts/gessi316_54145_linework.py"), "--output", output], { cwd: root, stdout: "pipe", stderr: "pipe" });
  if (run.exitCode !== 0) throw new Error(run.stderr.toString());
  const linework = JSON.parse(readFileSync(output, "utf8"));
  expect(linework.source_kind).toBe(sourceKind);
  expect(linework.official_sources.native_dwg.sha256).toBe(dwgHash);
  expect(linework.sheet).toMatchObject({ eligible_entity_count: 1290, all_eligible_entities_partitioned_once: true });
  expect(linework.nominal_dimension_cross_check).toMatchObject({ maximum_ifc_projection_delta_mm: 0.181335, pass: true });
  expect(linework.identity_gates.g001_variant_cad_used).toBeFalse();
  for (const [view, count] of Object.entries(pathCounts)) expect(linework.views[view].paths_mm).toHaveLength(count);
  rmSync(temporary, { recursive: true, force: true });
}, 60_000);

test("single-instance 54145 review has exact G000 blue DWG paths in all three views", () => {
  const manifest = JSON.parse(readFileSync(join(product, "manifest.json"), "utf8"));
  const candidate = JSON.parse(readFileSync(join(product, "candidate-representations.json"), "utf8"));
  const profile = JSON.parse(readFileSync(join(product, "profile.json"), "utf8"));
  const inventory = JSON.parse(readFileSync(join(root, "pipeline/decisions/highpoly-product-inventory.json"), "utf8"));
  expect(inventory.products.find((item: any) => item.type_name === "Gessi316 54145")).toMatchObject({
    representative_global_id: "3jT4sCgpHC98VSIUdGUNYH", folder: "output/review/highpoly-types/gessi316-54145",
    status: "review_ready_pending_approval", source_strategy: "exact_public_official_54145_g000_native_dwg_and_technical_pdf_archived",
  });
  expect(manifest).toMatchObject({ representative_global_id: "3jT4sCgpHC98VSIUdGUNYH", registered_instance_global_ids: ["3jT4sCgpHC98VSIUdGUNYH"], geometry_product_count: 1, whole_model_render: false, formal_ifc_bytes_unchanged: true, article_number: "54145" });
  expect(profile.profiles["gessi316-54145"].official_reference).toMatchObject({ article_number: "54145", configuration: "G000" });
  expect(candidate).toMatchObject({ profile_key: "gessi316-54145", article_number: "54145", source_kind: sourceKind, source_label_zh: sourceLabelZh, source_dwg_sha256: dwgHash, official_cad_used: true, third_party_cad_used: false });
  expect(manifest.views.map((view: any) => view.silhouette_path_count)).toEqual([1, 3, 5]);
  for (const view of manifest.views) {
    expect(view.drawing_line_source_kind).toBe(sourceKind);
    expect(view.official_cad_path_count).toBe(pathCounts[view.view as keyof typeof pathCounts]);
    expect(view.blue_line_present).toBeTrue();
    expect(view.white_mask_present).toBeTrue();
    expect(view.ifc_dwg_compatibility.pass).toBeTrue();
    expect(Math.max(...view.ifc_dwg_compatibility.absolute_delta_mm)).toBeLessThanOrEqual(0.5);
    expect(candidate.views[view.view].official_native_dwg_paths_mm).toHaveLength(pathCounts[view.view as keyof typeof pathCounts]);
    const svg = readFileSync(join(root, view.svg), "utf8");
    expect(svg).toContain('class="simplified-proxy-silhouette geometry-derived"');
    expect(svg).toContain('class="official-reference-mask"');
    expect(svg).toContain('class="official-reference native-dwg"');
    expect(svg).toContain('data-source-kind="native_dwg"');
    expect(svg).toContain("#1677c8");
  }
});

test("54145 project plan and R17 elevations retain context under native-DWG overlays", () => {
  const context = JSON.parse(readFileSync(join(product, "project-context-manifest.json"), "utf8"));
  expect(context).toMatchObject({ source_kind: sourceKind, source_label_zh: sourceLabelZh, official_cad_used: true, third_party_cad_used: false, project_context_retained: true, walls_and_surrounding_project_elements_retained: true, overlay_top_layer_with_white_mask: true, blue_line_present: true, pass: true });
  expect(context.views.map((item: any) => [item.view, item.candidate_view])).toEqual([["plan", "plan"], ["front", "front"], ["side", "side"]]);
  for (const view of context.views) {
    expect(view.overlay.source_dwg_sha256).toBe(dwgHash);
    expect(view.overlay.path_count).toBe(pathCounts[view.candidate_view as keyof typeof pathCounts]);
    expect(view.overlay.fit).toMatchObject({ scale_svg_units_per_mm: 0.02, uniform_scale_preserved: true, transformation_mode: "axis_swap_rigid_reflection_and_translation_only", pass: true });
    expect(Math.max(...view.overlay.fit.bbox_absolute_delta_svg_units)).toBeLessThanOrEqual(view.overlay.fit.bbox_tolerance_svg_units);
    expect(sha256(join(root, view.review_preview))).toBe(view.review_preview_sha256);
    const svg = readFileSync(join(root, view.output), "utf8");
    expect(svg).toContain("IfcWall");
    expect(svg).toContain("IfcSanitaryTerminal");
    if (view.view === "plan") expect(svg).toContain("IfcFurniture");
    expect(svg).toContain('class="official-reference-envelope-mask"');
    expect(svg).toContain('class="official-reference-mask"');
    expect(svg).toContain('class="official-reference native-dwg project-context-overlay"');
  }
});

test("54145 Bonsai evidence contains four actual renders from one isolated IFC Body", () => {
  const manifest = JSON.parse(readFileSync(join(product, "manifest.json"), "utf8"));
  const evidence = JSON.parse(readFileSync(join(product, "bonsai-review-manifest.json"), "utf8"));
  expect(manifest.bonsai_review.manifest_sha256).toBe(sha256(join(product, "bonsai-review-manifest.json")));
  expect(evidence).toMatchObject({ mode: "actual_bonsai_ifc_body_camera_render", representative_global_id: "3jT4sCgpHC98VSIUdGUNYH", geometry_product_count: 1, whole_model_render: false, formal_ifc_bytes_unchanged: true });
  expect(evidence.bonsai_session.saved_active_representation).toBe("Body");
  expect(evidence.bonsai_session.saved_camera_count).toBe(4);
  expect(evidence.renders.map((render: any) => render.view)).toEqual(["plan", "front", "side", "iso"]);
  for (const render of evidence.renders) {
    expect(render.camera_type).toBe("ORTHO");
    expect(render.render_operation).toBe("bpy.ops.render.render(write_still=True)");
    expect(sha256(join(root, render.path))).toBe(render.sha256);
  }
  expect(sha256(join(root, evidence.bonsai_session.path))).toBe(evidence.bonsai_session.sha256);
  expect(sha256(join(root, evidence.isolated_ifc))).toBe(evidence.isolated_ifc_sha256);
});

test("pending 54145 approval leaves formal IFC byte-identical and creates no derived IFC", () => {
  const candidate = JSON.parse(readFileSync(join(product, "candidate-representations.json"), "utf8"));
  const approval = JSON.parse(readFileSync(join(root, "pipeline/decisions/gessi316-54145-drawing-approval.json"), "utf8"));
  expect(candidate).toMatchObject({ formal_ifc_write_allowed: false, review_status: "visual_review_pending" });
  expect(approval).toMatchObject({ status: "pending", approved_views: [], derived_ifc_write_allowed: false, formal_authoritative_ifc_write_allowed: false, scope });
  expect(approval.candidate_manifest_sha256).toBe(sha256(join(product, "manifest.json")));
  expect(existsSync(join(product, "Gessi316-54145-derived-drawing.ifc"))).toBeFalse();
  expect(sha256(formal)).toBe(formalHash);
});

test("54145 writer rejects missing apply and the pending approval record", () => {
  const temporary = mkdtempSync(join(tmpdir(), "gessi316-54145-gate-"));
  const output = join(temporary, "forbidden.ifc");
  const script = join(root, "pipeline/scripts/gessi316_54145_drawing_ifc.py");
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

test("temporary scoped approval writes exact 54145 G000 representations and source links", () => {
  const manifest = join(product, "manifest.json");
  const temporary = mkdtempSync(join(tmpdir(), "gessi316-54145-approved-"));
  const approval = join(temporary, "approval.json");
  const output = join(temporary, "derived.ifc");
  const report = join(temporary, "report.json");
  writeFileSync(approval, JSON.stringify({
    schema_version: 1, profile_key: "gessi316-54145", ifc_type_name: "Gessi316 54145",
    candidate_manifest_sha256: sha256(manifest), status: "approved", reviewer: "automated gate fixture",
    review_date: "2026-08-23", approved_views: ["plan", "front", "side"], derived_ifc_write_allowed: true,
    formal_authoritative_ifc_write_allowed: false, scope, approval_evidence: "temporary automated writer verification only",
  }, null, 2));
  const run = Bun.spawnSync([
    "python3", join(root, "pipeline/scripts/gessi316_54145_drawing_ifc.py"), "--input", formal,
    "--manifest", manifest, "--approval", approval, "--output", output, "--report", report, "--apply",
  ], { cwd: root, stdout: "pipe", stderr: "pipe" });
  if (run.exitCode !== 0) throw new Error(run.stderr.toString());
  const result = JSON.parse(readFileSync(report, "utf8"));
  expect(result).toMatchObject({
    pass: true, formal_ifc_bytes_unchanged: true,
    representations: { plan: "Gessi54145Plan", front: "Gessi54145Front", side: "Gessi54145Side" },
    representation_path_counts: pathCounts, source_kind: sourceKind, source_label_zh: sourceLabelZh,
    official_cad_geometry_included: true, source_property_set: "Pset_Gessi31654145DrawingSource",
  });
  expect(result.representation_geometry_source).toBe("official_native_dwg_paths_mm extracted from GPF5414500000G000_3.dwg");
  expect(result.source_document_associations).toHaveLength(6);
  expect(result.source_document_associations).toContain("GESSI316-54145-OFFICIAL-NATIVE-DWG-ZIP");
  expect(result.source_document_associations).toContain("GESSI316-54145-OFFICIAL-TECHNICAL-PDF");
  const inspect = Bun.spawnSync(["python3", "-c", [
    "import ifcopenshell, ifcopenshell.util.element, json, sys", "m=ifcopenshell.open(sys.argv[1])",
    "p=m.by_guid('3jT4sCgpHC98VSIUdGUNYH')", "print(json.dumps(ifcopenshell.util.element.get_pset(p, 'Pset_Gessi31654145DrawingSource')))"
  ].join(";"), output], { cwd: root, stdout: "pipe", stderr: "pipe" });
  if (inspect.exitCode !== 0) throw new Error(inspect.stderr.toString());
  expect(JSON.parse(inspect.stdout.toString())).toMatchObject({
    ArticleNumber: "54145", Configuration: "G000", TechnicalDrawingNumber: "GPF5414500000G000",
    SourceDwgSha256: dwgHash, OfficialProductCadStatus: "public_official_api_exact_54145_g000_native_dwg_acquired",
    OfficialCadUsed: "true", OfficialVectorEvidenceUsedAsCadGeometry: "false", ExcludedVariant: "G001",
    DimensionCrossCheckPass: "true", EvidenceScope: scope,
  });
  expect(sha256(formal)).toBe(formalHash);
  rmSync(temporary, { recursive: true, force: true });
}, 120_000);
