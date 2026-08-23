import { expect, test } from "bun:test";
import { createHash } from "node:crypto";
import { existsSync, mkdtempSync, readFileSync, rmSync, writeFileSync } from "node:fs";
import { tmpdir } from "node:os";
import { join, resolve } from "node:path";

const root = resolve(import.meta.dir, "../..");
const product = join(root, "output/review/highpoly-types/gessi316-54146");
const sourceDir = join(product, "official-source");
const formal = join(root, "2504 GBTB Yanlord Zhuhai.ifc");
const formalHash = "7a50b87e8f48a7c2bfbaa9f1dabdfca7155fa684b2325c66aa1ae15c8ab0a25c";
const sourceKind = "native_dwg";
const sourceLabelZh = "Gessi 官方精确型号 54146 G000 原生 DWG 图纸表达";
const scope = "official Gessi exact 54146 G000 family reference; not a project shop drawing";
const dwgHash = "c8ddb90f61565d5273a32574f777812bbb1d9ef33e2b98479f1e7c381e69950d";
const zipHash = "aaf9964d57674a03b7e7a8f9eb647df84af74bab1b575c9bad306d83bd40237e";
const pdfHash = "ee7909886b1181306902d0a0e146b73f938bf96e3e12e8f1ba27ded996b2b8ce";
const pathCounts = { plan: 16, front: 607, side: 660 };

function sha256(path: string) {
  return createHash("sha256").update(readFileSync(path)).digest("hex");
}

test("exact public Gessi 54146 G000 DWG and PDF are archived and cross-verified", () => {
  const access = JSON.parse(readFileSync(join(sourceDir, "source-access-record.json"), "utf8"));
  const revalidation = JSON.parse(readFileSync(join(sourceDir, "official-source-revalidation.json"), "utf8"));
  const pdfVerification = JSON.parse(readFileSync(join(sourceDir, "official-pdf-verification.json"), "utf8"));
  expect(access).toMatchObject({
    manufacturer: "Gessi",
    family: "Gessi316 Meccanica",
    project_ifc_type_name: "Gessi316 54146",
    resolved_article_number: "54146",
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
    official_api_width_depth_height_mm: [300, 300, 276],
    native_dwg_plan_envelope_mm: [300, 300],
    native_dwg_front_envelope_mm: [300, 280.25],
    native_dwg_side_envelope_mm: [300, 280.25],
    project_ifc_body_local_xyz_mm: [300, 299.944916, 279.371241],
    maximum_ifc_projection_delta_mm: 0.878759,
    tolerance_mm: 1,
    pass: true,
    geometry_stretched: false,
  });
  expect(access.drawing_geometry_source).toMatchObject({
    source_kind: sourceKind,
    source_label_zh: sourceLabelZh,
    source_dwg_sha256: dwgHash,
    official_cad_used: true,
    third_party_cad_used: false,
  });
  expect(access.identity_and_geometry_policy).toMatchObject({ "54146_g000_native_dwg_used": true, g001_variant_cad_used: false, wall_mounted_54145_cad_used: false, third_party_cad_used: false });
  expect(sha256(join(sourceDir, "GPF5414600000G000_3.dwg"))).toBe(dwgHash);
  expect(sha256(join(sourceDir, "GPF5414600000G000_arc.zip"))).toBe(zipHash);
  expect(sha256(join(sourceDir, "GPF5414600000G000_1.pdf"))).toBe(pdfHash);
  expect(revalidation.public_access.exact_54146_g000_native_dwg_publicly_downloadable).toBeTrue();
  expect(revalidation.drawing_geometry_policy.wall_mounted_54145_cad_used).toBeFalse();
  expect(pdfVerification.mechanical_cross_check.pass).toBeTrue();
  expect(pdfVerification.visual_review.status).toBe("codex_visual_qa_pass");
  for (const evidence of access.official_identity_sources) {
    expect(sha256(join(root, evidence.local_path))).toBe(evidence.sha256);
  }
});

test("54146 G000 native DWG extraction partitions all paths repeatably", () => {
  const temporary = mkdtempSync(join(tmpdir(), "gessi316-54146-linework-"));
  const output = join(temporary, "linework.json");
  const run = Bun.spawnSync(["python3", join(root, "pipeline/scripts/gessi316_54146_linework.py"), "--output", output], { cwd: root, stdout: "pipe", stderr: "pipe" });
  if (run.exitCode !== 0) throw new Error(run.stderr.toString());
  const linework = JSON.parse(readFileSync(output, "utf8"));
  expect(linework).toMatchObject({ source_kind: sourceKind, sheet: { eligible_entity_count: 1283, all_eligible_entities_partitioned_once: true } });
  expect(linework.official_sources.native_dwg.sha256).toBe(dwgHash);
  expect(linework.nominal_dimension_cross_check).toMatchObject({ maximum_ifc_projection_delta_mm: 0.878759, pass: true });
  expect(linework.identity_gates.wall_mounted_54145_cad_used).toBeFalse();
  for (const [view, count] of Object.entries(pathCounts)) expect(linework.views[view].paths_mm).toHaveLength(count);
  rmSync(temporary, { recursive: true, force: true });
}, 60_000);

test("one Gessi 54146 representative supplies three exact G000 native-DWG views", () => {
  const manifest = JSON.parse(readFileSync(join(product, "manifest.json"), "utf8"));
  const candidate = JSON.parse(readFileSync(join(product, "candidate-representations.json"), "utf8"));
  const inventory = JSON.parse(readFileSync(join(root, "pipeline/decisions/highpoly-product-inventory.json"), "utf8"));
  const entry = inventory.products.find((item: any) => item.type_name === "Gessi316 54146");
  expect(entry).toMatchObject({
    folder: "output/review/highpoly-types/gessi316-54146",
    status: "review_ready_pending_approval",
    representative_global_id: "04DLh1Jk9Dcu9ibcaE0id8",
    source_strategy: "exact_public_official_54146_g000_native_dwg_and_technical_pdf_archived",
  });
  expect(inventory.review_ready_pending_approval_product_count).toBe(38);
  expect(manifest.representative_global_id).toBe("04DLh1Jk9Dcu9ibcaE0id8");
  expect(manifest.registered_instance_global_ids).toEqual(["04DLh1Jk9Dcu9ibcaE0id8"]);
  expect(manifest.geometry_product_count).toBe(1);
  expect(manifest.whole_model_render).toBeFalse();
  expect(manifest.bounds_mm.size).toEqual([300, 299.944916, 279.371241]);
  expect(manifest.views.map((item: any) => item.silhouette_path_count)).toEqual([1, 20, 23]);
  expect(candidate).toMatchObject({ source_kind: sourceKind, source_label_zh: sourceLabelZh, source_dwg_sha256: dwgHash, official_cad_used: true, third_party_cad_used: false });
  for (const view of manifest.views) {
    expect(view.official_cad_path_count).toBe(pathCounts[view.view as keyof typeof pathCounts]);
    expect(view.blue_line_present).toBeTrue();
    expect(view.white_mask_present).toBeTrue();
    expect(view.ifc_dwg_compatibility.pass).toBeTrue();
    expect(Math.max(...view.ifc_dwg_compatibility.absolute_delta_mm)).toBeLessThanOrEqual(1);
    expect(candidate.views[view.view].official_native_dwg_paths_mm).toHaveLength(pathCounts[view.view as keyof typeof pathCounts]);
    const svg = readFileSync(join(root, view.svg), "utf8");
    expect(svg).toContain('class="simplified-proxy-silhouette geometry-derived"');
    expect(svg).toContain('class="official-reference-mask"');
    expect(svg).toContain('class="official-reference native-dwg"');
    expect(svg).toContain("#1677c8");
  }
});

test("project plan and elevations retain walls, furniture and surrounding fixtures", () => {
  const context = JSON.parse(readFileSync(join(product, "project-context-manifest.json"), "utf8"));
  expect(context).toMatchObject({
    project_context_retained: true,
    walls_and_surrounding_project_elements_retained: true,
    overlay_top_layer_with_white_mask: true,
    source_kind: sourceKind,
    source_label_zh: sourceLabelZh,
    official_cad_used: true,
    blue_line_present: true,
    pass: true,
  });
  expect(context.views.map((item: any) => [item.view, item.candidate_view])).toEqual([
    ["plan", "plan"], ["front", "front"], ["side", "side"],
  ]);
  for (const view of context.views) {
    expect(view.overlay.fit.pass).toBeTrue();
    expect(view.overlay.source_dwg_sha256).toBe(dwgHash);
    expect(view.overlay.path_count).toBe(pathCounts[view.candidate_view as keyof typeof pathCounts]);
    expect(Math.max(...view.overlay.fit.bbox_absolute_delta_svg_units)).toBeLessThanOrEqual(0.08);
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

test("actual Bonsai IFC Body render saves four orthographic cameras", () => {
  const evidence = JSON.parse(readFileSync(join(product, "bonsai-review-manifest.json"), "utf8"));
  expect(evidence.mode).toBe("actual_bonsai_ifc_body_camera_render");
  expect(evidence.representative_global_id).toBe("04DLh1Jk9Dcu9ibcaE0id8");
  expect(evidence.geometry_product_count).toBe(1);
  expect(evidence.whole_model_render).toBeFalse();
  expect(evidence.bonsai_session.saved_active_representation).toBe("Body");
  expect(evidence.bonsai_session.saved_camera_count).toBe(4);
  expect(evidence.renders.map((item: any) => item.view)).toEqual(["plan", "front", "side", "iso"]);
  for (const render of evidence.renders) {
    expect(render.camera_type).toBe("ORTHO");
    expect(render.render_operation).toBe("bpy.ops.render.render(write_still=True)");
    expect(sha256(join(root, render.path))).toBe(render.sha256);
  }
  expect(sha256(join(root, evidence.bonsai_session.path))).toBe(evidence.bonsai_session.sha256);
  expect(sha256(join(root, evidence.isolated_ifc))).toBe(evidence.isolated_ifc_sha256);
});

test("pending Gessi 54146 review cannot write an IFC", () => {
  const approval = JSON.parse(readFileSync(join(root, "pipeline/decisions/gessi316-54146-drawing-approval.json"), "utf8"));
  const manifest = join(product, "manifest.json");
  expect(approval.status).toBe("pending");
  expect(approval.derived_ifc_write_allowed).toBeFalse();
  expect(approval.formal_authoritative_ifc_write_allowed).toBeFalse();
  expect(approval.candidate_manifest_sha256).toBe(sha256(manifest));
  expect(existsSync(join(product, "Gessi316-54146-derived-drawing.ifc"))).toBeFalse();
  expect(sha256(formal)).toBe(formalHash);
  const temporary = mkdtempSync(join(tmpdir(), "gessi316-54146-gate-"));
  const output = join(temporary, "forbidden.ifc");
  const script = join(root, "pipeline/scripts/gessi316_54146_drawing_ifc.py");
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

test("temporary approval writes verified representations and source associations", () => {
  const manifest = join(product, "manifest.json");
  const temporary = mkdtempSync(join(tmpdir(), "gessi316-54146-approved-"));
  const approval = join(temporary, "approval.json");
  const output = join(temporary, "derived.ifc");
  const report = join(temporary, "report.json");
  writeFileSync(approval, JSON.stringify({
    schema_version: 1,
    profile_key: "gessi316-54146",
    ifc_type_name: "Gessi316 54146",
    candidate_manifest_sha256: sha256(manifest),
    status: "approved",
    reviewer: "automated gate fixture",
    review_date: "2026-08-24",
    approved_views: ["plan", "front", "side"],
    derived_ifc_write_allowed: true,
    formal_authoritative_ifc_write_allowed: false,
    scope,
    approval_evidence: "temporary automated writer verification only",
  }, null, 2));
  const run = Bun.spawnSync([
    "python3", join(root, "pipeline/scripts/gessi316_54146_drawing_ifc.py"),
    "--input", formal, "--manifest", manifest, "--approval", approval,
    "--output", output, "--report", report, "--apply",
  ], { cwd: root, stdout: "pipe", stderr: "pipe" });
  expect(run.exitCode).toBe(0);
  if (run.exitCode !== 0) throw new Error(run.stderr.toString());
  const result = JSON.parse(readFileSync(report, "utf8"));
  expect(result.pass).toBeTrue();
  expect(result.formal_ifc_bytes_unchanged).toBeTrue();
  expect(result.representations).toEqual({
    plan: "Gessi54146Plan", front: "Gessi54146Front", side: "Gessi54146Side",
  });
  expect(result.representation_path_counts).toEqual(pathCounts);
  expect(result.source_kind).toBe(sourceKind);
  expect(result.source_label_zh).toBe(sourceLabelZh);
  expect(result.official_cad_geometry_included).toBeTrue();
  expect(result.representation_geometry_source).toBe("official_native_dwg_paths_mm extracted from GPF5414600000G000_3.dwg");
  expect(result.source_document_associations).toEqual([
    "GESSI316-54146-OFFICIAL-CATALOGUE",
    "GESSI316-54146-OFFICIAL-NATIVE-DWG-ZIP",
    "GESSI316-54146-OFFICIAL-PRODUCT-API",
    "GESSI316-54146-OFFICIAL-PRODUCT-ROUTE",
    "GESSI316-54146-OFFICIAL-TECHNICAL-PDF",
    "GESSI316-54146-SOURCE-ACCESS-RECORD",
  ]);
  const inspect = Bun.spawnSync(["python3", "-c", [
    "import ifcopenshell, ifcopenshell.util.element, json, sys", "m=ifcopenshell.open(sys.argv[1])",
    "p=m.by_guid('04DLh1Jk9Dcu9ibcaE0id8')", "print(json.dumps(ifcopenshell.util.element.get_pset(p, 'Pset_Gessi31654146DrawingSource')))"
  ].join(";"), output], { cwd: root, stdout: "pipe", stderr: "pipe" });
  if (inspect.exitCode !== 0) throw new Error(inspect.stderr.toString());
  expect(JSON.parse(inspect.stdout.toString())).toMatchObject({
    ArticleNumber: "54146", Configuration: "G000", TechnicalDrawingNumber: "GPF5414600000G000",
    SourceDwgSha256: dwgHash, OfficialProductCadStatus: "public_official_api_exact_54146_g000_native_dwg_acquired",
    OfficialCadUsed: "true", OfficialVectorEvidenceUsedAsCadGeometry: "false", ExcludedVariant: "G001",
    DimensionCrossCheckPass: "true", EvidenceScope: scope,
  });
  expect(sha256(formal)).toBe(formalHash);
  rmSync(temporary, { recursive: true, force: true });
}, 120_000);
