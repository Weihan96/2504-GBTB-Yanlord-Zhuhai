import { expect, test } from "bun:test";
import { createHash } from "node:crypto";
import { existsSync, mkdtempSync, readFileSync, rmSync, writeFileSync } from "node:fs";
import { tmpdir } from "node:os";
import { dirname, join, resolve } from "node:path";

const root = resolve(import.meta.dir, "../..");
const product = join(root, "output/review/highpoly-types/gessi316-54294");
const sourceDir = join(product, "official-source");
const formal = join(root, "2504 GBTB Yanlord Zhuhai.ifc");
const formalHash = "7a50b87e8f48a7c2bfbaa9f1dabdfca7155fa684b2325c66aa1ae15c8ab0a25c";
const sourceKind = "native_dwg";
const sourceLabelZh = "Gessi 官方精确型号 54294 原生 DWG 图纸表达";
const dwgHash = "dfe8bcf7e100c14187c387b75a35e3fe7f718c61595fadf66adfe92ca7648ff4";
const zipHash = "fad0d98e94a83bb488b5f0703472862d628759f21700c1a87c3af343f4c1bdb9";
const pdfHash = "82058bb27752f27e6fc92029783d67c15fd443befcd393f7d572fd0147a581e3";
const pathCounts = { plan: 1481, front: 1099, side: 745 };
const scope = "official Gessi exact 54294 family reference and 45089_54294 article combination; not a project shop drawing";

function sha256(path: string) {
  return createHash("sha256").update(readFileSync(path)).digest("hex");
}

test("exact public Gessi 54294 DWG and PDF sources are archived and cross-verified", () => {
  const access = JSON.parse(readFileSync(join(sourceDir, "source-access-record.json"), "utf8"));
  const revalidation = JSON.parse(readFileSync(join(sourceDir, "official-source-revalidation.json"), "utf8"));
  const pdfVerification = JSON.parse(readFileSync(join(sourceDir, "official-pdf-verification.json"), "utf8"));
  expect(access).toMatchObject({
    manufacturer: "Gessi",
    project_ifc_type_name: "Gessi316 54294",
    resolved_article_number: "45089_54294",
    external_visible_product_code: "54294",
    companion_built_in_product_code: "45089",
    scope,
    pass: true,
  });
  expect(access.official_product_cad).toMatchObject({
    authentication_required: false,
    acquired: true,
    exact_project_configuration_match: true,
  });
  expect(access.drawing_geometry_source).toMatchObject({
    source_kind: sourceKind,
    source_label_zh: sourceLabelZh,
    official_cad_used: true,
    third_party_cad_used: false,
    adjacent_product_cad_used: false,
    source_dwg_sha256: dwgHash,
  });
  expect(access.identity_and_geometry_policy).toMatchObject({
    "54294_native_dwg_used_for_three_views": true,
    "45089_companion_dwg_used_as_54294_geometry": false,
    adjacent_gessi_product_cad_used: false,
    third_party_cad_used: false,
  });
  expect(sha256(join(sourceDir, "GPF5429400000G000_3.dwg"))).toBe(dwgHash);
  expect(sha256(join(sourceDir, "GPF5429400000G000_arc.zip"))).toBe(zipHash);
  expect(sha256(join(sourceDir, "GPF5429400000G000_1.pdf"))).toBe(pdfHash);
  expect(revalidation).toMatchObject({
    all_downloaded_bytes_match_local_archive: true,
    all_zip_members_match_extracted_files: true,
    pass: true,
  });
  expect(revalidation.zip_members["54294"].member_sha256).toBe(dwgHash);
  expect(revalidation.drawing_geometry_policy.companion_45089_dwg_used_as_54294_geometry).toBeFalse();
  expect(pdfVerification.pass).toBeTrue();
  expect(pdfVerification.mechanical_cross_check.pass).toBeTrue();
  expect(pdfVerification.visual_review.status).toBe("codex_visual_qa_pass");
});

test("native DWG linework extraction is repeatable and dimensionally exact", () => {
  const temporary = mkdtempSync(join(tmpdir(), "gessi316-54294-linework-"));
  const output = join(temporary, "linework.json");
  const run = Bun.spawnSync([
    "python3", join(root, "pipeline/scripts/gessi316_54294_linework.py"), "--output", output,
  ], { cwd: root, stdout: "pipe", stderr: "pipe" });
  if (run.exitCode !== 0) throw new Error(run.stderr.toString());
  const linework = JSON.parse(readFileSync(output, "utf8"));
  expect(linework.source_kind).toBe(sourceKind);
  expect(linework.official_sources.native_dwg.sha256).toBe(dwgHash);
  expect(linework.official_sources.native_dwg_zip.sha256).toBe(zipHash);
  expect(linework.official_sources.technical_vector_pdf.sha256).toBe(pdfHash);
  expect(linework.identity_gates).toMatchObject({
    "54294_native_dwg_used": true,
    "45089_companion_dwg_used_as_54294_geometry": false,
    adjacent_gessi_product_cad_used: false,
    third_party_cad_used: false,
  });
  expect(linework.nominal_dimension_cross_check).toMatchObject({
    native_dwg_plan_envelope_mm: [361.7, 209.812412],
    absolute_delta_mm: [0.3, 0.187588],
    tolerance_mm: 0.5,
    pass: true,
  });
  for (const [view, count] of Object.entries(pathCounts)) {
    expect(linework.views[view].source_dwg_sha256).toBe(dwgHash);
    expect(linework.views[view].paths_mm).toHaveLength(count);
  }
  rmSync(temporary, { recursive: true, force: true });
}, 60_000);

test("single-instance three-view SVGs place exact native DWG blue linework above IFC comparisons", () => {
  const manifest = JSON.parse(readFileSync(join(product, "manifest.json"), "utf8"));
  const candidate = JSON.parse(readFileSync(join(product, "candidate-representations.json"), "utf8"));
  expect(manifest.representative_global_id).toBe("2iKOL78$H0N9Yd9$ky3pW4");
  expect(manifest.registered_instance_global_ids).toEqual(["2iKOL78$H0N9Yd9$ky3pW4"]);
  expect(manifest.geometry_product_count).toBe(1);
  expect(manifest.whole_model_render).toBeFalse();
  expect(manifest.formal_ifc_bytes_unchanged).toBeTrue();
  expect(manifest.drawing_source).toMatchObject({ source_kind: sourceKind, source_dwg_sha256: dwgHash, official_cad_used: true, third_party_cad_used: false });
  expect(manifest.companion_45089_dwg_used_as_54294_geometry).toBeFalse();
  expect(candidate).toMatchObject({ source_kind: sourceKind, source_dwg_sha256: dwgHash, official_cad_used: true, third_party_cad_used: false, review_status: "visual_review_pending" });
  for (const view of manifest.views) {
    expect(view.drawing_line_source_kind).toBe(sourceKind);
    expect(view.official_cad_path_count).toBe(pathCounts[view.view as keyof typeof pathCounts]);
    expect(view.blue_line_present).toBeTrue();
    expect(view.white_mask_present).toBeTrue();
    expect(view.ifc_dwg_compatibility.pass).toBeTrue();
    expect(candidate.views[view.view].official_native_dwg_paths_mm).toHaveLength(pathCounts[view.view as keyof typeof pathCounts]);
    const svg = readFileSync(join(root, view.svg), "utf8");
    expect(svg).toContain('class="original-highpoly"');
    expect(svg).toContain('class="simplified-proxy-silhouette geometry-derived"');
    expect(svg).toContain('class="official-reference-mask"');
    expect(svg).toContain('class="official-reference native-dwg"');
    expect(svg).toContain("#1677c8");
    expect(svg.indexOf('class="official-reference-mask"')).toBeLessThan(svg.indexOf('class="official-reference native-dwg"'));
  }
});

test("project drawings retain walls and furniture while exact blue Gessi CAD stays topmost", () => {
  const context = JSON.parse(readFileSync(join(product, "project-context-manifest.json"), "utf8"));
  expect(context).toMatchObject({ source_kind: sourceKind, official_cad_used: true, third_party_cad_used: false, project_context_retained: true, walls_and_surrounding_project_elements_retained: true, overlay_top_layer_with_white_mask: true, blue_line_present: true, pass: true });
  expect(context.views.map((view: any) => [view.view, view.candidate_view])).toEqual([["plan", "plan"], ["front", "side"], ["side", "front"]]);
  for (const view of context.views) {
    expect(view.overlay.source_dwg_sha256).toBe(dwgHash);
    expect(view.overlay.path_count).toBe(pathCounts[view.candidate_view as keyof typeof pathCounts]);
    expect(view.overlay.fit.uniform_scale_preserved).toBeTrue();
    expect(Math.max(...view.overlay.fit.bbox_absolute_delta_svg_units)).toBeLessThanOrEqual(view.overlay.fit.bbox_tolerance_svg_units);
    expect(sha256(join(root, view.review_preview))).toBe(view.review_preview_sha256);
    const svgPath = join(root, view.output);
    const svg = readFileSync(svgPath, "utf8");
    expect(svg).toContain("IfcWall");
    expect(svg).toContain("IfcSanitaryTerminal");
    expect(svg).toContain('class="official-reference-envelope-mask"');
    expect(svg).toContain('class="official-reference-mask"');
    expect(svg).toContain('class="official-reference native-dwg project-context-overlay"');
    const overlay = svg.slice(svg.indexOf('class="official-reference native-dwg project-context-overlay"'));
    expect(overlay.indexOf('class="official-reference-envelope-mask"')).toBeLessThan(overlay.indexOf('class="official-reference-mask"'));
    expect(overlay.indexOf('class="official-reference-mask"')).toBeLessThan(overlay.indexOf('class="official-reference native-dwg"'));
    for (const match of svg.matchAll(/(?:href|xlink:href)="([^"#]+)"/g)) {
      if (/^(?:https?:|data:)/.test(match[1])) continue;
      expect(existsSync(resolve(dirname(svgPath), match[1]))).toBeTrue();
    }
  }
});

test("Bonsai evidence is actual isolated IFC Body rendering from four saved cameras", () => {
  const evidence = JSON.parse(readFileSync(join(product, "bonsai-review-manifest.json"), "utf8"));
  expect(evidence).toMatchObject({ mode: "actual_bonsai_ifc_body_camera_render", representative_global_id: "2iKOL78$H0N9Yd9$ky3pW4", geometry_product_count: 1, whole_model_render: false, formal_ifc_bytes_unchanged: true });
  expect(evidence.bonsai_session.saved_active_representation).toBe("Body");
  expect(evidence.bonsai_session.saved_camera_count).toBe(4);
  expect(evidence.renders.map((render: any) => render.view)).toEqual(["plan", "front", "side", "iso"]);
  for (const render of evidence.renders) {
    expect(render.camera_type).toBe("ORTHO");
    expect(render.render_operation).toBe("bpy.ops.render.render(write_still=True)");
    expect(sha256(join(root, render.path))).toBe(render.sha256);
  }
  expect(sha256(join(root, evidence.isolated_ifc))).toBe(evidence.isolated_ifc_sha256);
});

test("pending review cannot write an IFC and formal IFC remains byte-identical", () => {
  const temporary = mkdtempSync(join(tmpdir(), "gessi316-54294-gate-"));
  const output = join(temporary, "forbidden.ifc");
  const script = join(root, "pipeline/scripts/gessi316_54294_drawing_ifc.py");
  const approval = JSON.parse(readFileSync(join(root, "pipeline/decisions/gessi316-54294-drawing-approval.json"), "utf8"));
  const manifest = join(product, "manifest.json");
  expect(approval).toMatchObject({ status: "pending", derived_ifc_write_allowed: false, formal_authoritative_ifc_write_allowed: false, scope });
  expect(approval.candidate_manifest_sha256).toBe(sha256(manifest));
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

test("temporary scoped approval writes open native-DWG representations and source associations", () => {
  const manifest = join(product, "manifest.json");
  const temporary = mkdtempSync(join(tmpdir(), "gessi316-54294-approved-"));
  const approvalPath = join(temporary, "approval.json");
  const output = join(temporary, "derived.ifc");
  const report = join(temporary, "report.json");
  writeFileSync(approvalPath, JSON.stringify({
    schema_version: 1, profile_key: "gessi316-54294", ifc_type_name: "Gessi316 54294",
    candidate_manifest_sha256: sha256(manifest), status: "approved", reviewer: "automated gate fixture",
    review_date: "2026-08-23", approved_views: ["plan", "front", "side"], derived_ifc_write_allowed: true,
    formal_authoritative_ifc_write_allowed: false, scope, approval_evidence: "temporary automated writer verification only",
  }, null, 2));
  const run = Bun.spawnSync(["python3", join(root, "pipeline/scripts/gessi316_54294_drawing_ifc.py"), "--input", formal, "--manifest", manifest, "--approval", approvalPath, "--output", output, "--report", report, "--apply"], { cwd: root, stdout: "pipe", stderr: "pipe" });
  if (run.exitCode !== 0) throw new Error(run.stderr.toString());
  const result = JSON.parse(readFileSync(report, "utf8"));
  expect(result.pass).toBeTrue();
  expect(result.formal_ifc_bytes_unchanged).toBeTrue();
  expect(result.representation_path_counts).toEqual(pathCounts);
  expect(result.source_kind).toBe(sourceKind);
  expect(result.source_label_zh).toBe(sourceLabelZh);
  expect(result.official_cad_geometry_included).toBeTrue();
  expect(result.representation_geometry_source).toBe("official_native_dwg_paths_mm from GPF5429400000G000_3.dwg");
  expect(result.source_document_associations).toContain("GESSI316-45089-54294-OFFICIAL-NATIVE-DWG");
  expect(result.source_document_associations).toContain("GESSI316-45089-54294-OFFICIAL-TECHNICAL-PDF");
  expect(result.source_document_associations).toContain("GESSI316-45089-54294-SOURCE-ACCESS-RECORD");
  expect(sha256(formal)).toBe(formalHash);
  rmSync(temporary, { recursive: true, force: true });
}, 120_000);
