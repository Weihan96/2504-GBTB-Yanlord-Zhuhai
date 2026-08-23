import { expect, test } from "bun:test";
import { createHash } from "node:crypto";
import { existsSync, mkdtempSync, readFileSync, rmSync, writeFileSync } from "node:fs";
import { tmpdir } from "node:os";
import { join, resolve } from "node:path";

const root = resolve(import.meta.dir, "../..");
const product = join(root, "output/review/highpoly-types/geberit-154-446-ks-1");
const formal = join(root, "2504 GBTB Yanlord Zhuhai.ifc");
const formalHash = "7a50b87e8f48a7c2bfbaa9f1dabdfca7155fa684b2325c66aa1ae15c8ab0a25c";
const productPageArchiveHash = "cf4578d1f1d60078c683548281f22e30fb71608f7d4413a21c840bde212fa02f";
const sourceHashes: Record<string, string> = {
  A: "de09a6bbf99ecbe3c832cbee7ac9d398e95b927246f274643b513bfe960586c3",
  G: "bc41f8db7c7989de9d9089374e03c3b9deea94e2634bbb0f9418cad982b4b6e4",
  L: "dedb47964cc69310bd9af3982379c428c8b20494b7588296aedc5d07380e3421",
  P: "46e1d9fe9b7afda6f4a72970a5e35d7bc0935c2bdc890ccca232f0b0896a8e82",
};

function sha256(path: string) {
  return createHash("sha256").update(readFileSync(path)).digest("hex");
}

test("exact archived CleanLine50 DWGs retain the old article identity", () => {
  const profile = JSON.parse(readFileSync(join(product, "profile.json"), "utf8"));
  const reference = profile.profiles["geberit-154-446-ks-1"].official_reference;
  expect(reference.current_catalog_status).toContain("154.446.KS.2 replaces legacy 154.446.KS.1");
  expect(reference.replacement_cad_used).toBeFalse();
  expect(sha256(join(root, reference.product_page_archive))).toBe(productPageArchiveHash);
  const linework = JSON.parse(readFileSync(join(product, "official-native-dwg-linework.json"), "utf8"));
  expect(linework.source_kind).toBe("native_dwg");
  expect(linework.article_number).toBe("154.446.KS.1");
  expect(linework.replacement_article_number).toBe("154.446.KS.2");
  expect(linework.retirement_identity.replacement_cad_used).toBeFalse();
  expect(linework.ifc_type_name).toBe("Geberit 154.446.KS.1");
  expect(linework.pass).toBeTrue();
  expect(linework.three_view_code_mapping).toMatchObject({
    G: "Grundriss / plan",
    A: "Ansicht / front elevation",
    L: "left side elevation",
  });
  for (const [code, expected] of Object.entries(sourceHashes)) {
    const source = linework.official_sources[code];
    expect(source.sha256).toBe(expected);
    expect(source.url).toBe(`https://cdn.data.geberit.com/cad/154.446.KS.1_${code}.dwg`);
    expect(sha256(join(root, source.path))).toBe(expected);
  }
  const access = JSON.parse(
    readFileSync(join(product, "official-source/source-access-record.json"), "utf8"),
  );
  expect(access).toMatchObject({
    manufacturer: "Geberit",
    article_number: "154.446.KS.1",
    replacement_article_number: "154.446.KS.2",
    replacement_cad_used: false,
    source_kind: "native_dwg",
    official_cad_acquired: true,
    official_cad_exact_archived_article_match: true,
    official_cad_used: true,
    third_party_cad_used: false,
    blue_line_present: true,
    representation_geometry_source: "official_native_dwg_paths_mm",
    pass: true,
  });
  expect(access.drawing_view_mapping).toEqual({ plan: "G", front: "A", side: "L" });
  expect(access.official_3d_identity_only.used_as_plan_or_elevation_geometry).toBeFalse();
  expect(access.official_cdn_revalidation).toMatchObject({
    all_declared_urls_accessible: true,
    all_downloaded_bytes_match_local_archive: true,
    replacement_cad_used: false,
    pass: true,
  });
  expect(sha256(join(root, access.official_cdn_revalidation.path))).toBe(
    access.official_cdn_revalidation.sha256,
  );
  for (const [code, expected] of Object.entries(sourceHashes)) {
    expect(access.official_native_dwg[code].sha256).toBe(expected);
    expect(sha256(join(root, access.official_native_dwg[code].path))).toBe(expected);
  }
});

test("official CDN revalidation proves exact old-article bytes remain accessible", () => {
  const evidence = JSON.parse(
    readFileSync(join(product, "official-source/official-cdn-revalidation.json"), "utf8"),
  );
  expect(evidence).toMatchObject({
    manufacturer: "Geberit",
    article_number: "154.446.KS.1",
    replacement_article_number: "154.446.KS.2",
    replacement_cad_used: false,
    drawing_view_mapping: { plan: "G", front: "A", side: "L" },
    identity_only_code: "P",
    all_declared_urls_accessible: true,
    all_downloaded_bytes_match_local_archive: true,
    pass: true,
  });
  for (const [code, expected] of Object.entries(sourceHashes)) {
    expect(evidence.results[code]).toMatchObject({
      url: `https://cdn.data.geberit.com/cad/154.446.KS.1_${code}.dwg`,
      http_status: 200,
      expected_sha256: expected,
      downloaded_sha256: expected,
      local_sha256: expected,
      downloaded_bytes_match_local_archive: true,
      pass: true,
    });
  }
});

test("official G/A/L views mechanically cross-check the isolated IFC Body", () => {
  const manifest = JSON.parse(readFileSync(join(product, "manifest.json"), "utf8"));
  expect(manifest).toMatchObject({
    source_kind: "native_dwg",
    official_cad_acquired: true,
    official_cad_exact_archived_article_match: true,
    official_cad_used: true,
    third_party_cad_used: false,
    blue_line_present: true,
    representation_geometry_source: "official_native_dwg_paths_mm",
    simplified_proxy_comparison_included: true,
  });
  expect(sha256(join(root, manifest.official_source_access_record))).toBe(
    manifest.official_source_access_record_sha256,
  );
  expect(manifest.geometry_product_count).toBe(1);
  expect(manifest.whole_model_render).toBeFalse();
  expect(manifest.formal_ifc_bytes_unchanged).toBeTrue();
  expect(manifest.review_status).toBe("visual_review_pending");
  expect(manifest.approved_for_drawing_ifc).toBeFalse();
  expect(manifest.official_reference.source_kind).toBe("native_dwg");
  expect(manifest.official_reference.article_number).toBe("154.446.KS.1");
  expect(manifest.official_reference.replacement_cad_used).toBeFalse();
  expect(manifest.views.map((view: any) => view.blue_line_native_dwg_code)).toEqual(["G", "A", "L"]);
  expect(manifest.views.map((view: any) => view.official_reference_path_count)).toEqual([122, 118, 176]);
  for (const view of manifest.views) {
    expect(view.blue_line_source_kind).toBe("native_dwg");
    expect(view.blue_line_article_number).toBe("154.446.KS.1");
    expect(view.mechanical_cross_check.pass).toBeTrue();
    const svg = readFileSync(join(root, view.svg), "utf8");
    expect(svg).toContain('class="official-reference-mask"');
    expect(svg).toContain('class="official-reference native-dwg"');
    expect(svg).toContain('data-source-kind="native_dwg"');
    expect((svg.match(/<path class="official-reference native-dwg"/g) ?? []).length).toBe(1);
    const officialPath = svg.match(
      /<path class="official-reference native-dwg"[^>]* d="([^"]+)"/,
    );
    expect(officialPath).not.toBeNull();
    expect((officialPath?.[1].match(/\bM /g) ?? []).length).toBe(
      view.official_reference_path_count,
    );
    expect(svg.indexOf('class="official-reference-mask"')).toBeLessThan(svg.indexOf('class="official-reference native-dwg"'));
  }
  const plan = manifest.views[0].mechanical_cross_check;
  expect(Math.max(...plan.absolute_delta_mm)).toBeLessThanOrEqual(0.2);
  const front = manifest.views[1].mechanical_cross_check;
  expect(Math.max(...front.absolute_delta_mm)).toBeLessThanOrEqual(0.2);
  const side = manifest.views[2].mechanical_cross_check;
  expect(side.height_absolute_delta_mm).toBeLessThanOrEqual(0.2);
  expect(side.visible_width_difference_mm).toBeGreaterThan(20);
  expect(side.note).toContain("visible 53.4 mm channel component");
});

test("project drawings retain walls and surrounding elements with top-layer official overlays", () => {
  const context = JSON.parse(readFileSync(join(product, "project-context-manifest.json"), "utf8"));
  expect(context.source_kind).toBe("native_dwg");
  expect(context.article_number).toBe("154.446.KS.1");
  expect(context.replacement_cad_used).toBeFalse();
  expect(context.project_context_retained).toBeTrue();
  expect(context.walls_and_surrounding_project_elements_retained).toBeTrue();
  expect(context.blue_line_top_layer_with_white_mask).toBeTrue();
  expect(context.official_dwg_uniform_project_scale_preserved).toBeTrue();
  expect(context.project_scale_svg_units_per_mm).toBe(0.02);
  expect(context.pass).toBeTrue();
  expect(context.views.flatMap((view: any) => view.overlays)).toHaveLength(3);
  for (const view of context.views) {
    const svg = readFileSync(join(root, view.output), "utf8");
    expect(svg).toContain('class="official-native-dwg project-context-overlay"');
    expect(svg).toContain('data-ifc-guid="2S2c498tb7$gzdukjhCGVQ"');
    expect(svg).toContain('data-source-kind="native_dwg"');
    expect(sha256(join(root, view.output))).toBe(view.output_sha256);
    expect(view.overlays[0].fit.uniform_scale_preserved).toBeTrue();
    expect(view.overlays[0].fit.transformation_mode).toBe("axis_swap_rigid_reflection_and_translation_only");
    expect(sha256(join(root, view.review_preview))).toBe(view.review_preview_sha256);
  }
  const plan = readFileSync(join(product, "project-context-sanitary-plan.svg"), "utf8");
  expect(plan).toContain("IfcWall");
  expect(plan).toContain("IfcSanitaryTerminal");
  expect(plan).toContain("IfcCovering");
  expect(plan).toContain("IfcRailing");
});

test("Bonsai evidence is an actual local-axis camera render of the IFC Body", () => {
  const evidence = JSON.parse(readFileSync(join(product, "bonsai-review-manifest.json"), "utf8"));
  expect(evidence.mode).toBe("actual_bonsai_ifc_body_camera_render");
  expect(evidence.geometry_product_count).toBe(1);
  expect(evidence.whole_model_render).toBeFalse();
  expect(evidence.formal_ifc_bytes_unchanged).toBeTrue();
  expect(evidence.bonsai_session.saved_active_representation).toBe("Body");
  expect(evidence.bonsai_session.ifc_context_identifier).toBe("Body");
  expect(evidence.bonsai_session.saved_camera_count).toBe(4);
  expect(evidence.bonsai_session.local_axis_extents_m.x).toBeCloseTo(0.9, 5);
  expect(evidence.renders.map((render: any) => render.view)).toEqual(["plan", "front", "side", "iso"]);
  for (const render of evidence.renders) {
    expect(render.camera_type).toBe("ORTHO");
    expect(render.camera_axis_basis).toBe("IFC product local axes");
    expect(render.render_operation).toBe("bpy.ops.render.render(write_still=True)");
    expect(sha256(join(root, render.path))).toBe(render.sha256);
  }
  expect(sha256(join(root, evidence.bonsai_session.path))).toBe(evidence.bonsai_session.sha256);
  expect(sha256(join(root, evidence.isolated_ifc))).toBe(evidence.isolated_ifc_sha256);
});

test("pending CleanLine50 candidate cannot write any derived drawing IFC", () => {
  const candidate = JSON.parse(readFileSync(join(product, "candidate-representations.json"), "utf8"));
  const manifest = join(product, "manifest.json");
  const approval = JSON.parse(readFileSync(join(root, "pipeline/decisions/geberit-154-446-ks-1-drawing-approval.json"), "utf8"));
  expect(candidate.formal_ifc_write_allowed).toBeFalse();
  expect(candidate).toMatchObject({
    source_kind: "native_dwg",
    official_cad_used: true,
    third_party_cad_used: false,
    blue_line_present: true,
    representation_geometry_source: "official_native_dwg_paths_mm",
    simplified_proxy_comparison_included: true,
    simplified_proxy_comparison_source_kind: "geometry_derived_simplified_proxy",
  });
  expect(candidate.review_status).toBe("visual_review_pending");
  expect(approval.status).toBe("pending");
  expect(approval.derived_ifc_write_allowed).toBeFalse();
  expect(approval.formal_authoritative_ifc_write_allowed).toBeFalse();
  expect(approval.candidate_manifest_sha256).toBe(sha256(manifest));
  expect(existsSync(join(product, "Geberit-154-446-KS-1-derived-drawing.ifc"))).toBeFalse();
  expect(sha256(formal)).toBe(formalHash);
});

test("writer rejects both missing apply and the pending approval record", () => {
  const temporary = mkdtempSync(join(tmpdir(), "geberit-154446ks1-gate-"));
  const output = join(temporary, "forbidden.ifc");
  const script = join(root, "pipeline/scripts/geberit_154_446_ks_1_drawing_ifc.py");
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

test("a scoped temporary approval writes verified native-DWG representations and source associations", () => {
  const manifest = join(product, "manifest.json");
  const temporary = mkdtempSync(join(tmpdir(), "geberit-154446ks1-approved-"));
  const approvalPath = join(temporary, "approval.json");
  const output = join(temporary, "derived.ifc");
  const report = join(temporary, "report.json");
  writeFileSync(approvalPath, JSON.stringify({
    schema_version: 1,
    profile_key: "geberit-154-446-ks-1",
    article_number: "154.446.KS.1",
    candidate_manifest_sha256: sha256(manifest),
    status: "approved",
    reviewer: "automated gate fixture",
    review_date: "2026-08-23",
    approved_views: ["plan", "front", "side"],
    derived_ifc_write_allowed: true,
    formal_authoritative_ifc_write_allowed: false,
    scope: "official manufacturer exact archived article reference, not a project shop drawing",
    approval_evidence: "temporary automated writer verification only",
  }, null, 2));
  const run = Bun.spawnSync([
    "python3",
    join(root, "pipeline/scripts/geberit_154_446_ks_1_drawing_ifc.py"),
    "--input", formal,
    "--manifest", manifest,
    "--approval", approvalPath,
    "--output", output,
    "--report", report,
    "--apply",
  ], { cwd: root, stdout: "pipe", stderr: "pipe" });
  expect(run.exitCode).toBe(0);
  if (run.exitCode !== 0) throw new Error(run.stderr.toString());
  expect(existsSync(output)).toBeTrue();
  const result = JSON.parse(readFileSync(report, "utf8"));
  expect(result.pass).toBeTrue();
  expect(result.formal_ifc_bytes_unchanged).toBeTrue();
  expect(result.representations).toEqual({
    plan: "Geberit154446KS1Plan",
    front: "Geberit154446KS1Front",
    side: "Geberit154446KS1Side",
  });
  expect(result.representation_path_counts).toEqual({ plan: 122, front: 118, side: 176 });
  expect(result.representation_geometry_source).toBe("official_native_dwg_paths_mm");
  expect(result.proxy_geometry_included).toBeFalse();
  expect(result.source_document_associations).toEqual([
    "GEBERIT-154-446-KS-1-A-NATIVE-DWG",
    "GEBERIT-154-446-KS-1-G-NATIVE-DWG",
    "GEBERIT-154-446-KS-1-L-NATIVE-DWG",
    "GEBERIT-154-446-KS-1-OFFICIAL-PRODUCT-PAGE",
    "GEBERIT-154-446-KS-1-P-NATIVE-DWG-IDENTITY",
  ]);
  expect(sha256(formal)).toBe(formalHash);
  rmSync(temporary, { recursive: true, force: true });
}, 30_000);
