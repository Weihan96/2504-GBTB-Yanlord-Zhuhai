import { expect, test } from "bun:test";
import { createHash } from "node:crypto";
import { existsSync, mkdtempSync, readFileSync, rmSync, writeFileSync } from "node:fs";
import { tmpdir } from "node:os";
import { join, resolve } from "node:path";

const root = resolve(import.meta.dir, "../..");
const product = join(root, "output/review/highpoly-types/geberit-146-140");
const formal = join(root, "2504 GBTB Yanlord Zhuhai.ifc");
const formalHash = "7a50b87e8f48a7c2bfbaa9f1dabdfca7155fa684b2325c66aa1ae15c8ab0a25c";
const sourceHashes: Record<string, string> = {
  A: "8e7f2ec11311e77e8625c908d13586379fb8def19e6f65320ce61b15b2515172",
  G: "f6253efd2d736ef8fa0830455692f96c0f61b057a597f38ddc0544e3c9de5c0c",
  L: "e3c29b0f1778de7f3178a6dd6a2261dc3e7d46942e9c81e7a9943308204e9233",
  P: "229eea526831ff08de65125b3c10cbf3ea9380310a5fd674b74a77634cde7669",
};

function sha256(path: string) {
  return createHash("sha256").update(readFileSync(path)).digest("hex");
}

test("exact archived Geberit 146.140 native DWGs retain identity", () => {
  const linework = JSON.parse(readFileSync(join(product, "official-native-dwg-linework.json"), "utf8"));
  expect(linework.source_kind).toBe("native_dwg");
  expect(linework.article_number).toBe("146.140.11.1");
  expect(linework.ifc_type_name).toBe("Geberit 146.140");
  expect(linework.pass).toBeTrue();
  expect(linework.three_view_code_mapping).toMatchObject({
    G: "Grundriss / plan",
    A: "Ansicht / front elevation",
    L: "left side elevation",
  });
  for (const [code, expected] of Object.entries(sourceHashes)) {
    const source = linework.official_sources[code];
    expect(source.sha256).toBe(expected);
    expect(source.url).toBe(`https://cdn.data.geberit.com/cad/146.140.11.1_${code}.dwg`);
    expect(sha256(join(root, source.path))).toBe(expected);
  }
  expect(sha256(join(product, "official-source/geberit-aquaclean-sela-product-page.html"))).toBe(
    "8d2faaecf4e8b4e1aef578609cd21bf10fb44015cf8c916a318ebcad7bd0f157",
  );
  const access = JSON.parse(
    readFileSync(join(product, "official-source/source-access-record.json"), "utf8"),
  );
  const cdn = JSON.parse(
    readFileSync(join(product, "official-source/official-cdn-revalidation.json"), "utf8"),
  );
  expect(access).toMatchObject({
    manufacturer: "Geberit",
    article_number: "146.140.11.1",
    current_replacement_article_number: "146.220.11.1",
    replacement_cad_used: false,
    source_kind: "native_dwg",
    official_cad_acquired: true,
    official_cad_product_family_match: true,
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
    pass: true,
  });
  expect(sha256(join(root, access.official_cdn_revalidation.path))).toBe(
    access.official_cdn_revalidation.sha256,
  );
  expect(cdn).toMatchObject({
    manufacturer: "Geberit",
    article_number: "146.140.11.1",
    drawing_view_mapping: { plan: "G", front: "A", side: "L" },
    identity_only_code: "P",
    all_declared_urls_accessible: true,
    all_downloaded_bytes_match_local_archive: true,
    pass: true,
  });
  for (const [code, expected] of Object.entries(sourceHashes)) {
    expect(cdn.results[code]).toMatchObject({
      url: `https://cdn.data.geberit.com/cad/146.140.11.1_${code}.dwg`,
      http_status: 200,
      expected_sha256: expected,
      downloaded_sha256: expected,
      local_sha256: expected,
      downloaded_bytes_match_local_archive: true,
      pass: true,
    });
  }
  for (const [code, expected] of Object.entries(sourceHashes)) {
    expect(access.official_native_dwg[code].sha256).toBe(expected);
  }
});

test("official G/A/L views mechanically match the isolated IFC Body", () => {
  const manifest = JSON.parse(readFileSync(join(product, "manifest.json"), "utf8"));
  expect(manifest).toMatchObject({
    source_kind: "native_dwg",
    official_cad_acquired: true,
    official_cad_product_family_match: true,
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
  expect(manifest.profile_register).toBe("output/review/highpoly-types/geberit-146-140/profile.json");
  expect(sha256(join(root, manifest.profile_register))).toBe(manifest.profile_register_sha256);
  expect(manifest.project_context).toBe(
    "output/review/highpoly-types/geberit-146-140/project-context-manifest.json",
  );
  expect(manifest.bonsai_review).toBe(
    "output/review/highpoly-types/geberit-146-140/bonsai-review-manifest.json",
  );
  expect(manifest.official_reference.source_kind).toBe("native_dwg");
  expect(manifest.official_reference.article_number).toBe("146.140.11.1");
  expect(manifest.official_reference.current_catalog_status).toContain(
    "146.220.11.1 replaces legacy 146.140.11.1",
  );
  expect(sha256(join(root, manifest.official_reference.product_page_archive))).toBe(
    manifest.official_reference.product_page_archive_sha256,
  );
  expect(manifest.views.map((view: any) => view.blue_line_native_dwg_code)).toEqual(["G", "A", "L"]);
  expect(manifest.views.map((view: any) => view.official_reference_path_count)).toEqual([50, 82, 65]);
  for (const view of manifest.views) {
    expect(view.blue_line_source_kind).toBe("native_dwg");
    expect(view.blue_line_article_number).toBe("146.140.11.1");
    expect(view.dimension_cross_check.pass).toBeTrue();
    expect(Math.max(...view.dimension_cross_check.absolute_delta_mm)).toBeLessThanOrEqual(0.1);
    const svg = readFileSync(join(root, view.svg), "utf8");
    expect(svg).toContain('class="official-reference-mask"');
    expect(svg).toContain('class="official-reference native-dwg"');
    expect(svg).toContain('data-source-kind="native_dwg"');
    const officialPath = svg.match(
      /<path class="official-reference native-dwg"[^>]* d="([^"]+)"/,
    );
    expect((officialPath?.[1].match(/\bM /g) ?? []).length).toBe(
      view.official_reference_path_count,
    );
    expect(svg.indexOf('class="official-reference-mask"')).toBeLessThan(svg.indexOf('class="official-reference native-dwg"'));
  }
});

test("project context retains walls and surrounding elements with top-layer DWG overlays", () => {
  const context = JSON.parse(readFileSync(join(product, "project-context-manifest.json"), "utf8"));
  expect(context.source_kind).toBe("native_dwg");
  expect(context.article_number).toBe("146.140.11.1");
  expect(context.project_context_retained).toBeTrue();
  expect(context.walls_and_surrounding_project_elements_retained).toBeTrue();
  expect(context.blue_line_top_layer_with_white_mask).toBeTrue();
  expect(context.official_dwg_uniform_project_scale_preserved).toBeTrue();
  expect(context.project_scale_svg_units_per_mm).toBe(0.02);
  expect(context.pass).toBeTrue();
  expect(context.views.map((view: any) => view.view)).toEqual(["plan", "front", "side"]);
  expect(context.views.flatMap((view: any) => view.overlays)).toHaveLength(4);
  for (const overlay of context.views.flatMap((view: any) => view.overlays)) {
    expect(overlay.fit.scale_svg_units_per_mm).toBe(0.02);
    expect(overlay.fit.uniform_scale_preserved).toBeTrue();
    expect(Math.max(...overlay.fit.size_absolute_delta_svg_units)).toBeLessThanOrEqual(0.000303);
  }
  for (const crop of context.views[0].review_crops) {
    expect(sha256(join(root, crop.preview))).toBe(crop.preview_sha256);
  }
  for (const view of context.views.slice(1)) {
    expect(sha256(join(root, view.review_preview))).toBe(view.review_preview_sha256);
  }
  const plan = readFileSync(join(product, "project-context-sanitary-plan.svg"), "utf8");
  expect(plan).toContain("IfcWall");
  expect(plan).toContain("IfcSanitaryTerminal");
  expect(plan).toContain("IfcWindow");
  expect(plan).toContain('class="official-native-dwg project-context-overlay"');
  expect(plan).toContain('data-ifc-guid="0UtU7yPb10ku4gsbGoM_sp"');
  expect(plan).toContain('data-ifc-guid="1rhZG98PPCSxaLeMFLTYb9"');
});

test("Bonsai evidence is an actual camera render of the IFC Body", () => {
  const evidence = JSON.parse(readFileSync(join(product, "bonsai-review-manifest.json"), "utf8"));
  expect(evidence.mode).toBe("actual_bonsai_ifc_body_camera_render");
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
    expect(sha256(join(root, render.path))).toBe(render.sha256);
  }
  expect(sha256(join(root, evidence.bonsai_session.path))).toBe(evidence.bonsai_session.sha256);
  expect(sha256(join(root, evidence.isolated_ifc))).toBe(evidence.isolated_ifc_sha256);
});

test("unapproved Geberit candidate cannot write a derived drawing IFC", () => {
  const candidate = JSON.parse(readFileSync(join(product, "candidate-representations.json"), "utf8"));
  expect(candidate).toMatchObject({
    source_kind: "native_dwg",
    official_cad_used: true,
    third_party_cad_used: false,
    blue_line_present: true,
    representation_geometry_source: "official_native_dwg_paths_mm",
    simplified_proxy_comparison_included: true,
    simplified_proxy_comparison_source_kind: "geometry_derived_simplified_proxy",
  });
  expect(candidate.formal_ifc_write_allowed).toBeFalse();
  expect(candidate.review_status).toBe("visual_review_pending");
  expect(existsSync(join(product, "Geberit-146-140-derived-drawing.ifc"))).toBeFalse();
  expect(sha256(formal)).toBe(formalHash);
});

test("writer rejects both missing apply and the pending approval record", () => {
  const temporary = mkdtempSync(join(tmpdir(), "geberit-146140-gate-"));
  const output = join(temporary, "forbidden.ifc");
  const script = join(root, "pipeline/scripts/geberit_146_140_drawing_ifc.py");
  const approval = JSON.parse(
    readFileSync(join(root, "pipeline/decisions/geberit-146-140-drawing-approval.json"), "utf8"),
  );
  const manifest = join(product, "manifest.json");
  expect(approval.status).toBe("pending");
  expect(approval.derived_ifc_write_allowed).toBeFalse();
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

test("a scoped temporary approval writes verified native-DWG representations and source associations", () => {
  const manifest = join(product, "manifest.json");
  const temporary = mkdtempSync(join(tmpdir(), "geberit-146140-approved-"));
  const approvalPath = join(temporary, "approval.json");
  const output = join(temporary, "derived.ifc");
  const report = join(temporary, "report.json");
  writeFileSync(approvalPath, JSON.stringify({
    schema_version: 1,
    profile_key: "geberit-146-140",
    article_number: "146.140.11.1",
    candidate_manifest_sha256: sha256(manifest),
    status: "approved",
    reviewer: "automated gate fixture",
    review_date: "2026-08-23",
    approved_views: ["plan", "front", "side"],
    derived_ifc_write_allowed: true,
    scope: "official manufacturer product-family reference, not a project shop drawing",
    approval_evidence: "temporary automated writer verification only",
  }, null, 2));
  const run = Bun.spawnSync([
    "python3",
    join(root, "pipeline/scripts/geberit_146_140_drawing_ifc.py"),
    "--input", formal,
    "--manifest", manifest,
    "--approval", approvalPath,
    "--output", output,
    "--report", report,
    "--apply",
  ], { cwd: root, stdout: "pipe", stderr: "pipe" });
  expect(run.exitCode).toBe(0);
  expect(existsSync(output)).toBeTrue();
  const result = JSON.parse(readFileSync(report, "utf8"));
  expect(result.pass).toBeTrue();
  expect(result.formal_ifc_bytes_unchanged).toBeTrue();
  expect(result.representations).toEqual({
    plan: "Geberit146140Plan",
    front: "Geberit146140Front",
    side: "Geberit146140Side",
  });
  expect(result.representation_path_counts).toEqual({ plan: 50, front: 82, side: 65 });
  expect(result.representation_geometry_source).toBe("official_native_dwg_paths_mm");
  expect(result.proxy_geometry_included).toBeFalse();
  expect(result.source_document_associations).toEqual([
    "GEBERIT-146-140-A-NATIVE-DWG",
    "GEBERIT-146-140-G-NATIVE-DWG",
    "GEBERIT-146-140-L-NATIVE-DWG",
    "GEBERIT-146-140-OFFICIAL-PRODUCT-PAGE",
    "GEBERIT-146-140-P-NATIVE-DWG-IDENTITY",
  ]);
  expect(sha256(formal)).toBe(formalHash);
  rmSync(temporary, { recursive: true, force: true });
}, 30_000);
