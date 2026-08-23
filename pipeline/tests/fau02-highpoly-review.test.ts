import { expect, test } from "bun:test";
import { createHash } from "node:crypto";
import { existsSync, mkdtempSync, readFileSync, rmSync, writeFileSync } from "node:fs";
import { tmpdir } from "node:os";
import { join, resolve } from "node:path";

const root = resolve(import.meta.dir, "../..");
const product = join(root, "output/review/highpoly-types/fau02");
const formal = join(root, "2504 GBTB Yanlord Zhuhai.ifc");
const formalHash = "7a50b87e8f48a7c2bfbaa9f1dabdfca7155fa684b2325c66aa1ae15c8ab0a25c";
const sourceKind = "geometry_derived_simplified_proxy";
const sourceLabelZh = "基于原始高模几何生成的简化图纸表达";
const scope = "official Falper Cilindro GH2 nearest family candidate only; project FAU02 exact article is not mechanically established; GH2 CAD is not used as representation and is not a project shop drawing";

function sha256(path: string) {
  return createHash("sha256").update(readFileSync(path)).digest("hex");
}

test("official Falper GH2 files are archived as a rejected nearest-family match", () => {
  const source = join(product, "official-source");
  const access = JSON.parse(readFileSync(join(source, "source-access-record.json"), "utf8"));
  const analysis = JSON.parse(readFileSync(join(source, "cad-match-analysis.json"), "utf8"));
  expect(access.manufacturer).toBe("Falper");
  expect(access.project_ifc_type_name).toBe("FAU02");
  expect(access.nearest_official_family_candidate).toMatchObject({
    family: "Cilindro",
    article: "GH2",
    exact_project_configuration_match: false,
  });
  expect(access.official_product_cad).toMatchObject({
    acquired: true,
    nearest_candidate_article: "GH2",
    exact_project_configuration_match: false,
    used_as_representation: false,
    third_party_cad_used: false,
  });
  expect(access.drawing_geometry_source).toMatchObject({
    source_kind: sourceKind,
    source_label_zh: sourceLabelZh,
    official_cad_used: false,
    third_party_cad_used: false,
  });
  const expectedHashes: Record<string, string> = {
    "Falper-Cilindro-GH2-2D.dwg": "e3824d70a7fd0fd16decfc4960c7b5210098de01031f385776878955062adf98",
    "Falper-Cilindro-GH2-3D.dwg": "6e98520e1c0f503e5db2b5e28e758d6a983e390ef6bcfa3525b53dcd15cc2a13",
    "Falper-Cilindro-GH2-technical-sheet.pdf": "ff4a0bb0724621b62bf6d268744e8bc1320c683ca724684dfda0bc16d8618076",
    "Falper-Cilindro-GH2-native-dwg-reference.svg": "946a963a74e37c63b3f2b3278bd999678e887bddfaff8c5909d6697ad70615e0",
    "falper-cilindro-official-product-page.html": "b8069ff53c62d3c7224502a0c16601bdc90337ca996333827c67a344c3be02f8",
  };
  for (const [name, hash] of Object.entries(expectedHashes)) {
    expect(sha256(join(source, name))).toBe(hash);
  }
  expect(analysis.native_source_cross_checks.autocad_2d_equals_technical_dwg).toBeTrue();
  expect(analysis.exact_project_configuration_match).toBeFalse();
  expect(analysis.official_CAD_used_as_representation).toBeFalse();
  expect(Math.abs(analysis.mechanical_deltas.height_mm)).toBeGreaterThan(30);
  expect(Math.abs(analysis.mechanical_deltas.stem_diameter_mm)).toBeGreaterThan(2);
  expect(Math.abs(analysis.mechanical_deltas.floor_flange_diameter_mm)).toBeGreaterThan(16);
});

test("FAU02 has one representative, three Body-derived views and zero blue CAD paths", () => {
  const manifest = JSON.parse(readFileSync(join(product, "manifest.json"), "utf8"));
  const candidate = JSON.parse(readFileSync(join(product, "candidate-representations.json"), "utf8"));
  const inventory = JSON.parse(readFileSync(join(root, "pipeline/decisions/highpoly-product-inventory.json"), "utf8"));
  const entry = inventory.products.find((item: any) => item.type_name === "FAU02");
  expect(entry.folder).toBe("output/review/highpoly-types/fau02");
  expect(entry.status).toBe("review_ready_pending_approval");
  expect(manifest.representative_global_id).toBe("36ZX3QPyD7SvlXsDKMP8rY");
  expect(manifest.registered_instance_global_ids).toEqual(["36ZX3QPyD7SvlXsDKMP8rY"]);
  expect(manifest.geometry_product_count).toBe(1);
  expect(manifest.whole_model_render).toBeFalse();
  expect(manifest.formal_ifc_bytes_unchanged).toBeTrue();
  expect(manifest.official_cad_acquired).toBeTrue();
  expect(manifest.official_cad_exact_project_configuration_match).toBeFalse();
  expect(manifest.official_cad_used).toBeFalse();
  expect(manifest.blue_line_present).toBeFalse();
  expect(candidate.source_kind).toBe(sourceKind);
  expect(candidate.official_cad_used).toBeFalse();
  expect(candidate.third_party_cad_used).toBeFalse();
  expect(manifest.views.map((item: any) => item.silhouette_path_count)).toEqual([31, 1, 3]);
  for (const view of manifest.views) {
    expect(view.official_cad_path_count).toBe(0);
    expect(view.blue_line_present).toBeFalse();
    const svg = readFileSync(join(root, view.svg), "utf8");
    expect(svg).toContain('class="simplified-proxy-silhouette geometry-derived"');
    expect(svg).not.toContain('class="official-reference native-dwg"');
    expect(svg).not.toContain("#1677c8");
  }
});

test("FAU02 context drawings retain actual IFC occlusion, walls and furniture", () => {
  const context = JSON.parse(readFileSync(join(product, "project-context-manifest.json"), "utf8"));
  expect(context.source_kind).toBe(sourceKind);
  expect(context.source_label_zh).toBe(sourceLabelZh);
  expect(context.candidate_source_kind).toBe(sourceKind);
  expect(context.candidate_source_label_zh).toBe(sourceLabelZh);
  expect(context.project_context_retained).toBeTrue();
  expect(context.walls_and_surrounding_project_elements_retained).toBeTrue();
  expect(context.furniture_retained_in_elevations).toBeTrue();
  expect(context.existing_actual_ifc_projection_retained).toBeTrue();
  expect(context.top_layer_substitution_performed).toBeFalse();
  expect(context.occlusion_retained).toBeTrue();
  expect(context.official_cad_used).toBeFalse();
  expect(context.blue_line_present).toBeFalse();
  expect(context.views.map((item: any) => item.view)).toEqual(["plan", "front", "side"]);
  for (const view of context.views) {
    expect(view.existing_project_projection.ifc_guid).toBe("36ZX3QPyD7SvlXsDKMP8rY");
    expect(view.existing_project_projection.geometry_source).toBe("actual_project_IFC_representation");
    expect(view.existing_project_projection.occlusion_retained).toBeTrue();
    expect(sha256(join(root, view.output))).toBe(view.output_sha256);
    expect(sha256(join(root, view.review_preview))).toBe(view.review_preview_sha256);
    const svg = readFileSync(join(root, view.output), "utf8");
    expect(svg).toContain("IfcWall");
    expect(svg).toContain("IfcSanitaryTerminal");
    expect(svg).toContain("36ZX3QPyD7SvlXsDKMP8rY");
    if (view.view !== "plan") expect(svg).toContain("IfcFurniture");
    expect(svg).not.toContain('class="official-reference native-dwg"');
  }
});

test("FAU02 Bonsai evidence is rendered from one actual IFC Body with four cameras", () => {
  const evidence = JSON.parse(readFileSync(join(product, "bonsai-review-manifest.json"), "utf8"));
  expect(evidence.mode).toBe("actual_bonsai_ifc_body_camera_render");
  expect(evidence.representative_global_id).toBe("36ZX3QPyD7SvlXsDKMP8rY");
  expect(evidence.geometry_product_count).toBe(1);
  expect(evidence.whole_model_render).toBeFalse();
  expect(evidence.formal_ifc_bytes_unchanged).toBeTrue();
  expect(evidence.bonsai_session.saved_active_representation).toBe("Body");
  expect(evidence.bonsai_session.ifc_context_identifier).toBe("Body");
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

test("pending FAU02 review leaves the formal IFC unchanged and creates no derived IFC", () => {
  const approval = JSON.parse(readFileSync(join(root, "pipeline/decisions/fau02-drawing-approval.json"), "utf8"));
  const manifest = join(product, "manifest.json");
  expect(approval.status).toBe("pending");
  expect(approval.derived_ifc_write_allowed).toBeFalse();
  expect(approval.formal_authoritative_ifc_write_allowed).toBeFalse();
  expect(approval.candidate_manifest_sha256).toBe(sha256(manifest));
  expect(existsSync(join(product, "FAU02-derived-drawing.ifc"))).toBeFalse();
  expect(sha256(formal)).toBe(formalHash);
});

test("FAU02 writer rejects missing apply and the pending approval", () => {
  const temporary = mkdtempSync(join(tmpdir(), "fau02-gate-"));
  const output = join(temporary, "forbidden.ifc");
  const script = join(root, "pipeline/scripts/fau02_drawing_ifc.py");
  const withoutApply = Bun.spawnSync(["python3", script, "--input", formal, "--output", output], {
    cwd: root,
    stderr: "pipe",
  });
  expect(withoutApply.exitCode).not.toBe(0);
  expect(withoutApply.stderr.toString()).toContain("IFC write requires the explicit --apply flag");
  const pending = Bun.spawnSync(["python3", script, "--input", formal, "--output", output, "--apply"], {
    cwd: root,
    stderr: "pipe",
  });
  expect(pending.exitCode).not.toBe(0);
  expect(pending.stderr.toString()).toContain("approval gate rejected IFC write");
  expect(existsSync(output)).toBeFalse();
  expect(sha256(formal)).toBe(formalHash);
  rmSync(temporary, { recursive: true, force: true });
});

test("a scoped temporary approval writes Body-derived FAU02 views and source associations", () => {
  const manifest = join(product, "manifest.json");
  const temporary = mkdtempSync(join(tmpdir(), "fau02-approved-"));
  const approval = join(temporary, "approval.json");
  const output = join(temporary, "derived.ifc");
  const report = join(temporary, "report.json");
  writeFileSync(approval, JSON.stringify({
    schema_version: 1,
    profile_key: "fau02",
    ifc_type_name: "FAU02",
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
    "python3",
    join(root, "pipeline/scripts/fau02_drawing_ifc.py"),
    "--input", formal,
    "--manifest", manifest,
    "--approval", approval,
    "--output", output,
    "--report", report,
    "--apply",
  ], { cwd: root, stdout: "pipe", stderr: "pipe" });
  expect(run.exitCode).toBe(0);
  if (run.exitCode !== 0) throw new Error(run.stderr.toString());
  const result = JSON.parse(readFileSync(report, "utf8"));
  expect(result.pass).toBeTrue();
  expect(result.formal_ifc_bytes_unchanged).toBeTrue();
  expect(result.representations).toEqual({
    plan: "FAU02Plan",
    front: "FAU02Front",
    side: "FAU02Side",
  });
  expect(result.representation_path_counts).toEqual({ plan: 31, front: 1, side: 3 });
  expect(result.source_kind).toBe(sourceKind);
  expect(result.official_cad_geometry_included).toBeFalse();
  expect(result.source_document_associations).toContain("FALPER-FAU02-GH2-2D-DWG");
  expect(result.source_document_associations).toContain("FALPER-FAU02-CAD-MATCH-ANALYSIS");
  expect(result.source_document_associations.length).toBe(8);
  expect(sha256(formal)).toBe(formalHash);
  rmSync(temporary, { recursive: true, force: true });
}, 30_000);
