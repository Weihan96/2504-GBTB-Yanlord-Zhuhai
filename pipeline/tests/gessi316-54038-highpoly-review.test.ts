import { expect, test } from "bun:test";
import { createHash } from "node:crypto";
import { existsSync, mkdtempSync, readFileSync, rmSync, writeFileSync } from "node:fs";
import { tmpdir } from "node:os";
import { dirname, join, resolve } from "node:path";

const root = resolve(import.meta.dir, "../..");
const product = join(root, "output/review/highpoly-types/gessi316-54038");
const sourceDir = join(product, "official-source");
const formal = join(root, "2504 GBTB Yanlord Zhuhai.ifc");
const formalHash = "7a50b87e8f48a7c2bfbaa9f1dabdfca7155fa684b2325c66aa1ae15c8ab0a25c";
const sourceKind = "native_dwg";
const sourceLabelZh = "Gessi 官方精确型号 54038 原生 DWG 图纸表达";
const dwgHash = "2c56b532bbbb78e5668050d1a7b52a7d545153a85a05efc377c7a651890897f1";
const zipHash = "e6f112fddb2bf5de9143bfddca45c2736550fcda49ea5173fd8167d73e17e39c";
const pdfHash = "1afcd3d7bf26a249904507a350665e304d2787d3ed74b1c8aa4ecf99360c209c";
const pathCounts = { plan: 119, front: 328, side: 545 };
const scope = "official Gessi exact 54038 family reference and 54139_54038 article combination; not a project shop drawing";

function sha256(path: string) {
  return createHash("sha256").update(readFileSync(path)).digest("hex");
}

test("exact public Gessi 54038 DWG and PDF sources are archived and cross-verified", () => {
  const access = JSON.parse(readFileSync(join(product, "official-source/source-access-record.json"), "utf8"));
  const revalidation = JSON.parse(readFileSync(join(sourceDir, "official-source-revalidation.json"), "utf8"));
  const pdfVerification = JSON.parse(readFileSync(join(sourceDir, "official-pdf-verification.json"), "utf8"));
  expect(access.manufacturer).toBe("Gessi");
  expect(access.family).toBe("Gessi316");
  expect(access.project_ifc_type_name).toBe("Gessi316 54038");
  expect(access.resolved_article_number).toBe("54038");
  expect(access.scope).toBe(scope);
  expect(access.official_product_cad.authentication_required).toBeFalse();
  expect(access.official_product_cad.acquired).toBeTrue();
  expect(access.official_product_cad.exact_project_configuration_match).toBeTrue();
  expect(access.drawing_geometry_source).toMatchObject({
    source_kind: sourceKind,
    source_label_zh: sourceLabelZh,
    source_dwg_sha256: dwgHash,
    official_cad_used: true,
    third_party_cad_used: false,
  });
  expect(access.identity_and_geometry_policy).toMatchObject({ "54038_native_dwg_used": true, adjacent_gessi_product_cad_used: false, third_party_cad_used: false });
  expect(sha256(join(sourceDir, "GPF5403800000G000_3.dwg"))).toBe(dwgHash);
  expect(sha256(join(sourceDir, "GPF5403800000G000_arc.zip"))).toBe(zipHash);
  expect(sha256(join(sourceDir, "GPF5403800000G000_1.pdf"))).toBe(pdfHash);
  expect(revalidation.public_access.exact_54038_native_dwg_publicly_downloadable).toBeTrue();
  expect(revalidation.zip_member.member_sha256).toBe(dwgHash);
  expect(revalidation.drawing_geometry_policy.adjacent_gessi_product_cad_used).toBeFalse();
  expect(pdfVerification.pass).toBeTrue();
  expect(pdfVerification.mechanical_cross_check.pass).toBeTrue();
  expect(pdfVerification.visual_review.status).toBe("codex_visual_qa_pass");
  for (const source of access.official_identity_sources) {
    expect(sha256(join(root, source.local_path))).toBe(source.sha256);
  }
});

test("native DWG linework extraction is repeatable and dimensionally exact", () => {
  const temporary = mkdtempSync(join(tmpdir(), "gessi316-54038-linework-"));
  const output = join(temporary, "linework.json");
  const run = Bun.spawnSync(["python3", join(root, "pipeline/scripts/gessi316_54038_linework.py"), "--output", output], { cwd: root, stdout: "pipe", stderr: "pipe" });
  if (run.exitCode !== 0) throw new Error(run.stderr.toString());
  const linework = JSON.parse(readFileSync(output, "utf8"));
  expect(linework.source_kind).toBe(sourceKind);
  expect(linework.official_sources.native_dwg.sha256).toBe(dwgHash);
  expect(linework.official_sources.native_dwg_zip.sha256).toBe(zipHash);
  expect(linework.official_sources.technical_vector_pdf.sha256).toBe(pdfHash);
  expect(linework.identity_gates).toMatchObject({ "54038_native_dwg_used": true, adjacent_gessi_product_cad_used: false, third_party_cad_used: false });
  expect(linework.nominal_dimension_cross_check).toMatchObject({ native_dwg_plan_envelope_mm: [265, 101], absolute_delta_mm: [0, 0, 0.11], tolerance_mm: 0.5, pass: true });
  for (const [view, count] of Object.entries(pathCounts)) expect(linework.views[view].paths_mm).toHaveLength(count);
  rmSync(temporary, { recursive: true, force: true });
}, 60_000);

test("one representative covers both Gessi instances with three exact native-DWG blue views", () => {
  const manifest = JSON.parse(readFileSync(join(product, "manifest.json"), "utf8"));
  const candidate = JSON.parse(readFileSync(join(product, "candidate-representations.json"), "utf8"));
  const profile = JSON.parse(readFileSync(join(product, "profile.json"), "utf8"));
  const inventory = JSON.parse(readFileSync(join(root, "pipeline/decisions/highpoly-product-inventory.json"), "utf8"));
  const inventoryEntry = inventory.products.find((entry: any) => entry.type_name === "Gessi316 54038");
  expect(inventory.review_ready_pending_approval_product_count).toBe(38);
  expect(inventoryEntry.folder).toBe("output/review/highpoly-types/gessi316-54038");
  expect(inventoryEntry.status).toBe("review_ready_pending_approval");
  expect(manifest.representative_global_id).toBe("245NU$zZL0d9tYTVwwBdk$");
  expect(manifest.registered_instance_global_ids).toEqual([
    "1sgpvCENf4xQiaypKWW2JC",
    "245NU$zZL0d9tYTVwwBdk$",
  ]);
  expect(manifest.geometry_product_count).toBe(1);
  expect(manifest.whole_model_render).toBeFalse();
  expect(manifest.formal_ifc_bytes_unchanged).toBeTrue();
  expect(manifest.article_number).toBe("54139_54038");
  expect(manifest.profile_register).toBe("output/review/highpoly-types/gessi316-54038/profile.json");
  expect(profile.profiles["gessi316-54038"].official_reference.article_number).toBe("54038");
  expect(manifest.project_context.manifest_sha256).toBe(
    sha256(join(root, manifest.project_context.manifest)),
  );
  expect(manifest.bonsai_review.manifest_sha256).toBe(
    sha256(join(root, manifest.bonsai_review.manifest)),
  );
  expect(candidate.source_kind).toBe(sourceKind);
  expect(candidate.source_label_zh).toBe(sourceLabelZh);
  expect(candidate.source_dwg_sha256).toBe(dwgHash);
  expect(candidate.official_cad_used).toBeTrue();
  expect(candidate.third_party_cad_used).toBeFalse();
  expect(manifest.views.map((view: any) => view.silhouette_path_count)).toEqual([4, 3, 6]);
  for (const view of manifest.views) {
    expect(view.drawing_line_source_kind).toBe(sourceKind);
    expect(view.official_cad_path_count).toBe(pathCounts[view.view as keyof typeof pathCounts]);
    expect(view.blue_line_present).toBeTrue();
    expect(view.white_mask_present).toBeTrue();
    expect(view.ifc_dwg_compatibility.pass).toBeTrue();
    if (view.view !== "plan") expect(view.ifc_dwg_compatibility.full_vertical_envelope_comparison_applicable).toBeFalse();
    expect(candidate.views[view.view].official_native_dwg_paths_mm).toHaveLength(pathCounts[view.view as keyof typeof pathCounts]);
    const svg = readFileSync(join(root, view.svg), "utf8");
    expect(svg).toContain('class="simplified-proxy-silhouette geometry-derived"');
    expect(svg).toContain('class="official-reference-mask"');
    expect(svg).toContain('class="official-reference native-dwg"');
    expect(svg).toContain('data-source-kind="native_dwg"');
    expect(svg).toContain("#1677c8");
    expect(svg.indexOf('class="official-reference-mask"')).toBeLessThan(svg.indexOf('class="official-reference native-dwg"'));
  }
});

test("Gessi project drawings retain walls and sanitary context with masked top-layer native DWG", () => {
  const context = JSON.parse(readFileSync(join(product, "project-context-manifest.json"), "utf8"));
  expect(context.source_kind).toBe(sourceKind);
  expect(context.source_label_zh).toBe(sourceLabelZh);
  expect(context.official_cad_used).toBeTrue();
  expect(context.third_party_cad_used).toBeFalse();
  expect(context.project_context_retained).toBeTrue();
  expect(context.walls_and_surrounding_project_elements_retained).toBeTrue();
  expect(context.overlay_top_layer_with_white_mask).toBeTrue();
  expect(context.blue_line_present).toBeTrue();
  expect(context.pass).toBeTrue();
  expect(context.views.map((view: any) => [view.view, view.candidate_view])).toEqual([
    ["plan", "plan"],
    ["front", "side"],
    ["side", "front"],
  ]);
  expect(readFileSync(join(product, "project-context-ffl-plan.svg"), "utf8")).toContain(
    "IfcFurniture",
  );
  for (const view of context.views) {
    expect(view.overlay.fit.uniform_scale_preserved).toBeTrue();
    expect(view.overlay.fit.scale_svg_units_per_mm).toBe(0.02);
    expect(view.overlay.fit.transformation_mode).toBe(
      "axis_swap_rigid_reflection_and_translation_only",
    );
    expect(view.overlay.source_dwg_sha256).toBe(dwgHash);
    expect(view.overlay.path_count).toBe(pathCounts[view.candidate_view as keyof typeof pathCounts]);
    expect(sha256(join(root, view.review_preview))).toBe(view.review_preview_sha256);
    expect(Math.max(...view.overlay.fit.bbox_absolute_delta_svg_units)).toBeLessThanOrEqual(
      view.overlay.fit.bbox_tolerance_svg_units,
    );
    const svg = readFileSync(join(root, view.output), "utf8");
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
      expect(existsSync(resolve(dirname(join(root, view.output)), match[1]))).toBeTrue();
    }
  }
});

test("Gessi Bonsai evidence contains unclipped actual IFC Body camera renders", () => {
  const evidence = JSON.parse(readFileSync(join(product, "bonsai-review-manifest.json"), "utf8"));
  expect(evidence.mode).toBe("actual_bonsai_ifc_body_camera_render");
  expect(evidence.representative_global_id).toBe("245NU$zZL0d9tYTVwwBdk$");
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
    expect(render.camera_ortho_scale_m + 1e-8).toBeGreaterThanOrEqual(
      render.projected_width_m * render.framing_margin_factor,
    );
    expect(render.camera_ortho_scale_m + 1e-8).toBeGreaterThanOrEqual(
      render.projected_height_m * render.image_aspect * render.framing_margin_factor,
    );
    expect(sha256(join(root, render.path))).toBe(render.sha256);
  }
  expect(sha256(join(root, evidence.bonsai_session.path))).toBe(evidence.bonsai_session.sha256);
  expect(sha256(join(root, evidence.isolated_ifc))).toBe(evidence.isolated_ifc_sha256);
});

test("pending Gessi review leaves the formal IFC byte-identical and creates no derived IFC", () => {
  const candidate = JSON.parse(readFileSync(join(product, "candidate-representations.json"), "utf8"));
  expect(candidate.formal_ifc_write_allowed).toBeFalse();
  expect(candidate.review_status).toBe("visual_review_pending");
  expect(existsSync(join(product, "Gessi316-54038-derived-drawing.ifc"))).toBeFalse();
  expect(sha256(formal)).toBe(formalHash);
});

test("Gessi writer rejects missing apply and the pending human approval record", () => {
  const temporary = mkdtempSync(join(tmpdir(), "gessi316-54038-gate-"));
  const output = join(temporary, "forbidden.ifc");
  const script = join(root, "pipeline/scripts/gessi316_54038_drawing_ifc.py");
  const approval = JSON.parse(
    readFileSync(join(root, "pipeline/decisions/gessi316-54038-drawing-approval.json"), "utf8"),
  );
  const manifest = join(product, "manifest.json");
  expect(approval.status).toBe("pending");
  expect(approval.derived_ifc_write_allowed).toBeFalse();
  expect(approval.formal_authoritative_ifc_write_allowed).toBeFalse();
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

test("a scoped temporary approval writes verified open native-DWG representations and sources", () => {
  const manifest = join(product, "manifest.json");
  const temporary = mkdtempSync(join(tmpdir(), "gessi316-54038-approved-"));
  const approvalPath = join(temporary, "approval.json");
  const output = join(temporary, "derived.ifc");
  const report = join(temporary, "report.json");
  writeFileSync(approvalPath, JSON.stringify({
    schema_version: 1,
    profile_key: "gessi316-54038",
    ifc_type_name: "Gessi316 54038",
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
    join(root, "pipeline/scripts/gessi316_54038_drawing_ifc.py"),
    "--input", formal,
    "--manifest", manifest,
    "--approval", approvalPath,
    "--output", output,
    "--report", report,
    "--apply",
  ], { cwd: root, stdout: "pipe", stderr: "pipe" });
  if (run.exitCode !== 0) throw new Error(run.stderr.toString());
  expect(existsSync(output)).toBeTrue();
  const result = JSON.parse(readFileSync(report, "utf8"));
  expect(result.pass).toBeTrue();
  expect(result.formal_ifc_bytes_unchanged).toBeTrue();
  expect(result.representations).toEqual({
    plan: "Gessi54038Plan",
    front: "Gessi54038Front",
    side: "Gessi54038Side",
  });
  expect(result.representation_path_counts).toEqual(pathCounts);
  expect(result.source_kind).toBe(sourceKind);
  expect(result.source_label_zh).toBe(sourceLabelZh);
  expect(result.official_cad_geometry_included).toBeTrue();
  expect(result.representation_geometry_source).toBe("official_native_dwg_paths_mm extracted from GPF5403800000G000_3.dwg");
  expect(result.source_document_associations).toEqual([
    "GESSI316-54038-OFFICIAL-CATALOGUE",
    "GESSI316-54038-OFFICIAL-NATIVE-DWG-ZIP",
    "GESSI316-54038-OFFICIAL-PRODUCT-API",
    "GESSI316-54038-OFFICIAL-PRODUCT-PAGE",
    "GESSI316-54038-OFFICIAL-TECHNICAL-PDF",
    "GESSI316-54038-SOURCE-ACCESS-RECORD",
  ]);
  expect(sha256(formal)).toBe(formalHash);
  rmSync(temporary, { recursive: true, force: true });
}, 30_000);
