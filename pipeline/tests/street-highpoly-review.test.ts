import { expect, test } from "bun:test";
import { createHash } from "node:crypto";
import { existsSync, mkdtempSync, readFileSync, rmSync, writeFileSync } from "node:fs";
import { tmpdir } from "node:os";
import { join, resolve } from "node:path";

const root = resolve(import.meta.dir, "../..");
const product = join(root, "output/review/highpoly-types/street");
const formal = join(root, "2504 GBTB Yanlord Zhuhai.ifc");
const formalHash = "7a50b87e8f48a7c2bfbaa9f1dabdfca7155fa684b2325c66aa1ae15c8ab0a25c";
const sourceKind = "geometry_derived_simplified_proxy";
const sourceLabelZh = "基于原始高模几何生成的简化图纸表达";
const scope = "official antoniolupi Street family CAD evidence for a domestic-custom project top; no exact 1000 x 470 x 250 mm official configuration is present and no project shop drawing is claimed";

function sha256(path: string) {
  return createHash("sha256").update(readFileSync(path)).digest("hex");
}

test("official Street DXF proves that the domestic-custom STREET configuration has no exact CAD", () => {
  const access = JSON.parse(readFileSync(join(product, "official-source/source-access-record.json"), "utf8"));
  const audit = JSON.parse(readFileSync(join(product, "official-native-dxf-configuration-audit.json"), "utf8"));
  const evidencePath = join(product, "official-source/official-product-page-evidence.json");
  expect(access).toMatchObject({
    manufacturer: "antoniolupi",
    family: "Street",
    project_type_code: "STREET",
    ifc_type_description: "antoniolupi street240 prof. 40 + street4054 prof. 40",
    scope,
    pass: true,
  });
  expect(access.official_product_cad).toMatchObject({
    archive_sha256: "36161ac81bade86b7ec0419c0d9cf0c49d52ae0ecb9dced39006b9eddcb10e62",
    native_dxf_sha256: "72b0a893545c4d9ff2a1d3fc79b7618eb10c00c651533eed36071e4130277c7e",
    acquired: true,
    exact_family_located: true,
    exact_project_configuration_match: false,
    official_dxf_paths_used_as_project_representation: false,
  });
  expect(sha256(join(product, "official-source/ANTONIOLUPI-official-Street-2D-CAD.zip"))).toBe(
    access.official_product_cad.archive_sha256,
  );
  expect(sha256(join(product, "official-source/AL_Street.dxf"))).toBe(
    access.official_product_cad.native_dxf_sha256,
  );
  expect(access.official_identity_sources[0].sha256).toBe(sha256(evidencePath));
  expect(audit.ifc_description_cluster).toMatchObject({
    label: "street240 prof. 40 + street4054 prof. 40",
    configuration_bounds_mm: [1080, 400, 250],
    matches_project_body: false,
  });
  expect(audit.nearest_official_depth_cluster).toMatchObject({
    label: "street147 prof. 47 + street4754 prof. 40",
    configuration_bounds_mm: [1080, 470, 250],
    matches_project_body: false,
  });
  expect(audit).toMatchObject({
    project_body_bounds_mm: [1000, 470, 250],
    exact_1000_x_470_plan_cluster_present: false,
    exact_project_configuration_match: false,
    official_dxf_paths_used_as_project_representation: false,
    pass: true,
  });
  expect(audit.all_official_depth_470_plan_widths_mm).toContain(1080);
  expect(audit.all_official_depth_470_plan_widths_mm).not.toContain(1000);
  expect(access.configuration_cross_check).toMatchObject({
    project_ifc_body_bounds_mm: [1000, 470, 250],
    exact_1000_x_470_plan_cluster_present: false,
    geometry_scaled_to_official_cluster: false,
    pass: true,
  });
});

test("one STREET Body produces three black component views and zero blue parent paths", () => {
  const manifest = JSON.parse(readFileSync(join(product, "manifest.json"), "utf8"));
  const candidate = JSON.parse(readFileSync(join(product, "candidate-representations.json"), "utf8"));
  const register = JSON.parse(
    readFileSync(join(root, "pipeline/decisions/highpoly-drawing-profile-register.json"), "utf8"),
  );
  const inventory = JSON.parse(
    readFileSync(join(root, "pipeline/decisions/highpoly-product-inventory.json"), "utf8"),
  );
  expect(manifest.representative_global_id).toBe("1FgLPMw$5B4wBH2ySMkXE1");
  expect(manifest.registered_instance_global_ids).toEqual(["1FgLPMw$5B4wBH2ySMkXE1"]);
  expect(manifest.geometry_product_count).toBe(1);
  expect(manifest.whole_model_render).toBeFalse();
  expect(manifest.formal_ifc_bytes_unchanged).toBeTrue();
  expect(manifest.ifc_type_description).toBe("antoniolupi street240 prof. 40 + street4054 prof. 40");
  expect(manifest.bounds_mm.size).toEqual([1000, 470, 250]);
  expect(register.profiles["street"].representative_global_id).toBe("1FgLPMw$5B4wBH2ySMkXE1");
  expect(inventory.review_ready_pending_approval_product_count).toBe(38);
  expect(inventory.products.find((item: any) => item.type_name === "STREET")).toMatchObject({
    folder: "output/review/highpoly-types/street",
    status: "review_ready_pending_approval",
    source_strategy: "antoniolupi_Street_official_family_DXF_archived_domestic_custom_1000x470x250_exact_configuration_absent_geometry_derived_proxy",
  });
  expect(candidate.source_kind).toBe(sourceKind);
  expect(candidate.source_label_zh).toBe(sourceLabelZh);
  expect(candidate.official_cad_used).toBeFalse();
  expect(candidate.third_party_cad_used).toBeFalse();
  expect(manifest.views.map((view: any) => view.silhouette_path_count)).toEqual([3, 1, 2]);
  for (const view of manifest.views) {
    expect(view.official_cad_path_count).toBe(0);
    expect(view.blue_line_present).toBeFalse();
    const svg = readFileSync(join(root, view.svg), "utf8");
    expect(svg).toContain('class="simplified-proxy-silhouette geometry-derived"');
    expect(svg).toContain('data-source-kind="geometry_derived_simplified_proxy"');
    expect(svg).not.toContain('class="official-reference native-dwg"');
    expect(svg).not.toContain("#1677c8");
  }
});

test("STREET project context keeps the real sanitary plan without inventing elevations", () => {
  const context = JSON.parse(readFileSync(join(product, "project-context-manifest.json"), "utf8"));
  expect(context).toMatchObject({
    source_kind: sourceKind,
    source_label_zh: sourceLabelZh,
    official_cad_acquired: true,
    official_cad_used: false,
    official_family_cad_archived_as_rejected_reference: true,
    official_family_paths_used_as_project_representation: false,
    third_party_cad_used: false,
    project_context_retained: true,
    walls_and_surrounding_project_elements_retained: true,
    overlay_top_layer_with_white_mask: true,
    blue_line_present: false,
    pass: true,
  });
  expect(context.context_view_scope).toMatchObject({ included: ["plan"], excluded: ["front", "side"] });
  expect(context.context_view_scope.reason).toContain("project elevations are not invented");
  expect(context.views.map((view: any) => view.view)).toEqual(["plan"]);
  const view = context.views[0];
  expect(view.overlay.fit.uniform_scale_preserved).toBeTrue();
  expect(view.overlay.fit.scale_svg_units_per_mm).toBe(0.02);
  expect(view.overlay.fit.transformation_mode).toBe("axis_swap_rigid_reflection_and_translation_only");
  expect(Math.max(...view.overlay.fit.bbox_absolute_delta_svg_units)).toBeLessThanOrEqual(2.5);
  expect(sha256(join(root, view.review_preview))).toBe(view.review_preview_sha256);
  const svg = readFileSync(join(root, view.output), "utf8");
  expect(svg).toContain("IfcSanitaryTerminal");
  expect(svg).toContain('class="geometry-derived-simplified-proxy project-context-overlay"');
  expect(svg).toContain('class="geometry-derived-proxy-mask"');
  expect(svg.indexOf('class="geometry-derived-proxy-mask"')).toBeLessThan(
    svg.indexOf('class="geometry-derived-proxy"'),
  );
  expect(svg).not.toContain('class="official-reference native-dwg"');
});

test("STREET Bonsai evidence uses four cameras on one actual IFC Body", () => {
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

test("pending STREET approval blocks every IFC write and preserves the formal IFC", () => {
  const candidate = JSON.parse(readFileSync(join(product, "candidate-representations.json"), "utf8"));
  const approval = JSON.parse(readFileSync(join(root, "pipeline/decisions/street-drawing-approval.json"), "utf8"));
  expect(candidate.formal_ifc_write_allowed).toBeFalse();
  expect(candidate.review_status).toBe("visual_review_pending");
  expect(approval.status).toBe("pending");
  expect(approval.derived_ifc_write_allowed).toBeFalse();
  expect(approval.formal_authoritative_ifc_write_allowed).toBeFalse();
  expect(approval.candidate_manifest_sha256).toBe(sha256(join(product, "manifest.json")));
  expect(existsSync(join(product, "antoniolupi-Street-derived-drawing.ifc"))).toBeFalse();
  expect(sha256(formal)).toBe(formalHash);
});

test("STREET writer rejects missing apply and the pending human approval record", () => {
  const temporary = mkdtempSync(join(tmpdir(), "street-gate-"));
  const output = join(temporary, "forbidden.ifc");
  const script = join(root, "pipeline/scripts/street_drawing_ifc.py");
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

test("a scoped temporary approval writes verified STREET representations and source documents", () => {
  const manifest = join(product, "manifest.json");
  const temporary = mkdtempSync(join(tmpdir(), "street-approved-"));
  const approvalPath = join(temporary, "approval.json");
  const output = join(temporary, "derived.ifc");
  const report = join(temporary, "report.json");
  writeFileSync(approvalPath, JSON.stringify({
    schema_version: 1,
    profile_key: "street",
    ifc_type_name: "STREET",
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
    join(root, "pipeline/scripts/street_drawing_ifc.py"),
    "--input", formal,
    "--manifest", manifest,
    "--approval", approvalPath,
    "--output", output,
    "--report", report,
    "--apply",
  ], { cwd: root, stdout: "pipe", stderr: "pipe" });
  expect(run.exitCode).toBe(0);
  if (run.exitCode !== 0) throw new Error(run.stderr.toString());
  const result = JSON.parse(readFileSync(report, "utf8"));
  expect(result.pass).toBeTrue();
  expect(result.formal_ifc_bytes_unchanged).toBeTrue();
  expect(result.representations).toEqual({ plan: "StreetPlan", front: "StreetFront", side: "StreetSide" });
  expect(result.representation_path_counts).toEqual({ plan: 3, front: 1, side: 2 });
  expect(result.source_kind).toBe(sourceKind);
  expect(result.source_label_zh).toBe(sourceLabelZh);
  expect(result.official_cad_geometry_included).toBeFalse();
  expect(result.source_document_associations).toEqual([
    "ANTONIOLUPI-STREET-DXF-CONFIGURATION-AUDIT",
    "ANTONIOLUPI-STREET-OFFICIAL-CAD-ZIP",
    "ANTONIOLUPI-STREET-OFFICIAL-CATALOGUE-EXTRACT",
    "ANTONIOLUPI-STREET-OFFICIAL-NATIVE-DXF",
    "ANTONIOLUPI-STREET-OFFICIAL-PAGE-EVIDENCE",
    "ANTONIOLUPI-STREET-OFFICIAL-PRODUCT-PAGE",
    "ANTONIOLUPI-STREET-OFFICIAL-TECHNICAL-PDF",
    "ANTONIOLUPI-STREET-SOURCE-ACCESS-RECORD",
  ]);
  expect(sha256(formal)).toBe(formalHash);
  rmSync(temporary, { recursive: true, force: true });
}, 30_000);
