import { expect, test } from "bun:test";
import { createHash } from "node:crypto";
import { existsSync, mkdtempSync, readFileSync, rmSync, writeFileSync } from "node:fs";
import { tmpdir } from "node:os";
import { join, resolve } from "node:path";

const root = resolve(import.meta.dir, "../..");
const product = join(root, "output/review/highpoly-types/street-h");
const formal = join(root, "2504 GBTB Yanlord Zhuhai.ifc");
const formalHash = "7a50b87e8f48a7c2bfbaa9f1dabdfca7155fa684b2325c66aa1ae15c8ab0a25c";
const sourceKind = "geometry_derived_simplified_proxy";
const sourceLabelZh = "基于原始高模几何生成的简化图纸表达";
const scope = "official complete antoniolupi Street family-top CAD evidence only; the project STREET-H type is an isolated repeated sink-holder subcomponent, not the complete catalogue top and not a project shop drawing";

function sha256(path: string) {
  return createHash("sha256").update(readFileSync(path)).digest("hex");
}

test("official Street DXF and technical PDF prove a parent top, not STREET-H component CAD", () => {
  const access = JSON.parse(readFileSync(join(product, "official-source/source-access-record.json"), "utf8"));
  const linework = JSON.parse(readFileSync(join(product, "official-native-dxf-linework.json"), "utf8"));
  const evidencePath = join(product, "official-source/official-product-page-evidence.json");
  expect(access).toMatchObject({
    manufacturer: "antoniolupi",
    family: "Street",
    project_type_code: "STREET-H",
    ifc_type_description: "Sink holder",
    scope,
    pass: true,
  });
  expect(access.official_product_cad).toMatchObject({
    archive_sha256: "36161ac81bade86b7ec0419c0d9cf0c49d52ae0ecb9dced39006b9eddcb10e62",
    native_dxf_sha256: "72b0a893545c4d9ff2a1d3fc79b7618eb10c00c651533eed36071e4130277c7e",
    acquired: true,
    exact_parent_family_cluster_located: true,
    exact_project_component_match: false,
    official_parent_paths_used_as_project_representation: false,
  });
  expect(sha256(join(product, "official-source/ANTONIOLUPI-official-Street-2D-CAD.zip"))).toBe(
    access.official_product_cad.archive_sha256,
  );
  expect(sha256(join(product, "official-source/AL_Street.dxf"))).toBe(
    access.official_product_cad.native_dxf_sha256,
  );
  expect(access.official_identity_sources[0].sha256).toBe(sha256(evidencePath));
  expect(linework.selected_parent_family_cluster).toEqual({
    label_handle: "1388",
    label: "street240 prof. 40 + street4054 prof. 40",
    material_variant: "marble",
    views_present: ["plan", "front"],
    side_view_present: false,
  });
  expect(linework.views.plan).toMatchObject({ path_count: 7, bounds_mm: { size: [1080, 400] } });
  expect(linework.views.front).toMatchObject({ path_count: 1, bounds_mm: { size: [1080, 250] } });
  expect(linework.project_component_cross_check).toMatchObject({
    project_STREET_H_body_bounds_mm: [300, 150, 100],
    exact_STREET_H_component_match: false,
    official_parent_paths_used_as_STREET_H_representation: false,
    pass: true,
  });
  expect(access.cad_pdf_cross_check).toMatchObject({
    resolved_dxf_drawing_units_to_mm: 10,
    exact_bounds_match_after_documented_unit_resolution: true,
    side_view_present_in_selected_dxf_cluster: false,
    pass: true,
  });
});

test("one STREET-H Body produces three black component views and zero blue parent paths", () => {
  const manifest = JSON.parse(readFileSync(join(product, "manifest.json"), "utf8"));
  const candidate = JSON.parse(readFileSync(join(product, "candidate-representations.json"), "utf8"));
  const register = JSON.parse(
    readFileSync(join(root, "pipeline/decisions/highpoly-drawing-profile-register.json"), "utf8"),
  );
  const inventory = JSON.parse(
    readFileSync(join(root, "pipeline/decisions/highpoly-product-inventory.json"), "utf8"),
  );
  expect(manifest.representative_global_id).toBe("2ajpw0I9n1dBypfISg3ejX");
  expect(manifest.registered_instance_global_ids).toEqual([
    "2ajpw0I9n1dBypfISg3ejX",
    "3PRBtW00j1D8r2Or5PTFdh",
  ]);
  expect(manifest.geometry_product_count).toBe(1);
  expect(manifest.whole_model_render).toBeFalse();
  expect(manifest.formal_ifc_bytes_unchanged).toBeTrue();
  expect(manifest.ifc_type_description).toBe("Sink holder");
  expect(manifest.bounds_mm.size).toEqual([300, 150, 100]);
  expect(register.profiles["street-h"].representative_global_id).toBe("2ajpw0I9n1dBypfISg3ejX");
  expect(inventory.review_ready_pending_approval_product_count).toBe(38);
  expect(inventory.products.find((item: any) => item.type_name === "STREET-H")).toMatchObject({
    folder: "output/review/highpoly-types/street-h",
    status: "review_ready_pending_approval",
    source_strategy: "antoniolupi_Street_parent_family_official_DXF_archived_and_mechanically_verified_exact_STREET_H_component_match_rejected_geometry_derived_proxy",
  });
  expect(candidate.source_kind).toBe(sourceKind);
  expect(candidate.source_label_zh).toBe(sourceLabelZh);
  expect(candidate.official_cad_used).toBeFalse();
  expect(candidate.third_party_cad_used).toBeFalse();
  expect(manifest.views.map((view: any) => view.silhouette_path_count)).toEqual([1, 1, 5]);
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

test("STREET-H project context keeps the real sanitary plan without inventing elevations", () => {
  const context = JSON.parse(readFileSync(join(product, "project-context-manifest.json"), "utf8"));
  expect(context).toMatchObject({
    source_kind: sourceKind,
    source_label_zh: sourceLabelZh,
    official_cad_acquired: true,
    official_cad_used: false,
    official_parent_family_cad_archived_as_rejected_reference: true,
    official_parent_family_paths_used_as_project_representation: false,
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
  expect(Math.max(...view.overlay.fit.bbox_absolute_delta_svg_units)).toBeLessThanOrEqual(0.5);
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

test("STREET-H Bonsai evidence uses four cameras on one actual IFC Body", () => {
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

test("pending STREET-H approval blocks every IFC write and preserves the formal IFC", () => {
  const candidate = JSON.parse(readFileSync(join(product, "candidate-representations.json"), "utf8"));
  const approval = JSON.parse(readFileSync(join(root, "pipeline/decisions/street-h-drawing-approval.json"), "utf8"));
  expect(candidate.formal_ifc_write_allowed).toBeFalse();
  expect(candidate.review_status).toBe("visual_review_pending");
  expect(approval.status).toBe("pending");
  expect(approval.derived_ifc_write_allowed).toBeFalse();
  expect(approval.formal_authoritative_ifc_write_allowed).toBeFalse();
  expect(approval.candidate_manifest_sha256).toBe(sha256(join(product, "manifest.json")));
  expect(existsSync(join(product, "antoniolupi-Street-H-derived-drawing.ifc"))).toBeFalse();
  expect(sha256(formal)).toBe(formalHash);
});

test("STREET-H writer rejects missing apply and the pending human approval record", () => {
  const temporary = mkdtempSync(join(tmpdir(), "street-h-gate-"));
  const output = join(temporary, "forbidden.ifc");
  const script = join(root, "pipeline/scripts/street_h_drawing_ifc.py");
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

test("a scoped temporary approval writes verified STREET-H representations and source documents", () => {
  const manifest = join(product, "manifest.json");
  const temporary = mkdtempSync(join(tmpdir(), "street-h-approved-"));
  const approvalPath = join(temporary, "approval.json");
  const output = join(temporary, "derived.ifc");
  const report = join(temporary, "report.json");
  writeFileSync(approvalPath, JSON.stringify({
    schema_version: 1,
    profile_key: "street-h",
    ifc_type_name: "STREET-H",
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
    join(root, "pipeline/scripts/street_h_drawing_ifc.py"),
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
  expect(result.representations).toEqual({ plan: "StreetHPlan", front: "StreetHFront", side: "StreetHSide" });
  expect(result.representation_path_counts).toEqual({ plan: 1, front: 1, side: 5 });
  expect(result.source_kind).toBe(sourceKind);
  expect(result.source_label_zh).toBe(sourceLabelZh);
  expect(result.official_cad_geometry_included).toBeFalse();
  expect(result.source_document_associations).toEqual([
    "ANTONIOLUPI-STREET-H-DXF-LINEWORK-AUDIT",
    "ANTONIOLUPI-STREET-H-OFFICIAL-CAD-ZIP",
    "ANTONIOLUPI-STREET-H-OFFICIAL-CATALOGUE-EXTRACT",
    "ANTONIOLUPI-STREET-H-OFFICIAL-NATIVE-DXF",
    "ANTONIOLUPI-STREET-H-OFFICIAL-PAGE-EVIDENCE",
    "ANTONIOLUPI-STREET-H-OFFICIAL-PRODUCT-PAGE",
    "ANTONIOLUPI-STREET-H-OFFICIAL-TECHNICAL-PDF",
    "ANTONIOLUPI-STREET-H-SOURCE-ACCESS-RECORD",
  ]);
  expect(sha256(formal)).toBe(formalHash);
  rmSync(temporary, { recursive: true, force: true });
}, 30_000);
