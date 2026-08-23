import { expect, test } from "bun:test";
import { createHash } from "node:crypto";
import { existsSync, mkdtempSync, readFileSync, rmSync, writeFileSync } from "node:fs";
import { tmpdir } from "node:os";
import { join, resolve } from "node:path";

const root = resolve(import.meta.dir, "../..");
const product = join(root, "output/review/highpoly-types/tab02");
const formal = join(root, "2504 GBTB Yanlord Zhuhai.ifc");
const formalHash = "7a50b87e8f48a7c2bfbaa9f1dabdfca7155fa684b2325c66aa1ae15c8ab0a25c";
const sourceKind = "geometry_derived_simplified_proxy";
const sourceLabelZh = "基于原始高模几何生成的简化图纸表达";
const scope = "exact RODA BERNARDO 367 manufacturer identity and public-catalogue dimensions; official 2D/3D downloads require reserved-area authentication and were not acquired; not a project shop drawing and not official CAD geometry";

function sha256(path: string) {
  return createHash("sha256").update(readFileSync(path)).digest("hex");
}

test("exact official Bernardo 367 identity is archived without claiming CAD geometry", () => {
  const access = JSON.parse(readFileSync(join(product, "official-source/source-access-record.json"), "utf8"));
  expect(access).toMatchObject({
    manufacturer: "RODA",
    family: "Bernardo",
    model: "BERNARDO 367",
    designer: "Rodolfo Dordoni",
    project_type_code: "TAB02",
    identity_status: "exact_official_model_confirmed_by_type_description_surface_styles_and_dimensions",
    scope,
    pass: true,
  });
  expect(access.project_ifc_identity_evidence).toMatchObject({
    ifc_type_global_id: "07$ZVT0GfEAOX0ntYYymTr",
    ifc_type_description: "Bernardo 367",
    representative_global_id: "2eX84IyLr8_e34nuOoRPLQ",
    surface_style_names: ["Bernardo_serizzo", "Bernardo Metal milk"],
  });
  expect(access.official_product_cad).toMatchObject({
    authentication_required: true,
    acquired: false,
    exact_project_configuration_match: false,
    local_cad_files: [],
    third_party_cad_used: false,
  });
  expect(access.dimension_cross_check).toMatchObject({
    official_model_width_depth_height_mm: [500, 500, 670],
    official_base_diameter_height_mm: [300, 40],
    project_ifc_body_local_xyz_mm: [500, 500, 670],
    body_minus_official_mm: [0, 0, 0],
    maximum_absolute_delta_mm: 0,
    status: "exact_model_geometry_consistent_no_scaling_or_fitting",
  });
  expect(access.drawing_geometry_source).toEqual({
    source_kind: sourceKind,
    source_label_zh: sourceLabelZh,
    source_label_en: "simplified drawing representation derived from the original high-poly geometry",
    official_cad_used: false,
    third_party_cad_used: false,
  });
  const expectedHashes: Record<string, string> = {
    "RODA-official-Bernardo-367-catalogue-extract.pdf": "e4e58ccefd17c9fb50942b992c75d5018c0d30ef2e3f309dab4590c16aff32d1",
    "roda-bernardo-official-product-page.html": "8f04c1ddaa757b07a0d4cc62b419990fb10e22c073c62b889247a303d96c2318",
  };
  for (const [name, hash] of Object.entries(expectedHashes)) {
    expect(sha256(join(product, "official-source", name))).toBe(hash);
  }
  for (const source of access.official_identity_sources) {
    expect(sha256(join(root, source.local_path))).toBe(source.sha256);
  }
});

test("one TAB02 Body produces three geometry-derived views and zero blue CAD paths", () => {
  const manifest = JSON.parse(readFileSync(join(product, "manifest.json"), "utf8"));
  const candidate = JSON.parse(readFileSync(join(product, "candidate-representations.json"), "utf8"));
  expect(manifest.profile_register).toBe("output/review/highpoly-types/tab02/profile.json");
  expect(sha256(join(root, manifest.profile_register))).toBe(manifest.profile_register_sha256);
  expect(manifest.representative_global_id).toBe("2eX84IyLr8_e34nuOoRPLQ");
  expect(manifest.registered_instance_global_ids).toEqual(["2eX84IyLr8_e34nuOoRPLQ"]);
  expect(manifest.geometry_product_count).toBe(1);
  expect(manifest.whole_model_render).toBeFalse();
  expect(manifest.formal_ifc_bytes_unchanged).toBeTrue();
  expect(manifest.ifc_type_description).toBe("Bernardo 367");
  expect(manifest.article_number).toBe("RODA BERNARDO 367 / TAB02");
  expect(manifest.bounds_mm.size).toEqual([500, 500, 670]);
  expect(candidate.source_kind).toBe(sourceKind);
  expect(candidate.source_label_zh).toBe(sourceLabelZh);
  expect(candidate.official_cad_used).toBeFalse();
  expect(candidate.third_party_cad_used).toBeFalse();
  expect(manifest.views.map((view: any) => view.silhouette_path_count)).toEqual([1, 1, 1]);
  for (const view of manifest.views) {
    expect(view.drawing_line_source_kind).toBe(sourceKind);
    expect(view.official_cad_path_count).toBe(0);
    expect(view.blue_line_present).toBeFalse();
    const svg = readFileSync(join(root, view.svg), "utf8");
    expect(svg).toContain('class="simplified-proxy-silhouette geometry-derived"');
    expect(svg).toContain('data-source-kind="geometry_derived_simplified_proxy"');
    expect(svg).not.toContain('class="official-reference native-dwg"');
    expect(svg).not.toContain("#1677c8");
  }
});

test("TAB02 project plan keeps walls and furniture at exact uniform project scale", () => {
  const context = JSON.parse(readFileSync(join(product, "project-context-manifest.json"), "utf8"));
  expect(context).toMatchObject({
    source_kind: sourceKind,
    source_label_zh: sourceLabelZh,
    official_cad_used: false,
    third_party_cad_used: false,
    project_context_retained: true,
    walls_and_surrounding_project_elements_retained: true,
    overlay_top_layer_with_white_mask: true,
    blue_line_present: false,
    project_scale_svg_units_per_mm: 0.02,
    context_view_scope: {
      included: ["plan"],
      excluded: ["front", "side"],
    },
    pass: true,
  });
  expect(context.views.map((view: any) => [view.view, view.candidate_view])).toEqual([["plan", "plan"]]);
  const plan = context.views[0];
  expect(plan.overlay).toMatchObject({
    ifc_guid: "2eX84IyLr8_e34nuOoRPLQ",
    path_count: 1,
    fit: {
      rotate_quarter_turns: 0,
      alignment_mode: "centre",
      scale_svg_units_per_mm: 0.02,
      bbox_absolute_delta_svg_units: [0, 0],
      uniform_scale_preserved: true,
      pass: true,
    },
  });
  expect(sha256(join(root, plan.review_preview))).toBe(plan.review_preview_sha256);
  const svg = readFileSync(join(root, plan.output), "utf8");
  expect(svg).toContain("IfcSpace");
  expect(svg).toContain("IfcFurniture");
  expect(svg).toContain('class="geometry-derived-simplified-proxy project-context-overlay"');
  expect(svg).toContain('class="geometry-derived-proxy-mask"');
  expect(svg).toContain('class="geometry-derived-proxy"');
  expect(svg).not.toContain('class="official-reference native-dwg"');
});

test("TAB02 Bonsai evidence contains four actual IFC Body camera renders", () => {
  const evidence = JSON.parse(readFileSync(join(product, "bonsai-review-manifest.json"), "utf8"));
  expect(evidence.mode).toBe("actual_bonsai_ifc_body_camera_render");
  expect(evidence.representative_global_id).toBe("2eX84IyLr8_e34nuOoRPLQ");
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
    expect(render.camera_ortho_scale_m + 2e-7).toBeGreaterThanOrEqual(
      render.projected_width_m * render.framing_margin_factor,
    );
    expect(render.camera_ortho_scale_m + 2e-7).toBeGreaterThanOrEqual(
      render.projected_height_m * render.image_aspect * render.framing_margin_factor,
    );
    expect(sha256(join(root, render.path))).toBe(render.sha256);
  }
  expect(sha256(join(root, evidence.bonsai_session.path))).toBe(evidence.bonsai_session.sha256);
  expect(sha256(join(root, evidence.isolated_ifc))).toBe(evidence.isolated_ifc_sha256);
});

test("pending TAB02 approval cannot write any derived or formal IFC", () => {
  const candidate = JSON.parse(readFileSync(join(product, "candidate-representations.json"), "utf8"));
  const approvalPath = join(root, "pipeline/decisions/tab02-drawing-approval.json");
  const approval = JSON.parse(readFileSync(approvalPath, "utf8"));
  const manifest = join(product, "manifest.json");
  const temporary = mkdtempSync(join(tmpdir(), "tab02-gate-"));
  const output = join(temporary, "forbidden.ifc");
  const script = join(root, "pipeline/scripts/tab02_drawing_ifc.py");
  expect(candidate.formal_ifc_write_allowed).toBeFalse();
  expect(candidate.review_status).toBe("visual_review_pending");
  expect(approval.status).toBe("pending");
  expect(approval.derived_ifc_write_allowed).toBeFalse();
  expect(approval.formal_authoritative_ifc_write_allowed).toBeFalse();
  expect(approval.candidate_manifest_sha256).toBe(sha256(manifest));
  expect(existsSync(join(product, "RODA-Bernardo-367-TAB02-derived-drawing.ifc"))).toBeFalse();

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

test("scoped temporary approval writes verifiable TAB02 representations and source documents", () => {
  const manifest = join(product, "manifest.json");
  const temporary = mkdtempSync(join(tmpdir(), "tab02-approved-"));
  const approvalPath = join(temporary, "approval.json");
  const output = join(temporary, "derived.ifc");
  const report = join(temporary, "report.json");
  writeFileSync(approvalPath, JSON.stringify({
    schema_version: 1,
    profile_key: "tab02",
    ifc_type_name: "TAB02",
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
    join(root, "pipeline/scripts/tab02_drawing_ifc.py"),
    "--input", formal,
    "--manifest", manifest,
    "--approval", approvalPath,
    "--output", output,
    "--report", report,
    "--apply",
  ], { cwd: root, stdout: "pipe", stderr: "pipe" });
  if (run.exitCode !== 0) throw new Error(run.stderr.toString());
  expect(run.exitCode).toBe(0);
  expect(existsSync(output)).toBeTrue();
  const result = JSON.parse(readFileSync(report, "utf8"));
  expect(result.pass).toBeTrue();
  expect(result.formal_ifc_bytes_unchanged).toBeTrue();
  expect(result.representations).toEqual({
    plan: "Tab02Plan",
    front: "Tab02Front",
    side: "Tab02Side",
  });
  expect(result.representation_path_counts).toEqual({ plan: 1, front: 1, side: 1 });
  expect(result.source_kind).toBe(sourceKind);
  expect(result.source_label_zh).toBe(sourceLabelZh);
  expect(result.official_cad_geometry_included).toBeFalse();
  expect(result.source_document_associations).toEqual([
    "RODA-BERNARDO-367-TAB02-OFFICIAL-CATALOGUE",
    "RODA-BERNARDO-367-TAB02-OFFICIAL-CATALOGUE-EXTRACT",
    "RODA-BERNARDO-367-TAB02-OFFICIAL-PAGE-SNAPSHOT",
    "RODA-BERNARDO-367-TAB02-OFFICIAL-PRODUCT-PAGE",
    "RODA-BERNARDO-367-TAB02-SOURCE-ACCESS-RECORD",
  ]);
  expect(sha256(formal)).toBe(formalHash);
  rmSync(temporary, { recursive: true, force: true });
}, 30_000);
