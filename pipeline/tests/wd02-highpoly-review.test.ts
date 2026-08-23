import { expect, test } from "bun:test";
import { createHash } from "node:crypto";
import { existsSync, mkdtempSync, readFileSync, rmSync, writeFileSync } from "node:fs";
import { tmpdir } from "node:os";
import { join, resolve } from "node:path";

const root = resolve(import.meta.dir, "../..");
const product = join(root, "output/review/highpoly-types/wd02");
const formal = join(root, "2504 GBTB Yanlord Zhuhai.ifc");
const formalHash = "7a50b87e8f48a7c2bfbaa9f1dabdfca7155fa684b2325c66aa1ae15c8ab0a25c";
const sourceKind = "geometry_derived_simplified_proxy";
const sourceLabelZh = "基于原始高模几何生成的简化图纸表达";

function sha256(path: string) {
  return createHash("sha256").update(readFileSync(path)).digest("hex");
}

test("official Senzafine family identity is archived without claiming exact WD02 CAD", () => {
  const access = JSON.parse(readFileSync(join(product, "official-source/source-access-record.json"), "utf8"));
  const evidencePath = join(product, "official-source/official-product-page-evidence.json");
  const evidence = JSON.parse(readFileSync(evidencePath, "utf8"));
  expect(access).toMatchObject({
    manufacturer: "Poliform",
    family: "Senzafine Wardrobe",
    project_type_code: "WD02",
    ifc_type_description: "Poliform Glass Wardrobe",
    pass: true,
  });
  expect(access.official_product_identity).toMatchObject({
    product_name: "Senzafine Wardrobe",
    type: "Wardrobe",
    modular_system: true,
    customised_compositions_supported: true,
  });
  expect(access.official_product_cad).toMatchObject({
    registration_form_and_captcha_required: true,
    public_exact_native_dwg_url_located: false,
    acquired: false,
    local_cad_files: [],
    exact_project_configuration_match: false,
  });
  expect(access.dimension_cross_check).toMatchObject({
    project_ifc_body_bounds_mm: [600, 628.406189, 2392.752686],
    exact_dimension_match_claimed: false,
    geometry_scaled_to_catalogue_example: false,
  });
  expect(access.drawing_geometry_source).toMatchObject({
    source_kind: sourceKind,
    source_label_zh: sourceLabelZh,
    official_cad_used: false,
    third_party_cad_used: false,
  });
  expect(access.official_identity_sources).toHaveLength(1);
  expect(access.official_identity_sources[0].sha256).toBe(sha256(evidencePath));
  expect(access.project_registry_evidence).toEqual([
    {
      local_path: "pipeline/decisions/furniture-product-register.csv",
      sha256: sha256(join(root, "pipeline/decisions/furniture-product-register.csv")),
      selector: "type_name=WD02",
      recorded_identity: "Poliform / Senzafine / modular_wardrobe_system / confirmed",
    },
    {
      local_path: "pipeline/decisions/source-evidence-register.csv",
      sha256: sha256(join(root, "pipeline/decisions/source-evidence-register.csv")),
      selector: "source_id=FUR-SRC-WD02",
      recorded_identity: "verified_official_source / Poliform / Senzafine",
    },
  ]);
  expect(evidence.official_product_page_observations.identity.system_character).toContain("modular wardrobe system");
  expect(evidence.drawing_geometry_conclusion.official_native_dwg_acquired).toBeFalse();
});

test("one WD02 Body produces three geometry-derived views with zero blue CAD paths", () => {
  const manifest = JSON.parse(readFileSync(join(product, "manifest.json"), "utf8"));
  const candidate = JSON.parse(readFileSync(join(product, "candidate-representations.json"), "utf8"));
  const profile = JSON.parse(readFileSync(join(product, "profile.json"), "utf8"));
  const register = JSON.parse(
    readFileSync(join(root, "pipeline/decisions/highpoly-drawing-profile-register.json"), "utf8"),
  );
  const inventory = JSON.parse(
    readFileSync(join(root, "pipeline/decisions/highpoly-product-inventory.json"), "utf8"),
  );
  expect(manifest.representative_global_id).toBe("3yuXF4PtnBHgIXDlmXFJ$7");
  expect(manifest.registered_instance_global_ids).toEqual(["3yuXF4PtnBHgIXDlmXFJ$7"]);
  expect(manifest.geometry_product_count).toBe(1);
  expect(manifest.whole_model_render).toBeFalse();
  expect(manifest.formal_ifc_bytes_unchanged).toBeTrue();
  expect(manifest.ifc_type_description).toBe("Poliform Glass Wardrobe");
  expect(manifest.bounds_mm.size).toEqual([600, 628.406189, 2392.752686]);
  expect(profile.profiles.wd02.representative_global_id).toBe("3yuXF4PtnBHgIXDlmXFJ$7");
  expect(register.profiles.wd02.representative_global_id).toBe("3yuXF4PtnBHgIXDlmXFJ$7");
  expect(inventory.review_ready_pending_approval_product_count).toBe(38);
  expect(inventory.products.find((item: any) => item.type_name === "WD02")).toMatchObject({
    folder: "output/review/highpoly-types/wd02",
    status: "review_ready_pending_approval",
    source_strategy: "poliform_senzafine_glass_wardrobe_official_family_and_project_registry_identity_verified_gated_resources_exact_project_native_cad_not_acquired_geometry_derived_proxy",
  });
  expect(candidate.source_kind).toBe(sourceKind);
  expect(candidate.source_label_zh).toBe(sourceLabelZh);
  expect(candidate.official_cad_used).toBeFalse();
  expect(candidate.third_party_cad_used).toBeFalse();
  expect(manifest.views.map((view: any) => view.silhouette_path_count)).toEqual([4, 3, 3]);
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

test("WD02 project context keeps real walls and furniture in all three views", () => {
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
    pass: true,
  });
  expect(context.context_view_scope).toMatchObject({ included: ["plan", "front", "side"], excluded: [] });
  expect(context.context_view_scope.reason).toContain("actual projection group");
  expect(context.views.map((view: any) => view.view)).toEqual(["plan", "front", "side"]);
  expect(context.semantic_view_mapping).toEqual({
    plan: {
      candidate_axes: [0, 1],
      source: "drawings/Furniture Plan.svg",
      rotate_quarter_turns: 1,
      project_direction: "+Z",
    },
    front: {
      candidate_axes: [0, 2],
      source: "drawings/elevations/native/EL-09-34-R22-PX.svg",
      rotate_quarter_turns: 0,
      project_direction: "+X",
    },
    side: {
      candidate_axes: [1, 2],
      source: "drawings/elevations/native/EL-09-35-R22-NY.svg",
      rotate_quarter_turns: 0,
      project_direction: "-Y",
    },
  });
  for (const view of context.views) {
    expect(view.overlay.source_kind).toBe(sourceKind);
    expect(view.overlay.fit.uniform_scale_preserved).toBeTrue();
    expect(view.overlay.fit.scale_svg_units_per_mm).toBe(0.02);
    expect(view.overlay.fit.transformation_mode).toBe("axis_swap_rigid_reflection_and_translation_only");
    expect(Math.max(...view.overlay.fit.bbox_absolute_delta_svg_units)).toBeLessThanOrEqual(
      view.overlay.fit.bbox_tolerance_svg_units,
    );
    expect(sha256(join(root, view.review_preview))).toBe(view.review_preview_sha256);
    const svg = readFileSync(join(root, view.output), "utf8");
    expect(svg).toContain("IfcFurniture");
    expect(svg).toContain('class="geometry-derived-simplified-proxy project-context-overlay"');
    expect(svg).toContain('class="geometry-derived-proxy-mask"');
    expect(svg.indexOf('class="geometry-derived-proxy-mask"')).toBeLessThan(
      svg.indexOf('class="geometry-derived-proxy"'),
    );
  }
  expect(context.views[1].overlay.fit.bbox_absolute_delta_svg_units[1]).toBeLessThanOrEqual(0.08);
  expect(context.views[2].overlay.fit.bbox_absolute_delta_svg_units[1]).toBeLessThanOrEqual(0.08);
});

test("WD02 Bonsai evidence uses four cameras on the actual IFC Body", () => {
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

test("pending WD02 approval blocks every IFC write and preserves the formal IFC", () => {
  const candidate = JSON.parse(readFileSync(join(product, "candidate-representations.json"), "utf8"));
  const approval = JSON.parse(readFileSync(join(root, "pipeline/decisions/wd02-drawing-approval.json"), "utf8"));
  expect(candidate.formal_ifc_write_allowed).toBeFalse();
  expect(candidate.review_status).toBe("visual_review_pending");
  expect(approval.status).toBe("pending");
  expect(approval.derived_ifc_write_allowed).toBeFalse();
  expect(approval.formal_authoritative_ifc_write_allowed).toBeFalse();
  expect(approval.candidate_manifest_sha256).toBe(sha256(join(product, "manifest.json")));
  expect(existsSync(join(product, "Poliform-Senzafine-WD02-derived-drawing.ifc"))).toBeFalse();
  expect(sha256(formal)).toBe(formalHash);
});

test("WD02 writer rejects missing apply and the pending human approval record", () => {
  const temporary = mkdtempSync(join(tmpdir(), "wd02-gate-"));
  const output = join(temporary, "forbidden.ifc");
  const script = join(root, "pipeline/scripts/wd02_drawing_ifc.py");
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

test("a scoped temporary approval writes verified WD02 representations and source documents", () => {
  const manifest = join(product, "manifest.json");
  const temporary = mkdtempSync(join(tmpdir(), "wd02-approved-"));
  const approvalPath = join(temporary, "approval.json");
  const output = join(temporary, "derived.ifc");
  const report = join(temporary, "report.json");
  writeFileSync(approvalPath, JSON.stringify({
    schema_version: 1,
    profile_key: "wd02",
    ifc_type_name: "WD02",
    candidate_manifest_sha256: sha256(manifest),
    status: "approved",
    reviewer: "automated gate fixture",
    review_date: "2026-08-23",
    approved_views: ["plan", "front", "side"],
    derived_ifc_write_allowed: true,
    formal_authoritative_ifc_write_allowed: false,
    scope: "official manufacturer and Senzafine modular wardrobe-family identity only; not an exact project configuration, not official CAD geometry and not a project shop drawing",
    approval_evidence: "temporary automated writer verification only",
  }, null, 2));
  const run = Bun.spawnSync([
    "python3",
    join(root, "pipeline/scripts/wd02_drawing_ifc.py"),
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
  expect(result.representations).toEqual({ plan: "Wd02Plan", front: "Wd02Front", side: "Wd02Side" });
  expect(result.representation_path_counts).toEqual({ plan: 4, front: 3, side: 3 });
  expect(result.source_kind).toBe(sourceKind);
  expect(result.source_label_zh).toBe(sourceLabelZh);
  expect(result.official_cad_geometry_included).toBeFalse();
  expect(result.source_document_associations).toEqual([
    "POLIFORM-SENZAFINE-WD02-OFFICIAL-PAGE-EVIDENCE",
    "POLIFORM-SENZAFINE-WD02-OFFICIAL-PRODUCT-PAGE",
    "POLIFORM-SENZAFINE-WD02-OFFICIAL-TECHNICAL-CATALOG-ROUTE",
    "POLIFORM-SENZAFINE-WD02-SOURCE-ACCESS-RECORD",
  ]);
  expect(sha256(formal)).toBe(formalHash);
  rmSync(temporary, { recursive: true, force: true });
}, 30_000);
