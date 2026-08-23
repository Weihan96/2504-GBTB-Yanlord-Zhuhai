import { expect, test } from "bun:test";
import { createHash } from "node:crypto";
import { existsSync, mkdtempSync, readFileSync, rmSync, writeFileSync } from "node:fs";
import { tmpdir } from "node:os";
import { join, resolve } from "node:path";

const root = resolve(import.meta.dir, "../..");
const product = join(root, "output/review/highpoly-types/sis04");
const formal = join(root, "2504 GBTB Yanlord Zhuhai.ifc");
const formalHash = "7a50b87e8f48a7c2bfbaa9f1dabdfca7155fa684b2325c66aa1ae15c8ab0a25c";
const sourceKind = "geometry_derived_simplified_proxy";
const sourceLabelZh = "基于原始高模几何生成的简化图纸表达";

function sha256(path: string) {
  return createHash("sha256").update(readFileSync(path)).digest("hex");
}

test("official Molteni identity and dimensions do not masquerade as Wall Unit CAD", () => {
  const access = JSON.parse(readFileSync(join(product, "official-source/source-access-record.json"), "utf8"));
  const sources = Object.fromEntries(access.official_identity_sources.map((item: any) => [item.kind, item]));
  expect(access.manufacturer).toBe("Molteni&C");
  expect(access.product).toBe("Sistema 7 Wall Unit");
  expect(access.project_ifc_type_description).toBe("Sistema 7 Wall Unit 4 Doors");
  expect(access.official_product_cad.acquired).toBeFalse();
  expect(access.official_product_cad.local_cad_files).toEqual([]);
  expect(access.official_product_cad.near_name_sistema_7_doors_dwg_excluded).toBeTrue();
  expect(access.official_product_cad.third_party_cad_used).toBeFalse();
  expect(access.dimension_cross_check).toMatchObject({
    selected_official_width_height_depth_mm: [1960, 722, 331],
    project_ifc_body_local_xyz_mm: [1961.974548, 369, 721.991028],
    width_height_maximum_absolute_delta_mm: 1.974548,
    project_body_depth_excess_mm: 38,
    pass: true,
  });
  expect(sources.manufacturer_kitchen_catalogue.printed_pages_visually_checked).toEqual([272, 273, 274, 275]);
  expect(sources.manufacturer_kitchen_catalogue.pdf_file_pages_visually_checked).toEqual([156, 157]);
  for (const source of access.official_identity_sources) {
    expect(sha256(join(root, source.local_path))).toBe(source.sha256);
  }
  expect(access.drawing_geometry_source).toMatchObject({
    source_kind: sourceKind,
    source_label_zh: sourceLabelZh,
    official_cad_used: false,
    third_party_cad_used: false,
  });
});

test("single SIS04 instance has three geometry-derived views and no blue CAD line", () => {
  const manifest = JSON.parse(readFileSync(join(product, "manifest.json"), "utf8"));
  const candidate = JSON.parse(readFileSync(join(product, "candidate-representations.json"), "utf8"));
  const profile = JSON.parse(readFileSync(join(product, "profile.json"), "utf8"));
  expect(manifest.representative_global_id).toBe("1MzM8Ms2vFo8KEm503j9w2");
  expect(manifest.registered_instance_global_ids).toEqual(["1MzM8Ms2vFo8KEm503j9w2"]);
  expect(manifest.geometry_product_count).toBe(1);
  expect(manifest.whole_model_render).toBeFalse();
  expect(manifest.formal_ifc_bytes_unchanged).toBeTrue();
  expect(manifest.profile_register).toBe("output/review/highpoly-types/sis04/profile.json");
  expect(profile.profiles.sis04.representative_global_id).toBe("1MzM8Ms2vFo8KEm503j9w2");
  expect(manifest.bounds_mm.size).toEqual([1961.974548, 369, 721.991028]);
  expect(manifest.views.map((view: any) => view.silhouette_path_count)).toEqual([4, 4, 1]);
  expect(candidate.source_kind).toBe(sourceKind);
  expect(candidate.source_label_zh).toBe(sourceLabelZh);
  expect(candidate.official_cad_used).toBeFalse();
  expect(candidate.third_party_cad_used).toBeFalse();
  expect(candidate.views.plan.projection_axes).toEqual([0, 1]);
  expect(candidate.views.front.projection_axes).toEqual([0, 2]);
  expect(candidate.views.side.projection_axes).toEqual([1, 2]);
  for (const view of manifest.views) {
    expect(view.drawing_line_source_kind).toBe(sourceKind);
    expect(view.official_cad_path_count).toBe(0);
    expect(view.blue_line_present).toBeFalse();
    const svg = readFileSync(join(root, view.svg), "utf8");
    expect(svg).toContain('class="simplified-proxy-silhouette geometry-derived"');
    expect(svg).not.toContain('class="official-reference native-dwg"');
    expect(svg).not.toContain("#1677c8");
  }
});

test("project plan and R04 elevations retain walls and furniture with masked top-layer proxies", () => {
  const context = JSON.parse(readFileSync(join(product, "project-context-manifest.json"), "utf8"));
  expect(context.source_kind).toBe(sourceKind);
  expect(context.source_label_zh).toBe(sourceLabelZh);
  expect(context.official_cad_used).toBeFalse();
  expect(context.walls_and_surrounding_project_elements_retained).toBeTrue();
  expect(context.overlay_top_layer_with_white_mask).toBeTrue();
  expect(context.blue_line_present).toBeFalse();
  expect(context.pass).toBeTrue();
  expect(context.semantic_view_mapping).toMatchObject({
    plan: { candidate_axes: [0, 1], rotate_quarter_turns: 1 },
    front: { candidate_axes: [0, 2], project_direction: "+X" },
    side: { candidate_axes: [1, 2], project_direction: "+Y" },
  });
  expect(context.review_annotation_suppression).toMatchObject({
    plan: "official-elevation-anchor groups only",
    geometry_removed: false,
  });
  expect(context.review_preview_background).toMatchObject({
    plan: "opaque white review-only paper background",
    full_project_svg_remains_transparent: true,
    drawing_geometry_changed: false,
  });
  for (const view of context.views) {
    expect(view.overlay.fit.uniform_scale_preserved).toBeTrue();
    expect(view.overlay.fit.scale_svg_units_per_mm).toBe(0.02);
    expect(view.overlay.fit.transformation_mode).toBe("axis_swap_rigid_reflection_and_translation_only");
    expect(Math.max(...view.overlay.fit.bbox_absolute_delta_svg_units)).toBeLessThanOrEqual(
      view.overlay.fit.bbox_tolerance_svg_units,
    );
    expect(sha256(join(root, view.review_preview))).toBe(view.review_preview_sha256);
    const svg = readFileSync(join(root, view.output), "utf8");
    expect(svg).toContain("IfcFurniture");
    if (view.view !== "plan") expect(svg).toContain("IfcWall");
    expect(svg).toContain('class="geometry-derived-simplified-proxy project-context-overlay"');
    expect(svg).toContain('class="geometry-derived-proxy-mask"');
    expect(svg).toContain('class="geometry-derived-proxy"');
    expect(svg.indexOf('class="geometry-derived-proxy-mask"')).toBeLessThan(
      svg.indexOf('class="geometry-derived-proxy"'),
    );
  }
  const planSvg = readFileSync(join(product, "project-context-furniture-plan.svg"), "utf8");
  expect(planSvg).not.toContain('class="official-elevation-anchor"');
  const planReviewSvg = readFileSync(join(product, "project-context-furniture-plan-review.svg"), "utf8");
  expect(planReviewSvg).toContain('<svg style="background:#ffffff"');
  expect(planReviewSvg).toContain('class="review-paper-background"');
});

test("actual Bonsai file contains four camera renders from the isolated IFC Body", () => {
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
  expect(existsSync(join(product, "sis04-bonsai-review.blend1"))).toBeFalse();
});

test("pending SIS04 review cannot write an IFC and leaves the formal IFC byte-identical", () => {
  const candidate = JSON.parse(readFileSync(join(product, "candidate-representations.json"), "utf8"));
  const approval = JSON.parse(readFileSync(join(root, "pipeline/decisions/sis04-drawing-approval.json"), "utf8"));
  expect(candidate.formal_ifc_write_allowed).toBeFalse();
  expect(candidate.review_status).toBe("visual_review_pending");
  expect(approval.status).toBe("pending");
  expect(approval.derived_ifc_write_allowed).toBeFalse();
  expect(approval.formal_authoritative_ifc_write_allowed).toBeFalse();
  expect(approval.candidate_manifest_sha256).toBe(sha256(join(product, "manifest.json")));
  expect(existsSync(join(product, "Molteni-Sistema7-SIS04-derived-drawing.ifc"))).toBeFalse();
  expect(sha256(formal)).toBe(formalHash);
});

test("writer rejects missing apply and the pending SIS04 approval record", () => {
  const temporary = mkdtempSync(join(tmpdir(), "sis04-gate-"));
  const output = join(temporary, "forbidden.ifc");
  const script = join(root, "pipeline/scripts/sis04_drawing_ifc.py");
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

test("a scoped temporary approval writes verified SIS04 representations and source associations", () => {
  const manifest = join(product, "manifest.json");
  const temporary = mkdtempSync(join(tmpdir(), "sis04-approved-"));
  const approvalPath = join(temporary, "approval.json");
  const output = join(temporary, "derived.ifc");
  const report = join(temporary, "report.json");
  writeFileSync(approvalPath, JSON.stringify({
    schema_version: 1,
    profile_key: "sis04",
    ifc_type_name: "SIS04",
    candidate_manifest_sha256: sha256(manifest),
    status: "approved",
    reviewer: "automated gate fixture",
    review_date: "2026-08-23",
    approved_views: ["plan", "front", "side"],
    derived_ifc_write_allowed: true,
    formal_authoritative_ifc_write_allowed: false,
    scope: "exact Molteni&C Sistema 7 Wall Unit family identity and standard dimensions; project IFC configuration is labelled 4 Doors; exact official CAD not acquired; not a project shop drawing",
    approval_evidence: "temporary automated writer verification only",
  }, null, 2));
  const run = Bun.spawnSync([
    "python3",
    join(root, "pipeline/scripts/sis04_drawing_ifc.py"),
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
  expect(result.representations).toEqual({ plan: "Sis04Plan", front: "Sis04Front", side: "Sis04Side" });
  expect(result.representation_path_counts).toEqual({ plan: 4, front: 4, side: 1 });
  expect(result.source_kind).toBe(sourceKind);
  expect(result.source_label_zh).toBe(sourceLabelZh);
  expect(result.official_cad_geometry_included).toBeFalse();
  expect(result.source_document_associations).toEqual([
    "MOLTENI-SISTEMA7-WALL-UNIT-EXCLUDED-NEAR-NAME-DOWNLOAD",
    "MOLTENI-SISTEMA7-WALL-UNIT-OFFICIAL-KITCHEN-CATALOGUE",
    "MOLTENI-SISTEMA7-WALL-UNIT-OFFICIAL-PRODUCT-PAGE",
    "MOLTENI-SISTEMA7-WALL-UNIT-SOURCE-ACCESS-RECORD",
  ]);
  expect(sha256(formal)).toBe(formalHash);
  rmSync(temporary, { recursive: true, force: true });
}, 30_000);
