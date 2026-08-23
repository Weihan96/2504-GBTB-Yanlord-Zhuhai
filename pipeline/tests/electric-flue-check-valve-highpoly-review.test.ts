import { expect, test } from "bun:test";
import { createHash } from "node:crypto";
import { existsSync, mkdtempSync, readFileSync, rmSync, writeFileSync } from "node:fs";
import { tmpdir } from "node:os";
import { join, resolve } from "node:path";

const root = resolve(import.meta.dir, "../..");
const product = join(root, "output/review/highpoly-types/electric-flue-check-valve");
const formal = join(root, "2504 GBTB Yanlord Zhuhai.ifc");
const formalHash = "7a50b87e8f48a7c2bfbaa9f1dabdfca7155fa684b2325c66aa1ae15c8ab0a25c";
const sourceKind = "geometry_derived_simplified_proxy";
const sourceLabelZh = "基于原始高模几何生成的简化图纸表达";

function sha256(path: string) {
  return createHash("sha256").update(readFileSync(path)).digest("hex");
}

test("project evidence preserves unresolved identity and rejects near-product substitution", () => {
  const access = JSON.parse(readFileSync(join(product, "official-source/source-access-record.json"), "utf8"));
  expect(access.manufacturer).toBe("unresolved_from_project_IFC");
  expect(access.family).toBe("unresolved_from_project_IFC");
  expect(access.resolved_article_number).toBeNull();
  expect(access.external_research_boundary.exact_manufacturer_or_model_resolved).toBeFalse();
  expect(access.external_research_boundary.exact_official_CAD_resolved).toBeFalse();
  expect(access.external_research_boundary.rejected_near_candidates.map((item: any) => item.article_number ?? item.publication)).toEqual([
    "906271",
    "CN209977425U",
  ]);
  expect(access.official_product_cad.acquired).toBeFalse();
  expect(access.official_product_cad.local_cad_files).toEqual([]);
  expect(access.drawing_geometry_source).toMatchObject({
    source_kind: sourceKind,
    source_label_zh: sourceLabelZh,
    official_cad_used: false,
    third_party_cad_used: false,
  });
  for (const source of access.identity_evidence_sources) {
    expect(sha256(join(root, source.local_path))).toBe(source.sha256);
  }
  for (const source of access.project_projection_sources) {
    expect(sha256(join(root, source.local_path))).toBe(source.sha256);
  }
});

test("one untyped proxy produces three geometry-derived views and zero blue CAD paths", () => {
  const manifest = JSON.parse(readFileSync(join(product, "manifest.json"), "utf8"));
  const candidate = JSON.parse(readFileSync(join(product, "candidate-representations.json"), "utf8"));
  const profile = JSON.parse(readFileSync(join(product, "profile.json"), "utf8"));
  const inventory = JSON.parse(readFileSync(join(root, "pipeline/decisions/highpoly-product-inventory.json"), "utf8"));
  const inventoryEntry = inventory.products.find((entry: any) => entry.type_name === "Electric flue check valve");
  expect(inventory.review_ready_pending_approval_product_count).toBe(38);
  expect(inventoryEntry.status).toBe("review_ready_pending_approval");
  expect(inventoryEntry.folder).toBe("output/review/highpoly-types/electric-flue-check-valve");
  expect(manifest.representative_global_id).toBe("1faflkXXH6M9cnYPE9Liir");
  expect(manifest.registered_instance_global_ids).toEqual(["1faflkXXH6M9cnYPE9Liir"]);
  expect(manifest.geometry_product_count).toBe(1);
  expect(manifest.whole_model_render).toBeFalse();
  expect(manifest.mesh_vertex_count).toBe(766);
  expect(manifest.mesh_face_count).toBe(1400);
  expect(manifest.bounds_mm.size).toEqual([210, 240.000004, 247.500107]);
  expect(manifest.formal_ifc_bytes_unchanged).toBeTrue();
  expect(profile.profiles["electric-flue-check-valve"].official_reference.manufacturer).toBe("unresolved_from_project_IFC");
  expect(candidate.source_kind).toBe(sourceKind);
  expect(candidate.source_label_zh).toBe(sourceLabelZh);
  expect(candidate.official_cad_used).toBeFalse();
  expect(candidate.third_party_cad_used).toBeFalse();
  expect(manifest.views.map((view: any) => view.silhouette_path_count)).toEqual([2, 5, 10]);
  for (const view of manifest.views) {
    expect(view.official_cad_path_count).toBe(0);
    expect(view.blue_line_present).toBeFalse();
    const svg = readFileSync(join(root, view.svg), "utf8");
    expect(svg).toContain('data-source-kind="geometry_derived_simplified_proxy"');
    expect(svg).not.toContain('class="official-reference native-dwg"');
    expect(svg).not.toContain("#1677c8");
  }
});

test("two real R04 elevations retain project context and masked top-layer proxies", () => {
  const context = JSON.parse(readFileSync(join(product, "project-context-manifest.json"), "utf8"));
  expect(context.views.map((view: any) => [view.view, view.candidate_view])).toEqual([
    ["r04-px", "plan"],
    ["r04-ny", "front"],
  ]);
  expect(context.walls_and_surrounding_project_elements_retained).toBeTrue();
  expect(context.overlay_top_layer_with_white_mask).toBeTrue();
  expect(context.blue_line_present).toBeFalse();
  expect(context.pass).toBeTrue();
  for (const view of context.views) {
    expect(view.overlay.fit.uniform_scale_preserved).toBeTrue();
    expect(view.overlay.fit.scale_svg_units_per_mm).toBe(0.02);
    expect(view.overlay.fit.bbox_absolute_delta_svg_units).toEqual([0, 0]);
    expect(sha256(join(root, view.review_preview))).toBe(view.review_preview_sha256);
    const svg = readFileSync(join(root, view.output), "utf8");
    expect(svg).toContain("IfcWall");
    expect(svg).toContain('class="geometry-derived-simplified-proxy project-context-overlay"');
    expect(svg).toContain('class="geometry-derived-proxy-mask"');
    expect(svg).toContain('class="geometry-derived-proxy"');
    expect(svg).not.toContain('class="official-reference native-dwg"');
  }
});

test("actual Bonsai evidence contains four isolated IFC Body camera renders", () => {
  const evidence = JSON.parse(readFileSync(join(product, "bonsai-review-manifest.json"), "utf8"));
  expect(evidence.mode).toBe("actual_bonsai_ifc_body_camera_render");
  expect(evidence.representative_global_id).toBe("1faflkXXH6M9cnYPE9Liir");
  expect(evidence.geometry_product_count).toBe(1);
  expect(evidence.whole_model_render).toBeFalse();
  expect(evidence.formal_ifc_bytes_unchanged).toBeTrue();
  expect(evidence.bonsai_session.saved_active_representation).toBe("Body");
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

test("pending review leaves formal IFC byte-identical and creates no derived IFC", () => {
  const candidate = JSON.parse(readFileSync(join(product, "candidate-representations.json"), "utf8"));
  expect(candidate.formal_ifc_write_allowed).toBeFalse();
  expect(candidate.review_status).toBe("visual_review_pending");
  expect(existsSync(join(product, "Electric-flue-check-valve-derived-drawing.ifc"))).toBeFalse();
  expect(sha256(formal)).toBe(formalHash);
});

test("writer rejects missing apply and pending human approval", () => {
  const temporary = mkdtempSync(join(tmpdir(), "electric-flue-gate-"));
  const output = join(temporary, "forbidden.ifc");
  const script = join(root, "pipeline/scripts/electric_flue_check_valve_drawing_ifc.py");
  const approvalPath = join(root, "pipeline/decisions/electric-flue-check-valve-drawing-approval.json");
  const approval = JSON.parse(readFileSync(approvalPath, "utf8"));
  const manifest = join(product, "manifest.json");
  expect(approval.status).toBe("pending");
  expect(approval.derived_ifc_write_allowed).toBeFalse();
  expect(approval.formal_authoritative_ifc_write_allowed).toBeFalse();
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

test("scoped temporary approval writes verified proxy representations and project sources", () => {
  const temporary = mkdtempSync(join(tmpdir(), "electric-flue-approved-"));
  const manifest = join(product, "manifest.json");
  const approvalPath = join(temporary, "approval.json");
  const output = join(temporary, "derived.ifc");
  const report = join(temporary, "report.json");
  writeFileSync(approvalPath, JSON.stringify({
    schema_version: 1,
    profile_key: "electric-flue-check-valve",
    ifc_type_name: "Electric flue check valve",
    candidate_manifest_sha256: sha256(manifest),
    status: "approved",
    reviewer: "automated gate fixture",
    review_date: "2026-08-23",
    approved_views: ["plan", "front", "side"],
    derived_ifc_write_allowed: true,
    formal_authoritative_ifc_write_allowed: false,
    scope: "project IFC generic identity and observable geometry only; not manufacturer product identification, not official CAD geometry, not a system design and not a project shop drawing",
    approval_evidence: "temporary automated writer verification only",
  }, null, 2));
  const run = Bun.spawnSync([
    "python3",
    join(root, "pipeline/scripts/electric_flue_check_valve_drawing_ifc.py"),
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
  expect(result.representations).toEqual({
    plan: "ElectricFlueValvePlan",
    front: "ElectricFlueValveFront",
    side: "ElectricFlueValveSide",
  });
  expect(result.representation_path_counts).toEqual({ plan: 2, front: 5, side: 10 });
  expect(result.source_kind).toBe(sourceKind);
  expect(result.source_label_zh).toBe(sourceLabelZh);
  expect(result.official_cad_geometry_included).toBeFalse();
  expect(result.source_document_associations).toEqual([
    "ELECTRIC-FLUE-VALVE-ELEC-PROJECT-REVIEW",
    "ELECTRIC-FLUE-VALVE-M401-PROJECT-REVIEW",
    "ELECTRIC-FLUE-VALVE-R04-NY-PROJECTION",
    "ELECTRIC-FLUE-VALVE-R04-PX-PROJECTION",
    "ELECTRIC-FLUE-VALVE-RCP-PROJECT-REVIEW",
    "ELECTRIC-FLUE-VALVE-SOURCE-ACCESS-RECORD",
  ]);
  expect(sha256(formal)).toBe(formalHash);
  rmSync(temporary, { recursive: true, force: true });
}, 30_000);
