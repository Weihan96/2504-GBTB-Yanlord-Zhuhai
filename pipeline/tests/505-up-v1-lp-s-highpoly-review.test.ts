import { expect, test } from "bun:test";
import { createHash } from "node:crypto";
import { existsSync, mkdtempSync, readFileSync, rmSync, writeFileSync } from "node:fs";
import { tmpdir } from "node:os";
import { join, resolve } from "node:path";

const root = resolve(import.meta.dir, "../..");
const product = join(root, "output/review/highpoly-types/505-up-v1-lp-s");
const formal = join(root, "2504 GBTB Yanlord Zhuhai.ifc");
const formalHash = "7a50b87e8f48a7c2bfbaa9f1dabdfca7155fa684b2325c66aa1ae15c8ab0a25c";
const sourceKind = "geometry_derived_simplified_proxy";
const sourceLabelZh = "基于原始高模几何生成的简化图纸表达";
const scope = "exact manufacturer family and native family CAD; project-specific V1.LP.S composition is not mechanically matched to an official catalogue composition and is not a project shop drawing";

function sha256(path: string) {
  return createHash("sha256").update(readFileSync(path)).digest("hex");
}

test("official native 505 UP family DWGs are archived but do not masquerade as the project composition", () => {
  const access = JSON.parse(readFileSync(join(product, "official-source/source-access-record.json"), "utf8"));
  expect(access).toMatchObject({ manufacturer: "Molteni&C", family: "505 UP System", designer: "Nicola Gallizia", scope, pass: true });
  expect(access.official_product_cad).toMatchObject({ acquired: true, exact_project_configuration_match: false, official_cad_used_as_candidate_geometry: false });
  expect(access.configuration_cross_check).toMatchObject({ exact_match_count: 0, pass: true });
  expect(access.project_identity_evidence).toMatchObject({ family_identity_status: "confirmed", project_suffix_status: "unverified_project_or_asset_composition_suffix" });
  const dwgs = access.official_identity_sources.filter((source: any) => source.kind === "manufacturer_native_family_dwg");
  expect(dwgs.map((source: any) => source.sha256)).toEqual([
    "dfe4a6bcd655a3ed19343813a4aba8e139b847cafe8ea69b47623a6e1e71add3",
    "1e7988d5c7f00522736144fbecde1814032e5673af8cd13a47979e4347231646",
  ]);
  for (const source of dwgs) expect(sha256(join(root, source.local_path))).toBe(source.sha256);
  expect(access.native_dwg_conversion_audit.official_reference_svgs.map((item: any) => sha256(join(root, item.path)))).toEqual(
    access.native_dwg_conversion_audit.official_reference_svgs.map((item: any) => item.sha256),
  );
  expect(access.drawing_geometry_source).toMatchObject({ source_kind: sourceKind, source_label_zh: sourceLabelZh, official_cad_used: false, third_party_cad_used: false });
});

test("one 505 UP Body produces three geometry-derived candidates with zero blue CAD paths", () => {
  const manifest = JSON.parse(readFileSync(join(product, "manifest.json"), "utf8"));
  const candidate = JSON.parse(readFileSync(join(product, "candidate-representations.json"), "utf8"));
  expect(manifest).toMatchObject({ representative_global_id: "19MpdkWqXC7uhUNhLQgrce", geometry_product_count: 1, whole_model_render: false, formal_ifc_bytes_unchanged: true, review_status: "visual_review_pending", approved_for_drawing_ifc: false, pass: true });
  expect(manifest.project_context.manifest_sha256).toBe(sha256(join(root, manifest.project_context.manifest)));
  expect(manifest.bonsai_review.manifest_sha256).toBe(sha256(join(root, manifest.bonsai_review.manifest)));
  expect(candidate).toMatchObject({ article_number: "505 UP System / project V1.LP.S", source_kind: sourceKind, source_label_zh: sourceLabelZh, official_cad_used: false, third_party_cad_used: false, formal_ifc_write_allowed: false });
  expect(Object.fromEntries(Object.entries(candidate.views).map(([view, item]: any) => [view, item.proxy_paths_mm.length]))).toEqual({ plan: 19, front: 281, side: 93 });
  for (const view of ["plan", "front", "side"]) {
    expect(candidate.views[view].official_cad_paths_mm).toEqual([]);
    const svg = readFileSync(join(product, `${view}.svg`), "utf8");
    expect(svg).toContain(`data-source-kind="${sourceKind}"`);
    expect(svg).not.toContain('class="official-reference native-dwg"');
    expect(svg).not.toContain("#1677c8");
  }
});

test("project plan retains real walls and furniture while no project elevation is invented", async () => {
  const context = JSON.parse(readFileSync(join(product, "project-context-manifest.json"), "utf8"));
  expect(context.views.map((view: any) => view.view)).toEqual(["plan"]);
  expect(context.context_view_scope).toMatchObject({ included: ["plan"], excluded: ["front", "side"] });
  expect(context).toMatchObject({ project_context_retained: true, walls_and_surrounding_project_elements_retained: true, overlay_top_layer_with_white_mask: true, blue_line_present: false, pass: true });
  const view = context.views[0];
  expect(view.source).toBe("drawings/Furniture Plan.svg");
  expect(view.overlay.fit).toMatchObject({ rotate_quarter_turns: 0, uniform_scale_preserved: true, scale_svg_units_per_mm: 0.02, pass: true });
  expect(Math.max(...view.overlay.fit.bbox_absolute_delta_svg_units)).toBeLessThanOrEqual(view.overlay.fit.bbox_tolerance_svg_units);
  expect(sha256(join(root, view.review_preview))).toBe(view.review_preview_sha256);
  const svg = readFileSync(join(root, view.output), "utf8");
  expect(svg).toContain('class="geometry-derived-simplified-proxy project-context-overlay"');
  expect(svg.indexOf('class="geometry-derived-proxy-mask"')).toBeLessThan(svg.indexOf('class="geometry-derived-proxy"'));
  const nativeElevations = new Bun.Glob("drawings/elevations/native/*.svg");
  for await (const path of nativeElevations.scan(root)) expect(readFileSync(join(root, path), "utf8")).not.toContain("19MpdkWqXC7uhUNhLQgrce");
});

test("actual Bonsai session contains four camera renders of the isolated MODEL_VIEW Body", () => {
  const evidence = JSON.parse(readFileSync(join(product, "bonsai-review-manifest.json"), "utf8"));
  expect(evidence).toMatchObject({ mode: "actual_bonsai_ifc_body_camera_render", representative_global_id: "19MpdkWqXC7uhUNhLQgrce", geometry_product_count: 1, whole_model_render: false, formal_ifc_bytes_unchanged: true, pass: true });
  expect(evidence.bonsai_session).toMatchObject({ saved_active_representation: "Body", ifc_context_identifier: "Body", ifc_target_view: "MODEL_VIEW", saved_camera_count: 4, front_camera_local_y_sign: 1 });
  expect(evidence.renders.map((render: any) => render.view)).toEqual(["plan", "front", "side", "iso"]);
  for (const render of evidence.renders) {
    expect(render.render_operation).toBe("bpy.ops.render.render(write_still=True)");
    expect(sha256(join(root, render.path))).toBe(render.sha256);
  }
  expect(sha256(join(root, evidence.bonsai_session.path))).toBe(evidence.bonsai_session.sha256);
  expect(existsSync(join(product, "Molteni-505-UP-V1-LP-S-bonsai-review.blend1"))).toBeFalse();
});

test("pending approval rejects IFC writes and preserves the formal IFC bytes", () => {
  const temporary = mkdtempSync(join(tmpdir(), "molteni-505-pending-"));
  const output = join(temporary, "forbidden.ifc");
  const script = join(root, "pipeline/scripts/505_up_v1_lp_s_drawing_ifc.py");
  const approval = JSON.parse(readFileSync(join(root, "pipeline/decisions/505-up-v1-lp-s-drawing-approval.json"), "utf8"));
  expect(approval).toMatchObject({ status: "pending", derived_ifc_write_allowed: false, formal_authoritative_ifc_write_allowed: false });
  expect(approval.candidate_manifest_sha256).toBe(sha256(join(product, "manifest.json")));
  const missingApply = Bun.spawnSync(["python3", script, "--input", formal, "--output", output], { cwd: root, stderr: "pipe" });
  expect(missingApply.stderr.toString()).toContain("IFC write requires the explicit --apply flag");
  const pending = Bun.spawnSync(["python3", script, "--input", formal, "--output", output, "--apply"], { cwd: root, stderr: "pipe" });
  expect(pending.stderr.toString()).toContain("approval gate rejected IFC write");
  expect(existsSync(output)).toBeFalse();
  expect(sha256(formal)).toBe(formalHash);
  rmSync(temporary, { recursive: true, force: true });
});

test("scoped temporary approval writes verified representations and native-DWG source associations", () => {
  const manifest = join(product, "manifest.json");
  const temporary = mkdtempSync(join(tmpdir(), "molteni-505-approved-"));
  const approvalPath = join(temporary, "approval.json");
  const output = join(temporary, "derived.ifc");
  const report = join(temporary, "report.json");
  writeFileSync(approvalPath, JSON.stringify({ schema_version: 1, profile_key: "505-up-v1-lp-s", ifc_type_name: "505 UP V1.LP.S", candidate_manifest_sha256: sha256(manifest), status: "approved", reviewer: "automated gate fixture", review_date: "2026-08-23", approved_views: ["plan", "front", "side"], derived_ifc_write_allowed: true, formal_authoritative_ifc_write_allowed: false, scope, approval_evidence: "temporary automated writer verification only" }, null, 2));
  const run = Bun.spawnSync(["python3", join(root, "pipeline/scripts/505_up_v1_lp_s_drawing_ifc.py"), "--input", formal, "--manifest", manifest, "--approval", approvalPath, "--output", output, "--report", report, "--apply"], { cwd: root, stdout: "pipe", stderr: "pipe" });
  if (run.exitCode !== 0) throw new Error(run.stderr.toString());
  const result = JSON.parse(readFileSync(report, "utf8"));
  expect(result).toMatchObject({ pass: true, formal_ifc_bytes_unchanged: true, representation_path_counts: { plan: 19, front: 281, side: 93 }, official_cad_geometry_included: false, source_kind: sourceKind, source_label_zh: sourceLabelZh });
  expect(result.representations).toEqual({ plan: "Molteni505UpPlan", front: "Molteni505UpFront", side: "Molteni505UpSide" });
  expect(result.source_document_associations).toContain("MOLTENI-505-UP-OFFICIAL-TECHNICAL-DWG");
  expect(result.source_document_associations).toContain("MOLTENI-505-UP-OFFICIAL-INSPIRING-DWG");
  expect(sha256(formal)).toBe(formalHash);
  rmSync(temporary, { recursive: true, force: true });
}, 30_000);
