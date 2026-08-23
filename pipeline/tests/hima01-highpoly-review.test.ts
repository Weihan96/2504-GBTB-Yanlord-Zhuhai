import { expect, test } from "bun:test";
import { createHash } from "node:crypto";
import { existsSync, mkdtempSync, readFileSync, rmSync, writeFileSync } from "node:fs";
import { tmpdir } from "node:os";
import { join, resolve } from "node:path";

const root = resolve(import.meta.dir, "../..");
const product = join(root, "output/review/highpoly-types/hima01");
const formal = join(root, "2504 GBTB Yanlord Zhuhai.ifc");
const formalHash = "7a50b87e8f48a7c2bfbaa9f1dabdfca7155fa684b2325c66aa1ae15c8ab0a25c";
const sourceKind = "geometry_derived_simplified_proxy";
const sourceLabelZh = "基于原始高模几何生成的简化图纸表达";

function sha256(path: string) {
  return createHash("sha256").update(readFileSync(path)).digest("hex");
}

test("official Hima identity does not masquerade as acquired CAD geometry", () => {
  const access = JSON.parse(readFileSync(join(product, "official-source/source-access-record.json"), "utf8"));
  const globalProfile = JSON.parse(
    readFileSync(join(root, "pipeline/decisions/highpoly-drawing-profile-register.json"), "utf8"),
  ).profiles.hima01;
  const evidence = JSON.parse(
    readFileSync(join(product, "official-source/official-product-page-evidence.json"), "utf8"),
  );
  expect(access.manufacturer).toBe("Poliform");
  expect(access.family).toBe("Hima");
  expect(access.official_2d_dwg.published_on_product_page).toBeTrue();
  expect(access.official_2d_dwg.acquired).toBeFalse();
  expect(access.official_2d_dwg.local_path).toBeNull();
  expect(access.official_2d_dwg.sha256).toBeNull();
  expect(access.official_2d_dwg.access_status).toBe(
    "published_registration_form_and_captcha_required_not_acquired",
  );
  expect(access.official_2d_dwg.public_exact_native_dwg_url_located).toBeFalse();
  expect(access.official_2d_dwg.download_link_behavior).toBe(
    "same-page registration modal; no public native-DWG asset URL exposed",
  );
  expect(access.official_technical_sheet).toBe(
    "https://www.poliform.it/assets/pdf/201250-hima-poliform-en-us.pdf",
  );
  expect(access.official_technical_sheet_archive).toMatchObject({
    acquired: false,
    local_path: null,
    sha256: null,
    direct_https_status: 403,
  });
  expect(sha256(join(root, access.official_product_page_evidence.local_path))).toBe(
    access.official_product_page_evidence.sha256,
  );
  expect(evidence.official_product_page_observations.download_listing).toMatchObject({
    label: "2D (DWG)",
    published_size: "2 MB",
    listed_on_official_page: true,
  });
  expect(evidence.official_product_page_observations.download_access_form).toMatchObject({
    submitted: false,
    personal_data_entered: false,
    captcha_bypassed: false,
  });
  expect(evidence.official_product_page_observations.download_listing.link_behavior_observed).toBe(
    "same-page registration modal; no public native-DWG asset URL exposed",
  );
  expect(evidence.official_technical_publications[0]).toMatchObject({
    url: "https://www.poliform.it/assets/pdf/201250-hima-poliform-en-us.pdf",
    observed_page_count: 4,
    local_archive_acquired: false,
    direct_https_status: 403,
  });
  expect(access.official_dimension_cross_check).toMatchObject({
    published_nominal_height_mm: 1000,
    project_ifc_body_height_mm: 999,
    pass: true,
  });
  expect(access.drawing_geometry_source.source_kind).toBe(sourceKind);
  expect(access.drawing_geometry_source.source_label_zh).toBe(sourceLabelZh);
  expect(access.drawing_geometry_source.official_cad_used).toBeFalse();
  expect(access.drawing_geometry_source.third_party_cad_used).toBeFalse();
  expect(globalProfile.drawing_source).toEqual(access.drawing_geometry_source);
  expect(globalProfile.official_reference).toMatchObject({
    source_kind: "manufacturer_product_page_and_technical_publication_identity_only",
    technical_sheet: "https://www.poliform.it/assets/pdf/201250-hima-poliform-en-us.pdf",
    technical_sheet_legacy: "https://www.poliform.it/assets/pdf/200467-hima-poliform-en.pdf",
    official_2d_dwg_status: "published_registration_form_and_captcha_required_not_acquired",
    official_2d_dwg_used: false,
  });
});

test("single-instance review has three transparent geometry-derived views and no blue CAD line", () => {
  const manifest = JSON.parse(readFileSync(join(product, "manifest.json"), "utf8"));
  const candidate = JSON.parse(readFileSync(join(product, "candidate-representations.json"), "utf8"));
  const profile = JSON.parse(readFileSync(join(product, "profile.json"), "utf8"));
  expect(manifest.representative_global_id).toBe("2xmcLzu1rDTeMzRuNxPDyE");
  expect(manifest.registered_instance_global_ids).toEqual(["2xmcLzu1rDTeMzRuNxPDyE"]);
  expect(manifest.geometry_product_count).toBe(1);
  expect(manifest.whole_model_render).toBeFalse();
  expect(manifest.formal_ifc_bytes_unchanged).toBeTrue();
  expect(manifest.source_kind).toBe(sourceKind);
  expect(manifest.source_label_zh).toBe(sourceLabelZh);
  expect(manifest.official_cad_used).toBeFalse();
  expect(manifest.third_party_cad_used).toBeFalse();
  expect(manifest.blue_line_present).toBeFalse();
  expect(manifest.profile_register).toBe("output/review/highpoly-types/hima01/profile.json");
  expect(profile.profiles.hima01.representative_global_id).toBe("2xmcLzu1rDTeMzRuNxPDyE");
  expect(manifest.project_context.manifest_sha256).toBe(
    sha256(join(root, manifest.project_context.manifest)),
  );
  expect(manifest.bonsai_review.manifest_sha256).toBe(
    sha256(join(root, manifest.bonsai_review.manifest)),
  );
  expect(manifest.identity_checks.screen_element_count).toMatchObject({ expected: 3, actual: 3, pass: true });
  expect(manifest.identity_checks.nominal_height.pass).toBeTrue();
  expect(candidate.source_kind).toBe(sourceKind);
  expect(candidate.source_label_zh).toBe(sourceLabelZh);
  expect(candidate.official_cad_used).toBeFalse();
  expect(candidate.third_party_cad_used).toBeFalse();
  expect(manifest.views.map((view: any) => view.silhouette_path_count)).toEqual([8, 17, 27]);
  for (const view of manifest.views) {
    expect(view.drawing_line_source_kind).toBe(sourceKind);
    expect(view.drawing_line_source_label_zh).toBe(sourceLabelZh);
    expect(view.official_cad_path_count).toBe(0);
    expect(view.blue_line_present).toBeFalse();
    const svg = readFileSync(join(root, view.svg), "utf8");
    expect(svg).toContain('class="simplified-proxy-silhouette geometry-derived"');
    expect(svg).toContain('data-source-kind="geometry_derived_simplified_proxy"');
    expect(svg).not.toContain('class="official-reference native-dwg"');
    expect(svg).not.toContain("#1677c8");
  }
});

test("project drawings retain walls and furniture with masked black top-layer proxies", () => {
  const context = JSON.parse(readFileSync(join(product, "project-context-manifest.json"), "utf8"));
  expect(context.source_kind).toBe(sourceKind);
  expect(context.source_label_zh).toBe(sourceLabelZh);
  expect(context.official_cad_used).toBeFalse();
  expect(context.third_party_cad_used).toBeFalse();
  expect(context.project_context_retained).toBeTrue();
  expect(context.walls_and_surrounding_project_elements_retained).toBeTrue();
  expect(context.overlay_top_layer_with_white_mask).toBeTrue();
  expect(context.blue_line_present).toBeFalse();
  expect(context.pass).toBeTrue();
  for (const view of context.views) {
    expect(view.overlay.source_kind).toBe(sourceKind);
    expect(view.overlay.official_cad_used).toBeFalse();
    expect(view.overlay.fit.uniform_scale_preserved).toBeTrue();
    expect(view.overlay.fit.scale_svg_units_per_mm).toBe(0.02);
    expect(view.overlay.fit.transformation_mode).toBe(
      "axis_swap_rigid_reflection_and_translation_only",
    );
    expect(sha256(join(root, view.review_preview))).toBe(view.review_preview_sha256);
    expect(Math.max(...view.overlay.fit.bbox_absolute_delta_svg_units)).toBeLessThanOrEqual(
      view.overlay.fit.bbox_tolerance_svg_units,
    );
    const svg = readFileSync(join(root, view.output), "utf8");
    expect(svg).toContain("IfcFurniture");
    if (view.view !== "plan") expect(svg).toContain("IfcWall");
    expect(svg).toContain('class="geometry-derived-simplified-proxy project-context-overlay"');
    expect(svg).toContain('class="geometry-derived-proxy-mask"');
    expect(svg).toContain('class="geometry-derived-proxy"');
    expect(svg.indexOf('class="geometry-derived-proxy-mask"')).toBeLessThan(
      svg.indexOf('class="geometry-derived-proxy"'),
    );
    expect(svg).not.toContain('class="official-reference native-dwg"');
  }
});

test("Bonsai evidence contains unclipped actual IFC Body camera renders", () => {
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
});

test("pending Hima review leaves the formal IFC byte-identical and creates no derived IFC", () => {
  const candidate = JSON.parse(readFileSync(join(product, "candidate-representations.json"), "utf8"));
  expect(candidate.formal_ifc_write_allowed).toBeFalse();
  expect(candidate.review_status).toBe("visual_review_pending");
  expect(existsSync(join(product, "Poliform-Hima-HIMA01-derived-drawing.ifc"))).toBeFalse();
  expect(sha256(formal)).toBe(formalHash);
});

test("writer rejects missing apply and the pending human approval record", () => {
  const temporary = mkdtempSync(join(tmpdir(), "hima01-gate-"));
  const output = join(temporary, "forbidden.ifc");
  const script = join(root, "pipeline/scripts/hima01_drawing_ifc.py");
  const approval = JSON.parse(readFileSync(join(root, "pipeline/decisions/hima01-drawing-approval.json"), "utf8"));
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

test("a scoped temporary approval writes verified Hima representations and source associations", () => {
  const manifest = join(product, "manifest.json");
  const temporary = mkdtempSync(join(tmpdir(), "hima01-approved-"));
  const approvalPath = join(temporary, "approval.json");
  const output = join(temporary, "derived.ifc");
  const report = join(temporary, "report.json");
  writeFileSync(approvalPath, JSON.stringify({
    schema_version: 1,
    profile_key: "hima01",
    ifc_type_name: "HIMA01",
    candidate_manifest_sha256: sha256(manifest),
    status: "approved",
    reviewer: "automated gate fixture",
    review_date: "2026-08-23",
    approved_views: ["plan", "front", "side"],
    derived_ifc_write_allowed: true,
    formal_authoritative_ifc_write_allowed: false,
    scope: "manufacturer family identity and nominal dimensions only; not official CAD geometry and not a project shop drawing",
    approval_evidence: "temporary automated writer verification only",
  }, null, 2));
  const run = Bun.spawnSync([
    "python3",
    join(root, "pipeline/scripts/hima01_drawing_ifc.py"),
    "--input", formal,
    "--manifest", manifest,
    "--approval", approvalPath,
    "--output", output,
    "--report", report,
    "--apply",
  ], { cwd: root, stdout: "pipe", stderr: "pipe" });
  expect(run.exitCode).toBe(0);
  if (run.exitCode !== 0) throw new Error(run.stderr.toString());
  expect(existsSync(output)).toBeTrue();
  const result = JSON.parse(readFileSync(report, "utf8"));
  expect(result.pass).toBeTrue();
  expect(result.formal_ifc_bytes_unchanged).toBeTrue();
  expect(result.representations).toEqual({
    plan: "Hima01Plan",
    front: "Hima01Front",
    side: "Hima01Side",
  });
  expect(result.representation_path_counts).toEqual({ plan: 8, front: 17, side: 27 });
  expect(result.source_kind).toBe(sourceKind);
  expect(result.source_label_zh).toBe(sourceLabelZh);
  expect(result.official_cad_geometry_included).toBeFalse();
  expect(result.source_document_associations).toEqual([
    "POLIFORM-HIMA-OFFICIAL-NEWS-2022-TECHNICAL-DATA",
    "POLIFORM-HIMA-OFFICIAL-PAGE-EVIDENCE",
    "POLIFORM-HIMA-OFFICIAL-PRODUCT-PAGE",
    "POLIFORM-HIMA-OFFICIAL-TECHNICAL-SHEET",
    "POLIFORM-HIMA-SOURCE-ACCESS-RECORD",
  ]);
  expect(sha256(formal)).toBe(formalHash);
  rmSync(temporary, { recursive: true, force: true });
}, 30_000);
