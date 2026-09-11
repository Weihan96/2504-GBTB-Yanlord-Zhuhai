import { expect, test } from "bun:test";
import { createHash } from "node:crypto";
import { existsSync, mkdtempSync, readFileSync, rmSync, writeFileSync } from "node:fs";
import { tmpdir } from "node:os";
import { join, resolve } from "node:path";

const root = resolve(import.meta.dir, "../..");
const product = join(root, "output/review/highpoly-types/geberit-duofix-sigma-224-212");
const formal = join(root, "2504 GBTB Yanlord Zhuhai.ifc");
const formalHash = "7a50b87e8f48a7c2bfbaa9f1dabdfca7155fa684b2325c66aa1ae15c8ab0a25c";
const sourceHashes: Record<string, string> = {
  A: "f06952ac7bf91f0a0037ae85644df3df0938d549322cbfb59b7dfd9be01939d4",
  G: "841dd3c72d594cfd2c7f921b84f10a2305ab9957b531771662fb7f03eff97207",
  L: "e6a26c660ae9cd4fc8cb4b1effadcea38a4794b66523f0fea2c02c5e02b77f25",
  P: "5c92a663de35357a4bf0ff7eb0838bc5e178030afe7d769ead5632f50821652f",
};

function sha256(path: string) {
  return createHash("sha256").update(readFileSync(path)).digest("hex");
}

test("exact project-selected 224.212.00.2 official DWGs retain identity", () => {
  const profile = JSON.parse(readFileSync(join(product, "profile.json"), "utf8"));
  const reference = profile.profiles["geberit-duofix-sigma-224-212"].official_reference;
  expect(reference.catalogue_status).toContain("B=50 cm, H=112 cm, T=12 cm");
  expect(sha256(join(root, reference.local_catalogue_extract))).toBe(reference.local_catalogue_extract_sha256);
  expect(sha256(join(root, reference.project_received_box_label))).toBe(reference.project_received_box_label_sha256);
  const inventory = JSON.parse(
    readFileSync(join(root, "pipeline/decisions/highpoly-product-inventory.json"), "utf8"),
  );
  const inventoryEntry = inventory.products.find(
    (entry: any) => entry.type_name === "Geberit Duofix Sigma",
  );
  expect(inventoryEntry.folder).toBe(
    "output/review/highpoly-types/geberit-duofix-sigma-224-212",
  );
  expect(inventoryEntry.status).toBe("review_ready_pending_approval");

  const linework = JSON.parse(readFileSync(join(product, "official-native-dwg-linework.json"), "utf8"));
  expect(linework.source_kind).toBe("native_dwg");
  expect(linework.article_number).toBe("224.212.00.2");
  expect(linework.ifc_type_name).toBe("Geberit Duofix Sigma");
  expect(linework.pass).toBeTrue();
  expect(linework.nominal_dimensions_mm).toEqual({ width: 500, height: 1120, cistern_depth: 120 });
  expect(linework.three_view_code_mapping).toMatchObject({
    G: "Grundriss / plan",
    A: "Ansicht / front elevation",
    L: "left side elevation",
  });
  for (const [code, expected] of Object.entries(sourceHashes)) {
    const source = linework.official_sources[code];
    expect(source.sha256).toBe(expected);
    expect(source.url).toBe(`https://cdn.data.geberit.com/cad/224.212.00.2_${code}.dwg`);
    expect(sha256(join(root, source.path))).toBe(expected);
  }
  const access = JSON.parse(
    readFileSync(join(product, "official-source/source-access-record.json"), "utf8"),
  );
  expect(access).toMatchObject({
    manufacturer: "Geberit",
    article_number: "224.212.00.2",
    source_kind: "native_dwg",
    official_cad_acquired: true,
    official_cad_exact_project_configuration_match: true,
    official_cad_used: true,
    third_party_cad_used: false,
    blue_line_present: true,
    representation_geometry_source: "official_native_dwg_paths_mm",
    pass: true,
  });
  expect(access.drawing_view_mapping).toEqual({ plan: "G", front: "A", side: "L" });
  expect(access.official_3d_identity_only.used_as_plan_or_elevation_geometry).toBeFalse();
  expect(access.source_label_zh).toBe("基于 Geberit 精确型号原生 DWG 的官方图纸表达");
  expect(access.official_source_revalidation).toMatchObject({
    all_declared_dwg_urls_accessible: true,
    all_downloaded_dwg_bytes_match_local_archive: true,
    catalogue_page_render_matches_local_extract: true,
    pass: true,
  });
  expect(sha256(join(root, access.official_source_revalidation.path))).toBe(
    access.official_source_revalidation.sha256,
  );
  for (const [code, expected] of Object.entries(sourceHashes)) {
    expect(access.official_native_dwg[code].sha256).toBe(expected);
    expect(sha256(join(root, access.official_native_dwg[code].path))).toBe(expected);
  }
  expect(sha256(join(root, linework.local_official_catalogue_extract.path))).toBe(
    linework.local_official_catalogue_extract.sha256,
  );
  expect(sha256(join(root, linework.received_box_label.path))).toBe(linework.received_box_label.sha256);
});

test("live official sources preserve exact DWG bytes and catalogue article identity", () => {
  const evidence = JSON.parse(
    readFileSync(join(product, "official-source/official-source-revalidation.json"), "utf8"),
  );
  expect(evidence).toMatchObject({
    manufacturer: "Geberit",
    article_number: "224.212.00.2",
    drawing_view_mapping: { plan: "G", front: "A", side: "L" },
    identity_only_code: "P",
    all_declared_dwg_urls_accessible: true,
    all_downloaded_dwg_bytes_match_local_archive: true,
    pass: true,
  });
  expect(evidence.catalogue).toMatchObject({
    http_status: 200,
    content_type: "application/pdf",
    article_page_pdf_index_one_based: 10,
    printed_page_number: 11,
    expected_local_extract_sha256: "8c3a17569e6fa53d9a769ab584fed79a865118c4a27018facd619b4a5cc5cdd3",
    expected_article_page_render_sha256: "757450d7fb7a85d61755fe61dce34e77d14ecc33e7522da056d8f97db3fee660",
    downloaded_article_page_render_sha256: "757450d7fb7a85d61755fe61dce34e77d14ecc33e7522da056d8f97db3fee660",
    local_extract_render_sha256: "757450d7fb7a85d61755fe61dce34e77d14ecc33e7522da056d8f97db3fee660",
    downloaded_article_page_render_matches_local_extract: true,
    pass: true,
  });
  expect(evidence.catalogue.identity_statement).toContain(
    "224.212.00.2 with B=50 cm, H=112 cm and T=12 cm",
  );
  for (const [code, expected] of Object.entries(sourceHashes)) {
    expect(evidence.dwg_results[code]).toMatchObject({
      url: `https://cdn.data.geberit.com/cad/224.212.00.2_${code}.dwg`,
      http_status: 200,
      expected_sha256: expected,
      downloaded_sha256: expected,
      local_sha256: expected,
      downloaded_bytes_match_local_archive: true,
      pass: true,
    });
  }
});

test("official G/A/L views mechanically cross-check one isolated IFC Body", () => {
  const manifest = JSON.parse(readFileSync(join(product, "manifest.json"), "utf8"));
  expect(manifest).toMatchObject({
    source_kind: "native_dwg",
    official_cad_acquired: true,
    official_cad_exact_project_configuration_match: true,
    official_cad_used: true,
    third_party_cad_used: false,
    blue_line_present: true,
    representation_geometry_source: "official_native_dwg_paths_mm",
    simplified_proxy_comparison_included: true,
  });
  expect(sha256(join(root, manifest.official_source_access_record))).toBe(
    manifest.official_source_access_record_sha256,
  );
  expect(manifest.geometry_product_count).toBe(1);
  expect(manifest.whole_model_render).toBeFalse();
  expect(manifest.formal_ifc_bytes_unchanged).toBeTrue();
  expect(manifest.registered_instance_global_ids).toEqual([
    "3hgNkx97vCTOC2eewCpMNk",
    "3pvAlH5C14v8uVEJ1LmK8M",
  ]);
  expect(manifest.views.map((view: any) => view.blue_line_native_dwg_code)).toEqual(["G", "A", "L"]);
  expect(manifest.views.map((view: any) => view.official_reference_path_count)).toEqual([253, 790, 225]);
  for (const view of manifest.views) {
    expect(view.blue_line_source_kind).toBe("native_dwg");
    expect(view.blue_line_article_number).toBe("224.212.00.2");
    expect(view.mechanical_cross_check.pass).toBeTrue();
    const svg = readFileSync(join(root, view.svg), "utf8");
    expect(svg).toContain('class="official-reference-mask"');
    expect(svg).toContain('class="official-reference native-dwg"');
    expect(svg).toContain('data-source-kind="native_dwg"');
    expect(svg).toContain('data-article-number="224.212.00.2"');
    expect(svg.indexOf('class="official-reference-mask"')).toBeLessThan(
      svg.indexOf('class="official-reference native-dwg"'),
    );
  }
  expect(manifest.views[0].mechanical_cross_check.absolute_delta_mm[1]).toBeLessThanOrEqual(20);
  expect(manifest.views[1].mechanical_cross_check.absolute_delta_mm[0]).toBeLessThanOrEqual(1);
  expect(manifest.views[1].mechanical_cross_check.absolute_delta_mm[1]).toBeLessThanOrEqual(1);
  expect(manifest.views[2].mechanical_cross_check.absolute_delta_mm[0]).toBeLessThanOrEqual(20);
  expect(manifest.views[2]).toMatchObject({
    blue_line_native_dwg_code: "L",
    view_direction: "left_side_local_negative_x",
    mechanical_cross_check: {
      handedness: {
        comparison: "asymmetric_path_centroids_share_left_view_handedness",
        view_direction: "left_side_local_negative_x",
        same_horizontal_side: true,
        pass: true,
      },
    },
  });
});

test("project context preserves uniform official scale, walls, and top-layer blue overlays", () => {
  const context = JSON.parse(readFileSync(join(product, "project-context-manifest.json"), "utf8"));
  expect(context.source_kind).toBe("native_dwg");
  expect(context.article_number).toBe("224.212.00.2");
  expect(context.project_context_retained).toBeTrue();
  expect(context.walls_and_surrounding_project_elements_retained).toBeTrue();
  expect(context.blue_line_top_layer_with_white_mask).toBeTrue();
  expect(context.official_dwg_uniform_project_scale_preserved).toBeTrue();
  expect(context.project_scale_svg_units_per_mm).toBe(0.02);
  expect(context.views.flatMap((view: any) => view.overlays)).toHaveLength(4);
  for (const view of context.views) {
    const svg = readFileSync(join(root, view.output), "utf8");
    expect(svg).toContain("IfcWall");
    expect(svg).toContain("IfcSanitaryTerminal");
    expect(svg).toContain('class="official-native-dwg project-context-overlay"');
    expect(svg).toContain('data-article-number="224.212.00.2"');
    for (const overlay of view.overlays) {
      expect(overlay.fit.uniform_scale_preserved).toBeTrue();
      expect(overlay.fit.scale_svg_units_per_mm).toBe(0.02);
      expect(overlay.fit.transformation_mode).toBe("axis_swap_rigid_reflection_and_translation_only");
    }
    if (view.view === "plan") {
      for (const crop of view.review_crops) {
        expect(sha256(join(root, crop.preview))).toBe(crop.preview_sha256);
      }
    } else {
      expect(sha256(join(root, view.review_preview))).toBe(view.review_preview_sha256);
    }
  }
});

test("Bonsai evidence is an unclipped actual camera render of the isolated IFC Body", () => {
  const evidence = JSON.parse(readFileSync(join(product, "bonsai-review-manifest.json"), "utf8"));
  expect(evidence.mode).toBe("actual_bonsai_ifc_body_camera_render");
  expect(evidence.geometry_product_count).toBe(1);
  expect(evidence.whole_model_render).toBeFalse();
  expect(evidence.formal_ifc_bytes_unchanged).toBeTrue();
  expect(evidence.bonsai_session.saved_active_representation).toBe("Body");
  expect(evidence.bonsai_session.ifc_context_identifier).toBe("Body");
  expect(evidence.bonsai_session.saved_camera_count).toBe(4);
  expect(evidence.bonsai_session.side_camera_local_x_sign).toBe(-1);
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

test("approved Duofix review still leaves the formal IFC byte-identical and creates no derived IFC", () => {
  const candidate = JSON.parse(readFileSync(join(product, "candidate-representations.json"), "utf8"));
  expect(candidate).toMatchObject({
    source_kind: "native_dwg",
    official_cad_used: true,
    third_party_cad_used: false,
    blue_line_present: true,
    representation_geometry_source: "official_native_dwg_paths_mm",
    simplified_proxy_comparison_included: true,
    simplified_proxy_comparison_source_kind: "geometry_derived_simplified_proxy",
  });
  expect(candidate.formal_ifc_write_allowed).toBeFalse();
  expect(candidate.review_status).toBe("visual_review_pending");
  expect(existsSync(join(product, "Geberit-Duofix-Sigma-224-212-derived-drawing.ifc"))).toBeFalse();
  expect(sha256(formal)).toBe(formalHash);
});

test("approved writer gate still requires an explicit apply operation", () => {
  const temporary = mkdtempSync(join(tmpdir(), "duofix-224212-gate-"));
  const output = join(temporary, "forbidden.ifc");
  const script = join(root, "pipeline/scripts/geberit_duofix_sigma_224_212_drawing_ifc.py");
  const approval = JSON.parse(
    readFileSync(join(root, "pipeline/decisions/geberit-duofix-sigma-224-212-drawing-approval.json"), "utf8"),
  );
  const manifest = join(product, "manifest.json");
  expect(approval.status).toBe("approved");
  expect(approval.derived_ifc_write_allowed).toBeTrue();
  expect(approval.formal_authoritative_ifc_write_allowed).toBeFalse();
  expect(approval.approved_views).toEqual(["plan", "front", "side"]);
  expect(approval.candidate_manifest_sha256).toBe(sha256(manifest));

  const withoutApply = Bun.spawnSync(["python3", script, "--input", formal, "--output", output], {
    cwd: root,
    stderr: "pipe",
  });
  expect(withoutApply.exitCode).not.toBe(0);
  expect(withoutApply.stderr.toString()).toContain("IFC write requires the explicit --apply flag");
  expect(existsSync(output)).toBeFalse();

  expect(existsSync(output)).toBeFalse();
  expect(sha256(formal)).toBe(formalHash);
  rmSync(temporary, { recursive: true, force: true });
});

test("a scoped temporary approval writes verified native-DWG representations and source associations", () => {
  const manifest = join(product, "manifest.json");
  const temporary = mkdtempSync(join(tmpdir(), "duofix-224212-approved-"));
  const approvalPath = join(temporary, "approval.json");
  const output = join(temporary, "derived.ifc");
  const report = join(temporary, "report.json");
  writeFileSync(approvalPath, JSON.stringify({
    schema_version: 1,
    profile_key: "geberit-duofix-sigma-224-212",
    article_number: "224.212.00.2",
    candidate_manifest_sha256: sha256(manifest),
    status: "approved",
    reviewer: "automated gate fixture",
    review_date: "2026-08-23",
    approved_views: ["plan", "front", "side"],
    derived_ifc_write_allowed: true,
    formal_authoritative_ifc_write_allowed: false,
    scope: "official manufacturer exact article family reference, not a project shop drawing",
    approval_evidence: "temporary automated writer verification only",
  }, null, 2));
  const run = Bun.spawnSync([
    "python3",
    join(root, "pipeline/scripts/geberit_duofix_sigma_224_212_drawing_ifc.py"),
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
    plan: "Duofix224212Plan",
    front: "Duofix224212Front",
    side: "Duofix224212Side",
  });
  expect(result.representation_path_counts).toEqual({ plan: 253, front: 790, side: 225 });
  expect(result.representation_geometry_source).toBe("official_native_dwg_paths_mm");
  expect(result.proxy_geometry_included).toBeFalse();
  expect(result.source_document_associations).toEqual([
    "GEBERIT-224-212-A-NATIVE-DWG",
    "GEBERIT-224-212-G-NATIVE-DWG",
    "GEBERIT-224-212-L-NATIVE-DWG",
    "GEBERIT-224-212-OFFICIAL-CATALOGUE-EXTRACT",
    "GEBERIT-224-212-P-NATIVE-DWG-IDENTITY",
    "GEBERIT-224-212-PROJECT-RECEIVED-BOX-LABEL",
  ]);
  expect(sha256(formal)).toBe(formalHash);
  rmSync(temporary, { recursive: true, force: true });
}, 30_000);
