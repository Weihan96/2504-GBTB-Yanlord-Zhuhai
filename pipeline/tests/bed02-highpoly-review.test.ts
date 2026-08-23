import { expect, test } from "bun:test";
import { createHash } from "node:crypto";
import { existsSync, mkdtempSync, readFileSync, rmSync } from "node:fs";
import { tmpdir } from "node:os";
import { join, resolve } from "node:path";

const root = resolve(import.meta.dir, "../..");
const product = join(root, "output/review/highpoly-types/bed02");
const formal = join(root, "2504 GBTB Yanlord Zhuhai.ifc");
const formalHash = "7a50b87e8f48a7c2bfbaa9f1dabdfca7155fa684b2325c66aa1ae15c8ab0a25c";
const sourceKind = "geometry_derived_simplified_proxy";
const sourceLabelZh = "基于原始高模几何生成的简化图纸表达";

function sha256(path: string) {
  return createHash("sha256").update(readFileSync(path)).digest("hex");
}

test("official Baxter identity is archived while CAD and dimension differences remain explicit", () => {
  const access = JSON.parse(readFileSync(join(product, "official-source/source-access-record.json"), "utf8"));
  expect(access.manufacturer).toBe("Baxter");
  expect(access.family).toBe("Viktor");
  expect(access.designer).toBe("Draga & Aurel");
  expect(access.project_type_code).toBe("BED02");
  expect(access.official_product_cad.authentication_required).toBeTrue();
  expect(access.official_product_cad.acquired).toBeFalse();
  expect(access.official_product_cad.local_cad_files).toEqual([]);
  expect(access.drawing_geometry_source).toMatchObject({
    source_kind: sourceKind,
    source_label_zh: sourceLabelZh,
    official_cad_used: false,
    third_party_cad_used: false,
  });
  expect(access.dimension_cross_check).toMatchObject({
    official_160x200_variant_overall_mm: [1720, 2340, 1060],
    project_ifc_description_mm: [1620, 2340, 1060],
    project_ifc_body_bounds_mm: [1675.710938, 2404.0177, 1107.145996],
    body_minus_official_mm: [-44.289062, 64.0177, 47.145996],
    status: "review_required_not_scaled_or_corrected",
  });
  const expectedHashes: Record<string, string> = {
    "baxter-viktor-official-product-page.html": "c2d4ae34a31626853063c8e7f445029f0b55ea0c1a2ffea5d215f3165ee0f9f8",
    "Baxter_Viktor_technical-sheet.pdf": "b762451b409c057df396735d2bdb90f3ef52b2642cff628b322cf20a927ae7f5",
    "Baxter_Viktor_assembly-instructions.pdf": "d5aeab11bcd2edeeb4274e9459c97af4044ec4049be87c7028ed6f4e5f72a385",
  };
  for (const [name, hash] of Object.entries(expectedHashes)) {
    expect(sha256(join(product, "official-source", name))).toBe(hash);
  }
  for (const source of access.official_identity_sources) {
    expect(sha256(join(root, source.local_path))).toBe(source.sha256);
  }
});

test("single-instance BED02 review has three geometry-derived views and zero blue CAD paths", () => {
  const manifest = JSON.parse(readFileSync(join(product, "manifest.json"), "utf8"));
  const candidate = JSON.parse(readFileSync(join(product, "candidate-representations.json"), "utf8"));
  const inventory = JSON.parse(readFileSync(join(root, "pipeline/decisions/highpoly-product-inventory.json"), "utf8"));
  const inventoryEntry = inventory.products.find((entry: any) => entry.type_name === "BED02");
  expect(inventoryEntry.folder).toBe("output/review/highpoly-types/bed02");
  expect(inventoryEntry.status).toBe("review_ready_pending_approval");
  expect(manifest.representative_global_id).toBe("1i_pqgLv9A7uuV7MjaArBW");
  expect(manifest.registered_instance_global_ids).toEqual(["1i_pqgLv9A7uuV7MjaArBW"]);
  expect(manifest.geometry_product_count).toBe(1);
  expect(manifest.whole_model_render).toBeFalse();
  expect(manifest.formal_ifc_bytes_unchanged).toBeTrue();
  expect(manifest.article_number).toBe("Baxter Viktor / BED02");
  expect(candidate.source_kind).toBe(sourceKind);
  expect(candidate.source_label_zh).toBe(sourceLabelZh);
  expect(candidate.official_cad_used).toBeFalse();
  expect(candidate.third_party_cad_used).toBeFalse();
  expect(manifest.views.map((view: any) => view.silhouette_path_count)).toEqual([1, 1, 4]);
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

test("BED02 project furniture plan retains room and furniture context with a uniform unforced overlay", () => {
  const context = JSON.parse(readFileSync(join(product, "project-context-manifest.json"), "utf8"));
  expect(context.source_kind).toBe(sourceKind);
  expect(context.source_label_zh).toBe(sourceLabelZh);
  expect(context.official_cad_used).toBeFalse();
  expect(context.third_party_cad_used).toBeFalse();
  expect(context.project_context_retained).toBeTrue();
  expect(context.walls_and_surrounding_project_elements_retained).toBeTrue();
  expect(context.overlay_top_layer_with_white_mask).toBeTrue();
  expect(context.blue_line_present).toBeFalse();
  expect(context.project_scale_svg_units_per_mm).toBe(0.02);
  expect(context.pass).toBeTrue();
  expect(context.views.map((view: any) => [view.view, view.candidate_view])).toEqual([["plan", "plan"]]);
  const view = context.views[0];
  expect(view.overlay.fit.alignment_mode).toBe("centre");
  expect(view.overlay.fit.scale_svg_units_per_mm).toBe(0.02);
  expect(Math.max(...view.overlay.fit.bbox_absolute_delta_svg_units)).toBeLessThanOrEqual(
    view.overlay.fit.bbox_tolerance_svg_units,
  );
  const svg = readFileSync(join(root, view.output), "utf8");
  expect(svg).toContain("IfcSpace");
  expect(svg).toContain("IfcFurniture");
  expect(svg).toContain('class="geometry-derived-simplified-proxy project-context-overlay"');
  expect(svg).toContain('class="geometry-derived-proxy-mask"');
  expect(svg).toContain('class="geometry-derived-proxy"');
  expect(svg.indexOf('class="geometry-derived-proxy-mask"')).toBeLessThan(
    svg.indexOf('class="geometry-derived-proxy"'),
  );
  expect(svg).not.toContain('class="official-reference native-dwg"');
});

test("BED02 Bonsai evidence contains unclipped actual IFC Body camera renders", () => {
  const evidence = JSON.parse(readFileSync(join(product, "bonsai-review-manifest.json"), "utf8"));
  expect(evidence.mode).toBe("actual_bonsai_ifc_body_camera_render");
  expect(evidence.representative_global_id).toBe("1i_pqgLv9A7uuV7MjaArBW");
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
    expect(render.camera_ortho_scale_m + 1e-7).toBeGreaterThanOrEqual(
      render.projected_width_m * render.framing_margin_factor,
    );
    expect(render.camera_ortho_scale_m + 1e-7).toBeGreaterThanOrEqual(
      render.projected_height_m * render.image_aspect * render.framing_margin_factor,
    );
    expect(sha256(join(root, render.path))).toBe(render.sha256);
  }
  expect(sha256(join(root, evidence.bonsai_session.path))).toBe(evidence.bonsai_session.sha256);
  expect(sha256(join(root, evidence.isolated_ifc))).toBe(evidence.isolated_ifc_sha256);
});

test("pending BED02 review leaves the formal IFC byte-identical and creates no derived IFC", () => {
  const candidate = JSON.parse(readFileSync(join(product, "candidate-representations.json"), "utf8"));
  expect(candidate.formal_ifc_write_allowed).toBeFalse();
  expect(candidate.review_status).toBe("visual_review_pending");
  expect(existsSync(join(product, "Baxter-Viktor-BED02-derived-drawing.ifc"))).toBeFalse();
  expect(sha256(formal)).toBe(formalHash);
});

test("BED02 writer rejects missing apply and the pending human approval record", () => {
  const temporary = mkdtempSync(join(tmpdir(), "bed02-gate-"));
  const output = join(temporary, "forbidden.ifc");
  const script = join(root, "pipeline/scripts/bed02_drawing_ifc.py");
  const approval = JSON.parse(readFileSync(join(root, "pipeline/decisions/bed02-drawing-approval.json"), "utf8"));
  const manifest = join(product, "manifest.json");
  expect(approval.status).toBe("pending");
  expect(approval.derived_ifc_write_allowed).toBeFalse();
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
