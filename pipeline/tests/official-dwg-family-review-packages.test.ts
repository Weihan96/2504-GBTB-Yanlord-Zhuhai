import { expect, test } from "bun:test";
import { createHash } from "node:crypto";
import { existsSync, readFileSync } from "node:fs";
import { join, resolve } from "node:path";

const root = resolve(import.meta.dir, "../..");
const formal = join(root, "2504 GBTB Yanlord Zhuhai.ifc");
const formalHash = "7a50b87e8f48a7c2bfbaa9f1dabdfca7155fa684b2325c66aa1ae15c8ab0a25c";
const sourceLabelZh = "基于原始高模几何生成的简化图纸表达";

const expected = {
  bed01: {
    source: "ab697db9c27448c62a9a77537d7cc4b6286577c335328113d703184be28d9c4a",
    kinds: [
      "native_dwg_2d_linework",
      "native_dwg_2d_linework",
      "native_dwg_2d_linework",
    ],
    exact: false,
  },
  bed02: {
    source: "3f0676a004f3e744779093188d75d58238153332d9d1810f28821cafc79a9bbd",
    kinds: ["native_dwg_2d_linework", "native_dwg_2d_linework", "native_dwg_2d_linework"],
    exact: false,
  },
  sis04: {
    source: "77f240f65f1ab68a168b363b18f8294122ec913259f93745a01d19d534678263",
    kinds: [
      "native_dwg_2d_linework",
      "native_dwg_dimension_envelope_not_dedicated_2d_view",
      "native_dwg_dimension_envelope_not_dedicated_2d_view",
    ],
    exact: true,
  },
  hima01: {
    source: "f0e1b980aa102f3e06fe602f9e7a9b40db45d76e45200c08e4b8b4a0eec1d523",
    kinds: [
      "native_dwg_2d_linework",
      "native_dwg_2d_linework",
      "native_dwg_dimension_envelope_not_dedicated_2d_view",
    ],
    exact: false,
  },
} as const;

function sha256(path: string) {
  return createHash("sha256").update(readFileSync(path)).digest("hex");
}

for (const [slug, gate] of Object.entries(expected)) {
  test(`${slug} official-DWG supplement preserves the primary candidate and IFC write gate`, () => {
    const folder = join(root, "output/review/highpoly-types", slug);
    const reference = JSON.parse(readFileSync(join(folder, "official-dwg-review-reference.json"), "utf8"));
    const review = JSON.parse(readFileSync(join(folder, "official-dwg-review-manifest.json"), "utf8"));
    const candidate = JSON.parse(readFileSync(join(folder, "candidate-representations.json"), "utf8"));

    expect(reference.source_dwg_sha256).toBe(gate.source);
    expect(sha256(join(root, reference.source_dwg))).toBe(gate.source);
    expect(reference.source_geometry_scaled).toBeFalse();
    expect(reference.source_geometry_anisotropically_fitted).toBeFalse();
    expect(reference.third_party_cad_used).toBeFalse();
    expect(reference.primary_candidate_replaced).toBeFalse();
    expect(reference.formal_ifc_bytes_unchanged).toBeTrue();
    expect(reference.official_download_url).toMatch(/^https:\/\//);
    expect(Object.keys(reference.views)).toEqual(["plan", "front", "side"]);
    expect(Object.values(reference.views).map((view: any) => view.geometry_kind)).toEqual(gate.kinds);
    for (const view of Object.values(reference.views) as any[]) {
      expect(view.path_count).toBeGreaterThan(0);
      expect(view.paths_mm.length).toBe(view.path_count);
      expect(view.source_scaled).toBeFalse();
    }

    expect(candidate.source_kind).toBe("geometry_derived_simplified_proxy");
    expect(candidate.source_label_zh).toBe(sourceLabelZh);
    expect(candidate.formal_ifc_write_allowed).toBeFalse();
    expect(review.primary_candidate_source_kind).toBe("geometry_derived_simplified_proxy");
    expect(review.primary_candidate_source_label_zh).toBe(sourceLabelZh);
    expect(review.primary_candidate_replaced).toBeFalse();
    expect(review.official_reference_used_as_primary_candidate).toBeFalse();
    expect(review.derived_ifc_write_performed).toBeFalse();
    expect(review.formal_authoritative_ifc_write_performed).toBeFalse();
    expect(review.formal_ifc_bytes_unchanged).toBeTrue();
    expect(review.views.map((view: any) => view.view)).toEqual(["plan", "front", "side"]);
    expect(review.views.map((view: any) => view.official_reference_geometry_kind)).toEqual(gate.kinds);
    for (const view of review.views) {
      expect(view.alignment.uniform_scale).toBe(1);
      expect(view.alignment.anisotropic_scale_used).toBeFalse();
      expect(view.alignment.source_geometry_deformed).toBeFalse();
      expect(view.primary_candidate_replaced).toBeFalse();
      expect(sha256(join(root, view.svg))).toBe(view.svg_sha256);
      const svg = readFileSync(join(root, view.svg), "utf8");
      expect(svg).toContain('class="actual-ifc-body"');
      expect(svg).toContain('class="geometry-derived-simplified-proxy"');
      expect(svg).toContain('class="official-reference native-dwg"');
      expect(svg).toContain('data-source-role="review_reference_only"');
      expect(svg).toContain('data-source-scaled="false"');
      expect(svg).toContain("Primary candidate unchanged");
      expect(svg).toContain("IFC write: none");
    }
    expect(sha256(join(root, review.official_dwg_review_reference))).toBe(
      review.official_dwg_review_reference_sha256,
    );
    expect(sha256(join(root, review.contact_sheet_png))).toBe(review.contact_sheet_png_sha256);
    expect(existsSync(join(root, review.index))).toBeTrue();
    expect(sha256(formal)).toBe(formalHash);
  });
}

test("product-specific source semantics remain explicit and non-interchangeable", () => {
  const bed01 = JSON.parse(
    readFileSync(join(root, "output/review/highpoly-types/bed01/official-dwg-review-reference.json"), "utf8"),
  );
  const bed02 = JSON.parse(
    readFileSync(join(root, "output/review/highpoly-types/bed02/official-dwg-review-reference.json"), "utf8"),
  );
  const sis04 = JSON.parse(
    readFileSync(join(root, "output/review/highpoly-types/sis04/official-dwg-review-reference.json"), "utf8"),
  );
  const hima01 = JSON.parse(
    readFileSync(join(root, "output/review/highpoly-types/hima01/official-dwg-review-reference.json"), "utf8"),
  );
  expect(bed01.note).toContain("independent native 2D file Casablanca_Letto.dwg");
  expect(bed01.selected_geometry_layers).toEqual(["_ARREDO", "_MATERASSO", "_PIEDINI", "_CUCITURE"]);
  expect(bed01.excluded_layers._QUOTE).toContain("excluded from approval linework");
  expect(bed01.views.plan.path_count).toBe(25);
  expect(bed01.views.front.path_count).toBe(22);
  expect(bed01.views.side.path_count).toBe(36);
  expect(bed01.views.plan.blue_stroke_style).toBe("solid");
  expect(bed01.superseded_error_candidates).toHaveLength(1);
  expect(bed01.superseded_error_candidates[0].status).toBe(
    "superseded_error_candidate_excluded_from_approval",
  );
  expect(bed01.superseded_error_candidates[0].geometry_kind).toBe(
    "native_dwg_3d_solid_projected_envelope",
  );
  expect(bed01.autocad_source_screenshots).toHaveLength(2);
  for (const screenshot of bed01.autocad_source_screenshots) {
    expect(existsSync(join(root, screenshot))).toBeTrue();
  }
  expect(sha256(join(root, bed01.source_archive))).toBe(
    "acf67fddcfaa5e5bb50eb90f78a5eac262bfa318a767a938fc29f73931bb5970",
  );
  expect(bed01.source_archive_member_sha256).toBe(bed01.source_dwg_sha256);
  expect(sha256(join(root, bed01.source_archive_inventory))).toBe(
    bed01.source_archive_inventory_sha256,
  );
  const bed01Inventory = JSON.parse(readFileSync(join(root, bed01.source_archive_inventory), "utf8"));
  expect(bed01Inventory.member_count).toBe(13);
  expect(bed01Inventory.directory_count).toBe(3);
  expect(bed01Inventory.dwg_members).toEqual(["CASABLANCA/2D/Casablanca_Letto.dwg"]);
  expect(bed01.local_byte_identical_archive_copies).toHaveLength(3);
  expect(bed02.note).toContain("six family variants");
  expect(bed02.note).toContain("non-uniformly changed");
  expect(sis04.configuration).toContain("1962 x 370 x 722");
  expect(sis04.note).toContain("331 mm is internal");
  expect(hima01.configuration).toContain("PVA11");
  expect(hima01.note).toContain("folded to a different Plan envelope");
  expect(sha256(formal)).toBe(formalHash);
});
