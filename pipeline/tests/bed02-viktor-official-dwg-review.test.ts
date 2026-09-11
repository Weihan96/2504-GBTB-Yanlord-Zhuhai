import { describe, expect, test } from "bun:test";
import { createHash } from "node:crypto";
import { readFileSync } from "node:fs";
import { join } from "node:path";

const root = join(import.meta.dir, "../..");
const folder = join(root, "output/review/highpoly-types/bed02");
const formalIfc = join(root, "2504 GBTB Yanlord Zhuhai.ifc");
const sourceDwg = join(folder, "official-source/official-download/Viktor_Letto.dwg");
const formalSha = "7a50b87e8f48a7c2bfbaa9f1dabdfca7155fa684b2325c66aa1ae15c8ab0a25c";
const sourceSha = "3f0676a004f3e744779093188d75d58238153332d9d1810f28821cafc79a9bbd";

function json(name: string) {
  return JSON.parse(readFileSync(join(folder, name), "utf8"));
}

function sha256(path: string) {
  return createHash("sha256").update(readFileSync(path)).digest("hex");
}

function closeTo(actual: number[], expected: number[], tolerance = 0.001) {
  expect(actual).toHaveLength(expected.length);
  actual.forEach((value, index) => expect(Math.abs(value - expected[index])).toBeLessThanOrEqual(tolerance));
}

describe("BED02 corrected Viktor native-DWG review", () => {
  test("keeps the immutable source DWG and formal IFC", () => {
    expect(sha256(sourceDwg)).toBe(sourceSha);
    expect(sha256(formalIfc)).toBe(formalSha);
  });

  test("selects geometry from ARREDO and layer 0 while excluding dimensions and text", () => {
    const reference = json("official-dwg-review-reference.json");
    expect(reference.generator).toBe("pipeline/scripts/repair_bed02_viktor_official_dwg_review.py");
    expect(reference.source_dwg_sha256).toBe(sourceSha);
    expect(reference.selected_geometry_layers).toEqual(["ARREDO", "0"]);
    expect(reference.excluded_layers._QUOTE).toContain("dimensions");
    expect(reference.excluded_entity_types.MTEXT).toContain("false Front path");
    expect(reference.source_geometry_scaled).toBe(false);
    expect(reference.source_geometry_anisotropically_fitted).toBe(false);
    expect(reference.primary_candidate_replaced).toBe(false);
  });

  test("restores the complete 160x200 Plan Front and Side paths", () => {
    const reference = json("official-dwg-review-reference.json");
    const expected = {
      plan: { count: 51, layers: { ARREDO: 28, "0": 23 }, size: [1723.392281, 2340.708974] },
      front: { count: 99, layers: { ARREDO: 91, "0": 8 }, size: [1733.666161, 1060] },
      side: { count: 84, layers: { ARREDO: 77, "0": 7 }, size: [2351.640392, 1065.443243] },
    } as const;
    for (const view of ["plan", "front", "side"] as const) {
      const record = reference.views[view];
      expect(record.path_count).toBe(expected[view].count);
      expect(record.path_count_by_source_layer).toEqual(expected[view].layers);
      closeTo(record.bounds_mm.size, [...expected[view].size]);
      expect(record.blue_stroke_style).toBe("solid");
      expect(record.source_scaled).toBe(false);
      expect(record.paths_mm.length).toBe(expected[view].count);
    }
  });

  test("records the ARREDO-only SVGs as superseded extraction, not a DWG error", () => {
    const reference = json("official-dwg-review-reference.json");
    const old = reference.superseded_extractions[0];
    expect(old.status).toBe("superseded_incomplete_extraction_excluded_from_approval");
    expect(old.source_dwg_valid).toBe(true);
    expect(old.selected_geometry_layers).toEqual(["ARREDO"]);
    expect(old.reason).toContain("DWG is correct");
    expect(old.views.plan.path_count).toBe(28);
    closeTo(old.views.plan.bounds_mm.size, [1715.424395, 2169.757421]);
  });

  test("renders solid blue native paths with only rigid alignment", () => {
    const reference = json("official-dwg-review-reference.json");
    const packageManifest = json("official-dwg-review-manifest.json");
    const expectedCounts = { plan: 51, front: 99, side: 84 };
    expect(packageManifest.selected_geometry_layers).toEqual(["ARREDO", "0"]);
    for (const record of packageManifest.views) {
      const view = record.view as keyof typeof expectedCounts;
      expect(record.official_reference_path_count).toBe(expectedCounts[view]);
      expect(record.blue_stroke_style).toBe("solid");
      expect(record.alignment.uniform_scale).toBe(1);
      expect(record.alignment.anisotropic_scale_used).toBe(false);
      expect(record.alignment.source_geometry_deformed).toBe(false);
      expect(record.alignment.view_direction_reflection_x).toBe(view === "side");
      closeTo(record.alignment.source_bounds_before_alignment_mm.size, reference.views[view].bounds_mm.size);

      const svg = readFileSync(join(root, record.svg), "utf8");
      expect(svg).toContain('class="official-reference native-dwg"');
      expect(svg).toContain('stroke="#1677c8"');
      expect(svg).toContain('data-source-scaled="false"');
      expect(svg).not.toContain('stroke-dasharray="11 7"');
      expect(svg.length).toBeGreaterThan(10_000);
      expect(sha256(join(root, record.svg))).toBe(record.svg_sha256);
    }
  });

  test("keeps the black high-poly proxy primary and all write gates closed", () => {
    const manifest = json("manifest.json");
    const candidate = json("candidate-representations.json");
    const packageManifest = json("official-dwg-review-manifest.json");
    expect(manifest.blue_line_present).toBe(true);
    expect(manifest.blue_line_role).toContain("review_reference_only");
    expect(manifest.source_kind).toBe("geometry_derived_simplified_proxy");
    expect(manifest.official_cad_used).toBe(false);
    expect(manifest.approved_for_drawing_ifc).toBe(false);
    expect(manifest.derived_ifc_write_allowed).toBe(false);
    expect(manifest.formal_ifc_write).toBe("not performed");
    expect(candidate.source_kind).toBe("geometry_derived_simplified_proxy");
    expect(candidate.formal_ifc_write_allowed).toBe(false);
    expect(packageManifest.derived_ifc_write_performed).toBe(false);
    expect(packageManifest.formal_authoritative_ifc_write_performed).toBe(false);
    expect(packageManifest.formal_ifc_sha256).toBe(formalSha);
    expect(packageManifest.formal_ifc_bytes_unchanged).toBe(true);
  });
});
