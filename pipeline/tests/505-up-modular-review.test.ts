import { describe, expect, test } from "bun:test";
import { readFileSync } from "node:fs";
import { join } from "node:path";

const root = join(import.meta.dir, "../..");
const out = join(root, "output/review/highpoly-types/505-up-v1-lp-s");

describe("505 UP modular review evidence", () => {
  const audit = JSON.parse(readFileSync(join(out, "505-up-modular-rule-comparison.json"), "utf8"));
  const manifest = JSON.parse(readFileSync(join(out, "505-up-modular-review-manifest.json"), "utf8"));

  test("mechanically decomposes the project width and height", () => {
    expect(audit.checks.width_module_sum.status).toBe("pass");
    expect(audit.checks.width_module_sum.delta_mm).toBe(0);
    expect(audit.checks.body_height_modules.status).toBe("pass");
  });

  test("does not hide the outstanding compliance failures", () => {
    expect(audit.checks.floor_clearance.status).toBe("fail");
    expect(audit.checks.right_display_depth.status).toBe("fail");
    expect(audit.checks.left_grille_identity.official_blue_allowed).toBe(false);
  });

  test("uses the native component without stretching", () => {
    expect(audit.native_dwg_atom.native_size_mm[0]).toBeCloseTo(608, 6);
    expect(audit.native_dwg_atom.native_size_mm[1]).toBeCloseTo(768, 6);
    expect(audit.native_dwg_atom.scale_x).toBe(1);
    expect(audit.native_dwg_atom.scale_y).toBe(1);
    expect(audit.checks.native_atom_transform.placement).toBe("translation_only");
  });

  test("keeps formal IFC outside the prototype boundary", () => {
    expect(manifest.formal_ifc_sha256).toBe("7a50b87e8f48a7c2bfbaa9f1dabdfca7155fa684b2325c66aa1ae15c8ab0a25c");
    expect(manifest.formal_ifc_unchanged).toBe(true);
    expect(manifest.derived_ifc_written).toBe(false);
  });

  test("emits non-empty review SVGs with explicit source semantics", () => {
    const overlay = readFileSync(join(out, "505-up-modular-official-component-overlay-front.svg"), "utf8");
    const atom = readFileSync(join(out, "official-dwg-atom-w608-h768.svg"), "utf8");
    expect(overlay).toContain("actual IFC Body");
    expect(overlay).toContain("native DWG atom");
    expect(overlay.length).toBeGreaterThan(10000);
    expect(atom).toContain("scale 1:1");
  });
});
