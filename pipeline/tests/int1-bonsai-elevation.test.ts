import { expect, test } from "bun:test";
import { existsSync, readFileSync } from "node:fs";
import { join } from "node:path";

const root = process.cwd();
const manifestPath = join(root, "drawings/elevations/native/manifest.json");
const manifest = JSON.parse(readFileSync(manifestPath, "utf8"));

test("native Bonsai elevation manifest covers all registered views", () => {
  expect(manifest.pass).toBe(true);
  expect(manifest.view_count).toBe(36);
  expect(manifest.sheet_count).toBe(9);
  expect(manifest.demolish_wall_count).toBe(0);
  expect(manifest.linework_mode_counts).toEqual({ OPENCASCADE: 36 });
  expect(manifest.complexity_exclusion_occurrences).toBe(0);
  expect(manifest.complexity_exclusion_unique_global_ids).toEqual([]);
  // Eleven lightweight objects are visible in the original 36 views; the
  // twelfth is scoped only by the public-space P01/P02 batch.
  expect(manifest.lightweight_elevation_unique_global_ids).toHaveLength(11);
  expect(new Set(manifest.views.map((view: { view_id: string }) => view.view_id)).size).toBe(36);
});

test("every native view is vector-only and linked from the IFC manifest", () => {
  for (const view of manifest.views) {
    const svgPath = join(root, view.svg);
    expect(existsSync(svgPath)).toBe(true);
    const svg = readFileSync(svgPath, "utf8");
    expect(svg).toContain('data-scale="1:50"');
    expect(svg).toContain('id="noninteger-highlights"');
    expect(svg).not.toMatch(/<image\b/i);
    expect(view.target_view).toBe("ELEVATION_VIEW");
    expect(view.scale).toBe("1/50");
    expect(view.has_underlay).toBe(false);
    expect(view.demolish_wall_count).toBe(0);
  }
});

test("nine review sheets and phone PNGs are present", () => {
  for (let number = 1; number <= 9; number += 1) {
    const sheet = `EL-${String(number).padStart(2, "0")}`;
    const sheetSvg = join(root, `drawings/elevations/${sheet}-bonsai-native.svg`);
    const phonePng = join(root, `output/images/elevations/${sheet}-bonsai-native.png`);
    expect(existsSync(sheetSvg)).toBe(true);
    expect(existsSync(phonePng)).toBe(true);
    const svg = readFileSync(sheetSvg, "utf8");
    expect(svg).toContain("native Bonsai Drawing SVG composition");
    expect(svg).not.toMatch(/<image\b/i);
  }
});
