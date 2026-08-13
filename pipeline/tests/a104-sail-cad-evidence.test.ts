import { expect, test } from "bun:test";
import { resolve } from "node:path";

const root = resolve(import.meta.dir, "../..");
const auditPath = resolve(root, "drawings/evidence/RIMADESIO-Sail-monorotaia-mechanical-audit.json");

test("Sail official CAD audit keeps generic examples out of project nominal dimensions", async () => {
  const audit = await Bun.file(auditPath).json();
  expect(audit.source.official_dwg_sha256).toBe(
    "4eeaeb27b0668bb22b9b4fb191091e07eabbb2b259d90cedcdc45584a1f68e52",
  );
  expect(audit.cad_inventory.entity_count).toBe(1338);
  expect(audit.cad_inventory.generic_dimensions).toEqual({
    panel_width_mm: [1000],
    opening_width_mm: [976, 989, 1978],
    rail_width_mm: [2011, 2037, 4022],
    panel_height_mm: [2670],
    opening_height_mm: [2666, 2678],
  });
  expect(audit.ifc_group.members.sort()).toEqual([
    "0zjVS5FBbBewgUkk0fdfiv",
    "2D5BPoo2XFSvhTdfPenCh7",
  ]);
  for (const occurrence of Object.values(audit.ifc_occurrences) as any[]) {
    expect(occurrence.overall_width_mm).toBeNull();
    expect(occurrence.overall_height_mm).toBeNull();
  }
  expect(audit.comparisons.m05_panel_width.numeric_match).toBe(true);
  expect(audit.comparisons.m05_panel_width.project_nominal_dimension_proven).toBe(false);
  expect(audit.comparisons.m05_panel_height.numeric_match).toBe(false);
  expect(audit.comparisons.m06_rail_long_axis.numeric_match).toBe(false);
  expect(audit.evidence_boundary.formal_ifc_write_allowed).toBe(false);
  expect(audit.gates.mechanical_pass).toBe(true);
});
