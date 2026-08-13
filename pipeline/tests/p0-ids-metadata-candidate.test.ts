import { expect, test } from "bun:test";
import { readFileSync } from "node:fs";
import { resolve } from "node:path";

const root = resolve(import.meta.dir, "../..");
const script = resolve(root, "pipeline/scripts/p0_ids_metadata_candidate.py");

test("P0 metadata candidate is limited to geometry-neutral, evidenced fields", () => {
  const source = readFileSync(script, "utf8");
  expect(source).toContain('properties={"Length": record["length_mm"]}');
  expect(source).toContain('PredefinedType = "DOOR"');
  expect(source).toContain('PredefinedType = "WINDOW"');
  expect(source).toContain('"wall Tag"');
  expect(source).toContain('"M05-M07 OverallWidth/OverallHeight"');
    expect(source).not.toMatch(/\.OperationType\s*=(?!=)/);
  expect(source).not.toContain("OverallWidth =");
  expect(source).not.toContain("OverallHeight =");
});

test("P0 metadata candidate is idempotent after the wall-tag IDS batch", () => {
  const source = readFileSync(script, "utf8");
  expect(source).toContain("EXPECTED_SOURCE_IDS_PASS = {392, 499, 587}");
  expect(source).toContain("MINIMUM_CANDIDATE_IDS_PASS = 499");
  expect(source).toContain("max(");
  expect(source).toContain('"protected_products_geometry_exact"');
  expect(source).toContain('"fills_voids_relationships_unchanged"');
});
