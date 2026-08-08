import { expect, test } from "bun:test";

const source = await Bun.file("pipeline/scripts/p0_semantics_batch_postwrite_audit.py").text();

test("P0 semantics postwrite audit compares the formal IFC to the approved candidate", () => {
  expect(source).toContain('classes=("IfcDoor", "IfcWindow", "IfcOpeningElement", "IfcCovering", "IfcSlab")');
  expect(source).toContain('"root_ids_equal"');
  expect(source).toContain('"fills_relationship_ids_equal"');
  expect(source).toContain('"maximum_placement_matrix_delta"');
  expect(source).toContain('"geometry_over_tolerance"');
  expect(source).not.toContain("apply_semantics(");
});
