import { expect, test } from "bun:test";
import { readFileSync } from "node:fs";
import { resolve } from "node:path";

const root = resolve(import.meta.dir, "../..");
const source = readFileSync(resolve(root, "pipeline/scripts/a103_wall_tag_candidate.py"), "utf8");
const review = readFileSync(resolve(root, "pipeline/scripts/a103_wall_tag_review.py"), "utf8");

test("A-103 wall tags remain a review-only geometry-neutral candidate", () => {
  expect(source).toContain('"existing": "EW01-EW84"');
  expect(source).toContain('"new": "NW01-NW04"');
  expect(source).toContain("north_to_south_then_west_to_east_then_global_id");
  expect(source).toContain('"formal_ifc_write_allowed": False');
  expect(source).toContain("protected_products_geometry_exact");
  expect(source).toContain("EXPECTED_CANDIDATE_IDS_PASS = 587");
  expect(source).toContain("EXPECTED_SOURCE_IDS_PASS = {499, 587}");
  expect(source).not.toContain("FORMAL_IFC.write");
  expect(review).toContain("label_collision_count");
  expect(review).toContain('"formal_ifc_write_allowed": False');
  expect(review).toContain("A-103 墙编号一次性审核候选");
});
