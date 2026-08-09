import { expect, test } from "bun:test";
import { readFileSync } from "node:fs";
import { resolve } from "node:path";

const source = readFileSync(resolve(import.meta.dir, "../scripts/mep_positioning_blender_review.py"), "utf8");

test("MEP Blender review uses clear depth-aware discipline groups", () => {
  expect(source).toContain('ROOT_COLLECTION = "MEP_POSITIONING_REVIEW"');
  expect(source).toContain('"01_PLUM_BLUE"');
  expect(source).toContain('"02_ELEC_PURPLE"');
  expect(source).toContain('"03_REVIEW_RED"');
  expect(source).toContain('obj.show_in_front = False');
  expect(source).toContain('area.spaces.active.shading.type = "SOLID"');
  expect(source).toContain('area.spaces.active.shading.show_xray = False');
  expect(source).toContain('area.spaces.active.overlay.show_wireframes = False');
  expect(source).toContain("FURNITURE_YELLOW");
  expect(source).not.toContain("model.write(");
  expect(source).not.toContain("bpy.ops.wm.save_as_mainfile");
});
