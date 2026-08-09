import { expect, test } from "bun:test";
import { readFileSync } from "node:fs";
import { resolve } from "node:path";

const script = resolve(import.meta.dir, "../scripts/int1_handoff_blender_review.py");

test("INT1 Blender review uses three depth-aware groups without IFC writes", () => {
  const source = readFileSync(script, "utf8");
  expect(source).toContain('ROOT_COLLECTION = "INT1_HANDOFF_REVIEW"');
  expect(source).toContain('"01_NAMED_GREEN"');
  expect(source).toContain('"02_GEOMETRY_ORANGE"');
  expect(source).toContain('"03_LEGACY_RED"');
  expect(source).toContain('obj.show_in_front = False');
  expect(source).toContain('area.spaces.active.shading.type = "SOLID"');
  expect(source).toContain('area.spaces.active.shading.show_xray = False');
  expect(source).toContain('area.spaces.active.overlay.show_wireframes = False');
  expect(source).toContain('FURNITURE_YELLOW');
  expect(source).not.toContain("model.write(");
  expect(source).not.toContain("bpy.ops.wm.save_as_mainfile");
});
