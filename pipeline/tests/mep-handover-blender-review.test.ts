import { expect, test } from "bun:test";
import { readFileSync } from "node:fs";
import { resolve } from "node:path";

const source = readFileSync(resolve(import.meta.dir, "../scripts/mep_handover_blender_review.py"), "utf8");

test("developer handover Blender review is depth-aware and cannot write IFC", () => {
  expect(source).toContain('ROOT_COLLECTION = "DEVELOPER_HANDOVER_MEP_REFERENCE"');
  expect(source).toContain('"01_ELECTRIC_PURPLE"');
  expect(source).toContain('"02_SWITCH_ORANGE"');
  expect(source).toContain('"03_PLUM_BLUE"');
  expect(source).toContain('"04_RELOCATION_REVIEW_RED"');
  expect(source).toContain("obj.show_in_front = False");
  expect(source).toContain('space.shading.type = "SOLID"');
  expect(source).toContain('largest.type = "VIEW_3D"');
  expect(source).toContain('window.workspace = layout_workspace');
  expect(source).toContain('if area.type == "CONSOLE"');
  expect(source).toContain("FURNITURE_YELLOW");
  expect(source).toContain("IFC_CONTEXT_GREY");
  expect(source).not.toContain("model.write(");
  expect(source).not.toContain("bpy.ops.wm.save_as_mainfile");
});
