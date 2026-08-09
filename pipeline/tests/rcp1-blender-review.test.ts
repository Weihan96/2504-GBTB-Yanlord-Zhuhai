import { expect, test } from "bun:test";
import { resolve } from "node:path";

const root = resolve(import.meta.dir, "../..");
const scriptPath = resolve(root, "pipeline/scripts/rcp1_blender_review.py");

test("RCP1 Blender review separates demolition references from built walls", async () => {
  const source = await Bun.file(scriptPath).text();
  expect(source).toContain('DEMOLITION_COLLECTION = "A102_DEMOLITION_REFERENCE"');
  expect(source).toContain('DEMOLITION_STATUSES = {"DEMOLISH", "DEMOLISHED"}');
  expect(source).toContain("source.hide_set(True)");
  expect(source).toContain("collection.hide_viewport = True");
  expect(source).toContain('reference.display_type = "SOLID"');
  expect(source).toContain("reference.show_in_front = False");
  expect(source).toContain('space.shading.type = "SOLID"');
  expect(source).toContain("space.shading.show_xray = False");
  expect(source).toContain("space.overlay.show_wireframes = False");
  expect(source).toContain("RCP1_OPTION_A");
  expect(source).toContain("RCP1_OPTION_B");
  expect(source).toContain("RCP1_A06_固定机位");
  expect(source).toContain("R20 客厅：送回风待深化");
  expect(source).not.toContain('display_type = "WIRE"');
});
