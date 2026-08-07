import { expect, test } from "bun:test";
import { resolve } from "node:path";

const root = resolve(import.meta.dir, "../..");
const scriptPath = resolve(root, "pipeline/scripts/p0_blender_review.py");

test("combined Blender review preserves true depth and clear category colors", async () => {
  const source = await Bun.file(scriptPath).text();
  expect(source).toContain('obj.show_in_front = False');
  expect(source).toContain('space.shading.show_xray = False');
  expect(source).toContain('space.shading.type = "SOLID"');
  expect(source).toContain('"a104_review"');
  expect(source).toContain('"wet_floor"');
  expect(source).toContain('"material_pending"');
  expect(source).toContain('bpy.ops.object.select_all(action="DESELECT")');
  expect(source).toContain('bpy.ops.object.mode_set(mode="OBJECT")');
  expect(source).toContain('obj.display_type = "SOLID"');
  expect(source).toContain('add_bbox_outline(collection, obj)');
  expect(source).not.toContain('display_type = "WIRE"');
});
