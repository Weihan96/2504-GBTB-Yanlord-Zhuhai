import { expect, test } from "bun:test";
import { resolve } from "node:path";

const root = resolve(import.meta.dir, "../..");
const script = resolve(root, "pipeline/scripts/rcp1_legacy_base_audit.py");

test("RCP1 legacy audit is read-only and covers the exact legacy base", async () => {
  const source = await Bun.file(script).text();
  for (const globalId of [
    "0Ik2RcgGbFOhdYTJPgh5AQ",
    "0f2ZLauDH8lRnYj6oervDm",
    "0hHnbLj0X4jPDz4o3QJo1l",
    "10Wm8ivdX7dAVfz4cV8l5Q",
    "1hZRB0eOX8OA8rcjke67P0",
    "1QBdVekDnBsOleyo9PM6rT",
    "1yW7DASIz8qA$2j8z9tdl2",
  ]) {
    expect(source).toContain(globalId);
  }
  expect(source).toContain("BL|RPIZ-22FSN6QD_curve_.001");
  expect(source).toContain('"formal_ifc_write_allowed": False');
  expect(source).toContain('"legacy_blend_write_allowed": False');
  expect(source).toContain('"final_remodel_hvac_design_inferred": False');
  expect(source).not.toContain("model.write(");
  expect(source).not.toContain("bpy.ops.wm.save");
  expect(source).not.toContain("bpy.ops.wm.save_as_mainfile");
});
