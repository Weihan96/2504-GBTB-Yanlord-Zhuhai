import { expect, test } from "bun:test";
import { resolve } from "node:path";

const root = resolve(import.meta.dir, "../..");
const scriptPath = resolve(root, "pipeline/scripts/rcp1_route_blender_review.py");

test("RCP1 route review shows confirmed constraint skeletons with true depth", async () => {
  const source = await Bun.file(scriptPath).text();
  expect(source).toContain('COLLECTION_NAME = "RCP1_ROUTE_CONSTRAINT_REVIEW"');
  expect(source).toContain("A02→H03");
  expect(source).toContain("A03→H04→H02");
  expect(source).toContain("H01 室外机接口位置");
  expect(source).toContain("H07 冷凝水排放接口");
  expect(source).toContain("不是最终风管、冷媒管或厂家接口");
  expect(source).toContain('PAIRING_REVIEW_HIDDEN_PREFIXES = ("RCP1_PIPE_", "RCP1_SERVICE_ARROW_")');
  expect(source).toContain("obj.hide_set(True)");
  expect(source).toContain("obj.show_in_front = False");
  expect(source).toContain("space.shading.show_xray = False");
  expect(source).toContain("space.overlay.show_wireframes = False");
  expect(source).not.toContain('display_type = "WIRE"');
});
