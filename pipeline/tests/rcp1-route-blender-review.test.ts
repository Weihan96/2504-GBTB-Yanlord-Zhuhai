import { expect, test } from "bun:test";
import { resolve } from "node:path";

const root = resolve(import.meta.dir, "../..");
const scriptPath = resolve(root, "pipeline/scripts/rcp1_route_blender_review.py");

test("RCP1 route review shows pairings without presenting them as routes", async () => {
  const source = await Bun.file(scriptPath).text();
  expect(source).toContain('COLLECTION_NAME = "RCP1_ROUTE_READINESS_REVIEW"');
  expect(source).toContain("不是风管/冷媒管");
  expect(source).toContain("H01 未配对 / 外部接口链待确认");
  expect(source).toContain("H06→H01？仅穿墙链候选");
  expect(source).toContain('PAIRING_REVIEW_HIDDEN_PREFIXES = ("RCP1_PIPE_", "RCP1_SERVICE_ARROW_")');
  expect(source).toContain("obj.hide_set(True)");
  expect(source).toContain("obj.show_in_front = False");
  expect(source).toContain("space.shading.show_xray = False");
  expect(source).toContain("space.overlay.show_wireframes = False");
  expect(source).not.toContain('display_type = "WIRE"');
});
