import { expect, test } from "bun:test";
import { spawnSync } from "node:child_process";

test("PVC110 controlled exception requires exact byte identity for zero change", () => {
  const result = spawnSync(
    "python3",
    [
      "-c",
      `import json,sys; sys.path.insert(0,"pipeline/scripts"); from pvc110_controlled_exception_audit import byte_identity_world_change_mm; print(json.dumps({"same":byte_identity_world_change_mm("abc","abc"),"different":byte_identity_world_change_mm("abc","def")}))`,
    ],
    { cwd: process.cwd(), encoding: "utf8" },
  );
  expect(result.status).toBe(0);
  expect(JSON.parse(result.stdout)).toEqual({ same: 0, different: null });
});

test("PVC110 controlled exception is exact and never authorizes formal write", async () => {
  const source = await Bun.file(
    "pipeline/scripts/pvc110_controlled_exception_audit.py",
  ).text();
  expect(source).toContain('DECISION_ID = "COORD-FLOW-C003-PVC110-BUNDLE"');
  expect(source).toContain('"formal_write_allowed": False');
  expect(source).toContain('"six_read_only_branches_recognized"');
  expect(source).toContain("targeted world-mesh comparison in C003 postwrite audit");
  const result = spawnSync(
    "python3",
    [
      "-c",
      `import json,sys; sys.path.insert(0,"pipeline/scripts"); from pvc110_controlled_exception_audit import TARGET_IDS; print(json.dumps(TARGET_IDS))`,
    ],
    { cwd: process.cwd(), encoding: "utf8" },
  );
  expect(result.status).toBe(0);
  expect(JSON.parse(result.stdout)).toEqual([
    "178mqyyzzFowLcbXcH6prO",
    "0bfVg4Ys1CevZs$qxhkXTo",
  ]);
});
