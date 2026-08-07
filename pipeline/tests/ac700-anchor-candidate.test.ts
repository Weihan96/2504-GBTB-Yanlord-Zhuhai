import { expect, test } from "bun:test";
import { resolve } from "node:path";

const root = resolve(import.meta.dir, "../..");

test("AC700 candidate scope and type evidence are exact", () => {
  const source = `
import json,sys
sys.path.insert(0,"pipeline/scripts")
from ac700_anchor_candidate import AC700_GLOBAL_ID,AC700_TYPE_NAME,AC700_DESCRIPTION
print(json.dumps([AC700_GLOBAL_ID,AC700_TYPE_NAME,AC700_DESCRIPTION]))
`;
  const result = Bun.spawnSync(["python3", "-c", source], { cwd: root });
  expect(result.exitCode).toBe(0);
  expect(JSON.parse(result.stdout.toString())).toEqual([
    "1PUCikoaP5fgiYt8sJd8$6",
    "AC700",
    "RPIZ-22FSLN5QD/P 700x447x192",
  ]);
});

test("AC700 face evidence rejects unsupported semantic sides", () => {
  const source = `
import sys
sys.path.insert(0,"pipeline/scripts")
from ac700_anchor_candidate import face_anchor_evidence
try:
  face_anchor_evidence(None,"centre",0.01)
except ValueError as error:
  print(str(error))
`;
  const result = Bun.spawnSync(["python3", "-c", source], { cwd: root });
  expect(result.exitCode).toBe(0);
  expect(result.stdout.toString()).toContain("unsupported AC700 face side");
});
