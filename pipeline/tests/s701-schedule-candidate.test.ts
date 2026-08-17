import { expect, test } from "bun:test";
import { createHash } from "node:crypto";
import { mkdtempSync, readFileSync } from "node:fs";
import { tmpdir } from "node:os";
import { join, resolve } from "node:path";

const root = resolve(import.meta.dir, "../..");
const script = resolve(root, "pipeline/scripts/s701_schedule_candidate.py");
const ifcPath = resolve(root, "2504 GBTB Yanlord Zhuhai.ifc");
const currentIfcHash = createHash("sha256").update(readFileSync(ifcPath)).digest("hex");

test("S-701 compiles a current read-only evidence schedule and PNG", () => {
  const temporary = mkdtempSync(join(tmpdir(), "s701-"));
  const report = join(temporary, "s701-report.json");
  const result = Bun.spawnSync(["python3", script, "--root", root, "--input-ifc", ifcPath,
    "--expected-ifc-sha256", currentIfcHash, "--review-csv", join(temporary, "review.csv"),
    "--output-svg", join(temporary, "candidate.svg"), "--proof-png", join(temporary, "proof.png"),
    "--report", report], { cwd: root });
  expect(result.exitCode).toBe(0);
  const payload = JSON.parse(readFileSync(report, "utf8"));
  expect(payload.source_ifc_sha256).toBe(currentIfcHash);
  expect(payload.summary.record_count).toBe(53);
  expect(payload.summary.section_counts).toEqual({
    "家具产品身份": 9, "家电与移动厨电": 18, "门窗与五金": 19, "暖通设备类型": 3, "洁具与排水": 1, "墙面材料系统": 3,
  });
  expect(payload.gates.confirmed_candidate_unresolved_split).toBe(true);
  expect(payload.gates.automatic_ifc_write_allowed).toBe(false);
  expect(payload.gates.construction_release_ready).toBe(false);
  expect(readFileSync(join(temporary, "proof.png")).subarray(0, 8)).toEqual(Buffer.from("89504e470d0a1a0a", "hex"));
}, 20_000);

test("S-701 source has no IFC write path", () => {
  const source = readFileSync(script, "utf8");
  expect(source).not.toContain("model.write(");
  expect(source).not.toContain("ifcopenshell.api.run(");
});
