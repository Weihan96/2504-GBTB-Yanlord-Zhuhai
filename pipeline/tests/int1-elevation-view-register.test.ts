import { expect, test } from "bun:test";
import { mkdtempSync, readFileSync, writeFileSync } from "node:fs";
import { tmpdir } from "node:os";
import { join, resolve } from "node:path";

const root = resolve(import.meta.dir, "../..");
const script = join(root, "pipeline/scripts/int1_elevation_view_register.py");
const planDwg = join(root, "../图纸/矩阵纵横/01 成品房(D1户型-115)平面系统图.dwg");
const planDxf = join(root, "tmp/dwg/d1-handover-plan.dxf");
const elevationDwg = join(root, "../图纸/矩阵纵横/02 成品房(D1户型-115)立面图.dwg");

const requiredColumns = [
  "view_id", "sheet_id", "official_title", "anchor_id", "direction",
  "source_handle", "source_viewport_handle", "source_scale",
  "source_anchor_handle", "source_x_mm", "source_y_mm", "ifc_x_mm",
  "ifc_y_mm", "space_reference", "space_global_id", "review_status",
  "source_locator",
];

function run(output: string, sources: string[] = [planDwg, planDxf, elevationDwg]) {
  return Bun.spawnSync([
    "python3", script,
    "--plan-dwg", sources[0],
    "--plan-dxf", sources[1],
    "--elevation-dwg", sources[2],
    "--output", output,
  ], { cwd: root, stdout: "pipe", stderr: "pipe" });
}

function parseCsv(path: string): Record<string, string>[] {
  const lines = readFileSync(path, "utf8").trimEnd().split("\n");
  const headers = lines[0].split(",");
  return lines.slice(1).map((line) => Object.fromEntries(
    line.split(",").map((value, index) => [headers[index], value]),
  ));
}

test("compiles the 36 official elevation views in deterministic view order", () => {
  const temp = mkdtempSync(join(tmpdir(), "int1-elevation-"));
  const output = join(temp, "register.csv");
  const result = run(output);
  expect(result.exitCode).toBe(0);

  const report = JSON.parse(result.stdout.toString());
  const rows = parseCsv(output);
  const headers = readFileSync(output, "utf8").split("\n", 1)[0].split(",");
  expect(requiredColumns.every((column) => headers.includes(column))).toBe(true);
  expect(rows).toHaveLength(36);
  expect(rows.map((row) => row.view_id)).toEqual(
    Array.from({ length: 36 }, (_, index) => String(index + 1).padStart(2, "0")),
  );
  expect(new Set(rows.map((row) => row.sheet_id))).toEqual(
    new Set(Array.from({ length: 9 }, (_, index) => `EL-${String(index + 1).padStart(2, "0")}`)),
  );
  expect(new Set(rows.map((row) => row.direction))).toEqual(new Set(["+Y", "+X", "-Y", "-X"]));
  expect(rows.every((row) => row.source_viewport_handle === "243CB4")).toBe(true);
  expect(rows.every((row) => row.source_scale === "1:50@A2")).toBe(true);
  expect(rows.every((row) => row.anchor_method === "resolved_leader_or_arrow_target_not_marker_circle")).toBe(true);
  expect(report.view_count).toBe(36);
  expect(report.sheet_count).toBe(9);
  expect(report.anchor_policy).toContain("marker circle insertion is not an anchor");
  expect(report.unindexed_scopes).toEqual([
    "I-503/玄关：01 DWG立面索引未发现官方视图；不得伪称来自DWG",
  ]);
});

test("locks source hashes and records exact duplicate removal for views 02 through 05", () => {
  const output = join(mkdtempSync(join(tmpdir(), "int1-elevation-hash-")), "register.csv");
  const result = run(output);
  expect(result.exitCode).toBe(0);
  const report = JSON.parse(result.stdout.toString());
  expect(report.hashes).toEqual({
    plan_dwg: "ba355f6a90732ad07f843d59e8bab5e1da9daffe7aed74889d1d265b2ce22d7e",
    plan_dxf: "706483de83d90e526d7d7cf4f0b097902a4f97868c16aba5f6fac9d9e942f0f0",
    elevation_dwg: "80c0e01a3e713941f729c3b01750181c18ea4d0ff4a5ae8d65800686f19646dc",
  });
  expect(report.duplicates_removed).toEqual({
    "02": ["24415B"],
    "03": ["24416A"],
    "04": ["244179"],
    "05": ["244188"],
  });
  const rows = parseCsv(output);
  expect(rows.find((row) => row.view_id === "02")?.source_handle).toBe("243C78");
  expect(rows.find((row) => row.view_id === "05")?.duplicate_source_handles).toBe("244188");
});

test("preserves the resolved CAD and IFC anchors with current Space references", () => {
  const output = join(mkdtempSync(join(tmpdir(), "int1-elevation-anchor-")), "register.csv");
  const result = run(output);
  expect(result.exitCode).toBe(0);
  const byId = Object.fromEntries(parseCsv(output).map((row) => [row.view_id, row]));
  expect(byId["01"]).toMatchObject({
    anchor_id: "A1", source_anchor_handle: "243CB7",
    source_x_mm: "126967.662", source_y_mm: "-49873.101",
    ifc_x_mm: "-2285.579", ifc_y_mm: "127.511",
    space_reference: "R14 次卧", space_global_id: "0WyQ2Z9pX5qgwfdTOZwAkw",
  });
  expect(byId["02"]).toMatchObject({
    anchor_id: "A2", source_anchor_handle: "243C70",
    ifc_x_mm: "1820.962", ifc_y_mm: "626.980",
    space_reference: "R20 客厅", space_global_id: "2wgBPVUpv2DvcZCfbe6fdv",
  });
  expect(byId["06"]).toMatchObject({
    anchor_id: "A3", source_anchor_handle: "243CC6",
    ifc_x_mm: "1820.962", ifc_y_mm: "-4240.404",
    space_reference: "R04 中厨", space_global_id: "2fhEbDfK1EkhJwlPikNm$b",
  });
  expect(byId["36"]).toMatchObject({
    anchor_id: "A12", source_anchor_handle: "25F70A",
    ifc_x_mm: "5023.876", ifc_y_mm: "1722.710",
    space_reference: "R22 书房", space_global_id: "0XmeOOraz9tP2_CetH7YYh",
  });
});

test("rejects a source file whose locked hash changes before writing a register", () => {
  const temp = mkdtempSync(join(tmpdir(), "int1-elevation-tamper-"));
  const changedPlan = join(temp, "changed-plan.dwg");
  writeFileSync(changedPlan, Buffer.concat([readFileSync(planDwg), Buffer.from("changed")]));
  const output = join(temp, "should-not-exist.csv");
  const result = run(output, [changedPlan, planDxf, elevationDwg]);
  expect(result.exitCode).toBe(2);
  expect(result.stderr.toString()).toContain("01 plan DWG SHA-256 mismatch");
  expect(() => readFileSync(output)).toThrow();
});
