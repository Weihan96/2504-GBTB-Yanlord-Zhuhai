import { expect, test } from "bun:test";
import { resolve } from "node:path";

const root = resolve(import.meta.dir, "../..");
const script = resolve(root, "pipeline/scripts/wfin_wall_finish_candidate.py");

test("WFIN generator records the confirmed room rule and stays read-only", async () => {
  const source = await Bun.file(script).text();
  expect(source).toContain('TADELAKT_SPACE_NAMES = {');
  expect(source).toContain('"次卧"');
  expect(source).toContain('"主卫干区"');
  expect(source).toContain('"客卫"');
  expect(source).toContain('"automatic_ifc_write_allowed": False');
  expect(source).not.toContain("model.write(");
});

test("WFIN candidate segments CLADDING at mixed finish boundaries", async () => {
  const temp = resolve(root, "build/test-wfin");
  const result = Bun.spawnSync([
    "python3", script,
    "--input", resolve(root, "2504 GBTB Yanlord Zhuhai.ifc"),
    "--source-svg", resolve(root, "drawings/Wall Finish Plan.svg"),
    "--output-svg", resolve(temp, "candidate.svg"),
    "--register", resolve(temp, "register.csv"),
    "--segment-register", resolve(temp, "segments.csv"),
    "--report", resolve(temp, "report.json"),
    "--containment-tolerance-mm", "0.5",
  ], { cwd: root });
  expect(result.exitCode).toBe(0);
  const report = await Bun.file(resolve(temp, "report.json")).json();
  expect(report.cladding_count).toBe(51);
  expect(report.single_finish_object_count).toBe(48);
  expect(report.mixed_finish_object_count).toBe(3);
  expect(report.finish_segment_count).toBe(54);
  expect(report.tadelakt_segment_count).toBe(26);
  expect(report.white_wall_segment_count).toBe(28);
  expect(report.main_bath_dry_segment_count).toBe(3);
  expect(report.mixed_finish_object_ids.sort()).toEqual([
    "06wFwLoDD6ie5iCTnc_yad",
    "0e0XOb$L18ZBVYJiQJrQ1p",
    "3bMoS7bIT8wBmjqRvdPunu",
  ]);
  expect(report.source_svg_object_count + report.fallback_object_count).toBe(51);
  expect(report.qa.formal_ifc_unchanged).toBe(true);
  const segments = await Bun.file(resolve(temp, "segments.csv")).text();
  expect((segments.match(/,TADELAKT,/g) ?? []).length).toBe(26);
  expect((segments.match(/,WHITE_WALL,/g) ?? []).length).toBe(28);
});
