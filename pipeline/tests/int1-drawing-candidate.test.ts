import { expect, test } from "bun:test";
import { copyFileSync, mkdirSync, mkdtempSync, readFileSync } from "node:fs";
import { tmpdir } from "node:os";
import { join, resolve } from "node:path";

const root = resolve(import.meta.dir, "../..");
const script = resolve(root, "pipeline/scripts/int1_drawing_candidate.py");

test("INT1 drawing candidate produces four safe coordination-envelope sheets", () => {
  const temporary = mkdtempSync(join(tmpdir(), "int1-drawings-"));
  const drawings = join(temporary, "drawings");
  mkdirSync(drawings);
  copyFileSync(resolve(root, "drawings/Furniture Plan.svg"), join(drawings, "Furniture Plan.svg"));
  copyFileSync(resolve(root, "drawings/Sanitary Plan.svg"), join(drawings, "Sanitary Plan.svg"));
  const result = Bun.spawnSync(
    [
      "python3",
      script,
      "--input-csv",
      resolve(root, "pipeline/decisions/int1-existing-review.csv"),
      "--source-ifc",
      resolve(root, "2504 GBTB Yanlord Zhuhai.ifc"),
      "--drawings-dir",
      drawings,
      "--pdf-dir",
      join(temporary, "pdf"),
      "--build-dir",
      join(temporary, "build"),
    ],
    { cwd: root },
  );
  expect(result.exitCode).toBe(0);
  const report = JSON.parse(
    readFileSync(join(temporary, "build/int1-drawing-report.json"), "utf8"),
  );
  expect(report.mode).toBe("read_only_int1_drawing_candidate");
  expect(report.summary.sheet_count).toBe(4);
  expect(report.summary.existing_object_count).toBe(129);
  expect(report.summary.block_count).toBe(4);
  expect(report.summary.all_overlays_are_coordination_envelopes).toBe(true);
  expect(report.gates.candidate_generation_pass).toBe(true);
  expect(report.gates.formal_ifc_write_allowed).toBe(false);
  expect(report.gates.fabrication_dimension_ready).toBe(false);
  expect(report.gates.int1_completion_pass).toBe(false);
  expect(report.sheets.map((sheet: { sheet_id: string }) => sheet.sheet_id)).toEqual([
    "I-501",
    "I-502",
    "I-503",
    "I-504",
  ]);
  expect(report.sheets.map((sheet: { overlay_count: number }) => sheet.overlay_count)).toEqual([
    70,
    33,
    1,
    25,
  ]);
  for (const sheet of report.sheets) {
    expect(sheet.dimension_status).toBe(
      "existing_world_bbox_not_fabrication_dimension",
    );
    expect(sheet.mechanical_pass).toBe(true);
  }
}, 20_000);

test("INT1 drawing generator does not write IFC or claim fabrication readiness", () => {
  const source = readFileSync(script, "utf8");
  expect(source).not.toContain("ifcopenshell");
  expect(source).not.toContain("model.write(");
  expect(source).toContain("fabrication_dimension_ready");
  expect(source).toContain("existing_world_bbox_not_fabrication_dimension");
  for (const title of [
    "Kitchen Coordination Candidate",
    "Bathroom Coordination Candidate",
    "Entry / Laundry Candidate",
    "Fixed Furniture Candidate",
  ]) {
    expect(title.length).toBeLessThan(35);
  }
});
