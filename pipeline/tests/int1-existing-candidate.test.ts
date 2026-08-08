import { expect, test } from "bun:test";
import { mkdtempSync, readFileSync } from "node:fs";
import { tmpdir } from "node:os";
import { join, resolve } from "node:path";

const root = resolve(import.meta.dir, "../..");
const script = resolve(root, "pipeline/scripts/int1_existing_candidate.py");
const ifc = resolve(root, "2504 GBTB Yanlord Zhuhai.ifc");
const roles = resolve(root, "pipeline/decisions/furniture-installation-role.csv");
const products = resolve(root, "pipeline/decisions/furniture-product-register.csv");
const constructionAudit = resolve(root, "build/construction-surfaces/audit.json");

test("INT1 candidate preserves role split and blocks unverified fabrication inputs", () => {
  const temporary = mkdtempSync(join(tmpdir(), "int1-existing-"));
  const output = join(temporary, "build");
  const decisionCsv = join(temporary, "int1-existing-review.csv");
  const result = Bun.spawnSync(
    [
      "python3",
      script,
      "--input",
      ifc,
      "--output-dir",
      output,
      "--decision-csv",
      decisionCsv,
      "--role-decisions",
      roles,
      "--product-register",
      products,
      "--construction-audit",
      constructionAudit,
    ],
    { cwd: root },
  );
  expect(result.exitCode).toBe(0);

  const report = JSON.parse(
    readFileSync(join(output, "int1-existing-report.json"), "utf8"),
  );
  expect(report.mode).toBe("read_only_existing_object_candidate");
  expect(report.summary.furniture_total).toBe(89);
  expect(report.summary.fixed_furniture).toBe(70);
  expect(report.summary.loose_furniture).toBe(19);
  expect(report.summary.furniture_by_sheet).toEqual({
    "I-501": 63,
    "I-502": 1,
    "I-504": 25,
  });
  expect(report.summary.sheet_ids).toEqual(["I-501", "I-502", "I-503", "I-504"]);
  expect(report.summary.block_count).toBe(4);
  expect(report.summary.geberit_flush_plate_semantic_corrections).toBe(2);
  expect(report.gates.candidate_generation_pass).toBe(true);
  expect(report.gates.ifc_write_allowed).toBe(false);
  expect(report.gates.fabrication_dimensions_ready).toBe(false);
  expect(report.gates.int1_completion_pass).toBe(false);

  const scopes = report.blockers.map((item: { scope: string }) => item.scope);
  expect(scopes).toEqual([
    "equipment_installation_drawings",
    "gas_meter_site_survey",
    "sanitary_rough_in_drawings",
    "hardware_motion_envelopes",
  ]);
  expect(report.records.every((item: { dimension_status: string }) =>
    item.dimension_status === "existing_world_bbox_not_fabrication_dimension"
  )).toBe(true);

  const byId = new Map(
    report.records.map((item: { global_id: string }) => [item.global_id, item]),
  );
  expect(byId.get("3JAkt8PsX7vPfGKWLK5EKp")?.sheet_id).toBe("I-503");
  expect(byId.get("3JAkt8PsX7vPfGKWLK5EKp")?.type_description).toBe(
    "FN23BQH W600D660H1655",
  );
  expect(byId.get("0UOnmuAdP1MPy6p3olwiEU")?.type_name).toBe("WD01");
  expect(byId.get("3PQOXKxgj6IftqWXFXMQXG")?.type_name).toBe("OV01");
  expect(byId.get("288GLY62v8kPPydA1lAK8W")?.type_name).toBe("HD01");
  for (const globalId of [
    "2gFgcOYEXEaQWAzcKulTFt",
    "2lDPsdQevFSfeOThtjSlPG",
  ]) {
    expect(byId.get(globalId)?.review_status).toBe(
      "semantic_correction_required_flush_plate_not_wcseat",
    );
  }

  const csv = readFileSync(decisionCsv, "utf8");
  expect((csv.match(/,BLOCK,/g) ?? []).length).toBe(4);
  expect(csv).toContain("existing_world_bbox_not_fabrication_dimension");
}, 20_000);

test("INT1 script contains no formal IFC write path", () => {
  const source = readFileSync(script, "utf8");
  expect(source).not.toContain("model.write(");
  expect(source).not.toContain("ifcopenshell.api.run(");
  expect(source).not.toContain("copyfile(");
});
