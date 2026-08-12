import { expect, test } from "bun:test";
import { mkdtempSync, readFileSync } from "node:fs";
import { tmpdir } from "node:os";
import { join, resolve } from "node:path";

const root = resolve(import.meta.dir, "../..");
const script = resolve(root, "pipeline/scripts/int1_existing_candidate.py");
const ifc = resolve(root, "2504 GBTB Yanlord Zhuhai.ifc");
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
  expect(report.summary.equipment_ssot_linked_object_records).toBe(123);
  expect(report.summary.equipment_ssot_unlinked_scoped_records).toBe(0);
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

  const kitchenProducts = new Map(
    report.kitchen_product_installation_requirements.map((item: any) => [item.equipment_id, item]),
  );
  expect(kitchenProducts.get("APP-011")?.model).toContain("HB754G2B1W");
  expect(kitchenProducts.get("APP-012")?.model).toContain("ER9EPA33MP/01");
  expect(kitchenProducts.get("APP-013")?.model).toContain("LS33R6VB9W/01");
  for (const equipmentId of ["APP-009", "APP-010"]) {
    const dishwasher = kitchenProducts.get(equipmentId);
    expect(dishwasher?.sheet_id).toBe("I-501");
    expect(dishwasher?.use_location).toBe("西厨岛台");
    expect(dishwasher?.procurement_status).toBe("candidate");
    expect(dishwasher?.final_product_confirmed).toBe(false);
    expect(dishwasher?.project_interface_status).toBe("unlocated");
    expect(dishwasher?.official_source_ids).toEqual([
      "APP-DW-INSTALL-001",
      "APP-DW-SPEC-001",
    ]);
    expect(dishwasher.requirements.every((item: any) =>
      item.status === "confirmed" &&
      ["official_exact_model", "user_input"].includes(item.value_origin)
    )).toBe(true);
  }
  const ovenRequirements = new Map(
    kitchenProducts.get("APP-011").requirements.map((item: any) => [item.parameter_key, item]),
  );
  expect(ovenRequirements.get("rated_current")?.value).toBe("16");
  expect(ovenRequirements.get("niche_width_min")?.source_id).toBe("APP-011-OFFICIAL-001");

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
  expect(csv).not.toContain("\r\n");
  expect((csv.match(/,BLOCK,/g) ?? []).length).toBe(4);
  expect(csv).toContain("existing_world_bbox_not_fabrication_dimension");
}, 20_000);

test("INT1 script contains no formal IFC write path", () => {
  const source = readFileSync(script, "utf8");
  expect(source).not.toContain("model.write(");
  expect(source).not.toContain("ifcopenshell.api.run(");
  expect(source).not.toContain("copyfile(");
});
