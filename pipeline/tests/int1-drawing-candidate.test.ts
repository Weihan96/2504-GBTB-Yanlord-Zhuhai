import { expect, test } from "bun:test";
import { copyFileSync, mkdirSync, mkdtempSync, readFileSync, writeFileSync } from "node:fs";
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
  const register = readFileSync(
    resolve(root, "pipeline/decisions/int1-existing-review.csv"),
    "utf8",
  );
  const sourceHash = register.trim().split("\n")[1].split(",").at(-1);
  const existingReport = join(temporary, "int1-existing-report.json");
  const requirements = [
    ["niche_height_min", "780", "mm"],
    ["niche_height_max", "835", "mm"],
    ["niche_width_min", "600", "mm"],
    ["niche_width_max", "608", "mm"],
    ["niche_depth_min", "550", "mm"],
    ["water_connection", "G3/4 cold water", ""],
    ["drain_connection_od", "38", "mm"],
  ].map(([parameter_key, value, unit], index) => ({
    requirement_id: `REQ-T${index}`,
    discipline: "MULTI",
    parameter_key,
    value,
    unit,
    value_origin: "official_exact_model",
    status: "confirmed",
    source_id: "APP-DW-INSTALL-001",
    blocks_release: false,
  }));
  writeFileSync(existingReport, JSON.stringify({
    source: { ifc_sha256: sourceHash },
    kitchen_product_installation_requirements: ["APP-009", "APP-010"].map(
      (equipment_id) => ({
        equipment_id,
        sheet_id: "I-501",
        item_name: "岛台洗碗机",
        manufacturer: "Siemens",
        model: "Siemens SJ85ZX26MC",
        procurement_status: "candidate",
        decision_status: "partial",
        use_location: "西厨岛台",
        official_source_ids: ["APP-DW-INSTALL-001", "APP-DW-SPEC-001"],
        requirements,
        product_constraints_available: true,
        project_interface_status: "unlocated",
        final_product_confirmed: false,
        release_blockers: ["project interface unlocated", "final procurement unconfirmed"],
      }),
    ),
  }));
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
      "--existing-report",
      existingReport,
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
  expect(report.summary.installation_interface_row_count).toBe(2);
  expect(report.summary.unlocated_interface_count).toBe(2);
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
  expect(report.sheets.map((sheet: { installation_interface_row_count: number }) =>
    sheet.installation_interface_row_count
  )).toEqual([2, 0, 0, 0]);
  const kitchenSvg = readFileSync(join(drawings, "I-501-kitchen-existing-candidate.svg"), "utf8");
  expect(kitchenSvg).toContain('data-equipment-id="APP-009"');
  expect(kitchenSvg).toContain('data-equipment-id="APP-010"');
  expect(kitchenSvg).toContain('data-project-interface-status="unlocated"');
  expect(kitchenSvg).toContain("CANDIDATE / NOT PURCHASED");
  expect(kitchenSvg).toContain("G3/4 cold water · drain Ø38 mm");
  expect(kitchenSvg).toContain("Niche H780–835 mm · W600–608 mm · D≥550 mm");
  expect(kitchenSvg).toContain("rough-in XYZ / valves / hose path / opening position");
  expect((kitchenSvg.match(/data-equipment-id=/g) ?? []).length).toBe(2);
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
