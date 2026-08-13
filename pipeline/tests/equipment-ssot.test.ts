import { expect, test } from "bun:test";
import { readFileSync } from "node:fs";
import { resolve } from "node:path";

const root = resolve(import.meta.dir, "../..");
const script = resolve(root, "pipeline/scripts/equipment_ssot.py");

test("equipment SSOT validates projections and covers every scoped IFC object", () => {
  const run = Bun.spawnSync(["python3", script, "all"], { cwd: root });
  expect(run.exitCode).toBe(0);
  const report = JSON.parse(run.stdout.toString());
  expect(report.validate).toMatchObject({
    master_count: 153,
    requirement_count: 578,
    source_count: 118,
    schema_version: "1.0.0",
  });
  expect(report.projections).toMatchObject({
    appliance_rows: 18,
    furniture_rows: 9,
    elec_evidence_rows: 24,
    hvac_evidence_rows: 4,
    furniture_role_rows: 51,
    mode: "check",
  });
  expect(report.ifc_coverage).toMatchObject({
    scope_count: 162,
    covered_count: 162,
    missing_count: 0,
    duplicate_count: 0,
  });
  expect(report.ifc_coverage.class_counts).toEqual({
    IfcElectricAppliance: 19,
    IfcFurniture: 89,
    IfcSanitaryTerminal: 27,
    IfcWasteTerminal: 3,
    IfcSensor: 1,
    IfcDoor: 8,
    IfcWindow: 11,
    IfcElementAssembly: 3,
    IfcBuildingElementProxy: 1,
  });
}, 180_000);

test("equipment SSOT never writes the formal IFC", () => {
  const source = readFileSync(script, "utf8");
  expect(source).not.toContain("model.write(");
  expect(source).not.toContain("ifc.write(");
});

test("equipment SSOT can audit an IFC evidence hash refresh without writing", () => {
  const sourcePath = resolve(root, "pipeline/decisions/source-evidence-register.csv");
  const before = readFileSync(sourcePath);
  const run = Bun.spawnSync(["python3", script, "refresh-ifc-source-hashes", "--dry-run"], { cwd: root });
  expect(run.exitCode, run.stderr.toString()).toBe(0);
  const report = JSON.parse(run.stdout.toString());
  expect(report.scoped_ifc_object_count).toBe(162);
  expect(report.target_source_count).toBeGreaterThan(1);
  expect(report.dry_run).toBe(true);
  expect(readFileSync(sourcePath)).toEqual(before);
}, 30_000);

test("Geberit receipt shortfalls and installation unknowns block release mechanically", () => {
  const requirementsPath = resolve(root, "pipeline/decisions/equipment-installation-requirements.csv");
  const extract = String.raw`
import csv,json,sys
targets={"SAN-003","SAN-014","SAN-001","DRAIN-GEB-001","DRAIN-GEB-002","DRAIN-GEB-004"}
with open(sys.argv[1],encoding="utf-8-sig",newline="") as stream:
    rows=[row for row in csv.DictReader(stream) if row["equipment_id"] in targets]
print(json.dumps(rows,ensure_ascii=False))
`;
  const run = Bun.spawnSync(["python3", "-c", extract, requirementsPath], { cwd: root });
  expect(run.exitCode, run.stderr.toString()).toBe(0);
  const rows = JSON.parse(run.stdout.toString()) as Array<Record<string, string>>;
  const get = (equipmentId: string, key: string) =>
    rows.find((row) => row.equipment_id === equipmentId && row.parameter_key === key);

  for (const equipmentId of ["SAN-003", "SAN-014"]) {
    expect(get(equipmentId, "required_quantity")?.value_number).toBe("2");
    expect(get(equipmentId, "received_quantity")?.value_number).toBe("1");
    expect(get(equipmentId, "purchase_shortfall")).toMatchObject({
      value_number: "1",
      status: "pending",
      blocks_release: "yes",
    });
  }

  const blockerKeys: Record<string, string[]> = {
    "SAN-014": ["completed_wc_height", "wall_anchor_design", "boxing_finished_thickness", "nuna_hidden_interface_coordination"],
    "SAN-003": ["panel_centre_height", "finished_surface_thickness", "service_access_detail"],
    "SAN-001": ["actual_cut_length", "placement_mode", "drain_centre_coordinate", "tile_layout_coordination"],
    "DRAIN-GEB-002": ["final_installation_level", "waterproofing_system", "site_pipe_orientation"],
    "DRAIN-GEB-001": ["final_assembly_orientation"],
    "DRAIN-GEB-004": ["washing_machine_hose_connection", "installation_height", "floor_slope", "waterproofing_detail"],
  };
  for (const [equipmentId, keys] of Object.entries(blockerKeys)) {
    for (const key of keys) {
      expect(get(equipmentId, key), `${equipmentId}.${key}`).toMatchObject({
        status: "pending",
        blocks_release: "yes",
      });
    }
  }
});

test("custom sinks and the guest-bath mirror candidate keep fabrication unknowns explicit", () => {
  const equipmentPath = resolve(root, "pipeline/decisions/equipment-register.csv");
  const requirementsPath = resolve(root, "pipeline/decisions/equipment-installation-requirements.csv");
  const extract = String.raw`
import csv,json,sys
equipment_ids={"SAN-KIT-001","SAN-022","SAN-021","SAN-018","FURIFC-003","FUR-BATHG-MIRROR-001"}
def rows(path):
    with open(path,encoding="utf-8-sig",newline="") as stream:
        return list(csv.DictReader(stream))
print(json.dumps({
  "masters":[row for row in rows(sys.argv[1]) if row["equipment_id"] in equipment_ids],
  "requirements":[row for row in rows(sys.argv[2]) if row["equipment_id"] in equipment_ids],
},ensure_ascii=False))
`;
  const run = Bun.spawnSync(["python3", "-c", extract, equipmentPath, requirementsPath], { cwd: root });
  expect(run.exitCode, run.stderr.toString()).toBe(0);
  const data = JSON.parse(run.stdout.toString()) as {
    masters: Array<Record<string, string>>;
    requirements: Array<Record<string, string>>;
  };
  const master = (equipmentId: string) => data.masters.find((row) => row.equipment_id === equipmentId);
  const requirement = (equipmentId: string, key: string) =>
    data.requirements.find((row) => row.equipment_id === equipmentId && row.parameter_key === key);

  expect(master("SAN-KIT-001")).toMatchObject({
    model: "1014 850",
    selector_kind: "logical_input",
    decision_status: "confirmed",
  });
  expect(requirement("SAN-KIT-001", "product_width")?.value_number).toBe("838");
  expect(requirement("SAN-KIT-001", "cutout_width")?.value_number).toBe("790");
  expect(requirement("SAN-KIT-001", "formal_ifc_placement")).toMatchObject({ status: "pending", blocks_release: "yes" });

  for (const [equipmentId, key] of [
    ["SAN-022", "support_cutout_and_shop_drawing"],
    ["FURIFC-003", "fabrication_shop_drawing"],
    ["SAN-021", "exact_custom_article_code"],
    ["SAN-018", "exact_custom_configuration"],
    ["FUR-BATHG-MIRROR-001", "finished_width"],
    ["FUR-BATHG-MIRROR-001", "lighting_integration_decision"],
  ]) {
    expect(requirement(equipmentId, key), `${equipmentId}.${key}`).toMatchObject({
      status: "pending",
      blocks_release: "yes",
    });
  }
  expect(master("FUR-BATHG-MIRROR-001")).toMatchObject({
    procurement_status: "candidate",
    decision_status: "candidate",
  });
});

test("Rimadesio Sail CAD proves the single-track family without closing project dimensions", () => {
  const equipmentPath = resolve(root, "pipeline/decisions/equipment-register.csv");
  const requirementsPath = resolve(root, "pipeline/decisions/equipment-installation-requirements.csv");
  const sourcesPath = resolve(root, "pipeline/decisions/source-evidence-register.csv");
  const extract = String.raw`
import csv,hashlib,json,pathlib,sys
def rows(path):
    with open(path,encoding="utf-8-sig",newline="") as stream:
        return list(csv.DictReader(stream))
equipment=rows(sys.argv[1]); requirements=rows(sys.argv[2]); sources=rows(sys.argv[3])
source=next(row for row in sources if row["source_id"]=="RIMADESIO-SAIL-MONOROTAIA-001")
cad=pathlib.Path(sys.argv[4])
print(json.dumps({
  "masters":[row for row in equipment if row["equipment_id"] in {"DW-M05","DW-M06"}],
  "requirements":[row for row in requirements if row["equipment_id"] in {"DW-M05","DW-M06"}],
  "source":source,
  "actual_sha256":hashlib.sha256(cad.read_bytes()).hexdigest(),
},ensure_ascii=False))
`;
  const cadPath = resolve(root, "drawings/evidence/RIMADESIO-official-Sail-monorotaia.dwg");
  const run = Bun.spawnSync(["python3", "-c", extract, equipmentPath, requirementsPath, sourcesPath, cadPath], { cwd: root });
  expect(run.exitCode, run.stderr.toString()).toBe(0);
  const data = JSON.parse(run.stdout.toString()) as {
    masters: Array<Record<string, string>>;
    requirements: Array<Record<string, string>>;
    source: Record<string, string>;
    actual_sha256: string;
  };
  expect(data.source).toMatchObject({
    sha256: "4eeaeb27b0668bb22b9b4fb191091e07eabbb2b259d90cedcdc45584a1f68e52",
    status: "verified_official_family_cad",
    formal_ifc_write_allowed: "no",
  });
  expect(data.actual_sha256).toBe(data.source.sha256);
  for (const equipmentId of ["DW-M05", "DW-M06"]) {
    expect(data.masters.find((row) => row.equipment_id === equipmentId)?.source_ids).toContain("RIMADESIO-SAIL-MONOROTAIA-001");
    expect(data.requirements.find((row) => row.equipment_id === equipmentId && row.parameter_key === "operation_type")).toMatchObject({
      status: "candidate",
      source_id: "RIMADESIO-SAIL-MONOROTAIA-001",
      blocks_release: "yes",
    });
    for (const key of ["nominal_width", "nominal_height"]) {
      expect(data.requirements.find((row) => row.equipment_id === equipmentId && row.parameter_key === key)).toMatchObject({
        status: "pending",
        blocks_release: "yes",
      });
    }
  }
});
