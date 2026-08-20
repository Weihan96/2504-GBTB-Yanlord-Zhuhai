import { expect, test } from "bun:test";
import { createHash } from "node:crypto";
import { readFileSync } from "node:fs";
import { resolve } from "node:path";

const root = resolve(import.meta.dir, "../..");
const script = resolve(root, "pipeline/scripts/equipment_ssot.py");

function csvRecordCount(path: string): number {
  const program = [
    "import csv,sys",
    "with open(sys.argv[1],encoding='utf-8-sig',newline='') as stream:",
    " print(sum(1 for _ in csv.DictReader(stream)))",
  ].join("\n");
  const run = Bun.spawnSync(["python3", "-c", program, path], { cwd: root });
  expect(run.exitCode, run.stderr.toString()).toBe(0);
  return Number(run.stdout.toString().trim());
}

function fileSha256(path: string): string {
  return createHash("sha256").update(readFileSync(path)).digest("hex");
}

test("equipment SSOT validates projections and covers every scoped IFC object", () => {
  const run = Bun.spawnSync(["python3", script, "all"], { cwd: root });
  expect(run.exitCode).toBe(0);
  const report = JSON.parse(run.stdout.toString());
  expect(report.validate).toMatchObject({
    master_count: csvRecordCount(resolve(root, "pipeline/decisions/equipment-register.csv")),
    requirement_count: csvRecordCount(resolve(root, "pipeline/decisions/equipment-installation-requirements.csv")),
    source_count: csvRecordCount(resolve(root, "pipeline/decisions/source-evidence-register.csv")),
    schema_version: "1.1.0",
  });
  expect(report.projections).toMatchObject({
    appliance_rows: 19,
    furniture_rows: 9,
    elec_evidence_rows: 69,
    hvac_evidence_rows: 4,
    furniture_role_rows: 51,
    mode: "check",
  });
  expect(report.ifc_coverage).toMatchObject({
    scope_count: 165,
    covered_count: 165,
    missing_count: 0,
    duplicate_count: 0,
  });
  expect(report.ifc_coverage.class_counts).toEqual({
    IfcElectricAppliance: 21,
    IfcUnitaryEquipment: 1,
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

test("local evidence uses an exact lowercase SHA-256 instead of a truncated digest", () => {
  const extract = String.raw`
import csv,json,sys
with open(sys.argv[1],encoding="utf-8-sig",newline="") as stream:
    print(json.dumps(list(csv.DictReader(stream)),ensure_ascii=False))
`;
  const register = resolve(root, "pipeline/decisions/source-evidence-register.csv");
  const run = Bun.spawnSync(["python3", "-c", extract, register], { cwd: root });
  expect(run.exitCode, run.stderr.toString()).toBe(0);
  const rows = JSON.parse(run.stdout.toString()) as Array<Record<string, string>>;
  for (const row of rows) {
    if (row.sha256.startsWith("not_applicable_")) continue;
    expect(row.sha256).toMatch(/^[0-9a-f]{64}$/);
    if (!row.local_path) continue;
    const localPath = resolve(root, row.local_path);
    expect(fileSha256(localPath)).toBe(row.sha256);
  }
}, 30_000);

test("owner inputs preserve aliases, candidates, unknowns, and provenance", () => {
  const extract = String.raw`
import csv,json,sys
def rows(path):
    with open(path,encoding="utf-8-sig",newline="") as stream:
        return list(csv.DictReader(stream))
masters=rows(sys.argv[1]); requirements=rows(sys.argv[2]); sources=rows(sys.argv[3]); rules=rows(sys.argv[4])
ids={"APP-004","APP-005","APP-011","APP-014","APP-015","APP-016","APP-017","APP-019","SAN-023","SAN-024","SAN-025","CTRL-ENTRY-A","NET-AP-R09","NET-AP-R14","SENSOR-GAS-R04","SENSOR-001"}
print(json.dumps({
  "masters":[row for row in masters if row["equipment_id"] in ids],
  "requirements":[row for row in requirements if row["equipment_id"] in ids],
  "sources":[row for row in sources if row["source_id"] in {"OWNER-INBOX-20260814-001","OWNER-INPUT-APP019-20260815","ZOYLIGHT-D1-20250815-001"} or row["local_path"].startswith("drawings/evidence/owner-input-20260814/")],
  "rules":[row for row in rules if row["rule_id"] in {"ELEC-DES-020","ELEC-DES-022","ELEC-DES-033","ELEC-DES-035","ELEC-DES-037","ELEC-DES-040","ELEC-DES-041"}],
},ensure_ascii=False))
`;
  const run = Bun.spawnSync([
    "python3", "-c", extract,
    resolve(root, "pipeline/decisions/equipment-register.csv"),
    resolve(root, "pipeline/decisions/equipment-installation-requirements.csv"),
    resolve(root, "pipeline/decisions/source-evidence-register.csv"),
    resolve(root, "pipeline/decisions/elec-design-rules.csv"),
  ], { cwd: root });
  expect(run.exitCode, run.stderr.toString()).toBe(0);
  const data = JSON.parse(run.stdout.toString()) as {
    masters: Array<Record<string, string>>;
    requirements: Array<Record<string, string>>;
    sources: Array<Record<string, string>>;
    rules: Array<Record<string, string>>;
  };
  const master = (id: string) => data.masters.find((row) => row.equipment_id === id);
  const requirement = (id: string, key: string) =>
    data.requirements.find((row) => row.equipment_id === id && row.parameter_key === key);
  const rule = (id: string) => data.rules.find((row) => row.rule_id === id);

  expect(master("APP-014")).toMatchObject({
    model: "WS7FSB0C1C",
    quantity: "1",
    procurement_status: "candidate",
    decision_status: "partial",
  });
  expect(requirement("APP-014", "functional_alias")).toMatchObject({
    value_text: "APP-015",
    status: "confirmed",
  });
  expect(requirement("APP-014", "installation_manual")).toMatchObject({
    value_text: "8001325304_D",
    status: "confirmed",
    blocks_release: "no",
  });
  expect(requirement("APP-014", "price_comparison_alternative")).toMatchObject({
    value_text: "WS7060BC1C",
    status: "confirmed",
  });
  expect(master("APP-015")).toMatchObject({
    quantity: "0",
    procurement_status: "not_applicable",
    decision_status: "superseded",
    schedule_included: "no",
  });
  expect(requirement("APP-015", "alias_of")).toMatchObject({ value_text: "APP-014", status: "confirmed" });
  expect(master("APP-017")).toMatchObject({
    manufacturer: "Siemens",
    model: "WG54M7D20W + WQ55M7U20W",
    procurement_status: "candidate",
  });
  expect(requirement("APP-017", "coordination_clear_width")?.value_number).toBe("650");
  expect(requirement("APP-017", "coordination_clear_height")?.value_number).toBe("1900");
  expect(requirement("APP-017", "coordination_clear_depth")?.value_number).toBe("800");
  expect(requirement("APP-017", "independent_power_plug_count")?.value_number).toBe("2");
  expect(requirement("APP-017", "stack_connector_order_number")?.value_text).toBe("17008829");
  expect(requirement("APP-015", "rated_power")).toMatchObject({
    value_number: "0",
    status: "confirmed",
    blocks_release: "no",
  });
  for (const key of ["water_required", "drain_required", "gas_required", "ventilation_required"]) {
    expect(requirement("APP-015", key), `APP-015.${key}`).toMatchObject({
      value_text: "not_applicable_alias",
      status: "confirmed",
      blocks_release: "no",
    });
  }

  expect(requirement("APP-004", "rated_power")).toMatchObject({ status: "pending", blocks_release: "no" });
  expect(requirement("APP-004", "design_reserve_power_w")).toMatchObject({
    value_number: "2600",
    value_origin: "research_conclusion",
    status: "confirmed",
  });
  expect(requirement("APP-011", "appliance_plug_rating_a")?.value_number).toBe("16");
  expect(requirement("APP-011", "wall_socket_rating_a")?.value_number).toBe("16");
  expect(requirement("APP-017", "wall_socket_rating_a")?.value_number).toBe("10");
  expect(requirement("APP-017", "wall_socket_quantity")?.value_number).toBe("2");
  expect(requirement("APP-017", "shared_branch_circuit_permission")?.value_text).toBe("yes");
  expect(requirement("APP-017", "branch_breaker_rating_a")?.value_number).toBe("16");
  expect(requirement("APP-004", "candidate_power_gs3")).toMatchObject({ value_number: "2120", status: "candidate" });
  expect(requirement("APP-004", "candidate_power_e1_prima_exp")).toMatchObject({
    value_text: "1600_or_2600_by_version",
    status: "candidate",
  });
  expect(requirement("APP-005", "rated_power")).toMatchObject({ value_number: "1200", status: "confirmed" });
  expect(master("APP-016")).toMatchObject({ procurement_status: "candidate", decision_status: "partial" });
  for (const [equipmentId, referenceLabel] of [
    ["SAN-023", "MF287"],
    ["SAN-024", "CZ356"],
    ["SAN-025", "CZ028"],
  ] as const) {
    expect(master(equipmentId)).toMatchObject({
      procurement_status: "purchased_delivery_unverified",
      decision_status: "partial",
      quantity: "1",
    });
    expect(master(equipmentId).model).toContain(referenceLabel);
    expect(data.requirements.some((row) => row.equipment_id === equipmentId && row.parameter_key === "exact_manufacturer_article_and_sku" && row.status === "pending" && row.blocks_release === "yes")).toBe(true);
  }
  expect(requirement("APP-016", "rated_power")).toMatchObject({ status: "pending", blocks_release: "yes" });
  expect(master("APP-019")).toMatchObject({
    item_name: "冰淇淋机",
    quantity: "1",
    model: "Mini Lussino 4080 / CREAMi 220 V（候选）",
    procurement_status: "candidate",
    decision_status: "partial",
  });
  for (const key of [
    "rated_power", "water_required", "drain_required", "appliance_plug_rating_a",
    "wall_socket_rating_a", "branch_breaker_rating_a", "dedicated_branch_circuit",
    "interface_center_coordinates",
  ]) {
    expect(requirement("APP-019", key), `APP-019.${key}`).toMatchObject({
      status: "pending",
      blocks_release: "yes",
    });
  }
  expect(requirement("APP-019", "appliance_plug_rating_a")?.value_text).toBe("unknown");
  expect(requirement("APP-019", "interface_center_coordinates")?.value_text).toBe("unknown");
  expect(requirement("APP-019", "project_clear_width")).toMatchObject({ value_number: "500", status: "confirmed", blocks_release: "no" });
  expect(requirement("APP-019", "project_clear_depth")).toMatchObject({ value_number: "450", status: "confirmed", blocks_release: "no" });
  expect(requirement("APP-019", "project_clear_height")).toMatchObject({ value_number: "450", status: "confirmed", blocks_release: "no" });
  expect(requirement("APP-019", "minimum_support_load")).toMatchObject({ value_number: "20", status: "confirmed", blocks_release: "no" });
  const iceCreamSource = data.sources.find((row) => row.source_id === "OWNER-INPUT-APP019-20260815");
  expect(iceCreamSource).toMatchObject({
    local_path: "drawings/evidence/OWNER-INPUT-APP-019-20260815.md",
    status: "confirmed_owner_input_identity_only",
    formal_ifc_write_allowed: "no",
  });

  expect(master("CTRL-ENTRY-A")).toMatchObject({ procurement_status: "candidate", decision_status: "candidate" });
  expect(master("CTRL-ENTRY-A")?.source_ids).toContain("ZOYLIGHT-D1-20250815-001");
  for (const key of ["exact_sku", "terminal_diagram", "minimum_backbox_clear_depth"]) {
    expect(requirement("CTRL-ENTRY-A", key), `CTRL-ENTRY-A.${key}`).toMatchObject({
      status: "pending",
      blocks_release: "yes",
    });
  }
  for (const id of ["NET-AP-R09", "NET-AP-R14"]) {
    expect(requirement(id, "final_model"), `${id}.final_model`).toMatchObject({ status: "pending", blocks_release: "yes" });
    expect(requirement(id, "cable_continuity_test"), `${id}.cable_continuity_test`).toMatchObject({ status: "pending" });
  }
  expect(master("SENSOR-GAS-R04")).toMatchObject({ model: "", decision_status: "pending" });
  expect(requirement("SENSOR-GAS-R04", "gas_company_approval")).toMatchObject({ status: "pending", blocks_release: "yes" });
  expect(requirement("SENSOR-GAS-R04", "consultation_candidates")).toMatchObject({ status: "candidate" });
  expect(requirement("SENSOR-001", "final_detector_type_and_model")).toMatchObject({ status: "pending" });

  expect(rule("ELEC-DES-020")?.status).toBe("superseded_by_device_first_selection");
  expect(rule("ELEC-DES-022")?.status).toBe("superseded_by_gas_authority_gate");
  for (const id of ["ELEC-DES-033", "ELEC-DES-035"]) {
    expect(rule(id)?.status).toBe("research_conclusion");
  }
  expect(rule("ELEC-DES-037")?.value).toContain("actual_length_x_W_per_m");
  expect(rule("ELEC-DES-040")?.value).toContain("minimum_2_per_usable_side");
  expect(rule("ELEC-DES-041")?.status).toBe("confirmed_process_exact_model_unknown");

  const snapshot = data.sources.find((row) => row.source_id === "OWNER-INBOX-20260814-001");
  expect(snapshot).toMatchObject({
    local_path: "drawings/evidence/OWNER-INPUT-20260814-SYNCED.md",
    sha256: "83b1fc425ec8632c997455a7924aa8669eb2dc06856c6ca4de0eea7da2d89053",
    status: "verified_owner_input_snapshot",
    formal_ifc_write_allowed: "no",
  });
  const snapshotText = readFileSync(resolve(root, snapshot!.local_path), "utf8");
  expect(snapshotText).toContain("record_status: synced_archive");
  expect(snapshotText).toContain("final_reconciliation_completed_at: 2026-08-14T11:01:15+08:00");
  expect(snapshotText).toContain("source_payload_sha256: 0dde421f925a8f16a0a52706db29b544d115cebd596b850cb9f18dea26da27ec");
  expect(data.sources.filter((row) => row.local_path.startsWith("drawings/evidence/owner-input-20260814/"))).toHaveLength(18);
  const lighting = data.sources.find((row) => row.source_id === "ZOYLIGHT-D1-20250815-001");
  expect(lighting).toMatchObject({
    local_path: "drawings/evidence/lighting/ZOYLIGHT-D1-20250815.dwg",
    sha256: "f857849f015fcec58f6f8d52a842ce341f6512dc1feabfbb9f465d1bdfcb4522",
    status: "verified_owner_directed_design_source_pending_circuit_extraction",
    formal_ifc_write_allowed: "no",
  });
  const gasTemplate = readFileSync(resolve(root, "drawings/evidence/A106-燃气公司咨询模板.md"), "utf8");
  expect(gasTemplate).toContain("仅在贵司确认后选用");
  expect(gasTemplate).toContain("如均不适用，请填准入品牌和准确型号");
  expect(gasTemplate).toContain("不代表燃气公司已批准");
});

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
  expect(report.scoped_ifc_object_count).toBe(165);
  expect(report.target_source_count).toBeGreaterThan(1);
  expect(report.dry_run).toBe(true);
  expect(readFileSync(sourcePath)).toEqual(before);
}, 30_000);

test("Geberit receipt shortfalls remain blockers while withdrawn drain components are not current-scheme blockers", () => {
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

  const currentBlockerKeys: Record<string, string[]> = {
    "SAN-014": ["completed_wc_height", "wall_anchor_design", "boxing_finished_thickness", "nuna_hidden_interface_coordination"],
    "SAN-003": ["panel_centre_height", "finished_surface_thickness", "service_access_detail"],
  };
  for (const [equipmentId, keys] of Object.entries(currentBlockerKeys)) {
    for (const key of keys) {
      expect(get(equipmentId, key), `${equipmentId}.${key}`).toMatchObject({
        status: "pending",
        blocks_release: "yes",
      });
    }
  }

  const withdrawnKeys: Record<string, string[]> = {
    "SAN-001": ["actual_cut_length", "placement_mode", "drain_centre_coordinate", "tile_layout_coordination"],
    "DRAIN-GEB-002": ["final_installation_level", "waterproofing_system", "site_pipe_orientation"],
    "DRAIN-GEB-001": ["final_assembly_orientation"],
    "DRAIN-GEB-004": ["washing_machine_hose_connection", "installation_height", "floor_slope", "waterproofing_detail"],
  };
  for (const [equipmentId, keys] of Object.entries(withdrawnKeys)) {
    for (const key of keys) {
      expect(get(equipmentId, key), `${equipmentId}.${key}`).toMatchObject({
        value_text: "not_applicable_current_scheme",
        status: "not_applicable",
        blocks_release: "no",
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
dxf=pathlib.Path(sys.argv[5])
print(json.dumps({
  "masters":[row for row in equipment if row["equipment_id"] in {"DW-M05","DW-M06"}],
  "requirements":[row for row in requirements if row["equipment_id"] in {"DW-M05","DW-M06"}],
  "source":source,
  "actual_sha256":hashlib.sha256(cad.read_bytes()).hexdigest(),
  "dxf_sha256":hashlib.sha256(dxf.read_bytes()).hexdigest(),
},ensure_ascii=False))
`;
  const cadPath = resolve(root, "drawings/evidence/RIMADESIO-official-Sail-monorotaia.dwg");
  const dxfPath = resolve(root, "drawings/evidence/RIMADESIO-Sail-monorotaia.dxf");
  const run = Bun.spawnSync(["python3", "-c", extract, equipmentPath, requirementsPath, sourcesPath, cadPath, dxfPath], { cwd: root });
  expect(run.exitCode, run.stderr.toString()).toBe(0);
  const data = JSON.parse(run.stdout.toString()) as {
    masters: Array<Record<string, string>>;
    requirements: Array<Record<string, string>>;
    source: Record<string, string>;
    actual_sha256: string;
    dxf_sha256: string;
  };
  expect(data.source).toMatchObject({
    sha256: "4eeaeb27b0668bb22b9b4fb191091e07eabbb2b259d90cedcdc45584a1f68e52",
    status: "verified_official_family_cad",
    formal_ifc_write_allowed: "no",
  });
  expect(data.source.locator).toContain("76b40456832e6c3e217c8bf58dc58ee49214beca710269660dc11a087209adb4");
  expect(data.dxf_sha256).toBe("76b40456832e6c3e217c8bf58dc58ee49214beca710269660dc11a087209adb4");
  expect(data.source.evidence).toContain("轨道长 2011/2037/4022 mm");
  expect(data.source.notes).toContain("RIMADESIO-Sail-monorotaia-mechanical-audit.json");
  expect(data.actual_sha256).toBe(data.source.sha256);
  for (const equipmentId of ["DW-M05", "DW-M06"]) {
    expect(data.masters.find((row) => row.equipment_id === equipmentId)?.source_ids).toContain("RIMADESIO-SAIL-MONOROTAIA-001");
  }
  const requirement = (equipmentId: string, key: string) =>
    data.requirements.find((row) => row.equipment_id === equipmentId && row.parameter_key === key);
  expect(requirement("DW-M05", "operation_type")).toMatchObject({ status: "confirmed", blocks_release: "no" });
  expect(requirement("DW-M05", "nominal_width")).toMatchObject({ value_number: "1000", status: "confirmed", blocks_release: "no" });
  expect(requirement("DW-M05", "nominal_height")).toMatchObject({ value_number: "2400", status: "confirmed", blocks_release: "no" });
  expect(requirement("DW-M06", "operation_type")).toMatchObject({ status: "candidate", blocks_release: "yes" });
  expect(requirement("DW-M06", "nominal_width")).toMatchObject({ value_number: "2000", status: "confirmed", blocks_release: "no" });
  expect(requirement("DW-M06", "nominal_height")).toMatchObject({ status: "pending", blocks_release: "yes" });
});
