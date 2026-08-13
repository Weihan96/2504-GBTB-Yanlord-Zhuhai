import { expect, test } from "bun:test";
import { resolve } from "node:path";

const root = resolve(import.meta.dir, "../..");
const script = resolve(root, "pipeline/scripts/rcp1_route_readiness_candidate.py");
const output = resolve(root, "build/rcp1/route-readiness-candidate.test.json");

test("RCP1 route graph accepts shared openings and ordered multi-hop paths", async () => {
  const source = await Bun.file(script).text();
  expect(source).toContain('"blender_preview_is_source_of_truth": False');
  expect(source).toContain('"equipment_positions_may_move": False');
  expect(source).not.toContain("itertools.permutations");
  expect(source).not.toContain("model.write(");

  const process = Bun.spawn(
    [
      "python3",
      script,
      "--input",
      resolve(root, "2504 GBTB Yanlord Zhuhai.ifc"),
      "--hvac-report",
      resolve(root, "build/rcp1/hvac-remodel-candidate.json"),
      "--legacy-report",
      resolve(root, "build/rcp1/legacy-base-audit.json"),
      "--delivery-dwg",
      resolve(root, "../图纸/矩阵纵横/D户型交付竣工图.dwg"),
      "--routes",
      resolve(root, "pipeline/decisions/rcp1-hvac-route-register.csv"),
      "--waypoints",
      resolve(root, "pipeline/decisions/rcp1-hvac-route-waypoints.csv"),
      "--equipment-register",
      resolve(root, "pipeline/decisions/equipment-register.csv"),
      "--requirements",
      resolve(root, "pipeline/decisions/equipment-installation-requirements.csv"),
      "--output",
      output,
    ],
    { cwd: root, stdout: "pipe", stderr: "pipe" },
  );
  const [exitCode, stderr] = await Promise.all([
    process.exited,
    new Response(process.stderr).text(),
  ]);
  expect(stderr).toBe("");
  expect(exitCode).toBe(0);

  const report = await Bun.file(output).json();
  expect(report.equipment).toHaveLength(6);
  expect(report.openings).toHaveLength(7);
  expect(report.airside_readiness.a02_and_a06_both_serve_dining).toBe(true);
  const routes = Object.fromEntries(report.confirmed_route_graph.map((item: any) => [item.route_id, item]));
  expect(routes["RCP1-SERVICE-A02"].waypoints.map((item: any) => item.anchor_id)).toEqual(["A02", "H03"]);
  expect(routes["RCP1-SERVICE-A03"].waypoints.map((item: any) => item.anchor_id)).toEqual(["A03", "H04", "H02"]);
  expect(report.interface_roles.H01.role).toBe("outdoor_unit_interface_at_existing_opening");
  expect(report.interface_roles.H07.role).toBe("condensate_discharge_interface_at_existing_opening");
  expect(report.interface_roles.H01.ifc_entity_remains).toBe("IfcOpeningElement");
  expect(report.interface_roles.H07.ifc_entity_remains).toBe("IfcOpeningElement");
  expect(report.equipment_opening_mapping).toHaveLength(6);
  expect(report.interface_coverage).toHaveLength(7);
  const mappingByEquipment = Object.fromEntries(
    report.equipment_opening_mapping.map((item: any) => [item.equipment_id, item]),
  );
  expect(mappingByEquipment.A01).toMatchObject({
    opening_id: "H05", status: "geometry_candidate", opening_rank_by_clearance: 1,
    formal_ifc_write_allowed: false,
  });
  expect(mappingByEquipment.A02).toMatchObject({
    opening_id: "H03", status: "user_confirmed", formal_ifc_write_allowed: false,
  });
  expect(mappingByEquipment.A03).toMatchObject({
    opening_id: "H04", status: "user_confirmed", formal_ifc_write_allowed: false,
  });
  expect(mappingByEquipment.A04).toMatchObject({
    opening_id: "H06", status: "geometry_candidate", opening_rank_by_clearance: 1,
  });
  expect(mappingByEquipment.A05).toMatchObject({
    opening_id: "H07", status: "geometry_candidate", opening_rank_by_clearance: 1,
  });
  expect(mappingByEquipment.A06).toMatchObject({
    opening_id: "H03", status: "geometry_candidate", opening_rank_by_clearance: 1,
  });
  expect(mappingByEquipment.A01.minimum_clearance_candidate_mm).toBeCloseTo(249.562, 2);
  expect(mappingByEquipment.A04.minimum_clearance_candidate_mm).toBeCloseTo(297.868, 2);
  expect(mappingByEquipment.A05.minimum_clearance_candidate_mm).toBeCloseTo(620.958, 2);
  expect(mappingByEquipment.A06.minimum_clearance_candidate_mm).toBeCloseTo(343.452, 2);
  const interfaceById = Object.fromEntries(
    report.interface_coverage.map((item: any) => [item.opening_id, item]),
  );
  expect(interfaceById.H01).toMatchObject({
    role: "user_confirmed_outdoor_unit_interface",
    coverage_status: "user_confirmed_route_or_endpoint",
  });
  expect(interfaceById.H02.route_ids).toContain("RCP1-SERVICE-A03");
  expect(interfaceById.H03.equipment_mappings).toEqual([
    { equipment_id: "A02", status: "user_confirmed" },
    { equipment_id: "A06", status: "geometry_candidate" },
  ]);
  expect(interfaceById.H07).toMatchObject({
    role: "user_confirmed_condensate_endpoint_with_A05_candidate",
    coverage_status: "user_confirmed_route_or_endpoint",
  });
  expect(report.legacy_reference.component_topology_matches).toBe(true);
  expect(report.legacy_reference.maximum_component_dimension_difference_mm).toBeLessThan(0.016);
  expect(report.gates.confirmed_shared_and_multihop_graph_ready).toBe(true);
  expect(report.gates.six_equipment_opening_mappings_registered).toBe(true);
  expect(report.gates.seven_interfaces_have_controlled_roles).toBe(true);
  expect(report.gates.confirmed_and_candidate_mapping_split_preserved).toBe(true);
  expect(report.gates.manufacturer_constraints_consumed_from_ssot).toBe(true);
  expect(report.gates.a01_a04_refrigerant_and_condensate_nominals_available).toBe(true);
  const a02Requirements = Object.fromEntries(report.manufacturer_interface_inputs.by_equipment.A02.requirements
    .map((item: any) => [item.requirement_name, item.value]));
  expect(a02Requirements.gas_pipe_od).toBe(12.7);
  expect(a02Requirements.liquid_pipe_od).toBe(6.35);
  expect(a02Requirements.drain_pipe_od).toBe(32);
  const a05Requirements = Object.fromEntries(report.manufacturer_interface_inputs.by_equipment.A05.requirements
    .map((item: any) => [item.requirement_name, item.value]));
  expect(a05Requirements.service_access_width_min).toBe(450);
  expect(a05Requirements.supply_duct_length_min).toBe(1000);
  expect(report.manufacturer_interface_inputs.by_equipment.A06.release_blocking_requirements)
    .toContain("nominal_body_width");
  expect(report.gates.formal_ifc_write_allowed).toBe(false);
}, 30_000);
