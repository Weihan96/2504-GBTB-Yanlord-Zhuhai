import { expect, test } from "bun:test";
import { resolve } from "node:path";

const root = resolve(import.meta.dir, "../..");
const modulePath = resolve(root, "pipeline/scripts/normalization_exit_audit.py");

function runPython(source: string): string {
  const result = Bun.spawnSync(["python3", "-c", source], { cwd: root });
  expect(result.exitCode).toBe(0);
  return result.stdout.toString().trim();
}

test("normalization exit audit parses one canonical 15-step table", () => {
  const output = runPython(`
import importlib.util, json, pathlib, sys
sys.path.insert(0, str(pathlib.Path(${JSON.stringify(modulePath)}).parent))
spec = importlib.util.spec_from_file_location("normalization_exit_audit", ${JSON.stringify(modulePath)})
module = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = module
spec.loader.exec_module(module)
text = "\\n".join(f"| {i} | Stage {i} | Operation {i} | Gate {i} |" for i in range(1, 16))
print(json.dumps(module.parse_standard_steps(text)))
`);
  const parsed = JSON.parse(output);
  expect(parsed).toHaveLength(15);
  expect(parsed[0]).toMatchObject({ step: 1, stage: "Stage 1" });
  expect(parsed[14]).toMatchObject({ step: 15, acceptance: "Gate 15" });
});

test("normalization exit audit scopes unresolved decisions without absorbing drawing PM", () => {
  const output = runPython(`
import importlib.util, json, pathlib, sys
sys.path.insert(0, str(pathlib.Path(${JSON.stringify(modulePath)}).parent))
spec = importlib.util.spec_from_file_location("normalization_exit_audit", ${JSON.stringify(modulePath)})
module = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = module
spec.loader.exec_module(module)
rows = [{"decision_id": "COORD-X"}, {"decision_id": "FLOOR-X"}, {"decision_id": "A103-X"}]
print(json.dumps([module.is_normalization_decision(row) for row in rows]))
`);
  expect(JSON.parse(output)).toEqual([true, true, false]);
});

test("furniture origin exceptions require an exact implemented decision boundary", () => {
  const output = runPython(`
import importlib.util, json, pathlib, sys
sys.path.insert(0, str(pathlib.Path(${JSON.stringify(modulePath)}).parent))
spec = importlib.util.spec_from_file_location("normalization_exit_audit", ${JSON.stringify(modulePath)})
module = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = module
spec.loader.exec_module(module)
records = [
    {"global_id": "A", "normalization_disposition": "controlled_exception_keep_complex_origin", "review_required": False, "automatic_write_allowed": False},
    {"global_id": "B", "normalization_disposition": "controlled_exception_keep_complex_origin", "review_required": False, "automatic_write_allowed": False},
]
exact = [{"decision_id": "FURNITURE-ANCHOR-C003-LOOSE", "status": "implemented", "object_guid": "A; B"}]
partial = [{"decision_id": "FURNITURE-ANCHOR-C003-LOOSE", "status": "implemented", "object_guid": "A"}]
print(json.dumps({
    "exact": sorted(module.approved_furniture_origin_exceptions(exact, records)),
    "partial": sorted(module.approved_furniture_origin_exceptions(partial, records)),
}))
`);
  expect(JSON.parse(output)).toEqual({ exact: ["A", "B"], partial: [] });
});

test("direct host-relative Opening origins require the implemented delegation rule", () => {
  const output = runPython(`
import importlib.util, json, pathlib, sys
sys.path.insert(0, str(pathlib.Path(${JSON.stringify(modulePath)}).parent))
spec = importlib.util.spec_from_file_location("normalization_exit_audit", ${JSON.stringify(modulePath)})
module = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = module
spec.loader.exec_module(module)
records = [
    {"global_id": "A", "normalization_disposition": "delegated_to_host_placement", "review_required": False, "automatic_write_allowed": False},
    {"global_id": "B", "normalization_disposition": "independent_opening_anchor_review", "review_required": True, "automatic_write_allowed": False},
]
implemented = [{"decision_id": "COORD-OPENING-C003-HOSTREL", "scope": "opening-origin", "status": "implemented", "object_guid": "A"}]
pending = [{"decision_id": "COORD-OPENING-C003-HOSTREL", "scope": "opening-origin", "status": "pending", "object_guid": "A"}]
print(json.dumps({
    "implemented": sorted(module.approved_delegated_opening_origins(implemented, records)),
    "pending": sorted(module.approved_delegated_opening_origins(pending, records)),
}))
`);
  expect(JSON.parse(output)).toEqual({ implemented: ["A"], pending: [] });
});

test("independent Opening queue requires the exact pending decision boundary", () => {
  const output = runPython(`
import importlib.util, json, pathlib, sys
sys.path.insert(0, str(pathlib.Path(${JSON.stringify(modulePath)}).parent))
spec = importlib.util.spec_from_file_location("normalization_exit_audit", ${JSON.stringify(modulePath)})
module = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = module
spec.loader.exec_module(module)
records = [
    {"global_id": "A", "normalization_disposition": "independent_opening_anchor_review", "review_required": True, "automatic_write_allowed": False},
    {"global_id": "B", "normalization_disposition": "delegated_to_host_placement", "review_required": False, "automatic_write_allowed": False},
]
exact = [{"decision_id": "COORD-OPENING-C003-INDEPENDENT", "scope": "opening-origin", "status": "pending", "review_required": "no", "object_guid": "A"}]
wrong = [{"decision_id": "COORD-OPENING-C003-INDEPENDENT", "scope": "opening-origin", "status": "pending", "review_required": "no", "object_guid": "B"}]
print(json.dumps({
    "exact": sorted(module.recorded_independent_opening_origins(exact, records)),
    "wrong": sorted(module.recorded_independent_opening_origins(wrong, records)),
}))
`);
  expect(JSON.parse(output)).toEqual({ exact: ["A"], wrong: [] });
});

test("implemented independent Opening batch requires its exact stable GlobalId boundary", () => {
  const output = runPython(`
import importlib.util, json, pathlib, sys
sys.path.insert(0, str(pathlib.Path(${JSON.stringify(modulePath)}).parent))
spec = importlib.util.spec_from_file_location("normalization_exit_audit", ${JSON.stringify(modulePath)})
module = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = module
spec.loader.exec_module(module)
exact_ids = "0_eLd9WVfFLvk6I8goe8it; 0vxak9xlD8Z97XrPYRGunC; 2LCZDONiLBShgK9CpOImZJ"
exact = [{"decision_id": "COORD-OPENING-C003-INDEPENDENT", "scope": "opening-origin", "status": "implemented", "review_required": "no", "object_guid": exact_ids}]
wrong_ids = [{"decision_id": "COORD-OPENING-C003-INDEPENDENT", "scope": "opening-origin", "status": "implemented", "review_required": "no", "object_guid": "0_eLd9WVfFLvk6I8goe8it"}]
pending = [{"decision_id": "COORD-OPENING-C003-INDEPENDENT", "scope": "opening-origin", "status": "pending", "review_required": "no", "object_guid": exact_ids}]
print(json.dumps({
    "exact": module.implemented_independent_opening_batch(exact),
    "wrong_ids": module.implemented_independent_opening_batch(wrong_ids),
    "pending": module.implemented_independent_opening_batch(pending),
}))
`);
  expect(JSON.parse(output)).toEqual({ exact: true, wrong_ids: false, pending: false });
});

test("Direction exceptions come from one exact implemented decision row", () => {
  const output = runPython(`
import importlib.util, json, pathlib, sys
sys.path.insert(0, str(pathlib.Path(${JSON.stringify(modulePath)}).parent))
spec = importlib.util.spec_from_file_location("normalization_exit_audit", ${JSON.stringify(modulePath)})
module = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = module
spec.loader.exec_module(module)
exact = [{"decision_id": "COORD-DIRECTION-C003-001", "scope": "direction-noise-exception", "status": "implemented", "review_required": "no", "object_guid": "#43042; #2081575"}]
wrong_scope = [{"decision_id": "COORD-DIRECTION-C003-001", "scope": "direction-noise", "status": "implemented", "review_required": "no", "object_guid": "#43042; #2081575"}]
print(json.dumps({
    "exact": sorted(module.approved_direction_exception_ids(exact)),
    "wrong_scope": sorted(module.approved_direction_exception_ids(wrong_scope)),
}))
`);
  expect(JSON.parse(output)).toEqual({ exact: [43042, 2081575], wrong_scope: [] });
});

test("controlled exception register requires every canonical decision", () => {
  const output = runPython(`
import importlib.util, json, pathlib, sys
sys.path.insert(0, str(pathlib.Path(${JSON.stringify(modulePath)}).parent))
spec = importlib.util.spec_from_file_location("normalization_exit_audit", ${JSON.stringify(modulePath)})
module = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = module
spec.loader.exec_module(module)
rows = [
    {"decision_id": decision_id, "scope": scope, "status": "implemented"}
    for decision_id, scope in module.REQUIRED_CONTROLLED_EXCEPTION_DECISIONS.items()
]
complete = module.controlled_exception_decision_gate(rows)
incomplete = module.controlled_exception_decision_gate(rows[:-1])
print(json.dumps({"complete": complete, "incomplete": incomplete}))
`);
  const parsed = JSON.parse(output);
  expect(parsed.complete.pass).toBe(true);
  expect(parsed.complete.recorded).toHaveLength(9);
  expect(parsed.incomplete.pass).toBe(false);
  expect(parsed.incomplete.missing_or_invalid).toHaveLength(1);
});

test("project decision register contains every controlled exception", () => {
  const output = runPython(`
import importlib.util, json, pathlib, sys
sys.path.insert(0, str(pathlib.Path(${JSON.stringify(modulePath)}).parent))
spec = importlib.util.spec_from_file_location("normalization_exit_audit", ${JSON.stringify(modulePath)})
module = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = module
spec.loader.exec_module(module)
rows = module.read_decisions(pathlib.Path("pipeline/decisions/p0-review.csv"))
print(json.dumps(module.controlled_exception_decision_gate(rows)))
`);
  expect(JSON.parse(output)).toEqual({
    pass: true,
    recorded: [
      "COORD-DIRECTION-C003-001",
      "COORD-FLOW-C003-PVC110-BUNDLE",
      "COORD-MATL-C003-GB01",
      "COORD-SOCKET-C003-CONTROLLED",
      "COORD-WALL-C003-G3",
      "COORD-WALL-C005-CURVE",
      "FLOOR-TILE-A105-DRAIN-001",
      "FLOOR-TILE-A105-GROUT-001",
      "FURNITURE-ANCHOR-C003-LOOSE",
    ],
    missing_or_invalid: [],
  });
});

test("socket origin exceptions require exact implemented appliance records", () => {
  const output = runPython(`
import importlib.util, json, pathlib, sys
sys.path.insert(0, str(pathlib.Path(${JSON.stringify(modulePath)}).parent))
spec = importlib.util.spec_from_file_location("normalization_exit_audit", ${JSON.stringify(modulePath)})
module = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = module
spec.loader.exec_module(module)
ids = sorted(module.SOCKET_CONTROLLED_EXCEPTION_IDS)
records = [
    *[{"global_id": value, "within_review_tolerance": False, "ifc_class": "IfcElectricAppliance", "review_queue": "service_installation_anchor_review"} for value in ids],
    {"global_id": "B", "within_review_tolerance": False, "ifc_class": "IfcWall", "review_queue": "service_installation_anchor_review"},
]
exact = [{"decision_id": "COORD-SOCKET-C003-CONTROLLED", "scope": "electric-appliance-origin-exception", "status": "implemented", "review_required": "no", "object_guid": "; ".join(ids)}]
wrong_scope = [{"decision_id": "COORD-SOCKET-C003-CONTROLLED", "scope": "object-placement", "status": "implemented", "review_required": "no", "object_guid": "; ".join(ids)}]
wrong_class = [{"decision_id": "COORD-SOCKET-C003-CONTROLLED", "scope": "electric-appliance-origin-exception", "status": "implemented", "review_required": "no", "object_guid": "B"}]
print(json.dumps({
    "exact": sorted(module.approved_socket_origin_exceptions(exact, records)),
    "wrong_scope": sorted(module.approved_socket_origin_exceptions(wrong_scope, records)),
    "wrong_class": sorted(module.approved_socket_origin_exceptions(wrong_class, records)),
}))
`);
  expect(JSON.parse(output)).toEqual({
    exact: [
      "0laejMoxn8Lu_X3FZaCXmi",
      "27MTenki57DQsfMryX_1U0",
      "2OOjqQDMHDjRcXQCniWXnp",
      "3KXtmVvejA78j_iAS$kydj",
    ],
    wrong_scope: [],
    wrong_class: [],
  });
});

test("PVC110 origin exceptions require the exact bundled-product decision", () => {
  const output = runPython(`
import importlib.util, json, pathlib, sys
sys.path.insert(0, str(pathlib.Path(${JSON.stringify(modulePath)}).parent))
spec = importlib.util.spec_from_file_location("normalization_exit_audit", ${JSON.stringify(modulePath)})
module = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = module
spec.loader.exec_module(module)
ids = sorted(module.PVC110_CONTROLLED_EXCEPTION_IDS)
records = [
    *[{"global_id": value, "within_review_tolerance": False, "ifc_class": "IfcFlowSegment", "review_queue": "service_installation_anchor_review"} for value in ids],
]
exact = [{"decision_id": "COORD-FLOW-C003-PVC110-BUNDLE", "scope": "flow-segment-bundle-exception", "status": "implemented", "review_required": "no", "object_guid": "; ".join(ids)}]
wrong = [{"decision_id": "COORD-FLOW-C003-PVC110-BUNDLE", "scope": "flow-segment-bundle-exception", "status": "implemented", "review_required": "no", "object_guid": ids[0]}]
print(json.dumps({
    "exact": sorted(module.approved_pvc110_origin_exceptions(exact, records)),
    "wrong": sorted(module.approved_pvc110_origin_exceptions(wrong, records)),
}))
`);
  expect(JSON.parse(output)).toEqual({
    exact: ["0bfVg4Ys1CevZs$qxhkXTo", "178mqyyzzFowLcbXcH6prO"],
    wrong: [],
  });
});

test("C003 downstream handoffs require exact non-overlapping responsibility rows", () => {
  const output = runPython(`
import importlib.util, json, pathlib, sys
sys.path.insert(0, str(pathlib.Path(${JSON.stringify(modulePath)}).parent))
spec = importlib.util.spec_from_file_location("normalization_exit_audit", ${JSON.stringify(modulePath)})
module = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = module
spec.loader.exec_module(module)
eligible = set()
rows = []
for index, (decision_id, owner) in enumerate(module.ORIGIN_HANDOFF_DECISIONS.items()):
    global_id = f"G{index}"
    eligible.add(global_id)
    rows.append({
        "decision_id": decision_id,
        "scope": "c003-origin-responsibility",
        "status": "delegated",
        "review_required": "yes",
        "object_guid": global_id,
        "proposed_value": f"{owner} 负责后续确认",
    })
exact = module.approved_origin_handoffs(rows, eligible)
overlap = [dict(row) for row in rows]
overlap[1]["object_guid"] = overlap[0]["object_guid"]
wrong_scope = [dict(row) for row in rows]
wrong_scope[0]["scope"] = "object-placement"
print(json.dumps({
    "exact": {key: sorted(value) for key, value in exact.items()},
    "overlap": module.approved_origin_handoffs(overlap, eligible),
    "wrong_scope": module.approved_origin_handoffs(wrong_scope, eligible),
}))
`);
  const parsed = JSON.parse(output);
  expect(Object.keys(parsed.exact).sort()).toEqual(
    ["A104", "A105", "WFIN", "PLUM", "ELEC", "RCP1", "INT1", "DET1"].sort(),
  );
  expect(parsed.overlap).toEqual({});
  expect(parsed.wrong_scope).toEqual({});
});
