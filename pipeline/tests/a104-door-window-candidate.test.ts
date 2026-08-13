import { expect, test } from "bun:test";
import { createHash } from "node:crypto";
import { resolve } from "node:path";

const root = resolve(import.meta.dir, "../..");
const modulePath = resolve(root, "pipeline/scripts/a104_door_window_candidate.py");

function runPython(body: string) {
  return Bun.spawnSync(["python3", "-c", `
import importlib.util, json, pathlib, sys
sys.path.insert(0, str(pathlib.Path(${JSON.stringify(modulePath)}).parent))
spec = importlib.util.spec_from_file_location("a104_door_window_candidate", ${JSON.stringify(modulePath)})
module = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = module
spec.loader.exec_module(module)
${body}
`], { cwd: root });
}

test("A104 opening interval audit distinguishes overlap and gap", () => {
  const result = runPython(`print(json.dumps([
    module.interval_relationship((0, 1180), (580, 1760)),
    module.interval_relationship((0, 1465), (1465.05, 2930.05)),
  ]))`);
  expect(result.exitCode).toBe(0);
  const values = JSON.parse(result.stdout.toString());
  expect(values[0].kind).toBe("OVERLAP");
  expect(values[0].overlap_mm).toBe(600);
  expect(values[1].kind).toBe("GAP_OR_TOUCH");
  expect(values[1].gap_mm).toBeCloseTo(0.05, 6);
});

test("A104 operation labels remain human-readable without inventing types", () => {
  const result = runPython(`print(json.dumps([
    module.format_operation("SINGLE_SWING_LEFT"),
    module.format_operation("TRIPLE_PANEL_RIGHT"),
    module.format_operation(""),
  ], ensure_ascii=False))`);
  expect(result.exitCode).toBe(0);
  expect(JSON.parse(result.stdout.toString())).toEqual(["单扇左开", "三扇右分格", "待确认"]);
});

test("A104 known semantic conflict and unhosted pair require review", () => {
  const result = runPython(`print(json.dumps([
    module.review_metadata({"global_id": module.MASTER_BEDROOM_DOOR_CANDIDATE}),
    module.review_metadata({"global_id": "2D5BPoo2XFSvhTdfPenCh7"}),
  ], ensure_ascii=False))`);
  expect(result.exitCode).toBe(0);
  const values = JSON.parse(result.stdout.toString());
  expect(values[0].review_group).toBe("A104-R01");
  expect(values[0].review_required).toBe("yes");
  expect(values[1].review_group).toBe("A104-R02");
});

test("A104 keeps M07 swing review open after its identity is confirmed", () => {
  const result = runPython(`
records = [{"global_id": module.MASTER_BEDROOM_DOOR, "review_group": "", "review_required": "no", "review_question": "", "confidence": 1.0}]
module.close_confirmed_reviews(records, [])
print(json.dumps(records[0], ensure_ascii=False))
`);
  expect(result.exitCode).toBe(0);
  const value = JSON.parse(result.stdout.toString());
  expect(value.review_group).toBe("A104-R03");
  expect(value.review_required).toBe("yes");
  expect(value.review_question).toContain("OperationType=NOTDEFINED");
  expect(value.review_question).toContain("Master A 仅为暂定优选");
  expect(value.review_question).toContain("仍待取证");
});

test("A104 uses the current IFC hash by default and rejects a mismatched caller freeze", async () => {
  const ifcPath = resolve(root, "2504 GBTB Yanlord Zhuhai.ifc");
  const currentSha = createHash("sha256")
    .update(new Uint8Array(await Bun.file(ifcPath).arrayBuffer()))
    .digest("hex");
  const result = runPython(`
from pathlib import Path
values = [
    module.validate_source_sha(Path(${JSON.stringify(ifcPath)})),
    module.validate_source_sha(Path(${JSON.stringify(ifcPath)}), ${JSON.stringify(currentSha)}),
]
try:
    module.validate_source_sha(Path(${JSON.stringify(ifcPath)}), "0" * 64)
except RuntimeError as error:
    mismatch = str(error)
else:
    raise AssertionError("mismatched caller freeze was accepted")
print(json.dumps({"values": values, "mismatch": mismatch}))
`);
  expect(result.exitCode).toBe(0);
  const value = JSON.parse(result.stdout.toString());
  expect(value.values).toEqual([currentSha, currentSha]);
  expect(value.mismatch).toContain("formal IFC SHA-256 mismatch");
});

test("A104 candidate identifiers sort deterministically by plan position", () => {
  const result = runPython(`
rows = [
  {"global_id":"B", "bbox":{"centre_mm":[100, 500, 0]}},
  {"global_id":"A", "bbox":{"centre_mm":[-100, 500, 0]}},
  {"global_id":"C", "bbox":{"centre_mm":[0, 100, 0]}},
]
print(json.dumps([row["global_id"] for row in sorted(rows, key=module.plan_sort_key)]))
`);
  expect(result.exitCode).toBe(0);
  expect(JSON.parse(result.stdout.toString())).toEqual(["A", "B", "C"]);
});

test("A104 reserves D identifiers for A102 demolition walls", async () => {
  const source = await Bun.file(modulePath).text();
  expect(source).toContain('record["candidate_id"] = f"M{index:02d}"');
  expect(source).toContain("D prefix is reserved by A-102 demolition walls");
});

test("A104 closes superseded human reviews when formal tags names and groups exist", () => {
  const result = runPython(`
model = module.ifcopenshell.open("2504 GBTB Yanlord Zhuhai.ifc")
records = [
  {"global_id": product.GlobalId, "current_tag": str(product.Tag or ""), "candidate_id": str(product.Tag or "")}
  for product in [*model.by_type("IfcDoor"), *model.by_type("IfcWindow")]
]
print(json.dumps(module.formal_semantics_state(model, records), ensure_ascii=False))
`);
  expect(result.exitCode).toBe(0);
  expect(JSON.parse(result.stdout.toString())).toEqual({
    tags_pass: true,
    names_pass: true,
    groups_pass: true,
    complete: true,
  });
});

test("A104 consumes the Sail CAD audit without expanding its IFC write boundary", () => {
  const auditPath = resolve(root, "drawings/evidence/RIMADESIO-Sail-monorotaia-mechanical-audit.json");
  const result = runPython(`
audit = module.load_sail_cad_audit(
    pathlib.Path(${JSON.stringify(auditPath)}),
    "446a9fc8447d38e16dacad7dcce3df097b03f8a457c6bc80fca03d3560891632",
)
print(json.dumps({
    "entity_count": audit["cad_inventory"]["entity_count"],
    "rail_widths": audit["cad_inventory"]["generic_dimensions"]["rail_width_mm"],
    "formal_ifc_write_allowed": audit["evidence_boundary"]["formal_ifc_write_allowed"],
}))
`);
  expect(result.exitCode, result.stderr.toString()).toBe(0);
  expect(JSON.parse(result.stdout.toString())).toEqual({
    entity_count: 1338,
    rail_widths: [2011, 2037, 4022],
    formal_ifc_write_allowed: false,
  });
});
