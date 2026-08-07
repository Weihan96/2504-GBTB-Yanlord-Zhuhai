import { expect, test } from "bun:test";
import { resolve } from "node:path";

const root = resolve(import.meta.dir, "../..");
const script = resolve(root, "pipeline/scripts/structural_origin_reset_candidate.py");
const targets = resolve(root, "pipeline/decisions/c003-covering-origin-reset-targets.csv");

test("covering origin targets are ten unique integer geometry anchors", () => {
  const code = `
import importlib.util, json, pathlib, sys
sys.path.insert(0, str(pathlib.Path(${JSON.stringify(script)}).parent))
spec = importlib.util.spec_from_file_location("structural", ${JSON.stringify(script)})
module = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = module
spec.loader.exec_module(module)
rows = module.read_targets(pathlib.Path(${JSON.stringify(targets)}))
print(json.dumps([{"class": row["expected_class"], "global_id": row["global_id"], "target": row["target_mm"].tolist()} for row in rows]))
`;
  const result = Bun.spawnSync(["python3", "-c", code], { cwd: root });
  expect(result.exitCode).toBe(0);
  const rows = JSON.parse(result.stdout.toString());
  expect(rows).toHaveLength(10);
  expect(new Set(rows.map((row: any) => row.global_id)).size).toBe(10);
  expect(rows.every((row: any) => row.class === "IfcCovering")).toBe(true);
  expect(rows.every((row: any) => row.target.every(Number.isInteger))).toBe(true);
});

test("new OwnerHistory is accepted only when attached to allowed roots", () => {
  const code = `
import importlib.util, json, pathlib, sys
sys.path.insert(0, str(pathlib.Path(${JSON.stringify(script)}).parent))
spec = importlib.util.spec_from_file_location("structural", ${JSON.stringify(script)})
module = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = module
spec.loader.exec_module(module)
source = module.ifcopenshell.file(schema="IFC4")
candidate = module.ifcopenshell.file(schema="IFC4")
person = candidate.create_entity("IfcPerson")
organization = candidate.create_entity("IfcOrganization", Name="Org")
person_org = candidate.create_entity("IfcPersonAndOrganization", ThePerson=person, TheOrganization=organization)
application = candidate.create_entity("IfcApplication", ApplicationDeveloper=organization, Version="1", ApplicationFullName="Test", ApplicationIdentifier="TEST")
history = candidate.create_entity("IfcOwnerHistory", OwningUser=person_org, OwningApplication=application, ChangeAction="ADDED", CreationDate=1)
project = candidate.create_entity("IfcProject", GlobalId="0abcdefghijklmnopqrstuv", OwnerHistory=history, Name="P")
accepted = module.scoped_owner_history_additions(source, candidate, {project.GlobalId})
rejected = module.scoped_owner_history_additions(source, candidate, set())
print(json.dumps({"accepted": accepted, "rejected": rejected}))
`;
  const result = Bun.spawnSync(["python3", "-c", code], { cwd: root });
  expect(result.exitCode).toBe(0);
  const parsed = JSON.parse(result.stdout.toString());
  expect(parsed.accepted.scoped).toBe(true);
  expect(parsed.accepted.count).toBe(1);
  expect(parsed.rejected.scoped).toBe(false);
});
