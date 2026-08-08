import { expect, test } from "bun:test";
import { resolve } from "node:path";

const root = resolve(import.meta.dir, "../..");
const modulePath = resolve(root, "pipeline/scripts/ids_validate.py");

test("IDS compact report keeps failed GlobalIds without STEP payloads", () => {
  const script = `
import importlib.util, json
spec = importlib.util.spec_from_file_location("ids_validate", ${JSON.stringify(modulePath)})
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
raw = {
  "status": False,
  "total_specifications": 1,
  "total_specifications_pass": 0,
  "total_specifications_fail": 1,
  "total_checks": 2,
  "total_checks_pass": 1,
  "total_checks_fail": 1,
  "specifications": [{
    "name": "Door",
    "status": False,
    "total_applicable": 1,
    "total_applicable_pass": 0,
    "total_applicable_fail": 1,
    "total_checks": 2,
    "total_checks_pass": 1,
    "total_checks_fail": 1,
    "requirements": [{
      "description": "Tag required",
      "status": False,
      "total_applicable": 1,
      "total_pass": 0,
      "total_fail": 1,
      "failed_entities": [{"global_id": "ABC", "reason": "empty", "element": "STEP"}]
    }]
  }]
}
print(json.dumps(module.compact_report(raw)))
`;
  const result = Bun.spawnSync(["python3", "-c", script], { cwd: root });
  expect(result.exitCode).toBe(0);
  const report = JSON.parse(result.stdout.toString());
  expect(report.specifications[0].requirements[0].failed_global_ids).toEqual(["ABC"]);
  expect(JSON.stringify(report)).not.toContain("STEP");
});

test("Space Reference follows the IFC4 Pset template data type", async () => {
  const ids = await Bun.file(resolve(root, "pipeline/ids/p0-construction-information.ids")).text();
  const reference = ids.match(/<property dataType="([^"]+)"[^>]*>[\s\S]*?<propertySet><simpleValue>Pset_SpaceCommon<\/simpleValue><\/propertySet>[\s\S]*?<baseName><simpleValue>Reference<\/simpleValue><\/baseName>/);
  expect(reference?.[1]).toBe("IFCIDENTIFIER");
});
