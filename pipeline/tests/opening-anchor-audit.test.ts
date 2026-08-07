import { expect, test } from "bun:test";
import { resolve } from "node:path";

const root = resolve(import.meta.dir, "../..");
const modulePath = resolve(root, "pipeline/scripts/opening_anchor_audit.py");

test("only one unfilled direct child Opening delegates its origin to the host", () => {
  const source = `
import importlib.util, json, pathlib, sys
sys.path.insert(0, str(pathlib.Path(${JSON.stringify(modulePath)}).parent))
spec = importlib.util.spec_from_file_location("opening_anchor_audit", ${JSON.stringify(modulePath)})
module = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = module
spec.loader.exec_module(module)
cases = [(1, 0, True), (1, 0, False), (1, 1, True), (2, 0, True), (0, 0, False)]
print(json.dumps([module.classify_opening_relationship(*case) for case in cases]))
`;
  const result = Bun.spawnSync(["python3", "-c", source], { cwd: root });
  expect(result.exitCode).toBe(0);
  const records = JSON.parse(result.stdout.toString());
  expect(records[0]).toMatchObject({
    normalization_disposition: "delegated_to_host_placement",
    review_required: false,
    automatic_write_allowed: false,
  });
  for (const record of records.slice(1)) {
    expect(record.normalization_disposition).toBe(
      "independent_opening_anchor_review",
    );
    expect(record.review_required).toBe(true);
    expect(record.automatic_write_allowed).toBe(false);
  }
});
