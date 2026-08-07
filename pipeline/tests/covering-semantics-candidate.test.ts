import { expect, test } from "bun:test";
import { resolve } from "node:path";

const root = resolve(import.meta.dir, "../..");
const modulePath = resolve(root, "pipeline/scripts/covering_semantics_candidate.py");

test("baseboard semantics require explicit naming and matching geometry", () => {
  const source = `
import importlib.util, json, pathlib, sys
sys.path.insert(0, str(pathlib.Path(${JSON.stringify(modulePath)}).parent))
spec = importlib.util.spec_from_file_location("covering_semantics_candidate", ${JSON.stringify(modulePath)})
module = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = module
spec.loader.exec_module(module)
check = module.baseboard_evidence
print(json.dumps({
  "accepted": check("Baseboard.001", None, [20, 1200, 100]),
  "wrong_name": check("Covering", None, [20, 1200, 100]),
  "too_high": check("Baseboard", None, [20, 1200, 150]),
  "already_typed": check("Baseboard", "FLOORING", [20, 1200, 100]),
}))
`;
  const result = Bun.spawnSync(["python3", "-c", source], { cwd: root });
  expect(result.exitCode).toBe(0);
  const parsed = JSON.parse(result.stdout.toString());
  expect(parsed.accepted[0]).toBe(true);
  expect(parsed.wrong_name[0]).toBe(false);
  expect(parsed.too_high[0]).toBe(false);
  expect(parsed.already_typed[0]).toBe(false);
});
