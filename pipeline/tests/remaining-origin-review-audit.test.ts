import { expect, test } from "bun:test";
import { resolve } from "node:path";

const root = resolve(import.meta.dir, "../..");
const modulePath = resolve(root, "pipeline/scripts/remaining_origin_review_audit.py");

test("remaining origin queues never authorize a generic write", () => {
  const source = `
import importlib.util, pathlib, sys
sys.path.insert(0, str(pathlib.Path(${JSON.stringify(modulePath)}).parent))
spec = importlib.util.spec_from_file_location("remaining_origin_review_audit", ${JSON.stringify(modulePath)})
module = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = module
spec.loader.exec_module(module)
print(module.SURFACE_CLASSES)
print(module.SERVICE_CLASSES)
print(module.SECONDARY_BUILDING_CLASSES)
`;
  const result = Bun.spawnSync(["python3", "-c", source], { cwd: root });
  expect(result.exitCode).toBe(0);
  expect(result.stdout.toString()).toContain("IfcCovering");
  expect(result.stdout.toString()).toContain("IfcLightFixture");
  expect(result.stdout.toString()).toContain("IfcBuildingElementProxy");
  expect(result.stdout.toString()).not.toContain("automatic_write_allowed=True");
});

test("remaining origin audit accepts a configurable positive tolerance", () => {
  const result = Bun.spawnSync([
    "python3",
    modulePath,
    "--help",
  ], { cwd: root });
  expect(result.exitCode).toBe(0);
  expect(result.stdout.toString()).toContain("--tolerance-mm");
});

test("service subgroups remain descriptive anchor-review labels", () => {
  const source = `
import importlib.util, json, pathlib, sys
sys.path.insert(0, str(pathlib.Path(${JSON.stringify(modulePath)}).parent))
spec = importlib.util.spec_from_file_location("remaining_origin_review_audit", ${JSON.stringify(modulePath)})
module = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = module
spec.loader.exec_module(module)
class FakeProduct:
    def __init__(self, kind):
        self.kind = kind
    def is_a(self, query=None):
        return self.kind if query is None else self.kind == query
product = FakeProduct("IfcLightFixture")
queue = module.review_queue(product)
subgroup = module.review_subgroup(product, queue[0])
print(json.dumps({"queue": queue, "subgroup": subgroup}))
`;
  const result = Bun.spawnSync(["python3", "-c", source], { cwd: root });
  expect(result.exitCode).toBe(0);
  const parsed = JSON.parse(result.stdout.toString());
  expect(parsed.queue[0]).toBe("service_installation_anchor_review");
  expect(parsed.subgroup[0]).toBe("light_fixture_anchor_review");
  expect(parsed.subgroup.join(" ")).not.toContain("integer snap");
});

test("assembly subgroups follow represented child roles", () => {
  const source = `
import importlib.util, json, pathlib, sys
sys.path.insert(0, str(pathlib.Path(${JSON.stringify(modulePath)}).parent))
spec = importlib.util.spec_from_file_location("remaining_origin_review_audit", ${JSON.stringify(modulePath)})
module = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = module
spec.loader.exec_module(module)
class Child:
    def __init__(self, kind): self.kind=kind
    def is_a(self, query=None): return self.kind if query is None else self.kind == query
class Relation:
    def __init__(self, kinds): self.RelatedObjects=[Child(kind) for kind in kinds]
class Assembly:
    Name="BayWindow"
    IsDecomposedBy=[Relation(["IfcWall","IfcWall"])]
    def is_a(self, query=None): return "IfcElementAssembly" if query is None else query == "IfcElementAssembly"
module.ifcopenshell.util.element.get_container=lambda product: None
print(json.dumps(module.review_subgroup(Assembly(), "secondary_building_element_anchor_review")))
`;
  const result = Bun.spawnSync(["python3", "-c", source], { cwd: root });
  expect(result.exitCode).toBe(0);
  expect(JSON.parse(result.stdout.toString())[0]).toBe(
    "wall_assembly_anchor_review",
  );
});
