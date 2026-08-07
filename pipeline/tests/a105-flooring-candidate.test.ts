import { expect, test } from "bun:test";
import { resolve } from "node:path";

const root = resolve(import.meta.dir, "../..");
const modulePath = resolve(root, "pipeline/scripts/a105_flooring_candidate.py");

test("A105 tile level classification separates current flooring from lower references", () => {
  const source = `
import importlib.util, json, pathlib, sys
sys.path.insert(0, str(pathlib.Path(${JSON.stringify(modulePath)}).parent))
spec = importlib.util.spec_from_file_location("a105_flooring_candidate", ${JSON.stringify(modulePath)})
module = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = module
spec.loader.exec_module(module)
print(json.dumps([
    module.tile_level(-40.0, 5.7),
    module.tile_level(-1020.0, -1000.0),
    module.tile_level(580.0, 600.0),
]))
`;
  const result = Bun.spawnSync(["python3", "-c", source], { cwd: root });
  expect(result.exitCode).toBe(0);
  expect(JSON.parse(result.stdout.toString())).toEqual([
    "near_ffl_flooring",
    "below_ffl_reference",
    "outside_a105_scope",
  ]);
});

test("A105 lower deletion boundary includes each hosted drain Opening", () => {
  const source = `
import importlib.util, json, pathlib, sys
sys.path.insert(0, str(pathlib.Path(${JSON.stringify(modulePath)}).parent))
spec = importlib.util.spec_from_file_location("a105_flooring_candidate", ${JSON.stringify(modulePath)})
module = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = module
spec.loader.exec_module(module)
print(json.dumps({
    "keep": sorted(module.expected_removed_global_ids("keep")),
    "delete": sorted(module.expected_removed_global_ids("delete")),
    "tiles": sorted(module.LOWER_TILE_GLOBAL_IDS),
    "openings": sorted(module.LOWER_DEPENDENT_OPENING_HOSTS),
    "delete_roots": sorted(module.expected_removed_root_global_ids("delete", [
        {"void_relationship_global_ids": ["REL-A"]},
        {"void_relationship_global_ids": ["REL-B"]},
    ])),
}))
`;
  const result = Bun.spawnSync(["python3", "-c", source], { cwd: root });
  expect(result.exitCode).toBe(0);
  const payload = JSON.parse(result.stdout.toString());
  expect(payload.keep).toEqual([]);
  expect(payload.tiles).toHaveLength(4);
  expect(payload.openings).toHaveLength(4);
  expect(payload.delete).toEqual(
    [...payload.tiles, ...payload.openings].sort(),
  );
  expect(payload.delete_roots).toEqual(
    [...payload.delete, "REL-A", "REL-B"].sort(),
  );
});

test("A105 drain slot is owned by eight tile Boolean openings", () => {
  const source = `
import importlib.util, json, pathlib, sys
sys.path.insert(0, str(pathlib.Path(${JSON.stringify(modulePath)}).parent))
spec = importlib.util.spec_from_file_location("a105_flooring_candidate", ${JSON.stringify(modulePath)})
module = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = module
spec.loader.exec_module(module)
print(json.dumps({
    name: {
        "tiles": sorted(group["low_y_tile_ids"] | group["high_y_tile_ids"]),
        "openings": sorted(group["swept_opening_ids"] | group["tessellated_opening_ids"]),
        "swept": len(group["swept_opening_ids"]),
        "tessellated": len(group["tessellated_opening_ids"]),
    }
    for name, group in module.DRAIN_GAP_GROUPS.items()
}))
`;
  const result = Bun.spawnSync(["python3", "-c", source], { cwd: root });
  expect(result.exitCode).toBe(0);
  const payload = JSON.parse(result.stdout.toString());
  expect(Object.keys(payload).sort()).toEqual([
    "guest_bathroom",
    "main_bathroom",
  ]);
  for (const group of Object.values(payload) as Array<{
    tiles: string[];
    openings: string[];
    swept: number;
    tessellated: number;
  }>) {
    expect(group.tiles).toHaveLength(4);
    expect(group.openings).toHaveLength(4);
    expect(group.swept).toBe(2);
    expect(group.tessellated).toBe(2);
  }
});
