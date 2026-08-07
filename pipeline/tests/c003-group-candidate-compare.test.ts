import { expect, test } from "bun:test";
import { resolve } from "node:path";

const repositoryRoot = resolve(import.meta.dir, "../..");
const modulePath = resolve(repositoryRoot, "pipeline/scripts/c003_group_candidate_compare.py");

function runPython(source: string): string {
  const result = Bun.spawnSync(["python3", "-c", source], { cwd: repositoryRoot });
  expect(result.exitCode).toBe(0);
  return result.stdout.toString().trim();
}

test("C003 group plan parses approved typed transforms", () => {
  const output = runPython(`
import importlib.util, json, pathlib, sys, tempfile
sys.path.insert(0, str(pathlib.Path(${JSON.stringify(modulePath)}).parent))
spec = importlib.util.spec_from_file_location("c003_group_candidate_compare", ${JSON.stringify(modulePath)})
module = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = module
spec.loader.exec_module(module)
with tempfile.TemporaryDirectory() as directory:
    path = pathlib.Path(directory) / "plan.csv"
    path.write_text("group_id,ifc_class,global_id,cardinalize,translation_x_mm,translation_y_mm,translation_z_mm,snap_bbox_axis,snap_bbox_side,snap_bbox_target_mm,basis,confidence,human_review,status\\nC003-G1,IfcWall,W1,true,1,0,0,,,,test,1.00,no,approved\\n", encoding="utf-8")
    print(json.dumps(module.read_plan(path, "C003-G1")))
`);
  const rows = JSON.parse(output) as Array<Record<string, unknown>>;
  expect(rows).toHaveLength(1);
  expect(rows[0]).toMatchObject({
    global_id: "W1",
    cardinalize: true,
    translation_mm: [1, 0, 0],
    snap_bbox_target_mm: null,
    confidence: 1,
  });
});

test("C003 group plan parses copied Profile control points", () => {
  const output = runPython(`
import importlib.util, json, pathlib, sys, tempfile
sys.path.insert(0, str(pathlib.Path(${JSON.stringify(modulePath)}).parent))
spec = importlib.util.spec_from_file_location("c003_group_candidate_compare", ${JSON.stringify(modulePath)})
module = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = module
spec.loader.exec_module(module)
with tempfile.TemporaryDirectory() as directory:
    path = pathlib.Path(directory) / "plan.csv"
    path.write_text("group_id,ifc_class,global_id,cardinalize,translation_x_mm,translation_y_mm,translation_z_mm,snap_bbox_axis,snap_bbox_side,snap_bbox_target_mm,copy_shared_profile,profile_coordinate_axis,profile_point_indices,profile_coordinate_targets_mm,profile_snap_increment_mm,basis,confidence,human_review,status\\nC003-G2,IfcWall,W2,false,0,0,0,,,,true,x,0;2;5,,1,test,1.00,yes,approved\\n", encoding="utf-8")
    print(json.dumps(module.read_plan(path, "C003-G2")))
`);
  const rows = JSON.parse(output) as Array<Record<string, unknown>>;
  expect(rows[0]).toMatchObject({
    global_id: "W2",
    copy_shared_profile: true,
    profile_coordinate_axis: "x",
    profile_point_indices: [0, 2, 5],
    profile_coordinate_targets_mm: [],
    profile_snap_increment_mm: 1,
  });
});

test("C003 group plan parses explicit Profile targets", () => {
  const output = runPython(`
import importlib.util, json, pathlib, sys, tempfile
sys.path.insert(0, str(pathlib.Path(${JSON.stringify(modulePath)}).parent))
spec = importlib.util.spec_from_file_location("c003_group_candidate_compare", ${JSON.stringify(modulePath)})
module = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = module
spec.loader.exec_module(module)
with tempfile.TemporaryDirectory() as directory:
    path = pathlib.Path(directory) / "plan.csv"
    path.write_text("group_id,ifc_class,global_id,cardinalize,translation_x_mm,translation_y_mm,translation_z_mm,snap_bbox_axis,snap_bbox_side,snap_bbox_target_mm,copy_shared_profile,profile_coordinate_axis,profile_point_indices,profile_coordinate_targets_mm,profile_snap_increment_mm,basis,confidence,human_review,status\\nC003-G3,IfcOpeningElement,O1,false,0,0,0,,,,true,y,0;1,-300.4;-300.4,,test,1.00,no,approved\\n", encoding="utf-8")
    print(json.dumps(module.read_plan(path, "C003-G3")))
`);
  const rows = JSON.parse(output) as Array<Record<string, unknown>>;
  expect(rows[0]).toMatchObject({
    global_id: "O1",
    profile_coordinate_targets_mm: [-300.4, -300.4],
    profile_snap_increment_mm: null,
  });
});

test("Profile control point snapping rejects shifts over the safety limit", () => {
  const output = runPython(`
import importlib.util, json, pathlib, sys
sys.path.insert(0, str(pathlib.Path(${JSON.stringify(modulePath)}).parent))
spec = importlib.util.spec_from_file_location("c003_group_candidate_compare", ${JSON.stringify(modulePath)})
module = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = module
spec.loader.exec_module(module)
accepted = module.snap_profile_coordinate(-800.300537109375, 1.0, 0.5)
try:
    module.snap_profile_coordinate(12.49, 10.0, 0.5)
except RuntimeError as error:
    rejected = str(error)
print(json.dumps({"accepted": accepted, "rejected": rejected}))
`);
  const result = JSON.parse(output) as {
    accepted: [number, number];
    rejected: string;
  };
  expect(result.accepted[0]).toBe(-800);
  expect(result.accepted[1]).toBeCloseTo(0.300537109375, 10);
  expect(result.rejected).toContain("exceeds 0.500000 mm");
});

test("Extrusion depth edits reject shifts over the safety limit", () => {
  const output = runPython(`
import importlib.util, json, pathlib, sys
sys.path.insert(0, str(pathlib.Path(${JSON.stringify(modulePath)}).parent))
spec = importlib.util.spec_from_file_location("c003_group_candidate_compare", ${JSON.stringify(modulePath)})
module = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = module
spec.loader.exec_module(module)
class Solid:
    id = lambda self: 12
    Depth = 100.49
class Product:
    GlobalId = "P1"
    is_a = lambda self, value=None: "IfcCovering" if value is None else value == "IfcCovering"
class Model:
    by_guid = lambda self, value: Product()
module.body_extruded_solid = lambda product: Solid()
accepted = module.apply_extrusion_depth_row(Model(), {"global_id": "P1", "ifc_class": "IfcCovering", "extrusion_depth_target_mm": 100.0})
try:
    module.apply_extrusion_depth_row(Model(), {"global_id": "P1", "ifc_class": "IfcCovering", "extrusion_depth_target_mm": 101.5})
except RuntimeError as error:
    rejected = str(error)
print(json.dumps({"accepted": accepted, "rejected": rejected}))
`);
  const result = JSON.parse(output) as {
    accepted: { shift_mm: number };
    rejected: string;
  };
  expect(result.accepted.shift_mm).toBeCloseTo(-0.49, 10);
  expect(result.rejected).toContain("exceeds 1.000000 mm");
});
