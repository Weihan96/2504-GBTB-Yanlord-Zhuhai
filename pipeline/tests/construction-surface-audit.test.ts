import { expect, test } from "bun:test";
import { resolve } from "node:path";

const root = resolve(import.meta.dir, "../..");
const modulePath = resolve(root, "pipeline/scripts/construction_surface_audit.py");

test("construction surface orientation does not imply semantics", () => {
  const source = `
import importlib.util, json, pathlib, sys
sys.path.insert(0, str(pathlib.Path(${JSON.stringify(modulePath)}).parent))
spec = importlib.util.spec_from_file_location("construction_surface_audit", ${JSON.stringify(modulePath)})
module = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = module
spec.loader.exec_module(module)
print(json.dumps({
  "floor": module.geometry_orientation([3000, 4000, 20], 100),
  "wall": module.geometry_orientation([20, 4000, 2700], 100),
  "other": module.geometry_orientation([400, 400, 400], 100),
}))
`;
  const result = Bun.spawnSync(["python3", "-c", source], { cwd: root });
  expect(result.exitCode).toBe(0);
  expect(JSON.parse(result.stdout.toString())).toEqual({
    floor: "horizontal_thin",
    wall: "vertical_thin",
    other: "other",
  });
});

test("explicit semantics alone never authorize an anchor write", () => {
  const source = `
import importlib.util, json, pathlib, sys
sys.path.insert(0, str(pathlib.Path(${JSON.stringify(modulePath)}).parent))
spec = importlib.util.spec_from_file_location("construction_surface_audit", ${JSON.stringify(modulePath)})
module = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = module
spec.loader.exec_module(module)
print(json.dumps({
  "explicit_integer_vertex": module.automatic_anchor_write_allowed("explicit_ifc_semantics", True, True),
  "candidate_integer_vertex": module.automatic_anchor_write_allowed("orientation_candidate_requires_review", True, True),
}))
`;
  const result = Bun.spawnSync(["python3", "-c", source], { cwd: root });
  expect(result.exitCode).toBe(0);
  const parsed = JSON.parse(result.stdout.toString());
  expect(parsed.explicit_integer_vertex[0]).toBe(false);
  expect(parsed.candidate_integer_vertex[0]).toBe(false);
  expect(parsed.explicit_integer_vertex[1]).toContain("class-specific confirmed construction anchor");
});

test("construction surface semantics inherit an explicit type PredefinedType", () => {
  const source = `
import importlib.util, json, pathlib, sys
sys.path.insert(0, str(pathlib.Path(${JSON.stringify(modulePath)}).parent))
spec = importlib.util.spec_from_file_location("construction_surface_audit", ${JSON.stringify(modulePath)})
module = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = module
spec.loader.exec_module(module)

class TypeObject:
    PredefinedType = "CLADDING"

class Product:
    PredefinedType = None
    Name = "Covering"
    def is_a(self, name=None):
        return "IfcCovering" if name is None else name == "IfcCovering"

module.ifcopenshell.util.element.get_type = lambda product: TypeObject()
print(json.dumps(module.semantic_status(Product(), "vertical_thin")))
`;
  const result = Bun.spawnSync(["python3", "-c", source], { cwd: root });
  expect(result.exitCode).toBe(0);
  expect(JSON.parse(result.stdout.toString())).toEqual([
    "explicit_ifc_semantics",
    "effective PredefinedType=CLADDING from type",
  ]);
});

test("construction surfaces are grouped for review without inferring semantics", () => {
  const source = `
import importlib.util, json, pathlib, sys
sys.path.insert(0, str(pathlib.Path(${JSON.stringify(modulePath)}).parent))
spec = importlib.util.spec_from_file_location("construction_surface_audit", ${JSON.stringify(modulePath)})
module = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = module
spec.loader.exec_module(module)
batch = module.review_batch
print(json.dumps({
  "explicit": batch("IfcCovering", "Plane", "FLOORING", None, "horizontal_thin", [0,0,0], [1000,1000,0]),
  "baseboard": batch("IfcCovering", "Baseboard.001", None, None, "vertical_thin", [0,0,0], [2000,20,100]),
  "tile_below": batch("IfcCovering", "Covering", None, "TerrazzoMosaicTile", "horizontal_thin", [0,0,-1020], [800,800,-1000]),
  "tile_ffl": batch("IfcCovering", "Covering", None, "TerrazzoMosaicTile", "horizontal_thin", [0,0,-40], [800,800,-10]),
  "lower_slab": batch("IfcSlab", "Slab", None, None, "other", [0,0,-230], [3000,3000,0]),
  "upper_slab": batch("IfcSlab", "Slab", None, None, "other", [0,0,2770], [3000,3000,3000]),
}))
`;
  const result = Bun.spawnSync(["python3", "-c", source], { cwd: root });
  expect(result.exitCode).toBe(0);
  expect(JSON.parse(result.stdout.toString())).toEqual({
    explicit: "covering_existing_explicit_semantics",
    baseboard: "covering_baseboard_name_candidates",
    tile_below: "covering_tile_instances_below_ffl",
    tile_ffl: "covering_tile_instances_near_ffl",
    lower_slab: "slab_below_or_at_ffl_candidates",
    upper_slab: "slab_above_ceiling_candidates",
  });
});

test("review subgroups describe evidence without assigning construction semantics", () => {
  const source = `
import importlib.util, json, pathlib, sys
sys.path.insert(0, str(pathlib.Path(${JSON.stringify(modulePath)}).parent))
spec = importlib.util.spec_from_file_location("construction_surface_audit", ${JSON.stringify(modulePath)})
module = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = module
spec.loader.exec_module(module)
subgroup = module.review_subgroup
print(json.dumps({
  "wall_surface": subgroup("IfcCovering", "FFL", "covering_vertical_candidates", [0,0,100], [20,3000,2700], [20,3000,2600]),
  "high_strip": subgroup("IfcCovering", "STUD", "covering_vertical_candidates", [0,0,2420], [9,3500,2760], [9,3500,340]),
  "raised_surface": subgroup("IfcCovering", "FFL", "covering_horizontal_candidates", [0,0,580], [810,620,600], [810,620,20]),
  "upper_slab": subgroup("IfcSlab", "RFL", "slab_above_ceiling_candidates", [0,0,2770], [3000,3000,3000], [3000,3000,230]),
  "service_furniture": subgroup("IfcFurniture", "KITCHEN", "furniture_fixed_or_loose_review", [0,0,0], [600,600,900], [600,600,900]),
  "bath_fixture": subgroup("IfcSanitaryTerminal", "BATHM", "sanitary_installation_anchor_review", [0,0,0], [600,600,900], [600,600,900]),
}))
`;
  const result = Bun.spawnSync(["python3", "-c", source], { cwd: root });
  expect(result.exitCode).toBe(0);
  const parsed = JSON.parse(result.stdout.toString());
  expect(parsed.wall_surface[0]).toBe("covering_vertical_full_height_surfaces");
  expect(parsed.high_strip[0]).toBe("covering_high_level_vertical_strips");
  expect(parsed.raised_surface[0]).toBe("covering_raised_horizontal_surfaces");
  expect(parsed.upper_slab[0]).toBe("slab_rfl_geometry_review");
  expect(parsed.service_furniture[0]).toBe("furniture_service_zone_fixed_or_loose_review");
  expect(parsed.bath_fixture[0]).toBe("sanitary_bathroom_anchor_review");
  expect(Object.values(parsed).flat().join(" ")).not.toContain("FLOORING");
  expect(Object.values(parsed).flat().join(" ")).not.toContain("CLADDING");
});

test("horizontal footprints are associated by geometric overlap", () => {
  const source = `
import importlib.util, json, pathlib, sys
import numpy as np
from shapely.geometry import box
sys.path.insert(0, str(pathlib.Path(${JSON.stringify(modulePath)}).parent))
spec = importlib.util.spec_from_file_location("construction_surface_audit", ${JSON.stringify(modulePath)})
module = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = module
spec.loader.exec_module(module)
vertices=np.array([[0,0,0],[10,0,0],[10,10,0],[0,10,0]],dtype=float)
faces=np.array([[0,1,2],[0,2,3]],dtype=int)
footprint=module.projected_footprint(vertices,faces)
spaces=[
 {"global_id":"A","long_name":"Room A","geometry":box(0,0,6,10)},
 {"global_id":"B","long_name":"Room B","geometry":box(6,0,10,10)},
]
print(json.dumps(module.footprint_space_overlaps(footprint,spaces)))
`;
  const result = Bun.spawnSync(["python3", "-c", source], { cwd: root });
  expect(result.exitCode).toBe(0);
  const parsed = JSON.parse(result.stdout.toString());
  expect(parsed.map((item: { global_id: string }) => item.global_id)).toEqual(["A", "B"]);
  expect(parsed[0].covering_overlap_ratio).toBeCloseTo(0.6);
  expect(parsed[1].covering_overlap_ratio).toBeCloseTo(0.4);
});
