import { expect, test } from "bun:test";
import { mkdtempSync, readFileSync, rmSync, writeFileSync } from "node:fs";
import { tmpdir } from "node:os";
import { join, resolve } from "node:path";

const root = resolve(import.meta.dir, "../..");
const script = resolve(root, "pipeline/scripts/int1_elevation_blender_render.py");
const source = readFileSync(script, "utf8");

test("INT1 elevation renderer is render-only and restores every mutated state family", () => {
  expect(source).toContain("EXPECTED_VIEW_COUNT = 36");
  expect(source).toContain('OUTPUT_DIR_REL = Path("build/int1/elevations/raw")');
  expect(source).toContain('MANIFEST_REL = Path("build/int1/elevation-render-manifest.json")');
  expect(source).toContain("def capture_scene_state(");
  expect(source).toContain("def restore_scene_state(");
  expect(source).toContain('"render_state"');
  expect(source).toContain('"camera_state"');
  expect(source).toContain('"object_visibility_state"');
  expect(source).toContain('"viewport_shading_state"');
  expect(source).toContain('obj.show_in_front = False');
  expect(source).toContain('space.shading.type = "SOLID"');
  expect(source).toContain('space.shading.show_xray = False');
  expect(source).toContain("def create_sanitized_render_copies(");
  expect(source).toContain("def create_isolated_render_scene(");
  expect(source).toContain("dependency graph contains no Bonsai source meshes");
  expect(source).toContain("visible_width = visible_height * aspect");
  expect(source).toContain("frame_u_min = lateral_center - visible_width / 2.0");
  expect(source).toContain("mesh.from_pydata(vertices, [], faces)");
  expect(source).toContain("Originals are never edited.");
  expect(source).toContain('bpy.ops.render.render(write_still=True, scene=render_scene.name)');
  expect(source).not.toContain("model.write(");
  expect(source).not.toContain("ifc.write(");
  expect(source).not.toContain("save_ifc_file");
  expect(source).not.toContain("bpy.ops.wm.save_as_mainfile");
  expect(source).not.toContain("bpy.ops.wm.save_mainfile");
});

test("INT1 elevation renderer excludes CAD semantics and all registered DEMOLISH walls", () => {
  expect(source).toContain('A102_WALL_REVIEW_REL = Path("pipeline/decisions/a102-wall-review.csv")');
  expect(source).toContain('A102_DEMOLITION_REVIEW_REL = Path("pipeline/decisions/a102-demolition-review.csv")');
  expect(source).toContain('EXCLUDED_IFC_CLASSES = {"IfcSpace", "IfcGrid", "IfcAnnotation", "IfcOpeningElement"}');
  expect(source).toContain("load_demolition_ids(demolition_register)");
  expect(source).toContain('if str(getattr(entity, "GlobalId", "")) in demolition_ids:');
  expect(source).toContain("def is_typed_ceiling(");
  expect(source).toContain('return "CEILING" in {occurrence_predefined, type_predefined}');
  expect(source).toContain('entity_container_name(entity).upper() == "DCL" and "ceiling" in name');
  expect(source).toContain("bpy.ops.render.render(write_still=True, scene=render_scene.name)");
  expect(source).toContain("remove_sanitized_render_copies(bpy, collection, copies)");
});

test("manifest schema validator accepts 36 strict views and rejects visible demolition", () => {
  const python = `
import importlib.util, json, sys
spec = importlib.util.spec_from_file_location("int1_elevation_blender_render", ${JSON.stringify(script)})
module = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = module
spec.loader.exec_module(module)
projected = {"u_min_mm": 0, "u_max_mm": 1, "z_min_mm": 0, "z_max_mm": 1, "depth_min_mm": 0, "depth_max_mm": 1}
obj = {"global_id": "G", "ifc_class": "IfcWall", "name": "W", "container": "DCL", "bbox_min_mm": [0,0,0], "bbox_max_mm": [1,1,1], "projected": projected, "source": "world_bbox"}
view = {"view_id": "01", "sheet_id": "EL-01", "direction": "+Y", "anchor_id": "A1", "space_reference": "R01", "space_global_id": "S", "png": {"path": "x.png", "sha256": "h", "width_px": 1800, "height_px": 1200}, "camera": {"location_mm": [0,0,0], "rotation_euler_rad": [0,0,0], "ortho_scale_mm": 3000, "clip_start_mm": 10, "clip_end_mm": 5000, "resolution_px": [1800,1200]}, "frame": {"u_min_mm": 0, "u_max_mm": 1, "z_min_mm": 0, "z_max_mm": 1}, "objects": [obj], "demolish_visible_count": 0}
manifest = {"source_ifc_sha256": "i", "view_register_sha256": "r", "generated_at": "t", "views": [dict(view, view_id=str(index).zfill(2)) for index in range(1,37)], "demolish_visible_count": 0}
module.validate_manifest_schema(manifest)
manifest["views"][0]["demolish_visible_count"] = 1
try:
    module.validate_manifest_schema(manifest)
except ValueError as error:
    print(json.dumps({"valid_count": 36, "rejected": "visible DEMOLISH walls" in str(error)}))
else:
    raise RuntimeError("schema validator accepted visible DEMOLISH wall")
`;
  const result = Bun.spawnSync(["python3", "-c", python], { cwd: root });
  expect(result.exitCode).toBe(0);
  expect(JSON.parse(result.stdout.toString())).toEqual({ valid_count: 36, rejected: true });
});

test("INT1 elevation manifest retains camera and world/projected geometry evidence", () => {
  for (const field of [
    '"source_ifc_sha256"',
    '"view_register_sha256"',
    '"generated_at"',
    '"view_id"',
    '"sheet_id"',
    '"direction"',
    '"anchor_id"',
    '"space_reference"',
    '"space_global_id"',
    '"png"',
    '"camera"',
    '"frame"',
    '"objects"',
    '"location_mm"',
    '"rotation_euler_rad"',
    '"ortho_scale_mm"',
    '"clip_start_mm"',
    '"clip_end_mm"',
    '"resolution_px"',
    '"u_min_mm"',
    '"u_max_mm"',
    '"z_min_mm"',
    '"z_max_mm"',
    '"depth_min_mm"',
    '"depth_max_mm"',
    '"container"',
    '"bbox_min_mm"',
    '"bbox_max_mm"',
    '"source": "world_bbox"',
    '"demolish_visible_count": 0',
    '"state_restored"',
    '"formal_ifc_write": False',
    '"blend_save": False',
  ]) {
    expect(source).toContain(field);
  }
  expect(source).toContain('for engine in ("BLENDER_WORKBENCH_NEXT", "BLENDER_WORKBENCH")');
  expect(source).toContain("def validate_manifest_schema(");
  expect(source).toContain("validate_manifest_schema(manifest)");
});

test("36-row fixture accepts aliases, sorts views and normalizes four directions", () => {
  const fixtureDir = mkdtempSync(join(tmpdir(), "int1-elevation-render-"));
  const fixture = join(fixtureDir, "views.csv");
  try {
    const directions = ["N", "E", "S", "W"];
    const rows = ["view_number,el_sheet,anchor_id,ifc_anchor_x_mm,ifc_anchor_y_mm,view_direction,space_guid"];
    for (let index = 36; index >= 1; index -= 1) {
      rows.push(`${index},EL-${String(Math.ceil(index / 4)).padStart(2, "0")},A${Math.ceil(index / 4)},${index * 10},${-index * 10},${directions[(index - 1) % 4]},SPACE-${index}`);
    }
    writeFileSync(fixture, `${rows.join("\n")}\n`);
    const python = `
import importlib.util, json, pathlib, sys
spec = importlib.util.spec_from_file_location("int1_elevation_blender_render", ${JSON.stringify(script)})
module = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = module
spec.loader.exec_module(module)
rows = module.load_view_rows(pathlib.Path(${JSON.stringify(fixture)}))
print(json.dumps({"count": len(rows), "first": rows[0]["view_id"], "last": rows[-1]["view_id"], "directions": sorted(set(row["direction"] for row in rows))}))
`;
    const result = Bun.spawnSync(["python3", "-c", python], { cwd: root });
    expect(result.exitCode).toBe(0);
    expect(JSON.parse(result.stdout.toString())).toEqual({
      count: 36,
      first: "01",
      last: "36",
      directions: ["+X", "+Y", "-X", "-Y"],
    });
  } finally {
    rmSync(fixtureDir, { recursive: true, force: true });
  }
});
