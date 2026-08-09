import { expect, test } from "bun:test";
import { resolve } from "node:path";

const root = resolve(import.meta.dir, "../..");
const scriptPath = resolve(root, "pipeline/scripts/rcp1_hvac_constraint_preview.py");

test("HVAC authoring preview updates geometry without writing IFC or blend", async () => {
  const source = await Bun.file(scriptPath).text();
  expect(source).toContain('COLLECTION_NAME = "RCP1_HVAC_CONSTRAINT_AUTHORING"');
  expect(source).toContain('NODE_GROUP_NAME = "RCP1_HVAC_CONSTRAINT_GEOMETRY"');
  expect(source).toContain('nodes.new("GeometryNodeFilletCurve")');
  expect(source).toContain('nodes.new("GeometryNodeCurveToMesh")');
  expect(source).toContain('SERVICE_EXIT_M = 0.15');
  expect(source).toContain('"equipment_local_positive_x_service_face"');
  expect(source).toContain("orthogonal_world_points");
  expect(source).toContain('"orthogonal_points_world_mm"');
  expect(source).toContain('"fillet_radius_m": FILLET_RADIUS_M');
  expect(source).toContain("automatic_preview_handler");
  expect(source).toContain('bl_idname = "rcp1_hvac.add_bend"');
  expect(source).toContain('bl_idname = "rcp1_hvac.rebuild_preview"');
  expect(source).toContain('bl_idname = "rcp1_hvac.export_candidate"');
  expect(source).toContain('"formal_ifc_write_allowed": False');
  expect(source).toContain('"source_of_truth": "IFC plus approved route and waypoint registers"');
  expect(source).toContain('bl_category = "RCP1 HVAC"');
  expect(source).toContain('"RCP1_HVAC_REMODEL_REVIEW"');
  expect(source).toContain('anchor_id == "H01"');
  expect(source).toContain('anchor_id == "H07"');
  expect(source).toContain('"endpoint_anchors"');
  expect(source).toContain('"ifc_write": False');
  expect(source).toContain('"blend_save": False');
  expect(source).not.toContain("save_and_load_ifc");
  expect(source).not.toContain("bpy.ops.wm.save_as_mainfile");
  expect(source).not.toContain("model.write(");
  expect(source).toContain("obj.show_in_front = False");
  expect(source).toContain("show_xray = False");
  expect(source).toContain("show_wireframes = False");
});
