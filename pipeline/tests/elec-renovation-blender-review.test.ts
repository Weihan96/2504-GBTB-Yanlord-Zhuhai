import { expect, test } from "bun:test";
import { readFileSync } from "node:fs";
import { resolve } from "node:path";

const source = readFileSync(resolve(import.meta.dir, "../scripts/elec_renovation_blender_review.py"), "utf8");

test("renovation electrical Blender review uses true depth and separated scopes", () => {
  expect(source).toContain('ROOT_COLLECTION = "RENOVATION_ELEC_ROUND1"');
  expect(source).toContain('"01_BEDSIDE_LIGHT_CYAN"');
  expect(source).toContain('"02_NEW_SOCKET_GREEN"');
  expect(source).toContain('"03_CABINET_POWER_AMBER"');
  expect(source).toContain('"04_KITCHEN_SOCKET_RECHECK_RED"');
  expect(source).toContain('"05_DOORWAY_CONTROL_MAGENTA"');
  expect(source).toContain('"06_BEDROOM_AP_BLUE"');
  expect(source).toContain('"07_LIVING_STUDY_ROUTER_VIOLET"');
  expect(source).toContain('REVIEW_LABEL_COLLECTION = "08_CURRENT_REVIEW_CALLOUTS"');
  expect(source).toContain('"CTRL-ENTRY": ("E302  CTRL-ENTRY"');
  expect(source).toContain('"NS-01": ("E303  NS-01"');
  expect(source).toContain('"A106-AP-R09": ("E304  AP-R09"');
  expect(source).toContain('"NET-ROUTER-LIVING-STUDY": ("E304  ROUTER-R20/R22"');
  expect(source).toContain('obj["review_overlay_only"] = True');
  expect(source).toContain('"Callouts are review overlays at Z=3.15 m; source markers retain their true installation depth."');
  expect(source).toContain("obj.show_in_front = False");
  expect(source).toContain('space.shading.type = "SOLID"');
  expect(source).toContain("space.shading.show_xray = False");
  expect(source).toContain("space.overlay.show_wireframes = False");
  expect(source).toContain("FURNITURE_YELLOW");
  expect(source).toContain("hide_demolition_walls");
  expect(source).toContain('row["ifc_status_candidate"] == "DEMOLISH"');
  expect(source).toContain("obj.hide_set(True)");
  expect(source).toContain("electrical review contains a partial DEMOLISH set");
  expect(source).toContain("loaded DEMOLISH walls remain visible");
  expect(source).not.toContain("model.write(");
  expect(source).not.toContain("bpy.ops.wm.save_as_mainfile");
});
