import { expect, test } from "bun:test";
import { mkdtempSync, mkdirSync, readFileSync, writeFileSync } from "node:fs";
import { tmpdir } from "node:os";
import { join, resolve } from "node:path";
import { spawnSync } from "node:child_process";
import { createHash } from "node:crypto";

const root = resolve(import.meta.dir, "../..");
const script = join(root, "pipeline/scripts/int1_elevation_candidate.py");
const source = readFileSync(script, "utf8");

test("elevation compiler preserves non-integer evidence and never writes IFC", () => {
  expect(source).toContain('FLOAT_TOLERANCE_MM = 0.01');
  expect(source).toContain('HIGH_TOLERANCE_MM = 0.1');
  expect(source).toContain('automatic_ifc_write_allowed": False');
  expect(source).toContain('space_top_used_as_ceiling": False');
  expect(source).not.toContain("ifcopenshell.api");
  expect(source).not.toContain("model.write(");
});

test("elevation compiler assembles 36 views into nine official sheets", () => {
  const temporary = mkdtempSync(join(tmpdir(), "int1-elevation-"));
  mkdirSync(join(temporary, "pipeline/decisions"), { recursive: true });
  mkdirSync(join(temporary, "build/int1/elevations/raw"), { recursive: true });
  const ifc = join(temporary, "formal.ifc");
  writeFileSync(ifc, "IFC fixture\n");
  const ifcHash = createHash("sha256").update(readFileSync(ifc)).digest("hex");
  const counts: Record<string, number> = {"EL-01":3,"EL-02":2,"EL-03":4,"EL-04":4,"EL-05":4,"EL-06":5,"EL-07":4,"EL-08":6,"EL-09":4};
  const directions = ["+Y", "+X", "-Y", "-X"];
  const header = "view_id,sheet_id,official_title,anchor_id,direction,source_handle,source_viewport_handle,source_scale,source_anchor_handle,source_x_mm,source_y_mm,ifc_x_mm,ifc_y_mm,space_reference,space_global_id,review_status,source_locator";
  const rows: string[] = [header];
  const views: any[] = [];
  const png = Buffer.from("iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAQAAAC1HAwCAAAAC0lEQVR42mNk+A8AAQUBAScY42YAAAAASUVORK5CYII=", "base64");
  let number = 1;
  for (const [sheet, count] of Object.entries(counts)) {
    for (let index = 0; index < count; index++, number++) {
      const id = String(number).padStart(2, "0");
      const pngPath = join(temporary, `build/int1/elevations/raw/${id}.png`);
      writeFileSync(pngPath, png);
      const pngHash = createHash("sha256").update(png).digest("hex");
      rows.push(`${id},${sheet},Fixture Elevation,A${number},${directions[(number-1)%4]},H${id},243CB4,1:50@A2,C${id},0,0,0,0,R01,SPACE,candidate,fixture`);
      views.push({
        view_id:id,sheet_id:sheet,direction:directions[(number-1)%4],anchor_id:`A${number}`,space_reference:"R01",space_global_id:"SPACE",
        png:{path:`build/int1/elevations/raw/${id}.png`,sha256:pngHash,width_px:1,height_px:1},
        camera:{location_mm:[0,0,1400],rotation_euler_rad:[0,0,0],ortho_scale_mm:3000,clip_start_mm:1,clip_end_mm:3000,resolution_px:[1,1]},
        frame:{u_min_mm:0,u_max_mm:3000.25,z_min_mm:0,z_max_mm:2700},demolish_visible_count:0,
        objects:[{global_id:`G${id}`,ifc_class:"IfcFurniture",name:"Fixture",container:"FUR",bbox_min_mm:[0,0,0],bbox_max_mm:[1000.25,500,900],source:"world_bbox",projected:{u_min_mm:0,u_max_mm:1000.25,z_min_mm:0,z_max_mm:900,depth_min_mm:0,depth_max_mm:500}}]
      });
    }
  }
  const register = join(temporary, "pipeline/decisions/int1-elevation-view-register.csv");
  writeFileSync(register, rows.join("\n") + "\n");
  const registerHash = createHash("sha256").update(readFileSync(register)).digest("hex");
  const manifest = join(temporary, "build/int1/elevation-render-manifest.json");
  writeFileSync(manifest, JSON.stringify({source_ifc_sha256:ifcHash,view_register_sha256:registerHash,views}));
  const result = spawnSync("python3", [script,"--root",temporary,"--input-ifc",ifc,"--views",register,"--manifest",manifest,"--skip-render-png"], {encoding:"utf8"});
  expect(result.status, result.stderr).toBe(0);
  const report = JSON.parse(readFileSync(join(temporary,"build/int1/elevation-sheet-candidate.json"),"utf8"));
  expect(report.summary.sheet_count).toBe(9);
  expect(report.summary.view_count).toBe(36);
  expect(report.summary.magenta_object_count).toBe(36);
  expect(report.gates.automatic_ifc_write_allowed).toBe(false);
  expect(report.gates.construction_release_ready).toBe(false);
  expect(readFileSync(join(temporary,"drawings/elevations/EL-08-ifc-elevation-candidate.svg"),"utf8")).toContain('data-sheet-id="EL-08"');
});
