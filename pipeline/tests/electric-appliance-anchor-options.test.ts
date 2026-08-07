import { expect, test } from "bun:test";
import { resolve } from "node:path";

const root = resolve(import.meta.dir, "../..");

function runPython(source: string): string {
  const result = Bun.spawnSync(["python3", "-c", source], { cwd: root });
  expect(result.exitCode).toBe(0);
  return result.stdout.toString().trim();
}

test("electric appliance option scope is an exact four-object set", () => {
  const output = runPython(`
import json,sys
sys.path.insert(0,"pipeline/scripts")
from electric_appliance_anchor_options import TARGET_TYPES
print(json.dumps(TARGET_TYPES,sort_keys=True))
`);
  expect(JSON.parse(output)).toEqual({
    "0UOnmuAdP1MPy6p3olwiEU": "WD01",
    "1PUCikoaP5fgiYt8sJd8$6": "AC700",
    "288GLY62v8kPPydA1lAK8W": "HD01",
    "3PQOXKxgj6IftqWXFXMQXG": "OV01",
  });
});

test("bbox feature evidence identifies an edge midpoint but rejects an arbitrary origin", () => {
  const output = runPython(`
import json,sys,numpy as np
sys.path.insert(0,"pipeline/scripts")
from electric_appliance_anchor_options import current_origin_feature
print(json.dumps({
  "edge": current_origin_feature(np.array([-447.5,-345.0,0.0]),np.array([447.5,0.0,873.0]),0.01),
  "arbitrary": current_origin_feature(np.array([-410.9,-238.9,-190.5]),np.array([436.9,238.9,1.5]),0.01),
}))
`);
  expect(JSON.parse(output)).toEqual({
    edge: "bbox_edge_midpoint",
    arbitrary: null,
  });
});

test("option reports serialize numpy evidence without losing numeric values", () => {
  const output = runPython(`
import json,sys,numpy as np
sys.path.insert(0,"pipeline/scripts")
from electric_appliance_anchor_options import json_default
print(json.dumps({"array":np.array([1.0,2.0]),"scalar":np.float64(0.1)},default=json_default))
`);
  expect(JSON.parse(output)).toEqual({ array: [1, 2], scalar: 0.1 });
});
