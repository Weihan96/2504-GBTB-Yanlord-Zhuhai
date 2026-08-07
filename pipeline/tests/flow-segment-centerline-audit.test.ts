import { expect, test } from "bun:test";
import { spawnSync } from "node:child_process";

function runPython(expression: string) {
  return spawnSync("python3", ["-c", expression], {
    cwd: process.cwd(),
    encoding: "utf8",
  });
}

test("flow segment audit separates disconnected components", () => {
  const result = runPython(
    `import json,sys; sys.path.insert(0,"pipeline/scripts"); from flow_segment_centerline_audit import linked_component_groups; import numpy as np; a=np.array([[0.,0.,0.],[1.,0.,0.]]); b=np.array([[1.,0.,0.],[2.,0.,0.]]); c=np.array([[10.,0.,0.],[11.,0.,0.]]); print(json.dumps(linked_component_groups([a,b,c],0.1)))`,
  );
  expect(result.status).toBe(0);
  const [groups, pairs] = JSON.parse(result.stdout);
  expect(groups).toEqual([[0, 1], [2]]);
  expect(pairs.filter((pair: any) => pair.linked_within_tolerance)).toHaveLength(1);
});

test("straight circular prism exposes evidence but never write authority", () => {
  const result = runPython(
    `import json,sys; sys.path.insert(0,"pipeline/scripts"); from flow_segment_centerline_audit import principal_axis_record; import numpy as np; angles=np.linspace(0,2*np.pi,24,endpoint=False); p=np.array([[r*np.cos(a),r*np.sin(a),z] for z in (0.,1000.) for a in angles for r in (50.,)]); print(json.dumps(principal_axis_record(p)))`,
  );
  expect(result.status).toBe(0);
  const record = JSON.parse(result.stdout);
  expect(record.straight_circular_prism_candidate).toBe(true);
  expect(record.axis_extent_mm).toBeCloseTo(1000);
  expect(record.derived_centreline_is_write_authority).toBe(false);
});

test("flow segment audit requires a configurable positive connection tolerance", () => {
  const result = spawnSync(
    "python3",
    ["pipeline/scripts/flow_segment_centerline_audit.py", "--help"],
    { cwd: process.cwd(), encoding: "utf8" },
  );
  expect(result.status).toBe(0);
  expect(result.stdout).toContain("--connection-tolerance-mm");
});

test("flow segment family includes the IFC4 pipe subtype queue", () => {
  const source = Bun.file(
    "pipeline/scripts/flow_segment_centerline_audit.py",
  ).text();
  return expect(source).resolves.toContain('"IfcPipeSegment"');
});

test("typed multibody fittings are not mistaken for untyped pipe bundles", () => {
  const result = runPython(
    `import json,sys; sys.path.insert(0,"pipeline/scripts"); from flow_segment_centerline_audit import classify_product; print(json.dumps({"typed":classify_product(2,2,True),"untyped":classify_product(3,3,False),"controlled":classify_product(3,3,False,True)}))`,
  );
  expect(result.status).toBe(0);
  const classification = JSON.parse(result.stdout);
  expect(classification.typed).toEqual([
    "typed_multibody_product",
    "preserve_as_one_typed_product_and_audit_its_installation_or_connector_datum",
    false,
  ]);
  expect(classification.untyped[0]).toBe("disconnected_geometry_bundle");
  expect(classification.untyped[2]).toBe(true);
  expect(classification.controlled).toEqual([
    "controlled_disconnected_bundle",
    "preserve_existing_ifc_product_and_audit_each_geometry_branch_individually",
    false,
  ]);
});
