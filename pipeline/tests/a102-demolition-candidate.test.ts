import { expect, test } from "bun:test";
import { resolve } from "node:path";

const root = resolve(import.meta.dir, "../..");
const modulePath = resolve(root, "pipeline/scripts/a102_demolition_candidate.py");

function runPython(body: string) {
  return Bun.spawnSync(["python3", "-c", `
import importlib.util, json, pathlib, sys
spec = importlib.util.spec_from_file_location("a102_demolition_candidate", ${JSON.stringify(modulePath)})
module = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = module
spec.loader.exec_module(module)
${body}
`], { cwd: root });
}

test("A102 source boundary contains 11 already-removed and 2 planned segments", () => {
  const result = runPython(`
from collections import Counter
print(json.dumps({
  "count": len(module.CANONICAL_RECTS),
  "status": Counter(rect[1] for rect in module.CANONICAL_RECTS),
  "nominal": len(module.NOMINAL_BBOXES),
  "confidence": len(module.CONFIDENCE),
}))
`);
  expect(result.exitCode).toBe(0);
  expect(JSON.parse(result.stdout.toString())).toEqual({
    count: 13,
    status: { ALREADY_REMOVED: 11, PLANNED_DEMOLITION: 2 },
    nominal: 13,
    confidence: 13,
  });
});

test("A102 grid fit stays within 0.5 mm and candidates stay on 50 mm grid", () => {
  const result = runPython(`
sx, bx, xr = module.linear_fit(module.X_ANCHORS)
sy, by, yr = module.linear_fit(module.Y_ANCHORS)
records = module.build_records(module.CANONICAL_RECTS, (sx, bx), (sy, by), "pdf", "dwg")
print(json.dumps({
  "residual": max(max(abs(value) for value in xr), max(abs(value) for value in yr)),
  "grid": max(abs(float(record[key]) / 50 - round(float(record[key]) / 50)) for record in records for key in (
    "candidate_x_min_mm", "candidate_y_min_mm", "candidate_x_max_mm", "candidate_y_max_mm"
  )),
  "review": all(record["review_required"] == "yes" and record["formal_ifc_write_allowed"] == "no" for record in records),
}))
`);
  expect(result.exitCode).toBe(0);
  const value = JSON.parse(result.stdout.toString());
  expect(value.residual).toBeLessThanOrEqual(0.5);
  expect(value.grid).toBe(0);
  expect(value.review).toBe(true);
});

test("A102 protected adjusted door opening is explicit in D13 evidence", () => {
  const result = runPython(`
sx, bx, _ = module.linear_fit(module.X_ANCHORS)
sy, by, _ = module.linear_fit(module.Y_ANCHORS)
record = module.build_records(module.CANONICAL_RECTS, (sx, bx), (sy, by), "pdf", "dwg")[-1]
print(json.dumps(record))
`);
  expect(result.exitCode).toBe(0);
  const record = JSON.parse(result.stdout.toString());
  expect(record.candidate_id).toBe("D13");
  expect(record.basis).toContain("0hKdvAZkn1TejLgJhK_vDp");
  expect(record.basis).toContain("1YxMx6s0r3ZPPohkRKXWbl");
  expect(record.review_status).toBe("pending");
});
