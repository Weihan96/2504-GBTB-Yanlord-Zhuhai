import { expect, test } from "bun:test";
import { readFileSync } from "node:fs";
import { resolve } from "node:path";

const targets = resolve(
  import.meta.dir,
  "../decisions/c003-construction-surface-origin-reset-targets.csv",
);

test("construction surface targets are explicit integer points on rendered faces", () => {
  const lines = readFileSync(targets, "utf8").trim().split("\n");
  const rows = lines.slice(1).map((line) => line.split(","));
  const classes = Object.fromEntries(
    ["IfcCovering", "IfcSlab"].map((ifcClass) => [
      ifcClass,
      rows.filter((row) => row[0] === ifcClass).length,
    ]),
  );
  expect(rows).toHaveLength(54);
  expect(new Set(rows.map((row) => row[1])).size).toBe(54);
  expect(classes).toEqual({ IfcCovering: 45, IfcSlab: 9 });
  expect(rows.every((row) => row.slice(2, 5).every((value) => Number.isInteger(Number(value))))).toBe(true);
  expect(new Set(rows.map((row) => row[5]))).toEqual(
    new Set(["integer_point_on_axis_aligned_surface"]),
  );
});
