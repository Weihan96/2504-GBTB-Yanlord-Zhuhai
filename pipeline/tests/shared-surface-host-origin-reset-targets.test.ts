import { expect, test } from "bun:test";
import { readFileSync } from "node:fs";
import { resolve } from "node:path";

import { splitCsvLine } from "../src/cli";

test("shared surface host targets are the exact five coupled Coverings", () => {
  const path = resolve(
    process.cwd(),
    "pipeline/decisions/c003-shared-surface-host-origin-reset-targets.csv",
  );
  const [header, ...values] = readFileSync(path, "utf8")
    .trim()
    .split(/\r?\n/)
    .map(splitCsvLine);
  const rows = values.map((value) =>
    Object.fromEntries(header.map((key, index) => [key, value[index]])),
  );
  expect(rows).toHaveLength(5);
  expect(new Set(rows.map((row) => row.global_id))).toEqual(
    new Set([
      "3UpLXAml9Bufs6wAjUjPJT",
      "06wFwLoDD6ie5iCTnc_yad",
      "0rwOcjvjn5WBdoq4bxiHSk",
      "2roMaBy7v9OvqUBSpWl2Q2",
      "37DESmCqj8qAHVsS$K2iyb",
    ]),
  );
  expect(rows.every((row) => row.expected_class === "IfcCovering")).toBe(true);
  expect(
    rows.every((row) =>
      ["target_x_mm", "target_y_mm", "target_z_mm"].every(
        (key) => Number.isInteger(Number(row[key])),
      ),
    ),
  ).toBe(true);
  expect(
    rows.every(
      (row) => row.anchor_kind === "integer_point_on_axis_aligned_surface",
    ),
  ).toBe(true);
});
