import { expect, test } from "bun:test";
import { splitCsvLine } from "../src/cli";

test("CSV parser preserves quoted commas and doubled quotes", () => {
  expect(splitCsvLine('A,"basis, with comma","a ""quote"""')).toEqual(["A", "basis, with comma", 'a "quote"']);
});
