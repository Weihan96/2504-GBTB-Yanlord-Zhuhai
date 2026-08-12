import { expect, test } from "bun:test";
import { readFileSync } from "node:fs";

const script = readFileSync(
  "pipeline/scripts/int1_lightweight_elevation_representation.py",
  "utf8",
);

test("lightweight elevations are curve-only and preserve detailed bodies", () => {
  expect(script).toContain('"Body", "ELEVATION_VIEW"');
  expect(script).toContain('"Body", "Curve3D"');
  expect(script).toContain("createIfcPolyline");
  expect(script).not.toContain("createIfcCircle");
  expect(script).not.toContain("createIfcFace");
  expect(script).not.toContain("createIfcExtrudedAreaSolid");
  expect(script).toContain(
    "definition.Representations = tuple(definition.Representations) + (representation,)",
  );
});

test("all twelve controlled exclusions receive a lightweight profile", () => {
  const ids = [
    "0V9CnYT0n3GvkLmXhxvv$u",
    "0YgqsahOH4FvZ_nLwl_Tlf",
    "2wViHGrwX2hOYCb3MDNb$9",
    "2S2c498tb7$gzdukjhCGVQ",
    "2iKOL78$H0N9Yd9$ky3pW4",
    "3hgNkx97vCTOC2eewCpMNk",
    "3pvAlH5C14v8uVEJ1LmK8M",
    "2xmcLzu1rDTeMzRuNxPDyE",
    "3IQBEqO5vDI8Z9k1Ltge_N",
    "1MzM8Ms2vFo8KEm503j9w2",
    "1O9JRXCI56VRUbpuLJy86Z",
    "1i_pqgLv9A7uuV7MjaArBW",
  ];
  for (const id of ids) expect(script).toContain(id);
});

test("Bonsai elevation compiler prioritises lightweight target-view bodies", () => {
  const compiler = readFileSync(
    "pipeline/scripts/int1_bonsai_elevation.py",
    "utf8",
  );
  expect(compiler).toContain(
    'element, "Model", "Body", "ELEVATION_VIEW"',
  );
  expect(compiler).toContain("lightweight_elevation_global_ids");
  expect(compiler).toContain("INT1_BONSAI_OUTPUT_DIR");
});
