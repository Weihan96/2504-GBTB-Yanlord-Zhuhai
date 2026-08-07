import { expect, test } from "bun:test";
import { analyzeIfcText } from "../src/ifc-step";
import { buildQaReport } from "../src/qa";
import { DECISION_STATUSES } from "../src/cli";

test("delegated is a valid decision status for a scoped downstream handoff", () => {
  expect(DECISION_STATUSES.has("delegated")).toBe(true);
});

test("planned is a valid decision status for recorded future work", () => {
  expect(DECISION_STATUSES.has("planned")).toBe(true);
});

test("release gate blocks known incomplete IFC data", () => {
  const snapshot = analyzeIfcText(`ISO-10303-21;
HEADER;
FILE_SCHEMA(('IFC4'));
ENDSEC;
DATA;
#1=IFCSPACE('0AAAAAAAAAAAAAAAAAAAAA',$,'SPACE-01',$,'PROVISIONAL_GRID_CELL',#2,#3,'Room',$,.INTERNAL.,$);
ENDSEC;
END-ISO-10303-21;`);
  const report = buildQaReport(snapshot, [], { ifcSchema: "IFC4", baselineCounts: { IFCSPACE: 1 } });
  expect(report.summary.releasable).toBe(false);
  expect(report.summary.block).toBeGreaterThan(0);
  expect(report.gates.find((gate) => gate.id === "SPACE-REVIEW")?.status).toBe("block");
});
