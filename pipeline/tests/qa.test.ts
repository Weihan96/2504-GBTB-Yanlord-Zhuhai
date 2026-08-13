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

test("implemented Space review evidence supersedes legacy provisional ObjectType", () => {
  const snapshot = analyzeIfcText(`ISO-10303-21;
HEADER;
FILE_SCHEMA(('IFC4'));
ENDSEC;
DATA;
#1=IFCSPACE('0AAAAAAAAAAAAAAAAAAAAA',$,'SPACE-01',$,'PROVISIONAL_GRID_CELL',#2,#3,'Room',$,.INTERNAL.,$);
#4=IFCPROPERTYSINGLEVALUE('Reference',$,IFCIDENTIFIER('R01'),$);
#5=IFCPROPERTYSET('0BBBBBBBBBBBBBBBBBBBBB',$,'Pset_SpaceCommon',$,(#4));
#6=IFCRELDEFINESBYPROPERTIES('0CCCCCCCCCCCCCCCCCCCCC',$,$,$,(#1),#5);
ENDSEC;
END-ISO-10303-21;`);
  const report = buildQaReport(
    snapshot,
    [],
    { ifcSchema: "IFC4", baselineCounts: { IFCSPACE: 1 } },
    [],
    { path: "space-review.csv", total: 1, implemented: 1, ready: true, errors: [] },
  );
  expect(report.gates.find((gate) => gate.id === "SPACE-REVIEW")?.status).toBe("pass");
});

test("annotation composition and paired elevation names are baseline evidence", () => {
  const snapshot = analyzeIfcText(`ISO-10303-21;
HEADER;
FILE_SCHEMA(('IFC4'));
ENDSEC;
DATA;
#1=IFCANNOTATION('0AAAAAAAAAAAAAAAAAAAAA',$,'EL-01',$,'DRAWING',#2,#3);
#4=IFCANNOTATION('0BBBBBBBBBBBBBBBBBBBBB',$,'EL-01',$,'ELEVATION',#5,$);
ENDSEC;
END-ISO-10303-21;`);
  const report = buildQaReport(snapshot, [], {
    ifcSchema: "IFC4",
    baselineCounts: { IFCANNOTATION: 2 },
    baselineAnnotationTypes: { DRAWING: 1, ELEVATION: 1 },
    baselineElevationDrawingPairs: 1,
  });
  expect(report.gates.find((gate) => gate.id === "BASELINE-COUNTS")?.status).toBe("pass");
});

test("concentrated human review blocks QA while controlled root items remain open", () => {
  const snapshot = analyzeIfcText(`ISO-10303-21;
HEADER;
FILE_SCHEMA(('IFC4'));
ENDSEC;
DATA;
ENDSEC;
END-ISO-10303-21;`);
  const report = buildQaReport(
    snapshot,
    [],
    { ifcSchema: "IFC4", baselineCounts: {} },
    [],
    undefined,
    {
      path: "build/release/concentrated-human-review-current.json",
      sourceIfcSha256: snapshot.source.sha256,
      reviewItemCount: 7,
      openRootReviewItemCount: 2,
      unmappedBlockerCount: 0,
      mappingComplete: true,
      constructionReleaseReady: false,
      errors: [],
    },
  );
  const gate = report.gates.find((item) => item.id === "CONCENTRATED-HUMAN-REVIEW");
  expect(gate?.status).toBe("block");
  expect(gate?.summary).toContain("2 concentrated root review item(s) remain open");
});

test("concentrated human review fails closed on unmapped blockers", () => {
  const snapshot = analyzeIfcText(`ISO-10303-21;
HEADER;
FILE_SCHEMA(('IFC4'));
ENDSEC;
DATA;
ENDSEC;
END-ISO-10303-21;`);
  const report = buildQaReport(
    snapshot,
    [],
    { ifcSchema: "IFC4", baselineCounts: {} },
    [],
    undefined,
    {
      path: "build/release/concentrated-human-review-current.json",
      sourceIfcSha256: snapshot.source.sha256,
      reviewItemCount: 1,
      openRootReviewItemCount: 1,
      unmappedBlockerCount: 1,
      mappingComplete: false,
      constructionReleaseReady: false,
      errors: [],
    },
  );
  const gate = report.gates.find((item) => item.id === "CONCENTRATED-HUMAN-REVIEW");
  expect(gate?.status).toBe("block");
  expect(gate?.summary).toContain("1 release blocker(s) lack a controlled root review item");
});

test("concentrated human review passes when mapping is complete and no root items remain", () => {
  const snapshot = analyzeIfcText(`ISO-10303-21;
HEADER;
FILE_SCHEMA(('IFC4'));
ENDSEC;
DATA;
ENDSEC;
END-ISO-10303-21;`);
  const report = buildQaReport(
    snapshot,
    [],
    { ifcSchema: "IFC4", baselineCounts: {} },
    [],
    undefined,
    {
      path: "build/release/concentrated-human-review-current.json",
      sourceIfcSha256: snapshot.source.sha256,
      reviewItemCount: 0,
      openRootReviewItemCount: 0,
      unmappedBlockerCount: 0,
      mappingComplete: true,
      constructionReleaseReady: false,
      errors: [],
    },
  );
  expect(report.gates.find((item) => item.id === "CONCENTRATED-HUMAN-REVIEW")?.status).toBe("pass");
});
