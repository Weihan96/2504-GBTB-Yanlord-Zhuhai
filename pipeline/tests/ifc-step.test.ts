import { describe, expect, test } from "bun:test";
import { analyzeIfcText, decodeStepString, splitStepArgs } from "../src/ifc-step";

const SAMPLE_IFC = `ISO-10303-21;
HEADER;
FILE_NAME('sample.ifc','2026-08-04T00:00:00+08:00',('Author'),('Org'),'IfcOpenShell','Bonsai','');
FILE_SCHEMA(('IFC4'));
ENDSEC;
DATA;
#1=IFCPROJECT('0AAAAAAAAAAAAAAAAAAAAA',$,'Project',$,$,$,$,$,$);
#2=IFCBUILDINGSTOREY('0BBBBBBBBBBBBBBBBBBBBB',$,'FFL',$,$,$,$,$,$,$);
#3=IFCSPACE('0CCCCCCCCCCCCCCCCCCCCC',$,'SPACE-01','Description','PROVISIONAL_GRID_CELL',#10,#11,'\\X2\\4E3B5367\\X0\\',$,.INTERNAL.,$);
#4=IFCDOOR('0DDDDDDDDDDDDDDDDDDDDD',$,'Door',$,$,#12,#13,$,2000.,800.,.DOOR.,.SINGLE_SWING_LEFT.,$);
#5=IFCWINDOW('0EEEEEEEEEEEEEEEEEEEEE',$,'Window',$,$,#14,#15,'W01',1800.,1200.,.WINDOW.,.SINGLE_PANEL.,$);
#6=IFCANNOTATION('0FFFFFFFFFFFFFFFFFFFFF',$,'Wall Plan',$,'DRAWING',#16,#17);
#7=IFCGRIDAXIS('08',#18,.T.);
#8=IFCRELCONTAINEDINSPATIALSTRUCTURE('0GGGGGGGGGGGGGGGGGGGGG',$,$,$,(#3),#2);
#9=IFCRELFILLSELEMENT('0HHHHHHHHHHHHHHHHHHHHH',$,$,$,#20,#4);
ENDSEC;
END-ISO-10303-21;`;

describe("STEP parser", () => {
  test("splits nested arguments without breaking quoted commas", () => {
    expect(splitStepArgs("'A,B',(#1,#2),.T.")).toEqual(["'A,B'", "(#1,#2)", ".T."]);
  });

  test("decodes IFC X2 unicode strings", () => {
    expect(decodeStepString("'\\X2\\4E3B5367\\X0\\'")).toBe("主卧");
  });

  test("extracts P0 identity, relationships, and opening data", () => {
    const snapshot = analyzeIfcText(SAMPLE_IFC);
    expect(snapshot.schema).toBe("IFC4");
    expect(snapshot.projectName).toBe("Project");
    expect(snapshot.spaces.total).toBe(1);
    expect(snapshot.spaces.records[0].longName).toBe("主卧");
    expect(snapshot.spaces.provisional).toBe(1);
    expect(snapshot.spaces.containedInStorey).toBe(1);
    expect(snapshot.spaces.aggregatedUnderStorey).toBe(0);
    expect(snapshot.doors.records[0].fillsOpening).toBe(true);
    expect(snapshot.doors.missingTag).toBe(1);
    expect(snapshot.windows.missingTag).toBe(0);
    expect(snapshot.annotations.byObjectType.DRAWING).toBe(1);
    expect(snapshot.grids.axes).toEqual(["08"]);
    expect(snapshot.integrity.duplicateStepIds).toEqual([]);
    expect(snapshot.integrity.duplicateGlobalIds).toEqual([]);
    expect(snapshot.integrity.hasValidFooter).toBe(true);
  });
});
