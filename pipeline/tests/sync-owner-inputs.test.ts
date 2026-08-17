import { expect, test } from "bun:test";
import { copyFileSync, existsSync, mkdtempSync, readFileSync, writeFileSync } from "node:fs";
import { tmpdir } from "node:os";
import { join, resolve } from "node:path";

const root = resolve(import.meta.dir, "../..");

test("owner input sync is dry-run by default and never writes IFC", () => {
  const temp = mkdtempSync(join(tmpdir(), "owner-inputs-"));
  const report = join(temp, "report.json");
  const openItems = join(temp, "open.md");
  const decisions = join(root, "pipeline/decisions/owner-input-register.csv");
  const before = readFileSync(decisions, "utf8");
  const result = Bun.spawnSync([
    "python3", "pipeline/scripts/sync_owner_inputs.py",
    "--input", decisions,
    "--report", report,
    "--open-items", openItems,
  ], { cwd: root });
  expect(result.exitCode).toBe(0);
  expect(readFileSync(decisions, "utf8")).toBe(before);
  const parsed = JSON.parse(readFileSync(report, "utf8"));
  expect(parsed.mode).toBe("dry-run");
  expect(parsed.formal_ifc_write).toBe(false);
  expect(parsed.summary.release_blocking_open_count).toBeGreaterThan(0);
  expect(parsed.summary.automatic_close_open_count).toBe(0);
  expect(parsed.summary.human_or_external_closeout_open_count).toBe(
    parsed.summary.release_blocking_open_count,
  );
  expect(readFileSync(openItems, "utf8")).toContain("必须由人审、现场、厂家或主管方关闭");
  const entry = parsed.normalized_inputs.decisions.find(
    (row: { input_id: string }) => row.input_id === "E302-ENTRY-SIDE",
  );
  expect(entry.candidate_value).not.toBe("");
  expect(entry.effective_value).toBe("");
  const adopted = parsed.normalized_inputs.decisions.find(
    (row: { input_id: string }) => row.input_id === "E303-NS01-FORM",
  );
  expect(adopted).toMatchObject({
    decision_status: "采用候选",
    closeout_status: "decision_confirmed_evidence_pending",
  });
  expect(parsed.summary.release_blocking_open_ids).toContain("E303-NS01-FORM");
  const wallTags = parsed.normalized_inputs.decisions.find(
    (row: { input_id: string }) => row.input_id === "A103-WALL-TAG-SCHEME",
  );
  expect(wallTags.closeout_status).toBe("closed_by_human_confirmation");
});

test("default dry-run writes only stdout", () => {
  const temp = mkdtempSync(join(tmpdir(), "owner-inputs-stdout-"));
  const result = Bun.spawnSync([
    "python3", join(root, "pipeline/scripts/sync_owner_inputs.py"),
  ], { cwd: temp });
  expect(result.exitCode).toBe(0);
  expect(JSON.parse(result.stdout.toString()).mode).toBe("dry-run");
  expect(existsSync(join(temp, "build"))).toBe(false);
});

test("owner input sync rejects a custom confirmation without a value", () => {
  const temp = mkdtempSync(join(tmpdir(), "owner-inputs-invalid-"));
  const source = readFileSync(join(root, "pipeline/decisions/owner-input-register.csv"), "utf8");
  const invalid = source.replace(
    ",Entry A,,待填写,,门口控制应位于真实墙面",
    ",,,自定义确认,,门口控制应位于真实墙面",
  );
  expect(invalid).not.toBe(source);
  const input = join(temp, "invalid.csv");
  writeFileSync(input, invalid);
  const result = Bun.spawnSync([
    "python3", "pipeline/scripts/sync_owner_inputs.py",
    "--input", input,
    "--report", join(temp, "report.json"),
    "--open-items", join(temp, "open.md"),
  ], { cwd: root });
  expect(result.exitCode).toBe(2);
  expect(result.stderr.toString()).toContain("custom confirmation requires user_value");
});

test("owner input sync rejects workbook-owned candidate changes", () => {
  const temp = mkdtempSync(join(tmpdir(), "owner-inputs-protected-"));
  const source = readFileSync(join(root, "pipeline/decisions/owner-input-register.csv"), "utf8");
  const invalid = source.replace(
    "Entry A（靠进入后顺手侧且避开门扇/门套）",
    "Entry B",
  );
  const input = join(temp, "invalid.csv");
  writeFileSync(input, invalid);
  const result = Bun.spawnSync([
    "python3", "pipeline/scripts/sync_owner_inputs.py", "--input", input,
  ], { cwd: root });
  expect(result.exitCode).toBe(2);
  expect(result.stderr.toString()).toContain("protected fields changed: candidate_value");
});

test("owner input sync rejects incomplete closeout ownership coverage", () => {
  const temp = mkdtempSync(join(tmpdir(), "owner-inputs-closeout-"));
  const source = readFileSync(
    join(root, "pipeline/decisions/owner-input-closeout-rules.csv"),
    "utf8",
  );
  const incomplete = source
    .split("\n")
    .filter((line) => !line.startsWith("E304-AP-POWER,"))
    .join("\n");
  const closeout = join(temp, "closeout.csv");
  writeFileSync(closeout, incomplete);
  const result = Bun.spawnSync([
    "python3", "pipeline/scripts/sync_owner_inputs.py",
    "--closeout-rules", closeout,
  ], { cwd: root });
  expect(result.exitCode).toBe(2);
  expect(result.stderr.toString()).toContain("closeout rules: IDs changed");
  expect(result.stderr.toString()).toContain("E304-AP-POWER");
});

test("pending appliance fields never become effective design inputs", () => {
  const result = Bun.spawnSync([
    "python3", "pipeline/scripts/sync_owner_inputs.py",
  ], { cwd: root });
  expect(result.exitCode).toBe(0);
  const parsed = JSON.parse(result.stdout.toString());
  const pending = parsed.normalized_inputs.appliances.find(
    (row: { status: string }) => row.status === "待填写",
  );
  expect(pending.candidate_use_location).not.toBe("");
  expect(Object.values(pending.effective).every((value) => value === "")).toBe(true);
});

test("zero electrical load is accepted only for gas appliances", () => {
  const accepted = Bun.spawnSync([
    "python3", "pipeline/scripts/sync_owner_inputs.py",
    "--input", "pipeline/decisions/owner-input-register.csv",
  ], { cwd: root });
  expect(accepted.exitCode).toBe(0);

  const temp = mkdtempSync(join(tmpdir(), "owner-inputs-zero-load-"));
  const source = readFileSync(join(root, "pipeline/decisions/appliance-input-register.csv"), "utf8");
  const invalid = source.replace(
    "APP-011,烤箱,固定厨房电器,西厨高柜,西厨高柜,西厨高柜,西厨高柜,1,3400,",
    "APP-011,烤箱,固定厨房电器,西厨高柜,西厨高柜,西厨高柜,西厨高柜,1,0,",
  );
  const input = join(temp, "invalid-appliances.csv");
  writeFileSync(input, invalid);
  const rejected = Bun.spawnSync([
    "python3", "pipeline/scripts/sync_owner_inputs.py",
    "--input", "pipeline/decisions/owner-input-register.csv",
    "--appliances", input,
  ], { cwd: root });
  expect(rejected.exitCode).toBe(2);
  expect(rejected.stderr.toString()).toContain("zero is allowed only for a gas appliance");
});

test("readonly SSOT sheets report full parity and apply rebuilds them without CSV writeback", () => {
  const temp = mkdtempSync(join(tmpdir(), "owner-inputs-parity-"));
  const workbook = join(temp, "owner-inputs.xlsx");
  const decisions = join(temp, "owner-input-register.csv");
  const appliances = join(temp, "appliance-input-register.csv");
  const closeout = join(temp, "owner-input-closeout-rules.csv");
  const pm = join(temp, "pm.md");
  const equipment = join(temp, "equipment-register.csv");
  const requirements = join(temp, "equipment-installation-requirements.csv");
  const evidence = join(temp, "source-evidence-register.csv");
  for (const [source, target] of [
    ["output/forms/滨海湾施工输入清单.xlsx", workbook],
    ["pipeline/decisions/owner-input-register.csv", decisions],
    ["pipeline/decisions/appliance-input-register.csv", appliances],
    ["pipeline/decisions/owner-input-closeout-rules.csv", closeout],
    ["drawings/滨海湾装修施工图深化工作管理.md", pm],
    ["pipeline/decisions/equipment-register.csv", equipment],
    ["pipeline/decisions/equipment-installation-requirements.csv", requirements],
    ["pipeline/decisions/source-evidence-register.csv", evidence],
  ] as const) {
    copyFileSync(join(root, source), target);
  }

  const equipmentText = readFileSync(equipment, "utf8");
  const changedEquipment = equipmentText.replace(
    "APP-001,APPLIANCE,移动厨电,火锅电器",
    "APP-001,APPLIANCE,移动厨电,火锅设备（parity test）",
  );
  expect(changedEquipment).not.toBe(equipmentText);
  writeFileSync(equipment, changedEquipment);
  writeFileSync(
    requirements,
    readFileSync(requirements, "utf8").trimEnd()
      + "\nREQ-PARITY-TEST,APP-001,ELEC,parity_test,yes,,,,project_candidate,candidate,,,yes,test only\n",
  );
  const evidenceText = readFileSync(evidence, "utf8");
  const changedEvidence = evidenceText.replace(
    "强弱电箱位于过道高柜",
    "强弱电箱位于过道高柜（parity test）",
  );
  expect(changedEvidence).not.toBe(evidenceText);
  writeFileSync(evidence, changedEvidence);

  const args = [
    "--input", workbook,
    "--decisions", decisions,
    "--appliances", appliances,
    "--closeout-rules", closeout,
    "--pm", pm,
    "--equipment-register", equipment,
    "--installation-requirements", requirements,
    "--evidence-register", evidence,
  ];
  const dryRun = Bun.spawnSync([
    "python3", "pipeline/scripts/sync_owner_inputs.py", ...args,
  ], { cwd: root });
  expect(dryRun.exitCode, dryRun.stderr.toString()).toBe(0);
  const before = JSON.parse(dryRun.stdout.toString());
  expect(before.readonly_view_parity.workbook_structure_match).toBe(true);
  expect(before.readonly_view_parity.views["设备主表"]).toMatchObject({
    row_count_match: true,
    columns_match: true,
    keys_match: true,
    values_match: false,
    all_match: false,
  });
  expect(before.readonly_view_parity.views["安装条件"]).toMatchObject({
    row_count_match: false,
    columns_match: true,
    keys_match: false,
    all_match: false,
  });
  expect(before.readonly_view_parity.views["证据索引"]).toMatchObject({
    columns_match: true,
    values_match: false,
    all_match: false,
  });

  const canonicalBefore = new Map([
    [equipment, readFileSync(equipment)],
    [requirements, readFileSync(requirements)],
    [evidence, readFileSync(evidence)],
  ]);
  const applied = Bun.spawnSync([
    "python3", "pipeline/scripts/sync_owner_inputs.py", ...args, "--apply",
  ], { cwd: root });
  expect(applied.exitCode, applied.stderr.toString()).toBe(0);
  const after = JSON.parse(applied.stdout.toString());
  expect(after.readonly_view_parity.all_match).toBe(false);
  expect(after.readonly_view_parity_after_apply.all_match).toBe(true);
  const tableAudit = Bun.spawnSync([
    "python3", "-c", String.raw`
import json,sys,zipfile
from xml.etree import ElementTree as ET
ns="http://schemas.openxmlformats.org/spreadsheetml/2006/main"
with zipfile.ZipFile(sys.argv[1]) as archive:
    assert archive.testzip() is None
    result={}
    for name in archive.namelist():
        if not name.startswith("xl/tables/table") or not name.endswith(".xml"):
            continue
        root=ET.fromstring(archive.read(name))
        result[root.attrib["name"]]=root.attrib["ref"]
print(json.dumps(result,ensure_ascii=False))
`, workbook,
  ], { cwd: root });
  expect(tableAudit.exitCode, tableAudit.stderr.toString()).toBe(0);
  const tableRefs = JSON.parse(tableAudit.stdout.toString());
  expect(tableRefs["OwnerDecisionInputs"]).toBe("A1:M41");
  expect(tableRefs["设备主表Table"]).toBe("A2:AB165");
  expect(tableRefs["安装条件Table"]).toBe("A2:N852");
  expect(tableRefs["证据索引Table"]).toBe("A2:X190");
  const summaryAudit = Bun.spawnSync([
    "python3", "-c", String.raw`
import json,sys,zipfile
from xml.etree import ElementTree as ET
ns="http://schemas.openxmlformats.org/spreadsheetml/2006/main"
with zipfile.ZipFile(sys.argv[1]) as archive:
    root=ET.fromstring(archive.read("xl/worksheets/sheet1.xml"))
    result={}
    for reference in ("B13","B14","B15","B16"):
        cell=root.find(f".//{{{ns}}}c[@r='{reference}']")
        value=cell.find(f"{{{ns}}}v") if cell is not None else None
        result[reference]=value.text if value is not None else None
print(json.dumps(result,ensure_ascii=False))
`, workbook,
  ], { cwd: root });
  expect(summaryAudit.exitCode, summaryAudit.stderr.toString()).toBe(0);
  expect(JSON.parse(summaryAudit.stdout.toString())).toEqual({
    B13: "40",
    B14: "33",
    B15: "19",
    B16: "3",
  });
  for (const [path, contents] of canonicalBefore) {
    expect(readFileSync(path)).toEqual(contents);
  }

  const verified = Bun.spawnSync([
    "python3", "pipeline/scripts/sync_owner_inputs.py", ...args,
  ], { cwd: root });
  expect(verified.exitCode, verified.stderr.toString()).toBe(0);
  const verifiedReport = JSON.parse(verified.stdout.toString());
  expect(verifiedReport.readonly_view_parity.all_match).toBe(true);
  for (const view of Object.values(verifiedReport.readonly_view_parity.views) as Array<Record<string, unknown>>) {
    expect(view).toMatchObject({
      row_count_match: true,
      columns_match: true,
      keys_match: true,
      values_match: true,
      all_match: true,
    });
  }
}, 30_000);
