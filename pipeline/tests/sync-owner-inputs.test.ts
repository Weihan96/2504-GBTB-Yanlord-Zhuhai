import { expect, test } from "bun:test";
import { existsSync, mkdtempSync, readFileSync, writeFileSync } from "node:fs";
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
  const workbook = join(temp, "input.xlsx");
  const script = `
from pathlib import Path
import csv, zipfile
from xml.sax.saxutils import escape

root = Path(${JSON.stringify(root)})
rows = list(csv.reader(Path(${JSON.stringify(input)}).open(encoding="utf-8-sig")))
decision_rows = list(csv.reader((root / "pipeline/decisions/owner-input-register.csv").open(encoding="utf-8-sig")))
display = {
  "decisions": ${JSON.stringify([
    "ID", "专业", "优先级", "阻塞无保留发布", "需要你确认", "常见候选（不等于确认）",
    "你的确认值", "单位", "状态", "证据/链接", "现有依据", "同步目标", "备注",
  ])},
  "appliances": ${JSON.stringify([
    "设备ID", "设备名称", "类别", "候选存放位置", "候选使用位置", "确认存放位置",
    "确认使用位置", "数量", "额定功率(W)", "同时使用组", "需给水", "需排水",
    "需燃气", "需通风", "型号", "证据/链接", "状态", "备注",
  ])},
}
def sheet_xml(data, header):
    values = [header] + data[1:]
    body = []
    for r, row in enumerate(values, 1):
        cells = []
        for c, value in enumerate(row, 1):
            n, letters = c, ""
            while n:
                n, rem = divmod(n - 1, 26); letters = chr(65 + rem) + letters
            cells.append(f'<c r="{letters}{r}" t="inlineStr"><is><t>{escape(value)}</t></is></c>')
        body.append(f'<row r="{r}">' + ''.join(cells) + '</row>')
    return '<?xml version="1.0" encoding="UTF-8" standalone="yes"?><worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main"><sheetData>' + ''.join(body) + '</sheetData></worksheet>'
with zipfile.ZipFile(Path(${JSON.stringify(workbook)}), "w") as z:
    z.writestr("xl/workbook.xml", '<?xml version="1.0"?><workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships"><sheets><sheet name="设计决策" sheetId="1" r:id="rId1"/><sheet name="家电清单" sheetId="2" r:id="rId2"/></sheets></workbook>')
    z.writestr("xl/_rels/workbook.xml.rels", '<?xml version="1.0"?><Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships"><Relationship Id="rId1" Target="worksheets/sheet1.xml"/><Relationship Id="rId2" Target="worksheets/sheet2.xml"/></Relationships>')
    z.writestr("xl/worksheets/sheet1.xml", sheet_xml(decision_rows, display["decisions"]))
    z.writestr("xl/worksheets/sheet2.xml", sheet_xml(rows, display["appliances"]))
`;
  const built = Bun.spawnSync(["python3", "-c", script], { cwd: root });
  expect(built.exitCode).toBe(0);
  const rejected = Bun.spawnSync([
    "python3", "pipeline/scripts/sync_owner_inputs.py", "--input", workbook,
  ], { cwd: root });
  expect(rejected.exitCode).toBe(2);
  expect(rejected.stderr.toString()).toContain("zero is allowed only for a gas appliance");
});
