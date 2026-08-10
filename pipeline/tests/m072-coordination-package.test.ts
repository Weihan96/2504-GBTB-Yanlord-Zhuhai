import { expect, test } from "bun:test";

const source = await Bun.file("pipeline/scripts/m072_coordination_package.py").text();

test("M072 package requires the complete drawing register and same-IFC candidates", () => {
  expect(source).toContain('"A-001", "A-101", "A-102", "A-103", "A-104", "A-105"');
  expect(source).toContain('"source_report_same_ifc_count"');
  expect(source).toContain('"legacy_blocks_disclosed"');
  expect(source).toContain('"candidate_pdf_pass_count"');
  expect(source).toContain('"candidate_output_pass_count"');
  expect(source).toContain('EXPECTED_PLANNED_SHEETS = ("A-101",)');
  expect(source).toContain('"m072_mechanical_pass"');
});

test("M072 package treats delegated work as disclosed rather than resolved", () => {
  expect(source).toContain('"delegated": "已移交"');
  expect(source).toContain('“已移交/已计划”不等于“已解决”');
});

test("M072 language covers mixed candidate outputs and the 2026-08-19 milestone", () => {
  expect(source).toContain("# 滨海湾 2026-08-19 待复核施工候选包");
  expect(source).toContain('[["图号", "图名", "输出文件", "状态"]]');
  expect(source).toContain("任一已登记候选图不能从当前 IFC 与受控决策证据重现");
  expect(source).not.toContain("# 滨海湾 M072 72小时协调候选包");
});

test("M072 PDF audit checks page size and rendered ink", () => {
  expect(source).toContain('abs(info["width_mm"] - 500.0) <= 0.2');
  expect(source).toContain('raster["ink_ratio"] >= 0.002');
  expect(source).toContain('raster["ink_bbox_page_ratio"] >= 0.35');
  expect(source).toContain('"---" for _ in rows[0]');
});

test("M072 candidate output audit accepts only validated PDF and SVG artifacts", () => {
  expect(source).toContain('target.suffix.lower() == ".pdf"');
  expect(source).toContain('target.suffix.lower() == ".svg"');
  expect(source).toContain('record["byte_count"] >= 1000');
  expect(source).toContain('record["viewbox_present"]');
});
