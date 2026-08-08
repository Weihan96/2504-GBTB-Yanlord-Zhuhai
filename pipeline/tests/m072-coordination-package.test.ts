import { expect, test } from "bun:test";

const source = await Bun.file("pipeline/scripts/m072_coordination_package.py").text();

test("M072 package requires the complete drawing register and same-IFC candidates", () => {
  expect(source).toContain('"A-001", "A-101", "A-102", "A-103", "A-104", "A-105"');
  expect(source).toContain('"source_report_same_ifc_count"');
  expect(source).toContain('"legacy_blocks_disclosed"');
  expect(source).toContain('"candidate_pdf_pass_count"');
  expect(source).toContain('"m072_mechanical_pass"');
});

test("M072 package treats delegated work as disclosed rather than resolved", () => {
  expect(source).toContain('"delegated": "已移交"');
  expect(source).toContain('“已移交/已计划”不等于“已解决”');
});

test("M072 PDF audit checks page size and rendered ink", () => {
  expect(source).toContain('abs(info["width_mm"] - 500.0) <= 0.2');
  expect(source).toContain('raster["ink_ratio"] >= 0.002');
  expect(source).toContain('raster["ink_bbox_page_ratio"] >= 0.35');
  expect(source).toContain('"---" for _ in rows[0]');
});
