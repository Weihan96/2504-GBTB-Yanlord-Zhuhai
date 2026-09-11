import { expect, test } from 'bun:test';
import { readFileSync, existsSync } from 'node:fs';
import { createHash } from 'node:crypto';
import { join, resolve } from 'node:path';

const root = resolve(import.meta.dir, '../..');
const dir = join(root, 'output/review/highpoly-types/bst02');
const read = (path: string) => JSON.parse(readFileSync(path, 'utf8'));
const hash = (path: string) => createHash('sha256').update(readFileSync(path)).digest('hex');

test('BST02 classifies actual connected components and keeps the approved input geometry', () => {
  const audit = read(join(dir, 'bst02-semantic-segmentation.json'));
  expect(audit.representation_item_count).toBe(2);
  expect(audit.connected_component_count).toBe(3);
  expect(audit.components.map((c: any) => c.semantic).sort()).toEqual([
    'cylindrical_cabinet_body', 'stepped_top_cap', 'vertical_handle',
  ]);
  expect(audit.mesh_face_count).toBe(1418);
  expect(audit.actual_height_mm).toBe(250);
  expect(audit.height_discrepancy_mm).toBe(100);
  expect(audit.unsupported_invented_features).toEqual([]);
  expect(hash(join(root, audit.isolated_ifc))).toBe(audit.isolated_ifc_sha256);
  expect(hash(join(root, '2504 GBTB Yanlord Zhuhai.ifc'))).toBe('7a50b87e8f48a7c2bfbaa9f1dabdfca7155fa684b2325c66aa1ae15c8ab0a25c');
});

test('top view hides triangulation; elevations keep the actual recessed cap seam', () => {
  const audit = read(join(dir, 'bst02-semantic-segmentation.json'));
  const plan = audit.views.plan;
  expect(plan).toHaveLength(2);
  expect(plan.filter((p: any) => p.component === 'cylindrical_cabinet_body')).toEqual([]);
  const top = plan.find((p: any) => p.component === 'stepped_top_cap');
  expect(top.closed).toBeTrue();
  for (const [x, y] of top.path_mm) expect(Math.abs(Math.hypot(x, y) - 280)).toBeLessThan(0.01);
  const handle = plan.find((p: any) => p.component === 'vertical_handle');
  for (const [x, y] of handle.path_mm) expect(Math.hypot(x, y)).toBeGreaterThan(279.9);
  for (const v of ['front', 'side']) {
    const seam = audit.views[v].find((p: any) => p.kind === 'visible_depth_discontinuity');
    expect(seam.path_mm).toEqual([[-266, 230], [266, 230]]);
  }
  expect(audit.meaningless_coplanar_internal_lines).toBe(0);
});

test('candidate artifacts and context hashes are current, and IFC gates remain closed', () => {
  const manifest = read(join(dir, 'manifest.json'));
  const candidate = read(join(dir, 'candidate-representations.json'));
  const audit = read(join(dir, 'bst02-semantic-segmentation.json'));
  expect(hash(join(dir, 'candidate-representations.json'))).toBe(manifest.candidate_representations_sha256);
  expect(hash(join(dir, 'bst02-semantic-segmentation.json'))).toBe(manifest.semantic_segmentation.sha256);
  for (const entry of manifest.views) {
    expect(hash(join(root, entry.svg))).toBe(entry.svg_sha256);
    expect(candidate.views[entry.view].proxy_paths_mm).toEqual(audit.views[entry.view].map((p: any) => p.path_mm));
    expect(candidate.views[entry.view].official_cad_paths_mm).toEqual([]);
    const svg = readFileSync(join(root, entry.svg), 'utf8');
    expect(svg).toContain('基于原始高模几何生成的简化图纸表达');
    expect(svg).toContain('stroke="#000"');
    expect(svg).not.toContain('original-highpoly');
  }
  const context = read(join(dir, 'project-context-manifest.json'));
  expect(hash(join(dir, 'project-context-manifest.json'))).toBe(manifest.project_context.manifest_sha256);
  for (const v of context.views) {
    expect(hash(join(root, v.review_crop))).toBe(v.review_crop_sha256);
    expect(hash(join(root, v.review_preview))).toBe(v.review_preview_sha256);
  }
  const approval = read(join(root, 'pipeline/decisions/bst02-drawing-approval.json'));
  expect(approval.status).toBe('pending');
  expect(approval.candidate_manifest_sha256).toBe(hash(join(dir, 'manifest.json')));
  expect(approval.derived_ifc_write_allowed).toBeFalse();
  expect(approval.formal_authoritative_ifc_write_allowed).toBeFalse();
  expect(existsSync(join(dir, 'Baxter-Beside-BST02-derived-drawing.ifc'))).toBeFalse();
});
