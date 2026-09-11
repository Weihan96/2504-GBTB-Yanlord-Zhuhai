import { expect, test } from 'bun:test';
import { readFileSync } from 'node:fs';
import { resolve } from 'node:path';
const root=resolve(import.meta.dir,'../..');
const audit=JSON.parse(readFileSync(resolve(root,'output/review/highpoly-types/wd02/wd02-semantic-segmentation.json'),'utf8'));
test('WD02 classifies every actual Body item and separates closed glass from opaque panels',()=>{
  expect(audit.representation_item_count).toBe(35);
  expect(new Set(audit.components.map((c:any)=>c.item_id)).size).toBe(35);
  expect(audit.components.find((c:any)=>c.semantic==='door_glass').size_mm[1]).toBeCloseTo(7.001,2);
  const frontGlass=audit.visibility.front.find((c:any)=>c.semantic==='door_glass');
  expect(frontGlass.transparent).toBeTrue();
  for(const semantic of ['upper_shelf','clothes_rail','middle_drawer_front_frame','lower_drawer_front_frame']) {
    expect(audit.visibility.front.find((c:any)=>c.semantic===semantic).visible_area_mm2).toBeGreaterThan(0);
  }
});
test('WD02 top and side occlusion remove covered internal furniture',()=>{
  for(const view of ['plan','side']) {
    expect(audit.visibility[view].find((c:any)=>c.semantic==='clothes_rail').visible_area_mm2).toBeLessThan(0.01);
    expect(audit.visibility[view].find((c:any)=>c.semantic==='middle_drawer_front_frame').visible_area_mm2).toBeLessThan(0.01);
  }
  expect(audit.visibility.plan.find((c:any)=>c.semantic==='top_deck').visible_area_mm2).toBeGreaterThan(290000);
  expect(audit.visibility.side.find((c:any)=>c.semantic==='left_side_panel').visible_area_mm2).toBeGreaterThan(1300000);
});
test('WD02 candidate component rings are closed, bounded and finite',()=>{
  for(const view of Object.values(audit.views) as any[]) for(const item of view.paths) {
    const path=item.path;
    expect(path[0]).toEqual(path.at(-1));
    expect(path.length).toBeGreaterThanOrEqual(4);
    expect(path.flat().every(Number.isFinite)).toBeTrue();
    expect(audit.components.some((c:any)=>c.item_id===item.item_id)).toBeTrue();
  }
});
