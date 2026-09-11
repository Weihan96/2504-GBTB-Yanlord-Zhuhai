// Publish linework previews from verified native Bonsai output, never from camera PNGs.
import fs from 'node:fs/promises';
import path from 'node:path';
import os from 'node:os';
import { createHash } from 'node:crypto';
import { execFileSync } from 'node:child_process';
import sharp from '/Users/jiaxinchen/.cache/codex-runtimes/codex-primary-runtime/dependencies/node/node_modules/sharp';

const root = path.resolve(import.meta.dir, '../..');
const out = path.join(root, 'output/review/approved-product-library');
const runtime = path.join(out, 'runtime');
const read = async p => JSON.parse(await fs.readFile(p, 'utf8'));
const hash = bytes => createHash('sha256').update(bytes).digest('hex');
const digest = async p => hash(await fs.readFile(p));
const write = (p, value) => fs.writeFile(p, JSON.stringify(value, null, 2) + '\n');
const index = () => execFileSync('git', ['ls-files', '--stage', '-z'], { cwd: root });
const beforeIndex = hash(index());
const formal = path.join(root, '2504 GBTB Yanlord Zhuhai.ifc');
const beforeFormal = await digest(formal);
const cleanup = await read(path.join(out, 'portable-cleanup-plan.json'));
if (beforeIndex !== cleanup.protected_index_sha256 || beforeFormal !== cleanup.formal_ifc_sha256)
  throw Error('Protected baseline changed; inspect before publishing');
const catalog = await read(path.join(runtime, 'catalog.json'));
const original = await read(path.join(out, 'all-review-catalog.json'));
const staging = await fs.mkdtemp(path.join(os.tmpdir(), 'ifc-linework-previews-'));
const copies = [];
const products = [];
const roles = ['plan', 'front', 'side'];
for (const p of catalog.products) {
  const sourceDir = path.join(out, 'portable-validation', p.id);
  const proof = await read(path.join(sourceDir, 'result.json'));
  if (proof.verdict !== 'technical_pass' || proof.source_sha256 !== p.ifc_sha256 ||
      await digest(path.join(runtime, p.ifc_path)) !== p.ifc_sha256)
    throw Error(`${p.id}: native drawing evidence is stale`);
  const old = original.products.find(q => q.id === p.id);
  const evidence = {};
  for (const view of roles) {
    const drawing = proof.drawings.find(d => d.role === view[0].toUpperCase() + view.slice(1));
    const bytes = await fs.readFile(path.join(sourceDir, `${view}.svg`));
    if (hash(bytes) !== drawing.sha256 || drawing.target_geometry_nodes <= 0 ||
        !drawing.selected.length || drawing.selected.some(s => s.role !== drawing.role || s.old_body_or_annotation_emitted))
      throw Error(`${p.id}/${view}: invalid direction-specific linework`);
    const svg = bytes.toString();
    // Keep the runtime portable: linework may reference local SVG defs, not external resources.
    if (/<image\b/i.test(svg) || /(?:href\s*=\s*["'](?!#)|url\(\s*["']?(?!#)[a-z])/i.test(svg))
      throw Error(`${p.id}/${view}: drawing contains an external/image dependency`);
    const svgRelative = `previews/${p.id}/${view}.svg`;
    const pngRelative = `previews/${p.id}/${view}.png`;
    const stageSvg = path.join(staging, `${p.id}-${view}.svg`);
    const stagePng = path.join(staging, `${p.id}-${view}.png`);
    await fs.writeFile(stageSvg, bytes);
    await sharp(bytes, { density: 180 }).resize({ width: 1600, height: 1200, fit: 'contain', background: 'white' })
      .flatten({ background: 'white' }).png().toFile(stagePng);
    const stats = await sharp(stagePng).stats();
    if (stats.channels.every(c => c.min === c.max)) throw Error(`${p.id}/${view}: blank raster`);
    // Retain the pre-existing camera/review image at its historical source, before replacement.
    const previous = path.resolve(out, old.previews[view]);
    const previousHash = await digest(previous);
    if (!p.preview_evidence && await digest(path.join(runtime, p.previews[view])) !== previousHash)
      throw Error(`${p.id}/${view}: unknown existing preview; refusing overwrite`);
    evidence[view] = { kind: 'bonsai_body_linework', role: drawing.role, source_ifc_sha256: p.ifc_sha256,
      svg: svgRelative, svg_sha256: hash(bytes), png_sha256: await digest(stagePng),
      generator: proof.generator, width: 1600, height: 1200 };
    copies.push([stageSvg, path.join(runtime, svgRelative)], [stagePng, path.join(runtime, pngRelative)]);
    p.previews[view] = pngRelative;
    products.push({ id: p.id, view, previous_image_preserved: path.relative(root, previous),
      previous_sha256: previousHash, ...evidence[view] });
  }
  evidence.iso = { ...p.preview_evidence?.iso, kind: 'camera_render_3d', png_sha256: await digest(path.join(runtime, p.previews.iso)) };
  p.preview_evidence = evidence;
}
// Validate every product before publishing any new preview or catalog.
for (const [source, target] of copies) await fs.copyFile(source, target);
catalog.preview_contract = 'plan_front_side_bonsai_linework_iso_camera_render';
await write(path.join(runtime, 'catalog.json'), catalog);
for (const name of await fs.readdir(path.join(root, 'pipeline/addons/highpoly_review_library'))) {
  if (name.endsWith('.py')) await fs.copyFile(path.join(root, 'pipeline/addons/highpoly_review_library', name),
    path.join(runtime, 'addon/highpoly_review_library', name));
}
if (hash(index()) !== beforeIndex || await digest(formal) !== beforeFormal) throw Error('Protected files changed');
await write(path.join(out, 'drawing-preview-validation.json'), {
  status: 'packaged_pending_ui_verification', version: '0.7.1', product_count: catalog.products.length,
  linework_count: products.length, preview_contract: catalog.preview_contract, products,
  formal_sha256: beforeFormal, protected_index_sha256: beforeIndex,
  boundary: 'Only preview rasters, SVG copies, runtime catalog and addon mirror changed. IFC, approvals and 3D previews unchanged.',
  courseEvidence: { mode: 'embedded-course-index', lesson: '085', fact: 'Drawing generates SVG; viewing output is separate from camera rendering.' },
  staging_directory: staging,
});
console.log(JSON.stringify({ products: catalog.products.length, linework_previews: products.length, staging }));
