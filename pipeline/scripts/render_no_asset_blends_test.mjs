import fs from 'node:fs/promises';
import path from 'node:path';
import sharp from '/Users/jiaxinchen/.cache/codex-runtimes/codex-primary-runtime/dependencies/node/node_modules/sharp';
const out = path.resolve('output/review/approved-product-library');
const reportPath = path.join(out, 'no-asset-blends-test.json');
const report = JSON.parse(await fs.readFile(reportPath, 'utf8'));
const destination = path.join(out, 'no-asset-blends-test');
await fs.mkdir(destination, {recursive: true});
const parts = [];
for (const [row, result] of report.results.entries()) {
  const label = Buffer.from(`<svg width="1500" height="35"><rect width="1500" height="35" fill="#eee"/><text x="10" y="24" font-family="sans-serif" font-size="20">${result.product} | NO ASSET BLENDS | Plan / Front / Side</text></svg>`);
  parts.push({input: label, left: 0, top: row * 435});
  for (const [col, drawing] of result.drawings.entries()) {
    const file = path.join(destination, result.product + '-' + drawing.role.toLowerCase() + '.svg');
    await fs.copyFile(drawing.svg, file);
    drawing.persisted_svg = file;
    const png = await sharp(file, {density: 140}).resize({width: 500, height: 400, fit: 'contain', background: 'white'}).flatten({background: 'white'}).png().toBuffer();
    await fs.writeFile(file.replace(/\.svg$/, '.png'), png);
    parts.push({input: png, left: col * 500, top: row * 435 + 35});
  }
}
const contact = path.join(destination, 'contact.png');
await sharp({create: {width: 1500, height: report.results.length * 435, channels: 3, background: 'white'}}).composite(parts).png().toFile(contact);
report.visual = {status: 'rendered_pending_inspection', images: 6, contact};
await fs.writeFile(reportPath, JSON.stringify(report, null, 2) + '\n');
console.log(contact);
