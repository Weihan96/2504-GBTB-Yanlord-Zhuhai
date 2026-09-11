// Rasterize actual Bonsai SVGs for the final visual gate; no source IFC edits.
import fs from 'node:fs/promises';
import path from 'node:path';
import sharp from '/Users/jiaxinchen/.cache/codex-runtimes/codex-primary-runtime/dependencies/node/node_modules/sharp';
const out=path.resolve('output/review/approved-product-library');
const validation=JSON.parse(await fs.readFile(path.join(out,'portable-isolation-validation.json'),'utf8'));
if(validation.status!=='structured_pass'||validation.svg_count!==117)throw new Error('All 117 native drawings must pass first');
const catalog=JSON.parse(await fs.readFile(path.join(out,'runtime/catalog.json'),'utf8'));
const sheets=[];
let parts=[];
for(const [i,p] of catalog.products.entries()){
  const directory=path.join(out,'portable-validation',p.id);
  const row=i%6;
  const label=Buffer.from(`<svg width="1500" height="30"><rect width="1500" height="30" fill="#eee"/><text x="12" y="22" font-size="18" font-family="sans-serif">${p.id} — Plan / Front / Side</text></svg>`);
  parts.push({input:label,left:0,top:row*310});
  for(const [j,view] of ['plan','front','side'].entries()){
    const svg=path.join(directory,view+'.svg');
    const png=await sharp(svg,{density:140}).resize({width:500,height:280,fit:'contain',background:'white'}).flatten({background:'white'}).png().toBuffer();
    await fs.writeFile(path.join(directory,view+'.png'),png);
    parts.push({input:png,left:j*500,top:row*310+30});
  }
  if(row===5||i===catalog.products.length-1){
    const file=path.join(out,'portable-validation',`contact-${sheets.length+1}.png`);
    await sharp({create:{width:1500,height:(row+1)*310,channels:3,background:'white'}}).composite(parts).png().toFile(file);
    sheets.push(file);parts=[];
  }
}
await fs.writeFile(path.join(out,'portable-raster-validation.json'),JSON.stringify({status:'rendered_pending_human_inspection',images:117,sheets},null,2)+'\n');
console.log(JSON.stringify({images:117,sheets}));
