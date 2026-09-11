// Reduce only asset-card rasters. Preserve source renders and all drawing files.
import fs from 'node:fs/promises';
import path from 'node:path';
import os from 'node:os';
import {createHash} from 'node:crypto';
import {execFileSync} from 'node:child_process';
import sharp from '/Users/jiaxinchen/.cache/codex-runtimes/codex-primary-runtime/dependencies/node/node_modules/sharp';

const root=path.resolve(import.meta.dir,'../..');
const out=path.join(root,'output/review/approved-product-library');
const runtime=path.join(out,'runtime');
const read=async p=>JSON.parse(await fs.readFile(p,'utf8'));
const hash=b=>createHash('sha256').update(b).digest('hex');
const digest=async p=>hash(await fs.readFile(p));
const index=()=>hash(execFileSync('git',['ls-files','--stage','-z'],{cwd:root}));
const plan=await read(path.join(out,'portable-cleanup-plan.json'));
if(index()!==plan.protected_index_sha256 || await digest(path.join(root,'2504 GBTB Yanlord Zhuhai.ifc'))!==plan.formal_ifc_sha256)
  throw Error('Protected baseline changed');
const catalog=await read(path.join(runtime,'catalog.json'));
const original=await read(path.join(out,'all-review-catalog.json'));
const staging=await fs.mkdtemp(path.join(os.tmpdir(),'ifc-asset-thumbs-'));
const files=[];
for(const p of catalog.products){
  const target=path.join(runtime,p.previews.iso);
  if(await digest(target)!==p.preview_evidence.iso.png_sha256)throw Error(`${p.id}: unknown thumbnail change`);
  const source=path.resolve(out,original.products.find(q=>q.id===p.id).previews.iso);
  const sourceHash=await digest(source);
  const temp=path.join(staging,p.id+'.png');
  const old=await sharp(target).metadata();
  await sharp(source).resize({width:384,height:384,fit:'inside',withoutEnlargement:true})
    .png({palette:true,colours:128,dither:.8,compressionLevel:9,effort:10}).toFile(temp);
  const meta=await sharp(temp).metadata();
  if(Math.max(meta.width,meta.height)>384)throw Error('Oversized thumbnail');
  const record={id:p.id,path:p.previews.iso,source:path.relative(root,source),source_sha256:sourceHash,
    before_bytes:(await fs.stat(target)).size,after_bytes:(await fs.stat(temp)).size,
    before_dimensions:[old.width,old.height],dimensions:[meta.width,meta.height],sha256:await digest(temp)};
  files.push(record);
  p.preview_evidence.iso={kind:'camera_render_3d',purpose:'asset_thumbnail',png_sha256:record.sha256,
    width:meta.width,height:meta.height,max_dimension:384,source_sha256:sourceHash};
}
// Preflight all images before replacing any generated runtime raster.
for(const f of files)await fs.copyFile(path.join(staging,f.id+'.png'),path.join(runtime,f.path));
await fs.writeFile(path.join(runtime,'catalog.json'),JSON.stringify(catalog,null,2)+'\n');
const report={status:'compressed_pending_live_verification',version:'0.8.1',files,
  before_bytes:files.reduce((s,f)=>s+f.before_bytes,0),after_bytes:files.reduce((s,f)=>s+f.after_bytes,0),
  original_renders_preserved:true,protected_index_sha256:plan.protected_index_sha256,staging};
await fs.writeFile(path.join(out,'thumbnail-compression-validation.json'),JSON.stringify(report,null,2)+'\n');
if(index()!==plan.protected_index_sha256)throw Error('Index changed');
console.log(JSON.stringify({images:files.length,before_bytes:report.before_bytes,after_bytes:report.after_bytes}));
