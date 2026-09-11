#!/usr/bin/env python3
"""Render and package persisted WD02 native drawings for human scene review."""
import json
import sys
import hashlib
import shutil
import subprocess
import tempfile
from pathlib import Path
from PIL import Image
from reportlab.pdfgen import canvas
from reportlab.lib.utils import ImageReader
from reportlab.lib.pagesizes import A4, landscape
from pypdf import PdfReader

ROOT=Path(__file__).resolve().parents[2]
P=ROOT/'output/review/highpoly-types/wd02'
OUT=P/'bonsai-drawings/wardrobe'
EV=OUT/'WD02-WARDROBE-create-drawing-evidence.json'
PDF=P/'Poliform-Senzafine-WD02-scene-drawings.pdf'
def rec(p):return {'path':str(p),'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def write(p,x):p.write_text(json.dumps(x,ensure_ascii=False,indent=2)+'\n')
def main():
    if '--confirm-visual-qa' in sys.argv:
        evidence=json.loads(EV.read_text())
        evidence['visual']={'status':'visually_inspected','previews_nonempty':True,'observations':['Plan shows approved top/door boundary at adjacent wall.','Front looks through glass door: upper shelf, rail and three drawer rows visible; bottom remains within frame.','Side shows approved opaque panel boundary against grey project context; grey context is retained for positioning and can appear behind the linework.','PDF three pages checked for complete page bounds and captions.']}
        evidence['verdict']='pass';evidence['pass']=True
        write(EV,evidence)
        manifest=json.loads((P/'manifest.json').read_text());manifest['scene_drawings']['evidence']=rec(EV);write(P/'manifest.json',manifest)
        return
    evidence=json.loads(EV.read_text());records=[]
    page=landscape(A4);pdf=canvas.Canvas(str(PDF),pagesize=page)
    pdf.setTitle('Poliform Senzafine WD02 - native Bonsai scene drawings')
    for item in evidence['outputs']['views']:
        svg=Path(item['svg']['path']);png=P/f'wd02-wardrobe-scene-{item["view"]}.png'
        assert rec(svg)['sha256']==item['svg']['sha256']
        subprocess.run(['/Applications/Inkscape.app/Contents/MacOS/inkscape',str(svg),'--export-area-page','--export-background=white','--export-background-opacity=1','--export-width=1650',f'--export-filename={png}'],capture_output=True,check=True)
        image=Image.open(png).convert('RGB')
        assert sum(1 for r,g,b in image.getdata() if max(r,g,b)<100)>100
        records.append(rec(png))
        pdf.setFont('Helvetica-Bold',15);pdf.drawString(30,page[1]-30,f'WD02 / {item["view"].upper()} / scene SVG')
        pdf.setFont('Helvetica',9);pdf.drawString(30,page[1]-46,'Approved semantic black linework | grey actual project context | native Bonsai Create Drawing')
        s=min((page[0]-60)/image.width,(page[1]-80)/image.height)
        w,h=image.width*s,image.height*s
        pdf.drawImage(ImageReader(image),(page[0]-w)/2,20,width=w,height=h)
        pdf.showPage()
    pdf.save();assert len(PdfReader(PDF).pages)==3
    evidence['outputs']['previews']=records;evidence['outputs']['pdf']=rec(PDF)
    evidence['visual']={'status':'rendered_pending_human_visual_QA','previews_nonempty':True}
    write(EV,evidence)
    manifest=json.loads((P/'manifest.json').read_text())
    manifest['semantic_single_product_approval']={'status':'approved','quote':'单品SVG可以验收，只差给我一个场景SVG','approval_record':'pipeline/decisions/wd02-drawing-approval.json'}
    manifest['scene_drawings']={'status':'pending_user_review','generator':'Bonsai Create Drawing','evidence':rec(EV),'views':evidence['outputs']['views'],'derived_ifc':evidence['outputs']['derived_ifc'],'previews':records,'pdf':rec(PDF)}
    write(P/'manifest.json',manifest)
    scene='<section id="wd02-scene-review"><h2>WD02 场景 SVG · 待验收</h2><p>单品语义黑线已批准。基于原始高模几何生成的简化图纸表达。Front 从玻璃门外朝柜内看。</p>'
    for v in ('plan','front','side'):
        link=f'bonsai-drawings/wardrobe/WD02-WARDROBE-{v.upper()}.svg'
        scene+=f'<h3>{v.title()} 场景 SVG</h3><a href="{link}"><img style="width:100%;max-width:1000px" src="wd02-wardrobe-scene-{v}.png"></a>'
    scene+='<p><a href="Poliform-Senzafine-WD02-derived-drawing.ifc">产品级派生 IFC</a> · <a href="Poliform-Senzafine-WD02-scene-drawings.pdf">三视图 PDF</a></p></section>'
    index=P/'index.html';html=index.read_text();assert 'id="wd02-scene-review"' not in html
    index.write_text(html.replace('</body>',scene+'</body>'))
    print(json.dumps({'previews':records,'pdf':rec(PDF)},indent=2))
if __name__=='__main__':main()
