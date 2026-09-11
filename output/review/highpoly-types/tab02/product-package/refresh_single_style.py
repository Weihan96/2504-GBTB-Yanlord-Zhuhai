"""Display-only thumbnail line weight; scene SVG and all coordinates unchanged."""
import json,xml.etree.ElementTree as ET
import migrate_tab02 as task
from verify_tab02 import render
r=json.loads(task.REPORT.read_text());m=json.loads((task.OUT/'manifest.json').read_text())
for row in r['independent_validation']['view_checks']:
    v=row['view'];p=task.OUT/f'TAB02-SINGLE-{v.upper()}.svg';tree=ET.parse(p)
    for e in tree.getroot().iter():
        if e.tag.rsplit('}',1)[-1]=='line':e.set('style','stroke:#111820;stroke-width:0.045;fill:none;stroke-linecap:round;stroke-linejoin:round')
    tree.write(p,encoding='utf-8',xml_declaration=True);render(p,p.with_suffix('.png'),1000)
    row['single_product_svg']=task.pkg.record(p);row['single_product_thumbnail']=task.pkg.record(p.with_suffix('.png'))
    m['single_product_views'][v]=row['single_product_svg'];m['library_previews'][v]=row['single_product_thumbnail']
m['scene_outputs']=r['independent_validation']['view_checks']
m['preview_provenance']['2d']='Exact target LINEWORK cropped from verified real Bonsai scene SVG; display-only thinner rounded stroke, no changed coordinates.'
task.write(task.REPORT,r);task.write(task.OUT/'manifest.json',m)
