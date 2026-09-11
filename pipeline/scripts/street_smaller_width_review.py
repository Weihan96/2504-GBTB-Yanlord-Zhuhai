"""Review smaller Street module evidence without inventing native CAD."""
import json
import xml.etree.ElementTree as ET
from pathlib import Path
from street_linework import normalised_insert
from street_h_linework import records
from falper_sorgente_linework import sha256

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'output/review/highpoly-types/street'
SRC = OUT / 'official-source/AL_Street.dxf'

def main():
    entities, blocks = records(SRC)
    plan = normalised_insert(entities, blocks, '1425', '*U37')
    front = normalised_insert(entities, blocks, '1430', '*U25')
    assert plan['bounds_mm']['size'] == [1080.0, 470.0]
    assert front['bounds_mm']['size'] == [1080.0, 250.0]
    # Preserve the original grey Body rendering and black review silhouette.
    old = ET.parse(OUT / 'plan.svg').getroot()
    layers = ''.join(ET.tostring(x, encoding='unicode') for x in old if x.tag.endswith('path'))
    paths = []
    for path in plan['paths_mm']:
        points = [(120.175 + (x-40)*.85965, 455.789+(470-y)*.85965) for x,y in path]
        paths.append('M '+' L '.join(f'{x:.3f} {y:.3f}' for x,y in points))
    svg = '<svg xmlns="http://www.w3.org/2000/svg" width="1200" height="1080" viewBox="0 0 1200 1080"><rect width="1200" height="1080" fill="white"/>'
    for y,t in [(45,'Street centered sink / closed sides: size evidence'),(80,'Blue solid: native DXF 1080 x 470 / Grey Body + black: project 1000 x 470'),(115,'Common anchor: width center + rear outer edge. No scaling; official width +80 mm.'),(155,'Smaller inferred option: 270 + 450 + 270 = 990 mm (technical PDF p14).'),(188,'990 is a dimensional inference, NOT an existing native DXF view or confirmed order size.'),(221,'PDF p6 instead states minimum 1080: manufacturer confirmation remains necessary.')]:
        svg += f'<text x="40" y="{y}" font-family="Arial" font-size="20" fill="#233446">{t}</text>'
    svg += layers + '<path d="'+' '.join(paths)+'" fill="none" stroke="#087fd6" stroke-width="2.4"/>'
    svg += '<text x="80" y="965" font-family="Arial" font-size="18">No native 990 mm view found. The 450 mm basin is smaller than the 540 mm basin.</text>'
    svg += '<text x="80" y="1000" font-family="Arial" font-size="18">Source label: STREET147 prof.47 + STREET4754 prof.47 / handles 1425, 1430, 1439</text></svg>'
    (OUT/'street-smaller-width-evidence.svg').write_text(svg)
    audit={'schema_version':1,'review_status':'pending_user_review','project_global_id':'1FgLPMw$5B4wBH2ySMkXE1','project_bounds_mm':[1000,470,250], 'native_smaller_complete_view_found':False,'smallest_stored_native_plan_width_mm':1080,'inferred_smaller_centered_configuration':{'top':'STREET147 prof.47','basin':'STREET4745 prof.47','dimension_chain_mm':[270,450,270],'inferred_width_mm':990,'depth_mm':470,'height_mm':250,'difference_from_project_mm':[-10,0,0],'source_pdf_page':14,'status':'dimensional_inference_requires_manufacturer_confirmation','native_dxf_view_available':False,'caveat':'PDF page 6 off-center variant explicitly gives min 1080 mm; do not generalize 990 to every variant.'},'native_reference':{'label':'STREET147 prof.47 + STREET4754 prof.47','plan':plan,'front':front,'side':None,'anchor':'width center + rear outer edge','width_delta_mm':80,'scale':1},'dxf_sha256':sha256(SRC),'dxf_path':str(SRC.relative_to(ROOT)),'technical_pdf_sha256':sha256(OUT/'official-source/ANTONIOLUPI-official-Street-technical.pdf'),'technical_pdf_url':'https://www.antoniolupi.it/uploads/2025/11/27/street_st_web_1764249905846.pdf','official_product_page':'https://www.antoniolupi.it/en/products/sinks/street','derived_ifc_write_allowed':False,'formal_ifc_write_allowed':False}
    (OUT/'street-smaller-width-audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2)+'\n')

if __name__ == '__main__': main()
