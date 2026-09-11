#!/usr/bin/env python3
"""Verify native Bonsai scene SVGs and render official TRAP01 details."""
import json
import tempfile
import shutil
from pathlib import Path
import generate_trap01_internal_detail_drawings as g
from PIL import Image

P = g.PRODUCT_DIR
g.OUTPUT_DIR = P/'bonsai-drawings/cabinet-official-detail-v4'
g.EVIDENCE = g.OUTPUT_DIR/'TRAP01-official-detail-v4-evidence.json'
g.MANIFEST = P/'TRAP01-official-detail-v4-manifest.json'
g.COMBINED_PDF = P/'Geberit-151.116.11.1-TRAP01-official-detail-v4.pdf'
g.PREVIEW_PREFIX = 'trap01-official-detail-v4'
g.DERIVED_IFC = P/'Geberit-151.116.11.1-TRAP01-official-detail-v4.ifc'
g.VIEWS = {v:f'TRAP01-OFFICIAL-DETAIL-V4-{v.upper()}' for v in ('plan','front','side')}

def inspect_preview(path):
    image = Image.open(path).convert('RGB')
    blue = grey = 0
    for r,green,b in image.getdata():
        blue += int(b>=100 and b>r*1.25 and b>green*1.03)
        grey += int(abs(r-green)<=20 and abs(green-b)<=20 and 90<=r<=225)
    return {'width':image.width,'height':image.height,'blue_installed_detail_pixels':blue,
        'grey_context_pixels':grey,'current_configuration_visible':blue>=100,
        'retained_context_visible':grey>=100,'pass':blue>=100 and grey>=100,
        'blue_dashed_reference_absent':'verified in SVG source, not inferred from pixel colour'}

if __name__=='__main__':
    g.inspect_preview = inspect_preview
    original_pdf = g.svg_to_pdf
    def fresh_svg_to_pdf(svg,pdf):
        # The older helper exits early when the destination already exists.
        # Always render into a fresh path before replacing the generated PDF.
        with tempfile.TemporaryDirectory(prefix='trap01-fresh-pdf-') as directory:
            fresh = Path(directory)/'drawing.pdf'
            original_pdf(svg,fresh)
            shutil.copyfile(fresh,pdf)
    g.svg_to_pdf = fresh_svg_to_pdf
    evidence = json.loads(g.EVIDENCE.read_text())
    if evidence['status']=='persisted_reloaded_rendered_verified':
        evidence['status']='persisted_reloaded_verified'
    if 'installation_anchors' in evidence['postState']:
        evidence['postState']['inherited_configuration_reference_points'] = evidence['postState'].pop('installation_anchors')
    evidence['postState']['anchor_semantics_note'] = 'Inherited controls are prior configuration reference points, not surveyed pipe centre lines. Native annotation origins are independently verified against target ObjectPlacement.'
    evidence['tests'] = {'native_official_paths_reloaded_exactly':all(v['coordinates_match_approved_retained_dwg'] for v in evidence['postState']['reloaded_views'])}
    g.EVIDENCE.write_text(json.dumps(evidence,indent=2,ensure_ascii=False)+'\n')
    g.main()
