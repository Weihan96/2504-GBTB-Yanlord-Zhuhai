#!/usr/bin/env python3
"""Expose restored native curves and unresolved proof requirements without approval."""
import json
import math
from pathlib import Path

from falper_sorgente_linework import ROOT, sha256, write_json

PRODUCT = ROOT / "output/review/highpoly-types/marilyn-01"


def main():
    linework = json.loads((PRODUCT / "official-native-dwg-linework.json").read_text())
    views = {}
    elements = ['<svg xmlns="http://www.w3.org/2000/svg" width="1500" height="850" viewBox="0 0 1500 850">', '<rect width="1500" height="850" fill="white"/>', '<g font-family="Arial" fill="#20252b"><text x="30" y="40" font-size="25">Marilyn 01 — restored native curve audit, pending AutoCAD close-up verification</text><text x="30" y="73" font-size="16">Black: previously delivered native lines · Red: restored rational SPLINE · Orange: fit-point curve approximation requiring verification</text></g>']
    for vi, (name, view) in enumerate(linework["views"].items()):
        paths = view["paths_mm"]
        ids = view["native_entity_indices"]
        rational = set(view["rational_spline_entity_indices"])
        provisional = set(view["fit_point_spline_entity_indices"])
        b = view["bounds_mm"]
        lo, hi = b["minimum"], b["maximum"]
        scale = min(425 / (hi[0]-lo[0]), 585 / (hi[1]-lo[1]))
        x0, y0 = 30 + 500 * vi, 135
        elements.append(f'<text x="{x0}" y="115" font-family="Arial" font-size="21">{name.upper()}: {len(paths)} native entities</text>')
        endpoint_inventory = []
        for idx, path in zip(ids, paths):
            colour = '#d17b00' if idx in provisional else '#cc263d' if idx in rational else '#24272b'
            coords = [f'{x0+(p[0]-lo[0])*scale:.3f},{y0+(hi[1]-p[1])*scale:.3f}' for p in path]
            elements.append(f'<path d="M {" L ".join(coords)}" stroke="{colour}" stroke-width="1.8" fill="none"><title>Native DWG entity {idx}</title></path>')
            for end_name, point in [('start',path[0]),('end',path[-1])]:
                nearest = min(math.dist(point, q) for other, op in zip(ids, paths) if other != idx for q in (op[0],op[-1]))
                endpoint_inventory.append({'native_entity_index':idx,'endpoint':end_name,'nearest_other_endpoint_mm':round(nearest,6)})
        views[name] = {'path_count':len(paths),'restored_rational_indices':sorted(rational),'fit_point_approximation_indices':sorted(provisional),'endpoint_inventory':endpoint_inventory,'closure_proven':False,'closure_note':'Individual manufacturer seam/fold paths can legitimately be open; endpoint distance alone does not prove a missing contour. AutoCAD close-ups are still required.'}
    elements += ['<text x="30" y="780" font-family="Arial" font-size="17">Counts now 44/96/78 after restoring negative-Z OCS arcs. Generated diagnostic, not an AutoCAD screenshot.</text>', '<text x="30" y="812" font-family="Arial" font-size="17">Support outline cycle verified separately. No synthetic joins. Two Front fit curves still await verification.</text>', '</svg>']
    target = PRODUCT / "marilyn-restored-curves-diagnostic.svg"
    target.write_text(''.join(elements))
    write_json(PRODUCT / "marilyn-unfinished-review-audit.json", {
        'schema_version':1,'status':'pending_autocad_closeup_and_visual_approval',
        'source_dwg_sha256':linework['source_dwg_sha256'],
        'source_dwg_hash_verified':sha256(PRODUCT/'official-source/Marilyn_Abaco.dwg')==linework['source_dwg_sha256'],
        'source_linework_sha256':sha256(PRODUCT/'official-native-dwg-linework.json'),
        'diagnostic_svg':str(target.relative_to(ROOT)), 'diagnostic_svg_sha256':sha256(target),
        'views':views,
        'remaining_work':[
            'Capture readable exact 860 x 1000 x 940 variant Plan/Front/Side in AutoCAD; existing full-modelspace image is too small for closure proof.',
            'Resolve two Front fit-only native SPLINE curves against AutoCAD evaluation. Existing generator uses quadratic interpolation, not a verified native cubic curve.',
            'Verify manufacturer outline joins and central component semantics against the readable screenshots before claiming full visual completeness.'
        ],
        'autocad_session_observation':'2026-09-06: File > Open Recent opens exact archived Marilyn_Abaco.dwg. Modelspace viewport is constrained to an approximately 95-pixel strip; Zoom Extents works within that strip, so readable close-ups are still unavailable. Source DWG was not saved.',
        'side_support_closure_audit':'output/review/highpoly-types/marilyn-01/marilyn-side-support-closure-audit.json',
        'derived_ifc_write_allowed':False,'formal_ifc_write_allowed':False,
    })
    print(json.dumps({'diagnostic':str(target),'status':'pending_autocad_closeup_and_visual_approval','counts':{k:v['path_count'] for k,v in views.items()}}))


if __name__ == '__main__':
    main()
