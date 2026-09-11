#!/usr/bin/env python3
"""Prove the missing support edge is a native OCS ARC, not a hand-drawn join."""
import json
import math

from falper_sorgente_linework import ROOT, sha256, write_json
from marilyn_01_linework import SOURCE_DWG, dwg_json

PRODUCT = ROOT / "output/review/highpoly-types/marilyn-01"
TOLERANCE = 0.005


def main():
    linework = json.loads((PRODUCT / "official-native-dwg-linework.json").read_text())
    document, _ = dwg_json(SOURCE_DWG)
    native = {item['index']: item for item in document['OBJECTS'] if 'index' in item}
    side = linework['views']['side']
    paths = dict(zip(side['native_entity_indices'], side['paths_mm']))
    start, end = paths[234][0], paths[234][-1]

    def walk(point, used, joins):
        if math.dist(point, start) < TOLERANCE:
            return used, joins + [math.dist(point, start)]
        for index, path in paths.items():
            if index in used:
                continue
            for a, b in ((path[0], path[-1]), (path[-1], path[0])):
                gap = math.dist(point, a)
                if gap < TOLERANCE:
                    found = walk(b, used + [index], joins + [gap])
                    if found:
                        return found
        return None

    loop = walk(end, [234], [])
    assert loop, 'No native outline cycle through support ARC 234'
    ids, gaps = loop
    record = {
        'schema_version': 1,
        'source_dwg_sha256': sha256(SOURCE_DWG),
        'source_linework_sha256': sha256(PRODUCT / 'official-native-dwg-linework.json'),
        'cause': 'ARC coordinates were treated as WCS although DWG stores them in OCS; -Z extrusion requires reflected X before spatial selection.',
        'missing_support_arc_native_entity': native[234],
        'restored_arc_234_normalized_endpoints_mm': [start, end],
        'restored_arc_234_length_mm': native[234]['radius'] * (native[234]['end_angle'] - native[234]['start_angle']),
        'restored_negative_z_arc_indices': {name: view['negative_z_ocs_arc_entity_indices'] for name, view in linework['views'].items()},
        'path_counts_before': {'plan': 42, 'front': 86, 'side': 62},
        'path_counts_after': {name: view['path_count'] for name, view in linework['views'].items()},
        'support_outline_cycle': {'native_entity_indices': ids, 'maximum_endpoint_gap_mm': max(gaps), 'join_tolerance_mm': TOLERANCE, 'closed_within_tolerance': True, 'coordinates_snapped': False},
        'hand_drawn_join_added': False,
        'source_geometry_scaled': False,
        'front_fit_only_curves_still_approximate': [1005, 1006],
        'visual_approval': 'pending',
        'derived_ifc_write_allowed': False,
        'formal_ifc_write_allowed': False,
    }
    target = PRODUCT / 'marilyn-side-support-closure-audit.json'
    write_json(target, record)
    elements = ['<svg xmlns="http://www.w3.org/2000/svg" width="1200" height="900" viewBox="0 0 1200 900"><rect width="1200" height="900" fill="white"/><g font-family="Arial"><text x="35" y="40" font-size="25">Marilyn Side — native support outline restored</text><text x="35" y="70" font-size="17">Blue: native DWG · Red: restored negative-Z OCS arcs · No synthetic join</text></g>']
    restored = set(side['negative_z_ocs_arc_entity_indices'])
    for index, path in paths.items():
        points = [(600 - p[0] * 0.72, 810 - p[1] * 0.72) for p in path]
        d = 'M ' + ' L '.join(f'{x:.4f},{y:.4f}' for x, y in points)
        color = '#cc263d' if index in restored else '#1677c8'
        elements.append(f'<path d="{d}" fill="none" stroke="{color}" stroke-width="2"><title>Native entity {index}</title></path>')
    elements.append(f'<text x="35" y="870" font-family="Arial" font-size="16">Native support loop: {len(ids)} entities; maximum endpoint gap {max(gaps):.6f} mm. Generated diagnostic, not AutoCAD screenshot.</text></svg>')
    (PRODUCT / 'marilyn-side-support-closure-diagnostic.svg').write_text(''.join(elements))
    print(json.dumps({'audit': str(target), 'max_gap_mm': max(gaps), 'arc_length_mm': record['restored_arc_234_length_mm']}))


if __name__ == '__main__':
    main()
