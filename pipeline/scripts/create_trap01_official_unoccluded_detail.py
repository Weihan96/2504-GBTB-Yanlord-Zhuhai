#!/usr/bin/env python3
"""Render approved TRAP01 native 2D installation lines in unoccluded Drawings.

Runs in a private Blender process using the installed public bridge handlers;
does not start a server or alter the unrelated interactive Blender session.
"""
from pathlib import Path
import json
import sys
import copy
import hashlib
import traceback
import math

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'pipeline/scripts'))
PRODUCT = ROOT / 'output/review/highpoly-types/trap01'
OUT = PRODUCT / 'bonsai-drawings/cabinet-official-detail-v4'
IFC = PRODUCT / 'Geberit-151.116.11.1-TRAP01-official-detail-v4.ifc'
INPUT_IFC = PRODUCT / 'Geberit-151.116.11.1-TRAP01-derived-drawing.ifc'
EVIDENCE = OUT / 'TRAP01-official-detail-v4-evidence.json'

def configure_components(view, paths, alignment, dx, dz):
    """Move intact official joints; shorten only source straight pipe spans.

    Source path indices are pinned to the hashed 151.116.11.1 family files.
    A/L elbow tangent: Z=336.5. Horizontal pipe: X=136..391.
    A/L upper assembly moves as one part; its far-end inlet additionally moves X.
    """
    result, audit = [], []
    for index, path in enumerate(paths):
        transformed = []
        for x,y in path:
            if view=='plan':
                if index in (8,9,10,11,14,15) or index>=92:
                    x += dx
                elif index in (12,13) and x>300:
                    x += dx
            elif view=='front':
                if index<=23:
                    y += dz
                elif index in (24,25) and y>300:
                    y += dz
                if 8<=index<=19:
                    x += dx
                elif index in (20,21) and x>300:
                    x += dx
            elif view=='side':
                if index<=24 and index not in (20,21):
                    y += dz
                elif index in (20,21) and y>300:
                    y += dz
            transformed.append([x+alignment[0],y+alignment[1]])
        straight_adjustable = (view=='plan' and index in (12,13)) or (view=='front' and index in (20,21,24,25)) or (view=='side' and index in (20,21))
        deltas = [(q[0]-p[0],q[1]-p[1]) for p,q in zip(path,transformed)]
        rigid_error = max(max(abs(a-b) for a,b in zip(d,deltas[0])) for d in deltas)
        assert straight_adjustable or rigid_error<1e-8
        if straight_adjustable:
            assert len(path)==2
            assert abs(path[0][0]-path[1][0])<1e-8 or abs(path[0][1]-path[1][1])<1e-8
        result.append(transformed)
        audit.append({'source_path_index':index,'vertex_count':len(path),'straight_adjustable_span':straight_adjustable,
            'fixed_path_rigid_translation_error_mm':rigid_error if not straight_adjustable else None})
    return result,audit

def main():
    import bpy
    import addon_utils
    import ifcopenshell.util.placement
    import ifcopenshell.util.element
    from bonsai import tool
    import bonsai_bridge as bridge
    import create_trap01_internal_mapped_detail_drawings as mapped
    import create_trap01_internal_detail_drawings as base
    import trap01_review as review
    import trap01_linework as native
    OUT.mkdir(parents=True, exist_ok=True)
    for code, expected_hash in native.EXPECTED.items():
        assert base.sha256(PRODUCT/'official-source'/f'151.116.11.1_{code}.dwg')==expected_hash
    # Each run starts from the approved source. No previous evidence is deleted.
    verify_only = '--verify-only' in sys.argv
    bpy.ops.bim.load_project(filepath=str(IFC if verify_only else INPUT_IFC), should_start_fresh_session=True)
    pre = bridge._h_get_scene_info({})
    save_as = {'resumed_saved_output':str(IFC)} if verify_only else bridge._h_save_ifc_file({'output_path':str(IFC), 'overwrite':True, 'reload':False})
    # Save-as does not update Bonsai's document path until reload.
    if not verify_only:
        bridge._reload_ifc_project(str(IFC))
    candidate = json.loads((PRODUCT/'candidate-representations.json').read_text())
    source = json.loads((PRODUCT/'official-native-dwg-linework.json').read_text())
    # The legacy extractor excluded CIRCLE from its contour entity whitelist.
    # Restore all three exact native G-view inlet circles, without altering the
    # historical extraction JSON or changing existing path indices.
    g_payload = native.shared.dwg_json(PRODUCT/'official-source/151.116.11.1_G.dwg')
    layers = {native.shared.handle_value(e['handle']):e['name'] for e in g_payload['OBJECTS'] if e.get('object')=='LAYER'}
    contour = [e for e in g_payload['OBJECTS'] if e.get('entity') and layers.get(native.shared.handle_value(e.get('layer',[0])))==native.shared.CONTOUR_LAYER]
    old_entities = [e for e in contour if e['entity'] in {'LINE','SPLINE','ARC','ELLIPSE'}]
    oriented = [[[-p[1],p[0]] for p in native.entity_path(e)] for e in old_entities]
    minimum = [min(p[a] for path in oriented for p in path) for a in (0,1)]
    circles = [e for e in contour if e['entity']=='CIRCLE']
    assert len(circles)==3
    for circle in circles:
        cx,cy = circle['center'][:2]
        radius = circle['radius']
        path = []
        for i in range(73):
            angle = math.tau*i/72
            x,y = cx+radius*math.cos(angle),cy+radius*math.sin(angle)
            path.append([-y-minimum[0],x-minimum[1]])
        source['views']['plan']['paths_mm'].append(path)
    approved_audit = json.loads((PRODUCT/'adjustable-overlay-audit.json').read_text())
    actual = copy.deepcopy(candidate)
    official_audit = {}
    dx = -approved_audit['views']['front']['clipped_extension_mm']['horizontal']
    dz = -approved_audit['views']['front']['clipped_extension_mm']['vertical']
    for view in ('plan','front','side'):
        aligned, retained, clipped, audit = review.align_and_segment_official(
            view, candidate['views'][view]['proxy_paths_mm'], source['views'][view]['paths_mm'])
        configured, path_audit = configure_components(view,source['views'][view]['paths_mm'],audit['translation_mm'],dx,dz)
        actual['views'][view]['proxy_paths_mm'] = configured
        official_audit[view] = dict(audit, source_dwg_sha256=source['views'][view]['source_dwg_sha256'],
            configured_paths_mm=configured, configured_path_count=len(configured),
            upper_assembly_translation_z_mm=dz, distal_joint_translation_x_mm=dx,
            path_transform_audit=path_audit, reference_only_excess_paths_mm=clipped,
            note='Old retained crop omitted the whole upper assembly. Current paths preserve every official source path with rigid joint translations and straight pipe shortening only.')
        if view=='plan':
            official_audit[view]['restored_native_circles'] = [{'radius_mm':e['radius'],'source_handle':e.get('handle'),'centre_native_mm':e['center']} for e in circles]
    actual_path = OUT/'native-retained-installation-input.json'
    actual_path.write_text(json.dumps(actual,indent=2)+'\n')
    mapped.CANDIDATE = actual_path
    mapped.DERIVED_IFC = IFC
    mapped.OUTPUT_DIR = OUT
    mapped.PREPERSIST = OUT/'TRAP01-official-detail-v4-prepersist.json'
    base.DERIVED_IFC = IFC
    base.SOURCE_KIND = 'native_dwg_configured_installation'
    base.SOURCE_LABEL_ZH = '基于官方原生DWG固定接头刚性移位及直管长度调整的二维安装细节表达'
    base.BLACK = '#1677c8'
    base.EXPECTED_PATH_COUNTS = {v:len(actual['views'][v]['proxy_paths_mm']) for v in actual['views']}
    for view, definition in mapped.VIEW_DEFINITIONS.items():
        definition['drawing_name'] = f'TRAP01-OFFICIAL-DETAIL-V4-{view.upper()}'
    base.VIEW_DEFINITIONS = mapped.VIEW_DEFINITIONS
    original_curve = base.curve_representation
    def official_curve(model, context, view, paths):
        rep = original_curve(model,context,view,paths)
        for item in rep.Items:
            for styled in item.StyledByItem:
                for style in styled.Styles:
                    style.Name = 'TRAP01 official configured installation blue'
                    style.CurveColour.Red = 22/255
                    style.CurveColour.Green = 119/255
                    style.CurveColour.Blue = 200/255
        return rep
    base.curve_representation = official_curve
    if not verify_only:
        mapped.main()
    data = json.loads(mapped.PREPERSIST.read_text())
    model = tool.Ifc.get()
    target = model.by_guid(base.TARGET_GLOBAL_ID)
    before_placement = ifcopenshell.util.placement.get_local_placement(target.ObjectPlacement).tolist()
    before_body_ids = [r.id() for r in target.Representation.Representations if r.RepresentationIdentifier=='Body']
    # Mark old contour-only Drawing state superseded without deleting its files.
    superseded = []
    for drawing in model.by_type('IfcAnnotation'):
        if drawing.ObjectType=='DRAWING' and 'TRAP01-CABINET-INTERNAL-DETAIL' in (drawing.Name or ''):
            drawing.Description = 'SUPERSEDED by TRAP01-OFFICIAL-DETAIL-V4: old proxy-only outline omitted native DWG component details.'
            superseded.append(drawing.GlobalId)
    for rec in data['outputs']['views']:
        ann = model.by_guid(rec['annotation_global_id'])
        ann.Name = f"TRAP01 configured official 2D details / {rec['view']}"
        ann.Description = base.SOURCE_LABEL_ZH + '; fixed geometry is not scaled; adjustable excess is reference-only and omitted from this scene.'
    persisted = {'saved':True,'reloaded':True,'method':'loaded previously persisted IFC for independent verification'} if verify_only else bridge._h_save_ifc_file({'output_path':str(IFC),'overwrite':True,'reload':True})
    model = tool.Ifc.get()
    target = model.by_guid(base.TARGET_GLOBAL_ID)
    reloaded = []
    for rec in data['outputs']['views']:
        ann = model.by_guid(rec['annotation_global_id'])
        paths = [p for rep in ann.Representation.Representations for item in rep.Items if item.is_a('IfcGeometricCurveSet') for p in item.Elements]
        assert len(paths)==rec['persisted_path_count_expected']
        expected = actual['views'][rec['source_product_view']]['proxy_paths_mm']
        expected3d = [[base.coordinates_mm(rec['source_product_view'],*p) for p in path] for path in expected]
        measured = [[tuple(p.Coordinates) for p in path.Points] for path in paths]
        assert [len(p) for p in measured] == [len(p) for p in expected3d]
        coordinate_error = max(abs(a-b) for p,q in zip(measured,expected3d) for point,expect in zip(p,q) for a,b in zip(point,expect))
        assert coordinate_error < 1e-8, coordinate_error
        placement = ifcopenshell.util.placement.get_local_placement(ann.ObjectPlacement).tolist()
        placement_error = max(abs(a-b) for row,expected_row in zip(placement,before_placement) for a,b in zip(row,expected_row))
        assert placement_error < 0.001, placement_error
        reloaded.append({'view':rec['view'],'source_product_view':rec['source_product_view'],
            'persisted_paths':len(paths),'coordinates_match_approved_retained_dwg':True,'maximum_coordinate_roundtrip_error_mm':coordinate_error,
            'annotation_placement_matches_target_within_0_001_mm':True,'maximum_placement_error_mm':placement_error,'placement_mm':placement})
    assert ifcopenshell.util.placement.get_local_placement(target.ObjectPlacement).tolist()==before_placement
    assert [r.id() for r in target.Representation.Representations if r.RepresentationIdentifier=='Body']==before_body_ids
    original = __import__('ifcopenshell').open(str(INPUT_IFC))
    original_target = original.by_guid(base.TARGET_GLOBAL_ID)
    assert ifcopenshell.util.placement.get_local_placement(original_target.ObjectPlacement).tolist()==before_placement
    def body_digest(file, element):
        body = [r for r in element.Representation.Representations if r.RepresentationIdentifier=='Body']
        entities = {e.id():str(e) for r in body for e in file.traverse(r)}
        return hashlib.sha256('\n'.join(entities[i] for i in sorted(entities)).encode()).hexdigest()
    input_body_hash = body_digest(original,original_target)
    assert input_body_hash == body_digest(model,target)
    anchors = {}
    for role, item in data['adjustable_components'].items():
        pset = item['pset']
        points = {}
        for key in ('InterfacePointMm','InstalledEndpointMm'):
            local = json.loads(pset[key])
            world = [sum(before_placement[i][j]*local[j] for j in range(3))+before_placement[i][3] for i in range(3)]
            points[key] = {'local_mm':local,'world_mm':world}
        anchors[role] = points
    blend = PRODUCT/'Geberit-151.116.11.1-TRAP01-official-detail-v4.blend'
    bpy.ops.wm.save_as_mainfile(filepath=str(blend), check_existing=False)
    data.update({'task':'TRAP01 approved official 2D installed configuration visible through cabinet occluders',
        'status':'persisted_reloaded_verified', 'verdict':'pending_rendered_verification',
        'preState':{'provider_scene':pre,'input_ifc':str(INPUT_IFC),'input_ifc_sha256':base.sha256(INPUT_IFC)},
        'execution':{'adapter':'installed public bonsai_bridge handlers in private Blender process; no new MCP server',
            'provider_file':bridge.__file__,'provider_sha256':base.sha256(Path(bridge.__file__)),
            'provider_version':bridge.bl_info.get('version'),'generator':'bpy.ops.bim.create_drawing','blender':bpy.app.version_string},
        'persistence':{'initial_save_as':save_as,'save_and_reload':persisted,'derived_ifc_sha256':base.sha256(IFC),'blend':str(blend)},
        'official_linework_audit':official_audit,
        'postState':{'reloaded_views':reloaded,'target_body_unchanged':True,'target_placement_unchanged':True,
            'target_body_graph_sha256':input_body_hash,'inherited_configuration_reference_points':anchors,
            'anchor_semantics_note':'Inherited controls are prior drawing-configuration reference points, not independently surveyed pipe centre lines. Annotation origin is verified against actual target ObjectPlacement; source shape is component-adjusted native DWG.',
            'old_outline_drawings_superseded':superseded,'formal_ifc_sha256':base.sha256(base.FORMAL_IFC)},
        'line_semantics':{'blue_solid':'complete official native 2D paths; rigid joints repositioned and straight pipes shortened to configured envelope',
            'blue_dashed':'official adjustable excess reference only; kept in audit, omitted from installed scene',
            'grey':'retained project context'},
        'user_scene_approval':'pending'})
    EVIDENCE.write_text(json.dumps(data,indent=2,ensure_ascii=False)+'\n')
    print('TRAP01_V4_SUCCESS '+str(EVIDENCE))

if __name__=='__main__':
    try:
        main()
    except Exception:
        traceback.print_exc()
        raise
