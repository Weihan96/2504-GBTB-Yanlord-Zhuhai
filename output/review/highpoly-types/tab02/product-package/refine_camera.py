"""Keep the pending front camera inside the terrace, away from the door wall."""
import json
import migrate_tab02 as task
from pure_product_package import graph_from_json,graph_to_json
recipe=json.loads(task.RECIPE.read_text())
graph=graph_from_json(recipe['graph'])
camera=graph.by_guid(recipe['views']['front']['drawing_guid'])
camera.ObjectPlacement.RelativePlacement.Location.Coordinates=(1750.0,3600.0,650.0)
for view,dims,loc in [('front',(2200.,1400.,1600.),(1750.,3600.,650.)),('side',(1600.,1400.,2300.),(900.,3950.,650.))]:
    d=graph.by_guid(recipe['views'][view]['drawing_guid']);d.ObjectPlacement.RelativePlacement.Location.Coordinates=loc
    block=next(x for x in graph.traverse(d.Representation) if x.is_a('IfcBlock'))
    block.XLength,block.YLength,block.ZLength=dims
    block.Position.Location.Coordinates=(-dims[0]/2,-dims[1]/2,-dims[2])
recipe['graph']=graph_to_json(graph)
recipe['svg_display_style']={'script':'migrate_tab02.py','function':'style_svg','operation':'Style only after real CreateDrawing; no geometry changes','background_projection':'fill:none;stroke:#aaa;stroke-width:0.06','background_cut':'fill:#eee;stroke:#999;stroke-width:0.10','approved_target':'stroke:#111820;stroke-width:0.18;fill:none'}
recipe['camera_policy']='Pending living-room north terrace views; Front at y3600 beyond y3200..3400 door wall and 50 mm import tolerance; Plan at z1700 above table top z680. No ceiling slab intersects the table XY. World +Z is screen up for elevations.'
task.write(task.RECIPE,recipe)
r=json.loads(task.REPORT.read_text())
for v in r['views']:
    v['camera_matrix']=task.pkg.placement(graph.by_guid(v['drawing_guid'])).tolist()
r['first_camera_attempt_rejected']='Front y3400 coincided with wall edge and Bonsai drawing tolerance cut the wall; moved only camera to y3600.'
task.write(task.REPORT,r)
