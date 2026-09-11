"""Verify the complete external runtime adapter twice, before Bonsai import."""
import json,ifcopenshell
import migrate_tab02 as task
from pure_product_package import attach_recipe
from repair_runtime_metadata import repair
formal=ifcopenshell.open(str(task.FORMAL));project=ifcopenshell.open(str(task.FORMAL));pure=ifcopenshell.open(str(task.SINGLE));recipe=json.loads(task.RECIPE.read_text())
attach_recipe(project,pure,recipe);repair(project,formal)
count=len(list(project));attach_recipe(project,pure,recipe);second_repairs=repair(project,formal)
assert len(list(project))==count
task.write(task.OUT/'runtime-idempotence.json',{'runtime':'attach_recipe + recipe runtime_metadata_adapter','first_then_second_entity_delta':0,'second_repairs':second_repairs,'formal_sha256':task.pkg.sha256(task.FORMAL),'pure_sha256':task.pkg.sha256(task.SINGLE),'scene_recipe_sha256':task.pkg.sha256(task.RECIPE),'save_boundary':'in-memory original formal copy only','bonsai_camera_round_trip_normalization':'Bonsai may normalize orthographic camera block dimensions during import/save; saved scene is validated against real SVG separately.'})
print('Runtime adapter second attachment: 0 new entities')
