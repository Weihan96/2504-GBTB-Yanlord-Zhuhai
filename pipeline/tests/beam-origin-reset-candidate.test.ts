import { expect, test } from "bun:test";
import { spawnSync } from "node:child_process";

function runPython(expression: string) {
  return spawnSync("python3", ["-c", expression], {
    cwd: process.cwd(),
    encoding: "utf8",
  });
}

test("beam origin reset transforms local coordinates into the new placement", () => {
  const result = runPython(
    `import json,sys,numpy as np; sys.path.insert(0,"pipeline/scripts"); from beam_origin_reset_candidate import transform_point; m=np.identity(4); m[:3,3]=[10,-20,30]; print(json.dumps(transform_point((1,2,3),m)))`,
  );
  expect(result.status).toBe(0);
  expect(JSON.parse(result.stdout)).toEqual([11, -18, 33]);
});

test("beam origin reset measures points on triangle faces and edges", () => {
  const result = runPython(
    `import json,sys; sys.path.insert(0,"pipeline/scripts"); from beam_origin_reset_candidate import point_triangle_distance; print(json.dumps([point_triangle_distance((0.25,0.25,0),(0,0,0),(1,0,0),(0,1,0)),point_triangle_distance((2,0,0),(0,0,0),(1,0,0),(0,1,0)),point_triangle_distance((0.25,0.25,0.5),(0,0,0),(1,0,0),(0,1,0))]))`,
  );
  expect(result.status).toBe(0);
  const distances = JSON.parse(result.stdout);
  expect(distances[0]).toBeCloseTo(0, 9);
  expect(distances[1]).toBeCloseTo(1, 9);
  expect(distances[2]).toBeCloseTo(0.5, 9);
});

test("beam origin reset accepts an omitted swept-solid position", () => {
  const result = runPython(
    `import sys,numpy as np,ifcopenshell; sys.path.insert(0,"pipeline/scripts"); from beam_origin_reset_candidate import transform_body_item; f=ifcopenshell.file(schema="IFC4"); d=f.createIfcDirection((0.,0.,1.)); p=f.createIfcRectangleProfileDef("AREA",None,None,100.,200.); s=f.createIfcExtrudedAreaSolid(p,None,d,300.); m=np.identity(4); m[:3,3]=[10,20,30]; print(transform_body_item(f,s,m),s.Position.Location.Coordinates,len(list(f)))`,
  );
  expect(result.status).toBe(0);
  expect(result.stdout.trim()).toContain("2 (10.0, 20.0, 30.0)");
});

test("beam origin reset counter-transforms tessellated and box representations", () => {
  const result = runPython(
    `import sys,numpy as np,ifcopenshell; sys.path.insert(0,"pipeline/scripts"); from beam_origin_reset_candidate import transform_body_item; f=ifcopenshell.file(schema="IFC4"); pts=f.createIfcCartesianPointList3D(((0.,0.,0.),(1.,0.,0.),(0.,1.,0.))); face=f.createIfcIndexedPolygonalFace((1,2,3)); mesh=f.createIfcPolygonalFaceSet(pts,False,(face,),None); corner=f.createIfcCartesianPoint((0.,0.,0.)); box=f.createIfcBoundingBox(corner,1.,2.,3.); m=np.identity(4); m[:3,3]=[10,20,30]; print(transform_body_item(f,mesh,m),pts.CoordList); print(transform_body_item(f,box,m),corner.Coordinates)`,
  );
  expect(result.status).toBe(0);
  expect(result.stdout).toContain("((10.0, 20.0, 30.0), (11.0, 20.0, 30.0)");
  expect(result.stdout).toContain("(10.0, 20.0, 30.0)");
});

test("beam origin reset counter-transforms an exclusive mapped item", () => {
  const result = runPython(
    `import sys,numpy as np,ifcopenshell; sys.path.insert(0,"pipeline/scripts"); from beam_origin_reset_candidate import transform_body_item; f=ifcopenshell.file(schema="IFC4"); source_origin=f.createIfcAxis2Placement3D(f.createIfcCartesianPoint((0.,0.,0.))); context=f.createIfcGeometricRepresentationContext(None,"Model",3,1e-5,source_origin,None); mapped_rep=f.createIfcShapeRepresentation(context,"Body","Tessellation",()); source=f.createIfcRepresentationMap(source_origin,mapped_rep); local=f.createIfcCartesianPoint((1.,2.,3.)); target=f.createIfcCartesianTransformationOperator3D(None,None,local,1.,None); item=f.createIfcMappedItem(source,target); holder=f.createIfcShapeRepresentation(context,"Body","MappedRepresentation",(item,)); m=np.identity(4); m[:3,3]=[10,20,30]; print(transform_body_item(f,item,m),local.Coordinates,len(f.get_inverse(target)),len(f.get_inverse(local)))`,
  );
  expect(result.status).toBe(0);
  expect(result.stdout.trim()).toBe("0 (11.0, 22.0, 33.0) 1 1");
});

test("beam origin reset counter-transforms direct Curve3D representations", () => {
  const result = runPython(
    `import sys,numpy as np,ifcopenshell; sys.path.insert(0,"pipeline/scripts"); from beam_origin_reset_candidate import transform_body_item; f=ifcopenshell.file(schema="IFC4"); pts=f.createIfcCartesianPointList3D(((0.,0.,0.),(1.,2.,3.))); curve=f.createIfcIndexedPolyCurve(pts,None,False); curves=f.createIfcGeometricCurveSet((curve,)); m=np.identity(4); m[:3,3]=[10,20,30]; print(transform_body_item(f,curves,m),pts.CoordList)`,
  );
  expect(result.status).toBe(0);
  expect(result.stdout.trim()).toBe(
    "0 ((10.0, 20.0, 30.0), (11.0, 22.0, 33.0))",
  );
});

test("beam origin reset transforms an in-product shared coordinate list once", () => {
  const result = runPython(
    `import sys,numpy as np,ifcopenshell; sys.path.insert(0,"pipeline/scripts"); from beam_origin_reset_candidate import transform_body_item; f=ifcopenshell.file(schema="IFC4"); pts=f.createIfcCartesianPointList3D(((0.,0.,0.),(1.,0.,0.),(0.,1.,0.))); face=f.createIfcIndexedPolygonalFace((1,2,3)); a=f.createIfcPolygonalFaceSet(pts,False,(face,),None); b=f.createIfcPolygonalFaceSet(pts,False,(face,),None); m=np.identity(4); m[:3,3]=[10,20,30]; done=set(); allowed={a.id(),b.id(),pts.id(),face.id()}; print(transform_body_item(f,a,m,done,allowed),transform_body_item(f,b,m,done,allowed),pts.CoordList)`,
  );
  expect(result.status).toBe(0);
  expect(result.stdout.trim()).toBe(
    "0 0 ((10.0, 20.0, 30.0), (11.0, 20.0, 30.0), (10.0, 21.0, 30.0))",
  );
});
