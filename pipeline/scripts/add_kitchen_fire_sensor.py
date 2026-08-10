#!/usr/bin/env python3
"""Add the confirmed R04 kitchen fire-sensor location to the formal IFC."""

from pathlib import Path
import hashlib
import os
import tempfile

import ifcopenshell
from ifcopenshell.guid import new as new_guid


ROOT = Path(__file__).resolve().parents[2]
IFC = ROOT / "2504 GBTB Yanlord Zhuhai.ifc"
EXPECTED = "7521c09991f3d0c7b7d91ca2324fd55ad961d8e32e9e3e9a9777a4cc19b06e81"
SPACE_GUID = "2fhEbDfK1EkhJwlPikNm$b"
STOREY_NAME = "KITCHEN"
POSITION_MM = (1800.0, -4576.0, 2400.0)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def add_sensor(model):
    space = model.by_guid(SPACE_GUID)
    if space is None or space.is_a() != "IfcSpace":
        raise RuntimeError("R04 space not found")
    storey = next((s for s in model.by_type("IfcBuildingStorey") if s.Name == STOREY_NAME), None)
    if storey is None:
        raise RuntimeError("KITCHEN storey not found")

    existing = [s for s in model.by_type("IfcSensor") if s.Name == "A106-FIRE-R04"]
    if existing:
        raise RuntimeError("A106-FIRE-R04 already exists")

    owner_history = model.by_type("IfcOwnerHistory")[0]
    sensor = model.create_entity(
        "IfcSensor",
        GlobalId=new_guid(),
        OwnerHistory=owner_history,
        Name="A106-FIRE-R04",
        Description=(
            "厨房火灾探测器已确认点位；最终感温/感烟/复合类型、产品、数量、"
            "供电通信及厂家安装条件待消防/设备方确认。"
        ),
        ObjectType="Confirmed location only; product pending",
        ObjectPlacement=None,
        Representation=None,
        PredefinedType="FIRESENSOR",
    )

    # IFC project units are millimetres; the placement is intentionally body-free.
    point = model.create_entity("IfcCartesianPoint", Coordinates=POSITION_MM)
    axis = model.create_entity("IfcAxis2Placement3D", Location=point)
    sensor.ObjectPlacement = model.create_entity(
        "IfcLocalPlacement", PlacementRelTo=storey.ObjectPlacement, RelativePlacement=axis
    )

    storey_rel = next(
        (r for r in model.by_type("IfcRelContainedInSpatialStructure") if r.RelatingStructure == storey),
        None,
    )
    if storey_rel is None:
        raise RuntimeError("KITCHEN containment relation not found")
    storey_rel.RelatedElements = tuple(storey_rel.RelatedElements) + (sensor,)

    model.create_entity(
        "IfcRelReferencedInSpatialStructure",
        GlobalId=new_guid(),
        OwnerHistory=owner_history,
        RelatedElements=(sensor,),
        RelatingStructure=space,
    )
    return sensor


def main():
    before = sha256(IFC)
    if before != EXPECTED:
        raise RuntimeError(f"formal IFC hash changed before write: {before}")
    model = ifcopenshell.open(IFC)
    before_counts = {t: len(model.by_type(t)) for t in ("IfcProduct", "IfcSensor", "IfcProductDefinitionShape")}
    sensor = add_sensor(model)
    fd, temp_name = tempfile.mkstemp(prefix="formal-ifc-", suffix=".ifc", dir=IFC.parent)
    os.close(fd)
    temp = Path(temp_name)
    try:
        model.write(temp)
        os.replace(temp, IFC)
    finally:
        if temp.exists():
            temp.unlink()
    after = sha256(IFC)
    reloaded = ifcopenshell.open(IFC)
    check = reloaded.by_guid(sensor.GlobalId)
    if check is None or check.is_a() != "IfcSensor" or check.PredefinedType != "FIRESENSOR":
        raise RuntimeError("saved sensor failed semantic reload check")
    if check.Representation is not None:
        raise RuntimeError("sensor unexpectedly has a body representation")
    print({"guid": sensor.GlobalId, "hash_before": before, "hash_after": after,
           "counts_before": before_counts,
           "counts_after": {t: len(reloaded.by_type(t)) for t in before_counts},
           "position_mm": POSITION_MM, "space": SPACE_GUID})


if __name__ == "__main__":
    main()
