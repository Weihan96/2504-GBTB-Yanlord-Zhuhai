"""Product Body view roles survive context reuse, insertion and IFC reload.

IfcShapeAspect identifies each directional Body, not another annotation object.
The camera normal is in product-local coordinates, never thumbnail/world space.
"""
import json
import math
import numpy as np

ROLES = ("Plan", "Front", "Side")
MARKER = "review_library_body_view"


def views(product):
    result = {}
    shape = product.Representation
    if not shape:
        return result
    for aspect in getattr(shape, "HasShapeAspects", ()):
        try:
            data = json.loads(aspect.Description or "{}")
        except (ValueError, TypeError):
            continue
        if not isinstance(data, dict) or data.get(MARKER) != 1:
            continue
        role = aspect.Name
        if role not in ROLES or role in result or len(aspect.ShapeRepresentations) != 1:
            raise ValueError("单品 Body 视向记录重复或损坏")
        rep = aspect.ShapeRepresentations[0]
        if rep not in shape.Representations or rep.RepresentationIdentifier != "Body" or rep.ContextOfItems.ContextIdentifier != "Body":
            raise ValueError("单品视向必须指向本产品的 Body 表示")
        normal = np.asarray(data["camera_normal_local"], dtype=float)
        if normal.shape != (3,) or not np.isfinite(normal).all() or not np.isclose(np.linalg.norm(normal), 1):
            raise ValueError("无效的产品局部相机法向")
        result[role] = (rep, normal)
    if result and set(result) != set(ROLES):
        raise ValueError("单品 Body 缺少 Plan、Front 或 Side")
    return result


def validate(product):
    directional = views(product)
    assert set(directional) == set(ROLES)
    reps = product.Representation.Representations
    assert len(reps) == 4
    assert all(r.RepresentationIdentifier == "Body" and r.ContextOfItems.ContextIdentifier == "Body" for r in reps)
    assert len([r for r in reps if r.ContextOfItems.TargetView == "MODEL_VIEW"]) == 1
    assert directional["Plan"][0].ContextOfItems.TargetView == "PLAN_VIEW"
    assert all(directional[role][0].ContextOfItems.TargetView == "ELEVATION_VIEW" for role in ("Front", "Side"))
    return directional


def choose_role(directional, camera_normal_local, target_view, tolerance_degrees=5):
    """Canonical views only; oblique/section/model views use the 3D Body.

    Opposite viewing directions use the same geometric projection plane (with
    normal camera mirroring); no additional unreviewed back elevation is claimed.
    """
    if target_view in ("PLAN_VIEW", "REFLECTED_PLAN_VIEW"):
        roles = ("Plan",)
    elif target_view == "ELEVATION_VIEW":
        roles = ("Front", "Side")
    else:
        return None
    normal = np.asarray(camera_normal_local, dtype=float)
    if normal.shape != (3,) or not np.isfinite(normal).all() or np.linalg.norm(normal) < 1e-8:
        raise ValueError("无效的相机法向")
    normal = normal / np.linalg.norm(normal)
    ranked = sorted(((abs(float(np.dot(normal, directional[role][1]))), role) for role in roles), reverse=True)
    return ranked[0][1] if ranked[0][0] >= math.cos(math.radians(tolerance_degrees)) else None
