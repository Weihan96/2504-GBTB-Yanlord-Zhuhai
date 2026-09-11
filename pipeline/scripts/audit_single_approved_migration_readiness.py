"""Read-only readiness audit; this never creates or changes an IFC."""
import json
from pathlib import Path
import ifcopenshell
import numpy as np
import review_product_package as pkg

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "output/review/approved-product-library"
FORMAL = ROOT / "2504 GBTB Yanlord Zhuhai.ifc"
EXPECTED = "7a50b87e8f48a7c2bfbaa9f1dabdfca7155fa684b2325c66aa1ae15c8ab0a25c"


def main():
    assert pkg.sha256(FORMAL) == EXPECTED
    formal = ifcopenshell.open(str(FORMAL))
    scope = json.loads((OUT / "acceptance-scope.json").read_text())
    rows = []
    for key in scope["single_product_accepted_scene_not_final"]:
        directory = ROOT / "output/review/highpoly-types" / key
        approval_path = ROOT / "pipeline/decisions" / f"{key}-drawing-approval.json"
        approval = json.loads(approval_path.read_text())
        paths = list(directory.glob("*bonsai-isolated.ifc"))
        assert len(paths) == 1, (key, paths)
        source_path = paths[0]
        source = ifcopenshell.open(str(source_path))
        products = [e for e in source.by_type("IfcElement") if e.Representation]
        assert len(products) == 1, (key, [e.GlobalId for e in products])
        target = products[0]
        baseline = formal.by_guid(target.GlobalId)
        reps = list(target.Representation.Representations)
        identifiers = [r.RepresentationIdentifier for r in reps]
        views = [name for name in identifiers if any(v in (name or "").lower() for v in ("plan", "front", "side"))]
        body_matches = pkg.body_fingerprint(target) == pkg.body_fingerprint(baseline)
        placement_error = float(np.max(np.abs(pkg.placement(target) - pkg.placement(baseline))))
        rows.append({"id": key, "global_id": target.GlobalId, "approval": pkg.record(approval_path),
            "single_views_approved": approval["approved_views"],
            "recorded_derived_ifc_write_allowed": approval.get("derived_ifc_write_allowed"),
            "approval_scope": approval.get("scope"), "approval_note": approval.get("note"),
            "source_ifc": pkg.record(source_path), "representation_identifiers": identifiers,
            "existing_view_representations": views,
            "body_matches_formal": body_matches, "placement_matrix_max_error_project_units": placement_error,
            "unit_scale_matches_formal": pkg.unit_util.calculate_unit_scale(source) == pkg.unit_util.calculate_unit_scale(formal),
            "work_if_included": "Verify and package existing view geometry, create and validate external scene recipe"
                if len(views) == 3 else "First write approved 2D geometry into a new pure IFC, then create and validate external scene recipe",
            "scope_confirmation_pending": True, "ifc_written": False})
    result = {"status": "read_only_readiness_audit", "products": rows,
        "single_approval_not_promoted_to_scene_approval": True,
        "formal_ifc_sha256_before": EXPECTED, "formal_ifc_sha256_after": pkg.sha256(FORMAL),
        "ifc_writes_performed": 0,
        "question": "Should the seven single-SVG-approved products also be migrated and included with scene approval explicitly pending?"}
    assert result["formal_ifc_sha256_after"] == EXPECTED
    (OUT / "single-approved-readiness.json").write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n")
    print(json.dumps([{k: r[k] for k in ("id", "representation_identifiers", "body_matches_formal", "placement_matrix_max_error_project_units", "recorded_derived_ifc_write_allowed")} for r in rows], ensure_ascii=False))


if __name__ == "__main__":
    main()
