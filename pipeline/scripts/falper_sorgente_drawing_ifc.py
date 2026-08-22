#!/usr/bin/env python3
"""Approval-gated writer for Falper Sorgente drawing representations.

The pending approval record shipped with the repository is intentionally
rejected.  After explicit user approval, a reviewer must bind that approval
to the exact review manifest hash and approve plan/front/side.  The command
then writes a separate derived IFC; it never accepts the formal IFC path as
its output.
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import ifcopenshell
import ifcopenshell.guid
import ifcopenshell.util.element
import ifcopenshell.util.representation

from falper_sorgente_linework import (
    EXPECTED,
    PRODUCT_URL,
    ROOT,
    SCOPE,
    load_json,
    relative,
    sha256,
    write_json,
)


REPRESENTATIVE_GLOBAL_ID = "350tdaubr8QP3Cu2YMQZIN"
REQUIRED_VIEWS = {"plan", "front", "side"}
DEFAULT_MANIFEST = ROOT / "output/review/highpoly-types/falper-sorgente/manifest.json"
DEFAULT_CANDIDATE = ROOT / "output/review/highpoly-types/falper-sorgente/candidate-representations.json"
DEFAULT_APPROVAL = ROOT / "pipeline/decisions/falper-sorgente-drawing-approval.json"
DEFAULT_OFFICIAL_LINEWORK = (
    ROOT / "pipeline/decisions/falper-sorgente-official-dwg-linework.json"
)
EXPECTED_NATIVE_DWG_PATH_COUNTS = {"plan": 5, "front": 4, "side": 4}


def require_approval(approval: dict, manifest_path: Path) -> None:
    errors = []
    if approval.get("status") != "approved":
        errors.append("status must be approved")
    if approval.get("formal_ifc_write_allowed") is not True:
        errors.append("formal_ifc_write_allowed must be true")
    if set(approval.get("approved_views", [])) != REQUIRED_VIEWS:
        errors.append("approved_views must be exactly plan, front and side")
    if not str(approval.get("reviewer") or "").strip():
        errors.append("reviewer is required")
    if not re.fullmatch(r"\d{4}-\d{2}-\d{2}", str(approval.get("review_date") or "")):
        errors.append("review_date must be YYYY-MM-DD")
    if approval.get("candidate_manifest_sha256") != sha256(manifest_path):
        errors.append("candidate_manifest_sha256 does not match the reviewed manifest")
    if approval.get("profile_key") != "falper-sorgente" or approval.get("model_code") != "WFB":
        errors.append("approval identity must be falper-sorgente / WFB")
    if approval.get("scope") != SCOPE:
        errors.append("approval scope must remain family_reference_not_project_shop_drawing")
    if errors:
        raise RuntimeError("approval gate rejected IFC write: " + "; ".join(errors))


def representation_context(model, identifier: str, target_view: str):
    existing = ifcopenshell.util.representation.get_context(
        model, "Model", identifier, target_view
    )
    if existing is not None:
        return existing
    parent = next(
        context
        for context in model.by_type("IfcGeometricRepresentationContext", include_subtypes=False)
        if context.ContextType == "Model"
    )
    return model.create_entity(
        "IfcGeometricRepresentationSubContext",
        ContextIdentifier=identifier,
        ContextType="Model",
        ParentContext=parent,
        TargetScale=None,
        TargetView=target_view,
        UserDefinedTargetView=None,
    )


def curve_representation(model, context, identifier: str, view: str, paths):
    polylines = []
    for path in paths:
        if len(path) < 2:
            continue
        points = []
        for first, second in path:
            if view == "plan":
                coordinates = (float(first), float(second), 0.0)
            elif view == "front":
                coordinates = (float(first), 0.0, float(second))
            else:
                coordinates = (0.0, float(first), float(second))
            points.append(model.create_entity("IfcCartesianPoint", Coordinates=coordinates))
        polylines.append(model.create_entity("IfcPolyline", Points=points))
    curve_set = model.create_entity("IfcGeometricCurveSet", Elements=polylines)
    return model.create_entity(
        "IfcShapeRepresentation",
        ContextOfItems=context,
        RepresentationIdentifier=identifier,
        RepresentationType="GeometricCurveSet",
        Items=[curve_set],
    )


def official_native_dwg_view_paths(manifest: dict, linework: dict) -> dict:
    """Place official WFB DWG paths in the reviewed representative coordinates."""
    reference = manifest.get("official_reference", {})
    if (
        reference.get("source_kind") != "native_dwg"
        or reference.get("model_code") != "WFB"
        or reference.get("source_dwg_sha256") != EXPECTED["wfb_2d"]
        or linework.get("source_kind") != "native_dwg"
    ):
        raise RuntimeError("official native WFB DWG identity gate failed")
    variant = linework.get("variants", {}).get("WFB", {})
    if (
        variant.get("source_kind") != "native_dwg"
        or variant.get("source_dwg_sha256") != EXPECTED["wfb_2d"]
        or variant.get("scope") != SCOPE
    ):
        raise RuntimeError("official WFB linework register identity gate failed")
    bounds = manifest.get("local_bounds_mm", {})
    minimum, maximum = bounds.get("minimum"), bounds.get("maximum")
    if not minimum or not maximum or len(minimum) != 3 or len(maximum) != 3:
        raise RuntimeError("reviewed representative bounds are missing")
    center_x = (float(minimum[0]) + float(maximum[0])) / 2.0
    center_y = (float(minimum[1]) + float(maximum[1])) / 2.0
    base_z = float(minimum[2])
    plan = [
        [[center_x + float(point[0]), center_y + float(point[1])] for point in path]
        for path in variant["views"]["plan"]["paths_mm"]
    ]
    elevation = variant["views"]["elevation"]["paths_mm"]
    paths = {
        "plan": plan,
        "front": [
            [[center_x + float(point[0]), base_z + float(point[1])] for point in path]
            for path in elevation
        ],
        "side": [
            [[center_y + float(point[0]), base_z + float(point[1])] for point in path]
            for path in elevation
        ],
    }
    counts = {view: len(view_paths) for view, view_paths in paths.items()}
    if counts != EXPECTED_NATIVE_DWG_PATH_COUNTS:
        raise RuntimeError(f"official native DWG path count gate failed: {counts}")
    if variant["views"]["plan"].get("native_entity_types") != ["CIRCLE"] * 5:
        raise RuntimeError("official native DWG plan entity identity drifted")
    if variant["views"]["elevation"].get("native_entity_types") != [
        "SPLINE", "SPLINE", "LINE", "LINE"
    ]:
        raise RuntimeError("official native DWG elevation entity identity drifted")
    return paths


def representation_path_count(representation) -> int:
    return sum(len(curve_set.Elements) for curve_set in representation.Items)


def is_bonsai_drawing_body_representation(representation, target_view: str) -> bool:
    context = representation.ContextOfItems
    return (
        representation.RepresentationIdentifier == "Body"
        and context.ContextType == "Model"
        and context.ContextIdentifier == "Body"
        and context.TargetView == target_view
    )


def add_source_pset(
    model,
    product,
    product_type,
    approval,
    manifest_hash: str,
    linework_register_path: Path,
) -> None:
    values = {
        "SourceKind": "native_dwg",
        "Manufacturer": "Falper",
        "Family": "Sorgente",
        "ModelCode": "WFB",
        "SourceProductPage": PRODUCT_URL,
        "Source2DArchive": "https://falper.it/wp-content/uploads/2025/11/Lavabi-Freestanding-Autocad-2D.zip",
        "SourceTechnicalDwgArchive": "https://falper.it/wp-content/uploads/2025/11/Lavabi-Freestanding-Scheda-tecnica-DWG.zip",
        "SourceDwgPath": "drawings/evidence/FALPER-official-Sorgente-WFB-2D.dwg",
        "SourceDwgSha256": EXPECTED["wfb_2d"],
        "SourcePdfPath": "drawings/evidence/FALPER-official-Sorgente-WFA-WFB.pdf",
        "SourcePdfSha256": EXPECTED["pdf"],
        "EvidenceScope": SCOPE,
        "ReviewStatus": "APPROVED",
        "Reviewer": approval["reviewer"],
        "ReviewDate": approval["review_date"],
        "ApprovedCandidateManifestSha256": manifest_hash,
        "ApprovedRepresentationIdentifiers": "FalperWFBPlan;FalperWFBFront;FalperWFBSide",
        "BonsaiDrawingRepresentationIdentifiers": "Body/PLAN_VIEW;Body/ELEVATION_VIEW",
        "BonsaiPlanGeometrySource": "FalperWFBPlan official native DWG paths",
        "BonsaiElevationGeometrySource": "FalperWFBFront official native DWG paths",
        "OfficialLineColourHex": "#1677c8",
        "RepresentationSourceMapping": (
            "FalperWFBPlan=plan;FalperWFBFront=elevation;FalperWFBSide=elevation"
        ),
        "RepresentationGeometrySource": "official_native_dwg_paths_mm",
        "ProxyGeometryIncluded": "false",
        "NativeDwgPlanPathCount": "5",
        "NativeDwgFrontPathCount": "4",
        "NativeDwgSidePathCount": "4",
        "NativeDwgLineworkRegister": relative(linework_register_path),
        "NativeDwgLineworkRegisterSha256": sha256(linework_register_path),
    }
    properties = [
        model.create_entity(
            "IfcPropertySingleValue",
            Name=name,
            Description=None,
            NominalValue=model.create_entity("IfcText", str(value)),
            Unit=None,
        )
        for name, value in values.items()
    ]
    pset = model.create_entity(
        "IfcPropertySet",
        GlobalId=ifcopenshell.guid.new(),
        OwnerHistory=product.OwnerHistory,
        Name="Pset_FalperSorgenteDrawingSource",
        Description="Mechanically verifiable source and human approval for derived drawing representations",
        HasProperties=properties,
    )
    model.create_entity(
        "IfcRelDefinesByProperties",
        GlobalId=ifcopenshell.guid.new(),
        OwnerHistory=product.OwnerHistory,
        Name="Falper Sorgente drawing source properties",
        Description=None,
        RelatedObjects=[product, product_type],
        RelatingPropertyDefinition=pset,
    )


def add_document_association(model, product, product_type) -> None:
    reference = model.create_entity(
        "IfcDocumentReference",
        Location="drawings/evidence/FALPER-official-Sorgente-WFB-2D.dwg",
        Identification="FALPER-SORGENTE-WFB-NATIVE-DWG",
        Name="Falper official Sorgente WFB native DWG",
        Description=(
            f"SHA-256 {EXPECTED['wfb_2d']}; {SCOPE}; "
            "official product-family reference, not a project shop drawing"
        ),
        ReferencedDocument=None,
    )
    model.create_entity(
        "IfcRelAssociatesDocument",
        GlobalId=ifcopenshell.guid.new(),
        OwnerHistory=product.OwnerHistory,
        Name="Falper Sorgente WFB source association",
        Description="Official native DWG source for approved derived drawing representations",
        RelatedObjects=[product, product_type],
        RelatingDocument=reference,
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--candidate", type=Path, default=DEFAULT_CANDIDATE)
    parser.add_argument("--approval", type=Path, default=DEFAULT_APPROVAL)
    parser.add_argument("--official-linework", type=Path, default=DEFAULT_OFFICIAL_LINEWORK)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--report", type=Path)
    parser.add_argument("--apply", action="store_true")
    args = parser.parse_args()
    formal = args.input.resolve()
    output = args.output.resolve()
    manifest_path = args.manifest.resolve()
    linework_path = args.official_linework.resolve()
    if not args.apply:
        raise RuntimeError("IFC write requires the explicit --apply flag")
    if output == formal:
        raise RuntimeError("formal IFC cannot be the output; write a separate derived IFC")
    if output.exists():
        raise RuntimeError(f"refusing to overwrite existing output: {output}")
    formal_hash = sha256(formal)
    if formal_hash != EXPECTED["formal_ifc"]:
        raise RuntimeError("formal IFC hash mismatch")
    manifest = load_json(manifest_path)
    candidate = load_json(args.candidate.resolve())
    official_linework = load_json(linework_path)
    approval = load_json(args.approval.resolve())
    require_approval(approval, manifest_path)
    if (
        manifest.get("formal_ifc_sha256") != formal_hash
        or manifest.get("approved_for_drawing_ifc") is not False
        or manifest.get("review_status") != "visual_review_pending"
    ):
        raise RuntimeError("candidate manifest is not the expected pending review artifact")
    if candidate.get("formal_ifc_write_allowed") is not False or set(candidate["views"]) != REQUIRED_VIEWS:
        raise RuntimeError("candidate representation identity drifted")
    if manifest.get("official_reference", {}).get("linework_register_sha256") != sha256(linework_path):
        raise RuntimeError("official linework register hash does not match reviewed manifest")
    official_paths = official_native_dwg_view_paths(manifest, official_linework)

    model = ifcopenshell.open(formal)
    product = model.by_guid(REPRESENTATIVE_GLOBAL_ID)
    if product is None:
        raise RuntimeError("Falper representative product is missing")
    product_type = next(relation.RelatingType for relation in product.IsTypedBy)
    if product_type.Name != "BS01":
        raise RuntimeError("Falper representative type identity drifted")
    contexts = {
        "plan": representation_context(model, "FalperWFBPlan", "PLAN_VIEW"),
        "front": representation_context(model, "FalperWFBFront", "ELEVATION_VIEW"),
        "side": representation_context(model, "FalperWFBSide", "ELEVATION_VIEW"),
    }
    bonsai_contexts = {
        "plan": representation_context(model, "Body", "PLAN_VIEW"),
        "front": representation_context(model, "Body", "ELEVATION_VIEW"),
    }
    # Bonsai Drawing only considers Model/Body subcontexts. Replace the legacy
    # PLAN_VIEW Body projection and add an ELEVATION_VIEW Body projection so
    # the real IFC Drawing cameras select the approved native-DWG linework.
    representations = [
        representation
        for representation in product.Representation.Representations
        if not is_bonsai_drawing_body_representation(representation, "PLAN_VIEW")
        and not is_bonsai_drawing_body_representation(representation, "ELEVATION_VIEW")
    ]
    for view in ("plan", "front", "side"):
        representations.append(
            curve_representation(
                model,
                contexts[view],
                contexts[view].ContextIdentifier,
                view,
                official_paths[view],
            )
        )
    representations.extend(
        [
            curve_representation(
                model,
                bonsai_contexts["plan"],
                "Body",
                "plan",
                official_paths["plan"],
            ),
            curve_representation(
                model,
                bonsai_contexts["front"],
                "Body",
                "front",
                official_paths["front"],
            ),
        ]
    )
    product.Representation.Representations = representations
    add_document_association(model, product, product_type)
    add_source_pset(
        model,
        product,
        product_type,
        approval,
        sha256(manifest_path),
        linework_path,
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    model.write(output)
    if sha256(formal) != formal_hash:
        output.unlink(missing_ok=True)
        raise RuntimeError("formal IFC bytes changed during derived write")
    derived = ifcopenshell.open(output)
    derived_product = derived.by_guid(REPRESENTATIVE_GLOBAL_ID)
    derived_type = next(relation.RelatingType for relation in derived_product.IsTypedBy)
    identifiers = {
        representation.RepresentationIdentifier
        for representation in derived_product.Representation.Representations
    }
    if not {"FalperWFBPlan", "FalperWFBFront", "FalperWFBSide"}.issubset(identifiers):
        output.unlink(missing_ok=True)
        raise RuntimeError("derived IFC representation post-write gate failed")
    bonsai_standard_representations = {
        view: [
            representation
            for representation in derived_product.Representation.Representations
            if is_bonsai_drawing_body_representation(representation, target_view)
        ]
        for view, target_view in {
            "plan": "PLAN_VIEW",
            "front": "ELEVATION_VIEW",
        }.items()
    }
    if any(len(items) != 1 for items in bonsai_standard_representations.values()):
        output.unlink(missing_ok=True)
        raise RuntimeError("derived IFC Bonsai Drawing representation gate failed")
    if {
        view: representation_path_count(items[0])
        for view, items in bonsai_standard_representations.items()
    } != {"plan": 5, "front": 4}:
        output.unlink(missing_ok=True)
        raise RuntimeError("derived IFC Bonsai Drawing path count gate failed")
    representation_views = {
        "FalperWFBPlan": "plan",
        "FalperWFBFront": "front",
        "FalperWFBSide": "side",
    }
    native_dwg_path_counts = {}
    for identifier, view in representation_views.items():
        representation = next(
            item
            for item in derived_product.Representation.Representations
            if item.RepresentationIdentifier == identifier
        )
        count = representation_path_count(representation)
        native_dwg_path_counts[view] = count
        if count != EXPECTED_NATIVE_DWG_PATH_COUNTS[view]:
            output.unlink(missing_ok=True)
            raise RuntimeError(f"derived IFC native DWG path count gate failed: {view}={count}")
    document_relations = [
        inverse
        for inverse in derived.get_inverse(derived_product)
        if inverse.is_a("IfcRelAssociatesDocument")
        and inverse.RelatingDocument.is_a("IfcDocumentReference")
        and inverse.RelatingDocument.Identification == "FALPER-SORGENTE-WFB-NATIVE-DWG"
    ]
    if (
        len(document_relations) != 1
        or derived_type not in document_relations[0].RelatedObjects
        or EXPECTED["wfb_2d"] not in (document_relations[0].RelatingDocument.Description or "")
    ):
        output.unlink(missing_ok=True)
        raise RuntimeError("derived IFC document association post-write gate failed")
    source_pset = ifcopenshell.util.element.get_pset(
        derived_product, "Pset_FalperSorgenteDrawingSource"
    ) or {}
    expected_properties = {
        "SourceKind": "native_dwg",
        "ModelCode": "WFB",
        "SourceDwgSha256": EXPECTED["wfb_2d"],
        "SourcePdfSha256": EXPECTED["pdf"],
        "EvidenceScope": SCOPE,
        "ReviewStatus": "APPROVED",
        "Reviewer": approval["reviewer"],
        "ReviewDate": approval["review_date"],
        "ApprovedCandidateManifestSha256": sha256(manifest_path),
        "ApprovedRepresentationIdentifiers": "FalperWFBPlan;FalperWFBFront;FalperWFBSide",
        "BonsaiDrawingRepresentationIdentifiers": "Body/PLAN_VIEW;Body/ELEVATION_VIEW",
        "BonsaiPlanGeometrySource": "FalperWFBPlan official native DWG paths",
        "BonsaiElevationGeometrySource": "FalperWFBFront official native DWG paths",
        "OfficialLineColourHex": "#1677c8",
        "RepresentationSourceMapping": (
            "FalperWFBPlan=plan;FalperWFBFront=elevation;FalperWFBSide=elevation"
        ),
        "RepresentationGeometrySource": "official_native_dwg_paths_mm",
        "ProxyGeometryIncluded": "false",
        "NativeDwgPlanPathCount": "5",
        "NativeDwgFrontPathCount": "4",
        "NativeDwgSidePathCount": "4",
        "NativeDwgLineworkRegister": relative(linework_path),
        "NativeDwgLineworkRegisterSha256": sha256(linework_path),
    }
    if any(source_pset.get(name) != value for name, value in expected_properties.items()):
        output.unlink(missing_ok=True)
        raise RuntimeError("derived IFC source property post-write gate failed")
    result = {
        "schema_version": 1,
        "profile_key": "falper-sorgente",
        "model_code": "WFB",
        "representative_global_id": REPRESENTATIVE_GLOBAL_ID,
        "formal_ifc": relative(formal),
        "formal_ifc_sha256": formal_hash,
        "formal_ifc_bytes_unchanged": True,
        "derived_ifc": relative(output),
        "derived_ifc_sha256": sha256(output),
        "candidate_manifest": relative(manifest_path),
        "candidate_manifest_sha256": sha256(manifest_path),
        "approval_record": relative(args.approval.resolve()),
        "approval_record_sha256": sha256(args.approval.resolve()),
        "approval": {
            "status": approval["status"],
            "reviewer": approval["reviewer"],
            "review_date": approval["review_date"],
            "approved_views": approval["approved_views"],
            "scope": approval["scope"],
        },
        "representations": sorted(identifiers),
        "representation_path_counts": native_dwg_path_counts,
        "bonsai_drawing_representation_path_counts": {"plan": 5, "front": 4},
        "bonsai_drawing_representation_selection": "Model/Body PLAN_VIEW and ELEVATION_VIEW",
        "representation_source_mapping": source_pset["RepresentationSourceMapping"],
        "representation_geometry_source": source_pset["RepresentationGeometrySource"],
        "proxy_geometry_included": False,
        "official_linework_register": relative(linework_path),
        "official_linework_register_sha256": sha256(linework_path),
        "source_document_association": "IfcDocumentReference/IfcRelAssociatesDocument",
        "source_property_set": "Pset_FalperSorgenteDrawingSource",
        "source_kind": source_pset["SourceKind"],
        "source_dwg_sha256": source_pset["SourceDwgSha256"],
        "source_pdf_sha256": source_pset["SourcePdfSha256"],
        "pass": True,
    }
    if args.report:
        write_json(args.report.resolve(), result)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
