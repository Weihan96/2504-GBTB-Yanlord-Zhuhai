#!/usr/bin/env python3
"""Approval-gated writer for geometry-derived Hima drawing representations.

The official Hima product page and technical sheet establish manufacturer and
family identity only. The checked-in candidate linework is explicitly derived
from the single representative IFC Body because the published official 2D DWG
has not been acquired. The formal authoritative IFC is never a valid output.
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import ifcopenshell
import ifcopenshell.guid
import ifcopenshell.util.element

from falper_sorgente_linework import ROOT, load_json, relative, sha256, write_json


REPRESENTATIVE_GLOBAL_ID = "2xmcLzu1rDTeMzRuNxPDyE"
IFC_TYPE_NAME = "HIMA01"
PROFILE_KEY = "hima01"
SOURCE_KIND = "geometry_derived_simplified_proxy"
SOURCE_LABEL_ZH = "基于原始高模几何生成的简化图纸表达"
PRODUCT_PAGE = "https://www.poliform.it/en/products/hima/"
TECHNICAL_SHEET = "https://www.poliform.it/assets/pdf/200467-hima-poliform-en.pdf"
NEWS_TECHNICAL_PUBLICATION = "https://www.poliform.it/assets/2022/05/Poliform_News_2022-2.pdf"
SCOPE = "manufacturer family identity and nominal dimensions only; not official CAD geometry and not a project shop drawing"
REQUIRED_VIEWS = {"plan", "front", "side"}
REPRESENTATIONS = {
    "plan": ("Hima01Plan", "PLAN_VIEW"),
    "front": ("Hima01Front", "ELEVATION_VIEW"),
    "side": ("Hima01Side", "ELEVATION_VIEW"),
}
EXPECTED_PATH_COUNTS = {"plan": 8, "front": 17, "side": 27}
FORMAL_SHA256 = "7a50b87e8f48a7c2bfbaa9f1dabdfca7155fa684b2325c66aa1ae15c8ab0a25c"
PRODUCT_DIR = ROOT / "output/review/highpoly-types/hima01"
DEFAULT_MANIFEST = PRODUCT_DIR / "manifest.json"
DEFAULT_CANDIDATE = PRODUCT_DIR / "candidate-representations.json"
DEFAULT_ACCESS_RECORD = PRODUCT_DIR / "official-source/source-access-record.json"
PAGE_EVIDENCE = PRODUCT_DIR / "official-source/official-product-page-evidence.json"
DEFAULT_APPROVAL = ROOT / "pipeline/decisions/hima01-drawing-approval.json"
PSET_NAME = "Pset_Hima01DrawingSource"
DOCUMENT_ID_PREFIX = "POLIFORM-HIMA-"
CLOSE_REPRESENTATION_PATHS = True
OFFICIAL_CAD_USED = False
OFFICIAL_CAD_GEOMETRY_INCLUDED = False


def require_approval(approval: dict, manifest_path: Path) -> None:
    errors = []
    if approval.get("status") != "approved":
        errors.append("status must be approved")
    if approval.get("derived_ifc_write_allowed") is not True:
        errors.append("derived_ifc_write_allowed must be true")
    if approval.get("formal_authoritative_ifc_write_allowed") is not False:
        errors.append("formal_authoritative_ifc_write_allowed must remain false")
    if set(approval.get("approved_views", [])) != REQUIRED_VIEWS:
        errors.append("approved_views must be exactly plan, front and side")
    if not str(approval.get("reviewer") or "").strip():
        errors.append("reviewer is required")
    if not re.fullmatch(r"\d{4}-\d{2}-\d{2}", str(approval.get("review_date") or "")):
        errors.append("review_date must be YYYY-MM-DD")
    if approval.get("candidate_manifest_sha256") != sha256(manifest_path):
        errors.append("candidate_manifest_sha256 does not match the reviewed manifest")
    if approval.get("profile_key") != PROFILE_KEY or approval.get("ifc_type_name") != IFC_TYPE_NAME:
        errors.append(f"approval identity must be {PROFILE_KEY} / {IFC_TYPE_NAME}")
    if approval.get("scope") != SCOPE:
        errors.append("approval scope must preserve the geometry-derived and non-shop-drawing limitations")
    if not str(approval.get("approval_evidence") or "").strip():
        errors.append("approval_evidence is required")
    if errors:
        raise RuntimeError("approval gate rejected IFC write: " + "; ".join(errors))


def representation_context(model, identifier: str, target_view: str):
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
        closed = list(path)
        if CLOSE_REPRESENTATION_PATHS and closed[0] != closed[-1]:
            closed.append(closed[0])
        points = []
        for first, second in closed:
            if view == "plan":
                coordinates = (float(first), float(second), 0.0)
            elif view == "front":
                coordinates = (float(first), 0.0, float(second))
            else:
                coordinates = (0.0, float(first), float(second))
            points.append(model.create_entity("IfcCartesianPoint", Coordinates=coordinates))
        polylines.append(model.create_entity("IfcPolyline", Points=points))
    if len(polylines) != EXPECTED_PATH_COUNTS[view]:
        raise RuntimeError(f"approved {view} geometry-derived path count drifted")
    curve_set = model.create_entity("IfcGeometricCurveSet", Elements=polylines)
    return model.create_entity(
        "IfcShapeRepresentation",
        ContextOfItems=context,
        RepresentationIdentifier=identifier,
        RepresentationType="GeometricCurveSet",
        Items=[curve_set],
    )


def candidate_paths(candidate: dict, access: dict) -> dict:
    if (
        candidate.get("profile_key") != PROFILE_KEY
        or candidate.get("representative_global_id") != REPRESENTATIVE_GLOBAL_ID
        or candidate.get("source_kind") != SOURCE_KIND
        or candidate.get("source_label_zh") != SOURCE_LABEL_ZH
        or candidate.get("official_cad_used") is not False
        or candidate.get("third_party_cad_used") is not False
        or candidate.get("formal_ifc_write_allowed") is not False
        or candidate.get("review_status") != "visual_review_pending"
    ):
        raise RuntimeError("pending Hima candidate source gate failed")
    drawing_source = access.get("drawing_geometry_source", {})
    if (
        access.get("official_2d_dwg", {}).get("acquired") is not False
        or drawing_source.get("source_kind") != SOURCE_KIND
        or drawing_source.get("official_cad_used") is not False
        or drawing_source.get("third_party_cad_used") is not False
        or access.get("scope") != SCOPE
        or access.get("official_product_page_evidence", {}).get("sha256") != sha256(PAGE_EVIDENCE)
        or access.get("official_dimension_cross_check", {}).get("pass") is not True
    ):
        raise RuntimeError("Hima official-source access record gate failed")
    paths = {}
    for view in REQUIRED_VIEWS:
        item = candidate.get("views", {}).get(view, {})
        if item.get("source_kind") != SOURCE_KIND or item.get("official_cad_paths_mm") != []:
            raise RuntimeError(f"Hima {view} linework must remain geometry-derived with no official CAD paths")
        paths[view] = item.get("proxy_paths_mm", [])
        if len(paths[view]) != EXPECTED_PATH_COUNTS[view]:
            raise RuntimeError(f"Hima {view} geometry-derived path-count gate failed")
    return paths


def add_document_associations(model, product, product_type, access_path):
    documents = (
        (
            "POLIFORM-HIMA-OFFICIAL-PRODUCT-PAGE",
            PRODUCT_PAGE,
            "Poliform official Hima product page",
            "Manufacturer and family identity evidence only; not drawing geometry",
        ),
        (
            "POLIFORM-HIMA-OFFICIAL-TECHNICAL-SHEET",
            TECHNICAL_SHEET,
            "Poliform official Hima technical sheet",
            "Nominal dimensions and family identity evidence only; not the source of the geometry-derived linework",
        ),
        (
            "POLIFORM-HIMA-OFFICIAL-NEWS-2022-TECHNICAL-DATA",
            NEWS_TECHNICAL_PUBLICATION,
            "Poliform official News 2022 Hima technical data",
            "Pages 126-127 establish configuration and nominal dimensions only; not the source of linework",
        ),
        (
            "POLIFORM-HIMA-OFFICIAL-PAGE-EVIDENCE",
            relative(PAGE_EVIDENCE),
            "Poliform Hima official-page access evidence",
            f"SHA-256 {sha256(PAGE_EVIDENCE)}; records the official DWG listing and registration/CAPTCHA boundary",
        ),
        (
            "POLIFORM-HIMA-SOURCE-ACCESS-RECORD",
            relative(access_path),
            "Hima official CAD access record",
            f"SHA-256 {sha256(access_path)}; official 2D DWG published but not acquired; no third-party CAD used",
        ),
    )
    identifiers = []
    for identification, location, name, description in documents:
        reference = model.create_entity(
            "IfcDocumentReference",
            Location=location,
            Identification=identification,
            Name=name,
            Description=description,
            ReferencedDocument=None,
        )
        model.create_entity(
            "IfcRelAssociatesDocument",
            GlobalId=ifcopenshell.guid.new(),
            OwnerHistory=product.OwnerHistory,
            Name=f"{name} association",
            Description=SCOPE,
            RelatedObjects=[product, product_type],
            RelatingDocument=reference,
        )
        identifiers.append(identification)
    return sorted(identifiers)


def add_source_pset(model, product, product_type, approval, manifest_hash, candidate_path, access_path):
    values = {
        "SourceKind": SOURCE_KIND,
        "SourceLabelZh": SOURCE_LABEL_ZH,
        "Manufacturer": "Poliform",
        "Family": "Hima",
        "IFCTypeName": IFC_TYPE_NAME,
        "SourceProductPage": PRODUCT_PAGE,
        "SourceTechnicalSheet": TECHNICAL_SHEET,
        "SourceNews2022TechnicalPublication": NEWS_TECHNICAL_PUBLICATION,
        "Official2DDwgStatus": "published_registration_form_and_captcha_required_not_acquired",
        "OfficialCadUsed": str(OFFICIAL_CAD_USED).lower(),
        "ThirdPartyCadUsed": "false",
        "EvidenceScope": SCOPE,
        "ReviewStatus": "APPROVED",
        "Reviewer": approval["reviewer"],
        "ReviewDate": approval["review_date"],
        "ApprovalEvidence": approval["approval_evidence"],
        "ApprovedCandidateManifestSha256": manifest_hash,
        "ApprovedRepresentationIdentifiers": "Hima01Plan;Hima01Front;Hima01Side",
        "RepresentationGeometrySource": "proxy_paths_mm derived from the isolated representative IFC Body",
        "OfficialCadGeometryIncluded": str(OFFICIAL_CAD_GEOMETRY_INCLUDED).lower(),
        "CandidateRepresentations": relative(candidate_path),
        "CandidateRepresentationsSha256": sha256(candidate_path),
        "SourceAccessRecord": relative(access_path),
        "SourceAccessRecordSha256": sha256(access_path),
        "OfficialProductPageEvidence": relative(PAGE_EVIDENCE),
        "OfficialProductPageEvidenceSha256": sha256(PAGE_EVIDENCE),
        "PlanPathCount": str(EXPECTED_PATH_COUNTS["plan"]),
        "FrontPathCount": str(EXPECTED_PATH_COUNTS["front"]),
        "SidePathCount": str(EXPECTED_PATH_COUNTS["side"]),
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
        Name=PSET_NAME,
        Description="Mechanically verifiable geometry-derived drawing source and human approval",
        HasProperties=properties,
    )
    model.create_entity(
        "IfcRelDefinesByProperties",
        GlobalId=ifcopenshell.guid.new(),
        OwnerHistory=product.OwnerHistory,
        Name="HIMA01 drawing source properties",
        Description=None,
        RelatedObjects=[product, product_type],
        RelatingPropertyDefinition=pset,
    )


def representation_path_count(representation) -> int:
    return sum(len(curve_set.Elements) for curve_set in representation.Items)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--candidate", type=Path, default=DEFAULT_CANDIDATE)
    parser.add_argument("--source-access-record", type=Path, default=DEFAULT_ACCESS_RECORD)
    parser.add_argument("--approval", type=Path, default=DEFAULT_APPROVAL)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--report", type=Path)
    parser.add_argument("--apply", action="store_true")
    args = parser.parse_args()
    formal = args.input.resolve()
    output = args.output.resolve()
    manifest_path = args.manifest.resolve()
    candidate_path = args.candidate.resolve()
    access_path = args.source_access_record.resolve()
    if not args.apply:
        raise RuntimeError("IFC write requires the explicit --apply flag")
    if output == formal:
        raise RuntimeError("formal IFC cannot be the output; write a separate derived IFC")
    if output.exists():
        raise RuntimeError(f"refusing to overwrite existing output: {output}")
    formal_hash = sha256(formal)
    if formal_hash != FORMAL_SHA256:
        raise RuntimeError("formal IFC hash mismatch")
    manifest = load_json(manifest_path)
    candidate = load_json(candidate_path)
    access = load_json(access_path)
    approval = load_json(args.approval.resolve())
    require_approval(approval, manifest_path)
    if (
        manifest.get("formal_ifc_sha256") != formal_hash
        or manifest.get("formal_ifc_bytes_unchanged") is not True
        or manifest.get("review_status") != "visual_review_pending"
        or manifest.get("approved_for_drawing_ifc") is not False
        or manifest.get("representative_global_id") != REPRESENTATIVE_GLOBAL_ID
        or manifest.get("drawing_source", {}).get("source_kind") != SOURCE_KIND
        or manifest.get("project_context", {}).get("pass") is not True
        or manifest.get("project_context", {}).get("manifest_sha256") != sha256(PRODUCT_DIR / "project-context-manifest.json")
        or manifest.get("bonsai_review", {}).get("mode") != "actual_bonsai_ifc_body_camera_render"
        or manifest.get("bonsai_review", {}).get("manifest_sha256") != sha256(PRODUCT_DIR / "bonsai-review-manifest.json")
    ):
        raise RuntimeError("review manifest is not the expected pending Hima candidate")
    paths = candidate_paths(candidate, access)
    model = ifcopenshell.open(formal)
    product = model.by_guid(REPRESENTATIVE_GLOBAL_ID)
    if product is None:
        raise RuntimeError("Hima representative product is missing")
    product_type = next((relation.RelatingType for relation in product.IsTypedBy), None)
    actual_identity = product_type.Name if product_type is not None else product.Name
    if actual_identity != IFC_TYPE_NAME:
        raise RuntimeError("Hima representative type identity drifted")
    existing_identifiers = {item.RepresentationIdentifier for item in product.Representation.Representations}
    if any(identifier in existing_identifiers for identifier, _ in REPRESENTATIONS.values()):
        raise RuntimeError("Hima derived drawing representations already exist")
    representations = list(product.Representation.Representations)
    for view, (identifier, target_view) in REPRESENTATIONS.items():
        context = representation_context(model, identifier, target_view)
        representations.append(curve_representation(model, context, identifier, view, paths[view]))
    product.Representation.Representations = representations
    document_ids = add_document_associations(model, product, product_type, access_path)
    add_source_pset(model, product, product_type, approval, sha256(manifest_path), candidate_path, access_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    model.write(output)
    if sha256(formal) != formal_hash:
        output.unlink(missing_ok=True)
        raise RuntimeError("formal IFC bytes changed during derived write")
    derived = ifcopenshell.open(output)
    derived_product = derived.by_guid(REPRESENTATIVE_GLOBAL_ID)
    derived_type = next((relation.RelatingType for relation in derived_product.IsTypedBy), None)
    representation_counts = {}
    for view, (identifier, _) in REPRESENTATIONS.items():
        representation = next(
            (item for item in derived_product.Representation.Representations if item.RepresentationIdentifier == identifier),
            None,
        )
        if representation is None:
            output.unlink(missing_ok=True)
            raise RuntimeError(f"derived IFC lost {identifier}")
        representation_counts[view] = representation_path_count(representation)
    if representation_counts != EXPECTED_PATH_COUNTS:
        output.unlink(missing_ok=True)
        raise RuntimeError("derived IFC geometry-derived path counts drifted")
    document_relations = [
        relation
        for relation in derived.get_inverse(derived_product)
        if relation.is_a("IfcRelAssociatesDocument")
        and relation.RelatingDocument.is_a("IfcDocumentReference")
        and (relation.RelatingDocument.Identification or "").startswith(DOCUMENT_ID_PREFIX)
    ]
    actual_document_ids = sorted(relation.RelatingDocument.Identification for relation in document_relations)
    if actual_document_ids != document_ids or (
        derived_type is not None
        and any(derived_type not in relation.RelatedObjects for relation in document_relations)
    ):
        output.unlink(missing_ok=True)
        raise RuntimeError("derived IFC Hima document associations failed post-write verification")
    source_pset = ifcopenshell.util.element.get_pset(derived_product, PSET_NAME) or {}
    expected_properties = {
        "SourceKind": SOURCE_KIND,
        "SourceLabelZh": SOURCE_LABEL_ZH,
        "OfficialCadUsed": str(OFFICIAL_CAD_USED).lower(),
        "ThirdPartyCadUsed": "false",
        "EvidenceScope": SCOPE,
        "ReviewStatus": "APPROVED",
        "Reviewer": approval["reviewer"],
        "ReviewDate": approval["review_date"],
        "ApprovedCandidateManifestSha256": sha256(manifest_path),
        "OfficialCadGeometryIncluded": str(OFFICIAL_CAD_GEOMETRY_INCLUDED).lower(),
    }
    if any(source_pset.get(name) != value for name, value in expected_properties.items()):
        output.unlink(missing_ok=True)
        raise RuntimeError("derived IFC Hima source property post-write gate failed")
    result = {
        "schema_version": 1,
        "profile_key": PROFILE_KEY,
        "ifc_type_name": IFC_TYPE_NAME,
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
        "representations": {view: identifier for view, (identifier, _) in REPRESENTATIONS.items()},
        "representation_path_counts": representation_counts,
        "representation_geometry_source": source_pset["RepresentationGeometrySource"],
        "official_cad_geometry_included": OFFICIAL_CAD_GEOMETRY_INCLUDED,
        "source_document_associations": document_ids,
        "source_property_set": PSET_NAME,
        "source_kind": source_pset["SourceKind"],
        "source_label_zh": source_pset["SourceLabelZh"],
        "pass": True,
    }
    if args.report:
        write_json(args.report.resolve(), result)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
