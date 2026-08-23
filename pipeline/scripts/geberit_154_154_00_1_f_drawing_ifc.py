#!/usr/bin/env python3
"""Approval-gated derived IFC writer for the project 154.154.00.1.F flange."""

import json
from pathlib import Path

import ifcopenshell
import ifcopenshell.guid

from falper_sorgente_linework import ROOT, relative, sha256
import hima01_drawing_ifc as shared


SCOPE = "project-model flange subcomponent of the exact Geberit 154.154.00.1 parent installation set; .F is a project authoring decomposition label, not a separately published manufacturer article; parent EPS is not component CAD geometry or a project shop drawing"
PRODUCT_PAGE = "https://catalog.geberit.co.uk/en-GB/product/PRO_199058"
PRODUCT_DIR = ROOT / "output/review/highpoly-types/geberit-154-154-00-1-f"
EXPECTED_PATH_COUNTS = {"plan": 2, "front": 7, "side": 6}
OFFICIAL_FILES = {
    "PARENT-PRODUCT-PAGE-ARCHIVE": "official-source/geberit-154-154-00-1-product-page.html",
    "PARENT-VECTOR-EPS-TOP": "official-source/DAS_199904-top-view.eps",
    "PARENT-VECTOR-EPS-FRONT": "official-source/DAS_199902-front-view.eps",
    "PARENT-VECTOR-EPS-PERSPECTIVE": "official-source/DAS_199900-perspective.eps",
}

shared.REPRESENTATIVE_GLOBAL_ID = "2jjNIn9gHBYwNSWwlM5T_i"
shared.IFC_TYPE_NAME = "Geberit 154.154.00.1.F"
shared.PROFILE_KEY = "geberit-154-154-00-1-f"
shared.SCOPE = SCOPE
shared.REPRESENTATIONS = {
    "plan": ("Geberit154154FlangePlan", "PLAN_VIEW"),
    "front": ("Geberit154154FlangeFront", "ELEVATION_VIEW"),
    "side": ("Geberit154154FlangeSide", "ELEVATION_VIEW"),
}
shared.EXPECTED_PATH_COUNTS = EXPECTED_PATH_COUNTS
shared.PRODUCT_DIR = PRODUCT_DIR
shared.DEFAULT_MANIFEST = PRODUCT_DIR / "manifest.json"
shared.DEFAULT_CANDIDATE = PRODUCT_DIR / "candidate-representations.json"
shared.DEFAULT_ACCESS_RECORD = PRODUCT_DIR / "official-source/source-access-record.json"
shared.DEFAULT_APPROVAL = ROOT / "pipeline/decisions/geberit-154-154-00-1-f-drawing-approval.json"
shared.PSET_NAME = "Pset_Geberit154154FlangeDrawingSource"
shared.DOCUMENT_ID_PREFIX = "GEBERIT-154-154-00-1-F-"


def candidate_paths(candidate: dict, access: dict) -> dict:
    if (
        candidate.get("profile_key") != shared.PROFILE_KEY
        or candidate.get("representative_global_id") != shared.REPRESENTATIVE_GLOBAL_ID
        or candidate.get("article_number") != "154.154.00.1 / project component .F"
        or candidate.get("source_kind") != shared.SOURCE_KIND
        or candidate.get("source_label_zh") != shared.SOURCE_LABEL_ZH
        or candidate.get("official_cad_used") is not False
        or candidate.get("third_party_cad_used") is not False
        or candidate.get("formal_ifc_write_allowed") is not False
        or candidate.get("review_status") != "visual_review_pending"
    ):
        raise RuntimeError("pending Geberit .F component candidate source gate failed")
    component_cad = access.get("official_component_cad", {})
    project_geometry = access.get("project_component_geometry", {})
    drawing_source = access.get("drawing_geometry_source", {})
    if (
        access.get("parent_manufacturer_article_number") != "154.154.00.1"
        or access.get("project_component_code") != "154.154.00.1.F"
        or access.get("project_component_role") != "Flange"
        or access.get("scope") != SCOPE
        or access.get("pass") is not True
        or component_cad.get("published_as_separate_manufacturer_article") is not False
        or component_cad.get("acquired") is not False
        or component_cad.get("exact_project_component_match") is not False
        or component_cad.get("local_cad_files") != []
        or component_cad.get("official_parent_vector_eps_archived") is not True
        or component_cad.get("official_parent_vector_eps_used_as_component_cad_geometry") is not False
        or set(component_cad.get("parent_native_dwg_http_statuses", {}).values()) != {404}
        or project_geometry.get("body_local_xyz_mm") != [88.0, 210.0, 104.999998]
        or project_geometry.get("official_standalone_component_dimensions_available") is not False
        or drawing_source.get("source_kind") != shared.SOURCE_KIND
        or drawing_source.get("official_cad_used") is not False
        or drawing_source.get("third_party_cad_used") is not False
    ):
        raise RuntimeError("Geberit .F component official-source boundary gate failed")
    paths = {}
    for view in shared.REQUIRED_VIEWS:
        item = candidate.get("views", {}).get(view, {})
        if item.get("source_kind") != shared.SOURCE_KIND or item.get("official_cad_paths_mm") != []:
            raise RuntimeError(f"Geberit .F {view} must contain no official component CAD paths")
        paths[view] = item.get("proxy_paths_mm", [])
        if len(paths[view]) != EXPECTED_PATH_COUNTS[view]:
            raise RuntimeError(f"Geberit .F {view} Body-derived path count drifted")
    return paths


def add_document_associations(model, product, product_type, access_path: Path):
    documents = [
        (
            "GEBERIT-154-154-00-1-F-OFFICIAL-PARENT-PRODUCT-PAGE",
            PRODUCT_PAGE,
            "Geberit official parent product page for 154.154.00.1",
            "Parent article and component relationship evidence only; .F is not a separately published article",
        ),
        (
            "GEBERIT-154-154-00-1-F-SOURCE-ACCESS-RECORD",
            relative(access_path),
            "Geberit 154.154.00.1.F source access record",
            f"SHA-256 {sha256(access_path)}; no standalone component CAD; no third-party CAD used",
        ),
    ]
    labels = {
        "PARENT-PRODUCT-PAGE-ARCHIVE": "Archived Geberit official parent product page",
        "PARENT-VECTOR-EPS-TOP": "Geberit official parent top-view EPS",
        "PARENT-VECTOR-EPS-FRONT": "Geberit official parent front-view EPS",
        "PARENT-VECTOR-EPS-PERSPECTIVE": "Geberit official parent perspective EPS",
    }
    for suffix, relative_path in OFFICIAL_FILES.items():
        path = PRODUCT_DIR / relative_path
        documents.append((
            f"GEBERIT-154-154-00-1-F-{suffix}",
            relative(path),
            labels[suffix],
            f"SHA-256 {sha256(path)}; complete-parent evidence only; never component representation linework",
        ))
    identifiers = []
    for identification, location, name, description in documents:
        reference = model.create_entity(
            "IfcDocumentReference", Location=location, Identification=identification,
            Name=name, Description=description, ReferencedDocument=None,
        )
        model.create_entity(
            "IfcRelAssociatesDocument", GlobalId=ifcopenshell.guid.new(),
            OwnerHistory=product.OwnerHistory, Name=f"{name} association", Description=SCOPE,
            RelatedObjects=[product, product_type], RelatingDocument=reference,
        )
        identifiers.append(identification)
    return sorted(identifiers)


def add_source_pset(model, product, product_type, approval, manifest_hash, candidate_path, access_path):
    source_hashes = {suffix: sha256(PRODUCT_DIR / path) for suffix, path in OFFICIAL_FILES.items()}
    values = {
        "SourceKind": shared.SOURCE_KIND,
        "SourceLabelZh": shared.SOURCE_LABEL_ZH,
        "Manufacturer": "Geberit",
        "Family": "CleanLine shower channel installation set",
        "ParentArticleNumber": "154.154.00.1",
        "ProjectComponentCode": "154.154.00.1.F",
        "ProjectComponentRole": "Flange",
        "ComponentCodeIsManufacturerArticle": "false",
        "ProjectIFCTypeName": shared.IFC_TYPE_NAME,
        "SourceParentProductPage": PRODUCT_PAGE,
        "OfficialComponentCadStatus": "standalone component not published; parent cadDrawings undefined; parent native A/G/L/P DWGs returned 404",
        "OfficialCadUsed": "false",
        "OfficialComponentCadUsed": "false",
        "ParentVectorEPSUsedAsComponentGeometry": "false",
        "ThirdPartyCadUsed": "false",
        "EvidenceScope": SCOPE,
        "ReviewStatus": "APPROVED",
        "Reviewer": approval["reviewer"],
        "ReviewDate": approval["review_date"],
        "ApprovalEvidence": approval["approval_evidence"],
        "ApprovedCandidateManifestSha256": manifest_hash,
        "ApprovedRepresentationIdentifiers": "Geberit154154FlangePlan;Geberit154154FlangeFront;Geberit154154FlangeSide",
        "RepresentationGeometrySource": "proxy_paths_mm derived from the isolated project .F IFC Body",
        "OfficialCadGeometryIncluded": "false",
        "CandidateRepresentations": relative(candidate_path),
        "CandidateRepresentationsSha256": sha256(candidate_path),
        "SourceAccessRecord": relative(access_path),
        "SourceAccessRecordSha256": sha256(access_path),
        "OfficialParentSourceFileSha256": json.dumps(source_hashes, sort_keys=True),
        "PlanPathCount": "2",
        "FrontPathCount": "7",
        "SidePathCount": "6",
        "ProjectIFCBodyLocalXYZMm": "88.0;210.0;104.999998",
        "OfficialStandaloneComponentDimensionsAvailable": "false",
        "NativeProjectElevationClippingRecorded": "true",
    }
    properties = [
        model.create_entity(
            "IfcPropertySingleValue", Name=name, Description=None,
            NominalValue=model.create_entity("IfcText", str(value)), Unit=None,
        )
        for name, value in values.items()
    ]
    pset = model.create_entity(
        "IfcPropertySet", GlobalId=ifcopenshell.guid.new(), OwnerHistory=product.OwnerHistory,
        Name=shared.PSET_NAME,
        Description="Mechanically verifiable project-component boundary, Body-derived source and human approval",
        HasProperties=properties,
    )
    model.create_entity(
        "IfcRelDefinesByProperties", GlobalId=ifcopenshell.guid.new(), OwnerHistory=product.OwnerHistory,
        Name="Geberit 154.154.00.1.F drawing source properties", Description=None,
        RelatedObjects=[product, product_type], RelatingPropertyDefinition=pset,
    )


shared.candidate_paths = candidate_paths
shared.add_document_associations = add_document_associations
shared.add_source_pset = add_source_pset


if __name__ == "__main__":
    shared.main()
