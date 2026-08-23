#!/usr/bin/env python3
"""Approval-gated derived IFC writer for the STREET-H subcomponent."""

import json

import ifcopenshell
import ifcopenshell.guid

from falper_sorgente_linework import ROOT, relative, sha256
import hima01_drawing_ifc as shared


SCOPE = "official complete antoniolupi Street family-top CAD evidence only; the project STREET-H type is an isolated repeated sink-holder subcomponent, not the complete catalogue top and not a project shop drawing"
PRODUCT_PAGE = "https://www.antoniolupi.it/en/products/sinks/street"
PRODUCT_DIR = ROOT / "output/review/highpoly-types/street-h"
SOURCE_DIR = PRODUCT_DIR / "official-source"
PAGE_EVIDENCE = SOURCE_DIR / "official-product-page-evidence.json"
CAD_ZIP = SOURCE_DIR / "ANTONIOLUPI-official-Street-2D-CAD.zip"
NATIVE_DXF = SOURCE_DIR / "AL_Street.dxf"
TECHNICAL_PDF = SOURCE_DIR / "ANTONIOLUPI-official-Street-technical.pdf"
CATALOGUE_EXTRACT = SOURCE_DIR / "ANTONIOLUPI-official-Street-LevantoRed-extract.pdf"
DXF_LINEWORK = PRODUCT_DIR / "official-native-dxf-linework.json"
EXPECTED_PATH_COUNTS = {"plan": 1, "front": 1, "side": 5}

shared.REPRESENTATIVE_GLOBAL_ID = "2ajpw0I9n1dBypfISg3ejX"
shared.IFC_TYPE_NAME = "STREET-H"
shared.PROFILE_KEY = "street-h"
shared.SCOPE = SCOPE
shared.REPRESENTATIONS = {
    "plan": ("StreetHPlan", "PLAN_VIEW"),
    "front": ("StreetHFront", "ELEVATION_VIEW"),
    "side": ("StreetHSide", "ELEVATION_VIEW"),
}
shared.EXPECTED_PATH_COUNTS = EXPECTED_PATH_COUNTS
shared.PRODUCT_DIR = PRODUCT_DIR
shared.DEFAULT_MANIFEST = PRODUCT_DIR / "manifest.json"
shared.DEFAULT_CANDIDATE = PRODUCT_DIR / "candidate-representations.json"
shared.DEFAULT_ACCESS_RECORD = SOURCE_DIR / "source-access-record.json"
shared.DEFAULT_APPROVAL = ROOT / "pipeline/decisions/street-h-drawing-approval.json"
shared.PSET_NAME = "Pset_StreetHDrawingSource"
shared.DOCUMENT_ID_PREFIX = "ANTONIOLUPI-STREET-H-"


def candidate_paths(candidate: dict, access: dict) -> dict:
    if (
        candidate.get("profile_key") != shared.PROFILE_KEY
        or candidate.get("representative_global_id") != shared.REPRESENTATIVE_GLOBAL_ID
        or candidate.get("source_kind") != shared.SOURCE_KIND
        or candidate.get("source_label_zh") != shared.SOURCE_LABEL_ZH
        or candidate.get("official_cad_used") is not False
        or candidate.get("third_party_cad_used") is not False
        or candidate.get("formal_ifc_write_allowed") is not False
        or candidate.get("review_status") != "visual_review_pending"
    ):
        raise RuntimeError("pending STREET-H candidate source gate failed")
    cad = access.get("official_product_cad", {})
    drawing_source = access.get("drawing_geometry_source", {})
    identity = access.get("project_identity_cross_check", {})
    cross_check = access.get("cad_pdf_cross_check", {})
    if (
        cad.get("acquired") is not True
        or cad.get("exact_parent_family_cluster_located") is not True
        or cad.get("exact_project_component_match") is not False
        or cad.get("official_parent_paths_used_as_project_representation") is not False
        or cad.get("native_dxf_sha256") != sha256(NATIVE_DXF)
        or cad.get("native_dxf_linework_sha256") != sha256(DXF_LINEWORK)
        or drawing_source.get("source_kind") != shared.SOURCE_KIND
        or drawing_source.get("official_cad_used") is not False
        or drawing_source.get("third_party_cad_used") is not False
        or identity.get("pass") is not True
        or cross_check.get("pass") is not True
        or access.get("scope") != SCOPE
        or sha256(CAD_ZIP) != "36161ac81bade86b7ec0419c0d9cf0c49d52ae0ecb9dced39006b9eddcb10e62"
        or sha256(TECHNICAL_PDF) != "0b883bbf6df617331cc2408da4860b31160af67a6b6924b3caaee0ec389d2946"
    ):
        raise RuntimeError("STREET-H official parent-family source gate failed")
    paths = {}
    for view in shared.REQUIRED_VIEWS:
        item = candidate.get("views", {}).get(view, {})
        if item.get("source_kind") != shared.SOURCE_KIND or item.get("official_cad_paths_mm") != []:
            raise RuntimeError(f"STREET-H {view} linework must remain geometry-derived with no official CAD paths")
        paths[view] = item.get("proxy_paths_mm", [])
        if len(paths[view]) != EXPECTED_PATH_COUNTS[view]:
            raise RuntimeError(f"STREET-H {view} geometry-derived path-count gate failed")
    return paths


def add_document_associations(model, product, product_type, access_path):
    documents = (
        (
            "ANTONIOLUPI-STREET-H-OFFICIAL-PRODUCT-PAGE",
            PRODUCT_PAGE,
            "antoniolupi official Street product page",
            "Manufacturer, Street family and download-surface evidence only; not STREET-H representation geometry",
        ),
        (
            "ANTONIOLUPI-STREET-H-OFFICIAL-PAGE-EVIDENCE",
            relative(PAGE_EVIDENCE),
            "antoniolupi Street official-page observation record",
            f"SHA-256 {sha256(PAGE_EVIDENCE)}; family identity and download labels",
        ),
        (
            "ANTONIOLUPI-STREET-H-OFFICIAL-CAD-ZIP",
            relative(CAD_ZIP),
            "antoniolupi official Street 2D CAD archive",
            f"SHA-256 {sha256(CAD_ZIP)}; complete parent-family top CAD, rejected as STREET-H representation geometry",
        ),
        (
            "ANTONIOLUPI-STREET-H-OFFICIAL-NATIVE-DXF",
            relative(NATIVE_DXF),
            "antoniolupi official AL_Street native DXF",
            f"SHA-256 {sha256(NATIVE_DXF)}; selected parent-family cluster STREET240 + STREET4054",
        ),
        (
            "ANTONIOLUPI-STREET-H-OFFICIAL-TECHNICAL-PDF",
            relative(TECHNICAL_PDF),
            "antoniolupi official Street technical data sheet",
            f"SHA-256 {sha256(TECHNICAL_PDF)}; page 1 proves complete-top dimensions and component mismatch",
        ),
        (
            "ANTONIOLUPI-STREET-H-OFFICIAL-CATALOGUE-EXTRACT",
            relative(CATALOGUE_EXTRACT),
            "antoniolupi official Street catalogue extract",
            f"SHA-256 {sha256(CATALOGUE_EXTRACT)}; Street and Rosso Levanto family evidence",
        ),
        (
            "ANTONIOLUPI-STREET-H-DXF-LINEWORK-AUDIT",
            relative(DXF_LINEWORK),
            "STREET-H parent-family native DXF mechanical audit",
            f"SHA-256 {sha256(DXF_LINEWORK)}; parent top paths are not used as the STREET-H representation",
        ),
        (
            "ANTONIOLUPI-STREET-H-SOURCE-ACCESS-RECORD",
            relative(access_path),
            "STREET-H manufacturer-source access record",
            f"SHA-256 {sha256(access_path)}; exact component CAD rejected and no third-party CAD used",
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
    access = json.loads(access_path.read_text(encoding="utf-8"))
    values = {
        "SourceKind": shared.SOURCE_KIND,
        "SourceLabelZh": shared.SOURCE_LABEL_ZH,
        "Manufacturer": "antoniolupi",
        "Family": "Street",
        "ProjectTypeCode": "STREET-H",
        "IFCTypeDescription": "Sink holder",
        "SourceProductPage": PRODUCT_PAGE,
        "OfficialCadAcquired": "true",
        "OfficialNativeDxfSha256": sha256(NATIVE_DXF),
        "OfficialCadUsed": "false",
        "ThirdPartyCadUsed": "false",
        "ExactProjectComponentCadMatch": "false",
        "RejectedParentFamilyCluster": "STREET240 prof. 40 + STREET4054 prof. 40",
        "RejectedParentFamilyBoundsMm": "1080 x 400 x 250",
        "ProjectStreetHBodyBoundsMm": "300 x 150 x 100",
        "EvidenceScope": SCOPE,
        "ReviewStatus": "APPROVED",
        "Reviewer": approval["reviewer"],
        "ReviewDate": approval["review_date"],
        "ApprovalEvidence": approval["approval_evidence"],
        "ApprovedCandidateManifestSha256": manifest_hash,
        "ApprovedRepresentationIdentifiers": "StreetHPlan;StreetHFront;StreetHSide",
        "RepresentationGeometrySource": "proxy_paths_mm derived from one isolated representative IFC Body",
        "OfficialCadGeometryIncluded": "false",
        "CandidateRepresentations": relative(candidate_path),
        "CandidateRepresentationsSha256": sha256(candidate_path),
        "SourceAccessRecord": relative(access_path),
        "SourceAccessRecordSha256": sha256(access_path),
        "OfficialDxfLineworkAudit": relative(DXF_LINEWORK),
        "OfficialDxfLineworkAuditSha256": sha256(DXF_LINEWORK),
        "CadPdfCrossCheckStatus": str(access["cad_pdf_cross_check"]["pass"]).lower(),
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
        Name=shared.PSET_NAME,
        Description="Mechanically verifiable geometry-derived STREET-H drawing source and human approval",
        HasProperties=properties,
    )
    model.create_entity(
        "IfcRelDefinesByProperties",
        GlobalId=ifcopenshell.guid.new(),
        OwnerHistory=product.OwnerHistory,
        Name="STREET-H drawing source properties",
        Description=None,
        RelatedObjects=[product, product_type],
        RelatingPropertyDefinition=pset,
    )


shared.candidate_paths = candidate_paths
shared.add_document_associations = add_document_associations
shared.add_source_pset = add_source_pset


if __name__ == "__main__":
    shared.main()
