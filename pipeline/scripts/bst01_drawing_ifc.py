#!/usr/bin/env python3
"""Approval-gated derived IFC writer for Baxter Ninfea / project BST01."""

import json
from pathlib import Path

import ifcopenshell
import ifcopenshell.guid

from falper_sorgente_linework import ROOT, relative, sha256
import hima01_drawing_ifc as shared


SCOPE = "exact Baxter Ninfea bedside-table family and nominal 42 x 42 x 45 cm dimensions; official right- and left-opening public vectors both archived because the project type does not record opening side; native 2D/3D/BIM requires login and was not acquired; public PDF/SVG are not representation CAD geometry or a project shop drawing"
PRODUCT_PAGE = "https://www.baxter.it/en/products/ninfea-tables-and-coffee-tables"
TECHNICAL_SHEET = "https://productsbook.baxter.it/product-pdf/Ninfea_tavoli-e-tavolini_TechnicalSheet.pdf?code=NINF&kind=indoor&lang=eng&sector=tavoli-e-tavolini"
RIGHT_OPENING_SVG = "https://productsbook.baxter.it/models/measurements/NINFCLE40D.svg"
LEFT_OPENING_SVG = "https://productsbook.baxter.it/models/measurements/NINFCLE40S.svg"
PRODUCT_DIR = ROOT / "output/review/highpoly-types/bst01"
PRODUCT_PAGE_ARCHIVE = PRODUCT_DIR / "official-source/baxter-ninfea-product-page.html"
TECHNICAL_SHEET_ARCHIVE = PRODUCT_DIR / "official-source/Baxter_Ninfea_TechnicalSheet.pdf"
RIGHT_OPENING_SVG_ARCHIVE = PRODUCT_DIR / "official-source/NINFCLE40D-right-opening.svg"
LEFT_OPENING_SVG_ARCHIVE = PRODUCT_DIR / "official-source/NINFCLE40S-left-opening.svg"
EXPECTED_PATH_COUNTS = {"plan": 21, "front": 1, "side": 5}

shared.REPRESENTATIVE_GLOBAL_ID = "3eic1dzkn5heTIn4PhF37v"
shared.IFC_TYPE_NAME = "BST01"
shared.PROFILE_KEY = "bst01"
shared.SCOPE = SCOPE
shared.REPRESENTATIONS = {
    "plan": ("Bst01Plan", "PLAN_VIEW"),
    "front": ("Bst01Front", "ELEVATION_VIEW"),
    "side": ("Bst01Side", "ELEVATION_VIEW"),
}
shared.EXPECTED_PATH_COUNTS = EXPECTED_PATH_COUNTS
shared.PRODUCT_DIR = PRODUCT_DIR
shared.DEFAULT_MANIFEST = PRODUCT_DIR / "manifest.json"
shared.DEFAULT_CANDIDATE = PRODUCT_DIR / "candidate-representations.json"
shared.DEFAULT_ACCESS_RECORD = PRODUCT_DIR / "official-source/source-access-record.json"
shared.DEFAULT_APPROVAL = ROOT / "pipeline/decisions/bst01-drawing-approval.json"
shared.PSET_NAME = "Pset_Bst01DrawingSource"
shared.DOCUMENT_ID_PREFIX = "BAXTER-NINFEA-BST01-"


def candidate_paths(candidate: dict, access: dict) -> dict:
    if (
        candidate.get("profile_key") != shared.PROFILE_KEY
        or candidate.get("representative_global_id") != shared.REPRESENTATIVE_GLOBAL_ID
        or candidate.get("article_number") != "Baxter Ninfea bedside table opening side unresolved / BST01"
        or candidate.get("source_kind") != shared.SOURCE_KIND
        or candidate.get("source_label_zh") != shared.SOURCE_LABEL_ZH
        or candidate.get("official_cad_used") is not False
        or candidate.get("third_party_cad_used") is not False
        or candidate.get("formal_ifc_write_allowed") is not False
        or candidate.get("review_status") != "visual_review_pending"
    ):
        raise RuntimeError("pending BST01 candidate source gate failed")
    official_cad = access.get("official_product_cad", {})
    drawing_source = access.get("drawing_geometry_source", {})
    dimensions = access.get("dimension_cross_check", {})
    if (
        access.get("resolved_variant") != "Ninfea bedside-table family, 42 x 42 x 45 cm; opening side unresolved"
        or access.get("scope") != SCOPE
        or access.get("pass") is not True
        or official_cad.get("authentication_required") is not True
        or official_cad.get("acquired") is not False
        or official_cad.get("exact_project_configuration_match") is not False
        or official_cad.get("project_opening_side_resolved") is not False
        or official_cad.get("local_cad_files") != []
        or official_cad.get("official_vector_pdf_archived") is not True
        or official_cad.get("official_left_and_right_measurement_svgs_archived") is not True
        or official_cad.get("official_vector_evidence_used_as_cad_geometry") is not False
        or drawing_source.get("source_kind") != shared.SOURCE_KIND
        or drawing_source.get("official_cad_used") is not False
        or drawing_source.get("official_public_vector_evidence_used_as_cad_geometry") is not False
        or drawing_source.get("third_party_cad_used") is not False
        or dimensions.get("pass_by_axis") != [True, True, True]
        or dimensions.get("pass") is not True
        or dimensions.get("absolute_delta_mm") != [20.969086, 20.968094, 0.0]
    ):
        raise RuntimeError("BST01 official-source, unresolved-opening-side or dimension gate failed")
    paths = {}
    for view in shared.REQUIRED_VIEWS:
        item = candidate.get("views", {}).get(view, {})
        if item.get("source_kind") != shared.SOURCE_KIND or item.get("official_cad_paths_mm") != []:
            raise RuntimeError(f"BST01 {view} must contain no official CAD paths")
        paths[view] = item.get("proxy_paths_mm", [])
        if len(paths[view]) != EXPECTED_PATH_COUNTS[view]:
            raise RuntimeError(f"BST01 {view} geometry-derived path count drifted")
    return paths


def add_document_associations(model, product, product_type, access_path: Path):
    documents = (
        ("BAXTER-NINFEA-BST01-OFFICIAL-PRODUCT-PAGE", PRODUCT_PAGE, "Baxter official Ninfea product page", "Exact family, designer and nominal dimensions; not drawing geometry"),
        ("BAXTER-NINFEA-BST01-OFFICIAL-PRODUCT-PAGE-ARCHIVE", relative(PRODUCT_PAGE_ARCHIVE), "Archived Baxter Ninfea product page", f"SHA-256 {sha256(PRODUCT_PAGE_ARCHIVE)}; identity evidence only"),
        ("BAXTER-NINFEA-BST01-OFFICIAL-TECHNICAL-SHEET", TECHNICAL_SHEET, "Baxter official Ninfea technical sheet", "Page 7 publishes both opening variants; evidence only, not CAD representation linework"),
        ("BAXTER-NINFEA-BST01-OFFICIAL-TECHNICAL-SHEET-ARCHIVE", relative(TECHNICAL_SHEET_ARCHIVE), "Archived Baxter Ninfea technical sheet", f"SHA-256 {sha256(TECHNICAL_SHEET_ARCHIVE)}; vector PDF not used as CAD geometry"),
        ("BAXTER-NINFEA-BST01-OFFICIAL-RIGHT-OPENING-SVG", RIGHT_OPENING_SVG, "Baxter official Ninfea right-opening public SVG", "Right-opening family evidence only; project opening side is unresolved"),
        ("BAXTER-NINFEA-BST01-OFFICIAL-RIGHT-OPENING-SVG-ARCHIVE", relative(RIGHT_OPENING_SVG_ARCHIVE), "Archived Baxter Ninfea right-opening SVG", f"SHA-256 {sha256(RIGHT_OPENING_SVG_ARCHIVE)}; not used as representation geometry"),
        ("BAXTER-NINFEA-BST01-OFFICIAL-LEFT-OPENING-SVG", LEFT_OPENING_SVG, "Baxter official Ninfea left-opening public SVG", "Left-opening family evidence only; project opening side is unresolved"),
        ("BAXTER-NINFEA-BST01-OFFICIAL-LEFT-OPENING-SVG-ARCHIVE", relative(LEFT_OPENING_SVG_ARCHIVE), "Archived Baxter Ninfea left-opening SVG", f"SHA-256 {sha256(LEFT_OPENING_SVG_ARCHIVE)}; not used as representation geometry"),
        ("BAXTER-NINFEA-BST01-SOURCE-ACCESS-RECORD", relative(access_path), "Baxter Ninfea BST01 source access record", f"SHA-256 {sha256(access_path)}; native CAD unavailable and project opening side unresolved"),
    )
    identifiers = []
    for identification, location, name, description in documents:
        reference = model.create_entity("IfcDocumentReference", Location=location, Identification=identification, Name=name, Description=description, ReferencedDocument=None)
        model.create_entity("IfcRelAssociatesDocument", GlobalId=ifcopenshell.guid.new(), OwnerHistory=product.OwnerHistory, Name=f"{name} association", Description=SCOPE, RelatedObjects=[product, product_type], RelatingDocument=reference)
        identifiers.append(identification)
    return sorted(identifiers)


def add_source_pset(model, product, product_type, approval, manifest_hash, candidate_path, access_path):
    access = json.loads(access_path.read_text(encoding="utf-8"))
    dimensions = access["dimension_cross_check"]
    values = {
        "SourceKind": shared.SOURCE_KIND,
        "SourceLabelZh": shared.SOURCE_LABEL_ZH,
        "Manufacturer": "Baxter",
        "Family": "Ninfea",
        "Designer": "Pietro Russo",
        "ProjectTypeCode": "BST01",
        "Variant": "Ninfea bedside-table family, 42 x 42 x 45 cm; opening side unresolved",
        "IFCTypeDescription": "Nifea Comodino diam42xh45",
        "SourceProductPage": PRODUCT_PAGE,
        "SourceTechnicalSheet": TECHNICAL_SHEET,
        "SourceRightOpeningSvg": RIGHT_OPENING_SVG,
        "SourceLeftOpeningSvg": LEFT_OPENING_SVG,
        "Official2D3DBimStatus": "manufacturer_login_required_not_acquired",
        "OfficialCadUsed": "false",
        "OfficialOpeningSideResolved": "false",
        "OfficialRightOpeningVectorEvidence": "archived_not_representation_geometry",
        "OfficialLeftOpeningVectorEvidence": "archived_not_representation_geometry",
        "OfficialPublicVectorRole": "identity, nominal-dimension and opening-variant evidence only; not representation linework",
        "ThirdPartyCadUsed": "false",
        "EvidenceScope": SCOPE,
        "ReviewStatus": "APPROVED",
        "Reviewer": approval["reviewer"],
        "ReviewDate": approval["review_date"],
        "ApprovalEvidence": approval["approval_evidence"],
        "ApprovedCandidateManifestSha256": manifest_hash,
        "ApprovedRepresentationIdentifiers": "Bst01Plan;Bst01Front;Bst01Side",
        "RepresentationGeometrySource": "proxy_paths_mm derived from the isolated representative IFC MODEL_VIEW Body",
        "OfficialCadGeometryIncluded": "false",
        "CandidateRepresentations": relative(candidate_path),
        "CandidateRepresentationsSha256": sha256(candidate_path),
        "SourceAccessRecord": relative(access_path),
        "SourceAccessRecordSha256": sha256(access_path),
        "OfficialProductPageArchiveSha256": sha256(PRODUCT_PAGE_ARCHIVE),
        "OfficialTechnicalSheetArchiveSha256": sha256(TECHNICAL_SHEET_ARCHIVE),
        "OfficialRightOpeningSvgArchiveSha256": sha256(RIGHT_OPENING_SVG_ARCHIVE),
        "OfficialLeftOpeningSvgArchiveSha256": sha256(LEFT_OPENING_SVG_ARCHIVE),
        "PlanPathCount": str(EXPECTED_PATH_COUNTS["plan"]),
        "FrontPathCount": str(EXPECTED_PATH_COUNTS["front"]),
        "SidePathCount": str(EXPECTED_PATH_COUNTS["side"]),
        "OfficialVariantOverallMm": json.dumps(dimensions["official_variant_overall_mm"]),
        "ProjectIfcBodyLocalXYZMm": json.dumps(dimensions["project_ifc_body_local_xyz_mm"]),
        "BodyMinusOfficialAbsoluteMm": json.dumps(dimensions["absolute_delta_mm"]),
        "DimensionCrossCheckPass": "true",
    }
    properties = [model.create_entity("IfcPropertySingleValue", Name=name, Description=None, NominalValue=model.create_entity("IfcText", str(value)), Unit=None) for name, value in values.items()]
    pset = model.create_entity("IfcPropertySet", GlobalId=ifcopenshell.guid.new(), OwnerHistory=product.OwnerHistory, Name=shared.PSET_NAME, Description="Mechanically verifiable geometry-derived source, unresolved opening side and human approval", HasProperties=properties)
    model.create_entity("IfcRelDefinesByProperties", GlobalId=ifcopenshell.guid.new(), OwnerHistory=product.OwnerHistory, Name="BST01 drawing source and unresolved opening-side properties", Description=None, RelatedObjects=[product, product_type], RelatingPropertyDefinition=pset)


shared.candidate_paths = candidate_paths
shared.add_document_associations = add_document_associations
shared.add_source_pset = add_source_pset


if __name__ == "__main__":
    shared.main()
