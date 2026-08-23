#!/usr/bin/env python3
"""Approval-gated derived IFC writer for Baxter Casablanca / project BED01."""

import json
from pathlib import Path

import ifcopenshell
import ifcopenshell.guid

from falper_sorgente_linework import ROOT, relative, sha256
import hima01_drawing_ifc as shared


SCOPE = "manufacturer family and 180 x 200 cm variant identity only; manufacturer dimension conflict requires review; not a project shop drawing and not drawing geometry"
PRODUCT_PAGE = "https://www.baxter.it/en/products/casablanca-beds"
TECHNICAL_SHEET = "https://dam.baxter.it/m/394e8cfc516bb3a/original/Baxter_Casablanca_letto_Indoor.pdf"
PRODUCT_DIR = ROOT / "output/review/highpoly-types/bed01"
EXPECTED_PATH_COUNTS = {"plan": 3, "front": 7, "side": 4}

shared.REPRESENTATIVE_GLOBAL_ID = "3IQBEqO5vDI8Z9k1Ltge_N"
shared.IFC_TYPE_NAME = "BED01"
shared.PROFILE_KEY = "bed01"
shared.SCOPE = SCOPE
shared.REPRESENTATIONS = {
    "plan": ("Bed01Plan", "PLAN_VIEW"),
    "front": ("Bed01Front", "ELEVATION_VIEW"),
    "side": ("Bed01Side", "ELEVATION_VIEW"),
}
shared.EXPECTED_PATH_COUNTS = EXPECTED_PATH_COUNTS
shared.PRODUCT_DIR = PRODUCT_DIR
shared.DEFAULT_MANIFEST = PRODUCT_DIR / "manifest.json"
shared.DEFAULT_CANDIDATE = PRODUCT_DIR / "candidate-representations.json"
shared.DEFAULT_ACCESS_RECORD = PRODUCT_DIR / "official-source/source-access-record.json"
shared.DEFAULT_APPROVAL = ROOT / "pipeline/decisions/bed01-drawing-approval.json"
shared.PSET_NAME = "Pset_Bed01DrawingSource"
shared.DOCUMENT_ID_PREFIX = "BAXTER-CASABLANCA-BED01-"


def candidate_paths(candidate: dict, access: dict) -> dict:
    if (
        candidate.get("profile_key") != shared.PROFILE_KEY
        or candidate.get("representative_global_id") != shared.REPRESENTATIVE_GLOBAL_ID
        or candidate.get("article_number") != "Baxter Casablanca 180 / BED01"
        or candidate.get("source_kind") != shared.SOURCE_KIND
        or candidate.get("source_label_zh") != shared.SOURCE_LABEL_ZH
        or candidate.get("official_cad_used") is not False
        or candidate.get("third_party_cad_used") is not False
        or candidate.get("formal_ifc_write_allowed") is not False
        or candidate.get("review_status") != "visual_review_pending"
    ):
        raise RuntimeError("pending BED01 candidate source gate failed")
    drawing_source = access.get("drawing_geometry_source", {})
    official_cad = access.get("official_product_cad", {})
    dimension_check = access.get("dimension_cross_check", {})
    if (
        official_cad.get("authentication_required") is not True
        or official_cad.get("acquired") is not False
        or official_cad.get("local_cad_files") != []
        or drawing_source.get("source_kind") != shared.SOURCE_KIND
        or drawing_source.get("official_cad_used") is not False
        or drawing_source.get("third_party_cad_used") is not False
        or dimension_check.get("status") != "manufacturer_sources_disagree_review_required_not_scaled_or_corrected"
        or access.get("scope") != SCOPE
    ):
        raise RuntimeError("BED01 official-source access record gate failed")
    paths = {}
    for view in shared.REQUIRED_VIEWS:
        item = candidate.get("views", {}).get(view, {})
        if item.get("source_kind") != shared.SOURCE_KIND or item.get("official_cad_paths_mm") != []:
            raise RuntimeError(f"BED01 {view} linework must remain geometry-derived with no official CAD paths")
        paths[view] = item.get("proxy_paths_mm", [])
        if len(paths[view]) != EXPECTED_PATH_COUNTS[view]:
            raise RuntimeError(f"BED01 {view} geometry-derived path-count gate failed")
    return paths


def add_document_associations(model, product, product_type, access_path: Path):
    documents = (
        (
            "BAXTER-CASABLANCA-BED01-OFFICIAL-PRODUCT-PAGE",
            PRODUCT_PAGE,
            "Baxter official Casablanca product page",
            "Manufacturer and variant identity evidence; currently lists 2100 x 2520 x 900 mm; not drawing geometry",
        ),
        (
            "BAXTER-CASABLANCA-BED01-OFFICIAL-TECHNICAL-SHEET",
            TECHNICAL_SHEET,
            "Baxter official Casablanca technical sheet",
            "Lists the 180 x 200 cm version as 2200 x 2520 x 900 mm; dimension evidence only, not the source of linework",
        ),
        (
            "BAXTER-CASABLANCA-BED01-SOURCE-ACCESS-RECORD",
            relative(access_path),
            "Baxter Casablanca CAD access and dimension-conflict record",
            f"SHA-256 {sha256(access_path)}; manufacturer login required; exact CAD not acquired; source conflict retained",
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
    dimensions = access["dimension_cross_check"]
    values = {
        "SourceKind": shared.SOURCE_KIND,
        "SourceLabelZh": shared.SOURCE_LABEL_ZH,
        "Manufacturer": "Baxter",
        "Family": "Casablanca",
        "Designer": "Paola Navone",
        "ProjectTypeCode": "BED01",
        "IFCTypeDescription": "Casablanca 180",
        "SourceProductPage": PRODUCT_PAGE,
        "SourceTechnicalSheet": TECHNICAL_SHEET,
        "Official2D3DBimStatus": "manufacturer_login_required_not_acquired",
        "OfficialCadUsed": "false",
        "ThirdPartyCadUsed": "false",
        "EvidenceScope": SCOPE,
        "ReviewStatus": "APPROVED",
        "Reviewer": approval["reviewer"],
        "ReviewDate": approval["review_date"],
        "ApprovalEvidence": approval["approval_evidence"],
        "ApprovedCandidateManifestSha256": manifest_hash,
        "ApprovedRepresentationIdentifiers": "Bed01Plan;Bed01Front;Bed01Side",
        "RepresentationGeometrySource": "proxy_paths_mm derived from the isolated representative IFC Body",
        "OfficialCadGeometryIncluded": "false",
        "CandidateRepresentations": relative(candidate_path),
        "CandidateRepresentationsSha256": sha256(candidate_path),
        "SourceAccessRecord": relative(access_path),
        "SourceAccessRecordSha256": sha256(access_path),
        "PlanPathCount": str(EXPECTED_PATH_COUNTS["plan"]),
        "FrontPathCount": str(EXPECTED_PATH_COUNTS["front"]),
        "SidePathCount": str(EXPECTED_PATH_COUNTS["side"]),
        "OfficialProductPageOverallMm": json.dumps(dimensions["official_product_page_180x200_variant_overall_mm"]),
        "OfficialTechnicalSheetOverallMm": json.dumps(dimensions["official_technical_sheet_180x200_variant_overall_mm"]),
        "ProjectIfcBodyBoundsMm": json.dumps(dimensions["project_ifc_body_bounds_mm"]),
        "BodyMinusProductPageMm": json.dumps(dimensions["body_minus_product_page_mm"]),
        "BodyMinusTechnicalSheetMm": json.dumps(dimensions["body_minus_technical_sheet_mm"]),
        "DimensionReviewStatus": dimensions["status"],
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
        Description="Mechanically verifiable geometry-derived source, manufacturer dimension conflict and human approval",
        HasProperties=properties,
    )
    model.create_entity(
        "IfcRelDefinesByProperties",
        GlobalId=ifcopenshell.guid.new(),
        OwnerHistory=product.OwnerHistory,
        Name="BED01 drawing source and dimension-conflict properties",
        Description=None,
        RelatedObjects=[product, product_type],
        RelatingPropertyDefinition=pset,
    )


shared.candidate_paths = candidate_paths
shared.add_document_associations = add_document_associations
shared.add_source_pset = add_source_pset


if __name__ == "__main__":
    shared.main()
