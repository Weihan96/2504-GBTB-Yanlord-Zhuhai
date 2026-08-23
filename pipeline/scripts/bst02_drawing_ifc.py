#!/usr/bin/env python3
"""Approval-gated derived IFC writer for Baxter Beside / project BST02."""

import json
from pathlib import Path

import ifcopenshell
import ifcopenshell.guid

from falper_sorgente_linework import ROOT, relative, sha256
import hima01_drawing_ifc as shared


SCOPE = "exact Baxter Beside 55 x 58 x 35 cm manufacturer identity and public measurement evidence; native 2D/3D/BIM requires login and was not acquired; public PDF/SVG are not representation CAD geometry or a project shop drawing; the 100 mm project Body height discrepancy is retained"
PRODUCT_PAGE = "https://www.baxter.it/en/products/beside-tables-and-coffee-tables"
TECHNICAL_SHEET = "https://productsbook.baxter.it/product-pdf/Beside_tavoli-e-tavolini_TechnicalSheet.pdf?code=BESI&kind=indoor&lang=eng&sector=tavoli-e-tavolini"
MEASUREMENT_SVG = "https://productsbook.baxter.it/models/measurements/BESICOBS55.svg"
PRODUCT_DIR = ROOT / "output/review/highpoly-types/bst02"
PRODUCT_PAGE_ARCHIVE = PRODUCT_DIR / "official-source/baxter-beside-product-page.html"
TECHNICAL_SHEET_ARCHIVE = PRODUCT_DIR / "official-source/Baxter_Beside_TechnicalSheet.pdf"
MEASUREMENT_SVG_ARCHIVE = PRODUCT_DIR / "official-source/BESICOBS55.svg"
EXPECTED_PATH_COUNTS = {"plan": 6, "front": 1, "side": 1}

shared.REPRESENTATIVE_GLOBAL_ID = "3hA0vKpcn44u4Tsx4tqiUz"
shared.IFC_TYPE_NAME = "BST02"
shared.PROFILE_KEY = "bst02"
shared.SCOPE = SCOPE
shared.REPRESENTATIONS = {
    "plan": ("Bst02Plan", "PLAN_VIEW"),
    "front": ("Bst02Front", "ELEVATION_VIEW"),
    "side": ("Bst02Side", "ELEVATION_VIEW"),
}
shared.EXPECTED_PATH_COUNTS = EXPECTED_PATH_COUNTS
shared.PRODUCT_DIR = PRODUCT_DIR
shared.DEFAULT_MANIFEST = PRODUCT_DIR / "manifest.json"
shared.DEFAULT_CANDIDATE = PRODUCT_DIR / "candidate-representations.json"
shared.DEFAULT_ACCESS_RECORD = PRODUCT_DIR / "official-source/source-access-record.json"
shared.DEFAULT_APPROVAL = ROOT / "pipeline/decisions/bst02-drawing-approval.json"
shared.PSET_NAME = "Pset_Bst02DrawingSource"
shared.DOCUMENT_ID_PREFIX = "BAXTER-BESIDE-BST02-"


def candidate_paths(candidate: dict, access: dict) -> dict:
    if (
        candidate.get("profile_key") != shared.PROFILE_KEY
        or candidate.get("representative_global_id") != shared.REPRESENTATIVE_GLOBAL_ID
        or candidate.get("article_number") != "Baxter Beside 55x58xh35 / BST02"
        or candidate.get("source_kind") != shared.SOURCE_KIND
        or candidate.get("source_label_zh") != shared.SOURCE_LABEL_ZH
        or candidate.get("official_cad_used") is not False
        or candidate.get("third_party_cad_used") is not False
        or candidate.get("formal_ifc_write_allowed") is not False
        or candidate.get("review_status") != "visual_review_pending"
    ):
        raise RuntimeError("pending BST02 candidate source gate failed")
    official_cad = access.get("official_product_cad", {})
    drawing_source = access.get("drawing_geometry_source", {})
    dimensions = access.get("dimension_cross_check", {})
    if (
        access.get("resolved_variant") != "Beside bedside table, 55 x 58 x 35 cm"
        or access.get("scope") != SCOPE
        or access.get("pass") is not True
        or official_cad.get("authentication_required") is not True
        or official_cad.get("acquired") is not False
        or official_cad.get("local_cad_files") != []
        or official_cad.get("official_vector_pdf_archived") is not True
        or official_cad.get("official_measurement_svg_archived") is not True
        or official_cad.get("official_vector_evidence_used_as_cad_geometry") is not False
        or drawing_source.get("source_kind") != shared.SOURCE_KIND
        or drawing_source.get("official_cad_used") is not False
        or drawing_source.get("official_public_vector_evidence_used_as_cad_geometry") is not False
        or drawing_source.get("third_party_cad_used") is not False
        or dimensions.get("pass_by_axis") != [True, True, False]
        or dimensions.get("pass") is not False
        or dimensions.get("review_required") is not True
        or dimensions.get("absolute_delta_mm") != [10.0, 2.110748, 100.0]
    ):
        raise RuntimeError("BST02 official-source and height-discrepancy gate failed")
    paths = {}
    for view in shared.REQUIRED_VIEWS:
        item = candidate.get("views", {}).get(view, {})
        if item.get("source_kind") != shared.SOURCE_KIND or item.get("official_cad_paths_mm") != []:
            raise RuntimeError(f"BST02 {view} must contain no official CAD paths")
        paths[view] = item.get("proxy_paths_mm", [])
        if len(paths[view]) != EXPECTED_PATH_COUNTS[view]:
            raise RuntimeError(f"BST02 {view} geometry-derived path count drifted")
    return paths


def add_document_associations(model, product, product_type, access_path: Path):
    documents = (
        ("BAXTER-BESIDE-BST02-OFFICIAL-PRODUCT-PAGE", PRODUCT_PAGE, "Baxter official Beside product page", "Exact family, designer and nominal dimensions; not drawing geometry"),
        ("BAXTER-BESIDE-BST02-OFFICIAL-PRODUCT-PAGE-ARCHIVE", relative(PRODUCT_PAGE_ARCHIVE), "Archived Baxter Beside product page", f"SHA-256 {sha256(PRODUCT_PAGE_ARCHIVE)}; identity evidence only"),
        ("BAXTER-BESIDE-BST02-OFFICIAL-TECHNICAL-SHEET", TECHNICAL_SHEET, "Baxter official Beside technical sheet", "Page 7 publishes vector views and dimensions; evidence only, not CAD representation linework"),
        ("BAXTER-BESIDE-BST02-OFFICIAL-TECHNICAL-SHEET-ARCHIVE", relative(TECHNICAL_SHEET_ARCHIVE), "Archived Baxter Beside technical sheet", f"SHA-256 {sha256(TECHNICAL_SHEET_ARCHIVE)}; vector PDF not used as CAD geometry"),
        ("BAXTER-BESIDE-BST02-OFFICIAL-MEASUREMENT-SVG", MEASUREMENT_SVG, "Baxter official Beside public measurement SVG", "Public dimensioned side, front and plan evidence; not native 2D CAD representation geometry"),
        ("BAXTER-BESIDE-BST02-OFFICIAL-MEASUREMENT-SVG-ARCHIVE", relative(MEASUREMENT_SVG_ARCHIVE), "Archived Baxter Beside measurement SVG", f"SHA-256 {sha256(MEASUREMENT_SVG_ARCHIVE)}; official vector evidence not used as CAD geometry"),
        ("BAXTER-BESIDE-BST02-SOURCE-ACCESS-RECORD", relative(access_path), "Baxter Beside BST02 source access record", f"SHA-256 {sha256(access_path)}; native 2D/3D/BIM login required; 100 mm Body height discrepancy retained"),
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
        "Family": "Beside",
        "Designer": "Studiopepe",
        "ProjectTypeCode": "BST02",
        "Variant": "Beside bedside table, 55 x 58 x 35 cm",
        "IFCTypeDescription": "Beside 55x58xh35",
        "SourceProductPage": PRODUCT_PAGE,
        "SourceTechnicalSheet": TECHNICAL_SHEET,
        "SourceMeasurementSvg": MEASUREMENT_SVG,
        "Official2D3DBimStatus": "manufacturer_login_required_not_acquired",
        "OfficialPublicVectorRole": "identity, nominal-dimension and view evidence only; not representation linework",
        "OfficialCadUsed": "false",
        "ThirdPartyCadUsed": "false",
        "EvidenceScope": SCOPE,
        "ReviewStatus": "APPROVED",
        "Reviewer": approval["reviewer"],
        "ReviewDate": approval["review_date"],
        "ApprovalEvidence": approval["approval_evidence"],
        "ApprovedCandidateManifestSha256": manifest_hash,
        "ApprovedRepresentationIdentifiers": "Bst02Plan;Bst02Front;Bst02Side",
        "RepresentationGeometrySource": "proxy_paths_mm derived from the isolated representative IFC MODEL_VIEW Body",
        "OfficialCadGeometryIncluded": "false",
        "CandidateRepresentations": relative(candidate_path),
        "CandidateRepresentationsSha256": sha256(candidate_path),
        "SourceAccessRecord": relative(access_path),
        "SourceAccessRecordSha256": sha256(access_path),
        "OfficialProductPageArchiveSha256": sha256(PRODUCT_PAGE_ARCHIVE),
        "OfficialTechnicalSheetArchiveSha256": sha256(TECHNICAL_SHEET_ARCHIVE),
        "OfficialMeasurementSvgArchiveSha256": sha256(MEASUREMENT_SVG_ARCHIVE),
        "PlanPathCount": str(EXPECTED_PATH_COUNTS["plan"]),
        "FrontPathCount": str(EXPECTED_PATH_COUNTS["front"]),
        "SidePathCount": str(EXPECTED_PATH_COUNTS["side"]),
        "OfficialVariantOverallMm": json.dumps(dimensions["official_variant_overall_mm"]),
        "ProjectIfcBodyLocalXYZMm": json.dumps(dimensions["project_ifc_body_local_xyz_mm"]),
        "BodyMinusOfficialAbsoluteMm": json.dumps(dimensions["absolute_delta_mm"]),
        "DimensionCrossCheckPass": "false",
        "HeightDiscrepancyReviewRequired": "true",
    }
    properties = [model.create_entity("IfcPropertySingleValue", Name=name, Description=None, NominalValue=model.create_entity("IfcText", str(value)), Unit=None) for name, value in values.items()]
    pset = model.create_entity("IfcPropertySet", GlobalId=ifcopenshell.guid.new(), OwnerHistory=product.OwnerHistory, Name=shared.PSET_NAME, Description="Mechanically verifiable geometry-derived source, official public evidence, retained height discrepancy and human approval", HasProperties=properties)
    model.create_entity("IfcRelDefinesByProperties", GlobalId=ifcopenshell.guid.new(), OwnerHistory=product.OwnerHistory, Name="BST02 drawing source and discrepancy properties", Description=None, RelatedObjects=[product, product_type], RelatingPropertyDefinition=pset)


shared.candidate_paths = candidate_paths
shared.add_document_associations = add_document_associations
shared.add_source_pset = add_source_pset


if __name__ == "__main__":
    shared.main()
