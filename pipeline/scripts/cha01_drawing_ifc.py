#!/usr/bin/env python3
"""Approval-gated derived IFC writer for Baxter Colette / CHA01."""

import json

import ifcopenshell
import ifcopenshell.guid

from falper_sorgente_linework import ROOT, relative, sha256
import hima01_drawing_ifc as shared


SCOPE = "exact Baxter Colette armchair 57 x 60 x 73 cm manufacturer identity; authenticated 2D/3D/BIM not acquired; technical PDFs and measurement SVG are identity evidence only; not a project shop drawing or official CAD geometry"
PRODUCT_PAGE = "https://www.baxter.it/en/products/colette-chairs"
CURRENT_TECHNICAL_SHEET = "https://productsbook.baxter.it/product-pdf/Colette_sedie_TechnicalSheet.pdf?code=COLE&kind=indoor&lang=eng&sector=sedie"
STABLE_TECHNICAL_SHEET = "https://dam.baxter.it/m/3b8ed18c6647c201/original/Baxter_Colette_sedia_Indoor.pdf"
MEASUREMENT_SVG = "https://productsbook.baxter.it/models/measurements/COLEPOCO57.svg"
LOGIN_PAGE = "https://www.baxter.it/en/user/login"
PRODUCT_DIR = ROOT / "output/review/highpoly-types/cha01"
SOURCE_DIR = PRODUCT_DIR / "official-source"
EXPECTED_PATH_COUNTS = {"plan": 1, "front": 6, "side": 4}

shared.REPRESENTATIVE_GLOBAL_ID = "1luHljRzDAhPNXxTDNu7qB"
shared.IFC_TYPE_NAME = "CHA01"
shared.PROFILE_KEY = "cha01"
shared.SCOPE = SCOPE
shared.PRODUCT_PAGE = PRODUCT_PAGE
shared.TECHNICAL_SHEET = CURRENT_TECHNICAL_SHEET
shared.REPRESENTATIONS = {
    "plan": ("Cha01Plan", "PLAN_VIEW"),
    "front": ("Cha01Front", "ELEVATION_VIEW"),
    "side": ("Cha01Side", "ELEVATION_VIEW"),
}
shared.EXPECTED_PATH_COUNTS = EXPECTED_PATH_COUNTS
shared.PRODUCT_DIR = PRODUCT_DIR
shared.DEFAULT_MANIFEST = PRODUCT_DIR / "manifest.json"
shared.DEFAULT_CANDIDATE = PRODUCT_DIR / "candidate-representations.json"
shared.DEFAULT_ACCESS_RECORD = SOURCE_DIR / "source-access-record.json"
shared.DEFAULT_APPROVAL = ROOT / "pipeline/decisions/cha01-drawing-approval.json"
shared.PSET_NAME = "Pset_Cha01DrawingSource"
shared.DOCUMENT_ID_PREFIX = "BAXTER-COLETTE-CHA01-"


def candidate_paths(candidate: dict, access: dict) -> dict:
    if (
        candidate.get("profile_key") != shared.PROFILE_KEY
        or candidate.get("representative_global_id") != shared.REPRESENTATIVE_GLOBAL_ID
        or candidate.get("article_number") != "Baxter Colette armchair 57 x 60 x 73 cm / CHA01"
        or candidate.get("source_kind") != shared.SOURCE_KIND
        or candidate.get("source_label_zh") != shared.SOURCE_LABEL_ZH
        or candidate.get("official_cad_used") is not False
        or candidate.get("third_party_cad_used") is not False
        or candidate.get("formal_ifc_write_allowed") is not False
        or candidate.get("review_status") != "visual_review_pending"
    ):
        raise RuntimeError("pending CHA01 candidate source gate failed")
    drawing_source = access.get("drawing_geometry_source", {})
    official_cad = access.get("official_product_cad", {})
    dimension_check = access.get("dimension_cross_check", {})
    if (
        official_cad.get("authentication_required") is not True
        or official_cad.get("acquired") is not False
        or official_cad.get("local_cad_files") != []
        or official_cad.get("third_party_cad_used") is not False
        or drawing_source.get("source_kind") != shared.SOURCE_KIND
        or drawing_source.get("official_cad_used") is not False
        or drawing_source.get("third_party_cad_used") is not False
        or dimension_check.get("status") != "exact_model_geometry_consistent_no_scaling_or_fitting"
        or access.get("scope") != SCOPE
        or access.get("pass") is not True
    ):
        raise RuntimeError("CHA01 official-source access record gate failed")
    paths = {}
    for view in shared.REQUIRED_VIEWS:
        item = candidate.get("views", {}).get(view, {})
        if item.get("source_kind") != shared.SOURCE_KIND or item.get("official_cad_paths_mm") != []:
            raise RuntimeError(f"CHA01 {view} linework must remain geometry-derived with no official CAD paths")
        paths[view] = item.get("proxy_paths_mm", [])
        if len(paths[view]) != EXPECTED_PATH_COUNTS[view]:
            raise RuntimeError(f"CHA01 {view} geometry-derived path-count gate failed")
    return paths


def add_document_associations(model, product, product_type, access_path):
    product_archive = SOURCE_DIR / "baxter-colette-product-page.html"
    current_sheet = SOURCE_DIR / "Baxter_Colette_current-technical-sheet.pdf"
    stable_sheet = SOURCE_DIR / "Baxter_Colette_sedia_Indoor.pdf"
    measurement = SOURCE_DIR / "COLEPOCO57.svg"
    login_archive = SOURCE_DIR / "baxter-login-page.html"
    documents = (
        ("BAXTER-COLETTE-CHA01-OFFICIAL-PRODUCT-PAGE", PRODUCT_PAGE, "Baxter official Colette product page", "Exact Colette armchair identity and dimensions; not drawing geometry"),
        ("BAXTER-COLETTE-CHA01-PRODUCT-PAGE-ARCHIVE", relative(product_archive), "Archived Baxter Colette product page", f"SHA-256 {sha256(product_archive)}; identity and authenticated-download evidence only"),
        ("BAXTER-COLETTE-CHA01-CURRENT-TECHNICAL-SHEET", CURRENT_TECHNICAL_SHEET, "Baxter current Colette technical sheet", "Page 37 identifies the 57 x 60 x 73 cm armchair; miniature vectors are not representation geometry"),
        ("BAXTER-COLETTE-CHA01-CURRENT-TECHNICAL-SHEET-ARCHIVE", relative(current_sheet), "Archived Baxter current Colette technical sheet", f"SHA-256 {sha256(current_sheet)}; identity and dimension evidence only"),
        ("BAXTER-COLETTE-CHA01-STABLE-TECHNICAL-SHEET", STABLE_TECHNICAL_SHEET, "Baxter stable Colette technical sheet", "Page 6 identifies the 57 x 60 x 73 cm little armchair; not drawing geometry"),
        ("BAXTER-COLETTE-CHA01-STABLE-TECHNICAL-SHEET-ARCHIVE", relative(stable_sheet), "Archived Baxter stable Colette technical sheet", f"SHA-256 {sha256(stable_sheet)}; identity and dimension evidence only"),
        ("BAXTER-COLETTE-CHA01-MEASUREMENT-SVG", MEASUREMENT_SVG, "Baxter COLEPOCO57 measurement SVG", "Exact manufacturer variant identifier and nominal dimensions; not CAD representation geometry"),
        ("BAXTER-COLETTE-CHA01-MEASUREMENT-SVG-ARCHIVE", relative(measurement), "Archived Baxter COLEPOCO57 measurement SVG", f"SHA-256 {sha256(measurement)}; identity evidence only"),
        ("BAXTER-COLETTE-CHA01-LOGIN-PAGE", LOGIN_PAGE, "Baxter authenticated download page", "Records the access boundary for manufacturer 2D, 3D and BIM downloads"),
        ("BAXTER-COLETTE-CHA01-LOGIN-PAGE-ARCHIVE", relative(login_archive), "Archived Baxter login page", f"SHA-256 {sha256(login_archive)}; no authenticated CAD was acquired"),
        ("BAXTER-COLETTE-CHA01-SOURCE-ACCESS-RECORD", relative(access_path), "Baxter Colette CAD access record", f"SHA-256 {sha256(access_path)}; official CAD not acquired; no third-party CAD used"),
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
        "Family": "Colette",
        "Model": "Armchair / Poltroncina 57 x 60 x 73 cm",
        "Designer": "Roberto Lazzeroni",
        "ProjectTypeCode": "CHA01",
        "IFCTypeDescription": "Baxter Colette",
        "SourceProductPage": PRODUCT_PAGE,
        "SourceCurrentTechnicalSheet": CURRENT_TECHNICAL_SHEET,
        "SourceStableTechnicalSheet": STABLE_TECHNICAL_SHEET,
        "SourceMeasurementSvg": MEASUREMENT_SVG,
        "Official2D3DBimStatus": "baxter_authentication_required_not_acquired",
        "OfficialCadUsed": "false",
        "ThirdPartyCadUsed": "false",
        "EvidenceScope": SCOPE,
        "ReviewStatus": "APPROVED",
        "Reviewer": approval["reviewer"],
        "ReviewDate": approval["review_date"],
        "ApprovalEvidence": approval["approval_evidence"],
        "ApprovedCandidateManifestSha256": manifest_hash,
        "ApprovedRepresentationIdentifiers": "Cha01Plan;Cha01Front;Cha01Side",
        "RepresentationGeometrySource": "proxy_paths_mm derived from the isolated representative IFC Body",
        "OfficialCadGeometryIncluded": "false",
        "CandidateRepresentations": relative(candidate_path),
        "CandidateRepresentationsSha256": sha256(candidate_path),
        "SourceAccessRecord": relative(access_path),
        "SourceAccessRecordSha256": sha256(access_path),
        "PlanPathCount": str(EXPECTED_PATH_COUNTS["plan"]),
        "FrontPathCount": str(EXPECTED_PATH_COUNTS["front"]),
        "SidePathCount": str(EXPECTED_PATH_COUNTS["side"]),
        "OfficialModelDimensionsMm": json.dumps(dimensions["official_model_width_depth_height_mm"]),
        "ProjectIfcBodyDimensionsMm": json.dumps(dimensions["project_ifc_body_local_xyz_mm"]),
        "BodyMinusOfficialMm": json.dumps(dimensions["body_minus_official_mm"]),
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
        Description="Mechanically verifiable geometry-derived drawing source and human approval",
        HasProperties=properties,
    )
    model.create_entity(
        "IfcRelDefinesByProperties",
        GlobalId=ifcopenshell.guid.new(),
        OwnerHistory=product.OwnerHistory,
        Name="CHA01 drawing source properties",
        Description=None,
        RelatedObjects=[product, product_type],
        RelatingPropertyDefinition=pset,
    )


shared.candidate_paths = candidate_paths
shared.add_document_associations = add_document_associations
shared.add_source_pset = add_source_pset


if __name__ == "__main__":
    shared.main()
