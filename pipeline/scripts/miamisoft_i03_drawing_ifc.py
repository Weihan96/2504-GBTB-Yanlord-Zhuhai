#!/usr/bin/env python3
"""Approval-gated derived IFC writer for Baxter Miami Soft I03."""

import json
from pathlib import Path

import ifcopenshell
import ifcopenshell.guid

from falper_sorgente_linework import ROOT, relative, sha256
import hima01_drawing_ifc as shared


SCOPE = "exact manufacturer family, I03 roll model and nominal dimensions; project IFC Description incorrectly states the I01 108 cm length while type code and Body geometry match I03; authenticated native 2D/3D/BIM not acquired; not a project shop drawing and not official CAD geometry"
PRODUCT_PAGE = "https://www.baxter.it/gb/prodotti/miami-soft-divani-e-poltrone"
TECHNICAL_SHEET = "https://dam.baxter.it/m/1e3eaec668311b07/original/Baxter_MiamiSoft_divano_Indoor.pdf"
MEASUREMENT_SVG = "https://productsbook.baxter.it/models/measurements/MIAMSORUSI03.svg"
PRODUCT_DIR = ROOT / "output/review/highpoly-types/miamisoft-i03"
PRODUCT_PAGE_ARCHIVE = PRODUCT_DIR / "official-source/baxter-miami-soft-product-page.html"
TECHNICAL_SHEET_ARCHIVE = PRODUCT_DIR / "official-source/Baxter_MiamiSoft_divano_Indoor.pdf"
MEASUREMENT_SVG_ARCHIVE = PRODUCT_DIR / "official-source/MIAMSORUSI03.svg"
EXPECTED_PATH_COUNTS = {"plan": 1, "front": 1, "side": 1}

shared.REPRESENTATIVE_GLOBAL_ID = "1hdNtfQRPCGfHLfqzdp6w0"
shared.IFC_TYPE_NAME = "MiamiSoft I03"
shared.PROFILE_KEY = "miamisoft-i03"
shared.SCOPE = SCOPE
shared.REPRESENTATIONS = {
    "plan": ("MiamiSoftI03Plan", "PLAN_VIEW"),
    "front": ("MiamiSoftI03Front", "ELEVATION_VIEW"),
    "side": ("MiamiSoftI03Side", "ELEVATION_VIEW"),
}
shared.EXPECTED_PATH_COUNTS = EXPECTED_PATH_COUNTS
shared.PRODUCT_DIR = PRODUCT_DIR
shared.DEFAULT_MANIFEST = PRODUCT_DIR / "manifest.json"
shared.DEFAULT_CANDIDATE = PRODUCT_DIR / "candidate-representations.json"
shared.DEFAULT_ACCESS_RECORD = PRODUCT_DIR / "official-source/source-access-record.json"
shared.DEFAULT_APPROVAL = ROOT / "pipeline/decisions/miamisoft-i03-drawing-approval.json"
shared.PSET_NAME = "Pset_MiamiSoftI03DrawingSource"
shared.DOCUMENT_ID_PREFIX = "BAXTER-MIAMI-SOFT-I03-"


def candidate_paths(candidate: dict, access: dict) -> dict:
    if (
        candidate.get("profile_key") != shared.PROFILE_KEY
        or candidate.get("representative_global_id") != shared.REPRESENTATIVE_GLOBAL_ID
        or candidate.get("article_number") != "Baxter Miami Soft I03 roll cushion"
        or candidate.get("source_kind") != shared.SOURCE_KIND
        or candidate.get("source_label_zh") != shared.SOURCE_LABEL_ZH
        or candidate.get("official_cad_used") is not False
        or candidate.get("third_party_cad_used") is not False
        or candidate.get("formal_ifc_write_allowed") is not False
        or candidate.get("review_status") != "visual_review_pending"
    ):
        raise RuntimeError("pending Miami Soft I03 candidate source gate failed")
    official_cad = access.get("official_product_cad", {})
    drawing_source = access.get("drawing_geometry_source", {})
    identity = access.get("project_identity_evidence", {})
    if (
        access.get("resolved_variant") != "I03 - roll cushion - 170 x Ø27 cm"
        or access.get("scope") != SCOPE
        or access.get("pass") is not True
        or official_cad.get("authentication_required") is not True
        or official_cad.get("acquired") is not False
        or official_cad.get("local_cad_files") != []
        or official_cad.get("official_vector_pdf_archived") is not True
        or official_cad.get("official_measurement_svg_archived") is not True
        or official_cad.get("official_vector_references_used_as_cad_geometry") is not False
        or drawing_source.get("source_kind") != shared.SOURCE_KIND
        or drawing_source.get("official_cad_used") is not False
        or drawing_source.get("third_party_cad_used") is not False
        or identity.get("exact_model_code_status") != "confirmed"
        or identity.get("module_status") != "confirmed_I03_roll_170cm"
        or identity.get("description_status") != "conflicts_with_official_I03_and_matches_adjacent_I01_length"
        or identity.get("body_geometry_status") != "matches_official_I03_within_15mm"
        or access.get("dimension_cross_check", {}).get("pass") is not True
    ):
        raise RuntimeError("Miami Soft I03 official-source access record gate failed")
    paths = {}
    for view in shared.REQUIRED_VIEWS:
        item = candidate.get("views", {}).get(view, {})
        if item.get("source_kind") != shared.SOURCE_KIND or item.get("official_cad_paths_mm") != []:
            raise RuntimeError(f"Miami Soft I03 {view} must contain no official CAD paths")
        paths[view] = item.get("proxy_paths_mm", [])
        if len(paths[view]) != EXPECTED_PATH_COUNTS[view]:
            raise RuntimeError(f"Miami Soft I03 {view} geometry-derived path count drifted")
    return paths


def add_document_associations(model, product, product_type, access_path: Path):
    documents = (
        ("BAXTER-MIAMI-SOFT-I03-OFFICIAL-PRODUCT-PAGE", PRODUCT_PAGE, "Baxter official Miami Soft product page", "Exact family, designer, I03 roll model and authenticated native-file access boundary; not drawing geometry"),
        ("BAXTER-MIAMI-SOFT-I03-OFFICIAL-PRODUCT-PAGE-ARCHIVE", relative(PRODUCT_PAGE_ARCHIVE), "Archived Baxter Miami Soft product page", f"SHA-256 {sha256(PRODUCT_PAGE_ARCHIVE)}; identity and access evidence only"),
        ("BAXTER-MIAMI-SOFT-I03-OFFICIAL-TECHNICAL-SHEET", TECHNICAL_SHEET, "Baxter official Miami Soft vector technical sheet", "Page 8 confirms I03 roll and 170 x diameter 27 cm; evidence only, not representation linework"),
        ("BAXTER-MIAMI-SOFT-I03-OFFICIAL-TECHNICAL-SHEET-ARCHIVE", relative(TECHNICAL_SHEET_ARCHIVE), "Archived Baxter Miami Soft technical sheet", f"SHA-256 {sha256(TECHNICAL_SHEET_ARCHIVE)}; vector PDF not used as CAD geometry"),
        ("BAXTER-MIAMI-SOFT-I03-OFFICIAL-MEASUREMENT-SVG", MEASUREMENT_SVG, "Baxter official Miami Soft I03 measurement SVG", "Exact I03 three-view dimension graphic; evidence only, not native DWG/DXF representation geometry"),
        ("BAXTER-MIAMI-SOFT-I03-OFFICIAL-MEASUREMENT-SVG-ARCHIVE", relative(MEASUREMENT_SVG_ARCHIVE), "Archived Baxter Miami Soft I03 measurement SVG", f"SHA-256 {sha256(MEASUREMENT_SVG_ARCHIVE)}; official vector reference not used as CAD geometry"),
        ("BAXTER-MIAMI-SOFT-I03-SOURCE-ACCESS-RECORD", relative(access_path), "Baxter Miami Soft I03 source access record", f"SHA-256 {sha256(access_path)}; native 2D/3D/BIM login required; no third-party CAD used"),
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
        "Family": "Miami Soft",
        "Designer": "Paola Navone",
        "ModelCode": "I03",
        "Module": "I03 roll cushion",
        "IFCTypeDescription": "Roll 108 x Ø27 cm",
        "ResolvedOfficialVariant": "I03 - roll cushion - 170 x Ø27 cm",
        "IfcDescriptionConflict": "Project IFC Description uses the adjacent official I01 108 cm length; type code and Body geometry identify official I03 170 cm",
        "SourceProductPage": PRODUCT_PAGE,
        "SourceTechnicalSheet": TECHNICAL_SHEET,
        "SourceMeasurementSvg": MEASUREMENT_SVG,
        "Official2D3DBimStatus": "manufacturer_login_required_not_acquired",
        "OfficialVectorPdfArchived": "true",
        "OfficialMeasurementSvgArchived": "true",
        "OfficialVectorReferencesUsedAsCadGeometry": "false",
        "OfficialCadUsed": "false",
        "ThirdPartyCadUsed": "false",
        "EvidenceScope": SCOPE,
        "ReviewStatus": "APPROVED",
        "Reviewer": approval["reviewer"],
        "ReviewDate": approval["review_date"],
        "ApprovalEvidence": approval["approval_evidence"],
        "ApprovedCandidateManifestSha256": manifest_hash,
        "ApprovedRepresentationIdentifiers": "MiamiSoftI03Plan;MiamiSoftI03Front;MiamiSoftI03Side",
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
        "OfficialNominalLengthDiameterMm": json.dumps(dimensions["official_nominal_length_diameter_mm"]),
        "ProjectIfcBodyLocalXYZMm": json.dumps(dimensions["project_ifc_body_local_xyz_mm"]),
        "BodyMinusOfficialLengthYDiameterZDiameterMm": json.dumps(dimensions["body_minus_official_length_y_diameter_z_diameter_mm"]),
        "DimensionCrossCheckPass": str(dimensions["pass"]).lower(),
    }
    properties = [
        model.create_entity("IfcPropertySingleValue", Name=name, Description=None, NominalValue=model.create_entity("IfcText", str(value)), Unit=None)
        for name, value in values.items()
    ]
    pset = model.create_entity(
        "IfcPropertySet",
        GlobalId=ifcopenshell.guid.new(),
        OwnerHistory=product.OwnerHistory,
        Name=shared.PSET_NAME,
        Description="Mechanically verifiable geometry-derived source, exact I03 identity evidence and human approval",
        HasProperties=properties,
    )
    model.create_entity(
        "IfcRelDefinesByProperties",
        GlobalId=ifcopenshell.guid.new(),
        OwnerHistory=product.OwnerHistory,
        Name="Miami Soft I03 drawing source and model properties",
        Description=None,
        RelatedObjects=[product, product_type],
        RelatingPropertyDefinition=pset,
    )


shared.candidate_paths = candidate_paths
shared.add_document_associations = add_document_associations
shared.add_source_pset = add_source_pset


if __name__ == "__main__":
    shared.main()
