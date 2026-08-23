#!/usr/bin/env python3
"""Approval-gated derived IFC writer for Baxter Stone / project BST03."""

import json
from pathlib import Path

import ifcopenshell
import ifcopenshell.guid

from falper_sorgente_linework import ROOT, relative, sha256
import hima01_drawing_ifc as shared


SCOPE = "exact manufacturer family, left-drawer variant and nominal dimensions; not a project shop drawing and not official CAD geometry"
PRODUCT_PAGE = "https://www.baxter.it/en/products/stone-beds"
TECHNICAL_SHEET = "https://dam.baxter.it/m/e4f4056bc7c6088/original/Baxter_Stone_letto_Indoor.pdf"
PRODUCT_DIR = ROOT / "output/review/highpoly-types/bst03"
PRODUCT_PAGE_ARCHIVE = PRODUCT_DIR / "official-source/baxter-stone-product-page.html"
TECHNICAL_SHEET_ARCHIVE = PRODUCT_DIR / "official-source/Baxter_Stone_letto_Indoor.pdf"
EXPECTED_PATH_COUNTS = {"plan": 1, "front": 4, "side": 2}

shared.REPRESENTATIVE_GLOBAL_ID = "1HZoxe$df4cBb4UXXH5J2S"
shared.IFC_TYPE_NAME = "BST03"
shared.PROFILE_KEY = "bst03"
shared.SCOPE = SCOPE
shared.REPRESENTATIONS = {
    "plan": ("Bst03Plan", "PLAN_VIEW"),
    "front": ("Bst03Front", "ELEVATION_VIEW"),
    "side": ("Bst03Side", "ELEVATION_VIEW"),
}
shared.EXPECTED_PATH_COUNTS = EXPECTED_PATH_COUNTS
shared.PRODUCT_DIR = PRODUCT_DIR
shared.DEFAULT_MANIFEST = PRODUCT_DIR / "manifest.json"
shared.DEFAULT_CANDIDATE = PRODUCT_DIR / "candidate-representations.json"
shared.DEFAULT_ACCESS_RECORD = PRODUCT_DIR / "official-source/source-access-record.json"
shared.DEFAULT_APPROVAL = ROOT / "pipeline/decisions/bst03-drawing-approval.json"
shared.PSET_NAME = "Pset_Bst03DrawingSource"
shared.DOCUMENT_ID_PREFIX = "BAXTER-STONE-BST03-"


def candidate_paths(candidate: dict, access: dict) -> dict:
    if (
        candidate.get("profile_key") != shared.PROFILE_KEY
        or candidate.get("representative_global_id") != shared.REPRESENTATIVE_GLOBAL_ID
        or candidate.get("article_number") != "Baxter Stone L drawer 45 / BST03"
        or candidate.get("source_kind") != shared.SOURCE_KIND
        or candidate.get("source_label_zh") != shared.SOURCE_LABEL_ZH
        or candidate.get("official_cad_used") is not False
        or candidate.get("third_party_cad_used") is not False
        or candidate.get("formal_ifc_write_allowed") is not False
        or candidate.get("review_status") != "visual_review_pending"
    ):
        raise RuntimeError("pending BST03 candidate source gate failed")
    official_cad = access.get("official_product_cad", {})
    drawing_source = access.get("drawing_geometry_source", {})
    project_identity = access.get("project_identity_evidence", {})
    if (
        access.get("resolved_variant") != "freestanding bedside table with L drawer, 45 x 45 x 46 cm"
        or access.get("scope") != SCOPE
        or access.get("pass") is not True
        or official_cad.get("authentication_required") is not True
        or official_cad.get("acquired") is not False
        or official_cad.get("local_cad_files") != []
        or official_cad.get("official_vector_pdf_archived") is not True
        or official_cad.get("official_vector_pdf_used_as_cad_geometry") is not False
        or drawing_source.get("source_kind") != shared.SOURCE_KIND
        or drawing_source.get("official_cad_used") is not False
        or drawing_source.get("third_party_cad_used") is not False
        or project_identity.get("confirmed_variant") != "freestanding_bedside_table_with_L_drawer"
        or project_identity.get("register_status") != "confirmed"
        or access.get("dimension_cross_check", {}).get("pass") is not True
    ):
        raise RuntimeError("BST03 official-source access record gate failed")
    paths = {}
    for view in shared.REQUIRED_VIEWS:
        item = candidate.get("views", {}).get(view, {})
        if item.get("source_kind") != shared.SOURCE_KIND or item.get("official_cad_paths_mm") != []:
            raise RuntimeError(f"BST03 {view} must contain no official CAD paths")
        paths[view] = item.get("proxy_paths_mm", [])
        if len(paths[view]) != EXPECTED_PATH_COUNTS[view]:
            raise RuntimeError(f"BST03 {view} geometry-derived path count drifted")
    return paths


def add_document_associations(model, product, product_type, access_path: Path):
    documents = (
        (
            "BAXTER-STONE-BST03-OFFICIAL-PRODUCT-PAGE",
            PRODUCT_PAGE,
            "Baxter official Stone product page",
            "Exact family, designer, left-drawer variant and nominal dimensions; not drawing geometry",
        ),
        (
            "BAXTER-STONE-BST03-OFFICIAL-PRODUCT-PAGE-ARCHIVE",
            relative(PRODUCT_PAGE_ARCHIVE),
            "Archived Baxter Stone product page",
            f"SHA-256 {sha256(PRODUCT_PAGE_ARCHIVE)}; identity evidence only",
        ),
        (
            "BAXTER-STONE-BST03-OFFICIAL-TECHNICAL-SHEET",
            TECHNICAL_SHEET,
            "Baxter official Stone technical sheet",
            "Page 10 publishes vector views and dimensions for left/right 45 cm variants; evidence only, not CAD representation linework",
        ),
        (
            "BAXTER-STONE-BST03-OFFICIAL-TECHNICAL-SHEET-ARCHIVE",
            relative(TECHNICAL_SHEET_ARCHIVE),
            "Archived Baxter Stone technical sheet",
            f"SHA-256 {sha256(TECHNICAL_SHEET_ARCHIVE)}; vector PDF not used as CAD geometry",
        ),
        (
            "BAXTER-STONE-BST03-SOURCE-ACCESS-RECORD",
            relative(access_path),
            "Baxter Stone BST03 source access record",
            f"SHA-256 {sha256(access_path)}; native 2D/3D/BIM login required; no third-party CAD used",
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
        "Family": "Stone",
        "Designer": "Federico Peri",
        "ProjectTypeCode": "BST03",
        "Variant": "freestanding bedside table with L drawer, 45 x 45 x 46 cm",
        "IFCTypeDescription": "Stone Bedside Table with Drawer D45",
        "SourceProductPage": PRODUCT_PAGE,
        "SourceTechnicalSheet": TECHNICAL_SHEET,
        "Official2D3DBimStatus": "manufacturer_login_required_not_acquired",
        "OfficialVectorPdfRole": "identity, variant and nominal-dimension evidence only; not representation linework",
        "OfficialCadUsed": "false",
        "ThirdPartyCadUsed": "false",
        "EvidenceScope": SCOPE,
        "ReviewStatus": "APPROVED",
        "Reviewer": approval["reviewer"],
        "ReviewDate": approval["review_date"],
        "ApprovalEvidence": approval["approval_evidence"],
        "ApprovedCandidateManifestSha256": manifest_hash,
        "ApprovedRepresentationIdentifiers": "Bst03Plan;Bst03Front;Bst03Side",
        "RepresentationGeometrySource": "proxy_paths_mm derived from the isolated representative IFC MODEL_VIEW Body",
        "OfficialCadGeometryIncluded": "false",
        "CandidateRepresentations": relative(candidate_path),
        "CandidateRepresentationsSha256": sha256(candidate_path),
        "SourceAccessRecord": relative(access_path),
        "SourceAccessRecordSha256": sha256(access_path),
        "OfficialProductPageArchiveSha256": sha256(PRODUCT_PAGE_ARCHIVE),
        "OfficialTechnicalSheetArchiveSha256": sha256(TECHNICAL_SHEET_ARCHIVE),
        "PlanPathCount": str(EXPECTED_PATH_COUNTS["plan"]),
        "FrontPathCount": str(EXPECTED_PATH_COUNTS["front"]),
        "SidePathCount": str(EXPECTED_PATH_COUNTS["side"]),
        "OfficialVariantOverallMm": json.dumps(dimensions["official_variant_overall_mm"]),
        "ProjectIfcBodyLocalXYZMm": json.dumps(dimensions["project_ifc_body_local_xyz_mm"]),
        "BodyMinusOfficialAbsoluteMm": json.dumps(dimensions["absolute_delta_mm"]),
        "DimensionCrossCheckPass": str(dimensions["pass"]).lower(),
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
        Description="Mechanically verifiable geometry-derived source, exact variant evidence and human approval",
        HasProperties=properties,
    )
    model.create_entity(
        "IfcRelDefinesByProperties",
        GlobalId=ifcopenshell.guid.new(),
        OwnerHistory=product.OwnerHistory,
        Name="BST03 drawing source and variant properties",
        Description=None,
        RelatedObjects=[product, product_type],
        RelatingPropertyDefinition=pset,
    )


shared.candidate_paths = candidate_paths
shared.add_document_associations = add_document_associations
shared.add_source_pset = add_source_pset


if __name__ == "__main__":
    shared.main()
