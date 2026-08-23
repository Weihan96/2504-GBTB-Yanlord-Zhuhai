#!/usr/bin/env python3
"""Approval-gated derived IFC writer for Baxter Viktor / project BED02."""

import json
from pathlib import Path

import ifcopenshell
import ifcopenshell.guid

from falper_sorgente_linework import ROOT, relative, sha256
import hima01_drawing_ifc as shared


SCOPE = "manufacturer family identity and nominal dimensions only; dimension difference requires review; not a project shop drawing and not drawing geometry"
PRODUCT_PAGE = "https://www.baxter.it/en/products/viktor-beds"
TECHNICAL_SHEET = "https://dam.baxter.it/m/e9f40f5d27de0f5/original/Baxter_Viktor_letto_Indoor.pdf"
ASSEMBLY_INSTRUCTIONS = "https://dam.baxter.it/m/55b97eeff80c31a9/original/Baxter_AssemblyInstructions_Viktor_bed.pdf"
PRODUCT_DIR = ROOT / "output/review/highpoly-types/bed02"
EXPECTED_PATH_COUNTS = {"plan": 1, "front": 1, "side": 4}

shared.REPRESENTATIVE_GLOBAL_ID = "1i_pqgLv9A7uuV7MjaArBW"
shared.IFC_TYPE_NAME = "BED02"
shared.PROFILE_KEY = "bed02"
shared.SCOPE = SCOPE
shared.REPRESENTATIONS = {
    "plan": ("Bed02Plan", "PLAN_VIEW"),
    "front": ("Bed02Front", "ELEVATION_VIEW"),
    "side": ("Bed02Side", "ELEVATION_VIEW"),
}
shared.EXPECTED_PATH_COUNTS = EXPECTED_PATH_COUNTS
shared.PRODUCT_DIR = PRODUCT_DIR
shared.DEFAULT_MANIFEST = PRODUCT_DIR / "manifest.json"
shared.DEFAULT_CANDIDATE = PRODUCT_DIR / "candidate-representations.json"
shared.DEFAULT_ACCESS_RECORD = PRODUCT_DIR / "official-source/source-access-record.json"
shared.DEFAULT_APPROVAL = ROOT / "pipeline/decisions/bed02-drawing-approval.json"
shared.PSET_NAME = "Pset_Bed02DrawingSource"
shared.DOCUMENT_ID_PREFIX = "BAXTER-VIKTOR-BED02-"


def candidate_paths(candidate: dict, access: dict) -> dict:
    if (
        candidate.get("profile_key") != shared.PROFILE_KEY
        or candidate.get("representative_global_id") != shared.REPRESENTATIVE_GLOBAL_ID
        or candidate.get("article_number") != "Baxter Viktor / BED02"
        or candidate.get("source_kind") != shared.SOURCE_KIND
        or candidate.get("source_label_zh") != shared.SOURCE_LABEL_ZH
        or candidate.get("official_cad_used") is not False
        or candidate.get("third_party_cad_used") is not False
        or candidate.get("formal_ifc_write_allowed") is not False
        or candidate.get("review_status") != "visual_review_pending"
    ):
        raise RuntimeError("pending BED02 candidate source gate failed")
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
        or dimension_check.get("status") != "review_required_not_scaled_or_corrected"
        or access.get("scope") != SCOPE
    ):
        raise RuntimeError("BED02 official-source access record gate failed")
    paths = {}
    for view in shared.REQUIRED_VIEWS:
        item = candidate.get("views", {}).get(view, {})
        if item.get("source_kind") != shared.SOURCE_KIND or item.get("official_cad_paths_mm") != []:
            raise RuntimeError(f"BED02 {view} linework must remain geometry-derived with no official CAD paths")
        paths[view] = item.get("proxy_paths_mm", [])
        if len(paths[view]) != EXPECTED_PATH_COUNTS[view]:
            raise RuntimeError(f"BED02 {view} geometry-derived path-count gate failed")
    return paths


def add_document_associations(model, product, product_type, access_path: Path):
    documents = (
        (
            "BAXTER-VIKTOR-BED02-OFFICIAL-PRODUCT-PAGE",
            PRODUCT_PAGE,
            "Baxter official Viktor product page",
            "Manufacturer, family and nominal-dimension identity evidence only; not drawing geometry",
        ),
        (
            "BAXTER-VIKTOR-BED02-OFFICIAL-TECHNICAL-SHEET",
            TECHNICAL_SHEET,
            "Baxter official Viktor technical sheet",
            "The 1720 x 2340 x 1060 mm variant is dimension evidence only; not the source of linework",
        ),
        (
            "BAXTER-VIKTOR-BED02-OFFICIAL-ASSEMBLY-INSTRUCTIONS",
            ASSEMBLY_INSTRUCTIONS,
            "Baxter official Viktor assembly instructions",
            "Manufacturer identity and assembly evidence only; not the source of linework",
        ),
        (
            "BAXTER-VIKTOR-BED02-SOURCE-ACCESS-RECORD",
            relative(access_path),
            "Baxter Viktor CAD access and dimension record",
            f"SHA-256 {sha256(access_path)}; manufacturer login required; exact CAD not acquired; dimension difference retained",
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
        "Family": "Viktor",
        "Designer": "Draga & Aurel",
        "ProjectTypeCode": "BED02",
        "IFCTypeDescription": "VIKTOR 162x234xh106",
        "SourceProductPage": PRODUCT_PAGE,
        "SourceTechnicalSheet": TECHNICAL_SHEET,
        "SourceAssemblyInstructions": ASSEMBLY_INSTRUCTIONS,
        "Official2D3DBimStatus": "manufacturer_login_required_not_acquired",
        "OfficialCadUsed": "false",
        "ThirdPartyCadUsed": "false",
        "EvidenceScope": SCOPE,
        "ReviewStatus": "APPROVED",
        "Reviewer": approval["reviewer"],
        "ReviewDate": approval["review_date"],
        "ApprovalEvidence": approval["approval_evidence"],
        "ApprovedCandidateManifestSha256": manifest_hash,
        "ApprovedRepresentationIdentifiers": "Bed02Plan;Bed02Front;Bed02Side",
        "RepresentationGeometrySource": "proxy_paths_mm derived from the isolated representative IFC Body",
        "OfficialCadGeometryIncluded": "false",
        "CandidateRepresentations": relative(candidate_path),
        "CandidateRepresentationsSha256": sha256(candidate_path),
        "SourceAccessRecord": relative(access_path),
        "SourceAccessRecordSha256": sha256(access_path),
        "PlanPathCount": str(EXPECTED_PATH_COUNTS["plan"]),
        "FrontPathCount": str(EXPECTED_PATH_COUNTS["front"]),
        "SidePathCount": str(EXPECTED_PATH_COUNTS["side"]),
        "OfficialVariantOverallMm": json.dumps(dimensions["official_160x200_variant_overall_mm"]),
        "ProjectIfcDescriptionMm": json.dumps(dimensions["project_ifc_description_mm"]),
        "ProjectIfcBodyBoundsMm": json.dumps(dimensions["project_ifc_body_bounds_mm"]),
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
        Description="Mechanically verifiable geometry-derived drawing source, dimension difference and human approval",
        HasProperties=properties,
    )
    model.create_entity(
        "IfcRelDefinesByProperties",
        GlobalId=ifcopenshell.guid.new(),
        OwnerHistory=product.OwnerHistory,
        Name="BED02 drawing source and dimension-review properties",
        Description=None,
        RelatedObjects=[product, product_type],
        RelatingPropertyDefinition=pset,
    )


shared.candidate_paths = candidate_paths
shared.add_document_associations = add_document_associations
shared.add_source_pset = add_source_pset


if __name__ == "__main__":
    shared.main()
