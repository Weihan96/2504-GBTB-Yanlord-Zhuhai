#!/usr/bin/env python3
"""Approval-gated derived IFC writer for Geberit CleanLine 154.154.00.1."""

import json
from pathlib import Path

import ifcopenshell
import ifcopenshell.guid

from falper_sorgente_linework import ROOT, relative, sha256
import hima01_drawing_ifc as shared


SCOPE = "exact manufacturer article identity and standard dimensions; not a project shop drawing and not official CAD geometry"
PRODUCT_PAGE = "https://catalog.geberit.co.uk/en-GB/product/PRO_199058"
PRODUCT_DIR = ROOT / "output/review/highpoly-types/geberit-154-154-00"
EXPECTED_PATH_COUNTS = {"plan": 1, "front": 3, "side": 1}
OFFICIAL_FILES = {
    "PRODUCT-PAGE-ARCHIVE": "official-source/geberit-154-154-00-1-product-page.html",
    "VECTOR-EPS-TOP": "official-source/DAS_199904-top-view.eps",
    "VECTOR-EPS-FRONT": "official-source/DAS_199902-front-view.eps",
    "VECTOR-EPS-PERSPECTIVE": "official-source/DAS_199900-perspective.eps",
}

shared.REPRESENTATIVE_GLOBAL_ID = "14EazrLgP8whYZWY_yCuKy"
shared.IFC_TYPE_NAME = "Geberit 154.154.00"
shared.PROFILE_KEY = "geberit-154-154-00"
shared.SCOPE = SCOPE
shared.REPRESENTATIONS = {
    "plan": ("Geberit154154Plan", "PLAN_VIEW"),
    "front": ("Geberit154154Front", "ELEVATION_VIEW"),
    "side": ("Geberit154154Side", "ELEVATION_VIEW"),
}
shared.EXPECTED_PATH_COUNTS = EXPECTED_PATH_COUNTS
shared.PRODUCT_DIR = PRODUCT_DIR
shared.DEFAULT_MANIFEST = PRODUCT_DIR / "manifest.json"
shared.DEFAULT_CANDIDATE = PRODUCT_DIR / "candidate-representations.json"
shared.DEFAULT_ACCESS_RECORD = PRODUCT_DIR / "official-source/source-access-record.json"
shared.DEFAULT_APPROVAL = ROOT / "pipeline/decisions/geberit-154-154-00-drawing-approval.json"
shared.PSET_NAME = "Pset_Geberit154154DrawingSource"
shared.DOCUMENT_ID_PREFIX = "GEBERIT-154-154-00-1-"


def candidate_paths(candidate: dict, access: dict) -> dict:
    if (
        candidate.get("profile_key") != shared.PROFILE_KEY
        or candidate.get("representative_global_id") != shared.REPRESENTATIVE_GLOBAL_ID
        or candidate.get("article_number") != "154.154.00.1"
        or candidate.get("source_kind") != shared.SOURCE_KIND
        or candidate.get("source_label_zh") != shared.SOURCE_LABEL_ZH
        or candidate.get("official_cad_used") is not False
        or candidate.get("third_party_cad_used") is not False
        or candidate.get("formal_ifc_write_allowed") is not False
        or candidate.get("review_status") != "visual_review_pending"
    ):
        raise RuntimeError("pending Geberit 154.154.00.1 candidate source gate failed")
    drawing_source = access.get("drawing_geometry_source", {})
    official_cad = access.get("official_product_cad", {})
    if (
        access.get("resolved_article_number") != "154.154.00.1"
        or access.get("scope") != SCOPE
        or access.get("pass") is not True
        or official_cad.get("cad_drawings_field") != "undefined"
        or official_cad.get("acquired") is not False
        or official_cad.get("local_cad_files") != []
        or official_cad.get("official_vector_eps_archived") is not True
        or official_cad.get("official_vector_eps_used_as_cad_geometry") is not False
        or set(official_cad.get("native_dwg_http_statuses", {}).values()) != {404}
        or drawing_source.get("source_kind") != shared.SOURCE_KIND
        or drawing_source.get("official_cad_used") is not False
        or drawing_source.get("third_party_cad_used") is not False
        or access.get("dimension_cross_check", {}).get("pass") is not True
    ):
        raise RuntimeError("Geberit 154.154.00.1 official-source access record gate failed")
    paths = {}
    for view in shared.REQUIRED_VIEWS:
        item = candidate.get("views", {}).get(view, {})
        if item.get("source_kind") != shared.SOURCE_KIND or item.get("official_cad_paths_mm") != []:
            raise RuntimeError(f"Geberit 154.154.00.1 {view} must contain no official CAD paths")
        paths[view] = item.get("proxy_paths_mm", [])
        if len(paths[view]) != EXPECTED_PATH_COUNTS[view]:
            raise RuntimeError(f"Geberit 154.154.00.1 {view} geometry-derived path count drifted")
    return paths


def add_document_associations(model, product, product_type, access_path: Path):
    documents = [
        (
            "GEBERIT-154-154-00-1-OFFICIAL-PRODUCT-PAGE",
            PRODUCT_PAGE,
            "Geberit official product page for article 154.154.00.1",
            "Exact article identity and standard dimensions; the official page publishes no CAD drawing payload",
        ),
        (
            "GEBERIT-154-154-00-1-SOURCE-ACCESS-RECORD",
            relative(access_path),
            "Geberit 154.154.00.1 source access record",
            f"SHA-256 {sha256(access_path)}; native A/G/L/P DWG URLs returned 404; no third-party CAD used",
        ),
    ]
    labels = {
        "PRODUCT-PAGE-ARCHIVE": "Archived Geberit official product page",
        "VECTOR-EPS-TOP": "Geberit official top-view EPS",
        "VECTOR-EPS-FRONT": "Geberit official front-view EPS",
        "VECTOR-EPS-PERSPECTIVE": "Geberit official perspective EPS",
    }
    for suffix, relative_path in OFFICIAL_FILES.items():
        path = PRODUCT_DIR / relative_path
        documents.append((
            f"GEBERIT-154-154-00-1-{suffix}",
            relative(path),
            labels[suffix],
            f"SHA-256 {sha256(path)}; identity and dimension evidence only; not representation linework",
        ))
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
    source_hashes = {
        suffix: sha256(PRODUCT_DIR / path)
        for suffix, path in OFFICIAL_FILES.items()
    }
    values = {
        "SourceKind": shared.SOURCE_KIND,
        "SourceLabelZh": shared.SOURCE_LABEL_ZH,
        "Manufacturer": "Geberit",
        "Family": "CleanLine shower channel installation set",
        "ArticleNumber": "154.154.00.1",
        "ProjectIFCTypeName": shared.IFC_TYPE_NAME,
        "SourceProductPage": PRODUCT_PAGE,
        "OfficialCadStatus": "cadDrawings undefined; native A/G/L/P DWG URLs returned 404",
        "OfficialCadUsed": "false",
        "ThirdPartyCadUsed": "false",
        "OfficialVectorEPSRole": "identity and standard-dimension evidence only; not representation linework",
        "EvidenceScope": SCOPE,
        "ReviewStatus": "APPROVED",
        "Reviewer": approval["reviewer"],
        "ReviewDate": approval["review_date"],
        "ApprovalEvidence": approval["approval_evidence"],
        "ApprovedCandidateManifestSha256": manifest_hash,
        "ApprovedRepresentationIdentifiers": "Geberit154154Plan;Geberit154154Front;Geberit154154Side",
        "RepresentationGeometrySource": "proxy_paths_mm derived from the isolated representative IFC Body",
        "OfficialCadGeometryIncluded": "false",
        "CandidateRepresentations": relative(candidate_path),
        "CandidateRepresentationsSha256": sha256(candidate_path),
        "SourceAccessRecord": relative(access_path),
        "SourceAccessRecordSha256": sha256(access_path),
        "OfficialSourceFileSha256": json.dumps(source_hashes, sort_keys=True),
        "PlanPathCount": str(EXPECTED_PATH_COUNTS["plan"]),
        "FrontPathCount": str(EXPECTED_PATH_COUNTS["front"]),
        "SidePathCount": str(EXPECTED_PATH_COUNTS["side"]),
        "OfficialDimensionsMm": json.dumps(dimensions["official_dimensions_mm"]),
        "ProjectIFCBodyLocalXYZMm": json.dumps(dimensions["project_ifc_body_local_xyz_mm"]),
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
        Description="Mechanically verifiable geometry-derived source, official evidence and human approval",
        HasProperties=properties,
    )
    model.create_entity(
        "IfcRelDefinesByProperties",
        GlobalId=ifcopenshell.guid.new(),
        OwnerHistory=product.OwnerHistory,
        Name="Geberit 154.154.00.1 drawing source properties",
        Description=None,
        RelatedObjects=[product, product_type],
        RelatingPropertyDefinition=pset,
    )


shared.candidate_paths = candidate_paths
shared.add_document_associations = add_document_associations
shared.add_source_pset = add_source_pset


if __name__ == "__main__":
    shared.main()
