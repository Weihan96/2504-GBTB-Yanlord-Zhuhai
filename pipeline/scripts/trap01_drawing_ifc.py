#!/usr/bin/env python3
"""Approval-gated derived-IFC writer for configured Geberit TRAP01 views."""

import json
from pathlib import Path

import ifcopenshell
import ifcopenshell.guid

from falper_sorgente_linework import ROOT, relative, sha256
import hima01_drawing_ifc as shared


SOURCE_KIND = "geometry_derived_simplified_proxy"
SOURCE_LABEL_ZH = "基于原始高模几何生成的简化图纸表达"
SCOPE = "exact Geberit 151.116.11.1 adjustable family reference; project instance is a shortened installation configuration; not a project shop drawing"
PRODUCT_PAGE = "https://catalog.geberit.us/en-US/product/PRO_185224"
PRODUCT_DIR = ROOT / "output/review/highpoly-types/trap01"
EXPECTED_PATH_COUNTS = {"plan": 16, "front": 1, "side": 8}
DWG_FILES = {
    "A": "151.116.11.1_A.dwg",
    "G": "151.116.11.1_G.dwg",
    "L": "151.116.11.1_L.dwg",
    "P": "151.116.11.1_P.dwg",
}
EVIDENCE_FILES = {
    "PRODUCT-PAGE-ARCHIVE": "geberit-PRO_185224-product-page.html",
    "PRODUCT-DATA-SHEET": "Geberit-PRO_185224-product-data-sheet.pdf",
    "INSTALLATION-INSTRUCTIONS": "966.798.00.0-installation-instructions.pdf",
    "MAINTENANCE-MANUAL": "969.459.00.0-maintenance-manual.pdf",
}

shared.REPRESENTATIVE_GLOBAL_ID = "2Ak2ma0lvBEA49UpplzUqi"
shared.IFC_TYPE_NAME = "TRAP01"
shared.PROFILE_KEY = "trap01"
shared.SOURCE_KIND = SOURCE_KIND
shared.SOURCE_LABEL_ZH = SOURCE_LABEL_ZH
shared.SCOPE = SCOPE
shared.REPRESENTATIONS = {
    "plan": ("Trap01Plan", "PLAN_VIEW"),
    "front": ("Trap01Front", "ELEVATION_VIEW"),
    "side": ("Trap01Side", "ELEVATION_VIEW"),
}
shared.EXPECTED_PATH_COUNTS = EXPECTED_PATH_COUNTS
shared.PRODUCT_DIR = PRODUCT_DIR
shared.DEFAULT_MANIFEST = PRODUCT_DIR / "manifest.json"
shared.DEFAULT_CANDIDATE = PRODUCT_DIR / "candidate-representations.json"
shared.DEFAULT_ACCESS_RECORD = PRODUCT_DIR / "official-source/source-access-record.json"
shared.DEFAULT_APPROVAL = ROOT / "pipeline/decisions/trap01-drawing-approval.json"
shared.PSET_NAME = "Pset_Trap01DrawingSource"
shared.DOCUMENT_ID_PREFIX = "GEBERIT-151-116-11-1-TRAP01-"


def candidate_paths(candidate: dict, access: dict) -> dict:
    if (
        candidate.get("profile_key") != shared.PROFILE_KEY
        or candidate.get("representative_global_id") != shared.REPRESENTATIVE_GLOBAL_ID
        or candidate.get("article_number") != "151.116.11.1"
        or candidate.get("source_kind") != SOURCE_KIND
        or candidate.get("source_label_zh") != SOURCE_LABEL_ZH
        or candidate.get("official_cad_acquired") is not True
        or candidate.get("official_cad_used") is not False
        or candidate.get("official_cad_used_as_representation") is not False
        or candidate.get("third_party_cad_used") is not False
        or candidate.get("formal_ifc_write_allowed") is not False
        or candidate.get("review_status") != "visual_review_pending"
    ):
        raise RuntimeError("pending TRAP01 candidate source gate failed")
    native = access.get("native_cad_selection", {})
    configuration = access.get("configuration_cross_check", {})
    article = access.get("article_resolution", {})
    if (
        access.get("resolved_article") != "151.116.11.1"
        or access.get("scope") != SCOPE
        or access.get("pass") is not True
        or native.get("source_kind") != "native_dwg"
        or native.get("official_cad_acquired") is not True
        or native.get("third_party_cad_used") is not False
        or native.get("pass") is not True
        or configuration.get("project_instance_is_shortened_configuration") is not True
        or configuration.get("official_default_family_paths_used_as_project_representation") is not False
        or configuration.get("geometry_scaled_or_stretched_to_match") is not False
        or article.get("selected_article_pass") is not True
    ):
        raise RuntimeError("TRAP01 official-source/configuration record gate failed")
    paths = {}
    for view in shared.REQUIRED_VIEWS:
        item = candidate.get("views", {}).get(view, {})
        if (
            item.get("source_kind") != SOURCE_KIND
            or item.get("official_cad_paths_mm") != []
            or item.get("official_family_paths_used_as_project_representation") is not False
            or item.get("geometry_scaled_or_stretched") is not False
        ):
            raise RuntimeError(f"TRAP01 {view} must use only configured Body-derived proxy paths")
        paths[view] = item.get("proxy_paths_mm", [])
        if len(paths[view]) != EXPECTED_PATH_COUNTS[view]:
            raise RuntimeError(f"TRAP01 {view} configured proxy path count drifted")
    return paths


def add_document_associations(model, product, product_type, access_path: Path):
    source = PRODUCT_DIR / "official-source"
    documents = [
        ("OFFICIAL-PRODUCT-PAGE", PRODUCT_PAGE, "Geberit official product page", "Exact family and article identity, attributes and native CAD links"),
        ("SOURCE-ACCESS-RECORD", relative(access_path), "TRAP01 official-source access record", f"SHA-256 {sha256(access_path)}; article selection and configured-instance limitation"),
        ("NATIVE-DWG-LINEWORK", relative(PRODUCT_DIR / "official-native-dwg-linework.json"), "TRAP01 native-DWG linework register", f"SHA-256 {sha256(PRODUCT_DIR / 'official-native-dwg-linework.json')}; family reference only"),
    ]
    for suffix, filename in DWG_FILES.items():
        path = source / filename
        documents.append((f"OFFICIAL-{suffix}-DWG", relative(path), f"Geberit official {suffix} native DWG", f"SHA-256 {sha256(path)}; adjustable-family reference only; not representation geometry"))
    for suffix, filename in EVIDENCE_FILES.items():
        path = source / filename
        documents.append((suffix, relative(path), f"Geberit official {filename}", f"SHA-256 {sha256(path)}; identity and configuration evidence only"))
    identifiers = []
    for suffix, location, name, description in documents:
        identification = f"GEBERIT-151-116-11-1-TRAP01-{suffix}"
        reference = model.create_entity("IfcDocumentReference", Location=location, Identification=identification, Name=name, Description=f"{description}; {SCOPE}", ReferencedDocument=None)
        model.create_entity("IfcRelAssociatesDocument", GlobalId=ifcopenshell.guid.new(), OwnerHistory=product.OwnerHistory, Name=f"{name} association", Description=SCOPE, RelatedObjects=[product, product_type], RelatingDocument=reference)
        identifiers.append(identification)
    return sorted(identifiers)


def add_source_pset(model, product, product_type, approval, manifest_hash, candidate_path, access_path):
    access = json.loads(access_path.read_text(encoding="utf-8"))
    values = {
        "SourceKind": SOURCE_KIND,
        "SourceLabelZh": SOURCE_LABEL_ZH,
        "Manufacturer": "Geberit",
        "Family": access["family"],
        "ArticleNumber": "151.116.11.1",
        "ProjectIFCTypeName": shared.IFC_TYPE_NAME,
        "SourceProductPage": PRODUCT_PAGE,
        "OfficialCadAcquired": "true",
        "OfficialCadUsed": "false",
        "OfficialCadRole": "exact article identity, adjustable-family envelope and fixed-width cross-check only",
        "ProjectConfiguration": "shortened installed configuration represented from isolated IFC Body",
        "GeometryScaledOrStretched": "false",
        "ThirdPartyCadUsed": "false",
        "EvidenceScope": SCOPE,
        "ReviewStatus": "APPROVED",
        "Reviewer": approval["reviewer"],
        "ReviewDate": approval["review_date"],
        "ApprovalEvidence": approval["approval_evidence"],
        "ApprovedCandidateManifestSha256": manifest_hash,
        "ApprovedRepresentationIdentifiers": "Trap01Plan;Trap01Front;Trap01Side",
        "RepresentationGeometrySource": "proxy_paths_mm derived from the configured isolated representative IFC Body",
        "OfficialCadGeometryIncluded": "false",
        "CandidateRepresentations": relative(candidate_path),
        "CandidateRepresentationsSha256": sha256(candidate_path),
        "SourceAccessRecord": relative(access_path),
        "SourceAccessRecordSha256": sha256(access_path),
        "OfficialNativeDwgSha256": json.dumps(access["native_cad_selection"]["sha256"], sort_keys=True),
        "PlanPathCount": str(EXPECTED_PATH_COUNTS["plan"]),
        "FrontPathCount": str(EXPECTED_PATH_COUNTS["front"]),
        "SidePathCount": str(EXPECTED_PATH_COUNTS["side"]),
        "ProjectIFCBodyLocalXYZMm": json.dumps(access["article_resolution"]["project_ifc_body_local_xyz_mm"]),
        "SelectedArticleFixedWidthDeltaMm": str(access["article_resolution"]["selected_fixed_width_absolute_delta_mm"]),
        "ExcludedAlternativeArticle": access["article_resolution"]["excluded_alternative"],
    }
    properties = [model.create_entity("IfcPropertySingleValue", Name=name, Description=None, NominalValue=model.create_entity("IfcText", str(value)), Unit=None) for name, value in values.items()]
    pset = model.create_entity("IfcPropertySet", GlobalId=ifcopenshell.guid.new(), OwnerHistory=product.OwnerHistory, Name=shared.PSET_NAME, Description="Mechanically verifiable configured Body drawing source, official CAD evidence and human approval", HasProperties=properties)
    model.create_entity("IfcRelDefinesByProperties", GlobalId=ifcopenshell.guid.new(), OwnerHistory=product.OwnerHistory, Name="TRAP01 drawing source properties", Description=None, RelatedObjects=[product, product_type], RelatingPropertyDefinition=pset)


shared.candidate_paths = candidate_paths
shared.add_document_associations = add_document_associations
shared.add_source_pset = add_source_pset


if __name__ == "__main__":
    shared.main()
