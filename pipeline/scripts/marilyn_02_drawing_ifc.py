#!/usr/bin/env python3
"""Approval-gated derived-IFC writer for Baxter Marilyn 02 pouf views."""

from __future__ import annotations

import json
from pathlib import Path

import marilyn_01_drawing_ifc as shared

from falper_sorgente_linework import ROOT, load_json, relative, sha256


shared.REPRESENTATIVE_GLOBAL_ID = "1THxa7p7n97w$wtLn4THjz"
shared.IFC_TYPE_NAME = "Marilyn 02"
shared.IFC_TYPE_DESCRIPTION = "Pouf with swivel base W80D62H45"
shared.PROFILE_KEY = "marilyn-02"
shared.SCOPE = "exact Baxter Marilyn pouf with swivel base 80 x 62 x 45 cm family CAD reference; not a project shop drawing"
shared.EXPECTED_PATH_COUNTS = {"plan": 10, "front": 52, "side": 56}
shared.REQUIRED_VIEWS = set(shared.EXPECTED_PATH_COUNTS)
shared.REPRESENTATIONS = {
    "plan": ("Marilyn02Plan", "PLAN_VIEW"),
    "front": ("Marilyn02Front", "ELEVATION_VIEW"),
    "side": ("Marilyn02Side", "ELEVATION_VIEW"),
}
shared.MATTE_SVG = "https://productsbook.baxter.it/models/measurements/MARIPFMN80.svg"
shared.GLOSSY_SVG = "https://productsbook.baxter.it/models/measurements/MARIPFML80.svg"
shared.PRODUCT_DIR = ROOT / "output/review/highpoly-types/marilyn-02"
shared.DEFAULT_MANIFEST = shared.PRODUCT_DIR / "manifest.json"
shared.DEFAULT_CANDIDATE = shared.PRODUCT_DIR / "candidate-representations.json"
shared.DEFAULT_LINEWORK = shared.PRODUCT_DIR / "official-native-dwg-linework.json"
shared.DEFAULT_ACCESS = shared.PRODUCT_DIR / "official-source/source-access-record.json"
shared.DEFAULT_APPROVAL = ROOT / "pipeline/decisions/marilyn-02-drawing-approval.json"
shared.PSET_NAME = "Pset_Marilyn02DrawingSource"
shared.DOCUMENT_ID_PREFIX = "BAXTER-MARILYN-02-"


def document_records(access_path: Path, linework_path: Path):
    source = shared.PRODUCT_DIR / "official-source"
    return (
        ("BAXTER-MARILYN-02-OFFICIAL-PRODUCT-PAGE", shared.PRODUCT_PAGE, "Baxter official Marilyn product page", "Current pouf identity, designer, dimensions and native package link"),
        ("BAXTER-MARILYN-02-OFFICIAL-PRODUCT-PAGE-ARCHIVE", relative(source / "baxter-marilyn-product-page.html"), "Archived Baxter Marilyn product page", f"SHA-256 {sha256(source / 'baxter-marilyn-product-page.html')}"),
        ("BAXTER-MARILYN-02-OFFICIAL-NATIVE-ZIP", shared.NATIVE_ZIP_URL, "Baxter official Marilyn native 2D/3D package", "Manufacturer package URL"),
        ("BAXTER-MARILYN-02-OFFICIAL-NATIVE-ZIP-ARCHIVE", relative(source / "Baxter_Marilyn_Armchair_2D_3D.zip"), "Archived Baxter Marilyn native 2D/3D package", f"SHA-256 {sha256(source / 'Baxter_Marilyn_Armchair_2D_3D.zip')}"),
        ("BAXTER-MARILYN-02-NATIVE-DWG", relative(source / "Marilyn_Abaco.dwg"), "Baxter official Marilyn_Abaco.dwg", f"SHA-256 {shared.DWG_SHA256}; authoritative exact-pouf Plan/Front/Side linework"),
        ("BAXTER-MARILYN-02-EXACT-3DS", relative(source / "Marilyn_pouf_80x62xh45.3ds"), "Baxter exact Marilyn pouf 3DS", f"SHA-256 {sha256(source / 'Marilyn_pouf_80x62xh45.3ds')}; exact variant identity evidence only"),
        ("BAXTER-MARILYN-02-OFFICIAL-TECHNICAL-SHEET", shared.TECHNICAL_SHEET, "Baxter official Marilyn technical sheet", "Current manufacturer family publication and exact pouf identity evidence"),
        ("BAXTER-MARILYN-02-OFFICIAL-TECHNICAL-SHEET-ARCHIVE", relative(source / "Baxter_Marilyn_current-technical-sheet.pdf"), "Archived Baxter Marilyn technical sheet", f"SHA-256 {sha256(source / 'Baxter_Marilyn_current-technical-sheet.pdf')}"),
        ("BAXTER-MARILYN-02-MATTE-MEASUREMENT-SVG", shared.MATTE_SVG, "Baxter official matte-frame pouf measurement SVG", f"Archived SHA-256 {sha256(source / 'MARIPFMN80.svg')}"),
        ("BAXTER-MARILYN-02-GLOSSY-MEASUREMENT-SVG", shared.GLOSSY_SVG, "Baxter official glossy-frame pouf measurement SVG", f"Archived SHA-256 {sha256(source / 'MARIPFML80.svg')}"),
        ("BAXTER-MARILYN-02-SOURCE-ACCESS-RECORD", relative(access_path), "Baxter Marilyn 02 source access record", f"SHA-256 {sha256(access_path)}"),
        ("BAXTER-MARILYN-02-NATIVE-LINEWORK-REGISTER", relative(linework_path), "Baxter Marilyn 02 native-DWG linework register", f"SHA-256 {sha256(linework_path)}; path counts 10/52/56"),
    )


def add_source_pset(model, product, product_type, approval: dict, manifest_path: Path, candidate_path: Path, access_path: Path, linework_path: Path):
    access = load_json(access_path)
    values = {
        "SourceKind": shared.SOURCE_KIND,
        "SourceLabelZh": shared.SOURCE_LABEL_ZH,
        "Manufacturer": "Baxter",
        "Family": "Marilyn",
        "Designer": "Draga & Aurel",
        "ModelCode": "Marilyn 02",
        "ResolvedOfficialVariant": access["resolved_variant"],
        "IFCTypeDescription": shared.IFC_TYPE_DESCRIPTION,
        "SourceProductPage": shared.PRODUCT_PAGE,
        "SourceNativePackage": shared.NATIVE_ZIP_URL,
        "SourceTechnicalSheet": shared.TECHNICAL_SHEET,
        "SourceDwgPath": relative(shared.PRODUCT_DIR / "official-source/Marilyn_Abaco.dwg"),
        "SourceDwgSha256": shared.DWG_SHA256,
        "ExactModel3dsPath": relative(shared.PRODUCT_DIR / "official-source/Marilyn_pouf_80x62xh45.3ds"),
        "ExactModel3dsSha256": access["native_cad_selection"]["exact_model_3ds_sha256"],
        "EvidenceScope": shared.SCOPE,
        "ReviewStatus": "APPROVED",
        "Reviewer": approval["reviewer"],
        "ReviewDate": approval["review_date"],
        "ApprovalEvidence": approval["approval_evidence"],
        "ApprovedCandidateManifestSha256": sha256(manifest_path),
        "ApprovedRepresentationIdentifiers": "Marilyn02Plan;Marilyn02Front;Marilyn02Side",
        "RepresentationGeometrySource": "official_native_dwg_paths_mm",
        "OfficialCadGeometryIncluded": "true",
        "ProxyGeometryIncluded": "false",
        "ThirdPartyCadUsed": "false",
        "CandidateRepresentations": relative(candidate_path),
        "CandidateRepresentationsSha256": sha256(candidate_path),
        "SourceAccessRecord": relative(access_path),
        "SourceAccessRecordSha256": sha256(access_path),
        "NativeDwgLineworkRegister": relative(linework_path),
        "NativeDwgLineworkRegisterSha256": sha256(linework_path),
        "PlanPathCount": str(shared.EXPECTED_PATH_COUNTS["plan"]),
        "FrontPathCount": str(shared.EXPECTED_PATH_COUNTS["front"]),
        "SidePathCount": str(shared.EXPECTED_PATH_COUNTS["side"]),
        "SideOrientationTransform": "mirror_x_for_ifc_yz_direction",
        "GeometryScaledToMatchIfc": "false",
        "OfficialNominalWidthDepthHeightMm": json.dumps(access["dimension_cross_check"]["official_nominal_width_depth_height_mm"]),
        "ProjectIfcBodyLocalXYZMm": json.dumps(access["dimension_cross_check"]["project_ifc_body_local_xyz_mm"]),
    }
    properties = [model.create_entity("IfcPropertySingleValue", Name=name, Description=None, NominalValue=model.create_entity("IfcText", str(value)), Unit=None) for name, value in values.items()]
    pset = model.create_entity("IfcPropertySet", GlobalId=shared.ifcopenshell.guid.new(), OwnerHistory=product.OwnerHistory, Name=shared.PSET_NAME, Description="Mechanically verifiable Baxter exact-pouf native-DWG source and human approval", HasProperties=properties)
    model.create_entity("IfcRelDefinesByProperties", GlobalId=shared.ifcopenshell.guid.new(), OwnerHistory=product.OwnerHistory, Name="Marilyn 02 drawing source properties", Description=None, RelatedObjects=[product, product_type], RelatingPropertyDefinition=pset)


shared.document_records = document_records
shared.add_source_pset = add_source_pset


if __name__ == "__main__":
    shared.main()
