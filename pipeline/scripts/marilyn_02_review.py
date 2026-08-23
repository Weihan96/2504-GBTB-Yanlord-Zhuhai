#!/usr/bin/env python3
"""Generate the Baxter Marilyn 02 pouf native-DWG/high-poly three-view review."""

from falper_sorgente_linework import ROOT
import marilyn_01_review as shared


shared.PROFILE_KEY = "marilyn-02"
shared.OUTPUT_DIR = ROOT / "output/review/highpoly-types/marilyn-02"
shared.REGISTER = shared.OUTPUT_DIR / "profile.json"
shared.ACCESS_RECORD = shared.OUTPUT_DIR / "official-source/source-access-record.json"
shared.SOURCE_REVALIDATION = shared.OUTPUT_DIR / "official-source/official-source-revalidation.json"
shared.SOURCE_REVALIDATION_TYPE = "Marilyn 02"
shared.ADJACENT_VARIANT_GEOMETRY_FIELD = "marilyn_01_geometry_used"
shared.EXACT_3DS_EVIDENCE_FIELD = "pouf_3ds_member_sha256"
shared.LINEWORK = shared.OUTPUT_DIR / "official-native-dwg-linework.json"
shared.EXPECTED_DESCRIPTION = "Pouf with swivel base W80D62H45"
shared.EXPECTED_PATH_COUNTS = {"plan": 10, "front": 52, "side": 56}
shared.SCOPE = "exact Baxter Marilyn pouf with swivel base 80 x 62 x 45 cm family CAD reference; not a project shop drawing"
shared.ARTICLE_LABEL = "Baxter Marilyn pouf with swivel base 80 x 62 x 45 cm"
shared.COMPARISON_TOLERANCE_MM = 35.0
shared.GENERATOR = "pipeline/scripts/marilyn_02_review.py"
shared.INDEX_TITLE = "Baxter Marilyn pouf / project Marilyn 02"
shared.INDEX_DESCRIPTION = "Blue line = exact Baxter native DWG pouf cluster. Grey = actual IFC Body. Black = simplified proxy. Official product page, exact native 3DS filename, nominal dimensions and project Body all identify the 80 x 62 x 45 cm pouf."
shared.DISCLOSURE_LINES = (
    "Exact official pouf cluster: 800 / 620 / 450 mm.",
    "No scaling or adjacent Marilyn variant substitution.",
)


if __name__ == "__main__":
    shared.main()
