# Molteni&C 505 UP V1.LP.S modular audit

This is a read-only product review. It does not approve or write any IFC.

## Phase 1 - rules versus current geometry

- The official 2026 English datasheet defines horizontal module centre widths `320 / 480 / 640 / 960 / 1280 / 1440 / 1920 mm`, a vertical pitch of `384 mm`, `32 mm` structural panels, and structure depths `320 / 400 mm`.
- The project width `2912 mm` decomposes exactly as `640 + 1280 + 960 + 32 end panel`. Clear widths are `608 / 1248 / 928 mm`.
- The Body height `2380 mm` decomposes exactly as `76 mm plinth + 6 × 384 mm`. The datasheet adds a `5 mm` floor clearance, so a floor-standing system is `2385 mm` from FFL. The IFC placement starts at world Z `0`; therefore the body is 5 mm too low unless another project datum supplies that clearance.
- The central opening is a valid `W1280` / clear `1248` / `H3` clear `1120` construction.
- The right protruding open element is approximately `610.000 × 766.000 mm` in front, close to the official `W608 × H768` atom. Its actual depth is `420.462 mm`, not the official open DISPLAY `D500`; the deficit is `79.538 mm`.
- The left grille is 14 separate slats, each `22 × 25 × 1528 mm`, at a `44 mm` pitch. An exact native-DWG Stripe/Harry's Bar component has not been mechanically matched, so these lines remain geometry-derived and must not be shown as verified official blue geometry.
- The model includes 12 mm central dividers at rule-aligned positions. External geometry cannot prove that required shelves are internally reinforced; this remains a specification check.

Conclusion: the overall module grid is compliant, but the current design is **not yet fully compliant** because of the DISPLAY depth and the 5 mm floor-clearance placement. The grille identity is unresolved.

## Phase 2 - atomic decomposition

| Atom | Project position / size | Official evidence | Result |
| --- | --- | --- | --- |
| End/side panels | 32 deep-panel thickness at x 0, 640, 1920, 2880 | Datasheet pp. 2, 5-6 | exact |
| Left column | W640 / clear 608 / D320 | Datasheet p. 2 | exact shell |
| Left grille | 14 × 22 × 25 × 1528, pitch 44 | Stripe/Harry's Bar family reference only | unresolved; no blue substitution |
| Centre column | W1280 / clear 1248 / H3 clear 1120 / D320 | Datasheet pp. 10, 18 | exact rule geometry |
| Right column | W960 / clear 928 / D320 | Datasheet pp. 2, 9 | exact shell |
| Right open DISPLAY | project ≈ 610 × 766 × 420.462; target 608 × 768 × 500 | Datasheet p. 23 + native DWG block `154738` | front close, depth fails |

The extracted DWG atom is placed with translation only. It is not stretched to the project mesh. The composite is labelled `module-composed official-component candidate`, not an official complete shop drawing.

## Phase 3 - parameter tool boundary

`505-up-parametric-schema.json` defines the allowed module widths, 384 mm vertical pitch, panel thicknesses, depth families, placement rules, and source provenance. `505-up-v1-lp-s-modular-composition.json` is the current data instance. `pipeline/scripts/generate_505_up_modular_review.py` validates it and reproduces this audit and the SVGs. This prototype intentionally stops before IFC authoring.
