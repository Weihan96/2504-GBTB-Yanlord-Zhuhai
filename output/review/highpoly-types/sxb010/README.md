# sxb010 Venetian-blind review

This folder is a visual-review candidate for the single IFC object `1O9JRXCI56VRUbpuLJy86Z`. It does not modify the formal IFC and contains no derived drawing IFC while approval is pending.

- `plan.svg`, `front.svg`, `side.svg`: grey actual IFC Body projection plus black geometry-derived simplified linework. No blue official-CAD layer is present.
- `project-context-ffl-plan.svg`: black semantic Plan proxy on the complete FFL Plan, with a white mask and the original walls, furniture and grids retained. The `official-elevation-anchor` annotation groups are suppressed in this review output because the EL-04 marker cluster obscured the blind; no project geometry is removed.
- `project-context-r07-front-elevation.svg`, `project-context-r07-side-elevation.svg`: black semantic Front and Side proxies on the native R07 project elevations. Their positions come from the formal IFC Body world bbox projected through each SVG's native `ifc:plane` and `ifc:matrix3`; there is no non-uniform fitting or blank-space guessing.
- `bonsai-camera-*.png`, `sxb010-bonsai-review.blend`: actual isolated IFC Body loaded by Bonsai, with four orthographic/isometric cameras and executed Blender Workbench renders.
- `official-source/`: archived Hunter Douglas official product page, official PDFs and the mechanically filtered official-download response.

The legacy proxy axes are X=width, Y=height and Z=depth. Semantic views are therefore Plan=XZ, Front=XY and Side=ZY. The generic Bonsai filenames are mapped explicitly in `bonsai-review-manifest.json`; the full blind-face render is semantic Front, not Plan.

Owner evidence confirms the Hunter Douglas 25 mm aluminium Venetian-blind family and matte graphite/grey-black direction only. Exact SKU, control system, finish code, operating side and field-measured shop dimensions remain pending. The official 1000 mm brochure example is not used as project drawing geometry. The brochure's printed pages 2-3 are PDF file pages 4-5; the source record preserves both page-number systems.
