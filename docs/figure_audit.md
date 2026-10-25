# Publication-readiness audit of the mouse-comparison figures

Scope: the eight figures in `results/figures/`, audited 2026-10-25 against the requirements a high-impact
journal applies at final print size (Nature and Cell Press figure guidelines: double-column width 180-183 mm,
no text below 5 pt at final size, vector text, 300 dpi minimum for raster content).

## What was wrong

| # | finding | severity |
|---|---|---|
| 1 | No PDF existed for any of the eight figures; PNG only. | blocking |
| 2 | Four figures were authored 278-315 mm wide (m06 278, m06a 291, m09 311, m10 315). Scaled to a 183 mm column, their 6-8 pt text would have landed at 4.3-5.0 pt, under the 5 pt floor. | blocking |
| 3 | `m04`, `m05`, `m07`, `m08` wrote PNGs at **200 dpi**, below the 300 dpi floor. | blocking |
| 4 | `savefig.bbox="tight"` cropped each page to its content, so saved pages drifted off the authored size (179-195 mm) and would not land at 183 mm in a journal template. | major |
| 5 | `m07` named populations in six-entry legends set at **4.5 pt**. | major |
| 6 | Panel letters collided with centred panel titles once panels were narrowed (`m06`, `m09`, `m10`: "CPathway activity"). | major |
| 7 | `m10` panel d drew two group annotations over each other and over the title, plus a legend that repeated the x-axis arm names on top of the data points. | major |
| 8 | `m06` panel a placed a seven-entry legend over the embedding it described; `m09` panel a placed its legend over the point cloud. | major |
| 9 | Machine labels reached the page: `AC_like`, `T_NK`, `T_cell_activation`, `uncertainty_knn`, `Liu2022_TableS2a_OPC-like-1`. | minor |
| 10 | `m05` panel d truncated gene strings mid-symbol (`...;HIST1H1B;CC`); `m04` panel b collided all 19 chromosome labels. | minor |
| 11 | Panel titles carried `n` and statistics (`mouse myeloid (n=11621)`, `inferred CNV, 300 random cells per group (sorted by clone_r)`), which belong in the caption. | minor |

## What changed

All eight figures are now authored at the final print size and saved through
`src/mouse_cmp/figio.save_figure()`, which writes a vector PDF plus a 300 dpi PNG and runs a geometric
legibility check on every save. The check measures the rendered text boxes and reports text that overlaps
other text or a spine, text outside the canvas, text below 5 pt, and any figure wider than the page. It
ignores two cases that a box test over-reports: angled tick labels, whose boxes overlap while their glyphs
stay clear, and colourbar labels, which legitimately sit against their own axes.

Populations are now named inside the plot (`figio.annotate_clusters`) instead of in legends in `m06` panel a
and `m07` panels a-c. Labels anchor on the densest bin of a 24x24 histogram rather than the median, which
falls outside a scattered population, and are then pushed off the spines and off each other by measuring the
rendered text box - reserving a fixed fraction of the panel is not sufficient, because a centred label is half
its own width wide.

`m09` and `m10` gained a `--refigure` flag that redraws from the saved tables without recomputing, so a
layout change costs seconds rather than the 15 minutes `m09`'s marker-gene step takes.

## Final state

| figure | page (mm) | before (mm wide) | text | raster content |
|---|---|---|---|---|
| m04_tumour_calls | 183 x 180 | 184 | vector, Type 42 | CNV heatmaps, 300 dpi |
| m05_state_comparison | 183 x 218 | 201 | vector, Type 42 | heatmaps, 300 dpi |
| m06_reference_projection | 183 x 168 | 278 | vector, Type 42 | embeddings, 300 dpi |
| m06a_joint_embedding | 183 x 183 | 291 | vector, Type 42 | embeddings, 300 dpi |
| m07_tme_comparison | 183 x 132 | 184 | vector, Type 42 | embeddings, 300 dpi |
| m08_treatment_effects | 183 x 168 | 165 | vector, Type 42 | none (fully vector) |
| m09_subpopulation_evidence | 183 x 163 | 311 | vector, Type 42 | heatmaps, 300 dpi |
| m10_treatment_programmes | 183 x 157 | 315 | vector, Type 42 | heatmaps, 300 dpi |

Every page is 182.9 mm (7.2 in) wide - one double-column width, no rescaling needed - and every height is
within the 247 mm type area. All eight report "no legibility issues". Point clouds and heatmaps are
rasterised at 300 dpi inside an otherwise vector page, which is the standard compromise for single-cell
figures: 311,618 human plus 67,313 mouse points as vector paths would produce an unopenable file.

## Not assessed

Colour-vision safety was not re-derived here; the palettes are unchanged from the earlier figure pass.
CMYK conversion is left to the journal, as is the final caption text.
