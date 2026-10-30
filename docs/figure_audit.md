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

## Colour system (added 2026-10-28)

The palettes had been declared separately inside `m06a`, `m07`, `m08`, `m10`, `m11` and `m12`, and `m05` and `m09`
carried a sixth and seventh set, so colours had drifted. The measured defects were:

| defect | evidence |
|---|---|
| the same colour meant three different things | `#de2d26` was OPC-like (m06a), the RT arm (m10) and the malignant lineage (m12); `#2b6cb0` was FUS and myeloid; `#805ad5` was FUS+RT and oligodendrocyte |
| the same entity changed colour between figures | OPC-like `#de2d26` in m06a vs `#b2182b` in m11, and **blue** `#3b6fb6` in m05; AC-like `#fd8d3c` vs `#e08214` |
| the same treatment arm changed colour between adjacent figures | in m08 RT was amber and FUS+RT red; in m10 RT was red and FUS+RT violet |
| one palette failed colour-vision testing | in m11, OPC-like and OC-like were **0.045** apart under deuteranopia simulation - indistinguishable |
| three panels coloured by palette index | `tab10`/`tab20` by position in whatever category list the script happened to build, so colours move when the category set changes |
| m09 panel a used matplotlib's default cycle | AC-like was blue there and orange everywhere else |

`src/mouse_cmp/palette.py` now holds the single source of truth. Its rules:

1. Colours come from published colour-vision-deficiency-safe qualitative sets - Okabe-Ito and Paul Tol's bright,
   muted and vibrant schemes - not from hand-picked hex or a `tab10` index.
2. Colour is not asked to carry eleven categories. Solving the eleven-colour packing exactly (binary search on the
   distance threshold plus maximum clique over the candidate sets) tops out at a worst pair of 0.183, which is not
   comfortable. So the palette has two tiers: the five malignant subpopulations are focal, saturated and separated
   as far as the space allows (worst pair 0.405), and the six microenvironment lineages are context - muted, lower
   chroma, and always directly labelled in plot. They read as a different class of thing, which is the point.
3. `cycling` takes neutral dark rather than a hue, because it is a cell-cycle state and not a lineage. That is
   semantically right and frees a hue for the four lineage-like subpopulations.
4. Assignment within a tier follows adjacency in the embedding, not preference: OPC-like, AC-like and OC-like
   border one another in figure 2a, so the two red-family colours go to OPC-like and MES-like, which do not touch.
5. Within a compartment, hue carries the compartment and lightness separates the members. `shades()` ramps over
   ABSOLUTE lightness (0.30 to 0.72), not a multiple of the base: multiplying meant a dark base could only produce
   dark variants, and the six myeloid subclusters came out as six nearly identical indigos.
6. In-plot labels are neutral dark `#1A1A1A` with a white halo, not the colour of the population they name. These
   labels sit ON their own points, so colouring the text by series fights the background - "OC-like" in sky blue
   over the red OPC-like territory was unreadable however much the text was darkened. Position already ties a label
   to its cluster. Bold is reserved for panel letters.
7. `report()` checks each family AND the unions that share a panel.

The last rule is the one that matters most, because its absence is what let three rounds of review pass over an
unreadable figure. The first palette validated SUBTYPE and LINEAGE separately and reported zero risky pairs in each
- but both are drawn in figure 2a, and across that union OC-like sat 0.063 from vascular and MES-like 0.077 from
T/NK. Validating families one at a time when they share a panel is not a weak check, it is the wrong check. The
union now reports a worst pair of 0.196, 3.1x better, and `report()` prints it on every run.

Distances throughout are the minimum Euclidean RGB distance over normal vision plus deuteranopia, protanopia and
tritanopia simulation - a pair is only as distinguishable as its worst viewing condition. Scoring the simulations
alone, as the first version did, can rate a pair safe when it is the colour-normal reader who cannot separate it.
