"""One colour system for every figure in this project.

Why this module exists: the palettes used to be declared inside m06a, m07, m10, m11 and m12 separately, and they
drifted. OPC-like was #de2d26 in one figure and #b2182b in another; #de2d26 meant OPC-like, the RT arm and the
malignant lineage depending on which figure you opened. A reader should be able to learn the colour code once.

Why it was rebuilt a second time: the first version validated each family on its own and reported no risky pairs,
but SUBTYPE and LINEAGE are drawn in the SAME panel (figure 2a), and across that union the closest pairs were
OC-like/vascular at 0.063 and MES-like/T/NK at 0.077 - indistinguishable. Validating families separately when they
share a panel is the bug that let an unreadable figure through three rounds of review.

The rules:
  * Colours come from published colour-vision-deficiency-safe qualitative sets - Okabe-Ito, and Paul Tol's bright,
    muted and vibrant schemes - not from hand-picked hex or a `tab10` index.
  * Eleven mutually distinguishable categorical colours do not exist. Solving the 11-colour packing exactly (binary
    search on the distance threshold + maximum clique) tops out at 0.183, below what is comfortable. So colour is
    not asked to carry eleven categories. There are two TIERS:
      - tier 1, the malignant subpopulations, are the focal entities: saturated, and separated as far as the
        colour space allows (worst pair 0.405, against 0.207 before).
      - tier 2, the microenvironment lineages, are context: muted, lower chroma, and always directly labelled in
        plot. They read as a different class of thing, which is the point (figure-style 4.2).
    Across the union of both tiers - what figure 2a actually shows - the worst pair is 0.196, 3.1x better than the
    0.063 it replaces.
  * `cycling` is a cell-cycle state, not a lineage, so it takes neutral dark rather than a hue. That is both
    semantically right and frees one hue for the four lineage-like subpopulations.
  * Within a compartment, hue carries the compartment and lightness separates the members (`shades`), because
    lightness is the channel that survives every form of colour-vision deficiency.
  * `report()` checks each family AND the unions that share a panel, so this class of bug cannot recur silently.

Distances are the minimum Euclidean RGB distance over normal vision plus deuteranopia, protanopia and tritanopia
simulation - a pair is only as distinguishable as its worst viewing condition.
"""
from __future__ import annotations

import itertools

# --- tier 1: malignant subpopulations - focal, saturated ----------------------------------------------------------
# Sourced from Tol vibrant (#CC3311), Okabe-Ito (#0072B2, #56B4E9) and Tol bright (#EE6677); chosen by exhaustive
# search over those sets for the largest worst-pair distance. OPC-like keeps the strongest colour because it is the
# one subpopulation this study finds conserved; MES-like takes the maximally contrasting blue because OPC-like
# against MES-like is the comparison the figures are read for.
SUBTYPE = {
    "OPC_like": "#CC3311",     # red
    "AC_like": "#0072B2",      # blue
    "MES_like": "#EE6677",     # rose
    "OC_like": "#56B4E9",      # sky
    "cycling": "#222222",      # neutral dark - a state, not a lineage
    "unassigned": "#D9D9D9",
}

# --- tier 2: microenvironment lineages - context, muted -----------------------------------------------------------
# Tol muted, plus a dark brown and a grey, chosen against tier 1 rather than on their own. Oligodendrocyte is
# deliberately far from OC-like (0.567): a malignant oligodendrocyte-like cell and a normal oligodendrocyte are the
# pair a reader must never conflate, so they are the pair pushed furthest apart.
LINEAGE = {
    "malignant": "#CC3311",        # the same red as OPC-like/tier 1 - one meaning, deliberately shared
    "myeloid": "#332288",          # indigo
    "T_NK": "#999933",             # olive
    "astrocyte": "#44AA99",        # teal
    "oligodendrocyte": "#543005",  # dark brown, never the OC-like sky
    "vascular": "#882255",         # wine
    "neuron": "#9E9E9E",           # grey: ambient context population
    "not projected": "#E8E8E8",
}

# --- treatment arms -----------------------------------------------------------------------------------------------
# Arms never share a panel with tier 1 or tier 2, so they are optimised for separation among themselves and read as
# their own system. Control is neutral grey because it is the reference level, not a treatment.
ARM = {
    "Control": "#4D4D4D",
    "FUS": "#88CCEE",
    "RT": "#E69F00",
    "FUS+RT": "#1B9E77",
}

# --- compartment base hues for subcluster shading ------------------------------------------------------------------
COMPARTMENT = {"myeloid": "#332288", "T_NK": "#999933", "vascular": "#882255", "tumour": "#CC3311",
               "other": "#737373"}

# author labels in the deposit -> the compartment whose hue family they are shaded within
AUTHOR_COMPARTMENT = {
    "OPC-like tumor cells": "tumour", "Neuronal-like tumor": "tumour", "Purkinje-like tumor": "tumour",
    "Astrocyte-like tumor cells": "tumour",
    "Mono/Macro": "myeloid", "Microglia": "myeloid", "pDC": "myeloid", "Neutrophils": "myeloid",
    "T cells": "T_NK", "NK cells": "T_NK", "B/Plasma": "T_NK",
    "Vascular cells": "vascular", "Meningeal fibroblasts (VLMC)": "vascular",
    "Choroid plexus epithelial cells": "other",
}


def _luminance(hex_colour: str) -> float:
    """WCAG relative luminance, used to keep categorical colours dark enough to read on white."""
    import matplotlib.colors as mc
    r, g, b = mc.to_rgb(hex_colour)
    f = lambda c: c / 12.92 if c <= 0.04045 else ((c + 0.055) / 1.055) ** 2.4
    return 0.2126 * f(r) + 0.7152 * f(g) + 0.0722 * f(b)


# Paul Tol's qualitative sets, which are built to stay separable under the common colour-vision deficiencies.
# The light members (#88CCEE sky, #DDCC77 sand) are kept in the list and filtered out by luminance rather than
# dropped by hand, so the rule that excludes them is visible.
_TOL_QUALITATIVE = ["#332288", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#CC6677", "#882255",
                    "#AA4499", "#6699CC", "#994455", "#4477AA", "#EE6677", "#228833", "#AA3377", "#EE7733",
                    "#0077BB", "#CC3311", "#009988", "#7570B3", "#1B9E77", "#D95F02", "#66A61E", "#E7298A",
                    "#543005", "#8C510A", "#01665E", "#35978F", "#762A83", "#5AAE61", "#B2182B", "#9970AB",
                    "#C51B7D", "#4D9221", "#BF812D", "#80CDC1", "#D6604D", "#878787"]
MAX_POINT_LUMINANCE = 0.40


def denovo(n: int):
    """Colours for de novo mouse clusters - a fifth family, kept clear of the four that carry meaning.

    Cluster IDs are nominal, so the colours must be distinguishable without implying an order or a relationship
    to a human subtype, a lineage or a treatment arm. Two constraints are applied, in this order:

    1. Readability on white. These colours are used for single cells plotted as small dots on a white page, so a
       colour must be dark enough to see. An earlier version sampled viridis, which satisfied the colour-vision
       criterion below but returned its yellow endpoint (#fde725, luminance 0.84) - invisible as points and
       unreadable as a label. Candidates above MAX_POINT_LUMINANCE are therefore rejected outright.
    2. Separation under colour-vision deficiency. Among the survivors, colours are chosen greedily to maximise the
       smallest pairwise distance across deuteranopia, protanopia and tritanopia simulations.
    """
    # Kept clear of both the subtype and the arm colours. The de novo clusters share a figure with the arms, and
    # although they never share a PANEL with the subtypes, a reader who has learned "red = OPC-like" over nine
    # figures should not meet red again as "M2" in the tenth. The pool below carries enough spare hues
    # (purple, magenta, wine, brown, teal) that this exclusion does not starve it.
    # Two exclusion radii, because the two risks are not the same size. The de novo clusters share a FIGURE with
    # the treatment arms (figure 11b), so an arm colour reappearing as a cluster is a within-figure clash: radius
    # 0.15. They never share a panel with the subtypes, so reusing a subtype hue only risks a reader's memory
    # across figures: radius 0.10. A single 0.15 radius against both is infeasible - nine reserved colours leave
    # almost nothing in a CVD-safe space, and it returned an all-green set.
    grid = [c for c in _TOL_QUALITATIVE
            if _luminance(c) <= MAX_POINT_LUMINANCE
            and min(distance(c, o) for o in ARM.values()) >= 0.15
            and min(distance(c, o) for o in SUBTYPE.values()) >= 0.10]
    if n > len(grid):
        raise ValueError(f"denovo() has {len(grid)} readable colours, {n} requested")

    # Exhaustive, not greedy: greedy seeding from the most separated pair returns a set whose worst pair is
    # 0.129, where an exhaustive search over the same pool reaches 0.195. The pool is small enough that there
    # is no reason to approximate.
    chosen = max(itertools.combinations(grid, n),
                 key=lambda s: min(distance(a, b) for a, b in itertools.combinations(s, 2)))
    # Order matters as well as membership. Clusters are handed colours in index order, and neighbouring indices
    # often sit next to each other in the embedding (M2 and M3 are the two largest and share a border), so the
    # sequence is permuted to maximise the smallest distance between CONSECUTIVE entries.
    return list(max(itertools.permutations(chosen),
                    key=lambda s: min(distance(s[i], s[i + 1]) for i in range(len(s) - 1))))


def shades(base: str, n: int, lo: float = 0.30, hi: float = 0.72):
    """n variants of one hue spread over ABSOLUTE lightness - the channel every colour-vision deficiency keeps.

    `lo` and `hi` are target lightnesses, not multipliers of the base. Multiplying the base lightness made the
    spread depend on how dark the base happened to be: the myeloid base (#332288, lightness 0.33) could only
    produce shades between 0.15 and 0.48, so six myeloid subclusters came out as six nearly identical indigos.
    An absolute ramp gives every compartment the same usable spread. Saturation is eased off at the light end,
    where a fully saturated pale colour looks washed out rather than light.
    """
    import colorsys
    import matplotlib.colors as mc
    h, l, s = colorsys.rgb_to_hls(*mc.to_rgb(base))
    if n == 1:
        return [base]
    out = []
    for i in range(n):
        f = i / (n - 1)
        li = lo + (hi - lo) * f
        out.append(mc.to_hex(colorsys.hls_to_rgb(h, li, s * (1.0 - 0.25 * f))))
    return out


def by_compartment(labels, compartment_of, order=None):
    """Colour a set of subclusters so that hue says which compartment and lightness says which member."""
    labels = list(order) if order else sorted(labels)
    groups = {}
    for l in labels:
        groups.setdefault(compartment_of.get(l, "other"), []).append(l)
    pal = {}
    for comp, members in groups.items():
        for l, c in zip(members, shades(COMPARTMENT.get(comp, "#737373"), len(members))):
            pal[l] = c
    return pal


# A fixed 14-colour sequence for the mouse-of-origin panels. Keyed by sorted sample name rather than by position in
# whatever category list a script happens to build, so a mouse keeps its colour when the cell set changes.
# Deliberately muted earth and slate tones: the mouse of origin is a nuisance variable, and it should not borrow the
# saturated colours that carry meaning elsewhere (a cyan mouse read as the oligodendrocyte lineage, a red one as
# malignant). Muted also keeps these panels visually subordinate to the ones that carry the result.
_MOUSE_SEQ = ["#8d7b68", "#b5a642", "#6b7b8c", "#a08679", "#7d8471", "#9b8aa6", "#c0a080",
              "#5f6b6d", "#ab9b8e", "#737d6b", "#8a7f9e", "#a89a7b", "#6e7f8d", "#947c6d"]


def mouse_palette(names):
    """Stable colour per mouse: sorted name -> fixed sequence, so it does not shift between figures."""
    names = sorted(set(names))
    return {n: _MOUSE_SEQ[i % len(_MOUSE_SEQ)] for i, n in enumerate(names)}


def _simulate(hex_colour, kind):
    """Vienot-style linear CVD simulation, good enough to catch pairs that collapse onto each other."""
    import numpy as np
    import matplotlib.colors as mc
    M = {"deuteranopia": [[0.625, 0.375, 0.0], [0.70, 0.30, 0.0], [0.0, 0.30, 0.70]],
         "protanopia": [[0.567, 0.433, 0.0], [0.558, 0.442, 0.0], [0.0, 0.242, 0.758]],
         "tritanopia": [[0.95, 0.05, 0.0], [0.0, 0.433, 0.567], [0.0, 0.475, 0.525]]}[kind]
    return np.clip(np.array(M) @ np.array(mc.to_rgb(hex_colour)), 0, 1)


def distance(a: str, b: str) -> float:
    """Worst-case distance between two colours over normal vision and the three dichromacies.

    Normal vision is included deliberately. An earlier version scored only the CVD simulations, which can rate a
    pair as safe when it is the colour-normal reader who cannot separate them.
    """
    import numpy as np
    import matplotlib.colors as mc
    d = [float(np.linalg.norm(np.array(mc.to_rgb(a)) - np.array(mc.to_rgb(b))))]
    d += [float(np.linalg.norm(_simulate(a, k) - _simulate(b, k)))
          for k in ("deuteranopia", "protanopia", "tritanopia")]
    return min(d)


def check(palette: dict, min_distance: float = 0.15, skip=("unassigned", "not projected")):
    """Pairs in `palette` closer than `min_distance` under the worst viewing condition."""
    items = [(k, v) for k, v in palette.items() if k not in skip]
    bad = [{"a": k1, "b": k2, "distance": round(distance(v1, v2), 3)}
           for (k1, v1), (k2, v2) in itertools.combinations(items, 2)
           if distance(v1, v2) < min_distance]
    return sorted(bad, key=lambda r: r["distance"])


# Families whose colours appear in the SAME panel somewhere in the figure set, and so have to be checked as a
# union rather than one at a time. SUBTYPE and LINEAGE share figure 2a; checking them separately is what let
# OC-like/vascular (0.063) and MES-like/T-NK (0.077) ship.
CO_OCCURRING = {"SUBTYPE + LINEAGE (figure 2a)": ("SUBTYPE", "LINEAGE")}


def report(min_distance: float = 0.15):
    """Printable audit of every palette, and of every union that shares a panel."""
    fams = {"SUBTYPE": SUBTYPE, "LINEAGE": LINEAGE, "ARM": ARM, "COMPARTMENT": COMPARTMENT}
    lines = []
    for name, pal in fams.items():
        bad = check(pal, min_distance)
        worst = min((distance(a, b) for a, b in itertools.combinations(
            [v for k, v in pal.items() if k not in ("unassigned", "not projected")], 2)), default=1.0)
        lines.append(f"{name}: {len(pal)} colours, worst pair {worst:.3f}, {len(bad)} below {min_distance}"
                     + ("" if not bad else " -> " + "; ".join(f"{b['a']}/{b['b']} {b['distance']}" for b in bad)))
    for label, names in CO_OCCURRING.items():
        union = {}
        for n in names:
            for k, v in fams[n].items():
                union[f"{n[:3].lower()}:{k}"] = v
        # 'malignant' in LINEAGE is the same red as OPC-like on purpose; it is one meaning, not a clash.
        union = {k: v for k, v in union.items() if not k.endswith(":malignant")}
        bad = check(union, min_distance, skip=())
        bad = [b for b in bad if not (b["a"].endswith("unassigned") or b["b"].endswith("unassigned")
                                      or b["a"].endswith("not projected") or b["b"].endswith("not projected"))]
        worst = min(distance(a, b) for a, b in itertools.combinations(
            [v for k, v in union.items() if "unassigned" not in k and "not projected" not in k], 2))
        lines.append(f"{label}: {len(union)} colours, worst pair {worst:.3f}, {len(bad)} below {min_distance}"
                     + ("" if not bad else " -> " + "; ".join(f"{b['a']}/{b['b']} {b['distance']}" for b in bad)))
    return "\n".join(lines)


if __name__ == "__main__":
    print(report())
