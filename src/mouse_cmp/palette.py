"""One colour system for every figure in this project.

Why this module exists: the palettes used to be declared inside m06a, m07, m10, m11 and m12 separately, and they
drifted. OPC-like was #de2d26 in one figure and #b2182b in another; #de2d26 meant OPC-like, the RT arm and the
malignant lineage depending on which figure you opened; and the m11 subtype palette put OPC-like and OC-like 0.045
apart under deuteranopia simulation, which is indistinguishable. A reader should be able to learn the colour code once.

The rules:
  * three variable families, three visually separate colour ranges - malignant subtypes are warm, microenvironment
    lineages are cool, treatment arms are neutral-to-violet. No colour appears in two families.
  * within a compartment, hue carries the compartment and lightness separates the members, because lightness is the
    channel that survives every form of colour-vision deficiency.
  * nothing is coloured by `tab10`/`tab20` index: an index palette silently re-colours everything when the category
    set changes.
  * `check(palette)` reports the closest pair under deuteranopia, protanopia and tritanopia simulation, so a future
    edit cannot quietly reintroduce an indistinguishable pair.
"""
from __future__ import annotations

import itertools

# --- malignant subpopulations: warm, plus black for cycling -------------------------------------------------------
SUBTYPE = {
    "OPC_like": "#a50f15",     # dark red
    "AC_like": "#f16913",      # orange
    "OC_like": "#a6761d",      # dark ochre: a pale tan was invisible as points on a white panel
    "MES_like": "#762a83",     # violet
    "cycling": "#111111",      # near-black
    "unassigned": "#d9d9d9",
}

# --- microenvironment lineages: cool ------------------------------------------------------------------------------
LINEAGE = {
    "malignant": "#a50f15",        # same red as the malignant family above - one meaning, deliberately shared
    "myeloid": "#2171b5",          # blue
    "T_NK": "#41ab5d",             # green
    "astrocyte": "#0f766e",        # teal, never the AC-like orange
    "oligodendrocyte": "#00b8d4",  # cyan, never the OC-like tan
    "neuron": "#b3b3b3",           # light grey: a darker grey collapsed onto the T/NK green for a deuteranope
    "vascular": "#ce1256",         # magenta
    "not projected": "#e8e8e8",
}

# --- treatment arms: neutral to violet, no overlap with the two families above -------------------------------------
ARM = {
    "Control": "#525252",
    "FUS": "#74a9cf",
    "RT": "#fdae61",
    "FUS+RT": "#8c6bb1",
}

# --- compartment base hues for subcluster shading ------------------------------------------------------------------
COMPARTMENT = {"myeloid": "#2171b5", "T_NK": "#41ab5d", "vascular": "#ce1256", "tumour": "#a50f15",
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


def denovo(n: int):
    """Colours for de novo mouse clusters - a fifth family, kept clear of the four that carry meaning.

    Cluster IDs are nominal, so they take positions along a perceptually uniform ramp (viridis) rather than a hue
    that would imply a relationship to a human subtype, a lineage or a treatment arm. The positions are not evenly
    spaced: evenly spaced samples of viridis put adjacent clusters 0.07 apart under deuteranopia simulation, so the
    positions are chosen greedily to maximise the smallest pairwise distance across all three simulations.
    """
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.colors as mc
    cm = plt.get_cmap("viridis")
    grid = [mc.to_hex(cm(x)) for x in np.linspace(0.0, 1.0, 41)]

    def d(a, b):
        return min(float(np.linalg.norm(_simulate(a, k) - _simulate(b, k)))
                   for k in ("deuteranopia", "protanopia", "tritanopia"))

    picked = [grid[0], grid[-1]]
    while len(picked) < n:
        best = max((min(d(c, q) for q in picked), c) for c in grid if c not in picked)
        picked.append(best[1])
    return picked[:n]


def shades(base: str, n: int, lo: float = 0.45, hi: float = 1.45):
    """n variants of one hue, spread over lightness - the channel every colour-vision deficiency keeps."""
    import colorsys
    import matplotlib.colors as mc
    h, l, s = colorsys.rgb_to_hls(*mc.to_rgb(base))
    if n == 1:
        return [base]
    out = []
    for i in range(n):
        f = lo + (hi - lo) * i / (n - 1)
        out.append(mc.to_hex(colorsys.hls_to_rgb(h, min(0.92, max(0.12, l * f)), s)))
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


def check(palette: dict, min_distance: float = 0.15, skip=("unassigned", "not projected")):
    """Closest pair of a palette under each CVD simulation; returns rows below `min_distance`."""
    import numpy as np
    items = [(k, v) for k, v in palette.items() if k not in skip]
    bad = []
    for kind in ("deuteranopia", "protanopia", "tritanopia"):
        for (k1, v1), (k2, v2) in itertools.combinations(items, 2):
            d = float(np.linalg.norm(_simulate(v1, kind) - _simulate(v2, kind)))
            if d < min_distance:
                bad.append({"vision": kind, "a": k1, "b": k2, "distance": round(d, 3)})
    return bad


def report():
    """Printable audit of every palette in this module."""
    lines = []
    for name, pal in (("SUBTYPE", SUBTYPE), ("LINEAGE", LINEAGE), ("ARM", ARM)):
        bad = check(pal)
        lines.append(f"{name}: {len(pal)} colours, {len(bad)} risky pair(s)"
                     + ("" if not bad else " -> " + "; ".join(f"{b['a']}/{b['b']} {b['vision'][:4]} {b['distance']}"
                                                              for b in bad)))
    return "\n".join(lines)


if __name__ == "__main__":
    print(report())
