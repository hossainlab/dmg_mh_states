"""Build the README hero figure: one UMAP, the compartments as this study defines them.

Why this exists: the README needs a single image that states the project's finding without a caption to lean on.
The deposit's own embedding is used, so the comparison is the depositors' coordinates judged by two different
criteria; and the two clusters deposited as `Neuronal-like tumor` and `Purkinje-like tumor` are named here for
what m02/m04/m06 established them to be - cerebellar neurons carrying no clone, failing the human mitochondrial
criterion at 98 % and 90 %, and transferred to astrocyte by the human reference model. Colouring by the deposited
label instead would put a refuted claim in the largest type on the most visible page in the repository.

Fourteen deposited labels collapse to four named compartments plus grey, because a hero figure with seven
unlabelled colours cannot be read, and one saturated colour is reserved for the clone-positive set so the eye
lands on the finding rather than on the largest cluster.

Inputs : data/mouse/DIPG1_allcell_metadata.csv (author UMAP + Cell_type_v2), data/processed/mouse_tumour_calls.csv
Outputs: docs/hero_tumour_compartment.png
Run: python scripts/m16_readme_hero.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.patheffects as pe  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
META = ROOT / "data/mouse/DIPG1_allcell_metadata.csv"
CALLS = ROOT / "data/processed/mouse_tumour_calls.csv"
OUT = ROOT / "docs/hero_tumour_compartment.png"

# Deposited labels, grouped. The two neuronal sets are not tumour (m04: zero clone-positive cells).
NEURON = ["Neuronal-like tumor", "Purkinje-like tumor"]
MYELOID = ["Microglia", "Mono/Macro", "pDC", "Neutrophils"]
LYMPHOID = ["T cells", "NK cells", "B/Plasma"]

# Muted everywhere except the clone-positive set, which is the one thing the figure is for.
DOT = {"neuron": "#DCC9A6", "other": "#D9D9D9", "myeloid": "#9DBBD8",
       "lymphoid": "#BBA2CE", "tumour": "#08306B"}
# Label text takes the same hue at higher saturation: the pale marker colours are illegible as type.
TXT = {"neuron": "#8A6A33", "myeloid": "#2B6CA3", "lymphoid": "#64428A", "tumour": "#08306B"}
HALO = [pe.withStroke(linewidth=2.4, foreground="white")]


def compartments() -> pd.DataFrame:
    md = pd.read_csv(META, usecols=["cell_id", "Cell_type_v2", "UMAP_1", "UMAP_2"])
    calls = pd.read_csv(CALLS, usecols=["cell_id", "tumour_call"])
    d = md.merge(calls, on="cell_id", validate="1:1")
    comp = pd.Series("other", index=d.index)
    comp[d.Cell_type_v2.isin(LYMPHOID)] = "lymphoid"
    comp[d.Cell_type_v2.isin(MYELOID)] = "myeloid"
    comp[d.Cell_type_v2.isin(NEURON)] = "neuron"
    comp[d.tumour_call.eq("tumour")] = "tumour"      # clone beats annotation, in both directions
    return d.assign(comp=comp)


def main() -> int:
    d = compartments()
    n = d.comp.value_counts()

    # Equal aspect on a portrait embedding: size the canvas to the data so there are no dead bands.
    xr = d.UMAP_1.max() - d.UMAP_1.min()
    yr = d.UMAP_2.max() - d.UMAP_2.min()
    height = 5.6
    fig, ax = plt.subplots(figsize=(height * xr / yr, height))
    for comp, size in [("neuron", 1.3), ("other", 1.3), ("myeloid", 1.3), ("lymphoid", 1.3), ("tumour", 3.4)]:
        k = d[d.comp == comp]
        ax.scatter(k.UMAP_1, k.UMAP_2, s=size, c=DOT[comp], lw=0, rasterized=True)

    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_aspect("equal")
    ax.margins(0.03)

    # Direct labels rather than a legend, so nothing has to be read twice. One leader, for the
    # cluster that sits under the overlay and cannot carry text.
    ax.text(-3.4, 6.2, "cerebellar neurons", color=TXT["neuron"], fontsize=9.5,
            ha="center", va="center", path_effects=HALO, zorder=5)
    ax.annotate("tumour", xy=(-4.3, -3.3), xytext=(-8.8, -6.2), color=TXT["tumour"], fontsize=9.5,
                ha="left", va="center", path_effects=HALO, zorder=5,
                arrowprops=dict(arrowstyle="-", lw=0.7, color=TXT["tumour"], shrinkA=2, shrinkB=2))
    ax.text(3.8, -8.8, "myeloid", color=TXT["myeloid"], fontsize=9.5,
            ha="left", va="center", path_effects=HALO, zorder=5)
    ax.text(10.9, 1.0, "T / NK / B", color=TXT["lymphoid"], fontsize=9.5,
            ha="center", va="center", path_effects=HALO, zorder=5)

    fig.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.01)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, dpi=400)

    renderer = fig.canvas.get_renderer()
    boxes = [t.get_window_extent(renderer) for t in fig.findobj(matplotlib.text.Text)
             if t.get_text().strip() and t.get_visible()]
    overlaps = sum(a.overlaps(b) for i, a in enumerate(boxes) for b in boxes[i + 1:])
    if overlaps:
        print(f"LABEL COLLISION: {overlaps} overlapping text boxes")
        return 1

    print(f"{OUT.relative_to(ROOT)}: {len(d):,} cells - "
          + ", ".join(f"{k} {v:,}" for k, v in n.items())
          + f"; {OUT.stat().st_size / 1e6:.1f} MB")
    return 0


if __name__ == "__main__":
    sys.exit(main())
