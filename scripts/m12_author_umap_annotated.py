#!/usr/bin/env python
"""m12 - the authors' own DIPG1 UMAP, annotated with this study's evidence.

The deposit ships the authors' embedding (obs UMAP_1/UMAP_2, UMAP_source 'umap'), so their published overview of all
67,313 cells can be reproduced exactly and then overlaid with what the copy-number and reference-projection analyses
say about the same cells. The purpose is to place one question on their figure: which of the clusters labelled "tumor"
actually carry the tumour clone.

Inputs : data/processed/mouse_all_hs.h5ad (author UMAP + labels), data/processed/mouse_tumour_calls.csv (CNV clone),
         results/tables/m06_projection_per_cell.csv.gz (lineage transferred from the human reference)
Outputs: results/tables/m12_author_umap_summary.csv, results/figures/m12_author_umap_annotated.png (+ .pdf)
Run: nice -n 10 python scripts/m12_author_umap_annotated.py [--force]
"""
import sys, argparse
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
import os as _os  # noqa: E402
_ATLAS = _os.environ.get("DMG_ATLAS_ROOT") or (ROOT.parent / "dmg_atlas" if (ROOT.parent / "dmg_atlas" / "src").is_dir() else ROOT.parent)
sys.path[:0] = [str(ROOT / "src"), str(Path(_ATLAS) / "src")]
import dmg_atlas.utils as U  # noqa: E402
import numpy as np, pandas as pd, h5py  # noqa: E402
from mouse_cmp import PROCESSED, RESULTS, LOGS  # noqa: E402
from mouse_cmp.figio import PAGE_W, save_figure, panel, pretty, annotate_clusters  # noqa: E402
from mouse_cmp import palette as PAL  # noqa: E402
from mouse_cmp.common import provenance, skip_existing  # noqa: E402

LOG = U.StepLogger(LOGS / "m12_author_umap_annotated.log", context="m12_author_umap_annotated")
IN = PROCESSED / "mouse_all_hs.h5ad"
CALLS = PROCESSED / "mouse_tumour_calls.csv"
PROJ = RESULTS / "tables/m06_projection_per_cell.csv.gz"
OUT = {"summary": RESULTS / "tables/m12_author_umap_summary.csv",
       "fig": RESULTS / "figures/m12_author_umap_annotated.png"}
LIN_COL = PAL.LINEAGE


def read_obs():
    """Author UMAP and labels straight from the deposit, without loading the expression matrix."""
    with h5py.File(IN, "r") as f:
        o = f["obs"]

        def cat(name):
            g = o[name]
            if isinstance(g, h5py.Group):
                levels = [c.decode() for c in g["categories"][:]]
                return pd.Categorical.from_codes(g["codes"][:], levels).astype(str)
            v = g[:]
            return np.array([x.decode() if isinstance(x, bytes) else x for x in v])  # h5py returns bytes

        return pd.DataFrame({"cell_id": cat("cell_id"), "author_label": cat("author_label"),
                             "umap_1": o["UMAP_1"][:], "umap_2": o["UMAP_2"][:]})


def main(force=False):
    if skip_existing(OUT.values(), force, LOG):
        return
    with LOG.step("assemble the authors' embedding with our calls"):
        d = read_obs()
        calls = pd.read_csv(CALLS).set_index("cell_id")
        proj = pd.read_csv(PROJ).set_index("cell_id")
        d = d.join(calls[["clone_r", "qc_pass", "tumour_call", "author_label"]].rename(
            columns={"author_label": "al2"}), on="cell_id")
        d = d.join(proj[["predicted_lineage", "prob_malignant"]], on="cell_id")
        d["defensible_tumour"] = (d.tumour_call == "tumour") & d.qc_pass.fillna(False)
        LOG.write(f"{len(d):,} cells, {d.author_label.nunique()} author labels, "
                  f"{int(d.defensible_tumour.sum())} clone-defined QC-pass tumour cells")

    with LOG.step("per-label summary"):
        s = (d.groupby("author_label")
             .agg(n_cells=("cell_id", "size"), median_clone_r=("clone_r", "median"),
                  n_defensible_tumour=("defensible_tumour", "sum"),
                  frac_predicted_malignant=("predicted_lineage", lambda x: (x == "malignant").mean()),
                  top_transferred_lineage=("predicted_lineage", lambda x: x.value_counts().index[0]
                                           if x.notna().any() else "n/a"))
             .sort_values("n_cells", ascending=False).reset_index())
        s["frac_of_defensible_tumour"] = s.n_defensible_tumour / max(1, int(d.defensible_tumour.sum()))
        s.round(4).to_csv(OUT["summary"], index=False)
        LOG.write("author labels against our evidence:\n" + s.round(3).to_string(index=False))
    figure(d, s)
    provenance(list(OUT.values()), [IN, CALLS, PROJ], "scripts/m12_author_umap_annotated.py",
               extra={"n_cells": int(len(d)), "n_defensible_tumour": int(d.defensible_tumour.sum())})


def figure(d, s):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    sys.path.insert(0, str(Path(_ATLAS) / "src"))
    from dmg_atlas.figstyle import apply_figure_style
    apply_figure_style()
    fig, axes = plt.subplots(2, 2, figsize=(PAGE_W, 7.0))
    xy = d[["umap_1", "umap_2"]].values

    ax = axes[0, 0]
    labs = sorted(d.author_label.unique())
    # hue says which compartment the authors assigned, lightness separates labels inside it
    acol = PAL.by_compartment(labs, PAL.AUTHOR_COMPARTMENT)
    for l in labs:
        g = d[d.author_label == l]
        ax.scatter(g.umap_1, g.umap_2, s=0.6, lw=0, color=acol[l], rasterized=True)
    ax.set(xticks=[], yticks=[])
    annotate_clusters(ax, xy, d.author_label.values, acol, fontsize=5.5, min_label=250)
    panel(ax, "a", "Author annotation, as published")

    ax = axes[0, 1]
    v = d.clone_r.values
    ok = np.isfinite(v)
    ax.scatter(xy[~ok, 0], xy[~ok, 1], s=0.5, lw=0, color="0.9", rasterized=True)
    sc_ = ax.scatter(xy[ok, 0], xy[ok, 1], c=v[ok], s=0.6, lw=0, cmap="RdBu_r", vmin=-0.4, vmax=0.4,
                     rasterized=True)
    ax.set(xticks=[], yticks=[])
    cb = plt.colorbar(sc_, ax=ax, shrink=0.7)
    cb.set_label("clone correlation", fontsize=6)
    panel(ax, "b", "CNV clone signal")

    ax = axes[1, 0]
    ax.scatter(xy[:, 0], xy[:, 1], s=0.5, lw=0, color="0.88", rasterized=True)
    g = d[d.defensible_tumour]
    ax.scatter(g.umap_1, g.umap_2, s=3.5, lw=0.15, edgecolor="white", color=PAL.LINEAGE["malignant"], rasterized=True)
    ax.set(xticks=[], yticks=[])
    ax.text(0.02, 0.02, f"{len(g)} clone-defined,\nQC-pass tumour cells", transform=ax.transAxes,
            fontsize=5.5, va="bottom", ha="left",
            bbox=dict(fc="white", ec="0.8", lw=0.4, pad=1.8))
    panel(ax, "c", "Cells that pass the clone test")

    ax = axes[1, 1]
    lin = d.predicted_lineage.fillna("not projected")
    for l in ["not projected"] + sorted(set(lin) - {"not projected"}):
        g = d[lin == l]
        ax.scatter(g.umap_1, g.umap_2, s=0.6, lw=0, rasterized=True,
                   color=LIN_COL.get(l, "0.9"))
    ax.set(xticks=[], yticks=[])
    annotate_clusters(ax, xy, lin.values, LIN_COL, fontsize=5.5, min_label=300,
                      background=("not projected",), names={k: pretty(k) for k in LIN_COL})
    panel(ax, "d", "Lineage transferred from the human atlas")

    fig.subplots_adjust(left=0.05, right=0.95, top=0.945, bottom=0.02, wspace=0.08, hspace=0.08)
    save_figure(fig, OUT["fig"], LOG)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--force", action="store_true")
    main(p.parse_args().force)
