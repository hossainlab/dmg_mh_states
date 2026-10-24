#!/usr/bin/env python
"""m06a — the human reference annotation, and mouse and human annotation in one embedding.

Answers the collaborator's "show human annotation first": panel (a) is the human atlas as annotated, with the four
malignant subpopulations (OPC-like, AC-like, OC-like, MES-like) and the cycling state shown as such rather than as one
"malignant" block; panels (b-d) put the mouse cells on the identical axes, coloured by the author annotation, by the
lineage the human reference model transfers to them, and - for tumour cells only - by which human subpopulation they
score highest for.

Inputs : ../data/processed/dmg_atlas_v1.h5ad (X_umap, lineage/cell-type annotation)
         data/human/malignant_state_assignment.csv.gz (human per-cell subpopulation; 167,012 malignant cells)
         data/human/metaprogram_families.csv, data/human/metaprograms.csv (which metaprograms map to each family)
         results/tables/m06_projection_per_cell.csv.gz (mouse kNN placement + transferred lineage)
         results/tables/m05_state_scores.csv.gz (mouse per-cell human-state scores)
         dmg_atlas.states.STATE_GENE_SETS (the marker sets that define the subpopulations)
Outputs: results/tables/m06a_human_subpopulations.csv   the human reference panel: cells, share, samples, markers, MPs
         results/tables/m06a_joint_coordinates.csv.gz   per cell of both species: annotation + embedding coordinates
         results/figures/m06a_joint_embedding.png
Run: nice -n 10 python scripts/m06a_joint_embedding.py [--force]

The mouse coordinates are the m06 kNN placement (mean X_umap of the 15 nearest reference cells in the scANVI latent):
a placement for display on the human axes, not a joint re-embedding, so mouse cells cannot form their own territory.
Human subpopulation shares are given under both assignment rules: 'state_raw' (argmax of the raw marker scores) and
'state' (argmax of the within-assay-class z-scores, the rule the atlas uses), because the two differ substantially.
"""
import sys, argparse
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "src"), str(ROOT.parent / "src")]
import dmg_atlas.utils as U  # noqa: E402
import numpy as np, pandas as pd, anndata as ad  # noqa: E402
from dmg_atlas.states import STATE_GENE_SETS  # noqa: E402
from mouse_cmp import RESULTS, LOGS, HUMAN, ATLAS_ROOT  # noqa: E402
from mouse_cmp.common import provenance, skip_existing  # noqa: E402

LOG = U.StepLogger(LOGS / "m06a_joint_embedding.log", context="m06a_joint_embedding")
STATES = ["OPC_like", "AC_like", "OC_like", "MES_like", "cycling"]
OUT = {"panel": RESULTS / "tables/m06a_human_subpopulations.csv",
       "coords": RESULTS / "tables/m06a_joint_coordinates.csv.gz",
       "fig": RESULTS / "figures/m06a_joint_embedding.png"}


def main(force):
    if skip_existing(OUT.values(), force, LOG):
        return
    with LOG.step("human annotation + subpopulations"):
        H = ad.read_h5ad(ATLAS_ROOT / "data/processed/dmg_atlas_v1.h5ad", backed="r")
        ho = H.obs[["lineage_label", "cell_type_atlas", "sample_id", "patient_id", "malignant_selected"]].copy()
        ho[["umap_1", "umap_2"]] = H.obsm["X_umap"][:]
        st = pd.read_csv(HUMAN / "malignant_state_assignment.csv.gz", index_col=0)
        ho["state_raw"] = st["state_raw"].reindex(ho.index)
        ho["state"] = st["state"].reindex(ho.index)
        # the annotation shown in panel (a): subpopulation for malignant cells, lineage for everything else
        ho["annotation"] = np.where(ho["state_raw"].notna(), ho["state_raw"].astype(str), ho["lineage_label"].astype(str))

        fam = pd.read_csv(HUMAN / "metaprogram_families.csv")
        mps = pd.read_csv(HUMAN / "metaprograms.csv")[["metaprogram", "n_samples", "replicated", "top_genes"]]
        fam = fam.merge(mps, on="metaprogram", how="left", suffixes=("", "_mp"))
        rows = []
        for s in STATES:
            m = ho["state_raw"].eq(s)
            mz = ho["state"].eq(s)
            f = s.replace("_like", "_like") if s != "cycling" else "cycling"
            members = fam[fam.family.astype(str) == f]
            rows.append({
                "subpopulation": s,
                "marker_genes_n": len(STATE_GENE_SETS[s]),
                "marker_genes": ";".join(STATE_GENE_SETS[s]),
                "n_cells_raw_rule": int(m.sum()),
                "pct_of_malignant_raw_rule": round(100 * m.sum() / ho["state_raw"].notna().sum(), 2),
                "n_cells_z_rule": int(mz.sum()),
                "pct_of_malignant_z_rule": round(100 * mz.sum() / ho["state"].notna().sum(), 2),
                "n_samples": int(ho.loc[m, "sample_id"].nunique()),
                "n_patients": int(ho.loc[m, "patient_id"].nunique()),
                "metaprograms_in_family": ";".join(members.metaprogram.astype(str)),
                "replicated_metaprograms": ";".join(members.loc[members.replicated.astype(str) == "True", "metaprogram"].astype(str)),
            })
        panel = pd.DataFrame(rows)
        panel.to_csv(OUT["panel"], index=False)
        LOG.write("human subpopulations:\n" + panel[["subpopulation", "n_cells_raw_rule", "pct_of_malignant_raw_rule",
                                                     "n_cells_z_rule", "pct_of_malignant_z_rule", "n_samples",
                                                     "metaprograms_in_family"]].to_string(index=False))

    with LOG.step("mouse placement + per-cell coordinates"):
        mo = pd.read_csv(RESULTS / "tables/m06_projection_per_cell.csv.gz", index_col=0)
        ms = pd.read_csv(RESULTS / "tables/m05_state_scores.csv.gz")
        ms = ms[ms.definition == "primary"].set_index("cell_id")
        mo["mouse_state_raw"] = ms["state_raw"].reindex(mo.index)
        coords = pd.concat([
            pd.DataFrame({"species": "human", "annotation": ho["annotation"].values,
                          "detail": ho["cell_type_atlas"].astype(str).values,
                          "umap_1": ho["umap_1"].values, "umap_2": ho["umap_2"].values}, index=ho.index),
            pd.DataFrame({"species": "mouse", "annotation": mo["author_label"].astype(str).values,
                          "detail": mo["predicted_lineage"].astype(str).values + "|" + mo["mouse_state_raw"].astype(str).values,
                          "umap_1": mo["umap_ref_1"].values, "umap_2": mo["umap_ref_2"].values}, index=mo.index)])
        coords.index.name = "cell_id"
        coords.to_csv(OUT["coords"])
        LOG.write(f"coordinates: {len(ho):,} human + {len(mo):,} mouse cells")

    with LOG.step("figure"):
        figure(ho, mo)
    provenance(list(OUT.values()),
               [ATLAS_ROOT / "data/processed/dmg_atlas_v1.h5ad", HUMAN / "malignant_state_assignment.csv.gz",
                RESULTS / "tables/m06_projection_per_cell.csv.gz", RESULTS / "tables/m05_state_scores.csv.gz"],
               "m06a_joint_embedding.py",
               extra={"mouse_coordinates": "m06 kNN placement (mean X_umap of 15 nearest reference cells), display only"})
    LOG.write(f"done; peak RSS {U.peak_rss_gb():.2f} GB")


def figure(ho, mo):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    sys.path.insert(0, str(ATLAS_ROOT / "src"))
    from dmg_atlas.figstyle import apply_figure_style, letter
    apply_figure_style()

    # Legends sit BELOW each panel with space reserved explicitly (fig.subplots_adjust, no bbox_inches='tight'),
    # labels are short and carry no counts - the counts live in the caption and in m06a_human_subpopulations.csv.
    SHORT = {"Neuronal-like tumor": "Neuronal-like", "Purkinje-like tumor": "Purkinje-like",
             "OPC-like tumor cells": "OPC-like", "Astrocyte-like tumor cells": "Astrocyte-like",
             "Meningeal fibroblasts (VLMC)": "Meningeal fibro.", "Choroid plexus epithelial cells": "Choroid plexus",
             "Mono/Macro": "Mono/Macro", "B/Plasma": "B/Plasma", "Vascular cells": "Vascular",
             "T cells": "T cells", "NK cells": "NK cells", "Microglia": "Microglia",
             "Neutrophils": "Neutrophils", "pDC": "pDC"}
    STATE_COL = {"OPC_like": "#de2d26", "AC_like": "#fd8d3c", "OC_like": "#8c510a", "MES_like": "#762a83",
                 "cycling": "#111111", "unassigned": "#cccccc"}
    # Malignant subpopulations keep the warm/black family (STATE_COL); microenvironment lineages use cool colours, so
    # that normal astrocyte is not confusable with AC-like and oligodendrocyte is not confusable with MES-like - the two
    # pairs this report contrasts.
    LIN_COL = {"myeloid": "#2b6cb0", "T_NK": "#38a169", "oligodendrocyte": "#00b8d4", "astrocyte": "#0f766e",
               "neuron": "#718096", "vascular": "#d53f8c", "Unknown": "#e2e8f0", "malignant": "#9aa5b1"}
    PRETTY = {"OPC_like": "OPC-like", "AC_like": "AC-like", "OC_like": "OC-like", "MES_like": "MES-like",
              "cycling": "cycling", "unassigned": "no state call", "malignant": "malignant, no state call",
              "T_NK": "T/NK", "Unknown": "unannotated"}

    def legend_below(ax, ncol):
        ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.015), ncol=ncol, fontsize=7, markerscale=9,
                  frameon=False, handletextpad=0.25, columnspacing=0.9, labelspacing=0.35, borderpad=0.0)

    xy = ho[["umap_1", "umap_2"]].values
    fig, axes = plt.subplots(2, 2, figsize=(11.5, 12.2))

    # (a) human annotation: malignant subpopulations drawn last so they sit on top
    ax = axes[0, 0]
    cats = list(ho["annotation"].value_counts().index)
    for lab in sorted(cats, key=lambda c: c in STATE_COL):
        g = ho[ho["annotation"].eq(lab)]
        ax.scatter(g.umap_1, g.umap_2, s=0.35, lw=0, rasterized=True,
                   color=STATE_COL.get(lab, LIN_COL.get(lab, "0.7")), label=PRETTY.get(lab, lab))
    ax.set(xticks=[], yticks=[])
    legend_below(ax, 4)
    letter(ax, "a", title="Human atlas annotation")

    # (b) mouse author annotation
    ax = axes[0, 1]
    ax.scatter(xy[:, 0], xy[:, 1], s=0.3, c="0.90", lw=0, rasterized=True)
    labs = mo["author_label"].astype(str)
    cats = list(labs.value_counts().index)
    pal = dict(zip(cats, (list(plt.cm.tab20.colors) * 2)[:len(cats)]))
    for lab in cats:
        g = mo[labs.eq(lab)]
        ax.scatter(g.umap_ref_1, g.umap_ref_2, s=1.1 if len(g) > 2000 else 3.0, color=pal[lab], lw=0,
                   rasterized=True, label=SHORT.get(lab, lab))
    ax.set(xticks=[], yticks=[])
    legend_below(ax, 4)
    letter(ax, "b", title="Mouse, author annotation")

    # (c) transferred lineage
    ax = axes[1, 0]
    ax.scatter(xy[:, 0], xy[:, 1], s=0.3, c="0.90", lw=0, rasterized=True)
    pl = mo["predicted_lineage"].astype(str)
    for lab in list(pl.value_counts().index):
        g = mo[pl.eq(lab)]
        ax.scatter(g.umap_ref_1, g.umap_ref_2, s=1.1, color=LIN_COL.get(lab, "0.5"), lw=0, rasterized=True,
                   label=PRETTY.get(lab, lab))
    ax.set(xticks=[], yticks=[])
    legend_below(ax, 4)
    letter(ax, "c", title="Mouse, lineage transferred from the human reference")

    # (d) mouse tumour cells by best-scoring human subpopulation
    ax = axes[1, 1]
    ax.scatter(xy[:, 0], xy[:, 1], s=0.3, c="0.90", lw=0, rasterized=True)
    st = mo["mouse_state_raw"].astype(str)
    for lab in [c for c in ["OPC_like", "cycling", "AC_like", "MES_like", "OC_like", "unassigned"] if (st == c).any()]:
        g = mo[st.eq(lab)]
        ax.scatter(g.umap_ref_1, g.umap_ref_2, s=4.5, color=STATE_COL.get(lab, "0.5"), lw=0, rasterized=True,
                   label=PRETTY.get(lab, lab))
    ax.set(xticks=[], yticks=[])
    legend_below(ax, 3)
    letter(ax, "d", title="Mouse tumour cells, best-scoring human subpopulation")

    fig.subplots_adjust(left=0.02, right=0.98, top=0.955, bottom=0.055, wspace=0.06, hspace=0.30)
    fig.savefig(OUT["fig"], dpi=300)
    plt.close(fig)


if __name__ == "__main__":
    p = argparse.ArgumentParser(); p.add_argument("--force", action="store_true")
    main(p.parse_args().force)
