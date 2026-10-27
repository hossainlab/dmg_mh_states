#!/usr/bin/env python
"""m11 - do the transferred human subtype labels mark distinct clusters of mouse tumour cells, and do their DE genes
validate them?

The subtype label is NOT a classifier transfer: m05 assigns it by scoring the five human subpopulation marker sets on
each mouse cell and taking the argmax, which returns a label for every cell whether or not any state actually fits.
Differential expression on those labels therefore recovers the defining markers by construction, so this script
separates the circular part of the answer from the part that carries information:

  embedding    UMAP of the 891 QC-pass tumour cells alone, coloured by subtype and, as the confound control, by mouse
               of origin - one animal contributes 45 % of the set, so apparent structure must be checked against it
  separation   adjusted Rand index and silhouette for subtype and for mouse, against a label-permutation baseline
  hold-out     the decisive test: every gene belonging to any of the five defining marker sets is removed, and the
               embedding, clustering and silhouette are recomputed. If the subtypes still separate on the remaining
               genes, the labels capture transcriptional structure beyond the genes that assigned them; if separation
               collapses to the permutation baseline, the labels are a restatement of the scoring step
  DE           Wilcoxon one-vs-rest per subtype, with every gene flagged as a defining marker of that subtype, of
               another subtype, or of none, so the validated part of each signature can be read off directly

Inputs : data/processed/mouse_all_hs.h5ad, data/processed/mouse_tumour_calls.csv, results/tables/m05_state_scores.csv.gz
Outputs: results/tables/m11_subtype_de.csv, m11_separation_metrics.csv, m11_cell_embedding.csv.gz,
         results/figures/m11_subtype_umap_de.png (+ .pdf)
Run: nice -n 10 python scripts/m11_subtype_umap_de.py [--force]
"""
import sys, argparse
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
import os as _os  # noqa: E402
_ATLAS = _os.environ.get("DMG_ATLAS_ROOT") or (ROOT.parent / "dmg_atlas" if (ROOT.parent / "dmg_atlas" / "src").is_dir() else ROOT.parent)
sys.path[:0] = [str(ROOT / "src"), str(Path(_ATLAS) / "src")]
import dmg_atlas.utils as U  # noqa: E402
import numpy as np, pandas as pd, anndata as ad, scanpy as sc  # noqa: E402
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score, silhouette_score  # noqa: E402
from dmg_atlas.states import STATE_GENE_SETS  # noqa: E402
from mouse_cmp import PROCESSED, RESULTS, LOGS  # noqa: E402
from mouse_cmp.figio import PAGE_W, save_figure, panel, pretty, annotate_clusters  # noqa: E402
from mouse_cmp import palette as PAL  # noqa: E402
from mouse_cmp.common import provenance, skip_existing  # noqa: E402

LOG = U.StepLogger(LOGS / "m11_subtype_umap_de.log", context="m11_subtype_umap_de")
STATES = ["OPC_like", "AC_like", "OC_like", "MES_like", "cycling"]
COL = PAL.SUBTYPE
SEED = 0
IN = PROCESSED / "mouse_all_hs.h5ad"
CALLS = PROCESSED / "mouse_tumour_calls.csv"
SCORES = RESULTS / "tables/m05_state_scores.csv.gz"
OUT = {"de": RESULTS / "tables/m11_subtype_de.csv",
       "qc": RESULTS / "tables/m11_subtype_qc.csv",
       "sep": RESULTS / "tables/m11_separation_metrics.csv",
       "emb": RESULTS / "tables/m11_cell_embedding.csv.gz",
       "fig": RESULTS / "figures/m11_subtype_umap_de.png"}


def embed(a, genes_drop=(), tag=""):
    """HVG -> PCA -> neighbours -> UMAP -> Leiden on the cells given, optionally with a gene set removed."""
    b = a[:, ~a.var_names.isin(set(genes_drop))].copy() if len(genes_drop) else a.copy()
    sc.pp.highly_variable_genes(b, n_top_genes=2000, flavor="seurat")
    b = b[:, b.var.highly_variable].copy()
    sc.pp.scale(b, max_value=10)
    sc.tl.pca(b, n_comps=20, svd_solver="arpack", random_state=SEED)
    sc.pp.neighbors(b, n_neighbors=15, random_state=SEED)
    sc.tl.umap(b, random_state=SEED)
    sc.tl.leiden(b, resolution=0.5, key_added="leiden", flavor="igraph", n_iterations=2, directed=False,
                 random_state=SEED)
    LOG.write(f"{tag or 'full'}: {b.n_vars} HVGs of {b.shape[0]} cells, {b.obs.leiden.nunique()} Leiden clusters")
    return b


def metrics(b, labels, sample, rng):
    """Agreement of the embedding's own clusters with the labels, and how separated the labels are in PCA space."""
    keep = labels != "unassigned"
    lab, smp, lei = labels[keep], sample[keep], b.obs.leiden.values[keep]
    pcs = b.obsm["X_pca"][keep]
    perm = [silhouette_score(pcs, rng.permutation(lab)) for _ in range(20)]
    return {"ari_leiden_vs_subtype": adjusted_rand_score(lab, lei),
            "ari_leiden_vs_mouse": adjusted_rand_score(smp, lei),
            "nmi_leiden_vs_subtype": normalized_mutual_info_score(lab, lei),
            "nmi_leiden_vs_mouse": normalized_mutual_info_score(smp, lei),
            "silhouette_subtype": silhouette_score(pcs, lab),
            "silhouette_mouse": silhouette_score(pcs, smp),
            "silhouette_subtype_permuted_mean": float(np.mean(perm)),
            "silhouette_subtype_permuted_sd": float(np.std(perm)),
            "n_cells": int(keep.sum()), "n_hvg": int(b.n_vars)}


def main(force=False):
    if skip_existing(OUT.values(), force, LOG):
        return
    rng = np.random.default_rng(SEED)
    with LOG.step("load the 891 QC-pass tumour cells"):
        calls = pd.read_csv(CALLS).set_index("cell_id")
        sco = pd.read_csv(SCORES)
        sco = sco[sco.definition == "primary"].set_index("cell_id")
        cells = sco.index[sco.index.isin(calls.index)]
        a = ad.read_h5ad(IN)
        a = a[a.obs_names.isin(cells)].copy()
        a.obs["subtype"] = sco.loc[a.obs_names, "state_raw"].values       # raw-score argmax, as in m06a panel d
        a.obs["subtype_z"] = sco.loc[a.obs_names, "state"].values          # within-assay-class z rule
        a.obs["mouse"] = calls.loc[a.obs_names, "sample"].values
        a.obs["arm"] = calls.loc[a.obs_names, "arm"].values
        a.X = a.layers["counts"].copy()
        sc.pp.normalize_total(a, target_sum=1e4)
        sc.pp.log1p(a)
        LOG.write(f"{a.shape[0]} cells x {a.shape[1]} orthologs; subtypes "
                  f"{a.obs.subtype.value_counts().to_dict()}; {a.obs.mouse.nunique()} mice, "
                  f"largest {a.obs.mouse.value_counts().iloc[0] / a.n_obs:.1%} of cells")

    with LOG.step("library complexity per subtype"):
        C = a.layers["counts"]
        a.obs["total_counts"] = np.asarray(C.sum(1)).ravel()
        a.obs["n_genes"] = np.asarray((C > 0).sum(1)).ravel()
        qc = (a.obs.groupby("subtype", observed=True)
              .agg(n_cells=("n_genes", "size"), median_genes=("n_genes", "median"),
                   median_counts=("total_counts", "median"),
                   top_mouse_share=("mouse", lambda s: s.value_counts(normalize=True).iloc[0]),
                   top_mouse=("mouse", lambda s: s.value_counts().index[0])).reset_index())
        import scipy.stats as _ss
        groups = [g.n_genes.values for _, g in a.obs[a.obs.subtype != "unassigned"].groupby("subtype", observed=True)]
        H, pk = _ss.kruskal(*groups)
        qc["kruskal_H_genes_per_cell"] = round(float(H), 2)
        qc["kruskal_p_genes_per_cell"] = pk
        qc.to_csv(OUT["qc"], index=False)
        LOG.write("library complexity by subtype (a one-vs-rest DE is confounded if these differ):\n" +
                  qc.to_string(index=False))

    with LOG.step("embedding, clustering and the marker hold-out"):
        defining = sorted({g for s in STATES for g in STATE_GENE_SETS[s]})
        present = [g for g in defining if g in a.var_names]
        full = embed(a, tag="all genes")
        held = embed(a, genes_drop=present, tag="defining markers removed")
        rows = [{"embedding": "all genes", "n_defining_genes_removed": 0, **metrics(full, a.obs.subtype.values, a.obs.mouse.values, rng)},
                {"embedding": "defining markers removed", "n_defining_genes_removed": len(present),
                 **metrics(held, a.obs.subtype.values, a.obs.mouse.values, rng)}]
        sep = pd.DataFrame(rows)
        sep.round(4).to_csv(OUT["sep"], index=False)
        LOG.write("separation:\n" + sep.round(3).to_string(index=False))

    with LOG.step("differential expression per subtype"):
        d = a[a.obs.subtype != "unassigned"].copy()
        d.obs["subtype"] = d.obs.subtype.astype("category")
        sc.tl.rank_genes_groups(d, "subtype", method="wilcoxon", tie_correct=True)
        de = sc.get.rank_genes_groups_df(d, None).rename(columns={"group": "subtype", "names": "gene"})
        own = {s: set(STATE_GENE_SETS[s]) for s in STATES}
        any_marker = set().union(*own.values())
        de["marker_of_this_subtype"] = [g in own.get(s, set()) for g, s in zip(de.gene, de.subtype)]
        de["marker_of_another_subtype"] = [(g in any_marker) and (g not in own.get(s, set()))
                                           for g, s in zip(de.gene, de.subtype)]
        de = de.sort_values(["subtype", "pvals_adj", "scores"], ascending=[True, True, False])
        de.to_csv(OUT["de"], index=False)
        sig = de[(de.pvals_adj < 0.05) & (de.logfoldchanges > 0)]
        LOG.write("up-regulated at padj<0.05 per subtype (defining markers of that subtype in brackets):\n" +
                  "\n".join(f"  {s:9} {len(g):4d} genes [{int(g.marker_of_this_subtype.sum())} defining]"
                            for s, g in sig.groupby("subtype", observed=True)))

    with LOG.step("figure + embedding table"):
        emb = pd.DataFrame({"cell_id": a.obs_names, "subtype": a.obs.subtype.values,
                            "subtype_z": a.obs.subtype_z.values, "mouse": a.obs.mouse.values,
                            "arm": a.obs.arm.values, "leiden": full.obs.leiden.values,
                            "umap_1": full.obsm["X_umap"][:, 0], "umap_2": full.obsm["X_umap"][:, 1],
                            "umap_holdout_1": held.obsm["X_umap"][:, 0],
                            "umap_holdout_2": held.obsm["X_umap"][:, 1]})
        emb.to_csv(OUT["emb"], index=False)
        figure(emb, sep, de, sig, qc, a.obs)
    provenance(list(OUT.values()), [IN, CALLS, SCORES], "scripts/m11_subtype_umap_de.py",
               extra={"n_cells": int(a.n_obs), "n_defining_markers_present": len(present), "seed": SEED})


def figure(emb, sep, de, sig, qc, obs):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    sys.path.insert(0, str(Path(_ATLAS) / "src"))
    from dmg_atlas.figstyle import apply_figure_style
    apply_figure_style()
    fig, axes = plt.subplots(2, 3, figsize=(PAGE_W, 6.3))

    ax = axes[0, 0]
    for s in ["unassigned"] + STATES:
        g = emb[emb.subtype == s]
        ax.scatter(g.umap_1, g.umap_2, s=7, lw=0, color=COL[s], alpha=0.9)
    ax.set(xticks=[], yticks=[])
    annotate_clusters(ax, emb[["umap_1", "umap_2"]].values, emb.subtype.values, COL, fontsize=6,
                      background=("unassigned",), names={s: pretty(s) for s in STATES})
    panel(ax, "a", "Transferred subtype")

    ax = axes[0, 1]
    mice = sorted(emb.mouse.unique())
    mcol = PAL.mouse_palette(mice)
    for m in mice:
        g = emb[emb.mouse == m]
        ax.scatter(g.umap_1, g.umap_2, s=7, lw=0, color=mcol[m], alpha=0.9)
    ax.set(xticks=[], yticks=[])
    annotate_clusters(ax, emb[["umap_1", "umap_2"]].values, emb.mouse.values, mcol, fontsize=6, min_label=40)
    panel(ax, "b", "Mouse of origin")

    ax = axes[0, 2]
    top = (sig[~sig.marker_of_this_subtype].groupby("subtype", observed=True).head(3)
           .sort_values(["subtype", "scores"], ascending=[True, False]))
    genes = list(dict.fromkeys(top.gene))
    M = (de[de.gene.isin(genes)].pivot_table(index="gene", columns="subtype", values="scores")
         .reindex(index=genes, columns=STATES))
    im = ax.imshow(M.values, cmap="RdBu_r", vmin=-np.nanmax(abs(M.values)), vmax=np.nanmax(abs(M.values)),
                   aspect="auto")
    ax.set(yticks=range(len(genes)), yticklabels=genes,
           xticks=range(len(STATES)), xticklabels=[pretty(s) for s in STATES])
    ax.tick_params(axis="y", labelsize=5)
    ax.tick_params(axis="x", rotation=45)
    cb = plt.colorbar(im, ax=ax, shrink=0.75)
    cb.set_label("Wilcoxon z", fontsize=6)
    panel(ax, "c", "Top genes outside the defining sets")

    ax = axes[1, 0]
    x = np.arange(2)
    w = 0.26
    s0 = sep.set_index("embedding")
    order = ["all genes", "defining markers removed"]
    ax.bar(x - w, s0.loc[order, "silhouette_subtype"], w, color=PAL.SUBTYPE["OPC_like"], label="subtype")
    ax.bar(x, s0.loc[order, "silhouette_subtype_permuted_mean"], w, color="0.65", label="subtype, labels permuted")
    ax.bar(x + w, s0.loc[order, "silhouette_mouse"], w, color=PAL.LINEAGE["myeloid"], label="mouse of origin")
    ax.axhline(0, color="k", lw=0.6)
    ax.set(xticks=x, xticklabels=["all genes", "markers\nremoved"])
    ax.set_ylabel("silhouette in PCA space", fontsize=6)
    ax.legend(fontsize=5.5, frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=3,
              handlelength=1.0, columnspacing=0.9)
    panel(ax, "d", "Is the separation real?")

    ax = axes[1, 1]
    order = [s for s in STATES if s in set(obs.subtype)]
    ax.boxplot([obs.n_genes[obs.subtype == s].values for s in order], showfliers=False, widths=0.6,
               medianprops=dict(color=PAL.SUBTYPE["OPC_like"]))
    ax.set(xticks=range(1, len(order) + 1), xticklabels=[pretty(s) for s in order], ylabel="genes detected per cell")
    ax.tick_params(axis="x", rotation=45, labelsize=5.5)
    panel(ax, "e", "Library complexity")

    ax = axes[1, 2]
    comp = (obs[obs.subtype != "unassigned"].groupby(["mouse", "subtype"], observed=True).size()
            .unstack(fill_value=0))
    comp = comp.div(comp.sum(1), axis=0).reindex(columns=order).sort_values(order[0], ascending=False)
    left = np.zeros(len(comp))
    for s in order:
        ax.barh(range(len(comp)), comp[s].values, left=left, color=COL[s], height=0.78)
        left += comp[s].values
    ax.set(yticks=range(len(comp)), yticklabels=comp.index, xlabel="share of tumour cells", xlim=(0, 1))
    ax.tick_params(axis="y", labelsize=5)
    ax.invert_yaxis()
    panel(ax, "f", "Subtype mix per mouse")

    fig.subplots_adjust(left=0.095, right=0.945, top=0.945, bottom=0.135, wspace=0.46, hspace=0.50)
    save_figure(fig, OUT["fig"], LOG)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--force", action="store_true")
    main(p.parse_args().force)
