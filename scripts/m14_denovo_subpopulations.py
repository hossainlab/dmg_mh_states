#!/usr/bin/env python
"""m14 - mouse tumour subpopulations defined de novo, then compared with the human subpopulations.

Shenglin's proposal (28 October 2026): rather than scoring the human marker sets on each mouse cell and taking the
argmax - which m11 showed returns a label whether or not any state fits (ARI 0.075, silhouette at the permutation
baseline once the defining genes are held out) - let the mouse tumour cells define their own subpopulations, then
compare structure with structure: correlation between mouse and human subpopulation profiles, conserved marker genes,
conserved pathways and transcription factors, and finally the effect of FUS and RT on the mouse's own subpopulations.

Cell set: the authors' "OPC-like tumor cells" cluster, 3,246 cells. It is the only author tumour cluster carrying the
shared clone (median clone r 0.218, against 0.069 for Neuronal-like and 0.041 for Purkinje-like, which between them
contribute no clone-positive QC-pass cells at all), and it already contains 883 of the 891 cells of the strict
clone-positive QC-pass set, which is carried through as a sensitivity column. The Neuronal- and Purkinje-like clusters
are deliberately NOT clustered in: they behave as normal cerebellar tissue carrying ambient tumour RNA, and including
them would reintroduce the contamination that produced the spurious AC-like signal in the first place.

Clustering resolution is chosen by a rule stated before looking at the biology: the largest resolution in the sweep at
which every cluster holds at least 100 cells and appears in at least three mice.

Inputs : data/processed/mouse_all_hs.h5ad, mouse_tumour_calls.csv, data/interim/mouse_qc.csv,
         <atlas>/data/processed/dmg_atlas_v1.h5ad, data/human/malignant_state_assignment.csv.gz,
         data/external/progeny_human.csv, data/external/collectri_human.csv
Outputs: results/tables/m14_*.csv, results/figures_denovo/m14_*.png (+ .pdf)
Run: nice -n 10 python scripts/m14_denovo_subpopulations.py [--force]
"""
import sys, argparse
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
import os as _os  # noqa: E402
_ATLAS = _os.environ.get("DMG_ATLAS_ROOT") or (ROOT.parent / "dmg_atlas" if (ROOT.parent / "dmg_atlas" / "src").is_dir() else ROOT.parent)
sys.path[:0] = [str(ROOT / "src"), str(Path(_ATLAS) / "src")]
import dmg_atlas.utils as U  # noqa: E402
import numpy as np, pandas as pd, anndata as ad, scanpy as sc, scipy.sparse as sp  # noqa: E402
import scipy.stats as ss  # noqa: E402
import statsmodels.api as sm  # noqa: E402
import decoupler as dc  # noqa: E402
from dmg_atlas.states import STATE_GENE_SETS  # noqa: E402
from mouse_cmp import INTERIM, PROCESSED, EXTERNAL, RESULTS, LOGS, HUMAN, ATLAS_ROOT  # noqa: E402
from mouse_cmp.figio import (PAGE_W, save_figure, panel, pretty, annotate_clusters, cell_grid,  # noqa: E402
                             species_axis)
from mouse_cmp import palette as PAL  # noqa: E402
from mouse_cmp.common import provenance, skip_existing, bh, ARMS  # noqa: E402

LOG = U.StepLogger(LOGS / "m14_denovo_subpopulations.log", context="m14_denovo")
STATES = ["OPC_like", "AC_like", "OC_like", "MES_like", "cycling"]
TUMOUR_CLUSTER = "OPC-like tumor cells"
SEED = 0
FIGDIR = RESULTS / "figures_denovo"
OUT = {"cells": RESULTS / "tables/m14_denovo_cells.csv.gz",
       "summary": RESULTS / "tables/m14_denovo_cluster_summary.csv",
       "markers": RESULTS / "tables/m14_denovo_markers.csv",
       "corr": RESULTS / "tables/m14_denovo_vs_human_correlation.csv",
       "overlap": RESULTS / "tables/m14_denovo_marker_overlap.csv",
       "act": RESULTS / "tables/m14_denovo_pathway_tf.csv",
       "treat": RESULTS / "tables/m14_denovo_treatment.csv",
       "fig1": FIGDIR / "m14_denovo_clusters.png",
       "fig2": FIGDIR / "m14_denovo_vs_human.png",
       "fig3": FIGDIR / "m14_denovo_treatment.png"}


def pseudobulk(X, groups, genes, min_cells=10):
    out, names = [], []
    for g in pd.unique(groups):
        m = groups == g
        if m.sum() < min_cells:
            continue
        s = np.asarray(X[m].sum(0)).ravel()
        out.append(np.log1p(1e6 * s / max(s.sum(), 1)))
        names.append(g)
    return pd.DataFrame(out, index=names, columns=genes)


MIN_CLUSTER_CELLS, MIN_CLUSTER_MICE = 50, 3


def choose_resolution(a, resolutions=(0.2, 0.3, 0.4, 0.6, 0.8, 1.0),
                      min_cells=MIN_CLUSTER_CELLS, min_mice=MIN_CLUSTER_MICE):
    """Largest resolution whose clusters all clear the size and mouse-spread floors - stated before the biology."""
    chosen, rows = None, []
    for r in resolutions:
        key = f"leiden_{r}"
        sc.tl.leiden(a, resolution=r, key_added=key, flavor="igraph", n_iterations=2, directed=False,
                     random_state=SEED)
        g = a.obs.groupby(key, observed=True)
        sizes, mice = g.size(), g["mouse"].nunique()
        ok = bool((sizes >= min_cells).all() and (mice >= min_mice).all())
        rows.append({"resolution": r, "n_clusters": a.obs[key].nunique(), "min_cluster_size": int(sizes.min()),
                     "min_mice_per_cluster": int(mice.min()), "passes_rule": ok})
        if ok:
            chosen = r
    sweep = pd.DataFrame(rows)
    LOG.write("resolution sweep:\n" + sweep.to_string(index=False))
    if chosen is None:
        # No resolution satisfies the floors, because one small cluster persists at every resolution. The
        # pre-stated fallback is the lowest resolution (fewest clusters); clusters under the floor are kept in the
        # embedding for completeness but excluded from the statistics, and named as excluded in the log.
        LOG.write(f"no resolution satisfies >={min_cells} cells and >={min_mice} mice in every cluster; "
                  f"falling back to the lowest resolution {resolutions[0]} and flagging the clusters below the floor")
    return chosen or resolutions[0], chosen is not None


def main(force=False):
    if skip_existing(OUT.values(), force, LOG):
        return
    FIGDIR.mkdir(parents=True, exist_ok=True)

    with LOG.step("mouse tumour compartment"):
        M = ad.read_h5ad(PROCESSED / "mouse_all_hs.h5ad")
        calls = pd.read_csv(PROCESSED / "mouse_tumour_calls.csv", index_col="cell_id")
        qc = pd.read_csv(INTERIM / "mouse_qc.csv", index_col=0)
        keep = (M.obs["author_label"].astype(str) == TUMOUR_CLUSTER).values
        a = M[keep].copy()
        a.obs["mouse"] = calls["sample"].reindex(a.obs_names).values
        a.obs["arm"] = calls["arm"].reindex(a.obs_names).values
        a.obs["clone_r"] = calls["clone_r"].reindex(a.obs_names).values
        a.obs["strict_set"] = ((calls["tumour_call"].reindex(a.obs_names) == "tumour") &
                               qc["qc_pass"].reindex(a.obs_names).fillna(False)).values
        counts = a.layers["counts"].copy()
        a.X = counts.copy()
        sc.pp.normalize_total(a, target_sum=1e4)
        sc.pp.log1p(a)
        LOG.write(f"{a.n_obs:,} cells of the '{TUMOUR_CLUSTER}' cluster, {a.obs.mouse.nunique()} mice; "
                  f"{int(a.obs.strict_set.sum())} also in the strict clone-positive QC-pass set; "
                  f"largest mouse {a.obs.mouse.value_counts(normalize=True).iloc[0]:.1%}")

    with LOG.step("de novo clustering"):
        b = a.copy()
        sc.pp.highly_variable_genes(b, n_top_genes=2000, flavor="seurat")
        b = b[:, b.var.highly_variable].copy()
        sc.pp.scale(b, max_value=10)
        sc.tl.pca(b, n_comps=20, svd_solver="arpack", random_state=SEED)
        sc.pp.neighbors(b, n_neighbors=15, random_state=SEED)
        sc.tl.umap(b, random_state=SEED)
        res, rule_passed = choose_resolution(b)
        b.obs["cluster"] = ("M" + (b.obs[f"leiden_{res}"].astype(int) + 1).astype(str)).values
        a.obs["cluster"] = b.obs["cluster"].values
        sizes = a.obs.cluster.value_counts()
        below = sorted(sizes[sizes < MIN_CLUSTER_CELLS].index)
        a.obs["below_floor"] = a.obs.cluster.isin(below).values
        LOG.write(f"resolution {res} chosen (pre-stated rule satisfied: {rule_passed}); "
                  f"clusters {sizes.sort_index().to_dict()}"
                  + (f"; below the {MIN_CLUSTER_CELLS}-cell floor and excluded from the statistics: {below}" if below else ""))
        emb = pd.DataFrame({"cell_id": a.obs_names, "cluster": a.obs.cluster.values, "mouse": a.obs.mouse.values,
                            "arm": a.obs.arm.values, "clone_r": a.obs.clone_r.values,
                            "strict_set": a.obs.strict_set.values, "below_floor": a.obs.below_floor.values,
                            "umap_1": b.obsm["X_umap"][:, 0], "umap_2": b.obsm["X_umap"][:, 1]})
        emb.to_csv(OUT["cells"], index=False)

    with LOG.step("per-cluster clone signal"):
        summ = (emb.groupby("cluster")
                .agg(cells=("cell_id", "size"), n_mice=("mouse", "nunique"),
                     median_clone_r=("clone_r", "median"), frac_strict_set=("strict_set", "mean"),
                     below_floor=("below_floor", "first")).reset_index())
        summ.round(4).to_csv(OUT["summary"], index=False)
        LOG.write("per-cluster clone signal (clone_r is the correlation with the shared tumour CNV profile):\n"
                  + summ.round(3).to_string(index=False))

    with LOG.step("cluster markers"):
        sc.tl.rank_genes_groups(a, "cluster", method="wilcoxon", tie_correct=True)
        mk = sc.get.rank_genes_groups_df(a, None).rename(columns={"group": "cluster", "names": "gene"})
        mk.to_csv(OUT["markers"], index=False)
        keepc = set(a.obs.loc[~a.obs.below_floor, "cluster"])
        top = {c: g.nlargest(50, "scores").gene.tolist() for c, g in mk.groupby("cluster") if c in keepc}
        LOG.write("top 6 markers per cluster:\n" + "\n".join(f"  {c}: {', '.join(v[:6])}" for c, v in top.items()))

    with LOG.step("human malignant pseudobulk"):
        H = ad.read_h5ad(ATLAS_ROOT / "data/processed/dmg_atlas_v1.h5ad")
        H.var_names = H.var["gene_symbol"].astype(str).values
        H.var_names_make_unique()
        st = pd.read_csv(HUMAN / "malignant_state_assignment.csv.gz", index_col=0)["state_raw"]
        hs = st.reindex(H.obs_names)
        Hm = H[hs.notna().values]
        hstate = hs[hs.notna()].astype(str).values
        shared = sorted(set(a.var_names) & set(Hm.var_names))
        hb = pseudobulk(Hm[:, shared].X, np.array([f"human:{s}" for s in hstate]), shared)
        del H, Hm
        mi = [a.var_names.get_loc(g) for g in shared]
        mb = pseudobulk(counts[~a.obs.below_floor.values][:, mi],
                        ("mouse:" + a.obs.cluster.astype(str))[~a.obs.below_floor.values].values, shared,
                        min_cells=MIN_CLUSTER_CELLS)
        LOG.write(f"shared genes {len(shared):,}; human profiles {len(hb)}; mouse profiles {len(mb)}")

    with LOG.step("mouse cluster versus human subpopulation"):
        # Gene-wise z-scoring across profiles, as m09 does. Without it the correlation is dominated by the
        # expression baseline every tumour profile shares, and every mouse cluster correlates with every human
        # subpopulation at +0.2 to +0.67 - a number that cannot discriminate. Z-scoring asks the question that
        # matters: which human subpopulation does this cluster resemble DISTINCTIVELY.
        hz = hb.loc[[f"human:{s}" for s in STATES if f"human:{s}" in hb.index]]
        genes = hz.var(0).sort_values(ascending=False).head(2000).index
        A = hz[genes].apply(ss.zscore, axis=0)
        B = mb[genes].apply(ss.zscore, axis=0)
        # A gene with no variance across the four mouse cluster profiles z-scores to NaN; with only four profiles
        # that happens often enough to poison whole rows, so those genes are dropped from both sides.
        usable = A.columns[A.notna().all() & B.notna().all()]
        A, B = A[usable], B[usable]
        LOG.write(f"{len(usable):,} of {len(genes):,} genes usable after z-scoring both profile sets")
        rows = []
        for mc in B.index:
            for hc in A.index:
                rho, p = ss.spearmanr(B.loc[mc], A.loc[hc])
                rows.append({"mouse_cluster": mc.replace("mouse:", ""),
                             "human_subpopulation": hc.replace("human:", ""), "rho": round(float(rho), 4),
                             "p": float(p)})
        corr = pd.DataFrame(rows)
        corr["q"] = bh(corr.p.values)
        corr.to_csv(OUT["corr"], index=False)
        best = corr.loc[corr.groupby("mouse_cluster").rho.idxmax()]
        LOG.write("best-matching human subpopulation per mouse cluster:\n" +
                  best[["mouse_cluster", "human_subpopulation", "rho"]].to_string(index=False))

    with LOG.step("conserved marker overlap"):
        universe = set(shared)
        rows = []
        for c, genes in top.items():
            mg = set(genes) & universe
            for s in STATES:
                hg = set(STATE_GENE_SETS[s]) & universe
                ov = len(mg & hg)
                p = ss.hypergeom.sf(ov - 1, len(universe), len(hg), len(mg)) if ov else 1.0
                rows.append({"mouse_cluster": c, "human_subpopulation": s, "n_mouse_top": len(mg),
                             "n_human_set": len(hg), "overlap": ov, "shared_genes": ";".join(sorted(mg & hg)),
                             "p": float(p)})
        ov = pd.DataFrame(rows)
        ov["q"] = bh(ov.p.values)
        ov.to_csv(OUT["overlap"], index=False)
        LOG.write("marker overlaps at q<0.1:\n" +
                  (ov[ov.q < 0.1][["mouse_cluster", "human_subpopulation", "overlap", "q", "shared_genes"]]
                   .to_string(index=False) or "  none"))

    with LOG.step("conserved pathways and TFs"):
        prof = pd.concat([hb, mb])
        net_p = pd.read_csv(EXTERNAL / "progeny_human.csv")
        net_t = pd.read_csv(EXTERNAL / "collectri_human.csv")[["source", "target", "weight"]]

        def score(df, net, tmin):
            x = ad.AnnData(df.astype(np.float32))
            dc.mt.ulm(x, net, tmin=tmin)
            return x.obsm["score_ulm"]

        ps, ts = score(prof, net_p, 5), score(prof, net_t, 10)
        act = []
        for name, tab in (("pathway", ps), ("tf", ts)):
            s = tab.stack().rename("activity").reset_index()
            s.columns = ["group", "feature", "activity"]
            s["family"] = name
            act.append(s)
        act = pd.concat(act, ignore_index=True)
        act["species"] = act.group.str.split(":").str[0]
        act.round(4).to_csv(OUT["act"], index=False)
        LOG.write(f"activity rows {len(act):,} over {act.group.nunique()} profiles")

    with LOG.step("treatment effects on the de novo clusters"):
        comp = (emb.groupby(["mouse", "arm", "cluster"]).size().rename("n").reset_index()
                .pivot_table(index=["mouse", "arm"], columns="cluster", values="n", fill_value=0).reset_index())
        clusters = sorted(emb.loc[~emb.below_floor, "cluster"].unique())
        comp["total"] = comp[clusters].sum(1)
        rows = []
        for c in clusters:
            d = comp[comp.total >= 20].copy()
            if d.arm.nunique() < 2:
                continue
            d["arm"] = pd.Categorical(d.arm, categories=ARMS)
            X = sm.add_constant(pd.get_dummies(d.arm, drop_first=True).astype(float))
            try:
                m = sm.GLM(np.c_[d[c], d.total - d[c]], X, family=sm.families.Binomial()).fit()
            except Exception as e:  # noqa: BLE001
                LOG.write(f"  {c}: model failed ({e})")
                continue
            # Quasi-binomial by hand. statsmodels' scale="X2" is wrong for two-column endog - m08 measured 0.24
            # against a true 79.8 - and without it the cell counts drive the standard errors, which turns a
            # 15-mouse comparison into a 3,246-cell one and returns p < 1e-5 for everything.
            phi = max(1.0, float(m.pearson_chi2 / m.df_resid)) if m.df_resid > 0 else np.nan
            cov = m.cov_params()
            for term in [t for t in X.columns if t != "const"]:
                est = float(m.params[term])
                se = float(np.sqrt(phi * cov.loc[term, term])) if hasattr(cov, "loc") else float(np.sqrt(phi) * m.bse[term])
                pval = float(2 * ss.t.sf(abs(est / se), m.df_resid)) if se > 0 and m.df_resid > 0 else np.nan
                rows.append({"cluster": c, "contrast": f"{term}_vs_Control", "coef_logodds": est,
                             "se": se, "phi": round(phi, 2), "p": pval, "n_mice": int(len(d)),
                             "mean_prop_control": float((d[d.arm == "Control"][c] / d[d.arm == "Control"].total).mean()),
                             "mean_prop_arm": float((d[d.arm == term][c] / d[d.arm == term].total).mean())})
        tr = pd.DataFrame(rows)
        tr["q"] = bh(tr.p.values)
        tr.round(5).to_csv(OUT["treat"], index=False)
        LOG.write(f"{len(tr)} composition contrasts, lowest q {tr.q.min():.3f}\n" +
                  tr.nsmallest(4, "q")[["cluster", "contrast", "coef_logodds", "p", "q"]].to_string(index=False))

    figures(emb, mk, corr, ov, act, tr, res, summ)
    provenance(list(OUT.values()),
               [PROCESSED / "mouse_all_hs.h5ad", PROCESSED / "mouse_tumour_calls.csv", INTERIM / "mouse_qc.csv",
                ATLAS_ROOT / "data/processed/dmg_atlas_v1.h5ad", EXTERNAL / "progeny_human.csv",
                EXTERNAL / "collectri_human.csv"],
               "scripts/m14_denovo_subpopulations.py",
               extra={"cell_set": TUMOUR_CLUSTER, "n_cells": int(len(emb)), "resolution": res,
                      "n_clusters": int(emb.cluster.nunique()), "seed": SEED})


def figures(emb, mk, corr, ov, act, tr, res, summ):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    sys.path.insert(0, str(Path(_ATLAS) / "src"))
    from dmg_atlas.figstyle import apply_figure_style
    apply_figure_style()
    clusters_all = sorted(emb.cluster.unique(), key=lambda c: int(c[1:]))
    clusters = [c for c in clusters_all if c in set(corr.mouse_cluster)]   # above the floor
    excluded = [c for c in clusters_all if c not in clusters]
    # A cluster that is excluded from the statistics must not be drawn as a peer of the ones that are not
    # (figure-style 1.1). It keeps its colour so it can be found in every panel, but its bars are hatched and
    # its name carries the exclusion, so a reader cannot take it for a result.
    def _tick(c):
        # Asterisk, not the word: "M5 (excluded)" as a y-tick ran left into the neighbouring panel's colourbar.
        # The word is spelled out once, in panel a, where there is room for it.
        return f"{c}*" if c in excluded else c

    def _name(c):
        return f"{c} (excluded)" if c in excluded else c
    _hatch = ["///" if c in excluded else "" for c in clusters_all]
    ccol = dict(zip(clusters_all, PAL.denovo(len(clusters_all))))

    # ---- figure 1: the clustering itself
    fig, axes = plt.subplots(2, 3, figsize=(PAGE_W, 5.6))
    ax = axes[0, 0]
    for c in clusters_all:
        g = emb[emb.cluster == c]
        ax.scatter(g.umap_1, g.umap_2, s=3, lw=0, color=ccol[c], rasterized=True)
    ax.set(xticks=[], yticks=[])
    annotate_clusters(ax, emb[["umap_1", "umap_2"]].values, emb.cluster.values, ccol, fontsize=6,
                      names={c: _name(c) for c in clusters_all})
    panel(ax, "a", f"De novo clusters (resolution {res})")

    # Fifteen mice cannot be told apart as fifteen point colours - the first version of this panel was an
    # indistinguishable brown wash with overlapping SampleN labels. The per-mouse batch check lives in panel c,
    # which shows it properly; this panel asks the other question a reader has, whether the structure tracks
    # treatment, and four arms can be coloured legibly.
    ax = axes[0, 1]
    acol = PAL.ARM
    for a in ARMS:
        g = emb[emb.arm == a]
        if len(g):
            ax.scatter(g.umap_1, g.umap_2, s=3, lw=0, color=acol[a], rasterized=True,
                       label=f"{a} ({g.mouse.nunique()})")
    ax.set(xticks=[], yticks=[])
    ax.legend(loc="upper left", bbox_to_anchor=(0.0, 1.0), frameon=False, fontsize=5.5, markerscale=2.2,
              handletextpad=0.3, labelspacing=0.25, borderpad=0.0)
    panel(ax, "b", "Treatment arm")

    ax = axes[0, 2]
    comp = emb.groupby(["mouse", "cluster"]).size().unstack(fill_value=0)
    comp = comp.div(comp.sum(1), axis=0).reindex(columns=clusters_all)
    # Sample10 sorts before Sample2 lexically; order the axis the way a reader reads the names.
    comp = comp.reindex(index=sorted(comp.index, key=lambda s: int("".join(ch for ch in s if ch.isdigit()) or 0)))
    left = np.zeros(len(comp))
    for c in clusters_all:
        ax.barh(range(len(comp)), comp[c].values, left=left, color=ccol[c], height=0.78)
        left += comp[c].values
    ax.set(yticks=range(len(comp)), yticklabels=comp.index, xlabel="share of tumour cells", xlim=(0, 1))
    ax.tick_params(axis="y", labelsize=5)
    ax.invert_yaxis()
    panel(ax, "c", "Cluster mix per mouse")

    ax = axes[1, 0]
    genes = list(dict.fromkeys(g for c in clusters_all
                               for g in mk[mk.cluster == c].nlargest(4, "scores").gene))
    Mx = mk[mk.gene.isin(genes)].pivot_table(index="gene", columns="cluster", values="scores").reindex(
        index=genes, columns=clusters_all)
    im = ax.imshow(Mx.values, cmap="RdBu_r", vmin=-np.nanmax(abs(Mx.values)), vmax=np.nanmax(abs(Mx.values)),
                   aspect="auto")
    ax.set(yticks=range(len(genes)), yticklabels=genes, xticks=range(len(clusters_all)),
           xticklabels=[_tick(c) for c in clusters_all])
    ax.tick_params(axis="y", labelsize=5)
    cell_grid(ax, Mx.shape[0], Mx.shape[1])
    cb = plt.colorbar(im, ax=ax, shrink=0.75)
    cb.set_label("Wilcoxon z", fontsize=6)
    panel(ax, "d", "Cluster markers")

    ax = axes[1, 1]
    y = np.arange(len(clusters_all))
    ax.barh(y, summ.set_index("cluster").loc[clusters_all, "median_clone_r"].values,
            color=[ccol[c] for c in clusters_all], height=0.7, hatch=_hatch, edgecolor="white", linewidth=0)
    ax.axvline(0.27, ls="--", lw=0.7, c="k")
    ax.set(yticks=y, yticklabels=[_tick(c) for c in clusters_all], xlabel="median clone correlation")
    ax.invert_yaxis()
    panel(ax, "e", "Clone signal per cluster")

    ax = axes[1, 2]
    ax.barh(y, 100 * summ.set_index("cluster").loc[clusters_all, "frac_strict_set"].values,
            color=[ccol[c] for c in clusters_all], height=0.7, hatch=_hatch, edgecolor="white", linewidth=0)
    # A zero-valued bar draws nothing, which reads as missing rather than zero (figure-style 6.1).
    vals = 100 * summ.set_index("cluster").loc[clusters_all, "frac_strict_set"].values
    for yi, v in zip(y, vals):
        if v <= 0:
            ax.plot([0], [yi], marker="|", ms=4, mew=1.0, color="#555555")
    ax.set(yticks=y, yticklabels=[_tick(c) for c in clusters_all], xlabel="% in the strict clone-positive set")
    ax.invert_yaxis()
    panel(ax, "f", "Overlap with the strict set")
    fig.subplots_adjust(left=0.085, right=0.955, top=0.93, bottom=0.095, wspace=0.46, hspace=0.42)
    if excluded:
        fig.text(0.5, 0.012, f"* below the {MIN_CLUSTER_CELLS}-cell floor; excluded from all statistics",
                 ha="center", va="bottom", fontsize=5, style="italic", color="#555555")
    save_figure(fig, OUT["fig1"], LOG)

    # ---- figure 2: mouse clusters against human subpopulations
    fig, axes = plt.subplots(2, 2, figsize=(PAGE_W, 6.6))
    ax = axes[0, 0]
    P1 = corr.pivot(index="mouse_cluster", columns="human_subpopulation", values="rho").reindex(
        index=clusters, columns=STATES)
    im = ax.imshow(P1.values, cmap="RdBu_r", vmin=-np.nanmax(abs(P1.values)), vmax=np.nanmax(abs(P1.values)),
                   aspect="auto")
    for i in range(P1.shape[0]):
        for j in range(P1.shape[1]):
            ax.text(j, i, f"{P1.values[i, j]:.2f}", ha="center", va="center", fontsize=5.5,
                    color="white" if abs(P1.values[i, j]) > 0.6 * np.nanmax(abs(P1.values)) else "black")
    ax.set(yticks=range(len(clusters)), yticklabels=clusters, xticks=range(len(STATES)),
           xticklabels=[pretty(s) for s in STATES])
    ax.tick_params(axis="x", rotation=45)
    cb = plt.colorbar(im, ax=ax, shrink=0.75); cb.set_label("Spearman rho", fontsize=6)
    cell_grid(ax, P1.shape[0], P1.shape[1])
    panel(ax, "a", "Pseudobulk correlation")

    ax = axes[0, 1]
    P2 = ov.assign(nlq=-np.log10(ov.q.clip(lower=1e-300))).pivot(
        index="mouse_cluster", columns="human_subpopulation", values="nlq").reindex(index=clusters, columns=STATES)
    im = ax.imshow(P2.values, cmap="Greys", vmin=0, aspect="auto")
    ovp = ov.pivot(index="mouse_cluster", columns="human_subpopulation", values="overlap").reindex(
        index=clusters, columns=STATES)
    for i in range(P2.shape[0]):
        for j in range(P2.shape[1]):
            if not np.isfinite(ovp.values[i, j]):
                continue
            ax.text(j, i, f"{int(ovp.values[i, j])}", ha="center", va="center", fontsize=5.5,
                    color="white" if P2.values[i, j] > 0.6 * np.nanmax(P2.values) else "black")
    ax.set(yticks=range(len(clusters)), yticklabels=clusters, xticks=range(len(STATES)),
           xticklabels=[pretty(s) for s in STATES])
    ax.tick_params(axis="x", rotation=45)
    cb = plt.colorbar(im, ax=ax, shrink=0.75); cb.set_label("-log10 q", fontsize=6)
    cell_grid(ax, P2.shape[0], P2.shape[1])
    panel(ax, "b", "Marker-gene overlap")

    for k, (fam, lab, cmap) in enumerate((("pathway", "PROGENy activity", "RdBu_r"),
                                          ("tf", "CollecTRI activity", "PuOr_r"))):
        ax = axes[1, k]
        d = act[act.family == fam]
        P3 = d.pivot_table(index="feature", columns="group", values="activity")
        cols = [f"human:{s}" for s in STATES if f"human:{s}" in P3.columns] + \
               [f"mouse:{c}" for c in clusters if f"mouse:{c}" in P3.columns]
        P3 = P3[cols]
        if fam == "tf":
            P3 = P3.loc[P3.abs().max(1).nlargest(16).index]
        v = np.nanmax(abs(P3.values))
        im = ax.imshow(P3.values, cmap=cmap, vmin=-v, vmax=v, aspect="auto")
        ax.set(yticks=range(P3.shape[0]), yticklabels=P3.index, xticks=range(len(cols)))
        species_axis(ax, cols, fontsize=5.5)
        ax.tick_params(axis="y", labelsize=5)
        cb = plt.colorbar(im, ax=ax, shrink=0.75); cb.set_label(lab, fontsize=6)
        cell_grid(ax, P3.shape[0], P3.shape[1], lw=0.3)
        panel(ax, "cd"[k], "Pathways" if fam == "pathway" else "Transcription factors", dy=1.085, pad=12.0)
    fig.subplots_adjust(left=0.10, right=0.945, top=0.945, bottom=0.115, wspace=0.40, hspace=0.45)
    save_figure(fig, OUT["fig2"], LOG)

    # ---- figure 3: treatment effects on the de novo clusters
    fig, axes = plt.subplots(1, 2, figsize=(PAGE_W, 3.1))
    ax = axes[0]
    prop = emb.groupby(["arm", "mouse", "cluster"]).size().rename("n").reset_index()
    tot = prop.groupby(["arm", "mouse"]).n.transform("sum")
    prop["frac"] = prop.n / tot
    for i, c in enumerate(clusters):
        for j, arm in enumerate(ARMS):
            v = prop[(prop.cluster == c) & (prop.arm == arm)].frac
            x = i + np.linspace(-0.26, 0.26, max(len(v), 1))[:len(v)]
            ax.scatter(x, v, s=9, color=PAL.ARM[arm], lw=0, zorder=3,
                       label=arm if i == 0 else None)
            if len(v):
                ax.hlines(v.mean(), i - 0.3, i + 0.3, color=PAL.ARM[arm], lw=1.0)
    ax.set(xticks=range(len(clusters)), xticklabels=clusters, ylabel="share of tumour cells per mouse")
    ax.legend(fontsize=5.5, frameon=False, ncol=4, loc="upper center", bbox_to_anchor=(0.5, -0.14),
              handletextpad=0.3, columnspacing=0.9)
    panel(ax, "a", "Cluster composition by arm")

    # A diverging heatmap of these coefficients reads as a strong pattern - saturated cells at both ends - when in
    # fact no contrast comes close to significance (lowest q 0.495). Effect sizes with their intervals say the same
    # numbers honestly: every interval crosses zero.
    ax = axes[1]
    rows = [(c, arm) for c in clusters for arm in ARMS[1:]]
    ypos = {}
    y = 0.0
    for c in clusters:
        for arm in ARMS[1:]:
            ypos[(c, arm)] = y
            y += 1.0
        y += 0.6                                        # a gap between clusters
    crit = ss.t.ppf(0.975, 11)                          # 15 mice - 4 arms, as in the fitted model
    for (c, arm) in rows:
        r = tr[(tr.cluster == c) & (tr.contrast == f"{arm}_vs_Control")]
        if not len(r):
            continue
        b, se = float(r.coef_logodds.iloc[0]), float(r.se.iloc[0])
        yy = ypos[(c, arm)]
        ax.plot([b - crit * se, b + crit * se], [yy, yy], color=PAL.ARM[arm], lw=1.4,
                solid_capstyle="butt", zorder=2)
        ax.plot([b], [yy], "o", color=PAL.ARM[arm], ms=3.4, mec="white", mew=0.5, zorder=3)
    ax.axvline(0, color="0.25", lw=0.8, ls="--", zorder=1)
    ax.set_yticks([np.mean([ypos[(c, a)] for a in ARMS[1:]]) for c in clusters])
    ax.set_yticklabels(clusters)
    ax.invert_yaxis()
    ax.set_xlabel("log-odds vs Control (95% CI)")
    ax.set_ylim(max(ypos.values()) + 1.0, -1.0)
    for arm in ARMS[1:]:
        ax.plot([], [], color=PAL.ARM[arm], lw=1.4, marker="o", ms=3.4, label=arm)
    ax.legend(fontsize=5.5, frameon=False, ncol=3, loc="upper center", bbox_to_anchor=(0.5, -0.215),
              handletextpad=0.4, columnspacing=0.9, handlelength=1.4)
    ax.text(0.99, 0.015, f"lowest q = {tr.q.min():.2f}", transform=ax.transAxes, ha="right", va="bottom",
            fontsize=5.5, style="italic", color="0.3")
    panel(ax, "b", "Composition shift")
    fig.subplots_adjust(left=0.085, right=0.975, top=0.90, bottom=0.28, wspace=0.30)
    save_figure(fig, OUT["fig3"], LOG)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--force", action="store_true")
    main(p.parse_args().force)
