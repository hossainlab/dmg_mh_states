#!/usr/bin/env python
"""m09 — does the mouse model contain the four human DMG tumour subpopulations? Markers, pathways, TFs, correlation.

Four independent lines of evidence per subpopulation (OPC-like, AC-like, OC-like, MES-like, plus cycling), each
computed for the human atlas and for the mouse on the same footing:
  markers      per marker gene of each human subpopulation: mean log-CPM and detection rate in the human cells of that
               subpopulation and in the mouse tumour cells, so presence is judged gene by gene, not only by a score
  pathways     PROGENy (14 pathways) on pseudobulk profiles
  TF activity  CollecTRI (1,185 TFs, 42,990 edges) on the same pseudobulk profiles
  correlation  per-subpopulation signature correlation, full matrix, every mouse tumour cluster included
Both species are scored with the HUMAN resources applied to the ortholog-mapped matrices, so the regulon and pathway
weights are identical across species; the mouse side is restricted to the 12,759 one-to-one orthologs, which lowers
per-TF target coverage and is reported as n_targets_mouse.

Inputs : data/processed/mouse_all_hs.h5ad, data/processed/mouse_tumour_calls.csv, data/interim/mouse_qc.csv
         ../data/processed/dmg_atlas_v1.h5ad, data/human/malignant_state_assignment.csv.gz
         data/external/collectri_human.csv, data/external/progeny_human.csv (omnipathdb.org, see resources_provenance.json)
Outputs: results/tables/m09_marker_genes.csv, m09_pathway_activity.csv, m09_tf_activity.csv,
         m09_subpopulation_correlation.csv, results/figures/m09_subpopulation_evidence.png
Run: nice -n 10 python scripts/m09_subpopulation_evidence.py [--force]
"""
import sys, argparse
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "src"), str(ROOT.parent / "src")]
import dmg_atlas.utils as U  # noqa: E402
import numpy as np, pandas as pd, anndata as ad, scipy.sparse as sp  # noqa: E402
import scipy.stats as ss  # noqa: E402
from dmg_atlas.states import STATE_GENE_SETS  # noqa: E402
from mouse_cmp import INTERIM, PROCESSED, RESULTS, LOGS, HUMAN, EXTERNAL, ATLAS_ROOT  # noqa: E402
from mouse_cmp.common import provenance, skip_existing  # noqa: E402

LOG = U.StepLogger(LOGS / "m09_subpopulation_evidence.log", context="m09_subpopulation_evidence")
STATES = ["OPC_like", "AC_like", "OC_like", "MES_like", "cycling"]
OUT = {"markers": RESULTS / "tables/m09_marker_genes.csv",
       "paths":   RESULTS / "tables/m09_pathway_activity.csv",
       "tfs":     RESULTS / "tables/m09_tf_activity.csv",
       "corr":    RESULTS / "tables/m09_subpopulation_correlation.csv",
       "fig":     RESULTS / "figures/m09_subpopulation_evidence.png"}


def cpm_log(X):
    tot = np.asarray(X.sum(1)).ravel()
    tot[tot == 0] = 1
    return np.log1p(sp.diags(1e6 / tot) @ X if sp.issparse(X) else (X * (1e6 / tot[:, None])))


def pseudobulk(X, groups, genes):
    """summed counts per group -> log CPM, returned as groups x genes."""
    out, names = [], []
    for g in pd.unique(groups):
        m = groups == g
        if m.sum() < 10:
            continue
        s = np.asarray(X[m].sum(0)).ravel()
        out.append(np.log1p(1e6 * s / max(s.sum(), 1))); names.append(g)
    return pd.DataFrame(out, index=names, columns=genes)


def main(force):
    if skip_existing(OUT.values(), force, LOG):
        return
    with LOG.step("load mouse (human symbols) + groups"):
        M = ad.read_h5ad(PROCESSED / "mouse_all_hs.h5ad")
        calls = pd.read_csv(PROCESSED / "mouse_tumour_calls.csv", index_col="cell_id")
        qc = pd.read_csv(INTERIM / "mouse_qc.csv", index_col=0)
        prim = ((calls["tumour_call"].reindex(M.obs_names) == "tumour") &
                qc["qc_pass"].reindex(M.obs_names).fillna(False)).values
        lab = M.obs["author_label"].astype(str).values
        mgroup = np.where(prim, "tumour", np.where(np.isin(lab, ["Microglia", "Mono/Macro"]), "myeloid_normal",
                          np.where(np.isin(lab, ["T cells", "NK cells"]), "T_NK_normal",
                          np.where(lab == "Astrocyte-like tumor cells", "astrocyte_like_normal",
                          np.where(np.isin(lab, ["Neuronal-like tumor", "Purkinje-like tumor"]), "cerebellar_neuron", "other")))))
        Xm = M.layers["counts"] if "counts" in M.layers else M.X
        LOG.write(f"mouse {M.shape}; primary tumour {int(prim.sum()):,}; groups {pd.Series(mgroup).value_counts().to_dict()}")

    with LOG.step("load human malignant cells by subpopulation"):
        H = ad.read_h5ad(ATLAS_ROOT / "data/processed/dmg_atlas_v1.h5ad")
        H.var_names = H.var["gene_symbol"].astype(str).values
        H.var_names_make_unique()
        st = pd.read_csv(HUMAN / "malignant_state_assignment.csv.gz", index_col=0)["state_raw"]
        hs = st.reindex(H.obs_names)
        keep = hs.notna().values
        Hm = H[keep]
        hstate = hs[keep].astype(str).values
        hsample = Hm.obs["sample_id"].astype(str).values
        LOG.write(f"human malignant {Hm.shape}; states {pd.Series(hstate).value_counts().to_dict()}")

        shared = sorted(set(M.var_names) & set(Hm.var_names))
        LOG.write(f"shared ortholog genes for the comparison: {len(shared):,}")

    with LOG.step("marker gene table"):
        rows = []
        Mi = {g: i for i, g in enumerate(M.var_names)}
        Hi = {g: i for i, g in enumerate(Hm.var_names)}
        Lm = cpm_log(Xm[prim])
        for s in STATES:
            for g in STATE_GENE_SETS[s]:
                r = {"subpopulation": s, "gene": g, "in_mouse_orthologs": g in Mi, "in_human_atlas": g in Hi}
                if g in Mi:
                    v = np.asarray(Lm[:, Mi[g]].todense()).ravel() if sp.issparse(Lm) else Lm[:, Mi[g]]
                    r["mouse_tumour_mean_logcpm"] = round(float(v.mean()), 4)
                    r["mouse_tumour_detect_frac"] = round(float((v > 0).mean()), 4)
                if g in Hi:
                    m = hstate == s
                    col = Hm[:, Hi[g]].X
                    col = np.asarray(col.todense()).ravel() if sp.issparse(col) else np.asarray(col).ravel()
                    tot = np.asarray(Hm.X.sum(1)).ravel(); tot[tot == 0] = 1
                    v = np.log1p(1e6 * col / tot)
                    r["human_subpop_mean_logcpm"] = round(float(v[m].mean()), 4)
                    r["human_subpop_detect_frac"] = round(float((v[m] > 0).mean()), 4)
                rows.append(r)
        mk = pd.DataFrame(rows)
        mk.to_csv(OUT["markers"], index=False)
        summ = (mk.groupby("subpopulation")
                .agg(genes=("gene", "size"), with_ortholog=("in_mouse_orthologs", "sum"),
                     mouse_detect_median=("mouse_tumour_detect_frac", "median"),
                     human_detect_median=("human_subpop_detect_frac", "median")).round(3))
        LOG.write("marker detection:\n" + summ.to_string())

    with LOG.step("pseudobulk + PROGENy + CollecTRI"):
        import decoupler as dc
        hb = pseudobulk(Hm[:, shared].X, pd.Series([f"human:{a}" for a in hstate]).values, shared)
        hb_s = pseudobulk(Hm[:, shared].X, pd.Series([f"human:{a}|{b}" for a, b in zip(hstate, hsample)]).values, shared)
        mb = pseudobulk(Xm[:, [M.var_names.get_loc(g) for g in shared]], np.array([f"mouse:{g}" for g in mgroup]), shared)
        mb_s = pseudobulk(Xm[:, [M.var_names.get_loc(g) for g in shared]],
                          np.array([f"mouse:{g}|{s}" for g, s in zip(mgroup, M.obs['sample'].astype(str))]), shared)
        prof = pd.concat([hb, mb])
        prof_s = pd.concat([hb_s, mb_s])
        net_p = pd.read_csv(EXTERNAL / "progeny_human.csv")
        net_t = pd.read_csv(EXTERNAL / "collectri_human.csv")[["source", "target", "weight"]]

        def score(df, net, tmin):
            a = ad.AnnData(df.astype(np.float32))
            dc.mt.ulm(a, net, tmin=tmin)
            sc_ = a.obsm["score_ulm"]; pv = a.obsm["padj_ulm"]
            return sc_, pv

        ps, pp = score(prof, net_p, 5)
        ts, tp = score(prof, net_t, 10)
        ps_s, _ = score(prof_s, net_p, 5)
        pa = (ps.stack().rename("activity").reset_index().rename(columns={"level_0": "group", "level_1": "pathway"})
              .merge(pp.stack().rename("padj").reset_index().rename(columns={"level_0": "group", "level_1": "pathway"}),
                     on=["group", "pathway"]))
        pa["species"] = pa.group.str.split(":").str[0]
        pa.round(5).to_csv(OUT["paths"], index=False)
        ta = (ts.stack().rename("activity").reset_index().rename(columns={"level_0": "group", "level_1": "tf"})
              .merge(tp.stack().rename("padj").reset_index().rename(columns={"level_0": "group", "level_1": "tf"}),
                     on=["group", "tf"]))
        ta["species"] = ta.group.str.split(":").str[0]
        cov = net_t[net_t.target.isin(shared)].groupby("source").size().rename("n_targets_shared")
        ta = ta.merge(cov, left_on="tf", right_index=True, how="left")
        ta.round(5).to_csv(OUT["tfs"], index=False)
        LOG.write(f"pathway rows {len(pa):,}; TF rows {len(ta):,}; profiles {list(prof.index)}")

    with LOG.step("per-subpopulation cross-species correlation"):
        hz = hb.loc[[f"human:{s}" for s in STATES if f"human:{s}" in hb.index]]
        var = hz.var(0).sort_values(ascending=False)
        genes = var.head(2000).index
        A = hz[genes].apply(ss.zscore, axis=0)
        B = mb[genes].apply(ss.zscore, axis=0)
        rows = []
        for mi in B.index:
            for hi in A.index:
                rho, p = ss.spearmanr(B.loc[mi], A.loc[hi])
                rows.append({"mouse_group": mi, "human_subpopulation": hi, "rho": round(float(rho), 4), "p": p,
                             "n_genes": len(genes)})
        corr = pd.DataFrame(rows)
        corr.to_csv(OUT["corr"], index=False)
        LOG.write("correlation (mouse tumour row):\n" +
                  corr[corr.mouse_group == "mouse:tumour"][["human_subpopulation", "rho"]].to_string(index=False))

    with LOG.step("figure"):
        figure(mk, pa, ta, corr)
    provenance(list(OUT.values()),
               [PROCESSED / "mouse_all_hs.h5ad", PROCESSED / "mouse_tumour_calls.csv",
                ATLAS_ROOT / "data/processed/dmg_atlas_v1.h5ad", EXTERNAL / "collectri_human.csv",
                EXTERNAL / "progeny_human.csv"],
               "m09_subpopulation_evidence.py",
               extra={"resources": "human CollecTRI + PROGENy applied to both species on the ortholog-mapped matrices",
                      "shared_genes": len(shared)})
    LOG.write(f"done; peak RSS {U.peak_rss_gb():.2f} GB")


def figure(mk, pa, ta, corr):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    sys.path.insert(0, str(ATLAS_ROOT / "src"))
    from dmg_atlas.figstyle import apply_figure_style, letter
    apply_figure_style()
    fig, axes = plt.subplots(2, 2, figsize=(12.5, 10))

    ax = axes[0, 0]
    d = mk.dropna(subset=["mouse_tumour_detect_frac", "human_subpop_detect_frac"])
    for s, g in d.groupby("subpopulation"):
        ax.scatter(g.human_subpop_detect_frac, g.mouse_tumour_detect_frac, s=14, label=f"{s} ({len(g)})", alpha=0.85)
    ax.plot([0, 1], [0, 1], ls=":", c="0.5", lw=0.8)
    ax.set(xlabel="detected in human cells of that subpopulation", ylabel="detected in mouse tumour cells",
           xlim=(-0.02, 1.02), ylim=(-0.02, 1.02))
    ax.legend(fontsize=6, frameon=False, loc="lower right")
    letter(ax, "a", title="Marker genes, gene by gene")

    ax = axes[0, 1]
    piv = pa.pivot_table(index="pathway", columns="group", values="activity")
    order = [c for c in piv.columns if c.startswith("human:")] + [c for c in piv.columns if c.startswith("mouse:")]
    piv = piv[order]
    v = np.nanmax(np.abs(piv.values))
    im = ax.imshow(piv.values, cmap="RdBu_r", vmin=-v, vmax=v, aspect="auto")
    ax.set(yticks=range(piv.shape[0]), yticklabels=piv.index, xticks=range(piv.shape[1]))
    ax.set_xticklabels([c.replace("human:", "H ").replace("mouse:", "M ") for c in piv.columns], rotation=90, fontsize=6)
    ax.axvline(sum(c.startswith("human:") for c in piv.columns) - 0.5, color="k", lw=1)
    plt.colorbar(im, ax=ax, shrink=0.7, label="PROGENy activity")
    letter(ax, "b", title="Pathway activity, human and mouse")

    ax = axes[1, 0]
    sel = ta[(ta.species == "human") & ta.group.isin([f"human:{s}" for s in STATES])]
    top = (sel.assign(a=sel.activity.abs()).sort_values("a", ascending=False)
           .drop_duplicates("tf").head(18).tf.tolist())
    piv2 = ta[ta.tf.isin(top)].pivot_table(index="tf", columns="group", values="activity")
    piv2 = piv2[[c for c in order if c in piv2.columns]]
    v = np.nanmax(np.abs(piv2.values))
    im = ax.imshow(piv2.values, cmap="PuOr_r", vmin=-v, vmax=v, aspect="auto")
    ax.set(yticks=range(piv2.shape[0]), yticklabels=piv2.index, xticks=range(piv2.shape[1]))
    ax.set_xticklabels([c.replace("human:", "H ").replace("mouse:", "M ") for c in piv2.columns], rotation=90, fontsize=6)
    ax.axvline(sum(c.startswith("human:") for c in piv2.columns) - 0.5, color="k", lw=1)
    plt.colorbar(im, ax=ax, shrink=0.7, label="CollecTRI activity")
    letter(ax, "c", title="Transcription-factor activity")

    ax = axes[1, 1]
    piv3 = corr.pivot(index="mouse_group", columns="human_subpopulation", values="rho")
    piv3.columns = [c.replace("human:", "") for c in piv3.columns]
    piv3.index = [i.replace("mouse:", "") for i in piv3.index]
    im = ax.imshow(piv3.values, cmap="RdBu_r", vmin=-0.6, vmax=0.6, aspect="auto")
    ax.set(xticks=range(piv3.shape[1]), yticks=range(piv3.shape[0]), yticklabels=piv3.index)
    ax.set_xticklabels(piv3.columns, rotation=45, ha="right")
    for i in range(piv3.shape[0]):
        for j in range(piv3.shape[1]):
            ax.text(j, i, f"{piv3.values[i, j]:.2f}", ha="center", va="center", fontsize=6,
                    color="white" if abs(piv3.values[i, j]) > 0.35 else "black")
    plt.colorbar(im, ax=ax, shrink=0.7, label="Spearman rho")
    letter(ax, "d", title="Pseudobulk correlation per subpopulation")

    fig.tight_layout()
    fig.savefig(OUT["fig"], dpi=300)
    plt.close(fig)


if __name__ == "__main__":
    p = argparse.ArgumentParser(); p.add_argument("--force", action="store_true")
    main(p.parse_args().force)
