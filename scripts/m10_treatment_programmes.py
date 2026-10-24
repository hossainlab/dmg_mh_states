#!/usr/bin/env python
"""m10 — how FUS, RT and FUS+RT shift tumour-cell pathways and TFs, and the immune programmes, with the genes.

Extends m08 (composition, state scores, pseudobulk DE) with the three things the collaborator asked for that m08 did
not cover: pathway activity, transcription-factor activity, and explicitly named immune programmes.

  tumour    per mouse, pseudobulk of its QC-pass tumour cells -> PROGENy pathway activity and CollecTRI TF activity;
            each tested arm-vs-control with mice as the unit (OLS on the per-mouse activity, BH within family)
  immune    four programmes scored per cell with AUCell inside the myeloid and T/NK compartments, averaged per mouse,
            tested the same way: inflammation, immune suppression, T-cell dysfunction/exhaustion, T-cell activation
            (plus a regulatory-T programme). Gene membership is listed in PROGRAMMES below and written to the output
            table, so every programme is auditable gene by gene; these are canonical markers, not a published set.
  genes     per-compartment pseudobulk differential expression per arm (per-mouse log CPM, OLS, BH over expressed
            genes), reported for the programme genes specifically and as the top movers overall

Mice: tumour analyses use the mice with >= 20 QC-pass tumour cells (as m08); immune analyses use every mouse with
>= 50 cells in the compartment. Human resources (CollecTRI, PROGENy) are applied to the ortholog-mapped matrix, as in m09.

Inputs : data/processed/mouse_all_hs.h5ad, mouse_tumour_calls.csv, data/interim/mouse_qc.csv,
         data/external/{collectri_human.csv,progeny_human.csv}
Outputs: results/tables/m10_tumour_pathway_tf.csv, m10_immune_programmes.csv, m10_immune_programme_tests.csv,
         m10_immune_de_genes.csv.gz, results/figures/m10_treatment_programmes.png
Run: nice -n 10 python scripts/m10_treatment_programmes.py [--force]
"""
import sys, argparse
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "src"), str(ROOT.parent / "src")]
import dmg_atlas.utils as U  # noqa: E402
import numpy as np, pandas as pd, anndata as ad, scipy.sparse as sp  # noqa: E402
import statsmodels.api as sm  # noqa: E402
from mouse_cmp import INTERIM, PROCESSED, RESULTS, LOGS, EXTERNAL, ATLAS_ROOT  # noqa: E402
from mouse_cmp.common import provenance, skip_existing, bh  # noqa: E402

LOG = U.StepLogger(LOGS / "m10_treatment_programmes.log", context="m10_treatment_programmes")
ARMS = ["FUS", "RT", "FUS+RT"]
MIN_TUMOUR, MIN_IMMUNE = 20, 50

# Mouse-nomenclature programmes, scored on the deposit's own matrix (data/interim/mouse_all.h5ad, 15,304 genes).
# The ortholog-mapped matrix loses many immune genes to the one-to-one filter, so programme scoring uses mouse symbols;
# the tumour pathway/TF analysis below still uses the ortholog matrix, because CollecTRI and PROGENy are human resources.
# Coverage is written into the output: several canonical genes (Ifng, Gzmb, Foxp3, Ctla4, Cd69, Tbx21) are ABSENT FROM
# THE DEPOSIT itself, so the activation and regulatory-T programmes are under-measured and flagged as such.
PROGRAMMES = {
 "inflammation": ["Il1a","Il1b","Il6","Tnf","Cxcl1","Cxcl2","Ccl2","Ccl3","Ccl4","Ptgs2","Nlrp3","Nfkb1","S100a8","S100a9"],
 "interferon": ["Isg15","Ifit1","Ifit3","Irf7","Stat1","Stat2","Mx1","Oas1a","Cxcl9","Cxcl10","B2m","Tap1","Psmb8","Psmb9","H2-K1","H2-D1"],
 "immune_suppression": ["Il10","Tgfb1","Arg1","Ido1","Cd274","Pdcd1lg2","Lgals9","Entpd1","Nt5e","Mrc1","Trem2","Apoe","Spp1","Vsir"],
 "T_cell_dysfunction": ["Pdcd1","Lag3","Havcr2","Tigit","Ctla4","Tox","Entpd1","Cd160","Btla","Eomes"],
 "T_cell_activation": ["Ifng","Gzma","Gzmb","Gzmk","Prf1","Il2","Cd69","Tbx21","Icos","Cd28","Il2ra","Klrg1","Nkg7","Ccl5"],
 "regulatory_T": ["Foxp3","Il2ra","Ctla4","Ikzf2","Tnfrsf18","Tnfrsf4"],
}
OUT = {"tum": RESULTS / "tables/m10_tumour_pathway_tf.csv",
       "cover": RESULTS / "tables/m10_programme_coverage.csv",
       "prog": RESULTS / "tables/m10_immune_programmes.csv",
       "tests": RESULTS / "tables/m10_immune_programme_tests.csv",
       "de": RESULTS / "tables/m10_immune_de_genes.csv.gz",
       "fig": RESULTS / "figures/m10_treatment_programmes.png"}


def per_mouse_pseudobulk(X, mouse, genes, min_cells):
    out, names = [], []
    for m in pd.unique(mouse):
        k = mouse == m
        if k.sum() < min_cells:
            continue
        s = np.asarray(X[k].sum(0)).ravel()
        out.append(np.log1p(1e6 * s / max(s.sum(), 1))); names.append(m)
    return pd.DataFrame(out, index=names, columns=genes)


def arm_tests(values: pd.DataFrame, arm_of: dict, family: str, unit: str):
    """OLS of each column on arm (reference = Control), mice as the unit; BH within family."""
    rows = []
    arms = pd.Series({i: arm_of[i] for i in values.index})
    X = pd.get_dummies(arms, drop_first=False).astype(float)
    for col in values.columns:
        y = values[col].astype(float)
        des = sm.add_constant(X[[a for a in ARMS if a in X.columns]], has_constant="add")
        try:
            fit = sm.OLS(y.values, des.values.astype(float)).fit()
        except Exception:
            continue
        for j, a in enumerate([a for a in ARMS if a in X.columns], start=1):
            rows.append({"family": family, "unit": unit, "feature": col, "contrast": f"{a}_vs_Control",
                         "estimate": fit.params[j], "se": fit.bse[j], "p": fit.pvalues[j],
                         "n_mice": int(len(y)), "mean_control": float(y[arms.eq("Control").reindex(y.index).fillna(False)].mean())})
    df = pd.DataFrame(rows)
    if len(df):
        df["q"] = bh(df.p.values)
    return df


def main(force):
    if skip_existing(OUT.values(), force, LOG):
        return
    import decoupler as dc
    with LOG.step("load"):
        M = ad.read_h5ad(PROCESSED / "mouse_all_hs.h5ad")
        calls = pd.read_csv(PROCESSED / "mouse_tumour_calls.csv", index_col="cell_id")
        qc = pd.read_csv(INTERIM / "mouse_qc.csv", index_col=0)
        X = M.layers["counts"] if "counts" in M.layers else M.X
        obs = M.obs
        prim = ((calls["tumour_call"].reindex(M.obs_names) == "tumour") &
                qc["qc_pass"].reindex(M.obs_names).fillna(False)).values
        qcp = qc["qc_pass"].reindex(M.obs_names).fillna(False).values
        lab = obs["author_label"].astype(str).values
        mouse = obs["sample"].astype(str).values
        arm_of = obs.groupby(obs["sample"].astype(str))["arm"].first().astype(str).to_dict()
        net_p = pd.read_csv(EXTERNAL / "progeny_human.csv")
        net_t = pd.read_csv(EXTERNAL / "collectri_human.csv")[["source", "target", "weight"]]
        LOG.write(f"{M.shape}; tumour {int(prim.sum()):,}; arms {pd.Series(arm_of).value_counts().to_dict()}")

    with LOG.step("tumour pathway + TF activity per mouse"):
        pb = per_mouse_pseudobulk(X[prim], mouse[prim], M.var_names.values, MIN_TUMOUR)
        LOG.write(f"tumour pseudobulk mice: {len(pb)} ({[f'{i}:{arm_of[i]}' for i in pb.index]})")
        a = ad.AnnData(pb.astype(np.float32)); dc.mt.ulm(a, net_p, tmin=5)
        path = a.obsm["score_ulm"]
        a2 = ad.AnnData(pb.astype(np.float32)); dc.mt.ulm(a2, net_t, tmin=10)
        tf = a2.obsm["score_ulm"]
        t1 = arm_tests(path, arm_of, "tumour_pathway_progeny", "mouse")
        t2 = arm_tests(tf[tf.abs().mean().sort_values(ascending=False).head(60).index], arm_of, "tumour_tf_collectri", "mouse")
        tum = pd.concat([t1, t2], ignore_index=True)
        tum.round(5).to_csv(OUT["tum"], index=False)
        sig = tum[tum.q < 0.1].sort_values("q")
        LOG.write(f"tumour tests {len(tum)}; q<0.1 {len(sig)}\n" +
                  (sig.head(12)[["family","feature","contrast","estimate","p","q"]].to_string(index=False) if len(sig) else "none"))

    with LOG.step("immune programme scoring (mouse symbols)"):
        Mm = ad.read_h5ad(INTERIM / "mouse_all.h5ad")
        Xmm = Mm.layers["counts"] if "counts" in Mm.layers else Mm.X
        assert list(Mm.obs_names) == list(M.obs_names), "cell order differs between the mouse-symbol and ortholog objects"
        comp = np.where(np.isin(lab, ["Microglia", "Mono/Macro", "pDC", "Neutrophils"]), "myeloid",
               np.where(np.isin(lab, ["T cells", "NK cells"]), "T_NK", "other"))
        net_prog = pd.DataFrame([{"source": k, "target": g, "weight": 1.0} for k, gs in PROGRAMMES.items() for g in gs])
        present = set(Mm.var_names)
        cover = pd.DataFrame([{"programme": k, "n_genes": len(v), "n_present": sum(g in present for g in v),
                               "missing_from_deposit": ";".join(g for g in v if g not in present)}
                              for k, v in PROGRAMMES.items()])
        cover.to_csv(RESULTS / "tables/m10_programme_coverage.csv", index=False)
        net_prog = net_prog[net_prog.target.isin(present)]
        LOG.write("programme coverage:\n" + cover.to_string(index=False))
        keep = qcp & np.isin(comp, ["myeloid", "T_NK"])
        sub = ad.AnnData(Xmm[keep].astype(np.float32), obs=pd.DataFrame(
            {"mouse": mouse[keep], "arm": [arm_of[m] for m in mouse[keep]], "compartment": comp[keep],
             "label": lab[keep]}, index=M.obs_names[keep]), var=pd.DataFrame(index=Mm.var_names))
        sub.X = sp.csr_matrix(sub.X)
        tot = np.asarray(sub.X.sum(1)).ravel(); tot[tot == 0] = 1
        sub.X = sp.diags(1e4 / tot) @ sub.X; sub.X.data = np.log1p(sub.X.data)
        dc.mt.aucell(sub, net_prog, tmin=3)
        sc_ = sub.obsm["score_aucell"]
        prog = sc_.join(sub.obs[["mouse", "arm", "compartment", "label"]])
        prog.to_csv(OUT["prog"])
        tests = []
        for cmp_ in ["myeloid", "T_NK"]:
            d = prog[prog.compartment == cmp_]
            per = d.groupby("mouse")[list(sc_.columns)].mean()
            per = per.loc[[i for i in per.index if (d.mouse == i).sum() >= MIN_IMMUNE]]
            tests.append(arm_tests(per, arm_of, f"immune_programme_{cmp_}", "mouse"))
        tst = pd.concat(tests, ignore_index=True)
        tst.round(5).to_csv(OUT["tests"], index=False)
        s = tst[tst.q < 0.1].sort_values("q")
        LOG.write(f"programme tests {len(tst)}; q<0.1 {len(s)}\n" +
                  (s[["family","feature","contrast","estimate","mean_control","p","q"]].to_string(index=False) if len(s) else "none"))

    with LOG.step("gene-level DE within myeloid and T/NK"):
        de = []
        for cmp_ in ["myeloid", "T_NK"]:
            k = qcp & (comp == cmp_)
            pbc = per_mouse_pseudobulk(Xmm[k], mouse[k], Mm.var_names.values, MIN_IMMUNE)
            expressed = pbc.columns[(pbc > 0).sum(0) >= max(3, int(0.5 * len(pbc)))]
            d = arm_tests(pbc[expressed], arm_of, f"de_{cmp_}", "mouse")
            d["programme"] = d.feature.map({g: k2 for k2, gs in PROGRAMMES.items() for g in gs})
            de.append(d)
        de = pd.concat(de, ignore_index=True)
        de.round(5).to_csv(OUT["de"], index=False)
        s = de[de.q < 0.1]
        LOG.write(f"DE tests {len(de):,}; q<0.1 {len(s)}\n" +
                  (s.sort_values('q').head(15)[["family","feature","contrast","estimate","p","q","programme"]].to_string(index=False) if len(s) else "none"))

    with LOG.step("figure"):
        figure(tum, prog, tst, de)
    provenance(list(OUT.values()),
               [PROCESSED / "mouse_all_hs.h5ad", PROCESSED / "mouse_tumour_calls.csv",
                EXTERNAL / "collectri_human.csv", EXTERNAL / "progeny_human.csv"],
               "m10_treatment_programmes.py",
               extra={"programmes": {k: len(v) for k, v in PROGRAMMES.items()}, "min_tumour_cells": MIN_TUMOUR,
                      "min_immune_cells": MIN_IMMUNE})
    LOG.write(f"done; peak RSS {U.peak_rss_gb():.2f} GB")


def figure(tum, prog, tst, de):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    sys.path.insert(0, str(ATLAS_ROOT / "src"))
    from dmg_atlas.figstyle import apply_figure_style, letter
    apply_figure_style()
    fig, axes = plt.subplots(2, 2, figsize=(12.5, 9.5))
    COL = {"Control": "#4a5568", "FUS": "#2b6cb0", "RT": "#de2d26", "FUS+RT": "#805ad5"}

    ax = axes[0, 0]
    d = tum[tum.family == "tumour_pathway_progeny"]
    piv = d.pivot_table(index="feature", columns="contrast", values="estimate")
    im = ax.imshow(piv.values, cmap="RdBu_r", vmin=-np.nanmax(abs(piv.values)), vmax=np.nanmax(abs(piv.values)), aspect="auto")
    ax.set(yticks=range(piv.shape[0]), yticklabels=piv.index, xticks=range(piv.shape[1]))
    ax.set_xticklabels([c.replace("_vs_Control", "") for c in piv.columns], rotation=0)
    qp = d.pivot_table(index="feature", columns="contrast", values="q")
    for i in range(piv.shape[0]):
        for j in range(piv.shape[1]):
            if qp.values[i, j] < 0.1:
                ax.text(j, i, "*", ha="center", va="center", fontsize=9)
    plt.colorbar(im, ax=ax, shrink=0.75, label="activity change vs control")
    letter(ax, "a", title="Tumour pathway activity by arm")

    ax = axes[0, 1]
    pm = prog[prog.compartment == "T_NK"].groupby(["mouse", "arm"])[["T_cell_dysfunction", "T_cell_activation"]].mean().reset_index()
    for arm, g in pm.groupby("arm"):
        ax.scatter(g.T_cell_activation, g.T_cell_dysfunction, s=42, color=COL.get(arm, "0.5"), label=f"{arm} ({len(g)})")
    ax.set(xlabel="T-cell activation (AUCell)", ylabel="T-cell dysfunction (AUCell)")
    ax.legend(fontsize=6, frameon=False)
    letter(ax, "b", title="T/NK compartment, one point per mouse")

    ax = axes[1, 0]
    d2 = tst.copy()
    d2["lab"] = d2.family.str.replace("immune_programme_", "") + ": " + d2.feature
    piv2 = d2.pivot_table(index="lab", columns="contrast", values="estimate")
    q2 = d2.pivot_table(index="lab", columns="contrast", values="q")
    im = ax.imshow(piv2.values, cmap="PuOr_r", vmin=-np.nanmax(abs(piv2.values)), vmax=np.nanmax(abs(piv2.values)), aspect="auto")
    ax.set(yticks=range(piv2.shape[0]), yticklabels=piv2.index, xticks=range(piv2.shape[1]))
    ax.set_xticklabels([c.replace("_vs_Control", "") for c in piv2.columns], rotation=0)
    for i in range(piv2.shape[0]):
        for j in range(piv2.shape[1]):
            if q2.values[i, j] < 0.1:
                ax.text(j, i, "*", ha="center", va="center", fontsize=9)
    plt.colorbar(im, ax=ax, shrink=0.75, label="AUCell change vs control")
    letter(ax, "c", title="Immune programmes by arm")

    ax = axes[1, 1]
    # the two best-ranked programme tests of the 36: myeloid immune suppression (p 0.026) and T/NK T-cell
    # activation (p 0.101). regulatory_T ranks 4th but is unmeasurable here (3/6 genes, no Foxp3), so it is not shown.
    panels = [("myeloid", "immune_suppression"), ("T_NK", "T_cell_activation")]
    order = ["Control", "FUS", "RT", "FUS+RT"]
    for i, (cmp_, feat) in enumerate(panels):
        d = prog[prog.compartment == cmp_].groupby(["mouse", "arm"])[[feat]].mean().reset_index()
        for j, arm in enumerate(order):
            v = d.loc[d.arm == arm, feat]
            x = i * 5 + j + np.linspace(-0.18, 0.18, len(v))
            ax.scatter(x, v, s=34, color=COL[arm], label=arm if i == 0 else None)
            ax.hlines(v.mean(), i * 5 + j - 0.3, i * 5 + j + 0.3, color=COL[arm], lw=1.6)
    ax.set(xticks=[0, 1, 2, 3, 5, 6, 7, 8], xticklabels=order * 2, ylabel="AUCell score (mean per mouse)")
    ax.tick_params(axis="x", rotation=60, labelsize=6)
    ax.text(1.5, ax.get_ylim()[1] * 0.995, "myeloid: immune suppression (15 mice, p 0.026)", ha="center", fontsize=6)
    ax.text(6.5, ax.get_ylim()[1] * 0.995, "T/NK: T-cell activation (9 mice, p 0.101)", ha="center", fontsize=6)
    ax.legend(fontsize=6, frameon=False, loc="lower left")
    letter(ax, "d", title="The two best-ranked programme tests, per mouse")

    fig.tight_layout(); fig.savefig(OUT["fig"], dpi=300); plt.close(fig)


if __name__ == "__main__":
    p = argparse.ArgumentParser(); p.add_argument("--force", action="store_true")
    main(p.parse_args().force)
