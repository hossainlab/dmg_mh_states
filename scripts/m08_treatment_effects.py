#!/usr/bin/env python
"""m08 — treatment effects (FUS, RT, FUS+RT vs untreated Control) on tumour-cell states / programmes and on the TME.

Inputs : data/interim/mouse_all.h5ad (m01), data/interim/mouse_qc.csv (m02), data/external/mouse_human_orthologs.csv (m03),
         data/processed/mouse_tumour_calls.csv (m04), data/processed/m05_cnmf_usages.parquet + results/tables/m05_cnmf_programs.csv
         + results/tables/m05_state_scores.csv.gz (m05), data/processed/mouse_tme_labels.csv (m07),
         data/human/metaprogram_genes.csv + metaprogram_families.csv, ../results/tables/genesets/phase6_curated_gene_sets.csv
Outputs (never overwritten unless --force):
  results/tables/m08_per_mouse.csv           long: mouse x set x category -> k, n, proportion / mean score (every modelled value)
  results/tables/m08_composition.csv         quasi-binomial GLM per set x category x contrast: log-odds coef, SE, p, BH q
                                             (primary model ~ arm; sensitivity ~ arm + block), n mice per arm
  results/tables/m08_scores.csv              OLS per score (tumour programme usage, human MP / state scores, pseudobulk gene-set
                                             scores) x contrast, mouse-level, inverse-variance not used (equal weight per mouse)
  results/tables/m08_pseudobulk_de.csv.gz    limma-style moderated t per population x contrast x gene (logFC, AveExpr, t, p, q)
  results/tables/m08_de_summary.csv          DE genes (q < 0.1, |logFC| > 0.5) per population x contrast, top up / down genes
  results/figures/m08_treatment_effects.png
Run: nice -n 10 python scripts/m08_treatment_effects.py [--force]

Design: 15 mice (Control 4, FUS 4, RT 3, FUS+RT 4; no ICI arm and no Sample3 in the deposit). The mouse is the unit throughout.
  Contrasts: FUS, RT, FUS+RT each vs Control, and FUS+RT vs RT (does FUS add to RT). block = inferred processing set (Samples
  1-4, 5-8, 9-12, 13-16; common.sample_block) enters only the sensitivity model. A mouse enters a set if its denominator /
  cell count >= MIN_N (20); with few tumour cells per mouse the tumour analyses retain 11 mice (RT and FUS+RT 2 each) and are
  exploratory by construction.
Sets
  composition model: binomial GLM with quasi-binomial dispersion phi = max(1, Pearson X2 / df), t on df_resid.
  compartment   : all cells of the mouse, categories = tumour (m04 call), myeloid, T_NK, B_plasma, vascular, neuron, other
                  (sampling-sensitive: the immune share varies 5-77 % between mice and is not explained by arm or block)
  immune_type   : qc_pass normal immune cells by author label (Microglia, Mono/Macro, T, NK, B/Plasma, Neutrophils, pDC)
  myeloid_label / T_NK_label : m07 subcluster labels (human marker rule) within compartment
  tumour_state  : primary tumour cells (m05 definition A), human state raw argmax
  tumour_program: primary tumour cells, argmax of normalised pooled cNMF usage (m05 programmes p1-p9)
Scores (per mouse mean over its cells): normalised cNMF usage per programme, human state scores, UCell of human lineage MPs
  (primary tumour cells); pseudobulk gene-set scores per population = mean over set genes of the gene-wise z-scored log-CPM.
Pseudobulk DE: populations tumour (primary), Microglia, Mono/Macro, T cells, NK cells (qc_pass, normal call for non-tumour);
  sums of raw counts per mouse (mice with >= MIN_N cells); log-CPM = log2((count + 0.5) / (lib + 1) * 1e6); genes with CPM >= 1
  in >= 3 mice; OLS on ~ arm; empirical-Bayes variance moderation (Smyth 2004 moment estimator of d0, s0^2 with a lowess trend
  on AveExpr, i.e. limma-trend); moderated t on residual df + d0; BH within population x contrast.
"""
import sys, argparse
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "src"), str(ROOT.parent / "src")]
import dmg_atlas.utils as U  # noqa: E402
import numpy as np, pandas as pd, scipy.sparse as sp, anndata as ad  # noqa: E402
import scipy.stats as ss, scipy.special as spc  # noqa: E402
import statsmodels.api as sm  # noqa: E402
from statsmodels.nonparametric.smoothers_lowess import lowess  # noqa: E402
from mouse_cmp import INTERIM, PROCESSED, EXTERNAL, RESULTS, LOGS, HUMAN, ATLAS_ROOT  # noqa: E402
from mouse_cmp.common import provenance, skip_existing, bh, sample_block, ARMS  # noqa: E402
from dmg_atlas import metaprograms as MP  # noqa: E402
from dmg_atlas.figstyle import apply_figure_style  # noqa: E402
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt  # noqa: E402

ap = argparse.ArgumentParser(); ap.add_argument("--force", action="store_true"); args = ap.parse_args()
IN, QC, ORTH = INTERIM / "mouse_all.h5ad", INTERIM / "mouse_qc.csv", EXTERNAL / "mouse_human_orthologs.csv"
CALLS, USE, PROGS, STS, TME = (PROCESSED / "mouse_tumour_calls.csv", PROCESSED / "m05_cnmf_usages.parquet", RESULTS / "tables/m05_cnmf_programs.csv",
                               RESULTS / "tables/m05_state_scores.csv.gz", PROCESSED / "mouse_tme_labels.csv")
MPG, FAM, GS = HUMAN / "metaprogram_genes.csv", HUMAN / "metaprogram_families.csv", ATLAS_ROOT / "results/tables/genesets/phase6_curated_gene_sets.csv"
T, FIG = RESULTS / "tables", RESULTS / "figures"
OUTS = [T / "m08_per_mouse.csv", T / "m08_composition.csv", T / "m08_scores.csv", T / "m08_pseudobulk_de.csv.gz", T / "m08_de_summary.csv",
        FIG / "m08_treatment_effects.png"]
log = U.StepLogger(LOGS / "m08_treatment_effects.log", "m08")
if skip_existing(OUTS, args.force, log):
    sys.exit(0)
MIN_N, Q_DE, LFC_DE = 20, 0.1, 0.5
CONTRASTS = {"FUS_vs_Control": {"FUS": 1}, "RT_vs_Control": {"RT": 1}, "FUS+RT_vs_Control": {"FUS+RT": 1},
             "FUS+RT_vs_RT": {"FUS+RT": 1, "RT": -1}}
GSETS = ["Interferon", "Hypoxia", "Glycolysis", "Stress_HSP", "IEG_dissociation", "OXPHOS_nuclear", "Cholesterol_biosynthesis", "Cycling_markers"]

with log.step("load"):
    a = ad.read_h5ad(IN)
    o = a.obs.copy()
    qc = pd.read_csv(QC, index_col=0).reindex(a.obs_names); calls = pd.read_csv(CALLS, index_col=0).reindex(a.obs_names)
    tme = pd.read_csv(TME, index_col=0)
    o["qc_pass"] = qc.qc_pass.values; o["tumour"] = (calls.tumour_call == "tumour").values
    o["primary"] = o.tumour & o.qc_pass
    mice = (o.groupby("sample", observed=True).agg(arm=("arm", "first")).assign(block=lambda d: sample_block(pd.Series(d.index.astype(str), index=d.index)).values))
    mice.index = mice.index.astype(str); mice["arm"] = pd.Categorical(mice.arm.astype(str), categories=ARMS)
    o["sample"] = o["sample"].astype(str)
    C = sp.csr_matrix(a.layers["counts"]); genes = a.var_names
    orth = pd.read_csv(ORTH); h2m = dict(zip(orth.human_symbol.dropna(), orth.mouse_symbol[orth.human_symbol.notna()]))


def design(m: pd.DataFrame, block: bool):
    X = pd.get_dummies(m.arm, drop_first=True).astype(float)
    if block:
        X = X.join(pd.get_dummies(m.block, drop_first=True, prefix="block").astype(float))
    return sm.add_constant(X, has_constant="add")


def cvec(X: pd.DataFrame, con: dict):
    v = np.zeros(X.shape[1]); cols = list(X.columns)
    for k, w in con.items():
        v[cols.index(k)] = w
    return v


per_mouse, comp_rows, score_rows = [], [], []

with log.step("composition sets"):
    lab = o.author_label.astype(str)
    CMP = {"Microglia": "myeloid", "Mono/Macro": "myeloid", "Neutrophils": "myeloid", "pDC": "myeloid", "T cells": "T_NK", "NK cells": "T_NK",
           "B/Plasma": "B_plasma", "Vascular cells": "vascular", "Meningeal fibroblasts (VLMC)": "vascular", "Neuronal-like tumor": "neuron",
           "Purkinje-like tumor": "neuron"}
    comp = lab.map(CMP).fillna("other"); comp[o.tumour] = "tumour"
    imm = o.qc_pass & ~o.tumour & lab.isin(["Microglia", "Mono/Macro", "T cells", "NK cells", "B/Plasma", "Neutrophils", "pDC"])
    st = pd.read_csv(STS, index_col=0, usecols=["cell_id", "state_raw", "definition"])
    st = st[st.definition == "primary"].state_raw
    U_ = pd.read_parquet(USE); U_ = U_.div(U_.sum(1), axis=0); U_.columns = [c.split("|")[-1] for c in U_.columns]
    sets = {"compartment": comp, "immune_type": lab[imm],
            "myeloid_label": tme.label[tme.compartment == "myeloid"].str.replace("|low_confidence", "", regex=False),
            "T_NK_label": tme.label[tme.compartment == "T_NK"].str.replace("|low_confidence", "", regex=False),
            "tumour_state": st, "tumour_program": U_.idxmax(1)}
    for sname, lab_s in sets.items():
        ct = pd.crosstab(o.loc[lab_s.index, "sample"], lab_s)
        n = ct.sum(1); ct = ct[n >= MIN_N]; n = n[n >= MIN_N]
        m = mice.loc[ct.index]
        for cat in ct.columns:
            for i in ct.index:
                per_mouse.append(dict(set=sname, sample=i, arm=m.arm[i], category=cat, k=int(ct.loc[i, cat]), n=int(n[i]), value=ct.loc[i, cat] / n[i]))
            y = np.column_stack([ct[cat].values, (n - ct[cat]).values])
            for blk in (False, True):
                X = design(m, blk)
                if X.shape[1] >= len(m) or m.arm.nunique() < 2:
                    continue
                try:
                    fit = sm.GLM(y, X.values, family=sm.families.Binomial()).fit()
                except Exception:
                    continue
                # quasi-binomial: phi = Pearson X2 / df (statsmodels' scale='X2' is wrong for two-column endog: 0.24 vs 79.8 here)
                phi = max(1.0, fit.pearson_chi2 / fit.df_resid) if fit.df_resid > 0 else np.nan
                for cn, con in CONTRASTS.items():
                    if not all(k in X.columns for k in con):
                        continue
                    v = cvec(X, con); est = float(v @ fit.params); se = float(np.sqrt(phi * (v @ fit.cov_params() @ v)))
                    df_r = len(m) - X.shape[1]
                    comp_rows.append(dict(set=sname, category=cat, contrast=cn, model="arm+block" if blk else "arm", coef_logodds=est, se=se, phi=phi,
                                          p=2 * ss.t.sf(abs(est / se), df_r) if se > 0 else np.nan, df_resid=df_r,
                                          n_mice=len(m), **{f"n_{a_}": int((m.arm == a_).sum()) for a_ in ARMS},
                                          mean_prop_control=float((ct[cat] / n)[m.arm == "Control"].mean()),
                                          mean_prop_arm=float((ct[cat] / n)[m.arm == next(iter(con))].mean())))
    COMPR = pd.DataFrame(comp_rows)
    COMPR["q"] = np.nan
    for key, d in COMPR.groupby(["set", "contrast", "model"]):
        COMPR.loc[d.index, "q"] = bh(d.p)
    log.write("composition hits (q < 0.1, model arm):\n" + COMPR[(COMPR.q < 0.1) & (COMPR.model == "arm")].round(4).to_string(index=False))


def ols_scores(name, family, df):
    """df: index mouse, columns scores -> OLS per column per contrast (both models)."""
    m = mice.loc[df.index]
    for col in df.columns:
        yv = df[col].astype(float)
        for blk in (False, True):
            X = design(m, blk)
            if X.shape[1] >= len(m):
                continue
            fit = sm.OLS(yv.values, X.values).fit()
            for cn, con in CONTRASTS.items():
                if not all(k in X.columns for k in con):
                    continue
                v = cvec(X, con); t_ = fit.t_test(v)
                score_rows.append(dict(family=family, population=name, score=col, contrast=cn, model="arm+block" if blk else "arm",
                                       estimate=float(t_.effect.ravel()[0]), se=float(t_.sd.ravel()[0]), p=float(t_.pvalue), n_mice=len(m),
                                       sd_control=float(yv[m.arm == "Control"].std()) if (m.arm == "Control").sum() > 1 else np.nan))
        for i in df.index:
            per_mouse.append(dict(set=f"score:{family}:{name}", sample=i, arm=mice.arm[i], category=col, k=np.nan, n=np.nan, value=float(df.loc[i, col])))


with log.step("tumour scores"):
    prim = o.index[o.primary]
    ms = o.loc[prim, "sample"]; okm = ms.value_counts()[lambda x: x >= MIN_N].index
    keepc = prim[ms.isin(okm).values]
    ols_scores("tumour", "cnmf_usage", U_.loc[U_.index.intersection(keepc)].groupby(o.loc[U_.index.intersection(keepc), "sample"]).mean())
    S = pd.read_csv(STS, index_col=0); S = S[S.definition == "primary"]
    scols = [c for c in S.columns if c.startswith("score_")]
    ols_scores("tumour", "human_state_score", S.loc[S.index.intersection(keepc), scols].groupby(o.loc[S.index.intersection(keepc), "sample"]).mean())
    mpg = pd.read_csv(MPG); fam = pd.read_csv(FAM).set_index("metaprogram")
    lineage = [m_ for m_ in fam.index if isinstance(fam.family.get(m_), str) and pd.isna(fam.excluded_reason.get(m_))]
    mps = {f"{m_}_{fam.family[m_]}": [h2m[g] for g in v if g in h2m] for m_, v in mpg.groupby("metaprogram").gene if m_ in lineage}
    idx = a.obs_names.get_indexer(keepc)
    UC = MP.ucell_scores(C[idx], genes.values, mps); UC.index = keepc
    ols_scores("tumour", "human_MP_ucell", UC.groupby(o.loc[keepc, "sample"]).mean())
    log.write(f"tumour scores: {len(okm)} mice ({mice.loc[okm].arm.value_counts().to_dict()})")

with log.step("pseudobulk DE + gene-set scores"):
    gs = pd.read_csv(GS); gs = gs[gs.gene_set.isin(GSETS)]
    gsm = {k: [h2m[g] for g in v if g in h2m] for k, v in gs.groupby("gene_set").gene}
    POPS = {"tumour": o.primary, **{l_: (o.qc_pass & ~o.tumour & (lab == l_)) for l_ in ["Microglia", "Mono/Macro", "T cells", "NK cells"]}}
    de_parts, summ = [], []
    for pop, msk in POPS.items():
        cells = np.flatnonzero(msk.values); sm_ = o["sample"].values[cells]
        vc = pd.Series(sm_).value_counts(); use = vc.index[vc >= MIN_N]
        if len(use) < 6:
            log.write(f"{pop}: {len(use)} mice with >= {MIN_N} cells, skipped"); continue
        PB = pd.DataFrame({s_: np.asarray(C[cells[sm_ == s_]].sum(0)).ravel() for s_ in use}, index=genes).T
        m = mice.loc[PB.index]
        lib = PB.sum(1)
        cpm = PB.div(lib, axis=0) * 1e6
        keepg = (cpm >= 1).sum(0) >= 3
        Y = np.log2((PB.loc[:, keepg] + 0.5).div(lib + 1, axis=0) * 1e6)
        # gene-set scores (z over mice of log-CPM, mean over set genes)
        Zs = (Y - Y.mean()) / (Y.std(ddof=0) + 1e-9)
        ols_scores(pop, "pseudobulk_geneset", pd.DataFrame({k: Zs[[g for g in v if g in Zs.columns]].mean(1) for k, v in gsm.items()}))
        X = design(m, False); Xv = X.values
        XtXi = np.linalg.pinv(Xv.T @ Xv); B = XtXi @ Xv.T @ Y.values
        R = Y.values - Xv @ B; dfr = Xv.shape[0] - Xv.shape[1]
        s2 = (R ** 2).sum(0) / dfr; amean = Y.values.mean(0)
        # limma-trend: lowess trend of the log-variance on AveExpr, then Smyth moment estimates of d0 / s0^2 around the trend
        zl = np.log(np.clip(s2, 1e-12, None)); e = zl - spc.digamma(dfr / 2) + np.log(dfr / 2)
        trend = lowess(e, amean, frac=0.3, return_sorted=False)
        evar = np.mean((e - trend) ** 2) * len(e) / (len(e) - 1) - spc.polygamma(1, dfr / 2)
        if evar > 0:
            yv = 0.5 + 1 / evar
            for _ in range(50):   # trigamma inverse (Newton, as limma)
                tri = spc.polygamma(1, yv); dif = tri * (1 - tri / evar) / spc.polygamma(2, yv); yv = yv + dif
                if np.all(-dif / yv < 1e-8):
                    break
            d0 = 2 * yv; s02 = np.exp(trend + spc.digamma(d0 / 2) - np.log(d0 / 2))
            s2p = (d0 * s02 + dfr * s2) / (d0 + dfr)
        else:
            d0 = np.inf; s2p = np.exp(trend)
        dft = dfr + (d0 if np.isfinite(d0) else 1e6)
        for cn, con in CONTRASTS.items():
            if not all(k in X.columns for k in con):
                continue
            v = cvec(X, con); lfc = v @ B; se = np.sqrt(s2p * (v @ XtXi @ v)); t_ = lfc / se
            p = 2 * ss.t.sf(np.abs(t_), dft)
            d = pd.DataFrame({"population": pop, "contrast": cn, "gene": Y.columns, "logFC": lfc, "AveExpr": amean, "t": t_, "p": p,
                              "q": bh(p), "df_total": dft, "d0": d0, "n_mice": len(m)})
            de_parts.append(d)
            sig = d[(d.q < Q_DE) & (d.logFC.abs() > LFC_DE)]
            summ.append(dict(population=pop, contrast=cn, n_mice=len(m), n_genes_tested=len(d), d0=float(d0), n_up=int((sig.logFC > 0).sum()),
                             n_down=int((sig.logFC < 0).sum()), top_up=";".join(sig[sig.logFC > 0].sort_values("t", ascending=False).gene.head(15)),
                             top_down=";".join(sig[sig.logFC < 0].sort_values("t").gene.head(15))))
        log.write(f"{pop}: {len(m)} mice {m.arm.value_counts().to_dict()}, {keepg.sum()} genes, d0 {d0:.1f}")
    DE = pd.concat(de_parts); SUMM = pd.DataFrame(summ)
    SC = pd.DataFrame(score_rows); SC["q"] = np.nan
    for key, d in SC.groupby(["family", "population", "contrast", "model"]):
        SC.loc[d.index, "q"] = bh(d.p)
    log.write("DE summary:\n" + SUMM.drop(columns=["top_up", "top_down"]).to_string(index=False))
    log.write("score hits (p < 0.05, model arm):\n" + SC[(SC.p < 0.05) & (SC.model == "arm")].round(4).to_string(index=False))

with log.step("write + figure"):
    PM = pd.DataFrame(per_mouse)
    PM.to_csv(OUTS[0], index=False); COMPR.to_csv(OUTS[1], index=False); SC.to_csv(OUTS[2], index=False); DE.to_csv(OUTS[3], index=False)
    SUMM.to_csv(OUTS[4], index=False)
    apply_figure_style()
    acol = {"Control": "0.45", "FUS": "#3b6fb6", "RT": "#d98b2b", "FUS+RT": "#b8433a"}
    panels = [("compartment", "tumour", "tumour cells / all cells"), ("compartment", "T_NK", "T/NK / all cells"),
              ("compartment", "myeloid", "myeloid / all cells"), ("immune_type", "Mono/Macro", "Mono/Macro / immune"),
              ("myeloid_label", "interferon_myeloid", "IFN-myeloid / myeloid"), ("myeloid_label", "homeostatic_microglia", "homeostatic MG / myeloid"),
              ("T_NK_label", "NK", "NK / T+NK"), ("tumour_state", "cycling", "cycling / tumour (state)")]
    sc_panels = [("human_MP_ucell", "tumour", "MP4_hypoxia", "tumour: MP4 hypoxia UCell"), ("human_MP_ucell", "tumour", "MP8_MES_like", "tumour: MP8 MES-like UCell"),
                 ("pseudobulk_geneset", "tumour", "Interferon", "tumour IFN set score"), ("pseudobulk_geneset", "Microglia", "Interferon", "microglia IFN set score")]
    fig, axes = plt.subplots(3, 4, figsize=(7.4, 7.0), gridspec_kw=dict(hspace=0.95, wspace=0.55)); axes = axes.ravel()
    rng = np.random.default_rng(0)

    def dots(ax, d, title, stat):
        for j, a_ in enumerate(ARMS):
            v = d[d.arm == a_].value
            ax.scatter(j + rng.uniform(-0.12, 0.12, len(v)), v, s=10, color=acol[a_], zorder=3)
            if len(v):
                ax.hlines(v.mean(), j - 0.25, j + 0.25, color="k", lw=0.8)
        ax.set_xticks(range(4), ARMS, rotation=45, ha="right"); ax.set_title(f"{title}\n{stat}", fontsize=5.5, loc="left")

    for i, (s_, c_, t_) in enumerate(panels):
        d = PM[(PM.set == s_) & (PM.category == c_)]
        r = COMPR[(COMPR.set == s_) & (COMPR.category == c_) & (COMPR.model == "arm")].set_index("contrast")
        stat = "p: " + " | ".join(f"{k.split('_')[0].replace('FUS+RT', 'F+R')} {r.p[k]:.2g}" for k in ["FUS_vs_Control", "RT_vs_Control", "FUS+RT_vs_Control"] if k in r.index)
        dots(axes[i], d, t_, stat)
    for i, (f_, pop, sc_, t_) in enumerate(sc_panels):
        d = PM[(PM.set == f"score:{f_}:{pop}") & (PM.category == sc_)]
        r = SC[(SC.family == f_) & (SC.population == pop) & (SC.score == sc_) & (SC.model == "arm")].set_index("contrast")
        stat = "p: " + " | ".join(f"{k.split('_')[0].replace('FUS+RT', 'F+R')} {r.p[k]:.2g}" for k in ["FUS_vs_Control", "RT_vs_Control", "FUS+RT_vs_Control"] if k in r.index)
        dots(axes[8 + i], d, t_, stat)
    fig.suptitle("Per-mouse values (dots), arm means (bars); p from mouse-level models vs Control (exploratory; 2-4 mice per arm)", fontsize=6.5, y=0.995)
    fig.savefig(OUTS[5], dpi=200, bbox_inches="tight"); plt.close(fig)
    provenance(OUTS, [IN, QC, ORTH, CALLS, USE, PROGS, STS, TME, MPG, FAM, GS], "scripts/m08_treatment_effects.py")
