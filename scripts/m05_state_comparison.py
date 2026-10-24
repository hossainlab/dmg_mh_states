#!/usr/bin/env python
"""m05 — do mouse tumour cells recapitulate the human DMG malignant states and metaprograms?

Inputs : data/processed/mouse_all_hs.h5ad (m03; one-to-one orthologs, human atlas symbols, raw counts),
         data/processed/mouse_tumour_calls.csv (m04), data/interim/mouse_qc.csv (m02),
         data/human/: metaprogram_genes.csv, metaprogram_families.csv, metaprogram_family_composition_assayz.csv,
         malignant_state_assignment.csv.gz, c2l_signatures_mp.csv, reference_gene_sets_scripted.csv;
         dmg_atlas.states.STATE_GENE_SETS (the human state sets, as used by ../scripts/18_composition.py)
Outputs (never overwritten unless --force):
  results/tables/m05_state_scores.csv.gz        per mouse cell of every tumour definition: 5 state scores, raw-argmax state,
                                                within-mouse-z state (mice with >= MIN_MOUSE cells), set memberships
  results/tables/m05_state_proportions.csv      state proportions per tumour definition (mouse) and per assay class (human)
  results/tables/m05_mp_ucell_by_group.csv      per lineage MP: median UCell in mouse tumour vs each mouse normal group
  results/tables/m05_mp_family_proportions.csv  MP-family assignment proportions, mouse tumour vs human (assay-class z rule)
  results/tables/m05_xspecies_correlation.csv   Spearman rho, mouse cluster pseudobulk vs human class signature (gene z-scores)
  results/tables/m05_cnmf_programs.csv          mouse de novo cNMF programmes (top-50 genes, human symbols) per run
  results/tables/m05_cnmf_ksel.csv              k selection per run (dmg_atlas.metaprograms rule)
  results/tables/m05_cnmf_vs_human.csv          every mouse programme x human MP / Liu 2022 set: overlap, Jaccard, hypergeom q
  data/processed/m05_cnmf_usages.parquet        per cell usages of the pooled primary run
  results/figures/m05_state_comparison.png
Run: nice -n 10 python scripts/m05_state_comparison.py [--force]

Tumour definitions (A is primary; the others show how the answer depends on the definition)
  A primary     : tumour_call == 'tumour' (m04) AND qc_pass (m02)
  B all_tumour  : tumour_call == 'tumour'
  C author_glial: author label OPC-like or astrocyte-like tumour (no CNV evidence used)
  D author_all  : every author 'tumor' label (includes the granule- and Purkinje-like clusters m04 calls normal)
Methods
  States: log-normalised (1e4, log1p) human-symbol counts; scanpy score_genes per state set over the cells of the definition
    (control-gene pool = those cells), seed 0; 'state_raw' = argmax (unassigned if max <= 0); 'state' = argmax of within-mouse
    z-scores (the human rule), only for mice with >= MIN_MOUSE cells. Human proportions from malignant_state_assignment.csv.gz.
    Within-mouse z forces every mouse to spread over the states, so composition is read from state_raw; state is reported for
    parity with the human table only.
  MPs: UCell (dmg_atlas.metaprograms.ucell_scores, max_rank 1500) of every human MP top-50 set (orthologs present) on raw counts of
    all mouse cells; families from metaprogram_families.csv (technical / contamination MPs excluded, as in the human atlas);
    family score = mean of member-MP z-scores with z computed over the primary tumour cells pooled (the analogue of the human
    within-assay-class rule); argmax, 'unassigned' if <= 0.
  Cross-species correlation: pseudobulk = summed counts per group -> CP10k -> log1p; human c2l class signatures (per-class
    mean expression) -> CP10k -> log1p; genes = the 2,000 one-to-one genes most variable across the human classes; each gene
    z-scored across groups within species; Spearman over genes. Mouse normal labels are positive controls (microglia vs human
    myeloid, T/NK vs T_NK, vascular vs vascular, normal astrocytes vs astrocyte, granule/Purkinje vs neuron).
  De novo programmes: dmg_atlas.metaprograms.run_cnmf_sample (2,000 within-run HVGs, k 4-9, 10 restarts, the human k rule) on
    raw ortholog counts. Runs: pooled primary; pooled primary without Sample6 (54 % of tumour cells); each mouse with >= 100
    primary cells. Overlap of top-50 genes with each human MP top-50 and each Liu 2022 Table S2a set: hypergeometric test,
    universe = the one-to-one genes, BH over all pairs of a run; 'best' = smallest q.
"""
import sys, argparse
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "src"), str(ROOT.parent / "src")]
import dmg_atlas.utils as U  # noqa: E402
import numpy as np, pandas as pd, scipy.sparse as sp, anndata as ad, scanpy as sc  # noqa: E402
import scipy.stats as ss  # noqa: E402
from mouse_cmp import INTERIM, PROCESSED, RESULTS, LOGS, HUMAN  # noqa: E402
from mouse_cmp.common import provenance, skip_existing, bh  # noqa: E402
from dmg_atlas.states import STATE_GENE_SETS  # noqa: E402
from dmg_atlas import metaprograms as MP  # noqa: E402
from dmg_atlas.figstyle import apply_figure_style  # noqa: E402
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt  # noqa: E402

ap = argparse.ArgumentParser(); ap.add_argument("--force", action="store_true"); args = ap.parse_args()
IN, CALLS, QC = PROCESSED / "mouse_all_hs.h5ad", PROCESSED / "mouse_tumour_calls.csv", INTERIM / "mouse_qc.csv"
H = {k: HUMAN / v for k, v in dict(mpg="metaprogram_genes.csv", fam="metaprogram_families.csv", famcomp="metaprogram_family_composition_assayz.csv",
                                    states="malignant_state_assignment.csv.gz", c2l="c2l_signatures_mp.csv", ref="reference_gene_sets_scripted.csv").items()}
T, FIG = RESULTS / "tables", RESULTS / "figures"
OUTS = [T / "m05_state_scores.csv.gz", T / "m05_state_proportions.csv", T / "m05_mp_ucell_by_group.csv", T / "m05_mp_family_proportions.csv",
        T / "m05_xspecies_correlation.csv", T / "m05_cnmf_programs.csv", T / "m05_cnmf_ksel.csv", T / "m05_cnmf_vs_human.csv",
        PROCESSED / "m05_cnmf_usages.parquet", FIG / "m05_state_comparison.png"]
log = U.StepLogger(LOGS / "m05_state_comparison.log", "m05")
if skip_existing(OUTS, args.force, log):
    sys.exit(0)
SEED, MIN_MOUSE, MIN_CNMF, N_VAR = 0, 20, 100, 2000
STATES = list(STATE_GENE_SETS)

with log.step("load + tumour definitions"):
    a = ad.read_h5ad(IN)
    calls = pd.read_csv(CALLS, index_col=0).reindex(a.obs_names)
    qc = pd.read_csv(QC, index_col=0).reindex(a.obs_names)
    o = a.obs
    lab = o.author_label.astype(str)
    DEF = {"primary": ((calls.tumour_call == "tumour") & qc.qc_pass).values,
           "all_tumour": (calls.tumour_call == "tumour").values,
           "author_glial": lab.isin(["OPC-like tumor cells", "Astrocyte-like tumor cells"]).values,
           "author_all": o.author_tumour.values.astype(bool)}
    log.write({k: int(v.sum()) for k, v in DEF.items()}.__repr__())
    counts = sp.csr_matrix(a.layers["counts"])
    a.X = counts.copy(); sc.pp.normalize_total(a, target_sum=1e4); sc.pp.log1p(a)
    # mouse groups for MP / correlation comparisons (normal astrocytes = author astrocyte-like cells called normal by m04)
    grp = lab.map({"Neuronal-like tumor": "granule-like neuron", "Purkinje-like tumor": "Purkinje-like neuron",
                   "Astrocyte-like tumor cells": "astrocyte (normal)", "OPC-like tumor cells": "OPC-like (normal call)"}).fillna(lab)
    grp[DEF["all_tumour"]] = "tumour:" + calls["sub"][DEF["all_tumour"]].where(calls.tumour_basis[DEF["all_tumour"]] != "clone_only", "clone_only")
    o["group"] = grp.values

with log.step("A: human state scores"):
    parts, props = [], []
    for name, m in DEF.items():
        b = a[m].copy()
        for st, genes in STATE_GENE_SETS.items():
            sc.tl.score_genes(b, [g for g in genes if g in b.var_names], score_name=f"score_{st}", random_state=SEED)
        S = b.obs[[f"score_{s}" for s in STATES]].copy()
        raw = S.idxmax(axis=1).str.replace("score_", "", regex=False); raw[S.max(axis=1) <= 0] = "unassigned"
        ms = b.obs["sample"].astype(str)
        okm = ms.map(ms.value_counts()) >= MIN_MOUSE
        Z = S[okm.values].groupby(ms[okm].values).transform(lambda x: (x - x.mean()) / (x.std(ddof=0) if x.std(ddof=0) > 0 else 1.0))
        st = Z.idxmax(axis=1).str.replace("score_", "", regex=False); st[Z.max(axis=1) <= 0] = "unassigned"
        S["state_raw"] = raw; S["state"] = st.reindex(S.index); S["definition"] = name; S["sample"] = ms.values; S["arm"] = b.obs.arm.astype(str).values
        parts.append(S)
        for col in ("state_raw", "state"):
            vc = S[col].value_counts(normalize=True)
            props += [dict(species="mouse", set=name, rule=col, state=k, proportion=v, n=int(S[col].notna().sum())) for k, v in vc.items()]
        # per-mouse composition summary (equal weight per mouse) for state_raw
        pm = pd.crosstab(S["sample"], S.state_raw, normalize="index")[lambda d: d.index.isin(ms.value_counts()[lambda x: x >= MIN_MOUSE].index)]
        props += [dict(species="mouse", set=name, rule="state_raw_mouse_mean", state=k, proportion=v, n=len(pm)) for k, v in pm.mean().items()]
    SC = pd.concat(parts)
    hs = pd.read_csv(H["states"], index_col=0, usecols=["Unnamed: 0", "sample_id", "assay_class", "state_raw", "state"])
    for ac, d in [("all", hs)] + list(hs.groupby("assay_class")):
        for col in ("state_raw", "state"):
            props += [dict(species="human", set=f"human_{ac}", rule=col, state=k, proportion=v, n=len(d)) for k, v in d[col].value_counts(normalize=True).items()]
        pm = pd.crosstab(d.sample_id, d.state_raw, normalize="index")
        props += [dict(species="human", set=f"human_{ac}", rule="state_raw_mouse_mean", state=k, proportion=v, n=len(pm)) for k, v in pm.mean().items()]
    PROPS = pd.DataFrame(props)
    log.write(PROPS[PROPS.rule == "state_raw"].pivot_table(index="set", columns="state", values="proportion").round(3).to_string())

with log.step("B: metaprogram UCell + families"):
    mpg = pd.read_csv(H["mpg"]); fam = pd.read_csv(H["fam"]).set_index("metaprogram")
    sets = {k: [g for g in v if g in a.var_names] for k, v in mpg.groupby("metaprogram").gene.agg(list).items()}
    UC = MP.ucell_scores(counts, a.var_names.values, sets); UC.index = a.obs_names
    lineage = [m for m in fam.index if isinstance(fam.family.get(m), str) and pd.isna(fam.excluded_reason.get(m))]
    prim = DEF["primary"]
    norm_groups = ["Microglia", "Mono/Macro", "T cells", "NK cells", "Vascular cells", "Meningeal fibroblasts (VLMC)", "astrocyte (normal)",
                   "granule-like neuron", "Purkinje-like neuron", "OPC-like (normal call)"]
    rows = []
    for m_ in sets:
        r = dict(metaprogram=m_, family=fam.family.get(m_), label=fam.label.get(m_), lineage=m_ in lineage, n_genes_mouse=len(sets[m_]),
                 median_tumour=float(UC.loc[prim, m_].median()))
        for g_ in norm_groups:
            sel = (o.group == g_).values & ~DEF["all_tumour"]
            r[f"median_{g_}"] = float(UC.loc[sel, m_].median())
        normals = [r[f"median_{g_}"] for g_ in norm_groups]
        r["tumour_over_max_normal"] = r["median_tumour"] / max(max(normals), 1e-6)
        r["top_group"] = (["tumour"] + norm_groups)[int(np.argmax([r["median_tumour"]] + normals))]
        rows.append(r)
    MPU = pd.DataFrame(rows)
    Zt = UC.loc[prim, lineage]; Zt = (Zt - Zt.mean()) / Zt.std(ddof=0)
    fams = fam.loc[lineage, "family"]
    F = pd.DataFrame({f: Zt[fams.index[fams == f]].mean(1) for f in sorted(fams.unique())})
    fa = F.idxmax(1); fa[F.max(1) <= 0] = "unassigned"
    famrows = [dict(species="mouse", set="primary", family=k, proportion=v, n=int(prim.sum())) for k, v in fa.value_counts(normalize=True).items()]
    hc = pd.read_csv(H["famcomp"])
    fcols = [c for c in hc.columns if c not in ("sample_id", "study", "assay", "assay_class", "n_cells")]
    for ac, d in [("all", hc)] + list(hc.groupby("assay_class")):
        w = d.n_cells / d.n_cells.sum()
        famrows += [dict(species="human", set=f"human_{ac}", family=c, proportion=float((d[c] * w).sum()), n=int(d.n_cells.sum())) for c in fcols]
    FAMP = pd.DataFrame(famrows)
    log.write(FAMP.pivot_table(index="set", columns="family", values="proportion").round(3).to_string())

with log.step("C: cross-species pseudobulk correlation"):
    sig = pd.read_csv(H["c2l"], index_col=0)
    sig = sig[sig.index.isin(a.var_names)]
    hl = np.log1p(sig / sig.sum() * 1e4)
    genes = hl.var(1).sort_values(ascending=False).index[:N_VAR]
    groups = sorted(set(o.group) - {"Neutrophils", "pDC", "B/Plasma", "Choroid plexus epithelial cells"})
    groups = [g_ for g_ in groups if (o.group == g_).sum() >= 30]
    gi = a.var_names.get_indexer(genes)
    PB = pd.DataFrame({g_: np.asarray(counts[(o.group == g_).values].sum(0)).ravel() for g_ in groups}, index=a.var_names)
    ml = np.log1p(PB / PB.sum() * 1e4).iloc[gi]
    zh = hl.loc[genes].sub(hl.loc[genes].mean(1), axis=0).div(hl.loc[genes].std(1) + 1e-9, axis=0)
    zm = ml.sub(ml.mean(1), axis=0).div(ml.std(1) + 1e-9, axis=0)
    XC = pd.DataFrame([dict(mouse_group=g_, human_class=c, n_mouse_cells=int((o.group == g_).sum()), rho=ss.spearmanr(zm[g_], zh[c])[0])
                       for g_ in groups for c in zh.columns])
    log.write(XC.pivot(index="mouse_group", columns="human_class", values="rho").round(2).to_string())

with log.step("D: de novo cNMF + overlap with human programmes"):
    runs = {"pooled_primary": prim, "pooled_primary_noSample6": prim & (o["sample"].astype(str) != "Sample6").values}
    for s_, n_ in o["sample"][prim].astype(str).value_counts().items():
        if n_ >= MIN_CNMF:
            runs[f"mouse_{s_}"] = prim & (o["sample"].astype(str) == s_).values
    ref = pd.read_csv(H["ref"]); ref = ref[ref.gene_set.str.startswith("Liu2022_TableS2a")]
    hsets = {**{f"{k}|{fam.label.get(k, '')}": set(v) for k, v in mpg.groupby("metaprogram").gene},
             **{k: set(v) for k, v in ref.groupby("gene_set").gene}}
    universe = set(a.var_names); N = len(universe)
    progs, ksels, ovl, usages = [], [], [], None
    for rn, m in runs.items():
        try:
            P_, W_, K_ = MP.run_cnmf_sample(counts[m], a.var_names.values, a.obs_names[m].values, rn)
        except ValueError as e:   # seurat_v3 loess can be singular on ~100-cell runs; recorded, run skipped
            ksels.append(pd.DataFrame([dict(sample_id=rn, run=rn, n_cells=int(m.sum()), rule=f"failed: {str(e)[:80]}")]))
            log.write(f"{rn}: {int(m.sum())} cells, cNMF failed ({e})"); continue
        progs.append(P_.assign(run=rn, n_cells=int(m.sum()))); ksels.append(K_.assign(run=rn, n_cells=int(m.sum())))
        if rn == "pooled_primary":
            usages = W_
        rows = []
        for p, d in P_.groupby("program"):
            top = set(d.gene)
            for hk, hv in hsets.items():
                hv_u = hv & universe; k_ = len(top & hv_u)
                rows.append(dict(run=rn, program=p, human_set=hk, n_top=len(top), n_human_in_universe=len(hv_u), overlap=k_,
                                 jaccard=k_ / len(top | hv_u), p=ss.hypergeom.sf(k_ - 1, N, len(hv_u), len(top)),
                                 shared=";".join(sorted(top & hv_u))))
        R = pd.DataFrame(rows); R["q"] = bh(R.p); ovl.append(R)
        log.write(f"{rn}: {int(m.sum())} cells, k={P_.k.iloc[0]}")
    PROG = pd.concat(progs); KSEL = pd.concat(ksels); OVL = pd.concat(ovl)
    best = OVL.sort_values("q").groupby(["run", "program"]).head(1)
    tops = PROG.groupby(["run", "program"]).gene.agg(lambda x: ";".join(list(x)[:12]))
    log.write("best human match per mouse programme:\n" + best.set_index(["run", "program"]).join(tops)[["human_set", "overlap", "jaccard", "q", "gene"]]
              .round(4).to_string())

with log.step("write + figure"):
    for d in (T, FIG):
        d.mkdir(parents=True, exist_ok=True)
    memb = pd.DataFrame(DEF, index=a.obs_names)
    SC.join(memb).to_csv(OUTS[0]); PROPS.to_csv(OUTS[1], index=False); MPU.to_csv(OUTS[2], index=False); FAMP.to_csv(OUTS[3], index=False)
    XC.to_csv(OUTS[4], index=False); PROG.to_csv(OUTS[5], index=False); KSEL.to_csv(OUTS[6], index=False); OVL.to_csv(OUTS[7], index=False)
    usages.to_parquet(OUTS[8])

    apply_figure_style()
    fig = plt.figure(figsize=(7.4, 9.6)); gs = fig.add_gridspec(3, 2, height_ratios=[1.25, 1, 1.1], hspace=0.6, wspace=0.55)
    ax = fig.add_subplot(gs[0, :])
    hc_order = ["malignant_OPC_like", "malignant_OC_like", "malignant_AC_like", "malignant_MES_like", "malignant_cycling", "malignant_hypoxia",
                "malignant_novel_MP15", "astrocyte", "oligodendrocyte", "neuron", "myeloid", "T_NK", "vascular"]
    M_ = XC.pivot(index="mouse_group", columns="human_class", values="rho")[hc_order]
    M_ = M_.loc[sorted(M_.index, key=lambda x: (not x.startswith("tumour"), x))]
    im = ax.imshow(M_.values, cmap="RdBu_r", vmin=-0.6, vmax=0.6, aspect="auto")
    ax.set_xticks(range(M_.shape[1]), [c.replace("malignant_", "mal. ") for c in M_.columns], rotation=45, ha="right")
    ax.set_yticks(range(M_.shape[0]), [f"{i} (n={int((o.group == i).sum())})" for i in M_.index])
    for i in range(M_.shape[0]):
        for j in range(M_.shape[1]):
            ax.text(j, i, f"{M_.values[i, j]:.2f}", ha="center", va="center", fontsize=5)
    fig.colorbar(im, ax=ax, shrink=0.8, label="Spearman rho"); ax.set_title("a  mouse cluster vs human class (z-scored pseudobulk, 2,000 genes)", loc="left")
    ax = fig.add_subplot(gs[1, 0])
    sp_ = PROPS[PROPS.rule == "state_raw_mouse_mean"].pivot_table(index="set", columns="state", values="proportion").fillna(0)
    sp_ = sp_.loc[[s for s in ["primary", "all_tumour", "author_glial", "author_all", "human_cells", "human_nuclei", "human_all"] if s in sp_.index]]
    cols = {"OPC_like": "#3b6fb6", "OC_like": "#7fb3e0", "AC_like": "#d98b2b", "MES_like": "#b8433a", "cycling": "#4f9a55", "unassigned": "#bbbbbb"}
    left = np.zeros(len(sp_))
    for st in cols:
        if st in sp_:
            ax.barh(range(len(sp_)), sp_[st], left=left, color=cols[st], label=st); left += sp_[st].values
    ax.set_yticks(range(len(sp_)), sp_.index); ax.invert_yaxis(); ax.set_xlabel("mean proportion per mouse / patient sample")
    ax.legend(fontsize=5, ncol=3, bbox_to_anchor=(0, -0.35), loc="upper left", frameon=False); ax.set_title("b  human state sets, raw argmax", loc="left")
    ax = fig.add_subplot(gs[1, 1])
    lm = MPU[MPU.lineage].set_index("metaprogram")
    yy = np.arange(len(lm))
    ax.scatter(lm.median_tumour, yy, c="k", s=12, label="mouse tumour", zorder=3)
    for g_, c_ in [("astrocyte (normal)", "#d98b2b"), ("granule-like neuron", "#8a63b5"), ("Microglia", "#4f9a55"), ("OPC-like (normal call)", "#7fb3e0")]:
        ax.scatter(lm[f"median_{g_}"], yy, s=8, c=c_, label=g_, alpha=0.8)
    ax.set_yticks(yy, [f"{i} {lm.family[i]}" for i in lm.index]); ax.invert_yaxis(); ax.set_xlabel("median UCell (human MP top-50)")
    ax.legend(fontsize=5, frameon=False, loc="lower right"); ax.set_title("c  human lineage MPs in mouse cells", loc="left")
    ax = fig.add_subplot(gs[2, :])
    sub = OVL[OVL.run == "pooled_primary"].copy(); sub["mlq"] = -np.log10(sub.q.clip(lower=1e-30))
    hs_order = [k for k in hsets if not k.startswith("Liu")] + sorted(k for k in hsets if k.startswith("Liu"))
    Mq = sub.pivot(index="program", columns="human_set", values="mlq")[hs_order]
    im = ax.imshow(Mq.values, cmap="Greys", vmin=0, vmax=10, aspect="auto")
    ax.set_xticks(range(Mq.shape[1]), [h.split("|")[0] + (" " + h.split("|")[1].split(":")[-1][:18] if "|" in h else "").replace("Liu2022_TableS2a_", "Liu ") for h in Mq.columns],
                  rotation=60, ha="right", fontsize=5)
    ax.set_yticks(range(Mq.shape[0]), [f"{p.split('|')[-1]}: {tops.loc[('pooled_primary', p)][:40]}" for p in Mq.index], fontsize=5)
    fig.colorbar(im, ax=ax, shrink=0.8, label="-log10 q (hypergeom.)"); ax.set_title("d  mouse de novo programmes (pooled primary) vs human MPs / Liu 2022", loc="left")
    fig.savefig(OUTS[9], dpi=200, bbox_inches="tight"); plt.close(fig)
    provenance(OUTS, [IN, CALLS, QC, *H.values()], "scripts/m05_state_comparison.py", {"definitions": {k: int(v.sum()) for k, v in DEF.items()}})
