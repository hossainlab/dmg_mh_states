#!/usr/bin/env python
"""m07 — mouse tumour microenvironment vs the human TME subclusters (myeloid, T/NK, vascular).

Inputs : data/processed/mouse_all_hs.h5ad (m03), data/processed/mouse_tumour_calls.csv (m04), data/interim/mouse_qc.csv (m02),
         data/human/tme_subclusters.csv, tme_subcluster_markers.csv (top-25 Wilcoxon genes per human subcluster),
         tme_subcluster_labels.csv.gz (per human cell label, step 17 — the same clustering as
         tme_subclusters.csv / markers; the step-22 'clean' table relabels every myeloid cluster lipid_phagocytic_TAM),
         composition_per_sample.csv (human cell-type proportions);
         dmg_atlas.states.TME_SUBTYPE_MARKERS (the marker sets that labelled the human subclusters)
Outputs (never overwritten unless --force):
  data/processed/mouse_tme_labels.csv        per mouse TME cell: compartment, subcluster, label (human rule)
  results/tables/m07_tme_subclusters.csv     per mouse subcluster: size, mouse/arm breadth, author-label mix, label, best score,
                                             margin, top markers, best-overlapping human subcluster (hypergeom q)
  results/tables/m07_tme_marker_scores.csv   mean score_genes per mouse subcluster x human marker set
  results/tables/m07_tme_vs_human_markers.csv every mouse x human subcluster (same compartment): top-25 marker overlap, q
  results/tables/m07_tme_label_proportions.csv  label proportions within compartment: mouse (cells; per-mouse mean) vs human
  results/tables/m07_celltype_proportions.csv   compartment proportions of non-tumour cells, mouse vs human (per sample mean)
  results/figures/m07_tme_comparison.png
Run: nice -n 10 python scripts/m07_tme_comparison.py [--force]

Method (mirrors ../scripts/17_tme_subclusters.py where the data allow)
  Cells: qc_pass (m02) and tumour_call == 'normal' (m04). Compartments by author label: myeloid = Microglia, Mono/Macro,
    Neutrophils, pDC; T_NK = T cells, NK cells; vascular = Vascular cells, Meningeal fibroblasts (VLMC). B/Plasma and choroid
    plexus have no human compartment; the mouse has no oligodendroglial compartment (no Mbp/Plp1/Mog-high cluster exists).
  Per compartment: log-normalised counts, 2,000 seurat_v3 HVGs, PCA 30, Harmony on mouse (sample) in place of the human scVI
    latent, kNN 15, Leiden 0.5 (seed 0). DE: Wilcoxon one-vs-rest, top 25 (human symbols of one-to-one orthologs).
  Label: score_genes per human marker set (orthologs present); subcluster label = highest mean set; '|low_confidence' if that
    mean <= 0.05 or the margin to the runner-up < 0.05 (the human rule). Marker overlap: hypergeometric, universe = one-to-one
    genes, BH within compartment.
  Proportions: human labels come from tme_subcluster_labels.csv.gz (step 17) ('|low_confidence' stripped for the comparison);
    human cell-type proportions from composition_per_sample.csv set 'celltype' (malignant excluded, renormalised).
"""
import sys, argparse
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "src"), str(ROOT.parent / "src")]
import dmg_atlas.utils as U  # noqa: E402
import numpy as np, pandas as pd, anndata as ad, scanpy as sc  # noqa: E402
import scipy.stats as ss  # noqa: E402
from mouse_cmp import INTERIM, PROCESSED, RESULTS, LOGS, HUMAN  # noqa: E402
from mouse_cmp.common import provenance, skip_existing, bh  # noqa: E402
from dmg_atlas.states import TME_SUBTYPE_MARKERS  # noqa: E402
from dmg_atlas.figstyle import apply_figure_style  # noqa: E402
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt  # noqa: E402

ap = argparse.ArgumentParser(); ap.add_argument("--force", action="store_true"); args = ap.parse_args()
IN, CALLS, QC = PROCESSED / "mouse_all_hs.h5ad", PROCESSED / "mouse_tumour_calls.csv", INTERIM / "mouse_qc.csv"
H = {k: HUMAN / v for k, v in dict(sub="tme_subclusters.csv", mk="tme_subcluster_markers.csv", lab="tme_subcluster_labels.csv.gz",
                                    comp="composition_per_sample.csv").items()}
T, FIG = RESULTS / "tables", RESULTS / "figures"
OUTS = [PROCESSED / "mouse_tme_labels.csv", T / "m07_tme_subclusters.csv", T / "m07_tme_marker_scores.csv", T / "m07_tme_vs_human_markers.csv",
        T / "m07_tme_label_proportions.csv", T / "m07_celltype_proportions.csv", FIG / "m07_tme_comparison.png"]
log = U.StepLogger(LOGS / "m07_tme_comparison.log", "m07")
if skip_existing(OUTS, args.force, log):
    sys.exit(0)
COMP = {"myeloid": ["Microglia", "Mono/Macro", "Neutrophils", "pDC"], "T_NK": ["T cells", "NK cells"],
        "vascular": ["Vascular cells", "Meningeal fibroblasts (VLMC)"]}
RES, N_NEIGH, SEED, MIN_SCORE, MIN_MARGIN, TOP_DE, MIN_MOUSE = 0.5, 15, 0, 0.05, 0.05, 25, 20

with log.step("load"):
    a = ad.read_h5ad(IN)
    calls = pd.read_csv(CALLS, index_col=0).reindex(a.obs_names); qc = pd.read_csv(QC, index_col=0).reindex(a.obs_names)
    keep = (qc.qc_pass & (calls.tumour_call == "normal")).values
    hmk = pd.read_csv(H["mk"]); hsub = pd.read_csv(H["sub"]).set_index("subcluster")
    universe = set(a.var_names); N = len(universe)

labels, subrows, scorerows, ovl = [], [], [], []
with log.step("subcluster + label per compartment"):
    import harmonypy
    for comp, labs in COMP.items():
        b = a[keep & a.obs.author_label.isin(labs).values].copy()
        sc.pp.highly_variable_genes(b, n_top_genes=2000, flavor="seurat_v3", layer="counts")
        sc.pp.normalize_total(b, target_sum=1e4); sc.pp.log1p(b)
        sc.pp.pca(b, n_comps=30, mask_var="highly_variable", random_state=SEED)
        ho = harmonypy.run_harmony(b.obsm["X_pca"], b.obs[["sample"]].astype(str), "sample", random_state=SEED, verbose=False)
        Zh = np.asarray(ho.Z_corr); b.obsm["X_harmony"] = Zh if Zh.shape[0] == b.n_obs else Zh.T
        sc.pp.neighbors(b, n_neighbors=N_NEIGH, use_rep="X_harmony", random_state=SEED)
        sc.tl.leiden(b, resolution=RES, key_added="sub", flavor="igraph", n_iterations=2, directed=False, random_state=SEED)
        sc.tl.umap(b, random_state=SEED)
        b.obs["sub"] = comp + "_" + b.obs["sub"].astype(str)
        sets = {k: [g for g in v if g in b.var_names] for k, v in TME_SUBTYPE_MARKERS[comp].items()}
        for k, g in sets.items():
            sc.tl.score_genes(b, g, score_name=f"ms_{k}", random_state=SEED)
        ms = b.obs.groupby("sub")[[f"ms_{k}" for k in sets]].mean(); ms.columns = list(sets)
        sc.tl.rank_genes_groups(b, "sub", method="wilcoxon")
        rg = b.uns["rank_genes_groups"]["names"]
        hm = {s: set(d.sort_values("scores", ascending=False).names.head(TOP_DE)) & universe
              for s, d in hmk[hmk.subcluster.str.startswith(comp + "_")].groupby("subcluster")}
        lab_map = {}
        for s in ms.index:
            r = ms.loc[s].sort_values(ascending=False)
            low = (r.iloc[0] <= MIN_SCORE) or (r.iloc[0] - r.iloc[1] < MIN_MARGIN)
            lab_map[s] = r.index[0] + ("|low_confidence" if low else "")
            top = set(rg[s][:TOP_DE])
            for hs_, hv in hm.items():
                k_ = len(top & hv)
                ovl.append(dict(compartment=comp, mouse_subcluster=s, human_subcluster=hs_, human_label=hsub.label.get(hs_), overlap=k_,
                                p=ss.hypergeom.sf(k_ - 1, N, len(hv), len(top)), shared=";".join(sorted(top & hv))))
            ob = b.obs[b.obs["sub"] == s]
            subrows.append(dict(compartment=comp, subcluster=s, n_cells=len(ob), n_mice=ob["sample"].nunique(),
                                frac_top_mouse=ob["sample"].value_counts(normalize=True).iloc[0], label=lab_map[s], best_set=r.index[0],
                                best_score=r.iloc[0], runner_up=r.index[1], margin=r.iloc[0] - r.iloc[1],
                                author_labels=";".join(f"{k}:{v}" for k, v in ob.author_label.astype(str).value_counts().items()),
                                arms=";".join(f"{k}:{v}" for k, v in ob.arm.astype(str).value_counts().items()), top_markers=";".join(rg[s][:15])))
        scorerows.append(ms.assign(compartment=comp))
        b.obs["label"] = b.obs["sub"].map(lab_map)
        labels.append(pd.DataFrame({"compartment": comp, "subcluster": b.obs["sub"], "label": b.obs["label"],
                                    "umap1": b.obsm["X_umap"][:, 0], "umap2": b.obsm["X_umap"][:, 1]}, index=b.obs_names))
        log.write(f"{comp}: {b.n_obs} cells, {ms.shape[0]} subclusters: {pd.Series(lab_map).to_dict()}")
    LAB = pd.concat(labels); SUB = pd.DataFrame(subrows); OVL = pd.DataFrame(ovl)
    OVL["q"] = np.concatenate([bh(d.p) for _, d in OVL.groupby("compartment", sort=False)])
    best = OVL.sort_values("q").groupby("mouse_subcluster").head(1).set_index("mouse_subcluster")
    SUB = SUB.join(best[["human_subcluster", "human_label", "overlap", "q", "shared"]].rename(columns=lambda c: c if c.startswith("human_") else "human_" + c).add_prefix("best_"), on="subcluster")
    log.write(SUB[["subcluster", "n_cells", "n_mice", "label", "best_score", "margin", "best_human_subcluster", "best_human_label",
                   "best_human_overlap", "best_human_q"]].round(3).to_string(index=False))

with log.step("proportions vs human"):
    o = a.obs
    rows = []
    LAB["label_clean"] = LAB.label.str.replace("|low_confidence", "", regex=False)
    hl = pd.read_csv(H["lab"], usecols=["sample_id", "compartment", "label"])
    hl["label_clean"] = hl.label.str.replace("|low_confidence", "", regex=False)
    for comp in COMP:
        m = LAB[LAB.compartment == comp]
        ms_ = o.loc[m.index, "sample"].astype(str)
        rows += [dict(species="mouse", compartment=comp, label=k, rule="cells", proportion=v, n=len(m)) for k, v in m.label_clean.value_counts(normalize=True).items()]
        pm = pd.crosstab(ms_, m.label_clean, normalize="index")[lambda d: d.index.isin(ms_.value_counts()[lambda x: x >= MIN_MOUSE].index)]
        rows += [dict(species="mouse", compartment=comp, label=k, rule="sample_mean", proportion=v, n=len(pm)) for k, v in pm.mean().items()]
        h_ = hl[hl.compartment == comp]
        rows += [dict(species="human", compartment=comp, label=k, rule="cells", proportion=v, n=len(h_)) for k, v in h_.label_clean.value_counts(normalize=True).items()]
        ph = pd.crosstab(h_.sample_id, h_.label_clean, normalize="index")[lambda d: d.index.isin(h_.sample_id.value_counts()[lambda x: x >= MIN_MOUSE].index)]
        rows += [dict(species="human", compartment=comp, label=k, rule="sample_mean", proportion=v, n=len(ph)) for k, v in ph.mean().items()]
    PROP = pd.DataFrame(rows)
    # compartment proportions among non-tumour cells
    CT = {"Microglia": "myeloid", "Mono/Macro": "myeloid", "Neutrophils": "myeloid", "pDC": "myeloid", "T cells": "T_NK", "NK cells": "T_NK",
          "B/Plasma": "B_plasma", "Vascular cells": "vascular", "Meningeal fibroblasts (VLMC)": "vascular", "Choroid plexus epithelial cells": "choroid_plexus",
          "Neuronal-like tumor": "neuron", "Purkinje-like tumor": "neuron", "Astrocyte-like tumor cells": "astrocyte", "OPC-like tumor cells": "OPC_like_unresolved"}
    nt = (calls.tumour_call == "normal").values
    ct = o.author_label.astype(str)[nt].map(CT); sm = o["sample"].astype(str)[nt]
    crow = [dict(species="mouse", rule="cells_all", celltype=k, proportion=v) for k, v in ct.value_counts(normalize=True).items()]
    crow += [dict(species="mouse", rule="sample_mean_all", celltype=k, proportion=v) for k, v in pd.crosstab(sm, ct, normalize="index").mean().items()]
    qn = qc.qc_pass.values[nt]
    crow += [dict(species="mouse", rule="sample_mean_qc_pass", celltype=k, proportion=v) for k, v in pd.crosstab(sm[qn], ct[qn], normalize="index").mean().items()]
    hc = pd.read_csv(H["comp"]); hc = hc[(hc["set"] == "celltype") & (hc.category != "malignant")]
    hc["p"] = hc.k / hc.groupby("sample_id").k.transform("sum")
    for ac, d in [("all", hc)] + list(hc.groupby("assay_class")):
        crow += [dict(species="human", rule=f"sample_mean_{ac}", celltype=k, proportion=v) for k, v in d.pivot_table(index="sample_id", columns="category", values="p").fillna(0).mean().items()]
    CTP = pd.DataFrame(crow)
    log.write(PROP[PROP.rule == "sample_mean"].pivot_table(index=["compartment", "label"], columns="species", values="proportion").round(3).to_string())
    log.write(CTP.pivot_table(index="celltype", columns=["species", "rule"], values="proportion").round(3).to_string())

with log.step("write + figure"):
    LAB.drop(columns="label_clean").to_csv(OUTS[0]); SUB.to_csv(OUTS[1], index=False); pd.concat(scorerows).to_csv(OUTS[2]); OVL.to_csv(OUTS[3], index=False)
    PROP.to_csv(OUTS[4], index=False); CTP.to_csv(OUTS[5], index=False)
    apply_figure_style()
    fig, axes = plt.subplots(2, 3, figsize=(7.4, 5.6), gridspec_kw=dict(height_ratios=[1, 1.1], hspace=0.55, wspace=0.6))
    for j, comp in enumerate(COMP):
        ax = axes[0, j]; m = LAB[LAB.compartment == comp]
        cats = sorted(m.label_clean.unique()); cmap = plt.get_cmap("tab10")
        for i, c in enumerate(cats):
            s_ = m.label_clean == c
            ax.scatter(m.umap1[s_], m.umap2[s_], s=1, color=cmap(i % 10), label=c, rasterized=True)
        ax.set_xticks([]); ax.set_yticks([]); ax.set_title(f"{'abc'[j]}  mouse {comp} (n={len(m)})", loc="left")
        ax.legend(fontsize=4.5, markerscale=4, frameon=False, loc="upper left", bbox_to_anchor=(0, -0.02), ncol=1)
        ax = axes[1, j]
        P = PROP[(PROP.compartment == comp) & (PROP.rule == "sample_mean")].pivot_table(index="label", columns="species", values="proportion").fillna(0)
        P = P.sort_values("human", ascending=True)
        y = np.arange(len(P)); ax.barh(y - 0.2, P.get("human", 0), 0.4, color="0.55", label="human (sample mean)")
        ax.barh(y + 0.2, P.get("mouse", 0), 0.4, color="#3b6fb6", label="mouse (mouse mean)")
        ax.set_yticks(y, [l.replace("_", " ") for l in P.index], fontsize=5); ax.set_xlabel("proportion of compartment")
        ax.set_title(f"{'def'[j]}  {comp} labels", loc="left")
        if j == 0:
            ax.legend(fontsize=5, frameon=False, loc="lower right")
    fig.savefig(OUTS[6], dpi=200, bbox_inches="tight"); plt.close(fig)
    provenance(OUTS, [IN, CALLS, QC, *H.values()], "scripts/m07_tme_comparison.py")
