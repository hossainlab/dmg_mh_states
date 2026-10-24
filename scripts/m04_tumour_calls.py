#!/usr/bin/env python
"""m04 — tumour vs normal calls for the mouse cells: inferred CNV, a shared tumour-clone CNV profile, and subclusters.

Why: the author labels call 64 % of cells tumour ('Neuronal-like tumor', 'Purkinje-like tumor'), but those clusters carry the
cerebellar granule-cell (Cbln3, Gabra6, Neurod1, Zic1) and Purkinje-cell (Pcp2, Pvalb, Car8, Calb1) programmes with
45-52 % mitochondrial reads. The model has no reporter / transgene feature in the count matrix, so tumour identity is tested by
copy number. Author labels are kept verbatim; the call made here is a separate column.

Inputs : data/interim/mouse_all.h5ad (m01), data/interim/mouse_qc.csv (m02), data/external/mouse_human_orthologs.csv (m03:
         mouse gene positions, GRCm39)
Outputs (never overwritten unless --force):
  data/processed/mouse_tumour_calls.csv   per cell: cnv_score, best_ref, clone_r, clone_pos, sub (Leiden subcluster of the
                                          OPC-/astrocyte-like author labels), sub_clone_frac, tumour_call, tumour_basis
  data/processed/mouse_cnv.npz            X_cnv (cells x windows, float16), window chromosome labels, obs names
  results/tables/m04_clone_by_label.csv   per author label: n, median clone_r, fraction clone-positive, median cnv_score
  results/tables/m04_clone_by_sample.csv  per sample: tumour calls, clone-positive OPC-like cells, between-mouse profile correlation
  results/tables/m04_clone_profile.csv    per chromosome: mean CNV of the seed cells in each mouse (shared-clone evidence)
  results/tables/m04_subclusters.csv      per subcluster: n, author label mix, clone fraction, QC medians, top markers, call
  results/figures/m04_tumour_calls.png    (a) clone_r by author label, (b) chromosome CNV of seed cells per mouse,
                                          (c) CNV heatmap of sampled cells per call group
Run: nice -n 10 python scripts/m04_tumour_calls.py [--force]

Method
  CNV: log-normalised counts; genes detected in >= 5 % of cells with a GRCm39 autosomal position; for every cell, the deviation
    from each normal reference label (Microglia, Mono/Macro, T, NK, B/Plasma, Vascular, VLMC, Choroid plexus) is clipped to
    +-3, smoothed by infercnvpy's running mean (100 genes, step 10, per chromosome), median-centred, and the reference with the
    smallest mean |deviation| is kept (same 'best_match' rule as the human pipeline, dmg_atlas.cnv.run_infercnv); values with
    |x| < 1.5 SD are zeroed; cnv_score = mean |x|.
  Clone profile: seed cells = author 'OPC-like tumor cells' with cnv_score above the 99th percentile of the lymphoid/vascular/
    epithelial references. Each cell is correlated (Pearson, over windows) with the mean seed profile of the OTHER mice
    (leave-one-mouse-out, so no mouse defines its own reference); clone_pos = clone_r above the 99th percentile of the same normal
    references. An implanted / clonal model shares its CNV profile across mice; a normal cell has none.
  Subclusters: author OPC-like + astrocyte-like cells, 2,000 seurat_v3 HVGs, PCA 30, kNN 20, Leiden 0.6 (seed 0).
  Call: tumour = (clone_pos AND author label OPC-/astrocyte-like) OR member of a subcluster whose clone-positive fraction >= 0.5
    (catches tumour cells too shallow for a confident single-cell CNV call); every other cell = normal. clone_pos alone is not
    used outside the glial candidate labels: there its rate (0.2-2.2 %) equals the 1 % false-positive rate of the null, i.e. no
    lineage other than the OPC-like cluster carries the clone above chance. tumour_basis records which rule fired. The 0.5 rule and all
    thresholds are fixed here before any human comparison is run.
"""
import sys, argparse
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "src"), str(ROOT.parent / "src")]
import dmg_atlas.utils as U  # noqa: E402
import numpy as np, pandas as pd, scipy.sparse as sp, anndata as ad, scanpy as sc  # noqa: E402
from infercnvpy.tl._infercnv import _running_mean_by_chromosome  # noqa: E402
from mouse_cmp import INTERIM, PROCESSED, EXTERNAL, RESULTS, LOGS  # noqa: E402
from mouse_cmp.common import provenance, skip_existing  # noqa: E402
from dmg_atlas.figstyle import apply_figure_style  # noqa: E402
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt  # noqa: E402

ap = argparse.ArgumentParser(); ap.add_argument("--force", action="store_true"); args = ap.parse_args()
IN, QC, ORTH = INTERIM / "mouse_all.h5ad", INTERIM / "mouse_qc.csv", EXTERNAL / "mouse_human_orthologs.csv"
T = RESULTS / "tables"; FIG = RESULTS / "figures"
OUTS = [PROCESSED / "mouse_tumour_calls.csv", PROCESSED / "mouse_cnv.npz", T / "m04_clone_by_label.csv", T / "m04_clone_by_sample.csv",
        T / "m04_clone_profile.csv", T / "m04_subclusters.csv", FIG / "m04_tumour_calls.png"]
log = U.StepLogger(LOGS / "m04_tumour_calls.log", "m04")
if skip_existing(OUTS, args.force, log):
    sys.exit(0)
REF = ["Microglia", "Mono/Macro", "T cells", "NK cells", "B/Plasma", "Vascular cells", "Meningeal fibroblasts (VLMC)",
       "Choroid plexus epithelial cells"]
NULL = ["Mono/Macro", "T cells", "NK cells", "B/Plasma", "Vascular cells", "Choroid plexus epithelial cells"]   # microglia excluded:
SEED_LABEL, SUB_LABELS = "OPC-like tumor cells", ["OPC-like tumor cells", "Astrocyte-like tumor cells"]  # 4.5 % above the cnv_score
MIN_FRAC, WINDOW, STEP, NOISE_SD, Q_NULL, SUB_FRAC = 0.05, 100, 10, 1.5, 0.99, 0.5                      # null in exploration
SEED = 0

with log.step("CNV (best-match reference)"):
    a = ad.read_h5ad(IN)
    lab = a.obs.author_label.astype(str).values
    orth = pd.read_csv(ORTH).set_index("mouse_symbol")
    orth = orth[orth.chromosome.isin([str(i) for i in range(1, 20)])]
    sc.pp.normalize_total(a, target_sum=1e4); sc.pp.log1p(a)
    frac = np.asarray((a.X > 0).mean(0)).ravel()
    g = [x for x, f in zip(a.var_names, frac) if f >= MIN_FRAC and x in orth.index]
    var = pd.DataFrame({"chromosome": "chr" + orth.loc[g, "chromosome"].astype(int).astype(str).str.zfill(2),
                        "start": orth.loc[g, "start"].astype(int).values, "end": orth.loc[g, "end"].astype(int).values}, index=g)
    var = var.iloc[np.lexsort((var.start.values, var.chromosome.values))]
    X = sp.csr_matrix(a[:, var.index].X)
    means = np.vstack([np.asarray(X[lab == c].mean(0)).ravel() for c in REF]).astype(np.float32)
    out, best = [], []
    for i in range(0, X.shape[0], 4000):
        xd = X[i:i + 4000].toarray(); res = None
        for k in range(len(REF)):
            chr_pos, r, _ = _running_mean_by_chromosome(np.clip(xd - means[k], -3, 3), var, WINDOW, STEP, False)
            r = r - np.median(r, axis=1, keepdims=True); s = np.abs(r).mean(1)
            if res is None:
                res, score, b = r, s, np.zeros(len(xd), int)
            else:
                m = s < score; res[m] = r[m]; score[m] = s[m]; b[m] = k
        out.append(res.astype(np.float32)); best.append(b)
    Xc = np.vstack(out); best = np.concatenate(best)
    Xc[np.abs(Xc) < NOISE_SD * Xc.std()] = 0
    chs = sorted(chr_pos, key=chr_pos.get); bounds = [chr_pos[c] for c in chs] + [Xc.shape[1]]
    win_chr = np.concatenate([[c] * (bounds[i + 1] - bounds[i]) for i, c in enumerate(chs)])
    o = pd.DataFrame({"author_label": lab, "sample": a.obs["sample"].astype(str).values, "arm": a.obs.arm.astype(str).values,
                      "cnv_score": np.abs(Xc).mean(1), "best_ref": np.array(REF)[best]}, index=a.obs_names)
    log.write(f"{len(var)} genes, {Xc.shape[1]} windows")

with log.step("leave-one-mouse-out clone profile"):
    null = o.author_label.isin(NULL).values
    seed = ((o.author_label == SEED_LABEL) & (o.cnv_score > np.quantile(o.cnv_score[null], Q_NULL))).values
    Xz = Xc - Xc.mean(1, keepdims=True); nx = np.linalg.norm(Xz, axis=1); nx[nx == 0] = np.nan
    r = np.full(len(o), np.nan); prof_rows = []
    for s in sorted(o["sample"].unique()):
        m = (o["sample"] == s).values
        p = Xc[seed & ~m].mean(0); pz = p - p.mean()
        r[m] = (Xz[m] @ pz) / (nx[m] * np.linalg.norm(pz))
        own = Xc[seed & m]
        prof_rows.append({"sample": s, "n_seed": int(own.shape[0]),
                          **({c: float(own[:, win_chr == c].mean()) for c in chs} if len(own) else {})})
    o["clone_r"] = r
    thr_r = float(np.nanquantile(r[null], Q_NULL))
    o["clone_pos"] = o.clone_r > thr_r
    P = np.vstack([Xc[seed & (o["sample"] == s).values].mean(0) for s in sorted(o["sample"].unique()) if (seed & (o["sample"] == s).values).sum() >= 3])
    cc = np.corrcoef(P); off = cc[np.triu_indices_from(cc, 1)]
    log.write(f"seed cells {seed.sum()}; clone_r null p99 {thr_r:.3f}; between-mouse seed profile r median {np.median(off):.2f} "
              f"(range {off.min():.2f}-{off.max():.2f}, {len(P)} mice with >= 3 seed cells)")

with log.step("subclusters of OPC-/astrocyte-like cells"):
    b = ad.read_h5ad(IN)
    b = b[b.obs.author_label.isin(SUB_LABELS)].copy()
    sc.pp.highly_variable_genes(b, n_top_genes=2000, flavor="seurat_v3", layer="counts")
    sc.pp.normalize_total(b, target_sum=1e4); sc.pp.log1p(b)
    sc.pp.pca(b, n_comps=30, mask_var="highly_variable", random_state=SEED); sc.pp.neighbors(b, n_neighbors=20, random_state=SEED)
    sc.tl.leiden(b, resolution=0.6, key_added="sub", flavor="igraph", n_iterations=2, directed=False, random_state=SEED)
    sc.tl.rank_genes_groups(b, "sub", method="wilcoxon")
    o["sub"] = pd.Series("s" + b.obs["sub"].astype(str), index=b.obs_names).reindex(o.index).fillna("none")
    qc = pd.read_csv(QC, index_col=0).reindex(o.index)
    o["qc_pass"] = qc.qc_pass.values
    sf = o[o["sub"] != "none"].groupby("sub").clone_pos.mean()
    o["sub_clone_frac"] = o["sub"].map(sf)
    rg = b.uns["rank_genes_groups"]["names"]
    subs = (o[o["sub"] != "none"].groupby("sub").agg(n=("clone_r", "size"), n_samples=("sample", "nunique"), clone_frac=("clone_pos", "mean"),
            median_clone_r=("clone_r", "median"), frac_qc_pass=("qc_pass", "mean"))
            .join(qc.groupby(o["sub"]).agg(median_genes=("n_genes", "median"), median_pct_mito=("pct_mito", "median")))
            .join(pd.crosstab(o["sub"], o.author_label)))
    subs["top_markers"] = [";".join(rg[k.lstrip("s")][:15]) for k in subs.index]
    subs["tumour_subcluster"] = subs.clone_frac >= SUB_FRAC

with log.step("calls + tables + figure"):
    tsub = set(subs.index[subs.tumour_subcluster])
    in_sub = o["sub"].isin(tsub)
    cpos = o.clone_pos & o.author_label.isin(SUB_LABELS)
    o["tumour_call"] = np.where(cpos | in_sub, "tumour", "normal")
    o["tumour_basis"] = np.select([cpos & in_sub, cpos, in_sub], ["clone+subcluster", "clone_only", "subcluster_only"], "none")
    subs["call"] = np.where(subs.tumour_subcluster, "tumour", "normal-majority")
    PROCESSED.mkdir(parents=True, exist_ok=True); T.mkdir(parents=True, exist_ok=True); FIG.mkdir(parents=True, exist_ok=True)
    o.to_csv(OUTS[0])
    np.savez_compressed(OUTS[1], X_cnv=Xc.astype(np.float16), win_chr=win_chr, obs_names=o.index.values.astype(str))
    bylab = o.groupby("author_label").agg(n=("clone_r", "size"), median_clone_r=("clone_r", "median"), frac_clone_pos=("clone_pos", "mean"),
                                          median_cnv_score=("cnv_score", "median"), n_tumour_call=("tumour_call", lambda x: int((x == "tumour").sum())))
    bylab.to_csv(OUTS[2])
    opc = o.author_label == SEED_LABEL
    bys = (o.groupby("sample").agg(arm=("arm", "first"), n_cells=("clone_r", "size"),
                                   n_tumour=("tumour_call", lambda x: int((x == "tumour").sum())))
           .join(o[opc].groupby("sample").agg(n_opc_like=("clone_r", "size"), frac_opc_like_clone_pos=("clone_pos", "mean")))
           .join(o[o.tumour_call == "tumour"].groupby("sample").agg(n_tumour_qc_pass=("qc_pass", "sum"))))
    bys["frac_tumour"] = bys.n_tumour / bys.n_cells
    bys.to_csv(OUTS[3])
    pd.DataFrame(prof_rows).to_csv(OUTS[4], index=False)
    subs.to_csv(OUTS[5])

    apply_figure_style()
    fig = plt.figure(figsize=(7.2, 7.2)); gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.15], hspace=0.45, wspace=0.35)
    ax = fig.add_subplot(gs[0, 0])
    order = bylab.sort_values("median_clone_r").index
    ax.boxplot([o.clone_r[o.author_label == l].dropna() for l in order], vert=False, showfliers=False, widths=0.6)
    ax.set_yticks(range(1, len(order) + 1), order); ax.axvline(thr_r, ls="--", lw=0.7, c="k")
    ax.set_xlabel("correlation with shared tumour CNV profile\n(leave-one-mouse-out)"); ax.set_title("a  CNV clone signal by author label", loc="left")
    ax = fig.add_subplot(gs[0, 1])
    PR = pd.DataFrame(prof_rows).set_index("sample"); PR = PR[PR.n_seed >= 3].drop(columns="n_seed")
    im = ax.imshow(PR.values * 100, aspect="auto", cmap="RdBu_r", vmin=-8, vmax=8)
    ax.set_xticks(range(PR.shape[1]), [c.replace("chr0", "").replace("chr", "") for c in PR.columns]); ax.set_yticks(range(len(PR)), PR.index)
    ax.set_xlabel("chromosome"); ax.set_title("b  seed-cell CNV per mouse", loc="left"); fig.colorbar(im, ax=ax, shrink=0.7, label="mean CNV x100")
    ax = fig.add_subplot(gs[1, :])
    rng = np.random.default_rng(SEED); grp = []
    for name, m in [("tumour (call)", o.tumour_call == "tumour"), ("OPC-like, normal call", opc & (o.tumour_call == "normal")),
                    ("granule-like (author tumour)", o.author_label == "Neuronal-like tumor"), ("Purkinje-like (author tumour)", o.author_label == "Purkinje-like tumor"),
                    ("immune / stromal", o.author_label.isin(REF))]:
        idx = np.flatnonzero(m.values); idx = rng.choice(idx, min(300, len(idx)), replace=False)
        idx = idx[np.argsort(-o.clone_r.values[idx])]; grp.append((name, idx))
    H = np.vstack([Xc[i] for _, i in grp])
    ax.imshow(H, aspect="auto", cmap="RdBu_r", vmin=-0.25, vmax=0.25, interpolation="nearest")
    y = 0
    for name, i in grp:
        ax.axhline(y - 0.5, c="k", lw=0.5); ax.text(-3, y + len(i) / 2, name, ha="right", va="center", fontsize=6); y += len(i)
    for c in chs[1:]:
        ax.axvline(chr_pos[c] - 0.5, c="0.6", lw=0.4)
    ax.set_xticks([(chr_pos[c] + bounds[i + 1]) / 2 for i, c in enumerate(chs)], [c.replace("chr0", "").replace("chr", "") for c in chs])
    ax.set_yticks([]); ax.set_xlabel("chromosome (GRCm39 gene order)"); ax.set_title("c  inferred CNV, 300 random cells per group (sorted by clone_r)", loc="left")
    fig.savefig(OUTS[6], dpi=200, bbox_inches="tight"); plt.close(fig)
    provenance(OUTS, [IN, QC, ORTH], "scripts/m04_tumour_calls.py",
               {"clone_r_threshold": thr_r, "n_seed": int(seed.sum()), "tumour_subclusters": sorted(tsub)})
    log.write(f"tumour calls {(o.tumour_call == 'tumour').sum()} ({o.tumour_basis.value_counts().to_dict()}); qc_pass tumour "
              f"{int((o.qc_pass & (o.tumour_call == 'tumour')).sum())}\nby label:\n{bylab.round(3).to_string()}\nsubclusters:\n"
              f"{subs.drop(columns='top_markers').round(2).to_string()}\nby sample:\n{bys.round(3).to_string()}")
