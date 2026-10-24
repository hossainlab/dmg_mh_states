#!/usr/bin/env python
"""m06 — project the ortholog-mapped mouse cells onto the human scANVI reference (zero-shot).

Complements m05: m05 compares *cluster* pseudobulk profiles with the human classes, which needs no model; this step asks
the human reference model itself to label every mouse cell, and — more importantly — measures how far outside the human
manifold those cells sit, calibrated against the five human query cohorts that went through the identical code path
(../scripts/11g_remap_queries_zeroshot.py).

Inputs : data/processed/mouse_all_hs.h5ad (m03; one-to-one orthologs, human atlas symbols, raw counts in layers['counts'])
         data/processed/mouse_tumour_calls.csv (m04), data/interim/mouse_qc.csv (m02)
         ../models/scanvi/ (reference scANVI: model.pt, latent_f16.npy, obs_names.txt)
         ../results/tables/atlas_core_obs_labels.csv (reference scanvi_label per cell)
         ../results/tables/phase5_reference_knn_distances.csv (reference 15-NN distance distribution, 50,000 cells)
         ../data/processed/dmg_atlas_v1.h5ad (human X_umap, for placing mouse cells on the reference embedding)
         ../data/processed/query_mapping.h5ad (the five human queries, for calibration)
Outputs (never overwritten unless --force):
  results/tables/m06_projection_per_cell.csv.gz   per mouse cell: predicted lineage, posterior, kNN transfer, uncertainty, UMAP placement
  results/tables/m06_projection_by_label.csv      per author label: predicted composition, posterior, uncertainty, agreement
  results/tables/m06_calibration.csv              mouse vs each human query cohort on the same confidence statistics
  results/tables/m06_tumour_concordance.csv       m04 tumour call x predicted lineage
  results/figures/m06_reference_projection.png
Run: nice -n 10 python scripts/m06_reference_mapping.py [--force]

Method
  Zero-shot: SCANVI.prepare_query_anndata() aligns the mouse matrix to the reference's 3,000 HVG symbols (genes without a
  one-to-one ortholog are zero-filled) and load_query_data() aligns the registry; the model is NOT trained, so the encoder
  and classifier weights stay exactly the reference's. That identity is verified tensor-by-tensor against model.pt and the
  maximum absolute difference is written into the provenance sidecar (it must be 0).
  Registry columns are filled as for the human queries: lineage_label='Unknown', assay='scRNA' (whole-cell 10x),
  study='mouse_dipg1', sample_id=the mouse, pct_mito from m02.
  kNN transfer: 15 nearest reference cells in the 30-d latent (pynndescent), majority scanvi_label vote; uncertainty_knn is
  the percentile of the mouse cell's mean 15-NN distance within the reference's own 15-NN distance distribution, so 1.0
  means "further from the reference than any human reference cell".
  UMAP placement is the mean X_umap of those 15 neighbours — a placement for display only, not a UMAP transform.
"""
import sys, argparse, json, time
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "src"), str(ROOT.parent / "src")]
import dmg_atlas.utils as U  # noqa: E402
import numpy as np, pandas as pd, anndata as ad  # noqa: E402
from mouse_cmp import INTERIM, PROCESSED, RESULTS, LOGS, ATLAS_ROOT  # noqa: E402
from mouse_cmp.common import provenance, skip_existing  # noqa: E402

K = 15
REF_DIR = ATLAS_ROOT / "models/scanvi"
LOG = U.StepLogger(LOGS / "m06_reference_mapping.log", context="m06_reference_mapping")
OUT = {
    "per_cell": RESULTS / "tables/m06_projection_per_cell.csv.gz",
    "by_label": RESULTS / "tables/m06_projection_by_label.csv",
    "calib":    RESULTS / "tables/m06_calibration.csv",
    "conc":     RESULTS / "tables/m06_tumour_concordance.csv",
    "fig":      RESULTS / "figures/m06_reference_projection.png",
}


def encode(mouse):
    """Zero-shot encode: reference weights unchanged (verified)."""
    import torch, scvi
    from scvi.model import SCANVI
    torch.set_num_threads(8)
    scvi.settings.seed = 1
    scvi.settings.dl_num_workers = 0
    q = mouse.copy()
    q.X = q.layers["counts"]
    q.obs["lineage_label"] = "Unknown"
    q.obs["assay"] = "scRNA"
    q.obs["study"] = "mouse_dipg1"
    q.obs["sample_id"] = q.obs["sample"].astype(str).values
    q.obs["pct_mito"] = q.obs["percent.mt"].astype(np.float32).values
    for c in ["sample_id", "study", "assay", "lineage_label"]:
        q.obs[c] = q.obs[c].astype("category")
    ref_vars = list(torch.load(REF_DIR / "model.pt", map_location="cpu", weights_only=False)[
        "attr_dict"]["registry_"]["field_registries"]["X"]["state_registry"]["column_names"])
    n_missing = int(len(set(ref_vars) - set(map(str, q.var_names))))
    SCANVI.prepare_query_anndata(q, str(REF_DIR))
    model = SCANVI.load_query_data(q, str(REF_DIR), accelerator="cpu", device="auto")
    model.module.eval()
    model.is_trained = True
    latent = model.get_latent_representation(batch_size=4096)
    probs = model.predict(soft=True, batch_size=4096)
    ref_sd = torch.load(REF_DIR / "model.pt", map_location="cpu", weights_only=False)["model_state_dict"]
    qsd = model.module.state_dict()
    diffs = {k: float((qsd[k].float() - ref_sd[k].float()).abs().max()) for k in ref_sd
             if (k.startswith("z_encoder") or k.startswith("classifier")) and k in qsd and qsd[k].shape == ref_sd[k].shape}
    verify = {"tensors_compared": len(diffs), "max_abs_diff_vs_reference": max(diffs.values()) if diffs else None,
              "n_reference_hvg": len(ref_vars), "n_reference_hvg_absent_from_mouse_orthologs": n_missing,
              "frac_reference_hvg_present": round(1 - n_missing / len(ref_vars), 4), "scvi_tools": scvi.__version__}
    LOG.write(f"zero-shot verification: {verify}")
    return latent.astype(np.float32), probs, verify


def main(force: bool):
    if skip_existing(OUT.values(), force, LOG):
        return
    with LOG.step("load mouse ortholog object + m02/m04 annotations"):
        mouse = ad.read_h5ad(PROCESSED / "mouse_all_hs.h5ad")
        qc = pd.read_csv(INTERIM / "mouse_qc.csv", index_col=0)
        calls = pd.read_csv(PROCESSED / "mouse_tumour_calls.csv", index_col="cell_id")
        mouse.obs["qc_pass"] = qc["qc_pass"].reindex(mouse.obs_names).fillna(False).astype(bool).values
        mouse.obs["tumour_call"] = calls["tumour_call"].reindex(mouse.obs_names).fillna("normal").astype(str).values
        LOG.write(f"mouse {mouse.shape}, qc_pass {int(mouse.obs.qc_pass.sum()):,}, "
                  f"tumour {int((mouse.obs.tumour_call == 'tumour').sum()):,}")

    with LOG.step("zero-shot encode with the reference scANVI"):
        t0 = time.time()
        latent, probs, verify = encode(mouse)
        pred = probs.idxmax(axis=1).values.astype(str)
        maxprob = probs.max(axis=1).values.astype(np.float32)
        LOG.write(f"{mouse.n_obs:,} cells in {time.time()-t0:.0f}s; predicted {pd.Series(pred).value_counts().to_dict()}")

    with LOG.step("kNN transfer against the reference latent"):
        import pynndescent
        Zref = np.load(REF_DIR / "latent_f16.npy").astype(np.float32)
        ref_names = np.array((REF_DIR / "obs_names.txt").read_text().splitlines())
        ref_lab = (pd.read_csv(ATLAS_ROOT / "results/tables/atlas_core_obs_labels.csv", index_col=0)
                   .reindex(ref_names)["scanvi_label"].astype(str).values)
        index = pynndescent.NNDescent(Zref, n_neighbors=30, metric="euclidean", random_state=0, n_jobs=8)
        index.prepare()
        nn_idx, dist = index.query(latent, k=K)
        labs = sorted(set(ref_lab))
        votes = ref_lab[nn_idx]
        counts = np.stack([(votes == l).sum(1) for l in labs], 1)
        knn_lineage = np.array(labs)[counts.argmax(1)]
        knn_vote_frac = (counts.max(1) / K).astype(np.float32)
        knn_mean_dist = dist.mean(1).astype(np.float32)
        ref_d = np.sort(pd.read_csv(ATLAS_ROOT / "results/tables/phase5_reference_knn_distances.csv")["knn_mean_dist"].values)
        uncertainty_knn = (np.searchsorted(ref_d, knn_mean_dist) / len(ref_d)).astype(np.float32)
        H = ad.read_h5ad(ATLAS_ROOT / "data/processed/dmg_atlas_v1.h5ad", backed="r")
        hu = pd.DataFrame(H.obsm["X_umap"], index=H.obs_names).reindex(ref_names).values
        umap_ref = np.stack([hu[nn_idx, 0].mean(1), hu[nn_idx, 1].mean(1)], 1).astype(np.float32)
        LOG.write(f"kNN: classifier/kNN agreement {np.mean(knn_lineage == pred):.3f}; "
                  f"median uncertainty_knn {np.median(uncertainty_knn):.3f}")

    with LOG.step("assemble per-cell table"):
        o = mouse.obs
        per = pd.DataFrame({
            "cell_id": mouse.obs_names, "sample": o["sample"].astype(str).values, "arm": o["arm"].astype(str).values,
            "author_label": o["author_label"].astype(str).values, "author_compartment": o["author_compartment"].astype(str).values,
            "qc_pass": o["qc_pass"].values, "tumour_call": o["tumour_call"].values,
            "predicted_lineage": pred, "scanvi_maxprob": maxprob, "knn_lineage": knn_lineage,
            "knn_vote_frac": knn_vote_frac, "knn_mean_dist": knn_mean_dist, "uncertainty_knn": uncertainty_knn,
            "classifier_knn_agree": knn_lineage == pred,
            "umap_ref_1": umap_ref[:, 0], "umap_ref_2": umap_ref[:, 1]}).set_index("cell_id")
        for c in probs.columns:
            per[f"prob_{c}"] = probs[c].values.astype(np.float32)
        per.to_csv(OUT["per_cell"])

    with LOG.step("per-label summary, tumour concordance, calibration vs the human queries"):
        comp = pd.crosstab(per.author_label, per.predicted_lineage, normalize="index")
        by = per.groupby("author_label").agg(**{
            "n_cells": ("scanvi_maxprob", "size"), "qc_pass": ("qc_pass", "sum"),
            "median_maxprob": ("scanvi_maxprob", "median"),
            "frac_maxprob_ge_0.8": ("scanvi_maxprob", lambda s: (s >= 0.8).mean()),
            "median_uncertainty_knn": ("uncertainty_knn", "median"),
            "frac_classifier_knn_agree": ("classifier_knn_agree", "mean")})
        by.join(comp.add_prefix("pred_")).round(4).to_csv(OUT["by_label"])

        prim = per[per.qc_pass]
        conc = (per.assign(predicted_malignant=per.predicted_lineage.eq("malignant"))
                .groupby(["tumour_call", "qc_pass"])
                .agg(n=("predicted_lineage", "size"), frac_predicted_malignant=("predicted_malignant", "mean"),
                     median_maxprob=("scanvi_maxprob", "median"), median_uncertainty_knn=("uncertainty_knn", "median"))
                .round(4))
        conc.to_csv(OUT["conc"])

        Q = ad.read_h5ad(ATLAS_ROOT / "data/processed/query_mapping.h5ad", backed="r").obs
        rows = [{"cohort": f"mouse ({name})", "n_cells": len(d), "species": "mouse",
                 "median_maxprob": d.scanvi_maxprob.median(), "frac_maxprob_ge_0.8": (d.scanvi_maxprob >= 0.8).mean(),
                 "median_uncertainty_knn": d.uncertainty_knn.median(),
                 "frac_uncertainty_knn_gt_0.9": (d.uncertainty_knn > 0.9).mean(),
                 "frac_classifier_knn_agree": d.classifier_knn_agree.mean(),
                 "frac_predicted_malignant": d.predicted_lineage.eq("malignant").mean()}
                for name, d in [("all cells", per), ("QC-pass", prim),
                                ("m04 tumour, QC-pass", per[(per.tumour_call == "tumour") & per.qc_pass])]]
        for st, g in Q.groupby("study", observed=True):
            rows.append({"cohort": f"human query: {st}", "n_cells": len(g), "species": "human",
                         "median_maxprob": g.scanvi_maxprob.median(), "frac_maxprob_ge_0.8": (g.scanvi_maxprob >= 0.8).mean(),
                         "median_uncertainty_knn": g.uncertainty_knn.median(),
                         "frac_uncertainty_knn_gt_0.9": (g.uncertainty_knn > 0.9).mean(),
                         "frac_classifier_knn_agree": g.classifier_knn_agree.astype(str).eq("True").mean(),
                         "frac_predicted_malignant": g.predicted_lineage.astype(str).eq("malignant").mean()})
        calib = pd.DataFrame(rows).round(4)
        calib.to_csv(OUT["calib"], index=False)
        LOG.write("calibration:\n" + calib.to_string(index=False))

    with LOG.step("figure"):
        figure(per, calib, comp)
    provenance(list(OUT.values()),
               [PROCESSED / "mouse_all_hs.h5ad", PROCESSED / "mouse_tumour_calls.csv", INTERIM / "mouse_qc.csv",
                REF_DIR / "model.pt", ATLAS_ROOT / "data/processed/query_mapping.h5ad"],
               "m06_reference_mapping.py", extra={"zero_shot_verification": verify, "k_neighbours": K})
    LOG.write(f"done; peak RSS {U.peak_rss_gb():.2f} GB")


def figure(per, calib, comp):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    sys.path.insert(0, str(ATLAS_ROOT / "src"))
    from dmg_atlas.figstyle import apply_figure_style, letter
    apply_figure_style()
    fig, axes = plt.subplots(2, 2, figsize=(11, 9))

    ax = axes[0, 0]
    H = ad.read_h5ad(ATLAS_ROOT / "data/processed/dmg_atlas_v1.h5ad", backed="r")
    hu = H.obsm["X_umap"][:]
    ax.scatter(hu[:, 0], hu[:, 1], s=0.3, c="0.88", lw=0, rasterized=True)
    pal = {"malignant": "#de2d26", "myeloid": "#2b6cb0", "T_NK": "#38a169", "oligodendrocyte": "#805ad5",
           "astrocyte": "#dd6b20", "neuron": "#718096", "vascular": "#d53f8c"}
    for lab, g in per.groupby("predicted_lineage"):
        ax.scatter(g.umap_ref_1, g.umap_ref_2, s=1.2, lw=0, rasterized=True,
                   color=pal.get(lab, "0.4"), label=f"{lab} ({len(g):,})")
    ax.set(xticks=[], yticks=[])
    ax.legend(markerscale=6, fontsize=6, frameon=False, loc="upper left")
    letter(ax, "a", title="Mouse cells placed on the human reference")

    ax = axes[0, 1]

    def ecdf(ax, v, **kw):
        v = np.sort(np.asarray(v, float))
        ax.step(v, np.arange(1, len(v) + 1) / len(v), where="post", **kw)

    Q = ad.read_h5ad(ATLAS_ROOT / "data/processed/query_mapping.h5ad", backed="r").obs
    ecdf(ax, per.uncertainty_knn, color="#de2d26", lw=1.8, label=f"mouse, all ({len(per):,})")
    mt = per[(per.tumour_call == "tumour") & per.qc_pass]
    ecdf(ax, mt.uncertainty_knn, color="#742a2a", lw=1.8, ls="--", label=f"mouse, m04 tumour ({len(mt):,})")
    for st, g in Q.groupby("study", observed=True):
        ecdf(ax, g.uncertainty_knn, lw=0.9, label=f"human: {st} ({len(g):,})")
    ax.set(xlabel="uncertainty_knn (percentile of reference 15-NN distance)",
           ylabel="cumulative fraction of cells", xlim=(0, 1), ylim=(0, 1))
    ax.legend(fontsize=6, frameon=False, loc="upper left")
    letter(ax, "b", title="Distance from the human manifold")

    ax = axes[1, 0]
    sub = comp.loc[comp.index.isin(per.author_label.value_counts().head(12).index)]
    im = ax.imshow(sub.values, cmap="Blues", vmin=0, vmax=1, aspect="auto")
    ax.set(xticks=range(sub.shape[1]), yticks=range(sub.shape[0]), yticklabels=sub.index)
    ax.set_xticklabels(sub.columns, rotation=45, ha="right")
    for i in range(sub.shape[0]):
        for j in range(sub.shape[1]):
            if sub.values[i, j] > 0.02:
                ax.text(j, i, f"{sub.values[i, j]:.2f}", ha="center", va="center", fontsize=5.5,
                        color="white" if sub.values[i, j] > 0.5 else "black")
    plt.colorbar(im, ax=ax, shrink=0.7, label="fraction")
    letter(ax, "c", title="Author label vs transferred lineage")

    ax = axes[1, 1]
    c = calib.set_index("cohort")
    y = np.arange(len(c))
    col = np.where(c.species == "mouse", "#de2d26", "#2b6cb0")
    ax.barh(y - 0.2, c.median_uncertainty_knn, height=0.38, color=col, label="median uncertainty_knn")
    ax.barh(y + 0.2, c.frac_predicted_malignant, height=0.38, color=col, alpha=0.45,
            label="fraction predicted malignant")
    for i, (u, m) in enumerate(zip(c.median_uncertainty_knn, c.frac_predicted_malignant)):
        ax.text(u + 0.01, i - 0.2, f"{u:.2f}", va="center", fontsize=5.5)
        ax.text(m + 0.01, i + 0.2, f"{m:.2f}", va="center", fontsize=5.5)
    ax.set(yticks=y, yticklabels=c.index, xlim=(0, 1.15), xlabel="value")
    ax.tick_params(labelsize=6)
    ax.legend(fontsize=6, frameon=False, loc="lower right")
    ax.invert_yaxis()
    letter(ax, "d", title="Mouse against the human query cohorts")

    fig.tight_layout()
    fig.savefig(OUT["fig"], dpi=300)
    plt.close(fig)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--force", action="store_true")
    main(p.parse_args().force)
