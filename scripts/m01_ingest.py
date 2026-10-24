#!/usr/bin/env python
"""m01 — ingest the mouse DMG deposit (DIPG1 all-cell object) into data/interim with standardised design columns.

Replaces the planned m01_convert_rds.R: the deposit now arrives as an h5ad export of the standardized Seurat object
(uns['source'] = 'DIPG1_allcell.standardized_seurat.rds (raw RNA counts)') plus a meta.data CSV.
Inputs : data/mouse/DIPG1_allcell_raw_counts.h5ad, data/mouse/DIPG1_allcell_metadata.csv
Outputs (never overwritten unless --force):
  data/interim/mouse_all.h5ad          X = layers['counts'] = raw UMI counts (CSR, float32), obs = deposit meta.data verbatim
                                       + design columns: sample, arm (= Group), block (inferred, see common.sample_block),
                                       author_label (= Cell_type_v2, verbatim), author_tumour (label contains 'tumor'),
                                       author_compartment (tumour / myeloid / lymphoid / vascular_stromal / epithelial)
  results/qc/m01_deposit_checks.csv    integrity checks (counts integer, X == counts, nCount/percent.mt reproduced, meta CSV == obs)
  results/qc/m01_design.csv            sample x arm x block x n_cells x author compartment counts
Run: nice -n 10 python scripts/m01_ingest.py [--force]   (from mouse_comparison/)

Notes: the deposit is already filtered (min 200 genes/cell, scDblFinder singlets only) and has 15,304 genes. No ICI arm and
no Sample3 are present; arms are Control (4 mice), FUS (4), RT (3), FUS+RT (4).
"""
import sys, argparse
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "src"), str(ROOT.parent / "src")]
import dmg_atlas.utils as U  # noqa: E402  (thread caps first)
import numpy as np, pandas as pd, scipy.sparse as sp, anndata as ad  # noqa: E402
from mouse_cmp import INTERIM, RESULTS, LOGS  # noqa: E402
from mouse_cmp.common import provenance, skip_existing, sample_block, ARMS  # noqa: E402

ap = argparse.ArgumentParser(); ap.add_argument("--force", action="store_true"); args = ap.parse_args()
H5, META = ROOT / "data/mouse/DIPG1_allcell_raw_counts.h5ad", ROOT / "data/mouse/DIPG1_allcell_metadata.csv"
OUT, CHK, DES = INTERIM / "mouse_all.h5ad", RESULTS / "qc/m01_deposit_checks.csv", RESULTS / "qc/m01_design.csv"
log = U.StepLogger(LOGS / "m01_ingest.log", "m01")
if skip_existing([OUT, CHK, DES], args.force, log):
    sys.exit(0)
COMPARTMENT = {"Neuronal-like tumor": "tumour", "Purkinje-like tumor": "tumour", "OPC-like tumor cells": "tumour",
               "Astrocyte-like tumor cells": "tumour", "Mono/Macro": "myeloid", "Microglia": "myeloid", "Neutrophils": "myeloid",
               "pDC": "myeloid", "T cells": "lymphoid", "NK cells": "lymphoid", "B/Plasma": "lymphoid",
               "Meningeal fibroblasts (VLMC)": "vascular_stromal", "Vascular cells": "vascular_stromal",
               "Choroid plexus epithelial cells": "epithelial"}

with log.step("read + checks"):
    a = ad.read_h5ad(H5)
    meta = pd.read_csv(META, index_col=0)
    C = sp.csr_matrix(a.layers["counts"]).astype(np.float32)
    lib = np.asarray(C.sum(1)).ravel()
    mt = a.var_names.str.startswith("mt-")
    pmt = np.asarray(C[:, mt].sum(1)).ravel() / lib * 100
    num = [c for c in meta.columns if pd.api.types.is_numeric_dtype(meta[c])]
    chk = pd.DataFrame([
        ("n_cells", a.n_obs), ("n_genes", a.n_vars), ("n_mt_genes", int(mt.sum())),
        ("counts_integer", bool(np.allclose(C.data, np.round(C.data)))), ("X_equals_counts", (sp.csr_matrix(a.X) != C).nnz == 0),
        ("nCount_RNA_reproduced_max_abs_diff", float(np.abs(lib - a.obs.nCount_RNA.values).max())),
        ("percent_mt_reproduced_max_abs_diff", float(np.abs(pmt - a.obs["percent.mt"].values).max())),
        ("meta_csv_same_cells_same_order", bool((meta.index == a.obs_names).all())),
        ("meta_csv_numeric_cols_equal_obs", bool(all(np.allclose(meta[c].values, a.obs[c].astype(float).values, equal_nan=True) for c in num if c in a.obs))),
        ("all_scDblFinder_singlet", bool((a.obs["scDblFinder.class"] == "singlet").all())),
        ("min_nFeature_RNA", int(a.obs.nFeature_RNA.min())), ("duplicate_gene_names", int(a.var_names.duplicated().sum())),
        ("uns_source", str(a.uns.get("source", ""))),
    ], columns=["check", "value"])
    log.write("; ".join(f"{r.check}={r.value}" for r in chk.itertuples()))
    assert chk.set_index("check").loc[["counts_integer", "X_equals_counts", "meta_csv_same_cells_same_order"], "value"].all()

with log.step("design columns + write"):
    o = a.obs.copy()
    o["sample"] = o["Sample"].astype(str).astype("category")
    o["arm"] = pd.Categorical(o["Group"].astype(str), categories=ARMS)
    assert o.arm.notna().all(), "unexpected arm label"
    o["block"] = sample_block(o["sample"])
    o["author_label"] = o["Cell_type_v2"].astype(str).astype("category")
    o["author_tumour"] = o.author_label.str.contains("tumor").values
    o["author_compartment"] = o.author_label.map(COMPARTMENT).astype("category")
    assert o.author_compartment.notna().all(), set(o.author_label[o.author_compartment.isna()])
    out = ad.AnnData(X=C, obs=o, var=pd.DataFrame(index=a.var_names.astype(str)), obsm=dict(a.obsm), uns={"source": a.uns.get("source", "")})
    out.layers["counts"] = C
    INTERIM.mkdir(parents=True, exist_ok=True); out.write_h5ad(OUT, compression="gzip")
    des = (o.groupby(["sample", "arm", "block"], observed=True).size().rename("n_cells").to_frame()
           .join(pd.crosstab(o["sample"], o.author_compartment)).reset_index())
    des = des.sort_values(["arm", "block"])
    (RESULTS / "qc").mkdir(parents=True, exist_ok=True)
    chk.to_csv(CHK, index=False); des.to_csv(DES, index=False)
    provenance([OUT, CHK, DES], [H5, META], "scripts/m01_ingest.py")
    log.write("design:\n" + des.to_string(index=False))
