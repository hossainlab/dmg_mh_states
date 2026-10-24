#!/usr/bin/env python
"""m02 — QC flags for the mouse cells (flag, never remove), using the human atlas 'cells' thresholds (../configs/qc.yaml).

Inputs : data/interim/mouse_all.h5ad (m01), ../configs/qc.yaml
Outputs (never overwritten unless --force):
  data/interim/mouse_qc.csv          per cell: n_genes, n_counts, pct_mito (mouse 'mt-' prefix), pct_ribo (Rps/Rpl), pct_hb,
                                     fail_mito, fail_min_genes, fail_max_genes, qc_pass
  results/qc/m02_qc_summary.csv      per author label and per sample: n, medians, fraction failing each rule, qc_pass fraction
Run: nice -n 10 python scripts/m02_qc.py [--force]

Rules (same as the human whole-cell cohorts): pct_mito <= 20, n_genes >= 300, n_genes <= per-sample 0.995 quantile.
The deposit is pre-filtered (>= 200 genes, scDblFinder singlets), so no doublet calling is repeated. The granule- and
Purkinje-like clusters have median mito 45-52 % and fail the mito rule en bloc; they are damaged / ambient-rich normal
cerebellar neurons (m04) and are not analysed as tumour, so the human rule is kept rather than relaxed for them.
"""
import sys, argparse
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "src"), str(ROOT.parent / "src")]
import dmg_atlas.utils as U  # noqa: E402
import numpy as np, pandas as pd, scipy.sparse as sp, anndata as ad, yaml, re  # noqa: E402
from mouse_cmp import INTERIM, RESULTS, LOGS, ATLAS_ROOT  # noqa: E402
from mouse_cmp.common import provenance, skip_existing  # noqa: E402

ap = argparse.ArgumentParser(); ap.add_argument("--force", action="store_true"); args = ap.parse_args()
IN, CFG = INTERIM / "mouse_all.h5ad", ATLAS_ROOT / "configs/qc.yaml"
OUT, SUM = INTERIM / "mouse_qc.csv", RESULTS / "qc/m02_qc_summary.csv"
log = U.StepLogger(LOGS / "m02_qc.log", "m02")
if skip_existing([OUT, SUM], args.force, log):
    sys.exit(0)
cfg = yaml.safe_load(CFG.read_text())
TH = cfg["thresholds"]["cells"]; QMAX = cfg["max_genes_quantile"]

with log.step("metrics + flags"):
    a = ad.read_h5ad(IN)
    C = sp.csr_matrix(a.layers["counts"]); v = a.var_names
    lib = np.asarray(C.sum(1)).ravel()
    frac = lambda m: np.asarray(C[:, np.flatnonzero(m)].sum(1)).ravel() / lib * 100
    q = pd.DataFrame({"n_genes": np.diff(C.indptr), "n_counts": lib, "pct_mito": frac(v.str.startswith("mt-")),
                      "pct_ribo": frac(v.str.match(r"^Rp[sl]\d")), "pct_hb": frac(v.str.match(r"^Hb[ab]-"))}, index=a.obs_names)
    q["fail_mito"] = q.pct_mito > TH["pct_mito_max"]
    q["fail_min_genes"] = q.n_genes < TH["min_genes"]
    q["fail_max_genes"] = q.n_genes > q.n_genes.groupby(a.obs["sample"].values).transform(lambda x: x.quantile(QMAX))
    q["qc_pass"] = ~(q.fail_mito | q.fail_min_genes | q.fail_max_genes)
    q.to_csv(OUT)

with log.step("summary"):
    rows = []
    for col in ("author_label", "sample"):
        g = q.groupby(a.obs[col].astype(str).values)
        s = g.agg(n=("qc_pass", "size"), median_genes=("n_genes", "median"), median_counts=("n_counts", "median"),
                  median_pct_mito=("pct_mito", "median"), frac_fail_mito=("fail_mito", "mean"), frac_fail_min_genes=("fail_min_genes", "mean"),
                  frac_fail_max_genes=("fail_max_genes", "mean"), frac_qc_pass=("qc_pass", "mean"))
        rows.append(s.rename_axis("level").reset_index().assign(grouping=col))
    S = pd.concat(rows)[["grouping", "level", *s.columns]]
    S.to_csv(SUM, index=False)
    provenance([OUT, SUM], [IN, CFG], "scripts/m02_qc.py", {"thresholds": TH, "max_genes_quantile": QMAX})
    log.write(f"qc_pass {q.qc_pass.sum()}/{len(q)}\n" + S.round(3).to_string(index=False))
