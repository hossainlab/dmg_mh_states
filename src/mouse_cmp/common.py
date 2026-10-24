"""Shared helpers for the mouse_comparison scripts: input checksums, provenance sidecars, skip-unless-force, sample design.

Every results/ file gets a sidecar `<file>.provenance.json` with the producing script, UTC time and the sha256 of each input
(the main repo records input checksums per results file; the sidecar is the mouse sub-project's form of that rule).
Checksums of large inputs are cached in data/interim/.sha256_cache.json keyed by (path, size, mtime).
"""
from __future__ import annotations

import datetime, hashlib, json, sys
from pathlib import Path

import numpy as np
import pandas as pd

from mouse_cmp import ROOT, INTERIM

ARMS = ["Control", "FUS", "RT", "FUS+RT"]          # labels as in obs['Group']; reference arm first
_CACHE = INTERIM / ".sha256_cache.json"


def sha256(path: Path | str) -> str:
    p = Path(path).resolve()
    st = p.stat()
    key = f"{p}|{st.st_size}|{int(st.st_mtime)}"
    cache = json.loads(_CACHE.read_text()) if _CACHE.exists() else {}
    if key in cache:
        return cache[key]
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 24), b""):
            h.update(b)
    cache[key] = h.hexdigest()
    _CACHE.parent.mkdir(parents=True, exist_ok=True)
    _CACHE.write_text(json.dumps(cache, indent=1))
    return cache[key]


def provenance(outputs, inputs, script: str, extra: dict | None = None):
    """Write <output>.provenance.json next to every output file (inputs hashed once)."""
    sums = {str(Path(i).resolve().relative_to(ROOT.parent)) if str(Path(i).resolve()).startswith(str(ROOT.parent)) else str(i): sha256(i)
            for i in inputs if Path(i).is_file()}
    rec = {"script": script, "utc": datetime.datetime.utcnow().isoformat(timespec="seconds") + "Z", "inputs_sha256": sums,
           "argv": sys.argv[1:], **(extra or {})}
    for o in outputs:
        Path(str(o) + ".provenance.json").write_text(json.dumps(rec, indent=1))


def skip_existing(outs, force: bool, log=None) -> bool:
    if not force and all(Path(o).exists() for o in outs):
        msg = "all outputs exist; skipping (use --force to rebuild)"
        (log.write(msg) if log else print(msg))
        return True
    return False


def sample_block(sample: pd.Series) -> pd.Series:
    """Processing block inferred from sample numbering: Samples 1-4, 5-8, 9-12, 13-16 each hold one mouse per arm
    (Sample3 absent). Used only as an adjustment/sensitivity covariate; it is an inference, not recorded metadata."""
    n = sample.astype(str).str.replace("Sample", "", regex=False).astype(int)
    return ("B" + ((n - 1) // 4 + 1).astype(str)).astype("category")


def bh(p) -> np.ndarray:
    p = np.asarray(p, float); out = np.full_like(p, np.nan); ok = ~np.isnan(p)
    if ok.sum():
        from statsmodels.stats.multitest import multipletests
        out[ok] = multipletests(p[ok], method="fdr_bh")[1]
    return out
