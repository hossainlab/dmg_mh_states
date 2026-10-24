#!/usr/bin/env python
"""m03 — one-to-one mouse->human orthologs for the deposit genes, mouse gene positions, and the human-symbol counterpart object.

Inputs : data/interim/mouse_all.h5ad (m01); ../data/interim/gene_universe.csv (human atlas universe: GRCh38-2020-A symbol <->
         Ensembl ID); Ensembl 116 BioMart (jun2026.archive.ensembl.org; legacy archives <= 115 for mm10 are retired), raw
         query results cached under data/external/ensembl116_*.tsv (re-used unless --refetch).
Outputs (never overwritten unless --force):
  data/external/mouse_human_orthologs.csv  per deposit gene: mouse symbol, mouse Ensembl ID + how it was resolved
                                           (name / synonym / ambiguous / none), chromosome, start, end, biotype (GRCm39),
                                           human Ensembl ID, orthology type, confidence, human atlas symbol, one2one_in_atlas
  data/processed/mouse_all_hs.h5ad          same cells, var = human atlas symbols of one-to-one orthologs (counts summed if
                                           two deposit symbols resolve to one human gene; none expected), obs from m01
  results/qc/m03_ortholog_summary.csv       genes per resolution step; coverage of human comparator gene sets
                                           (state sets, MP top-50, TME subcluster markers) in the one-to-one space
Run: nice -n 10 python scripts/m03_orthologs.py [--force] [--refetch]

Method: the deposit carries mm10 2020-A symbols (e.g. H3f3a, Hist1h3c) and no Ensembl IDs. Symbol -> mouse gene ID by
current Ensembl name first; otherwise by Ensembl synonym when the synonym points to exactly one gene ID whose current name is
not itself a deposit symbol (rename rescue); otherwise dropped. Orthology: Ensembl Compara 'ortholog_one2one' only. Human gene
ID -> atlas symbol through the atlas universe table (GRCh38-2020-A IDs), so renamed human genes keep their atlas symbol.
Gene positions are GRCm39 (gene order is what inferCNV uses; mm10 vs mm39 order differences are negligible at this scale).
"""
import sys, argparse, io, time
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "src"), str(ROOT.parent / "src")]
import dmg_atlas.utils as U  # noqa: E402
import numpy as np, pandas as pd, scipy.sparse as sp, anndata as ad, requests  # noqa: E402
from mouse_cmp import INTERIM, PROCESSED, EXTERNAL, RESULTS, LOGS, HUMAN, ATLAS_ROOT  # noqa: E402
from mouse_cmp.common import provenance, skip_existing  # noqa: E402
from dmg_atlas.states import STATE_GENE_SETS, TME_SUBTYPE_MARKERS  # noqa: E402

ap = argparse.ArgumentParser(); ap.add_argument("--force", action="store_true"); ap.add_argument("--refetch", action="store_true")
args = ap.parse_args()
IN, UNI = INTERIM / "mouse_all.h5ad", ATLAS_ROOT / "data/interim/gene_universe.csv"
OUT_T, OUT_H, OUT_S = EXTERNAL / "mouse_human_orthologs.csv", PROCESSED / "mouse_all_hs.h5ad", RESULTS / "qc/m03_ortholog_summary.csv"
log = U.StepLogger(LOGS / "m03_orthologs.log", "m03")
if skip_existing([OUT_T, OUT_H, OUT_S], args.force, log):
    sys.exit(0)
MART = "https://jun2026.archive.ensembl.org/biomart/martservice"
QUERIES = {
    "genes": ["ensembl_gene_id", "external_gene_name", "chromosome_name", "start_position", "end_position", "gene_biotype"],
    "synonyms": ["ensembl_gene_id", "external_synonym"],
    "homologs": ["ensembl_gene_id", "hsapiens_homolog_ensembl_gene", "hsapiens_homolog_orthology_type",
                 "hsapiens_homolog_orthology_confidence"],
}


def biomart(name, attrs):
    path = EXTERNAL / f"ensembl116_mmusculus_{name}.tsv"
    if path.exists() and not args.refetch:
        return pd.read_csv(path, sep="\t", dtype=str), path
    xml = ('<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query><Query virtualSchemaName="default" formatter="TSV" header="0" '
           'uniqueRows="1" datasetConfigVersion="0.6"><Dataset name="mmusculus_gene_ensembl" interface="default">'
           + "".join(f'<Attribute name="{a}"/>' for a in attrs) + "</Dataset></Query>")
    for attempt in range(4):
        r = requests.post(MART, data={"query": xml}, timeout=900)
        if r.ok and "[success]" not in r.text[:200] and not r.text.lstrip().startswith("<"):
            break
        time.sleep(30 * (attempt + 1))
    else:
        raise RuntimeError(f"BioMart {name}: {r.status_code} {r.text[:300]}")
    df = pd.read_csv(io.StringIO(r.text), sep="\t", header=None, names=attrs, dtype=str)
    EXTERNAL.mkdir(parents=True, exist_ok=True); df.to_csv(path, sep="\t", index=False)
    return df, path


with log.step("BioMart"):
    q = {k: biomart(k, v) for k, v in QUERIES.items()}
    genes, syn, hom = (q[k][0] for k in QUERIES)
    log.write(f"genes {len(genes)}, synonyms {len(syn)}, homolog rows {len(hom)}")

with log.step("resolve symbols -> orthologs"):
    a = ad.read_h5ad(IN)
    dep = pd.Index(a.var_names)
    by_name = genes.groupby("external_gene_name").ensembl_gene_id.agg(list)
    syn = syn.dropna()
    by_syn = syn.groupby("external_synonym").ensembl_gene_id.agg(lambda x: sorted(set(x)))
    gname = genes.set_index("ensembl_gene_id").external_gene_name
    rows = []
    for s in dep:
        ids = by_name.get(s, [])
        if len(ids) == 1:
            rows.append((s, ids[0], "name")); continue
        if len(ids) > 1:
            rows.append((s, None, "ambiguous_name")); continue
        cand = [i for i in by_syn.get(s, []) if gname.get(i) not in set(dep)]
        rows.append((s, cand[0], "synonym") if len(cand) == 1 else (s, None, "ambiguous_synonym" if cand else "none"))
    T = pd.DataFrame(rows, columns=["mouse_symbol", "mouse_gene_id", "resolved_by"])
    T = T.merge(genes.rename(columns={"ensembl_gene_id": "mouse_gene_id", "external_gene_name": "mouse_name_e116",
                                      "chromosome_name": "chromosome", "start_position": "start", "end_position": "end",
                                      "gene_biotype": "biotype"}), on="mouse_gene_id", how="left")
    # a gene ID claimed by two deposit symbols (e.g. old+new name both present) is ambiguous: keep the exact-name match
    dupid = T.mouse_gene_id.notna() & T.duplicated("mouse_gene_id", keep=False)
    T.loc[dupid & (T.resolved_by == "synonym"), ["mouse_gene_id", "resolved_by"]] = [None, "synonym_collides"]
    h1 = hom[hom.hsapiens_homolog_orthology_type == "ortholog_one2one"].drop_duplicates("ensembl_gene_id")
    T = T.merge(h1.rename(columns={"ensembl_gene_id": "mouse_gene_id", "hsapiens_homolog_ensembl_gene": "human_gene_id",
                                   "hsapiens_homolog_orthology_type": "orthology_type",
                                   "hsapiens_homolog_orthology_confidence": "orthology_confidence"}), on="mouse_gene_id", how="left")
    uni = pd.read_csv(UNI)
    T["human_symbol"] = T.human_gene_id.map(dict(zip(uni.ensembl_id, uni.canonical_symbol)))
    T["one2one_in_atlas"] = T.human_symbol.notna()
    # two mouse genes -> one human gene cannot happen for one2one; guard anyway
    assert not T[T.one2one_in_atlas].human_symbol.duplicated().any()
    T.to_csv(OUT_T, index=False)
    log.write(T.resolved_by.value_counts().to_dict().__repr__() + f"; one2one in atlas universe: {T.one2one_in_atlas.sum()}")

with log.step("human-symbol object + coverage"):
    keep = T[T.one2one_in_atlas]
    idx = a.var_names.get_indexer(keep.mouse_symbol)
    C = a.layers["counts"][:, idx].tocsr()
    h = ad.AnnData(X=C, obs=a.obs.copy(), var=pd.DataFrame({"mouse_symbol": keep.mouse_symbol.values}, index=keep.human_symbol.values),
                   obsm=dict(a.obsm))
    h.layers["counts"] = C
    PROCESSED.mkdir(parents=True, exist_ok=True); h.write_h5ad(OUT_H, compression="gzip")
    hs = set(keep.human_symbol); atlas = set(uni.canonical_symbol)
    cov = [("deposit_genes", len(T), None), *[(f"resolved_{k}", v, None) for k, v in T.resolved_by.value_counts().items()],
           ("one2one_ortholog_any", int(T.human_gene_id.notna().sum()), None), ("one2one_in_atlas_universe", len(keep), None)]
    def setcov(name, genes):
        g = [x for x in dict.fromkeys(genes) if x in atlas]
        cov.append((name, sum(x in hs for x in g), len(g)))
    for k, v in STATE_GENE_SETS.items():
        setcov(f"state:{k}", v)
    mpg = pd.read_csv(HUMAN / "metaprogram_genes.csv")
    for k, v in mpg.groupby("metaprogram").gene:
        setcov(f"metaprogram:{k}", list(v))
    for comp, d in TME_SUBTYPE_MARKERS.items():
        for k, v in d.items():
            setcov(f"tme_marker:{comp}:{k}", v)
    S = pd.DataFrame(cov, columns=["item", "n_covered", "n_in_atlas"])
    S["frac"] = S.n_covered / S.n_in_atlas
    S.to_csv(OUT_S, index=False)
    provenance([OUT_T, OUT_H, OUT_S], [IN, UNI, HUMAN / "metaprogram_genes.csv"] + [p for _, p in q.values()], "scripts/m03_orthologs.py",
               {"biomart": MART})
    log.write("coverage:\n" + S.to_string(index=False))
