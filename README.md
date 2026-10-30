# Mouse DMG model vs. the human DMG atlas

Compares a mouse diffuse midline glioma (DMG) model with the human H3 K27-altered DMG atlas
(311,618 cells, 57 samples, 6 studies).

**Questions**

1. **Recapitulation** — do mouse tumour cells reproduce the human malignant states (OPC-, AC-, OC-, MES-like, cycling)
   and replicated cNMF metaprograms? Does the mouse microenvironment (TME) resemble the human TME subclusters?
2. **Treatment effects** — how do focused ultrasound (FUS), radiotherapy (RT) and FUS+RT shift tumour-cell states,
   pathways and TME composition relative to untreated control?

📄 **Full results: [`docs/REPORT.md`](docs/REPORT.md)** · step notes: [`docs/m01-m08_notes.md`](docs/m01-m08_notes.md)

---

## Quick start

This is a standalone repository that reads the human atlas from a sibling `dmg_atlas` checkout
(it imports `dmg_atlas/src/dmg_atlas` and reads human outputs from its `data`, `results` and `models`).
Set `DMG_ATLAS_ROOT` if the atlas lives somewhere other than `../dmg_atlas`.

```bash
nice -n 10 python scripts/m01_ingest.py            # run steps in order, m01 → m14
nice -n 10 python scripts/m01_ingest.py --force    # recompute an existing output
```

## Pipeline

| step | script | what it does |
|---|---|---|
| 1 | `m01_ingest.py` | integrity checks; design columns (sample, arm, block); author labels kept verbatim |
| 2 | `m02_qc.py` | human whole-cell QC rule (mito ≤ 20 %, genes ≥ 300) — flags only, no filtering |
| 3 | `m03_orthologs.py` | Ensembl 116 one-to-one orthologs → human atlas symbols (12,759 genes) |
| 4 | `m04_tumour_calls.py` | inferred CNV + shared-clone score → tumour calls (1,479 cells; 891 QC-pass) |
| 5 | `m05_state_comparison.py` | human states / metaprograms scored on mouse tumour cells; pseudobulk correlation; mouse cNMF vs. human MPs |
| 6 | `m06_reference_mapping.py` | zero-shot projection onto the human scANVI reference; lineage transfer, calibrated uncertainty |
| 6a | `m06a_joint_embedding.py` | human annotation and mouse cells shown on the same UMAP axes |
| 7 | `m07_tme_comparison.py` | myeloid / T-NK / vascular subclusters with the human labelling rule; proportions vs. human |
| 8 | `m08_treatment_effects.py` | per-mouse composition (quasi-binomial), state scores (OLS), pseudobulk DE (limma-trend) |
| 9 | `m09_subpopulation_evidence.py` | per subpopulation: markers, PROGENy pathways, CollecTRI TFs, correlation — human and mouse alike |
| 10 | `m10_treatment_programmes.py` | treatment effects on tumour pathways/TFs and named immune programmes (AUCell) |
| 11 | `m11_subtype_umap_de.py` | do the transferred subtype labels mark clusters? silhouette vs. permutation, per-subtype DE |
| 12 | `m12_author_umap_annotated.py` | the authors' own published UMAP annotated with clone, strict tumour set and transferred lineage |
| 13 | `m13_figure_index.py` | per-panel figure index: what each panel shows, how to read it, key numbers, caveats |
| 14 | `m14_denovo_subpopulations.py` | de novo clustering of the authors' OPC-like tumour cluster; cross-species markers/pathways/TFs; treatment |

Each script writes to `results/tables/`, `results/figures/`, `results/figures_denovo/` or `results/qc/`.
Shared helpers: `src/mouse_cmp/common.py` (pipeline), `src/mouse_cmp/figio.py` (figure sizing and legibility
checks), `src/mouse_cmp/palette.py` (the single colour system, colour-vision checked).

## Data

### Mouse input (`data/mouse/`, not tracked)

| file | sha256 |
|---|---|
| `DIPG1_allcell_raw_counts.h5ad` — raw counts, 67,313 cells × 15,304 genes | `fb7b0f0a…fce61a` |
| `DIPG1_allcell_metadata.csv` — Seurat `meta.data` (= `obs`) | `86f17a3d…244715` |

Full hashes in `configs/mouse.yaml`.

**Design:** 15 mice — Control 4, FUS 4, RT 3, FUS+RT 4.

> ⚠️ The deposit has **no immune-checkpoint (ICI) arm** and **no Sample3**. Several T-cell genes
> (`Ifng`, `Gzmb`, `Foxp3`, `Ctla4`, …) are also missing; see the report.

### Human comparators (`data/human/`, read-only symlinks)

| files | used for |
|---|---|
| `dmg_atlas_v1.h5ad`, `atlas_core.h5ad` | reference expression, labels, embeddings |
| `malignant_calls.csv`, `malignant_state_assignment.csv.gz` | malignant selection, per-cell state scores |
| `metaprograms.csv`, `metaprogram_genes.csv`, `metaprogram_members.csv`, `metaprogram_families.csv`, `replication_matrix.csv` | metaprogram gene sets, families, replication status |
| `tme_subclusters.csv`, `tme_subcluster_markers.csv`, `tme_subcluster_marker_scores.csv`, `tme_subcluster_labels*.csv.gz` | TME subcluster comparison |
| `phase6_curated_gene_sets.csv`, `reference_gene_sets_scripted.csv` | Filbin/Liu states, curated reference sets |
| `c2l_signatures_mp.csv`, `composition_per_sample.csv`, `sample_manifest.csv` | class signatures, human composition, sample metadata |
| `scanvi/`, `scvi_malignant/` | trained models for projecting mouse cells |

### External resources (`data/external/`)

Ensembl 116 mouse genes / homologs / synonyms, the derived ortholog table, and OmniPath CollecTRI + PROGENy
(sources in `resources_provenance.json`).

## Layout

```
configs/mouse.yaml     input paths, metadata columns, treatment arms, ortholog policy, human comparators
src/mouse_cmp/         path constants and shared helpers
scripts/               numbered entry points (m01 … m10)
docs/                  REPORT.md (results) and per-step notes
data/                  mouse/ interim/ processed/ external/ human/   (large files not tracked)
results/               tables/ figures/ qc/
logs/                  step logs (not tracked)
```

## Conventions

- Numbered CLI scripts; each step logs time and memory (`StepLogger`) to `logs/`.
- Every `results/` file has a `.provenance.json` sidecar with input sha256 checksums.
- Existing outputs are never overwritten without `--force`.
- Human atlas outputs are never modified from here.
