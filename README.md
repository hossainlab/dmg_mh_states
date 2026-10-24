# mouse_comparison — mouse DMG model vs the human DMG atlas

Self-contained sub-project for two questions:

1. **Recapitulation** — how well do mouse DMG tumour cells reproduce the human malignant states (OPC-/AC-/OC-/MES-like,
   cycling) and the replicated cNMF metaprograms of the human atlas, and does the mouse microenvironment resemble the
   human TME subclusters?
2. **Treatment effects** — how do focused ultrasound (FUS), radiotherapy (RT), FUS+RT and immune checkpoint inhibition
   (ICI) shift tumour-cell states and TME composition relative to untreated control?

Everything here follows the main repo conventions (`../CLAUDE.md`): numbered CLI scripts, `StepLogger` timing/RSS logs,
inputs' checksums recorded with every `results/` file, earlier outputs never overwritten without `--force`.
Human atlas outputs are never modified from here.

## Layout

```
mouse_comparison/
  README.md            this file
  configs/mouse.yaml   mouse object path, metadata column names, treatment arms, ortholog policy, human comparators
  src/mouse_cmp/       path constants (ROOT, RAW, INTERIM, PROCESSED, HUMAN, EXTERNAL, RESULTS, CONFIGS, LOGS)
  scripts/             numbered entry points, m01_... (see plan below)
  docs/                analysis notes per step
  data/mouse/          mouse deposit as received (h5ad + meta CSV) (not tracked)
  data/interim/        ingested h5ad, QC flags              (not tracked)
  data/processed/      ortholog-mapped / annotated objects  (not tracked)
  data/external/       ortholog table, mouse gene sets      (not tracked)
  data/human/          symlinks to human atlas outputs      (not tracked)
  results/tables|figures|qc
  logs/                                                     (not tracked)
```

## Inputs

Mouse: `data/mouse/DIPG1_allcell_raw_counts.h5ad` (h5ad export of the standardized Seurat object, raw counts;
sha256 `fb7b0f0a60993e98d9e2c3295bba13c4808d4cb0f4f0ef4734fda53de5fce61a`) and `data/mouse/DIPG1_allcell_metadata.csv`
(meta.data, identical to obs; sha256 `86f17a3d9264fae610e59a4381deb59089e19ea9b98c72751653c8776c244715`). Replaces the
truncated RDS. 15 mice: Control 4, FUS 4, RT 3, FUS+RT 4. **No ICI arm and no Sample3 in the deposit.**

Human (symlinks in `data/human/`, pointing into `../data/processed`, `../results/tables`, `../models`):

| link | used for |
|---|---|
| `dmg_atlas_v1.h5ad`, `atlas_core.h5ad` | reference expression, labels, embeddings |
| `malignant_calls.csv`, `malignant_state_assignment.csv.gz` | malignant selection, per-cell state scores |
| `metaprograms.csv`, `metaprogram_genes.csv`, `metaprogram_members.csv`, `replication_matrix.csv`, `metaprogram_scores.parquet` | metaprogram gene sets and replication status |
| `tme_subclusters.csv`, `tme_subcluster_markers.csv`, `tme_subcluster_marker_scores.csv`, `tme_subcluster_labels_clean.csv.gz` | TME subcluster comparison |
| `phase6_curated_gene_sets.csv`, `reference_gene_sets_scripted.csv` | Filbin/Liu state and curated reference sets |
| `c2l_signatures_mp.csv` | per-class reference signatures |
| `scanvi`, `scvi_malignant` | optional projection of ortholog-mapped mouse cells |
| `sample_manifest.csv` | human treatment / timing metadata |

## Steps (run from mouse_comparison/: `nice -n 10 python scripts/<script> [--force]`)

| script | status | what |
|---|---|---|
| `m01_ingest.py` | done | integrity checks, design columns (sample, arm, inferred block), author labels verbatim |
| `m02_qc.py` | done | human whole-cell QC rule (mito <= 20 %, genes >= 300), flag only |
| `m03_orthologs.py` | done | Ensembl 116 one-to-one orthologs -> human atlas symbols (12,759 genes); mouse gene positions |
| `m04_tumour_calls.py` | done | inferred CNV + shared-clone score + subclusters -> tumour calls (1,479; 891 QC-pass) |
| `m05_state_comparison.py` | done | human states / MPs on mouse tumour cells, cross-species pseudobulk correlation, mouse cNMF vs human MPs |
| `m06_reference_mapping.py` | done | zero-shot projection onto the human scANVI reference (weights verified unchanged); per-cell lineage, kNN transfer, uncertainty calibrated against the five human query cohorts |
| `m07_tme_comparison.py` | done | myeloid / T-NK / vascular subclusters with the human labelling rule; proportions vs human |
| `m08_treatment_effects.py` | done | mouse-level composition (quasi-binomial), scores (OLS), pseudobulk DE (limma-trend) |

Findings and caveats: `docs/m01-m08_notes.md`. Every `results/` file has a `.provenance.json` sidecar with input sha256s.
Shared helpers: `src/mouse_cmp/common.py`.
