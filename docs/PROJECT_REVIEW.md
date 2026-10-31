# Project review - `dmg_mh_states`

Full review, 1 November 2026, against the state at commit `797b023`. Scope: pipeline integrity, statistical method,
reproducibility, figures, documentation and release readiness. Every number below was read from the repository or
recomputed, not taken from prose.

## 1. Release blocker

**Deposit-derived per-cell data is already public.** `data/processed/mouse_tme_labels.csv` (2.0 MB, 18,303 rows:
cell barcode, compartment, subcluster, label, UMAP coordinates) is in commit `2d24e4d`, which is an ancestor of
`origin/main` on **github.com/hossainlab/dmg_mh_states**, a **public** repository. Three provenance sidecars for
mouse-derived files are published with it.

The data belongs to Shenglin and Dr Mei and is unpublished. `.gitignore` correctly excludes the raw `.h5ad`, so the
count matrix is not exposed, but a derived per-cell annotation table is. This is a decision for the user and the
depositors, not a change to make unilaterally: purging it requires rewriting published history and force-pushing.
Making the repository private is the reversible first step.

## 2. What is sound

- **Reproducible.** A full `--force` recompute of `m04`-`m12` takes about 2.5 minutes and reproduces **every results
  table byte-for-byte**; only gzip internal timestamps and provenance sidecars differ. This was verified twice.
- **Provenance.** 66 sidecars record the producing script and the sha256 of every input. Exactly one output lacks
  one: `results/tables/m13_figure_index.csv`, which is an index of other outputs rather than an analysis result.
- **Unit of replication.** Treatment inference uses the mouse, not the cell, throughout `m08`, `m10` and `m14`.
  Quasi-binomial dispersion is computed by hand as `phi = max(1, Pearson X2 / df_resid)` with a t-test on
  `df_resid`, because `statsmodels`' `scale="X2"` is wrong for two-column endog - it returned 0.24 against a true
  79.8 and produced cell-level p-values.
- **Multiple testing.** Corrected within each family in `m05`, `m07`, `m08`, `m10` and `m14`, and test counts are
  reported alongside the survivors (12 of 245,308 pseudobulk DE tests at q<0.1; 10 of 528 state/MP tests; 12
  composition contrasts with lowest q 0.495).
- **Pre-registered rules, and honest failure.** The clone threshold is fixed before any human comparison. When the
  `m14` resolution rule could not be satisfied at any resolution in 0.2-1.0, the script says so, falls back
  explicitly, and flags the offending cluster rather than merging it away.
- **Cross-method checking.** The RT-hypoxia result (`m08` panel i, q 0.0062) was tested with an independent method
  and did not replicate (PROGENy, +0.31, p 0.877). It is reported as a single-method hypothesis everywhere it
  appears. A project that overturns its own headline is doing the right thing.
- **Colour and legibility.** All 13 figures render at 183 mm with no geometric legibility failures. The palette is
  drawn from published colour-vision-safe sets; worst pair within the malignant subpopulations 0.405, and 0.196
  across the union that shares figure 2a.

## 3. Faults found and fixed in this review

| # | fault | status |
|---|---|---|
| 1 | Figure 1a told the reader "a box left of the dashed line carries no clone" - but **every** label's median sits left of it, including OPC-like (0.218 vs threshold 0.267), the one label that does carry the clone and supplies all 1,471 tumour calls. The line is a per-cell threshold; a median below it only means fewer than half the cells clear it. | Fixed: the panel now prints the fraction of each label's cells clearing the line (OPC-like 40%, all others 0-2%). |
| 2 | Figure 11 drew M5 as a peer of M1-M4 in four panels while excluding it from every statistic. | Fixed: hatched bars, asterisked ticks, "(excluded)" spelled out in panel a, one footnote giving the reason. |
| 3 | All 13 figure links in `REPORT.md` were **broken inside the delivered zip** - a recipient opened the report to broken images. | Fixed: `m15_package.py` rewrites links to the package layout, ships the PNGs, and fails the build if any link does not resolve. |
| 4 | The report embedded only 9 of 13 figures; Figure 10 was cited in the argument but had no section and was never shown. | Fixed: Figure 10 now has its own section, and Figures 11-13 are embedded in the `m14` section. |
| 5 | Figure 11f's 0% bar for M5 drew nothing, reading as missing rather than zero. | Fixed: zero values get a baseline stub. |

## 4. Open weaknesses

These are judgements about the work, not defects to patch.

1. **The treatment question is underpowered and will stay that way.** 15 mice, 3-4 per arm. Four independent
   parameterisations - author labels, transferred subtypes, state scores, de novo clusters - all return nothing, and
   the single tumour-transcriptome result (12 of 245,308 tests) is the only thing that clears correction. The
   honest reading is that this design can exclude large effects and nothing else. The report says so.
2. **No environment pinning.** There is no `environment.yml`, `requirements.txt` or lockfile. The analysis currently
   runs on scanpy 1.11.5, anndata 0.12.19, numpy 2.4.6, pandas 2.3.3, scipy 1.17.1, scikit-learn 1.7.2,
   statsmodels 0.15.0, decoupler 2.2.0, matplotlib 3.11.2. Byte-identical reproduction has only been demonstrated
   on this machine with these versions. This is the cheapest remaining improvement.
3. **Cell-level DE treats cells as independent.** `m11` and `m14` use `scanpy`'s Wilcoxon test on cells, which is
   conventional for marker discovery but inflates significance, because cells within a mouse are not independent.
   It is used only to name clusters, never to support a treatment claim, and the `m11` analysis exists precisely to
   show that the resulting labels do not mark clusters. Worth one sentence in the methods.
4. **`m09_subpopulation_correlation.csv` carries an uncorrected `p` column.** Thirty tests, no BH column. No claim
   in the report rests on it - the correlations are reported as effect sizes - but the column should be labelled or
   corrected before anyone reuses the table.
5. **The Q1 verdict is a convergence judgement, not a test.** "Two of five" rests on four descriptive lines of
   evidence agreeing in direction, with no formal statistic and no confidence interval. The report is explicit
   about this and tells the reader to read the ordering rather than the absolute values. Defensible, and the
   limitation should survive into any manuscript.
6. **One animal dominates the strict tumour set.** Sample6 contributes 45% of the 891 cells and 61% of the
   OPC-like calls. Every per-subtype statement inherits that imbalance; `m11` panels b and f now show it directly.

## 5. Suggested order of work

1. Decide the GitHub exposure (section 1) - everything else can wait behind it.
2. Add a pinned `environment.yml`; it is half an hour and it is what makes the byte-identical reproduction a claim
   anyone else can check.
3. Label or correct the `m09` p column.
4. Add the pseudoreplication sentence to the methods.
