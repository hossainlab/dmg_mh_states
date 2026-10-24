# Mouse DMG model vs human DMG atlas — m01–m08 notes (2026-10-22)

Deposit: `data/mouse/DIPG1_allcell_raw_counts.h5ad` (sha256 `fb7b0f0a…c244715`) + `DIPG1_allcell_metadata.csv`
(`86f17a3d…de5fce61a`). 67,313 cells × 15,304 mm10 genes, raw UMI counts, pre-filtered (≥ 200 genes, scDblFinder singlets).
15 mice: Control 4 (S1, S5, S9, S13), FUS 4 (S2, S6, S10, S14), RT 3 (S7, S11, S15), FUS+RT 4 (S4, S8, S12, S16).
**No immune-checkpoint (ICI) arm and no Sample3 are in the deposit**; every treatment result below covers FUS / RT / FUS+RT only.

## 1. Which cells are tumour (m04) — changes the starting point

The author labels call 46,532 cells (69 %) tumour, but:

| author label | n | programme | median mito % | CNV clone-positive |
|---|---|---|---|---|
| Neuronal-like tumor | 29,915 | cerebellar granule cell (Cbln3 42 %, Gabra6, Neurod1, Zic1) | 52 | 2.2 % |
| Purkinje-like tumor | 12,999 | Purkinje cell (Pcp2 87 %, Pvalb 85 %, Car8, Calb1) | 45 | 1.4 % |
| Astrocyte-like tumor cells | 372 | astrocyte (Aqp4, Gja1, Slc1a3) | 12 | 2.2 % |
| OPC-like tumor cells | 3,246 | OPC/OC + cycling | 12 | **40 %** |
| immune / vascular / epithelial references | — | — | — | 0.2–1.8 % (null = 1 %) |

CNV (best-match inferCNV against the immune/vascular/epithelial labels, mouse gene order) finds one shared clone:
+chr1, +chr8, +chr11, +chr15/16, +chr17, −chr18, −chr7/−chr12. It appears in the seed cells of 13/14 mice (between-mouse
profile r median 0.63, up to 0.92), as expected for an implanted clonal line. Correlation of each cell with the leave-one-mouse-out
clone profile separates the OPC-like cluster (40 % positive) from every other label, which sits at the 1 % false-positive rate.
**The granule- and Purkinje-like "tumour" clusters are normal cerebellar neurons**, damaged and ambient-rich (98 % / 90 % fail the
human whole-cell mito ≤ 20 % rule). The sampled tissue includes cerebellum; there is no oligodendrocyte cluster at all.

Tumour call: clone-positive OPC-/astrocyte-like cells, plus subclusters of those labels with ≥ 50 % clone-positive cells
(s1 OPC/OC + cycling, s2, s12 interferon-high). **1,479 tumour cells, 891 QC-pass** (the primary set). They are very unevenly
spread: Sample6 holds 54 %, Sample9 has none, and most mice have 16–123.

## 2. Recapitulation of human malignant states and programmes (m05)

**Correlation of whole-cluster profiles** (z-scored pseudobulk, 2,000 one-to-one genes most variable across human classes).
Positive controls work: T ↔ human T_NK ρ 0.68, Mono/Macro and microglia ↔ myeloid 0.63–0.64, vascular 0.55, granule neurons ↔ neuron 0.48.

| mouse tumour cluster | best human classes (ρ) | AC-like | MES-like | neuron / astro |
|---|---|---|---|---|
| s1 (OPC/OC, cycling) | cycling 0.52, OC-like 0.50, OPC-like 0.37 | 0.14 | 0.14 | 0.11 / −0.31 |
| s2 | OPC-like 0.47, cycling 0.47, OC-like 0.37 | 0.13 | 0.03 | −0.09 / −0.31 |
| clone-only cells | cycling 0.45, OC 0.42, OPC 0.40 | 0.20 | 0.13 | 0.13 / −0.23 |
| s12 (IFN/MHC-I high) | MES-like 0.25, hypoxia 0.23 (weak); cycling/OC/OPC collapse to 0.12–0.13, so s12 is excluded from the 0.37–0.52 lineage range quoted for the other three clusters | 0.16 | 0.25 | −0.09 / −0.03 |

**Human metaprograms in mouse cells (UCell).**
- Enriched in tumour: cycling MP1 and MP11 (1.4–1.6× the best normal group), OC-like MP5 (1.1×).
- Higher in normal mouse astrocytes than in tumour, 5–10×: every AC-like MP (MP7, MP10, MP14, MP15).
- OPC-like MP3/MP6 are not tumour-specific. They are the nuclei-dominated, long-gene programmes of the human atlas, and this is a whole-cell assay.

**Mouse de novo cNMF** on the pooled primary set (k = 9), with hypergeometric overlap against the human MPs and Liu 2022 sets:

| programme | genes | best human match |
|---|---|---|
| p1 | Hmgb2, Top2a, Cdk1, Birc5 | **MP1 cycling**: 25/50 genes, Jaccard 0.37, q 6e-50 (also Liu G2M, S). Sample6 alone 21/50 (J 0.29); pooled without Sample6 only 8/50 (J 0.09, q 1e-10) — present but much weaker |
| p2 | Olig1, Fabp7, Bcas1, Cdkn2a, Stmn3 | Liu OPC-like-1, weak (4 genes) — OPC/OC-lineage |
| p3, p4 | Isg15, Ifit, Rsad2; Psmb8/9, B2m, Cxcl10 | interferon / MHC-I. No human malignant counterpart: best hits are the technical IEG MP13 and Liu MES-like (5 genes). Carried by true tumour cells (s12, Sample11) |
| p7 | — | AC-like MP7 / MP10, but in only 8 cells, mostly in the normal-astrocyte subcluster: astrocyte contamination, not a tumour AC-like state |
| p5, p6, p8 | — | technical: neuronal ambient, nuclear-transcript, mitochondrial |
| p9 | — | ER/UPR (Hspa5, Pdia6, Calr) |

The human state gene sets with raw argmax (mean over mice, primary set): cycling 32 % (human 6–11 %), OPC-like 26 %
(human 37–53 %), AC-like 24 %, MES-like 11 %, OC-like 5.5 %. The AC-like share comes from broad genes (Apoe, Clu, Cst3, Sparc)
in a tumour-only background and is contradicted by the MP and correlation results; treat state proportions as relative.

**Bottom line:** the model reproduces the human OPC/OC-lineage, highly proliferative malignant compartment (more cycling than
human). It does not reproduce the human AC-like or hypoxic/MES-like axes. Its one extra axis is an interferon/antigen-presentation
programme with no human malignant equivalent. The recapitulation claim is limited to ~900 QC-pass tumour cells, half from one mouse.

## 3. Microenvironment (m07)

Cells are QC-pass and called normal. Clustering is per compartment (Harmony by mouse in place of the human scVI latent), with the human marker-set labelling rule.
- **Myeloid** (11,621):
  - lipid-phagocytic TAM (Apoe, Trem2, C1q): 57 %, human 53 %;
  - homeostatic microglia (P2ry12, Tmem119): 20 %, human 0 %;
  - interferon-high, Ccr2⁺ Plac8⁺ monocyte-derived cells: 15 %, human 0 %;
  - human BAM (23 %) and proliferating myeloid (7 %) are absent.
- **T/NK** (6,378) is far more abundant than in human: 19 % of QC-pass normal cells vs 2–10 %.
  - NK cells are 34 % of T/NK, vs 15 % in human.
  - Tregs are 8 %; no human cluster carries that label.
  - A Lag3⁺/Pdcd1⁺ T cluster (1,208 cells, 7 mice) is the counterpart of the human exhausted-T label. It narrowly lost to CD8_T under the human margin rule.
- **Cross-species marker overlap** is modest (1–9 of the top-25 markers). Matches are cell-type-level: T_NK_1 ↔ human CD4_T, T_NK_6 ↔ proliferating T/NK, mouse capillary EC ↔ human capillary EC.
- **Composition differs grossly from human tumours:** no oligodendrocytes, very few astrocytes, a large cerebellar neuron fraction, and immune shares that swing 5–77 % between mice.

## 4. Treatment effects (m08) — exploratory

Unit = mouse; tumour analyses keep the 11 mice with ≥ 20 primary tumour cells (Control 3, FUS 4, RT 2, FUS+RT 2).
Composition models are quasi-binomial (φ up to 150, i.e. very overdispersed). Pseudobulk DE is limma-trend.

- **One composition row reaches q < 0.1 under the primary model** — OC-like tumour-state share, FUS+RT vs RT, q = 0.0986 (coefficient 2.46 log-odds, p = 0.016, 2 vs 2 mice, no overdispersion estimated) — and nothing else in any composition set does, but refitting the same 160 primary contrasts with the inferred processing block gives 9 rows at q < 0.1 (see the sensitivity bullet), so the flatness of composition is model-dependent; pseudobulk DE finds essentially no genes in microglia, Mono/Macro, T or NK. The tumour population has 1–8 genes (Vim up; Ptprz1, Tnr, Igsf9b down with FUS+RT).
- **Most consistent signal: RT moves tumour cells toward hypoxia/MES-like.**
  - Human MP4 hypoxia UCell: RT q = 0.006, FUS+RT p = 0.009.
  - MP8 MES-like: RT q = 0.04.
  - Human MES-like state score: RT q = 0.054.
  - OPC-like state proportion falls with RT (p = 0.03).
  - It rests on 2 RT mice vs 3 Controls, and the two RT mice also have higher mito.
- **With the inferred processing block as a covariate:** regulatory T cells rise with every arm (FUS q 0.047, RT q 0.03); monocytes rise with FUS (q 0.02); pDC/DC fall. Not significant in the primary model, so report only as sensitivity.
- **FUS:** no robust effect. The tumour fraction is higher (p = 0.06), driven by Sample6.

## 4b. Projection onto the human reference (m06)

The step listed below as "not run" was written and run on 2026-10-22. It asks the human reference model itself to label every
mouse cell, zero-shot: `SCANVI.prepare_query_anndata` aligns the mouse matrix to the reference's 3,000 HVG symbols
(**2,363 present, 637 zero-filled**, i.e. 78.8 % — scvi-tools warns below 80 %) and the reference weights are used unchanged.
That identity is verified tensor-by-tensor against `model.pt`: 34 encoder/classifier tensors, **max |difference| = 0.0**.

- **The CNV-based tumour calls of m04 are corroborated by a method that shares none of their assumptions.** Of the 891 QC-pass
  m04 tumour cells, **97.5 % are predicted malignant** by the human classifier, at median posterior 0.9986 and 98.9 %
  classifier/kNN agreement; of the 21,859 QC-pass cells called normal, 9.9 % are. The author's `OPC-like tumor cells` label is
  82.5 % malignant; `Astrocyte-like tumor cells` is 16.7 % malignant and 83.3 % astrocyte.
- **The granule- and Purkinje-like "tumour" clusters again do not behave like tumour.** Neuronal-like tumor (29,915 cells) maps
  to astrocyte 67.8 % / malignant 22.0 % at median posterior 0.732; Purkinje-like tumor (12,999) to astrocyte 75.1 % /
  malignant 24.5 % at median posterior 0.895. Their classifier/kNN agreement (0.829 and 0.737) is the lowest of any
  label above 1,000 cells; across all 14 labels only Neutrophils is lower (250 cells, 0.680, median posterior 0.541),
  the human atlas having no neutrophil class. The human reference has no cerebellar
  granule or Purkinje class, so these cells have nowhere correct to go — their low posterior is the signal, not their label.
- **Positive controls transfer as expected:** Microglia → myeloid 95.2 %, Mono/Macro → myeloid 91.8 %, T cells → T_NK 93.0 %,
  NK cells → T_NK 86.6 %, Vascular cells → vascular 70.0 %. pDC → myeloid 95.7 % and B/Plasma → myeloid 94.9 % are reference
  gaps (no dendritic or B-cell class), not errors.
- **Calibration against the five human query cohorts that used this identical code path** (`m06_calibration.csv`) is the result
  that licenses the comparison: median `uncertainty_knn` — the percentile of a cell's mean 15-NN distance within the
  reference's own distribution — is **0.795 for mouse tumour cells**, against 0.998 (Filbin), 0.992 (Liu), 0.977 (Braun),
  0.920 (LaBelle) and 0.574 (Siletti). Mouse DMG cells therefore sit *closer* to the human manifold than all three published human
  DMG cohorts (Filbin, Liu, LaBelle) and than the fetal-brain reference (Braun); only the adult-brain reference
  (Siletti, 0.574) is closer. The statistic tracks assay and platform more than species: Siletti (adult human brain, snRNA) is closest of
  all, because its neurons and oligodendrocytes land in dense normal-brain regions of the reference.
- **What the projection cannot do:** 31,354 of 67,313 cells (46.6 %) are assigned `astrocyte`, mostly the cerebellar neurons
  above. Roughly half the mouse object has no counterpart class in a reference built from human pons and thalamus, and a
  classifier always returns a label. Read the posterior and `uncertainty_knn` columns, never the label alone.

Conclusion: the projection and the m05 cluster-correlation agree, by independent routes — the mouse model reproduces the
human OPC/OC-lineage proliferative malignant compartment, and the comparison is as well-calibrated for these mouse cells as
for the human query cohorts already reported in the atlas paper.

## 6. Human subpopulations, pathways, TFs and treatment programmes (m06a, m09, m10; 23 Oct 2026)

**m06a — human annotation first.** The human reference panel (`m06a_human_subpopulations.csv`), from the 167,012
malignant atlas cells: OPC-like 84,979 cells (50.9 % under the raw-score rule; 25.6 % under the within-assay-class z
rule, 53 samples, MP3/MP6), AC-like 44,693 (26.8 % / 20.7 %, MP7/MP10/MP14/MP15), MES-like 25,960 (15.5 % / 19.1 %,
MP8), cycling 9,766 (5.9 % / 11.7 %, MP1/MP11), OC-like 1,105 (0.7 % / 17.9 %, MP5). `m06a_joint_embedding.png` shows
this annotation and places the mouse cells on the same axes (kNN placement, display only).

**m09 — four lines of evidence per subpopulation.** Human CollecTRI (1,185 TFs) and PROGENy (14 pathways) applied to
both species on the ortholog-mapped matrices, 12,687 shared genes.
- Pseudobulk correlation, mouse tumour vs each human subpopulation: cycling **+0.47**, OPC-like **+0.35**, OC-like
  −0.13, AC-like **−0.36**, MES-like **−0.40**.
- Marker detection (median fraction of cells, mouse tumour / human cells of that subpopulation): OPC-like 0.146/0.313,
  cycling 0.118/0.506, AC-like 0.084/0.397, OC-like 0.064/0.120, MES-like 0.040/0.292.
- TF activity: OLIG2 2.76 in mouse (human OPC-like 0.14, cycling 2.43), SOX10 1.08 (0.23), E2F4 5.72 (6.97), MYC 12.81
  (8.65), E2F1 1.67, FOXM1 0.76; HIF1A −1.2 against human MES-like +3.0; ASCL1 −0.27 and SOX2 −0.38 against 0.49–1.89
  and 1.32–2.04 across human subpopulations.
- Pathways: hypoxia −14.4 in mouse tumour against +9.4 in human MES-like; JAK-STAT +5.5 against −9.0 to −5.1 in every
  human subpopulation; STAT1 TF +1.6 against −3.3 (OPC-like) to +0.5 (MES-like).
Verdict: cycling and OPC-like present, AC-like and MES-like absent, OC-like equivocal, plus a mouse-only
interferon/JAK-STAT axis.

**m10 — treatment effects on pathways, TFs and immune programmes. Nothing survives correction.**
- Tumour pathway + TF activity, 11 mice: 222 tests, lowest q 0.228. Best nominal: Androgen pathway up with RT
  (+1.65, p 0.005), NR3C1 TF up with FUS+RT (+1.86, p 0.004).
- **The m08 hypoxia finding does not replicate under an independent method**: PROGENy hypoxia in the same tumour cells
  of the same mice gives +0.31, p 0.877 for RT vs control, against MP4 hypoxia q 0.006 in m08. Report the RT to
  hypoxia/MES shift as a single-method hypothesis.
- Immune programmes (AUCell, mouse nomenclature, myeloid 15 mice / T-NK 9 mice): 36 tests, lowest q 0.459. Myeloid
  immune suppression rises in all three arms (FUS+RT +0.042 p 0.026, FUS +0.035 p 0.058, RT +0.022 p 0.235, control mean
  0.164); T-cell activation rises with RT (+0.025 on 0.218, p 0.101); T/NK inflammation falls with FUS+RT (−0.018 on
  0.036, p 0.109).
- Gene level: 82,632 tests, lowest q 0.568. Nominal programme-gene movers: `Il1a` down with RT in myeloid (p 0.013),
  `Trem2` up with FUS+RT in T/NK (p 0.027), `Nt5e` up with FUS (p 0.041), `Mrc1` down with FUS+RT (p 0.049), `Nfkb1`
  down with FUS+RT (p 0.051), `Il1b` down with FUS+RT (p 0.079).
- **Coverage caveat (`m10_programme_coverage.csv`): the deposit's 15,304 genes do not include `Ifng`, `Gzma`, `Gzmb`,
  `Cd69`, `Tbx21`, `Icos`, `Klrg1`, `Foxp3`, `Ctla4`.** T-cell activation keeps 7 of 14 genes and regulatory-T 3 of 6,
  so those two programmes are under-measured; Tregs cannot be assessed by Foxp3 at all.

## 5. Caveats

- Treatment effects are underpowered (2–4 mice per arm; tumour cells concentrated in one FUS mouse). Use them as hypotheses.
- Tissue sampling differs between mice (immune 5–77 %, cerebellar neurons up to ~80 %). It drives compartment-level composition more than treatment does.
- The processing block (Samples 1–4, 5–8, …) is inferred from numbering, not recorded.
- Only one-to-one orthologs are used (12,759 genes). Human MP coverage is 82–98 %, except the ciliary MP9 at 44 %.
- m06 (zero-shot projection onto the human scANVI reference) was run on 2026-10-22; see section 4b. Only 78.8 % of the reference HVGs have one-to-one orthologs, so 637 genes enter the encoder as zeros; the m05 correlation remains the assumption-light comparison and the two agree.
- **ICI could not be assessed.** Its data need to be deposited before it can be analysed.
