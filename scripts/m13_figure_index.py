#!/usr/bin/env python
"""m13 - build the figure index: one row per panel, with the mechanical facts read from the files themselves.

The per-panel descriptions are curated (a figure's meaning cannot be derived from its bytes), but everything that can
be checked is read rather than typed: the page size comes from each PDF's MediaBox, the producing script and the input
checksums come from the provenance sidecar written next to every figure, and the file sizes and modification times come
from the filesystem. If a figure is re-rendered, re-running this script keeps the index honest.

Outputs: results/tables/m13_figure_index.csv
Run: python scripts/m13_figure_index.py
"""
import sys, re, json, csv
from pathlib import Path
from datetime import datetime, timezone
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from mouse_cmp import RESULTS  # noqa: E402

OUT = RESULTS / "tables/m13_figure_index.csv"
FIG = RESULTS / "figures"

# figure number -> (file stem, short title, report section)
FIGURES = [
    (1, "m04_tumour_calls", "Which mouse cells carry the tumour clone", "Which mouse cells are really tumour"),
    (2, "m06a_joint_embedding", "Human annotation first, then the mouse on the same axes", "The human reference"),
    (3, "m09_subpopulation_evidence", "Four independent lines of evidence per subpopulation", "Question 1"),
    (4, "m05_state_comparison", "State scores, human metaprograms, and the mouse's own programmes", "Question 1"),
    (5, "m06_reference_projection", "The human reference model applied to the mouse cells", "Question 1"),
    (6, "m07_tme_comparison", "Microenvironment, mouse against patient", "The microenvironment"),
    (7, "m08_treatment_effects", "Treatment effects on composition and state scores", "Question 2"),
    (8, "m10_treatment_programmes", "Treatment effects on pathways, TFs and immune programmes", "Question 2"),
    (9, "m11_subtype_umap_de", "The transferred subtype labels tested on the mouse's own embedding", "Follow-up"),
    (10, "m12_author_umap_annotated", "The authors' own DIPG1 UMAP annotated with this study's evidence", "Follow-up"),
    (11, "m14_denovo_clusters", "De novo clusters inside the authors' OPC-like tumour cluster", "Follow-up"),
    (12, "m14_denovo_vs_human", "Each de novo cluster tested against the human subpopulations", "Follow-up"),
    (13, "m14_denovo_treatment", "Treatment effects on the de novo clusters", "Follow-up"),
]

# Figures 11-13 are written to a separate directory so the de novo analysis stays visibly distinct from the
# label-transfer figures it is meant to be an alternative to.
FIGDIR = {11: RESULTS / "figures_denovo", 12: RESULTS / "figures_denovo", 13: RESULTS / "figures_denovo"}

# figure -> panel -> (panel title, what it shows, how to read it, key numbers, caveat)
PANELS = {
 1: {
  "a": ("CNV clone signal by author label",
        "Distribution of each cell's correlation with the shared tumour CNV profile, computed leave-one-mouse-out, one box per author label.",
        "The dashed line is the clone-positive threshold. A box sitting left of it means that label carries no clone.",
        "OPC-like tumor cells median 0.218; every other label <= 0.069.",
        "Leave-one-mouse-out avoids a cell contributing to the profile it is scored against."),
  "b": ("Seed-cell CNV per mouse",
        "Mean inferred copy number x100 per chromosome for each mouse with at least three seed cells.",
        "Red is gain, blue is loss. Read down a column to see whether a chromosome moves in the same direction across mice.",
        "One shared clone: +1, +8, +11, +15/16, +17, -18, -7/-12, present in 13 of 14 mice.",
        "Alternate chromosome labels are blanked; 19 labels collide at any legible size on a half-width axis."),
  "c": ("Inferred CNV by cell group",
        "Per-cell inferred copy number along the GRCm39 gene order, 300 randomly sampled cells per group.",
        "Cells are ordered within each group by clone correlation, so the clone-bearing cells sit at the top of their block.",
        "The tumour-call block shows the clone; the granule- and Purkinje-like blocks do not.",
        "A sample of 300 cells per group, not the full set."),
 },
 2: {
  "a": ("Human atlas annotation",
        "The human DMG atlas: 311,618 cells, malignant subpopulations in the warm palette, microenvironment lineages in the cool palette.",
        "This is the reference every mouse result is judged against; learn the colours here.",
        "167,012 malignant cells across 57 samples and 6 studies.",
        ""),
  "b": ("Mouse, author annotation",
        "Mouse cells on the same axes, coloured by the labels the depositing authors assigned.",
        "Hue is the compartment the authors assigned, lightness separates labels inside it.",
        "67,313 mouse cells, 14 author labels.",
        ""),
  "c": ("Mouse, lineage transferred from the human reference",
        "The same mouse cells coloured by the lineage the human reference model assigns.",
        "Compare with panel b: where the two disagree is where the authors' label and the human reference part company.",
        "31,354 mouse cells are absorbed into astrocyte.",
        "The human atlas has no cerebellar granule or Purkinje class, so those cells have no correct answer available and land in astrocyte at low posterior. This is reference coverage, not a statement about the mouse."),
  "d": ("Mouse tumour cells on the human malignant territory",
        "The 891 clone-defined tumour cells placed on the human malignant region, coloured by their best-scoring human subpopulation.",
        "Cropped to the human malignant territory; counts are in the box.",
        "OPC-like 230, cycling 285, AC-like 215, MES-like 102, OC-like 49, unassigned 10.",
        "The subtype label is a marker-score argmax, not a classifier transfer. Figure 9 tests whether it means anything."),
 },
 3: {
  "a": ("Marker genes, gene by gene",
        "For every marker gene of each human subpopulation: the fraction of mouse tumour cells and of human cells of that subpopulation in which it is detected.",
        "Points below the diagonal are detected less often in mouse than in human. Read the ordering within the mouse axis, not the absolute gap.",
        "Median mouse/human detection: OPC-like 0.146/0.330, cycling 0.118/0.503, AC-like 0.084/0.395, OC-like 0.064/0.130, MES-like 0.040/0.281.",
        "The absolute gap is expected: the human atlas is 88% nuclei and the mouse is whole cells."),
  "b": ("Pathway activity, human and mouse",
        "PROGENy activity for 14 pathways across the human subpopulations (left of the line) and the mouse populations (right).",
        "Both species are scored with the same human weights on ortholog-mapped matrices, so columns are directly comparable.",
        "Hypoxia -14.4 in mouse tumour against +9.4 in human MES-like; JAK-STAT +5.5 in mouse against -5.1 to -9.0 in every human subpopulation.",
        "Mouse-specific regulation is invisible to a human resource by construction."),
  "c": ("Transcription-factor activity",
        "CollecTRI activity for the 18 most variable TFs on the same profiles.",
        "Same layout as panel b: human left, mouse right.",
        "OLIG2 2.76 and SOX10 1.08 in mouse tumour; E2F4 5.72, MYC 12.81; HIF1A -1.20 against +2.98 in human MES-like; ASCL1 -0.27 and SOX2 -0.38 against positive in all human subpopulations.",
        "Per-TF target coverage is reduced by the ortholog restriction."),
  "d": ("Pseudobulk correlation per subpopulation",
        "Spearman correlation between each mouse group's pseudobulk profile and each human subpopulation's, over the 2,000 most variable genes.",
        "The tumour row is the one that answers Question 1; the other rows are normal mouse populations, shown as controls.",
        "Mouse tumour: cycling +0.47, OPC-like +0.35, OC-like -0.13, AC-like -0.36, MES-like -0.40.",
        ""),
 },
 4: {
  "a": ("Mouse cluster versus human class",
        "Spearman correlation between every mouse cluster and every human class on z-scored pseudobulk over 2,000 shared genes.",
        "Cluster sizes are in the row labels; tumour clusters are listed first.",
        "Clone-bearing tumour clusters correlate 0.37-0.52 with cycling and OPC-like.",
        "This is a class-level correlation and ranks OC-like higher than the subpopulation-level analysis in Figure 3d does."),
  "b": ("Human state sets, raw argmax",
        "Mean proportion of cells assigned to each human state, per mouse or per patient sample, under seven different cell-set definitions.",
        "Compare the mouse rows against the human rows at the bottom.",
        "Proportions shift substantially with the cell-set definition, which is why the definition is named on every row.",
        ""),
  "c": ("Human lineage MPs in mouse cells",
        "Median UCell score of each human metaprogram (top 50 genes) in mouse tumour cells and in each normal mouse population.",
        "Read across a row: if a normal population scores higher than tumour, the programme is not a tumour state here.",
        "Every AC-like metaprogram scores 4.9-8.7x higher in normal mouse astrocytes than in tumour cells.",
        ""),
  "d": ("Mouse de novo programmes versus human MPs",
        "The nine cNMF programmes factorised from the 891-cell primary tumour set against the human metaprograms and the Liu 2022 signatures.",
        "Colour is the hypergeometric -log10 q of the gene overlap; rows are labelled by each programme's top genes.",
        "p1 recovers MP1 cycling: 25 of 50 genes, Jaccard 0.368, q 5.5e-50.",
        "p5 (neuronal), p6 (nuclear) and p8 (mitochondrial) are technical programmes arising from ambient contamination, not tumour states."),
 },
 5: {
  "a": ("Mouse cells on the human reference",
        "All mouse cells placed on the human reference embedding by the zero-shot projection, coloured by transferred lineage.",
        "Grey is the human atlas underneath; lineages are named in the plot.",
        "",
        ""),
  "b": ("Distance from the human manifold",
        "Cumulative distribution of each cell's kNN distance from the human reference, for the mouse and for five human query cohorts.",
        "A curve to the left means cells sit closer to the reference. The mouse curve lying left of the human DMG cohorts is the point.",
        "Median kNN uncertainty: mouse tumour 0.795, against Filbin 0.998, Liu 0.992, LaBelle 0.920, fetal Braun 0.977, adult Siletti 0.574.",
        ""),
  "c": ("Author label vs transferred lineage",
        "For the 12 largest author labels, the fraction of cells assigned to each transferred lineage.",
        "Rows sum to 1. A row concentrated in one column means the two annotations agree.",
        "Granule-like and Purkinje-like go to astrocyte at 67.8% and 75.1%.",
        ""),
  "d": ("Mouse against the human query cohorts",
        "Median kNN uncertainty (solid) and the fraction of cells predicted malignant (pale) for the mouse sets and each human cohort.",
        "Compare the mouse bars against the human DMG cohorts rather than against zero.",
        "97.5% of the 891 clone-defined tumour cells are predicted malignant, at median posterior 0.9986.",
        ""),
 },
 6: {
  "a": ("Mouse myeloid",
        "The mouse myeloid compartment, clustered on its own, populations named in the plot.",
        "Hue is the compartment, lightness separates the populations within it.",
        "11,621 cells.",
        ""),
  "b": ("Mouse T/NK",
        "The mouse T and NK compartment, clustered on its own.",
        "Same shading rule as panel a.",
        "6,378 cells.",
        ""),
  "c": ("Mouse vascular",
        "The mouse vascular compartment, clustered on its own.",
        "Same shading rule as panel a.",
        "304 cells - the smallest compartment, and the least reliable.",
        "At 304 cells the vascular proportions are indicative only."),
  "d": ("Myeloid composition",
        "Mouse mean per mouse against human mean per patient sample, under the shared labelling rule.",
        "A bar present on one side only is a population with no counterpart in that species.",
        "Border-associated macrophage: human 23.3%, mouse 0. Homeostatic microglia: mouse 20.2%, human 0. Interferon myeloid: mouse 15.2%, human 0.",
        ""),
  "e": ("T/NK composition",
        "Same comparison for the T and NK compartment.",
        "Same reading rule as panel d.",
        "Exhausted T: human 20.0%, mouse 0. NK: mouse 33.5% against human 14.9%. Regulatory T: mouse 7.9%, no human counterpart.",
        "Mouse regulatory T cells are called without Foxp3, which is absent from the deposit."),
  "f": ("Vascular composition",
        "Same comparison for the vascular compartment.",
        "Same reading rule as panel d.",
        "Pericytes and tip/angiogenic endothelium are present in human and absent in mouse.",
        "Based on 304 mouse cells."),
 },
 7: {  # replaced at runtime by _figure7_panels(); kept so the dict shape stays uniform
  "_": ("Composition panels (a-d)",
        "Tumour cells, T/NK and myeloid as a share of all cells, and Mono/Macro as a share of immune cells.",
        "Dots are per-mouse values, bars are arm means; p-values above each panel are mouse-level models against Control.",
        "Nothing in composition reaches q < 0.1 under the primary model except one contrast (OC-like share, FUS+RT vs RT, q 0.0986).",
        "Two to four mice per arm; these p-values are exploratory."),
  "e": ("Myeloid and T/NK fractions (e-h)",
        "Interferon-myeloid and homeostatic microglia as a share of myeloid, NK as a share of T and NK, cycling as a share of tumour cells.",
        "Same dots-and-bars convention.",
        "No contrast in this row survives correction.",
        ""),
  "i": ("Tumour scores (i-l)",
        "Tumour MP4 hypoxia and MP8 MES-like UCell scores, and the tumour and microglia interferon set scores.",
        "Panel i is the figure's one nominally significant tumour-state result.",
        "MP4 hypoxia with RT +0.0153, q 0.0062; MP8 MES-like +0.0195, q 0.041 (11 mice).",
        "Panel i does NOT replicate under an independent method: PROGENy hypoxia in the same cells gives +0.31, p 0.877 (Figure 8a)."),
 },
 8: {
  "a": ("Tumour pathway activity by arm",
        "PROGENy pathway activity change against Control for each treatment arm, in tumour cells.",
        "Stars would mark q < 0.1; there are none. Blocks of colour here are not significant.",
        "222 pathway and TF tests, lowest q 0.228 (Androgen up with RT, +1.65, p 0.0054).",
        "This panel contains the independent hypoxia readout that contradicts Figure 7i."),
  "b": ("T/NK compartment, per mouse",
        "T-cell activation against T-cell dysfunction, one point per mouse, coloured by arm.",
        "Spread within an arm is as large as the difference between arms.",
        "9 mice with at least 50 T/NK cells.",
        ""),
  "c": ("Immune programmes by arm",
        "AUCell change against Control for six immune programmes in the myeloid and T/NK compartments.",
        "Same reading as panel a: colour without stars is not significance.",
        "36 programme tests, lowest q 0.459. Myeloid immune suppression is up in all three arms (FUS+RT +0.042 on a control mean of 0.164, p 0.026).",
        "Gene coverage is incomplete: T-cell activation 7/14 genes, regulatory T 3/6 with no Foxp3."),
  "d": ("Best-ranked programme tests",
        "Per-mouse values for the two best-ranked of the 36 programme tests.",
        "Dots are mice, bars are arm means; the x axis names the arms, so no legend is needed.",
        "Myeloid immune suppression (15 mice, p 0.026) and T/NK T-cell activation (9 mice, p 0.101).",
        "Neither survives correction."),
 },
 9: {
  "a": ("Transferred subtype",
        "The 891 QC-pass tumour cells on their own UMAP over 2,000 highly variable orthologs, coloured by transferred subtype.",
        "If the subtypes were discrete states they would occupy separate territory here. They do not.",
        "ARI between the embedding's own Leiden clusters and the subtype labels: 0.075.",
        "The subtype label is a marker-score argmax, which returns a label for every cell whether or not a state fits."),
  "b": ("Mouse of origin",
        "The same cells coloured by which mouse they came from - the confound control.",
        "Muted colours deliberately: this is a nuisance variable, not a result.",
        "One mouse contributes 45% of the cells; ARI with Leiden 0.064, silhouette -0.180.",
        ""),
  "c": ("Top genes outside the defining sets",
        "The strongest up-regulated genes per subtype that belong to no defining marker set, as Wilcoxon z across all five groups.",
        "Genes inside the defining sets are excluded because recovering them is circular.",
        "cycling validates (STMN1, HIST1H1B, TUBA1B); AC-like returns astrocyte genes (SLC6A11, ACSBG1, GRIA1); MES-like returns interferon (CXCL10, OAS3, USP18).",
        ""),
  "d": ("Is the separation real?",
        "Silhouette of the subtype labels in PCA space, against a 20-draw label permutation and against the mouse of origin, before and after removing the 137 defining marker genes.",
        "The decisive panel: if the red bar is inside the grey band, the labels carry no structure.",
        "Subtype -0.031 against a permuted baseline of -0.046 +/- 0.021; with markers removed, -0.045 against -0.045.",
        ""),
  "e": ("Library complexity",
        "Genes detected per cell by subtype.",
        "Explains why the OPC-like one-vs-rest DE returns a quarter of the transcriptome.",
        "OPC-like median 1,065 genes per cell against MES-like 484; Kruskal-Wallis H = 60.0, p = 2.9e-12.",
        "A one-vs-rest DE across groups differing two-fold in complexity is not interpretable."),
  "f": ("Subtype mix per mouse",
        "Subtype composition of each mouse's tumour cells.",
        "Shows how unevenly the subtype calls are distributed across animals.",
        "OPC-like is 61% Sample6.",
        ""),
 },
 10: {
  "a": ("Author annotation, as published",
        "The depositing authors' own UMAP of all 67,313 cells, reproduced from the coordinates shipped in the deposit.",
        "Hue is the compartment the authors assigned, lightness separates labels inside it.",
        "14 author labels; the coordinates match the authors' published figure exactly.",
        ""),
  "b": ("CNV clone signal",
        "Each cell's correlation with the shared tumour CNV profile, on the authors' own embedding.",
        "Dark red is clone-bearing. The signal is confined to one lobe.",
        "OPC-like tumor cells median 0.218; Neuronal-like 0.069; Purkinje-like 0.041.",
        ""),
  "c": ("Cells that pass the clone test",
        "The 891 clone-defined, QC-passing tumour cells highlighted against all other cells.",
        "Compare with panel a: the tumour is the small lobe, not the large mass labelled tumour.",
        "883 of the 891 (99.1%) fall in the OPC-like tumor cells cluster; 8 in Astrocyte-like tumor cells.",
        ""),
  "d": ("Lineage transferred from the human atlas",
        "The same cells coloured by the lineage the human reference assigns.",
        "Informative about malignant versus immune, not about fine cell types the reference never saw.",
        "Neuronal-like and Purkinje-like transfer to astrocyte, 21.9% and 24.5% predicted malignant; OPC-like tumor cells 82.5%.",
        "Choroid plexus epithelial cells come back 87.5% malignant and meningeal fibroblasts 53% - reference-coverage artefacts, since the human atlas has neither class."),
 },
}


# which of Dr Mei's questions each figure answers, quoted from her list
QUESTION = {
 1: "prerequisite for both questions: which mouse cells are really tumour",
 2: "'Show human annotation first' and 'for m06a, could you also show both human and mouse annotation in UMAP embedding'",
 3: "'Do mice have 4 human DMG tumor cell subpopulations (markers, pathways, TF activities, expression correlations)'",
 4: "'M05 d, do p1-p9 represent the NMF programs on tumor cells?' (panel d), plus state scores and metaprograms",
 5: "independent check on Q1: does the human reference model itself place the mouse cells where the marker evidence says",
 6: "'How do FUS, RT, and FUS+RT affect ... the tumor microenvironment' - the baseline mouse-vs-patient comparison",
 7: "'Check how FUS, RT, FUS+RT affect the tumor cell states' - composition, state scores, and DE genes",
 8: "'You can check DE genes, check pathway, check NMF programs' and 'how FUS, RT, FUS+RT affect inflammation, immune suppression, T cell dysfunction, T cell activation'",
 9: "'Generate a UMAP of the mouse tumor cells coloured by transferred subtype ... and examine the DE genes for each subtype to validate their identities'",
 10: "follow-up on Dr Mei's own published UMAP: which of the clusters labelled tumour carry the clone",
 11: "Dr Shenglin's proposal: 'cluster the mouse tumour cells de novo rather than transferring human labels'",
 12: "Dr Shenglin's proposal: 'then correlate the de novo clusters with the human subpopulations and compare markers, pathways and TFs'",
 13: "Dr Shenglin's proposal: 'then test FUS/RT effects on the de novo clusters'",
}


PANELS.update({
 11: {
  "a": ("De novo clusters",
        "UMAP of the 3,246 cells the authors label OPC-like tumor cells, re-embedded and Leiden-clustered without using any human label.",
        "Five clusters at resolution 0.2. This is the figure Dr Shenglin asked for: the mouse's own structure, found before any human annotation is consulted.",
        "M2 1,733 cells; M4 890; M3 500; M1 95; M5 28.",
        "No resolution between 0.2 and 1.0 satisfied the pre-stated floor of >=50 cells and >=3 mice in every cluster; the lowest resolution is used as the fallback and M5 (28 cells) is excluded from all statistics."),
  "b": ("Mouse of origin",
        "The same embedding coloured by which of the 15 mice each cell came from.",
        "A cluster occupying one colour would be a per-animal batch effect rather than a cell population.",
        "M2, M3 and M4 each draw on all 15 mice; M1 on 12; M5 on 10.",
        ""),
  "c": ("Cluster mix per mouse",
        "Share of each mouse's cells falling in each de novo cluster.",
        "Read across a row to see whether one animal dominates a cluster.",
        "",
        ""),
  "d": ("Cluster markers",
        "Top differentially expressed genes per de novo cluster, Wilcoxon, one column block per cluster.",
        "These names are what each cluster actually is, independent of any human label.",
        "M1 immune (ARHGDIB, AIF1, LCP1, CORO1A, DOCK2); M2 ribosomal (RPS12, RPL19, RPS18); M3 CSPG4/PTPRZ1/SOX6/VCAN-bearing; M4 neuronal (SYNE1, PCLO, RYR2, RBFOX3); M5 ciliated (CFAP44, SPEF2, DNAH6).",
        ""),
  "e": ("Clone signal per cluster",
        "Each cell's correlation with the shared tumour CNV profile, one box per de novo cluster.",
        "This is the panel that decides which clusters are tumour. The clone is genetic evidence and owes nothing to the human atlas.",
        "Median clone correlation M3 0.449, M2 0.259, M4 0.053, M1 0.047, M5 0.038.",
        "The clone profile is built leave-one-mouse-out."),
  "f": ("Overlap with the strict set",
        "Share of each de novo cluster that is also in the 891 clone-defined QC-pass tumour cells.",
        "An independent confirmation of panel e using the study's own strict tumour call.",
        "M2 38.1%, M3 36.4%, M1 7.4%, M4 3.7%, M5 0%.",
        ""),
 },
 12: {
  "a": ("Pseudobulk correlation",
        "Spearman correlation between each de novo mouse cluster and each human subpopulation, over the 1,995 most variable genes, z-scored across profiles within each species.",
        "Z-scoring removes the expression baseline every tumour profile shares, so a value here means distinctive resemblance, not overall similarity. Read a row for which human subpopulation a cluster leans towards.",
        "M4 OC-like +0.38 and cycling -0.47; M1 cycling +0.27; M2 cycling +0.21; M3 OC-like +0.17, OPC-like +0.14.",
        "With only four mouse profiles the z-score forces each row towards zero mean, so the magnitudes are small and the sign pattern, not the value, is the readable part. These numbers are not comparable with figure 3d, which uses six mouse groups."),
  "b": ("Marker-gene overlap",
        "Number of genes shared between each cluster's top 50 markers and each human subpopulation's defining set, with a hypergeometric test.",
        "This is the discriminating test in the figure. A dark cell is a conserved marker programme, not just a correlated profile.",
        "M3 x OPC-like: 4 genes (CSPG4, PTPRZ1, SOX6, VCAN), q 8.2e-5. M1 x MES-like: 2 genes (LGALS3, TGFBI), q 0.050. Nothing else below q 0.1.",
        "LGALS3 and TGFBI are also macrophage genes, and M1 is an immune cluster; the human MES-like set overlaps myeloid biology, so this hit is more likely contamination of the human set than a conserved mesenchymal programme."),
  "c": ("Pathways",
        "PROGENy pathway activity for the human subpopulations and the de novo mouse clusters side by side.",
        "Columns left of the rule are human, right are mouse. Compare a row across the rule.",
        "",
        "Human and mouse activities are computed on ortholog-mapped matrices but in separate runs; absolute values are not directly comparable across the rule."),
  "d": ("Transcription factors",
        "CollecTRI transcription factor activity, the 16 factors with the largest absolute activity.",
        "Same layout as panel c.",
        "",
        "Same caveat as panel c."),
 },
 13: {
  "a": ("Cluster composition by arm",
        "Share of each mouse's OPC-like-cluster cells falling in each de novo cluster, one dot per mouse, grouped by treatment arm.",
        "Spread within an arm is the noise level any treatment effect has to beat.",
        "15 mice: Control 4, FUS 4, RT 3, FUS+RT 4.",
        ""),
  "b": ("Composition shift",
        "Quasi-binomial log-odds coefficient per cluster and contrast against Control, with Benjamini-Hochberg q values.",
        "Nothing in this panel is significant. It is shown because the absence is the result.",
        "12 contrasts, lowest q 0.495; largest raw effects M3 with FUS +0.99 log-odds (p 0.053) and M4 with FUS -1.35 (p 0.146).",
        "Dispersion is estimated per cluster from the Pearson chi-square and the test uses the residual degrees of freedom, so the unit of replication is the mouse, not the cell."),
 },
})


def _figure7_panels():
    """The twelve panels of figure 7, in the order m08 draws them, with each panel's own p-values from the tables."""
    import pandas as pd
    comp = pd.read_csv(RESULTS / "tables/m08_composition.csv")
    sco = pd.read_csv(RESULTS / "tables/m08_scores.csv")
    comp_panels = [("compartment", "tumour", "tumour cells / all cells"), ("compartment", "T_NK", "T/NK / all cells"),
                   ("compartment", "myeloid", "myeloid / all cells"), ("immune_type", "Mono/Macro", "Mono/Macro / immune"),
                   ("myeloid_label", "interferon_myeloid", "IFN-myeloid / myeloid"),
                   ("myeloid_label", "homeostatic_microglia", "homeostatic MG / myeloid"),
                   ("T_NK_label", "NK", "NK / T+NK"), ("tumour_state", "cycling", "cycling / tumour (state)")]
    score_panels = [("human_MP_ucell", "tumour", "MP4_hypoxia", "tumour: MP4 hypoxia UCell"),
                    ("human_MP_ucell", "tumour", "MP8_MES_like", "tumour: MP8 MES-like UCell"),
                    ("pseudobulk_geneset", "tumour", "Interferon", "tumour IFN set score"),
                    ("pseudobulk_geneset", "Microglia", "Interferon", "microglia IFN set score")]

    def stat(df):
        df = df[df.model == "arm"].set_index("contrast")
        bits = []
        for c, lab in (("FUS_vs_Control", "FUS"), ("RT_vs_Control", "RT"), ("FUS+RT_vs_Control", "FUS+RT")):
            if c in df.index:
                bits.append(f"{lab} p {df.loc[c, 'p']:.3g} (q {df.loc[c, 'q']:.3g})")
        return "; ".join(bits)

    out = {}
    for i, (s_, c_, title) in enumerate(comp_panels):
        out["abcdefgh"[i]] = (
            title,
            f"Per-mouse proportion for {title}, one dot per mouse, horizontal bar at the arm mean.",
            "Compare the spread within an arm against the gap between arms before reading any difference.",
            stat(comp[(comp.set == s_) & (comp.category == c_)]),
            "Two to four mice per arm; p-values are exploratory and the composition result depends on whether the inferred processing block is modelled.")
    for i, (f_, pop, sc_, title) in enumerate(score_panels):
        out["ijkl"[i]] = (
            title,
            f"Per-mouse mean of {title}, one dot per mouse, horizontal bar at the arm mean.",
            "Scores are on the tumour cells of the 11 mice with at least 20 QC-pass tumour cells.",
            stat(sco[(sco.family == f_) & (sco.population == pop) & (sco.score == sc_)]),
            "Panel i does not replicate under an independent method: PROGENy hypoxia in the same cells gives +0.31, p 0.877 (figure 8a)."
            if sc_ == "MP4_hypoxia" else "")
    return out


def page_mm(pdf):
    d = pdf.read_bytes()
    m = re.search(rb"MediaBox\s*\[\s*([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)", d)
    return f"{float(m.group(3)) / 72 * 25.4:.0f} x {float(m.group(4)) / 72 * 25.4:.0f}" if m else ""


def main():
    rows = []
    for num, stem, title, section in FIGURES:
        d = FIGDIR.get(num, FIG)
        png, pdf = d / f"{stem}.png", d / f"{stem}.pdf"
        prov_path = d / f"{stem}.png.provenance.json"
        prov = json.loads(prov_path.read_text()) if prov_path.exists() else {}
        inputs = "; ".join(Path(k).name for k in prov.get("inputs_sha256", {}))
        rendered = datetime.fromtimestamp(png.stat().st_mtime, timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
        panels = _figure7_panels() if num == 7 else PANELS[num]
        for panel, (ptitle, shows, how, key, caveat) in panels.items():
            rows.append({
                "panel_tag": f"Fig {num}{panel}",
                "figure": num, "figure_title": title, "report_section": section,
                "question_answered": QUESTION[num],
                "panel": panel, "panel_title": ptitle,
                "what_it_shows": shows, "how_to_read_it": how, "key_numbers": key, "caveat": caveat,
                "png": png.name, "pdf": pdf.name, "page_mm": page_mm(pdf),
                "script": prov.get("script", ""), "inputs": inputs, "rendered_utc": rendered,
                "png_kb": round(png.stat().st_size / 1024), "pdf_kb": round(pdf.stat().st_size / 1024),
            })
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print(f"{OUT.name}: {len(rows)} panel rows across {len(FIGURES)} figures")
    missing = [r["figure"] for r in rows if not r["page_mm"] or not r["script"]]
    print("rows missing page size or script:", sorted(set(missing)) or "none")


if __name__ == "__main__":
    main()
