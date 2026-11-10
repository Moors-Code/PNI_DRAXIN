# Tumour-Conditioned Nerves Establish a DRAXIN+ Niche Associated with Metabolic Adaptation in Pancreatic Cancer

This repository contains the computational workflows and statistical analyses used in the study:

**“Tumour-Conditioned Nerves Establish a DRAXIN+ Niche Associated with Metabolic Adaptation in Pancreatic Cancer”**

Nerve categories used throughout: **InvN**, invaded nerves; **NI-TB**,
non-invaded nerves of a tumour-bearing section; **NI-Adj**, non-invaded nerves
of an adjacent, tumour-free section.

## Data Availability

| Dataset | Source | Accession |
|---|---|---|
| Spatial transcriptomics, raw (FASTQ) and processed (h5, spatial coordinates) | This study | [GSE333985](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE333985) |
| Bulk RNA-seq, raw (FASTQ) | This study | [GSE333991](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE333991) |
| CRISPR-Cas12 knockout screen sequencing, raw (FASTQ) | This study | [GSE333994](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE333994) |
| Human PDAC snRNA-seq, processed (h5) | Chen et al., 2025 | [GSE278688](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE278688) |
| Human PDAC scRNA-seq, processed (h5) | Chen et al., 2025 | [GSE278689](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE278689) |

## Repository Structure

### 1. `bulk_RNA_seq`

Scripts for bulk RNA-seq count processing, normalisation, batch-aware differential expression analysis, exploratory quality control, and reusable RNA-seq plotting and modelling functions.

| Script | Description |
|---|---|
| `count_merging.R` | Merges RNA-seq count matrices and extracts gene-length information for downstream RPKM normalisation and transcript-level quantification. |
| `DEG_modeling_experiment_libraryprep_2026_02_26.R` | Batch-corrected edgeR quasi-likelihood differential expression analysis comparing ctrl, DRAXIN, and diluted DRAXIN conditions across independent experiments. Includes TMM normalisation, dispersion estimation, PCA, volcano plots, and differential expression testing. |
| `functions.R` | Custom utility functions for RNA-seq analysis, including edgeR/DESeq2 differential expression workflows, normalisation, PCA, volcano plotting, correlation analysis, GSEA visualization, and sample-matching utilities. |

---

### 2. `Cas12_screen`

Workflows for Cas12a pooled-screen quality control, guide-count assembly, MAGeCK-compatible filtering and
reporting for the Cas12a pooled screen (DRAXIN-treated versus control).
Listed in the order they are run.

| Script | Description |
|---|---|
| `01_libraryQC.Rmd` | Quality control of the Cas12a guide-count libraries. Library-size normalisation (CPM and log2(CPM+1)), Spearman correlation between replicates, guide-level variability (CV versus abundance), PCA on variance-stabilized counts (DESeq2 VST), per-replicate distribution metrics including Gini index, and comparison against the plasmid library stock. |
| `02a_merging_guide_count.R` | Assembles the analysis count matrix, combining the per-condition screen count files with the library QC counts into a single guide × sample table. Guides are matched by identifier rather than by row order. |
| `02_ScreenQC.Rmd` | Quality control of the merged screen count matrix. Detected sgRNAs and skew ratio per sample, CPM and TMM-normalised logCPM, density and histogram distributions, PCA, MDS, Pearson correlation heatmap, pairwise scatterplots, MD plots, and control checks comparing common-essential and non-targeting guides against the remaining library. |
| `03_MAGECK_pre_filtering.Rmd` | Generates MAGeCK input matrices. Guides are retained on mean voom logCPM (cutoff 2.5) computed in a chosen set of reference samples; three strategies are run (library input + T0, T0 only, T0 + ctrl), each writing raw counts for MAGeCK together with the per-guide filtering decisions. |
| `run_mageck_all_variants.sh` | Runs MAGeCK on the screen: a paired test of the two DRAXIN-treated against the two control replicates, repeated across seven combinations of input matrix (unfiltered, T0-filtered, T0 + ctrl filtered) and normalisation method (non-targeting control, median, total). The analysis reported in the paper is `T0_PLUS_CONTROL_PAIRED_TOTAL`; the remaining variants are included to show that the reported hits do not depend on the filtering or normalisation choice. Run on an HPC cluster (SLURM) with MAGeCK 0.5.9.5 (RRA 0.5.9). |
| `04_MAGECK_figures_output.R` | Draws the published screen figures from the MAGeCK gene summary: ranked RRA score distributions for enriched and depleted arrays, and the volcano of −log10 FDR against log2 fold change. MAGeCK is not re-run, so the plotted values are those of the original analysis. |
| `05_MAGECK.report.Rmd` | Standard MAGeCK report for the T0 + ctrl paired comparison, giving sgRNA- and gene-level statistics, RRA and p-value distributions, per-guide read counts for top genes, and essential-gene dropout by GSEA. Adapted from the report template distributed with MAGeCK (W. Li). |

Paths in `run_mageck_all_variants.sh` are placeholders and must be set before
running. Steps `04` and `05` read the `T0_PLUS_CONTROL_PAIRED_TOTAL` output.

---

### 3. `Chen_reference_DRAXIN`

DRAXIN expression across cell types in the published sc/snRNA-seq reference
(Chen et al., 2025), as independent corroboration that the transcript is found
preferentially in Schwann cells. Restricted to tumour tissue of the PNI-high
patients, because adjacent-normal samples are not paired with tumour in that
cohort.

| Script | Description |
|---|---|
| `01_DRAXIN_chen_reference.R` | Takes the patients, not the cells, as replicates. Counts are summed per patient and cell type and tested with limma-voom, patient as a blocking factor, each cell type against the average of the others. Because Schwann cells are rare after dissociation, the primary analysis for cell identity is the within-patient detection fraction instead: the percentage of cells with DRAXIN detected in a cell type against the percentage among all other cells of the same patient, compared across patients by a paired test, with the same rule applied to every cell type. Descriptive values are computed per patient and averaged with equal weight, so a patient with many cells does not dominate, and the per-patient values are written out. |

Dissociated data carry no spatial context, so this addresses cell-type identity
only; it says nothing about invaded versus non-invaded nerves.

---

### 4. `enrichment_analysis`

Pathway and gene-set analysis across bulk RNA-seq, the CRISPR screen and the
spatial nerve data. Each analysis is split into a compute step, which saves
result objects, and a plot step, which reads them — so figures can be
regenerated without re-running GSEA.

| Script | Description |
|---|---|
| `00_annotate_gene_symbols.R` | Adds HGNC symbols to the edgeR differential expression tables by mapping Ensembl IDs with `org.Hs.eg.db`. |
| `01a_GSEA_compute.R`, `01b_GSEA_plots.R` | GSEA on the bulk RNA-seq comparisons (DRAXIN vs ctrl, diluted vs ctrl, DRAXIN vs diluted), ranked by signed −log10 p. Covers GO (BP, CC, MF), Reactome, Disease Ontology, KEGG and MSigDB Hallmark, followed by a pass across all MSigDB collections. The plot step draws dotplots split by NES sign and writes gene-symbol result tables. |
| `01c_ssGSEA_compute.R`, `01d_ssGSEA_plots.R` | ssGSEA pathway scoring of batch-corrected logCPM across 0, 17 and 68 nM DRAXIN, for a fixed panel of Hallmark, GO and KEGG gene sets. Scores are standardised both across all samples and relative to the control group; the plot step fits limma on the raw scores to test the three dose contrasts and annotates the panels with FDR. |
| `01e_heatmap_compute.R`, `01f_heatmap_plots.R` | Leading-edge gene heatmaps for selected pathways. The compute step takes the union of leading-edge genes across the three comparisons and maps them onto the batch-corrected expression matrix; the plot step draws one row-z-scored heatmap per pathway, split by dose. |
| `02a_CRISPR_GSEA_compute.R`, `02b_CRISPR_GSEA_plots.R` | GSEA on the MAGeCK gene summary, ranked by LFC × −log10 p, with the resistance and sensitisation arms ranked separately, across Hallmark, C2, C5 and Reactome. Multi-target guides, dedup-collapsed gene pairs and pooled control pseudo-genes are excluded from ranking and exported separately. The plot step produces the paper's NES comparison figure and the full results table. |
| `02c_CRISPR_GSEA_exploration_supplementary.R` | Broader exploration of the CRISPR GSEA results: both filtering modes, all collections, per-result dotplots and a paired resistance-versus-sensitisation NES scatter. Not part of the reported figures. |
| `03a_Figure2_GSEA_compute.R`, `03b_Figure2_GSEA_plots.R` | GSEA on the invaded versus non-invaded nerve differential expression, ranked by signed −log10 FDR. Covers GO, Reactome, WikiPathways, Disease Ontology, Network of Cancer Genes, KEGG, all MSigDB collections, and a custom cell-state signature set supplied as a GMT. |
| `03c_Figure2_ssGSEA_compute.R`, `03d_Figure2_ssGSEA_plots.R` | Pathway activity across NI-Adj, NI-TB and InvN nerves. Spot-level ssGSEA scores for an immune and complement pathway panel are aggregated to the nerve level and z-scored across nerves; the plot step draws the two Figure 2 panels. |
| `03e_Figure2_ssGSEA_stats_supplementary.R` | Repeated-measures ANOVA across nerve categories on the patient-aggregated pathway scores. Not used in the paper. |
| `combine_GSEA_results.Rmd` | Concatenates the per-collection GSEA result tables from the invaded-nerve analysis into one merged file for downstream pathway selection. |

`03c_Figure2_ssGSEA_compute.R` reads the spatial `SpatialExperiment` object
and its annotation from a separate project directory (`EXTERNAL_SPATIAL_DIR`
at the top of the script), not from this repository's `input/` folder, as the
spatial data are large and shared across projects. Set that path to wherever
the spatial data are located, or to the deposited dataset, before running it.

---

### 5. `IF_analysis`

Quantification of multiplexed immunofluorescence and RNAscope HiPlex staining in
the validation cohort (five patients, six FFPE blocks). Nerves are compared
within a tissue section throughout, so that staining and imaging are held
constant: NI-TB (non-invaded nerves in a tumour-bearing section) against InvN
(invaded nerves of the same section). NI-Adj nerves lie on a separate section
and are therefore shown descriptively, without statistics.

Three segmentation pipelines feed this folder — cell-level Horizon exports,
region-level QuPath exports, and the COMET immune panel — each with its own
export format. Within the Horizon pipeline the scripts are listed in the order
they are run; `00_settings.R` holds the paths and is sourced by the rest.

| Script | Description |
|---|---|
| `00_settings.R` | Shared settings for the Horizon pipeline: which cleaned export the analysis scripts use, the per-block threshold file, patient colours matching the QuPath and COMET figures, and the SuperPlot function used by several figures (small dots = nerves or ROIs coloured by patient, large dots = block medians, bar = median of block medians). |
| `01_clean_cancer.R` | Cleans the tumour-ROI exports. Parses patient, block, ROI and ROI type out of the annotation path, merges the duplicated per-ROI columns into one column per measure, and combines all exports of the same image — so an ROI re-exported later replaces the incomplete earlier copy. Writes per-image and per-ROI cleaning reports plus an issues report listing missing analyses, gaps in the ROI numbering, and ROIs present in more than one export. |
| `02_clean_nerve.R` | The same cleaning and combining for the nerve-ROI exports, with ROI names resolved to invaded or uninvaded. |
| `03_cancer_analysis.R` | PanCK⁺ tumour cells in perineural versus distal ROIs, for five RNAscope probes (mean dots per cell) and five proteins (mean intensity). ROIs are the unit of analysis and patient a random intercept, so the estimate is the within-tumour difference; reported as fold changes with 95% confidence intervals, Holm-adjusted, with leave-one-patient-out estimates and a patient-as-replicate model alongside. CD3 and CD45 serve as same-channel background controls. The Visium discovery patient is analysed separately and is not part of the cohort statistics. |
| `04_nerve_DRAXIN_analysis.R` | DRAXIN transcript density per nerve (dots per 100 µm², total dots over total cell area) against markers of the repair-like Schwann cell state: GAP43 density, and the share of p75NTR⁺, S100B⁺ and NESTIN⁺ cells. Linear mixed model of log2 DRAXIN with patient and block as random intercepts, Holm-adjusted over the four markers, with per-block Spearman correlations, adjustment for invasion status and an uninvaded-only analysis as robustness checks. Also compares DRAXIN in InvN against NI-TB nerves of the same sections, with random-slope, block-level and leave-one-block-out versions. |
| `05_nerve_NESTIN_ITGB4_analysis.R` | NESTIN and ITGB4 in InvN against NI-TB nerves of the same blocks, as the share of positive cells (binomial mixed model on positive and negative cell counts per nerve, with a nerve-level term for overdispersion) and as median cell intensity (linear mixed model on log2 intensity), the latter being threshold-free. Every model is refitted with a second optimiser and the convergence notes are written to the results table rather than the console. |
| `06_representative_examples.R` | Suggests which ROIs and nerves to show as example images. For each finding the block whose effect is closest to the typical effect across blocks is chosen, and within it the ROIs or nerves closest to their group median; for the DRAXIN correlation, nerves that lie close to the trend and far apart along it. Run after `03`, `04` and `05`. |
| `07_distance_to_nearest_nerve.R` | Distance from every annotated tumour ROI to the nearest annotated nerve, using the whole-image cell coordinates shared by the two cleaned exports. Distal ROIs were selected visually as lying away from any viable nerve; this states that criterion as a number and allows the imaging definition of "distal" to be compared with the 1.3 mm threshold used for Tum-DistN in the Visium analysis. Perineural ROIs are measured the same way as a positive control. |
| `fig2E_stats.R` | CD45⁺ cell density in the perineural area across NI-Adj, NI-TB and InvN, from the QuPath region-level export. Negative-binomial mixed model of the CD45⁺ count with nerve group as fixed effect, patient and image as random intercepts and log perineural area as an offset; run both on the full perineural area and after excluding tumour-occupied area, with all three pairwise comparisons Holm-adjusted within each model. |
| `fig2E_plot.R` | Draws the Figure 2E SuperPlot from the per-nerve table and p values written by `fig2E_stats.R`. |
| `markers_stats.R` | GAP43 and p75NTR in the nerve, p75NTR in the perineural ring, and tenascin-C in the 0–50 µm ring. Each intensity is normalised to the median intensity of the same marker in the normalisation rings of the same section, so that 1 equals section background; InvN is then compared against NI-TB with a linear mixed model on log intensity and patient as a random intercept, with a paired test on patient medians as a patient-level check. Because the comparison is within a section, the normalisation factor cancels out of it. |
| `markers_plot.R` | Marker SuperPlots from the outputs of `markers_stats.R`: NI-TB against InvN with the model p value as the main figures, and NI-Adj, NI-TB and InvN together as descriptive supplementary figures. |
| `markers_representative_nerves.R` | Ranks patients by how closely their InvN / NI-TB ratios match the overall effect, and within a patient the nerves closest to their group median across all markers at once, so that one multi-channel image panel can illustrate several findings. Very small nerves are skipped. |
| `plotting_relative_toNI_TB.R` | Brings the Horizon and QuPath markers onto a single scale: each nerve is expressed relative to the geometric mean of the NI-TB nerves of its own section, so NI-TB centres on 1 and InvN shows the within-section fold difference, drawn on a log2 axis. Nothing is re-analysed; the p values are read from the two statistics outputs and Holm-adjusted across the panels of the figure. |
| `Classification_Plots.Rmd` | Presence of a surrounding p75NTR⁺ perineurial layer in invaded against non-invaded nerves. Three annotators scored each nerve independently and the reported call is their majority vote, with Fleiss' kappa given alongside. All annotated nerves are included; an earlier version restricted to nerves without annotator comments was rejected because comments were recorded far more often for invaded nerves, and that restricted set is retained only as a sensitivity analysis. |
| `processing_ROIs.R` | Combines the raw per-ROI exports of the COMET immune panel into one per-cell table. The rule-based classification columns are identified by their values rather than their names, because the column names differ between exports, and the nuclear marker intensities and dot counts are matched by pattern for the same reason. |
| `immune_cells_0_200um.Rmd` | Immune composition and inflammatory milieu in the ROIs lying within 200 µm of a nerve, the window reported in the paper. Each ROI contributes one value per measure: the share of immune cells in each population, the CD45⁺ fraction of all cells, and the share of CD45⁺ cells positive for IFNγ, CXCL9, IL1B or PD-L1 — so the four inflammatory markers share one denominator and can be read on the same axis. Distance is used only to select the window; no distance modelling is performed. |

The Horizon scripts read the raw exports from `input/Cancer_analysis/`,
`input/Cancer_analysis_with_background/`, `input/Nerve_analysis/` and
`input/Nerve_analysis_with_control/`, together with `input/thresholds.csv`, which
holds one row per block and the per-marker intensity thresholds. `Image_Name` in
the threshold file must be written exactly as in the data; the scripts print both
lists when they do not match. A block without data for a marker, or without a
threshold for it, is left out of that marker only, and this is reported.

The QuPath scripts read `input/plotting_data.csv`. `plotting_relative_toNI_TB.R`
additionally reads `Per_nerve_values.csv` and `Statistics.csv` from the output of
`05_nerve_NESTIN_ITGB4_analysis.R`, which belongs to the Horizon project rather
than to this folder; set `NESTIN_DIR` at the top of the script to that output
directory before running it.

`Classification_Plots.Rmd` expects `input/nerve_classification.xlsx`, the invaded
or non-invaded call per image, and `input/nerve_annotator_1.xlsx` to `_3.xlsx`,
one per annotator with image, label and comments. `processing_ROIs.R` needs
`RAW_EXPORT_DIR` set to the folder holding the raw ROI exports. The rule-based
labels it writes were corrected manually before use, so
`immune_cells_0_200um.Rmd` reads
`input/combined_selected_columns_corrected_manually.csv` together with
`input/ROIs_immune_cells.csv`, which gives each ROI's distance to the nearest
nerve.

---

### 6. `K_sweep`

K-sensitivity analysis of the sparse NMF decomposition, performed during revision
to test the choice of K = 30 used in `spatial_data_processing`. The sweep refits
the NMF at K = 15, 20, 25, 30, 35 and 40 and asks whether the chosen K is
supported by reconstruction quality and by correspondence with RCTD-derived
cell-type fractions. Listed in the order they are run.

| Script | Description |
|---|---|
| `run_hvg_matrix_container.sh`, `make_hvg_matrix_container.R` | Builds the highly variable gene matrix the sweep is run on: per-sample gene variance (`modelGeneVar`) combined across samples (`combineVar`), top 5000 genes retained, saved as their log-count matrix. Submitted as a SLURM job. |
| `run_ksweep_container.sh` | SLURM array job fitting a sparse NMF to that matrix at each K, one task per K, and saving each fit together with its relative Frobenius reconstruction error. The NMF settings reported in the Methods are set here. |
| `04_Ksweep_analysis.R` | Combines the per-K reconstruction errors into one table, and for every K identifies the factor most positively correlated with the RCTD-derived Schwann-cell fraction and the one most positively correlated with the malignant-cell fraction, alongside the summed activity of the top three correlated factors. Cosine similarity of gene loadings is then used to ask whether the selected factor recovers a similar programme at different K, within and across fits. Produces the K-sensitivity supplementary figure. |

The first two scripts run on an HPC cluster (SLURM) inside a
Singularity/Apptainer container built from the Docker image
`bioconductor/bioconductor_docker:RELEASE_3_17` (Bioconductor 3.17, R 4.3), with
`RcppML` 0.5.5 installed; set `CONTAINER` in each shell script to that image
before submitting. `04_Ksweep_analysis.R` then runs locally and reads the sweep
outputs from `input/ksweep/`, keeping the cluster job's filenames, together with
`input/RCTD_weights.Rds`.

Because RCTD is estimated from the same spatial expression data as the NMF, this
is a reference-based consistency assessment rather than independent validation.

---

### 7. `LR_analysis`

Ligand-receptor co-expression at the tumour-nerve interface, from spatially
weighted bivariate scores through region-level statistics to over-representation
analysis of the prioritised interactions. Listed in the order they are run.

| Script | Description |
|---|---|
| `01_LR_analysis.ipynb` | Computes local ligand-receptor co-expression per spot with LIANA, using a spatial neighbour graph with an exponential kernel narrow enough to include mainly the directly neighbouring spots, and aggregates the scores to one value per region (one nerve or one tumour region per section). Each pair is then tested across invaded nerves, non-invaded nerves and tumour regions by Kruskal-Wallis with pairwise Mann-Whitney tests and FDR correction, and the direction of effect is checked per patient. Writes the master table of all tested pairs, the Tier 1 and Tier 2 tables, and the query and background gene lists the over-representation analysis reads. Produces Figures 3F and S10A-C. |
| `02_LR_patient_stratified_sensitivity.ipynb` | Sensitivity analysis for the patient structure, read-only: re-tests every pair with patient taken into account, by a van Elteren stratified rank test and by Stouffer combination of per-patient one-sided tests, and reports how many Tier 1 pairs remain significant and how many are higher in invaded nerves in all three patients. Nothing upstream is re-run and no existing table or figure is changed; results are written to a separate `sensitivity` folder. |
| `03_LR_interface_ORA.Rmd` | Over-representation analysis of the genes in the Tier 1 interactions against the genes of all tested pairs as background, and the enrichment network of the resulting terms. Produces Figure S10D. |

Paths are relative to the project folder, the one containing `input/`. Run the
notebooks and the script from there.

`01` reads `input/clean_adata.h5ad`, exported by
`05c_deconvolution_composition_figures.Rmd` in `spatial_data_processing`, and
writes the pseudobulk object and tables to `output/data/`, which `02` and `03`
then read. Tier 1 is the set used for the over-representation analysis and for
the pairs shown in the figures. Outputs named `exploratory_` are not in the
manuscript.

---

### 8. `pathologist_scoring`

Independent pathologist review of the computationally called nerves, and the
agreement between the two annotations. Listed in the order they are run.

| Script | Description |
|---|---|
| `01_pathologist_nerve_annotation_export.Rmd` | Prepares the material for review: one row per nerve with its identifier, section, patient, invasion status, spot count and centroid and bounding box, given in both the object's own coordinates and the embedded image's pixel units; and one labelled overview image per section, H&E with the spots coloured by area and each nerve identifier printed at its location, together with mirrored versions with and without labels for matching against the raw image files. |
| `02_pathologist_agreement_score.R` | Agreement between the computational and the pathologist annotation, reported as three separate things rather than one number: the computational detections the pathologist judged not to be nerve; invasion-status agreement on the nerves both sides gave a valid call for, as a fixed confusion matrix with raw and strict agreement, Cohen's kappa and sensitivity, specificity and predictive values; and the subtype breakdown, which shows where the disagreements concentrate. Labels from both files are mapped onto a fixed vocabulary, and nerves with an ambiguous, unexpected or missing call are counted and listed rather than dropped. Read-only on both input files. |

The pathologist's completed spreadsheet goes into `input/`; `02` reads the
export written by `01`.

---

### 9. `screening_candidates_in_vitro`

The 15 candidate genes taken into the in vitro screen: what is known about each
one from the spatial data, and their expression as a function of distance from a
nerve. Listed in the order they are run.

| Script | Description |
|---|---|
| `01_screening_candidates_summary_step1.R` | Assembles the candidate summary table: spatial association with invaded nerves, differential expression in invaded versus non-invaded nerves, ligand-receptor support in the revised LR analysis, mean expression per nerve category and patient, and how each candidate was nominated. The descriptive summary column is generated from these values rather than written by hand, so it always matches the numbers. Records the symbol correction for placental growth factor, nominated as PIGF and reported as PGF. |
| `02_screening_candidates_decay_step2.R` | One distance-decay panel per candidate, ordered by nerve association, each with the curves of the three nerve groups and, as the control, the geometry-matched relocated shapes of the same nerves. Distances, bins and controls are those of the plasma-cell analysis; the background line is the mean of tumour-bearing spots that are neither nerve nor in contact with one, the same background as the spatial test. Curves are averaged per nerve, then per patient, and shown as the mean and standard error across patients. |

`01` reads `lee_sham_all_candidates.csv` and `02` reads
`figure_bundle_for_real_vs_relocated.rds`, both written in
`spatial_enrichment_prioritisation`; place a copy of each in `input/`, together
with a copy of `00_decay_shared.R`, so the zones and bins are identical.

---

### 10. `spatial_data_processing`

Visium spatial transcriptomics workflow, from SpatialExperiment construction and
quality control through NMF decomposition, nerve identification and annotation,
pseudo-bulk differential expression, sc/snRNA-seq reference preprocessing and
RCTD deconvolution, to the DRAXIN and gene-panel analyses. Listed in the order
they are run.

Several steps are split into a compute script and a figures script. The compute
step writes a handoff file and is not meant to be re-run; the figures script
reads that file and never recomputes it, so figures can be regenerated freely.

| Script | Description |
|---|---|
| `00_Making_SPE.R` | Builds the integrated SpatialExperiment object from the Space Ranger output of the six Visium sections: raw counts subset to the filtered probe set, spatial coordinates and tissue images. Samples are ordered explicitly, because column order can affect downstream results and the pipeline's wildcards are unordered. |
| `01_QC_filtering.Rmd` | Quality control and filtering. Spots are retained on library size (3 median absolute deviations below the median, in log space, within each section) and on detected genes (> 500); genes with fewer than 100 counts across the dataset are dropped; the remaining spots are library-size normalised. Mitochondrial and haemoglobin fractions are inspected but not used for filtering. A few areas of adipose tissue are set in-tissue by hand after review of the H&E images. Produces Figure S1A-G. |
| `02a_nmf_fitting_euler.R` | Fits the sparse NMF model. Highly variable genes are selected per section and combined across sections, and the top 5000 are factorised with K = 30 and L1 = 0.5. Computationally heavy, so run on an HPC cluster; the visualisation script only ever reads its output. |
| `02b_nmf_visualization_and_export.Rmd` | Interprets and exports the fitted model: factor representation across tissue blocks, mean activity per spot, top gene loadings for all and for the five highlighted factors, and the cosine similarity between factors. Collapses the 30 factors into broader biological classes using a manually curated factor-to-class mapping, and writes the annotated object every later script reads. Produces Figures 1C, 1D, 1E, S3 and S4A-C. |
| `03a_nerve_identification_actual_compute.Rmd` | Calls nerves and classifies them. Spots with more than 20% nerve-factor signal are grouped into individual nerves by single-linkage clustering of their coordinates, nerves of fewer than two spots are dropped, and a nerve is called invaded when its mean tumour-factor contribution exceeds 1%. Also identifies the spots in direct contact with each nerve and the wider nerve-proximal zone. Writes the nerve annotation table. |
| `03b_nerve_identification_figures.Rmd` | Figures for the step above, read from the frozen annotation table rather than recomputed. Produces the nerve- and tumour-signal panels of Figure 1A and Figures S5A-C. |
| `04a_DE_analysis_COMPUTE.Rmd` | Defines the four spatial groups (invaded nerve, non-invaded nerve, nerve-adjacent tumour, nerve-remote tumour), aggregates spots into patch-level pseudo-bulk profiles, and runs the differential expression for all four contrasts with limma, modelling the dependence between pseudo-bulks from the same section. Prioritises invasion-associated and tumour-proximity genes by pairing each contrast with a control contrast. Writes the spot annotation table every later script reads. |
| `04b_DE_analysis_figures.Rmd` | Figures and full result tables for the step above. Produces Figures 1F, 1G, 1H, 2A, 3A, 3C, 6G and S17A-B. |
| `05a_scRNA_snRNA_preprocessing.Rmd` | Preprocesses and integrates the published sc/snRNA-seq reference (Chen et al., 2025) for use in deconvolution, following the authors' own parameters, with their cell-type annotation mapped on by barcode and ductal cells split into malignant and normal by tissue of origin. |
| `05b_RCTD_compute.Rmd` | Builds the RCTD reference and query objects and processes the output into row-normalised cell-type weights per spot. The deconvolution itself is heavy and was run on an HPC cluster. Writes the weight table every later script reads. |
| `05c_deconvolution_composition_figures.Rmd` | Cell-type composition from the deconvolution: per nerve, by location relative to a nerve, and as an enrichment heatmap across spatial categories; and four individual immune populations at spot and nerve level. Produces Figure 2F. Also exports the object for the Python analyses. |
| `05d_deconvolution_immune_macrophage_figures.Rmd` | Immune and macrophage composition across the three nerve categories, at spot and nerve level, with mixed-effects models of the adaptive-immune and macrophage-programme signal (patient and nerve as random effects). Produces Figures S6F-G. |
| `06_plasma_and_adipocyte_signal_plotting.Rmd` | Spatial maps of the plasma-cell deconvolution weight and of the adipocyte NMF programme, each with nerve outlines coloured by invasion status. Produces Figures S6A and S7A. |
| `07_gene_expression_heatmaps.Rmd` | Expression of curated gene panels across the three nerve categories, averaged within each patient and then across patients and z-scored per gene. Produces Figure 2G. |
| `08_draxin_analysis.Rmd` | DRAXIN across nerve categories, with the nerve rather than the spot as the unit of analysis. Correlates DRAXIN with the nerve and tumour NMF programmes, with the MSigDB Hallmark TGFb signature, and with the Chen et al. Schwann-cell subtype signatures and every refined subtype in that reference. Because nerves are nested within patients, each correlation is reported pooled, within section, and as the number of patients with a positive correlation; gene-set overlap between the signatures is quantified so the multiple-testing correction can be judged. Produces Figures 6D and S16A-C. |
| `general_functions.R` | Shared utility and plotting functions used across the spatial scripts. |

Paths are relative to the project folder, the one containing `input/`. Run the
scripts from there, or from `scripts/`, which they detect and step up from.

Four files are passed between scripts and treated as fixed inputs rather than
per-analysis results, because many scripts read each of them:
`Integrated_NMF.rds` (`02b`), `nerveID.csv` (`03a`), `DE_SPE_annotaion.csv`
(`04a`) and `RCTD_weights.Rds` (`05b`). The nerve and tumour-distal groupings in
`03a` and `04a` are assigned by clustering, which is deterministic in principle
but not guaranteed identical across package versions, and the nerve numbering in
`nerveID.csv` corresponds to the images that went for pathologist review. These
files are therefore the authoritative versions and the compute scripts that
produce them are not meant to be re-run.

Figure file names give the manuscript figure each output corresponds to. Outputs
named `exploratory_` are not in the manuscript; in `08_draxin_analysis.Rmd`,
outputs named `robustness_` are the second z-scoring scope, run as a robustness
check.

---

### 11. `spatial_enrichment_prioritisation`

Revision analyses that ask whether a signal is genuinely associated with the
nerves rather than with the tissue around them. Every comparison is made against
a geometry-matched control: each nerve's own shape, relocated elsewhere in the
same section, so a nerve is compared with itself displaced. Listed in the order
they are run.

| Script | Description |
|---|---|
| `00_decay_shared.R` | Shared definitions for the near-nerve analyses, sourced by the scripts below: the near and far zones, the curve bins, the three nerve groups, the single decay function every script uses, the pairing of each nerve with its own relocated controls, and the patient-level tests. The settings are stored in the handoff bundle and the later scripts refuse to run on a bundle built with different ones. |
| `00_scale_um_pixel_calibration.R` | Calibrates the pixel-to-micrometre conversion for the coordinates stored in `Integrated_NMF.rds`, using the fixed 100 µm Visium spot spacing: nearest-neighbour spacing is measured within each section and the median of the per-section medians gives the conversion every distance-based script reads. An optional cross-check against the Space Ranger scale factors confirms the coordinate scale is uniform across sections; set `SPACERANGER_DIR` to run it. |
| `01_plasma_decay_analysis_step1.Rmd` | Plasma-cell deconvolution weight as a function of distance from a nerve, for the three nerve categories, against the relocated controls. The three relocated replicates of a nerve are averaged into one control for that nerve, and patients are the replicates: the reported effect is the mean of the patient means with a patient-level paired test, with a nerve-level model blocked on patient as a sensitivity analysis. Writes the distance tables and the bundle the later decay scripts read. |
| `01_plasma_decay_analysis_supplement_step2.R` | True-to-scale zoom of one nerve per category beside its own relocated geometry, each with its distance profile. The example is the nerve whose whole profile is closest to its category's median, chosen among the most typical nerves to be similar in size across categories, so it is typical rather than clearest; a table records each example's percentile within its category. |
| `01_plasma_cell_decay_and_all_cells_step3.Rmd` | The same near-versus-far contrast for every deconvolved cell type, run separately per nerve category, with the same shared function and the same distance tables; a check confirms that the plasma-cell row reproduces the result of step 1 exactly. Correction within each nerve category. Exploratory screen. |
| `02_spatial_association_with_nerve_step1_lee_sham.R` | Spatial association of candidate genes with the invaded nerves: Lee's L between expression and the invaded-nerve mask, computed within each section, calibrated against 200 sham masks in which every invaded nerve's shape is relocated within that section. Real and sham masks are built identically, each with its own contact ring excluded. Per gene and patient the real value is ranked among its shams, and the rank scores are combined across patients by a test enumerated exactly under exchangeability, so strong evidence requires consistency across patients. Defines the core signature, the genes that are also above their shams in all three patients, which the later steps read. Includes a null audit on held-out sham masks. |
| `02_spatial_association_with_nerve_step2_decay_examples.R` | Expression of example genes as a function of distance from the invaded nerves, against the relocated shapes of the same nerves, with a true-to-scale zoom of one typical invaded nerve. Shows the shape of the profile; the test for all genes is step 1. |
| `02_spatial_association_with_nerve_step3_prioritisation_plot.R` | Places every gene upregulated in invaded nerves by its association with the nerve mask and with the interface ring, each as the difference from the mean of its relocated masks, and marks the core signature. Top genes of the tumour and nerve NMF factors are shown for reference and not used for selection. |
| `02_spatial_association_with_nerve_step4_celltype_origin.R` | Likely cellular source of each signature gene in the Chen et al. reference: the cell type with the highest patient-averaged fraction of expressing cells, tested by pseudo-bulk limma-voom with patient as a blocking factor, against the average of the other cell types and against the next-highest one. Genes detected in fewer than 5% of cells in every cell type are reported as not assignable. |
| `02_spatial_association_with_nerve_step5_ORA.R` | Over-representation analysis of the core signature, queried as one list, against two backgrounds: the genes upregulated in invaded nerves, the set the signature was selected from, and all genes tested in the differential expression. |
| `02_spatial_association_with_nerve_step6_robustness_InvN_vs_NITB.R` | Robustness check, not used in the manuscript. The manuscript reports invaded nerves against all non-invaded nerves; this script repeats the pseudo-bulk differential expression against non-invaded nerves of the tumour-bearing sections only, so that both groups lie in the same sections and the contrast cannot be carried by the difference between the tumour-bearing and the adjacent block. The test is genome-wide and the prioritised genes are extracted from the finished result, so no gene set is redefined. |

Paths are relative to the project folder, the one containing `input/`. Run the
scripts from there, or from `scripts/`, which they detect and step up from.

`step6` reads `DE_analysis_results.rds`, written by
`04a_DE_analysis_COMPUTE.Rmd` in `spatial_data_processing`; place a copy in
`input/`.

Outputs named `robustness_` are not in the manuscript. Figure file names
describe what each figure shows rather than its number in the paper.

---

### 12. `statistics`

Statistical analysis scripts for spatial transcriptomics, deconvolution outputs, DRAXIN-associated phenotypes, cell growth, wound closure, migration, attachment, flow cytometry, EdU incorporation, ROS assays, and related functional experiments.

| Script | Description |
|---|---|
| `08_Draxin_DE_analysis_2026_04_28_stats.Rmd` | Statistical analysis of nerve-level DRAXIN expression and associated transcriptional programs across NI and InvN categories using patient-adjusted linear models and estimated marginal means. |
| `22_stats_figure2.Rmd` | Statistical analysis of RCTD-derived cell-type proportions and NMF programme activities across NI and InvN nerve categories using patient-adjusted linear modelling. |
| `Cell_growth_BXPC3_tecan.Rmd` | Statistical analysis of DRAXIN-associated growth effects in BXPC3 cells using linear mixed-effects models and estimated marginal means. |
| `Cell_growth_other_PDAC_lines_tecan.Rmd` | Cell-growth analysis across multiple PDAC cell lines under DRAXIN treatment using mixed-effects modelling. |
| `Coating_preincubation_tecan.Rmd` | Analysis of coating and pre-incubation effects on cell-growth phenotypes using linear modelling approaches. |
| `Commerical_boiled_self_made.Rmd` | Comparative statistical analysis of recombinant DRAXIN preparations across experimental conditions. |
| `Cumulative_absorbance_22aa_peptide.Rmd` | Linear mixed-effects modelling of cumulative absorbance assays testing DRAXIN and 22-aa peptide interactions. |
| `Cumulative_absorbance_Draxin_cell_lines.Rmd` | Analysis of DRAXIN-associated cumulative absorbance measurements across multiple PDAC cell lines. |
| `Cumulative_absorbance_Draxin.Rmd` | Dose-response analysis of cumulative crystal violet absorbance under DRAXIN treatment conditions. |
| `Flow_apoptosis_necrosis_day3.Rmd` | Flow-cytometry-based analysis of apoptosis and necrosis measurements under DRAXIN treatment conditions. |
| `Red_dot2_percentage_tecan.Rmd` | Statistical analysis of RedDot2-positive cell proportions across treatment conditions. |
| `Single_cell_speed_coating.Rmd` | Analysis of single-cell migration speed across coating conditions using linear modelling and post hoc comparisons. |
| `Timing_of_addition_tecan.Rmd` | Statistical evaluation of treatment timing effects on cell-growth phenotypes. |
| `Wound_closure_15_screen.Rmd` | Linear modelling analysis of wound-healing assays across recombinant screening conditions. |
| `Wound_closure_DRAXIN_gradient.Rmd` | Dose-response analysis of wound closure velocity under increasing DRAXIN concentrations. |
| `Wound_closure_mitomcyin.Rmd` | Analysis of wound-healing assays incorporating Mitomycin treatment and DRAXIN interaction effects. |
