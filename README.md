# Tumour-Conditioned Nerves Establish a DRAXIN+ Niche Associated with Metabolic Adaptation in Pancreatic Cancer

This repository contains the computational workflows and statistical analyses used in the study:

**“Tumour-Conditioned Nerves Establish a DRAXIN+ Niche Associated with Metabolic Adaptation in Pancreatic Cancer”**

## Repository Structure

### 1. `bulk_RNA_seq`

Scripts for bulk RNA-seq count processing, normalization, batch-aware differential expression analysis, exploratory quality control, and reusable RNA-seq plotting and modelling functions.

| Script | Description |
|---|---|
| `count_merging.R` | Merges RNA-seq count matrices and extracts gene-length information for downstream RPKM normalization and transcript-level quantification. |
| `DEG_modeling_experiment_libraryprep_2026_02_26.R` | Batch-corrected edgeR quasi-likelihood differential expression analysis comparing ctrl, DRAXIN, and diluted DRAXIN conditions across independent experiments. Includes TMM normalization, dispersion estimation, PCA, volcano plots, and differential expression testing. |
| `functions.R` | Custom utility functions for RNA-seq analysis, including edgeR/DESeq2 differential expression workflows, normalization, PCA, volcano plotting, correlation analysis, GSEA visualization, and sample-matching utilities. |

---

### 2. `Cas12_screen`

Workflows for Cas12a pooled-screen quality control, guide-count assembly, MAGeCK-compatible filtering and
reporting for the Cas12a pooled screen (DRAXIN-treated versus control).
Listed in the order they are run.

| Script | Description |
|---|---|
| `01_libraryQC.Rmd` | Quality control of the Cas12a guide-count libraries. Library-size normalization (CPM and log2(CPM+1)), Spearman correlation between replicates, guide-level variability (CV versus abundance), PCA on variance-stabilized counts (DESeq2 VST), per-replicate distribution metrics including Gini index, and comparison against the plasmid library stock. |
| `02a_merging_guide_count.R` | Assembles the analysis count matrix, combining the per-condition screen count files with the library QC counts into a single guide × sample table. Guides are matched by identifier rather than by row order, which corrects the merging error present in the earlier version of this step. |
| `02_ScreenQC.Rmd` | Quality control of the merged screen count matrix. Detected sgRNAs and skew ratio per sample, CPM and TMM-normalized logCPM, density and histogram distributions, PCA, MDS, Pearson correlation heatmap, pairwise scatterplots, MD plots, and control checks comparing common-essential and non-targeting guides against the remaining library. |
| `03_MAGECK_pre_filtering.Rmd` | Generates MAGeCK input matrices. Guides are retained on mean voom logCPM (cutoff 2.5) computed in a chosen set of reference samples; three strategies are run (library input + T0, T0 only, T0 + ctrl), each writing raw counts for MAGeCK together with the per-guide filtering decisions. |
| `run_mageck_all_variants.sh` | Runs MAGeCK on the screen: a paired test of the two DRAXIN-treated against the two control replicates, repeated across seven combinations of input matrix (unfiltered, T0-filtered, T0 + ctrl filtered) and normalization method (non-targeting control, median, total). The analysis reported in the paper is `T0_PLUS_CONTROL_PAIRED_TOTAL`; the remaining variants are included to show that the reported hits do not depend on the filtering or normalization choice. Run on an HPC cluster (SLURM) with MAGeCK version `[VERSION]`. |
| `04_MAGECK_figures_output.R` | Draws the published screen figures from the MAGeCK gene summary: ranked RRA score distributions for enriched and depleted arrays, and the volcano of −log10 FDR against log2 fold change. MAGeCK is not re-run, so the plotted values are those of the original analysis. |
| `05_MAGECK.report.Rmd` | Standard MAGeCK report for the T0 + ctrl paired comparison, giving sgRNA- and gene-level statistics, RRA and p-value distributions, per-guide read counts for top genes, and essential-gene dropout by GSEA. Adapted from the report template distributed with MAGeCK (W. Li). |

Paths in `run_mageck_all_variants.sh` are placeholders and must be set before
running. Steps `04` and `05` read the `T0_PLUS_CONTROL_PAIRED_TOTAL` output.
---

### 3. `enrichment_analysis`

Pathway and gene-set analysis across bulk RNA-seq, the CRISPR screen and the
spatial nerve data. Each analysis is split into a compute step, which saves
result objects, and a plot step, which reads them — so figures can be
regenerated without re-running GSEA.

| Script | Description |
|---|---|
| `00_annotate_gene_symbols.R` | Adds HGNC symbols to the edgeR differential expression tables by mapping Ensembl IDs with `org.Hs.eg.db`. |
| `01a_GSEA_compute.R`, `01b_GSEA_plots.R` | GSEA on the bulk RNA-seq comparisons (DRAXIN vs ctrl, diluted vs ctrl, DRAXIN vs diluted), ranked by signed −log10 p. Covers GO (BP, CC, MF), Reactome, Disease Ontology, KEGG and MSigDB Hallmark, followed by a pass across all MSigDB collections. The plot step draws dotplots split by NES sign and writes gene-symbol result tables. |
| `01c_ssGSEA_compute.R`, `01d_ssGSEA_plots.R` | ssGSEA pathway scoring of batch-corrected logCPM across 0, 17 and 68 nM DRAXIN, for a fixed panel of Hallmark, GO and KEGG gene sets. Scores are standardized both across all samples and relative to the control group; the plot step fits limma on the raw scores to test the three dose contrasts and annotates the panels with FDR. |
| `01e_heatmap_compute.R`, `01f_heatmap_plots.R` | Leading-edge gene heatmaps for selected pathways. The compute step takes the union of leading-edge genes across the three comparisons and maps them onto the batch-corrected expression matrix; the plot step draws one row-z-scored heatmap per pathway, split by dose. |
| `02a_CRISPR_GSEA_compute.R`, `02b_CRISPR_GSEA_plots.R` | GSEA on the MAGeCK gene summary, ranked by LFC × −log10 p, with the resistance and sensitisation arms ranked separately, across Hallmark, C2, C5 and Reactome. Multi-target guides, dedup-collapsed gene pairs and pooled control pseudo-genes are excluded from ranking and exported separately. The plot step produces the paper's NES comparison figure and the full results table. |
| `02c_CRISPR_GSEA_exploration_supplementary.R` | Broader exploration of the CRISPR GSEA results: both filtering modes, all collections, per-result dotplots and a paired resistance-versus-sensitisation NES scatter. Not part of the reported figures. |
| `03a_Figure2_GSEA_compute.R`, `03b_Figure2_GSEA_plots.R` | GSEA on the invaded versus non-invaded nerve differential expression, ranked by signed −log10 FDR. Covers GO, Reactome, WikiPathways, Disease Ontology, Network of Cancer Genes, KEGG, all MSigDB collections, and a custom cell-state signature set supplied as a GMT. |
| `03c_Figure2_ssGSEA_compute.R`, `03d_Figure2_ssGSEA_plots.R` | Pathway activity across NI-Adj, NI-TB and InvN nerves. Spot-level ssGSEA scores for an immune and complement pathway panel are aggregated to the nerve level and z-scored across nerves; the plot step draws the two Figure 2 panels. |
| `03e_Figure2_ssGSEA_stats_supplementary.R` | Repeated-measures ANOVA across nerve categories on the patient-aggregated pathway scores. Not used in the paper. |
| `combine_GSEA_results.Rmd` | Concatenates the per-collection GSEA result tables from the invaded-nerve analysis into one merged file for downstream pathway selection. |

`03c_Figure2_ssGSEA_compute.R` reads the spatial `SpatialExperiment` object
and its annotation from a separate project directory (`external_spatial_dir`
at the top of the script), not from this repository's `input/` folder, as the
spatial data are large and shared across projects. Set that path to wherever
the spatial data are located, or to the deposited dataset, before running it.

---

### 4. `IF_analysis`

Image-analysis workflows for manual annotation consensus, immune-cell spatial profiling, RNAscope marker correlations, nerve-associated immune quantification, and publication-ready visualization.

| Script | Description |
|---|---|
| `Classification_Plots.Rmd` | Integration of manual nerve and tumour annotations from multiple reviewers to generate consensus classifications of surrounding P75NTR+ line structures across invaded and non-invaded nerves. |
| `immune_cell.Rmd` | Spatial immune profiling across CD45+ immune aggregates positioned at increasing distances from nerves to quantify immune composition, inflammatory cytokines, and immune checkpoint activity. |
| `P75NGF_vs_Draxin_RNAscope_correlation.R` | Correlation analysis between DRAXIN puncta per cell and proportions of P75NGFR+ or PS100B+ cells per nerve using Spearman statistics and scatterplot visualization. |

---

### 5. `spatial_data_processing`

Spatial transcriptomics workflows covering SpatialExperiment construction, quality control, NMF decomposition, nerve and tumour-region annotation, pseudo-bulk differential expression, sc/snRNA-seq reference preprocessing, RCTD deconvolution, DRAXIN-focused spatial analyses, LIANA ligand–receptor analysis, plasma-cell mapping, and selected-gene spatial visualization.

| Script | Description |
|---|---|
| `00_Making.SPE.R` | Generation of integrated SpatialExperiment objects from multiple 10x Visium datasets, including raw counts, spatial coordinates, and tissue images. |
| `01_Quality_Control.Rmd` | Spatial transcriptomics preprocessing and quality control across paired tumour-bearing tissue and adjacent healthy tissue samples, including spot filtering, QC metrics, and normalization. |
| `02_NMF_2026_03_02.Rmd` | Sparse non-negative matrix factorization analysis for identification of reference-free transcriptional programs across spatial transcriptomic spots. |
| `03_Nerve_ID_2026_03_02.Rmd` | Visualization and interpretation of NMF-derived transcriptional programs, including factor representation plots, gene loading plots, and spatial activity maps. |
| `04_Differenital_Expression_2026_03_12.Rmd` | Spatial differential expression analysis comparing InvN, NI, nerve-proximal tumour spots, and nerve-distal tumour spots using pseudo-bulk aggregation and limma-based modelling. |
| `05_Reference_preprocessing.Rmd` | Preprocessing and BBKNN integration of published scRNA-seq and snRNA-seq reference datasets for downstream spatial deconvolution. |
| `06_Deconvolution_RCTD.Rmd` | RCTD-based deconvolution of spatial transcriptomics data using integrated scRNA-seq/snRNA-seq references to estimate spot-level cell type composition. |
| `06_Deconvolution_RCTD_2026_04_27.Rmd` | Nerve-focused deconvolution analysis comparing RCTD-derived immune and macrophage-associated programs across three nerve categories. |
| `08_Draxin_DE_analysis_2026_02_22_modified.Rmd` | Analysis of DRAXIN expression and Schwann-associated transcriptional programs across invaded and non-invaded nerve categories. |
| `09_LR_analysis.ipynb` | Spatial ligand–receptor interaction analysis using LIANA to identify signalling interactions associated with nerve invasion based on local spatial covariation patterns. |
| `10_Gene_expression_for_figure2_and_3.Rmd` | Patient-averaged spatial gene expression profiling across nerve categories using curated immune, stromal, invasion-associated, and tumour-associated gene panels. |
| `12_Draxin_correlattion_with_NMF_nerve_tumor_TGFbeta.Rmd` | Correlation analysis between DRAXIN expression, nerve-associated and tumour-associated NMF programs, and TGFβ-response signatures across tumour-associated nerves. |
| `13_Plasma_cells.Rmd` | Spatial analysis of plasma cell enrichment using RCTD-derived plasma-cell weights across invaded and non-invaded nerve microenvironments. |
| `14_selected_genes_spatial.Rmd` | Spatial visualization of BTC, FGF5, and IL11 expression across tumour-bearing tissue and adjacent healthy tissue sections containing invaded and non-invaded nerves. |
| `general_functions.R` | Custom utility and visualization functions supporting clustering, mapping, dimensionality reduction, graph extraction, and publication-quality plotting workflows for spatial transcriptomics analyses. |

---

### 6. `statistics`

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
