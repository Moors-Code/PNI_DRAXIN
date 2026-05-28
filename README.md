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

Workflows for Cas12a pooled-screen quality control, guide-count normalization, MAGeCK-compatible filtering, and standard MAGeCK reporting for the selected T0 + ctrl filtering strategy.

| Script | Description |
|---|---|
| `01_libraryQC.Rmd` | Quality control and normalization of Cas12a guide-count libraries using CPM and variance-stabilized transformed counts. Includes replicate correlation analysis, PCA, guide variability assessment, and library composition summaries. |
| `02_ScreenQC.Rmd` | Quality control analysis of pooled CRISPR screen count matrices using raw and normalized sgRNA counts. Includes detected sgRNA metrics, skew ratios, PCA, MDS, density plots, correlation heatmaps, scatterplots, and MD plots. |
| `03_MAGECK_pre_filtering.Rmd` | Generation of MAGeCK-compatible filtered count matrices using voom-normalized logCPM thresholds across selected reference groups, including library input, T0, and ctrl samples. |
| `FILTERING_T0_PLUS_CONTROL_PAIRED_TOTAL.report.Rmd` | Standard MAGeCK report summarizing CRISPR screen enrichment results from the selected T0 + ctrl filtering strategy, including sgRNA- and gene-level enrichment statistics and quality assessment plots. |

---

### 3. `enrichment_analysis`

Pathway and gene-set analysis workflows integrating bulk RNA-seq, CRISPR screen outputs, invaded nerve signatures, ssGSEA scoring, GSEA, ORA, and pathway-level visualizations.

| Script | Description |
|---|---|
| `bulkRNAseq_01_GSEA.Rmd` | Automated edgeR-to-GSEA workflow performing GO, Reactome, KEGG, Disease Ontology, and MSigDB enrichment analyses across multiple differential expression comparisons. |
| `bulkRNAseq_02_ssGSEA.Rmd` | ssGSEA pathway scoring of batch-corrected logCPM data across ctrl, 17 nM DRAXIN, and 68 nM DRAXIN conditions with pathway-level statistical comparison using limma. |
| `bulkRNAseq_ENSEMBL_to_GENE_SYMBOL.R` | Annotation of differential expression tables by mapping Ensembl IDs to HGNC gene symbols using `org.Hs.eg.db`. |
| `CRISPR_GSEA_p_Value_based_scoring.Rmd` | GSEA analysis of MAGeCK gene-level CRISPR screen results using separate ranking strategies for resistance and sensitisation arms across Hallmark, C2, C5, and Reactome gene sets. |
| `CRISPR_GSEA_plotting_for_paper.Rmd` | Visualization and filtering of MAGeCK-GSEA results using Hallmark and C5 pathways with positive NES and adjusted p-value thresholds. |
| `figure2_GSEA_inv_vs_uninv.Rmd` | GSEA analysis comparing invaded nerves and non-invaded nerve differential expression profiles using ontology, pathway, and custom GMT-based enrichment workflows. |
| `figure2_ssGSEA_score_Rmd` | Pathway activity profiling across NI nerves from adjacent healthy tissue, NI nerves from tumour-bearing tissue, and invaded nerves using ssGSEA-derived pathway scores. |
| `figure3_enricher_342_invasion_genes.Rmd` | Over-representation analysis of invasion-associated gene signatures using a defined differential-expression background universe. |

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
