# This is the file for functions I created for Minkyoung's data
# Source this in your code source("functions.R") to be able to use the
# functions

# This function plots the whole data with a scatter plot.
plot_correlation_all <- function(data, highlight_data, x, y, 
                                 pearson_x = 0.3, pearson_y = 16,
                                 title = "") {
  ggplot(
    data = data,
    aes_string(
      x = x,
      y = y
    ),
  ) +
    # scale_x_continuous(
    #   breaks = c(5, 10, 15)
    # ) +
    # scale_y_continuous(
    #   breaks = c(5, 10, 15)
    # ) +
    geom_point(
      shape = 16,
      alpha = 0.4,
      color = "lightblue"
    ) +
    geom_smooth(method = lm) +
    stat_cor(method = "pearson",
             label.y = max(c(
               data[, x],
               data[, y]
             )) - 0.7
             ) +
    geom_point(
      data = highlight_data,
      aes_string(
        x = x,
        y = y
      ),
      color = "darkblue",
      size = 3,
      alpha = 0.7
    ) +
    geom_text_repel(data = highlight_data, 
                    aes(label=rownames(highlight_data)), size=5,
                    nudge_y = 0.2, force = 2
    ) +
    geom_text(
      x = pearson_x, y = pearson_y, label = "Pearson",
      check_overlap = T
    ) +
    theme_classic() +
    theme(
      axis.title = element_text(size = 20),
      axis.text = element_text(size = 17)
    ) +
    ggtitle(title)
}

# This function puts limit to x and y axes to zoom into certain parts
# CAUTION: It will filter out the other data points, therefore, your
# correlation values will be different than plotting all data points
plot_correlation_zoomed <- function(data, highlight_data, x, y, low_limit = 4,
                                    up_limit = 15, pearson_x = 4.3, pearson_y = 15.2) {
  ggplot(
    data = data,
    aes_string(
      x = x,
      y = y
    ),
  ) +
    scale_x_continuous(
      breaks = c(5, 10, 15),
      limits = c(low_limit, up_limit)
    ) +
    scale_y_continuous(
      breaks = c(5, 10, 15),
      limits = c(low_limit, up_limit)
    ) +
    geom_point(
      shape = 16,
      alpha = 0.4,
      color = "lightblue"
    ) +
    geom_smooth(method = lm) +
    stat_cor(method = "pearson",
             label.y = max(c(
               data[, x],
               data[, y]
             )) - 0.7
             ) +
    geom_point(
      data = highlight_data,
      aes_string(
        x = x,
        y = y
      ),
      color = "darkblue",
      size = 3,
      alpha = 0.7
    ) +
    geom_text(
      data = highlight_data,
      aes(label = rownames(highlight_data)),
      hjust = 0.5, vjust = -0.6,
      cex = 6
    ) +
    geom_text(
      x = pearson_x, y = pearson_y, label = "Pearson",
      check_overlap = T
    ) +
    theme_classic() +
    theme(
      axis.title = element_text(size = 20),
      axis.text = element_text(size = 17)
    )
}

# This plot is same with the first, but it doesn't have a highlight.
plot_correlation_no_highlight <- function(data, x, y, pearson_x = 0.3, pearson_y = 16) {
  ggplot(
    data = data,
    aes_string(
      x = x,
      y = y
    ),
  ) +
    # scale_x_continuous(
    #   breaks = c(5, 10, 15)
    # ) +
    # scale_y_continuous(
    #   breaks = c(5, 10, 15)
    # ) +
    geom_point(
      shape = 16,
      alpha = 0.4,
      color = "lightblue"
    ) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
    geom_smooth(method = lm) +
    stat_cor(
      method = "pearson",
      label.y = max(c(
        data[, x],
        data[, y]
      )) - 0.7
    ) +
    geom_text(
      x = pearson_x, y = pearson_y, label = "Pearson",
      check_overlap = T
    ) +
    theme_classic() +
    theme(
      axis.title = element_text(size = 20),
      axis.text = element_text(size = 17)
    )
}

# This plot is the same with the second, but without highlighting
plot_correlation_zoomed_no_highlight <- function(data, x, y, low_limit = 4,
                                                 up_limit = 15, pearson_x = 4.3, pearson_y = 15.2) {
  ggplot(
    data = data,
    aes_string(
      x = x,
      y = y
    ),
  ) +
    scale_x_continuous(
      breaks = c(5, 10, 15),
      limits = c(low_limit, up_limit)
    ) +
    scale_y_continuous(
      breaks = c(5, 10, 15),
      limits = c(low_limit, up_limit)
    ) +
    geom_point(
      shape = 16,
      alpha = 0.4,
      color = "lightblue"
    ) +
    geom_smooth(method = lm) +
    stat_cor(
      method = "pearson",
      label.y = max(c(
        data[, x],
        data[, y]
      )) - 0.7
    ) +
    geom_text(
      x = pearson_x, y = pearson_y, label = "Pearson",
      check_overlap = T
    ) +
    theme_classic() +
    theme(
      axis.title = element_text(size = 20),
      axis.text = element_text(size = 17)
    )
}

#For plotting samples against each other, you sometimes need to plot all
#samples in pairs. This function loops over the list, and returns a data frame
#with two columns. You can use this data frame to loop over the columns
#in the original data frame to quickly plot everything
make_pairs <- function(list_of_samples) {
  comparisons <- data.frame(matrix(nrow = 0, ncol = 2))
  k <- 1
  while (k <= length(list_of_samples)) {
    i <- k + 1
    while (i <= length(list_of_samples)) {
      comparisons <- rbind(comparisons, 
                           data.frame(X1 = list_of_samples[k], 
                                      X2 = list_of_samples[i]))
      i <- i + 1
    }
    k <- k + 1
  }
  return(comparisons)
}

# For extracting matching orgs and ctrls from the dataframe
extract_matching_samples <- function(dataframe, sample_ID) {
  orgs <- colnames(dataframe) %>%
    grep(pattern = paste0(sample_ID,
      "_org[0-9]{1,2}"),
      value = T
    )
  orgctrls <- colnames(dataframe) %>% 
    grep(pattern = paste0(sample_ID,
                           "_orgctrl[0-9]{1,2}"), 
                           value = T)
  orgctrls <- orgctrls[gsub("ctrl", "", orgctrls) %in% orgs]
  orgs <- orgs[orgs %in% gsub("ctrl", "", orgctrls)]
  matching_frame <- data.frame(orgctrls = sort(orgctrls), 
                               orgs = sort(orgs))
  return(matching_frame)
}

#Creates gene list with reshaped data frame for easy dot / jitter
#plotting. Accepts a vector of gene names aside from the data frame
#of the expression values. If there is a control sample, the other
#condition will be divided by the median of the control to normalize
#the data. This can be disabled by setting median_normalize_first_label
#to false
prepare_dot_plot_data_case_control <- function(dataframe, gene_vector, 
                                               sample_ID,
                                  condition_label=c("Control","H2O2"),
                                  median_normalize_first_label = T) {
  collection_of_genes <- data.frame(matrix(
    nrow = 0,
    ncol = 4,
    dimnames = list(
      NULL,
      c(
        "variable",
        "value",
        "Condition",
        "Gene"
      )
    )
  ))

  for (i in gene_vector) {
    only_gene_df <- dataframe[i, ]
    list_of_samples <- extract_matching_samples(only_gene_df, sample_ID)
    selected_sample_gene <- only_gene_df[, colnames(only_gene_df) %in%
      as.vector(unlist(list_of_samples))]
    melted_table <- melt(selected_sample_gene) %>%
      dplyr::mutate(Condition = ifelse(grepl("ctrl", variable),
        condition_label[1], condition_label[2]
      ))
    
    # Calculate control median for the plot
    if(median_normalize_first_label){
    the_median <- median(
      melted_table[melted_table$Condition == condition_label[1], 
                   "value"])
    }
    
    # Divide the values by median
    melted_table$value <- melted_table$value / the_median
    melted_table$Gene <- i
    collection_of_genes <- rbind(collection_of_genes, melted_table)
  }
  
  return(collection_of_genes)
}



# This is for median of ratios normalization, which is one of the
# best performing normalization method for RNASeq data for between
# sample comparison. This is explained very well in DESeq2's website,
# but this website had this function ready:
# https://scienceparkstudygroup.github.io/rna-seq-lesson/median_of_ratios_manual_normalization/index.html
mor_normalization = function(data){
  library(dplyr)
  library(tibble)
  
  # take the log
  log_data = log(data) 
  
  # find the psuedo-references per sample by taking the geometric mean
  log_data = log_data %>% 
    rownames_to_column('gene') %>% 
    mutate (gene_averages = rowMeans(log_data)) %>% 
    filter(gene_averages != "-Inf")
  
  # the last columns is the pseudo-reference column 
  pseudo_column = ncol(log_data)
  
  # where to stop before the pseudo column 
  before_pseduo = pseudo_column - 1
  
  # find the ratio of the log data to the pseudo-reference
  ratios = sweep(log_data[,2:before_pseduo], 1, log_data[,pseudo_column], "-")
  
  # find the median of the ratios
  sample_medians = apply(ratios, 2, median)
  
  # convert the median to a scaling factor
  scaling_factors = exp(sample_medians)
  
  # use scaling factors to scale the original data
  manually_normalized = sweep(data, 2, scaling_factors, "/")
  return(manually_normalized)
}

# This function simply gets the un-normalized count data for RNA Seq
# and applies DESeq analysis on it. This can be found in DESeq2 vignette.
get_DEG <- function(count_table, batch_info, design, result_name){
  library("DESeq2")
  dds <- DESeqDataSetFromMatrix(countData = count_table,
                                colData = batch_info,
                                design= design)
  dds <- DESeq(dds)
  res <- results(dds, name=result_name)
  return(res)
}


# This function is to use get_DEG, and extract the up-regulated genes
# compared to a base level condition. In my case, I compared conditions
# with GFP since GFP should be expressed everywhere, and the other
# conditions should be focused in one region (and expressed more).
get_upreg_genes <- function(selected_data_df,
                            grep_pattern, 
                            batch_df,
                            design,
                            base_condition,
                            result_name_DE,
                            only_padj=T){
  
  #I select samples
  dpp_vs_gfp_selected <- selected_data_df[,
                                          grep(grep_pattern,
                                               colnames(selected_data_df))]
  
  #And get the relevant batch info
  dpp_gfp_batch_selected <- batch_df[colnames(dpp_vs_gfp_selected),]
  
  #Add the column for comparison
  dpp_gfp_batch_selected$compare <- dpp_gfp_batch_selected$samples %>%
    stringr::str_split(pattern =  "_" ) %>% purrr::map(,.f = 1) %>%
    unlist() %>% as.factor()
  
  
  #I make the GFP 1 so that we see log change as if it is high in DPP or not
  dpp_gfp_batch_selected$compare <- relevel(dpp_gfp_batch_selected$compare, 
                                            base_condition)
  
  #Test between all GFP and DPP4
  dpp_gfp_selected_results <- get_DEG(count_table = dpp_vs_gfp_selected,
                                      batch_info = dpp_gfp_batch_selected,
                                      design = design,
                                      result_name = result_name_DE)
  
  #Get the results in data frame (deletes the meta data)
  dpp_gfp_selected_results <- data.frame(dpp_gfp_selected_results)
  
  #Get only the upregulated ones
  if(only_padj){
  upreg_dpp_gfp_selected <- dpp_gfp_selected_results[
    # dpp_gfp_selected_results$log2FoldChange>0 & 
      dpp_gfp_selected_results$padj < 0.05,
  ]
  
  #Remove NAs
  upreg_dpp_gfp_selected <- upreg_dpp_gfp_selected[
    !is.na(upreg_dpp_gfp_selected$log2FoldChange),]
  
  return(upreg_dpp_gfp_selected)
  } else {
    dpp_gfp_selected_results <- dpp_gfp_selected_results[
      !is.na(dpp_gfp_selected_results$log2FoldChange),]
    return(dpp_gfp_selected_results)
  }
  

}


# This is a function to create a Venn diagram to see the intersection
# between two data frames with genes.
compare_with_a_list <- function(df_of_genes,
                                df_to_compare,
                                venn_out_path,
                                merge_col_x = 0,
                                merge_col_y = 3){

  #Compare with Andreas' results
  andreas_intersect <- df_to_compare[df_to_compare$ext_gene %in% 
                                      rownames(df_of_genes),
  ]
  
  # # I want to convert Minkyoung's list to ensemble gene ID to make sure
  # # gene symbols are not a problem
  # 
  # ensembl = useMart("ensembl",dataset="mmusculus_gene_ensembl")
  # 
  # #I get the list of ensembl ids from biomart
  # ensembl_mk_list <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
  #       filters = "external_gene_name",
  #       values = rownames(upreg_dpp_gfp_selected),
  #       mart = ensembl)
  # 
  # #I see not all genes have been converted... So I check which genes are missing
  # missing_genes <- rownames(upreg_dpp_gfp_selected)[!(
  #   rownames(upreg_dpp_gfp_selected) %in% 
  #     ensembl_mk_list$external_gene_name)]
  # 
  # #I get the ensembl ids anyway to see if we have more than 15 intersection
  
  #Ensembl ID conversion made no difference so I continue with gene symbols
  
  venn.diagram(list(df_to_compare$ext_gene, rownames(df_of_genes)),
               filename = venn_out_path,
               category.names = c(paste0(
                 "Andreas (", nrow(df_to_compare),")") , paste0(
                 "Minkyoung (", nrow(df_of_genes),")")),
               output = TRUE ,
               disable.logging = T,
               imagetype="png",
               height = 1080 ,
               width = 1080 ,
               resolution = 600,
               compression = "lzw",
               lwd = 1,
               col=c("#440154ff", '#21908dff'),
               fill = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3)),
               cex = 0.5,
               fontfamily = "sans",
               cat.cex = 0.3,
               cat.default.pos = "outer",
               cat.pos = c(-20, 12),
               cat.dist = c(0.055, 0.055),
               cat.fontfamily = "sans",
               cat.col = c("#440154ff", '#21908dff')
               # rotation = 1
  )
  df_of_genes <- merge(df_of_genes, 
                       andreas_intersect, by.x=merge_col_x,
                       by.y = merge_col_y, all.x =T)
  
 return(df_of_genes)
}


# This is a function to get PCA and PCA plot done quickly
plot_pca <- function(data, colour_column = group, ctrl=F){

library(patchwork)
library(stringr)
data_t <- as.data.frame(t(data))
chosen_pca_obj <- prcomp(data_t, scale. = T, center = T, rank. = 20)
chosen_pca <- as.data.frame(chosen_pca_obj$x)
# And I add the group (as in GFP, DPP4, INTB4 e tc.)
chosen_pca$group <- str_split(rownames(chosen_pca),"\\.",simplify = T)[,2]
if(ctrl == T){
  controls <- grep(pattern = "c$",x =rownames(chosen_pca))
  chosen_pca$group[controls] <- paste0(chosen_pca$group[controls],"-c")
  nb.cols <- length(unique(chosen_pca$group))
  mycolors <- colorRampPalette(brewer.pal(24, "Dark2"))(nb.cols)
} else {
  nb.cols <- length(unique(chosen_pca$group))
  mycolors <- colorRampPalette(brewer.pal(24, "Dark2"))(nb.cols)
}

expl_variance <- fviz_eig(chosen_pca_obj,
                          ggtheme = theme(axis.title = element_text(size=24),
                                          axis.text = element_text(size=20),
                                          legend.text = element_text(size=20)))
labels_text <- paste0(str_split(
  rownames(chosen_pca),
  "\\.", simplify = T)[,1],"-" ,
  str_split(
    rownames(chosen_pca),
    "\\.", simplify = T)[,2],"-",
  str_split(
    rownames(chosen_pca),
    "\\.", simplify = T)[,3])
pca_plot1 <- ggplot(data = chosen_pca, aes(x=PC1, y=PC2, 
                              colour=group,
                              label = labels_text)) +
  geom_point(size=4) + geom_text_repel(size =7, show.legend = F) + 
  theme_light() + scale_color_manual(values = mycolors,
                                     name="Conditions") +
  theme(axis.title = element_text(size=24),
        axis.text = element_text(size=20),
        legend.text = element_text(size=20),
        legend.title = element_text(size=24))

# Normally, scale_color_manual(values = mycolors)
# But I edit for Minkyoung's request

pca_plot2 <- ggplot(data = chosen_pca, aes(x=PC1, y=PC3, 
                                           colour=group,
                                           label = labels_text)) +
  geom_point(size=4) + geom_text_repel(size =7, show.legend = F) + 
  theme_light() + 
  #scale_color_manual(values = c(
   # "#df6020", "#58a77c" , "#c39b3c" , "#ff0040","#6d32cd"
  #),
  scale_color_manual(values = mycolors,
name="Conditions") +
  theme(axis.title = element_text(size=24),
        axis.text = element_text(size=20),
        legend.text = element_text(size=20))
 
pca_plot3 <- ggplot(data = chosen_pca, aes(x=PC2, y=PC3, 
                                           colour=group,
                                           label = labels_text)) +
  geom_point(size=4) + geom_text_repel(size =7, show.legend = F) + 
  theme_light() + scale_color_manual(values = mycolors,
                                     name="Conditions")+
  theme(axis.title = element_text(size=24),
        axis.text = element_text(size=20),
        legend.text = element_text(size=20))

pca_plot <- pca_plot1 + pca_plot2 + pca_plot3

return(list(chosen_pca,expl_variance,pca_plot,
            pca_plot1,pca_plot2,pca_plot3))
}

# This is a function to get PCA and PCA plot done quickly
plot_pca_merged <- function(data, colour_column = group, ctrl=F){
  
  library(patchwork)
  library(stringr)
  data_t <- as.data.frame(t(data))
  chosen_pca_obj <- prcomp(data_t, scale. = T, center = T, rank. = 20)
  chosen_pca <- as.data.frame(chosen_pca_obj$x)
  # And I add the group (as in GFP, DPP4, INTB4 e tc.)
  groups <- str_split(rownames(chosen_pca),"_",simplify = T)
  chosen_pca$group <- paste(groups[,1], groups[,ncol(groups)],sep = "_")
  if(ctrl == T){
    controls <- grep(pattern = "c_",x =rownames(chosen_pca))
    chosen_pca$group[controls] <- paste0(chosen_pca$group[controls],"-c")
    nb.cols <- length(unique(chosen_pca$group))
    mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols)
  } else {
    nb.cols <- length(unique(chosen_pca$group))
    mycolors <- colorRampPalette(brewer.pal(8, "Dark2"))(nb.cols)
  }
  
  expl_variance <- fviz_eig(chosen_pca_obj,
                            ggtheme = theme(axis.title = element_text(size=24),
                                            axis.text = element_text(size=20),
                                            legend.text = element_text(size=20)))
  labels_text <- paste0(str_split(
    rownames(chosen_pca),
    "_", simplify = T)[,1],"-" ,str_split(
      rownames(chosen_pca),
      "_", simplify = T)[,2])
  pca_plot1 <- ggplot(data = chosen_pca, aes(x=PC1, y=PC2, 
                                             colour=group,
                                             label = labels_text)) +
    geom_point(size=4) + geom_text_repel(size =7, show.legend = F,
                                         vjust=-0.35, hjust=-0.7) + 
    theme_light() + scale_color_manual(values = mycolors,
                                       name="Conditions") +
    theme(axis.title = element_text(size=24),
          axis.text = element_text(size=20),
          legend.text = element_text(size=20),
          legend.title = element_text(size=24))
  
  # Normally, scale_color_manual(values = mycolors)
  # But I edit for Minkyoung's request
  
  pca_plot2 <- ggplot(data = chosen_pca, aes(x=PC1, y=PC3, 
                                             colour=group,
                                             label = labels_text)) +
    geom_point(size=4) + geom_text_repel(size =7, show.legend = F,
                                         vjust=-0.35, hjust=-0.7) + 
    theme_light() + 
    #scale_color_manual(values = c(
    # "#df6020", "#58a77c" , "#c39b3c" , "#ff0040","#6d32cd"
    #),
    scale_color_manual(values = mycolors,
                       name="Conditions") +
    theme(axis.title = element_text(size=24),
          axis.text = element_text(size=20),
          legend.text = element_text(size=20),
          legend.title = element_text(size=24))
  
  pca_plot3 <- ggplot(data = chosen_pca, aes(x=PC2, y=PC3, 
                                             colour=group,
                                             label = labels_text)) +
    geom_point(size=4) + geom_text_repel(size =7, show.legend = F,
                                         vjust=-0.35, hjust=-0.7) + 
    theme_light() + scale_color_manual(values = mycolors,
                                       name="Conditions")+
    theme(axis.title = element_text(size=24),
          axis.text = element_text(size=20),
          legend.text = element_text(size=20),
          legend.title = element_text(size=24))
  
  pca_plot <- pca_plot1 + pca_plot2 + pca_plot3
  
  return(list(chosen_pca,expl_variance,pca_plot,
              pca_plot1,pca_plot2,pca_plot3))
}




# This is a function to run differential gene expression on edgeR
do_DE_edgeR <- function(data, grep_pattern, base_level="GFP",
                        controls_in=F, read_thr=60){
  library(edgeR)
  library(stringr)
  clean_data <- data[,grep(grep_pattern, 
                                 colnames(data))] 
  if(controls_in){
  clean_data <- clean_data[,-grep("ctrl", colnames(clean_data))] 
  }
  conditions <- as.factor(str_split(colnames(clean_data),"_",simplify = T)[,1])
  conditions <- relevel(conditions, base_level)
  clean_dge <- DGEList(counts = clean_data, group = conditions)
  # I set the min count to 60 here based on how many samples we are interested
  # in. If we get dpp4, intb4, myonb and gfp, we end up with 58 samples. So,
  # 1 read per gene across all sample should be around 60.
  clean_dge <- clean_dge[filterByExpr(
    clean_dge,min.total.count=read_thr), 
    ,
    keep.lib.sizes = F]
  clean_dge <- calcNormFactors(clean_dge)

  design <- model.matrix(~conditions)
  clean_dge <- estimateDisp(clean_dge, design)
  result_deg <- exactTest(clean_dge)
  result_table <- result_deg$table
  result_table$fdr_pvalue <- p.adjust(result_table$PValue, method = "fdr")
  return(list(result_table,clean_dge))
}


# This is a function to run differential gene expression on edgeR
do_DE_edgeR_control_vs_treat <- function(data, grep_pattern,
                                         read_thr=60){
  library(edgeR)
  library(stringr)
  clean_data <- data[,grep(grep_pattern, 
                           colnames(data))] 
  conditions <- rep("treated", dim(clean_data)[2])
  conditions[grep("c$",colnames(clean_data))] <- "control"
  conditions <- relevel(as.factor(conditions), "control")
  clean_dge <- DGEList(counts = clean_data, group = conditions)
  # I set the min count to 60 here based on how many samples we are interested
  # in. If we get dpp4, intb4, myonb and gfp, we end up with 58 samples. So,
  # 1 read per gene across all sample should be around 60.
  clean_dge <- clean_dge[filterByExpr(
    clean_dge,min.count=read_thr), 
    ,
    keep.lib.sizes = F]
  clean_dge <- calcNormFactors(clean_dge)
  
  design <- model.matrix(~conditions)
  clean_dge <- estimateDisp(clean_dge, design)
  result_deg <- exactTest(clean_dge)
  result_table <- result_deg$table
  result_table$fdr_pvalue <- p.adjust(result_table$PValue, method = "fdr")
  return(result_table)
}

# This is a simple function to get rid of the clutter caused by the
# custom dot plot. It is pretty hard coded, so you are warned.
do_dot_plot <- function(data, title, top){
  myPalette <- colorRampPalette(brewer.pal(2, "Dark2"))
  if(!missing(top)){
  data <- data[seq(1,top),]
  }
  if(nrow(data)>20){
    plt <- data[1:20,]
    gg <- ggplot(plt) +
      geom_bar(stat="identity",aes(y=reorder(Description, -pvalue), 
                                   x=pvalue, fill=NES)) +
      theme_light() + 
      theme(axis.title.y = element_blank(),
            axis.text = element_text(size = 16),
            axis.title = element_text(size = 18),
            plot.title = element_text(size=24),
            legend.text = element_text(size=16),
            legend.title = element_text(size=18)) +
      xlab("P value") +
      ggtitle(title) +
      scale_fill_gradientn(colours = myPalette(100))
  } else {
    gg <- ggplot(data) +
      geom_bar(stat="identity",aes(y=reorder(Description, -pvalue), 
                                   x=pvalue, fill=NES)) +
      theme_light() + 
      theme(axis.title.y = element_blank(),
            axis.text = element_text(size = 16),
            axis.title = element_text(size = 18),
            plot.title = element_text(size=24),
            legend.text = element_text(size=16),
            legend.title = element_text(size=18)) +
      xlab("P value") +
      ggtitle(title) +
      scale_fill_gradientn(colours = myPalette(100))
  }
  return(gg)
}

# And this is for plotting two of them. It returns two ggplot
# variables in a list, so these should be used with patchwork
do_dot_plot_two_conditions <- function(data1,data2,title1,title2){
  if(nrow(data1)>50){
    g1 <- ggplot(data1[1:50,]) +
      geom_point(aes(y=reorder(Description, NES), 
                     x=NES, size=setSize), color=pvalue) +
      theme_light() + 
      theme(axis.title.y = element_blank()) +
      xlab("Normalized Enrichment Score") +
      ggtitle(title1)
  } else {
    g1 <- ggplot(data1) +
      geom_point(aes(y=reorder(Description, NES), 
                     x=NES, size=setSize), color=pvalue) +
      theme_light() + 
      theme(axis.title.y = element_blank()) +
      xlab("Normalized Enrichment Score") +
      ggtitle(title1)
    
  }
  if(nrow(data2) >50 ){
    g2 <- ggplot(data2[1:50,]) +
      geom_point(aes(y=reorder(Description, NES), 
                     x=NES, size=setSize), color=qvalues) +
      theme_light() + 
      theme(axis.title.y = element_blank()) +
      xlab("Normalized Enrichment Score") +
      ggtitle(title2)
    
  } else {
    
    g2 <- ggplot(data2) +
      geom_point(aes(y=reorder(Description, NES), 
                     x=NES, size=setSize), color=qvalues) +
      theme_light() + 
      theme(axis.title.y = element_blank()) +
      xlab("Normalized Enrichment Score") +
      ggtitle(title2)
    
  }
  return(list(g1,g2))
}



# This function is to create a volcano plot. The original one
# plots normally; from both sides. I will also put a one-sided version
# since sometimes we only care about one direction of the log-fold change.
do_volcano_plot<- function(data, up_label="UP", down_label="DOWN", 
                           pval_thr=0.05, logFC_thr=0.5, 
                           colors=c("blue", "red", "black")){
  
  data$diffexpressed <- "NO"
  # if log2Foldchange > 0.5 and pvalue < 0.05, set as "UP" 
  data$diffexpressed[data$logFC > logFC_thr & 
                       data$PValue < pval_thr] <- up_label
  # if log2Foldchange < -0.5 and pvalue < 0.05, set as "DOWN"
  data$diffexpressed[data$logFC < -logFC_thr &
                       data$PValue < pval_thr] <- down_label
  
  # Add the gene names
  data$delabel <- NA
  data$delabel[data$diffexpressed != "NO"] <- rownames(data)[
    data$diffexpressed != "NO"]
  
  # Now plot them
  gg <- ggplot(data=data, aes(x=logFC, y=-log10(PValue), 
                              col=diffexpressed, label=delabel)) +
    geom_point() + 
    theme_minimal() +
    geom_text_repel() +
    scale_color_manual(values=colors) +
    # scale_color_manual(values=c("red", "black")) +
    geom_vline(xintercept=c(-logFC_thr, logFC_thr), col="red") +
    geom_hline(yintercept=-log10(pval_thr), col="red")
    # ggtitle("INTB4 vs GFP") +
    # theme(legend.position = "none")

  return(gg)
}

do_volcano_plot_one_side<- function(data, up_label="UP", 
                           pval_thr=0.05, logFC_thr=0.5, 
                           colors=c("red", "black")){
  
  data$diffexpressed <- "NO"
  # if log2Foldchange > 0.5 and pvalue < 0.05, set as "UP" 
  data$diffexpressed[data$logFC > logFC_thr & 
                       data$PValue < pval_thr] <- up_label
  
  # Add the gene names
  data$delabel <- NA
  data$delabel[data$diffexpressed != "NO"] <- rownames(data)[
    data$diffexpressed != "NO"]
  
  # Now plot them
  gg <- ggplot(data=data, aes(x=logFC, y=-log10(PValue), 
                        col=diffexpressed, label=delabel)) +
    geom_point() + 
    theme_minimal() +
    geom_text_repel() +
    scale_color_manual(values=colors) +
    # scale_color_manual(values=c("red", "black")) +
    geom_vline(xintercept=c(-logFC_thr, logFC_thr), col="red") +
    geom_hline(yintercept=-log10(pval_thr), col="red")
  # ggtitle("INTB4 vs GFP") +
  # theme(legend.position = "none")
  
  return(gg)
}


# This is a function to make volcano plots
do_volcano <- function(df, right, left, lfc_col, p_col, 
                       logFC=2, threshold = 0.05, maxoverlaptext=10){
  df$diffexpressed <- "NO"
  # if log2Foldchange > 0.5 and pvalue < 0.05, set as "UP" 
  df$diffexpressed[df[,lfc_col] > 2 & 
                     df[,p_col] < 0.05] <- right
  # if log2Foldchange < -0.5 and pvalue < 0.05, set as "DOWN"
  df$diffexpressed[df[,lfc_col] < -2 &
                     df[,p_col] < 0.05] <- left
  
  # Add the gene names
  df$delabel <- NA
  df$delabel[df$diffexpressed != "NO"] <- rownames(df)[
    df$diffexpressed != "NO"]
  
  p <- df %>%
    ggplot(aes(x=df[,lfc_col], y=-log10(df[,p_col]), 
               col=diffexpressed,
               label=delabel
    )) +
    geom_point() + 
    geom_text_repel(max.overlaps = maxoverlaptext) +
    theme_classic() +
    scale_color_manual(values=c("red", "darkblue", "black"),
                       name="Differential in",
                       breaks=c(left, right, "NO"),
                       labels=c(left, right, "No significant difference")) +
    geom_vline(xintercept=c(-2, 2), col="red") +
    geom_hline(yintercept=-log10(0.05), col="red") +
    ggtitle(paste0(left, " vs ", right)) +
    scale_y_continuous(name="-log10(pValue)") +
    scale_x_continuous(name= "Log2FoldChange") +
    theme(
      axis.text = element_text(size = 16),
      axis.title = element_text(size = 19),
      text = element_text(size = 24)
    )
  return(p)
}
















