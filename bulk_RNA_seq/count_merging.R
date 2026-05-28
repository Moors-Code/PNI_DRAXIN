# This is the script to merge the data and get gene lengths for
# RPKM normalization

library(stringr)
library(openxlsx)
library(dplyr)
library(ggplot2)
library(edgeR)

# Get all the files
file_list <- 
  list.files(c("results/"),
             full.names = T)

file_list_dedup <- file_list[grep("dedup", file_list)]
file_list <- file_list[grep("dedup", file_list, invert = T)]

# Merge all the files together and put the sample name to columns
for (file in file_list){
  
  # if the merged dataset doesn't exist, create it
  if (!exists("dataset")){
    name <- gsub(x = basename(file), "_S[0-9]{1,2}\\S*", "")
    dataset <- read.table(file, header=FALSE, sep = "\t", row.names = 1)
    colnames(dataset) <- name
  } 
  else {
    name <- gsub(x = basename(file), "_S[0-9]{1,2}\\S*", "")
    temp_dataset <- read.table(file, header=FALSE, 
                               sep = "\t", row.names = 1)
    dataset<-cbind(dataset, temp_dataset)
    colnames(dataset) <- gsub("^V[0-9]*",name, colnames(dataset))
    rm(temp_dataset)
  }
}


# Merge all the files together and put the sample name to columns
for (file in file_list_dedup){
  
  # if the merged dataset doesn't exist, create it
  if (!exists("dataset_dedup")){
    name <- gsub(x = basename(file), "_S[0-9]{1,2}\\S*", "")
    dataset_dedup <- read.table(file, header=FALSE, sep = "\t", row.names = 1)
    colnames(dataset_dedup) <- name
  } 
  else {
    name <- gsub(x = basename(file), "_S[0-9]{1,2}\\S*", "")
    temp_dataset <- read.table(file, header=FALSE, 
                               sep = "\t", row.names = 1)
    dataset_dedup<-cbind(dataset_dedup, temp_dataset)
    colnames(dataset_dedup) <- gsub("^V[0-9]*",name, colnames(dataset_dedup))
    rm(temp_dataset)
  }
}


# I delete the rows with counts that show ambiguous or multimap
dataset <-dataset[-grep("__",rownames(dataset)),]
dataset_dedup <-dataset_dedup[-grep("__",rownames(dataset_dedup)),]

names(dataset) <- make.unique(names(dataset))
names(dataset_dedup) <- make.unique(names(dataset_dedup))

pdf("plots/read_counts.pdf", width = 16, height = 9)
# Convert the data to numeric
dataset <- mutate_all(dataset, function(x) as.numeric(x))
# I want to plot the counts so I make a read per million
pl_df <- as.data.frame(colSums(dataset)) / 1000000
colnames(pl_df) <- c("sum")
pl_df$samples <- rownames(pl_df)
pl_df$group <- str_split(row.names(pl_df), pattern = "_", simplify = T)[,1]

pl_df$date <- ifelse(
  grepl(pattern = "DRAXIN_rep[1-2]_1_2000", x = rownames(pl_df)),
  yes = paste(str_split(rownames(pl_df), "_", simplify = T)[,5],
              str_split(rownames(pl_df), "_", simplify = T)[,6],
              str_split(rownames(pl_df), "_", simplify = T)[,7]),
  no = paste(str_split(rownames(pl_df), "_", simplify = T)[,3],
             str_split(rownames(pl_df), "_", simplify = T)[,4],
             str_split(rownames(pl_df), "_", simplify = T)[,5]))

ggplot(data=pl_df, aes(x=reorder(samples, sum), y=sum, fill=date)) +
  geom_bar(stat = "identity", width = 0.8) +  # Adjust bar width
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 20),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 24),
    axis.text.y = element_text(size = 20),
    legend.position = "right",
    plot.margin = margin(20, 20, 20, 20)  # Add margins to the plot
  ) +
  ylab("Total mapped read counts (in millions)") +
  ggtitle("UMI extracted, no dedup")



# Convert the data to numeric
dataset_dedup <- mutate_all(dataset_dedup, function(x) as.numeric(x))
# I want to plot the counts so I make a read per million
pl_df <- as.data.frame(colSums(dataset_dedup)) / 1000000
colnames(pl_df) <- c("sum")
pl_df$samples <- rownames(pl_df)
pl_df$group <- str_split(row.names(pl_df), pattern = "_", simplify = T)[,1]

pl_df$date <- ifelse(
  grepl(pattern = "DRAXIN_rep[1-2]_1_2000", x = rownames(pl_df)),
  yes = paste(str_split(rownames(pl_df), "_", simplify = T)[,5],
              str_split(rownames(pl_df), "_", simplify = T)[,6],
              str_split(rownames(pl_df), "_", simplify = T)[,7]),
  no = paste(str_split(rownames(pl_df), "_", simplify = T)[,3],
             str_split(rownames(pl_df), "_", simplify = T)[,4],
             str_split(rownames(pl_df), "_", simplify = T)[,5]))

ggplot(data=pl_df, aes(x=reorder(samples, sum), y=sum, fill=date)) +
  geom_bar(stat = "identity", width = 0.8) +  # Adjust bar width
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 20),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 24),
    axis.text.y = element_text(size = 20),
    legend.position = "right",
    plot.margin = margin(20, 20, 20, 20)  # Add margins to the plot
  ) +
  ylab("Total mapped read counts (in millions)") + 
  ggtitle("UMI dedup")

dev.off()

# I want to check the correlation between all samples
cor_matrix <- cor(dataset)
library(pheatmap)
library(corrplot)

pdf ("plots/correlation.pdf", width = 9, height = 9)
corrplot(cor_matrix, method="color", type="upper", order="hclust", 
         tl.col="black", tl.srt=45, number=TRUE, number.cex=0.7, number.digits=2)

print(pheatmap(mat = cor_matrix, display_numbers = T, na_col="white", cluster_rows = T, cluster_cols = T))
dev.off()

cor_matrix_dedup <- cor(dataset_dedup)
pdf ("plots/correlation_dedup.pdf", width = 9, height = 9)
corrplot(cor_matrix_dedup, method="color", type="upper", order="hclust", 
         tl.col="black", tl.srt=45, number=TRUE, number.cex=0.7, number.digits=2)

print(pheatmap(mat = cor_matrix_dedup, display_numbers = T, na_col="white", cluster_rows = T, cluster_cols = T))

dev.off()

write.csv(dataset, "data/counts.csv")
write.csv(dataset_dedup, "data/counts_dedup.csv")
