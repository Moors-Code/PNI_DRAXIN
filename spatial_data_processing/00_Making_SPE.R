# Building SPE

# ---- Parsing ----

library(optparse)

parser <- OptionParser()
parser <- add_option(parser, c("-i", "--input"), type="character", 
		     help="Path to the cellranger bam output")
parser <- add_option(parser, c("-o", "--output"), type="character", 
		     help="Path to the SPE")
opt <- parse_args(parser)

# ---- Define Input ----

sample_dirs <- unlist(strsplit(opt$input, split=",", fixed=TRUE))
sample_dirs <- gsub("possorted_genome_bam.bam","",sample_dirs)
samples <- unlist(lapply(strsplit(sample_dirs, split="/", fixed=TRUE), function(X) X[length(X)-1]))

### Important ####
# Snakemake wildcards are in random order, hence the samples/sample_dirs order might vary from run to run 
# To ensure downstream reproducibility (which can be affected by column order) we manually order them here based on the first run we had
ordr <- c("54703bl28_highPNI","54703bl38_noPNI","6249bl23_noPNI","6249bl18_noPNIandhighPNI","70865bl15_highPNI","70865bl5_noPNI")
sample_dirs <- sample_dirs[match(ordr,samples)]
samples <- samples[match(ordr,samples)]

# ---- Build SPE ----

library(SpatialExperiment)
# Have to adjust the tissue_positions_list to tissue_positions that's why we have the manual function otherwise could use read10xVisium
read10xFun <- function(direct,smp) {
# we read in the raw counts to contain all spots but need to subset to only filtered probes.
  fnm <- file.path(direct, "filtered_feature_bc_matrix/features.tsv.gz")
  feats <- read.table(fnm)[,1]

  fnm <- file.path(direct, "raw_feature_bc_matrix")
  sce <- DropletUtils::read10xCounts(fnm)
  sce <- sce[feats,]
  colnames(sce) <- sce$Barcode <-  paste0(smp,"_",sce$Barcode)
  # Subset to same features

# read in image data
       
  img <- readImgData(
      path = file.path(direct, "spatial"),
      imageSources = file.path(direct, "spatial", "tissue_hires_image.png"),
      sample_id = smp, load=TRUE)

# read in spatial coordinates
  fnm <- file.path(direct, "spatial", "tissue_positions.csv")
  xyz <- read.csv(fnm, header = TRUE)
  rownames(xyz) <- xyz$barcode <- paste0(smp,"_",xyz$barcode)
  xyz <- xyz[colnames(sce),]
#      col.names = c(
#          "barcode", "in_tissue", "array_row", "array_col",
#          "pxl_row_in_fullres", "pxl_col_in_fullres"))

# construct observation & feature metadata
  rd <- S4Vectors::DataFrame(
      symbol = rowData(sce)$Symbol)
  spe <- SpatialExperiment(
    assays = list(counts = assay(sce)),
    rowData = rd, 
    colData = DataFrame(xyz), 
    spatialCoordsNames = c("pxl_col_in_fullres", "pxl_row_in_fullres"),
    imgData = img,
    sample_id = smp)
}

spelist <- lapply(1:length(samples), function(IDX){
  tmp <- read10xFun(sample_dirs[IDX], samples[IDX])
  return(tmp)
})
spe <- Reduce(cbind,spelist)
saveRDS(spe,opt$output)
