# this script runs normalization on L2/3 unnormalized counts, in preparation of DESC clustering
# these normalized counts are also used for the linear classifier analysis to identify EDGs

# normalization steps: cell normalization, natural log transform
  #cell-level normalization: the UMI count for each gene in each cell is divided by the total number of UMIs in the cell, multiplied by 10,000


library(tidyverse)
library(matrixStats)
library(Seurat)



########### read in 10X data into a Seurat object, do QC filtering and normalization using Seurat functions
#### Part 1: generate unnormalized counts with proper QC filters
## read in 10X counts matrix
data <- Read10X(data.dir = "~/Downloads/RNAseq/data/memonet_data/combined_cellranger_no-normalization/outs/filtered_feature_bc_matrix")
  #this is a dgCMatrix

## make seurat object
data_obj <- CreateSeuratObject(counts = data, min.cells = 0, min.features = 0)
length(colnames(data_obj))


## preprocessing - starting from the whole dataset without QC having been performed
# first subset cells
data_obj[["percent.mt"]] <- PercentageFeatureSet(data_obj, pattern = "^mt-")
data_obj_cellQC <- subset(data_obj, subset = nFeature_RNA > 200 & nCount_RNA < 30000 & nCount_RNA > 800 & percent.mt < 1)

# next subset genes - use min.cells parameter which is used when reading data into a Seurat object
# extract counts from object to make a new object using min.cells parameter
data_cellQC <- data_obj_cellQC@assays$RNA@counts
  #this generates a dgCMatrix which Seurat can read into an object
# read back into obj: filter genes found in < 3 cells
data_obj_QC_fin <- CreateSeuratObject(counts = data_cellQC, min.cells = 3, min.features = 0)
rm(data_obj,data_obj_cellQC,data_cellQC)


# save the cells and genes (for whole dataset of all cell types) remaining after QC filtering 
QC_genes <- rownames(data_obj_QC_fin)
QC_cells <- colnames(data_obj_QC_fin)
write.csv(QC_cells, "~/Downloads/RNAseq/QC/cells_after_QC.csv", quote=FALSE, row.names=FALSE)
write.csv(QC_genes, "~/Downloads/RNAseq/QC/genes_after_QC.csv", quote=FALSE, row.names=FALSE)


# how many cells per mouse? 
QC_cells <- as.data.frame(QC_cells)
names(QC_cells)[1] <- 'barcode'
QC_cells$sample <- sub("[[:print:]]*-", "",QC_cells$barcode)
QC_cells$sample <- as.character(QC_cells$sample)
table <- as.data.frame(table(QC_cells$sample))
colnames(table) <- c('mouse','n_cells')
write.csv(table, "~/Downloads/RNAseq/QC/cells_per_mouse.csv",row.names=FALSE)


## subset to L2/3 cells
# read in barcode file
cutoff <- '0'  #choose 0.3 or 0 
barcodes <- read.csv(paste0("~/Downloads/RNAseq/AIBSmapping/OA/barcode_files/L23barcodes-fromAIBS_",cutoff,".csv"))
names(barcodes)[1] <- 'barcode'
# subset
L23_data_obj <- data_obj_QC_fin[,colnames(data_obj_QC_fin) %in% barcodes$barcode]
rm(data_obj_QC_fin)
unnormalized_counts <- as.data.frame(L23_data_obj@assays$RNA@counts)


## save unnormalized counts to use in later analyses
write.csv(unnormalized_counts, paste0("~/Downloads/RNAseq/AIBSmapping/OA/count_matrices/unnormalized_counts_L23_",cutoff,".csv"), quote=FALSE)


#### Part 2: normalize 
L23_data_obj_norm <- NormalizeData(L23_data_obj, normalization.method = "LogNormalize",
                                   scale.factor = 10000)
normalized_counts <- as.data.frame(L23_data_obj_norm@assays$RNA@data)


## save counts as df 
if(cutoff == 0){
  # if the cutoff is 0, the normalized counts are for use outside of DESC, so saving as 15 decimals is fine (the default with write.csv)
  print("saving 0 cutoff")
  write.csv(normalized_counts, paste0("~/Downloads/RNAseq/AIBSmapping/OA/count_matrices/DESCnormalized_counts_L23_",cutoff,".csv"), quote=FALSE)
} else {
  # if cutoff is 0.3, save the file with the full 22 decimals for DESC
  print("saving 0.3 cutoff")
  normalized_counts <- format(normalized_counts, digits = 22)
  write.table(normalized_counts, paste0("~/Downloads/RNAseq/AIBSmapping/OA/count_matrices/DESCnormalized_counts_L23_",cutoff,".csv"), quote=FALSE, sep=',')
}



## For the 0.3 cutoff only: save a file of barcodes and stimulation corresponding to the normalized 0.3 counts file for use in the linear classifier analysis
cells <- as.data.frame(colnames(normalized_counts))
names(cells)[1] <- 'barcode'
cells$sample <- sub("[[:print:]]*-", "",cells$barcode)
cells$sample <- as.character(cells$sample)
cells <- cells %>%
  mutate(stim = recode(sample,
                       "1" = "control",
                       "2" = "train",
                       "3" = "control",
                       "4" = "train",
                       "5" = "control",
                       "6" = "train"))
write.csv(cells, paste0("~/Downloads/RNAseq/AIBSmapping/OA/count_matrices/sampleIDs_",cutoff,".csv"), row.names = FALSE, quote=FALSE)
cells_t <- as.data.frame(t(cells))
write.csv(cells, paste0("~/Downloads/RNAseq/AIBSmapping/OA/count_matrices/sampleIDs_",cutoff,"_transpose.csv"), row.names = FALSE, quote=FALSE)

# make sure the order is identical
cells_barcode <- cells$barcode
mtx_baarcode <- colnames(normalized_counts)
identical(cells_barcode,mtx_baarcode)


