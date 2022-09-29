# this script runs DESC-based normalization on L2/3 unnormalized counts
# the normalized counts will then be used to perform classifier analysis to identify DGs

#follow DESC normalization: cell normalization, natural log transform
  #cell-level normalization: the UMI count for each gene in each cell is divided by the total number of UMIs in the cell, multiplied by 10,000


library(tidyverse)
library(matrixStats)
library(Seurat)



########### Method 1: read in 10X data, do QC filtering in Seurat, then normalize in seurat
#### Part 1: generate unnormalized file with proper QC filters
## read in 10X matrix
data <- Read10X(data.dir = "~/Downloads/RNAseq/data/memonet_data/combined_cellranger_no-normalization/outs/filtered_feature_bc_matrix")
  #this is a dgCMatrix

## make seurat object
data_obj <- CreateSeuratObject(counts = data, min.cells = 0, min.features = 0)
length(colnames(data_obj))


## preprocessing - you're starting from the whole dataset without the QC steps that were performed
# since the original QC subset genes after cells, you have to subset cells here (even tho you will be subsetting to L2/3 cells)
data_obj[["percent.mt"]] <- PercentageFeatureSet(data_obj, pattern = "^mt-")
data_obj_cellQC <- subset(data_obj, subset = nFeature_RNA > 200 & nCount_RNA < 30000 & nCount_RNA > 800 & percent.mt < 1)
# remove counts from object to make a new object using min.cells parameter
data_cellQC <- data_obj_cellQC@assays$RNA@counts
  #this generates a dgCMatrix which seurat can read into an object
# read back into obj: filter genes found in < 3 cells
data_obj_QC_fin <- CreateSeuratObject(counts = data_cellQC, min.cells = 3, min.features = 0)
rm(data_obj,data_obj_cellQC,data_cellQC)

# save unnormalized counts for whole dataset after QC filtering
  # won't work on laptop, vector memory exhausted
  # instead save the colnames and rownames
QC_genes <- rownames(data_obj_QC_fin)
QC_cells <- colnames(data_obj_QC_fin)
write.csv(QC_cells, "~/Downloads/RNAseq/QC/cells_after_QC.csv", quote=FALSE, row.names=FALSE)
write.csv(QC_genes, "~/Downloads/RNAseq/QC/genes_after_QC.csv", quote=FALSE, row.names=FALSE)


## subset to L2/3 cells
# read in barcode file
barcodes <- read.csv("~/Downloads/RNAseq/AIBSmapping/OA/barcode_files/L23barcodes-fromAIBS_0.3.csv")
names(barcodes)[1] <- 'barcode'
# subset
L23_data_obj <- data_obj_QC_fin[,colnames(data_obj_QC_fin) %in% barcodes$barcode]
rm(data_obj_QC_fin)
unnormalized_counts <- as.data.frame(L23_data_obj@assays$RNA@counts)


## save unnormalized counts to use in other analyses
write.csv(unnormalized_counts, "~/Downloads/RNAseq/AIBSmapping/OA/count_matrices/unnormalized_counts_L23_0.3.csv", quote=FALSE)


#### Part 2: normalize 
L23_data_obj_norm <- NormalizeData(L23_data_obj, normalization.method = "LogNormalize",
                                   scale.factor = 10000)
normalized_counts <- as.data.frame(L23_data_obj_norm@assays$RNA@data)


## save counts as df to send to Nathan for classifier step
#DESC expects cells in rows and genes in columns
# save as full 22 digit file; using write.csv will round to 15
normalized_counts <- format(normalized_counts, digits = 22)
write.table(normalized_counts, "~/Downloads/RNAseq/AIBSmapping/OA/count_matrices/DESCnormalized_counts_L23_0.3.csv", quote=FALSE, sep=',')



## save file of barcodes and stim corresponding to the normalized counts file to reference during classifier training
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
write.csv(cells, "~/Downloads/RNAseq/AIBSmapping/OA/count_matrices/sampleIDs9.10.22.csv", row.names = FALSE, quote=FALSE)
cells_t <- as.data.frame(t(cells))
write.csv(cells, "~/Downloads/RNAseq/AIBSmapping/OA/count_matrices/sampleIDs9.10.22_transpose.csv", row.names = FALSE, quote=FALSE)

# make sure the order is identical
cells_barcode <- cells$barcode
mtx_barcode <- colnames(normalized_counts)
identical(cells_barcode,mtx_barcode)
