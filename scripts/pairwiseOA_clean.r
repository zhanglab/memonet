# this script is for mapping our dataset and transfering cell type labels: our data is query, reference is 10X AIBS dataset
# use this for chapter 1. Instead of using the barcode list from DESC QC filtering, do the filtering here.
# 'OA' denotes the mapping of Our dataset to AIBS dataset

library(Seurat)
library(scCATCH)
library(ggplot2)
library(patchwork)
library(cowplot)
library(tidyverse)
library(Matrix)
library(RColorBrewer) #for brewer.pal()

options(warn=1)
  #this will print error messages to slurm file
options(future.globals.maxSize = 55000 * 1024^2) # This is important because preparing SCT normalzied data for integration needs a large amount of memory
#the integration step will give an error such as "total size of the 6 globals that need to be exported for the future expression (‘FUN()’) is 49.31 GiB. This exceeds the maximum allowed size of 48.83 GiB" so just change the number being multiplied by 1024 to encompass the total size of the globals


############# read in datasets ##############
#####  Read in our data - QUERY #####
our_data.data = Read10X(data.dir = "/data/zhanglab/kdunton/6samples_cluster/deepseq_3_clustering/scrattch/combined_cellranger_no-normalization/outs/filtered_feature_bc_matrix/")
print("initial data:")
dim(our_data.data)

## make Seurat object
our_data = CreateSeuratObject(counts = our_data.data, min.cells = 0, project = "our_data")
  #don't filter genes here (with min.cells) bc you want to do this after filtering the cells
our_data$dataset <- 'our_data'
print("seurat obj initial:")
dim(our_data)

## calculate mitochondrial ratio
our_data[["percent.mt"]] <- PercentageFeatureSet(our_data, pattern = "^mt-")

## perform QC filtering
our_data <- subset(our_data, subset = nFeature_RNA > 200 & nCount_RNA > 800 & nCount_RNA < 30000 & percent.mt <1)
print("dims after cell filtering:")
dim(our_data)
#get list of the barcodes
barcodes <- our_data@meta.data
barcodes <- rownames(barcodes)

## subset initial data by the filtered cells
our_data.data <- our_data.data[,barcodes]
print("initial data, filtered cells:")
dim(our_data.data)

## read into seurat object again, this time with the filtered cells, and use min.cells=3 to filter genes. Want to filter genes after filtering cells, only because this is how pairwiseOA.r did it
our_data = CreateSeuratObject(counts = our_data.data, min.cells = 3, project = "our_data")
our_data$dataset <- 'our_data'
#add percent.mt column again 
our_data[["percent.mt"]] <- PercentageFeatureSet(our_data, pattern = "^mt-")
print("dims after gene filtering:")
dim(our_data)


print("fin1")
rm(our_data.data)




##### read in the 3 data files for 10X_sn_aibs dataset - REFERENCE #####
coo_aibs <- read.csv("/data/zhanglab/kdunton/6samples_cluster/deepseq_3_clustering/BICCN_integration/data_BICCN/10X-v3_sn_AIBS/aibs_matrix.mtx", header=FALSE)
  #the matrix is in COO format
genes <- read.csv("/data/zhanglab/kdunton/6samples_cluster/deepseq_3_clustering/BICCN_integration/data_BICCN/10X-v3_sn_AIBS/aibs_genes.tsv", check.names=FALSE, row.names=NULL)
  #since there are duplicate gene names, use row.names=NULL
metadata_sn_10X_aibs <- read.csv("/data/zhanglab/kdunton/6samples_cluster/deepseq_3_clustering/BICCN_integration/data_BICCN/10X-v3_sn_AIBS/aibs_barcodes.tsv")

## convert COO to dgcmatrix
counts <- Matrix::sparseMatrix(i = coo_aibs$V2, j = coo_aibs$V1, x = coo_aibs$V3)
#add row and col names
rownames(counts) <- genes$row.names
#the cell barcodes overlap with our_data, so add suffixes to these barcodes in order to differentiate
metadata_sn_10X_aibs$X <- paste("aibs",metadata_sn_10X_aibs$X,sep='_')
colnames(counts) <- metadata_sn_10X_aibs$X
print("aibs counts:")
head(counts[,1:5], 10)

## make into seurat obj
sn_10X_AIBS = CreateSeuratObject(counts = counts, min.cells = 3, project = "sn_10X_AIBS")
#add information to identify dataset of origin, when merging objects later
sn_10X_AIBS$dataset <- 'sn_10X_AIBS'

## add metadata to seurat obj
#subset metadata to cells, cell type, and cluster columns
metadata_sn_10X_aibs <- metadata_sn_10X_aibs[,c("X","cluster_label","subclass_label","class_label")]
#make barcodes the rownames, this will allow metadata columns to transfer properly to seurat meta.data slot
rownames(metadata_sn_10X_aibs) <- metadata_sn_10X_aibs$X
#once barcodes are transferred to rownames, delete the barcode (X) column
metadata_sn_10X_aibs$X <- NULL
print("aibs metadata:")
head(metadata_sn_10X_aibs,10)
#add metadata to seurat obj
sn_10X_AIBS <- AddMetaData(object = sn_10X_AIBS, metadata = metadata_sn_10X_aibs)

print("fin2")
#remove data you don't need to free memory
rm(counts)
rm(metadata_sn_10X_aibs)
rm(coo_aibs)
rm(genes)




########### Normalization ############
#specify which dataset in data.list is the reference
data.reference <- sn_10X_AIBS
rm(sn_10X_AIBS)
data.query <- our_data
rm(our_data)

data.query <- NormalizeData(data.query, verbose = FALSE)
data.query <- FindVariableFeatures(data.query, selection.method = "vst", nfeatures = 2000,
                                             verbose = FALSE)
data.reference <- NormalizeData(data.reference, verbose = FALSE)
data.reference <- FindVariableFeatures(data.reference, selection.method = "vst", nfeatures = 2000,
                                             verbose = FALSE)



############# scale, dimension reduction, visualization of reference ###########
# Run the standard workflow for visualization and clustering
data.reference <- ScaleData(data.reference, verbose = FALSE)
data.reference <- RunPCA(data.reference, npcs = 30, verbose = FALSE)
data.reference <- RunUMAP(data.reference, reduction = "pca", dims = 1:30, verbose = FALSE)

#original version with cluster_label labels on the umap space:
#bitmap("umap_referenceOA-celltypes.png", width = 14, height = 11, units = 'in', res = 300)
#DimPlot(data.reference, reduction = "umap", group.by = "cluster_label", label = TRUE, label.color='black',repel = TRUE, label.size = 5) + NoLegend() + ggtitle("AIBS dataset") + 
#  theme(axis.text = element_text(size = 15)) + theme(axis.title = element_text(size = 20)) + theme(plot.title = element_text(size = 20)) 
#ggsave(filename="umap_referenceOA-celltypes.svg",width = 14, height = 11)

# version with legend instead of labels on the umap space, and using subclass_label instead for broader labels:
bitmap("umap_AIBS_subclassLabelLegend.png", width = 14, height = 11, units = 'in', res = 300)
DimPlot(data.reference, reduction = "umap", group.by = "subclass_label") + ggtitle("AIBS dataset") +
  theme(axis.text = element_text(size = 15)) + theme(axis.title = element_text(size = 20)) + theme(plot.title = element_text(size = 20)) +
  guides(color = guide_legend(override.aes = list(size=4), ncol=2))
ggsave(filename="umap_AIBS_subclassLabelLegend.svg",width = 14, height = 11)
dev.off()



############# Assign reference data labels to query ############# 
# the purpose of this is to assign the reference cluster labels onto our cells
data.anchors <- FindTransferAnchors(reference = data.reference, query = data.query,
                                        dims = 1:30, reference.reduction = "pca")
predictions <- TransferData(anchorset = data.anchors, refdata = data.reference$cluster_label,
                            dims = 1:30)
data.query <- AddMetaData(data.query, metadata = predictions)
  #this step adds a column of predicted.id that comes from 'predictions' variable

# save seurat object and file of prediciton scores
#saveRDS(data.query, file = "OA-seuratObj.rds")
table <- data.frame(data.query@meta.data)
print("predictions table:")
print(table[1:5,1:10])
table <- subset(table, select=-c(orig.ident,nCount_RNA,nFeature_RNA,percent.mt,dataset))
table <- table %>% rename(
  AIBS_predicted_label=predicted.id
)
write.csv(table, "prediction_scores.csv")




############# Project query onto reference umap ############# 
data.reference <- RunUMAP(data.reference, dims = 1:30, reduction = "pca", return.model = TRUE)
#data.query <- MapQuery(anchorset = data.anchors, reference = data.reference, query = data.query,
 #                          refdata = list(celltype = "cluster_label"), reference.reduction = "pca", reduction.model = "umap")
data.query <- MapQuery(anchorset = data.anchors, reference = data.reference, query = data.query,
                           refdata = list(celltype = "subclass_label"), reference.reduction = "pca", reduction.model = "umap")

# original version with cluster_label labels on the umap space:
#bitmap("umap_OA_query-predictedLabels.png", width = 14, height = 11, units = 'in', res = 300)
#DimPlot(data.query, reduction = "ref.umap", group.by = "predicted.celltype", label = TRUE, label.size = 5, label.color='black', repel = TRUE) + 
 # ggtitle("Predicted AIBS labels") + labs(x = "UMAP_1",y = "UMAP_2") + 
  #NoLegend() + theme(axis.text = element_text(size = 15)) + theme(axis.title = element_text(size = 20)) + theme(plot.title = element_text(size = 20))
#ggsave(filename="umap_OA_query-predictedLabels.svg",width = 14, height = 11)

# version with legend instead of labels on the umap space, and using subclass_label for broader labels:
bitmap("umap_MEMONET_subclassLabelLegend.png", width = 14, height = 11, units = 'in', res = 300)
DimPlot(data.query, reduction = "ref.umap", group.by = "predicted.celltype") +
  ggtitle("MEMONET dataset") + labs(x = "UMAP_1",y = "UMAP_2") +
  theme(axis.text = element_text(size = 15)) + theme(axis.title = element_text(size = 20)) + theme(plot.title = element_text(size = 20)) +
  guides(color = guide_legend(override.aes = list(size=4), ncol=2))
ggsave(filename="umap_MEMONET_subclassLabelLegend.svg",width = 14, height = 11)
dev.off()


print("done!")

