#this script tests the accuracy of Seurat's mapping pipeline on the AIBS dataset
#run this in an array to lower runtime


library(Seurat)
#library(scCATCH)
library(ggplot2)
library(patchwork)
library(cowplot)
library(tidyverse)
library(Matrix)

args <- commandArgs(T)
jobid <- args[1]


options(warn=1)
  #this will print error messages to slurm file
options(future.globals.maxSize = 55000 * 1024^2) 


############# read in datasets ##############
##### read in the 3 data files for 10X_sn_aibs dataset #####
coo_aibs <- read.csv("~/Downloads/RNAseq/data/AIBS_data/aibs_matrix.mtx", header=FALSE)
#the matrix is in COO format
genes <- read.csv("~/Downloads/RNAseq/data/AIBS_data/aibs_genes.tsv", check.names=FALSE, row.names=NULL)
#since there are duplicate gene names, use row.names=NULL
metadata_sn_10X_aibs <- read.csv("~/Downloads/RNAseq/data/AIBS_data/10X-v3_sn_AIBS/aibs_barcodes.tsv")

## convert COO to dgcmatrix
counts <- Matrix::sparseMatrix(i = coo_aibs$V2, j = coo_aibs$V1, x = coo_aibs$V3)
#add row and col names
rownames(counts) <- genes$row.names
#the cell barcodes overlap with memonet data, so add suffixes to these barcodes in order to differentiate
metadata_sn_10X_aibs$X <- paste("aibs",metadata_sn_10X_aibs$X,sep='_')
colnames(counts) <- metadata_sn_10X_aibs$X
print("aibs counts:")
dim(counts)
head(counts[,1:5], 10)


## generate test query dataset- randomly sample 25% of cells from each cluster_label
cluster_labels <- unique(metadata_sn_10X_aibs$cluster_label)
dflist <- list()
  #initiate empty list to append loop df to
for(cl in cluster_labels){
  sampled <- metadata_sn_10X_aibs[sample(which(metadata_sn_10X_aibs$cluster_label==cl),round(0.25*length(which(metadata_sn_10X_aibs$cluster_label==cl)))),]
  dflist[[cl]] <- sampled
}
sampled_metadata <- do.call(rbind, dflist)
  #combine all df together from the dflist
print("query metadata:")
dim(sampled_metadata)
#subset the counts matrix by the sampled metadata
query_counts <- counts[,sampled_metadata$X]
print("query counts:")
dim(query_counts)

## subset original metadata and counts so that the sampled cells are not present in reference
ref_metadata <- subset(metadata_sn_10X_aibs, !X %in% sampled_metadata$X)
print("ref metadata:")
dim(ref_metadata)
counts <- counts[,ref_metadata$X]
print("ref counts:")
dim(counts)

## make ref into seurat obj
sn_10X_aibs = CreateSeuratObject(counts = counts, min.cells = 3, project = "sn_10X_aibs")
sn_10X_aibs$dataset <- 'sn_10X_aibs'

## add metadata to ref obj
ref_metadata <- ref_metadata[,c("X","cluster_label","subclass_label","class_label")]
rownames(ref_metadata) <- ref_metadata$X
ref_metadata$X <- NULL
print("aibs metadata:")
head(ref_metadata,10)
sn_10X_aibs <- AddMetaData(object = sn_10X_aibs, metadata = ref_metadata)

## make query into seurat obj
query = CreateSeuratObject(counts = query_counts, min.cells = 3, project = "query")
query$dataset <- 'query'

## add metadata to query obj
sampled_metadata <- sampled_metadata[,c("X","cluster_label","subclass_label","class_label")]
rownames(sampled_metadata) <- sampled_metadata$X
sampled_metadata$X <- NULL
print("aibs metadata:")
head(sampled_metadata,10)
query <- AddMetaData(object = query, metadata = sampled_metadata)


print("fin1")
rm(counts, metadata_sn_10X_aibs, genes, sampled_metadata, query_counts, ref_metadata, coo_aibs)





########### Normalization ############
data.reference <- sn_10X_aibs
rm(sn_10X_aibs)
data.query <- query
rm(query)

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





############# Project reference data labels onto query (our_data) ############# 
# the purpose of this is to assign the reference cluster labels onto query cells
data.anchors <- FindTransferAnchors(reference = data.reference, query = data.query,
                                         dims = 1:30, reference.reduction = "pca")
predictions <- TransferData(anchorset = data.anchors, refdata = data.reference$cluster_label,
                             dims = 1:30)
data.query <- AddMetaData(data.query, metadata = predictions)
 

table <- data.frame(data.query@meta.data)
table <- subset(table, select=-c(orig.ident,nCount_RNA,nFeature_RNA,dataset))
table <- table %>% rename(
  original_label=cluster_label,
  AIBS_predicted_label=predicted.id
)
write.csv(table, paste0("prediction_scores_",jobid,".csv"))


print("done!")

