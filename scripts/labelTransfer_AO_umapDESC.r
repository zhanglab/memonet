#this script is for performing label transfer of AIBS L2/3 cells onto our L2/3 classifier clusters, using the umap space from DESC

library(Seurat)
#library(scCATCH)
library(ggplot2)
#library(patchwork)
#library(cowplot)
library(tidyverse)
library(Matrix)
#library(RColorBrewer) #for brewer.pal()
print(sessionInfo)
sessionInfo()

options(warn=1)
  #this will print error messages to slurm file
options(future.globals.maxSize = 55000 * 1024^2) # This is important because preparing SCT normalzied data for integration needs a large amount of memory
#the integration step will give an error such as "total size of the 6 globals that need to be exported for the future expression (‘FUN()’) is 49.31 GiB. This exceeds the maximum allowed size of 48.83 GiB" so just change the number being multiplied by 1024 to encompass the total size of the globals


############# read in datasets ##############
#####  Read in our data L2/3 cells- REFERENCE #####
mat <- read.csv("/work/pi_yingzhang_uri_edu/kdunton/RNAseq/AIBSmapping/OA/count_matrices/unnormalized_counts_L23_0.3.csv", check.names=FALSE)
rownames(mat) <- mat[,1]
mat[,1] <- NULL
mtx <- as.matrix(mat)
dim(mtx)


## prepare metadata- read in the cluster file
path <- "/work/pi_yingzhang_uri_edu/kdunton/RNAseq/cluster_by_genes/0.3cutoff/DESC/clusters_n12.L0.85.csv"
print(paste0("input file:", path))
#print the path to the input file so you know where the analysis data came from
metadata <- read.csv(path)

colnames(metadata)[2] <- 'cluster_ID'
metadata$cluster_ID <- as.character(metadata$cluster_ID)
# add column for sample ID
metadata$sample <- sub("[[:print:]]*-", "",metadata$X)
#this takes the barcodes column and substitutes all the characters before the '-' with a blank, leaving just the sample number (1-6)
# add column for condition
#first have to coerce sample type (integer) to character
metadata$sample <- as.character(metadata$sample)
metadata <- metadata %>%
  mutate(stim = recode(sample,
                       "1" = "control",
                       "2" = "train",
                       "3" = "control",
                       "4" = "train",
                       "5" = "control",
                       "6" = "train"))



#make seurat obj
our_data = CreateSeuratObject(counts = mtx, min.cells = 0, project = "our_data")
our_data$dataset <- 'our_data'

#add metadata to seurat obj
#move barcodes to rownames
rownames(metadata) <- metadata$X
metadata$X <- NULL
our_data <- AddMetaData(object = our_data, metadata = metadata)

print("fin1")
rm(mtx, mat, metadata)




##### read in the AIBS dataset - QUERY #####
coo_aibs <- read.csv("/work/pi_yingzhang_uri_edu/kdunton/RNAseq/data/AIBS_data/aibs_matrix.mtx", header=FALSE)
  #the matrix is in COO format
genes <- read.csv("/work/pi_yingzhang_uri_edu/kdunton/RNAseq/data/AIBS_data/aibs_genes.tsv", check.names=FALSE, row.names=NULL)
  #since there are duplicate gene names, use row.names=NULL
metadata_sn_10X_aibs <- read.csv("/work/pi_yingzhang_uri_edu/kdunton/RNAseq/data/AIBS_data/aibs_barcodes.tsv")


## convert COO to dgcmatrix
counts <- Matrix::sparseMatrix(i = coo_aibs$V2, j = coo_aibs$V1, x = coo_aibs$V3)
#add row and col names
rownames(counts) <- genes$row.names
#the cell barcodes overlap with our_data, so add suffixes to these barcodes in order to differentiate
metadata_sn_10X_aibs$X <- paste("aibs",metadata_sn_10X_aibs$X,sep='_')
colnames(counts) <- metadata_sn_10X_aibs$X
print("aibs counts:")
head(counts[,1:5], 10)


## subset to L2/3 cells
metadata_sn_10X_aibs <- subset(metadata_sn_10X_aibs, subclass_label %in% 'L2/3 IT')
  #subclass_label is broadly L2/3 IT it doesn't have the _1,_2,_3 -- that's in cluster_label
print("dim:")
dim(counts)
counts <- counts[,metadata_sn_10X_aibs$X]
print("dim:")
dim(counts)


## make into seurat obj
sn_10X_AIBS = CreateSeuratObject(counts = counts, min.cells = 0, project = "sn_10X_AIBS")
#add information to identify dataset of origin, when merging objects later
sn_10X_AIBS$dataset <- 'sn_10X_AIBS'


## add metadata to seurat obj
#subset metadata to cells, cell type, and cluster columns
metadata <- metadata_sn_10X_aibs[,c("X","cluster_label","subclass_label","class_label")]
#make barcodes the rownames, this will allow metadata columns to transfer properly to seurat meta.data slot
rownames(metadata) <- metadata$X
#once barcodes are transferred to rownames, delete the barcode (X) column
metadata$X <- NULL
print("aibs metadata:")
head(metadata,10)
#add metadata to seurat obj
sn_10X_AIBS <- AddMetaData(object = sn_10X_AIBS, metadata = metadata)

print("fin2")
#remove data you don't need to free memory
rm(counts)
rm(metadata_sn_10X_aibs)
rm(metadata)
rm(coo_aibs)
rm(genes)



########### Normalization ############
data.reference <- our_data
rm(our_data)
data.query <- sn_10X_AIBS
rm(sn_10X_AIBS)

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
#data.reference <- RunUMAP(data.reference, reduction = "pca", dims = 1:30, verbose = FALSE, return.model=TRUE)
  #this step is only needed for computing 2D umap coordinates and getting a umap model to use for projecting cells
  #I specify the umap model parameters DESC used, and upload the DESC umap coordinates, so RunUMAP isn't needed


### in replacement of RunUMAP:
## load DESC umap data
umap <- read.csv("/work/pi_yingzhang_uri_edu/kdunton/RNAseq/cluster_by_genes/0.3cutoff/DESC/umap_n12.L0.85.csv")
#load cluster data in order to get the barcodes for the umap coordinates
clusters <- read.csv("/work/pi_yingzhang_uri_edu/kdunton/RNAseq/cluster_by_genes/0.3cutoff/DESC/clusters_n12.L0.85.csv")
#the umap file is indexed based on the ordering of the clusters file, so just add the barcode column of the clusters file to the umap file
umap$barcode <- clusters$X
umap$X <- NULL
  #this is the indexing for cells, don't need
umap <- rename(umap, umap_1=X0, umap_2=X1)
#format for adding to seurat obj
rownames(umap) <- umap$barcode
umap$barcode <- NULL
#make into a matrix for seurat
umap <- data.matrix(umap)

## add to object
data.reference[["umap.desc"]] <- CreateDimReducObject(embeddings = umap, key="UMAPdesc_",  assay = DefaultAssay(data.reference))
  #key designates the names of the columns, so for example if you plot this, the x and y axis will be named UMAPdesc_1 and UMAPdesc_2 for x and y respectively
  #umap.desc is the name of the slot the umap coordinates are added to

## add umap model parameters (the parameters used by DESC)
print("setup umap model with desc umap parameters")
# set UMAP model
umap.desc.model <- list()
umap.desc.model$n_epochs <- 300
umap.desc.model$alpha <-1
umap.desc.model$method <- "umap"
umap.desc.model$negative_sample_rate <- 5
umap.desc.model$gamma <- 1
umap.desc.model$approx_pow <- 0
umap.desc.model$metric$euclidean <- list()
umap.desc.model$embedding <- data.reference[["umap.desc"]]@cell.embeddings
ab_param <- uwot:::find_ab_params(spread = 1, min_dist = 0.5)
umap.desc.model$a <- ab_param["a"]
umap.desc.model$b <- ab_param["b"]
# add umap model to the slot that has DESC umap coordinates
data.reference[["umap.desc"]]@misc$model <- umap.desc.model


## plot ref using DESC umap- should look the same as DESC plot
bitmap("umap_clusterRef.png", width = 13, height = 11, units = 'in', res = 300)
DimPlot(data.reference, reduction = "umap.desc", group.by = "cluster_ID", label = TRUE, label.color='black',repel = TRUE) #+ NoLegend()
dev.off()




############# Project reference data labels onto query (our_data) ############# 
# the purpose of this is to assign the reference cluster labels onto our cells
data.anchors <- FindTransferAnchors(reference = data.reference, query = data.query,
                                        dims = 1:30, reference.reduction = "pca")
predictions <- TransferData(anchorset = data.anchors, refdata = data.reference$cluster_ID,
                            dims = 1:30)
  #this function isn't used in the ref, but I'm keeping it bc the purpose is to document the labels so I can save to file
data.query <- AddMetaData(data.query, metadata = predictions)

#save seurat object and file of prediciton scores
#saveRDS(data.query, file = "OA-seuratObj.rds")
table <- data.frame(data.query@meta.data)
table <- subset(table, select=-c(orig.ident,nCount_RNA,nFeature_RNA,dataset))
write.csv(table, "prediction_scores.csv")




############# Project query onto reference umap ############# 
## project aibs data into classifier DESC umap space
data.query <- MapQuery(anchorset = data.anchors, reference = data.reference, query = data.query,
                           refdata = list(celltype = "cluster_ID"), reference.reduction = "pca", reduction.model = "umap.desc")
  #this function uses pca info, but also uses the umap model "umap.desc" so that our cells project to the DESC umap space

#saveRDS(data.query, file = "aibs_Obj.rds")
#saveRDS(data.reference, file = "our_Obj.rds")


## color umap by predicted label
# to color cells with prediction scores below 0.2 as gray:
metadata <- data.query@meta.data
#metadata$predicted.celltype[metadata$prediction.score.max <= 0.4] <- "NA"
#data.query <- AddMetaData(object = data.query, metadata = metadata)
#colors <- c('brown1','darkgoldenrod1','darkkhaki','chartreuse2','darkturquoise','deepskyblue','cornflowerblue','darkorchid','darkgray')
#colors <- c('brown1','darkgoldenrod1','darkkhaki','chartreuse2','darkturquoise','deepskyblue','cornflowerblue','darkgray')#,'darkorchid','darkgray')

bitmap("umap_aibsL23_colored_by_cluster.png", width = 13, height = 11, units = 'in', res = 300)
DimPlot(data.query, reduction = "ref.umap", group.by = "predicted.celltype", label = TRUE,
              label.size = 3, label.color='black', repel = TRUE) + ggtitle("Predicted cluster labels") #+ NoLegend(), cols=colors
dev.off()

## color umap by aibs L2/3 subclasses
  #NA values automatically are gray
bitmap("umap_aibsL23_colored_by_aibsL23subclass.png", width = 13, height = 11, units = 'in', res = 300)
DimPlot(data.query, reduction = "ref.umap", group.by = "cluster_label", label = TRUE, label.size = 3, label.color='black', repel = TRUE) + ggtitle("Predicted cluster labels, colored by AIBS L2/3 IT subclass")
dev.off()



print("done!")