# this script denotes how many cells are removed at each filtering step for our dataset

library(Seurat)
#library(scCATCH)
library(ggplot2)
#library(patchwork)
#library(cowplot)
library(tidyverse)
library(Matrix)
#library(RColorBrewer) #for brewer.pal()

options(warn=1)
  #this will print error messages to slurm file
options(future.globals.maxSize = 55000 * 1024^2) # This is important because preparing SCT normalzied data for integration needs a large amount of memory
#the integration step will give an error such as "total size of the 6 globals that need to be exported for the future expression (‘FUN()’) is 49.31 GiB. This exceeds the maximum allowed size of 48.83 GiB" so just change the number being multiplied by 1024 to encompass the total size of the globals


#####  Read in MEMONET data - QUERY #####
our_data.data = Read10X(data.dir = "/Volumes/Katie_storage/memonet_backups/RNAseq_6.25.23/data/memonet_data/combined_cellranger_no-normalization/outs/filtered_feature_bc_matrix/")
print("initial data:")
dim(our_data.data)

## make Seurat object
our_data = CreateSeuratObject(counts = our_data.data, min.cells = 0, project = "our_data")
our_data$dataset <- 'our_data'
print("seurat obj initial:")
dim(our_data)

## calculate mitochondrial ratio
our_data[["percent.mt"]] <- PercentageFeatureSet(our_data, pattern = "^mt-")

## perform QC filtering
#our_data <- subset(our_data, subset = nFeature_RNA > 200 & nCount_RNA > 800 & nCount_RNA < 30000 & percent.mt <1)
print("original dims:")
dim(our_data)
our_data <- subset(our_data, subset = nFeature_RNA > 200)
print("dims after gene filtering:")
dim(our_data)
our_data <- subset(our_data, subset = nCount_RNA > 800)
print("dims after UMI lower cutoff filtering:")
dim(our_data)
our_data <- subset(our_data, subset = nCount_RNA < 30000)
print("dims after UMI upper cutoff filtering:")
dim(our_data)
our_data <- subset(our_data, subset = percent.mt <1)
print("dims after mt percentage filtering:")
dim(our_data)

#get list of the barcodes
barcodes <- our_data@meta.data
barcodes <- rownames(barcodes)

## subset initial data by the filtered cells
our_data.data <- our_data.data[,barcodes]
print("initial data, filtered cells:")
dim(our_data.data)

## read into seurat object again, this time with the filtered cells, and use min.cells=3 to filter genes
our_data = CreateSeuratObject(counts = our_data.data, min.cells = 3, project = "our_data")
our_data$dataset <- 'our_data'
#add percent.mt column again 
our_data[["percent.mt"]] <- PercentageFeatureSet(our_data, pattern = "^mt-")
print("dims after gene filtering:")
dim(our_data)


print("done!")


