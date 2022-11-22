#this script is the precursor for running DESC clustering. It subsets the counts matrix to the 3000 EDGs identified by the linear classifier
# input: normalized counts file from DESCnormalization.r
  # cell-level normalization: the UMI count for each gene in each cell is divided by the total number of UMIs in the cell, multiplied by 10,000
  # natural log transform



library(tidyverse)
library(matrixStats)
library(Seurat)



########### read in genelist to subset dataset to ########### 

### read in ranked gene list from linear classifier ###
genelist <- read.csv("~/Downloads/RNAseq/cluster_by_genes/0.3cutoff/PredictionGenesDescending0.3.csv", header=FALSE) 
genelist <- na.omit(genelist)
# subset to top 3000
genelist <- genelist[1:3000,]



########### read in L2/3 normalized counts file
countsi <- read.csv("~/Downloads/RNAseq/AIBSmapping/OA/count_matrices/DESCnormalized_counts_L23_0.3.csv",check.names=FALSE)
  #use check.names=FALSE so the dashes in cell names don't change to periods
counts <- countsi


#subset counts to just the EDGs
counts <- counts[genelist,] 
counts <- as.data.frame(t(counts))



########## save counts as a df to load into DESC
#DESC expects cells in rows and genes in columns

# save 22 digits
counts <- format(counts, digits = 22)
write.table(counts, "~/Downloads/RNAseq/cluster_by_genes/0.3cutoff/L23_0.3_EDGmtx.csv", quote=FALSE, sep=',')

