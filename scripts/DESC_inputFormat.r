#this script is the precursor for input to DESC for L2/3 clustering based on highly ranked genes from classifier
# input: normalized counts file from DESCnormalization.r
  # cell-level normalization: the UMI count for each gene in each cell is divided by the total number of UMIs in the cell, multiplied by 10,000
  # natural log transform
#subset to DGs


library(tidyverse)
library(matrixStats)
library(Seurat)



########### read in genelist to subset dataset to ########### 

### read in DG list ###
genelist <- read.csv("~/Downloads/RNAseq/cluster_by_genes/0.3cutoff/PredictionGenesDescending0.3.csv", header=FALSE) 
genelist <- na.omit(genelist)
# we want to subset to top 3000
genelist <- genelist[1:3000,]



########### read in L2/3 normalized counts file
countsi <- read.csv("~/Downloads/RNAseq/AIBSmapping/OA/count_matrices/DESCnormalized_counts_L23_0.3.csv",check.names=FALSE)
  #use check.names=FALSE so the dashes in cell names don't change to periods
counts <- countsi

#subset counts to just the high ranked genes
counts <- counts[genelist,] 
counts <- as.data.frame(t(counts))



########## save counts 
#DESC expects cells in rows and genes in columns

# save 22 digits
counts <- format(counts, digits = 22)
write.table(counts, "~/Downloads/RNAseq/cluster_by_genes/0.3cutoff/DESC/L23_0.3_DGmtx.csv", quote=FALSE, sep=',')



