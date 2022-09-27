# this script runs DESeq2 on DESC subclustering for one cluster vs another (using control and train cells0
# cluster X vs cluster Y- specify in code which cluster is the reference and it'll loop through and compare each of the other clusters to the ref
# for this, you want the input to be the counts matrix with genes in rows and cells in columns


library(tidyverse)
library(edgeR)
library(dplyr)
library(magrittr)
library(Matrix)
library(purrr)
library(reshape2)
library(S4Vectors)
library(tibble)
library(apeglm)
library(png)
library(DESeq2)

#for running in parallel (since we're using all the cells as samples rather than the aggregation of all cells in a mouse per cluster like the original DESeq2_DESC script does):
library("BiocParallel")
#register(MulticoreParam(10))


####### load MEX into R ########
mat <- read.csv("/work/pi_yingzhang_uri_edu/kdunton/RNAseq/AIBSmapping/OA/count_matrices/unnormalized_counts_L23_0.3.csv", check.names=FALSE)
rownames(mat) <- mat[,1]
mat[,1] <- NULL
mtx <- as.matrix(mat)
dim(mtx)




######## format metadata ########
#path <- "/work/pi_yingzhang_uri_edu/kdunton/RNAseq/cluster_by_genes/0.3cutoff/DESC/clusters_n12.L0.85.csv"
path <- "/work/pi_yingzhang_uri_edu/kdunton/RNAseq/cluster_by_genes/DESC_parameter_test/loop/clusters_dec22.n24.L0.85.csv"
# specify reference cluster
ref <- '2'

print(paste0("input file:", path))
#print the path to the input file so you know where the analysis data came from
clustersi <- read.csv(path)


## edit format of clusters
# rename column of barcodes from 'X' to 'barcodes'
colnames(clustersi)[1] <- 'barcodes'
colnames(clustersi)[2] <- 'cluster_label'
# add column for sample ID
clustersi$sample <- sub("[[:print:]]*-", "",clustersi$barcodes)
#this takes the barcodes column and substitutes all the characters before the '-' with a blank, leaving just the sample number (1-6)
# add column for condition
#first have to coerce sample type (integer) to character
clustersi$sample <- as.character(clustersi$sample)
clustersi <- clustersi %>%
  mutate(stim = recode(sample,
                       "1" = "control",
                       "2" = "train",
                       "3" = "control",
                       "4" = "train",
                       "5" = "control",
                       "6" = "train"))



#first create the overall directory for clusterX vs clusterY
dir_overall <- paste0("x_vs_", ref)
dir.create(dir_overall)
#then create a subdirectory for the all_cells scenario
dir <- paste0(dir_overall, "/all_cells")
#use dir variable when naming files
dir.create(dir)



## get list of cluster IDs
IDs <- levels(as_factor(clustersi$cluster_label))
# remove the ref cluster from the list
IDs <- IDs[!IDs %in% ref]


## loop through each cluster vs the ref cluster
print("starting loop")
for(subcluster in IDs){
  print(paste0("subcluster ", subcluster))
  
  # assign the clusters
  clusterX <- subcluster
  print(paste0("clusterX: ",clusterX))
  clusterY <- ref
  print(paste0("clusterY: ",clusterY))
  
  # subset to the pair of clusters
  clusters <- subset(clustersi, cluster_label %in% c(clusterX,clusterY))
  
  ## add comparison column to state clusterX and clusterY (could have done without this but wanted to keep format the same across scripts)
  metadata <- clusters %>% mutate(comparison = case_when(cluster_label == clusterX ~ "clusterX",
                                                         cluster_label == clusterY ~ "clusterY"))
  
  ######## subset gene count matrix #########
  metadata <- metadata %>% dplyr::mutate(barcodes = as.character(barcodes))
  #in order for subset to work correctly, need character type
  #need to ensure mutate function is used from dplyr package, bc there must be another package that has a mutate function too
  subset_matrix <- mtx[,metadata$barcodes]
  print("dim(subset_matrix:")
  print(dim(subset_matrix))
  print(subset_matrix[1:5,1:3])
  
  counts <- subset_matrix
  
  
  ########## make the barcodes column of metadata into rownames for correct input format to sce/dds object #########
  rownames(metadata) <- metadata[,1]
  #assign the first column (barcodes) to be rownames
  metadata[,1] <- NULL
  #then remove the barcodes column
  print(metadata[1:49,])
  #this shows that the rownames have indeed become the cell barcodes
  #printing 49 lines to show all the clusters to double check which cluster the script is evaluating
  
  #check that the cells are in the same order bt count matrix and metadata
  print(all(rownames(metadata) == colnames(counts)))
  #true means they are in the same order
  
  
  ######## run DESeq on all cells in clusterX to clusterY
  cluster_metadata <- metadata[,c("stim","comparison")]
  #only take the comparison column (barcodes are the rownames now), as that's all you need for this analysis
  #but only subsetting to 1 column makes it have no dimensions and that causes an error, so just add stim too
  
  dds <- DESeqDataSetFromMatrix(counts,
                                colData = cluster_metadata,
                                design = ~ comparison)
  
  #manually set sizeFactors 
  sizeFactors(dds)<- scran::calculateSumFactors(dds)
  normalization <- sizeFactors(dds)
  
  # Run DESeq2 differential expression analysis
  #first set contrast
  dds$comparison <- relevel(dds$comparison, ref="clusterY")
  dds <- DESeq(dds, parallel=TRUE, BPPARAM=MulticoreParam(8), useT=TRUE, minmu=1e-6, minReplicatesForReplace=Inf)
  
  print('check1')
  print("resultsNames(dds):")
  print(resultsNames(dds))
  #use this to figure out the name to use for coef for lfcShrink() 
  res <- results(dds,
                 alpha = 0.05)
  print('check2')
  print("resultsNames(dds):")
  print(resultsNames(dds))
  print(res)
  res <- lfcShrink(dds,
                   type="apeglm",
                   coef="comparison_clusterX_vs_clusterY")
  #coef necessary for apeglm method; can only use coef OR contrast but not both)
  print('check3')
  print("resultsNames(dds):")
  print(resultsNames(dds))
  print(res)
  
  # Set thresholds
  padj_cutoff <- 0.05
  
  # Turn the results object into a tibble for use with tidyverse functions
  res_tbl <- res %>%
    data.frame() %>%
    rownames_to_column(var="gene") %>%
    as_tibble()
  
  write.csv(res_tbl,
            paste0(dir, "/", clusterX, "_vs_", clusterY, "_all_genes.csv"),
            quote = FALSE,
            row.names = FALSE)
  
  # Subset the significant results
  sig_res <- dplyr::filter(res_tbl, padj < padj_cutoff) %>%
    dplyr::arrange(padj)
  
  write.csv(sig_res,
            paste0(dir, "/", clusterX, "_vs_", clusterY, "_sig_genes.csv"),
            quote = FALSE,
            row.names = FALSE)
}


warnings()
print("done!")


