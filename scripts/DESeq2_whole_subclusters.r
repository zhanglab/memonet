# this script runs DESeq2 on DESC subclustering for one cluster, using control and train cells
# it loops through all clusters
# it finds DE genes between one subcluster and all other subclusters
# for this, you want the input to be the counts matrix with genes in rows and cells in columns


#library(scater) #don't think this is necessary bc not using DESeq's default QC
#library(Seurat)  #I don't think this is necessary since not using seurat obj
library(tidyverse)
#library(cowplot)
#library(Matrix.utils)
library(edgeR)
library(dplyr)
library(magrittr)
library(Matrix)
library(purrr)
library(reshape2)
library(S4Vectors)
library(tibble)
#library(SingleCellExperiment)
#library(pheatmap)
library(apeglm)
library(png)
library(DESeq2)
#library(RColorBrewer)

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
path <- "/work/pi_yingzhang_uri_edu/kdunton/RNAseq/cluster_by_genes/0.3cutoff/DESC/clusters_n24.L0.7.csv"
print(paste0("input file:", path))
#print the path to the input file so you know where the analysis data came from
clusters <- read.csv(path)
#read in file


## edit format of clusters
colnames(clusters)[1] <- 'barcodes'
colnames(clusters)[2] <- 'cluster_label'
# add column for sample ID
clusters$sample <- sub("[[:print:]]*-", "",clusters$barcodes)
#this takes the barcodes column and substitutes all the characters before the '-' with a blank, leaving just the sample number (1-6)
#code came from scrattch_step3.r
# add column for condition
#first have to coerce sample type (integer) to character
clusters$sample <- as.character(clusters$sample)
clusters <- clusters %>%
  mutate(stim = recode(sample,
                       "1" = "control",
                       "2" = "train",
                       "3" = "control",
                       "4" = "train",
                       "5" = "control",
                       "6" = "train"))

dir.create("all_cells")



### start loop to run DESeq2 on each cluster compared to the others
#get list of clusters
IDs <- levels(as_factor(clusters$cluster_label))

print("starting loop")
for(subcluster in IDs){
  print(paste0("subcluster ", subcluster))

  ## add column to differentiate the cluster of interest from all others, since you want to find DE between a whole cluster and all others
  metadata <- clusters %>% mutate(comparison = case_when(cluster_label == subcluster ~ "compare",
                                                         cluster_label != subcluster ~ "others"))


  # no need to subset counts matrix since it's the unnormalized L2/3 file 
  counts <- mtx



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




  ######## run DESeq on all cells in one cluster compared to the others
  cluster_metadata <- metadata[,c("stim","comparison")]
    #only take the comparison column (barcodes are the rownames now), as that's all you need for this analysis
    #but only subsetting to 1 column makes it have no dimensions and that causes an error, so just add stim too

  dds <- DESeqDataSetFromMatrix(counts,
                                colData = cluster_metadata,
                                design = ~ comparison)

  #manually set sizeFactors 
  sizeFactors(dds)<- scran::calculateSumFactors(dds)
  normalization <- sizeFactors(dds)

  #save sizefactors only for the first iteration, as the file is the same over each iteration since no cell subsetting is happening (it's all L2/3)
  print("saving sizeFactors")
  if(file.exists("all_cells/normalized_sizeFactors_calculateSumFactors.csv")){
    print("sizeFactor file exists already")
  }
  else {
    write.csv(normalization,"all_cells/normalized_sizeFactors_calculateSumFactors.csv")
  }  #if the file exists, skip; if it doesn't exist, write to file https://www.datamentor.io/r-programming/if-else-statement/



  # Run DESeq2 differential expression analysis
  #first set contrast
  dds$comparison <- relevel(dds$comparison, ref="others")
  dds <- DESeq(dds, parallel=TRUE, BPPARAM=MulticoreParam(8), useT=TRUE, minmu=1e-6, minReplicatesForReplace=Inf)


  # Output results of Wald test for contrast for A vs B
  print('check1')
  print("resultsNames(dds):")
  print(resultsNames(dds))
    #use this to figure out the name to use for coef for lfcShrink() 
  res <- results(dds,
                 #contrast = contrast,
                 alpha = 0.05)
  print('check2')
  print("resultsNames(dds):")
  print(resultsNames(dds))
  print(res)
  res <- lfcShrink(dds,
                   #contrast =  contrast,
                   #res=res)
                   type="apeglm",
                   coef="comparison_compare_vs_others")
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
            paste0("all_cells/", subcluster, "_vs_", "others", "_all_genes.csv"),
            quote = FALSE,
            row.names = FALSE)

  # Subset the significant results
  sig_res <- dplyr::filter(res_tbl, padj < padj_cutoff) %>%
    dplyr::arrange(padj)

  write.csv(sig_res,
            paste0("all_cells/", subcluster, "_vs_", "others", "_sig_genes.csv"),
            quote = FALSE,
            row.names = FALSE)


  print("saving counts matrices")
  #save normalized counts to file, same method as for size factors
  normalized_counts <- counts(dds, 
                              normalized = TRUE)
  if(file.exists("all_cells/normalized_counts_from_dds.csv")){
    print("normalized counts file exists already")
  }
  else {
    write.csv(normalized_counts,"all_cells/normalized_counts_from_dds.csv")
  }

  #save UNnormalized counts to file
  unnormalized_counts <- counts(dds,
                              normalized = FALSE)
  if(file.exists("all_cells/unnormalized_counts_from_dds.csv")){
    print("unnormalized counts file exists already")
  }
  else {
    write.csv(unnormalized_counts,"all_cells/unnormalized_counts_from_dds.csv")
  }
}

warnings()
print("done!")


