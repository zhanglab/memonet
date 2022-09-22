# this script runs DESeq2 comparing ctrl vs train
# use for L2/3, glut, GABA


#library(Seurat) 
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


#for running in parallel:
library("BiocParallel")
#register(MulticoreParam(6))



####### Read in 10X data
matrix_dir = "/work/pi_yingzhang_uri_edu/kdunton/RNAseq/data/memonet_data/combined_cellranger_no-normalization/outs/filtered_feature_bc_matrix/"
barcode.path <- paste0(matrix_dir, "barcodes.tsv.gz")
features.path <- paste0(matrix_dir, "features.tsv.gz")
matrix.path <- paste0(matrix_dir, "matrix.mtx.gz")
mat <- readMM(file = matrix.path)
feature.names = read.delim(features.path,
                           header = FALSE,
                           stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path,
                           header = FALSE,
                           stringsAsFactors = FALSE)
colnames(mat) = barcode.names$V1
rownames(mat) = feature.names$V2
  #use V2 to get common gene names; V1 is for ensembl IDs
mtx <- as.matrix(mat)
print("whole dataset dim")
print(dim(mtx))
rm(mat)


### subset to genes and cells remaining after QC
# read in cells 
QC_cells <- read.csv("/work/pi_yingzhang_uri_edu/kdunton/RNAseq/QC/cells_after_QC.csv")
# make into vector
QC_cells <- QC_cells[,1]
## subset cells
counts <- mtx[,QC_cells]

## remove genes expressed in <3 cells
# do it this way instead of as above for cells bc some genes are duplicates in mtx so it doesn't work
# Output a logical vector for every gene on whether the more than zero counts per cell
nonzero <- counts > 0
# Sums all TRUE values and returns TRUE if more than 10 TRUE values per gene
keep_genes <- Matrix::rowSums(nonzero) >= 3
# Only keeping those genes expressed in more than 10 cells
counts <- counts[keep_genes, ]

print("counts after QC subset")
print(dim(counts))
rm(mtx)


####### subset to cells of interest (L2/3, glut, gaba)
# read in barcode file
barcodes <- read.csv("/work/pi_yingzhang_uri_edu/kdunton/RNAseq/AIBSmapping/OA/barcode_files/AIBS-defined_glut_barcodes.csv")
names(barcodes)[1] <- 'barcode'
# subset
counts <- counts[,barcodes$barcode]
print("counts after barcode list subset")
print(dim(counts))


######## format metadata ########
# add column for sample ID
barcodes$sample <- sub("[[:print:]]*-", "",barcodes$barcode)
  #this takes the barcodes column and substitutes all the characters before the '-' with a blank, leaving just the sample number (1-6)
# add column for condition
  #first have to coerce sample type (integer) to character
barcodes$sample <- as.character(barcodes$sample)
barcodes <- barcodes %>%
  mutate(stim = recode(sample,
                            "1" = "control",
                            "2" = "train",
                            "3" = "control",
                            "4" = "train",
                            "5" = "control",
                            "6" = "train"))
metadata <- barcodes



########## make the barcodes column of metadata into rownames for correct input format to sce/dds object #########
metadata <- metadata %>% dplyr::mutate(barcode = as.character(barcode))
metadata[1:5,]
rownames(metadata) <- metadata[,'barcode']
  #assign the first column (barcodes) to be rownames
metadata$barcode <- NULL
  #then remove the barcodes column
dim(metadata)
metadata[1:5,]
  #this shows that the rownames have indeed become the cell barcodes
dim(counts)
counts[1:5,1:3]

#check that the cells are in the same order bt count matrix and metadata
print(all(rownames(metadata) == colnames(counts)))
  #true means they are in the same order



######## run DESeq on all cells train vs control
dds <- DESeqDataSetFromMatrix(counts,
                              colData = metadata,
                              design = ~ stim)

#manually set sizeFactors 
sizeFactors(dds) <- scran::calculateSumFactors(dds)
normalization <- sizeFactors(dds) 
write.csv(normalization, "normalized_sizeFactors_calculateSumFactors.csv")
  #this saves the size factors that you would divide the raw matrix by to get the normalized counts


# Run DESeq2 differential expression analysis
#first set contrast
dds$stim <- relevel(as_factor(dds$stim), ref="control")
dds <- DESeq(dds, parallel=TRUE, BPPARAM=MulticoreParam(6), useT=TRUE, minmu=1e-6, minReplicatesForReplace=Inf)


# Output results of Wald test 
print("resultsNames(dds):")
print(resultsNames(dds))
  #use this to figure out the name to use for coef for lfcShrink()  
res <- results(dds,
               alpha = 0.05)
res <- lfcShrink(dds,
                 type="apeglm",
                 coef="stim_train_vs_control")
                    #coef necessary for apeglm method; can only use coef OR contrast but not both

# Set thresholds
padj_cutoff <- 0.05

# Turn the results object into a tibble for use with tidyverse functions
res_tbl <- res %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>%
  as_tibble()

write.csv(res_tbl,
          "train_vs_control_all_genes.csv",
          quote = FALSE,
          row.names = FALSE)

# Subset the significant results
sig_res <- dplyr::filter(res_tbl, padj < padj_cutoff) %>%
  dplyr::arrange(padj)

write.csv(sig_res,
          "train_vs_control_sig_genes.csv",
          quote = FALSE,
          row.names = FALSE)


print("saving counts matrices")
normalized_counts <- counts(dds, 
                            normalized = TRUE)
write.csv(normalized_counts, "normalized_counts_from_dds.csv")

unnormalized_counts <- counts(dds,
                              normalized = FALSE)
write.csv(unnormalized_counts, "unnormalized_counts_from_dds.csv")




print("done!")
