#this script subsets the DEG list to IEGs

library(tidyverse)
library(schoolmath) #for is.positive



## Read in DE results  
# GABA:
genes <- read.csv("~/Downloads/RNAseq/AIBSmapping/OA/DESeq2/GABA_tr_vs_ctrl/train_vs_control_sig_genes.csv") 

## make column to indicate whether gene is upregulated or downregulated
genes <- genes[!is.na(genes$log2FoldChange), ]
#first remove NA values from log2FoldChange column
genes$boolean <- is.positive(genes$log2FoldChange)
genes$boolean <- as.character(genes$boolean)
#need to convert from logical to character type for recode to work
genes <- genes %>% mutate(result = recode(boolean,
                                            'TRUE'='upregulated',
                                            'FALSE'='downregulated'))
#remove boolean column now that it's done its purpose
genes$boolean <- NULL

## make IEG list
genelist <- c('Arc','Fos','Fosb','Fos-l1','Fos-l2','Junb','Jund','Jun','Npas4','Nr4a3','Nr4a2','Egr2','Nptx2','Cebpb','Scg2','Nefm','Vgf','Syt4')

## subset DE results to IEGs
genes1 <- subset(genes, gene %in% genelist)
write.csv(IEGs, "~/Downloads/RNAseq/AIBSmapping/OA/DESeq2/GABA_tr_vs_ctrl/GABA_IEG_DEresults.csv", row.names = FALSE)



## Read in DE results  
# Glut:
genes <- read.csv("~/Downloads/RNAseq/AIBSmapping/OA/DESeq2/glut_tr_vs_ctrl/train_vs_control_sig_genes.csv") 

## make column to indicate whether gene is upregulated or downregulated
genes <- genes[!is.na(genes$log2FoldChange), ]
#first remove NA values from log2FoldChange column
genes$boolean <- is.positive(genes$log2FoldChange)
genes$boolean <- as.character(genes$boolean)
#need to convert from logical to character type for recode to work
genes <- genes %>% mutate(result = recode(boolean,
                                          'TRUE'='upregulated',
                                          'FALSE'='downregulated'))
#remove boolean column now that it's done its purpose
genes$boolean <- NULL

## make IEG list
genelist <- c('Arc','Fos','Fosb','Fos-l1','Fos-l2','Junb','Jund','Jun','Npas4','Nr4a3','Nr4a2','Egr2','Nptx2','Cebpb','Scg2','Nefm','Vgf','Syt4')

## subset DE results to IEGs
genes1 <- subset(genes, gene %in% genelist)
write.csv(IEGs, "~/Downloads/RNAseq/AIBSmapping/OA/DESeq2/glut_tr_vs_ctrl/glut_IEG_DEresults.csv", row.names = FALSE)



## Read in DE results  
# L2/3:
genes <- read.csv("~/Downloads/RNAseq/AIBSmapping/OA/DESeq2/L23_0.3_tr_vs_ctrl/train_vs_control_sig_genes.csv") 

## make column to indicate whether gene is upregulated or downregulated
genes <- genes[!is.na(genes$log2FoldChange), ]
#first remove NA values from log2FoldChange column
genes$boolean <- is.positive(genes$log2FoldChange)
genes$boolean <- as.character(genes$boolean)
#need to convert from logical to character type for recode to work
genes <- genes %>% mutate(result = recode(boolean,
                                          'TRUE'='upregulated',
                                          'FALSE'='downregulated'))
#remove boolean column now that it's done its purpose
genes$boolean <- NULL

## make IEG list
genelist <- c('Arc','Fos','Fosb','Fos-l1','Fos-l2','Junb','Jund','Jun','Npas4','Nr4a3','Nr4a2','Egr2','Nptx2','Cebpb','Scg2','Nefm','Vgf','Syt4')

## subset DE results to IEGs
genes1 <- subset(genes, gene %in% genelist)
write.csv(IEGs, "~/Downloads/RNAseq/AIBSmapping/OA/DESeq2/L23_0.3_tr_vs_ctrl/L23_IEG_DEresults.csv", row.names = FALSE)


