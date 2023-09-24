# this script compares the IEGs found to be differentially expressed across neuron subtypes

library(tidyverse)
library(schoolmath) #for is.positive


### L5
genes <- read.csv('~/Downloads/RNAseq/AIBSmapping/OA/DESeq2/L5_tr_vs_ctrl/train_vs_control_sig_genes_L5.csv')

## make column for upregulated or downregulated
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
## summarize IEGs
genelist_IEG <- c('Arc','Fos','Fosb','Fos-l1','Fos-l2','Junb','Jund','Jun','Npas4','Nr4a3','Nr4a2','Egr2','Nptx2','Cebpb','Scg2','Nefm','Vgf','Syt4','Bdnf','Ntrk2','Baz1a')
genelist <- genelist_IEG
genes <- subset(genes, padj < 0.05)
genes.L5 <- subset(genes, gene %in% genelist)
genes.L5 <- genes.L5[,c('gene','result')]
genes.L5$celltype <- 'L5'

### L2/3
genesL23 <- read.csv('~/Downloads/RNAseq/AIBSmapping/OA/DESeq2/L23_0.3_tr_vs_ctrl/train_vs_control_sig_genes.csv')
genesL23 <- genesL23[!is.na(genesL23$log2FoldChange), ]
genesL23$boolean <- is.positive(genesL23$log2FoldChange)
genesL23$boolean <- as.character(genesL23$boolean)
genesL23 <- genesL23 %>% mutate(result = recode(boolean,
                                                'TRUE'='upregulated',
                                                'FALSE'='downregulated'))
genesL23$boolean <- NULL
genesL23 <- subset(genesL23, padj < 0.05)
genes.L23 <- subset(genesL23, gene %in% genelist)
genes.L23 <- genes.L23[,c('gene','result')]
genes.L23$celltype <- 'L23'

## pvalb
genesPvalb <- read.csv('~/Downloads/RNAseq/AIBSmapping/OA/DESeq2/pvalb_tr_vs_ctrl/train_vs_control_sig_genes_pvalb.csv')
genesPvalb <- genesPvalb[!is.na(genesPvalb$log2FoldChange), ]
genesPvalb$boolean <- is.positive(genesPvalb$log2FoldChange)
genesPvalb$boolean <- as.character(genesPvalb$boolean)
genesPvalb <- genesPvalb %>% mutate(result = recode(boolean,
                                                    'TRUE'='upregulated',
                                                    'FALSE'='downregulated'))
genesPvalb$boolean <- NULL
genesPvalb <- subset(genesPvalb, padj < 0.05)
genes.pvalb <- subset(genesPvalb, gene %in% genelist)
genes.pvalb <- genes.pvalb[,c('gene','result')]
genes.pvalb$celltype <- 'pvalb'

## sst
genesSst <- read.csv('~/Downloads/RNAseq/AIBSmapping/OA/DESeq2/sst_tr_vs_ctrl/train_vs_control_sig_genes_sst.csv')
genesSst <- genesSst[!is.na(genesSst$log2FoldChange), ]
genesSst$boolean <- is.positive(genesSst$log2FoldChange)
genesSst$boolean <- as.character(genesSst$boolean)
genesSst <- genesSst %>% mutate(result = recode(boolean,
                                                'TRUE'='upregulated',
                                                'FALSE'='downregulated'))
genesSst$boolean <- NULL
genesSst <- subset(genesSst, padj < 0.05)
genes.sst <- subset(genesSst, gene %in% genelist)
genes.sst <- genes.sst[,c('gene','result')]
genes.sst$celltype <- 'sst'

## vip
genesVip <- read.csv('~/Downloads/RNAseq/AIBSmapping/OA/DESeq2/vip_tr_vs_ctrl/train_vs_control_sig_genes_vip.csv')
genesVip <- genesVip[!is.na(genesVip$log2FoldChange), ]
genesVip$boolean <- is.positive(genesVip$log2FoldChange)
genesVip$boolean <- as.character(genesVip$boolean)
genesVip <- genesVip %>% mutate(result = recode(boolean,
                                                'TRUE'='upregulated',
                                                'FALSE'='downregulated'))
genesVip$boolean <- NULL
genesVip <- subset(genesVip, padj < 0.05)
genes.vip <- subset(genesVip, gene %in% genelist)
genes.vip <- genes.vip[,c('gene','result')]
genes.vip$celltype <- 'vip'

## combine
neuron_celltypes <- rbind(genes.L23, genes.vip, genes.sst, genes.pvalb, genes.L5)
neuron_celltypes2 <- as.data.frame(pivot_wider(neuron_celltypes, names_from = celltype, values_from = result))
rownames(neuron_celltypes2) <- neuron_celltypes2$gene
neuron_celltypes2$gene <- NULL
write.csv(neuron_celltypes2, "~/Downloads/RNAseq/AIBSmapping/OA/DESeq2/IEGcomparison_neuronTypes2.csv", quote=FALSE)
  #this file holds data for table3



