### Visualize GO results: treeplots of individual clusters; heatmaps for specific GO terms
#https://yulab-smu.top/biomedical-knowledge-mining-book/enrichplot.html#enrichplot 



### install packages ###

#if (!require("BiocManager", quietly = TRUE))
 # install.packages("BiocManager")
#BiocManager::install("enrichplot")

#if (!require("BiocManager", quietly = TRUE))
 # install.packages("BiocManager")
#BiocManager::install("clusterProfiler")

#install DO.db, a dependency of DOSE and clusterProfiler
#if (!require("BiocManager", quietly = TRUE))
 # install.packages("BiocManager")
#BiocManager::install("DO.db")

#if (!require("BiocManager", quietly = TRUE))
 # install.packages("BiocManager")
#BiocManager::install("org.Mm.eg.db")



library(tidyverse)
library(schoolmath)
library(clusterProfiler) #runs GO
library(org.Mm.eg.db)
library(enrichplot) #visualize results


#--------------- Part A: GO analysis of EDGs ---------------#
genelist <- read.csv("~/Downloads/RNAseq/cluster_by_genes/0.3cutoff/PredictionGenesDescending0.3.csv", header=FALSE)
# subset to top 3000
genelist <- as.data.frame(genelist[1:3000,])
names(genelist)[1] <- 'gene'

#background list
background <- read.csv("~/Downloads/RNAseq/QC/genes_after_QC.csv")
names(background)[1] <- 'V1'

### convert symbol to entrez
DEG_entrez <- AnnotationDbi::select(org.Mm.eg.db, genelist$gene,
                                    "ENTREZID", "SYMBOL")
DEG_entrez <- DEG_entrez[,'ENTREZID']
DEG_entrez <- na.omit(DEG_entrez)

background_entrez <- AnnotationDbi::select(org.Mm.eg.db, background$V1,
                                           "ENTREZID", "SYMBOL")
background_entrez <- background_entrez[,'ENTREZID']
background_entrez <- na.omit(background_entrez)


### run GO
#https://yulab-smu.top/biomedical-knowledge-mining-book/clusterprofiler-go.html 
ego <- enrichGO(gene          = DEG_entrez,
                universe      = background_entrez,
                OrgDb         = org.Mm.eg.db,
                ont           = "BP",
                pAdjustMethod = "fdr",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.05,
                readable      = TRUE)
#readable will map entrezid to gene name
#gene has to be vector of entrez ids
#ont options: CC (cellular component), BP (bio process), MF (molecular function), ALL (for all 3)
# save results in table format
ego_save <- as.data.frame(ego)
write.csv(ego_save, "~/Downloads/RNAseq/cluster_by_genes/0.3cutoff/GO/EDG/EDG_GO.csv", quote=TRUE, row.names=FALSE)
#need quote=TRUE bc some values in the Description column have commas in them, and if they don't get quoted then R will interpret the comma to split into mulitple cells


### tree plot
ego2 <- pairwise_termsim(ego)
#https://rdrr.io/github/GuangchuangYu/enrichplot/man/treeplot.html
treeplot(ego2, offset=0.8, offset_tiplab=0.8, nWords=0, nCluster=5, showCategory=30, group_color=c()) 
  #offset is how close the category text is from the node text
  #offset_tiplab is how close the node text is from the end of the branch; use 0.3 for pos, 0.8 for neg
  #fontsize= changes just the category words not the GO terms
ggsave("~/Downloads/RNAseq/cluster_by_genes/0.3cutoff/GO/EDG/EDG_GO.png", width = 8, height=7.54)
ggsave("~/Downloads/RNAseq/cluster_by_genes/0.3cutoff/GO/EDG/EDG_GO.svg", width = 8, height=7.54)


#--------------- Part A-2: GO analysis of unique EDGs ---------------#
#run GO analysis on the unique EDGs (EDGs that don't show in the DE results of L2/3 train vs control)

# read in EDGs
EDGs <- read.csv("~/Downloads/RNAseq/cluster_by_genes/0.3cutoff/PredictionGenesDescending0.3.csv", header=FALSE) 
EDGs <- na.omit(EDGs)
# subset to top 3000
EDGs <- as.data.frame(EDGs[1:3000,])
colnames(EDGs)[1] <- 'gene'

# read in DEGs
DEG <- read.csv("~/Downloads/RNAseq/AIBSmapping/OA/DESeq2/L23_0.3_tr_vs_ctrl/train_vs_control_sig_genes.csv")

## get list of genes shared in EDG and DEG list - to have a column indicating with DEGs are also EDGs
DEG2 <- as.data.frame(DEG[,'gene'])
names(DEG2)[1] <- 'gene'
shared <- inner_join(EDGs, DEG2)
shared <- shared$gene
#make column in DEG list to indicate whether the gene is also in EDG list
DEG3 <- DEG
DEG3$Also_EDG <- ''
DEG3$Also_EDG <- ifelse(DEG3$gene %in% shared,"Yes", "No")
#https://stackoverflow.com/questions/16570302/how-to-add-a-factor-column-to-dataframe-based-on-a-conditional-statement-from-an
#make sure it worked and there are 1147 "yes"
table(DEG3$Also_EDG)
#save and add to the DEG excel sheet
write.csv(DEG3, "~/Downloads/RNAseq/AIBSmapping/OA/DESeq2/L23_0.3_tr_vs_ctrl/train_vs_control_sig_genes_EDGcolumn.csv", row.names = FALSE, quote=FALSE)

## get list of unique genes appearing in EDG but not DEG
unique_genes <- as.data.frame(setdiff(EDGs$gene,DEG$gene))
colnames(unique_genes)[1] <- 'gene'
#write.csv(unique_genes, "~/Downloads/RNAseq/cluster_by_genes/0.3cutoff/UniqueEDG_noDEG-allL23trvsctrl.csv", row.names=FALSE)

### run GO for the unique EDGs
genelist <- unique_genes
#genelist <- read.csv("~/Downloads/RNAseq/cluster_by_genes/0.3cutoff/UniqueEDG_noDEG-allL23trvsctrl.csv")

#background list
background <- read.csv("~/Downloads/RNAseq/QC/genes_after_QC.csv")
names(background)[1] <- 'V1'
# convert symbol to entrez
DEG_entrez <- AnnotationDbi::select(org.Mm.eg.db, genelist$gene,
                                    "ENTREZID", "SYMBOL")
DEG_entrez <- DEG_entrez[,'ENTREZID']
DEG_entrez <- na.omit(DEG_entrez)
background_entrez <- AnnotationDbi::select(org.Mm.eg.db, background$V1,
                                           "ENTREZID", "SYMBOL")
background_entrez <- background_entrez[,'ENTREZID']
background_entrez <- na.omit(background_entrez)
# run GO
#https://yulab-smu.top/biomedical-knowledge-mining-book/clusterprofiler-go.html 
ego2 <- enrichGO(gene          = DEG_entrez,
                 universe      = background_entrez,
                 OrgDb         = org.Mm.eg.db,
                 ont           = "BP",
                 pAdjustMethod = "fdr",
                 pvalueCutoff  = 0.05,
                 qvalueCutoff  = 0.05,
                 readable      = TRUE)
   #readable will map entrezid to gene name
   #gene has to be vector of entrez ids
   #ont options: CC (cellular component), BP (bio process), MF (molecular function), ALL (for all 3)
# save results in table format
ego_save2 <- as.data.frame(ego2@result) 
# the normal way of saving -- as.data.frame(ego2) -- doesn't work since nothing is significant
write.csv(ego_save2, "~/Downloads/RNAseq/cluster_by_genes/0.3cutoff/GO/EDG/uniqueEDG_GO.csv", quote=TRUE, row.names=FALSE)



#--------------- Part B. GO analysis of clusters: treeplot loop ---------------#
## saves a combined file of GO results for all clusters

### new analysis 0.3 cutoff
DEG <- read.csv("~/Downloads/RNAseq/cluster_by_genes/0.3cutoff/DESC/DESeq2/n25.L0.65/x_vs_0/all_cells/DEGstats_padj0.05.csv")

DEG <- subset(DEG, padj < 0.05)
regulation <- "up"  # choose 'up' or 'down'
DEG <- subset(DEG, result %in% paste0(regulation,"regulated"))

#background list
background <- read.csv("~/Downloads/RNAseq/QC/genes_after_QC.csv")
names(background)[1] <- 'V1'

### convert symbol to entrez
background_entrez <- AnnotationDbi::select(org.Mm.eg.db, background$V1,
                                           "ENTREZID", "SYMBOL")
background_entrez <- background_entrez[,'ENTREZID']
background_entrez <- na.omit(background_entrez)

# make empty df to append GO results to in loop
GOresults <- data.frame(matrix(0, ncol=9))
colnames(GOresults) <- c('cluster','Description','GeneRatio','BgRatio','pvalue','p.adjust','qvalue','geneID','Count')

### loop through clusters and run GO
clusters <- sort(unique(DEG$cluster))
for(cl in clusters){
  print(cl)
  DEG_cl <- subset(DEG, cluster %in% cl)

  ### convert symbol to entrez
  DEG_entrez <- AnnotationDbi::select(org.Mm.eg.db, DEG_cl$gene,
                                      "ENTREZID", "SYMBOL")
  DEG_entrez <- DEG_entrez[,'ENTREZID']
  DEG_entrez <- na.omit(DEG_entrez)


  ### run GO
  #https://yulab-smu.top/biomedical-knowledge-mining-book/clusterprofiler-go.html 
  ego <- clusterProfiler::enrichGO(gene          = DEG_entrez,
                  universe      = background_entrez,
                  OrgDb         = org.Mm.eg.db,
                  ont           = "BP",
                  pAdjustMethod = "fdr",
                  pvalueCutoff  = 0.05,
                  qvalueCutoff  = 0.05,
                  readable      = TRUE)
  #readable will map entrezid to gene name
  #gene has to be vector of entrez ids
  #ont options: CC (cellular component), BP (bio process), MF (molecular function), ALL (for all 3)

  ## save results
  ego_save <- as.data.frame(ego)
  ego_save$cluster <- cl
  ego_save <- ego_save[,c('cluster','Description','GeneRatio','BgRatio','pvalue','p.adjust','qvalue','geneID','Count')]
  GOresults <- rbind(GOresults,ego_save)

  ### tree plot
  ego2 <- enrichplot::pairwise_termsim(ego)
  p <- enrichplot::treeplot(ego2, offset=0.8, offset_tiplab=0.8, hexpand=0.2, nWords=0, nCluster=5, showCategory=30)
    #offset is how close the category text is from the node text
    #offset_tiplab is how close the node text is from the end of the branch
    #fontsize= changes just the category words not the GO terms
    #use plot zoom to screenshot
  ggsave(paste0("~/Downloads/RNAseq/cluster_by_genes/0.3cutoff/GO/clusters/n25.L0.65/x_vs_0/",cl,regulation,"_treeBio_n25.L0.65.png"), width = 10,height=7.54, plot=p)
  ggsave(file= paste0("~/Downloads/RNAseq/cluster_by_genes/0.3cutoff/GO/clusters/n25.L0.65/x_vs_0/",cl,regulation,"_treeBio_n25.L0.65.svg"),width = 10,height=7.54, plot= p)
  print("saved ggplot")
}
#remove first row of GOresults - it's just zeros from creation of the df
GOresults <- GOresults[-1,]
write.csv(GOresults, paste0("~/Downloads/RNAseq/cluster_by_genes/0.3cutoff/GO/clusters/n25.L0.65/x_vs_0/GOresults_n25.L0.65_Bio-",regulation,"regDEG.csv"), row.names = FALSE, quote=TRUE)
  #quote=true means any commas in the terms won't cause R to make a new cell
#------------------




#------------- Part C.1. Heatmaps and line plots of GO terms of interest -------------#
### compare terms of interest across clusters to show differences in learning-related terms
regulation <- 'upreg'
GO <- read.csv(paste0("~/Downloads/RNAseq/cluster_by_genes/0.3cutoff/GO/clusters/n25.L0.65/x_vs_0/GOresults_n25.L0.65_Bio-",regulation,"DEG.csv"), check.names=FALSE)
GOterms <- c('long-term synaptic potentiation','long-term synaptic depression','dendritic spine morphogenesis',
             'regulation of actin filament-based process','regulation of RNA splicing')
GO <- subset(GO, Description %in% GOterms)

# each cluster might have different genes contributing to a GO term, so combine genes across all clusters per term
for(term in GOterms){
  # split each term into a separate df
  x <- subset(GO, Description == term)
  # subset to the gene column
  x <- x[,'geneID']
  # replace / with ,
  x2 <- str_replace_all(x, c("/" = ","))
  # x2 is a list, where each sublist is the gene list from one cluster. Convert to string
  x2string <- toString(x2)
  # now all clusters are combined
  # remove duplicate genes that appear in multiple clusters
  x2string <- unique(unlist(strsplit(x2string, ",")))
  # becomes a character vector or something
  # put back to string
  x2string <- toString(x2string)
  assign(term,x2string)
}

# combine all terms into one df
df <- data.frame(matrix(0, nrow=length(GOterms), ncol=2))
colnames(df) <- c('term','genes')
df$term <- GOterms

for(term in GOterms){
  df$genes[df$term == term] <- get(term)
}

## expand df so each gene is in its own row but still labeled by GO term
df_separate <- separate_rows(df, genes, sep = ", ", convert = FALSE)
#need the space after the comma so the resulting cells don't leading white space
names(df_separate)[2] <- 'gene' 
# but some genes still have a leading white space:
genes_nowhite <- gsub("* ", "", df_separate$gene)
df_separate$gene <- genes_nowhite

## combine all GO term genes together and only keep unique ones - use to subset counts
# for generating mean gene expression per cluster, combine all strings together to get all genes together
#use `` for the df names bc there's spaces, and if you used '' that's like a string name so don't do that
allGenes <- paste(`long-term synaptic potentiation`,`long-term synaptic depression`,`dendritic spine morphogenesis`,
                  `regulation of actin filament-based process`,`regulation of RNA splicing`, sep=",", collapse=NULL)
split <- sort(unique(unlist(strsplit(allGenes, ","))))
# format genes
gene_df <- data.frame(matrix(0, nrow=length(split), ncol=1))
colnames(gene_df) <- c('gene')
i <- 1
for(gene in split){
  # remove leading white space that some genes have
  gene_new <- gsub(" ", "", gene)
  # add gene name to df
  gene_df[i,1] <- gene_new
  i <- i+1
}
# cast blank cells to NA
gene_df[gene_df==""] <- NA
gene_df <- na.omit(gene_df)
# remove duplicate gene names and save as vector
genes_unique <- as.data.frame(sort(unique(gene_df$gene)))


##### visualize GO term gene expressions per cluster  
# EJ: the line graphs are mean of z-scores across genes per cluster and s.e.m. 
# I did not re-normalize these z-scores across clusters here. 
# Re-normalization across clusters was done only for heat maps for visual illustration of the dominant clusters.
# Re-normalization: to show which one is the strongest, I did normalization per row one more time. 
# So if you have 5 values (y1,y2,y3,y4,y5) for a row, I normalized that row zi = (yi-mean(y1,y2,.. y5))/std(y1,y2,..y5). 
# That way, it was clearer to see which cluster has the highest expression for each gene.

### read in DESeq2 normalized values
countsi <- read.csv("~/Downloads/RNAseq/cluster_by_genes/0.3cutoff/DESC/DESeq2/n25.L0.65/all_cells/normalized_counts_from_dds.csv", check.names=FALSE)
  #/work/pi_yingzhang_uri_edu/kdunton/RNAseq/cluster_by_genes/0.3cutoff/DESC/DESeq2/n25.L0.65/all_cells/normalized_counts_from_dds.csv
counts <- countsi
names(counts)[1] <- 'gene'
rownames(counts) <- counts$gene
counts$gene <- NULL


### subset counts to GO genelist 
genelist <- genes_unique
counts <- counts[genelist$gene,]


### calculate z score per gene (row-wise)
# https://stackoverflow.com/questions/34707527/improving-my-r-code-to-calculate-z-score-of-dataframe
counts_z <- (counts-rowMeans(counts))/(rowSds(as.matrix(counts)))[row(counts)]

### read in cluster info
clusters <- read.csv("~/Downloads/RNAseq/cluster_by_genes/0.3cutoff/DESC/clusters_n25.L0.65.csv")
names(clusters)[1] <- 'barcode'
names(clusters)[2] <- 'cluster'

##  remove C4 since that is contaminated celltypes
clusters <- subset(clusters, cluster != '4')
counts_z <- counts_z[,clusters$barcode]

## add 'C' to the beginning of each cluster ID
clusters$cluster <- paste0('C',clusters$cluster)

## add cluster info to counts
counts_z <- as.data.frame(t(counts_z))
counts_z$barcode <- rownames(counts_z)
counts_z <- inner_join(clusters,counts_z, by='barcode')

## reformat for plotting
col_range <- colnames(counts_z)
counts_long <- gather(counts_z, gene, zscore, col_range[3]:col_range[length(col_range)], factor_key=TRUE)
# - data: Data object
# - key: Name of new key column (made from names of data columns)
# - value: Name of new value column
# - ...: Names of source columns that contain values
# - factor_key: Treat the new key column as a factor (instead of character vector)

## get gene average per cluster
counts_avg <- counts_long %>%
  group_by(cluster, gene) %>% 
  summarise_at(vars("zscore"), mean)


## add GO term label to each gene
# some genes appear in multiple terms so they will be repeated and have the same zscore
# https://stackoverflow.com/questions/65948556/r-merge-2-data-frames-one-of-them-with-repeated-measures-keeping-nas-where-ap
join <- merge(df_separate, counts_avg, sort=FALSE)



#### plots
join$cluster <- as.factor(join$cluster)

## LINE PLOT: plot average expression per term per cluster as a line graph (average all genes for a term together)
a1 <- list()
i <- 1
GOterms_order <- c('long-term synaptic potentiation','long-term synaptic depression','dendritic spine morphogenesis',
                   'regulation of RNA splicing','regulation of actin filament-based process')
for(term in GOterms_order){
  print(i)
  print(term)
  # subset to genes of one GO term
  join_sub <- join[join$term == term,]
  # calculate standard error for error bars
  stderror <- function(x) sd(x)/sqrt(length(x))
  std.error <- aggregate(join_sub$zscore, list(join_sub$cluster), FUN=stderror) 
  names(std.error)[1] <- 'cluster'
  names(std.error)[2] <- 'std.error'
  # calculate avg zscore per cluster
  means <- aggregate(join_sub$zscore, list(join_sub$cluster), FUN=mean) 
  names(means)[1] <- 'cluster'
  names(means)[2] <- 'avgZscore'
  # add st.error to the mean data
  data <- merge(means,std.error, by='cluster')
  a1[i] <- list(ggplot(data, aes(x=cluster,y=avgZscore, group=1)) + geom_line() + ggtitle(term) + xlab('') + ylab('Average Z-score') +
                  theme_classic() + geom_hline(yintercept=0, linetype=3, col = 'black') + 
                  geom_errorbar(aes(ymin=avgZscore-std.error, ymax=avgZscore+std.error), width=0.05))
  i <- i + 1
}
# arrange the lineplots together
z1 <- do.call(grid.arrange,c(a1, nrow=1))
# save as svg to get rid of extra legends
svg(filename="~/Downloads/RNAseq/cluster_by_genes/0.3cutoff/GO/clusters/n25.L0.65/x_vs_0/GO_line.svg", width = 12, height = 2)
plot(z1)
dev.off()



## HEATMAP: plot the gene exp per cluster per GO term as heatmap
# subset df by GO term and save the matrices into a list for plotting
mn <- list()
i <- 1
for(term in df$term){
  print(i)
  print(term)
  # subset to genes of one GO term
  join_sub <- join[join$term == term,]
  
  #make into matrix to do a clustered heatmap
  mtx <- join_sub[,c(1,3,4)]
  mtx <- reshape(mtx, idvar = "gene", timevar = "cluster", direction = "wide")
  rownames(mtx) <- mtx$gene
  mtx$gene <- NULL
  mtx <- mtx %>% 
    rename('C0' = zscore.C0,
           'C1' = zscore.C1,
           'C2' = zscore.C2,
           'C3' = zscore.C3,
           'C5' = zscore.C5)
  # reorder
  mtx <- mtx[,c('C0','C1','C2','C3','C5')]
  # normalize the rows (ie across clusters); this is basically doing z-score again to normalize each row of the heatmap to show the trend better
  cal_z_score <- function(x){
    (x - mean(x)) / sd(x)
  }
  mtx_norm <- t(apply(mtx, 1, cal_z_score))
  # add normalized values to list
  mn[i] <- list(mtx_norm)
  i <- i+1
}

# make a list to append heatmaps to https://stackoverflow.com/questions/51629421/combining-pheatmaps-in-r
a <- list()
colors <- viridis(n = 200)
for(i in 1:length(mn)){
  a[i] <- list(pheatmap(mn[[i]], cluster_cols=F, cluster_rows=T, treeheight_row = 0, color=colors, 
                        border_color = NA, cellwidth=15, cellheight=5, fontsize=6, main=GOterms[i]) [[4]])
  #want to cluster the rows so that genes group to show a trend
  #don't want to cluster columns bc that rearranges the clusters in each plot
  #treeheight_row=0 removes the dendrogram from the row clustering
}
# arrange the heatmaps together
z <- do.call(grid.arrange,c(a, nrow=1))
svg(filename="~/Downloads/RNAseq/cluster_by_genes/0.3cutoff/GO/clusters/n25.L0.65/x_vs_0/GO_heatmap.svg", width = 12, height = 12)
plot(z)
dev.off()



#------------- Part C.2. Heatmaps and line plots of GO terms of interest specific to C5 -------------#
### compare terms special to C5
regulation <- 'upreg'
GO <- read.csv(paste0("~/Downloads/RNAseq/cluster_by_genes/0.3cutoff/GO/clusters/n25.L0.65/x_vs_0/GOresults_n25.L0.65_Bio-",regulation,"DEG.csv"), check.names=FALSE)
GOterms <- c('ATP biosynthetic process')
GO <- subset(GO, Description %in% GOterms)

# each cluster might have different genes contributing to a GO term, so combine genes across all clusters per term
for(term in GOterms){
  # split each term into a separate df
  x <- subset(GO, Description == term)
  # subset to the gene column
  x <- x[,'geneID']
  # replace / with ,
  x2 <- str_replace_all(x, c("/" = ","))
  # x2 is a list, where each sublist is the gene list from one cluster. Convert to string
  x2string <- toString(x2)
  # now all clusters are combined
  # remove duplicate genes that appear in multiple clusters
  x2string <- unique(unlist(strsplit(x2string, ",")))
  # becomes a character vector or something
  # put back to string
  x2string <- toString(x2string)
  assign(term,x2string)
}

# combine all terms into one df
df <- data.frame(matrix(0, nrow=length(GOterms), ncol=2))
colnames(df) <- c('term','genes')
df$term <- GOterms

for(term in GOterms){
  df$genes[df$term == term] <- get(term)
}


## expand df so each gene is in its own row but still labeled by GO term
df_separate <- separate_rows(df, genes, sep = ", ", convert = FALSE)
#need the space after the comma so the resulting cells don't leading white space
names(df_separate)[2] <- 'gene' 
# but some genes still have a leading white space:
genes_nowhite <- gsub("* ", "", df_separate$gene)
df_separate$gene <- genes_nowhite


##### visualize GO term gene expressions per cluster  
### read in DESeq2 normalized values
countsi <- read.csv("~/Downloads/RNAseq/cluster_by_genes/0.3cutoff/DESC/DESeq2/n25.L0.65/all_cells/normalized_counts_from_dds.csv", check.names=FALSE)
counts <- countsi
names(counts)[1] <- 'gene'
rownames(counts) <- counts$gene
counts$gene <- NULL


### subset counts to C5-related terms' genes
# C5 terms:
atp <- subset(df_separate, term=='ATP biosynthetic process')
glycolysis <- data.frame(term = 'glycolysis', gene = c('Aldoa','Bpgm','Eno1','Eno2','Gapdh','Gpi1','Hk1','Hk2','Hkdc1','Pfkl','Pfkm','Pgam1','Pkg1','Pklr','Pkm','Tpl1','Aldoart1','Pgam2','Aldoart2'))
NE <- data.frame(term = 'NE receptors', gene = c('Adra1a','Adra1b','Adra1d','Adra2a','Adra2c','Adrb1','Adrb2','Adrb3'))
C5terms <- rbind(atp,glycolysis,NE)
genelist <- C5terms
C5terms_unique <- unique(C5terms$term)
# subset
counts <- counts[genelist$gene,]


### calculate z score per gene (row-wise)
# https://stackoverflow.com/questions/34707527/improving-my-r-code-to-calculate-z-score-of-dataframe
counts_z <- (counts-rowMeans(counts))/(rowSds(as.matrix(counts)))[row(counts)]

### read in cluster info
clusters <- read.csv("~/Downloads/RNAseq/cluster_by_genes/0.3cutoff/DESC/clusters_n25.L0.65.csv")
names(clusters)[1] <- 'barcode'
names(clusters)[2] <- 'cluster'

##  remove C4 since that is contaminated celltypes
clusters <- subset(clusters, cluster != '4')
counts_z <- counts_z[,clusters$barcode]

## add 'C' to the beginning of each cluster ID
clusters$cluster <- paste0('C',clusters$cluster)

## add cluster info to counts
counts_z <- as.data.frame(t(counts_z))
counts_z$barcode <- rownames(counts_z)
counts_z <- inner_join(clusters,counts_z, by='barcode')

## reformat for plotting
col_range <- colnames(counts_z)
counts_long <- gather(counts_z, gene, zscore, col_range[3]:col_range[length(col_range)], factor_key=TRUE)
# - data: Data object
# - key: Name of new key column (made from names of data columns)
# - value: Name of new value column
# - ...: Names of source columns that contain values
# - factor_key: Treat the new key column as a factor (instead of character vector)

## get gene average per cluster
counts_avg <- counts_long %>%
  group_by(cluster, gene) %>% 
  summarise_at(vars("zscore"), mean)


## add GO term label to each gene
# some genes appear in multiple terms so they will be repeated and have the same zscore
# https://stackoverflow.com/questions/65948556/r-merge-2-data-frames-one-of-them-with-repeated-measures-keeping-nas-where-ap
join <- merge(C5terms, counts_avg, sort=FALSE)


#### plots
join$cluster <- as.factor(join$cluster)

## LINE PLOT: plot average expression per term per cluster as a line graph (average all genes for a term together)
a1 <- list()
i <- 1
for(term in C5terms_unique){
  print(i)
  print(term)
  # subset to genes of one GO term
  join_sub <- join[join$term == term,]
  # calculate standard error for error bars
  stderror <- function(x) sd(x)/sqrt(length(x))
  std.error <- aggregate(join_sub$zscore, list(join_sub$cluster), FUN=stderror) 
  names(std.error)[1] <- 'cluster'
  names(std.error)[2] <- 'std.error'
  # calculate avg zscore per cluster
  means <- aggregate(join_sub$zscore, list(join_sub$cluster), FUN=mean) 
  names(means)[1] <- 'cluster'
  names(means)[2] <- 'avgZscore'
  # add st.error to the mean data
  data <- merge(means,std.error, by='cluster')
  a1[i] <- list(ggplot(data, aes(x=cluster,y=avgZscore, group=1)) + geom_line() + ggtitle(term) + xlab('') + ylab('Average Z-score') +
                  theme_classic() + geom_hline(yintercept=0, linetype=3, col = 'black') + 
                  geom_errorbar(aes(ymin=avgZscore-std.error, ymax=avgZscore+std.error), width=0.05))
  i <- i + 1
}
# arrange the lineplots together
z1 <- do.call(grid.arrange,c(a1, nrow=1))
svg(filename="~/Downloads/RNAseq/cluster_by_genes/0.3cutoff/GO/clusters/n25.L0.65/x_vs_0/C5_GOtermLinePlot_zscore.svg", width = 7, height = 2)
plot(z1)
dev.off()



## HEATMAP: plot the gene exp per cluster per GO term as heatmap
# subset df by GO term and save the matrices into a list for plotting
mn <- list()
i <- 1
for(term in C5terms_unique){
  print(i)
  print(term)
  # subset to genes of one GO term
  join_sub <- join[join$term == term,]
  
  #make into matrix to do a clustered heatmap
  mtx <- join_sub[,c(1,3,4)]
  mtx <- reshape(mtx, idvar = "gene", timevar = "cluster", direction = "wide")
  rownames(mtx) <- mtx$gene
  mtx$gene <- NULL
  mtx <- mtx %>% 
    rename('C0' = zscore.C0,
           'C1' = zscore.C1,
           'C2' = zscore.C2,
           'C3' = zscore.C3,
           'C5' = zscore.C5)
  # reorder
  mtx <- mtx[,c('C0','C1','C2','C3','C5')]
  # normalize the rows (ie across clusters)
  cal_z_score <- function(x){
    (x - mean(x)) / sd(x)
  }
  mtx_norm <- t(apply(mtx, 1, cal_z_score))
  # add normalized values to list
  mn[i] <- list(mtx_norm)
  i <- i+1
}

# make a list to append heatmaps to https://stackoverflow.com/questions/51629421/combining-pheatmaps-in-r
a <- list()
colors <- viridis(n = 200)
GOtermsC5 <- c('ATP biosynthetic process', 'glycolysis','NE receptors')
for(i in 1:length(mn)){
  a[i] <- list(pheatmap(mn[[i]], cluster_cols=F, cluster_rows=T, treeheight_row = 0, color=colors, 
                        border_color = NA, cellwidth=20, cellheight=5, fontsize=6, main=GOtermsC5[i]) [[4]])
  #want to cluster the rows so that genes group to show a trend
  #don't want to cluster columns bc that rearranges the clusters in each plot
  #treeheight_row=0 removes the dendrogram from the row clustering
}
# arrange the heatmaps together
z <- do.call(grid.arrange,c(a, nrow=1))
svg(filename="~/Downloads/RNAseq/cluster_by_genes/0.3cutoff/GO/clusters/n25.L0.65/x_vs_0/C5_GOtermHeatmap.svg", width = 10, height = 4)
plot(z)
dev.off()


