# visualize/summarize DE results: line plot of IEGs per cluster


library(tidyverse)
library(schoolmath) #for is.positive
library(RColorBrewer) #for colorblind palette 
library(viridis) # for heatmap colors
library(matrixStats) # for rowSds
library(pheatmap)
library(svglite)



##### IEG line plots across cluster ##### 

#--------- Part A. FIGURE 4a. DE of each cluster vs others. Plot z-score expression of IEGs across clusters, indicating DE significance result with open/closed dots ###
### read in DEG data
# ref others
DEG <- read.csv("~/Downloads/RNAseq/cluster_by_genes/0.3cutoff/DESC/DESeq2/n25.L0.65/all_cells/DEGstats_allGenes.csv")

# rename the results column values
DEG <- DEG %>% mutate(result = recode(result,
                                      'upregulated'= 'up-regulated',
                                      'downregulated' = 'down-regulated'))
# don't subset to padj < 0.05, bc you want to plot gene as a line across all clusters and those that aren't sig will be open circle

# subset to IEGs
genelist <- c('Arc','Cebpb','Egr2','Fos','Fosb','Jun','Junb','Jund','Nefm','Npas4','Nptx2','Nr4a2','Nr4a3','Scg2','Syt4','Vgf')
DEG <- subset(DEG, gene %in% genelist)

# add column to state if gene is significant or not
DEG <- DEG %>% mutate(significance = case_when(padj <= 0.05 ~ "sig", padj > 0.05 ~ "not_sig"))
# remove any NA values
DEG <- na.omit(DEG)
DEG$significance <- as.factor(DEG$significance)
DEG$cluster <- as.factor(DEG$cluster)


### read in DESeq2 normalized values
countsi <- read.csv("~/Downloads/RNAseq/cluster_by_genes/0.3cutoff/DESC/DESeq2/n25.L0.65/all_cells/normalized_counts_from_dds.csv", check.names=FALSE)
counts <- countsi
names(counts)[1] <- 'gene'
rownames(counts) <- counts$gene
counts$gene <- NULL

### subset counts to IEGs
genelist <- c('Arc','Cebpb','Egr2','Fos','Fosb','Jun','Junb','Jund','Nefm','Npas4','Nptx2','Nr4a2','Nr4a3','Scg2','Syt4','Vgf')
counts <- counts[genelist,]

### calculate z score per gene (row-wise)
# https://stackoverflow.com/questions/34707527/improving-my-r-code-to-calculate-z-score-of-dataframe
counts_z <- (counts-rowMeans(counts))/(rowSds(as.matrix(counts)))[row(counts)]

### read in cluster info
clusters <- read.csv("~/Downloads/RNAseq/cluster_by_genes/0.3cutoff/DESC/clusters_n25.L0.65.csv")
names(clusters)[1] <- 'barcode'
names(clusters)[2] <- 'cluster'

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

## add DEG result column
DEG2 <- DEG[,c(1,7,9)]
DEG2$significance <- as.character(DEG2$significance)
DEG2$cluster <- as.character(DEG2$cluster)
counts_avg$cluster <- as.character(counts_avg$cluster)
counts_avg <- inner_join(counts_avg, DEG2, by=c('cluster','gene'))

## plot
#print the names of the colors so you can remove the light yellow one
brewer.pal(n = 12, name = "Paired")
colors_colorblind <- c("#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C","#FDBF6F","#FF7F00","#CAB2D6","#6A3D9A","#B15928")
# add additional colors
colors <- c(colors_colorblind,'azure4','aquamarine','darkorchid','blue','deeppink')

ggplot(counts_avg, aes(x=cluster, y=zscore, color=gene, shape=significance, group=gene)) +
  #group=gene is what tells ggplot which dots to connect with a line
  geom_point(size=2) +
  scale_shape_manual(values = c(21,19)) +
  stat_summary(fun='identity', geom="line") +
  scale_color_manual(values=colors) +
  geom_hline(yintercept=0, linetype=2, col = 'black') +
  ylab('Average Z-score') + xlab('Cluster') +
  theme_classic() +
  ggtitle(paste0('DE results: IEGs')) +
  theme(axis.text = element_text(size = 18)) +
  theme(axis.title = element_text(size = 18)) +
  theme(plot.title = element_text(size = 20)) +
  theme(legend.text = element_text(size = 16)) +
  theme(legend.title = element_text(size = 18)) 
ggsave("~/Downloads/RNAseq/cluster_by_genes/0.3cutoff/DESC/DESeq2/n25.L0.65/all_cells/IEG_lineplot_zscore.svg")
ggsave("~/Downloads/RNAseq/cluster_by_genes/0.3cutoff/DESC/DESeq2/n25.L0.65/all_cells/IEG_lineplot_zscore.png", dpi=300)



#--------- Part B. FIGURE S4. DE of each cluster vs C0. IEGs plotted with LFC values 
### read in DEG data
# ref C0; use DEGstats file with no padj cutoff:
DEG <- read.csv("~/Downloads/RNAseq/cluster_by_genes/0.3cutoff/DESC/DESeq2/n25.L0.65/x_vs_0/all_cells/DEGstats_allGenes.csv")

# rename the results column values
DEG <- DEG %>% mutate(result = recode(result,
              'upregulated'= 'up-regulated',
              'downregulated' = 'down-regulated'))
# don't subset to padj < 0.05, bc you want to plot gene as a line across all clusters and those that aren't sig will be open circle

# subset to IEGs
genelist <- c('Arc','Cebpb','Egr2','Fos','Fosb','Jun','Junb','Jund','Nefm','Npas4','Nptx2','Nr4a2','Nr4a3','Scg2','Syt4','Vgf')
DEG <- subset(DEG, gene %in% genelist)

# add column to state if gene is significant or not
DEG <- DEG %>% mutate(significance = case_when(padj <= 0.05 ~ "sig", padj > 0.05 ~ "not_sig"))
# remove any NA values
DEG <- na.omit(DEG)
DEG$significance <- as.factor(DEG$significance)
DEG$cluster <- as.factor(DEG$cluster)



### plot ###

#print the names of the colors so you can remove the light yellow one (hard to see)
brewer.pal(n = 12, name = "Paired")
colors_colorblind <- c("#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C","#FDBF6F","#FF7F00","#CAB2D6","#6A3D9A","#B15928")
# add additional colors
colors <- c(colors_colorblind,'azure4','aquamarine','darkorchid','blue','deeppink')

# plot logFoldChange values -- LINE PLOT
# IEGs
ggplot(DEG, aes(x=cluster, y=log2FoldChange, color=gene, shape=significance, group=gene)) +
    #group=gene is what tells ggplot which dots to connect with a line
  geom_point(size=2) +
  scale_shape_manual(values = c(21,19)) +
  stat_summary(fun='identity', geom="line") +
    #geom_line won't work for factors, but clusters need to be factors or else cluster 4 will be on the plot (which doesn't make sensse when showing DE results in ref to 4)
  scale_color_manual(values=colors) +
  geom_hline(yintercept=0, linetype=2, col = 'black') +
  ylab('Log2FoldChange') + xlab('Cluster') +
  theme_classic() +
  ggtitle(paste0('DE results: IEGs')) +
  theme(axis.text = element_text(size = 18)) +
  theme(axis.title = element_text(size = 18)) +
  theme(plot.title = element_text(size = 20)) +
  theme(legend.text = element_text(size = 16)) +
  theme(legend.title = element_text(size = 18)) 
ggsave("~/Downloads/RNAseq/cluster_by_genes/0.3cutoff/DESC/DESeq2/n25.L0.65/x_vs_0/all_cells/IEG_lineplot_refC0.svg")
ggsave("~/Downloads/RNAseq/cluster_by_genes/0.3cutoff/DESC/DESeq2/n25.L0.65/x_vs_0/all_cells/IEG_lineplot_refC0.png", dpi=300)









