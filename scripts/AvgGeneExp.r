# this script calculates mean gene expression per cluster
# https://docs.google.com/document/d/1KnvVFjAHZKsvH2o1ejF_YeyLvNi82fYuZduwp92bzoI/edit#


library(tidyverse)


#### read in normalized counts
counts <- read.csv("/work/pi_yingzhang_uri_edu/kdunton/RNAseq/AIBSmapping/OA/count_matrices/DESCnormalized_counts_L23_0.3.csv", check.names=FALSE)
# need genes in a column for them to be retained after using gather()
counts$gene <- rownames(counts)
# rearrange so gene column is first
counts <- counts[,c(6434,1:6433)]

# format counts to long
col_range <- colnames(counts)
counts_long <- gather(counts, key=cell, value=expression, col_range[2]:col_range[length(col_range)], factor_key=TRUE)
# - data: Data object
# - key: Name of new key column (made from names of data columns)
# - value: Name of new value column
# - ...: Names of source columns that contain values
# - factor_key: Treat the new key column as a factor (instead of character vector)

# add cluster column
clusters <- read.csv("/work/pi_yingzhang_uri_edu/kdunton/RNAseq/cluster_by_genes/0.3cutoff/DESC/clusters_n12.L0.85.csv")
names(clusters)[1] <- 'cell'
names(clusters)[2] <- 'cluster'
counts_long2 <- left_join(counts_long, clusters, by='cell')


genes_unique <- as.data.frame(counts$gene)
colnames(genes_unique)[1] <- 'gene'

mean_df <- data.frame(matrix(0, nrow=0, ncol=3))
colnames(mean_df) <- c('cluster','Avg','gene')
for(gene in genes_unique$gene){
  print(gene)
  # subset counts to a single gene
  one_gene <- counts_long2[counts_long2$gene==gene,]
  # calculate mean per cluster
  means <- one_gene %>%
    group_by(cluster) %>%
    summarise(Avg = mean(expression))
  # add gene column so the loop iterations are discernable
  means$gene <- gene
  # append to df
  mean_df <- rbind(mean_df, means)
}

write.csv(mean_df,"AllGenes_meanPerCluster_n12.L0.85.csv",row.names=FALSE)

print("done")








