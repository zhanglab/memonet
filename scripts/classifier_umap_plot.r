# label each cluster on the plot rather than in legend


library(tidyverse)


########### load umap scatterpoints and cluster info #########
clusters <- read.csv("~/Downloads/RNAseq/cluster_by_genes/0.3cutoff/DESC/clusters_n25.L0.65.csv")
scatterpoints <- read.csv("~/Downloads/RNAseq/cluster_by_genes/0.3cutoff/DESC/umap_n25.L0.65.csv")


## edit format of scatterpoints
# rename columns so the x and y coordinates for umap are V1 and V2, respectively
scatterpoints <- rename(scatterpoints, V1 = X0, V2 = X1)
# remove first column 'X' (just adds up from 0)
scatterpoints <- subset(scatterpoints, select = -X)


## edit format of clusters
# rename column of barcodes from 'X' to 'barcodes'; also rename cluster column
names(clusters)[1] <- 'barcodes'
names(clusters)[2] <- 'cluster'
clusters$sample <- sub("[[:print:]]*-", "",clusters$barcodes)
clusters$sample <- as.character(clusters$sample)
clusters <- clusters %>%
  mutate(stim = recode(sample,
                       "1" = "control",
                       "2" = "train",
                       "3" = "control",
                       "4" = "train",
                       "5" = "control",
                       "6" = "train"))


## join scatterpoints and clusters
# join by row index since only clusters has the barcode info
# add rowname column to clusters and scatterpoints so the join can use those columns
clusters <- rownames_to_column(clusters)
scatterpoints <- rownames_to_column(scatterpoints)
# use an inner join since both columns are equal
metadata <- inner_join(clusters, scatterpoints, by = c("rowname" = "rowname"))
# remove rowname column
metadata$rowname <- NULL
metadata$cluster <- as.factor(metadata$cluster)


### Plot clusters colored by cluster
ggplot(metadata, aes(V1,V2, color = cluster)) +
  geom_point(size = 0.8) + #alpha=0.5
  guides(color = guide_legend(override.aes = list(size=2))) +
  labs(x='UMAP_1', y='UMAP_2', color = '') +
  theme(legend.text=element_text(size=12)) +
  theme_classic()
ggsave(file= "~/Downloads/RNAseq/cluster_by_genes/0.3cutoff/DESC/figures/cluster_umap.svg",width = 10,height=7.54)





