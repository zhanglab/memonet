# this script determines whether there are glial contaminants, and whether they are removed at a certain prediction score cutoff value

library(tidyverse)



#---------- Part A: choose a L2/3 prediction score cutoff -----------#
L23_all <- read.csv("~/Downloads/RNAseq/AIBSmapping/OA/barcode_files/L23barcodes-fromAIBS_0.csv")
names(L23_all)[1] <- 'X'

### read in DESC-normalized data for all L2/3 cells (6507)
countsi <- read.csv("~/Downloads/RNAseq/AIBSmapping/OA/count_matrices/DESCnormalized_counts_L23_0.csv", check.names=FALSE)
names(countsi)[1] <- 'gene'
rownames(countsi) <- countsi$gene
countsi$gene <- NULL
counts <- countsi 


### loop through cutoffs and add mean expression to df 
cutoffs <- seq(0,0.5,by=0.01)
cutoffs_exp <- data.frame(matrix(0, nrow=length(cutoffs), ncol=8))
colnames(cutoffs_exp) <- c('cutoff','Cux2','Otof','Rtn4rl1','Slc30a3','Cacna2d3','Mertk', 'n_cells_removed')
cutoffs_exp$cutoff <- cutoffs
for(cutoff in cutoffs){
  i <- cutoff
  print(i)
  # read in the cutoff barcode file
  L23_cutoff <- read.csv(paste0("~/Downloads/RNAseq/AIBSmapping/OA/barcode_files/L23barcodes-fromAIBS_",i,".csv"))
  names(L23_cutoff)[1] <- 'X'
  
  # get the cells that are removed with the cutoff
  removed <- subset(L23_all, !X %in% L23_cutoff$X)
  
  # cutoffs that are too small won't remove any cells, so for these you need to set the values to NA:
  if (length(removed$X) == 0) {
    cutoffs_exp$n_cells_removed[cutoffs_exp$cutoff == i] <- 0
    cutoffs_exp$Cux2[cutoffs_exp$cutoff == i] <- NA
    cutoffs_exp$Mertk[cutoffs_exp$cutoff == i] <- NA
    cutoffs_exp$Otof[cutoffs_exp$cutoff == i] <- NA
    cutoffs_exp$Rtn4rl1[cutoffs_exp$cutoff == i] <- NA
    cutoffs_exp$Slc30a3[cutoffs_exp$cutoff == i] <- NA
    cutoffs_exp$Cacna2d3[cutoffs_exp$cutoff == i] <- NA
    
  # or if only 1 cell is removed, you can't calculate the average because there is just one value for each gene:
  } else if (length(removed$X) == 1) {
    # subset to removed cells 
    counts_removed <- as.data.frame(counts[,removed$X]) 
    # it loses the rownames so add them again
    rownames(counts_removed) <- rownames(counts)
    ncells <- length(removed$X)
    
    # add each gene's single expression value to df
    cutoffs_exp$Cux2[cutoffs_exp$cutoff == i] <- counts_removed['Cux2',1]
    cutoffs_exp$Mertk[cutoffs_exp$cutoff == i] <- counts_removed['Mertk',1]
    cutoffs_exp$Otof[cutoffs_exp$cutoff == i] <- counts_removed['Otof',1]
    cutoffs_exp$Rtn4rl1[cutoffs_exp$cutoff == i] <- counts_removed['Rtn4rl1',1]
    cutoffs_exp$Slc30a3[cutoffs_exp$cutoff == i] <- counts_removed['Slc30a3',1]
    cutoffs_exp$Cacna2d3[cutoffs_exp$cutoff == i] <- counts_removed['Cacna2d3',1]
    cutoffs_exp$n_cells_removed[cutoffs_exp$cutoff == i] <- ncells
 
    # for all other scenarios, you can calculate the average:   
  } else {
    # subset to removed cells 
    counts_removed <- counts[,removed$X] 
    ncells <- length(removed$X)
    
    # get mean expression for the removed cells
    counts_removed$mean <- rowMeans(counts_removed)
    counts_removed_markers <- counts_removed[c('Cux2','Mertk','Otof','Rtn4rl1','Slc30a3','Cacna2d3'),]
    Cux2 <- counts_removed_markers['Cux2','mean']
    Mertk <- counts_removed_markers['Mertk','mean']
    Otof <- counts_removed_markers['Otof','mean']
    Rtn4rl1 <- counts_removed_markers['Rtn4rl1','mean']
    Slc30a3 <- counts_removed_markers['Slc30a3','mean']
    Cacna2d3 <- counts_removed_markers['Cacna2d3','mean']
    
    # add to df
    cutoffs_exp$Cux2[cutoffs_exp$cutoff == i] <- Cux2
    cutoffs_exp$Mertk[cutoffs_exp$cutoff == i] <- Mertk
    cutoffs_exp$Otof[cutoffs_exp$cutoff == i] <- Otof
    cutoffs_exp$Rtn4rl1[cutoffs_exp$cutoff == i] <- Rtn4rl1
    cutoffs_exp$Slc30a3[cutoffs_exp$cutoff == i] <- Slc30a3
    cutoffs_exp$Cacna2d3[cutoffs_exp$cutoff == i] <- Cacna2d3
    cutoffs_exp$n_cells_removed[cutoffs_exp$cutoff == i] <- ncells
  }
}
write.csv(cutoffs_exp, "~/Downloads/RNAseq/AIBSmapping/OA/barcode_files/markerGeneExp_predictionCutoffs.csv", row.names = FALSE)





#---------- Part B: investigate C4 cell type contamination -----------#
##### C4 has GO terms characterizing other cell types- is C4  glial contaminants? #####
## read in normalized expression counts
# DESC-normalized counts:
countsi <- read.csv("~/Downloads/RNAseq/AIBSmapping/OA/count_matrices/DESCnormalized_counts_L23_0.3.csv", check.names=FALSE)
counts <- countsi
counts <- as.data.frame(t(counts))
counts$barcode <- rownames(counts)

## read in cluster info
clusters <- read.csv("~/Downloads/RNAseq/cluster_by_genes/0.3cutoff/DESC/clusters_n25.L0.65.csv")
names(clusters)[1] <- 'barcode'
names(clusters)[2] <- 'cluster'

## join cluster info to counts
counts <- inner_join(clusters,counts,by='barcode')
rownames(counts) <- counts$barcode
counts$barcode <- NULL

## calculate avg expression of markers
marker_exp <- data.frame(matrix(0, nrow=length(unique(clusters$cluster)), ncol=1))
colnames(marker_exp) <- c('cluster')
marker_exp$cluster <- unique(clusters$cluster)
markers <- c('Cux2','Otof','Rtn4rl1','Slc30a3','Cacna2d3','Mertk')

for(marker in markers){
  avg <- counts %>%
    group_by(cluster) %>%
    summarise(Avg = mean(get(marker)))
  names(avg)[2] <- marker
  marker_exp <- inner_join(avg,marker_exp, by='cluster')
}
# order columns
marker_exp <- marker_exp[,c('cluster','Cux2','Otof','Rtn4rl1','Slc30a3','Cacna2d3','Mertk')]
write.csv(marker_exp, "~/Downloads/RNAseq/cluster_by_genes/0.3cutoff/DESC/L23_glial_marker_exp.csv", row.names = FALSE)


