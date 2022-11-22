# this script plots the train and control proportion of MEMONET cells in the 3 L2/3 subtypes and classifier clusters, respectively



library(tidyverse)
library(matrixStats) #for rowSds



#### read in data ####
# choose to calculate the L2/3 proportion or cluster proportion below

#------ for train/control proportion across the 3 L2/3 subtypes
clusters <- read.csv("~/Downloads/RNAseq/AIBSmapping/OA/prediction_scores.csv")
clusters <- clusters[,c('X','AIBS_predicted_label')]
# add sample and stim column
clusters$sample <- sub("[[:print:]]*-", "",clusters$X)
clusters$sample <- as.character(clusters$sample)
clusters <- clusters %>%
  mutate(stim = recode(sample,
                       "1" = "control",
                       "2" = "train",
                       "3" = "control",
                       "4" = "train",
                       "5" = "control",
                       "6" = "train"))
# in this case, the cluster column we want to look at is the aibs predicted label
clusters <- clusters %>% rename(cluster = AIBS_predicted_label)
# subset to L2/3 cells
L23 <- read.csv("~/Downloads/RNAseq/AIBSmapping/OA/barcode_files/L23barcodes-fromAIBS_0.3.csv")
names(L23)[1] <- 'X'
clusters <- subset(clusters, X %in% L23$X)
table(clusters$cluster,clusters$stim)

outfile <- 'L23subclass_'
#------


#------ for train/control proportion across L2/3 clusters
# clusters:
clusters <- read.csv("~/Downloads/RNAseq/cluster_by_genes/0.3cutoff/DESC/clusters_n25.L0.65.csv")
colnames(clusters)[2] <- 'cluster'

# add sample and stim column
clusters$sample <- sub("[[:print:]]*-", "",clusters$X)
clusters$sample <- as.character(clusters$sample)
clusters <- clusters %>%
  mutate(stim = recode(sample,
                       "1" = "control",
                       "2" = "train",
                       "3" = "control",
                       "4" = "train",
                       "5" = "control",
                       "6" = "train"))
outfile <- 'clusters_'
#-----



#### calculate the observed proportion of train and control ####
# here, cluster can refer to the L2/3 clusters OR the 3 L2/3 subtypes

# table gives the frequency of train and control cells in each cluster. All you need is to add a proportion column
table <- as.data.frame(table(clusters$cluster, clusters$stim))
colnames(table) <- c('cluster','stim','Freq')
# calculate the total train and control cells
table2 <- table %>%
  group_by(stim) %>%
  summarise(Total = sum(Freq))
# join df so the sum is attached to the combos of cluster and stim
table3 <- inner_join(table, table2, by = 'stim')
# calculate proportion
table4 <- table3 %>% mutate(Prop = Freq / Total)
# check that the prop of train,control adds to 1
check_table <- table4 %>%
  group_by(stim) %>%
  summarise(Sum = sum(Prop))
table4 <- table4[,c('cluster','stim','Prop')]

# pivot table4 so you can get the difference of train and control proportions (for calculating stats)
table4b <- table4 %>% pivot_wider(names_from = stim, values_from = Prop)
# calculate train-control proportions
table4_observed <- table4b %>% mutate(mean = train - control)
table4_observed <- table4_observed[,c('cluster','mean')]
## format with other columns that the iteration data will have
table4_observed$o_or_e <- 'Observed'
table4_observed$stError <- 'NA'
# reorder 
table4_observed <- table4_observed[,c(1,3,2,4)]



#### plot observed train and control proportions ####
# choose which to plot: L2/3 subclass OR cluster

## train/ctrl prop per L2/3 subclass
ggplot(data=table4, aes(x=cluster, y=Prop, fill=stim)) +
  geom_bar(stat='identity', position='dodge', width=0.5) +
  ylab('Proportion') + xlab('L2/3 subclass') +
  theme_classic() +
  scale_fill_manual(values=c("Red","Black"), name='') +
  ggtitle('Proportion of train and control per L2/3 subclass') +
  theme(axis.text = element_text(size = 12)) + 
  theme(legend.text = element_text(size = 12)) +
  theme(legend.title = element_text(size = 14)) +
  theme(axis.title = element_text(size = 14)) +
  theme(plot.title = element_text(size = 16))
ggsave(file= paste0("~/Downloads/RNAseq/AIBSmapping/OA/",outfile,"tr_ctrl_prop_hist.svg"))
ggsave(file= paste0("~/Downloads/RNAseq/AIBSmapping/OA/",outfile,"tr_ctrl_prop_hist.png"))


## train/ctrl prop per cluster
ggplot(data=table4, aes(x=cluster, y=Prop, fill=stim)) +
  geom_bar(stat='identity', position='dodge', width=0.5) +
  ylab('Proportion') + xlab('Cluster') +
  theme_classic() +
  scale_fill_manual(values=c("Red","Black"), name='') +
  ggtitle('Proportion of train and control per cluster') +
  theme(axis.text = element_text(size = 12)) + 
  theme(legend.text = element_text(size = 12)) +
  theme(legend.title = element_text(size = 14)) +
  theme(axis.title = element_text(size = 14)) +
  theme(plot.title = element_text(size = 16))
ggsave(file= paste0("~/Downloads/RNAseq/cluster_by_genes/0.3cutoff/DESC/figures/",outfile,"tr_ctrl_prop_hist.svg"))
ggsave(file= paste0("~/Downloads/RNAseq/cluster_by_genes/0.3cutoff/DESC/figures/",outfile,"tr_ctrl_prop_hist.png"))




#### now for stats - randomly shuffle the control and train labels ####
set.seed(12345)
#get vector of stim column
condition <- clusters$stim
#shuffle the vector - each cell gets a shuffled control/train label so the result will have a row for each cell (6506) and column for each iteration (1000)
n <- 1000
rand_labels <- data.frame(matrix(0, nrow=length(condition), ncol=n))
for(i in 1:n){
  rand_labels[,i] <- sample(condition)
}

# add cluster info to rand_labels
rand_labels$cluster <- clusters$cluster
# reorder so the cluster column is first and then the shuffles follow
rand_labels <- rand_labels[,c(1001,1:1000)]
## want to have a row for each cluster so need to group the cells of each cluster together
rand_labels_long <- gather(rand_labels, key=iteration, value=stim, X1:X1000, factor_key=TRUE)
# - data: Data object
# - key: Name of new key column (made from names of data columns)
# - value: Name of new value column
# - ...: Names of source columns that contain values
# - factor_key: Treat the new key column as a factor (instead of character vector)
# now convert the 'control' and 'train' labels to counts for each cluster per iteration
rand_labels_counts <- as.data.frame(table(rand_labels_long$cluster, rand_labels_long$iteration ,rand_labels_long$stim))
# split into control and train datasets and then remove the stim column
rand_train <- subset(rand_labels_counts, Var3 %in% 'train')
rand_train$Var3 <- NULL
rand_ctrl <- subset(rand_labels_counts, Var3 %in% 'control')
rand_ctrl$Var3 <- NULL
# finally change from long to wide
rand_train <- spread(data=rand_train,key=Var2,value=Freq)
rand_ctrl <- spread(data=rand_ctrl,key=Var2,value=Freq)
rand_train$Var1 <- NULL
rand_ctrl$Var1 <- NULL

## calculate proportions from the raw counts
#get total count of control cells
ccnt <- sum(table[table$stim=='control',]$Freq)
#get total count of trained cells
tcnt <- sum(table[table$stim=='train',]$Freq)
#divide each train count by the total train
rand_train_prop <- rand_train[,] / tcnt
#divide each control count by the total control
rand_ctrl_prop <- rand_ctrl[,] / ccnt
#check that each column adds to 1
check_prop_tr <- as.data.frame(colSums(rand_train_prop))
check_prop_ctrl <- as.data.frame(colSums(rand_ctrl_prop))

## calculate proportion difference: train - control
rand_prop_diff <- rand_train_prop - rand_ctrl_prop




#### calculate p-values ####
cluster_ids <- sort(unique(clusters$cluster))
# add observed values to the expected
#- first format expected data
rand_prop_diff2 <- rand_prop_diff
rand_prop_diff2$cluster <- cluster_ids
rand_prop_diff2$o_or_e <- 'Expected'
rand_prop_diff2 <- rand_prop_diff2[,c(1001,1002,1:1000)]
rand_prop_diff2 <- gather(rand_prop_diff2, key=iteration, value=PropDiff, X1:X1000, factor_key=TRUE)
rand_prop_diff2$iteration <- NULL
#- format observed data
table4_observed2 <- table4_observed[,1:3]
table4_observed2 <- table4_observed2 %>% rename(PropDiff = mean)
#- join
obs_expected2 <- rbind(table4_observed2, rand_prop_diff2)

# calculate p-values
p.values <- data.frame(matrix(0, nrow=length(cluster_ids), ncol=1))
colnames(p.values) <- 'pvalue'
i <- 0
  #counter for adding to rows of p.values
for(cl in cluster_ids){
  # subset by cluster
  subset_df <- subset(obs_expected2, cluster %in% cl)
  # take the absolute value
  subset_df$PropDiff <- abs(subset_df$PropDiff)
  # order p-values smallest to largest
  subset_df <- subset_df[order(subset_df$PropDiff),]
  # get row index for the observed value
  row_index <- which(subset_df$o_or_e=='Observed')
  
  # calculate p-value
    # p < (b+1)/(m+1)
    # b = the number of |shuffled value| > |observation|
    # m = number of shuffles
  m <- 1000
  b <- 1001 - row_index
  pvalue <- (b+1)/(m+1)
  
  # save to df
  i <- i+1
  cl_index <- i
  print(cl_index)
  p.values[cl_index,1] <- pvalue
}

# rename with the correct cluster numbers, since the loop started from 1 not 0
p.values$cluster <- cluster_ids
p.values <- p.values[,c(2,1)]






#### calculate mean and standard error of the null distribution (iteration data) ####
# since there are 1000 iterations, need to get average proportion difference
rand_stats <- data.frame(matrix(0, nrow=length(cluster_ids), ncol=0))
rand_stats$cluster <- cluster_ids
rand_stats$mean <- rowMeans(rand_prop_diff[,])
rand_stats$stError <- rowSds(as.matrix(rand_prop_diff[,]))
# add o_or_e column
rand_stats$o_or_e <- 'Expected' 
# reorder
rand_stats <- rand_stats[,c(1,4,2,3)]




### summary table 
summary <- rbind(table4_observed, rand_stats)
summary$stError <- NULL
summary <- spread(summary, key=o_or_e, value=mean)
summary <- summary %>% rename(ObservedPropDiff=Observed,
                              ExpectedPropDiff=Expected)
summary <- summary[,c(1,3,2)]
summary$standardError <- rand_stats$stError
p.values2 <- as.data.frame(p.values[,"pvalue"])
colnames(p.values2) <- "pvalue"
summary <- cbind(summary, p.values2)
summary <- summary[,c(1,5,2,3,4)]


# for L2/3 subclass:
write.csv(summary, paste0("~/Downloads/RNAseq/AIBSmapping/OA/",outfile,"tr_ctrl_prop_table.csv"), row.names = FALSE)

# for cluster:
write.csv(summary, paste0("~/Downloads/RNAseq/cluster_by_genes/0.3cutoff/DESC/figures/",outfile,"tr_ctrl_prop_table.csv"), row.names = FALSE)

