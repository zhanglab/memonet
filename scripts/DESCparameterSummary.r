# this script summarizes info for each DESC parameter setting so they can be compared

library(tidyverse)
library(caret)

setwd("~/Desktop/")



### summarize the avg silhouette score per cluster ###
filenames <- Sys.glob("~/Downloads/RNAseq/cluster_by_genes/DESC_parameter_test/test_loop/sil_scores_n*.csv")

overall_sil_score <- data.frame(matrix(0, nrow=length(filenames), ncol=3))
colnames(overall_sil_score) <- c('parameters','overall_silhouette_score','n_clusters')
i <- 1

for(file in filenames){
  print(i)
  ## read in file
  filename <- basename(file)
  #this extracts the actual filename from the path, so you can extract the parameter details
  print(filename)
  parameters <- sub("sil_scores_", "",filename)
  parameters <- sub(".csv", "",parameters)
  sil <- read.csv(file)
  names(sil)[1] <- 'barcode'
  names(sil)[2] <- 'cluster'
  
  ## calculate number of clusters that were generated
  n_clusters <- length(unique(sil$cluster))
  
  ## calculate average silhouette score per cluster
  mean_sil_cl <- sil %>%
    group_by(cluster) %>%
    summarise(Avg = mean(silhouette_samples))
  overall_sil <- mean(mean_sil_cl$Avg)
  
  overall_sil_score[i,1] <- parameters
  overall_sil_score[i,2] <- overall_sil
  overall_sil_score[i,3] <- n_clusters
  i <- i + 1
}


## heatmaps ##
overall_sil_score$n <- sub(".L[[:print:]]*", "",overall_sil_score$parameters)
#remove everything after .L
overall_sil_score$n <- sub("n", "",overall_sil_score$n)
#remove the n prefix
overall_sil_score$L <- sub("[[:print:]]*.L", "",overall_sil_score$parameters)
#keep everything after .L

# L by n, showing sil score
ggplot(overall_sil_score, aes(x=L, y=n, fill=overall_silhouette_score)) +
  geom_tile(color="gray",lwd=0.05, linetype=1) + #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  geom_text(aes(label = round(x=overall_silhouette_score,digits=3)), color = "black", size = 3) +
  scale_fill_gradient(low="white",high="blue") +
  xlab("Louvain resolution") + ylab("nearest neighbors")
ggsave(file="heatmap_silScore.png")
ggsave(file="heatmap_silScore.svg")

# L by n, showing n_clusters
ggplot(overall_sil_score, aes(x=L, y=n, fill=n_clusters)) +
  geom_tile(color="gray",lwd=0.05, linetype=1) + #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  geom_text(aes(label = n_clusters), color = "black", size = 3) +
  scale_fill_gradient(low="white",high="blue") +
  xlab("Louvain resolution") + ylab("nearest neighbors")
ggsave(file="heatmap_nclusters.png")
ggsave(file="heatmap_nclusters.svg")





### determine if there's a prominent cluster that most AIBS cells map to ###
filenames <- Sys.glob("~/Downloads/RNAseq/AIBSmapping/AO/parameter_loop/*prediction_scores.csv")
AOsummary <- data.frame(matrix(0, nrow=length(filenames), ncol=3))
colnames(AOsummary) <- c('parameters','largest_cluster_AIBS','second_largest_cluster_AIBS')
i <- 1
for(file in filenames){
  filename <- basename(file)
  #this extracts the actual filename from the path, so you can extract the parameter details
  print(filename)
  # get parameter setting
  parameter <- sub("_prediction_scores.csv","",filename)
  
  data <- read.csv(file)
  
  ## get observed counts of AIBS cells per cluster
  table <- as.data.frame(table(data$predicted.id))
  colnames(table) <- c('cluster','Freq')
  table$percent <- round(100 * table$Freq/sum(table$Freq), 2)
  
  largest <- max(table$percent)
  n <- length(table$percent)
  second_largest <- sort(table$percent,partial=n-1)[n-1]
  
  AOsummary[i,'parameters'] <- parameter
  AOsummary[i,'largest_cluster_AIBS'] <- largest
  AOsummary[i,'second_largest_cluster_AIBS'] <- second_largest
  
  i <- i+1
}

# combine AO results with silhouette scores
AOsummary <- inner_join(AOsummary, overall_sil_score, by="parameters")
# calculate difference between the 2 largest aibs-mapped clusters, to see which parameter settings have a prominent aibs-mapped cluster
AOsummary$diff_AIBS_mapping <- AOsummary$largest_cluster_AIBS - AOsummary$second_largest_cluster_AIBS



## heatmaps ##
#plot diff bt largest and second largest aibs-mapped cluster
ggplot(AOsummary, aes(x=L, y=n, fill=diff_AIBS_mapping)) +
  geom_tile(color="gray",lwd=0.05, linetype=1) + #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  geom_text(aes(label = round(x=diff_AIBS_mapping,digits=3)), color = "black", size = 3) +
  scale_fill_gradient(low="white",high="blue") +
  xlab("Louvain resolution") + ylab("nearest neighbors")
ggsave(file="heatmap_aibsPercentageDiff.png")
ggsave(file="heatmap_aibsPercentageDiff.svg")






### check how well the original clusters are conserved in each parameter clustering result
filenames <- Sys.glob("~/Downloads/RNAseq/cluster_by_genes/DESC_parameter_test/test_loop/clusters*.csv")
summary <- data.frame(matrix(0, nrow=length(filenames), ncol=9))
colnames(summary) <- c('parameters','C0_conserved_value','C1_conserved_value','C2_conserved_value','C3_conserved_value','C4_conserved_value','C5_conserved_value','C6_conserved_value','C7_conserved_value')
i <- 1
for(file in filenames){
  filename <- basename(file)
  #this extracts the actual filename from the path, so you can extract the parameter details
  print(filename)
  parameter <- sub("clusters_dec22.","",filename)
  parameter <- sub(".csv","",parameter)
  
  #### create confusion matrix ####
  x <- read.csv("~/Downloads/memonet/cluster_by_genes/NathanClassifierTrVsCtrl_3.8.22/clusters3.24.22_TrCtrlTogether.csv")
  colnames(x)[2] <- 'x_cluster'
  y <- read.csv(file)
  colnames(y)[2] <- 'y_cluster'
  
  #join files together; use inner_join to only keep cells that appear in both
  data <- inner_join(y,x,by='X')
  
  #testing across datasets can leave you with different labels in each column, so use union
  #it adds the missing labels to each dataset and gives them zeros
  u <- union(data$y_cluster,data$x_cluster)
  t <- table(factor(data$y_cluster,u), factor(data$x_cluster,u))
  confusion <- confusionMatrix(t)
  #view the matrix stored in confusion variable
  confusion_df2 <- as.data.frame(confusion$table)
  confusion_df2 <- confusion_df2 %>% rename(analysis_y=Var1, analysis_x=Var2)
  
  ## uncomment this if you want normalized by column
  #normalize by column (each column of confuison matrix should add up to 1)
  wide_confusion <- as.data.frame(pivot_wider(confusion_df2,names_from = analysis_x,values_from = Freq))
  #pivot so that it resembles a confusion matrix visualization, so that you can scale the values by column
  #set the analysis_y to rownames so all the values are numbers (in order for scale() to work)
  row.names(wide_confusion) <- wide_confusion$analysis_y
  wide_confusion$analysis_y <- NULL
  #get the sum of each column
  c1 <- colSums(wide_confusion)
  scaled_confusion <- as.data.frame(scale(wide_confusion, center=FALSE, scale=c1))
  #pass the colSums variable into scale so that each value is divided by the sum
  #change back to long format for ggplot
  #first have to change rownames (analysis_y) back into a column
  scaled_confusion <- rownames_to_column(scaled_confusion,var="analysis_y")
  scaled_confusion <- pivot_longer(scaled_confusion, cols=!analysis_y, names_to ="analysis_x", values_to="Freq")
  #use all columns except analysis_y to put into a column, which is set to "analysis_x"
  #the values get moved to one column called Freq
  scaled_confusion[is.na(scaled_confusion)] = 0
  
  
  #### summarize how each cluster is conserved
  # calculate max percentage value per cluster
  conserve <- scaled_confusion %>%
    group_by(analysis_x) %>%
    summarise(max = max(Freq))
  
  summary[i,'parameters'] <- parameter
  summary[i,'C0_conserved_value'] <- conserve$max[conserve$analysis_x=='0']
  summary[i,'C1_conserved_value'] <- conserve$max[conserve$analysis_x=='1']
  summary[i,'C2_conserved_value'] <- conserve$max[conserve$analysis_x=='2']
  summary[i,'C3_conserved_value'] <- conserve$max[conserve$analysis_x=='3']
  summary[i,'C4_conserved_value'] <- conserve$max[conserve$analysis_x=='4']
  summary[i,'C5_conserved_value'] <- conserve$max[conserve$analysis_x=='5']
  summary[i,'C6_conserved_value'] <- conserve$max[conserve$analysis_x=='6']
  summary[i,'C7_conserved_value'] <- conserve$max[conserve$analysis_x=='7']
  i <- i+1
}

# join with the other data: AO results, sil scores
allData <- inner_join(AOsummary, summary, by= 'parameters')
allData <- allData[,c(1,6,7,4,2,3,8,9:16)]
write.csv(allData, "parameterSettingsSummary.csv", row.names = FALSE)  



