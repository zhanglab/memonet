# this script summarizes info for each DESC parameter setting so they can be compared

library(tidyverse)
library(caret)

setwd("~/Desktop/")



### summarize the avg silhouette score per cluster ###
filenames <- Sys.glob("~/Downloads/RNAseq/cluster_by_genes/DESC_parameter_test/sil_scores_n*.csv")

overall_sil_score <- data.frame(matrix(0, nrow=length(filenames), ncol=2))
colnames(overall_sil_score) <- c('parameters','overall_silhouette_score')
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

  ## calculate average silhouette score per cluster
  mean_sil_cl <- sil %>%
    group_by(cluster) %>%
    summarise(Avg = mean(silhouette_samples))
  overall_sil <- mean(mean_sil_cl$Avg)
  
  overall_sil_score[i,1] <- parameters
  overall_sil_score[i,2] <- overall_sil
  i <- i + 1
}


## heatmap ##
# plot Louvain_res by n_neighbors, showing sil score

# get n_neighbors parameter value
overall_sil_score$n <- sub(".L[[:print:]]*", "",overall_sil_score$parameters)
  #remove everything after .L
overall_sil_score$n <- sub("n", "",overall_sil_score$n)
  #remove the n prefix
# get Louvain_res parameter value
overall_sil_score$L <- sub("[[:print:]]*.L", "",overall_sil_score$parameters)
  #keep everything after .L


ggplot(overall_sil_score, aes(x=L, y=n, fill=overall_silhouette_score)) +
  geom_tile(color="gray",lwd=0.05, linetype=1) +
  geom_text(aes(label = round(x=overall_silhouette_score,digits=3)), color = "black", size = 3) +
  scale_fill_gradient(low="white",high="blue") +
  xlab("Louvain resolution") + ylab("nearest neighbors") + labs(fill = "Silhouette score")
ggsave(file="~/Downloads/RNAseq/cluster_by_genes/DESC_parameter_test/heatmap_silScore.png")
ggsave(file="~/Downloads/RNAseq/cluster_by_genes/DESC_parameter_test/heatmap_silScore.svg")




