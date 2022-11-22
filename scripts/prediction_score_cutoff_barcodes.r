#this script makes barcode files for various prediction score cutoff values

library(tidyverse)


data <- read.csv("~/Downloads/RNAseq/AIBSmapping/OA/prediction_scores.csv")
data <- subset(data, AIBS_predicted_label %in% c("L2/3 IT_1","L2/3 IT_2","L2/3 IT_3"))
# calculate overall L2/3 prediction score: sum of scores for IT_1,_2,_3; we want to look at the overall confidence for L2/3
data$L23_score <- data$prediction.score.L2.3.IT_1 + data$prediction.score.L2.3.IT_2 + data$prediction.score.L2.3.IT_3


## loop through multiple cutoff values and subset the barcodes
cutoffs <- cutoffs <- seq(0,0.5,by=0.01)
for(x in cutoffs){
  cutoff <- x
  data1 <- subset(data, L23_score > cutoff)
  data1 <- data1[,1]
  write.csv(data1, paste0("~/Downloads/RNAseq/AIBSmapping/OA/barcode_files/L23barcodes-fromAIBS_",cutoff,".csv"), quote=FALSE, row.names=FALSE)
}
