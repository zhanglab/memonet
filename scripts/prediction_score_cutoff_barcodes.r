#this script makes separate barcode files for different prediction score cutoff values

data <- read.csv("~/Downloads/RNAseq/AIBSmapping/OA/prediction_scores.csv")
data <- subset(data, AIBS_predicted_label %in% c("L2/3 IT_1","L2/3 IT_2","L2/3 IT_3"))
#add column that combines the prediction scores for IT_1,_2,_3; we want to look at the overall confidence for L2/3
data$L23_score <- data$prediction.score.L2.3.IT_1 + data$prediction.score.L2.3.IT_2 + data$prediction.score.L2.3.IT_3


## loop through different cutoff values and subset the barcodes
for(x in c(0.25,0.30,0.35)){
  cutoff <- x
  data1 <- subset(data, L23_score > cutoff)
  data1 <- data1[,1]
  write.csv(data1, paste0("L23barcodes-fromAIBS_",cutoff,".csv"), quote=FALSE, row.names=FALSE)
}
