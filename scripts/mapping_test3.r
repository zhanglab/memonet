#this script visualizes the results of the AIBS mapping test; calculates the mean and median max.prediction.score for AIBS test mapping and for MEMONET data mapped to AIBS


library(ggplot2)
library(tidyverse)
library(caret)
library(svglite)



######## read in test scores ########

### read in the file containing all iteration runs ###
data <- read.csv("~/Downloads/RNAseq/AIBSmapping/test/prediction_cutoff.csv")

### calculate mis-classification percentage across the 100 runs ###
# add column for match
data$match <- data$original_label == data$AIBS_predicted_label
match_tbl <- as.data.frame(table(data$match))
false <- match_tbl[1,2]
total <- length(data$X)
print("mis-classification percentage:")
(false/total)*100


### get average and median max.prediction.score for AIBS whole dataset and for L2/3 cells ###
summary_scores <- data.frame(matrix(ncol = 4, nrow = 4))
colnames(summary_scores) <- c('dataset','mean prediction.score.max','median prediction.score.max','st.error')
summary_scores$dataset <- c('AIBS all cells','AIBS L2/3 cells','MEMONET all cells','MEMONET L2/3 cells')
# calculate for all cells in AIBS dataset
avg <- mean(data$prediction.score.max)
median <- median(data$prediction.score.max)
summary_scores$`mean prediction.score.max`[summary_scores$dataset == 'AIBS all cells' ] <- avg
summary_scores$`median prediction.score.max`[summary_scores$dataset == 'AIBS all cells' ] <- median
# calculate for L2/3 cells in AIBS dataset
data_L23 <- subset(data, subclass_label %in% "L2/3 IT")
avg <- mean(data_L23$prediction.score.max)
median <- median(data_L23$prediction.score.max)
summary_scores$`mean prediction.score.max`[summary_scores$dataset == 'AIBS L2/3 cells' ] <- avg
summary_scores$`median prediction.score.max`[summary_scores$dataset == 'AIBS L2/3 cells' ] <- median

### calculate standard error: BOOTSTRAP - sample with replacement 
# calculate for all cells
set.seed(13579)   # set a seed for consistency/reproducibility
n.obs <- length(data$X)  # the number of observations to sample in each bootstrap (the same as the number of values you have)
B <- 1000  # the number of bootstrap samples
# bootstrap the max prediction scores
boot_df <- matrix(sample(data$prediction.score.max, size= B*n.obs, 
                         replace=TRUE), ncol=B, nrow=n.obs)
# calculate median of each bootstrap: gives 1000 medians
boot_medians <- colMedians(boot_df)
# calculate sd of the bootstrap medians: this will be reported as the standard error
boot_se <- sd(boot_medians)
summary_scores$`st.error`[summary_scores$dataset == 'AIBS all cells' ] <- boot_se

# calculate for L2/3 cells
set.seed(13579)   # set a seed for consistency/reproducibility
n.obs <- length(data_L23$X)  # the number of observations to sample in each bootstrap (the same as the number of values you have)
B <- 1000  # the number of bootstrap samples
# bootstrap the max prediction scores
boot_df <- matrix(sample(data_L23$prediction.score.max, size= B*n.obs, 
                         replace=TRUE), ncol=B, nrow=n.obs)
# calculate median of each bootstrap: gives 1000 medians
boot_medians <- colMedians(boot_df)
# calculate stdev of the bootstrap medians: this will be reported as the standard error
boot_se <- sd(boot_medians)
summary_scores$`st.error`[summary_scores$dataset == 'AIBS L2/3 cells' ] <- boot_se

#** continue adding to summary_scores in the next section that reads in OA mapping results (MEMONET mapped to AIBS)



### plot confusion matrix showing original labels to predicted labels ###
data$AIBS_predicted_label <- as.factor(data$AIBS_predicted_label)
data$original_label <- as.factor(data$original_label)
confusion <- confusionMatrix(data=data$AIBS_predicted_label, reference=data$original_label)
  #data is the predicted label; reference is the true label aka the original label in the AIBS dataset
#view matrix 
confusion_df <- as.data.frame(confusion$table)
confusion_df <- confusion_df %>% rename(Original_label=Reference)

#normalize by column (each column of confusion matrix should add up to 1)
wide_confusion <- as.data.frame(pivot_wider(confusion_df,names_from = Original_label,values_from = Freq))
#pivot so that it resembles a confusion matrix visualization, so that you can scale the values by column
#set the Predictions to rownames so all the values are numbers (in order for scale() to work)
row.names(wide_confusion) <- wide_confusion$Prediction
wide_confusion$Prediction <- NULL
#get the sum of each column
c1 <- colSums(wide_confusion)
scaled_confusion <- as.data.frame(scale(wide_confusion, center=FALSE, scale=c1))
#pass the colSums variable into scale so that each value is divided by the sum
#change back to long format for ggplot
scaled_confusion <- rownames_to_column(scaled_confusion,var="Prediction")
#first have to change rownames (predictions) back into a column
scaled_confusion <- pivot_longer(scaled_confusion, cols=!Prediction, names_to ="Original_label", values_to="Freq")
#use all columns except Prediction to put into a column, which is set to "Original_label"
#the values get moved to one column called Freq
scaled_confusion[is.na(scaled_confusion)] = 0
  #set any NA values to 0

#visualize
ggplot(scaled_confusion, aes(x=Original_label, y=Prediction, fill=Freq)) + 
  geom_tile(color="gray",lwd=0.05, linetype=1) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  xlab('Original cell type') + ylab('Predicted cell type') + labs(fill = "Frequency") +
  scale_fill_gradient(low="white",high="blue") +
  #theme(axis.text = element_text(size = 12)) +
  theme(axis.title = element_text(size = 16)) +
  theme(plot.title = element_text(size = 20)) +
  theme(legend.text = element_text(size = 14)) +
  theme(legend.title = element_text(size = 16)) 
ggsave("~/Downloads/RNAseq/AIBSmapping/test/AIBStest_confusionMtx.png", width=10,height=10)
ggsave("~/Downloads/RNAseq/AIBSmapping/test/AIBStest_confusionMtx.svg", width=10,height=10)




#----------------
### read in prediction scores from OA mapping - Our (MEMONET) data onto AIBS ###
data2 <- read.csv("~/Downloads/RNAseq/AIBSmapping/OA/prediction_scores.csv")

### get average and median max.prediction.score for whole MEMONET dataset and for L2/3 cells ###
# calculate for all cells
avg <- mean(data2$prediction.score.max)
median <- median(data2$prediction.score.max)
summary_scores$`mean prediction.score.max`[summary_scores$dataset == 'MEMONET all cells' ] <- avg
summary_scores$`median prediction.score.max`[summary_scores$dataset == 'MEMONET all cells' ] <- median
# calculate for L2/3 cells
data2_L23 <- subset(data2, AIBS_predicted_label %in% c("L2/3 IT_1","L2/3 IT_2","L2/3 IT_3"))
avg <- mean(data2_L23$prediction.score.max)
median <- median(data2_L23$prediction.score.max)
summary_scores$`mean prediction.score.max`[summary_scores$dataset == 'MEMONET L2/3 cells' ] <- avg
summary_scores$`median prediction.score.max`[summary_scores$dataset == 'MEMONET L2/3 cells' ] <- median


### calculate standard error: BOOTSTRAP - sample with replacement 
# calculate for all cells
set.seed(13579)   # set a seed for consistency/reproducibility
n.obs <- length(data2$X)  # the number of observations to sample in each bootstrap (the same as the number of values you have)
B <- 1000  # the number of bootstrap samples
# bootstrap the max prediction scores
boot_df <- matrix(sample(data2$prediction.score.max, size= B*n.obs, 
                         replace=TRUE), ncol=B, nrow=n.obs)
# calculate median of each bootstrap: gives 1000 medians
boot_medians <- colMedians(boot_df)
# calculate stdev of the bootstrap medians: this will be reported as the standard error
boot_se <- sd(boot_medians)
summary_scores$`st.error`[summary_scores$dataset == 'MEMONET all cells' ] <- boot_se

# calculate for L2/3 cells
set.seed(13579)   # set a seed for consistency/reproducibility
n.obs <- length(data2_L23$X)  # the number of observations to sample in each bootstrap (the same as the number of values you have)
B <- 1000  # the number of bootstrap samples
# bootstrap the max prediction scores
boot_df <- matrix(sample(data2_L23$prediction.score.max, size= B*n.obs, 
                         replace=TRUE), ncol=B, nrow=n.obs)
# calculate median of each bootstrap: gives 1000 medians
boot_medians <- colMedians(boot_df)
# calculate stdev of the bootstrap medians: this will be reported as the standard error
boot_se <- sd(boot_medians)
summary_scores$`st.error`[summary_scores$dataset == 'MEMONET L2/3 cells' ] <- boot_se


### save summary_scores to file
write.csv(summary_scores, "~/Downloads/RNAseq/AIBSmapping/test/maxPredictionScores-AIBStest_and_OA.csv", row.names=FALSE)




