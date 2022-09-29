#this script concatenates a table of all iterations of testA-prediction_cutoff.r


library(tidyverse)
library(schoolmath) #for is.positive


filenames <- Sys.glob("~/Downloads/RNAseq/AIBSmapping/test/prediction_scores_*.csv")
  #Sys.glob is a base R function that finds all files with a specified pattern

#initiate empty list to append df to in loop
all_results <- list()

######## loop through each file ########
for(file in filenames) {
  df <- read.csv(file)
    #read in each prediction score iteration file
  filename <- basename(file)
    #this extracts the actual filename from the path, so you can extract the cluster ID
  print(filename)
  iteration <- sub("prediction_scores_", "",filename)
    #replace "prediction_scores_" with a blank
  iteration <- sub(".csv","",iteration)
    #replace ".csv" with a blank
    #now you have just the number of the array job aka the iteration number
  df$iteration <- iteration
  all_results <- do.call("rbind", list(all_results,df))
}


###### save to file #####
write.csv(all_results, "prediction_cutoff.csv", row.names=FALSE)

print("done!")




