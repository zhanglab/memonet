#this script makes a table of all the genes and their DESeq2 results from DESC
#you can filter to padj<0.1 (or any value) within the script, so you don't end up with a huge file
#run this in the DESeq2 output folder that you get the _all_genes.csv or _sig_genes.csv files from

library(tidyverse)
library(schoolmath) #for is.positive

#usage: module load r/4.2.0  **don't do this in interactive or it won't work
#       nohup srun Rscript /work/pi_yingzhang_uri_edu/kdunton/memonet_github_repo/memonet/scripts/DESC-DESeq2-table.r &


filenames <- Sys.glob("*_all_genes.csv")
  #this just looks at files in whatever directory you run it in
  #Sys.glob is a base R function that finds all files with a specified pattern


#initiate empty df to append to in loop; make sure the colnames match what's generated in loop
new_df <- data.frame(matrix(ncol=8,nrow=0))
x <- c("gene","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj","cluster")
colnames(new_df) <- x


######## loop through each file ########
for(file in filenames) {
  #read in the DESeq result file
  df <- read.csv(file)
  filename <- basename(file)
  print(filename)
    #this extracts the actual filename from the path, so you can extract the cluster ID
  clusterID <- sub("_[[:print:]]*", "",filename)
    #take everything after the first _ of the filename and replace with a blank
  
  subset_df <- df
  if(dim(subset_df)[1] == 0){
    next
    #this says if the row dimension (indicated with [1]) is 0, skip this iteration of the loop
  }
  
  #then add column for cluster
  subset_df$cluster <- clusterID
  
  #assign a name for the subset_df based on clusterID, so you can add all of them together outside the loop
  #assign(paste0("df_",clusterID), subset_df)
  new_df <- do.call("rbind", list(new_df,subset_df))
}



######### format df ########
## make column for upregulated or downregulated
new_df <- new_df[!is.na(new_df$log2FoldChange), ]
  #first remove NA values from log2FoldChange column
new_df$boolean <- is.positive(new_df$log2FoldChange)
new_df$boolean <- as.character(new_df$boolean)
  #need to convert from logical to character type for recode to work
new_df <- new_df %>% mutate(result = recode(boolean,
                                            'TRUE'='upregulated',
                                            'FALSE'='downregulated'))
#remove boolean column now that it's done its purpose
new_df$boolean <- NULL




###### save to file #####
## save all genes regardless of significance
write.csv(new_df, "DEGstats_allGenes.csv", row.names=FALSE)
## save only sig genes
new_df <- subset(new_df, padj < 0.05)
write.csv(new_df, "DEGstats_padj0.05.csv", row.names=FALSE)






