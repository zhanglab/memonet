# this script calculates md5 checksums for the fastq files, before uploading to GEO ncbi repository
library(tools)

filenames <- Sys.glob("/data/zhanglab/jingwang/brain/RNAseq/Takaki/Deep/fastqs/JB*")

checksums <- data.frame(matrix(ncol = 2, nrow = length(filenames)))
colnames(checksums) <- c('filename','checksum')
i <- 1
for(file in filenames){
  filename <- basename(file)
  checksums[i,1] <- filename
  checksum <- md5sum(file)
  checksums[i,2] <- checksum
  i <- i + 1
}
write.csv(checksums, "checksums.csv", row.names=FALSE)


