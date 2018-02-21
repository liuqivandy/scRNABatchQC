source("makeSCRNAdata.R")

path <- "D:/Work/scRNABatchQC" 
fileNames <- c("count1.csv", "count2.csv", "count3.csv")
sampleNames <- c("S1", "S2", "S3")

dat <- list()
  
for (n in 1:length(fileNames)) {
  counttable <- as.matrix(read.csv(file.path(path, fileNames[n]), row.names = 1))

  dat[[n]] <- makeSCRNAdata(counttable)
}
names(dat) <- sampleNames

