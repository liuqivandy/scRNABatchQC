source("makeSCRNAdata.R")

readInCountTables <- function(path, fileNames, sampleNames) {
  #path = YOUR_FILE_PATH, 
  #filenames = c("count1.csv","count2.csv", "count3.csv"),
  #sampleNames = c("S1","S2", "S3")
  
  stopifnot(length(fileNames) == length(sampleNames))
  
  dat <- list()
  
  for (n in 1:length(filenames)) {
    dat[[n]] <- makeSCRNAdata(fread(file.path(path, fileNames[n])))
  }
  names(dat) <- sampleNames

  return(dat)
}
