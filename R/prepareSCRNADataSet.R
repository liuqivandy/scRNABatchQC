##' prepareSCRNADataSet
##'
##' The function read multiple count table and prepare statistics information for each table
##'
##' @param sampleTable the sample table with first column as sample name, second column as file location and an optional third column as transform
##' @return a named list of makeSCRNAdata result
##' @export prepareSCRNADataSet
##' @examples 
##' #sampleTable<-data.frame(Sample=c("S1", "S2", "S3"),
##' #                        File=file.path(("Z:/shengq1/20180214_scRNABatchQC/", c("S1.csv", "S2.csv", "S3.csv")),
##' #                        Transform=c(1,1,1))
##' #sces<-prepareSCRNADataSet(sampleTable)
prepareSCRNADataSet <- function(sampleTable){
  result <- list()
  
  n <- 1
  for (n in 1:nrow(sampleTable)) {
    sampleName<-as.character(sampleTable[n,1])
    countFile<-as.character(sampleTable[n,2])
    cat(sampleName, "\n")
    counts<-as.matrix(read.csv(countFile, row.names=1, header=T))
    result[[n]] <- prepareSCRNAData(counts)
  }
  names(result) <- sampleTable[, 1]
  return(result)
}

DEBUG=TRUE
if(DEBUG){
  source("prepareSCRNAData.R")
  sampleTable<-data.frame(Sample=c("S1", "S2", "S3"),
                          File=file.path(paste0("Z:/shengq1/20180214_scRNABatchQC/", c("S1.csv", "S2.csv", "S3.csv"))),
                          stringsAsFactors = FALSE)
  sces<-prepareSCRNADataSet(sampleTable)
}
