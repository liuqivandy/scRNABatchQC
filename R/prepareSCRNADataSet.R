

##' prepareSCRNADataSet
##'
##' The function read multiple count table and prepare statistics information for each table
##'
##' @param sampleTable the sample table with first column as sample name, second column as file location and an optional third column as transform
##' @return a named list of makeSCRNAdata result
##' @export prepareSCRNADataSet
##' @examples 
##' #sampleTable<-data.frame(Sample=c("S1", "S2", "S3"),
##' #                        File=paste0("Z:/shengq1/20180214_scRNABatchQC/", c("qi_m1.csv", "qi_m2.csv", "s1_pan_qi_1.csv")),
##' #                        Transform=c(1,1,1))
##' #sces<-prepareSCRNADataSet(sampleTable)
prepareSCRNADataSet<-function(sampleTable){
  result<-list()
  
  n<-1
  for (n in 1:nrow(sampleTable)) {
    sampleName <- as.character(sampleTable[n,1])
    sampleFile <- as.character(sampleTable[n,2])
    sampleTransform <- ifelse(ncol(sampleTable) > 2, sampleTable[n,3], 0)
    cat(sampleName, "\n")
    counttable <- as.matrix(read.csv(sampleFile, row.names = 1))
    if(sampleTransform){
      counttable <- t(counttable)
    }
    result[[n]] <- prepareSCRNAData(counttable)
  }
  names(result) <- sampleTable[,1]
  return(result)
}

DEBUG=FALSE
if(DEBUG){
  source("prepareSCRNAData.R")
  sampleTable<-data.frame(Sample=c("S1", "S2", "S3"),
                          File=paste0("Z:/shengq1/20180214_scRNABatchQC/", c("qi_m1.csv", "qi_m2.csv", "s1_pan_qi_1.csv")),
                          Transform=c(1,1,1))
  sces<-prepareSCRNADataSet(sampleTable)
}