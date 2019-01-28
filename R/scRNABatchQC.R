scRNABatchQC<-function(inputfiles,samplenames=NULL, organism=c("hsapiens","mmusculus"),outputFile="report.html", sampleRatio=1,cache=FALSE, nHVGs=1000,nPC=10,sf=10000,logFC=1,FDR=0.01){
  organism<-match.arg(organism)
  if(!missing(cache) & cache){
    plotData<-prepareReportData(inputfiles=inputfiles,samplenames=samplenames, organism=organism,  sampleRatio=sampleRatio,nHVGs=nHVGs, nPC=nPC, sf=sf, logFC=logFC, FDR=FDR)
  }else{
    plotData<-prepareReportData(inputfiles=inputfiles,samplenames=samplenames, organism=organism,  sampleRatio=sampleRatio,nHVGs=nHVGs, nPC=nPC, sf=sf, logFC=logFC, FDR=FDR)
  }
  
  #reportRmd <- "E:/sqh/programs/scRNABatchQC/inst/report/scRNABatchQCreport.Rmd"
  reportRmd <- system.file("report/scRNABatchQCreport.Rmd", package="scRNABatchQC")
  
  outputFile <- getAbsolutePath(outputFile)
  output_dir = dirname(outputFile)
  output_file = basename(outputFile)
  
  cat("Output report to:", outputFile, "\n")
  rmarkdown::render(reportRmd,
                    output_dir = output_dir,
                    output_file = output_file,
                    params = list(data = plotData))
}
